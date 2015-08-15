#!/usr/bin/env python3

# Tool to extract method signatures from a header file and do stuff with
# that.

import sys
import os
import argparse
import re
import collections
import json

function_pattern = re.compile(r"""
    \b
    (?P<PREFIX> FLINT_DLL | ([A-Z0-9_]+)_INLINE)
    \s+
    (?P<FRTYPE> [\w\s*]+ )
    \b
    (?P<FNAME> \w+ )
    \s*
    \(
    (?P<ARGS> [\w\s*,]*)
    \)
""", re.VERBOSE | re.DOTALL)

argument_pattern = re.compile(r"""
    \s*
    (?P<ATYPE> [\w\s*]+ )
    \b
    (?P<ANAME> \w+ )
    \s*
""", re.VERBOSE | re.DOTALL)

type_component_pattern = re.compile(r"""
    \s*
    (?P<TCOMP>
        (\w+\b)
        |
        ([*])
    )
    \s*
""", re.VERBOSE | re.DOTALL)

def parse_type(st):
    """Extracts type components from the given string.
    The type is returned as a single space-separated string.
    """
    global type_component_pattern
    res = []
    pos = 0
    while pos < len(st):
        m = type_component_pattern.match(st, pos)
        if m is None:
            if st[pos:].isspace():
                break
            else:
                raise ValueError("invalid type component")
        res.append(m.group('TCOMP'))
        pos = m.end()
    return ' '.join(res)

def parse_args(st):
    """Extracts a list of arguments from the given string."""
    global argument_pattern
    res = []
    for argst in st.split(','):
        m = argument_pattern.match(argst)
        if m is None or m.end() != len(argst):
            raise ValueError("invalid argument string")
        atype = parse_type(m.group('ATYPE'))
        aname = m.group('ANAME')
        res.append(Arg(atype, aname))
    return res

def parse_func(st, cname=""):
    """Extracts a function from the given string."""
    global function_pattern
    m = function_pattern.match(st)
    if m is None or m.end() != len(st):
        raise ValueError("invalid function string")
    rtype = parse_type(m.group('FRTYPE'))
    args = parse_args(m.group('ARGS'))
    return Func(cname, m.group('FNAME'), rtype, args, m.group('PREFIX'))

def type_re(typestr):
    res = ''
    for comp in typestr.split():
        if comp == '*':
            res += r'[*]\s*'
        else:
            res += r'\b' + comp + r'\b\s*'
    return res

def func_re(fun):
    """Returns a regular expression matching the given function."""
    rawpat = (
        type_re(fun.rtype)
        + r'\b' + fun.name
        + r'\s*\(\s*'
        + r',\s*'.join(type_re(a.type) + r'\b' + a.name + r'\s*'
                       for a in fun.args)
        + r'\)\s*\{'
        )
    # print("RAWPAT for", str(fun))
    # print(rawpat)
    return re.compile(rawpat, re.DOTALL)

Arg = collections.namedtuple('Arg', ['type', 'name'])

class Func:
    def __init__(self, cname, name, rtype, args, prefix=""):
        self.cname = cname
        self.name = name
        self.rtype = rtype
        self.args = args
        self.prefix = prefix

    def find_code(self, s):
        """Searches for the code for this function in the given string."""
        return func_re(self).search(s)

    def is_hidden(self):
        return self.name.startswith('_')

    def short_name(self):
        ret = self.name.lstrip('_')
        if ret.startswith(self.cname):
            ret = ret[len(self.cname):]
        return ret.lstrip('_')

    def is_inline(self):
        if self.prefix.endswith('_INLINE'):
            return True
        elif self.prefix == 'FLINT_DLL':
            return False
        else:
            raise ValueError("invalid prefix: {}".format(self.prefix))

    def __str__(self):
        return (
            self.prefix + ' ' + self.rtype + ' ' + self.name + '('
            + ', '.join(t + ' ' + n for t,n in self.args) + ')'
            )

    def __eq__(self, other):
        return type(other) is Func and self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not (self == other)

def get_funcs(fn, cname):
    """Extracts all function signatures from the given file.

    Returns a list of all function names (in declaration order),
    and a dictionary mapping the names to their signature.
    """
    global function_pattern
    with open(fn) as fin:
        orig = fin.read()
    fnames, fsigs = [], {}
    for m in function_pattern.finditer(orig):
        f = parse_func(m.group(0), cname)
        fnames.append(f.short_name())
        fsigs[f.short_name()] = f
    return fnames, fsigs

def get_code(cname):
    """Goes through all .c files in specified dir and finds what functions
    are implemented."""
    res = []
    if not os.path.isdir(cname):
        print("ERROR: no directory", cname)
    else:
        for filename in os.listdir(cname):
            if filename.endswith('.c'):
                cname = filename[:-2]
                res.append(cname)
    return res

def get_tests(cname):
    """Finds test files in specified dir."""
    res = []
    tdir = os.path.join(cname, 'test')
    if not os.path.isdir(cname):
        print("ERROR: no directory", cname)
    else:
        for filename in os.listdir(tdir):
            if filename.startswith('t-') and filename.endswith('.c'):
                cname = filename[2:-2]
                res.append(cname)
    return res

def get_cname(fn):
    fbn = os.path.basename(fn)
    if '.' in fbn:
        return fbn[:fbn.find('.')]
    else:
        return fbn

def set_default_db(dbdir):
    global db_file
    db_file = os.path.join(dbdir, 'fun_db.json')

def open_db():
    global db_file, db
    db = {}
    try:
        with open(db_file) as dbin:
            db = json.load(dbin)
    except IOError:
        print("WARNING: {} not found; loading empty database".format(db_file))
    except OSError:
        print("WARNING: {} unreadable; loading empty database".format(db_file))
    except ValueError:
        print("WARNING: {} badly formatted; loading empty database".format(db_file))
    for cname in db:
        sigs = db[cname]['signatures']
        for fname in sigs:
            orig = sigs[fname]
            sigs[fname] = Func(
                cname,
                orig['name'],
                orig['rtype'],
                [Arg(t,n) for (t,n) in orig['args']],
                orig['prefix']
            )

def json_default(obj):
    if type(obj) is Func:
        res = dict(obj.__dict__)
        del res['cname']
        return res
    else:
        return obj

def save_db():
    global db_file, db, backup_db
    if backup_db and os.path.exists(db_file):
        os.rename(db_file, '.'.join(db_file, backup_db))
    with open(db_file,'w') as dbout:
        json.dump(db, dbout, indent=4, sort_keys=True, default=json_default)

def list_diffs(title, stored, found):
    stset = set(stored)
    fset = set(found)
    if stset == fset:
        if stored != found:
            print(title)
            print("INCONSISTENT: same entries, different order")
            print()
    else:
        print(title)
        for x in stored:
            if x not in fset:
                print("NOT FOUND:", x)
        for x in found:
            if x not in stset:
                print("MISSING:", x)
        print()

def dict_diffs(title, stored, found):
    if stored != found:
        print(title)
        for x in stored:
            if x not in found:
                print("NOT FOUND:", x)
            elif stored[x] != found[x]:
                print("INCONSISTENT:", x)
        for x in found:
            if x not in stored:
                print("MISSING:", x)
        print()


if __name__ == '__main__':
    commands = {'check', 'update', 'code_pop', 'code_gen', 'tests_pop', 'tests_gen'}

    set_default_db(os.path.dirname(os.path.abspath(sys.argv[0])))

    parser = argparse.ArgumentParser(
        description="Extracts method signatures from a header file and does some useful processing.",
        )
    
    parser.add_argument('command',
        help="The action to take ({}, default check)".format('|'.join(commands)),
        choices=commands, nargs='?', default='check'
        )

    parser.add_argument('header_file', 
        help="The header file to get signatures from (default: all files)",
        nargs='?',
        )

    parser.add_argument('--database', '-d',
        help="Database file to use (default {})".format(db_file),
        )

    parser.add_argument('--backup', '-b',
        help="Extension for database backup; (default: 'bak')",
        )
    
    parser.add_argument('--no-backup',
        help="Don't save any database backup",
        action='store_const', dest='backup', const=None,
        )

    parser.add_argument('--quick', '-q',
        help="Don't ask any questions",
        action='store_const', const=True, # FIXME better way?
        )

    args = parser.parse_args()
    command = args.command
    one_file = args.header_file is not None
    if args.database is not None:
        db_file = args.database
    backup_db = args.backup

    ask = not args.quick

    open_db()
    if one_file:
        os.chdir(os.path.dirname(os.path.abspath(args.header_file)))
        headers = [os.path.basename(args.header_file)]
    else:
        os.chdir(os.path.dirname(os.path.dirname(db_file)))
        headers = list(db)

    for hf in headers:
        if not hf.endswith('.h'):
            hf = hf + '.h'
        cname = get_cname(hf)
        print('Running', command, 'on', cname, '...')
        if command == 'check':
            if cname not in db:
                print("MISSING class name", cname)
                print()
                continue
            dbc = db[cname]

            fnames, fsigs = get_funcs(hf, cname)
            list_diffs("Header function names", dbc['functions'], fnames)
            dict_diffs("Header function signatures", dbc['signatures'], fsigs)

            codefs = get_code(cname)
            list_diffs("Code files", sorted(dbc['code']), sorted(codefs))
            for codef, flist in dbc['code'].items():
                if codef == 'inlines':
                    continue
                cfn = os.path.join(cname, codef+'.c')
                if os.path.exists(cfn):
                    with open(cfn) as cfin:
                        cf_all = cfin.read()
                    for fname in flist:
                        f = dbc['signatures'][fname]
                        if not f.find_code(cf_all):
                            print("INCONSISTENT: Code for {} not found in {}"
                                .format(fname, codef))
                else:
                    print("CODE FILE NOT FOUND:", codef)

            ft = get_tests(cname)
            list_diffs("Test names", sorted(dbc['tests']), sorted(ft))
            tested = set()
            for tname, flist in dbc['tests'].items():
                if not os.path.exists(
                        os.path.join(cname, 'test', 't-'+tname+'.c')):
                    print("TEST FILE NOT FOUND:", tname)
                tested.update(flist)
            list_diffs("Tested functions", sorted(dbc['functions']), sorted(tested))

        elif command == 'update':
            if cname not in db:
                db[cname] = {'code':{}, 'tests':{}}
            dbc = db[cname]
            
            fnames, fsigs = get_funcs(hf, cname)
            dbc['functions'] = fnames
            dbc['signatures'] = fsigs

            codefs = get_code(cname)
            oldcode = dbc['code']
            dbc['code'] = {}
            for codef in codefs:
                cfn = os.path.join(cname, codef+'.c')
                flist = []
                with open(cfn) as cfin:
                    cf_all = cfin.read()
                if codef in oldcode:
                    for fname in oldcode[codef]:
                        if fname in fsigs:
                            f = fsigs[fname]
                            if not f.is_inline() and f.find_code(cf_all):
                                flist.append(fname)
                elif codef in fsigs:
                    if fsigs[codef].find_code(cf_all):
                        flist.append(codef)
                if codef == 'inlines':
                    for fname, f in fsigs.items():
                        if f.is_inline():
                            flist.append(fname)
                elif not flist:
                    if ask:
                        while True:
                            funcs = set( 
                                input("What functions are coded in {}? " 
                                      .format(codef))
                                .split())
                            if all(f in fsigs for f in funcs):
                                flist = sorted(funcs)
                                break
                            print("ERROR: invalid function name(s) given.")
                    else:
                        print("WARNING: empty flist for", codef)
                dbc['code'][codef] = flist

            for tname in get_tests(cname):
                if tname not in dbc['tests']:
                    if tname in fsigs:
                        dbc['tests'][tname] = [tname]
                    else:
                        if ask:
                            while True:
                                funcs = set( 
                                    input("What functions are tested in {}? " 
                                          .format(tname))
                                    .split())
                                if all(f in fsigs for f in funcs):
                                    break
                                print("ERROR: invalid function name(s) given.")
                        else:
                            print("WARNING: don't know what's tested in", tname)
                            funcs = []
                        dbc['tests'][tname] = sorted(funcs)

        elif command == 'code_pop':
            for fname, f in dbc['signatures'].items():
                if fname not in dbc['code']:
                    if f.is_inline():
                        dbc['code'][fname] = 'inlines'
                    else:
                        dbc['code'][fname] = fname
                
        elif command == 'code_gen':
            raise NotImplementedError(command)
            pass #TODO
        elif command == 'tests_pop':
            raise NotImplementedError(command)
            pass #TODO
        elif command == 'tests_gen':
            raise NotImplementedError(command)
            pass #TODO
        else:
            raise ValueError("Command {} not supported".format(command))
        print()

    save_db()
