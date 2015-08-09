#!/usr/bin/env python3

# Tool to extract method signatures from a header file and do stuff with
# that.

import sys
import os
import argparse
import re
import collections
import json

class Func:
    fpat = re.compile(r"""
        \b
        (?P<PREFIX> FLINT_DLL | ([A-Z0-9_]+)_INLINE)
        \s+
        (?P<FRTYPE> [^(){};]+ )
        \s+
        (?P<FNAME> \w+ )
        \s*
        \(
        (?P<ARGS> [^)]*)
        \)
    """, re.VERBOSE | re.DOTALL)

    argpat = re.compile(r"""
        \s*
        (?P<ATYPE> [^(){};,]+ )
        \s+
        (?P<ANAME> \w+ )
        \s*
        (, | $)
    """, re.VERBOSE | re.DOTALL)

    Arg = collections.namedtuple('Arg', ['type', 'name'])

    def get_dict(self):
        res = {
            'prefix': self.prefix,
            'rtype': self.rtype,
            'name': self.name,
            'args': self.args,
            }
        return res

    def __init__(self, s, cname=""):
        if type(s) is dict:
            self.prefix = s['prefix']
            self.rtype = s['rtype']
            self.name = s['name']
            self.args = [self.Arg(t,n) for (t,n) in s['args']]
            return

        self.cname = cname
        m = self.fpat.match(s)
        if m is None or m.end() != len(s):
            raise ValueError("invalid function string")
        self.prefix = m.group('PREFIX')
        self.rtype = re.sub(r'\s+', ' ', m.group('FRTYPE'))
        self.name = m.group('FNAME')
        argstr = m.group('ARGS')
        pos = 0
        self.args = []
        while pos < len(argstr):
            am = self.argpat.match(argstr, pos)
            if am is None:
                break
            self.args.append(
                self.Arg(re.sub(r'\s+',' ',am.group('ATYPE')), am.group('ANAME')))
            pos = am.end()

    def find_code(self, s):
        """Searches for the code for this function in the given string."""
        rawpat = (
            r'\s+'.join(self.rtype.split()).replace('*',r'[*]')
            + r'\s+' + self.name
            + r'\s*\(\s*'
            + r',\s*'.join(r'\s+'.join(a.type.split()).replace('*',r'[*]') 
                           + r'\s+' + a.name + r'\s*'
                           for a in self.args)
            + r'\)\s*\{'
            )
        print("RAWPAT for", str(self))
        print(rawpat)
        codepat = re.compile(rawpat, re.DOTALL)
        return codepat.match(s)

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

def json_default(obj):
    return obj.get_dict()

def get_funcs(fn, cname):
    """Extracts all function signatures from the given file.

    Returns a list of all function names (in declaration order),
    and a dictionary mapping the names to their signature.
    """
    with open(fn) as fin:
        orig = fin.read()
    fnames, fsigs = [], {}
    for m in Func.fpat.finditer(orig):
        f = Func(m.group(0), cname)
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
        cname = get_cname(hf)
        print('Running', command, 'on', cname, '...')
        if command == 'check':
            if cname not in db:
                print("MISSING class name", cname)
                print()
                continue
            dbc = db['cname']

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
