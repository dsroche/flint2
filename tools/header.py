#!/usr/bin/env python3

# Tool to extract method signatures from a header file and do stuff with
# that.

import sys
import os
import argparse
import re
import collections
import json
import copy
import itertools

#### Globally defined regular expressions ####

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

#### Class definitions ####

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
        """Returns a UNIQUE short version.
        e.g. _fmpz_poly_add becomes add_"""
        ret = self.name
        if self.is_hidden():
            ret = ret.lstrip('_') + '_';
        if ret.startswith(self.cname):
            ret = ret[len(self.cname):]
        return ret.lstrip('_')
    
    def shorter_name(self):
        """Returns a NON-UNIQUE short version.
        e.g. _fmpz_poly_add becomes add"""
        return self.short_name().rstrip('_')

    def is_inline(self):
        if self.prefix.endswith('_INLINE'):
            return True
        elif self.prefix == 'FLINT_DLL':
            return False
        else:
            raise ValueError("invalid prefix: {}".format(self.prefix))

    def write_impl(self, outfile):
        print(self.rtype.strip(), ' ', self.name.strip(), '(', 
              file=outfile, sep='', end='')
        print(', '.join(t.strip()+' '+n.strip() for (t,n) in self.args), 
              file=outfile, end='')
        print(")\n{\n    /* TODO */", file=outfile);
        if self.rtype != 'void':
            print("    return 0; /* FIXME dummy return value */", file=outfile)
        print("}", file=outfile)

    def __str__(self):
        return (
            self.prefix + ' ' + self.rtype + ' ' + self.name + '('
            + ', '.join(t + ' ' + n for t,n in self.args) + ')'
            )

    def __eq__(self, other):
        return type(other) is Func and self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not (self == other)

#### PARSING FUNCTIONS ####

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
    """Generates a regular expression to search for the given type name."""
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
    return re.compile(rawpat, re.DOTALL)

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

#### FILE SEARCHING FUNCTIONS ####

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
    if not os.path.isdir(tdir):
        print("ERROR: no directory", tdir)
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

#### OTHER FILE STUFF ####

banner = None
def get_banner():
    """Loads the banner text from banner.txt."""
    global db_dir, banner
    if banner is None:
        banfile = os.path.join(db_dir, 'banner.txt')
        with open(banfile) as banin:
            banner = banin.read()
        authfile = os.path.join(db_dir, 'author.txt')
        if os.path.exists(authfile):
            with open(authfile) as authin:
                banner += authin.read()
        else:
            print("WARNING: no author file", authfile)
    return banner

def code_file(cname, short_name):
    return os.path.join(cname, short_name + '.c')

#### PERSISTENT DATABASE MANIPULATION ####

def set_default_db(dbdir):
    global db_dir, db_file
    db_dir = dbdir
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

#### HELPER FUNCTIONS FOR CHECK ####

def list_diffs(title, stored, found):
    stset = set(stored)
    fset = set(found)
    if stset == fset:
        if stored != found:
            print()
            print(title)
            print("INCONSISTENT: same entries, different order")
            print()
    else:
        print()
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
        print()
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

##### OTHER HELPER FUNCS #####

last_confirm = None

def confirm(text, level=2, remember=False):
    global ask_level, last_confirm
    if ask_level < level or (remember and last_confirm):
        print(text+"...")
        return True
    elif remember and (last_confirm is False):
        print(text + "? N")
        return False
    elif ask_level >= level:
        question = text + "? [y/n"
        if remember:
            question += "/Y/N"
        question += "] "
        while True:
            yn = input(question).strip()
            if not remember:
                yn = yn.lower()
            if yn == 'n':
                return False
            elif yn == 'y':
                return True
            elif yn == 'N':
                last_confirm = False
                return False
            elif yn == 'Y':
                last_confirm = True
                return True

def forget_confirm():
    global last_confirm
    last_confirm = None

def trim_list(lst, trimto):
    return [x for x in trimto if x in lst]

def trim_dict(dic, trimto):
    return {x:dic[x] for x in trimto if x in dic}


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

    parser.add_argument('functions',
        help="Function(s) to perform the operation on (default: all in header).",
        nargs='*',
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

    parser.add_argument('--confirm', '-c',
        help="Confirmation level. 0=no questions, default=1",
        type=int, default=1
        )

    args = parser.parse_args()
    command = args.command
    one_file = args.header_file is not None
    if args.database is not None:
        db_file = args.database
    backup_db = args.backup

    ask_level = args.confirm

    open_db()
    db_orig = copy.deepcopy(db)

    if one_file:
        os.chdir(os.path.dirname(os.path.abspath(args.header_file)))
        headers = [os.path.basename(args.header_file)]
    else:
        os.chdir(os.path.dirname(os.path.dirname(db_file)))
        headers = list(db)

    if args.functions:
        all_funcs = False
        do_funcs = args.functions
    else:
        all_funcs = True

    for hf in headers:
        if not hf.endswith('.h'):
            hf = hf + '.h'
        cname = get_cname(hf)
        print('Running', command, 'on', cname, '...')
        if not all_funcs:
            print("Only considering functions", ", ".join(do_funcs))

        ##### CHECK COMMAND #####
        if command == 'check':
            if cname not in db:
                print("MISSING class name", cname)
                print()
                continue
            dbc = db[cname]

            fnames, fsigs = get_funcs(hf, cname)
            funlist = dbc['functions']
            sigdict = dbc['signatures']
            if not all_funcs:
                fnames = trim_list(fnames, do_funcs)
                fsigs = trim_dict(fsigs, do_funcs)
                funlist = trim_list(funlist, do_funcs)
                sigdict = trim_dict(sigdict, do_funcs)

            list_diffs("Header function names", funlist, fnames)
            dict_diffs("Header function signatures", sigdict, fsigs)

            codefs = get_code(cname)
            stored_code = dbc['code']
            if all_funcs:
                list_diffs("Code files", sorted(stored_code), sorted(codefs))
            else:
                stored_code = {f:L for (f,L) in stored_code.items() 
                               if any(x in L for x in do_funcs)}
                for codef in stored_code:
                    if codef not in codefs:
                        print("Code file missing:", codef)
            coded = set()
            for codef, flist in stored_code.items():
                coded.update(flist)
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
                if not flist:
                    print("WARNING: No coded functions listed for file ", codef)
            if all_funcs:
                list_diffs("Coded functions", sorted(dbc['functions']), sorted(coded))
            else:
                for fun in do_funcs:
                    if fun not in coded:
                        print("FUNCTION NOT CODED:", fun)

            ft = get_tests(cname)
            stored_tests = dbc['tests']
            if all_funcs:
                list_diffs("Test names", sorted(stored_tests), sorted(ft))
            else:
                stored_tests = {f:L for (f,L) in stored_tests.items()
                                if any(x in L for x in do_funcs)}
                for testf in stored_tests:
                    if testf not in ft:
                        print("Test file missing:", testf)
            tested = set()
            for tname, flist in stored_tests.items():
                if not os.path.exists(
                        os.path.join(cname, 'test', 't-'+tname+'.c')):
                    print("TEST FILE NOT FOUND:", tname)
                tested.update(flist)
            if all_funcs:
                list_diffs("Tested functions", sorted(dbc['functions']), sorted(tested))
            else:
                for fun in do_funcs:
                    if fun not in tested:
                        print("FUNCTION NOT TESTED:", fun)
            for tfile, tnames in stored_tests.items():
                if not tnames:
                    print("WARNING: No tests listed for file", tfile)

        ##### UPDATE COMMAND #########
        elif command == 'update':
            if cname not in db:
                print("Creating new database entry for", cname)
                db[cname] = {'code':{}, 'tests':{}}
            dbc = db[cname]
            
            fnames, fsigs = get_funcs(hf, cname)
            oldfuncs, oldsigs = dbc['functions'], dbc['signatures']
            dbc['functions'] = fnames
            dbc['signatures'] = fsigs

            if not all_funcs:
                fnames = trim_list(fnames, do_funcs)
                fsigs = trim_dict(fsigs, do_funcs)
            
            for fun in fnames:
                if fun not in dbc['functions']:
                    if confirm("Add function " + fun, remember=True):
                        dbc['functions'].append(fun)
            for fun, sig in fsigs.items():
                if fun not in dbc['signatures'] or sig != dbc['signatures'][fun]:
                    if confirm("Update signature for " + fun, remember=True):
                        dbc['signatures'][fun] = sig
            forget_confirm()

            if all_funcs:
                for fun in dbc['functions']:
                    if (fun not in fnames 
                        and confirm('Remove function ' + fun, remember=True)
                       ):
                        dbc['functions'].remove(fun)
                for fun in dbc['signatures']:
                    if (fun not in fsigs 
                        and confirm('Remove signature for ' + fun, remember=True)
                       ):
                        del dbc['signatures'][fun]
            forget_confirm()

            codefs = get_code(cname)
            oldcode = dbc['code']
            dbc['code'] = {}
            for codef in codefs:
                hcodef = codef+'_'
                cfn = code_file(cname, codef)
                flist = []
                with open(cfn) as cfin:
                    cf_all = cfin.read()
                if codef in oldcode:
                    for fname in oldcode[codef]:
                        if fname in fsigs:
                            f = fsigs[fname]
                            if not f.is_inline() and f.find_code(cf_all):
                                flist.append(fname)
                elif codef in fsigs or hcodef in fsigs:
                    if codef in fsigs and fsigs[codef].find_code(cf_all):
                        flist.append(codef)
                    if hcodef in fsigs and fsigs[hcodef].find_code(cf_all):
                        flist.append(hcodef)
                if not all_funcs and not any(fun in flist for fun in do_funcs):
                    if codef in oldcode:
                        dbc['code'][codef] = oldcode[codef]
                    continue
                if codef == 'inlines':
                    for fname, f in fsigs.items():
                        if f.is_inline():
                            flist.append(fname)
                elif not flist:
                    if ask_level:
                        while True:
                            funcs = set( 
                                input("What functions are coded in {}? " 
                                      .format(codef))
                                .split())
                            if all(f in dbc['signatures'] for f in funcs):
                                flist = sorted(funcs)
                                break
                            print("ERROR: invalid function name(s) given.")
                    else:
                        print("WARNING: empty flist for", codef)
                if ((codef not in oldcode 
                     or sorted(oldcode[codef]) != sorted(flist))
                    and confirm("Update code file "+codef, remember=True)
                   ):
                    dbc['code'][codef] = flist
                elif codef in oldcode:
                    dbc['code'][codef] = oldcode[codef]
            for codef in oldcode:
                if (codef not in dbc['code'] 
                    and not confirm("Remove code for "+codef,remember=True)
                   ):
                    dbc['code'][codef] = oldcode[codef]
            forget_confirm()

            for tname in get_tests(cname):
                if tname not in dbc['tests']:
                    htname = tname+'_'
                    if tname in fsigs or htname in fsigs:
                        if confirm("Add tested functions to "+tname, remember=True):
                            dbc['tests'][tname] = []
                            if tname in fsigs:
                                dbc['tests'][tname].append(tname)
                            if htname in fsigs:
                                dbc['tests'][tname].append(htname)
                    else:
                        if ask_level:
                            while True:
                                funcs = set( 
                                    input("What functions are tested in {}? " 
                                          .format(tname))
                                    .split())
                                if all(f in dbc['signatures'] for f in funcs):
                                    break
                                print("ERROR: invalid function name(s) given.")
                        else:
                            print("WARNING: don't know what's tested in", tname)
                            funcs = []
                        if funcs or confirm("Add empty test list to "+tname, 
                                            remember=True):
                            dbc['tests'][tname] = sorted(funcs)
            forget_confirm()

        ##### CODE_POP COMMAND #####
        elif command == 'code_pop':
            if cname not in db:
                print("MISSING class name", cname)
                print()
                continue
            dbc = db[cname]

            coded = set(itertools.chain(*dbc['code'].values()))
            for fname, f in dbc['signatures'].items():
                if (fname not in coded 
                    and (all_funcs or fname in do_funcs)
                   ):
                    if f.is_inline():
                        if confirm("Populate inline code for " + fname, remember=True):
                            if 'inlines' not in dbc['code']:
                                dbc['code']['inlines'] = []
                            dbc['code']['inlines'].append(fname)
                    else:
                        if confirm("Populate code for " + fname, remember=True):
                            sname = f.shorter_name()
                            if sname not in dbc['code']:
                                dbc['code'][sname] = []
                            dbc['code'][sname].append(fname)
            forget_confirm()
            
        ##### CODE GEN COMMAND #####
        elif command == 'code_gen':
            if cname not in db:
                print("MISSING class name", cname)
                print()
                continue
            dbc = db[cname]

            for codef, funcs in dbc['code'].items():
                if not any(fun in do_funcs for fun in funcs):
                    continue
                cfn = code_file(cname, codef)
                write_fun = set(funcs)
                if os.path.exists(cfn):
                    with open(cfn) as cfin:
                        cf_all = cfin.read()
                    for fname in funcs:
                        sig = dbc['signatures'][fname]
                        if not sig.is_inline() and sig.find_code(cf_all):
                            write_fun.discard(fname)
                else:
                    if confirm("Create file " + cfn, remember=True):
                        banner = get_banner()
                        with open(cfn,'w') as cfout:
                            cfout.write(banner)
                            print("", file=cfout)
                            print('#include "{}.h"'.format(cname), file=cfout)
                    else:
                        continue
                if (write_fun and 
                    confirm("Write functions " + str(write_fun) + " to " + codef,
                            remember=True)
                   ):
                    with open(cfn,'a') as cfout:
                        for fname in sorted(write_fun):
                            print("", file=cfout)
                            sig = dbc['signatures'][fname]
                            sig.write_impl(cfout)
            forget_confirm()

        elif command == 'tests_pop':
            raise NotImplementedError(command)
            pass #TODO
        elif command == 'tests_gen':
            raise NotImplementedError(command)
            pass #TODO
        else:
            raise ValueError("Command {} not supported".format(command))
        print()

    if db != db_orig:
        if confirm("Save modified database"):
            save_db()
