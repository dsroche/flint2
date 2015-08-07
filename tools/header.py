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

    def __init__(self, s, cname=""):
        self.cname = cname
        m = self.fpat.fullmatch(s)
        if m is None:
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
            self.prefix + ' ' + self.frtype + ' ' + self.fname + '('
            + ', '.join(t + ' ' + n for t,n in self.args) + ')'
            )

def get_funcs(fn, cname):
    """Extracts all function signatures from the given file."""
    with open(fn) as fin:
        orig = fin.read()
    funcs = {}
    for m in Func.fpat.finditer(orig):
        funcs.append(Func(m.group(0), cname))
    return funcs

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
    except OSError:
        print("WARNING: {} unreadable; loading empty database".format(db_file))
    except ValueError:
        print("WARNING: {} badly formatted; loading empty database".format(db_file))

def save_db():
    global db_file, db, backup_db
    if backup_db and os.path.exists(db_file):
        os.rename(db_file, '.'.join(db_file, backup_db))
    with open(db_file,'w') as dbout:
        json.dump(db, dbout, indent=4, sort_keys=True)

if __name__ == '__main__':
    commands = {'check', 'update', 'impl', 'tests'}

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

    args = parser.parse_args()
    command = args.command
    one_file = args.header_file is not None
    if args.database is not None:
        db_file = args.database
    backup_db = args.backup

    open_db()
    if one_file:
        os.chdir(os.path.dirname(os.path.abspath(args.header_file)))
        headers = [os.path.basename(args.header_file)]
    else:
        os.chdir(os.path.dirname(os.path.dirname(db_file)))
        headers = list(db)

    if command == 'check':
        for hf in headers:
            print('Checking', hf, '...')
            cname = get_cname(hf)
            if cname not in db:
                print("MISSING class name", cname)
                print()
                continue
            dbc = db['cname']
            funcs = get_funcs(hf, cname)
            print_diffs(dbc['functions'], ) #TODO HERE

            print()

    header_path = os.path.abspath(args.header_file)
    os.chdir(os.path.dirname(header_path))
    header_file = os.path.basename(header_path)

    funcs = get_funcs(header_file)
    for f in funcs:
        print("Function", f.short_name(), "has", len(f.args), "args")
