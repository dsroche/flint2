#!/bin/bash --norc

set -e
set -u

srcdir=$(dirname "$(readlink -f "$0")")

function usage {
  echo "Usage: $0"
  echo "Removes all the extra files and gets ready for \"public view\" of fmpz_spoly"
  echo "WARNING: This script will self-destruct!"
  [[ $# -eq 1 ]] && exit $1
}

[[ $# -eq 0 ]] || usage 1

read -p "WARNING: This script will delete stuff! You really want to do it? [y/N]?" -n1 yn
[[ -z $yn ]] && yn=n || echo
if [[ ! $yn = y ]]; then
    echo "ABORTING"
    exit
fi

cd "$srcdir"

exec 4<<DELETEFILES
TODO.txt
ourtest
Makefile
cleanup.sh
DELETEFILES

while read -u4 filename; do
    if [[ $filename =~ ^# ]]; then
        echo "comment: $filename"
    elif [[ $filename =~ ^[/.] ]]; then
        echo "ILLEGAL FILENAME: $filename (skipping)"
    elif [[ -d $filename ]]; then
        rm -r "$filename/"
    elif [[ -e $filename ]]; then
        rm "$filename"
    fi
done

exec 4<&-

sed -i 's/^\#\( *fmpz_spoly\)/ \1/' ../Makefile.in
