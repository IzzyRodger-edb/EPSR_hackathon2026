#!/usr/bin/env bash
# Build Fortran reference binaries for readxyz/readpdb with an explicit INTERFACE for ctype
# (the original call passes 11 of 13 arguments without an interface and crashes with SIGSEGV).
# Other group-C programs (writexyz, bonds, alignmole) are the unmodified reference binaries.
# usage: build_patched_readxyz_readpdb.sh FREF_DIR OUT_DIR   (FREF_DIR = folder built by build_fortran_reference.sh)
set -euo pipefail
FREF=$(realpath "$1"); OUT=$(mkdir -p "$2" && realpath "$2")
HERE=$(cd "$(dirname "$0")" && pwd)
SRC=$(cd "$HERE/../../fortran_src" && pwd)
LISTS="$HERE/fortran_link_lists.txt"
cd "$OUT"; cp "$FREF"/*.inc "$FREF"/*.mod .
for p in readxyz readpdb; do
  cp "$SRC/$p.f" ${p}_patched.f
  python3 - "$p" <<'PY'
import sys
p=sys.argv[1]; f=p+'_patched.f'; s=open(f).read()
iface='''      INTERFACE
         subroutine ctype(a1,ats,iexch,eread,wread,amasst,qt,rqt,natomsph,rhosph,radminsph,radmaxsph,pseudofacsph)
         character*3 a1,ats
         integer iexch
         real eread,wread,amasst,qt,rqt
         real, intent(in),optional :: natomsph,rhosph,radminsph,radmaxsph,pseudofacsph
         end subroutine ctype
      END INTERFACE
'''
key='      subroutine %s(readfile,rho,ierr)\n\n      use molecom_mod\n'%p
assert key in s
open(f,'w').write(s.replace(key,key+iface,1))
PY
  gfortran -c -O2 -fno-backslash -ffixed-line-length-none -Wno-unused-label -fopenmp -w -fno-tree-vectorize ${p}_patched.f -o $p.o
  objs=$(grep "^$p:" "$LISTS" | cut -d: -f2)
  gfortran -fopenmp -o $p $(for o in $objs; do [ "$o" = "$p.o" ] && echo "$OUT/$o" || echo "$FREF/$o"; done)
done
for p in writexyz bonds alignmole; do ln -sf "$FREF/$p" "$p"; done
echo "built patched readxyz/readpdb in $OUT"
