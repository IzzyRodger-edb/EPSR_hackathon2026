#!/usr/bin/env bash
# Build the original Fortran EPSR programs (version 26) from ../../fortran_src
# into a separate build folder, for use as references against which the Rust
# port is validated.
#
# usage: build_fortran_reference.sh BUILD_DIR [--vectorize]
#
# Produces BUILD_DIR/epsr and BUILD_DIR/<program> for every program listed in
# tools/fortran_link_lists.txt (partials, coord, nwts, makeato, ...).  The
# sources are copied to BUILD_DIR; fortran_src is never modified.
#
# By default the code is compiled with `-O2 -fno-tree-vectorize`.  gfortran
# >= 12 auto-vectorises some loops at -O2 and then calls glibc's SIMD math
# library (libmvec, e.g. `_ZGVbN4v_sinf`) whose results differ from the
# scalar `sinf`/`expf` by a few ulp; older compilers (used for the official
# EPSR builds) produce scalar code.  The Rust port reproduces the scalar
# build bit for bit.  Pass --vectorize to get gfortran's default behaviour.
#
# Run a program as:  ulimit -s unlimited; OMP_NUM_THREADS=1 BUILD_DIR/<prog> <folder>/ <command> <input file>
set -euo pipefail
BUILD=$1
VEC=${2:-}
HERE=$(cd "$(dirname "$0")" && pwd)
SRC=$(cd "$HERE/../../fortran_src" && pwd)
LISTS="$HERE/fortran_link_lists.txt"
mkdir -p "$BUILD"
cp "$SRC"/*.inc "$BUILD"/
cd "$BUILD"
FLAGS="-c -O2 -fno-backslash -ffixed-line-length-none -Wno-unused-label -fopenmp -w"
if [ "$VEC" != "--vectorize" ]; then FLAGS="$FLAGS -fno-tree-vectorize"; fi
# Modules must be compiled first, in dependency order.
MODULES="lngamma_mod.f90 reallocation_routines.f90 molecom_mod.f90 stringutility.f90 threebody_mod.f90
coordconstraint_mod.f90 containerforces_mod.f90 cellneighbourlist.f ergpairreference.f ato_manipulation_rev4.f
FCreadline_windows.f90 pfroutines_mod.f90 molecule_coordinates_mod.f90 sharm_routines_rev1.f sdfcube.f dockato.f"
compile() {
  local s=$1
  cp "$SRC/$s" .
  # sharm_diffuse (ato_manipulation_rev4.f, called by epsr for the `diffuse`
  # keyword) has the same run-time FORMAT "(1h(i1,1h,,..." as sharm_rev2.f
  # (see below): insert the missing comma, otherwise epsr stops with "Missing
  # comma between descriptors".  sharm_diffuse is only called by epsr.
  if [ "$s" = ato_manipulation_rev4.f ]; then
    sed -i "s/formatout='(1h('\/\/il1fmt/formatout='(1h(,'\/\/il1fmt/" "$s"
    grep -q "(1h(,'//il1fmt" "$s"
  fi
  gfortran $FLAGS "$s"
}
for s in $MODULES; do compile "$s"; done
OBJS=$(grep -v '^#' "$LISTS" | cut -d: -f2 | tr ' ' '\n' | grep '\.o$' | sort -u)
for o in $OBJS; do
  [ -f "$o" ] && continue
  b=${o%.o}
  if [ -f "$SRC/$b.f90" ]; then compile "$b.f90"; else compile "$b.f"; fi
done
grep -v '^#' "$LISTS" | while IFS=: read -r prog objs; do
  name=$prog
  case $prog in sdf_rev2) name=sdf;; sharm_rev2) name=sharm;; chains_sp_criterion) name=chains;;
    rings_sp_criterion) name=rings;; mapgrquad) name=mapgr;; introtcluster) name=randomise;; bondsclus) name=bonds;; esac
  gfortran -fopenmp -o "$name" $objs
done
# readxyz.f and readpdb.f call ctype with 11 of its 13 (optional) arguments
# without an explicit interface and crash on every input.  Install builds with
# an explicit interface as the reference (originals kept as *.orig).
"$HERE/build_patched_readxyz_readpdb.sh" "$BUILD" "$BUILD/patched" > /dev/null
for p in readxyz readpdb; do mv "$p" "$p.orig"; cp "patched/$p" "$p"; done
# sharm_rev2.f and sdf_rev2.f build a run-time FORMAT "(1h(i1,1h,,..." with no
# comma after the Hollerith 1h( ; gfortran rejects it at run time, so the
# coefficient files are never written.  Reference = the same source with the
# comma inserted (originals kept as *.orig).
for p in sharm:sharm_rev2 sdf:sdf_rev2; do
  name=${p%%:*}; src=${p##*:}
  sed "s/formatout='(1h('\/\/il1fmt/formatout='(1h(,'\/\/il1fmt/" "$SRC/$src.f" > ${src}_patched.f
  grep -q "(1h(,'//il1fmt" ${src}_patched.f
  gfortran $FLAGS ${src}_patched.f -o ${src}_patched.o
  objs=$(grep "^$src:" "$LISTS" | cut -d: -f2 | sed "s/ $src.o / ${src}_patched.o /; s/^$src.o /${src}_patched.o /")
  mv "$name" "$name.orig"; gfortran -fopenmp -o "$name" ${src}_patched.o $(echo $objs | tr ' ' '\n' | grep -v "^$src.o$\|^${src}_patched.o$")
done
echo "built Fortran reference programs in $BUILD ($( [ "$VEC" = "--vectorize" ] && echo vectorised || echo scalar math ))"
