# EPSR-test-ex45 — readpdb, triclinic cell from stdin, no dihedrals

Program: `readpdb` (`fortran_src/readpdb.f`).

Input `cluster.pdb` (hand-built with a one-off helper script that was not kept; the file itself is the input): a valid triclinic CRYST1 record, followed by a second
CRYST1 record whose angle field cannot be read -> cell asked on stdin (21 23.5 26.25 / 80 96 110); a water,
a DMSO, a TIP4P-like water (dummy site) and a methanol, each in its own residue (no CONECT records).

stdin: cell, van der Waals factor `0` (invalid -> default 0.7), dihedrals `no`.

Expected output: `cluster.ato`.

**Fortran reference.** The original `readxyz`/`readpdb` binaries crash (SIGSEGV) at the first call of `ctype`: the routine has optional arguments but is called with 11 of its 13 arguments and no explicit interface, so gfortran passes the hidden character lengths (3, 3) where the addresses of `radmaxsph`/`pseudofacsph` are expected and `present()` is true. The reference (`fortran_reference.tar.gz`) was therefore produced with `readxyz.f`/`readpdb.f` compiled with an explicit `INTERFACE` block for `ctype` added to the `readxyz`/`readpdb` subroutine (no other change; absent arguments then mean `radmaxsphere = 0.5*sigcor`, `pseudofacmol = 1`), linked with the standard reference objects (build script: `rust/tools/build_patched_readxyz_readpdb.sh <fref> <out>`). Run the case with that build folder: `rust/tools/run_ex_test.sh <this folder> <patched build> <work>`.

Verification: `ALL IDENTICAL`.
