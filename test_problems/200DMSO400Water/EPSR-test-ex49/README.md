# EPSR-test-ex49 — readpdb with single-atom residues

Program: `readpdb` (`fortran_src/readpdb.f`).

Input `ions.pdb` (hand-made): three waters, a lone Na and a lone Cl residue (element names in columns 13-14,
`NA`/`CL` -> upper-case symbol match) and a `HETATM` record (ignored: only `ATOM` records are read).
stdin: van der Waals factor 0.75, no dihedrals.

Quirk of the Fortran reproduced here: for a residue with one atom `nal` is never set (it is only assigned
when a second atom is added), so it keeps 0; the loops `naf..nal` are empty and the `.ato` file contains the
ion molecules with `nato = 1` but no atom records (the ions are effectively lost).

Expected output: `ions.ato`.

**Fortran reference.** The original `readxyz`/`readpdb` binaries crash (SIGSEGV) at the first call of `ctype`: the routine has optional arguments but is called with 11 of its 13 arguments and no explicit interface, so gfortran passes the hidden character lengths (3, 3) where the addresses of `radmaxsph`/`pseudofacsph` are expected and `present()` is true. The reference (`fortran_reference.tar.gz`) was therefore produced with `readxyz.f`/`readpdb.f` compiled with an explicit `INTERFACE` block for `ctype` added to the `readxyz`/`readpdb` subroutine (no other change; absent arguments then mean `radmaxsphere = 0.5*sigcor`, `pseudofacmol = 1`), linked with the standard reference objects (build script: `rust/tools/build_patched_readxyz_readpdb.sh <fref> <out>`). Run the case with that build folder: `rust/tools/run_ex_test.sh <this folder> <patched build> <work>`.

Verification: `ALL IDENTICAL`.
