# EPSR-test-ex40 — readxyz, cell read from standard input

Program: `readxyz` (`fortran_src/readxyz.f`): converts an `.xyz` file into an `.ato` file in which every
atom is a separate one-atom molecule (non-cubic box header, default Lennard-Jones parameters per element,
one molecule type per element label).

Inputs: `DMSO.xyz`, `WaterMolecule.xyz` copied from EPSR-test. Their 2nd line ("Filename: ...") does not
contain `<word> a b c alpha beta gamma`, so the cell is read from standard input:

* `readxyz DMSO.xyz < readxyz_dmso.stdin` — triclinic cell 20 22 25 / 80 95 100 (exercises
  `convertalphabetagamma` and the general cell vectors);
* `readxyz WaterMolecule < readxyz_water.stdin` — no extension given; the first stdin line is not numeric
  (read error -> asked again), then the six values are spread over two lines (list-directed read), cubic cell 30 Å.

Expected outputs: `DMSO.ato`, `WaterMolecule.ato`.

**Fortran reference.** The original `readxyz`/`readpdb` binaries crash (SIGSEGV) at the first call of `ctype`: the routine has optional arguments but is called with 11 of its 13 arguments and no explicit interface, so gfortran passes the hidden character lengths (3, 3) where the addresses of `radmaxsph`/`pseudofacsph` are expected and `present()` is true. The reference (`fortran_reference.tar.gz`) was therefore produced with `readxyz.f`/`readpdb.f` compiled with an explicit `INTERFACE` block for `ctype` added to the `readxyz`/`readpdb` subroutine (no other change; absent arguments then mean `radmaxsphere = 0.5*sigcor`, `pseudofacmol = 1`), linked with the standard reference objects (build script: `rust/tools/build_patched_readxyz_readpdb.sh <fref> <out>`). Run the case with that build folder: `rust/tools/run_ex_test.sh <this folder> <patched build> <work>`.

Verification: `ALL IDENTICAL` (Rust `epsr-rs ./ readxyz ...` vs patched Fortran).
