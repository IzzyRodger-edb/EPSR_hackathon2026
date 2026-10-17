# EPSR-test-ex41 — readxyz, cell from the file, more elements, error paths

Program: `readxyz` (`fortran_src/readxyz.f`).

* `readxyz mixture.xyz` — hand-made file; line 2 `Cell 18.5 19.0 21.2 88.0 92.5 105.0 ...` gives the cell
  (words 2-7). 12 atoms of N, H, Na, Cl, C, O, K, Br, F (all branches of the default-parameter selection:
  H, O, N, C and "other"); the coordinate list is ended by a blank line (fewer than 4 words).
* `readxyz titled.xyz < titled.stdin` — line 2 has 7 words that are not numbers (read error) -> cell from
  stdin (hexagonal 12 13 14 / 90 90 120).
* `readxyz badlabel.xyz` — atom label `Xx` is not a chemical symbol -> "Problems with .xyz file", no output.
* `readxyz short.xyz` — fewer atoms than announced on line 1 -> no output.

Expected outputs: `mixture.ato`, `titled.ato` (nothing for the two error cases).

**Fortran reference.** The original `readxyz`/`readpdb` binaries crash (SIGSEGV) at the first call of `ctype`: the routine has optional arguments but is called with 11 of its 13 arguments and no explicit interface, so gfortran passes the hidden character lengths (3, 3) where the addresses of `radmaxsph`/`pseudofacsph` are expected and `present()` is true. The reference (`fortran_reference.tar.gz`) was therefore produced with `readxyz.f`/`readpdb.f` compiled with an explicit `INTERFACE` block for `ctype` added to the `readxyz`/`readpdb` subroutine (no other change; absent arguments then mean `radmaxsphere = 0.5*sigcor`, `pseudofacmol = 1`), linked with the standard reference objects (build script: `rust/tools/build_patched_readxyz_readpdb.sh <fref> <out>`). Run the case with that build folder: `rust/tools/run_ex_test.sh <this folder> <patched build> <work>`.

Verification: `ALL IDENTICAL`.
