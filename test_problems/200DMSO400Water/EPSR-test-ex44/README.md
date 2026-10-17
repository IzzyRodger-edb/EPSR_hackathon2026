# EPSR-test-ex44 — readpdb with CONECT merging and dihedral groups

Program: `readpdb` (`fortran_src/readpdb.f`): converts a PDB file into an `.ato` file (molecules from residue
numbers and CONECT records, bonds from van der Waals radii, bond-angle constraints, optional DIH groups).

Input `mixture.pdb` was generated with a one-off helper script (not kept; the file itself is the input) from the DMSO and water coordinates of
EPSR-test plus a hand-made methanol; atom serial numbers start at 101 (CONECT look-up of serials):

* residue 1 = DMSO heavy atoms, residue 2 = water, residue 3 = the six DMSO H atoms: CONECT C1-H and C2-H
  records merge residue 3 into residue 1 although they are not contiguous (atom re-ordering branch);
* a TIP4P-like water with an `MW` site (unknown element -> dummy `QW` atom, `qtype=1`);
* a complete second DMSO (CONECT within one residue: nothing to merge), methanol split over two contiguous
  residues (contiguous merge), eight more waters; a short `CONECT 9999` record (ignored).

stdin (`readpdb.stdin`): van der Waals factor 0.7, dihedrals `y`. Needs `vanderWaalsRadii.txt` (copied from
EPSR-test) in the current directory.

The Fortran calls `dihedral(im,...)` with absolute atom numbers although `dihedral` adds `naf(im)-1`; this is
only correct for molecule 1, so the DIH angles of the later molecules use unrelated atoms (the methanol
values are NaN). The translation reproduces this.

Expected output: `mixture.ato` (non-cubic header, DIH rotational groups with their angles).

**Fortran reference.** The original `readxyz`/`readpdb` binaries crash (SIGSEGV) at the first call of `ctype`: the routine has optional arguments but is called with 11 of its 13 arguments and no explicit interface, so gfortran passes the hidden character lengths (3, 3) where the addresses of `radmaxsph`/`pseudofacsph` are expected and `present()` is true. The reference (`fortran_reference.tar.gz`) was therefore produced with `readxyz.f`/`readpdb.f` compiled with an explicit `INTERFACE` block for `ctype` added to the `readxyz`/`readpdb` subroutine (no other change; absent arguments then mean `radmaxsphere = 0.5*sigcor`, `pseudofacmol = 1`), linked with the standard reference objects (build script: `rust/tools/build_patched_readxyz_readpdb.sh <fref> <out>`). Run the case with that build folder: `rust/tools/run_ex_test.sh <this folder> <patched build> <work>`.

Verification: `ALL IDENTICAL` (needed the NaN fix in `fortran/format.rs`: gfortran writes `NaN` unsigned).
