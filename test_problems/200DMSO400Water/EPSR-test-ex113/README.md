# EPSR-test-ex113 — three-body potentials with ireset 1, more than 9 angles

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.

`.pcof`: appended after the pair potentials
```
tb O1  H1  O1  8.0 12.0 0.2 1.2 2.6 1 0 0 180.0 1.0 175.0 0.9 ... 135.0 0.1 130.0
tb O   S   C   4.0 10.0 0.1 1.8 2.0 1 1 1 106.0 1.0
```
The first group has 10 angle/amplitude pairs plus an odd trailing word: `(nwords-12)/2 = 10` is capped at 9.
`.EPSR.inp`: `ireset 1` (the Fortran resets it to 0 on saving, so pass 2 runs with `ireset 0`).
`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 2`.

Fortran quirk reproduced: with `ireset 1` `read_pcof` regenerates the default pair potentials, which resets
`nsymbol` to `nsymboldata + npairs`, so the `tb` symbols are no longer counted.  When the `ndata` line of the
`.EPSR.inp` is then processed, the per-dataset symbol values (pairs) are moved but the `tb` values are not,
so at save time the symbol slot read for each group does not hold the `tb` text (`nwords <= 12`): the groups are
not re-read and the in-memory values are written — width 1.0 rad → `57.30`, cosines/piconv → angles such as
`-57.30`.  Pass 2 reads these values back (angles of -57.3° etc.).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex113 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in every pass (byte-for-byte after masking the `.out` date line and the fort.35 timings).
