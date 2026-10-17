# EPSR-test-ex112 — three-body potentials on DMSO, zero amplitudes, molecule shakes

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.

`.pcof`: appended after the pair potentials
```
tb S   C   O   0.0 10.0 0.1 2.0 2.0 1 1 1 107.0 1.0
tb C   S   C   5.0 10.0 0.1 2.0 2.0 1 1 1 97.0 0.0
tb C   H   O1  15.0 20.0 0.1 1.2 2.8 1 0 0 180.0 1.0 150.0 0.5
```
* group 1 has amplitude 0 and group 2 has all angle amplitudes 0 (`read_pcof` then sets the amplitude to 0):
  triplets are found but contribute nothing (`abs(threebamp) > 0` tests, "Energy for each triplet" = 0);
* group 3: methyl C–H···O(water) contacts; the H atoms move in the methyl internal rotations, so
  `save_threeb_energies` after internal rotations is exercised.

`.EPSR.inp`: `ntimes 6`, `intra 3` (several molecule shakes per pass: intramolecular moves do not update the
triplet energies, `fmolec1` is re-run after the shake).  `commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 2`.

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex112 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in every pass (byte-for-byte after masking the `.out` date line and the fort.35 timings).
