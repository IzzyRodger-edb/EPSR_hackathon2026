# EPSR-test-ex95 — data types nrtype 1, 3 and 4

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets); only the files named below differ.

The three datasets are read in the three non-GUDRUN formats of `runepsr17` ("read in data to be fitted", `epsr_standalone_rev3.f`):

* data 1 `D6DMSO_D2O.nrt1`, **nrtype 1**: list-directed `Q S(Q)` (errors set to 0), no header, point-mode Q (bin centres of the `.mint01` histogram). Note the Fortran reads the first record, then `backspace`s, so the first line is data; the last point is kept.
* data 2 `D6DMSO_HDO.nrt3`, **nrtype 3**: GENIE format `(i7,3g14.7)`, histogram Q (bin boundaries as in the `.mint01` file, converted to centres by EPSR); 4 header lines (the Fortran reads record 1 + 4 records and backspaces, so record 5 is the first data line); last point ignored.
* data 3 `D6DMSO_H2O.nrt4`, **nrtype 4**: GENIE format, point-mode Q, plus a deliberately rubbish last point (S=9.99) which EPSR ignores.

Inputs: `scripts/make_data.py` creates the three files from the base `.mint01` files (kept in the folder); `scripts/make_inp.sh` sets `datafile`/`nrtype` in the `.EPSR.inp` (via `scripts/setkw.py`). All three files give 594 points (log: `Read 594 points from file ...`), the same as the original `.mint01` data.

Rust: `epsr-core/src/datafile.rs` (`read_data_points`, incl. the fixed-format `(i7,3g14.7)` read with BLANK=NULL/PAD semantics and the iostat chain of the header reads).

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 3` (each pass reads back the files written by the previous one).

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex95 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).

## Verification
Rust port vs Fortran reference: **ALL IDENTICAL** in all 3 passes (byte-for-byte after masking the `.out` date line and the fort.35 timings).
