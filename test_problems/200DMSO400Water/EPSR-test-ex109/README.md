# EPSR-test-ex109 — starting from scratch: no `.pcof`, no accumulators; 7-value `.ato` step line; no "real" atoms

Program: `epsr`. Inputs made by `scripts/make_inputs.sh` from the EPSR-test files.
`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `epsr nosph.EPSR.inp`, `repeat 3`.

Note: the Fortran `epsr` command needs an existing `.EPSR.inp` (`setup_input_file` returns if it cannot be
opened; new input files are made by the interactive EPSR shell, not by `epsr`). The first thing a user runs is
therefore an `.EPSR.inp` + `.ato` + data/weights with nothing else:

* `200DMSO400Water`: only `.EPSR.inp`, `.ato`, data and weights (no `.pcof`, `.EPSR.inpa`, `.erg`, `.g01`,
  `.y01`, `.f01`, `.s01`). `read_pcof` cannot open the `.pcof` ("Problems with specified .pcof file"):
  `ncoeffp`, `roverlap`, `rminfac` are taken from the symbol table, `npairspcof = 0` forces
  `ireset = iinit = 1` and the default pair potentials (`rminex = max(rminfac*dcor, roverlap)`, zero EP);
  `inita` initialises all accumulators, `initg` finds no `.?01` files, the `.erg` is created. Every output
  file (incl. `.pcof`, `.inpa`) is created and read back in passes 2 and 3.
  Also: `ncoeffp 2000` in the `.EPSR.inp` -> `check_symbol` caps it to mcoeff = 1000 and rewrites the symbol
  (the Fortran message prints the value with i3, i.e. `***`); the `.ato` step-size line has 7 values, the 7th
  being `tempdihreq` (`minitc` `nwords.eq.7` branch; the file is rewritten with 10 values).
  Quirk: the default pair potentials do not set `efacp`, so the EP would not be refined in the first run
  (here ereq = 0 anyway, the `.pcof` defaults).
* `nosph` (`scripts/clone_run.py`, `scripts/zero_natomsphere.py`): same scratch start, with the
  coarse-graining value `natomsphere` (6th number of each atom-type line) set to 0 for every type, so
  `natomrcg = natomcg = 0` in `minitc`: the atom fractions `cfrac`/`cfraccg` are left un-normalised (counts)
  and several outputs contain NaN (e.g. `.erg` columns 6-8, `.r01`). This is what the Fortran does; it is
  reproduced, NaNs included.

## How to run
```
rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex109 <fortran build dir> [work dir] [--save-reference]
```
The Fortran output files after the last pass are stored in `fortran_reference.tar.gz`
(gfortran -O2 -fno-tree-vectorize reference build, `rust/tools/build_fortran_reference.sh`).
Inputs are generated from the base EPSR-test files by `scripts/` (see `scripts/make_inputs.*`).

## Verification
gcov: `npairspcof.le.0` branch of the missing-file case, default pair potentials, `nwords.eq.7` in `minitc`.
Rust vs Fortran: **ALL IDENTICAL** in all three passes.
