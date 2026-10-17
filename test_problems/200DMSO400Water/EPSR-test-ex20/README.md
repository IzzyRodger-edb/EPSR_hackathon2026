# EPSR-test-ex20 — `nwts`, single samples (the EPSR-test weights)

Program: `nwts` (`fortran_src/nwts.f`, keyword file `.NWTS.dat`).

Inputs: `200DMSO400Water.ato` and the three `DMSODWater{D,H,Null}.NWTS.dat` files of
`EPSR-test` (unchanged). Each is one sample (`nsamples 1`, `normtot 0`, per atom):
D6-DMSO in D2O (H and exchangeable H1 = deuterium), in H2O (H1 natural) and in the
HDO "null" mixture (H1 = 0.64 natural + 0.36 D, i.e. two isotopes; H1 exchangeable).

Exercises: reading the keyword file (the `fnameato` branch of `check_symbol` that
sets the number of components, labels and symbols from the `.ato` file), the
isotope abundance lists, `scatt_len` (Sears 1991 table), the averaging of
exchangeable H, the total weights file `<run>tot.wts`, rewriting of the `.ato`
file (`moloutc`) and of the input file (`save_files_others`).

Expected outputs: `DMSODWater{D,H,Null}.NWTStot.wts` — identical to the files
shipped in `EPSR-test`; the `.NWTS.dat` and `.ato` files are rewritten unchanged.

Verification: `rust/tools/run_ex_test.sh test_problems/200DMSO400Water/EPSR-test-ex20 <fortran build>` →
**ALL IDENTICAL** (reference in `fortran_reference.tar.gz`).
