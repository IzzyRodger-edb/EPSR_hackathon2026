# EPSR-test-ex31 — `mixato` (basic)

**Program:** `mixato` (`fortran_src/mixato.f`): mixes single-molecule-type `.ato` files
into one box.

**Inputs:** `DMSO.ato`, `WaterMolecule.ato` (copied from `EPSR-test`), `mixato.stdin`.
No `system_commands.txt`, so `search_for_files` lists nothing (the empty `system_ls`
command fails, `dirlist.txt` is created empty) and the file names are typed.

**Command:** `mixato < mixato.stdin` — 2 files (number asked on stdin), 20 DMSO + 40 water
(the second name given with its `.ato` extension), density 0.1 atoms/Å³, output `DMSOWaterMix`.

**Exercises:** repeated `minitc` reads accumulating atom/molecule types, copying of
bonds and rotational groups, template-file check, box size from density,
`moloutc` (molecule-scale factors reset to 1 for the added molecules).

**Expected outputs:** `DMSOWaterMix.ato`, `dirlist.txt` (empty).

**Verification:** ALL IDENTICAL.
