# EPSR-test-ex36 — `addato`

**Program:** `addato` (`fortran_src/addato.f`): inserts copies of molecules at random
positions/orientations into an existing box (box size unchanged).

**Inputs:** `DMSO.ato`, `WaterMolecule.ato`, `200DMSO400Water.ato` (from `EPSR-test`),
`addato1.stdin`, `addato2.stdin`.

**Commands:**
1. `addato 2 < addato1.stdin` — add 5 DMSO (name given with `.ato`) and 10 water to
   `200DMSO400Water` → `bigbox.ato` (615 molecules).
2. `addato < addato2.stdin` — number of files asked on stdin; a non-existent file to add
   (`NoFile`, retried) and a non-existent target (`NoFile2`, retried); add 4 water to
   `bigbox.ato` → `bigbox2.ato`.

**Notes:** as in the Fortran, the template of every added molecule type becomes the name
of the file it was read from (the existence test opens the blank-padded
`character*256` template name, which never exists).

**Expected outputs:** `bigbox.ato`, `bigbox2.ato`, `dirlist.txt`.

**Verification:** ALL IDENTICAL.
