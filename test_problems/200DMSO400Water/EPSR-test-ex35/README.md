# EPSR-test-ex35 — `randomise` (+ mixato, fmole)

**Program:** `randomise` (`fortran_src/introtcluster.f`): random position and random
whole-molecule rotation (`rotatm`) of every molecule, using the `.ato` random-number state.

**Inputs:** `DMSO.ato`, `WaterMolecule.ato`, `200DMSO400Water.ato` (from `EPSR-test`),
`mixato.stdin`, `randomise.stdin`.

**Commands:**
1. `mixato < mixato.stdin` — 15 DMSO + 30 water at 0.09 atoms/Å³ → `box.ato`.
2. `randomise box` — randomise the freshly mixed box.
3. `fmole box 20`.
4. `randomise < randomise.stdin` — file name typed on stdin.
5. `randomise 200DMSO400Water.ato` — the 600-molecule EPSR test configuration.

**Expected outputs:** `box.ato`, `200DMSO400Water.ato`, `dirlist.txt`.

**Verification:** ALL IDENTICAL.
