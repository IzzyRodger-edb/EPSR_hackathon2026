# EPSR-test-ex38 — `makelatticeato`

**Program:** `makelatticeato` (`fortran_src/makelatticeato.f`): replicates the contents
of a `.ato` box `na x nb x nc` times, optionally converting every atom to a tethered
one-atom molecule.

**Inputs:** `DMSO.ato`, `WaterMolecule.ato` (from `EPSR-test`), `mla3.stdin`, `mla4.stdin`.

**Commands:**
1. `makelatticeato DMSO 2 2 2 dmso222 n` — cubic box kept (dmax doubled), bonds and
   rotational groups renumbered.
2. `makelatticeato WaterMolecule.ato 3 2 1 waterlat y` — non-cubic result, conversion
   to tethered atoms (`T00000` flags, new molecule types `moltype01..`).
3. `makelatticeato < mla3.stdin` — everything asked on stdin, convert = `y`.
4. `makelatticeato DMSO.ato 1 1 2 dmso112 < mla4.stdin` — convert asked (`n`).

**Fortran quirk reproduced:** the convert answer is compared over the column range of
the first word of the previous `parse` (the Fortran parses an uninitialised variable).

**Expected outputs:** `dmso222.ato`, `waterlat.ato`, `dmso121.ato`, `dmso112.ato`, `dirlist.txt`.

**Verification:** ALL IDENTICAL.
