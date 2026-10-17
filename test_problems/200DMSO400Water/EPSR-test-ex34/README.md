# EPSR-test-ex34 — `fmole` (prompts, repulsive core, missing template)

**Program:** `fmole`.

**Inputs:** `DMSO.ato`; `DMSOcore.ato` = `DMSO.ato` with `ecore dcore = 0.3 2.5` (so the
repulsive WCA core potential of `generate_mp` acts on the non-bonded H...H pairs);
`Water.ato` = `WaterMolecule.ato` renamed (its template `WaterMolecule.ato` is absent);
`fmole1.stdin`, `fmole2.stdin`.

**Commands:**
1. `fmole < fmole1.stdin` — file name typed (no listing) and number of shakes (150) asked on stdin.
2. `fmole Water < fmole2.stdin` — 30 shakes; `update_ato` cannot open the template
   (bonds kept from the file).

**Expected outputs:** updated `DMSOcore.ato`, `Water.ato`, `dirlist.txt` (empty).

**Verification:** ALL IDENTICAL.
