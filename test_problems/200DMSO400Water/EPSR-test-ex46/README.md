# EPSR-test-ex46 — bonds (bondsclus.f)

Program: `bonds` (`fortran_src/bondsclus.f`): for each run of consecutive molecules of the same type writes
the matrix of mean intramolecular distances (upper triangle) and r.m.s. spreads (lower triangle) to `.bonds`.

* `bonds 200DMSO400Water.ato` — 200 DMSO + 400 water (two tables, averages over 200/400 molecules);
* `bonds DMSO`, `bonds WaterMolecule.ato` — single molecules (spread 0);
* `bonds alternating.ato` — 7 molecules D D W W D W D (randomly rotated and slightly distorted copies of
  DMSO/water, built with a one-off helper script that was not kept; molecule order DDWWDWD): a new table each time the type changes.

Expected outputs: `200DMSO400Water.bonds`, `DMSO.bonds`, `WaterMolecule.bonds`, `alternating.bonds`.

Verification: `ALL IDENTICAL` against the original Fortran reference binary.
