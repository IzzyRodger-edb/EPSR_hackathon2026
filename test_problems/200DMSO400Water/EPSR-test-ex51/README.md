# EPSR-test-ex51 — `partials` with a centre atom, rho(z) and rho(r)

Program: `partials`. Three inputs, run in sequence, whole list repeated twice:

* `centreO1.PARTIALS.dat`: `atom-c O1` (only pairs with O1 at the origin, rho = N1
  normalisation), nz 200, `rnormalise 0.8`, `nsumt -1` (no accumulation).
* `rhozS.PARTIALS.dat`: `atom-z S` — density profile |z| about S atoms (dorhoz), nsumt 0.
* `rhorO.PARTIALS.dat`: `atom-r O` — density in the xy-plane about O atoms (dorhor),
  `rnormalise 0.5`, nsumt 0.

Outputs compared: the three `.PARTIALS.g01` files, the rewritten `.dat` files, `fort.14`
(from the last run).

Verification: ALL IDENTICAL (passes 1 and 2).
