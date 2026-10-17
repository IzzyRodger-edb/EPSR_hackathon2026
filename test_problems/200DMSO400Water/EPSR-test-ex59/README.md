# EPSR-test-ex59 — `fluctuations` without accumulation, seeded random numbers

Program: `fluctuations`, input `single.FLUCTUATIONS.dat`: one sphere size (7.5 A),
1500 samples, `nsumt -1` (the nsumts < 0 branch: no accumulation, nsumt stays -1), initial
`ntotal`/`nsqtotal` values given (ignored since nsumt = 0 internally), repeat 2.
`randomno.dat` (copied from ex56) seeds the generator, so initran1 reads a saved state.

Outputs compared: `single.FLUCTUATIONS.n01`, `single.FLUCTUATIONS.dat`, `randomno.dat`.

Verification: ALL IDENTICAL (passes 1 and 2).
