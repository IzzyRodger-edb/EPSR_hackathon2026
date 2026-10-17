# EPSR-test-ex57 — `fluctuations`, three sphere sizes, three chained runs

Program: `fluctuations` (`fortran_src/fluctuations.f`), input `numfluct.FLUCTUATIONS.dat`:
3000 random sphere centres, radii 4, 6.5 and 10 A, nsize 60, nsumt 0, repeat 3.
No `randomno.dat` initially (generator reset), then chained through `randomno.dat`.

Outputs compared: `numfluct.FLUCTUATIONS.n01` (distribution of N - <N>, accumulated),
`numfluct.FLUCTUATIONS.dat` (nsumt and the running real*8 averages `ntotal`, `nsqtotal`,
`szero` written list-directed), `randomno.dat`.

Verification: ALL IDENTICAL (passes 1, 2 and 3).
