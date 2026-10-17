# EPSR-test-ex56 — `voids`, chained runs with accumulation

Program: `voids`, input `chain.VOIDS.dat` (radius 2.2 A, `include S O`), repeat 2.

`chain.VOIDS.n01` and `randomno.dat` were produced by one Fortran run of the same input with
`nsumt 0`; the input was then set to `nsumt 99`. Each pass therefore reads the previous
histograms (running average with nsum = 100) and continues the random-number sequence from
`randomno.dat`. (Because voids writes nsumt into the `npix` entry, nsumt itself stays 99 and
npix becomes 100 — the same 101-pixel starting grid — so the chain remains valid.)

Outputs compared: `chain.VOIDS.n01`, `chain.VOIDS.dat`, `fort.14`, `randomno.dat`.

Verification: ALL IDENTICAL (passes 1 and 2).
