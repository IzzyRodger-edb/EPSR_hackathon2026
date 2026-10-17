# EPSR-test-ex55 — `voids`, two distributions, fresh random numbers

Program: `voids` (`fortran_src/voids.f`), input `cavity.VOIDS.dat`, single run:

* distribution 1: radius 2.0 A, all atoms occupy pixels (blank `include`; note the Fortran
  writes the blank value using the word boundaries of the previous keyword — 3 blanks).
* distribution 2: radius 2.5 A, `include O1 H1 XX` (XX not in the .ato file and silently
  dropped): only water atoms occupy pixels and occupied/unoccupied are swapped.
* No `randomno.dat` present: the generator is reset (idum = 0), file created and saved.

Pixel grid is adjusted to the radius (99^3 pixels), voids found by repeated
neighbour-labelling, void-size histogram + void-pixel g(r) (+ inter-void part) from 100
random origins.

Outputs compared: `cavity.VOIDS.n01`, `cavity.VOIDS.dat` (note the Fortran stores nsumt and
ndist into the `npix`/`voidmax` entries), `fort.14` (void volume fraction), `randomno.dat`.

Verification: ALL IDENTICAL.
