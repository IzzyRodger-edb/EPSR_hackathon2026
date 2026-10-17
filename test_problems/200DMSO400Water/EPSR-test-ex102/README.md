# EPSR-test-ex102 — coarse graining with modified sphere parameters and sizefactor > 1

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets).

`.EPSR.inp`: `qinscalefc 3.0 0.5`, `ntimes 3`, `sizefactor 1.2 0.9 1.0E+06` (molecules shrunk; the threshold
is so large that the size factor declines every iteration: 1.2 -> 1.08 -> 1.0 (clamped), which also stops
intramolecular rotations while > 1).
`.ato`: DMSO `scalefacmol` = 1.00..1.40 in steps of 0.05, water 1.0.  Single-precision quirk: `int((1.4-1.0)/0.2)+1` = 2 (not 3) bins for the DMSO types.  Atom-type coarse-graining
columns: S `-0.08 0.3 1.9 2.0` (negative `natomsph` = a density, converted with `cgmolfac`), C `1.0 0.5 2.0 2.0`
(`radminsph` > 0), O `1.0 0.0 1.5 2.0`, H `1.0 0.0 0.0 2.0` (radius 0 -> density 0), O1 `-0.05 0.0 1.7 1.0`,
H1 `1.0 0.2 0.6 0.5` (`cgmolfac` < 1 is reset to 1 on reading).  `.pcof`: the `rminex` values of EPSR-test-ex1.

`commands.txt`: `epsr 200DMSO400Water.EPSR.inp`, `repeat 2` (the second pass reads back the `.ato` with its per-molecule scale factors, the `.EPSR.inp`, `.pcof`, `.inpa` and accumulators written by the first).

## How the inputs were made
Starting from the 16 base files of `EPSR-test`, `scripts/make_cg.py <folder> <mode>` edits the `.EPSR.inp`
keywords, the 10th word (`scalefacmol`) of every molecule header line of the `.ato` and, for some modes,
the coarse-graining columns (`natomsph radminsph radmaxsph cgmolfac`) of the atom-type lines.

## Coarse-graining code exercised (Fortran -> Rust `epsr-core/src/coarsegrain.rs`)
* `runepsr17`: `coarsegraining = qscaleinputfactor > 1.0` from `qinscalefc`.
* `inita` (`epsr_standalone_rev3.f`) `if(coarsegraining)` block: `scalefacmoluse = totalatoms**selfbroadpwr`,
  per-type scale-factor bins (`scalefacratio = 0.2`, `nsffbintype`, `assignatomscalefac`), spherical form
  factors `sphereformfactor` (`gaussianqr2`), normalisation by the average number of atoms per molecule,
  intramolecular form factors `intraformfactor` / `intraformfactorself`, and the list-directed files
  `formfactors.txt`, `intraformfactors.txt`, `intraformfactorsself.txt` (compared byte for byte).
* `en_mole_pair` / `en_change` (`molfit.f`): `scalefac = 0.5*(scalefacmol(im)+scalefacmol(jm))` scales the
  distance used for the reference-potential part of `en_atom_inter_scale` (affects MC acceptance, energies, `.erg`).
* `atogen` (`ato_manipulation_rev4.f`): same `scalefac` on the `rminex` close-contact count `potex`, and the
  `wavinter` histograms per form-factor bin (normalised by `cprod2`).
* NB (Fortran): the form factors and `wavinter` are computed but never used for S(Q): the coarse-graining code
  in `gensq` is commented out (it only fills the unused temporary `prodssff`), so `gensq` is unchanged.

## Verification
`rust/tools/run_ex_test.sh` (gfortran -O2 -fno-tree-vectorize reference, `fortran_reference.tar.gz` = Fortran
output after the last pass): **ALL IDENTICAL** in both passes, including the three form-factor files.
