# EPSR-test-ex101 — coarse graining with distributed molecule scale factors

Program: `epsr`. Based on EPSR-test (200 DMSO + 400 water, 3 neutron datasets).

`.EPSR.inp`: `qinscalefc 1.5 0.33` (coarse graining on, `selfbroadpwr` 0.33 -> `scalefacmoluse` =
(atoms per molecule)**0.33), `ntimes 3`.
`.ato`: `scalefacmol` = 1.0/1.25/1.5/1.7 for DMSO and 0.8/1.0/1.2 for water (cycled by molecule number), so
there are 4 (DMSO types) and 3 (water types) form-factor bins, bin averaging (`ntermssff > 1`) and type minima
below 1.0.  `.pcof`: the `rminex` values of EPSR-test-ex1 (S-S 4.0, two water pairs) so the scaled `potex`
count of `atogen` matters.  The pair scale factor changes the MC energies (first-pass energy/molecule
-30.18 vs -39.51 kJ/mol for the same input with `qinscalefc 1.0`).

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
