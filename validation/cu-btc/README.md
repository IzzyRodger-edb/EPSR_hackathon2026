# Cu-BTC (HKUST-1): a crystalline metal-organic framework from a PDB file

This case validates the Rust port and the Python wrapper on a crystal fitted to real
neutron data: the empty Cu-BTC framework measured on NIMROD, from the Dissolve
example in `../../experimental_data/cu-btc`. It differs from the water and silica
cases in three ways: the configuration is read from a PDB file instead of being
built with `mixato`, the force field is applied with `changeato`, and EPSR
calculates Bragg peaks.

| | |
|---|---|
| system | one cubic unit cell, a = 26.3336 Å, 624 atoms (48 Cu, 192 O, 288 C, 96 H), 300 K |
| data | `Empty-CuBTC.mint01` (NIMROD, Q = 0.25 to 50 Å⁻¹) |
| molecules in EPSR | 24 Cu₂ paddlewheels (2 atoms) and 32 benzene-1,3,5-tricarboxylate linkers (18 atoms) |
| atom types | `CU`, `O`, `C1` (carboxylate C), `C2` (ring C bonded to C1), `C3` (ring C-H), `H` |

## From the Dissolve input to EPSR inputs

**Configuration (`make_pdb.py` → `CuBTC.pdb` → `readpdb` → `CuBTC.ato`).**
The Dissolve species is the whole framework as one 624-atom molecule with periodic
Cu–O bonds. EPSR cannot hold a periodic network as one molecule, so
`make_pdb.py` splits it at the Cu–O bonds:

* each linker (C₉H₃O₆, held together by the Dissolve C1–O, C1–C2, C2–C3 and C3–H
  bonds) becomes one PDB residue, written with its atoms in the order O, C1, C2, C3, H;
* each Cu₂ paddlewheel (two Cu atoms 2.629 Å apart; the next Cu is 8 Å away) becomes
  one residue. Single-atom residues cannot be used: `readpdb` writes them without
  their atom record, in the Fortran and, faithfully, in the port
  (`test_problems/200DMSO400Water/EPSR-test-ex49`).

Other details that `readpdb` dictates: only `ATOM` records are read (not `HETATM`);
the element is taken from columns 13–14 and compared with upper-case symbols, so
copper is written as `CU` and keeps that label in EPSR; the box comes from the
`CRYST1` record. `readpdb` finds the bonds from van der Waals radii
(`vanderWaalsRadii.txt`, copied from `EPSR-test`) times the factor given on standard
input (`readpdb.stdin`: 0.68, no dihedral groups). With 0.68 the bonds are exactly
those of the Dissolve species (192 C1–O, 96 C1–C2, 192 C2–C3, 96 C3–H) plus the 24
Cu–Cu pairs, and no spurious O···O bond inside a carboxylate (2.26 Å against
0.68 × 3.1 = 2.11 Å). The 192 Cu–O coordination bonds become intermolecular and are
described only by the Lennard-Jones and Coulomb terms of the reference potential
(and by the empirical potential). That is the main approximation of this setup.
The crystal positions are kept: there is no `randomise` step.

**Force field (`changeato.stdin`).** `readpdb` gives every atom type default
parameters; `changeato` replaces them with the Dissolve `PairPotentials` block
(ε in kJ/mol, σ in Å, charge in e): Cu 0.02092 / 3.144 / +1.098, O 0.401664 /
3.033 / −0.665, C1 0.39748 / 3.473 / +0.778, C2 −0.092, C3 −0.014 (same ε, σ as C1),
H 0.06276 / 2.846 / +0.109, and sets the H mass to 1.008 (natural hydrogen; `readpdb`
assumes deuterium). The temperature is 300 K and the Monte Carlo step sizes are
those of the EPSR test case. The Dissolve bond, angle, torsion and improper terms
are not reproduced: EPSR holds the bond lengths and bond angles found by `readpdb`
with its own harmonic restraints.

**Neutron weights (`CuBTC.NWTS.dat` → `nwts`).** One sample, natural isotopes,
no exchangeable hydrogen (the H atoms are bonded to carbon).

**EPSR input (`CuBTC.EPSR.inp`).** The silica input with these changes:

| keyword | value | why |
|---|---|---|
| title, `fnameato`, `fnamepcof` | `CuBTC.EPSR`, `CuBTC.ato`, `CuBTC.pcof` | |
| `qwidthqmax` | `0.2 5.0` | Bragg peaks up to 5 Å⁻¹ (Dissolve `QMax 5`), Gaussian width 0.2 Å⁻¹ (see below) |
| `hklqmin` | `0.0 10.0 0.0` | all reflections, 10 Å truncation radius for the Bragg g(r) (as `EPSR-test-ex120`), no Debye–Waller factor |
| `mplicities` | `1 1 1` | the box is one unit cell |
| `qmin` | `0.05 0.0` | fit the whole data range |
| `rshmin` | `0.7` | below the shortest distance, C–H 0.93 Å |
| `datafile`, `wtsfile` | `Empty-CuBTC.mint01`, `CuBTC.NWTStot.wts` | |

The Bragg width was chosen with 5-iteration trial runs on the reference potential.
Dissolve's `BraggQBroadening 0.055` makes EPSR's calculated F(Q) overshoot the data
badly at the lowest Q (−5.4 against −0.25 at Q = 0.3 Å⁻¹); a width of 0.2 Å⁻¹ gives
the best fit of the values tried:

| width (Å⁻¹) | 0.055 | 0.1 | 0.15 | 0.2 | 0.3 | no Bragg peaks |
|---|---|---|---|---|---|---|
| R-factor after 5 iterations | 0.075 | 0.085 | 0.037 | 0.031 | 0.032 | 0.039 |

`ereq` is set to 10 kJ/mol, the EPSR test-case value; the Dissolve example has no
EPSR module and so no `EReq`. The schedule (`steps.txt`) is that of the other
cases: a scratch start, 99 iterations on the reference potential, then 2 × 100
iterations with the empirical potential.

## Results

`run_validation.sh validation/cu-btc <fortran build> <work>`: the Fortran programs
and the Rust CLI on separate copies of the folder, every file compared after every
step.

| step | | Fortran | Rust | |
|---|---|---|---|---|
| 1 | `readpdb CuBTC.pdb` | 0.0 s | 0.0 s | ALL IDENTICAL |
| 2 | `changeato CuBTC.ato` | 0.0 s | 0.0 s | ALL IDENTICAL |
| 3 | `nwts CuBTC.NWTS.dat` | 0.0 s | 0.0 s | ALL IDENTICAL |
| 4 | `epsr` (scratch start, 1 iteration) | 0.5 s | 0.6 s | ALL IDENTICAL |
| 6 | `epsr` (99 iterations, reference potential) | 35.0 s | 41.9 s | ALL IDENTICAL |
| 9 | `epsr` (100 iterations, `ereq` 10) | 40.1 s | 78.4 s | ALL IDENTICAL |
| 10 | `epsr` (100 more iterations) | 40.1 s | 77.9 s | ALL IDENTICAL |

**ALL STEPS IDENTICAL**: 36 files after the last step, including the Bragg files
`.bragg` (the reflection list: 80 726 reflections at 940 unique Q values) and
`.inten` (the Bragg intensity of every atom pair on the Q grid up to 5 Å⁻¹). The Python wrapper, run on the same steps
by `python/tests/test_cu_btc.py`, also writes files byte-identical to the Fortran
after every EPSR run; `fortran_reference.tar.gz` (made by
`make_fortran_reference.py`) holds the Fortran folder after each of the four runs
for that test.

[gr_agreement.png](gr_agreement.png) overlays the partial g(r) of four atom pairs
from the three codes (made by `../plot_gr_agreement.py`); the differences are zero.

State after each run (from `CuBTC.EPSR.erg`, identical in all three codes):

| after iteration | energy / molecule (kJ/mol) | EPSR pressure | EP amplitude | R-factor |
|---|---|---|---|---|
| 1 (scratch start) | −413.0 | 138.8 | 0 | 3.08×10⁻² |
| 100 (reference potential only) | −528.7 | 17.9 | 0 | 2.29×10⁻² |
| 200 (`ereq` 10) | −538.4 | 8.6 | 10.0 | 1.69×10⁻² |
| 300 | −558.0 | 2.6 | 10.0 | 1.64×10⁻² |

The fit is much poorer than for the liquids (R-factor 1.6×10⁻² against 2×10⁻⁴),
as expected from the approximations above. The framework stays intact: the
molecules keep their geometry (Cu–Cu 2.64 Å, C1–O 1.26 Å, C3–H 0.93 Å), each Cu
keeps about four O neighbours (4.3 within 2.8 Å) and the paddlewheels stay 7.95 Å
apart, but the intermolecular Cu–O distance relaxes from 1.95 Å in the crystal to
2.25 Å because it is held only by Lennard-Jones and Coulomb forces. A better model
would need Cu–O bonds, which means a different split of the framework.

## Files

| file | purpose |
|---|---|
| `make_pdb.py` | converts `experimental_data/cu-btc/input.txt` to `CuBTC.pdb` |
| `CuBTC.pdb`, `readpdb.stdin`, `vanderWaalsRadii.txt` | input and answers of `readpdb` |
| `changeato.stdin` | the Dissolve force field, as `changeato` commands |
| `CuBTC.NWTS.dat`, `CuBTC.EPSR.inp` | neutron-weights and EPSR input files |
| `Empty-CuBTC.mint01` | the data (a copy of `experimental_data/cu-btc/data/`) |
| `steps.txt` | the calculation, step by step |
| `make_fortran_reference.py`, `fortran_reference.tar.gz` | the Fortran results after every EPSR run, for `python/tests/test_cu_btc.py` |
| `gr_agreement.png` | partial g(r) from the Fortran, the Rust CLI and the Python wrapper, with their differences |
