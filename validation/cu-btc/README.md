# Cu-BTC (HKUST-1): a crystalline metal-organic framework from a PDB file

The empty Cu-BTC framework measured on NIMROD, from the Dissolve example in
`../../experimental_data/cu-btc`: one cubic unit cell (a = 26.3336 Å, 624 atoms:
48 Cu, 192 O, 288 C, 96 H) at 300 K, fitted to `Empty-CuBTC.mint01`
(Q = 0.25 to 50 Å⁻¹) with Bragg peaks calculated up to 5 Å⁻¹.

## Run

```bash
# Fortran and Rust CLI side by side, every file compared after every step
validation/run_validation.sh validation/cu-btc FORTRAN_BUILD_DIR WORK_DIR

# the same steps through the Python wrapper, compared with the archived Fortran results
cd python && uv run pytest tests/test_cu_btc.py
```

`steps.txt` is the calculation: `readpdb`, `changeato`, `nwts`, then EPSR runs of
1, 99, 100 and 100 iterations with `ereq` set to 10 kJ/mol before the last two.
`fortran_reference.tar.gz` holds the Fortran folder after each EPSR run
(made by `make_fortran_reference.py FORTRAN_BUILD_DIR fortran_reference.tar.gz`).

## How the inputs were made

* **`make_pdb.py` → `CuBTC.pdb`.** EPSR cannot hold a periodic network as one
  molecule, so the Dissolve species is split at its Cu–O bonds into 32
  benzene-tricarboxylate linkers (18 atoms) and 24 Cu₂ paddlewheels (two Cu atoms
  2.63 Å apart). Single Cu residues are not possible: `readpdb` writes single-atom
  residues without their atom record, in the Fortran and in the port. `readpdb`
  reads only `ATOM` records and takes the element from columns 13–14, so copper is
  written as `CU` and keeps that label.
* **`readpdb` (`readpdb.stdin`: van der Waals factor 0.68, no dihedrals).** With
  0.68 the bonds found are exactly the Dissolve ones (192 C1–O, 96 C1–C2, 192 C2–C3,
  96 C3–H) plus the 24 Cu–Cu pairs. The 192 Cu–O coordination bonds become
  intermolecular, described only by Lennard-Jones and Coulomb terms: the main
  approximation of this setup. The crystal positions are kept (no `randomise`).
* **`changeato` (`changeato.stdin`).** The Dissolve `PairPotentials` block (ε, σ,
  charge per atom type), H mass 1.008, 300 K, the EPSR test-case Monte Carlo step
  sizes. The Dissolve bond, angle, torsion and improper terms are not reproduced;
  EPSR holds the bond lengths and angles found by `readpdb` with its own restraints.
* **`CuBTC.NWTS.dat`:** natural isotopes, no exchangeable hydrogen.
* **`CuBTC.EPSR.inp`:** the silica input with `qwidthqmax 0.2 5.0` (Bragg peaks to
  5 Å⁻¹, Gaussian width 0.2 Å⁻¹: Dissolve's 0.055 makes EPSR overshoot the data at
  low Q; 0.2 fitted best in 5-iteration trials), `hklqmin 0.0 10.0 0.0`,
  `qmin 0.05`, `rshmin 0.7` (below the 0.93 Å C–H bond).

## Results

`run_validation.sh`: **ALL STEPS IDENTICAL**, 36 files after the last step,
including the Bragg files `.bragg` (80 726 reflections at 940 unique Q values) and
`.inten`. The Python wrapper's files are byte-identical to the Fortran after every
EPSR run too. [gr_agreement.png](gr_agreement.png) overlays the partial g(r) of
four pairs from the three codes (`../plot_gr_agreement.py`).

| after iteration | energy / molecule (kJ/mol) | pressure | EP amplitude | R-factor |
|---|---|---|---|---|
| 1 (scratch start) | −413.0 | 138.8 | 0 | 3.08×10⁻² |
| 100 (reference potential only) | −528.7 | 17.9 | 0 | 2.29×10⁻² |
| 200 (`ereq` 10) | −538.4 | 8.6 | 10.0 | 1.69×10⁻² |
| 300 | −558.0 | 2.6 | 10.0 | 1.64×10⁻² |

The fit is far poorer than for the liquids (R-factor 1.6×10⁻² against 2×10⁻⁴), as
expected from the approximations above. The framework stays intact (Cu–Cu 2.64 Å,
C1–O 1.26 Å, C3–H 0.93 Å, about four O around each Cu, paddlewheels 8 Å apart),
but the intermolecular Cu–O distance relaxes from 1.95 Å to 2.25 Å. A better model
would need Cu–O bonds, which means a different split of the framework.

Timing per EPSR run (Fortran / Rust CLI): 0.5 / 0.6 s for 1 iteration, 35 / 42 s for
99 on the reference potential, 40 / 78 s for 100 with the empirical potential.
