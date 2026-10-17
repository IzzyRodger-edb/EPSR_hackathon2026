# Example: liquid water (H₂O / D₂O / HDO)

This example analyses an EPSR refinement of liquid water against three neutron
diffraction data sets measured with different hydrogen isotopes. Run it with:

```bash
cd visualisation/python_visualisation_toolkit
uv sync --extra viewer                               # once
uv run python example_H2O/plot_H2O.py                # ~1 min; figures -> example_H2O/figures/
uv run python example_H2O/plot_H2O.py --movies       # also the animated GIFs (~5 min)
uv run python example_H2O/plot_H2O.py --interactive  # also open the 3D viewer (iteration slider)
```

`--data` and `--out` change the input and output folders.

## The calculation

| | |
|---|---|
| system | 1000 water molecules (3000 atoms) in a cubic box of 31.07 Å, 0.1 atoms/Å³ (0.997 g/cm³), 300 K |
| model | rigid-ish SPC/Fw-like water: O–H 0.976 Å, H–O–H 107.1°, charges −0.82/+0.41 e |
| data | SANDALS (ISIS) neutron data, 3 isotopic samples: D₂O (`SLS18502`), HDO 50:50 (`SLS18500`), H₂O (`SLS18498`) |
| refinement | 4 `epsr` runs, 300 iterations in total (below) |
| origin | built from the Dissolve water example; validated byte-for-byte against the Fortran EPSR (`validation/water`) |

| run | iterations | stage | `ereq` |
|---|---|---|---|
| 1 | 1 | scratch start from random positions | 0 |
| 2 | 2–100 | equilibrate on the reference potential only | 0 |
| 3 | 101–200 | empirical potential (EP) on | 3 kJ/mol |
| 4 | 201–300 | continue, restarting from run 3's files | 3 kJ/mol |

## Data

| path | content |
|---|---|
| `data/trajectory/` | one snapshot per EPSR iteration (300): F(Q) data and simulation, partial g(r)/S(Q), potentials, positions, convergence numbers. Format: [toolkit README §5](../README.md#5-trajectory-folder-format). Recorded by `visualisation/recorder` |
| `data/epsr_output/` | the normal EPSR files after iteration 300: `waterbox.EPSR.*` result tables (`.u01`, `.g01`, ... see `epsrvis.outputs.TABLES`), `.erg` history, `.out`, final configuration `waterbox.ato`, potentials `waterbox.pcof`, weights `*.NWTStot.wts` |

## What the script does

1. Loads the trajectory and prints an overview and the fit at the end of each
   run: energy, pressure, EP amplitude, and R and R² per data set.
2. `make_report`: the standard figures `01`–`09` and `summary.csv/.md`
   (see the [toolkit README §6](../README.md#6-figures-make_report)).
3. Two water-specific figures, built from toolkit functions (a template for
   your own analyses):
   - `10_hydrogen_bonds.png`: hydrogen bonds per molecule (O···H < 2.45 Å
     between molecules) and the tetrahedral order q of the O atoms, against
     iteration; the distribution of q before and after the EP is switched on.
   - `11_isotope_fit.png`: R-factor of each isotopic sample relative to its
     starting value.
4. `12_final_epsr_output.png`: an overview read from the plain EPSR files in
   `data/epsr_output/`. It works for any EPSR folder, without a recording.
5. 3D views (pyvista, off-screen): `config_iter001/100/300.png` show a 10 Å
   slab through the box next to the partial g(r) of the same iteration.
   `--movies` adds `config_movie.gif` and `fq_fit.gif`.

## Results

| | iter 100 (reference only) | iter 200 (EP) | iter 300 (EP) |
|---|---|---|---|
| intermolecular energy (kJ/mol) | −34.98 | −46.23 | −44.89 |
| pressure (kbar) | 5.49 | 3.54 | 3.87 |
| EP amplitude / ereq (kJ/mol) | 0 / 0 | 3.0 / 3 | 3.0 / 3 |
| R-factor (mean) | 6.14×10⁻⁴ | 2.10×10⁻⁴ | 2.22×10⁻⁴ |
| R²: D₂O / HDO / H₂O | 0.848 / 0.770 / 0.930 | 0.976 / 0.855 / 0.946 | 0.973 / 0.846 / 0.945 |
| O–O first peak; O neighbours (r < 3.45 Å) | 2.82 Å; 5.56 | 2.79 Å; 5.08 | 2.79 Å; 5.12 |
| O···H first peak (H-bond); H per O (r < 2.4 Å) | 1.89 Å; 1.68 | 1.86 Å; 1.75 | 1.83 Å; 1.76 |
| H-bonds per molecule; mean tetrahedral order q | 3.49; 0.46 | 3.58; 0.54 | 3.60; 0.54 |

How to read it:

- **The refinement has converged.** After the EP is switched on (iteration 101),
  its amplitude reaches `ereq` = 3 kJ/mol within ~40 iterations. The R-factor,
  R², energy and first-shell numbers are then flat over the last ~150 iterations
  (`01_convergence`, `07_first_shell`).
- **The EP is what fits the D₂O data.** With the reference potential alone, R²
  for D₂O is 0.85; with the EP it is 0.97. The HDO sample stays the hardest to
  fit (R² ≈ 0.85).
- **The remaining misfit is at low Q** (Q < 2 Å⁻¹) in all three samples
  (`02_fq_fit`, `04_residual_map`). Inelasticity effects in hydrogen-containing
  samples are largest there, and EPSR cannot fit them.
- **The EP makes the water more structured.** The O–O first peak sharpens, the
  O···H hydrogen-bond minimum deepens (`06_potentials`, bottom row), and the
  tetrahedral order rises from 0.46 to 0.54 (`10_hydrogen_bonds`). The final
  ~3.6 hydrogen bonds per molecule and ~5 O neighbours are typical of liquid
  water models.
- **Iterations 1–20 are a random start with overlapping molecules.** Energies
  of ~10⁹ kJ/mol, R² < 0 and "8 hydrogen bonds per molecule" there are
  artefacts. The plots scale their axes to the later values.

## Outputs

| file | content |
|---|---|
| `01_convergence.png` | energies, pressure, R-factor (◇ = EPSR's own values), R², EP amplitude vs ereq, density, MC acceptance |
| `02_fq_fit.png` | data vs EPSR F(Q) + difference, three samples × end of each run |
| `03_fq_evolution.png` | simulated F(Q) through the refinement |
| `04_residual_map.png` | data − EPSR over Q and iteration |
| `05_partials.png` | O–O, O–H, H–H g(r) and S(Q) at the middle and end of each run |
| `06_potentials.png` | EP per pair; g(r) with the reference potential only vs with the EP |
| `07_first_shell.png` | first-peak position and height, coordination numbers vs iteration |
| `08_local_structure.png` | O neighbours, H-bonds per O, H–O–H angle, O···O···O angle |
| `09_displacement.png` | RMS displacement of the atoms since iteration 1 |
| `10_hydrogen_bonds.png` | H-bonds per molecule and tetrahedral order (example-specific) |
| `11_isotope_fit.png` | relative improvement of the fit per sample (example-specific) |
| `12_final_epsr_output.png` | overview from the EPSR output files |
| `config_iter001/100/300.png` | 3D configuration (10 Å slab) + partial g(r) |
| `summary.csv`, `summary.md` | numbers at the end of each run |
| `fq_fit.gif`, `config_movie.gif` | animations (`--movies`) |
