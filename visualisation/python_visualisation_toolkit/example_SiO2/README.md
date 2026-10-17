# Example: amorphous silica (SiO₂)

This example analyses an EPSR refinement of amorphous silica against one
neutron diffraction data set, and follows how the tetrahedral Si–O network
forms during the refinement. Run it with:

```bash
cd visualisation/python_visualisation_toolkit
uv sync --extra viewer                                 # once
uv run python example_SiO2/plot_SiO2.py                # ~1–2 min; figures -> example_SiO2/figures/
uv run python example_SiO2/plot_SiO2.py --movies       # also the animated GIFs (~5 min)
uv run python example_SiO2/plot_SiO2.py --interactive  # also open the 3D viewer (iteration slider)
```

`--data` and `--out` change the input and output folders.

## The calculation

| | |
|---|---|
| system | 1500 Si + 3000 O atoms (each a one-atom "molecule") in a cubic box of 40.77 Å, 0.0664 atoms/Å³ (2.21 g/cm³), 298 K |
| model | Lennard-Jones + charges: Si ε 0.175 kJ/mol, σ 1.03 Å, +2 e; O ε 0.165 kJ/mol, σ 3.5 Å, −1 e |
| data | NIMROD (ISIS) neutron total scattering, `SiO2_NIMROD.mint01`, fitted from Q = 0.1 Å⁻¹ |
| refinement | 4 `epsr` runs, 300 iterations in total (below) |
| origin | built from the Dissolve silica example; validated byte-for-byte against the Fortran EPSR (`validation/silica_bulk`) |

| run | iterations | stage | `ereq` |
|---|---|---|---|
| 1 | 1 | scratch start from random positions | 0 |
| 2 | 2–100 | equilibrate on the reference potential only | 0 |
| 3 | 101–200 | empirical potential (EP) on | 50 kJ/mol |
| 4 | 201–300 | continue, restarting from run 3's files | 50 kJ/mol |

## Data

| path | content |
|---|---|
| `data/trajectory/` | one snapshot per EPSR iteration (300): F(Q) data and simulation, partial g(r)/S(Q), potentials, positions, convergence numbers. Format: [toolkit README §5](../README.md#5-trajectory-folder-format). Recorded by `visualisation/recorder` |
| `data/epsr_output/` | the normal EPSR files after iteration 300: `silica.EPSR.*` result tables, `.erg` history, `.out`, final configuration `silica.ato`, potentials `silica.pcof`, weights `SiO2.NWTStot.wts` |

## What the script does

1. Loads the trajectory and prints an overview and the fit at the end of each
   run.
2. **Convergence check**: the linear trend of the R-factor, pressure, energy and
   EP amplitude over the last 50 iterations. It warns if the EP amplitude has
   not reached `ereq`.
3. `make_report`: the standard figures `01`–`09` and `summary.csv/.md`
   (see the [toolkit README §6](../README.md#6-figures-make_report)).
4. A silica-specific figure, built from toolkit functions (a template for your
   own analyses). `10_network.png` shows, against iteration: the fraction of Si
   with 3, 4 and 5 O neighbours (Si–O < 2.0 Å); the fraction of bridging (2 Si)
   and non-bridging (1 Si) O; and the mean tetrahedral order q of Si.
5. `11_final_epsr_output.png`: an overview read from the plain EPSR files in
   `data/epsr_output/`.
6. 3D views (pyvista, off-screen) of a 12 Å slab next to the partial g(r):
   `config_iter001.png` (random start) and `config_iter300.png` are coloured by
   element. `config_iter100_coordination.png` and
   `config_iter300_coordination.png` colour each Si by its number of O
   neighbours, so under-coordinated Si (defects) stand out. `--movies` adds
   `config_movie.gif` and `fq_fit.gif`.

## Results

| | iter 100 (reference only) | iter 200 (EP) | iter 300 (EP) |
|---|---|---|---|
| intermolecular energy (kJ/mol) | −876.2 | −862.3 | −845.3 |
| pressure (kbar) | 53.7 | 88.0 | 108.7 |
| EP amplitude / ereq (kJ/mol) | 0 / 0 | 15.1 / 50 | 23.7 / 50 |
| R-factor | 1.04×10⁻³ | 5.23×10⁻⁴ | 2.30×10⁻⁴ |
| R² | 0.631 | 0.815 | 0.919 |
| Si–O first peak; O per Si (r < 2.22 Å) | 1.56 Å; 3.69 | 1.59 Å; 3.82 | 1.59 Å; 3.88 |
| 4-fold / 3-fold Si (r < 2.0 Å) | 60.1 % / 39.4 % | 74.7 % / 25.1 % | 80.7 % / 18.9 % |
| bridging O | 79.0 % | 85.9 % | 88.8 % |
| mean tetrahedral order q of Si | 0.838 | 0.875 | 0.881 |
| most probable O–Si–O / Si–O–Si angle | 107° / 141° | 107° / 155° | 107° / 155° |

How to read it:

- **The refinement is still improving; it has not converged.** The EP amplitude
  has reached only 23.7 of the requested 50 kJ/mol and is still rising (+6 per
  100 iterations). The pressure keeps increasing (+18 kbar per 100 iterations).
  R² is still going up (`01_convergence`). The script's convergence check flags
  this. This schedule is enough to compare the Fortran and Rust codes, but not a
  production refinement. Run more iterations before interpreting the structure
  quantitatively.
- **The EP builds the network.** The reference potential alone leaves 39 % of
  Si only 3-coordinated. With the EP the fraction of 4-fold Si rises steadily
  (60 → 81 %), and so do the bridging-O fraction (79 → 89 %) and the
  tetrahedral order (`10_network`, `08_local_structure`). The 3D views show the
  3-coordinated Si (teal) being replaced by 4-coordinated Si (green).
- **The geometry is that of a tetrahedral network.** Si–O peaks at 1.59 Å
  (`05_partials`), the O–Si–O angle at ≈107° (close to the tetrahedral 109.5°),
  and the most probable Si–O–Si angle shifts from ≈141° to ≈155° as the EP acts.
- **Averaging matters.** The R-factor of a single configuration (dotted line
  in `01_convergence`) fluctuates far more than the running average EPSR fits.
- **Iterations 1–20 are a random start with overlapping atoms.** Energies of
  ~10⁹ kJ/mol and R² < 0 there are artefacts.

## Outputs

| file | content |
|---|---|
| `01_convergence.png` | energies, pressure, R-factor (◇ = EPSR's own values), R², EP amplitude vs ereq, density, MC acceptance |
| `02_fq_fit.png` | data vs EPSR F(Q) + difference at the end of each run |
| `03_fq_evolution.png` | simulated F(Q) through the refinement |
| `04_residual_map.png` | data − EPSR over Q and iteration |
| `05_partials.png` | Si–Si, Si–O, O–O g(r) and S(Q) at the middle and end of each run |
| `06_potentials.png` | EP per pair; g(r) with the reference potential only vs with the EP |
| `07_first_shell.png` | first-peak position and height, coordination numbers vs iteration |
| `08_local_structure.png` | O around Si, Si around O, O–Si–O and Si–O–Si angle distributions |
| `09_displacement.png` | RMS displacement of the atoms since iteration 1 |
| `10_network.png` | Si coordination, bridging O, tetrahedral order vs iteration (example-specific) |
| `11_final_epsr_output.png` | overview from the EPSR output files |
| `config_iter*.png` | 3D configuration (12 Å slab) + partial g(r); `_coordination` = Si coloured by O neighbours |
| `summary.csv`, `summary.md` | numbers at the end of each run |
| `fq_fit.gif`, `config_movie.gif` | animations (`--movies`) |
