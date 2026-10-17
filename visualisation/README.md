# EPSR visualisation and analysis tools

Tools for following an EPSR refinement **iteration by iteration**:

* how well the simulated structure factors F(Q) fit each measured data set;
* the key convergence numbers: energy, pressure, R-factor, R², χ², EP amplitude,
  density, Monte Carlo acceptance;
* the partial g(r) / S(Q), the empirical potential, and first-shell structure
  (peak positions, coordination numbers);
* the atomic configuration in 3D (pyvista), next to the partial g(r) of the same
  iteration.

They work on the two validated cases in [`../validation`](../validation): water and
amorphous silica.

> **For the Python wrapper, use [`python_visualisation_toolkit/`](python_visualisation_toolkit/README.md).**
> It is the integration-ready version of these tools (package `epsrvis`). It adds
> a recorder that captures a running wrapper `Simulation` iteration by iteration,
> an in-memory `Trajectory`, tests, and documented H₂O and SiO₂ examples with
> their data.

## Why a recorder is needed

EPSR only writes its result files (`.u01`, `.g01`, ...) after the **last** iteration of
each `epsr` call, and overwrites them each time. It computes the R-factor only
then, too: in the `.erg` history the R-factor column just repeats the last value.
So the files alone cannot show how the fit develops.

[`recorder/`](recorder) is a small standalone Rust crate. Like
`validation/api_driver`, it runs a case's `steps.txt` through the `epsr-core`
library, with one difference. Instead of `sim.run()` it calls
`sim.iterate(1..=niter)` and then `sim.write_fort35()`, which is exactly what
`run()` does. Between iterations it copies the in-memory results to disk through
the public accessors. Those accessors only read the simulation, so the EPSR files
are unaffected. This was checked: all of them match the validated Rust run
(see *Results*). `epsr-core` is not modified.

Recorded per iteration (in `data/<case>/trajectory/`):

| file | shape per iteration | content |
|------|---------------------|---------|
| `history.tsv` | one line | EPSR `Summary` (energies, pressure, EP amplitude, ...), `ereq`, box, density, MC step sizes and move counters |
| `fq_data.f32`, `fq_err.f32` | [data set, Q] | measured F(Q) on the simulation Q grid, and its error |
| `fq_sim.f32` | [data set, Q] | simulated F(Q), running average — what EPSR fits |
| `fq_inst.f32` | [data set, Q] | simulated F(Q) of the current configuration only |
| `gr_inter.f32`, `gr_intra.f32` | [pair, r] | partial g(r), inter- and intramolecular |
| `sq_inter.f32`, `sq_intra.f32` | [pair, Q] | partial S(Q) |
| `ep.f32`, `ref_pot.f32` | [pair, r] | empirical and reference potential (kJ/mol) |
| `pos.f32` | [atom, 3] | atom positions (Å); box centred on the origin, molecules kept whole |
| `meta.json` | — | grids, labels, fitted Q range per data set, molecule of each atom |

`data/<case>/epsr_output/` holds the normal EPSR files after the last run: `.EPSR.*`,
`.ato`, `.pcof`, weights.

## Setup

Python packages are managed with [uv](https://docs.astral.sh/uv/) (Python 3.12,
numpy, scipy, matplotlib, pyvista, imageio):

```bash
cd visualisation
uv sync
```

To record a case again (≈10 min for water, ≈30 min for silica, single core):

```bash
visualisation/collect_data.sh water ~/epsr-work/rec-water ~/epsr-work/val-water/rust
visualisation/collect_data.sh silica_bulk ~/epsr-work/rec-silica ~/epsr-work/val-silica/rust
```

The third argument is optional. When given, the recorder's EPSR files are compared
with that validated folder using `rust/tools/compare_outputs.py`.

## Usage

All commands are run from `visualisation/`.

```bash
# every diagnostic figure + summary table (add --gif for an animated F(Q)/g(r) GIF)
uv run epsr-vis report data/water/trajectory --out figures/water --gif

# interactive: 3D configuration + partial g(r), slider over EPSR iterations
uv run epsr-vis view data/water/trajectory
uv run epsr-vis view data/silica_bulk/trajectory --color-by coordination --slab 12

# the same view, off-screen, as an image or a movie
uv run epsr-vis snapshot data/water/trajectory --step 100 --slab 10 --out figures/water/config_100.png
uv run epsr-vis movie data/silica_bulk/trajectory --every 5 --out figures/silica_bulk/config.gif

# overview of any finished EPSR folder from its own output files (no recorder needed)
uv run epsr-vis final data/water/epsr_output --out figures/water/final_outputs.png
```

Viewer options: `--gr inter|total` (intermolecular g(r) as in `.g01`, or with the
intramolecular part added), `--slab W` (only atoms with |z| < W/2, to see inside the
box), `--color-by coordination` (colour the centre atoms of the case's first
coordination analysis by neighbour count, e.g. Si by number of O), `--no-bonds`,
`--atom-scale`.

From Python:

```python
from epsr_vis import Trajectory
from epsr_vis import metrics, plots
from epsr_vis.viewer import ConfigurationViewer

t = Trajectory.load("data/water/trajectory")
fit = metrics.fit_metrics(t)          # rfactor, rfactor_inst, r2, chi2: [step, dataset]
t.fq_sim[-1, 0]                       # simulated F(Q), last iteration, first data set
t.gr_total(99)                        # partial g(r) at iteration 100, [pair, r]
plots.convergence(t).savefig("c.png")
ConfigurationViewer(t).show()
```

## Figures (`epsr-vis report`)

| figure | what it shows | what to look for |
|--------|---------------|------------------|
| `01_convergence` | energies, pressure, R-factor per data set, R², EP amplitude vs `ereq`, density, MC acceptance | plateaus in R and energy; the EP amplitude reaching `ereq`; acceptance staying in a sensible range |
| `02_fq_fit` | data (black) vs EPSR (colour) and the difference, per data set, at the end of each run | where in Q the misfit is; how much each stage improves it |
| `03_fq_evolution` | simulated F(Q) through the refinement (light → dark) | how quickly each feature is captured |
| `04_residual_map` | data − EPSR over Q and iteration | which Q regions the EP fixes, and which stay wrong (e.g. low-Q inelasticity) |
| `05_partials` | partial g(r) and S(Q) at the middle and end of each run | how the EP changes the local structure |
| `06_potentials` | empirical potential through the refinement; reference + EP at the end | whether the EP stays small next to the reference potential |
| `07_first_shell` | first-peak position, height, coordination number vs iteration | structural convergence: whether the local structure has settled |
| `08_local_structure` | coordination-number and bond-angle distributions from the coordinates | e.g. H-bonds per water, tetrahedrality; Si coordination defects |
| `09_displacement` | RMS atom displacement since the first iteration | how much the configuration still moves (whether it decorrelates) |
| `summary.md/.csv` | key numbers at the end of each run | — |

### Metric definitions

* **R-factor**: mean over the fitted Q points of (data − simulation)², with the
  running-average simulation. This is EPSR's own definition (`rfacd` in `.out`).
  The recomputed value matches EPSR's at the end of every run to 7 significant
  figures; the dashboard draws EPSR's values as ◇ markers.
  *Current configuration*: the same for one configuration only (EPSR's "raw"
  R-factor). The gap between the two shows the configurational noise.
* **R²** = 1 − Σ(data − sim)² / Σ(data − mean data)² over the fitted Q range.
* **χ²** = mean((data − sim)/σ)², with the data errors σ. Treat it as relative:
  σ from Gudrun is statistical only, so χ² is much larger than 1 everywhere.
* **Coordination number** N_ab(r) = ρ c_b ∫ 4πr² g_ab(r) dr (as `.z01`). The
  first-shell values integrate to the lowest point of the smoothed g(r) within
  1.2 Å after the first peak.
* **Acceptance**: from EPSR's MOLFIT counters, for the last MC cycle of each
  iteration. Rotation and translation of a molecule are one combined move.
* **Density**: constant in these runs (constant-volume MC). It is shown so that runs
  with pressure/box control (`control-p`) are covered too.

## Results on the validated cases

**The recording does not change the simulation.** For both cases every EPSR file
the recorder wrote is identical to the validated Rust run
(`compare_outputs.py`: ALL IDENTICAL). The R-factor recomputed at every iteration
equals EPSR's `rfacd` wherever EPSR computes it (the end of each run).

| | water (iter 100 → 200 → 300) | silica (iter 100 → 200 → 300) |
|---|---|---|
| stage | reference only → EP (ereq 3) → EP | reference only → EP (ereq 50) → EP |
| R-factor (mean of data sets) | 6.1×10⁻⁴ → 2.1×10⁻⁴ → 2.2×10⁻⁴ | 1.0×10⁻³ → 5.2×10⁻⁴ → 2.3×10⁻⁴ |
| R² | D₂O 0.85 → 0.98 → 0.97; HDO 0.77 → 0.85 → 0.85; H₂O 0.93 → 0.95 → 0.95 | 0.63 → 0.81 → 0.92 |
| energy (kJ/mol per molecule) | −35.0 → −46.2 → −44.9 | −876 → −862 → −845 |
| pressure (kbar) | 5.5 → 3.5 → 3.9 | 54 → 88 → 109 |
| EP amplitude / ereq (kJ/mol) | 0 → 3.0/3 → 3.0/3 | 0 → 15.1/50 → 23.7/50 |
| density | 0.997 g/cm³ (0.1 atoms/Å³) | 2.21 g/cm³ (0.0664 atoms/Å³) |
| first shell (final) | O–O 2.79 Å, 5.1 O within 3.45 Å; O···H 1.83 Å, 1.76 H-bonded H per O | Si–O 1.59 Å, 3.88 O per Si; O–O 2.55 Å; Si–Si 3.09 Å |
| from coordinates | ~50 % of O have 2 H-bonds; H–O–H ≈ 107° | 4-coordinated Si: 60 % → 75 % → 81 %; O–Si–O peak ≈ 108° |

What the figures show:

* **Water converges.** The EP reaches `ereq` = 3 within ~40 iterations of being
  switched on. R and R² are flat over the last 150 iterations. The first-shell
  numbers settle. Most of the remaining misfit is at Q < 2 Å⁻¹, in all three
  samples (`04_residual_map`). That is the region where inelasticity effects in
  hydrogenous data are largest, and EPSR cannot fit them. The HDO set is fitted worst
  (R² 0.85). The EP sharpens the O–O peak and deepens the O···H minimum
  (`06_potentials`, bottom row).
* **Silica has not converged.** The EP amplitude is still rising (23.7 of 50
  kJ/mol), the pressure keeps increasing and R² is still improving. The fraction of
  4-coordinated Si also keeps rising (`08_local_structure`; in 3D:
  `config_iter*_coordination.png`). The refinement needs more iterations.
* The gap between the running-average and the current-configuration R-factor
  (dotted line) is small for water but large for silica. A single silica
  configuration is noisy: the averaging matters.

## Layout

```
visualisation/
├── README.md
├── pyproject.toml, uv.lock     Python project (uv)
├── collect_data.sh             record a case and copy its data here
├── recorder/                   Rust: per-iteration recorder (uses epsr-core, does not modify it)
├── src/epsr_vis/
│   ├── io.py                   Trajectory loader; readers for .?01 tables and .erg
│   ├── metrics.py              R-factor, R², χ², acceptance, coordination, angles, MSD
│   ├── plots.py                matplotlib figures
│   ├── viewer.py               pyvista 3D configuration + linked g(r) chart
│   ├── cases.py                per-system settings (colours, bonds, analyses)
│   ├── style.py                palette and plot style
│   └── cli.py                  `epsr-vis` command
├── data/<case>/                trajectory/ (recorded), epsr_output/ (final EPSR files)
└── figures/<case>/             output of `epsr-vis report` etc.
```

To support another system, add an entry to `CASES` in `cases.py`: element of each
atom type, bonds to draw, coordination and angle analyses, and stage names.
