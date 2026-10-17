# epsrvis — EPSR analysis and visualisation toolkit

`epsrvis` analyses and plots EPSR refinements run with the Rust port of EPSR
(`rust/`). It follows a refinement **iteration by iteration**:

- how well the simulated F(Q) fits each data set;
- the convergence numbers: energy, pressure, R-factor, R², χ², EP amplitude,
  density, Monte Carlo acceptance;
- partial g(r) and S(Q), the empirical potential (EP), and coordination
  numbers and bond angles;
- the atomic configuration in 3D (pyvista), next to the partial g(r) of the
  same iteration.

It is written to be **dropped into the Python wrapper** of `epsr-core` (e.g. as
`epsr.vis`): it uses relative imports only and depends only on
numpy/scipy/matplotlib. pyvista is optional and only needed for the 3D viewer.
Until the wrapper exists, it works on trajectories recorded by the Rust
recorder in [`../recorder`](../recorder).

Two worked examples come with data:

| folder | system | data |
|---|---|---|
| [example_H2O/](example_H2O/README.md) | 1000 water molecules | H₂O, D₂O, HDO neutron data (SANDALS) |
| [example_SiO2/](example_SiO2/README.md) | 1500 Si + 3000 O, amorphous silica | one neutron data set (NIMROD) |

---

## 1. Installation

The project is managed with [uv](https://docs.astral.sh/uv/):

```bash
cd visualisation/python_visualisation_toolkit
uv sync --extra viewer          # numpy, scipy, matplotlib + pyvista, imageio (3D viewer), pytest
uv run pytest                   # 7 tests, ~6 s
```

Without the 3D viewer: `uv sync` (core dependencies only).

## 2. Quick start

```python
from epsrvis import Trajectory, make_report, metrics, plots

traj = Trajectory.load("example_H2O/data/trajectory")    # 300 iterations, memory-mapped
make_report(traj, "figures/")                            # 9 figures + summary.csv/.md

fit = metrics.fit_metrics(traj)                          # per iteration, per data set
fit["rfactor"][-1], fit["r2"][-1]                        # final R-factor and R² of each data set
traj.fq_sim[99, 0]                                       # simulated F(Q), iteration 100, first data set
traj.gr_total(99)                                        # partial g(r) (inter + intra), [pair, r]

fig = plots.convergence(traj)                            # any single figure
fig.savefig("convergence.png")

from epsrvis.viewer import ConfigurationViewer           # needs pyvista
ConfigurationViewer(traj, slab=10).show()                # 3D view + g(r), slider over iterations
```

The same from the command line:

```bash
uv run epsrvis report  example_H2O/data/trajectory --out figures/ --gif
uv run epsrvis view    example_SiO2/data/trajectory --color-by coordination --slab 12
uv run epsrvis snapshot example_H2O/data/trajectory --step 100 --out config.png
uv run epsrvis movie   example_H2O/data/trajectory --every 5 --out config.gif
uv run epsrvis final   example_H2O/data/epsr_output --out final.png   # plain EPSR output files
```

## 3. Package layout

| module | contents |
|---|---|
| `trajectory.py` | `Trajectory`: the per-iteration record; `load()`/`save()` of the folder format (§5) |
| `recorder.py` | `TrajectoryRecorder`, `record_run()`: record a running wrapper `Simulation` (§4) |
| `metrics.py` | `fit_metrics` (R-factor, R², χ²), `acceptance`, `mass_density`, `running_coordination`, `first_shell`, `shell_history`, `coordination_distribution`, `bond_angle_distribution`, `tetrahedral_order`, `mean_squared_displacement`, `neighbour_pairs` |
| `plots.py` | one matplotlib function per figure (table in §6), `fq_animation`, `final_outputs`, `summary_rows` |
| `report.py` | `make_report()`: every figure + `summary.csv`/`summary.md`; `FIGURES` lists them |
| `viewer.py` | `ConfigurationViewer`: pyvista 3D configuration + linked g(r) chart; `show`, `screenshot`, `movie` |
| `live.py` | `LiveMonitor`: redraw R-factor and F(Q) while the wrapper runs |
| `outputs.py` | `read_epsr_table` (`.?01` files), `read_erg`, `OutputFolder`, `TABLES` (what each file holds) |
| `systems.py` | per-system settings: elements, colours, bonds, coordination/angle analyses, stage names |
| `style.py` | colour palette and matplotlib style |
| `cli.py` | the `epsrvis` command |

## 4. Integration with the Python wrapper

### 4.1 Why a recorder

EPSR writes its result files (`.u01`, `.g01`, ...) only on the **last**
iteration of each `epsr` call, and overwrites them each time. It also computes
the R-factor only then. To follow a refinement, the wrapper must run one
iteration at a time and copy its in-memory results after each one. `record_run`
does exactly what `Epsr::run` does (`iterate(1..=niter)`, then `write_fort35()`),
so the simulation folder ends up the same as after `sim.run()`. The Rust recorder
does the same, and its files were checked to be identical to the validated runs.

```python
import epsr                                     # the PyO3 wrapper (to be written)
from epsrvis import TrajectoryRecorder, record_run, make_report
from epsrvis.live import LiveMonitor

rec = TrajectoryRecorder(case="water")          # case selects the settings in systems.py
for run in range(4):
    # edit .EPSR.inp / .pcof between runs as usual (niter, ereq, ...)
    sim = epsr.Simulation("work/water", "waterbox.EPSR.inp")             # Epsr::setup
    record_run(sim, rec, pcof="work/water/waterbox.pcof",                # == sim.run(), recorded
               callback=LiveMonitor(every=10))                           # optional live plot
traj = rec.trajectory()                         # in memory: analyse straight away
rec.save("work/water/trajectory")               # or keep it (same format as the Rust recorder)
make_report(traj, "work/water/figures")
```

`ereq` is stored in the `.pcof`, not in the simulation, so pass `pcof=` (or
`ereq=`) to label the stages.

### 4.2 What the wrapper must expose

The recorder is duck-typed, because the wrapper's API is not fixed yet. Each
item is looked up under the name in `rust/docs/PYTHON_WRAPPER_NOTES.md` §5 *or*
under the Rust accessor name. It may be a property or a method, and per-pair or
per-data-set quantities may be given by index (`gr_inter(ip)`) or as one 2D
array (`gr_inter()`).

| needed | accepted names | Rust source (`epsr-core`) |
|---|---|---|
| one iteration | `iterate(i)` | `Epsr::iterate` |
| iterations per run | `niter` | `sim.ctl.niter` |
| grids | `q` / `q_grid`, `r` / `r_grid` | `q_grid()`, `r_grid()` |
| data sets | `n_datasets`, `data_q_range(id)`, `data_fq(id)`, `simulated_fq(id)` | `api.rs` |
| partials | `gr_inter`, `gr_intra`, `sq_inter`, `sq_intra` | `api.rs` |
| potentials | `empirical_potential`, `reference_potential` | `api.rs` |
| labels | `pair_labels`, `atom_types`, `atom_labels`, `molecule_atoms` | `api.rs` |
| configuration | `positions` / `atom_positions`, `box_size` | `api.rs` |
| convergence | `summary` (dict or object) | `Summary`; keys `energd`/`energy`, `pressd`/`pressure`, `energabs`/`ep_amplitude`, `rfactot`/`r_factor`, `rfacd`/`r_factor_per_dataset`, `foldsv`, `foldsb`, `chisq`, `sentb`, `epcup`, `epccoeff` (`recorder.SUMMARY_ALIASES`) |
| *optional:* finish a run | `write_fort35()` | `Epsr::write_fort35` |
| *optional:* data set names | `get_parameter("datafile", id)` | `api.rs` |
| *optional:* data errors | `data_errors(id)` | `sim.errs` (no accessor yet) |
| *optional:* F(Q) of the current configuration | `simulated_fq_current(id)` | same sum as `simulated_fq` on `sim.qwind` (see `../recorder/src/main.rs`) |
| *optional:* MC counters | `mc_counters()` → 2×4 (attempts, rejects × intra, internal rotation, rotation, translation) | `st.mol.nmovei/ninter/nwhole/nmovem`, as in the `.out` file |
| *optional:* step sizes | `step_sizes()` → dict `sizefactor, stepm, stepr, stepmi, stepri` | `st.mol` |

Missing optional items are recorded as NaN or left out. The plots then say "not
recorded" instead of failing; `test_recorder_without_optional_quantities`
checks this. The three optional items without an accessor (data errors,
current-configuration F(Q), MC counters) are worth adding to `api.rs`.
`../recorder/src/main.rs` shows how each is computed from the public state.

`tests/conftest.py::ReplaySimulation` is a working mock of the wrapper. It
replays a recorded run through the API in this table, so the whole path
(wrapper → recorder → disk → analysis) is tested before the wrapper exists.

### 4.3 Embedding as a sub-package

Copy `src/epsrvis/` into the wrapper package (e.g. `python/epsr/vis/`). Then add
numpy/scipy/matplotlib to its dependencies, and pyvista/imageio as an optional
extra. Nothing imports pyvista unless `epsrvis.viewer` is imported.

## 5. Trajectory folder format

Written by `Trajectory.save`, `TrajectoryRecorder.save` and the Rust recorder;
read by `Trajectory.load`, which memory-maps the arrays.

| file | content |
|---|---|
| `meta.json` | `case`, sizes (`nq`, `nr`, `natom`, `nmol`, `npairs`, `ndata`), grids `q`, `r`, `atom_types`, `pairs`, `datasets`, `q_range` (fitted Q indices per data set), `atom_labels`, `molecule_of_atom`, `arrays` (shape of one step of each array) |
| `history.tsv` | one line per iteration; columns in `trajectory.HISTORY_COLUMNS` |
| `<name>.f32` | little-endian float32, one block per iteration; names and meaning in `trajectory.ARRAYS` |
| `pos_initial.f32` | configuration before the first iteration (optional) |

Per-iteration arrays: `fq_data`, `fq_err`, `fq_sim`, `fq_inst` [data set, Q];
`gr_inter`, `gr_intra`, `ep`, `ref_pot` [pair, r]; `sq_inter`, `sq_intra` [pair, Q];
`pos` [atom, 3] in Å, box centred on the origin, molecules kept whole.

## 6. Figures (`make_report`)

| file | content |
|---|---|
| `01_convergence` | energies, pressure, R-factor per data set (◇ = EPSR's own values), R², EP amplitude vs `ereq`, density, MC acceptance |
| `02_fq_fit` | data vs EPSR F(Q), and their difference, at the end of each run |
| `03_fq_evolution` | simulated F(Q) through the refinement (light → dark) |
| `04_residual_map` | data − EPSR over Q (x) and iteration (y) |
| `05_partials` | partial g(r) and S(Q) at the middle and end of each run |
| `06_potentials` | EP per pair through the refinement; g(r) before vs after the EP is switched on |
| `07_first_shell` | first-peak position and height, coordination number vs iteration |
| `08_local_structure` | coordination-number and bond-angle distributions from the coordinates |
| `09_displacement` | RMS atom displacement since the first iteration |
| `summary.csv/.md` | key numbers at the end of each run |
| `fq_fit.gif` | (`gif=True`) animated F(Q) fit and partial g(r) |

## 7. Metric definitions

- **R-factor**: mean over the fitted Q points of (data − simulation)², with the
  running-average simulation. This is EPSR's `rfacd`; the test suite checks that
  the two agree at every run end. *Current configuration*: the same, for one
  configuration only (EPSR's raw R-factor).
- **R²** = 1 − Σ(data − sim)² / Σ(data − mean)² over the fitted Q range.
- **χ²** = mean(((data − sim)/σ)²). σ holds statistical errors only, so treat χ²
  as relative.
- **Coordination number** N_ab(r) = ρ c_b ∫ 4πr² g_ab(r) dr (EPSR's `.z01`
  convention: intermolecular g(r)). `shell_history` integrates every iteration
  to the *final* first minimum, so a change means the structure changed, not
  the minimum.
- **Tetrahedral order** q = 1 − 3/8 Σ_{j<k} (cos ψ_jk + 1/3)² over the 4 nearest
  neighbours: 1 for a perfect tetrahedron, 0 for a random arrangement.
- **Acceptance**: from EPSR's MOLFIT counters in the last MC cycle of each
  iteration. Rotation and translation of a molecule are one combined move.

## 8. Adding a system

Add an entry to `CASES` in `systems.py`, keyed by the `case` name given to the
recorder. It holds the element of each atom type (colours, radii, masses), the
bonds to draw, the coordination and angle analyses, and the stage names. An
unknown case still works: elements are guessed from the atom-type labels, and
the structural analyses are skipped.

## 9. Relation to `visualisation/`

`visualisation/src/epsr_vis` is the first version of these tools. This toolkit
adds the recorder for the wrapper, the in-memory `Trajectory`, tests and the
examples. The example data are copies of `visualisation/data/` (recorded by
`visualisation/recorder`, built from the cases in `validation/`).
