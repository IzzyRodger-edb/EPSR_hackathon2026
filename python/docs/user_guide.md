# User guide

This guide walks through the EPSR workflow in Python. The examples use the
DMSO/water test case in `EPSR-test/` (200 D6-DMSO + 400 water molecules, three
neutron datasets). Atom types `S`, `C`, `O`, `H` belong to DMSO and `O1`, `H1`
to water.

## 1. How EPSR works, in one paragraph

EPSR builds an atomistic model of a liquid or amorphous material that agrees with
measured diffraction data. A Monte Carlo simulation runs with a *reference
potential* (Lennard-Jones plus Coulomb, defined in the `.ato` file). After every
`ntimes` Monte Carlo cycles the simulated structure factors are compared with the
data, and the difference is used to refine an *empirical potential* (EP) that
pushes the simulation towards the data. The EP is limited in size by the
parameter `ereq`. Repeating this many times gives a model whose partial pair
distribution functions g(r) can be analysed.

A simulation folder holds:

| file | contents |
|---|---|
| `<name>.EPSR.inp` | run parameters (`ntimes`, `niter`, data sets, ...) |
| `<name>.ato` | the configuration, atom types, reference-potential parameters, random-number state |
| `<name>.pcof` | EP control parameters (`ereq`, ...) and the EP coefficients of every atom pair |
| `*.wts`, data files | neutron or X-ray weights and the measured data |
| `<name>.EPSR.*` | results written by each run (`.g01` partial g(r), `.u01` simulated F(Q), `.erg` history, ...) |

## 2. A first run

```python
import epsr

proj = epsr.Project.from_template("EPSR-test", "work/dmso")
sim = proj.simulation()
summary = sim.run()
print(summary.r_factor, summary.energy, summary.pressure)
```

`Project.from_template` copies the regular files of a folder into a new one
(it refuses to overwrite a non-empty folder unless `overwrite=True`). **Always run
on a copy**: EPSR rewrites `.ato`, `.pcof`, `.EPSR.inp` and `.EPSR.inpa` at the
end of each run and appends a line to `.erg`, and creating a `Simulation` already
rewrites `.qdr`.

`proj.simulation()` reads every input (the Rust `Epsr::setup`). `run()` performs
`niter` iterations and writes the results. It returns a `Summary`:

| field | meaning |
|---|---|
| `r_factor`, `r_factor_per_dataset` | overall and per-dataset R-factor (fit quality) |
| `energy`, `pressure` | intermolecular energy per molecule (kJ/mol) and pressure (kbar) |
| `energy_inter`, `energy_intra` | running intermolecular and intramolecular energies per molecule |
| `ep_amplitude` | absolute energy of the empirical potential, at most `ereq` |
| `ep_coefficient`, `gr_roughness`, `chisq`, ... | the other quantities EPSR writes to `.out` |

`summary.fortran` holds the same values under the Fortran variable names.

### Iterating step by step

```python
sim = proj.simulation(overrides={"niter": 3})
for s in sim.iterations():          # yields a Summary after each iteration
    print(s.iteration, s.r_factor)
```

`iterations()` writes exactly the same files as `run()`. Iterations must run in
order and only up to `niter`; a finished `Simulation` cannot be run again. Create
a new one for the next run.

## 3. Inspecting results

All arrays are float32 numpy copies.

```python
q, r = sim.q, sim.r                 # Q grid (1/Å), r grid (Å)
sim.pair_labels                     # [('S', 'S'), ('S', 'C'), ..., ('H1', 'H1')]
sim.gr_inter()                      # (n_pairs, nr) intermolecular partial g(r)
sim.gr_inter("O1", "H1")            # one pair; label order does not matter
sim.gr_intra(3), sim.sq_inter(), sim.sq_intra()
sim.current_gr()                    # (2, n_pairs, nr) histograms of the current configuration only
```

The partials are running averages over the last `inter` iterations, exactly
what EPSR writes to `.g01`, `.y01`, `.f01` and `.s01`.

Data and fits, per dataset (0-based index, in the order of the `data n` blocks):

```python
for i, ds in enumerate(sim.datasets):
    f = sim.fit(i)                  # restricted to the fitted Q range
    print(ds.name, f.q[0], f.q[-1], abs(f.residual).max())
sim.data_fq(0), sim.simulated_fq(0), sim.data_errors(0), sim.data_q_range(0)
```

Configuration:

```python
sim.positions                       # (n_atoms, 3) Å
sim.atom_labels, sim.molecule_atoms # per-atom labels; atom index range of each molecule
sim.box_size, sim.density           # Å, atoms/Å³
sim.get_parameter("ntimes")         # current value of any input keyword, as text
sim.get_parameter("tweak", 1)       # a per-dataset keyword of dataset 1 (0-based)
```

Result files can also be read directly, without a `Simulation`:

```python
g = proj.read_table("g01")          # epsr.io.Table
g.x, g["O1-H1"], g.errors["O1-H1"]
g.to_dataframe()                    # with the pandas extra
proj.history()                      # the whole .erg file as a History
```

## 4. Changing parameters: overrides

Overrides change input keywords after the files are read, as if the lines had
been appended to the `.EPSR.inp` file. Keywords of the `.pcof` header (`ereq`,
`rminpt`, `Gaussian`, ...) work the same way.

```python
sim = proj.simulation(overrides={
    "ntimes": 10,
    "ereq": 15.0,                       # EP amplitude limit
    "potfac": [1.0, 1.0],               # lists become space-separated values
    "data": {2: {"tweak": 1.05}},       # dataset blocks use EPSR's 1-based numbers
})
sim = proj.simulation(overrides=["ntimes 10", "data 2", "tweak 1.05"])   # raw lines also work
```

Note that `data` numbers in overrides are 1-based like the input file, while
every other dataset index in the API is 0-based.

Overrides are written back to the input files at the end of the run, so they
persist. An unknown keyword raises `EpsrInputError`.

## 5. Refining: repeated runs

A real refinement is hundreds of EPSR runs. `epsr.refine` runs them in a row,
each with a fresh `Simulation`:

```python
hist = epsr.refine(
    proj, n_runs=200,
    overrides={"ntimes": 5},            # applied to every run
    initial_overrides={"ireset": 1},    # first run only
    callback=lambda k, sim: print(k, sim.summary.r_factor),  # return True to stop
    until=lambda h: h.r_factor[-1] < 5e-5,
    snapshot_every=50,                  # copies to <folder>.snapshots/run0050, ...
)
hist.r_factor, hist.energy, hist.ep_amplitude      # numpy arrays, one value per run
hist.r_factor_per_dataset                          # (n_runs, n_datasets)
hist.to_csv("history.csv"); hist.to_dataframe()
epsr.plot.history(hist)
```

A progress bar appears when tqdm is installed and stdout is a terminal
(`progress=True/False` forces it).

## 6. Perturbing the potential

`sim.potential` gives access to the potentials of a live simulation.

```python
pot = sim.potential
pot.empirical("O1", "H1")           # EP of one pair (kJ/mol on sim.r)
pot.reference(), pot.total()        # LJ + Coulomb; reference + empirical
pot.params                          # ereq, potfac, refpotfac, feedback, ncoeffp, pdstep, ...
pot.type_params                     # per atom type: epsilon, sigma, charge, mass, fraction
```

The EP of each pair is a sum of basis functions (Poisson or Gaussian, spacing
`pdstep`) with `ncoeffp - 1` coefficients; the engine builds the r-space tables
from them. Perturbations therefore act on the coefficients:

```python
c = pot.ep_coefficients()           # (n_pairs, ncoeffp - 1)
c[pot.pair_index("O1", "H1")] *= 0.5
pot.set_ep_coefficients(c)          # rebuilds the EP immediately
pot.perturb_ep(("O1", "H1"), 0.002)  # add to one pair's coefficients
pot.scale_ep(0.9)                   # scale every pair
pot.zero_ep()                       # like ireset 2
sim.run()                           # the edited EP drives the run and is saved to .pcof
```

**The ereq limit.** Whenever the EP is rebuilt, the engine measures its absolute
energy and scales all coefficients down so it does not exceed `ereq`. Scaling
down, zeroing and perturbing within the limit act as given. Scaling the EP beyond
`ereq` is undone. To allow a larger EP, raise the limit when you create the
simulation: `proj.simulation(overrides={"ereq": 20.0})`.

To run Monte Carlo with a fixed EP, switch refinement off with
`overrides={"potfac": [0.0, 0.0]}`; the coefficients are then saved unchanged.

The reference-potential parameters (`type_params`) are read-only; edit the `.ato`
file, for example with the `changeato` tool, to change them.

## 7. The setup and analysis tools

All 35 programs of the EPSR suite are available through `epsr.tools.run(folder,
command, args, answers)`, with the same arguments as on the EPSR command line.

```python
from epsr import tools
tools.commands()                                       # the 35 names

# keyword-file programs take their .dat file as argument
tools.nwts(proj.folder, "DMSODWaterD.NWTS.dat")        # neutron weights
tools.write_keyword_file(proj.folder / "w.PARTIALS.dat", "partials",
                         {"fnameato": "200DMSO400Water.ato", "atom-c": "O1"})
tools.partials(proj.folder, "w.PARTIALS.dat")
epsr.read_table(proj.folder / "w.PARTIALS.g01")

# interactive programs read their answers from standard input
tools.run(proj.folder, "writexyz", ["DMSO.ato"], answers=["Y", "1", "10 10", "-10 10", "10 20 30", "3 4"])
```

* `tools.keywords("nwts")` shows a program's keyword table, defaults and
  descriptions. `write_keyword_file` takes global keywords as a dict and repeated
  blocks (`component`, `distribution`, ...) as a list of dicts; missing keywords
  take the defaults.
* The answers each interactive program expects are easiest to copy from the
  `*.stdin` files of the cases in `test_problems/200DMSO400Water/EPSR-test-ex*`
  (their `README.md` explains each answer).
* Like the Fortran, `tools.run` silently does nothing when a keyword file is
  missing; the typed helpers (`tools.nwts`, `tools.partials`, ...) raise
  `FileNotFoundError` instead.
* Several programs use the current working directory (they read
  `system_commands.txt` and write `dirlist.txt` there). `tools.run` therefore
  switches to the folder for the duration of the call; tool runs are serialised
  by a lock.

## 8. Plotting

`epsr.plot` replaces the PGPLOT/gnuplot front ends, which are not ported. It needs
the `plot` extra.

```python
epsr.plot.fits(sim)                                 # data and fit for every dataset
epsr.plot.partials(sim, [("O1", "H1"), ("O1", "O1")])
epsr.plot.potentials(sim, which="empirical")
epsr.plot.history(hist, ["r_factor", "energy"])
```

Each function accepts `ax=` to draw into existing axes and returns the axes.

### Following a refinement iteration by iteration: `epsr.vis`

EPSR writes its result files, and computes the R-factor, only at the end of each
run, so `epsr.plot` shows the final state. `epsr.vis` records **every iteration**
and analyses the whole refinement. Replace `Simulation.run()` and `epsr.refine`
by their recorded versions. They run the same iterations
(`Simulation.iterations()`), so the folder ends up exactly as without recording
(`tests/test_vis.py` checks this file by file).

```python
from epsr import vis

recorder = vis.TrajectoryRecorder(case="water")      # case: settings in vis.systems
vis.run_recorded(project.simulation(), recorder)     # instead of simulation.run()
history = vis.refine_recorded(project, 2, recorder,  # instead of epsr.refine(...)
                              overrides={"niter": 100},
                              initial_overrides={"ereq": 3.0})
trajectory = recorder.trajectory()                   # arrays [iteration, ...]
recorder.save("work/trajectory")                     # vis.Trajectory.load(...) later

vis.make_report(trajectory, "figures/")              # 9 figures + summary.csv/.md
fit = vis.metrics.fit_metrics(trajectory)            # R-factor, R², chi² per iteration
vis.plots.partials_vs_reference(trajectory, vis.load_reference_partials({...}))
```

The figures show convergence (energy, pressure, R-factor, R², EP amplitude
against `ereq`, density, MC acceptance), the F(Q) fits and residuals, the
partials, the empirical potential, coordination numbers and bond angles. With
the `viewer` extra, `epsr.vis.viewer.ConfigurationViewer(trajectory).show()`
opens a 3D view of the configuration next to the partial g(r), with a slider
over the iterations. `screenshot()` and `movie()` (GIF or MP4) work off-screen.
`examples/05_experimental_water.py` and `examples/06_experimental_SiO2.py`
end with these figures.

`epsr.vis` needs the `plot` extra (matplotlib, scipy). It is the visualisation
toolkit of `visualisation/python_visualisation_toolkit` (standalone package
`epsrvis`, with documented examples). The shared modules are identical in both
places.

## 9. Parallel runs

**Inside one simulation.** The engine parallelises its expensive loops (energies,
pressure, structure factors, the S(Q) to g(r) transform, the potential fit) over
all CPU cores with rayon. Results are byte-identical at any number of threads. Set
`RAYON_NUM_THREADS` before starting Python to limit the cores used, for example
`RAYON_NUM_THREADS=4 uv run python script.py`; the pool size is fixed once the
first simulation runs.

**Several folders at once.** The engine releases the GIL while it runs, so
independent folders can also be refined from threads:

```python
projects = [epsr.Project.from_template("EPSR-test", f"work/run{i}") for i in range(4)]
histories = epsr.refine_many(projects, n_runs=50, overrides={"ntimes": 5})
```

Never use one folder from two threads at once. Tool runs change the working
directory and are therefore serialised. Concurrent simulations share the engine's
one pool of worker threads, so running several at once gives the same results
but does not make each of them faster; it helps most when `RAYON_NUM_THREADS` is
set low or the machine has many cores.

## 10. Errors

| exception | when | also a |
|---|---|---|
| `epsr.EpsrIOError` | a file is missing or unreadable | `OSError` |
| `epsr.EpsrInputError` | a malformed input file or an unknown override keyword | `ValueError` |
| `epsr.EpsrUnsupportedError` | an input on which the original Fortran itself crashes or hangs | `NotImplementedError` |

All three derive from `epsr.EpsrError`. The list of unsupported inputs is in
`rust/docs/FORTRAN_MAP.md` ("Still unsupported").

## 11. Troubleshooting

* *The progress messages do not appear in my notebook.* With `verbose=True` the
  engine prints from Rust straight to the process's standard output, which
  Jupyter does not capture. Watch the terminal that runs the kernel.
* *My edits to the EP had no effect.* Check whether the EP was already at `ereq`
  (see section 6), and whether refinement (`potfac`) changed it again during the run.
* *Results differ slightly between machines.* EPSR computes in single precision
  and uses the platform's libm; the last digits can differ outside Linux/glibc.
