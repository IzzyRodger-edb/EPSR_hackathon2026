# API reference

Everything listed here is importable from `epsr` unless a module is named.
Arrays are numpy `float32` copies; changing them never changes the simulation.
Dataset indices are 0-based (except in `overrides`, see `render_overrides`).
A *pair* argument accepts an index into `pair_labels`, two atom-type labels
(`"O1", "H1"`) or a label tuple (`("O1", "H1")`), in either order.

The Fortran variable name of each quantity is given in brackets where it helps
cross-referencing with `rust/docs/ARCHITECTURE.md`.

## Simulation

```python
Simulation(folder, inp_name=None, *, overrides=None, verbose=False)
```

One EPSR run on `folder`. `inp_name` is detected when the folder holds exactly
one `*.EPSR.inp`. Reading the inputs rewrites the `.qdr` file; the last iteration
rewrites `.ato`, `.pcof`, `.EPSR.inp`, `.EPSR.inpa` and all result files.

**Running**

| member | description |
|---|---|
| `run() -> Summary` | all `niter` iterations, then the output files (`Epsr::run`); releases the GIL; allowed once |
| `iterate(i) -> Summary` | iteration `i`; must be the next one, `1..niter`; the last writes the output files |
| `iterations()` | generator over the remaining iterations, yielding `Summary` |
| `last_iteration`, `completed` | progress of this run |
| `summary -> Summary` | results of the most recent iteration |

**Inputs and metadata**

| member | description |
|---|---|
| `folder`, `inp_name`, `run_name`, `overrides` | paths, run name (prefix of the outputs), override lines applied |
| `output_path(ext) -> Path` | e.g. `output_path("g01")` |
| `control -> Control` | run-control values (`ntimes`, `niter`, `nq`, `qstep`, ...) |
| `datasets -> list[Dataset]` | the `data n` blocks |
| `get_parameter(name, dataset=None) -> str or None` | current value of any keyword |
| `n_pairs`, `n_datasets`, `nq`, `nr`, `n_atoms`, `n_molecules` | sizes |
| `pair_labels`, `atom_types`, `pair_index(pair_or_type, other_type=None)` | pair and type labels; pair lookup |
| `weights_files` | weights-file paths per dataset |

**Results**

| member | shape | description |
|---|---|---|
| `q`, `r` | `(nq,)`, `(nr,)` | Q grid (Å⁻¹), r grid (Å) |
| `gr_inter(pair=None, other_type=None)` | `(n_pairs, nr)` or `(nr,)` | running-average intermolecular partial g(r) (`gorinter`) |
| `gr_intra(...)` | same | intramolecular partial g(r) (`gorintra`) |
| `sq_inter(...)`, `sq_intra(...)` | `(n_pairs, nq)` or `(nq,)` | partial S(Q) (`soqinter`, `soqintra`) |
| `current_gr()`, `current_sq()` | `(2, n_pairs, n)` | current configuration only, index 0 inter, 1 intra (`wav`, `qwind`) |
| `data_fq(i)`, `data_errors(i)` | `(nq,)` | measured F(Q) on the Q grid and its errors (`datai`, `errs`) |
| `simulated_fq(i)` | `(nq,)` | simulated total F(Q) |
| `data_q_range(i)` | `(lo, hi)` | fitted Q indices, inclusive |
| `fit(i) -> Fit` | | data, errors, simulation over the fitted range |
| `positions` | `(n_atoms, 3)` | atom positions (Å) |
| `atom_labels`, `molecule_atoms` | | per-atom labels; `range` of atom indices per molecule |
| `box_size`, `density` | | box side (Å), atoms/Å³ (`dmax`, `rho`) |
| `potential -> Potential` | | see below |

## Summary, Control, Dataset, Fit

Frozen dataclasses.

* `Summary`: `iteration` (`iter`), `energy_inter` (`foldsv`), `energy_intra`
  (`foldsb`), `energy` (`energd`), `pressure` (`pressd`), `ep_amplitude`
  (`energabs`), `r_factor` (`rfactot`), `r_factor_raw` (`rfacraw`),
  `r_factor_per_dataset` (`rfacd`), `ep_up` (`epcup`), `ep_coefficient`
  (`epccoeff`), `gr_roughness` (`sentb`), `chisq`; `fortran` (raw dict);
  `as_dict()`.
* `Control`: the fields of the Rust `Control` struct (`nq`, `qstep`, `ireset`,
  `iinit`, `ntimes`, `niter`, `intra`, `rotfreq`, `inter`, ..., `pressure`).
* `Dataset`: `infili` (data file as written in the input), `fnamedata`,
  `fnamewts`, `nrtype`, `rshmin`, `szeros`, `tweak`, `efilereq`, `porodadd`,
  `nqmind`, `nqmaxd`, `porodamp`; property `name` (data file stem).
* `Fit`: `name`, `q`, `data`, `errors`, `simulated`, `q_range`; property `residual`.

## render_overrides

```python
render_overrides(overrides) -> list[str]
```

Converts overrides to `.EPSR.inp` lines. A sequence of strings passes through.
A mapping gives `"key value"` lines; lists and tuples become space-separated
values, booleans `T`/`F`. The key `"data"` maps **1-based** EPSR dataset numbers
to keyword mappings and is emitted after all global keywords:
`{"ereq": 12.5, "data": {2: {"tweak": 1.05}}}` gives
`["ereq 12.5", "data 2", "tweak 1.05"]`.

## Potential

`sim.potential`; a view onto the live simulation.

| member | description |
|---|---|
| `r` | r grid of the tables (Å) |
| `empirical(pair=None, other_type=None)` | EP (kJ/mol, `xavs`), `(n_pairs, nr)` or `(nr,)` |
| `reference(...)` | Lennard-Jones + Coulomb (`ecoul`) |
| `total(...)` | reference + empirical (`xavtot`) |
| `ep_coefficients()` | `(n_pairs, ncoeffp - 1)` EP coefficients (`pcoefp(2:)`) |
| `ep_flags()` | `(n_pairs,)` refinement flags (`pcoefp(1)`) |
| `set_ep_coefficients(coef)` | replace the coefficients and rebuild the EP |
| `set_ep_flags(flags)` | set the refinement flags (scalar or per pair) |
| `perturb_ep(pair, delta)` | add `delta` (scalar or vector) to one pair's coefficients; `pair` is an index or a label tuple |
| `scale_ep(factor)` | multiply every coefficient (`multiply_pot_coeffs`) |
| `zero_ep()` | set the EP to zero |
| `shared_pairs()` | groups of pairs that share one set of `.pcof` coefficients |
| `params -> PotParams` | `ereq`, `potfac`, `potfacamp`, `refpotfac`, `feedback`, `ncoeffp`, `pdmin`, `pdstep`, `pdmax`, `psigma2`, `gaussian` |
| `type_params -> list[TypeParams]` | `label`, `epsilon` (kJ/mol), `sigma` (Å), `charge` (e), `mass` (amu), `fraction` |
| `pair_index(pair_or_type, other_type=None)` | as `Simulation.pair_index` |

After every write the engine rebuilds the EP tables (`generate_ep`) and applies
the `ereq` amplitude limit: if the EP's absolute energy exceeds `ereq`, all
coefficients are scaled down to it. Pairs listed by `shared_pairs()` must be
given identical rows; `perturb_ep` updates all of them.

## Project

```python
Project(folder, inp_name=None)
Project.from_template(template, work, *, inp_name=None, overwrite=False)
```

| member | description |
|---|---|
| `folder`, `inp_name`, `run_name` | |
| `simulation(overrides=None, *, verbose=False)` | a new `Simulation` on the folder |
| `snapshot(label, *, overwrite=False) -> Path` | copy to `<folder>.snapshots/<label>/` |
| `output(ext) -> Path` | `<folder>/<run_name>.<ext>` |
| `read_table(ext) -> Table` | a result table, merging `…02`, `…03` continuation files |
| `history() -> History` | the `.erg` file |
| `input_file() -> KeywordFile` | the parsed `.EPSR.inp` |
| `files() -> dict` | paths of the `inp`, `ato`, `pcof`, `data`, `weights`, `erg`, `out` files |

`copy_folder(src, dst, *, overwrite=False)` copies the regular files of a folder.

## refine, refine_many, History

```python
refine(project, n_runs, *, overrides=None, initial_overrides=None, callback=None,
       until=None, snapshot_every=None, verbose=False, progress=None) -> History
refine_many(projects, n_runs, *, workers=None, **refine_kwargs) -> list[History]
```

`refine` runs `n_runs` EPSR runs on one folder, each with a fresh `Simulation`.
`overrides` apply to every run, `initial_overrides` to the first only.
`callback(k, sim)` is called after run `k` (1-based); returning `True` stops.
`until(history)` stops when true. `snapshot_every=n` snapshots every `n` runs.
`refine_many` refines different folders in parallel threads.

`History` holds columns as numpy arrays, accessed as attributes (`hist.r_factor`)
or with `column(name)`. From `refine`: `run`, `iteration`, `energy`,
`energy_inter`, `energy_intra`, `pressure`, `ep_amplitude`, `r_factor`,
`r_factor_raw`, `ep_up`, `ep_coefficient`, `gr_roughness`, `chisq`,
`r_factor_per_dataset` (2-D); `summaries` keeps the `Summary` objects. From
`History.from_erg(path)`: `run` and the `.erg` columns `energy`, `pressure`,
`r_factor_raw`, `ep_amplitude`, `density`, `ep_up`, `ep_coefficient`,
`parab_a`, `parab_b`, `parab_c`. Methods: `append(summary)`,
`from_summaries(list)`, `to_csv(path)`, `to_dataframe()`, `len()`, indexing
(row as dict), iteration.

## epsr.tools

| function | description |
|---|---|
| `commands() -> list[str]` | the 35 program names |
| `run(folder, command, args=(), answers="", *, verbose=False)` | run a program; `answers` is a string or a list of lines for standard input; runs inside `folder` |
| `keywords(command) -> ToolKeywords` | keyword table: `extinp`, `dataprompt`, `symbols` (`ToolKeyword` entries with `name`, `value_type`, `flag`, `default`, `description`), `global_keywords`, `block_keywords`, `defaults()`, `descriptions()` |
| `write_keyword_file(path, command, values, blocks=None, *, title=None) -> Path` | write a keyword input file in EPSR's layout |
| `nwts`, `xwts`, `partials`, `coord`, `triangles`, `torangles` `(folder, dat_file)` | keyword-file programs; raise `FileNotFoundError` for a missing file |
| `writexyz(folder, ato_file, answers="")` | convert a configuration to `.xyz` |

## epsr.io

| function or class | description |
|---|---|
| `read_table(path) -> Table` | an EPSR column file; `Table.x`, `x_label`, `columns`, `errors`, `labels`, `table[label]`, `to_dataframe()` |
| `read_table_series(first) -> Table` | the same, merging continuation files `…02`, `…03` |
| `read_erg(path) -> ndarray` | `(n, 10)`; column names in `ERG_COLUMNS` |
| `read_keyword_file(path) -> KeywordFile` | `title`, `values`, `blocks`, `block_keyword`, `descriptions`, `trailer` (lines after `q`), `words(key)` |

Repeated labels in a table header get a `#2` suffix (`.z01` lists each pair twice).

## epsr.plot

`fits(sim, datasets=None, *, offset=1.0, ax=None)`,
`partials(sim, pairs=None, *, kind="gr_inter", offset=1.0, ax=None)`,
`potentials(sim, pairs=None, *, which="empirical", ax=None)`,
`history(hist, columns=("r_factor", "energy", "ep_amplitude"), *, axes=None)`.
Each returns the axes it drew on.

## epsr.vis

Analysis of a refinement iteration by iteration (`plot` extra; the 3D viewer
needs the `viewer` extra). See the user guide, section 8.

| name | purpose |
|---|---|
| `TrajectoryRecorder(case="", record_positions=True)` | collects one snapshot per iteration; `.trajectory()`, `.save(folder)`, `.nsteps` |
| `run_recorded(sim, recorder, *, iteration_callback=None) -> Summary` | `sim.run()`, recorded |
| `refine_recorded(project, n_runs, recorder, *, overrides, initial_overrides, callback, iteration_callback, verbose) -> History` | `epsr.refine`, recorded |
| `record_run(sim, recorder, ereq=None, pcof=None, callback=None)` | the lower-level loop used by `run_recorded` |
| `Trajectory.load(folder)` / `.save(folder)` | per-iteration arrays `fq_data, fq_err, fq_sim, fq_inst` [step, dataset, Q], `gr_inter, gr_intra, ep, ref_pot` [step, pair, r], `sq_inter, sq_intra` [step, pair, Q], `positions` [step, atom, 3], and `history` (convergence numbers per step) |
| `make_report(trajectory, folder, figures=None, gif=False)` | figures `01_convergence` ... `09_displacement`, `summary.csv`, `summary.md` |
| `metrics` | `fit_metrics` (R-factor, R², chi² per iteration), `acceptance`, `running_coordination`, `shell_history`, `coordination_distribution`, `bond_angle_distribution`, `tetrahedral_order`, `mean_squared_displacement`, `neighbour_pairs` |
| `plots` | one function per figure, plus `partials_vs_reference`, `fq_animation`, `final_outputs` |
| `read_isis_gr`, `load_reference_partials` | experimental partial g(r) files of the ISIS Disordered Materials Database |
| `viewer.ConfigurationViewer(trajectory, gr="inter", slab=None, color_by="type", reference=None, off_screen=False, window_size=(1600, 820))` | pyvista 3D configuration + partial g(r): `show()`, `screenshot(path, step)`, `movie(path, every, fps)` |
| `live.LiveMonitor(every=5)` | `iteration_callback` that redraws the R-factor and the F(Q) fits during a run |

`python -m epsr.vis report|view|snapshot|movie|final ...` is the command line.

## Errors

`EpsrError`; `EpsrIOError(EpsrError, OSError)`; `EpsrInputError(EpsrError,
ValueError)`; `EpsrUnsupportedError(EpsrError, NotImplementedError)`. Run-order
mistakes (`iterate` out of order, `run` twice) raise `ValueError` or
`RuntimeError`.

## epsr._core

The extension module mirrors the Rust API with Fortran names; `epsr/_core.pyi`
lists its signatures. Use it only when the Python layer lacks something; see
[architecture.md](architecture.md).
