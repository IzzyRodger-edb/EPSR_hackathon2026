# EPSR, reinvented

**EPSR** (Empirical Potential Structure Refinement, A. K. Soper) turns neutron and X-ray diffraction data into atomistic models of liquids, glasses and disordered materials. It is a standard analysis route for total-scattering data at the ISIS Neutron and Muon Source and other neutron and X-ray facilities. For 30 years it has been a Fortran program run through an interactive shell. It is hard to install and hard to script, its results are hard to reproduce, and it has no interface an agent can drive.

This repository is our entry to [**Reinvent the Wheel: A Tools for Science Hack**](https://luma.com/vdq0qw9u?tk=deOgNI) (Science Works × Inherent Laboratories, London, 16–18 October 2026).

## 1. Summary

We rebuilt EPSR as a modern, verified stack that a scientist can use from a terminal, a Python notebook, a web browser or an AI agent:

1. **A Rust translation of the whole Fortran EPSR suite**: the `epsr` engine and all 35 setup and analysis programs. It is a drop-in replacement, verified **byte for byte** against the original.
2. **A Python package, `epsr`**, that drives the Rust engine in-process. It returns NumPy arrays, runs refinement loops and draws publication figures through the `epsrvis` toolkit.
3. **Mellor**, an LLM agent that takes measured data and returns a finished, reproducible refinement with a PDF report.
4. **A web app** around Mellor: upload data, watch the run, get seaborn figures and a report, and push figures straight into an Overleaf paper.

> **Status**
>
> - ✅ Byte-identical to the Fortran on the supplied test case, 115 feature tests and 3 real-data refinements, at every step
> - ✅ Identical output at any thread count
> - ✅ 5 bugs in the original Fortran found and fixed, including one that silently corrupts three-body energies
> - ✅ 222 Python package tests · 51 web app tests (end to end on real water data)

Everything is MIT-licensed (see [LICENSE](LICENSE)).

---

## 2. What this means for the community

### A drop-in replacement

Nobody has to change how they work to benefit. `epsr-rs` takes **the same command-line arguments, reads the same input files and writes the same output files** as the Fortran programs:

```bash
# Fortran:   epsr  <folder> epsr 200DMSO400Water.EPSR.inp
epsr-rs EPSR-test/ epsr 200DMSO400Water.EPSR.inp
epsr-rs EPSR-test/ nwts DMSODWaterD.NWTS.dat      # the setup and analysis programs too
```

Because the output is byte-identical, an existing EPSR project folder (including one behind a published paper) can be continued, re-run or audited with the Rust engine and give exactly the same numbers. Groups can switch with no re-validation of their own and no change to their files, and the Python package, the agent and the web app all sit on top of that same verified engine.

### Before and after

| task | EPSR today | EPSR, reinvented |
|---|---|---|
| **Install** | build a Fortran suite and its shell by hand | `cargo build`, or `uv sync` for the Python package |
| **Run a refinement** | interactive EPSRshell session, typed commands | `epsr-rs` in a script, `sim.run()` in Python, or upload in a browser |
| **Analyse results** | parse fixed-format output files | partial g(r), F(Q), EPs and coordination numbers as NumPy arrays |
| **Set up a new system** | expert work: molecules, boxes, weights, schedules | Mellor proposes the system and writes the input files |
| **Reproduce a result** | lab notebook and memory | every run ships a standalone `run_<name>.py` that needs no LLM |
| **Check convergence** | read R-factors off the screen | per-iteration recording, convergence plots and a written verdict |
| **Make figures for a paper** | export and replot by hand | seaborn figures (PNG and PDF), pushed straight to Overleaf |
| **Use more cores** | OpenMP code with a data race | parallel Rust with output identical at any thread count |
| **Trust the code** | trust the binary | 115 feature cases with archived Fortran references; `cargo test` in seconds |

### Why it matters

* **Correctness of published science.** Translating line by line exposed real defects in the original ([FORTRAN_BUGS.md](FORTRAN_BUGS.md)). The worst one: after every save, the Fortran overwrites its three-body angles (stored as cosines) with radians, so from the second iteration of any multi-iteration run using three-body terms the energy is wrong. In one test case the Fortran reports +1954 kJ/mol per molecule where the correct value is −39, corrupting Monte Carlo acceptance and pressure. Anyone who has used three-body terms has been exposed to this. Five such bugs are now fixed in the Rust, and the rest are documented with proposed fixes.
* **Reproducibility.** A refinement becomes a script that anyone can re-run with open-source packages, instead of a sequence of typed shell commands.
* **A lower barrier to entry.** Setting up EPSR well takes expertise that students and facility users often don't have. Mellor and the web app let a non-specialist go from reduced data to a reviewed, documented refinement, while the expert can still read and edit every input file.
* **A maintainable future.** The method rests on a large legacy Fortran code. The Rust port is documented routine by routine ([rust/docs/FORTRAN_MAP.md](rust/docs/FORTRAN_MAP.md)) and pinned to the original's output by tests, so it can be extended without fear of silently changing results.
* **Part of the modern scientific toolchain.** EPSR can now run inside Jupyter, on HPC batch systems, in automated pipelines and under AI agents, alongside NumPy, pandas and the rest of the Python ecosystem.

---

## 3. Quick start

Tested on Linux (including WSL2 on Windows).

**Required**

| tool | why | install (Ubuntu / WSL) |
|---|---|---|
| Rust toolchain (stable, edition 2021) | builds the engine | `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \| sh` |
| C compiler and linker | needed by the Rust and Python builds | `sudo apt install build-essential` |
| Python ≥ 3.10 | Python package, Mellor, web app | usually present; `sudo apt install python3` |
| [uv](https://docs.astral.sh/uv/) | Python package management | `curl -LsSf https://astral.sh/uv/install.sh \| sh` |

**Optional**

| tool | why |
|---|---|
| `gfortran` | only to rebuild the original Fortran and re-run the comparison yourself |
| `git` | the web app's Overleaf publishing |
| a [Hugging Face token](https://huggingface.co/settings/tokens) | Mellor on a hosted model (a local model or the `scripted` backend need none) |

```bash
# the engine, as a drop-in replacement for the Fortran `epsr`
cd rust && cargo build --release && cargo test --release

# the Python package: refine water against real SANDALS data
cd ../python && uv sync --all-extras
uv run python examples/05_experimental_water.py --quick

# the web app with the agent, on a hosted model
cd ../webapp && uv sync --extra agent      # builds the Rust engine on first sync
export HF_TOKEN=hf_your_token_here
uv run epsrweb                             # http://127.0.0.1:8000
```

---

## 4. What we built

```
 browser ── webapp/ ── mellor/ (LLM agent) ── python/ (epsr package) ── rust/ (engine)
                │              │                       │
          accounts, jobs,   proposes the model,     epsr.vis / epsrvis
          seaborn, Overleaf  writes the report      (figures, per-iteration recording)
```

### The Rust engine

[`rust/`](rust/README.md) translates EPSR version 26 (`fortran_src/epsr_standalone_rev3.f` and everything it calls) into about 35,000 lines of Rust whose only dependency is `rayon` for parallel loops. The port covers every optional engine feature: X-ray data, non-cubic boxes, tethered molecules and dihedral groups, container forces, three-body terms, coordination constraints, coarse graining, pressure control, Bragg and diffuse scattering, and lattice refinement. All 35 setup and analysis programs (`makemole`, `mixato`, `fmole`, `nwts`, `partials`, `sdf`, `sharm`, …) are ported too, behind the single `epsr-rs` command. How the engine was parallelised without losing byte-identity is in [PARALLELISM_PLAN.md](PARALLELISM_PLAN.md).

### The Python package

[`python/`](python/README.md) is a PyO3 extension built with `uv sync`. It calls the Rust engine's public API unchanged, so its output is still byte-identical to the Fortran:

```python
import epsr
project = epsr.Project.from_template("../EPSR-test", "work/dmso")
sim = project.simulation(overrides={"ntimes": 5})
summary = sim.run()                       # one EPSR run in Rust, GIL released
sim.gr_inter("O1", "H1")                  # partial g(r) as a NumPy array
history = epsr.refine(project, n_runs=20) # refinement loop with convergence history
epsr.tools.run(folder, "nwts", ["sample.NWTS.dat"])   # any of the 35 tools
```

`epsr.vis`, built on the [`visualisation/`](visualisation/README.md) toolkit `epsrvis`, records every iteration and produces convergence plots, F(Q) fits and residual maps, partials, EPs, coordination numbers, bond angles and 3D configurations. Examples 05–07 refine water and silica against real data. The package has 222 tests.

### Mellor and the web app

The agent and the web app are described in sections [7](#7-mellor-an-ai-agent-for-structure-refinement) and [8](#8-the-web-app).

---

## 5. Validation against the Fortran

**Our standard was exact equivalence, not "close agreement".** Monte Carlo trajectories are chaotic: once any accept/reject decision comes out differently, the two runs drift apart for good. Byte-identical output after hundreds of coupled iterations is therefore a very strict test.

![Energy of the water refinement from the Fortran and the Rust over 300 steps, and the divergence caused by moving one molecule by a millionth of an ångström](docs/figures/rust_matches_fortran.png)

On the left, the Fortran and the Rust run the same 300-step water refinement and produce the same energy at every step. Every one of the 39 output files is identical too. On the right, both programs continue for 30 more steps. The Rust stays identical to the Fortran. A copy whose input differs only by moving one of the 1000 molecules by 10⁻⁶ Å (a millionth of the size of an atom) matches for 12 steps and then diverges. Nudges of 10⁻⁵ to 10⁻² Å diverge within one or two steps. (A 10⁻⁷ Å nudge, about the precision of the coordinate file, changed no Monte Carlo decision within 30 steps.) The run data and the script that draws these figures are in [docs/figures/](docs/figures/make_readme_figures.py).

The original Fortran was rebuilt as a reference (gfortran `-O2 -fno-tree-vectorize`), and the two codes were compared file by file at three levels:

| level | what was compared | result |
|---|---|---|
| **Supplied test case** ([`EPSR-test/`](EPSR-test), 200 DMSO + 400 water, three neutron data sets) | every output file over 12 chained EPSR runs, plus input variants; FNV hashes of every file are pinned in `cargo test` | byte-identical |
| **115 feature tests** ([`test_problems/200DMSO400Water/EPSR-test-ex*`](rust/docs/TEST_CASES.md)) | one case per feature and per tool: X-ray data, Gaussian/Morse potentials, triclinic boxes, `ereq` control, resets, every setup and analysis program, error paths and more; 232 chained passes, each with an archived Fortran reference | byte-identical (bar the 5 deliberate bug fixes below) |
| **Real experimental data** ([`validation/`](validation/README.md), see [section 6](#6-use-with-unseen-experimental-data)) | three systems built from scratch and refined for 300 iterations; Fortran, Rust CLI and Python wrapper run side by side, with **every file compared after every step** | all three cases: every step identical |

We also checked that the parallel Rust gives identical output at any thread count (`tests/thread_invariance.rs`). That check also fixes a data race present in the Fortran's OpenMP code.

**Fortran bugs found along the way.** [FORTRAN_BUGS.md](FORTRAN_BUGS.md) catalogues them. Five are fixed in the Rust as documented, deliberate divergences, the largest being the three-body corruption described in [section 2](#why-it-matters). The rest are reproduced on purpose so that the equivalence still holds, and each has a proposed fix and its blast radius. Inputs on which the Fortran crashes or reads uninitialised memory raise a clear "not supported by the Rust port" error instead.

To reproduce: [rust/docs/VALIDATION.md](rust/docs/VALIDATION.md) and [validation/README.md § How to reproduce](validation/README.md#how-to-reproduce).

---

## 6. Use with unseen experimental data

The port was developed against the supplied DMSO/water test case. To show that it works on data it was never tuned on, we built three new systems from scratch from the [Dissolve](https://www.projectdissolve.com/) examples and refined each against measured neutron data from the ISIS Neutron and Muon Source, running the Fortran, the Rust CLI and the Python wrapper side by side.

| case | system | data (ISIS Neutron and Muon Source instrument) | final R-factor |
|---|---|---|---|
| water | 1000 H₂O, 300 K | SANDALS: H₂O, D₂O, HDO | 2.3 × 10⁻⁴ (from 2.1 × 10⁻² at a random start) |
| silica | 1500 Si + 3000 O, amorphous | NIMROD | 2.3 × 10⁻⁴ |
| Cu-BTC (HKUST-1) | 624-atom MOF unit cell, with Bragg peaks | NIMROD | 1.6 × 10⁻² |

All three codes agree at every step. The overlay plots for [water](validation/water/gr_agreement.png), [silica](validation/silica_bulk/gr_agreement.png) and [Cu-BTC](validation/cu-btc/gr_agreement.png) show zero difference for every pair, because the three `.g01` files are literally the same file.

### What the refinement looks like

The figures below come from the water case in [`validation/water`](validation/water): 1000 water molecules fitted to three SANDALS measurements. The model starts as a random box of molecules, the simulated signal is pulled onto the measured one, and the structure that comes out is the known structure of water.

![Measured and simulated scattering signal of heavy water at the start, after 100 steps and after 300 steps](docs/figures/fit_before_after.png)

The first 100 steps use only the physics model (the reference potential), which gets the general shape right. Switching on the data-driven correction (the empirical potential) closes most of the remaining gap, most visibly around the main peak near Q = 2 Å⁻¹.

![Mismatch between simulation and experiment for the three water samples over 300 refinement steps](docs/figures/mismatch_history.png)

One model has to fit all three isotopic samples at once. That constraint is what makes the result trustworthy rather than a curve fit, and the mismatch for all three falls to about 2 × 10⁻⁴.

![Pair distribution functions of the refined model against the published reference partials of water](docs/figures/structure_vs_published.png)

The refined model shows the familiar features of liquid water: the O–O peak at 2.79 Å and the O···H hydrogen bond at 1.80 Å. Its pair distribution functions closely follow the published reference partials of Soper (2000). Two caveats apply. The reference is itself derived from water diffraction data by an EPSR-type refinement, so this is a consistency check rather than an independent measurement. Our 300-step schedule is also short, which is why the oxygen–oxygen curve is still noisy between 3 and 5 Å.

### Limitations and next steps

* **The agent needs validating on more unseen data.** The engine is verified against the Fortran, but Mellor's *choices* (the system it proposes, the potentials, the refinement schedule) have so far been checked only on the five systems in its validated library (water, silica, DMSO/water, benzene, Cu-BTC). The most important next step is to test Mellor's accuracy on additional experimental data it has never seen, from the ISIS Neutron and Muon Source and other facilities, against refinements by EPSR experts. Until then, any model outside the library is flagged as unvalidated in its report.
* **The LLM is an assistant, not an authority.** It never produces a number, but it can propose a poor starting model. Every proposal is saved as JSON and as a standalone script, so a scientist can review and correct it; the `scripted` backend runs the same pipeline with no model at all.
* **Some Fortran bugs are reproduced on purpose.** Fixing them would change every output and break byte-identity with the original. They are listed, with proposed fixes, in [FORTRAN_BUGS.md](FORTRAN_BUGS.md).
* **Not everything around the engine is ported.** The PGPLOT/gnuplot plotting front-ends are replaced by `epsrvis`, and `system_commands.txt` (shell configuration for the interactive EPSRshell) is not read.
* **Installation still needs a compiler.** There are no prebuilt wheels or binaries yet, so a Rust toolchain and a C compiler are required (see [section 3](#3-quick-start)). Testing so far has been on Linux and WSL2.
* **Web app extras.** Without the agent extra the web app runs in a clearly labelled demonstration mode, and direct Overleaf publishing needs an Overleaf plan with the Git bridge (otherwise a zip is offered).

---

## 7. Mellor: an AI agent for structure refinement

Setting up an EPSR refinement takes expertise: choosing molecules and force fields, building boxes, computing scattering weights, scheduling the reference-potential and EP stages, and reading the results. [**Mellor**](mellor/README.md) lowers that barrier. You give it Gudrun `.mint01` data files and, optionally, a sentence about the sample, and it returns a finished refinement:

1. **Inspect and propose.** A Hugging Face language model (hosted, e.g. Qwen3-235B, or local, e.g. Qwen3-4B) reads the data and proposes a *system specification* as JSON. When the sample matches one of its validated library systems (water, silica, DMSO/water, benzene, Cu-BTC), Mellor uses that system.
2. **Generate and run.** Ordinary code turns the specification into EPSR input files and a **standalone Python script**. It then runs the Rust engine through `epsr`: a scratch start, equilibration on the reference potential, then EP refinement.
3. **Analyse and report.** `epsrvis` makes the figures, the analysis gives fit quality per data set, a convergence verdict and first-shell distances and coordination numbers, and all of it goes into a PDF report.
4. **Discuss.** The model writes the discussion and answers follow-up questions. Its tools can query the results, search for comparable data (locally, the Dissolve examples, Zenodo) and launch new refinements.

**The model never produces a number.** It only proposes the specification and writes prose; every number comes from the verified engine. Every run folder contains `run_<name>.py`, which reproduces the whole calculation with only the open-source `epsr` and `epsrvis` packages and no LLM. Models outside the validated library are flagged as unvalidated in the report. A `scripted` backend runs the same pipeline with no model at all.

```bash
cd mellor && uv sync
uv run mellor run ../experimental_data/water/data/*.mint01 --out work/water \
    --backend hf-api --model Qwen/Qwen3-235B-A22B-Instruct-2507
uv run mellor chat work/water
```

Details: [mellor/docs/architecture.md](mellor/docs/architecture.md) and [mellor/docs/integration.md](mellor/docs/integration.md).

---

## 8. The web app

[`webapp/`](webapp/README.md) is a FastAPI application that puts Mellor in front of scientists who would never touch a Fortran shell:

| feature | how |
|---|---|
| **Upload data** | `.mint01` files are validated in the request, before a job is queued. The validator catches EPSR's silent traps: the 5000-point limit, data beyond `nq × qstep`, and file names with blanks. |
| **Agent runs EPSR** | one worker process per run calls `Mellor.run_pipeline`; progress is streamed live over Server-Sent Events; runs can be cancelled and survive server restarts |
| **Results** | seaborn figures drawn from Mellor's tidy CSV tables (PNG and vector PDF), the PDF report, the reproducible script and the EPSR folder; follow-up chat with the agent |
| **Publish to Overleaf** | tick figures and push them into the `figures/` folder of an Overleaf project through its Git bridge, with ready-made `\includegraphics` blocks. Without Overleaf premium, a zip plus `figures.tex` is offered instead. |
| **Accounts** | scrypt passwords, signed cookies, CSRF protection, per-user run folders as the trust boundary, Overleaf tokens encrypted at rest and never placed on a command line |
| **About and Instructions** | the theory and a user guide, at `/about` and `/instructions` |

```bash
cd webapp
uv sync --extra agent      # builds the Rust engine on first sync
uv run epsrweb             # http://127.0.0.1:8000
```

Without the agent extra the app runs in a clearly labelled **demonstration mode** with synthetic results. `/healthz` reports exactly which component is missing. The 51 tests run the HTTP layer end to end on real water data and pin every assumption the app makes about Mellor's API.

---

## 9. Supplementary: the theory, in brief

A diffraction experiment on a disordered material measures the total structure factor F(Q). It is a weighted sum of the partial structure factors of every pair of atom types α, β:

$$F(Q) = \sum_{\alpha \le \beta} (2-\delta_{\alpha\beta})\, c_\alpha c_\beta\, b_\alpha b_\beta \,\big[S_{\alpha\beta}(Q) - 1\big]$$

Here c are the atomic fractions and b the scattering lengths (form factors for X-rays). Each S_αβ(Q) is the Fourier transform of a partial pair distribution function g_αβ(r). A single measurement mixes all the partials together. Isotope substitution helps: H₂O, D₂O and HDO have the same structure but different weights b, so measuring several of them gives several independent combinations of the same partials.

EPSR inverts these data with a model:

* **Monte Carlo on a reference potential.** A periodic box of molecules (thousands of atoms) is simulated with a physically sensible potential: Lennard-Jones plus Coulomb, with intramolecular geometry held by harmonic constraints.
* **An empirical potential (EP).** Every few cycles, the simulated F(Q) of each data set is compared with the measurement. The residual, D(Q) − F_sim(Q), is inverted through the weights matrix into a correction for each atom pair. That correction is added to the potential as a sum of Poisson functions in r, which pushes the simulation towards the data.
* **A cap on the EP.** The EP's total amplitude is limited by `ereq`, so the data refine the model rather than overwrite the physics. A large EP is a sign that the reference potential is poor.
* **Ensemble averages.** Once the fit stabilises, EPSR accumulates partial g(r), coordination numbers, bond angles and spatial density functions over many configurations. All of them are consistent with every data set at once.

Fit quality is given by the R-factor, the mean squared difference between data and simulation over the fitted Q range. For water refined against good neutron data it falls to a few × 10⁻⁴.

Key references: Soper, *Chem. Phys.* **202**, 295 (1996); Soper, *Phys. Rev. B* **72**, 104204 (2005); Soper, *EPSRshell user's guide*, RAL-TR-2011-012. More references are in [mellor/docs/integration.md §7.3](mellor/docs/integration.md).

---

## 10. Repository map and acknowledgements

| path | contents |
|---|---|
| [`rust/`](rust/README.md) | the Rust engine (`epsr-core`), the drop-in CLI `epsr-rs`, the Fortran reference build and comparison tools, and the architecture, Fortran-map and validation docs |
| [`python/`](python/README.md) | the `epsr` Python package (PyO3), examples, tests, user guide and API reference |
| [`visualisation/`](visualisation/README.md) | the per-iteration recorder and the `epsrvis` analysis and plotting toolkit |
| [`validation/`](validation/README.md) | water, silica and Cu-BTC refined against measured data by the Fortran, the Rust CLI and the Python wrapper |
| [`mellor/`](mellor/README.md) | the LLM agent |
| [`webapp/`](webapp/README.md) | the web application around Mellor |
| [`test_problems/`](rust/docs/TEST_CASES.md) | the 115 `EPSR-test-ex*` verification cases with archived Fortran references |
| [`EPSR-test/`](EPSR-test) | the original 200 DMSO / 400 water test case |
| [`experimental_data/`](experimental_data) | the Dissolve example systems and measured neutron data |
| [`fortran_src/`](fortran_src) | the original Fortran, untouched |
| [`docs/figures/`](docs/figures/make_readme_figures.py) | the figures in this README, the script that draws them and the run data they use |
| [FORTRAN_BUGS.md](FORTRAN_BUGS.md) | bugs found in the Fortran and their status in the port |
| [PARALLELISM_PLAN.md](PARALLELISM_PLAN.md) | how the engine was parallelised without losing byte-identity |
| [DATA_ACQUISITION_FEASIBILITY.md](DATA_ACQUISITION_FEASIBILITY.md) | feasibility study for an agent that fetches reduced data from the ISIS Neutron and Muon Source, ORNL and other facilities |

### Acknowledgements

EPSR is the work of Alan Soper (ISIS Neutron and Muon Source, STFC). The example systems and data come from the [Dissolve](https://www.projectdissolve.com/) project (T. G. A. Youngs). Built at Reinvent the Wheel 2026, hosted by Science Works and Inherent Laboratories.
