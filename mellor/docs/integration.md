# Integrating Mellor into a web application

This guide is for the developer (or coding agent) building a web app around Mellor.
It says what Mellor is, how to install it next to your app, which Python calls to
make, what files a run produces, and what the app has to provide itself (workers,
storage, uploads). Mellor has no HTTP layer of its own; the app is a thin layer over
the Python API below.

Status: the package is being finished on 2026-10-17. The names and signatures in
this guide are the project's fixed contract (`docs/architecture.md`, `src/mellor/spec.py`
and the module docstrings). Two things are still landing today and are marked
"(planned)" below.

## 0. The brief, mapped onto Mellor

The web app is specified as: upload data; the LLM takes it and runs the Rust EPSR
code; results come back with seaborn plots and PDF summaries; selected figures can be
clicked to upload them to a `/figures` folder of an Overleaf project; users and
authentication; an About page (theory) and an Instructions page. This is how each
requirement maps onto what Mellor provides:

| requirement | Mellor side | the app's side |
|---|---|---|
| upload data | `.mint01` files (two to four per run is typical) and an optional free-text prompt are the inputs of `Session`/`run_pipeline` (§4.1); `mellor.data.read_mint01` validates a file | upload form, per-user storage, size limits (§6) |
| the LLM runs the Rust EPSR code | `Mellor.run_pipeline` does the whole thing: the model proposes the model system, code generates the inputs and runs the engine (§1, §4) | a job queue and a worker; progress via `on_event` (§4.3) |
| results: plots and PDF | `report.pdf`, `figures/*.png` (matplotlib, made by the repository's `epsrvis` toolkit), `analysis.json` (§5); tidy CSV tables for seaborn plots the app draws itself (§5.1) | results page, seaborn rendering (§5.1), downloads |
| click a figure to upload it to Overleaf | the figure files and their captions (`mellor.figures.FIGURE_CAPTIONS`) | the Overleaf Git bridge with a per-user token (§7.1) |
| users and authentication | nothing; Mellor is per run folder | accounts, per-user run folders and quotas, token storage (§7.2) |
| About and Instructions pages | the theory summary, references and figure explanations in §7.3 and the documents it points to | the pages |

## 1. What Mellor is

Mellor is a Python package (`mellor/` in the `EPSR_hackathon` repository) that turns
measured structure-factor data (Gudrun `.mint01` files) into a finished EPSR
refinement with a language model from Hugging Face in the loop:

1. the model inspects the data and proposes a *system specification* (JSON);
2. code generates the EPSR input files and a standalone Python script, runs the Rust
   EPSR engine through the `epsr` wrapper, records the refinement and analyses it;
3. figures come from the repository's `epsrvis` toolkit; a PDF report is assembled;
4. the model writes the discussion and then answers follow-up questions, with tools
   to query the results and to look for comparable data (local, Dissolve examples,
   Zenodo).

The model never runs anything itself: it only produces the specification and prose.
All numbers come from EPSR. Your app therefore never has to run EPSR, parse EPSR
files or talk to the model directly; it calls Mellor and serves the files Mellor writes.

## 2. Installation

Mellor depends on two sibling packages of the same repository by path
(`../python`, the EPSR wrapper with a Rust extension, and
`../visualisation/python_visualisation_toolkit`, the `epsrvis` toolkit), so the app
needs the repository checkout and a Rust toolchain for the one-off build.

```bash
# prerequisites: Python >= 3.10, uv, Rust (stable >= 1.85)
curl -LsSf https://astral.sh/uv/install.sh | sh
curl https://sh.rustup.rs -sSf | sh
export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"

git clone <EPSR_hackathon repository> && cd EPSR_hackathon/mellor
uv sync                   # core: epsr (compiles the engine, a few minutes), epsrvis, reportlab, huggingface_hub
uv sync --extra local     # only if the server has a GPU and should run the model locally (adds torch, transformers)
uv run pytest -m "not slow and not gpu and not network" -q   # sanity check, seconds
```

To use Mellor from your own project instead of from the `mellor/` folder, declare it
as a path dependency; with uv:

```toml
[project]
dependencies = ["mellor"]
[tool.uv.sources]
mellor = { path = "../EPSR_hackathon/mellor" }
```

Keep the checkout in place: `mellor` resolves its two sibling dependencies relative
to its own folder.

## 3. Configuration

| variable / option | meaning |
|---|---|
| `HF_TOKEN` | Hugging Face token (fine-grained, with "Make calls to Inference Providers"). Needed by the `hf-api` backend. |
| `MELLOR_BACKEND` | default backend when none is given: `auto` (local model if torch is installed, else `hf-api` if `HF_TOKEN` is set), `transformers`, `hf-api`, `scripted`. |
| `MELLOR_MODEL` | default model id. The local default is `Qwen/Qwen3-4B-Instruct-2507` (8 GB bf16). |
| `RAYON_NUM_THREADS` | caps the CPU threads of the EPSR engine (it uses every core by default). |

Backend strings accepted by `mellor.llm.backend.make_backend`:

* `"hf-api:Qwen/Qwen3-235B-A22B-Instruct-2507"`: hosted inference through Hugging
  Face Inference Providers, billed per token to the account of `HF_TOKEN`. This is
  the backend for a server without a GPU. Good tool-calling models available today:
  `Qwen/Qwen3-235B-A22B-Instruct-2507`, `openai/gpt-oss-120b`,
  `meta-llama/Llama-3.3-70B-Instruct`, `deepseek-ai/DeepSeek-V3.2`.
  (planned, today) `provider=` and `bill_to=` options so the app can pin a provider
  and bill an organisation.
* `"transformers[:model]"`: the model runs in the server process on the GPU (needs the
  `local` extra, about 8.4 GB of GPU memory for the default model).
* `"scripted"`: no model; replies come from a script. Use it in tests and for a
  "no LLM" mode (see §4.4).

## 4. The Python API

Everything below is importable without torch; the local backend imports torch only
when it is constructed.

### 4.1 One run, end to end

```python
from pathlib import Path
from mellor.llm.backend import make_backend
from mellor.llm.tools import Session
from mellor.llm.agent import Mellor

session = Session(
    work_dir=Path("/srv/mellor/runs/0042"),          # one fresh folder per run
    data_files=[Path("/srv/uploads/0042/SLS18498-H2O.mint01"), ...],
    quick=False,                                     # True: 10-iteration schedule for previews/tests
)
backend = make_backend("hf-api:Qwen/Qwen3-235B-A22B-Instruct-2507")

events = []
mellor = Mellor(backend, session, on_event=events.append)   # on_event: progress stream, see §4.3
outcome = mellor.run_pipeline(user_prompt="Liquid water at 25 C, three isotopic samples.")

outcome.report_path      # Path to <work_dir>/report.pdf
outcome.final_message    # short text for the user
outcome.spec_source      # "model", "library:water", "user" or "fallback:library:water"
session.analysis         # the analysis dict (also saved as analysis.json)
session.figures          # {"01_convergence": Path(...), ...}
```

`run_pipeline` blocks for the whole run: a few minutes with `quick=True`, roughly
ten minutes for the full water schedule on a 20-core machine (300 EPSR iterations),
plus the model's time. Run it in a worker (§6).

### 4.2 Follow-up chat

```python
answer = mellor.chat("What is the O-H hydrogen-bond distance and how good is the fit?")
```

`chat` keeps the conversation of the same `Mellor` object. To continue a run in a
later process (a new request, after a restart), rebuild the session from its folder:

```python
session = Session(work_dir=Path("/srv/mellor/runs/0042"))
session.load_existing()                      # reads spec.json, manifest.json, analysis.json, figures/, report.pdf
mellor = Mellor(make_backend("hf-api:..."), session)
answer = mellor.chat("Find comparable benzene data and tell me what it is.")
```

The earlier conversation is not restored into the model's context (the transcript is
on disk, `mellor_transcript.jsonl`, if the app wants to show it); the run's results
are, through the tools. The model may call `find_similar_data`, `download_data`
(into `<work_dir>/data`) and `build_and_run` (a new refinement in the same folder,
which overwrites the previous run's `epsr/` and results) during chat. If the app
wants each refinement kept, give every run its own `work_dir`.

### 4.3 Progress events

`on_event` receives dicts while the pipeline or a chat turn runs:

| `type` | fields | when |
|---|---|---|
| `phase` | `phase`: `setup`, `narrative`, `report` | the pipeline enters a phase |
| `assistant` | `content`, `tool_calls` | the model replied |
| `tool_call` | `name`, `arguments` | a tool is about to run |
| `tool_result` | `name`, `result` (truncated JSON) | a tool finished |
| `progress` | `line` | one line of the EPSR run log (setup programs, iterations) |
| `note` | `content` | fallbacks and warnings (e.g. library fallback used) |

They are suitable for Server-Sent Events or a WebSocket. The callback runs in the
worker's thread; keep it quick and non-blocking.

### 4.4 Without a model (fixed specification)

The pipeline functions are plain Python and can be called without any model, for a
"use the validated water model" button, for tests, or when no token is configured:

```python
from mellor import data, library, runner, analysis, figures, report

datasets = [data.read_mint01(p) for p in paths]
spec = library.make_spec("water", datasets)                 # or SystemSpec.load("spec.json")
result = runner.run_spec(spec, work_dir, quick=False, progress=print)
summary = analysis.analyse(result)
analysis.save_analysis(summary, work_dir / "analysis.json")
pngs = figures.make_figures(result)
pdf = report.build_pdf(result, summary, pngs, analysis.analysis_markdown(summary))
```

`library.KNOWN_SYSTEMS` lists the validated systems (`water`, `silica`,
`dmso_water`, `benzene`, `cu_btc`); `library.auto_spec(datasets)` picks one from the
data titles. `mellor.spec.SystemSpec` is the specification; `SystemSpec.from_dict`
and `validate()` give precise error lists, which an app can surface in a form.
`mellor.export.write_tables(work_dir)` adds the CSV tables of §5.1.

## 5. What a run folder contains (what to serve)

Every run is one self-describing folder (`session.work_dir`):

| path | content | serve as |
|---|---|---|
| `report.pdf` | the full report | download |
| `figures/*.png`, `figures/summary.csv`, `figures/summary.md` | the figures (captions in `mellor.figures.FIGURE_CAPTIONS`) and the summary table | gallery, table |
| `tables/*.csv` (after `mellor.export.write_tables`) | long-format tables of the fits, partials, potentials and history (§5.1) | the app's own seaborn plots |
| `analysis.json` | every computed number: runs, fit quality per data set, convergence verdict, first-shell table, flags | results page |
| `manifest.json` | per-run summaries written by the script (R-factor, energy, pressure, EP amplitude, timings, versions) | status |
| `spec.json` | the specification that was run | "model" tab, editable re-run |
| `run_<name>.py` | the standalone Python script that reproduces the run with only `epsr` and `epsrvis` | download ("check it yourself") |
| `epsr/` | the EPSR simulation folder (inputs, `.ato`, `.pcof`, every output file) | download as archive |
| `trajectory/` | per-iteration record in the `epsrvis` format (tens to a few hundred MB) | keep for chat; do not serve raw |
| `data/` | copies of the input `.mint01` files | – |
| `discussion.md` | the model's discussion text | results page |
| `mellor_transcript.jsonl` | every model message and tool call | debug / transparency view |
| `run.log` | console output of the script | live log |

Re-running `mellor report <folder>` (or `report.build_pdf`) regenerates the figures
and the PDF from these files without a model.

### 5.1 Tidy tables for seaborn and pandas

Mellor's own figures are matplotlib PNGs made by `epsrvis`. For the app's seaborn
plots, do not parse EPSR files: `mellor.export` gives the results as long-format
tables (lists of flat row dicts, or CSV files) that pandas reads directly.

```python
from mellor import export
paths = export.write_tables(run_dir)            # -> {"fits": run_dir/"tables"/"fits.csv", ...}
tables = export.tidy_tables(run_dir)            # the same as lists of dicts, no files
```

| table | columns | one row per |
|---|---|---|
| `fits` | `dataset, q, data, error, simulated, residual` | fitted Q point of each data set (final iteration) |
| `partials` | `pair, r, g_inter, g_intra, g_total` | r point of each atom pair (final partial g(r)) |
| `potentials` | `pair, r, empirical, reference` | r point of each pair (kJ/mol) |
| `history` | `step, run, iteration, ereq, energy, pressure, ep_amplitude, r_factor, r_factor[<dataset>]…, r2[<dataset>]…` | EPSR iteration |
| `runs` | `run, stage, niter, ereq, r_factor, energy, pressure, ep_amplitude, duration_s` | EPSR run |

Seaborn examples (add `pandas` and `seaborn` to the app's own dependencies):

```python
import pandas as pd, seaborn as sns
fits = pd.read_csv(run_dir / "tables" / "fits.csv")
g = sns.relplot(data=fits.melt(["dataset", "q"], ["data", "simulated"]), x="q", y="value",
                hue="variable", col="dataset", kind="line", facet_kws={"sharey": False})
g.set_axis_labels("Q (Å⁻¹)", "F(Q)")

partials = pd.read_csv(run_dir / "tables" / "partials.csv")
sns.lineplot(data=partials[partials.r < 8], x="r", y="g_total", hue="pair")

history = pd.read_csv(run_dir / "tables" / "history.csv")
sns.lineplot(data=history, x="step", y="r_factor").set(yscale="log")
```

The tables come from the recorded trajectory; without one (`--no-record`) they fall
back to the EPSR output files, with no errors and no per-iteration R-factor. The
first-shell numbers (peak positions, coordination numbers) are in `analysis.json`
under `shells`, and the per-run numbers under `runs`.

## 6. Process model and limits

* **Blocking calls.** `run_pipeline`, `chat` and `runner.run_spec` block. Put them in
  a worker (thread pool for the `hf-api` backend; a separate process per job is
  better because the EPSR run is CPU-bound and the model call is not).
* **EPSR runs in a subprocess** (the generated script), so the app's process never
  changes its working directory and a job can be killed by killing the worker.
  Cancellation inside Mellor is not implemented: kill the worker process and delete
  or mark the folder.
* **One folder, one run at a time.** Never run two pipelines on the same `work_dir`.
* **The local backend holds one model per process** (8.4 GB of GPU memory for the
  default model) and is not safe for concurrent `generate` calls from threads; use
  one worker process per GPU, or the `hf-api` backend, which is stateless.
* **CPU.** The engine uses every core; set `RAYON_NUM_THREADS` per worker when
  several jobs share a machine.
* **Token budget.** A full pipeline is a few tens of thousands of tokens of
  conversation; with the hosted models listed above that is well under 5 p per run.
* **Disk.** A water run is a few hundred MB with the trajectory.
* **Uploads.** `mellor.data.read_mint01` validates the file format and raises
  `ValueError` on anything that is not a three-column Gudrun table; impose your own
  size and count limits (typical files are 100 kB to 1 MB; two to four files per run).
* **What the model can do.** Through its tools the model can read files inside the run
  folder, download files from `http(s)` URLs into `<work_dir>/data` (50 MB limit) and
  start refinements. It cannot run shell commands or touch anything outside the run
  folder. Give each run its own folder under a quota and treat the folder as
  user-controlled data.

## 7. A minimal HTTP mapping

| endpoint | does |
|---|---|
| `POST /runs` (multipart: `.mint01` files, optional prompt, backend, quick) | create `work_dir`, save uploads, enqueue a job that builds `Session` + `Mellor` and calls `run_pipeline`; return the run id |
| `GET /runs/{id}` | status from the job queue; when finished, `analysis.json`, the figure list, `spec.json`, flags |
| `GET /runs/{id}/events` (SSE) | the `on_event` stream (and `run.log` tail) |
| `GET /runs/{id}/report.pdf`, `/figures/{name}.png`, `/script.py`, `/epsr.zip` | files from the run folder |
| `POST /runs/{id}/chat` (`{"message": ...}`) | `Session.load_existing()` (or a cached `Mellor` per run), `mellor.chat(message)`, return the answer and the new events |

A FastAPI sketch of the job:

```python
def job(run_dir: Path, files: list[Path], prompt: str | None, backend: str, quick: bool) -> dict:
    session = Session(work_dir=run_dir, data_files=files, quick=quick)
    agent = Mellor(make_backend(backend), session, on_event=lambda e: publish(run_dir, e))
    outcome = agent.run_pipeline(prompt)
    export.write_tables(run_dir)                               # CSV tables for the app's plots
    return {"report": str(outcome.report_path), "spec_source": outcome.spec_source,
            "message": outcome.final_message}
```

### 7.1 Uploading selected figures to an Overleaf project

Overleaf has no public upload API. The supported route is its Git bridge, which the
app can drive server-side:

* **Entitlement.** Git integration is an Overleaf premium feature (Overleaf Commons
  institutions, group subscriptions or individual subscriptions). Without it the
  button cannot work; give those users the fallback below.
* **Credentials.** Each user creates a Git authentication token in their Overleaf
  Account Settings ("Git integration"). Git operations use the username `git` with
  that token as the password; one token covers all of a user's projects, and
  collaborators can push, not only owners. The app stores the token encrypted per
  user and never logs it.
* **Project.** The remote is `https://git.overleaf.com/<project id>` (the id is the
  last path segment of the project's URL). Let the user paste the project URL once.
* **Flow per upload.** In a temporary directory: `git clone --depth 1
  https://git:<token>@git.overleaf.com/<id>`; copy the chosen PNGs into `figures/`
  (create it if absent; keep Mellor's file names, prefixed with the run name to avoid
  collisions, e.g. `waterbox_02_fq_fit.png`); `git add`, `git commit`, `git push`.
  If the push is rejected, `git pull --rebase` once and push again. Delete the clone.
  Show the user a ready-to-paste snippet:

  ```latex
  \begin{figure}\centering
  \includegraphics[width=\linewidth]{figures/waterbox_02_fq_fit.png}
  \caption{<caption from mellor.figures.FIGURE_CAPTIONS>}\label{fig:waterbox_02_fq_fit}
  \end{figure}
  ```
* **Fallback without premium.** Offer a zip of the selected figures plus a
  `figures.tex` file with the snippets, for drag-and-drop upload in the Overleaf
  editor. (Overleaf's GitHub synchronisation is a separate premium feature; if a
  user already links the project to GitHub, pushing to that repository's `figures/`
  folder and pulling in Overleaf works too.)
* **Formats.** Mellor writes PNG at 150 dpi (a few hundred kB each). For vector
  figures, regenerate with the `epsrvis` plotting functions and `fig.savefig("x.pdf")`
  from the trajectory; `mellor.figures.make_figures(result, dpi=300)` gives sharper PNGs.

### 7.2 Users, authentication and per-user state

Mellor keeps no user state: everything about a run is in its folder. For a
multi-user app:

* one run folder per run under a per-user root (`runs/<user id>/<run id>/`); never let
  a request address a folder outside the caller's root (the model's tools can read
  any file inside the run folder, so the folder is the trust boundary);
* the Hugging Face token can be one server-side key with per-user quotas on runs
  and chat turns (a run is well under 5 p of tokens), or a per-user key stored
  encrypted; set `HF_TOKEN` per worker accordingly;
* the transcript (`mellor_transcript.jsonl`), the prompt and the uploads are the
  user's data: scope them to the account, and delete the folder when the user does;
* cap concurrent jobs per user (one EPSR run uses every core for minutes) and total
  disk (a few hundred MB per run with the trajectory; `--no-record` keeps only the
  EPSR files but loses the per-iteration figures).

### 7.3 Material for the About (theory) and Instructions pages

**Theory, in short.** EPSR (Empirical Potential Structure Refinement, A. K. Soper)
builds an atomistic model of a liquid or amorphous material that agrees with
measured diffraction data. A Monte Carlo simulation of a periodic box of molecules
runs on a *reference potential* (Lennard-Jones plus Coulomb). After every few cycles
the simulated structure factor F(Q) is compared with each measured data set, and the
difference is used to refine an *empirical potential* (EP), pair by pair, that pushes
the simulation towards the data; its amplitude is capped by `ereq`. Several isotopic
samples of the same material (H₂O, D₂O and an H/D mixture, for instance) weight the
same partial structure factors differently, which is what makes the problem well
determined. The outputs are configurations consistent with the data and the partial
pair distribution functions g(r), from which peak positions, coordination numbers
and bond-angle distributions follow. Longer explanations in this repository:
`python/docs/user_guide.md` §1 and §2 (the workflow and the files), `validation/README.md`
(the three validated systems on real data and their numbers), and
`visualisation/python_visualisation_toolkit/README.md` §6 and §7 (what each figure
shows; definitions of the R-factor, R², χ², coordination number and tetrahedral order).

**References.**
* A. K. Soper, "Empirical potential Monte Carlo simulation of fluid structure",
  *Chem. Phys.* **202**, 295–306 (1996).
* A. K. Soper, "Partial structure factors from disordered materials diffraction data:
  an approach using empirical potential structure refinement", *Phys. Rev. B* **72**,
  104204 (2005).
* A. K. Soper, "Empirical Potential Structure Refinement: EPSRshell, a user's guide",
  Rutherford Appleton Laboratory Technical Report RAL-TR-2011-012 (2011): the manual
  of the original Fortran program that the Rust engine reproduces.
* A. K. Soper, "GudrunN and GudrunX: programs for correcting raw neutron and X-ray
  diffraction data to differential scattering cross section", RAL-TR-2011-013 (2011):
  the program that produces `.mint01` files.
* T. G. A. Youngs, "Dissolve: next generation software for the interrogation of total
  scattering data by empirical potential generation", *Mol. Phys.* **117**, 3464–3477
  (2019): the source of the example systems and data.

**What this software is.** The engine is a Rust translation of the Fortran EPSR
suite (`rust/`), verified byte for byte against the original on the test case and on
three real systems (water on SANDALS, amorphous silica and the Cu-BTC framework on
NIMROD). `epsr` (`python/`) is its Python interface; `epsrvis` (`visualisation/`) the
analysis and plotting toolkit; Mellor the agent on top. The language model only
proposes the model system and writes prose; every number comes from the engine.

**Instructions page: what users must know.**
* *Input.* Gudrun `.mint01` files: three columns, Q in Å⁻¹, differential cross
  section (F(Q)) in barn sr⁻¹ atom⁻¹, and its error; header lines start with `#` and
  the second one is the sample title, which Mellor reads to recognise the sample
  (H₂O/D₂O/HDO, SiO₂, benzene, DMSO, Cu-BTC). Upload all isotopic samples of one
  material together; add a sentence about the sample (composition, temperature,
  density) if the titles do not say it.
* *What happens.* Mellor inspects the data, picks or adapts a validated model from its
  library (water, silica, DMSO/water, benzene, Cu-BTC), generates the EPSR inputs and a
  Python script, and runs: a scratch start, an equilibration on the reference
  potential, then refinement runs with the empirical potential. The quick schedule
  (10 iterations per run) takes a minute or two and is for previews; the full schedule
  (99 + 2 × 100 iterations) takes about ten minutes for water on a 20-core machine.
* *Reading the results.* The R-factor is the mean squared difference between data and
  simulation over the fitted Q range: for a liquid refined against good neutron data
  it falls to a few 10⁻⁴ (water: 2.3 × 10⁻⁴); a crystal modelled with intermolecular
  bonds only stays near 10⁻². R² near 1 means the fit follows the data. "EP amplitude
  at ereq" means the empirical potential is as large as allowed: a larger `ereq` may
  improve the fit, but a large EP means the reference potential is poor. The pressure
  should be modest (a few kbar) for a liquid at its measured density. Partial g(r):
  first-peak positions are bond or neighbour distances (water O···H 1.8 Å, O···O 2.8 Å),
  coordination numbers count neighbours within the first minimum.
* *Checking it yourself.* Every run folder contains `run_<name>.py`, a standalone
  Python script with every input file embedded; `python run_<name>.py --out DIR`
  reproduces the run with the open-source wrapper and toolkit, no language model
  involved. The report's appendix lists the same files.
* *Limits.* Models outside the library are the language model's proposal and are
  marked unvalidated in the report; force-field parameters must then be checked. EPSR
  runs with EPSR's defaults (`nq 600`, `qstep 0.05`, `ntimes 5`, `feedback 0.9`).

## 8. Command-line equivalents (useful for smoke tests)

```bash
uv run mellor run data/*.mint01 --out runs/0042 --backend hf-api --model Qwen/Qwen3-235B-A22B-Instruct-2507 --no-chat
uv run mellor run data/*.mint01 --out runs/0043 --system water --quick --backend scripted --no-chat   # no model
uv run mellor chat runs/0042
uv run mellor report runs/0042
uv run mellor search "benzene neutron" --sources local,dissolve,zenodo
uv run mellor info
```

## 9. Where to read more

* `docs/architecture.md`: the design, the data flow, the run-folder and manifest
  contracts, the agent's phases.
* `src/mellor/spec.py`: the specification model and its JSON schema.
* `src/mellor/llm/tools.py`: the tools the model can call and the `Session` fields.
* `src/mellor/llm/backend.py`: the backends and the message format.
* `tests/test_pipeline.py`: a complete scripted run, the best template for an
  integration test of the app.
