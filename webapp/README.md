# EPSR web app

A web front end for the **Mellor** agent: upload measured structure-factor data, let the
agent set up and run the Rust EPSR engine, read the results as seaborn figures and a PDF
report, and publish chosen figures straight into the `figures/` folder of an Overleaf
project.

Built against the contract in [`../mellor/docs/integration.md`](../mellor/docs/integration.md).
Mellor has no HTTP layer of its own; this app is the thin layer over its Python API.

```
browser ── FastAPI ── job manager ── one process per run ── mellor.run_pipeline ── EPSR engine
                │                                                   │
            SQLite                                        the run folder (figures,
       (accounts, run states,                              tables, analysis.json,
        Overleaf push log)                                   report.pdf, epsr/)
```

## What it does

| the brief | how |
|---|---|
| upload your data | `.mint01` files, validated against the Gudrun contract *before* a job is queued |
| the agent runs the Rust EPSR code | `Mellor.run_pipeline` in a worker process, progress streamed by Server-Sent Events |
| results with seaborn plots | the app draws its own figures from the run's long-format CSV tables, as PNG **and** vector PDF |
| PDF summaries | Mellor's `report.pdf`, served per run |
| click figures to upload them to Overleaf | tick figures, push to `figures/` over Overleaf's Git bridge, with ready-made LaTeX blocks |
| users and authentication | accounts, signed session cookies, one run root per user, encrypted token storage |
| About and Instructions pages | the theory and the user guide, at `/about` and `/instructions` |

## Install and run

Requires Python 3.10 or newer and [uv](https://docs.astral.sh/uv/)
(`curl -LsSf https://astral.sh/uv/install.sh | sh`). uv creates `.venv/` and installs
from `uv.lock`; the dev group (pytest, httpx) is included by default.

```bash
cd webapp
uv sync
uv run epsrweb                     # http://127.0.0.1:8000
```

That gives you the whole app in **demonstration mode**. To make it real, three packages
from this repository have to be importable by the server process — the agent, the engine
and the plotting toolkit. The `agent` extra installs all three from the sibling folders
(`../mellor`, and through it `../python` and `../visualisation/python_visualisation_toolkit`):

```bash
cd webapp
uv sync --extra agent              # builds epsr, the Rust engine, on first sync
uv sync --extra agent --extra local   # ... plus torch/transformers for the local model
uv run epsrweb
```

Note that `uv sync` is exact: a later plain `uv sync` removes the agent again, so keep
passing `--extra agent` (or `--all-extras`).

Restart the server afterwards. `/healthz` reports `mellor=yes` when it worked, and the
banner at the top of every page disappears. Keep the checkout intact: Mellor resolves
`../python` and `../visualisation/python_visualisation_toolkit` relative to its own folder.

#### The engine build is the hard step

`epsr` is a Rust extension that uv builds with `maturin` during the sync. It needs a Rust
toolchain and a C compiler:

```bash
sudo apt install build-essential                                # or your distro's equivalent
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh   # then open a new shell
```

Until that build succeeds the app reports exactly which piece is missing — on its pages, in
`/healthz` and in `mellor_status()` — and keeps working in demonstration mode rather than
failing a job ten minutes in.

### Demonstration mode

Without Mellor the app still serves every page, validates uploads and exercises the whole
Overleaf path — but no simulation runs. Instead a run folder is written with the right
shape (same file names, same CSV columns, same `analysis.json` keys) and synthetic contents.
The measured F(Q) curves in it are your real uploaded data; the simulated curves,
potentials, g(r), R-factors and convergence history are not. Every such run is stamped
`demo` in the database, banner-marked on its results page and in its PDF, and its
`analysis.json` carries `"demo": true`. Nothing in it should be read as a result.

## Configuration

Everything is read from the environment once, at start-up ([`config.py`](src/epsrweb/config.py)).

| variable | default | meaning |
|---|---|---|
| `EPSRWEB_DATA_ROOT` | `webapp/var` | database, run folders and the generated secret live here |
| `EPSRWEB_SECRET_KEY` | generated into `var/secret.key` | signs session cookies and derives the key that encrypts stored tokens. **Set it in a deployment**; changing it logs everyone out and makes stored tokens unreadable |
| `EPSRWEB_WORKERS` | `1` | refinements running at once. One EPSR run uses every core, so raise this only together with `RAYON_NUM_THREADS` |
| `EPSRWEB_MAX_RUNS_PER_USER` | `1` | queued plus running jobs per account |
| `EPSRWEB_MAX_UPLOAD_BYTES` | 8 MiB | per file (a real `.mint01` is 100 kB to 1 MB) |
| `EPSRWEB_MAX_FILES` | `6` | files per run |
| `RAYON_NUM_THREADS` | unset | caps the engine's threads, set per worker process |
| `MELLOR_BACKEND` | `auto` | default backend: `auto`, `transformers`, `hf-api`, `scripted` |
| `MELLOR_MODEL` | `Qwen/Qwen3-235B-A22B-Instruct-2507` | default model id |
| `HF_TOKEN` | unset | server-wide Hugging Face token, used when a user has not stored their own |
| `EPSRWEB_ALLOW_REGISTRATION` | `true` | set false to close sign-ups |
| `EPSRWEB_COOKIE_SECURE` | `false` | set true behind HTTPS |
| `EPSRWEB_DISABLE_OVERLEAF_PUSH` | `false` | hides the push button, leaving the zip download |
| `EPSRWEB_GIT_TIMEOUT` | `120` | seconds per git invocation |
| `EPSRWEB_HOST`, `EPSRWEB_PORT` | `127.0.0.1`, `8000` | development server only |

## Publishing to Overleaf

Overleaf has no public upload API. The supported programmatic route is its **Git bridge**,
which the server drives: clone the project, copy the ticked figures into `figures/`, commit,
push; on a rejection, rebase once and retry. The user supplies a project URL and a Git
authentication token once, in Settings.

Two consequences are built into the UI rather than hidden:

- **Git integration is an Overleaf premium feature.** Accounts without it get a clear
  message and the fallback: a zip of the same figures plus a `figures.tex` of ready-made
  `\includegraphics` blocks, to drag into the editor. (A project already synchronised with
  GitHub can be fed through that repository instead.)
- **One token grants write access to every project in the account.** So it is encrypted at
  rest with a key derived from `EPSRWEB_SECRET_KEY`, never logged, never rendered back into
  a page, and never placed on a command line or in a stored remote: git receives it through
  a `GIT_ASKPASS` helper and the child process's environment only, and everything git
  prints is passed through `redact()` before it is stored or shown. Users can delete it from
  Settings, and are told to revoke it in Overleaf as well.

Figures are named `<run id>_<figure>.<ext>` in the project, so two runs of the same system
cannot overwrite each other, and PDF is preferred over PNG because these figures are meant
for a paper.

## Layout

| file | holds |
|---|---|
| [`app.py`](src/epsrweb/app.py) | the FastAPI factory, the pages, accounts and settings |
| [`runs_routes.py`](src/epsrweb/runs_routes.py) | everything under `/runs`: create, watch, read, publish, delete |
| [`mellor_adapter.py`](src/epsrweb/mellor_adapter.py) | **the only module that calls Mellor**, and the demo fallback switch |
| [`jobs.py`](src/epsrweb/jobs.py) | the queue: one `multiprocessing.Process` per run, with cancellation |
| [`worker.py`](src/epsrweb/worker.py) | what runs in that process |
| [`events.py`](src/epsrweb/events.py) | the progress stream, as a JSONL file the SSE endpoint tails |
| [`plots.py`](src/epsrweb/plots.py) | the app's seaborn figures, from `tables/*.csv` |
| [`overleaf.py`](src/epsrweb/overleaf.py) | the Git bridge and the zip fallback |
| [`storage.py`](src/epsrweb/storage.py) | run folders, upload validation, safe path handling |
| [`security.py`](src/epsrweb/security.py) | scrypt passwords, signed cookies, CSRF, token encryption |
| [`auth.py`](src/epsrweb/auth.py), [`db.py`](src/epsrweb/db.py) | accounts and SQLite |
| [`runviews.py`](src/epsrweb/runviews.py) | turning a run folder into what the results page needs |
| [`demo.py`](src/epsrweb/demo.py) | the synthetic run folder used without Mellor |

## Design decisions worth knowing

**One process per job, not a thread pool.** `run_pipeline` blocks, the engine saturates the
CPU, and Mellor has no cancellation: the documented way to stop a run is to kill its worker.
A `ProcessPoolExecutor` task cannot be interrupted, so the manager starts a process per job
and keeps the handle. "Stop" terminates it and keeps the partial folder. Runs that were in
flight when the server stopped are marked `interrupted` at the next start-up, not lost.

**Progress goes through a file.** The worker appends Mellor's `on_event` dicts to
`events.jsonl`; the SSE endpoint tails it. That crosses the process boundary, lets the
results page replay the whole log on reload, and survives a restart.

**Uploads are validated in the request, not the worker**, mirroring `mellor.data.read_mint01`
plus the EPSR-specific traps: the 5000-point `MQD` limit, data above `nq × qstep` that is
read and silently discarded, a missing zero sentinel row, and file names containing blanks
(EPSR truncates a data path at the first blank, so names are rewritten).

**The run folder is the trust boundary.** The model's tools can read any file inside the
folder they are given, so every run gets its own under `runs/<user id>/`, requests are
resolved with `resolve_within` and refused if they escape, and only the folders in
`runviews.SERVABLE_DIRS` are reachable. `trajectory/` is deliberately not served: it is the
per-iteration record and runs to hundreds of megabytes.

**Seaborn figures come from the CSV tables, never from EPSR files.** `mellor.export` writes
long-format tables for exactly this purpose; parsing `.u01`/`.g01` in the web layer would
duplicate the toolkit and rot.

## Tests

```bash
uv run pytest -q
```

51 tests, in four groups. **Units**: passwords, sessions, CSRF, token encryption, the
`.mint01` validator, path-escape refusals, Overleaf id parsing, redaction, LaTeX generation
and the zip. **The Mellor contract**: `test_mellor_source.py` parses the sibling package's
source and pins every assumption the adapter makes — signatures, return types, dataclass
fields, figure names, table names, the event types the log pane renders — so a rename fails
in a second instead of ten minutes into a run; `test_mellor_contract.py` checks the same
things by import where the engine is built. **The HTTP layer** end to end — register, sign in, upload two real water data sets from
`../experimental_data/`, wait for the job through the queue, and check the figures, the PDF,
the CSV tables, the zip with its LaTeX blocks, the chat endpoint, deletion, and that one
account cannot reach another's run. The suite runs in about two minutes and needs no
network, no GPU and no Overleaf account; nothing in it ever pushes to a remote.

## Deployment notes

- Put it behind a reverse proxy with TLS and set `EPSRWEB_COOKIE_SECURE=1`. For the SSE
  endpoint, disable response buffering (`X-Accel-Buffering: no` is already sent) and allow
  long-lived connections.
- Do **not** use `uvicorn --reload`: a reload would orphan a running refinement. The
  development server deliberately does not enable it.
- Size the machine for the engine, not the web layer: one run is minutes of every core and a
  few hundred megabytes of disk. `EPSRWEB_WORKERS=1` with `EPSRWEB_MAX_RUNS_PER_USER=1` is
  the honest default for a shared box.
- Hosted-model tokens are billed per token to whoever owns them. A full pipeline is a few
  tens of thousands of tokens — well under 5p per run with the models Mellor suggests — but
  set per-user quotas before opening registration to the world.
- The only state is `EPSRWEB_DATA_ROOT`: back up `epsrweb.db` with the `runs/` tree, or
  neither.
