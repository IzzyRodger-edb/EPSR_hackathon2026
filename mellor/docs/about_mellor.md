# Mellor: how it works and how to set it up

## How Mellor is set up

Mellor is a Python package that sits on top of our `epsr` Python wrapper, which runs the Rust engine, and the `epsrvis` plotting toolkit. Its guiding rule is that the model decides and the code executes.

- **The model writes only a specification.** It produces a JSON description of the system: the molecules, force-field parameters, density, temperature and isotopic make-up of each sample. It also writes prose for the report.
- **Code does everything else.** It generates every EPSR input file plus a standalone Python script that a researcher can read and rerun. Mellor runs that same script, so the checkable script is exactly what produced the results.
- **Numbers never come from the model.** Every value in the tables, figures and PDF report is computed from the engine's output.
- **The model works through twelve tools.** They let it inspect the data and start from a library of validated systems: water, silica, DMSO/water, benzene and Cu-BTC. It can then validate its specification and launch the run. Later tools query the results, search Zenodo and the Dissolve examples for comparable data, and rebuild the report.
- **Guardrails keep it honest.** Invalid specifications are rejected with clear errors. After repeated failures, Mellor falls back to the matching validated system and records why. Each phase is capped at 30 model calls, and every message and tool call is saved to a transcript.

## How the Hugging Face model is used

Mellor uses open models from [Hugging Face](https://huggingface.co), with two interchangeable backends:

- **Locally:** [Qwen3-4B-Instruct-2507](https://huggingface.co/Qwen/Qwen3-4B-Instruct-2507) runs on a laptop GPU with 12 GB of memory, through the `transformers` library. It's free and works offline once downloaded. Tool calls go through the model's own chat template and are parsed from its replies.
- **Hosted:** [Qwen3-235B-A22B-Instruct-2507](https://huggingface.co/Qwen/Qwen3-235B-A22B-Instruct-2507) runs on [Hugging Face Inference Providers](https://huggingface.co/docs/inference-providers), billed to your Hugging Face account through an access token. A full run uses a few tens of thousands of tokens.

There's also a no-model mode, which builds the specification directly from the validated library. That makes everything testable without a model.

## Setting up the Hugging Face model

You need [uv](https://docs.astral.sh/uv/) and a Hugging Face account. Run every command from the `mellor/` folder.

### Option 1: the local model, free, needs a GPU

1. Install the local-model extras:
   ```bash
   cd mellor
   uv sync --extra local
   ```
   You need an NVIDIA GPU with about 10 GB of memory. The model takes about 8 GB in bfloat16. It also runs on a CPU, but much more slowly and with about 16 GB of RAM.
2. Optionally, download the model in advance. Otherwise it downloads automatically on first use, into `~/.cache/huggingface/hub`.
   ```bash
   uv run hf download Qwen/Qwen3-4B-Instruct-2507
   ```
3. Check that Mellor can see the GPU and the model:
   ```bash
   uv run mellor info
   ```
4. Run a refinement. With the local extras installed, the default backend picks the local model:
   ```bash
   uv run mellor run ../experimental_data/water/data/*.mint01 --out work/water
   ```

### Option 2: the hosted model, paid, needs no GPU

1. On Hugging Face, create a [fine-grained access token](https://huggingface.co/settings/tokens) with the permission *Make calls to Inference Providers*.
2. Make the token available to Mellor:
   ```bash
   cd mellor
   uv sync
   export HF_TOKEN=hf_your_token_here
   ```
3. Run with the hosted backend:
   ```bash
   uv run mellor run ../experimental_data/water/data/*.mint01 --out work/water \
       --backend hf-api --model Qwen/Qwen3-235B-A22B-Instruct-2507
   ```

### In the web app

1. Install the agent:
   ```bash
   cd webapp
   uv sync --extra agent                  # hosted model only
   uv sync --extra agent --extra local    # also the local model
   ```
2. For the hosted model, either export `HF_TOKEN` before `uv run epsrweb`, or let each user paste their own token on the **Settings** page.
3. Choose the backend on the **New run** form. If you leave the model field blank, it uses Qwen3-4B locally and Qwen3-235B when hosted.

## A refined prompt to create Mellor

The same request, rewritten as a structured prompt for an AI coding agent. It folds in the design decisions above:

```text
Build "Mellor", an AI agent that turns processed neutron or X-ray total-scattering data into a
finished, independently checkable EPSR refinement.

Context
- This repository contains a Rust port of EPSR (rust/), a Python wrapper that drives it
  in-process (python/, package `epsr`) and a plotting toolkit (visualisation/, package
  `epsrvis`). Treat all three as read-only. Put Mellor in a new top-level `mellor/` package that
  depends on them.
- Input: Gudrun .mint01 files (columns Q, F(Q), error), usually several isotopic samples of one
  material. Example data is in experimental_data/water/data.

What Mellor must do
1. Read the .mint01 files and describe them to the language model.
2. Have the model propose the atomistic model as a JSON specification: molecules and geometry,
   Lennard-Jones and Coulomb parameters, density, temperature, the isotopic composition of each
   sample and the refinement schedule. Prefer adapting a library of validated systems (water,
   silica, ...) over inventing force-field values.
3. Validate the specification in code and return clear errors so the model can correct itself.
   After repeated failures, fall back to the closest validated system and record why.
4. Generate every EPSR input file from the specification in code, plus a standalone Python
   script that uses the `epsr` wrapper, so a researcher can read and rerun exactly what was done.
   Produce the results by running that script.
5. Analyse the results and plot them with `epsrvis`: fit quality per data set, convergence,
   partial g(r) and S(Q), empirical potentials and coordination numbers.
6. Build a PDF report with the specification, tables, figures and a short discussion written by
   the model.
7. Afterwards, answer follow-up questions about the run and search for comparable published
   data, for example on Zenodo, that the user could refine next.

Rules
- The model decides; the code executes. The model only writes the specification and prose, and
  acts only through tool calls.
- Every number in the tables, figures and report must come from EPSR output, never from the
  model. Label the model's text as such in the report.
- Log every model message and tool call to a transcript in the run folder.
- Keep one self-contained folder per run.

Model
- Use an open model from Hugging Face behind one interface with two backends: local, through
  `transformers` (default Qwen/Qwen3-4B-Instruct-2507, which fits a 12 GB GPU), and hosted,
  through Hugging Face Inference Providers with HF_TOKEN (for example
  Qwen/Qwen3-235B-A22B-Instruct-2507).
- Add a no-model mode that runs a validated library system directly, so the whole pipeline can
  be tested without a model.

Deliverables
- A Python API and a `mellor` command-line tool: `mellor run DATA.mint01 ...` and
  `mellor chat RUN_DIR`.
- Tests, including an end-to-end run on the water data.
- Documentation of the Python API and the run-folder layout, so that a web app can be built on
  top of Mellor.
```
