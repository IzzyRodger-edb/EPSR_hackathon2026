# Mellor

Mellor is an LLM agent (Hugging Face model) that turns measured structure-factor
data (`.mint01`) into a finished EPSR refinement: it proposes the model, generates the
EPSR inputs and an equivalent Python script, runs the Rust EPSR engine through the
`epsr` wrapper, plots with `epsrvis` and writes a PDF report, then answers questions.

The package is being completed today; this README will carry the install and run
instructions when it is. Until then:

* **Building a web app around Mellor?** Read [docs/integration.md](docs/integration.md):
  the Python API, the run-folder contract, tidy tables for your own plots, the
  Overleaf figure upload, users, and material for the About and Instructions pages.
* **Design:** [docs/architecture.md](docs/architecture.md).

## Example Usage

### 1. local model (Qwen3-4B on your GPU, free)
uv run mellor run ../experimental_data/water/data/*.mint01 --out work/water

### 2. hosted model on your Hugging Face credits
export HF_TOKEN=hf_...
uv run mellor run ../experimental_data/water/data/*.mint01 --out work/water --backend hf-api --model Qwen/Qwen3-235B-A22B-Instruct-2507

### 3. no model at all: the validated water setup, factual summary instead of a discussion
uv run mellor run ../experimental_data/water/data/*.mint01 --out work/water --system water --backend scripted
