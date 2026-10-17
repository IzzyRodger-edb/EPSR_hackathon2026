# Mellor: architecture

Mellor is an LLM agent that turns measured structure-factor data (Gudrun `.mint01`
files) into a finished EPSR refinement: it proposes the atomistic model, writes the
EPSR input files **and** an equivalent standalone Python script, runs the Rust EPSR
engine through the `epsr` Python wrapper, analyses and plots the results with the
`epsrvis` toolkit, and assembles a PDF report. Afterwards it answers questions about
the run and can look for comparable experimental data.

The language model comes from Hugging Face (default `Qwen/Qwen3-4B-Instruct-2507`,
run locally with `transformers`; the Hugging Face Inference API is the alternative).

## 1. Principles

1. **The model decides, the code executes.** The LLM only produces a *system
   specification* (`SystemSpec`, JSON) and prose. Every file EPSR reads is generated
   deterministically from the specification by `mellor.inputs`, so a researcher can
   check the specification, the generated files and the script, not a transcript.
2. **The Python script is the single source of truth.** `mellor.script` renders a
   self-contained `run_<name>.py` that embeds the generated input files verbatim and
   drives the wrapper (`epsr`, `epsrvis`) step by step. Mellor *executes that script*
   (`mellor.runner`), so what the researcher can re-run independently is exactly what
   produced the results.
3. **Numbers come from the engine, never from the model.** Tables, metrics and figures
   are computed by `mellor.analysis`/`mellor.figures` from the recorded trajectory and
   the EPSR output files. The LLM's narrative is clearly labelled as such in the report.
4. **Degrade gracefully.** Without a GPU, the Inference API or a scripted backend can be
   used; without any model, `mellor run --system water` builds the specification from
   the built-in library. If the model fails to produce a valid specification, Mellor
   falls back to the library and says so in the report.
5. **Additive.** Nothing in `python/`, `rust/` or `visualisation/` is modified
   (see the repository's rules); `mellor/` depends on them as path packages.

## 2. Data flow

```
.mint01 files ──► data.read_mint01 ──► Mint01 descriptions ─┐
                                                            ▼
            library (known systems) ◄──── LLM (Mellor) ──► SystemSpec (JSON)
                                                            │ spec.validate()
                                                            ▼
                           inputs.generate(spec) ──► EpsrInputs (.mol, mixato answers,
                                                      .NWTS.dat, .EPSR.inp, pdb answers)
                                                            │
                           script.render_script(spec, inputs) ──► run_<name>.py
                                                            │
                           runner.run_spec ── subprocess ──► epsr/ (simulation folder)
                                                            ├─► trajectory/ (epsrvis format)
                                                            └─► manifest.json
                                                            │
                   analysis.analyse(result) ──► analysis dict (fit quality, shells, ...)
                   figures.make_figures(result) ──► figures/*.png (epsrvis.make_report + extras)
                                                            │
                   LLM narrative (discussion markdown) ─────┤
                                                            ▼
                   report.build_pdf(...) ──► report.pdf
                                                            │
                   chat: Mellor + tools (get_results, find_similar_data, download_data, run again)
```

## 3. Run folder layout

`mellor run` writes everything into one folder (`--out`, default `work/<name>`):

| path | written by | content |
|---|---|---|
| `data/` | runner | copies of the input `.mint01` files |
| `spec.json` | runner | the `SystemSpec` that was run |
| `run_<name>.py` | script | the standalone Python script (the "Python input file") |
| `epsr/` | the script | the EPSR simulation folder: `.mol`, `.ato`, `.NWTS.dat`, `.wts`, `.EPSR.inp`, `.pcof`, all outputs |
| `trajectory/` | the script | per-iteration record in the `epsrvis` folder format |
| `manifest.json` | the script | per-run summaries (see §6), timings, versions |
| `run.log` | runner | console output of the script |
| `analysis.json` | analysis | the analysis dict |
| `figures/` | figures | PNG figures, `summary.csv`, `summary.md` |
| `report.pdf` | report | the PDF report |
| `mellor_transcript.jsonl` | agent | every model message and tool call |

## 4. Modules (`src/mellor/`)

| module | owner of | depends on |
|---|---|---|
| `spec.py` | `SystemSpec` and its parts, validation, JSON (+ JSON schema for the LLM) | – |
| `data.py` | `.mint01` reading and description, isotope hints from titles, raw-data figure | numpy, matplotlib |
| `library.py` | built-in systems (water, silica, DMSO/water, benzene, argon, Cu-BTC) → `SystemSpec` | spec, data |
| `inputs.py` | `SystemSpec` → EPSR input texts (`EpsrInputs`), `.mol` / `.NWTS.dat` / `.EPSR.inp` renderers | spec, epsr.tools (keyword tables) |
| `script.py` | `run_<name>.py` renderer | spec, inputs |
| `runner.py` | executes the script in a subprocess, `RunResult`, manifest loading | spec, inputs, script |
| `analysis.py` | numbers from the trajectory and output files | epsrvis.metrics, runner |
| `figures.py` | figures via `epsrvis.make_report`, `epsrvis.plots.final_outputs`, data overlay | epsrvis, data, runner |
| `report.py` | PDF (reportlab) | runner, analysis, figures |
| `sources.py` | similar-data search: local repository, Dissolve examples (GitHub), Zenodo; download | urllib |
| `llm/backend.py` | `ChatBackend` protocol; `TransformersBackend`, `HFInferenceBackend`, `ScriptedBackend`; `make_backend` | transformers/torch (optional), huggingface_hub |
| `llm/tools.py` | `Tool`, `ToolRegistry`, `Session`, the tool functions exposed to the model | everything above |
| `llm/prompts.py` | system prompts | – |
| `llm/agent.py` | `Mellor`: the tool-call loop, `run_pipeline`, `chat` | backend, tools |
| `cli.py` | `mellor run | chat | report | spec | search` | all |

## 5. The specification (`SystemSpec`)

See `spec.py` for the exact fields. In short:

* `atom_types`: EPSR atom types (label of at most 3 characters, element, Lennard-Jones
  ε [kJ/mol] and σ [Å], mass [amu], charge [e]).
* `molecules`: templates (name, atoms with coordinates and type labels, bonds with
  1-based indices, number of molecules in the box). Angles are every pair of bonded
  neighbours around an atom; values come from the coordinates unless given.
* `kind`: `molecular` (makemole → mixato → randomise → fmole), `atomic` (single-atom
  molecules: no fmole) or `pdb` (readpdb + changeato on a PDB file, crystal positions kept).
* `density` [atoms/Å³], `temperature` [K].
* `datasets`: one per `.mint01` file: isotopic composition of every atom type
  (`[mass_number, fraction]` pairs, mass number 0 = natural), exchangeable types, and
  the per-dataset `.EPSR.inp` keywords (`nrtype 5`, `rshmin`, `tweak`, ...).
* `schedule`: scratch run (1 iteration), equilibration on the reference potential
  (`equilibration_iterations`, `ereq` 0), then `refinement_runs` runs of
  `refinement_iterations` with `ereq` set; `ntimes`, `fmole_iterations`.
* `epsr_keywords`: overrides of the `.EPSR.inp` template (`qmin`, `fwhmq`, `qwidthqmax`, ...).

The `.EPSR.inp` template is `validation/water/waterbox.EPSR.inp` (the EPSR test-case
values: `nq 600`, `qstep 0.05`, `ntimes 5`, `feedback 0.9`, `fwhmq 0.02`, ...). The
`.mol` format and the `nwts` keyword file are documented in
`validation/PYTHON_WRAPPER_WORKFLOW.md` §2. `makemole` truncates labels to 3 characters
and reads the four trailing values of a `potential` line as `nsph rmin rmax cg`
(defaults `1 0 0.5σ 1`).

## 6. The manifest (`manifest.json`, written by the generated script)

```json
{
  "name": "waterbox", "title": "...", "kind": "molecular", "quick": false,
  "created": "2026-10-17T14:02:11", "python": "3.12.3",
  "versions": {"epsr": "0.1.0", "epsrvis": "0.1.0", "mellor": "0.1.0"},
  "epsr_folder": "epsr", "trajectory": "trajectory", "inp_name": "waterbox.EPSR.inp",
  "box_side": 31.0723, "n_atoms": 3000, "n_molecules": 1000,
  "atom_types": ["OW", "HW"], "pairs": ["OW-OW", "OW-HW", "HW-HW"],
  "datasets": [{"label": "D2O", "file": "SLS18502-D2O.mint01", "wts": "D2O.NWTStot.wts"}],
  "setup_steps": ["makemole Water", "mixato ...", "randomise waterbox", "fmole waterbox 100", "nwts ..."],
  "runs": [
    {"run": 1, "stage": "scratch", "niter": 1, "ereq": 0.0, "duration_s": 1.3,
     "summary": {"r_factor": 0.021, "r_factor_per_dataset": [..], "energy": -12.3,
                 "pressure": 1.2, "ep_amplitude": 0.0, "energy_inter": .., "energy_intra": .., "chisq": ..}},
    {"run": 2, "stage": "reference potential", ...},
    {"run": 3, "stage": "EP refinement", ...}
  ],
  "final": {"r_factor": ..., "...": "the last run's summary again"}
}
```

## 7. The agent

`Mellor.run_pipeline(data_files, prompt)`:

1. **Model setup (LLM).** The model gets the data descriptions and the library listing
   and must call `build_and_run(spec)` with a valid `SystemSpec`. `validate_spec`
   returns precise errors so it can correct itself. After `max_spec_attempts` failures
   Mellor uses `library.auto_spec` and records the fallback.
2. **Run (code).** `runner.run_spec` executes the generated script.
3. **Analysis and figures (code).**
4. **Narrative (LLM).** The model receives the analysis dict and writes the discussion
   (markdown). If it fails, `analysis.analysis_markdown` is used.
5. **Report (code).** `report.build_pdf`.

`Mellor.chat(message)` keeps the conversation and the same tools, plus
`get_results`, `find_similar_data`, `download_data`, `build_and_run` for a new run.

Tool calls use the OpenAI-style function-calling format in messages; the
`TransformersBackend` renders them with the model's chat template (`tools=`) and parses
`<tool_call>{...}</tool_call>` blocks (Qwen/Hermes format), with a fenced-JSON fallback.

## 8. Testing

* Unit tests for every module (`tests/`), running without a model or network.
* `tests/test_pipeline.py`: the whole pipeline on the water data with the
  `ScriptedBackend` and the quick schedule (marked `slow`, about a minute).
* Equivalence: the library's water specification regenerates the inputs of
  `validation/water` (parsed values equal), and the generated script reproduces the
  wrapper test workflow's EPSR output files.
* `gpu`/`network` marked tests exercise the real model and the data sources.
