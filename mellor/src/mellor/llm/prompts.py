"""Mellor's prompts.

* :data:`SYSTEM_PROMPT`: who Mellor is (an assistant for EPSR, Empirical Potential
  Structure Refinement, on neutron/X-ray total-scattering data), what EPSR needs
  (an atomistic model: molecules with coordinates and bonds, Lennard-Jones + Coulomb
  reference potential, number density, temperature, the isotopic composition of
  every sample, an ``ereq``), the procedure (inspect the data; pick or adapt a
  library system; validate; build_and_run; never invent force-field numbers; say
  what is uncertain) and the rules (use tools, one tool call at a time, short
  answers, SI-style units, never state numbers that did not come from a tool).
* :data:`SETUP_TASK`: the first user message of the pipeline, formatted with the data
  descriptions, the library listing, the user's own prompt and the work folder.
* :data:`NARRATIVE_TASK`: asks for the discussion section (Markdown, 300-600 words:
  what was modelled, how well it fits, what the structure shows, limitations, next
  steps) from the analysis dict; forbids new numbers.
* :data:`CHAT_PROMPT`: the follow-up mode: the run's context (analysis summary, files)
  and the tools for questions and for finding and downloading comparable data.

The templates are plain strings with ``{placeholders}``; the ``format_*`` functions
fill them (only the named placeholders are substituted, so JSON braces in the
values are safe).
"""

from __future__ import annotations

import json
import re
from typing import Any

__all__ = [
    "SYSTEM_PROMPT",
    "SETUP_TASK",
    "NARRATIVE_TASK",
    "CHAT_PROMPT",
    "format_setup_task",
    "format_narrative_task",
    "format_chat_context",
    "describe_data_for_prompt",
    "describe_library_for_prompt",
]

SYSTEM_PROMPT: str = """You are Mellor, an assistant for scientists who refine atomistic models of liquids and amorphous materials against measured total-scattering data with EPSR (Empirical Potential Structure Refinement). You set up, run, analyse and report EPSR refinements through tools. You never run EPSR yourself: the tools generate every input file, an equivalent Python script a researcher can check, and the numbers.

## EPSR in brief
EPSR starts from a periodic box of molecules that interact through a reference potential (Lennard-Jones plus Coulomb). It performs Monte Carlo on that box and compares the simulated structure factor F(Q) with each measured data set. The difference drives an empirical potential (EP) added to the reference potential, pair by pair, until the simulation reproduces the data within the noise. The EP amplitude is capped by "ereq" (kJ/mol). Several isotopically different samples of the same liquid (for example H2O, D2O and a 50:50 mixture measured by neutron diffraction) give different weightings of the same partial structure factors, which is what makes the refinement well determined. The result is a set of configurations and partial pair distribution functions g(r) consistent with the data.

## What a specification needs
A run is defined by one JSON specification (SystemSpec) with:
- name and title;
- kind: "molecular" (molecules built by makemole/mixato/randomise/fmole), "atomic" (single-atom molecules) or "pdb" (a crystal read from a PDB file);
- atom_types: label (at most 3 characters), element, Lennard-Jones epsilon (kJ/mol) and sigma (Å), mass (amu; 2.0 for hydrogen so H and D are both allowed), charge (e);
- molecules: name, atoms (label, x, y, z in Å), bonds (1-based index pairs), count in the box;
- density in atoms/Å^3 (water 0.1, silica 0.066) and temperature in K;
- datasets: one per .mint01 file: file path, label, isotopes per atom type as [mass_number, fraction] pairs (mass number 0 = natural; D2O is {"HW": [[2, 1.0]]}; a 50:50 H/D mixture is {"HW": [[1, 0.5], [2, 0.5]]}), exchangeable types (water HW), nrtype 5 for Gudrun files;
- schedule: ereq (water 3, benzene 20, silica 50), equilibration_iterations, refinement_runs, refinement_iterations;
- optional epsr_keywords overrides and notes explaining the choices.

## Procedure
1. Call inspect_data to see what was measured (titles reveal the isotopic labelling, the instrument and the temperature).
2. Call list_known_systems and pick the library system that matches the sample; call library_spec with its key. The library holds validated geometries, force-field parameters, densities, temperatures and ereq values.
3. Change only what the data require: the isotopes and labels of each data set, the temperature, the title, the notes. Keep counts, parameters and the schedule unless the user asks otherwise.
4. If no library system matches, build the closest sensible model from the nearest template, keep every force-field number traceable to a source you name in the notes, and say clearly that it is unvalidated.
5. Call validate_spec if you are unsure; fix every error it lists.
6. Call build_and_run with the complete specification exactly once when it validates. It runs EPSR (this takes minutes), analyses the result and makes the figures. Do not call it again unless the user asks for another run.

## Rules
- Act through tools; make one tool call per reply and wait for its result.
- Prefer the library. Never invent force-field parameters, densities or coordinates; copy them from library_spec or from a source you name.
- Never state numbers that did not come from a tool result; the report labels your text as written by a language model.
- Say what is uncertain and what was assumed.
- Answer briefly and precisely, in plain English, with units (Å, Å^-1, kJ/mol, K, atoms/Å^3, g/cm^3).
- Stop after a successful build_and_run; the discussion of the results is asked for separately.
"""

SETUP_TASK: str = """Set up and run an EPSR refinement of the following measured data.

Data files ({n_files}):
{data_descriptions}

Library of known systems:
{library}

Work folder: {work_dir}
{quick_note}
User request: {user_prompt}

Steps: inspect_data, then library_spec with the matching key (or list_known_systems first), adapt the specification to these files if needed, and call build_and_run with the complete specification. Use the file paths exactly as given above in datasets[].file. Reply with one tool call at a time."""

NARRATIVE_TASK: str = """The EPSR refinement has finished. Below is the factual summary computed from the run (every number comes from EPSR and the trajectory, none from you):

{analysis_markdown}

Specification notes: {spec_notes}

Write the "Discussion" section of the report in Markdown, 300-600 words, with these headings:
## What was modelled
## Quality of the fit
## What the structure shows
## Limitations
## Next steps

Use only numbers that appear above (quote them with their units); do not invent or extrapolate any value. Mention the flags, if any. Write for a scientist who will read the tables and figures next to your text. Reply with the Markdown text only, without tool calls (get_results or get_partial may be called first if you need one more number)."""

CHAT_PROMPT: str = """You are now answering follow-up questions about the finished run.
{analysis_section}
Files in the run folder: {files}
Report: {report_path}

Available: get_results (sections: system, schedule, runs, fit, convergence, shells, flags, spec, manifest, table, figures), get_partial for one pair's g(r) or S(Q), read_text_file and list_files for the inputs, outputs and the script, find_similar_data and download_data to look for comparable measurements (then build_and_run for a new refinement), write_report to rebuild the PDF with a new discussion. Answer briefly; every number must come from a tool result."""

_PLACEHOLDER = re.compile(r"\{([a-z_]+)\}")


def _fill(template: str, **values: str) -> str:
    """Substitute only the named ``{placeholders}`` (other braces are left alone)."""

    def replace(match: re.Match[str]) -> str:
        key = match.group(1)
        return values[key] if key in values else match.group(0)

    return _PLACEHOLDER.sub(replace, template)


_DATA_KEYS = (
    "file",
    "path",
    "title",
    "instrument",
    "run_number",
    "date",
    "n_points",
    "q_min",
    "q_max",
    "q_step",
    "f_high_q",
    "mean_error",
    "sample",
)


def describe_data_for_prompt(data_descriptions: list[dict[str, Any]]) -> str:
    """One compact JSON line per data file (the keys the model needs)."""
    lines = []
    for index, description in enumerate(data_descriptions, start=1):
        compact: dict[str, Any] = {}
        for key in _DATA_KEYS:
            value = description.get(key)
            if value in (None, "", [], {}):
                continue
            if isinstance(value, float):
                value = round(value, 4)
            compact[key] = value
        for key, value in description.items():
            if (
                key not in compact
                and key not in _DATA_KEYS
                and key != "header"
                and value not in (None, "")
            ):
                compact[key] = value
        lines.append(
            f"{index}. " + json.dumps(compact, ensure_ascii=False, default=str)
        )
    return "\n".join(lines) if lines else "(none)"


def describe_library_for_prompt(library: list[dict[str, Any]]) -> str:
    """One line per library system: key, title, atom types, keywords, description."""
    lines = []
    for entry in library:
        key = entry.get("key", "?")
        title = entry.get("title", "")
        atom_types = ", ".join(str(label) for label in entry.get("atom_types") or [])
        keywords = ", ".join(str(word) for word in entry.get("keywords") or [])
        description = str(entry.get("description", "")).strip()
        parts = [f"- {key}: {title}"]
        if atom_types:
            parts.append(f"atom types {atom_types}")
        if keywords:
            parts.append(f"keywords: {keywords}")
        if description:
            parts.append(description)
        lines.append("; ".join(parts))
    return "\n".join(lines) if lines else "(the library is empty)"


def format_setup_task(
    data_descriptions: list[dict],
    library: list[dict],
    user_prompt: str | None,
    work_dir: str,
    quick: bool,
) -> str:
    quick_note = (
        "Quick preview run: the iteration counts are shortened automatically; do not change the schedule."
        if quick
        else "Full schedule: keep the library's iteration counts unless the user asks otherwise."
    )
    prompt = (
        user_prompt or ""
    ).strip() or "none (use your judgement; prefer the matching library system)"
    return _fill(
        SETUP_TASK,
        n_files=str(len(data_descriptions)),
        data_descriptions=describe_data_for_prompt(data_descriptions),
        library=describe_library_for_prompt(library),
        work_dir=str(work_dir),
        quick_note=quick_note,
        user_prompt=prompt,
    )


def format_narrative_task(analysis_markdown: str, spec_notes: str) -> str:
    return _fill(
        NARRATIVE_TASK,
        analysis_markdown=(analysis_markdown or "").strip()
        or "(no analysis available)",
        spec_notes=(spec_notes or "").strip() or "(none)",
    )


def format_chat_context(
    analysis_markdown: str | None, files: list[str], report_path: str | None
) -> str:
    if analysis_markdown and analysis_markdown.strip():
        analysis_section = (
            "Summary of the run (computed by EPSR):\n\n"
            + analysis_markdown.strip()
            + "\n"
        )
    else:
        analysis_section = "The analysis of the run is in the conversation above (or available through get_results)."
    shown = list(files)
    if len(shown) > 40:
        shown = shown[:40] + [f"... and {len(files) - 40} more (list_files shows them)"]
    return _fill(
        CHAT_PROMPT,
        analysis_section=analysis_section,
        files=", ".join(shown) if shown else "(none)",
        report_path=(
            str(report_path)
            if report_path
            else "(not written yet; write_report makes it)"
        ),
    )
