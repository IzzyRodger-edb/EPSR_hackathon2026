//! The 35 auxiliary EPSR programs (`epsr_core::tools`) and their keyword tables.

use crate::convert::to_py;
use epsr_core::fortran::TextFile;
use epsr_core::symtab_tools;
use epsr_core::tools::{run_tool as core_run_tool, ToolContext, COMMANDS};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::path::PathBuf;

/// Run one auxiliary program in `folder`, as `epsr-rs <folder> <command> <args...>`.
///
/// `answers` is the text the Fortran program would read from standard input.
/// Files in `folder` are created or rewritten in place. Releases the GIL.
#[pyfunction]
#[pyo3(signature = (folder, command, args=Vec::new(), answers=String::new(), verbose=false))]
pub fn run_tool(py: Python<'_>, folder: PathBuf, command: &str, args: Vec<String>, answers: String, verbose: bool) -> PyResult<()> {
    let command = command.to_ascii_lowercase();
    if !COMMANDS.contains(&command.as_str()) {
        return Err(PyValueError::new_err(format!("unknown EPSR tool {command:?}; valid tools: {}", COMMANDS.join(", "))));
    }
    py.detach(move || {
        let stdin = TextFile::from_string(&answers);
        let mut ctx = ToolContext::new(&folder, args, stdin, verbose);
        core_run_tool(&command, &mut ctx)
    })
    .map_err(to_py)
}

/// Names of the auxiliary programs accepted by `run_tool`.
#[pyfunction]
pub fn tool_commands() -> Vec<&'static str> {
    COMMANDS.to_vec()
}

/// Names of the programs that have a keyword-file table (`tool_keywords`).
#[pyfunction]
pub fn tool_tables() -> Vec<&'static str> {
    symtab_tools::ALL.iter().map(|t| t.name).collect()
}

/// The keyword table of a keyword-file driven program (e.g. `nwts`, `partials`).
///
/// Returns a dict with `name`, `extinp`, `prompt`, `dataprompt`,
/// `dataprompttext`, `symbols` (list of `(symbol, type, iflag, default,
/// description)`) and `marks` (dict, e.g. `nsymboluptodata`: number of
/// keywords before the first repeated `dataprompt` block).
#[pyfunction]
pub fn tool_keywords<'py>(py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyDict>> {
    let t = symtab_tools::ALL
        .iter()
        .find(|t| t.name.eq_ignore_ascii_case(name))
        .ok_or_else(|| PyValueError::new_err(format!("no keyword table for {name:?}; tables exist for: {}", tool_tables().join(", "))))?;
    let d = PyDict::new(py);
    d.set_item("name", t.name)?;
    d.set_item("extinp", t.extinp)?;
    d.set_item("prompt", t.prompt)?;
    d.set_item("dataprompt", t.dataprompt)?;
    d.set_item("dataprompttext", t.dataprompttext)?;
    let syms = PyList::empty(py);
    for (sym, ty, iflag, default, descr) in t.symbols {
        syms.append((*sym, ty.to_string(), *iflag, *default, *descr))?;
    }
    d.set_item("symbols", syms)?;
    let marks = PyDict::new(py);
    for (k, v) in t.marks {
        marks.set_item(*k, *v)?;
    }
    d.set_item("marks", marks)?;
    Ok(d)
}
