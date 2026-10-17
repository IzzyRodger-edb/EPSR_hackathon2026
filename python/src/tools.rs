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
pub fn run_tool(
    py: Python<'_>,
    folder: PathBuf,
    command: &str,
    args: Vec<String>,
    answers: String,
    verbose: bool,
) -> PyResult<()> {
    let command = command.to_ascii_lowercase();
    if !COMMANDS.contains(&command.as_str()) {
        return Err(PyValueError::new_err(format!(
            "unknown EPSR tool {command:?}; valid tools: {}",
            COMMANDS.join(", ")
        )));
    }
    py.detach(move || {
        let standard_input = TextFile::from_string(&answers);
        let mut context = ToolContext::new(&folder, args, standard_input, verbose);
        core_run_tool(&command, &mut context)
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
    symtab_tools::ALL.iter().map(|table| table.name).collect()
}

/// The keyword table of a keyword-file driven program (e.g. `nwts`, `partials`).
///
/// Returns a dict with `name`, `extinp`, `prompt`, `dataprompt`,
/// `dataprompttext`, `symbols` (list of `(symbol, type, iflag, default,
/// description)`) and `marks` (dict, e.g. `nsymboluptodata`: number of
/// keywords before the first repeated `dataprompt` block).
#[pyfunction]
pub fn tool_keywords<'py>(py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyDict>> {
    let table = symtab_tools::ALL
        .iter()
        .find(|table| table.name.eq_ignore_ascii_case(name))
        .ok_or_else(|| {
            PyValueError::new_err(format!(
                "no keyword table for {name:?}; tables exist for: {}",
                tool_tables().join(", ")
            ))
        })?;
    let description = PyDict::new(py);
    description.set_item("name", table.name)?;
    description.set_item("extinp", table.extinp)?;
    description.set_item("prompt", table.prompt)?;
    description.set_item("dataprompt", table.dataprompt)?;
    description.set_item("dataprompttext", table.dataprompttext)?;
    let symbols = PyList::empty(py);
    for (symbol, value_type, flag, default, symbol_description) in table.symbols {
        symbols.append((
            *symbol,
            value_type.to_string(),
            *flag,
            *default,
            *symbol_description,
        ))?;
    }
    description.set_item("symbols", symbols)?;
    let marks = PyDict::new(py);
    for (mark_name, position) in table.marks {
        marks.set_item(*mark_name, *position)?;
    }
    description.set_item("marks", marks)?;
    Ok(description)
}
