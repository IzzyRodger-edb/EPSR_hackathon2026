//! `epsr._core`: PyO3 bindings for the Rust EPSR engine (`epsr-core`).
//!
//! This layer is deliberately thin: it mirrors the engine API with the Fortran
//! names and converts arrays to float32 numpy copies. The user-facing API is the
//! pure-Python package `epsr` (see `python/epsr/`).

use pyo3::prelude::*;

mod convert;
mod sim;
mod tools;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<sim::Simulation>()?;
    m.add_function(wrap_pyfunction!(tools::run_tool, m)?)?;
    m.add_function(wrap_pyfunction!(tools::tool_commands, m)?)?;
    m.add_function(wrap_pyfunction!(tools::tool_tables, m)?)?;
    m.add_function(wrap_pyfunction!(tools::tool_keywords, m)?)?;
    m.add("VERSION_DATE", epsr_core::VERSION_DATE)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
