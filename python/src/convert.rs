//! Conversions between engine types and Python objects.

use epsr_core::EpsrError;
use numpy::{PyArray1, PyArray2, PyArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Map an engine error to the Python exception classes defined in
/// `epsr/errors.py` (`EpsrIOError(OSError)`, `EpsrInputError(ValueError)`,
/// `EpsrUnsupportedError(NotImplementedError)`, all subclasses of `EpsrError`).
///
/// The classes live in Python so that they can use multiple inheritance with the
/// built-in exception types.
pub fn to_py(e: EpsrError) -> PyErr {
    let name = match &e {
        EpsrError::Io(_) => "EpsrIOError",
        EpsrError::Input(_) => "EpsrInputError",
        EpsrError::Unsupported(_) => "EpsrUnsupportedError",
    };
    let msg = e.to_string();
    Python::attach(|py| {
        match py.import("epsr.errors").and_then(|m| m.getattr(name)) {
            Ok(cls) => match cls.call1((msg.clone(),)) {
                Ok(inst) => PyErr::from_value(inst),
                Err(err) => err,
            },
            // The pure-Python package is missing (extension imported on its own).
            Err(_) => PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(msg),
        }
    })
}

/// A 1-D float32 array copied from a slice.
pub fn arr1<'py>(py: Python<'py>, s: &[f32]) -> Bound<'py, PyArray1<f32>> {
    PyArray1::from_slice(py, s)
}

/// A 2-D float32 array `(rows.len(), ncols)` copied from equally long rows.
/// An empty row list gives shape `(0, ncols)`.
pub fn arr2<'py>(py: Python<'py>, rows: Vec<&[f32]>, ncols: usize) -> PyResult<Bound<'py, PyArray2<f32>>> {
    let mut flat = Vec::with_capacity(rows.len() * ncols);
    for r in &rows {
        if r.len() != ncols {
            return Err(PyValueError::new_err(format!("internal error: row of length {} (expected {ncols})", r.len())));
        }
        flat.extend_from_slice(r);
    }
    let a = PyArray1::from_vec(py, flat);
    a.reshape([rows.len(), ncols])
}
