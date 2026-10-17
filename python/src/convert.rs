//! Conversions between engine types and Python objects.

use epsr_core::EpsrError;
use numpy::{PyArray1, PyArray2, PyArrayMethods};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;

/// Map an engine error to the Python exception classes defined in
/// `epsr/errors.py` (`EpsrIOError(OSError)`, `EpsrInputError(ValueError)`,
/// `EpsrUnsupportedError(NotImplementedError)`, all subclasses of `EpsrError`).
///
/// The classes live in Python so that they can use multiple inheritance with the
/// built-in exception types.
pub fn to_py(error: EpsrError) -> PyErr {
    let class_name = match &error {
        EpsrError::Io(_) => "EpsrIOError",
        EpsrError::Input(_) => "EpsrInputError",
        EpsrError::Unsupported(_) => "EpsrUnsupportedError",
    };
    let message = error.to_string();
    Python::attach(|py| {
        let exception_class = py
            .import("epsr.errors")
            .and_then(|errors_module| errors_module.getattr(class_name));
        match exception_class {
            Ok(exception_class) => match exception_class.call1((message.clone(),)) {
                Ok(exception) => PyErr::from_value(exception),
                Err(construction_error) => construction_error,
            },
            // The pure-Python package is missing (extension imported on its own).
            Err(_) => PyRuntimeError::new_err(message),
        }
    })
}

/// A 1-D float32 array copied from a slice.
pub fn array_1d<'py>(py: Python<'py>, values: &[f32]) -> Bound<'py, PyArray1<f32>> {
    PyArray1::from_slice(py, values)
}

/// A 2-D float32 array `(rows.len(), n_columns)` copied from equally long rows.
/// An empty row list gives shape `(0, n_columns)`.
pub fn array_2d<'py>(
    py: Python<'py>,
    rows: Vec<&[f32]>,
    n_columns: usize,
) -> PyResult<Bound<'py, PyArray2<f32>>> {
    let mut all_values = Vec::with_capacity(rows.len() * n_columns);
    for row in &rows {
        if row.len() != n_columns {
            return Err(PyValueError::new_err(format!(
                "internal error: row of length {} (expected {n_columns})",
                row.len()
            )));
        }
        all_values.extend_from_slice(row);
    }
    PyArray1::from_vec(py, all_values).reshape([rows.len(), n_columns])
}
