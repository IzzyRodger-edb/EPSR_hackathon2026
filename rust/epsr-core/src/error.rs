//! Error type for the EPSR engine.

use std::fmt;

#[derive(Debug, Clone)]
pub enum EpsrError {
    /// A file could not be opened/read/written.
    Io(String),
    /// An input file was malformed (the Fortran would report a READ error).
    Input(String),
    /// The input requests a feature of the Fortran code that has not been
    /// translated yet (see `docs/FORTRAN_MAP.md` for the list).
    Unsupported(String),
}

impl fmt::Display for EpsrError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EpsrError::Io(s) => write!(f, "I/O error: {s}"),
            EpsrError::Input(s) => write!(f, "input error: {s}"),
            EpsrError::Unsupported(s) => write!(f, "not supported by the Rust port: {s}"),
        }
    }
}

impl std::error::Error for EpsrError {}

impl From<std::io::Error> for EpsrError {
    fn from(e: std::io::Error) -> Self {
        EpsrError::Io(e.to_string())
    }
}

impl From<crate::fortran::IoErr> for EpsrError {
    fn from(e: crate::fortran::IoErr) -> Self {
        match e {
            crate::fortran::IoErr::Eof => EpsrError::Input("unexpected end of file".into()),
            crate::fortran::IoErr::Bad(s) => EpsrError::Input(s),
        }
    }
}

pub type Result<T> = std::result::Result<T, EpsrError>;

/// Shorthand for "this branch of the Fortran is not ported".
pub fn unsupported<T>(what: &str) -> Result<T> {
    Err(EpsrError::Unsupported(what.to_string()))
}
