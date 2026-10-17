//! Fortran compatibility layer.
//!
//! These modules contain no EPSR physics: they reproduce the parts of the
//! Fortran language and runtime (gfortran) whose behaviour EPSR's results
//! depend on, so that the physics modules can be translated statement by
//! statement and still produce byte-identical output files.
//!
//! | module        | purpose |
//! |---------------|---------|
//! | [`array`]     | 1-based, column-major `F1`/`F2`/`F3` arrays |
//! | [`fmath`]     | integer powers, `NINT`, `SIGN`, ... with gfortran's exact rounding |
//! | [`format`]    | FORMAT interpreter and list-directed output (`write(u,fmt)`, `write(u,*)`) |
//! | [`listread`]  | list-directed input, sequential text files, `backspace` |
//! | [`text`]      | EPSR's `parse`/`concatenate` word splitter and character helpers |
//! | [`output`]    | buffered output files (`open(status='unknown')` + `write`) |

pub mod array;
pub mod fmath;
pub mod format;
pub mod listread;
pub mod output;
pub mod text;

pub use array::{F1, F2, F3};
pub use format::{fd, fi, fmt, fr, ListWriter, FV};
pub use listread::{IoErr, IoRes, ListRead, TextFile};
pub use text::{concatenate, fstr, parse, Words};
