//! Buffered sequential output files.
//!
//! A Fortran `open(u, file=f, status='unknown')` followed by `write`s and a
//! `close` replaces the file contents; [`OutFile`] collects the records and
//! writes them when [`OutFile::close`] is called (or on drop).  [`OutFile::append`]
//! models a unit that was read to end-of-file and then written (as EPSR does
//! with the `.erg` file).

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

pub struct OutFile {
    path: PathBuf,
    buf: String,
    append: bool,
    closed: bool,
}

impl OutFile {
    /// Open for writing from the beginning (existing contents discarded).
    pub fn create(path: &Path) -> OutFile {
        OutFile { path: path.to_path_buf(), buf: String::new(), append: false, closed: false }
    }
    /// Open positioned after the existing contents.
    pub fn append(path: &Path) -> OutFile {
        OutFile { path: path.to_path_buf(), buf: String::new(), append: true, closed: false }
    }
    /// Write one record (a newline is added).  Multi-record strings (with
    /// embedded `'\n'`, e.g. from format reversion or `/`) are written as is.
    pub fn line(&mut self, s: &str) {
        self.buf.push_str(s);
        self.buf.push('\n');
    }
    pub fn path(&self) -> &Path {
        &self.path
    }
    pub fn close(mut self) -> std::io::Result<()> {
        self.flush_now()
    }
    fn flush_now(&mut self) -> std::io::Result<()> {
        if self.closed {
            return Ok(());
        }
        self.closed = true;
        if self.append {
            let mut f = fs::OpenOptions::new().create(true).append(true).open(&self.path)?;
            f.write_all(self.buf.as_bytes())
        } else {
            fs::write(&self.path, self.buf.as_bytes())
        }
    }
}

impl Drop for OutFile {
    fn drop(&mut self) {
        let _ = self.flush_now();
    }
}
