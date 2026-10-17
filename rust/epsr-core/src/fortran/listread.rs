//! Emulation of Fortran *input*: list-directed reads (`read(u,*) a,b,c`),
//! the few fixed formats EPSR uses (`'(a)'`, `'(a3)'`, `'(1x,a3,1x,a3)'`),
//! and sequential text files with `backspace`.
//!
//! List-directed semantics reproduced (gfortran behaviour):
//! * values are separated by blanks, commas or record ends, so a single READ
//!   may consume several lines;
//! * when a READ statement completes, the rest of the current line is skipped;
//! * a `/` terminates the statement leaving remaining items unchanged;
//! * `r*value` repeat counts and empty (null) values are honoured;
//! * reals accept `1.0`, `.5`, `5.`, `1e5`, `1.0D+01`, `1.0+01` forms and are
//!   converted with correct rounding directly to `f32` (gfortran uses `strtof`);
//! * reading an integer from a token such as `1.0` is an error, as in gfortran.

use std::fs;
use std::path::Path;

/// Error from a Fortran READ, mapped onto `iostat` values.
#[derive(Clone, Debug, PartialEq)]
pub enum IoErr {
    /// End of file / end of internal record (`iostat < 0`).
    Eof,
    /// Conversion error (`iostat > 0`).
    Bad(String),
}

impl IoErr {
    /// The `iostat` value gfortran would return (-1 for EOF, >0 otherwise).
    pub fn iostat(&self) -> i32 {
        match self {
            IoErr::Eof => -1,
            IoErr::Bad(_) => 5010,
        }
    }
}

pub type IoRes<T> = Result<T, IoErr>;

/// Convert an `IoRes` into the integer `iostat` code.
pub fn iostat<T>(r: &IoRes<T>) -> i32 {
    match r {
        Ok(_) => 0,
        Err(e) => e.iostat(),
    }
}

/// Parse a Fortran real constant token.
pub fn parse_real(tok: &str) -> IoRes<f32> {
    let t = normalise_real(tok)?;
    t.parse::<f32>().map_err(|_| IoErr::Bad(format!("bad real {tok:?}")))
}

/// Parse a Fortran double precision token.
pub fn parse_dble(tok: &str) -> IoRes<f64> {
    let t = normalise_real(tok)?;
    t.parse::<f64>().map_err(|_| IoErr::Bad(format!("bad real {tok:?}")))
}

fn normalise_real(tok: &str) -> IoRes<String> {
    let s = tok.trim();
    if s.is_empty() {
        return Err(IoErr::Bad("empty real".into()));
    }
    let low = s.to_ascii_lowercase();
    if matches!(low.trim_start_matches(['+', '-']), "inf" | "infinity" | "nan") {
        return Ok(low);
    }
    let b = s.as_bytes();
    let mut out = String::with_capacity(s.len() + 2);
    let mut i = 0;
    if b[0] == b'+' || b[0] == b'-' {
        out.push(b[0] as char);
        i = 1;
    }
    let mut nd = 0;
    let mut seen_dot = false;
    while i < b.len() && (b[i].is_ascii_digit() || (b[i] == b'.' && !seen_dot)) {
        if b[i] == b'.' {
            seen_dot = true;
        } else {
            nd += 1;
        }
        out.push(b[i] as char);
        i += 1;
    }
    if nd == 0 {
        return Err(IoErr::Bad(format!("bad real {tok:?}")));
    }
    if i < b.len() {
        let c = b[i].to_ascii_lowercase();
        let mut j = i;
        if c == b'e' || c == b'd' || c == b'q' {
            j += 1;
        } else if c != b'+' && c != b'-' {
            return Err(IoErr::Bad(format!("bad real {tok:?}")));
        }
        out.push('e');
        if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
            out.push(b[j] as char);
            j += 1;
        }
        let st = j;
        while j < b.len() && b[j].is_ascii_digit() {
            out.push(b[j] as char);
            j += 1;
        }
        if j != b.len() || j == st {
            return Err(IoErr::Bad(format!("bad real {tok:?}")));
        }
    }
    Ok(out)
}

/// Parse a Fortran integer token.
pub fn parse_int(tok: &str) -> IoRes<i32> {
    let s = tok.trim();
    let s2 = s.strip_prefix('+').unwrap_or(s);
    if s2.is_empty() || !s2.trim_start_matches('-').bytes().all(|c| c.is_ascii_digit()) || s2 == "-" {
        return Err(IoErr::Bad(format!("bad integer {tok:?}")));
    }
    s2.parse::<i32>().map_err(|_| IoErr::Bad(format!("bad integer {tok:?}")))
}

/// One list-directed value as seen by the reader.
#[derive(Clone, Debug)]
enum Tok {
    Val(String),
    Null,
    Slash,
}

/// Split one record into list-directed tokens.
fn tokenize(line: &str) -> Vec<Tok> {
    let b = line.as_bytes();
    let mut v = Vec::new();
    let mut i = 0;
    let mut last_was_comma = true; // a leading comma means a null value
    let mut pending_value = false;
    while i < b.len() {
        let c = b[i];
        if c == b' ' || c == b'\t' || c == b'\r' {
            i += 1;
            continue;
        }
        if c == b',' {
            if last_was_comma && !pending_value {
                v.push(Tok::Null);
            }
            last_was_comma = true;
            pending_value = false;
            i += 1;
            continue;
        }
        if c == b'/' {
            v.push(Tok::Slash);
            return v;
        }
        // quoted string
        if c == b'\'' || c == b'"' {
            let q = c;
            let mut s = String::new();
            i += 1;
            while i < b.len() {
                if b[i] == q {
                    if i + 1 < b.len() && b[i + 1] == q {
                        s.push(q as char);
                        i += 2;
                        continue;
                    }
                    i += 1;
                    break;
                }
                s.push(b[i] as char);
                i += 1;
            }
            v.push(Tok::Val(s));
            last_was_comma = false;
            pending_value = true;
            continue;
        }
        let st = i;
        while i < b.len() && !matches!(b[i], b' ' | b'\t' | b',' | b'/' | b'\r') {
            i += 1;
        }
        let tok = &line[st..i];
        // repeat count r*value or r*
        if let Some(p) = tok.find('*') {
            if let Ok(r) = tok[..p].parse::<usize>() {
                let val = &tok[p + 1..];
                for _ in 0..r {
                    if val.is_empty() {
                        v.push(Tok::Null);
                    } else {
                        v.push(Tok::Val(val.to_string()));
                    }
                }
                last_was_comma = false;
                pending_value = true;
                continue;
            }
        }
        v.push(Tok::Val(tok.to_string()));
        last_was_comma = false;
        pending_value = true;
    }
    v
}

/// A sequential formatted text file opened for reading, with a record pointer.
///
/// The lines are normally read completely when the file is opened.  A file
/// can also be backed by a stream ([`TextFile::from_reader`], e.g. standard
/// input for the interactive EPSR utilities); lines are then read on demand.
#[derive(Default)]
pub struct TextFile {
    pub lines: Vec<String>,
    pub pos: usize,
    reader: Option<Box<dyn std::io::BufRead + Send>>,
}

impl std::fmt::Debug for TextFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TextFile").field("lines", &self.lines.len()).field("pos", &self.pos).field("streaming", &self.reader.is_some()).finish()
    }
}

impl TextFile {
    /// `open(..., status='old')` + read everything.
    pub fn open(path: &Path) -> std::io::Result<TextFile> {
        let raw = fs::read(path)?;
        let text = String::from_utf8_lossy(&raw);
        let mut lines: Vec<String> = text.split('\n').map(|s| s.trim_end_matches('\r').to_string()).collect();
        if text.ends_with('\n') {
            lines.pop();
        }
        Ok(TextFile { lines, pos: 0, reader: None })
    }
    pub fn from_string(s: &str) -> TextFile {
        let mut lines: Vec<String> = s.split('\n').map(|s| s.trim_end_matches('\r').to_string()).collect();
        if s.ends_with('\n') {
            lines.pop();
        }
        TextFile { lines, pos: 0, reader: None }
    }
    /// A file whose records are read on demand from `r` (e.g. standard input).
    pub fn from_reader(r: Box<dyn std::io::BufRead + Send>) -> TextFile {
        TextFile { lines: Vec::new(), pos: 0, reader: Some(r) }
    }
    /// `read(u,'(a)') line`.
    pub fn read_line(&mut self) -> IoRes<String> {
        if self.pos >= self.lines.len() {
            if let Some(r) = self.reader.as_mut() {
                let mut l = String::new();
                match r.read_line(&mut l) {
                    Ok(0) | Err(_) => return Err(IoErr::Eof),
                    Ok(_) => {
                        let l = l.trim_end_matches('\n').trim_end_matches('\r').to_string();
                        self.lines.push(l);
                    }
                }
            } else {
                return Err(IoErr::Eof);
            }
        }
        self.pos += 1;
        Ok(self.lines[self.pos - 1].clone())
    }
    /// `read(u,*)` with no items: skip one record.
    pub fn skip_record(&mut self) -> IoRes<()> {
        self.read_line().map(|_| ())
    }
    /// `backspace(u)`.
    pub fn backspace(&mut self) {
        if self.pos > 0 {
            self.pos -= 1;
        }
    }
    /// Begin a list-directed READ statement.
    pub fn list(&mut self) -> ListRead<'_> {
        ListRead { src: Src::File(self), toks: Vec::new(), k: 0, started: false, terminated: false }
    }
}

enum Src<'a> {
    File(&'a mut TextFile),
    Str,
}

/// A list-directed READ statement in progress.  Obtain values in order with
/// [`ListRead::int`], [`ListRead::real`], ...; each returns `Ok(None)` if the
/// value was null or the statement was terminated by `/` (the Fortran variable
/// keeps its old value).  Dropping the reader ends the statement.
pub struct ListRead<'a> {
    src: Src<'a>,
    toks: Vec<Tok>,
    k: usize,
    started: bool,
    terminated: bool,
}

impl<'a> ListRead<'a> {
    /// An internal list-directed read from a character variable: `read(s,*) ...`.
    pub fn from_str(s: &str) -> ListRead<'static> {
        ListRead { src: Src::Str, toks: tokenize(s), k: 0, started: true, terminated: false }
    }

    fn next_tok(&mut self) -> IoRes<Tok> {
        if self.terminated {
            return Ok(Tok::Null);
        }
        loop {
            if self.k < self.toks.len() {
                let t = self.toks[self.k].clone();
                self.k += 1;
                if let Tok::Slash = t {
                    self.terminated = true;
                    return Ok(Tok::Null);
                }
                return Ok(t);
            }
            match &mut self.src {
                Src::Str => return Err(IoErr::Eof),
                Src::File(f) => {
                    let line = f.read_line()?;
                    self.started = true;
                    self.toks = tokenize(&line);
                    self.k = 0;
                }
            }
        }
    }
    fn raw(&mut self) -> IoRes<Option<String>> {
        match self.next_tok()? {
            Tok::Val(s) => Ok(Some(s)),
            _ => Ok(None),
        }
    }
    pub fn int(&mut self) -> IoRes<Option<i32>> {
        match self.raw()? {
            Some(s) => parse_int(&s).map(Some),
            None => Ok(None),
        }
    }
    pub fn real(&mut self) -> IoRes<Option<f32>> {
        match self.raw()? {
            Some(s) => parse_real(&s).map(Some),
            None => Ok(None),
        }
    }
    pub fn dble(&mut self) -> IoRes<Option<f64>> {
        match self.raw()? {
            Some(s) => parse_dble(&s).map(Some),
            None => Ok(None),
        }
    }
    /// A list-directed character item (unquoted word or quoted string).
    pub fn word(&mut self) -> IoRes<Option<String>> {
        self.raw()
    }
    /// Read into an existing variable, leaving it unchanged on a null value.
    pub fn get_i(&mut self, v: &mut i32) -> IoRes<()> {
        if let Some(x) = self.int()? {
            *v = x;
        }
        Ok(())
    }
    pub fn get_r(&mut self, v: &mut f32) -> IoRes<()> {
        if let Some(x) = self.real()? {
            *v = x;
        }
        Ok(())
    }
    pub fn get_d(&mut self, v: &mut f64) -> IoRes<()> {
        if let Some(x) = self.dble()? {
            *v = x;
        }
        Ok(())
    }
    /// Finish a READ with no items consumed: still consumes one record.
    pub fn finish(mut self) -> IoRes<()> {
        if !self.started {
            if let Src::File(f) = &mut self.src {
                f.read_line()?;
            }
        }
        Ok(())
    }
}

/// Internal list-directed read of a single real: `read(s,*) x`.
pub fn read_real(s: &str) -> IoRes<f32> {
    ListRead::from_str(s).real()?.ok_or(IoErr::Eof)
}
/// Internal list-directed read of a single integer: `read(s,*) i`.
pub fn read_int(s: &str) -> IoRes<i32> {
    ListRead::from_str(s).int()?.ok_or(IoErr::Eof)
}
/// Internal list-directed read of a character item.
pub fn read_word(s: &str) -> IoRes<String> {
    ListRead::from_str(s).word()?.ok_or(IoErr::Eof)
}

/// Formatted read `'(1x,a3,1x,a3)'` style: extract columns [from, from+len)
/// (0-based start) from a record, blank-padded, trailing blanks removed.
pub fn columns(line: &str, from: usize, len: usize) -> String {
    let b = line.as_bytes();
    let mut s = String::new();
    for k in from..from + len {
        s.push(if k < b.len() { b[k] as char } else { ' ' });
    }
    s.trim_end_matches(' ').to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reals() {
        assert_eq!(parse_real("1.0").unwrap(), 1.0);
        assert_eq!(parse_real("-.244E+01").unwrap(), -2.44);
        assert_eq!(parse_real("5.").unwrap(), 5.0);
        assert_eq!(parse_real("2").unwrap(), 2.0);
        assert_eq!(parse_real("1.0D+01").unwrap(), 10.0);
        assert_eq!(parse_real("1.0+01").unwrap(), 10.0);
        assert!(parse_real("abc").is_err());
        assert!(parse_int("1.0").is_err());
        assert_eq!(parse_int("-12").unwrap(), -12);
    }
    #[test]
    fn multi_line_record() {
        let mut f = TextFile::from_string("3 1 0.5 2\n 0.25 3 0.125 extra\nnext\n");
        let mut r = f.list();
        let n = r.int().unwrap().unwrap();
        let mut v = Vec::new();
        for _ in 0..n {
            v.push((r.int().unwrap().unwrap(), r.real().unwrap().unwrap()));
        }
        drop(r);
        assert_eq!(v, vec![(1, 0.5), (2, 0.25), (3, 0.125)]);
        assert_eq!(f.read_line().unwrap(), "next");
    }
    #[test]
    fn internal_reads() {
        assert_eq!(read_int("50 1").unwrap(), 50);
        assert_eq!(read_real("   8.80699679E-02").unwrap(), 8.80699679E-02);
        assert_eq!(read_word(" Gaussian  x").unwrap(), "Gaussian");
        let mut r = ListRead::from_str("1 2");
        assert_eq!(r.int().unwrap(), Some(1));
        assert_eq!(r.int().unwrap(), Some(2));
        assert_eq!(r.int(), Err(IoErr::Eof));
    }
}
