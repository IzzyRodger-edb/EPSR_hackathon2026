//! Emulation of Fortran formatted and list-directed *output*, byte-compatible
//! with gfortran.
//!
//! EPSR writes all of its results as ASCII text produced by Fortran `WRITE`
//! statements.  Down-stream tools (and the user's existing workflows) depend on
//! those exact layouts, and the validation of this port compares output files
//! byte-for-byte with the original program, so the Rust code reproduces
//! gfortran's formatting rules exactly:
//!
//! * [`fmt`] interprets a Fortran FORMAT string (e.g. `"(1x,i4,(5(1x,i4,1x,e9.3)))"`)
//!   against a list of values, including repeat counts, nested groups, `/`,
//!   `nX`, `nH` Hollerith constants, quoted strings and *format reversion*.
//!   Supported data edit descriptors: `Iw`, `Iw.m`, `I0`, `Fw.d`, `Ew.d`,
//!   `Ew.dEe`, `Aw`, `A`, `Lw`.
//! * [`ListWriter`] reproduces gfortran list-directed output (`write(u,*)`).
//!
//! All numeric conversion is done from the exact binary value with
//! round-half-even, which is what glibc `printf` (used internally by
//! libgfortran) does.

/// A value passed to a formatted write.
#[derive(Clone, Debug)]
pub enum FV<'a> {
    /// Default integer (`integer*4`) or `integer*8`.
    I(i64),
    /// `real*4`.
    R(f32),
    /// `real*8`.
    D(f64),
    /// A character variable: text plus its declared length (trailing blank
    /// padding up to that length is significant in Fortran output).
    Ch(&'a str, usize),
    /// A character expression whose length is exactly its text (e.g. a substring).
    S(&'a str),
    /// Logical.
    L(bool),
}

impl<'a> FV<'a> {
    fn as_f64(&self) -> Option<(f64, bool)> {
        match *self {
            FV::R(v) => Some((v as f64, false)),
            FV::D(v) => Some((v, true)),
            FV::I(v) => Some((v as f64, false)),
            _ => None,
        }
    }
}

/// Convenience constructors.
pub fn fi(v: i32) -> FV<'static> {
    FV::I(v as i64)
}
pub fn fr(v: f32) -> FV<'static> {
    FV::R(v)
}
pub fn fd(v: f64) -> FV<'static> {
    FV::D(v)
}

// ---------------------------------------------------------------------------
// Decimal conversion helpers
// ---------------------------------------------------------------------------

/// Return (digits, exp10) such that |v| ≈ 0.d1d2..dn × 10^exp10 with exactly `n`
/// significant digits (correctly rounded, ties-to-even).  For v == 0 returns n
/// zeros and exponent 0.
fn sig_digits(v: f64, n: usize) -> (String, i32) {
    let a = v.abs();
    if a == 0.0 {
        return ("0".repeat(n.max(1)), 0);
    }
    let s = format!("{:.*e}", n.saturating_sub(1), a);
    let (mant, exp) = s.split_once('e').unwrap();
    let e: i32 = exp.parse().unwrap();
    let digits: String = mant.chars().filter(|c| c.is_ascii_digit()).collect();
    (digits, e + 1)
}

fn stars(w: usize) -> String {
    "*".repeat(w)
}

/// gfortran's output of a non-finite value with `Fw.d`/`Ew.d` editing
/// (`write_infnan`): `NaN` (never signed) or `Inf`/`Infinity` with a `-`
/// sign when negative, right-justified; asterisks if it does not fit.
fn edit_infnan(v: f64, w: usize) -> String {
    if v.is_nan() {
        if w < 3 {
            return stars(w);
        }
        return format!("{:>w$}", "NaN", w = w);
    }
    let neg = v.is_sign_negative();
    if w < 3 || (neg && w < 4) {
        return stars(w);
    }
    let body = if (neg && w > 8) || (!neg && w > 7) { "Infinity" } else { "Inf" };
    let s = if neg { format!("-{body}") } else { body.to_string() };
    format!("{:>w$}", s, w = w)
}

/// Fortran `Ew.d` (or `Ew.dEe` when `e` is Some).
pub fn edit_e(v: f64, w: usize, d: usize, e: Option<usize>) -> String {
    if !v.is_finite() {
        return edit_infnan(v, w);
    }
    let neg = v.is_sign_negative();
    let (digits, mut ex) = sig_digits(v, d);
    if v == 0.0 {
        ex = 0;
    }
    let mut s = String::new();
    if neg {
        s.push('-');
    }
    s.push_str("0.");
    s.push_str(&digits[..d.min(digits.len())]);
    // exponent part
    let aex = ex.unsigned_abs();
    let sign = if ex < 0 { '-' } else { '+' };
    match e {
        None => {
            if aex <= 99 {
                s.push_str(&format!("E{}{:02}", sign, aex));
            } else if aex <= 999 {
                s.push_str(&format!("{}{:03}", sign, aex));
            } else {
                return stars(w);
            }
        }
        Some(ew) => {
            let es = format!("{:0width$}", aex, width = ew);
            if es.len() > ew {
                return stars(w);
            }
            s.push_str(&format!("E{}{}", sign, es));
        }
    }
    if s.len() > w {
        // drop the optional leading zero
        let alt = if neg { format!("-{}", &s[2..]) } else { s[1..].to_string() };
        if alt.len() > w {
            return stars(w);
        }
        s = alt;
    }
    format!("{:>w$}", s, w = w)
}

/// Fortran `Fw.d`.
pub fn edit_f(v: f64, w: usize, d: usize) -> String {
    if !v.is_finite() {
        return edit_infnan(v, w);
    }
    let neg = v.is_sign_negative();
    let mut body = format!("{:.*}", d, v.abs());
    if d == 0 {
        body.push('.');
    }
    let mut s = if neg { format!("-{}", body) } else { body.clone() };
    if s.len() > w || w == 0 {
        if body.starts_with("0.") {
            let b2 = &body[1..];
            s = if neg { format!("-{}", b2) } else { b2.to_string() };
        }
        if s.len() > w {
            return stars(w);
        }
    }
    format!("{:>w$}", s, w = w)
}

/// Fortran `Iw`, `Iw.m` and `I0`.
pub fn edit_i(v: i64, w: usize, m: Option<usize>) -> String {
    let mut digits = v.unsigned_abs().to_string();
    if let Some(m) = m {
        if m == 0 && v == 0 {
            digits = String::new();
        } else if digits.len() < m {
            digits = format!("{}{}", "0".repeat(m - digits.len()), digits);
        }
    }
    let s = if v < 0 { format!("-{}", digits) } else { digits };
    if w == 0 {
        return s;
    }
    if s.len() > w {
        return stars(w);
    }
    format!("{:>w$}", s, w = w)
}

// ---------------------------------------------------------------------------
// Format string parsing
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
enum Item {
    Group(u32, Vec<Item>),
    I(u32, usize, Option<usize>),
    F(u32, usize, usize),
    E(u32, usize, usize, Option<usize>),
    A(u32, Option<usize>),
    L(u32, usize),
    X(usize),
    Lit(String),
    Slash(u32),
    Colon,
    /// Scale factor `kP` (accepted, only `0P` semantics are implemented)
    #[allow(dead_code)]
    P(i32),
}

struct Parser<'a> {
    s: &'a [u8],
    p: usize,
}

impl<'a> Parser<'a> {
    fn peek(&self) -> Option<u8> {
        self.s.get(self.p).copied()
    }
    fn skip_ws(&mut self) {
        while let Some(c) = self.peek() {
            if c == b' ' || c == b',' {
                self.p += 1;
            } else {
                break;
            }
        }
    }
    fn num(&mut self) -> Option<usize> {
        let st = self.p;
        while let Some(c) = self.peek() {
            if c.is_ascii_digit() {
                self.p += 1;
            } else {
                break;
            }
        }
        if self.p > st {
            Some(std::str::from_utf8(&self.s[st..self.p]).unwrap().parse().unwrap())
        } else {
            None
        }
    }
    /// Parse items until the matching ')' (consumed).
    fn list(&mut self) -> Vec<Item> {
        let mut items = Vec::new();
        loop {
            self.skip_ws();
            let c = match self.peek() {
                None => break,
                Some(c) => c,
            };
            if c == b')' {
                self.p += 1;
                break;
            }
            if c == b'\'' || c == b'"' {
                items.push(Item::Lit(self.quoted(c)));
                continue;
            }
            if c == b'/' {
                self.p += 1;
                items.push(Item::Slash(1));
                continue;
            }
            if c == b':' {
                self.p += 1;
                items.push(Item::Colon);
                continue;
            }
            let neg = if c == b'-' {
                self.p += 1;
                true
            } else {
                false
            };
            let rep = self.num();
            self.skip_spaces_only();
            let c = self.peek().map(|c| c.to_ascii_uppercase());
            match c {
                Some(b'(') => {
                    self.p += 1;
                    let g = self.list();
                    items.push(Item::Group(rep.unwrap_or(1) as u32, g));
                }
                Some(b'H') => {
                    self.p += 1;
                    let n = rep.expect("Hollerith needs count");
                    let lit = std::str::from_utf8(&self.s[self.p..self.p + n]).unwrap().to_string();
                    self.p += n;
                    items.push(Item::Lit(lit));
                }
                Some(b'X') => {
                    self.p += 1;
                    items.push(Item::X(rep.unwrap_or(1)));
                }
                Some(b'P') => {
                    self.p += 1;
                    let k = rep.unwrap_or(0) as i32;
                    items.push(Item::P(if neg { -k } else { k }));
                }
                Some(b'/') => {
                    self.p += 1;
                    items.push(Item::Slash(rep.unwrap_or(1) as u32));
                }
                Some(b'I') => {
                    self.p += 1;
                    let w = self.num().unwrap_or(0);
                    let m = if self.peek() == Some(b'.') {
                        self.p += 1;
                        self.num()
                    } else {
                        None
                    };
                    items.push(Item::I(rep.unwrap_or(1) as u32, w, m));
                }
                Some(b'F') => {
                    self.p += 1;
                    let w = self.num().expect("Fw.d width");
                    assert_eq!(self.peek(), Some(b'.'));
                    self.p += 1;
                    let d = self.num().expect("Fw.d digits");
                    items.push(Item::F(rep.unwrap_or(1) as u32, w, d));
                }
                Some(b'E') => {
                    self.p += 1;
                    let w = self.num().expect("Ew.d width");
                    assert_eq!(self.peek(), Some(b'.'));
                    self.p += 1;
                    let d = self.num().expect("Ew.d digits");
                    let e = if matches!(self.peek(), Some(b'E') | Some(b'e')) {
                        self.p += 1;
                        self.num()
                    } else {
                        None
                    };
                    items.push(Item::E(rep.unwrap_or(1) as u32, w, d, e));
                }
                Some(b'A') => {
                    self.p += 1;
                    let w = self.num();
                    items.push(Item::A(rep.unwrap_or(1) as u32, w));
                }
                Some(b'L') => {
                    self.p += 1;
                    let w = self.num().unwrap_or(1);
                    items.push(Item::L(rep.unwrap_or(1) as u32, w));
                }
                Some(b'\'') | Some(b'"') => {
                    // e.g. 15x'text' with no comma is not used; treat repeat as error
                    let q = self.peek().unwrap();
                    items.push(Item::Lit(self.quoted(q)));
                }
                other => panic!(
                    "unsupported FORMAT descriptor {:?} in {:?}",
                    other.map(|c| c as char),
                    std::str::from_utf8(self.s).unwrap()
                ),
            }
        }
        items
    }
    fn skip_spaces_only(&mut self) {
        while self.peek() == Some(b' ') {
            self.p += 1;
        }
    }
    fn quoted(&mut self, q: u8) -> String {
        self.p += 1;
        let mut out = Vec::new();
        loop {
            let c = self.peek().expect("unterminated string in FORMAT");
            self.p += 1;
            if c == q {
                if self.peek() == Some(q) {
                    out.push(q);
                    self.p += 1;
                } else {
                    break;
                }
            } else {
                out.push(c);
            }
        }
        String::from_utf8(out).unwrap()
    }
}

fn parse_format(f: &str) -> Vec<Item> {
    let t = f.trim();
    let b = t.as_bytes();
    assert!(b.first() == Some(&b'('), "FORMAT must start with '(': {f:?}");
    let mut p = Parser { s: b, p: 1 };
    p.list()
}

// ---------------------------------------------------------------------------
// Format execution
// ---------------------------------------------------------------------------

struct Exec<'v, 'a> {
    vals: &'v [FV<'a>],
    next: usize,
    out: String,
    rec: String,
    pending: usize,
    done: bool,
}

impl<'v, 'a> Exec<'v, 'a> {
    fn emit(&mut self, s: &str) {
        if self.pending > 0 {
            for _ in 0..self.pending {
                self.rec.push(' ');
            }
            self.pending = 0;
        }
        self.rec.push_str(s);
    }
    fn newrec(&mut self) {
        self.out.push_str(&self.rec);
        self.out.push('\n');
        self.rec.clear();
        self.pending = 0;
    }
    fn take(&mut self) -> Option<FV<'a>> {
        if self.next < self.vals.len() {
            self.next += 1;
            Some(self.vals[self.next - 1].clone())
        } else {
            None
        }
    }
    /// Run a list of items once.  Returns false if output terminated (no more
    /// data for a data edit descriptor).
    fn run(&mut self, items: &[Item]) -> bool {
        for it in items {
            match it {
                Item::Group(n, g) => {
                    for _ in 0..*n {
                        if !self.run(g) {
                            return false;
                        }
                    }
                }
                Item::Lit(s) => {
                    let s = s.clone();
                    self.emit(&s);
                }
                Item::X(n) => self.pending += n,
                Item::Slash(n) => {
                    for _ in 0..*n {
                        self.newrec();
                    }
                }
                Item::Colon => {
                    if self.next >= self.vals.len() {
                        self.done = true;
                        return false;
                    }
                }
                Item::P(_) => {}
                Item::I(n, w, m) => {
                    for _ in 0..*n {
                        let v = match self.take() {
                            Some(v) => v,
                            None => {
                                self.done = true;
                                return false;
                            }
                        };
                        let iv = match v {
                            FV::I(i) => i,
                            other => panic!("I edit descriptor given {other:?}"),
                        };
                        let s = edit_i(iv, *w, *m);
                        self.emit(&s);
                    }
                }
                Item::F(n, w, d) => {
                    for _ in 0..*n {
                        let v = match self.take() {
                            Some(v) => v,
                            None => {
                                self.done = true;
                                return false;
                            }
                        };
                        let (x, _) = v.as_f64().expect("F edit needs a number");
                        let s = edit_f(x, *w, *d);
                        self.emit(&s);
                    }
                }
                Item::E(n, w, d, e) => {
                    for _ in 0..*n {
                        let v = match self.take() {
                            Some(v) => v,
                            None => {
                                self.done = true;
                                return false;
                            }
                        };
                        let (x, _) = v.as_f64().expect("E edit needs a number");
                        let s = edit_e(x, *w, *d, *e);
                        self.emit(&s);
                    }
                }
                Item::A(n, w) => {
                    for _ in 0..*n {
                        let v = match self.take() {
                            Some(v) => v,
                            None => {
                                self.done = true;
                                return false;
                            }
                        };
                        let (txt, len) = match v {
                            FV::Ch(s, l) => (s, l),
                            FV::S(s) => (s, s.len()),
                            other => panic!("A edit descriptor given {other:?}"),
                        };
                        let full = pad_to(txt, len);
                        let s = match w {
                            None => full,
                            Some(w) => {
                                if *w >= len {
                                    format!("{}{}", " ".repeat(*w - len), full)
                                } else {
                                    full[..*w].to_string()
                                }
                            }
                        };
                        self.emit(&s);
                    }
                }
                Item::L(n, w) => {
                    for _ in 0..*n {
                        let v = match self.take() {
                            Some(v) => v,
                            None => {
                                self.done = true;
                                return false;
                            }
                        };
                        let b = match v {
                            FV::L(b) => b,
                            other => panic!("L edit descriptor given {other:?}"),
                        };
                        let s = format!("{:>w$}", if b { "T" } else { "F" }, w = *w);
                        self.emit(&s);
                    }
                }
            }
        }
        true
    }
}

fn has_data_descriptor(items: &[Item]) -> bool {
    items.iter().any(|it| match it {
        Item::Group(_, g) => has_data_descriptor(g),
        Item::I(..) | Item::F(..) | Item::E(..) | Item::A(..) | Item::L(..) => true,
        _ => false,
    })
}

/// Pad (or truncate) `s` to exactly `len` characters, Fortran style.
pub fn pad_to(s: &str, len: usize) -> String {
    let n = s.chars().count();
    if n >= len {
        s.chars().take(len).collect()
    } else {
        let mut o = String::with_capacity(len);
        o.push_str(s);
        for _ in n..len {
            o.push(' ');
        }
        o
    }
}

/// Apply a Fortran FORMAT to `vals`, returning the produced records joined by
/// `'\n'` (no trailing newline).  Format reversion produces additional records.
pub fn fmt(format: &str, vals: &[FV]) -> String {
    let items = parse_format(format);
    let mut ex = Exec { vals, next: 0, out: String::new(), rec: String::new(), pending: 0, done: false };
    // index of the last top-level group, for reversion
    let rev_start = items.iter().rposition(|it| matches!(it, Item::Group(..))).unwrap_or(0);
    let mut first = true;
    loop {
        let slice = if first { &items[..] } else { &items[rev_start..] };
        let cont = ex.run(slice);
        first = false;
        if !cont || ex.done {
            break;
        }
        if ex.next >= vals.len() {
            break;
        }
        if !has_data_descriptor(&items) {
            break;
        }
        // reversion: new record
        ex.newrec();
    }
    let mut out = ex.out;
    out.push_str(&ex.rec);
    out
}

// ---------------------------------------------------------------------------
// List-directed output
// ---------------------------------------------------------------------------

/// gfortran list-directed `REAL(4)` field (16 characters, 9 significant digits).
pub fn list_real4(v: f32) -> String {
    list_real(v as f64, 9, 2)
}

/// gfortran list-directed `REAL(8)` field (25 characters, 17 significant digits).
pub fn list_real8(v: f64) -> String {
    list_real(v, 17, 3)
}

fn list_real(v: f64, prec: usize, expw: usize) -> String {
    // field width: 1 (sign) + 1 + 1 + (prec-1) + 1(E) + 1(sign) + expw  (+1 leading blank)
    let width = prec + expw + 5; // 16 for real4, 25 for real8
    let neg = v.is_sign_negative();
    let a = v.abs();
    let sgn = if neg { "-" } else { "" };
    if a == 0.0 {
        let s = format!("{}0.{}", sgn, "0".repeat(prec - 1));
        let fw = width - (expw + 2);
        return format!("{:>fw$}{}", s, " ".repeat(expw + 2), fw = fw);
    }
    if !a.is_finite() {
        let s = if a.is_nan() { "NaN".to_string() } else { format!("{}Infinity", sgn) };
        return format!("{:>w$}", s, w = width);
    }
    let s = format!("{:.*e}", prec - 1, a);
    let (mant, exp) = s.split_once('e').unwrap();
    let e: i32 = exp.parse().unwrap();
    let digits: String = mant.chars().filter(|c| c.is_ascii_digit()).collect();
    if e >= -1 && e < prec as i32 {
        // F form with `prec` significant digits
        let body = if e >= 0 {
            let ip = (e + 1) as usize;
            format!("{}.{}", &digits[..ip], &digits[ip..])
        } else {
            format!("0.{}", digits)
        };
        let s = format!("{}{}", sgn, body);
        let fw = width - (expw + 2);
        format!("{:>fw$}{}", s, " ".repeat(expw + 2), fw = fw)
    } else {
        let es = format!("{}{}.{}E{}{:0ew$}", sgn, &digits[..1], &digits[1..], if e < 0 { '-' } else { '+' }, e.unsigned_abs(), ew = expw);
        format!("{:>w$}", es, w = width)
    }
}

/// Builder for one record of gfortran list-directed output (`write(u,*) ...`).
///
/// ```
/// use epsr_core::fortran::format::ListWriter;
/// let mut w = ListWriter::new();
/// w.int(100).real(-40.6027222);
/// assert_eq!(w.finish(), "         100  -40.6027222    ");
/// ```
#[derive(Default)]
pub struct ListWriter {
    s: String,
    first: bool,
    last_char: bool,
}

impl ListWriter {
    pub fn new() -> Self {
        ListWriter { s: String::new(), first: true, last_char: false }
    }
    fn sep(&mut self, is_char: bool) {
        if self.first {
            self.s.push(' ');
            self.first = false;
        } else if !(is_char && self.last_char) {
            self.s.push(' ');
        }
        self.last_char = is_char;
    }
    /// Default integer (`integer*4`): width 11.
    pub fn int(&mut self, v: i32) -> &mut Self {
        self.sep(false);
        self.s.push_str(&format!("{:>11}", v));
        self
    }
    /// `integer*8`: width 20.
    pub fn int8(&mut self, v: i64) -> &mut Self {
        self.sep(false);
        self.s.push_str(&format!("{:>20}", v));
        self
    }
    pub fn real(&mut self, v: f32) -> &mut Self {
        self.sep(false);
        self.s.push_str(&list_real4(v));
        self
    }
    pub fn dble(&mut self, v: f64) -> &mut Self {
        self.sep(false);
        self.s.push_str(&list_real8(v));
        self
    }
    pub fn logical(&mut self, v: bool) -> &mut Self {
        self.sep(false);
        self.s.push(if v { 'T' } else { 'F' });
        self
    }
    /// Character item (written verbatim, including any trailing blanks given).
    pub fn chars(&mut self, v: &str) -> &mut Self {
        self.sep(true);
        self.s.push_str(v);
        self
    }
    pub fn finish(&mut self) -> String {
        if self.first {
            // empty write(*,*) produces an empty record
            return String::new();
        }
        std::mem::take(&mut self.s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn non_finite_match_gfortran() {
        let inf = f64::INFINITY;
        let got: Vec<String> = [3, 4, 7, 8, 9, 10].iter().map(|w| edit_f(inf, *w, 1)).collect();
        assert_eq!(got, ["Inf", " Inf", "    Inf", "Infinity", " Infinity", "  Infinity"]);
        let got: Vec<String> = [3, 4, 7, 8, 9, 10].iter().map(|w| edit_f(-inf, *w, 1)).collect();
        assert_eq!(got, ["***", "-Inf", "   -Inf", "    -Inf", "-Infinity", " -Infinity"]);
        assert_eq!(edit_f(f64::NAN, 2, 1), "**");
        assert_eq!(edit_f(-f64::NAN, 7, 2), "    NaN");
        assert_eq!(edit_e(f64::NAN, 12, 5, None), "         NaN");
        assert_eq!(edit_e(-inf, 12, 5, None), "   -Infinity");
    }
    #[test]
    fn list_reals_match_gfortran() {
        assert_eq!(list_real4(1.0), "  1.00000000    ");
        assert_eq!(list_real4(0.0), "  0.00000000    ");
        assert_eq!(list_real4(-1.0), " -1.00000000    ");
        assert_eq!(list_real4(0.1), " 0.100000001    ");
        assert_eq!(list_real4(0.0999999), "  9.99998972E-02");
        assert_eq!(list_real4(1.0e8), "  100000000.    ");
        assert_eq!(list_real4(1.0e9), "  1.00000000E+09");
        assert_eq!(list_real4(-0.000123), " -1.23000005E-04");
        assert_eq!(list_real4(-40.6027222), " -40.6027222    ");
        assert_eq!(list_real8(1.0), "  1.0000000000000000     ");
        assert_eq!(list_real8(0.1), " 0.10000000000000001     ");
        assert_eq!(list_real8(-3.3333333333333e-5), " -3.3333333333333003E-005");
        assert_eq!(list_real8(1.0e20), "  1.0000000000000000E+020");
    }
    #[test]
    fn list_record() {
        let mut w = ListWriter::new();
        w.int(1).real(2.0).int(-3).dble(4.0).int(5);
        assert_eq!(
            w.finish(),
            "           1   2.00000000              -3   4.0000000000000000                5"
        );
        let mut w = ListWriter::new();
        w.chars("abc").real(1.0).chars("de");
        assert_eq!(w.finish(), " abc   1.00000000     de");
    }
    #[test]
    fn edit_descriptors() {
        assert_eq!(
            fmt("(e14.7,1x,e13.6,1x,e12.5,1x,e9.3,1x,e9.3)", &[fr(1.0), fr(-0.123456789), fr(0.0), fr(1.82), fr(-2.44)]),
            " 0.1000000E+01 -0.123457E+00  0.00000E+00 0.182E+01 -.244E+01"
        );
        assert_eq!(
            fmt("(f10.5,1x,f8.3,1x,f6.3,1x,f14.8,1x,f10.2,f7.3)", &[fr(1.0), fr(-0.0005), fr(12.3456), fr(-0.1234567), fr(99999999.0), fr(-0.0)]),
            "   1.00000   -0.001 12.346    -0.12345670 ********** -0.000"
        );
        assert_eq!(fmt("(i0,1x,i5,1x,i2.2,1x,i3.3,1x,i4)", &[fi(42), fi(-7), fi(3), fi(5), fi(12345)]), "42    -7 03 005 ****");
        assert_eq!(fmt("(e13.6)", &[fr(0.9999999e10)]), " 0.100000E+11");
        let v: Vec<FV> = vec![
            fi(7), fi(1), fr(1.5), fi(2), fr(-2.5), fi(3), fr(3.5), fi(4), fr(4.5), fi(5), fr(5.5), fi(6), fr(6.5), fi(7), fr(7.5),
        ];
        assert_eq!(
            fmt("(1x,i4,(5(1x,i4,1x,e9.3)))", &v),
            "    7    1 0.150E+01    2 -.250E+01    3 0.350E+01    4 0.450E+01    5 0.550E+01\n    6 0.650E+01    7 0.750E+01"
        );
        assert_eq!(fmt("('(3x,a,1h-,a,',i2.2,'x,a,1h-,a,',i2.2,'x)')", &[fi(9), fi(5)]), "(3x,a,1h-,a,09x,a,1h-,a,05x)");
        assert_eq!(fmt("(3x,a,1h-,a,09x,a,1h-,a,05x)", &[FV::S("S"), FV::S("C"), FV::S("C"), FV::S("S")]), "   S-C         C-S");
        assert_eq!(fmt("(a11,1x,a)", &[FV::Ch("abc", 11), FV::S("x")]), "abc         x");
    }
}
