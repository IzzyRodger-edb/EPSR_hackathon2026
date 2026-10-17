//! EPSR's own text utilities (`parse.f`, `concatenate`, `indexlastdot`, ...)
//! and helpers for Fortran fixed-length `character` semantics.
//!
//! In the Fortran, every character variable is blank-padded to its declared
//! length and comparisons ignore trailing blanks.  The Rust port stores
//! character variables as `String`s *without* trailing blanks; [`fstr`]
//! applies the declared-length truncation that a Fortran assignment performs.
//! Where trailing blanks matter (formatted output with `A` editing) the
//! declared length is passed explicitly to the formatter (see
//! [`crate::fortran::format::FV::Ch`]).

/// Truncate to the declared Fortran length `len` and drop trailing blanks,
/// i.e. the value a `character*len` variable holds after assignment.
pub fn fstr(s: &str, len: usize) -> String {
    let t: String = s.chars().take(len).collect();
    t.trim_end_matches(' ').to_string()
}

/// Fortran `index(s, sub)`: 1-based position of the first occurrence, 0 if absent.
pub fn findex(s: &str, sub: &str) -> i32 {
    match s.find(sub) {
        Some(p) => p as i32 + 1,
        None => 0,
    }
}

/// `index(text, ' ') - 1` applied to a blank-padded variable whose stored value
/// has no trailing blanks: the length up to the first blank.
pub fn len_to_blank(s: &str) -> usize {
    match s.find(' ') {
        Some(p) => p,
        None => s.len(),
    }
}

/// Fortran substring `s(i:j)` (1-based, inclusive) of a blank-padded variable.
/// Positions beyond the stored text are blanks.
pub fn substr(s: &str, i: i32, j: i32) -> String {
    if j < i {
        return String::new();
    }
    let b = s.as_bytes();
    let mut o = String::with_capacity((j - i + 1) as usize);
    for k in i..=j {
        let c = if k >= 1 && (k as usize) <= b.len() { b[(k - 1) as usize] as char } else { ' ' };
        o.push(c);
    }
    o
}

/// `indexlastdot(text)` from `auxiliary_routines.f`: position of the last '.'
/// (1-based); if there is no dot, `len_trim(text)+1` is *not* returned — the
/// Fortran returns 0 in that case (the comment in the source is misleading).
pub fn indexlastdot(text: &str) -> i32 {
    let t = text.trim_end_matches(' ');
    match t.rfind('.') {
        Some(p) => p as i32 + 1,
        None => 0,
    }
}

/// Word boundaries produced by `parse` (1-based inclusive character positions,
/// exactly as `ncf(i)`, `ncl(i)` in the Fortran).
#[derive(Clone, Debug, Default)]
pub struct Words {
    pub ncf: Vec<i32>,
    pub ncl: Vec<i32>,
}

impl Words {
    pub fn nwords(&self) -> i32 {
        self.ncf.len() as i32
    }
    /// The text of word `iw` (1-based).
    pub fn word<'a>(&self, text: &'a str, iw: i32) -> &'a str {
        let a = self.ncf[(iw - 1) as usize] as usize;
        let b = self.ncl[(iw - 1) as usize] as usize;
        &text[a - 1..b.min(text.len())]
    }
    /// `text(ncf(i):ncl(j))` — span from the start of word i to the end of word j.
    pub fn span(&self, text: &str, i: i32, j: i32) -> String {
        substr(text, self.ncf[(i - 1) as usize], self.ncl[(j - 1) as usize])
    }
}

/// Maximum number of words kept (`nworddim` in commands.inc).
pub const NWORDDIM: usize = 201;
/// Length of the `sometext` buffer (`nstsize` in commands.inc).
pub const NSTSIZE: usize = 3020;

/// Translation of `subroutine parse(sometext,nwords,ncf,ncl,nspcword,nchartext,nworddim)`.
///
/// Splits `text` (conceptually blank-padded to `nchartext` characters) into
/// blank-separated words.  Scanning stops once more than `nspcword`
/// consecutive blanks have been seen, which is how EPSR ignores the
/// free-text description that follows each value in its input files.
pub fn parse_n(text: &str, nspcword: i32, nchartext: usize) -> Words {
    let b = text.as_bytes();
    let ch = |i: usize| -> u8 {
        if i >= 1 && i <= b.len() && i <= nchartext {
            b[i - 1]
        } else {
            b' '
        }
    };
    let mut w = Words::default();
    let n = nchartext;
    let mut nspace: i32 = 0;
    let mut i: usize = 1;
    while ch(i) == b' ' && i < n {
        i += 1;
    }
    let mut nfirst = i;
    let mut nlast = i;
    while i < n && nspace <= nspcword {
        i += 1;
        if ch(i) == b' ' {
            if nlast >= nfirst {
                if w.ncf.len() < NWORDDIM {
                    w.ncf.push(nfirst as i32);
                    w.ncl.push(nlast as i32);
                }
                nspace = 1;
            } else {
                nspace += 1;
            }
            nfirst = i + 1;
        } else {
            nlast = i;
        }
        // Optimisation (no behavioural change): once we are past the stored
        // text everything is blank, so the remaining iterations only count
        // blanks; finish them arithmetically.
        if i > b.len() && ch(i) == b' ' && nlast < nfirst {
            let remaining = n - i; // further iterations available
            let needed = (nspcword - nspace + 1).max(0) as usize; // until nspace > nspcword
            let steps = remaining.min(needed);
            i += steps;
            nspace += steps as i32;
            nfirst = i + 1;
            if nspace > nspcword || i >= n {
                break;
            }
        }
    }
    if ch(i) != b' ' && nlast >= nfirst && w.ncf.len() < NWORDDIM {
        w.ncf.push(nfirst as i32);
        w.ncl.push(nlast as i32);
    }
    w
}

/// `parse` with the usual `nchartext = nstsize`.
pub fn parse(text: &str, nspcword: i32) -> Words {
    parse_n(text, nspcword, NSTSIZE)
}

/// Translation of `concatenate(nfirst,nlast,return_string,ncharstring,nwords,ncf,ncl,input_string)`:
/// joins words `nfirst..=nlast` with single blanks (max 256 characters) and
/// truncates to `retlen` (the declared length of the destination).
pub fn concatenate(nfirst: i32, nlast: i32, w: &Words, input: &str, retlen: usize) -> (String, i32) {
    let ncharword: i32 = 256;
    let nwords = w.nwords();
    let nlast = nlast.min(nwords);
    let nfirst = nfirst.max(1);
    let mut out = String::new();
    let mut nchar: i32 = 0;
    let mut iw = nfirst;
    while iw <= nlast && nchar < ncharword {
        let cf = w.ncf[(iw - 1) as usize];
        let cl = w.ncl[(iw - 1) as usize];
        let mut ncharadd = cl - cf + 1;
        let mut nchartest = nchar + ncharadd;
        if iw > nfirst {
            nchartest += 1;
        }
        if nchartest > ncharword {
            ncharadd = ncharword - nchar;
            if iw > nfirst {
                ncharadd -= 1;
            }
        }
        let piece = substr(input, cf, cl.min(cf + ncharadd - 1));
        if iw == nfirst {
            out = piece;
            nchar = ncharadd;
        } else {
            if ncharadd > 0 {
                out = format!("{} {}", substr(&out, 1, nchar), piece);
            }
            nchar = nchar + ncharadd + 1;
        }
        iw += 1;
    }
    let nchar = nchar.min(retlen as i32).max(0);
    let res = substr(&out, 1, nchar);
    (res.trim_end_matches(' ').to_string(), nchar)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn parse_stops_after_long_gap() {
        let t = "feedback    0.9               Confidence factor - should be < 1. [0.8]";
        let w = parse(t, 12);
        assert_eq!(w.nwords(), 2);
        assert_eq!(w.word(t, 1), "feedback");
        assert_eq!(w.word(t, 2), "0.9");
        let w = parse(t, 20);
        assert!(w.nwords() > 2);
    }
    #[test]
    fn parse_leading_blanks() {
        let t = "   a  bb ccc";
        let w = parse(t, 5);
        assert_eq!(w.ncf, vec![4, 7, 10]);
        assert_eq!(w.ncl, vec![4, 8, 12]);
    }
    #[test]
    fn concatenate_words() {
        let t = "potfac      1.000000   1.000000";
        let w = parse(t, 12);
        let (s, n) = concatenate(2, w.nwords(), &w, t, 256);
        assert_eq!(s, "1.000000 1.000000");
        assert_eq!(n, 17);
    }
    #[test]
    fn lastdot() {
        assert_eq!(indexlastdot("a.b.c"), 4);
        assert_eq!(indexlastdot("abc"), 0);
    }
}
