//! Reading of the diffraction data files in `runepsr17`
//! (`epsr_standalone_rev3.f`, "read in data to be fitted"), for all data
//! types `nrtype`:
//!
//! * 1 — `Q, S(Q)` (list-directed, one set per line, errors set to 0);
//! * 2 — `Q, S(Q), E(Q)` (list-directed);
//! * 3, 4 — GENIE "SHOW DATA" format: `N, Q, S(Q), E(Q)` read with
//!   `format(i7,3g14.7)` after 4 header lines (3 = histogram, 4 = point mode);
//! * 5 — GUDRUN `.mint01`/`.dcs` format, `Q, S(Q), E(Q)` (list-directed),
//!   lines containing `#` are skipped (histogram).
//!
//! The record handling (first `read(8,*)`, header skips, `backspace(8)`)
//! is reproduced exactly: for nrtype 1 and 2 the first line of the file is
//! re-read as data, for nrtype 3 and 4 the 5th record is the first data line.

use crate::fortran::listread::{iostat, parse_int, parse_real, IoErr, IoRes, TextFile};
use crate::fortran::F1;
use crate::state::MQD;

/// One `g` (or `f`/`e`) input field of width `w` with `d` fractional digits
/// (`BLANK='NULL'`: blanks are ignored, an all-blank field is zero; without a
/// decimal point the value is scaled by `10**-d`).
fn fmt_real_field(field: &str, d: i32) -> IoRes<f32> {
    let s: String = field.chars().filter(|c| *c != ' ').collect();
    if s.is_empty() {
        return Ok(0.0);
    }
    if s.contains('.') {
        return parse_real(&s);
    }
    // implied decimal point: split mantissa / exponent
    let b = s.as_bytes();
    let mut k = 0;
    if !b.is_empty() && (b[0] == b'+' || b[0] == b'-') {
        k = 1;
    }
    while k < b.len() && b[k].is_ascii_digit() {
        k += 1;
    }
    let mant = &s[..k];
    let rest = &s[k..];
    let mut ex: i32 = 0;
    if !rest.is_empty() {
        let r = rest.trim_start_matches(['e', 'E', 'd', 'D', 'q', 'Q']);
        ex = parse_int(r)?;
    }
    if mant.trim_start_matches(['+', '-']).is_empty() {
        return Err(IoErr::Bad(format!("bad real {field:?}")));
    }
    parse_real(&format!("{}e{}", mant, ex - d))
}

/// One `i` input field.
fn fmt_int_field(field: &str) -> IoRes<i32> {
    let s: String = field.chars().filter(|c| *c != ' ').collect();
    if s.is_empty() {
        return Ok(0);
    }
    parse_int(&s)
}

/// `read(8,2224,iostat=ierr) n,qd,data,err` with `2224 format(i7,3g14.7)`.
fn read_genie_record(f: &mut TextFile) -> IoRes<(i32, f32, f32, f32)> {
    let line = f.read_line()?;
    // PAD='YES': the record is padded with blanks to the length required
    let mut b: Vec<u8> = line.into_bytes();
    if b.len() < 49 {
        b.resize(49, b' ');
    }
    let fld = |a: usize, w: usize| String::from_utf8_lossy(&b[a..a + w]).into_owned();
    let n = fmt_int_field(&fld(0, 7))?;
    let q = fmt_real_field(&fld(7, 14), 7)?;
    let s = fmt_real_field(&fld(21, 14), 7)?;
    let e = fmt_real_field(&fld(35, 14), 7)?;
    Ok((n, q, s, e))
}

/// The data-reading loop of `runepsr17` for one file (unit 8), from the first
/// `read(8,*,iostat=ierr)` up to and including "Ignore the last point which
/// may be rubbish".  Returns `nqd`.
pub fn read_data_points(f: &mut TextFile, nrtype: i32, qd: &mut F1<f32>, data: &mut F1<f32>, err: &mut F1<f32>) -> i32 {
    // Read a line to see if there are any problems
    let mut ierr = iostat(&f.skip_record());
    // for Genie style files skip the next 4 lines
    if nrtype > 2 && nrtype < 5 {
        for _ in 1..=4 {
            ierr = iostat(&f.skip_record());
        }
    }
    // for GUDRUN files, ignore any lines that have # in them
    if nrtype > 4 {
        let rd = |f: &mut TextFile| -> (i32, String) {
            let mut r = f.list();
            match r.word() {
                Ok(Some(w)) => (0, w),
                Ok(None) => (0, String::new()),
                Err(e) => (e.iostat(), String::new()),
            }
        };
        let (mut ie, mut hashsearch) = rd(f);
        while hashsearch.contains('#') && ie == 0 {
            let (a, b) = rd(f);
            ie = a;
            hashsearch = b;
        }
        ierr = ie;
    }
    // Backspace the data file 1 line to the first good line
    f.backspace();
    let mut nqd: i32 = 0;
    while (nqd as usize) < MQD && ierr == 0 {
        let nqd1 = nqd + 1;
        if nrtype == 1 {
            let rec = {
                let mut r = f.list();
                (|| -> IoRes<(Option<f32>, Option<f32>)> { Ok((r.real()?, r.real()?)) })()
            };
            match rec {
                Ok((a, b)) => {
                    if let Some(a) = a {
                        qd[nqd1] = a;
                    }
                    if let Some(b) = b {
                        data[nqd1] = b;
                    }
                    ierr = 0;
                }
                Err(e) => ierr = e.iostat(),
            }
            err[nqd1] = 0.0;
        } else if nrtype == 2 || nrtype >= 5 {
            let rec = {
                let mut r = f.list();
                (|| -> IoRes<(Option<f32>, Option<f32>, Option<f32>)> { Ok((r.real()?, r.real()?, r.real()?)) })()
            };
            match rec {
                Ok((a, b, c)) => {
                    if let Some(a) = a {
                        qd[nqd1] = a;
                    }
                    if let Some(b) = b {
                        data[nqd1] = b;
                    }
                    if let Some(c) = c {
                        err[nqd1] = c;
                    }
                    ierr = 0;
                }
                Err(e) => ierr = e.iostat(),
            }
        } else {
            match read_genie_record(f) {
                Ok((_n, a, b, c)) => {
                    qd[nqd1] = a;
                    data[nqd1] = b;
                    err[nqd1] = c;
                    ierr = 0;
                }
                Err(e) => ierr = e.iostat(),
            }
        }
        if ierr == 0 {
            nqd = nqd1;
        }
    }
    // Ignore the last point which may be rubbish
    if nrtype > 2 {
        nqd -= 1;
    }
    nqd
}
