//! Helpers shared by the structural-analysis programs (`partials`, `coord`,
//! `triangles`, `torangles`, `voids`, `fluctuations`, `local_density`).
//!
//! Translation of `initgofz`/`outgofz` (`read_and_save_files.f`),
//! `initran1`/`saveran1` (`auxiliary_routines.f`, with the file handling of
//! the callers) and of the Fortran idioms those programs share (reading
//! keyword values, the implicit `fort.14` unit written after
//! `save_files_others` has closed unit 14, ...).

use crate::error::{EpsrError, Result};
use crate::fortran::format::ListWriter;
use crate::fortran::listread::{read_int, read_real, ListRead, TextFile};
use crate::fortran::output::OutFile;
use crate::fortran::text::{fstr, indexlastdot};
use crate::fortran::{F1, F2};
use crate::output::{r_int, w_int};
use crate::rng::{Ran1, NTAB};
use crate::state::State;

/// `read(symbol_value(iref),*) i`.
pub fn sv_int(st: &State, iref: i32) -> Result<i32> {
    read_int(&st.sym.symbol_value[iref]).map_err(|e| EpsrError::Input(format!("reading integer from {:?}: {e:?}", st.sym.symbol_value[iref])))
}
/// `read(symbol_value(iref),*) x`.
pub fn sv_real(st: &State, iref: i32) -> Result<f32> {
    read_real(&st.sym.symbol_value[iref]).map_err(|e| EpsrError::Input(format!("reading real from {:?}: {e:?}", st.sym.symbol_value[iref])))
}
/// `read(symbol_value(iref),*) x8` (real*8).
pub fn sv_dble(st: &State, iref: i32) -> Result<f64> {
    let s = &st.sym.symbol_value[iref];
    ListRead::from_str(s)
        .dble()
        .ok()
        .flatten()
        .ok_or_else(|| EpsrError::Input(format!("reading real*8 from {s:?}")))
}
/// `read(symbol_value(iref),*) a, b` (two reals).
pub fn sv_reals2(st: &State, iref: i32) -> Result<(f32, f32)> {
    let s = &st.sym.symbol_value[iref];
    let mut r = ListRead::from_str(s);
    let a = r.real().ok().flatten();
    let b = r.real().ok().flatten();
    match (a, b) {
        (Some(a), Some(b)) => Ok((a, b)),
        _ => Err(EpsrError::Input(format!("reading two reals from {s:?}"))),
    }
}
/// `read(symbol_value(iref),*) i, j` (two integers).
pub fn sv_ints2(st: &State, iref: i32) -> Result<(i32, i32)> {
    let s = &st.sym.symbol_value[iref];
    let mut r = ListRead::from_str(s);
    let a = r.int().ok().flatten();
    let b = r.int().ok().flatten();
    match (a, b) {
        (Some(a), Some(b)) => Ok((a, b)),
        _ => Err(EpsrError::Input(format!("reading two integers from {s:?}"))),
    }
}
/// `write(symbol_value(iref),*) i`.
pub fn set_sv_int(st: &mut State, iref: i32, v: i32) {
    let mut w = ListWriter::new();
    w.int(v);
    st.sym.symbol_value[iref] = w.finish().trim_end().to_string();
}
/// `write(symbol_value(iref),*) x`.
pub fn set_sv_real(st: &mut State, iref: i32, v: f32) {
    let mut w = ListWriter::new();
    w.real(v);
    st.sym.symbol_value[iref] = w.finish().trim_end().to_string();
}
/// `write(symbol_value(iref),*) x8`.
pub fn set_sv_dble(st: &mut State, iref: i32, v: f64) {
    let mut w = ListWriter::new();
    w.dble(v);
    st.sym.symbol_value[iref] = w.finish().trim_end().to_string();
}

/// `read(symbol_value(iref),'(a)') atomc` with `character*3 atomc`, then the
/// search `do while (it.lt.ntype.and.ltypc.eq.0) ... if(atomc.eq.type(it))`.
/// Returns the matched type (0 if none); `ltypc` is only searched when 0.
pub fn match_centre(st: &State, iref: i32, ltypc: &mut i32) -> bool {
    let atom = fstr(&st.sym.symbol_value[iref], 3);
    let mut it = 0;
    let mut found = false;
    while it < st.mol.ntype && *ltypc == 0 {
        it += 1;
        if atom == st.mol.type_[it] {
            *ltypc = it;
            found = true;
        }
    }
    found
}

/// Fortran `type(it).eq.symbol_value(iref)` (blank-padded comparison).
pub fn type_eq(st: &State, it: i32, value: &str) -> bool {
    st.mol.type_[it] == value.trim_end_matches(' ')
}

/// `index(atom,' ')-1` substring of a `character*4` atom label.
pub fn upto_blank(s: &str) -> &str {
    match s.find(' ') {
        Some(p) => &s[..p],
        None => s,
    }
}

/// `gofz(iref,id) = gofz(iref,id) + v` with the Fortran's column-major address
/// arithmetic: a row index beyond the declared first dimension (e.g. a bin
/// index `nz+1` for an angle of exactly 180 degrees) addresses the next column,
/// exactly as the unchecked Fortran does.
pub fn add_flat(g: &mut F2<f32>, iref: i32, id: i32, v: f32) {
    let off = (iref as i64 - 1) + (id as i64 - 1) * g.n1 as i64;
    if off < 0 || off as usize >= g.data.len() {
        panic!("histogram index ({iref},{id}) outside the array");
    }
    g.data[off as usize] += v;
}

/// `subroutine outgofz(fnameinp,ext,nz,ndist,nsum,z,gofz,gzt,gzt2,titler,titledat,...)`:
/// running average of the histograms `gofz` with the previous results read
/// from `<fnameinp up to last dot>.<ext>01`, written back to that file.
#[allow(clippy::too_many_arguments)]
pub fn outgofz(
    st: &State,
    ext: &str,
    nz: i32,
    ndist: i32,
    nsum: &mut i32,
    z: &mut F1<f32>,
    gofz: &F2<f32>,
    gzt: &mut F2<f32>,
    gzt2: &mut F2<f32>,
    titler: &str,
    titledat: &[String],
) -> Result<()> {
    let fnameinp = &st.sym.fnameinp;
    let ilen = indexlastdot(fnameinp).max(0) as usize;
    let fname = format!("{}{}", &fnameinp[..ilen], &ext[..1]);
    let log = st.logger();
    let mut ierr = 0;
    let mut lptin = nz;
    if *nsum > 0 {
        match r_int(&fname, ndist, z, gzt, gzt2, &log)? {
            Some(l) => lptin = l,
            None => ierr = 1,
        }
    }
    if ierr != 0 {
        *nsum = -1;
    }
    if *nsum > 0 {
        for id in 1..=ndist {
            for j in 1..=nz {
                if gzt2[(j, id)] < 1.0e-20 {
                    gzt2[(j, id)] = 0.0;
                }
                if gzt[(j, id)].abs() < 1.0e-20 {
                    gzt[(j, id)] = 0.0;
                }
                if j <= lptin {
                    gzt2[(j, id)] = gzt2[(j, id)] * gzt2[(j, id)] + gzt[(j, id)] * gzt[(j, id)];
                } else {
                    gzt[(j, id)] = 0.0;
                    gzt2[(j, id)] = 0.0;
                }
            }
        }
    } else {
        for id in 1..=ndist {
            for j in 1..=nz {
                gzt[(j, id)] = 0.0;
                gzt2[(j, id)] = 0.0;
            }
        }
    }
    let (fc1, fc2): (f64, f64);
    if *nsum >= 0 {
        *nsum += 1;
        fc1 = (1.0f32 / *nsum as f32) as f64;
        fc2 = 1.0 - fc1;
    } else {
        fc1 = 1.0;
        fc2 = 0.0;
    }
    let lpt = nz;
    for id in 1..=ndist {
        for i in 1..=nz {
            let a2 = fc2 * gzt[(i, id)] as f64;
            let a1 = fc1 * gofz[(i, id)] as f64;
            let a = a1 + a2;
            let b2 = fc2 * gzt2[(i, id)] as f64;
            let g = gofz[(i, id)];
            let b1 = fc1 * (g * g) as f64;
            let b = b1 + b2;
            let mut error = b - a * a;
            if error > 0.0 {
                error = error.sqrt();
            } else {
                error = 0.0;
            }
            gzt[(i, id)] = a as f32;
            gzt2[(i, id)] = error as f32;
        }
    }
    w_int(lpt, &fname, ndist, z, gzt, gzt2, titler, titledat, &log)
}

/// Lines written to unit 14 after `save_files_others` has closed it: gfortran
/// then implicitly opens `fort.14` in the current directory.
pub fn write_fort14(ctx: &super::ToolContext, lines: &[String]) -> Result<()> {
    let mut out = OutFile::create(&ctx.path("fort.14"));
    for l in lines {
        out.line(l);
    }
    out.close().map_err(|e| EpsrError::Io(format!("fort.14: {e}")))
}

/// `open(14,file='randomno.dat',status='unknown'); call initran1(14); close(14)`.
pub fn initran1_file(ctx: &super::ToolContext, st: &mut State) -> Result<()> {
    let path = ctx.path("randomno.dat");
    let rec = match TextFile::open(&path) {
        Ok(mut f) => {
            let mut r = f.list();
            (|| -> crate::fortran::IoRes<(i32, i32, i32, [i32; NTAB])> {
                let mcall = r.int()?.unwrap_or(0);
                let idum = r.int()?.unwrap_or(0);
                let iy = r.int()?.unwrap_or(0);
                let mut iv = [0i32; NTAB];
                for v in iv.iter_mut() {
                    *v = r.int()?.unwrap_or(0);
                }
                Ok((mcall, idum, iy, iv))
            })()
            .ok()
        }
        Err(_) => {
            // status='unknown' creates the (empty) file
            std::fs::write(&path, b"").map_err(|e| EpsrError::Io(format!("randomno.dat: {e}")))?;
            None
        }
    };
    match rec {
        Some(v) => st.rng = Ran1::from_record(Some(v)),
        None => {
            st.rng.idum = 0;
            st.rng.iy = 0;
            st.rng.mcall = 0;
        }
    }
    Ok(())
}

/// `open(14,file='randomno.dat',status='unknown'); call saveran1(14); close(14)`.
pub fn saveran1_file(ctx: &super::ToolContext, st: &State) -> Result<()> {
    let r = &st.rng;
    let mut w = ListWriter::new();
    w.int(r.mcall).int(r.idum).int(r.iy);
    for v in r.iv.iter() {
        w.int(*v);
    }
    let mut out = OutFile::create(&ctx.path("randomno.dat"));
    out.line(&w.finish());
    out.close().map_err(|e| EpsrError::Io(format!("randomno.dat: {e}")))
}

// ---------------------------------------------------------------------------
// gfortran's intrinsic SIND / COSD for real(4)
// ---------------------------------------------------------------------------
//
// With gfortran >= 14 (Fortran 2023), a call of `sind`/`cosd` in a program
// unit that does not declare them EXTERNAL resolves to the *intrinsic*
// (`_gfortran_sind_r4`, `_gfortran_cosd_r4`, libgfortran `trigd.inc`), not to
// the functions of the same name in `auxiliary_routines.f`.  `triangles.f` is
// such a unit.  The routines below reproduce libgfortran's algorithm
// (reconstructed from the gfortran 15 library): reduction modulo 360, exact
// values at multiples of 30 degrees, octant reduction and a two-part
// conversion constant pi/180 = PIO180H + PIO180L applied with an FMA.

/// gfortran intrinsic SIND/COSD (moved to `fortran::fmath`, also used by `moloutc`).
pub use crate::fortran::fmath::{cosd_intrinsic, sind_intrinsic};
#[cfg(test)]
mod tests {
    /// Bit patterns (x, sind(x), cosd(x)) produced by gfortran 15 intrinsics.
    /// (The translation was checked against 801001 such values.)
    #[test]
    fn trig_intrinsics_match_gfortran() {
        let v: &[(u32, u32, u32)] = &[
            (0xC4090000, 0x3E0E8365, 0xBF7D8235),
            (0xC404C597, 0xBE1EA59E, 0xBF7CE8B2),
            (0xC388FFD3, 0x3F7F6078, 0x3D8ECFF2),
            (0xC2B6AAE6, 0xBF7FEE3E, 0xBCBEAF1A),
            (0x00000000, 0x00000000, 0x3F800000),
            (0x3F000000, 0x3C0EF9BE, 0x3F7FFD81),
            (0x43551C3F, 0xBF0BD73E, 0xBF566E47),
            (0x43CD8000, 0x3F46F309, 0x3F211B24),
            (0x3DCC985F, 0x3AE488F0, 0x3F7FFFE6),
            (0x3DCCCCCC, 0x3AE4C37F, 0x3F7FFFE6),
        ];
        for &(x, s, c) in v {
            let x = f32::from_bits(x);
            assert_eq!(super::sind_intrinsic(x).to_bits(), s, "sind({x})");
            assert_eq!(super::cosd_intrinsic(x).to_bits(), c, "cosd({x})");
        }
    }
}
