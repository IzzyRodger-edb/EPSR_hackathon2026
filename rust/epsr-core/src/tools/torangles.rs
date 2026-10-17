//! `torangles` — translation of `fortran_src/torangles.f` (`run_torangles`,
//! `find_torangles`).
//!
//! Torsion-angle distributions within molecules: for every molecule whose
//! first atom has type `atom1`, the dihedral angle about the axis between
//! atoms `ltyp1`-`ltyp2` (atom numbers within the molecule), using the
//! mid-points of the atom pairs `ltyp3` (from the first axis atom) and
//! `ltyp4` (from the second).  Distributions over -180..180 degrees are
//! accumulated in `<run>.TOR.c01`; the individual angles of the last run are
//! listed per molecule in `toranglesvalues.txt`.
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> torangles <input file>`
//! Keyword table: `crate::symtab_tools::TORANGLES`.

use super::partials::anal_common::*;
use super::ToolContext;
use crate::ato::{scapro, vecpro};
use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, FV};
use crate::fortran::output::OutFile;
use crate::fortran::{F1, F2};
#[allow(unused_imports)]
use crate::state::State;

const MDIST: usize = crate::state::MFILE;
const MZ: usize = 200;

struct Tor {
    nz: i32,
    ndist: i32,
    nfirst: Vec<i32>,
    ltyp: Vec<i32>,
    ltyp1: Vec<i32>,
    ltyp2: Vec<[i32; 2]>,
    ltyp3: Vec<[i32; 2]>,
}

/// Run the `torangles` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    if !ctx.open_input(&mut st, &crate::symtab_tools::TORANGLES)? {
        return Ok(());
    }
    let mut iref = 2;
    let nz = sv_int(&st, iref)?.min(MZ as i32);
    iref += 1;
    let mut nsumt = sv_int(&st, iref)?;
    let nfiles = st.sym.nfiles;
    let mut ndist = nfiles;
    if ndist > MDIST as i32 {
        ndist = MDIST as i32;
    }
    if nfiles > MDIST as i32 {
        return crate::error::unsupported("torangles: more than mfile distributions");
    }
    let icount = st.sym.search_for_symbol("atom1", 0);
    let nd = MDIST + 1;
    let mut t = Tor {
        nz,
        ndist,
        nfirst: vec![0; nd],
        ltyp: vec![0; nd],
        ltyp1: vec![0; nd],
        ltyp2: vec![[0; 2]; nd],
        ltyp3: vec![[0; 2]; nd],
    };
    let mut atom1 = vec![String::new(); nd];
    let mut atom2 = vec![String::new(); nd];
    let mut atom3 = vec![String::new(); nd];
    for id in 1..=nfiles {
        let mut iref = st.sym.ireff(icount, id);
        let mut ifound = 0;
        let mut it = 0;
        while it < st.mol.ntype && ifound == 0 {
            it += 1;
            if type_eq(&st, it, &st.sym.symbol_value[iref]) {
                ifound = it;
            }
        }
        if ifound == 0 {
            st.log(&fmt("(/'torangles> atom not matched to ','.ato file for distribution ',i3)", &[FV::I(id as i64)]));
            return Ok(());
        }
        let i = id as usize;
        t.nfirst[i] = ifound;
        atom1[i] = st.mol.type_[ifound].clone();
        iref += 1;
        t.ltyp[i] = sv_int(&st, iref)?;
        iref += 1;
        t.ltyp1[i] = sv_int(&st, iref)?;
        let m = &st.mol;
        let mut im = 0;
        let mut ifound = 0;
        while im < m.nmol && ifound == 0 {
            im += 1;
            if atom1[i] == m.type_[m.itype[m.naf[im]]] {
                ifound = im;
                atom2[i] = m.type_[m.itype[m.naf[im] + t.ltyp[i] - 1]].clone();
                atom3[i] = m.type_[m.itype[m.naf[im] + t.ltyp1[i] - 1]].clone();
            }
        }
        iref += 1;
        let (a, b) = sv_ints2(&st, iref)?;
        t.ltyp2[i] = [a, b];
        iref += 1;
        let (a, b) = sv_ints2(&st, iref)?;
        t.ltyp3[i] = [a, b];
    }
    let zstep = 360.0f32 / nz as f32;
    let zhalf = 0.5 * zstep;
    let mut z: F1<f32> = F1::new(MZ);
    for j in 1..=nz {
        z[j] = zstep * j as f32 - zhalf - 180.0;
    }
    let mut gofz: F2<f32> = F2::new(MZ, MDIST);
    st.normalise_coordinates();
    find_torangles(ctx, &st, &t, &mut gofz)?;

    let titler = "Phi[deg]";
    let ifmt = |v: i32| -> Result<&'static str> {
        if v < 10 {
            Ok("i1")
        } else if v < 100 {
            Ok("i2")
        } else if v < 1000 {
            Ok("i3")
        } else {
            crate::error::unsupported("torangles: atom number >= 1000 (format left undefined in the Fortran)")
        }
    };
    let mut titledat = Vec::new();
    for id in 1..=ndist as usize {
        let formatout = format!("(3x,a,2h:[,a,{},1h-,a,{},1h])", ifmt(t.ltyp[id])?, ifmt(t.ltyp1[id])?);
        let s = fmt(
            &formatout,
            &[
                FV::S(upto_blank(&atom1[id])),
                FV::S(upto_blank(&atom2[id])),
                FV::I(t.ltyp[id] as i64),
                FV::S(upto_blank(&atom3[id])),
                FV::I(t.ltyp1[id] as i64),
            ],
        );
        titledat.push(crate::fortran::fstr(&s, 29));
    }
    let mut gzt: F2<f32> = F2::new(MZ, MDIST);
    let mut gzt2: F2<f32> = F2::new(MZ, MDIST);
    outgofz(&st, "c", nz, ndist, &mut nsumt, &mut z, &gofz, &mut gzt, &mut gzt2, titler, &titledat)?;

    let icount = st.sym.search_for_symbol("nsize", 0);
    let mut iref = st.sym.ireff(icount, 1);
    set_sv_int(&mut st, iref, nz);
    iref += 1;
    set_sv_int(&mut st, iref, nsumt);
    iref += 1;
    set_sv_int(&mut st, iref, ndist);
    st.save_inp_file()
}

/// `subroutine find_torangles()`.
fn find_torangles(ctx: &ToolContext, st: &State, t: &Tor, gofz: &mut F2<f32>) -> Result<()> {
    let m = &st.mol;
    let pi = 3.141592653f32;
    let nz = t.nz;
    let ndist = t.ndist;
    let zstep = 360.0f32 / nz as f32;
    let mut ifindm = 0;
    let mut out = OutFile::create(&ctx.path("toranglesvalues.txt"));
    let formatout = format!("(1x,i6,{}(1x,e14.7))", ndist); // Fortran: (a7,i4,a11)
    let mut phisave = vec![0.0f32; ndist as usize + 1];
    for im1 in 1..=m.nmol {
        for v in phisave.iter_mut() {
            *v = 0.0;
        }
        for id in 1..=ndist as usize {
            if m.itype[m.naf[im1]] == t.nfirst[id] {
                ifindm += 1;
                let ia1 = m.naf[im1] + t.ltyp[id] - 1;
                let ia2 = m.naf[im1] + t.ltyp1[id] - 1;
                let (x1, y1, z1) = (m.xa[ia1], m.ya[ia1], m.za[ia1]);
                let (x2, y2, z2) = (m.xa[ia2], m.ya[ia2], m.za[ia2]);
                let mut cd = [x2 - x1, y2 - y1, z2 - z1];
                let dmag = scapro(&cd, &cd).sqrt();
                for v in cd.iter_mut() {
                    *v /= dmag;
                }
                let io1a = m.naf[im1] + t.ltyp2[id][0] - 1;
                let io1b = m.naf[im1] + t.ltyp2[id][1] - 1;
                let io2a = m.naf[im1] + t.ltyp3[id][0] - 1;
                let io2b = m.naf[im1] + t.ltyp3[id][1] - 1;
                let co1 = [
                    0.5 * (m.xa[io1a] + m.xa[io1b]) - x1,
                    0.5 * (m.ya[io1a] + m.ya[io1b]) - y1,
                    0.5 * (m.za[io1a] + m.za[io1b]) - z1,
                ];
                let co2 = [
                    0.5 * (m.xa[io2a] + m.xa[io2b]) - x2,
                    0.5 * (m.ya[io2a] + m.ya[io2b]) - y2,
                    0.5 * (m.za[io2a] + m.za[io2b]) - z2,
                ];
                let mut cv1 = vecpro(&co1, &cd);
                let mut cv2 = vecpro(&co2, &cd);
                let dmag = scapro(&cv1, &cv1).sqrt();
                for v in cv1.iter_mut() {
                    *v /= dmag;
                }
                let dmag = scapro(&cv2, &cv2).sqrt();
                for v in cv2.iter_mut() {
                    *v /= dmag;
                }
                let mut cosphi = scapro(&cv1, &cv2);
                if cosphi > 1.0 {
                    cosphi = 1.0;
                }
                if cosphi < -1.0 {
                    cosphi = -1.0;
                }
                let cpar = vecpro(&cv1, &cv2);
                let sinphi = scapro(&cpar, &cd);
                let phis = sinphi.atan2(cosphi) * 180.0 / pi;
                phisave[id] = phis;
                let iref = ((phis + 180.0) / zstep) as i32 + 1;
                add_flat(gofz, iref, id as i32, 1.0);
            }
        }
        let mut v: Vec<FV> = vec![FV::I(im1 as i64)];
        for id in 1..=ndist as usize {
            v.push(FV::R(phisave[id]));
        }
        out.line(&fmt(&formatout, &v));
    }
    out.close().map_err(|e| EpsrError::Io(format!("toranglesvalues.txt: {e}")))?;
    for id in 1..=ndist {
        let mut sum0 = 0.0f32;
        for i in 1..=nz {
            sum0 += gofz[(i, id)];
        }
        sum0 /= nz as f32;
        if sum0 != 0.0 {
            for i in 1..=nz {
                gofz[(i, id)] /= sum0;
            }
        }
    }
    st.log(&fmt("(1x,'No of molecules found = ',i5)", &[FV::I(ifindm as i64)]));
    Ok(())
}
