//! `partials` — translation of `fortran_src/partials.f` (`run_partials`,
//! `find_partials`).
//!
//! Site-site radial distribution functions (intermolecular only) between all
//! pairs of atom types of an `.ato` configuration, optionally restricted to a
//! centre atom type (`atom-c`), or computed as a density profile along z
//! (`atom-z`) or in the xy-plane (`atom-r`).  Results are accumulated over
//! runs (`nsumt`) in `<run>.PARTIALS.g01`; the list of pairs is written to
//! `fort.14` (unit 14 after `save_files_others` has closed it).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> partials <input file>`
//! Keyword table: `crate::symtab_tools::PARTIALS`.

use super::ToolContext;
use crate::error::Result;
use crate::fortran::format::{fmt, FV};
use crate::fortran::{F1, F2};
use crate::neighbour::minimagec;
use crate::state::{State, MPAIRS};

#[path = "anal_common.rs"]
pub mod anal_common;
use anal_common::*;

const MDIST: usize = MPAIRS;
const MZ: usize = 1000;

/// Run the `partials` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    if !ctx.open_input(&mut st, &crate::symtab_tools::PARTIALS)? {
        return Ok(());
    }
    let mut iref = 2;
    let mut nz = sv_int(&st, iref)?;
    nz = nz.min(MZ as i32);
    iref += 1;
    let mut rnormalise = sv_real(&st, iref)?;
    if rnormalise > 1.0 {
        rnormalise = 1.0;
    }
    iref += 1;
    let mut nsumt = sv_int(&st, iref)?;
    iref += 1;
    let mut ltypc = 0;
    match_centre(&st, iref, &mut ltypc);
    iref += 1;
    let dorhoz = match_centre(&st, iref, &mut ltypc);
    iref += 1;
    let dorhor = match_centre(&st, iref, &mut ltypc);
    st.log(&fmt("(l2,l2,1x,i4)", &[FV::L(dorhoz), FV::L(dorhor), FV::I(ltypc as i64)]));

    let ntype = st.mol.ntype;
    let mut ndist: i32;
    let mut ltyp1 = vec![0i32; MDIST + 1];
    let mut ltyp2 = vec![0i32; MDIST + 1];
    let mut atom1 = vec![String::new(); MDIST + 1];
    let mut atom2 = vec![String::new(); MDIST + 1];
    ndist = 0;
    for it in 1..=ntype {
        for jt in 1..=ntype {
            if (ndist as usize) < MPAIRS {
                if ltypc > 0 {
                    if it == ltypc {
                        ndist += 1;
                        let n = ndist as usize;
                        ltyp1[n] = it;
                        ltyp2[n] = jt;
                        atom1[n] = st.mol.type_[it].clone();
                        atom2[n] = st.mol.type_[jt].clone();
                    }
                } else if jt >= it {
                    ndist += 1;
                    let n = ndist as usize;
                    ltyp1[n] = it;
                    ltyp2[n] = jt;
                    atom1[n] = st.mol.type_[it].clone();
                    atom2[n] = st.mol.type_[jt].clone();
                }
            }
        }
    }

    let dmax = st.mol.dmax;
    if dorhoz {
        rnormalise = rnormalise * st.lat.clatz / dmax;
    }
    if dorhor {
        rnormalise = rnormalise * st.lat.alatx.min(st.lat.blaty) / dmax;
    }
    let zstep = 0.5 * rnormalise / nz as f32;
    let zhalf = 0.5 * zstep;
    let mut z: F1<f32> = F1::new(MZ);
    for j in 1..=nz {
        z[j] = (zstep * j as f32 - zhalf) * dmax;
    }
    let mut gofz: F2<f32> = F2::new(MZ, MDIST);
    st.normalise_coordinates();
    find_partials(&st, nz, ndist, zstep, ltypc, dorhoz, dorhor, &ltyp1, &ltyp2, &z, &mut gofz);

    let titler = "r";
    let mut titledat = Vec::with_capacity(ndist as usize);
    for id in 1..=ndist as usize {
        titledat.push(crate::fortran::fstr(&format!("   {}-{}", upto_blank(&atom1[id]), upto_blank(&atom2[id])), 29));
    }
    let mut gzt: F2<f32> = F2::new(MZ, MDIST);
    let mut gzt2: F2<f32> = F2::new(MZ, MDIST);
    outgofz(&st, "g", nz, ndist, &mut nsumt, &mut z, &gofz, &mut gzt, &mut gzt2, titler, &titledat)?;

    let icount = st.sym.search_for_symbol("nz", 0);
    let mut iref = st.sym.ireff(icount, 1);
    set_sv_int(&mut st, iref, nz);
    iref += 2;
    set_sv_int(&mut st, iref, nsumt);
    iref += 4;
    set_sv_int(&mut st, iref, ndist);
    st.save_inp_file()?;
    let mut lines = Vec::new();
    for id in 1..=ndist as usize {
        lines.push(fmt("(1x,i5,5x,a3,1x,a3)", &[FV::I(id as i64), FV::Ch(&atom1[id], 4), FV::Ch(&atom2[id], 4)]));
    }
    write_fort14(ctx, &lines)
}

/// `subroutine find_partials()`.
#[allow(clippy::too_many_arguments)]
fn find_partials(
    st: &State,
    nz: i32,
    ndist: i32,
    zstep: f32,
    ltypc: i32,
    dorhoz: bool,
    dorhor: bool,
    ltyp1: &[i32],
    ltyp2: &[i32],
    z: &F1<f32>,
    gofz: &mut F2<f32>,
) {
    let m = &st.mol;
    let l = &st.lat;
    let pi = std::f32::consts::PI;
    let dmax = m.dmax;
    let (gfac, vol): (f32, f32) = if dorhoz {
        (1.0, 2.0 * l.alatx * l.blaty * dmax * dmax)
    } else if dorhor {
        (pi, l.clatz * dmax)
    } else {
        (4.0 * pi / 3.0, l.alatx * l.blaty * l.clatz * dmax * dmax * dmax)
    };
    let natom = m.natom;
    let mut lsum1 = vec![0i32; ndist as usize + 1];
    let mut lsum2 = vec![0i32; ndist as usize + 1];
    for id in 1..=ndist as usize {
        lsum1[id] = 0;
        lsum2[id] = 0;
        for ia1 in 1..=natom {
            if m.itype[ia1] == ltyp2[id] {
                lsum2[id] += 1;
            }
            if m.itype[ia1] == ltyp1[id] {
                let im1 = m.mref[ia1];
                let xc1 = m.xadd[im1];
                let yc1 = m.yadd[im1];
                let zc1 = m.zadd[im1];
                let x1 = m.xa[ia1] + xc1;
                let y1 = m.ya[ia1] + yc1;
                let z1 = m.za[ia1] + zc1;
                lsum1[id] += 1;
                for ia2 in 1..=natom {
                    let im2 = m.mref[ia2];
                    if m.itype[ia2] == ltyp2[id] && im1 != im2 {
                        let x2 = m.xa[ia2] + m.xadd[im2];
                        let y2 = m.ya[ia2] + m.yadd[im2];
                        let z2 = m.za[ia2] + m.zadd[im2];
                        let mut xdif = x2 - x1;
                        let mut ydif = y2 - y1;
                        let mut zdif = z2 - z1;
                        minimagec(&mut xdif, &mut ydif, &mut zdif);
                        let d12 = if dorhoz {
                            zdif.abs()
                        } else if dorhor {
                            (xdif * xdif + ydif * ydif).sqrt()
                        } else {
                            (xdif * xdif + ydif * ydif + zdif * zdif).sqrt()
                        };
                        let iref = (d12 / zstep) as i32 + 1;
                        if iref > 0 && iref <= nz {
                            gofz[(iref, id as i32)] += 1.0;
                        }
                    }
                }
            }
        }
    }
    for id in 1..=ndist as usize {
        let rho = if dorhoz || dorhor {
            lsum1[id] as f32 * vol
        } else if ltypc > 0 {
            lsum1[id] as f32
        } else {
            lsum1[id] as f32 * lsum2[id] as f32 / vol
        };
        st.log(&fmt("(1x,i5,5(1x,e13.6))", &[FV::I(id as i64), FV::R(lsum1[id] as f32), FV::R(lsum2[id] as f32), FV::R(rho)]));
        let gfac1 = gfac * rho;
        let mut z0 = 0.0f32;
        let mut z03 = 0.0f32;
        for iz in 1..=nz {
            let zm = z[iz];
            let zdiff = zm - z0;
            let z1 = zm + zdiff;
            let z13 = if dorhoz {
                z1
            } else if dorhor {
                z1 * z1
            } else {
                z1 * z1 * z1
            };
            let rat = gfac1 * (z13 - z03);
            gofz[(iz, id as i32)] /= rat;
            z0 = z1;
            z03 = z13;
        }
    }
}
