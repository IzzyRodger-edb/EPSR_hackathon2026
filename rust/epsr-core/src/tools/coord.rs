//! `coord` — translation of `fortran_src/coord.f` (`run_coord`, `find_coord`).
//!
//! Coordination-number distributions: for each distribution, the histogram
//! of the number of `atom2` atoms within [rmin, rmax] of each `atom1` atom
//! (column 1) and the mean bond length for each coordination number
//! (column 2), accumulated over runs in `<run>.COORD.n01`.  Optionally only
//! `atom1` atoms within `rminrmaxc` of a centre atom type (`atom-c`, or
//! `atom-z`/`atom-r` for distances along z / in the xy-plane) are used.
//! The mean coordination numbers are written to `fort.14`.
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> coord <input file>`
//! Keyword table: `crate::symtab_tools::COORD`.

use super::partials::anal_common::*;
use super::ToolContext;
use crate::error::Result;
use crate::fortran::format::{fmt, FV};
use crate::fortran::{F1, F2};
use crate::neighbour::minimagec;
#[allow(unused_imports)]
use crate::state::State;

const MDIST: usize = 100;
const MZ: usize = 1000;

struct Coord {
    nz: i32,
    ndist: i32,
    ltypc: i32,
    dorhoz: bool,
    dorhor: bool,
    ltyp1: Vec<i32>,
    ltyp2: Vec<i32>,
    damin: Vec<f32>,
    damax: Vec<f32>,
    rminc: Vec<f32>,
    rmaxc: Vec<f32>,
}

/// Run the `coord` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    if !ctx.open_input(&mut st, &crate::symtab_tools::COORD)? {
        return Ok(());
    }
    let mut iref = 2;
    let nz = sv_int(&st, iref)?.min(MZ as i32);
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

    let mut ndist = st.sym.nfiles;
    if ndist > (MDIST / 2) as i32 {
        ndist = (MDIST / 2) as i32;
    }
    let ndist2 = 2 * ndist;
    let nd = ndist as usize + 1;
    let mut c = Coord {
        nz,
        ndist,
        ltypc,
        dorhoz,
        dorhor,
        ltyp1: vec![0; nd],
        ltyp2: vec![0; nd],
        damin: vec![0.0; nd],
        damax: vec![0.0; nd],
        rminc: vec![0.0; nd],
        rmaxc: vec![0.0; nd],
    };
    let mut atom1 = vec![String::new(); nd];
    let mut atom2 = vec![String::new(); nd];
    for id in 1..=ndist {
        let mut it = 0;
        let mut ifound1 = 0;
        let mut ifound2 = 0;
        let iref1 = st.sym.ireff(st.sym.nsymboluptodata + 1, id);
        let iref2 = iref1 + 1;
        let iref3 = iref2 + 1;
        let iref4 = iref3 + 1;
        while it < st.mol.ntype && (ifound1 == 0 || ifound2 == 0) {
            it += 1;
            if type_eq(&st, it, &st.sym.symbol_value[iref1]) {
                ifound1 = it;
            }
            if type_eq(&st, it, &st.sym.symbol_value[iref2]) {
                ifound2 = it;
            }
        }
        if ifound1 == 0 || ifound2 == 0 {
            return crate::error::unsupported("coord: atom type of a distribution not found in the .ato file (Fortran reads type(0))");
        }
        let i = id as usize;
        c.ltyp1[i] = ifound1;
        c.ltyp2[i] = ifound2;
        atom1[i] = st.mol.type_[ifound1].clone();
        atom2[i] = st.mol.type_[ifound2].clone();
        c.damin[i] = sv_real(&st, iref3)?;
        c.damax[i] = sv_real(&st, iref4)?;
        let (a, b) = sv_reals2(&st, iref4 + 1)?;
        c.rminc[i] = a;
        c.rmaxc[i] = b;
    }
    let mut z: F1<f32> = F1::new(MZ);
    for j in 1..=nz {
        z[j] = (j - 1) as f32;
    }
    let mut gofz: F2<f32> = F2::new(MZ, MDIST);
    st.normalise_coordinates();
    find_coord(&st, &c, &mut gofz);

    let titler = "Coordination_number";
    let mut titledat = vec![String::new(); MDIST];
    let mut idref2 = 0usize;
    for id in 1..=st.sym.nfiles as usize {
        let idref1 = idref2 + 1;
        idref2 = idref1 + 1;
        let t = fmt(
            "(1x,a,'-',a,'(',e9.3,',',e9.3,')')",
            &[FV::S(upto_blank(&atom1[id])), FV::S(upto_blank(&atom2[id])), FV::R(c.damin[id]), FV::R(c.damax[id])],
        );
        let t = crate::fortran::fstr(&t, 29);
        titledat[idref1 - 1] = t.clone();
        titledat[idref2 - 1] = t;
    }
    let mut gzt: F2<f32> = F2::new(MZ, MDIST);
    let mut gzt2: F2<f32> = F2::new(MZ, MDIST);
    outgofz(&st, "n", nz, ndist2, &mut nsumt, &mut z, &gofz, &mut gzt, &mut gzt2, titler, &titledat)?;

    let mut coordn = vec![(0.0f32, 0.0f32); nd];
    let mut idref2 = 0;
    for id in 1..=ndist as usize {
        let idref1 = idref2 + 1;
        idref2 = idref1 + 1;
        let mut sum1 = 0.0f32;
        let mut sum2 = 0.0f32;
        for i in 1..=nz {
            let term = z[i] * gzt[(i, idref1)];
            sum1 += term;
            sum2 += z[i] * term;
        }
        coordn[id].0 = sum1;
        let std = sum2 - sum1 * sum1;
        coordn[id].1 = if std > 0.0 { std.sqrt() } else { 0.0 };
    }

    let icount = st.sym.search_for_symbol("nsize", 0);
    let mut iref = st.sym.ireff(icount, 1);
    set_sv_int(&mut st, iref, nz);
    iref += 1;
    set_sv_int(&mut st, iref, nsumt);
    iref += 4;
    set_sv_int(&mut st, iref, ndist);
    st.save_inp_file()?;
    let mut lines = Vec::new();
    for id in 1..=ndist as usize {
        let a1 = FV::Ch(&atom1[id], 4);
        let a2 = FV::Ch(&atom2[id], 4);
        let l = if ltypc == 0 {
            fmt(
                "(1x,a3,1x,a3,4(1x,f10.5))",
                &[a1, a2, FV::R(c.damin[id]), FV::R(c.damax[id]), FV::R(coordn[id].0), FV::R(coordn[id].1)],
            )
        } else {
            fmt(
                "(1x,a3,1x,a3,4(1x,f10.5))",
                &[
                    a1,
                    a2,
                    FV::R(c.damin[id]),
                    FV::R(c.damax[id]),
                    FV::R(c.rminc[id]),
                    FV::R(c.rmaxc[id]),
                    FV::R(coordn[id].0),
                    FV::R(coordn[id].1),
                ],
            )
        };
        lines.push(l);
    }
    write_fort14(ctx, &lines)
}

/// `subroutine find_coord()`.
fn find_coord(st: &State, c: &Coord, gofz: &mut F2<f32>) {
    let m = &st.mol;
    let natom = m.natom;
    let dmax = m.dmax;
    let nz = c.nz;
    let mut mbond = vec![0i32; natom as usize + 1];
    let mut bav = vec![0.0f32; natom as usize + 1];
    for id in 1..=c.ndist as usize {
        let idref1 = 2 * id as i32 - 1;
        let idref2 = 2 * id as i32;
        let mut atom1list: Vec<i32> = Vec::new();
        let mut atom2list: Vec<i32> = Vec::new();
        for ia2 in 1..=natom {
            if m.itype[ia2] == c.ltyp1[id] {
                let mut carryon = true;
                if c.ltypc > 0 {
                    carryon = false;
                    let im2 = m.mref[ia2];
                    let x2 = m.xa[ia2] + m.xadd[im2];
                    let y2 = m.ya[ia2] + m.yadd[im2];
                    let z2 = m.za[ia2] + m.zadd[im2];
                    let mut ia4 = 0;
                    while ia4 < natom && !carryon {
                        ia4 += 1;
                        if m.itype[ia4] == c.ltypc {
                            let im4 = m.mref[ia4];
                            let mut xdif = x2 - (m.xa[ia4] + m.xadd[im4]);
                            let mut ydif = y2 - (m.ya[ia4] + m.yadd[im4]);
                            let mut zdif = z2 - (m.za[ia4] + m.zadd[im4]);
                            minimagec(&mut xdif, &mut ydif, &mut zdif);
                            let dsq12 = if c.dorhoz {
                                zdif * zdif
                            } else if c.dorhor {
                                xdif * xdif + ydif * ydif
                            } else {
                                xdif * xdif + ydif * ydif + zdif * zdif
                            };
                            let d12 = dsq12.sqrt() * dmax;
                            carryon = d12 >= c.rminc[id] && d12 < c.rmaxc[id];
                        }
                    }
                }
                if carryon {
                    atom1list.push(ia2);
                }
            }
            if m.itype[ia2] == c.ltyp2[id] {
                atom2list.push(ia2);
            }
        }
        for i in 1..=natom as usize {
            mbond[i] = 0;
            bav[i] = 0.0;
        }
        let mut momol = 0;
        let mut mbmol = 0;
        for &ia1 in &atom1list {
            let im1 = m.mref[ia1];
            momol += 1;
            let xa1 = m.xa[ia1] + m.xadd[im1];
            let ya1 = m.ya[ia1] + m.yadd[im1];
            let za1 = m.za[ia1] + m.zadd[im1];
            for &ia2 in &atom2list {
                let im2 = m.mref[ia2];
                let x2 = m.xa[ia2] + m.xadd[im2];
                let y2 = m.ya[ia2] + m.yadd[im2];
                let z2 = m.za[ia2] + m.zadd[im2];
                let mut xdif = x2 - xa1;
                let mut ydif = y2 - ya1;
                let mut zdif = z2 - za1;
                minimagec(&mut xdif, &mut ydif, &mut zdif);
                let dsq12 = xdif * xdif + ydif * ydif + zdif * zdif;
                let da12 = dsq12.sqrt() * dmax;
                if da12 >= c.damin[id] && da12 <= c.damax[id] {
                    mbond[ia1 as usize] += 1;
                    bav[ia1 as usize] += da12;
                    mbmol += 1;
                }
            }
        }
        for &ia1 in &atom1list {
            let ia1 = ia1 as usize;
            if mbond[ia1] > 0 {
                bav[ia1] /= mbond[ia1] as f32;
            } else {
                bav[ia1] = 0.0;
            }
            let iref = mbond[ia1] + 1;
            if iref > 0 && iref < nz {
                gofz[(iref, idref1)] += 1.0;
                gofz[(iref, idref2)] += bav[ia1];
            }
        }
        st.log(&format!(" {mbmol} {momol}"));
        let avbond;
        if momol > 0 {
            avbond = mbmol as f32 / momol as f32;
            for iref in 1..=nz {
                if gofz[(iref, idref1)] > 0.0 {
                    gofz[(iref, idref2)] /= gofz[(iref, idref1)];
                } else {
                    gofz[(iref, idref2)] = 0.0;
                }
                gofz[(iref, idref1)] /= momol as f32;
            }
        } else {
            avbond = 0.0;
        }
        st.log(&fmt(
            "(1x,'Average coord. number for type ',i2,' around type ',i2,' = ',f10.5)",
            &[FV::I(c.ltyp2[id] as i64), FV::I(c.ltyp1[id] as i64), FV::R(avbond)],
        ));
    }
}
