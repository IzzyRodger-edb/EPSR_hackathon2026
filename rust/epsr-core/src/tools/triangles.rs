//! `triangles` — translation of `fortran_src/triangles.f` (`run_triangles`,
//! `find_triangles`).
//!
//! Bond-angle distributions atom1-atom2-atom3 for atom pairs within
//! [rmin12, rmax12] and [rmin23, rmax23].  With `ltype 2` the third "atom"
//! is the average position of the atoms of type atom3 bonded to atom2.
//! Optionally the middle atom must lie within `rminrmaxc` of a centre atom
//! type (`atom-c`, `atom-z`, `atom-r`).  The distributions (normalised by
//! sin(theta) to unit area) are accumulated in `<run>.TRI.c01`; an order
//! parameter for each distribution is written to `fort.14`.
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> triangles <input file>`
//! Keyword table: `crate::symtab_tools::TRIANGLES`.

use super::partials::anal_common::*;
use super::ToolContext;
use crate::error::Result;
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::{F1, F2};
use crate::neighbour::minimagec;
#[allow(unused_imports)]
use crate::state::State;

const MDIST: usize = crate::state::MFILE;
const MZ: usize = 200;

struct Tri {
    nz: i32,
    ndist: i32,
    ltypc: i32,
    dorhoz: bool,
    dorhor: bool,
    ltyp: Vec<i32>,
    ltyp1: Vec<i32>,
    ltyp2: Vec<i32>,
    ltyp3: Vec<i32>,
    damin: Vec<[f32; 2]>,
    damax: Vec<[f32; 2]>,
    rminc: Vec<f32>,
    rmaxc: Vec<f32>,
}

/// Run the `triangles` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    if !ctx.open_input(&mut st, &crate::symtab_tools::TRIANGLES)? {
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

    let nfiles = st.sym.nfiles;
    let mut ndist = nfiles;
    if ndist > MDIST as i32 {
        ndist = MDIST as i32;
    }
    if nfiles > MDIST as i32 {
        return crate::error::unsupported("triangles: more than mfile distributions");
    }
    let icount1 = st.sym.search_for_symbol("atom1", 0);
    let icount2 = icount1 + 1;
    let icount3 = icount2 + 1;
    let nd = MDIST + 1;
    let mut t = Tri {
        nz,
        ndist,
        ltypc,
        dorhoz,
        dorhor,
        ltyp: vec![0; nd],
        ltyp1: vec![0; nd],
        ltyp2: vec![0; nd],
        ltyp3: vec![0; nd],
        damin: vec![[0.0; 2]; nd],
        damax: vec![[0.0; 2]; nd],
        rminc: vec![0.0; nd],
        rmaxc: vec![0.0; nd],
    };
    let mut atom1 = vec![String::new(); nd];
    let mut atom2 = vec![String::new(); nd];
    let mut atom3 = vec![String::new(); nd];
    let dmax = st.mol.dmax;
    for id in 1..=nfiles {
        let iref1 = st.sym.ireff(icount1, id);
        let iref2 = st.sym.ireff(icount2, id);
        let iref3 = st.sym.ireff(icount3, id);
        let (mut ifound1, mut ifound2, mut ifound3) = (0, 0, 0);
        let mut it = 0;
        while it < st.mol.ntype {
            it += 1;
            if type_eq(&st, it, &st.sym.symbol_value[iref1]) {
                ifound1 = it;
            }
            if type_eq(&st, it, &st.sym.symbol_value[iref2]) {
                ifound2 = it;
            }
            if type_eq(&st, it, &st.sym.symbol_value[iref3]) {
                ifound3 = it;
            }
        }
        if ifound1 == 0 || ifound2 == 0 || ifound3 == 0 {
            st.log(&fmt("(/'triangles> atoms not matched to ','.ato file for distribution ',i3)", &[FV::I(id as i64)]));
            return Ok(());
        }
        let i = id as usize;
        t.ltyp1[i] = ifound1;
        t.ltyp2[i] = ifound2;
        t.ltyp3[i] = ifound3;
        atom1[i] = st.mol.type_[ifound1].clone();
        atom2[i] = st.mol.type_[ifound2].clone();
        atom3[i] = st.mol.type_[ifound3].clone();
        let mut iref = iref3 + 1;
        t.ltyp[i] = sv_int(&st, iref)?;
        iref += 1;
        t.damin[i][0] = sv_real(&st, iref)?;
        iref += 1;
        t.damax[i][0] = sv_real(&st, iref)?;
        iref += 1;
        t.damin[i][1] = sv_real(&st, iref)?;
        iref += 1;
        t.damax[i][1] = sv_real(&st, iref)?;
        iref += 1;
        let (a, b) = sv_reals2(&st, iref)?;
        t.rminc[i] = a;
        t.rmaxc[i] = b;
        for k in 0..2 {
            t.damin[i][k] /= dmax;
            t.damax[i][k] /= dmax;
        }
        t.rminc[i] /= dmax;
        t.rmaxc[i] /= dmax;
    }
    let zstep = 180.0f32 / nz as f32;
    let zhalf = 0.5 * zstep;
    let mut z: F1<f32> = F1::new(MZ);
    for j in 1..=nz {
        z[j] = zstep * j as f32 - zhalf;
    }
    let mut gofz: F2<f32> = F2::new(MZ, MDIST);
    st.normalise_coordinates();
    find_triangles(&st, &mut t, &mut gofz);

    let titler = "theta[deg.]";
    let mut titledat = Vec::new();
    for id in 1..=ndist as usize {
        titledat.push(crate::fortran::fstr(
            &format!("   {}-{}-{}", upto_blank(&atom1[id]), upto_blank(&atom2[id]), upto_blank(&atom3[id])),
            29,
        ));
    }
    let mut gzt: F2<f32> = F2::new(MZ, MDIST);
    let mut gzt2: F2<f32> = F2::new(MZ, MDIST);
    outgofz(&st, "c", nz, ndist, &mut nsumt, &mut z, &gofz, &mut gzt, &mut gzt2, titler, &titledat)?;

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
        let mut sum1 = 0.0f32;
        let mut sum2 = 0.0f32;
        for iz in 1..=nz {
            let term = gzt[(iz, id as i32)] * sind_intrinsic(z[iz]);
            let mut q = cosd_intrinsic(z[iz]); // gfortran intrinsic COSD (see anal_common)
            if t.ltyp[id] == 1 {
                q += 1.0f32 / 3.0f32;
                q *= q;
            }
            sum1 += term;
            sum2 += q * term;
        }
        let rc = 0.5 * (t.rminc[id] + t.rmaxc[id]) * dmax;
        let mut w = ListWriter::new();
        w.int(id as i32).real(rc);
        if sum1 > 0.0 {
            if t.ltyp[id] == 1 {
                w.real(1.0 - (9.0 * sum2) / (4.0 * sum1));
            } else {
                w.real(sum2 / sum1);
            }
        } else {
            w.real(0.0);
        }
        lines.push(w.finish());
    }
    write_fort14(ctx, &lines)
}

/// `subroutine find_triangles()`.
fn find_triangles(st: &State, t: &mut Tri, gofz: &mut F2<f32>) {
    let m = &st.mol;
    let natom = m.natom;
    let nz = t.nz;
    let zstep = 180.0f32 / nz as f32;
    let piconv = 180.0f32 / 3.141592653f32;
    let ltypc = t.ltypc;
    for id in 1..=t.ndist as usize {
        let mut atom1list: Vec<i32> = Vec::new();
        let mut atom2list: Vec<i32> = Vec::new();
        let mut atom3list: Vec<i32> = Vec::new();
        for ia2 in 1..=natom {
            if m.itype[ia2] == t.ltyp1[id] {
                atom1list.push(ia2);
            }
            if m.itype[ia2] == t.ltyp3[id] {
                atom3list.push(ia2);
            }
            if m.itype[ia2] == t.ltyp2[id] {
                let mut carryon = true;
                if ltypc > 0 {
                    carryon = false;
                    let im2 = m.mref[ia2];
                    let x2 = m.xa[ia2] + m.xadd[im2];
                    let y2 = m.ya[ia2] + m.yadd[im2];
                    let z2 = m.za[ia2] + m.zadd[im2];
                    let mut ia4 = 0;
                    while ia4 < natom && !carryon {
                        ia4 += 1;
                        if m.itype[ia4] == ltypc {
                            let im4 = m.mref[ia4];
                            let mut xdif = x2 - (m.xa[ia4] + m.xadd[im4]);
                            let mut ydif = y2 - (m.ya[ia4] + m.yadd[im4]);
                            let mut zdif = z2 - (m.za[ia4] + m.zadd[im4]);
                            minimagec(&mut xdif, &mut ydif, &mut zdif);
                            let dsq12 = if t.dorhoz {
                                zdif * zdif
                            } else if t.dorhor {
                                xdif * xdif + ydif * ydif
                            } else {
                                xdif * xdif + ydif * ydif + zdif * zdif
                            };
                            let d12 = dsq12.sqrt();
                            carryon = d12 >= t.rminc[id] && d12 < t.rmaxc[id];
                        }
                    }
                }
                if carryon {
                    atom2list.push(ia2);
                }
            }
        }
        if t.ltyp1[id] == ltypc {
            t.damin[id][0] = t.rminc[id];
            t.damax[id][0] = t.rmaxc[id];
        }
        st.log(&format!(" {} {} {} {}", id, atom1list.len(), atom2list.len(), atom3list.len()));
        let centre1 = t.ltyp1[id] == ltypc;
        for &ia1 in &atom1list {
            let im1 = m.mref[ia1];
            let x1 = m.xa[ia1] + m.xadd[im1];
            let y1 = m.ya[ia1] + m.yadd[im1];
            let z1 = m.za[ia1] + m.zadd[im1];
            for &ia2 in &atom2list {
                let im2 = m.mref[ia2];
                let xc2 = m.xadd[im2];
                let yc2 = m.yadd[im2];
                let zc2 = m.zadd[im2];
                let x2 = m.xa[ia2] + xc2;
                let y2 = m.ya[ia2] + yc2;
                let z2 = m.za[ia2] + zc2;
                let mut xdif = x2 - x1;
                let mut ydif = y2 - y1;
                let mut zdif = z2 - z1;
                minimagec(&mut xdif, &mut ydif, &mut zdif);
                let dsq12 = if t.dorhoz && centre1 {
                    zdif * zdif
                } else if t.dorhor && centre1 {
                    xdif * xdif + ydif * ydif
                } else {
                    xdif * xdif + ydif * ydif + zdif * zdif
                };
                let d12 = dsq12.sqrt();
                if !(d12 >= t.damin[id][0] && d12 <= t.damax[id][0]) {
                    continue;
                }
                // closure computing the angle bin for a third point (x3,y3,z3)
                let add_angle = |x3: f32, y3: f32, z3: f32, dsq23: f32, d23: f32, gofz: &mut F2<f32>| {
                    let (mut xdif, mut ydif, mut zdif) = if t.dorhoz && centre1 {
                        (x3 - x2, y3 - y2, z3 - z1)
                    } else if t.dorhor && centre1 {
                        (x3 - x1, y3 - y1, z3 - z2)
                    } else {
                        (x3 - x1, y3 - y1, z3 - z1)
                    };
                    minimagec(&mut xdif, &mut ydif, &mut zdif);
                    let dsq13 = xdif * xdif + ydif * ydif + zdif * zdif;
                    let mut coss = (dsq12 + dsq23 - dsq13) / (2.0 * d12 * d23);
                    if coss > 1.0 {
                        coss = 1.0;
                    }
                    if coss < -1.0 {
                        coss = -1.0;
                    }
                    let deg = piconv * coss.acos();
                    let iref = (deg / zstep) as i32 + 1;
                    add_flat(gofz, iref, id as i32, 1.0);
                };
                if t.ltyp[id] == 1 {
                    for &ia3 in &atom3list {
                        if ia3 != ia2 && ia3 != ia1 {
                            let im3 = m.mref[ia3];
                            let x3 = m.xa[ia3] + m.xadd[im3];
                            let y3 = m.ya[ia3] + m.yadd[im3];
                            let z3 = m.za[ia3] + m.zadd[im3];
                            let mut xdif = x3 - x2;
                            let mut ydif = y3 - y2;
                            let mut zdif = z3 - z2;
                            minimagec(&mut xdif, &mut ydif, &mut zdif);
                            let dsq23 = xdif * xdif + ydif * ydif + zdif * zdif;
                            let d23 = dsq23.sqrt();
                            if d23 >= t.damin[id][1] && d23 <= t.damax[id][1] {
                                add_angle(x3, y3, z3, dsq23, d23, gofz);
                            }
                        }
                    }
                } else {
                    let mut x3 = 0.0f32;
                    let mut y3 = 0.0f32;
                    let mut z3 = 0.0f32;
                    let mut nf = 0;
                    for ib in 1..=m.nbond[ia2] {
                        let ia3 = m.iabond[(ia2, ib)];
                        if m.itype[ia3] == t.ltyp3[id] && ia3 != ia2 {
                            x3 += m.xa[ia3];
                            y3 += m.ya[ia3];
                            z3 += m.za[ia3];
                            nf += 1;
                        }
                    }
                    if nf > 0 {
                        x3 = x3 / nf as f32 + xc2;
                        y3 = y3 / nf as f32 + yc2;
                        z3 = z3 / nf as f32 + zc2;
                        let mut xdif = x3 - x2;
                        let mut ydif = y3 - y2;
                        let mut zdif = z3 - z2;
                        minimagec(&mut xdif, &mut ydif, &mut zdif);
                        let dsq23 = xdif * xdif + ydif * ydif + zdif * zdif;
                        let d23 = dsq23.sqrt();
                        add_angle(x3, y3, z3, dsq23, d23, gofz);
                    }
                }
            }
        }
    }
    for id in 1..=t.ndist {
        let mut sumg = 0.0f32;
        let mut beta = -0.5 * zstep;
        for iz in 1..=nz {
            beta += zstep;
            let betar = beta / piconv;
            gofz[(iz, id)] /= betar.sin();
            sumg += gofz[(iz, id)];
        }
        if sumg > 0.0 {
            for iz in 1..=nz {
                gofz[(iz, id)] /= sumg * zstep;
            }
        }
    }
}
