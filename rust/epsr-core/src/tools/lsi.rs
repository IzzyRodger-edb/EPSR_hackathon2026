//! `lsi` — translation of `fortran_src/lsi.f` (`run_lsi`, `find_lsi`): local
//! structure index (Shi, Russo & Tanaka) of bonded / non-bonded neighbours,
//! plus distributions of neighbour counts, r3..r6 and neighbour triangles.
//!
//! STATUS: translated, verified byte-identical (EPSR-test-ex66, ex67).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> lsi <input file>`
//! Keyword table: `crate::symtab_tools::LSI` (the `nsize` text and default are
//! computed from `mchain` as in `setup_symbols_lsi`).
//!
//! Outputs: `<run>.LSI.n01` (10 columns per distribution) and the rewritten
//! `<run>.LSI.dat` (with the running averages `nbonded`, `nnonbonded`,
//! `avbond`).  Quirks of the Fortran reproduced here: `gofz` is zeroed
//! (whole array) at the start of every distribution, so only the last
//! distribution keeps non-zero histograms; the running averages use
//! `fc = 1/nsumt` with the *input* `nsumt` (so runs 1 and 2 both use fc = 1);
//! a neighbour for which no atom pair matched has `iwmin = 0`, and the
//! Fortran then reads `dfs1(id,0)`/`dfs2(id,0)`, i.e. `dff2(id,mdist)` and
//! `dfs1(id,mdist)` (both 0 unless 30 criteria are given) — emulated.

use super::chains::{fmatch, read_criteria, set_nsize_default, upto_blank, Crit};
use super::clusters::{list_int, open_input_with, outgofz, rd_int};
use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::fmath::nint;
use crate::fortran::format::{fmt, pad_to, FV};
use crate::fortran::listread::ListRead;
use crate::fortran::text::{fstr, parse, substr};
use crate::fortran::{F1, F2};
use crate::neighbour::minimagec;
use crate::state::State;

// lsi.inc
const MDIST: usize = 30;
const MZ: usize = 200;
const MCHAIN: i32 = 200;
const MBOND: i32 = 20;

/// Running averages read from / written to the input file (per distribution).
#[derive(Default, Clone)]
struct Saves {
    navbondedsave: Vec<f32>,
    navnonbondedsave: Vec<f32>,
    avbondsave: Vec<f32>,
    navbondedsave2: Vec<f32>,
    navnonbondedsave2: Vec<f32>,
    avbondsave2: Vec<f32>,
}

/// `read(symbol_value(iref),*) a[,b]`
fn read_pair(s: &str, two: bool) -> Result<(f32, f32)> {
    let mut r = ListRead::from_str(s);
    let a = r.real()?.ok_or_else(|| EpsrError::Input(format!("cannot read {s:?}")))?;
    let b = if two { r.real()?.ok_or_else(|| EpsrError::Input(format!("cannot read {s:?}")))? } else { 0.0 };
    Ok((a, b))
}

/// `dfs1(id,iw)`, `dfs2(id,iw)` including the out-of-bounds `iw = 0` case
/// (column-major storage in COMMON /lsi_r/: dff1, dff2, dfs1, dfs2, ...).
fn dfs_range(c: &Crit, id: usize, iw: i32) -> (f32, f32) {
    if iw >= 1 {
        (c.dfs1[id][iw as usize], c.dfs2[id][iw as usize])
    } else {
        (c.dff2[id][MDIST], c.dfs1[id][MDIST])
    }
}

/// Run the `lsi` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    let fix = |st: &mut State| set_nsize_default(st, fmt("('Number of LSI values (max ',i3,')')", &[FV::I(MCHAIN as i64)]), MCHAIN);
    if !open_input_with(ctx, &mut st, &crate::symtab_tools::LSI, &fix)? {
        return Ok(());
    }
    let mut iref: i32 = 2;
    let mut nz = rd_int(&st.sym.symbol_value[iref])?;
    if nz > MZ as i32 {
        nz = MZ as i32;
    }
    iref += 1;
    let (lsimin, lsimax) = read_pair(&st.sym.symbol_value[iref], true)?;
    iref += 1;
    let mut nsumt = rd_int(&st.sym.symbol_value[iref])?;
    let mut ndist = st.sym.nfiles;
    if ndist > (MDIST / 6) as i32 {
        ndist = (MDIST / 6) as i32;
    }
    let ndist3 = 10 * ndist;
    if ndist3 as usize > MDIST {
        return crate::error::unsupported("lsi with more than 3 distributions (10*ndist columns overflow gofz(mz,mdist) in the Fortran)");
    }
    let nd = ndist as usize + 1;
    let mut sv = Saves {
        navbondedsave: vec![0.0; nd],
        navnonbondedsave: vec![0.0; nd],
        avbondsave: vec![0.0; nd],
        navbondedsave2: vec![0.0; nd],
        navnonbondedsave2: vec![0.0; nd],
        avbondsave2: vec![0.0; nd],
    };
    let c = match read_criteria(&st, ndist, MDIST, "lsi")? {
        Some(c) => c,
        None => return Ok(()),
    };
    let mut rtrimax = 0.0f32;
    for id in 1..=ndist {
        let idu = id as usize;
        let mut iref = st.sym.ireff(st.sym.nsymboluptodata + 12, id);
        let sometext = st.sym.symbol_value[iref].clone();
        let w = parse(&sometext, 5);
        if w.nwords() < 1 {
            return crate::error::unsupported("lsi: empty rtrimax value");
        }
        rtrimax = super::clusters::rd_real(w.word(&sometext, 1))?;
        let rd3 = |iref: i32, a: &mut f32, b: &mut f32| -> Result<()> {
            let v = st.sym.symbol_value[iref].clone();
            let w = parse(&v, 5);
            if w.nwords() >= 2 {
                let (x, y) = read_pair(&v, true)?;
                *a = x;
                *b = y * y;
            } else {
                let (x, _) = read_pair(&v, false)?;
                *a = x;
                *b = 0.0;
            }
            Ok(())
        };
        iref += 1;
        rd3(iref, &mut sv.navbondedsave[idu], &mut sv.navbondedsave2[idu])?;
        iref += 1;
        rd3(iref, &mut sv.navnonbondedsave[idu], &mut sv.navnonbondedsave2[idu])?;
        iref += 1;
        rd3(iref, &mut sv.avbondsave[idu], &mut sv.avbondsave2[idu])?;
    }
    let zstep = (lsimax - lsimin) / (nz - 1) as f32;
    let zhalf = 0.5 * zstep;
    let mut z = F1::<f32>::new(MZ);
    for j in 1..=nz {
        z[j] = lsimin + (j - 1) as f32 * zstep;
    }
    let mut gofz = F2::<f32>::new(MZ, MDIST);
    let mut gzt = F2::<f32>::new(MZ, MDIST);
    let mut gzt2 = F2::<f32>::new(MZ, MDIST);

    st.normalise_coordinates();
    let fc: f64 = if nsumt > 0 { 1.0 / nsumt as f64 } else { 1.0 };
    let fc1 = 1.0 - fc;

    let p = LsiPar { nz, ndist, lsimin, lsimax, zhalf, rtrimax };
    find_lsi(&st, &c, &p, fc, fc1, &mut sv, &mut gofz);

    let titler = "lsi";
    let mut titledat = vec![String::new(); MDIST];
    let mut idref10 = 0usize;
    for id in 1..=ndist as usize {
        let idref1 = idref10 + 1;
        idref10 = idref1 + 9;
        let a1 = upto_blank(&c.atom1[id][1]);
        let a2 = upto_blank(&c.atom2[id][1]);
        let a3 = upto_blank(&c.atom3[id][1]);
        let ilent1 = a1.len() as i32;
        let a3l1 = if crate::fortran::format::pad_to(&c.atom1[id][1], 4).contains(' ') { substr(&pad_to(&c.atom3[id][1], 4), 1, ilent1) } else { String::new() };
        let t = |s: String| fstr(&s, 29);
        titledat[idref1 - 1] = t(fmt("(4x,a,'-',a,'-',a)", &[FV::S(&a1), FV::S(&a2), FV::S(&a3)]));
        titledat[idref1] = titledat[idref1 - 1].clone();
        titledat[idref1 + 1] = titledat[idref1 - 1].clone();
        for k in 0..3 {
            let suffix = format!("{}_{}", a3, k + 1);
            titledat[idref1 + 2 + k] = t(fmt("(2x,a,'-',a,'-',a)", &[FV::S(&a3), FV::S(&a1), FV::S(&suffix)]));
        }
        for (k, r) in ["r3", "r4", "r5", "r6"].iter().enumerate() {
            titledat[idref1 + 5 + k] = t(fmt("(4x,a,'-',a,'-',a2)", &[FV::S(&a1), FV::S(&a3l1), FV::S(r)]));
        }
    }
    let fnameinp = st.sym.fnameinp.clone();
    outgofz(&st, &fnameinp, "n", nz, ndist3, &mut nsumt, &mut z, &gofz, &mut gzt, &mut gzt2, titler, &titledat)?;
    for id in 1..=ndist as usize {
        sv.navbondedsave2[id] = sv.navbondedsave2[id].sqrt();
        sv.navnonbondedsave2[id] = sv.navnonbondedsave2[id].sqrt();
        sv.avbondsave2[id] = sv.avbondsave2[id].sqrt();
    }
    let icount = st.sym.search_for_symbol("nsize", 0);
    let mut iref = st.sym.ireff(icount, 1);
    st.sym.symbol_value[iref] = list_int(nz);
    iref += 2;
    st.sym.symbol_value[iref] = list_int(nsumt);
    iref += 1;
    st.sym.symbol_value[iref] = list_int(ndist);
    for id in 1..=ndist {
        let idu = id as usize;
        let mut iref = st.sym.ireff(st.sym.nsymboluptodata + 13, id);
        st.sym.symbol_value[iref] = fmt("(2(1x,e12.5))", &[FV::R(sv.navbondedsave[idu]), FV::R(sv.navbondedsave2[idu])]);
        iref += 1;
        st.sym.symbol_value[iref] = fmt("(2(1x,e12.5))", &[FV::R(sv.navnonbondedsave[idu]), FV::R(sv.navnonbondedsave2[idu])]);
        iref += 1;
        st.sym.symbol_value[iref] = fmt("(2(1x,e12.5))", &[FV::R(sv.avbondsave[idu]), FV::R(sv.avbondsave2[idu])]);
    }
    st.save_inp_file()?;
    Ok(())
}

struct LsiPar {
    nz: i32,
    ndist: i32,
    lsimin: f32,
    lsimax: f32,
    zhalf: f32,
    rtrimax: f32,
}

/// `subroutine find_lsi(fc,fc1,navbondedsave,...)`
fn find_lsi(st: &State, c: &Crit, p: &LsiPar, fc: f64, fc1: f64, sv: &mut Saves, gofz: &mut F2<f32>) {
    let m = &st.mol;
    let nmol = m.nmol;
    let nz = p.nz;
    let lsimin = p.lsimin;
    let zhalf = p.zhalf;
    let rtrimax = p.rtrimax;
    let ty = |ia: i32| -> &str { &m.type_[m.itype[ia]] };
    let zstep = (p.lsimax - lsimin) / (nz - 1) as f32;
    let pi = std::f32::consts::PI;
    let dmax = m.dmax;
    let dmaxsqr = dmax * dmax;
    let nmu = nmol as usize + 1;
    let mb = MBOND as usize;
    // COMMON /lsi_i/, /lsi_r/ arrays (persist over distributions)
    let mut neighb = F2::<i32>::new(nmu, mb);
    let mut nneighb = F1::<i32>::new(nmu);
    let mut ia1min = F2::<i32>::new(nmu, mb);
    let mut ia2min = F2::<i32>::new(nmu, mb);
    let mut iwmin = F2::<i32>::new(nmu, mb);
    let mut dneighb = F2::<f32>::new(nmu, mb);
    let mut bneighb = F2::<bool>::new(nmu, mb);
    // local arrays
    let mut lsi = F1::<f32>::new(nmu);
    let mut r3 = F1::<f32>::new(nmu);
    let mut r4 = F1::<f32>::new(nmu);
    let mut r5 = F1::<f32>::new(nmu);
    let mut r6 = F1::<f32>::new(nmu);
    let mut avbondmol = F1::<f32>::new(nmu);
    let mut nneighbbond = F1::<f32>::new(nmu);
    let mut nneighbnonbond = F1::<f32>::new(nmu);
    let mut idref10 = 0;
    for id in 1..=p.ndist {
        let idu = id as usize;
        let idref1 = idref10 + 1;
        let idref2 = idref1 + 1;
        let idref3 = idref2 + 1;
        let idref4 = idref3 + 1;
        let idref5 = idref4 + 1;
        let idref6 = idref5 + 1;
        let idref7 = idref6 + 1;
        let idref8 = idref7 + 1;
        let idref9 = idref8 + 1;
        idref10 = idref9 + 1;
        gofz.fill(0.0);
        for i in 1..=nmol {
            nneighb[i] = 0;
            dneighb[(i, 1)] = 0.0;
            nneighbbond[i] = 0.0;
        }
        for im1 in 1..=nmol {
            let xc1 = m.xadd[im1];
            let yc1 = m.yadd[im1];
            let zc1 = m.zadd[im1];
            for im2 in 1..=nmol {
                if im1 == im2 {
                    continue;
                }
                let xc2 = m.xadd[im2];
                let yc2 = m.yadd[im2];
                let zc2 = m.zadd[im2];
                let mut mindist = dmax;
                let mut ia1save = 0;
                let mut ia2save = 0;
                let mut iwsave = 0;
                let mut iw = 0;
                while iw < c.nwordsave[idu] {
                    iw += 1;
                    let iwu = iw as usize;
                    let mut ia1 = m.naf[im1] - 1;
                    while ia1 < m.nal[im1] {
                        ia1 += 1;
                        if fmatch(ty(ia1), &c.atom1[idu][iwu]) {
                            let xa1 = m.xa[ia1] + xc1;
                            let ya1 = m.ya[ia1] + yc1;
                            let za1 = m.za[ia1] + zc1;
                            let mut ia2 = m.naf[im2] - 1;
                            while ia2 < m.nal[im2] {
                                ia2 += 1;
                                if fmatch(ty(ia2), &c.atom3[idu][iwu]) {
                                    let xa2 = m.xa[ia2] + xc2;
                                    let ya2 = m.ya[ia2] + yc2;
                                    let za2 = m.za[ia2] + zc2;
                                    let mut xdif = xa2 - xa1;
                                    let mut ydif = ya2 - ya1;
                                    let mut zdif = za2 - za1;
                                    minimagec(&mut xdif, &mut ydif, &mut zdif);
                                    let dsq13 = (xdif * xdif + ydif * ydif + zdif * zdif) * dmaxsqr;
                                    let da13 = dsq13.sqrt();
                                    if da13 < mindist {
                                        mindist = da13;
                                        ia1save = ia1;
                                        ia2save = ia2;
                                        iwsave = iw;
                                    }
                                }
                            }
                        }
                    }
                }
                let mut icount = 0;
                let mut foundneighb = false;
                while icount < nneighb[im1] && !foundneighb {
                    icount += 1;
                    if mindist < dneighb[(im1, icount)] {
                        let icountmax = MBOND.min(nneighb[im1] + 1);
                        let mut icount2 = icountmax;
                        while icount2 > icount {
                            let icount1 = icount2 - 1;
                            neighb[(im1, icount2)] = neighb[(im1, icount1)];
                            dneighb[(im1, icount2)] = dneighb[(im1, icount1)];
                            ia1min[(im1, icount2)] = ia1min[(im1, icount1)];
                            ia2min[(im1, icount2)] = ia2min[(im1, icount1)];
                            iwmin[(im1, icount2)] = iwmin[(im1, icount1)];
                            icount2 = icount1;
                        }
                        nneighb[im1] = icountmax;
                        neighb[(im1, icount)] = im2;
                        dneighb[(im1, icount)] = mindist;
                        ia1min[(im1, icount)] = ia1save;
                        ia2min[(im1, icount)] = ia2save;
                        iwmin[(im1, icount)] = iwsave;
                        foundneighb = true;
                    }
                }
                if !foundneighb && nneighb[im1] < MBOND {
                    nneighb[im1] += 1;
                    let n = nneighb[im1];
                    neighb[(im1, n)] = im2;
                    dneighb[(im1, n)] = mindist;
                    ia1min[(im1, n)] = ia1save;
                    ia2min[(im1, n)] = ia2save;
                    iwmin[(im1, n)] = iwsave;
                }
            }
            // which of the nearest neighbours of im1 are bonded to it
            for icount in 1..=nneighb[im1] {
                let iw = iwmin[(im1, icount)];
                let da13 = dneighb[(im1, icount)];
                let (d1, d2) = dfs_range(c, idu, iw);
                let found = da13 >= d1 && da13 <= d2;
                let mut foundia2 = false;
                if found {
                    let iwu = iw as usize;
                    let dsq13 = da13 * da13;
                    let ia1 = ia1min[(im1, icount)];
                    let im2 = neighb[(im1, icount)];
                    let ia3 = ia2min[(im1, icount)];
                    let xa3 = m.xadd[im2] + m.xa[ia3];
                    let ya3 = m.yadd[im2] + m.ya[ia3];
                    let za3 = m.zadd[im2] + m.za[ia3];
                    let mut ia2 = m.naf[im1] - 1;
                    while ia2 < m.nal[im1] && !foundia2 {
                        ia2 += 1;
                        foundia2 = fmatch(ty(ia2), &c.atom2[idu][iwu]);
                        if foundia2 {
                            let xdif = m.xa[ia2] - m.xa[ia1];
                            let ydif = m.ya[ia2] - m.ya[ia1];
                            let zdif = m.za[ia2] - m.za[ia1];
                            let dsq12 = (xdif * xdif + ydif * ydif + zdif * zdif) * dmaxsqr;
                            let da12 = dsq12.sqrt();
                            foundia2 = da12 >= c.dff1[idu][iwu] && da12 <= c.dff2[idu][iwu];
                            if foundia2 {
                                let xa2 = xc1 + m.xa[ia2];
                                let ya2 = yc1 + m.ya[ia2];
                                let za2 = zc1 + m.za[ia2];
                                let mut xdif = xa3 - xa2;
                                let mut ydif = ya3 - ya2;
                                let mut zdif = za3 - za2;
                                minimagec(&mut xdif, &mut ydif, &mut zdif);
                                let dsq23 = (xdif * xdif + ydif * ydif + zdif * zdif) * dmaxsqr;
                                let da23 = dsq23.sqrt();
                                foundia2 = da23 >= c.dft1[idu][iwu] && da23 <= c.dft2[idu][iwu];
                                if foundia2 {
                                    let anglefound = (dsq12 + dsq13 - dsq23) / (2.0 * da12 * da13);
                                    foundia2 = anglefound >= c.ang213min[idu][iwu] && anglefound <= c.ang213max[idu][iwu];
                                }
                            }
                        }
                    }
                }
                bneighb[(im1, icount)] = found && foundia2;
            }
        }
        // a non-bonded neighbour may be bonded the other way round
        for im1 in 1..=nmol {
            for icount1 in 1..=nneighb[im1] {
                if !bneighb[(im1, icount1)] {
                    let im2 = neighb[(im1, icount1)];
                    let mut icount2 = 0;
                    let mut found = false;
                    while icount2 < nneighb[im2] && !found {
                        icount2 += 1;
                        found = neighb[(im2, icount2)] == im1 && bneighb[(im2, icount2)];
                    }
                    bneighb[(im1, icount1)] = found;
                }
            }
            let mut dnonbondedmin = dmax;
            let mut dbondedav = 0.0f32;
            let mut nbonded = 0;
            for icount in 1..=nneighb[im1] {
                if bneighb[(im1, icount)] {
                    nbonded += 1;
                    dbondedav = ((nbonded - 1) as f32 * dbondedav + dneighb[(im1, icount)]) / nbonded as f32;
                } else if dneighb[(im1, icount)] < dnonbondedmin {
                    dnonbondedmin = dneighb[(im1, icount)];
                }
            }
            lsi[im1] = dnonbondedmin - dbondedav;
            r3[im1] = dneighb[(im1, 3)] - dbondedav;
            r4[im1] = dneighb[(im1, 4)] - dbondedav;
            r5[im1] = dneighb[(im1, 5)] - dbondedav;
            r6[im1] = dneighb[(im1, 6)] - dbondedav;
        }
        for im1 in 1..=nmol {
            nneighbbond[im1] = 0.0;
            nneighbnonbond[im1] = 0.0;
            avbondmol[im1] = 0.0;
            if nneighb[im1] >= 0 {
                for icount in 1..=nneighb[im1] {
                    if bneighb[(im1, icount)] {
                        nneighbbond[im1] += 1.0;
                        avbondmol[im1] += dneighb[(im1, icount)];
                    } else {
                        let iw = iwmin[(im1, icount)];
                        let da13 = dneighb[(im1, icount)];
                        let (d1, d2) = dfs_range(c, idu, iw);
                        if da13 >= d1 && da13 <= d2 {
                            nneighbnonbond[im1] += 1.0;
                        }
                    }
                }
            }
            if nneighbbond[im1] > 0.0 {
                avbondmol[im1] /= nneighbbond[im1];
            }
        }
        let mut navbonded = 0.0f32;
        let mut navnonbonded = 0.0f32;
        let mut avbond = 0.0f32;
        let mut lensum = 0;
        for im1 in 1..=nmol {
            navbonded += nneighbbond[im1];
            navnonbonded += nneighbnonbond[im1];
            if avbondmol[im1] > 0.0 {
                avbond += avbondmol[im1];
                lensum += 1;
            }
        }
        navbonded /= nmol as f32;
        navnonbonded /= nmol as f32;
        if lensum > 0 {
            avbond /= lensum as f32;
        }
        sv.navbondedsave[idu] = (fc * navbonded as f64 + fc1 * sv.navbondedsave[idu] as f64) as f32;
        sv.navnonbondedsave[idu] = (fc * navnonbonded as f64 + fc1 * sv.navnonbondedsave[idu] as f64) as f32;
        sv.avbondsave[idu] = (fc * avbond as f64 + fc1 * sv.avbondsave[idu] as f64) as f32;
        navbonded = 0.0;
        navnonbonded = 0.0;
        avbond = 0.0;
        lensum = 0;
        for im1 in 1..=nmol {
            let diff = nneighbbond[im1] - sv.navbondedsave[idu];
            navbonded += diff * diff;
            let diff = nneighbnonbond[im1] - sv.navnonbondedsave[idu];
            navnonbonded += diff * diff;
            if avbondmol[im1] > 0.0 {
                let diff = avbondmol[im1] - sv.avbondsave[idu];
                avbond += diff * diff;
                lensum += 1;
            }
        }
        navbonded /= nmol as f32;
        navnonbonded /= nmol as f32;
        if lensum > 0 {
            avbond /= lensum as f32;
        }
        sv.navbondedsave2[idu] = (fc * navbonded as f64 + fc1 * sv.navbondedsave2[idu] as f64) as f32;
        sv.navnonbondedsave2[idu] = (fc * navnonbonded as f64 + fc1 * sv.navnonbondedsave2[idu] as f64) as f32;
        sv.avbondsave2[idu] = (fc * avbond as f64 + fc1 * sv.avbondsave2[idu] as f64) as f32;

        let bin = |v: f32| -> i32 { nint((v - lsimin + zhalf) / zstep) + 1 };
        for im1 in 1..=nmol {
            let izref = 1 + nint(nneighbbond[im1]);
            gofz[(izref, idref1)] += 1.0;
            let izref = 1 + nint(nneighbnonbond[im1]);
            gofz[(izref, idref2)] += 1.0;
            for (v, idr) in [(lsi[im1], idref3), (r3[im1], idref7), (r4[im1], idref8), (r5[im1], idref9), (r6[im1], idref10)] {
                let izref = bin(v);
                if izref > 0 && izref < nz + 1 {
                    gofz[(izref, idr)] += 1.0;
                }
            }
        }
        // triangle distributions
        let astep = pi / nz as f32;
        let ahalf = 0.5 * astep;
        for im1 in 1..=nmol {
            for icount2 in 1..=3.min(nneighb[im1] - 1) {
                let da12 = dneighb[(im1, icount2)];
                let dsq12 = da12 * da12;
                let im2 = neighb[(im1, icount2)];
                let ia2 = ia2min[(im1, icount2)];
                let (qx, qy, qz) = atom_xyz(st, ia2);
                let xa2 = m.xadd[im2] + qx;
                let ya2 = m.yadd[im2] + qy;
                let za2 = m.zadd[im2] + qz;
                for icount3 in icount2 + 1..=4.min(nneighb[im1]) {
                    let da13 = dneighb[(im1, icount3)];
                    if (rtrimax > 0.0 && da12 <= rtrimax && da13 <= rtrimax) || rtrimax == 0.0 {
                        let dsq13 = da13 * da13;
                        let im3 = neighb[(im1, icount3)];
                        let ia3 = ia2min[(im1, icount3)];
                        let (qx, qy, qz) = atom_xyz(st, ia3);
                        let xa3 = m.xadd[im3] + qx;
                        let ya3 = m.yadd[im3] + qy;
                        let za3 = m.zadd[im3] + qz;
                        let mut xdif = xa3 - xa2;
                        let mut ydif = ya3 - ya2;
                        let mut zdif = za3 - za2;
                        minimagec(&mut xdif, &mut ydif, &mut zdif);
                        let dsq23 = (xdif * xdif + ydif * ydif + zdif * zdif) * dmaxsqr;
                        let _da23 = dsq23.sqrt(); // computed but unused in the Fortran
                        let mut anglefound = (dsq12 + dsq13 - dsq23) / (2.0 * da12 * da13);
                        anglefound = fmax(anglefound, -1.0);
                        anglefound = fmin(anglefound, 1.0);
                        let mut izref = nint((anglefound.acos() + ahalf) / astep);
                        izref = nz.min(izref);
                        izref = 1.max(izref);
                        if icount3 == 2 {
                            gofz[(izref, idref4)] += 1.0;
                        }
                        if icount3 == 3 {
                            gofz[(izref, idref5)] += 1.0;
                        }
                        if icount3 == 4 {
                            gofz[(izref, idref6)] += 1.0;
                        }
                    }
                }
            }
        }
        let mut sums = [0.0f32; 11];
        let idrefs = [0, idref1, idref2, idref3, idref4, idref5, idref6, idref7, idref8, idref9, idref10];
        for iz in 1..=nz {
            let betar = (pi * (iz as f32 - 0.5) / nz as f32).sin();
            gofz[(iz, idref4)] /= betar;
            gofz[(iz, idref5)] /= betar;
            gofz[(iz, idref6)] /= betar;
            for k in 1..=10 {
                sums[k] += gofz[(iz, idrefs[k])];
            }
        }
        sums[3] *= zstep;
        for k in 4..=6 {
            sums[k] = sums[k] * 180.0 / nz as f32;
        }
        for k in 7..=10 {
            sums[k] *= zstep;
        }
        for iz in 1..=nz {
            for k in 1..=10 {
                if sums[k] > 0.0 {
                    gofz[(iz, idrefs[k])] /= sums[k];
                }
            }
        }
    }
}

/// `xa(ia), ya(ia), za(ia)` including the out-of-bounds `ia = 0` read the
/// Fortran makes for neighbours without a matching atom pair (`ia2min = 0`):
/// in COMMON /molecom_r/ the word before `xa(1)` is `stepmi`, before `ya(1)`
/// is `xa(matom)` and before `za(1)` is `ya(matom)` (both 0).
fn atom_xyz(st: &State, ia: i32) -> (f32, f32, f32) {
    let m = &st.mol;
    if ia >= 1 {
        (m.xa[ia], m.ya[ia], m.za[ia])
    } else {
        (m.stepmi, 0.0, 0.0)
    }
}

/// gfortran `max(a,b)` for reals (maxss: returns b if a is NaN).
fn fmax(a: f32, b: f32) -> f32 {
    if a.is_nan() || b > a {
        b
    } else {
        a
    }
}
fn fmin(a: f32, b: f32) -> f32 {
    if a.is_nan() || b < a {
        b
    } else {
        a
    }
}
