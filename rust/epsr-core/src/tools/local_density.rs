//! `local_density` — translation of `fortran_src/local_density.f`
//! (`run_clusters`, `find_clusters` of that file).
//!
//! Despite its name this is the multiple-bond variant of `clusters`: it
//! uses the CLUSTERS keyword table (`.CLUSTERS.dat` input), but `atom1`,
//! `atom2`, `rmin` and `rmax` of each distribution may be lists, so that a
//! molecule-molecule "bond" exists if *any* of the listed atom pairs is
//! within its distance range.  For each distribution the cluster-size
//! histogram (with the percolating fraction in the first row) is
//! accumulated in `<run>.CLUSTERS.n01`; the fraction of molecules in clusters
//! of 2 or more is written to `fort.14`.
//!
//! Invocation: `epsr-rs <folder> local_density <input file>`
//! Keyword table: `crate::symtab_tools::CLUSTERS`.

use super::partials::anal_common::*;
use super::ToolContext;
use crate::error::Result;
use crate::fortran::format::{fmt, FV};
use crate::fortran::listread::{read_real, read_word};
use crate::fortran::text::{parse, parse_n};
use crate::fortran::{F1, F2};
use crate::neighbour::minimagec;
#[allow(unused_imports)]
use crate::state::State;

const MDIST: usize = crate::state::MFILE;
/// `mz` of `clusters_multibond.inc`.
const MZ: usize = 6000;
const MAXNCBONDS: usize = 10;

struct Def {
    ltyp: i32,
    ltyp1: i32,
    rmin: f32,
    rmax: f32,
}

/// Run the `local_density` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    if !ctx.open_input(&mut st, &crate::symtab_tools::CLUSTERS)? {
        return Ok(());
    }
    let mut iref = 2;
    let nz = sv_int(&st, iref)?;
    if nz as usize > MZ {
        return crate::error::unsupported("local_density: nsize larger than the Fortran array size 6000");
    }
    iref += 1;
    let mut nsumt = sv_int(&st, iref)?;
    let mut ndist = st.sym.nfiles;
    if ndist > MDIST as i32 {
        ndist = MDIST as i32;
    }
    let mut defs: Vec<Vec<Def>> = vec![Vec::new()];
    let mut atoms: Vec<Vec<(String, String)>> = vec![Vec::new()];
    for id in 1..=ndist {
        let iref1 = st.sym.ireff(st.sym.nsymboluptodata + 1, id);
        let (iref2, iref3, iref4) = (iref1 + 1, iref1 + 2, iref1 + 3);
        let atom1stext = st.sym.symbol_value[iref1].clone();
        let atom2stext = st.sym.symbol_value[iref2].clone();
        let rminstext = st.sym.symbol_value[iref3].clone();
        let rmaxstext = st.sym.symbol_value[iref4].clone();
        let w = parse(&atom1stext, 12);
        let atom1s: Vec<String> = (1..=w.nwords()).map(|iw| crate::fortran::fstr(w.word(&atom1stext, iw), 4)).collect();
        let w = parse(&atom2stext, 12);
        let mut atom2s = Vec::new();
        for iw in 1..=w.nwords() {
            let s = read_word(w.word(&atom2stext, iw)).map_err(|e| crate::error::EpsrError::Input(format!("atom2: {e:?}")))?;
            atom2s.push(crate::fortran::fstr(&s, 4));
        }
        let w = parse_n(&rminstext, 1, 256);
        let mut rmins = Vec::new();
        for iw in 1..=w.nwords() {
            rmins.push(read_real(w.word(&rminstext, iw)).map_err(|e| crate::error::EpsrError::Input(format!("rmin: {e:?}")))?);
        }
        let w = parse_n(&rmaxstext, 1, 256);
        let mut rmaxs = Vec::new();
        for iw in 1..=w.nwords() {
            rmaxs.push(read_real(w.word(&rmaxstext, iw)).map_err(|e| crate::error::EpsrError::Input(format!("rmax: {e:?}")))?);
        }
        let n = atom1s.len();
        if n != atom2s.len() || n != rmins.len() || n != rmaxs.len() {
            st.log(" Error: numbers of atoms and cutoffs in cluster bonds definitions not equal.");
            return Ok(());
        }
        if n > MAXNCBONDS {
            return crate::error::unsupported("local_density: more than 10 bond definitions in a distribution");
        }
        let mut d = Vec::new();
        let mut a = Vec::new();
        for ncbs in 0..n {
            let mut it = 0;
            let mut ifound1 = 0;
            let mut ifound2 = 0;
            while it < st.mol.ntype && (ifound1 == 0 || ifound2 == 0) {
                it += 1;
                if st.mol.type_[it] == atom1s[ncbs] {
                    ifound1 = it;
                }
                if st.mol.type_[it] == atom2s[ncbs] {
                    ifound2 = it;
                }
            }
            if ifound1 == 0 || ifound2 == 0 {
                return crate::error::unsupported("local_density: bond atom type not found in the .ato file (Fortran reads type(0))");
            }
            d.push(Def { ltyp: ifound1, ltyp1: ifound2, rmin: rmins[ncbs], rmax: rmaxs[ncbs] });
            a.push((st.mol.type_[ifound1].clone(), st.mol.type_[ifound2].clone()));
        }
        defs.push(d);
        atoms.push(a);
    }
    let mut z: F1<f32> = F1::new(MZ);
    for j in 1..=nz {
        z[j] = (j - 1) as f32;
    }
    let mut gofz: F2<f32> = F2::new(MZ, MDIST);
    st.normalise_coordinates();
    find_clusters(&st, nz, ndist, &defs, &mut gofz)?;

    let titler = "Cluster_size";
    let mut titledat = Vec::new();
    for id in 1..=ndist as usize {
        let mut t = vec![b' '; 29];
        let mut icharidx: usize = 1;
        for (a1, a2) in &atoms[id] {
            let s1 = upto_blank(a1);
            let s2 = upto_blank(a2);
            let nchars = s1.len() + s2.len() + 3;
            if icharidx + nchars > 29 {
                return crate::error::unsupported("local_density: column title longer than 29 characters (Fortran runtime error)");
            }
            let piece = crate::fortran::format::pad_to(&format!("{s1}-{s2},"), nchars + 1);
            t[icharidx - 1..icharidx + nchars].copy_from_slice(piece.as_bytes());
            icharidx += nchars;
        }
        if icharidx >= 3 {
            for c in t.iter_mut().skip(icharidx - 3) {
                *c = b' ';
            }
        } else {
            return crate::error::unsupported("local_density: distribution without bond definitions");
        }
        titledat.push(String::from_utf8(t).unwrap().trim_end().to_string());
    }
    let mut gzt: F2<f32> = F2::new(MZ, MDIST);
    let mut gzt2: F2<f32> = F2::new(MZ, MDIST);
    outgofz(&st, "n", nz, ndist, &mut nsumt, &mut z, &gofz, &mut gzt, &mut gzt2, titler, &titledat)?;

    let mut fract2 = vec![0.0f32; ndist as usize + 1];
    for id in 1..=ndist {
        let mut sum1 = 0.0f32;
        let mut sum2 = 0.0f32;
        for i in 1..=nz {
            sum1 += z[i] * gzt[(i, id)];
            if z[i] > 1.0 {
                sum2 += z[i] * gzt[(i, id)];
            }
        }
        fract2[id as usize] = sum2 / sum1;
    }
    let icount = st.sym.search_for_symbol("nsumt", 0);
    let mut iref = st.sym.ireff(icount, 1);
    set_sv_int(&mut st, iref, nsumt);
    iref += 1;
    set_sv_int(&mut st, iref, ndist);
    st.save_inp_file()?;
    let mut lines = Vec::new();
    for id in 1..=ndist as usize {
        lines.push(fmt("(/'Cluster ',i3,': fraction containing 2 or more molecules',f10.5)", &[FV::I(id as i64), FV::R(fract2[id])]));
    }
    write_fort14(ctx, &lines)
}

/// `subroutine find_clusters()` of `local_density.f`.
fn find_clusters(st: &State, nz: i32, ndist: i32, defs: &[Vec<Def>], gofz: &mut F2<f32>) -> Result<()> {
    let m = &st.mol;
    let nmol = m.nmol;
    let dmax = m.dmax;
    let nm = nmol as usize + 2;
    let mut mbond = vec![0i32; nm];
    let mut neighb: Vec<Vec<i32>> = vec![Vec::new(); nm];
    let mut labcluster = vec![0i32; nm];
    let mut lencluster = vec![0i32; nm];
    let mut ncllab = vec![0i32; nm];
    let damax = 0.0f32; // damax(id) of the common block is never set in local_density
    for id in 1..=ndist as usize {
        let d = &defs[id];
        for i in 1..=nmol as usize {
            mbond[i] = -1;
            neighb[i].clear();
        }
        for im1 in 1..=nmol - 1 {
            for ia1 in m.naf[im1]..=m.nal[im1] {
                for im2 in im1 + 1..=nmol {
                    for ia2 in m.naf[im2]..=m.nal[im2] {
                        for b in d {
                            if (m.itype[ia1] == b.ltyp && m.itype[ia2] == b.ltyp1) || (m.itype[ia2] == b.ltyp && m.itype[ia1] == b.ltyp1) {
                                let (u1, u2) = (im1 as usize, im2 as usize);
                                if mbond[u1] < 0 {
                                    mbond[u1] = 0;
                                }
                                if mbond[u2] < 0 {
                                    mbond[u2] = 0;
                                }
                                let xa1 = m.xa[ia1] + m.xadd[im1];
                                let ya1 = m.ya[ia1] + m.yadd[im1];
                                let za1 = m.za[ia1] + m.zadd[im1];
                                let xa2 = m.xa[ia2] + m.xadd[im2];
                                let ya2 = m.ya[ia2] + m.yadd[im2];
                                let za2 = m.za[ia2] + m.zadd[im2];
                                let mut xdif = xa2 - xa1;
                                let mut ydif = ya2 - ya1;
                                let mut zdif = za2 - za1;
                                minimagec(&mut xdif, &mut ydif, &mut zdif);
                                let dsq12 = xdif * xdif + ydif * ydif + zdif * zdif;
                                let da12 = dsq12.sqrt() * dmax;
                                if da12 >= b.rmin && da12 <= b.rmax {
                                    let bonded = neighb[u1][..mbond[u1] as usize].contains(&im2);
                                    if !bonded {
                                        if mbond[u1] as usize >= MDIST || mbond[u2] as usize >= MDIST {
                                            return crate::error::unsupported("local_density: a molecule bonded to more than 30 molecules (overflows neighb in the Fortran)");
                                        }
                                        mbond[u1] += 1;
                                        mbond[u2] += 1;
                                        neighb[u1].push(im2);
                                        neighb[u2].push(im1);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        let mut ncluster = 0;
        for im1 in 1..=nmol as usize {
            if mbond[im1] >= 0 {
                ncluster += 1;
                ncllab[im1] = ncluster;
                labcluster[ncluster as usize] = ncluster;
            } else {
                ncllab[im1] = nmol + 1;
            }
        }
        let mut nclusterold = 0;
        while nclusterold != ncluster {
            nclusterold = ncluster;
            for im1 in 1..=nmol as usize {
                if mbond[im1] >= 0 {
                    let mut nseek = ncllab[im1];
                    for k in 0..mbond[im1] as usize {
                        let im2 = neighb[im1][k] as usize;
                        if ncllab[im2] < nseek {
                            nseek = ncllab[im2];
                        }
                    }
                    ncllab[im1] = nseek;
                    for k in 0..mbond[im1] as usize {
                        let im2 = neighb[im1][k] as usize;
                        ncllab[im2] = nseek;
                    }
                }
            }
            ncluster = 0;
            for im1 in 1..=nmol as usize {
                if mbond[im1] >= 0 {
                    let mut iflag = 0;
                    let mut ic = 0;
                    while ic < ncluster && iflag == 0 {
                        ic += 1;
                        if ncllab[im1] == labcluster[ic as usize] {
                            iflag = ic;
                        }
                    }
                    if iflag == 0 {
                        ncluster += 1;
                        labcluster[ncluster as usize] = ncllab[im1];
                    }
                }
            }
            st.log(&format!(" {ncluster} {nclusterold}"));
        }
        for ic in 1..=ncluster as usize {
            let mut sum = 0.0f32;
            for im in 1..=nmol as usize {
                if ncllab[im] == labcluster[ic] {
                    sum += 1.0;
                }
            }
            lencluster[ic] = sum as i32;
        }
        let mut lensum = 0;
        for ic in 1..=ncluster as usize {
            lensum += lencluster[ic];
        }
        let avlen = if ncluster > 0 { lensum as f32 / ncluster as f32 } else { 0.0 };
        st.log(&fmt(
            "(/'For cluster ',i2,' there are ',i5,' clusters of average size = ',e12.5,' molecules')",
            &[FV::I(id as i64), FV::I(ncluster as i64), FV::R(avlen)],
        ));
        let mut npercolate = 0;
        for ic in 1..=ncluster as usize {
            let mut im1 = 0usize;
            let mut iflag = 0;
            while im1 < nmol as usize && iflag == 0 {
                im1 += 1;
                if ncllab[im1] == labcluster[ic] {
                    iflag = im1;
                }
            }
            let (mut xmin, mut ymin, mut zmin, mut xmax, mut ymax, mut zmax) = (0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32);
            for im2 in 1..=nmol as usize {
                if ncllab[im2] == labcluster[ic] {
                    let mut xdif = m.xadd[im2] - m.xadd[im1];
                    let mut ydif = m.yadd[im2] - m.yadd[im1];
                    let mut zdif = m.zadd[im2] - m.zadd[im1];
                    minimagec(&mut xdif, &mut ydif, &mut zdif);
                    if xdif > xmax {
                        xmax = xdif;
                    }
                    if xdif < xmin {
                        xmin = xdif;
                    }
                    if ydif > ymax {
                        ymax = ydif;
                    }
                    if ydif < ymin {
                        ymin = ydif;
                    }
                    if zdif > zmax {
                        zmax = zdif;
                    }
                    if zdif < zmin {
                        zmin = zdif;
                    }
                }
            }
            let xdif = xmax - xmin;
            let ydif = ymax - ymin;
            let zdif = zmax - zmin;
            let diftest = 1.0 - damax / dmax;
            if xdif >= diftest || ydif >= diftest || zdif >= diftest {
                npercolate += lencluster[ic];
            }
        }
        let fracpercolate = npercolate as f32 / lensum as f32;
        st.log(&fmt(
            "(/'Clusters> For cluster ',i2,' fraction of molecules which percolate is ',f10.5)",
            &[FV::I(id as i64), FV::R(fracpercolate)],
        ));
        let idi = id as i32;
        for ic in 1..=ncluster as usize {
            let iref = lencluster[ic] + 1;
            if iref <= nz {
                gofz[(iref, idi)] += 1.0;
            } else {
                gofz[(nz, idi)] += 1.0;
            }
        }
        for ic in 2..=nz {
            gofz[(ic, idi)] /= ncluster as f32;
        }
        gofz[(1, idi)] = fracpercolate;
    }
    Ok(())
}
