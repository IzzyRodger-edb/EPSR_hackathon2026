//! `rings` — translation of `fortran_src/rings_sp_criterion.f` (`run_rings`,
//! `find_rings`), using `match`/`shortest_path` of `chains_routines.f`
//! (in `chains.rs`).
//!
//! STATUS: translated, verified byte-identical (EPSR-test-ex64, ex65).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> rings <input file>`
//! Keyword table: `crate::symtab_tools::RINGS` (the `nsize` text and default
//! are computed from `mring` as in `setup_symbols_rings`).
//!
//! `run_rings` is identical to `run_chains` apart from the limit on `nsize`
//! (`mring`), the axis title and the search routine, so the driver is shared
//! (`chains::run_chains_or_rings`).  Outputs: `<run>.RINGS.n01` (histograms of
//! bonded neighbours and of ring sizes), the rewritten `<run>.RINGS.dat`, and
//! `fort.14` (means and spreads).

use super::chains::{find_bonds, shortest_path, Crit, Neigh, SpWork, MBOND, MMOL, MRING};
use super::ToolContext;
use crate::error::Result;
use crate::fortran::{F1, F2};
use crate::state::State;

/// Run the `rings` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    super::chains::run_chains_or_rings(ctx, true)
}

/// `subroutine find_rings()` (shortest-path ring criterion of Franzblau).
pub(crate) fn find_rings(st: &State, c: &Crit, nz: i32, ndist: i32, gofz: &mut F2<f32>) {
    let nmol = st.mol.nmol;
    let mut ng = Neigh::new(nmol, MBOND);
    let nmu = nmol as usize + 1;
    let mut spvalue = F1::<i32>::new(nmu);
    let mut carryon = F1::<bool>::new(nmu);
    let mut lenring = F1::<i32>::new(MMOL as usize);
    let mz = 200usize;
    let mut nrlab = vec![0i32; mz + 2];
    let mut rneighb = vec![0i32; mz + 2];
    let mut spseq = vec![0i32; 202];
    let _ = MRING;
    let mut idref2 = 0;
    for id in 1..=ndist {
        let idref1 = idref2 + 1;
        idref2 = idref1 + 1;
        let _nfsum = find_bonds(st, c, id, &mut ng);
        let nspmax = nz / 2 + 1;
        let mut w = SpWork::new(nmol, nspmax);
        let mut nring = 0;
        let mut im1 = 0;
        // common integers of chains.inc
        let mut spmin: i32;
        let mut spmax: i32;
        while im1 < nmol && nring < MMOL {
            im1 += 1;
            if ng.nneighb[im1] > 1 {
                for im in 1..=nmol {
                    carryon[im] = im >= im1;
                }
                shortest_path(im1, nspmax, nmol, &ng, &mut w);
                for im in 1..=nmol {
                    spvalue[im] = w.spvaluet[im];
                }
                nrlab[0] = im1;
                nrlab[1] = im1;
                spseq[0] = -1;
                spseq[1] = 0;
                rneighb[1] = 0;
                spmin = spvalue[ng.neighb[(im1, 1)]];
                for itest in 2..=ng.nneighb[im1] {
                    let v = spvalue[ng.neighb[(im1, itest)]];
                    if v < spmin {
                        spmin = v;
                    }
                }
                let mut icount: usize = 1;
                let mut ibase = 1;
                while icount > 0 && nring < MMOL {
                    let isub = if icount == 1 { 1 } else { 0 };
                    if rneighb[icount] >= ng.nneighb[nrlab[icount]] - isub {
                        if icount == 2 {
                            carryon[nrlab[icount]] = false;
                        }
                        icount -= 1;
                        if ibase > 1 {
                            ibase -= 1;
                        }
                        if icount > 1 {
                            spmin = nmol + 1;
                            let nr = nrlab[icount];
                            for itest in 1..=ng.nneighb[nr] {
                                let nb = ng.neighb[(nr, itest)];
                                if spvalue[nb] < spmin && nb != nrlab[icount - 1] {
                                    spmin = spvalue[nb];
                                }
                            }
                        } else {
                            spmin = 1;
                        }
                    } else {
                        rneighb[icount] += 1;
                        let im2 = ng.neighb[(nrlab[icount], rneighb[icount])];
                        if carryon[im2]
                            && spvalue[im2] == spmin
                            && spvalue[im2] < nmol
                            && ng.nneighb[im2] > 1
                            && im2 != nrlab[icount - 1]
                            && (icount as i32) < nz
                            && nring < MMOL
                        {
                            let mut ifound = 0;
                            let mut ic = 1;
                            while ifound == 0 && ic < icount {
                                ic += 1;
                                if im2 == nrlab[ic] {
                                    ifound = im2;
                                }
                            }
                            if ifound == 0 {
                                if im2 == im1 {
                                    let mut found = false;
                                    let ic_n = icount as i32;
                                    let mut ic1: i32 = 1;
                                    while ic1 < ic_n - 2 && !found {
                                        ic1 += 1;
                                        let mut ic2 = ic1 + 1;
                                        shortest_path(nrlab[ic1 as usize], nspmax, nmol, &ng, &mut w);
                                        while ic2 < ic_n && !found {
                                            ic2 += 1;
                                            spmin = (ic2 - ic1).min(ic_n + ic1 - ic2);
                                            spmax = (ic2 - ic1).max(ic_n + ic1 - ic2);
                                            let r2 = nrlab[ic2 as usize];
                                            found = spmin > w.spvaluet[r2];
                                            if !found && spmax > spmin {
                                                let mut icc = spmin;
                                                found = w.nspvaluet[(r2, icc)] > 1;
                                                while !found && icc < spmax - 1 && icc < nspmax {
                                                    icc += 1;
                                                    found = w.nspvaluet[(r2, icc)] > 0;
                                                }
                                            }
                                        }
                                    }
                                    if !found && nring < MMOL {
                                        nring += 1;
                                        lenring[nring] = icount as i32;
                                    }
                                } else {
                                    icount += 1;
                                    nrlab[icount] = im2;
                                    rneighb[icount] = 0;
                                    spseq[icount] = spvalue[im2];
                                    if icount > 1 {
                                        spmin = nmol + 1;
                                        for itest in 1..=ng.nneighb[im2] {
                                            let nb = ng.neighb[(im2, itest)];
                                            if spvalue[nb] < spmin && nb != nrlab[icount - 1] {
                                                spmin = spvalue[nb];
                                            }
                                        }
                                    } else {
                                        spmin = 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            carryon[im1] = false;
        }
        for im1 in 1..=nmol {
            let iref = ng.nneighb[im1] + 1;
            if iref <= nz {
                gofz[(iref, idref1)] += 1.0;
            }
        }
        for ir in 1..=nring {
            let iref = lenring[ir] + 1;
            if iref <= nz {
                gofz[(iref, idref2)] += 1.0;
            }
        }
    }
}
