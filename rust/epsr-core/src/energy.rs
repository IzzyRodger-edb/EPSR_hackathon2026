//! Potential energy evaluation.
//!
//! Translation of:
//! * `en_atom_inter_scale`, `en_mole_pair`, `en_atom_intra`, `en_mole_intra`,
//!   `en_mole_dihedral`, `en_change`, `en_change_intra`, `check_tethered`
//!                                   (`molfit.f`)
//! * `fmolec1`, `press`, `abs_energy_ep`, `adjust_intra_wts`
//!                                   (`ato_manipulation_rev4.f`)
//!
//! Energies are looked up from tabulated potentials on the r grid
//! (`xstar` in box-fraction units, spacing `cellst/dmax`) with linear
//! interpolation using the stored gradients:
//! * intermolecular: empirical potential `xavs`/`xavs1` + reference potential `ecoul`/`ecoul1`;
//! * intramolecular non-bonded: `xav`/`xav1`;
//! * bonded (bond / angle restraints): harmonic springs `amassr/tb * (tb - r)^2 * tempb * kT`.
//!
//! Three-body terms, coordination constraints, container forces, tethered
//! molecules and dihedral (`DIH`) groups are not part of the port; the input
//! readers reject files that use them, so the corresponding branches are
//! omitted here.

use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::F1;
use crate::neighbour::minimagec;
use crate::state::State;
use rayon::prelude::*;

/// `function en_atom_inter_scale(x1,y1,z1,x2,y2,z2,stepcl,xstar,icref,npt2,scalefac)`.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub fn en_atom_inter_scale(
    st: &State,
    x1: f32,
    y1: f32,
    z1: f32,
    x2: f32,
    y2: f32,
    z2: f32,
    stepcl: f32,
    icref: i32,
    npt2: i32,
    scalefac: f32,
) -> f32 {
    let p = &st.pot;
    let xstar = &st.dist.xstar;
    let mut sume: f32 = 0.0;
    let mut xdif = x2 - x1;
    let mut ydif = y2 - y1;
    let mut zdif = z2 - z1;
    minimagec(&mut xdif, &mut ydif, &mut zdif);
    let dsq = xdif * xdif + ydif * ydif + zdif * zdif;
    let mut dtest = dsq.sqrt();
    let rat = dtest / stepcl;
    let iz = rat as i32 + 2;
    if iz <= npt2 {
        let dtdif = dtest - xstar[iz];
        sume = p.xavs[(iz, icref)] + p.xavs1[(iz, icref)] * dtdif;
    }
    dtest *= scalefac;
    let rat = dtest / stepcl;
    let iz = rat as i32 + 2;
    if iz <= npt2 {
        let dtdif = dtest - xstar[iz];
        sume = sume + p.ecoul[(iz, icref)] + p.ecoul1[(iz, icref)] * dtdif;
    }
    sume
}

/// `function en_dihedral(co1,co2,cd,dihedralang)`: (cos φ − cos φ₀)² + (sin φ − sin φ₀)².
fn en_dihedral(co1: &[f32; 3], co2: &[f32; 3], cd: &[f32; 3], dihedralang: f32) -> f32 {
    use crate::ato::{scapro, vecpro};
    let mut cv1 = vecpro(co1, cd);
    let mut cv2 = vecpro(co2, cd);
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
    let sinphi = scapro(&cpar, cd);
    let cosdih = dihedralang.cos();
    let sindih = dihedralang.sin();
    let difc = cosdih - cosphi;
    let difs = sindih - sinphi;
    difc * difc + difs * difs
}

impl State {
    /// `function en_mole_pair(im, jm)`: intermolecular energy of a molecule pair.
    pub fn en_mole_pair(&self, im: i32, jm: i32) -> f32 {
        let m = &self.mol;
        let stepcl = self.dist.cellst / m.dmax;
        let npt2 = self.dist.npt2;
        let mut scalefac: f32 = 1.0;
        if m.coarsegraining {
            scalefac = 0.5 * (m.scalefacmol[im] + m.scalefacmol[jm]);
        }
        let bothtethered = m.tethered[im] && m.tethered[jm];
        let mut sume: f32 = 0.0;
        for ia in m.naf[im]..=m.nal[im] {
            let it = m.itype[ia];
            let x1 = m.xa[ia] + m.xadd[im];
            let y1 = m.ya[ia] + m.yadd[im];
            let z1 = m.za[ia] + m.zadd[im];
            let qtypeia = m.qtype[ia] == 1 && bothtethered;
            for ja in m.naf[jm]..=m.nal[jm] {
                if m.qtype[ja] == 1 && qtypeia {
                    continue;
                }
                let jt = m.itype[ja];
                let icref = m.npairref[(it, jt)];
                let x2 = m.xa[ja] + m.xadd[jm];
                let y2 = m.ya[ja] + m.yadd[jm];
                let z2 = m.za[ja] + m.zadd[jm];
                if self.cc.ncoordconstraints > 0 && self.cc.check_coord_bonded(ia, ja, m) {
                    sume += crate::coordconstraint::en_atom_inter_scale_coordextra(self, x1, y1, z1, x2, y2, z2, stepcl, icref, npt2, scalefac);
                } else {
                    sume += en_atom_inter_scale(self, x1, y1, z1, x2, y2, z2, stepcl, icref, npt2, scalefac);
                }
            }
        }
        sume
    }

    /// `function en_atom_intra(ia,ja,x1,y1,z1,x2,y2,z2,sumevib,sumeang)`:
    /// energy between two atoms of the same molecule.
    #[allow(clippy::too_many_arguments)]
    pub fn en_atom_intra(&self, ia: i32, ja: i32, x1: f32, y1: f32, z1: f32, x2: f32, y2: f32, z2: f32, sumevib: &mut f32, sumeang: &mut f32) -> f32 {
        let m = &self.mol;
        if ja == m.imageatom[ia] {
            return 0.0;
        }
        let mut sume: f32 = 0.0;
        let stepcl = self.dist.cellst / m.dmaxmol;
        let it = m.itype[ia];
        let jt = m.itype[ja];
        let icref = m.npairref[(it, jt)];
        let xdifs = x2 - x1;
        let ydifs = y2 - y1;
        let zdifs = z2 - z1;
        let mut xdif = xdifs;
        let mut ydif = ydifs;
        let mut zdif = zdifs;
        minimagec(&mut xdif, &mut ydif, &mut zdif);
        let mut jref = 0;
        let mut iflag = 0;
        let mut image = ia;
        while jref < m.nbond[image] && iflag == 0 {
            jref += 1;
            if ja == m.iabond[(image, jref)] {
                iflag = jref;
            }
        }
        if iflag == 0 {
            // bonds may be defined on the image atom of ia (molecules spanning the box)
            image = 0;
            if m.imageatom[ia] > 0 {
                image = m.imageatom[ia];
                jref = 0;
                while jref < m.nbond[image] && iflag == 0 {
                    jref += 1;
                    if ja == m.iabond[(image, jref)] {
                        iflag = jref;
                    }
                }
            }
        }
        if iflag > 0 && image > 0 {
            let dsq = xdif * xdif + ydif * ydif + zdif * zdif;
            let dtest = dsq.sqrt();
            let wtms = m.amassr[icref];
            let tbv = m.tbond[(image, jref)];
            let tb = stepcl.max(tbv.abs());
            let mut wtm = wtms / tb;
            if tbv < 0.0 {
                wtm *= m.tempang;
            }
            let ddif = tb - dtest;
            sume = wtm * ddif * ddif * m.tempb * m.simt;
            if tbv > 0.0 {
                *sumevib += sume;
            }
            if tbv < 0.0 {
                *sumeang += sume;
            }
            // refine_intra: the EP also acts on bonded pairs
            let iz = (dtest / stepcl) as i32 + 2;
            if m.refine_intra && iz <= self.dist.npt2 && m.qtype[ia] == 0 && m.qtype[ja] == 0 {
                let dtdif = dtest - self.dist.xstar[iz];
                sume += m.refine_intra_fac * (self.pot.xavs[(iz, icref)] + self.pot.xavs1[(iz, icref)] * dtdif);
            }
        } else {
            // atom pairs forming a dihedral of a DIH group are handled by en_mole_dihedral
            let im = m.mref[ia];
            let mut found = false;
            let mut ir = m.nrf[im] - 1;
            while ir < m.nrl[im] && !found {
                ir += 1;
                let mut kd = 0;
                while kd < m.ndihedralang[ir] && !found {
                    kd += 1;
                    let io1 = m.labrotdih[(1, kd, ir)];
                    let io2 = m.labrotdih[(2, kd, ir)];
                    found = (io1 == ia && io2 == ja) || (io1 == ja && io2 == ia);
                }
            }
            if found {
                return 0.0;
            }
            let xdif0 = xdifs * m.sizefactor;
            let ydif0 = ydifs * m.sizefactor;
            let zdif0 = zdifs * m.sizefactor;
            let (mut xdif, mut ydif, mut zdif) = (xdif0, ydif0, zdif0);
            minimagec(&mut xdif, &mut ydif, &mut zdif);
            let dsq = xdif * xdif + ydif * ydif + zdif * zdif;
            let dtest = dsq.sqrt();
            let rat = dtest / stepcl;
            let iz = rat as i32 + 2;
            if iz <= self.dist.npt2 {
                let dtdif = dtest - self.dist.xstar[iz];
                sume = self.pot.xav[(iz, icref)] + self.pot.xav1[(iz, icref)] * dtdif;
                if m.refine_intra && m.qtype[ia] == 0 && m.qtype[ja] == 0 {
                    sume += m.refine_intra_fac * (self.pot.xavs[(iz, icref)] + self.pot.xavs1[(iz, icref)] * dtdif);
                }
            }
        }
        sume
    }

    /// `function en_mole_intra(im, sumevib, sumeang, sumedih)`: total
    /// intramolecular energy of molecule `im`.
    pub fn en_mole_intra(&self, im: i32, sumevib: &mut f32, sumeang: &mut f32, sumedih: &mut f32) -> f32 {
        let m = &self.mol;
        let mut e: f32 = 0.0;
        for ia in m.naf[im]..m.nal[im] {
            for ja in ia + 1..=m.nal[im] {
                let sume = self.en_atom_intra(ia, ja, m.xa[ia], m.ya[ia], m.za[ia], m.xa[ja], m.ya[ja], m.za[ja], sumevib, sumeang);
                e += sume;
            }
        }
        if m.ngrot[im] > 0 && m.tempdih > 0.0 {
            let sume = m.tempdih * m.simt * self.en_mole_dihedral(im);
            e += sume;
            *sumedih += sume;
        }
        e
    }

    /// `function en_mole_dihedral(im)`: energy of the dihedral restraints of
    /// the `DIH` rotational groups of molecule `im`.
    pub fn en_mole_dihedral(&self, im: i32) -> f32 {
        let m = &self.mol;
        let mut e: f32 = 0.0;
        for ir in m.nrf[im]..=m.nrl[im] {
            if m.rotgrp[ir] == "DIH" {
                let ia1 = m.labaxi[(ir, 1)];
                let ia2 = m.labaxi[(ir, 2)];
                let (x1, y1, z1) = (m.xa[ia1], m.ya[ia1], m.za[ia1]);
                let (x2, y2, z2) = (m.xa[ia2], m.ya[ia2], m.za[ia2]);
                let mut cd = [x2 - x1, y2 - y1, z2 - z1];
                let dmag = crate::ato::scapro(&cd, &cd).sqrt();
                for v in cd.iter_mut() {
                    *v /= dmag;
                }
                let mut kd = 0;
                while kd < m.ndihedralang[ir] {
                    kd += 1;
                    let io1 = m.labrotdih[(1, kd, ir)];
                    let io2 = m.labrotdih[(2, kd, ir)];
                    let co1 = [m.xa[io1] - x1, m.ya[io1] - y1, m.za[io1] - z1];
                    let co2 = [m.xa[io2] - x2, m.ya[io2] - y2, m.za[io2] - z2];
                    e += en_dihedral(&co1, &co2, &cd, m.dihedralang[(kd, ir)]);
                }
            }
        }
        e
    }

    /// `function en_tethered(im, ia)`: squared distance of the tethered atom
    /// (or molecule centre) from its tethering point.
    pub fn en_tethered(&self, im: i32, ia: i32) -> f32 {
        let m = &self.mol;
        let mut xdif = m.xadd[im] - m.phix[im];
        let mut ydif = m.yadd[im] - m.phiy[im];
        let mut zdif = m.zadd[im] - m.phiz[im];
        if ia > 0 && ia <= m.nato[im] {
            let iaref = m.naf[im] + ia - 1;
            xdif += m.xa[iaref];
            ydif += m.ya[iaref];
            zdif += m.za[iaref];
        }
        minimagec(&mut xdif, &mut ydif, &mut zdif);
        xdif * xdif + ydif * ydif + zdif * zdif
    }

    /// `function en_change_tethered(im, ia)`: change of `en_tethered` relative
    /// to the saved centre (`xs`, `ys`, `zs`).
    pub fn en_change_tethered(&self, im: i32, ia: i32) -> f32 {
        let m = &self.mol;
        let mut xdif = m.xs - m.phix[im];
        let mut ydif = m.ys - m.phiy[im];
        let mut zdif = m.zs - m.phiz[im];
        if ia > 0 && ia <= m.nato[im] {
            let iaref = m.naf[im] + ia - 1;
            xdif += m.xa[iaref];
            ydif += m.ya[iaref];
            zdif += m.za[iaref];
        }
        minimagec(&mut xdif, &mut ydif, &mut zdif);
        let term = xdif * xdif + ydif * ydif + zdif * zdif;
        let mut xdif = m.xadd[im] - m.phix[im];
        let mut ydif = m.yadd[im] - m.phiy[im];
        let mut zdif = m.zadd[im] - m.phiz[im];
        if ia > 0 && ia <= m.nato[im] {
            let iaref = m.naf[im] + ia - 1;
            xdif += m.xa[iaref];
            ydif += m.ya[iaref];
            zdif += m.za[iaref];
        }
        minimagec(&mut xdif, &mut ydif, &mut zdif);
        xdif * xdif + ydif * ydif + zdif * zdif - term
    }

    /// `function check_tethered(im, ia)`: may a tethered molecule move here?
    /// (Always true for untethered molecules.)
    pub fn check_tethered(&self, im: i32, ia: i32) -> bool {
        let m = &self.mol;
        if !m.tethered[im] {
            return true;
        }
        let mut xdif = m.xadd[im] - m.phix[im];
        let mut ydif = m.yadd[im] - m.phiy[im];
        let mut zdif = m.zadd[im] - m.phiz[im];
        if ia > 0 && ia <= m.nato[im] {
            let iaref = m.naf[im] + ia - 1;
            xdif += m.xa[iaref];
            ydif += m.ya[iaref];
            zdif += m.za[iaref];
        }
        minimagec(&mut xdif, &mut ydif, &mut zdif);
        let dist = (xdif * xdif + ydif * ydif + zdif * zdif).sqrt();
        let ok = dist <= m.tol;
        if !ok {
            return self.en_change_tethered(im, ia) <= 0.0;
        }
        ok
    }

    /// `function en_change_intra(im, ...)`: change of the intramolecular
    /// energy of molecule `im` relative to the stored value (`ergold(1)`).
    pub fn en_change_intra(&mut self, im: i32, sumevib: &mut f32, sumeang: &mut f32, sumedih: &mut f32) -> f32 {
        let t0 = std::time::Instant::now();
        *sumevib = 0.0;
        *sumeang = 0.0;
        *sumedih = 0.0;
        let e = self.en_mole_intra(im, sumevib, sumeang, sumedih);
        self.erg.ergnew[1] = e;
        let r = self.erg.ergnew[1] - self.erg.ergold[1];
        self.timer.sumtime2 += t0.elapsed().as_millis() as i64 / 1000;
        r
    }

    /// `function en_change(im)`: change in intermolecular energy of molecule
    /// `im` with its stored neighbours (`nbref(2..nnbref)`), filling `ergnew`.
    pub fn en_change(&mut self, im: i32) -> f32 {
        let t0 = std::time::Instant::now();
        let stepcl = self.dist.cellst / self.mol.dmax;
        let nnbref = self.erg.nnbref;
        let npt2 = self.dist.npt2;
        for inn in 2..=nnbref {
            self.erg.ergnew[inn] = 0.0;
        }
        for inn in 2..=nnbref {
            let jm = self.erg.nbref[inn];
            let m = &self.mol;
            // Scale factor on intermolecular spacings
            let mut scalefac: f32 = 1.0;
            if m.coarsegraining {
                scalefac = 0.5 * (m.scalefacmol[im] + m.scalefacmol[jm]);
            }
            let bothtethered = m.tethered[im] && m.tethered[jm];
            let mut sume: f32 = 0.0;
            for ia in m.naf[im]..=m.nal[im] {
                let it = m.itype[ia];
                let x1 = m.xa[ia] + m.xadd[im];
                let y1 = m.ya[ia] + m.yadd[im];
                let z1 = m.za[ia] + m.zadd[im];
                let qtypeia = m.qtype[ia] == 1 && bothtethered;
                for ja in m.naf[jm]..=m.nal[jm] {
                    if m.qtype[ja] == 1 && qtypeia {
                        continue;
                    }
                    let jt = m.itype[ja];
                    let icref = m.npairref[(it, jt)];
                    let x2 = m.xa[ja] + m.xadd[jm];
                    let y2 = m.ya[ja] + m.yadd[jm];
                    let z2 = m.za[ja] + m.zadd[jm];
                    if self.cc.ncoordconstraints > 0 && self.cc.check_coord_bonded(ia, ja, m) {
                        sume += crate::coordconstraint::en_atom_inter_scale_coordextra(self, x1, y1, z1, x2, y2, z2, stepcl, icref, npt2, scalefac);
                    } else {
                        sume += en_atom_inter_scale(self, x1, y1, z1, x2, y2, z2, stepcl, icref, npt2, scalefac);
                    }
                }
            }
            self.erg.ergnew[inn] += sume;
        }
        let mut fdiff: f32 = 0.0;
        for inn in 2..=nnbref {
            fdiff = fdiff + self.erg.ergnew[inn] - self.erg.ergold[inn];
        }
        if self.tb.nthreeb > 0 {
            self.en_change_triplets(im, &mut fdiff);
        }
        if self.cont.numcontainerforces > 0 {
            let e = self.en_containerforce(im);
            self.cont.ergnew_container[im] = e;
            fdiff = fdiff + e - self.cont.ergold_container[im];
        }
        if self.cc.ncoordconstraints > 0 {
            let e = self.cc.en_mole_coord(im, &self.mol);
            self.cc.ergnew_coord[im] = e;
            fdiff = fdiff + e - self.cc.ergold_coord[im];
        }
        self.timer.sumtime2 += t0.elapsed().as_millis() as i64 / 1000;
        fdiff
    }

    /// `subroutine fmolec1(ierr, num_threads, ftest, ftestb, ftethered, fcoordconstraints)`:
    /// compute and store all intra- and intermolecular energies.  Returns
    /// `(ftest, ftestb)` = (intermolecular, intramolecular) energy per molecule.
    pub fn fmolec1(&mut self) -> (f32, f32) {
        self.initialise_ergpair();
        // Setup the triplets of atoms if three-body terms are requested
        if self.tb.nthreeb > 0 {
            self.log(&format!("fmolec1> No. of three-body terms = {:>3} {:>3}", self.tb.nthreeb, self.tb.nthreebangsmax));
            self.find_triplets();
        }
        if let Err(e) = self.setup_coord_constraints() {
            panic!("{e}");
        }
        let mut ftest: f32 = 0.0;
        let mut ftestb: f32 = 0.0;
        {
            let m = &mut self.mol;
            m.tempangwt = m.tempangreq;
            m.tempdihwt = m.tempdihreq;
            m.tempang = m.tempangwt;
            m.tempdih = m.tempdihwt * m.tempb;
        }
        let (mut sv, mut sa, mut sd) = (0f32, 0f32, 0f32);
        for im in 1..=self.mol.nmol {
            let energy = self.en_mole_intra(im, &mut sv, &mut sa, &mut sd);
            ftestb += energy;
            self.erg.set_ergpair(im, im, energy);
        }
        let mut sumetriplets: f32 = 0.0;
        if self.tb.nthreeb > 0 {
            for id in 1..=self.tb.nthreeb {
                self.tb.sume_triplet_dist[id] = 0.0;
            }
        }
        let mut sumecoordcon: f32 = 0.0;
        // Precompute the pair energies on the rayon pool: en_mole_pair is a
        // pure function of the (fixed at this point) positions, and each
        // (cell, molecule, neighbour) slot is written independently.  The
        // stored neighbour order, the set_ergpair insertion order and the
        // serial ftest sum below are unchanged, so the result is identical
        // at any thread count.
        let st: &State = self;
        let pair_energies: Vec<Vec<Vec<f32>>> = (1..st.cells.ncell + 1)
            .into_par_iter()
            .map(|ic| {
                st.sync_lattice();
                let mut out: Vec<Vec<f32>> = Vec::new();
                if st.cells.nmolsincell[ic] > 0 {
                    let mut nbref: F1<i32> = F1::new(st.erg.nbref.len());
                    let nnbref = st.cells.make_neighbour_list_from_cells(ic, &mut nbref);
                    out.reserve(st.cells.nmolsincell[ic] as usize);
                    for iref in 1..=st.cells.nmolsincell[ic] {
                        let im = st.cells.molsincell[(iref, ic)];
                        let mut e = vec![0.0f32; (nnbref + 1) as usize];
                        for inn in 1..=nnbref {
                            let jm = nbref[inn];
                            if jm > im {
                                e[inn as usize] = st.en_mole_pair(im, jm);
                            }
                        }
                        out.push(e);
                    }
                }
                out
            })
            .collect();
        let mut nbref = std::mem::take(&mut self.erg.nbref);
        let mut ergold = std::mem::take(&mut self.erg.ergold);
        for ic in 1..=self.cells.ncell {
            if self.cells.nmolsincell[ic] > 0 {
                let nnbref = self.cells.make_neighbour_list_from_cells(ic, &mut nbref);
                for iref in 1..=self.cells.nmolsincell[ic] {
                    let im = self.cells.molsincell[(iref, ic)];
                    let e = &pair_energies[(ic - 1) as usize][(iref - 1) as usize];
                    for inn in 1..=nnbref {
                        ergold[inn] = e[inn as usize];
                    }
                    // NB: the Fortran loop starts at inn=2 (the first candidate
                    // neighbour is never stored).  Reproduced as is.
                    let mut inn = 1;
                    while inn < nnbref {
                        inn += 1;
                        let jm = nbref[inn];
                        if ergold[inn] != 0.0 && jm > im {
                            self.erg.set_ergpair(im, jm, ergold[inn]);
                            ftest += ergold[inn];
                        }
                    }
                    // Get the triplet energy for this molecule
                    if self.tb.nthreeb > 0 {
                        let e = self.en_mole_triplets(im);
                        self.tb.ergold_triplet[im] = e;
                        sumetriplets += self.tb.ergold_triplet[im];
                        for id in 1..=self.tb.nthreeb {
                            self.tb.sume_triplet_dist[id] += self.tb.en_mole_triplet_dist[id];
                        }
                        ftest += self.tb.ergold_triplet[im];
                    }
                    if self.cc.ncoordconstraints > 0 {
                        let e = self.cc.en_mole_coord(im, &self.mol);
                        self.cc.ergold_coord[im] = e;
                        ftest += e;
                        sumecoordcon += e;
                    }
                    if self.cont.numcontainerforces > 0 {
                        let e = self.en_containerforce(im);
                        self.cont.ergold_container[im] = e;
                        ftest += e;
                    }
                }
                self.erg.nnbref = nnbref;
            }
        }
        self.erg.nbref = nbref;
        self.erg.ergold = ergold;
        let nmol = self.mol.nmol as f32;
        self.cc.fcoordconstraints = sumecoordcon / nmol;
        self.log(&format!(" {} {}", ftest / nmol, ftestb / nmol));
        let mut ftethered: f32 = 0.0;
        for im in 1..=self.mol.nmol {
            if self.mol.tethered[im] {
                ftethered += self.en_tethered(im, self.mol.tetheredatom[im]);
            }
        }
        self.ftethered = ftethered / nmol;
        let neighbsum = self.erg.get_neighbsum();
        self.log(&format!("\nfmolec1> Average no. of neighbours per molecule = {neighbsum}"));
        self.log(&fmt("('fmolec1> Average triplet energy per molecule = ',e12.5)", &[FV::R(sumetriplets / nmol)]));
        if self.tb.nthreeb > 0 {
            let mut w = ListWriter::new();
            w.chars("fmolec1> Total energy for each triplet: ");
            for id in 1..=self.tb.nthreeb {
                w.real(self.tb.sume_triplet_dist[id] / nmol);
            }
            let s = w.finish();
            self.log(&s);
        }
        if self.cc.ncoordconstraints > 0 {
            self.log(&fmt(
                "('fmolec1> Average coordination constraint energy per molecule = ',e12.5/)",
                &[FV::R(self.cc.fcoordconstraints)],
            ));
        }
        ftest /= nmol;
        ftestb /= nmol;
        (ftest, ftestb)
    }

    /// `subroutine adjust_intra_wts()`: recompute the bond-angle and dihedral
    /// weights (in this version simply reset to the requested values).
    pub fn adjust_intra_wts(&mut self) {
        let (mut sumevib, mut sumeang, mut sumedih) = (0f32, 0f32, 0f32);
        self.mol.tempang = 1.0;
        self.mol.tempdih = 10.0;
        for im in 1..=self.mol.nmol {
            let _ = self.en_mole_intra(im, &mut sumevib, &mut sumeang, &mut sumedih);
        }
        self.log(&format!(
            " adjust_intra_wts> Vibrational, bond angle and dihedral energies: {sumevib} {sumeang} {sumedih}"
        ));
        let m = &mut self.mol;
        m.tempangwt = m.tempangreq;
        m.tempang = m.tempangwt;
        m.tempdihwt = m.tempdihreq;
        m.tempdih = m.tempdihwt * m.tempb;
    }

    /// `real function abs_energy_ep()`: largest (max - min) range of the EP
    /// coefficients over all pairs — the "absolute energy" of the EP.
    pub fn abs_energy_ep(&self) -> f32 {
        let p = &self.pot;
        let mut energabs: f32 = 0.0;
        for ip in 1..=p.npairspcof {
            let mut mn: f32 = 0.0;
            let mut mx: f32 = 0.0;
            for i in 2..=p.ncoeffp {
                let term = p.pcoefp[(i, ip)];
                mn = mn.min(term);
                mx = mx.max(term);
            }
            energabs = energabs.max(mx - mn);
        }
        energabs
    }

    /// `subroutine press(energy, pressc, energf, pressf, energabs, energthreeb, prthreebout)`:
    /// intermolecular energy and pressure from the stored neighbour lists.
    /// Returns `(energy, pressc, energabs)`.
    pub fn press(&mut self) -> (f32, f32, f32) {
        let stepcl = self.dist.cellst / self.mol.dmax;
        let mut sumc: f32 = 0.0;
        let mut presst: f64 = 0.0;
        let energabs = self.abs_energy_ep();
        let npt1 = self.dist.npt1;
        // Each molecule's energy/virial contributions are a pure function of
        // the fixed positions and the stored neighbour structure, so they are
        // precomputed on the rayon pool into per-molecule buffers; the serial
        // pass below then loads them into the shared erg scratch and does the
        // sumc/presst reductions in the original element order, so the result
        // (and the final erg state) is identical at any thread count.
        let st: &State = self;
        let per_im: Vec<(Vec<f32>, Vec<f32>)> = (1..st.mol.nmol)
            .into_par_iter()
            .map(|im| {
                st.sync_lattice();
                let m = &st.mol;
                let p = &st.pot;
                let e = &st.erg;
                let xstar = &st.dist.xstar;
                let nnbref = e.nergpairm[im];
                let mut ergold = vec![0.0f32; (nnbref + 1).max(1) as usize];
                let mut ergnew = vec![0.0f32; (nnbref + 1).max(1) as usize];
                let xc1 = m.xadd[im];
                let yc1 = m.yadd[im];
                let zc1 = m.zadd[im];
                for jref in 1..=nnbref {
                    let jm = e.molergpair[(jref, im)];
                    if jm > im {
                        for ia in m.naf[im]..=m.nal[im] {
                            let it = m.itype[ia];
                            let x1 = m.xa[ia] + xc1;
                            let y1 = m.ya[ia] + yc1;
                            let z1 = m.za[ia] + zc1;
                            for ja in m.naf[jm]..=m.nal[jm] {
                                let xc2 = m.xadd[jm];
                                let yc2 = m.yadd[jm];
                                let zc2 = m.zadd[jm];
                                let jt = m.itype[ja];
                                let icref = m.npairref[(it, jt)];
                                let x2 = m.xa[ja] + xc2;
                                let y2 = m.ya[ja] + yc2;
                                let z2 = m.za[ja] + zc2;
                                let mut xdif = x2 - x1;
                                let mut ydif = y2 - y1;
                                let mut zdif = z2 - z1;
                                minimagec(&mut xdif, &mut ydif, &mut zdif);
                                let mut xcdif = xc2 - xc1;
                                let mut ycdif = yc2 - yc1;
                                let mut zcdif = zc2 - zc1;
                                minimagec(&mut xcdif, &mut ycdif, &mut zcdif);
                                let dsq = xdif * xdif + ydif * ydif + zdif * zdif;
                                let dtest = dsq.sqrt();
                                let rat = dtest / stepcl;
                                let iz = rat as i32 + 2;
                                if iz <= npt1 {
                                    let dtdif = dtest - xstar[iz];
                                    let ipz = iz;
                                    let prfac = if dtest > 0.0 { (xcdif * xdif + ycdif * ydif + zcdif * zdif) / dtest } else { 0.0 };
                                    let (enpot, prpot): (f32, f64) = if st.cc.ncoordconstraints > 0 && st.cc.check_coord_bonded(ia, ja, m) {
                                        (
                                            p.epot[(iz, icref)] + p.epot1[(iz, icref)] * dtdif,
                                            (0.5 * (p.epot1[(ipz, icref)] + p.epot1[(ipz - 1, icref)]) * prfac) as f64,
                                        )
                                    } else {
                                        (
                                            p.ecoul[(iz, icref)] + p.ecoul1[(iz, icref)] * dtdif,
                                            (0.5 * (p.ecoul1[(ipz, icref)] + p.ecoul1[(ipz - 1, icref)]) * prfac) as f64,
                                        )
                                    };
                                    let enwt = p.xavs[(iz, icref)] + p.xavs1[(iz, icref)] * dtdif;
                                    let prwt = 0.5 * (p.xavs1[(ipz, icref)] + p.xavs1[(ipz - 1, icref)]) * prfac;
                                    let j = jref as usize;
                                    ergold[j] = ergold[j] + enwt + enpot;
                                    ergnew[j] = (((ergnew[j] + prwt) as f64) + prpot) as f32;
                                }
                            }
                        }
                    }
                }
                (ergold, ergnew)
            })
            .collect();
        for im in 1..self.mol.nmol {
            self.erg.make_neighbour_list_from_mols(im);
            let nnbref = self.erg.nnbref;
            let (eo, en) = &per_im[(im - 1) as usize];
            for jref in 1..=nnbref {
                self.erg.ergold[jref] = eo[jref as usize];
                self.erg.ergnew[jref] = en[jref as usize];
            }
            for jref in 1..=nnbref {
                sumc += self.erg.ergold[jref];
                presst += self.erg.ergnew[jref] as f64;
            }
        }
        // Add on the triple body terms if present
        let mut energthreeb: f32 = 0.0;
        let mut pressthreeb: f64 = 0.0;
        if self.tb.nthreeb > 0 {
            self.find_triplets();
            for im in 1..=self.mol.nmol {
                energthreeb += self.en_mole_triplets(im);
                pressthreeb += self.pr_mole_triplets(im);
            }
            sumc += energthreeb;
            presst += pressthreeb;
            energthreeb /= self.mol.nmol as f32;
            self.tb.pressthreeb = (-((100.0 * (pressthreeb / 3.0)) / ((6.0222f32 * self.mol.volume) as f64))) as f32;
        }
        self.tb.energthreeb = energthreeb;
        if self.cc.ncoordconstraints > 0 {
            self.find_all_constraints();
            for im in 1..=self.mol.nmol {
                sumc += self.cc.en_mole_coord(im, &self.mol);
            }
        }
        let m = &self.mol;
        let nmol = m.nmol as f32;
        let energy = sumc / nmol;
        let pressc = ((((nmol * m.simt) as f64 - presst / 3.0) * 100.0) / ((6.0222f32 * m.volume) as f64)) as f32;
        (energy, pressc, energabs)
    }
}
