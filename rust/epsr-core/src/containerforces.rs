//! Container forces — translation of module `containerforces_mod`
//! (`containerforces_mod.f90`).
//!
//! A pair potential line in the `.pcof` file may give a maximum distance
//! `rmaxex` and amplitude `pcoefexmax` for a pair (Q-type pseudo atom,
//! real atom).  Every atom of the real type is then attracted towards the
//! *nearest* pseudo atom of that type once it is further than `rmaxex` away —
//! a way of keeping, e.g., solvent molecules inside a container.

use crate::fortran::{F1, F2};
use crate::neighbour::minimagec;
use crate::state::State;

#[derive(Clone, Debug, Default)]
pub struct ContainerForces {
    pub numcontainerforces: i32,
    pub numpairs: F1<i32>,
    pub numcentraltype: F1<i32>,
    pub numsecondarytype: F1<i32>,
    /// `pairrefstore(1|2, ic)`: central (Q) type and contained type.
    pub pairrefstore: F2<i32>,
    /// Atoms of the central type for each container force.
    pub centraltype: F2<i32>,
    pub ergold_container: F1<f32>,
    pub ergnew_container: F1<f32>,
}

impl State {
    /// `subroutine get_containerforce_list(nmol,ic,rmaxex,it,jt,...)`.
    pub fn get_containerforce_list(&mut self, ic: i32, rmaxex: f32, it: i32, jt: i32) {
        let m = &self.mol;
        let c = &mut self.cont;
        if ic == 1 {
            c.numcontainerforces = 0;
            c.ergold_container = F1::new(m.nmol.max(1) as usize);
            c.ergnew_container = F1::new(m.nmol.max(1) as usize);
        }
        c.numpairs.grow(ic as usize);
        c.numpairs[ic] = 0;
        if rmaxex <= 0.0 {
            return;
        }
        c.numcontainerforces += 1;
        let (mut typesecondary, mut typecentral) = (0, 0);
        for nt in 1..=m.ntype {
            if m.qtypea[nt] == 0 {
                if nt == it {
                    typesecondary = it;
                } else if nt == jt {
                    typesecondary = jt;
                }
            } else if nt == it {
                typecentral = it;
            } else if nt == jt {
                typecentral = jt;
            }
        }
        let mut numcentral = 0;
        let mut numsecondary = 0;
        for ia in 1..=m.natom {
            let t = m.itype[ia];
            if m.qtype[ia] == 1 && (t == it || t == jt) {
                numcentral += 1;
            } else if m.qtype[ia] == 0 && (t == it || t == jt) {
                numsecondary += 1;
            }
        }
        c.numpairs[ic] = numcentral * numsecondary;
        if c.numpairs[ic] > 0 {
            let n = c.numcontainerforces;
            c.pairrefstore.grow(2, n as usize);
            c.numcentraltype.grow(n as usize);
            c.numsecondarytype.grow(n as usize);
            c.centraltype.grow(numcentral as usize, n as usize);
            c.pairrefstore[(1, n)] = typecentral;
            c.pairrefstore[(2, n)] = typesecondary;
            c.numcentraltype[n] = numcentral;
            c.numsecondarytype[n] = numsecondary;
            let mut k = 0;
            for ia in 1..=m.natom {
                let t = m.itype[ia];
                if m.qtype[ia] == 1 && (t == it || t == jt) {
                    k += 1;
                    c.centraltype[(k, n)] = ia;
                }
            }
        } else {
            c.numpairs[ic] = 0;
            if c.numcontainerforces > 0 {
                c.numcontainerforces -= 1;
            }
        }
    }

    /// Squared distance (box units) from atom `ia` of molecule `im` to the
    /// nearest central atom of container force `ic`.
    fn nearest_central(&self, im: i32, ia: i32, ic: i32) -> f32 {
        let m = &self.mol;
        let c = &self.cont;
        let mut mindist: f32 = 0.0;
        let x1 = m.xadd[im] + m.xa[ia];
        let y1 = m.yadd[im] + m.ya[ia];
        let z1 = m.zadd[im] + m.za[ia];
        for jc in 1..=c.numcentraltype[ic] {
            let ja = c.centraltype[(jc, ic)];
            let jm = m.mref[ja];
            let mut xdif = m.xadd[jm] + m.xa[ja] - x1;
            let mut ydif = m.yadd[jm] + m.ya[ja] - y1;
            let mut zdif = m.zadd[jm] + m.za[ja] - z1;
            minimagec(&mut xdif, &mut ydif, &mut zdif);
            let dsq = xdif * xdif + ydif * ydif + zdif * zdif;
            if mindist == 0.0 || dsq < mindist {
                mindist = dsq;
            }
        }
        mindist
    }

    /// `function en_containerforce(im, ...)`: container energy of molecule `im`.
    pub fn en_containerforce(&self, im: i32) -> f32 {
        let c = &self.cont;
        if c.numcontainerforces == 0 {
            return 0.0;
        }
        let m = &self.mol;
        let p = &self.pot;
        let mut sume: f32 = 0.0;
        let mut gradlimit = self.inp.expecfpwr;
        let loggradlimit;
        if gradlimit > 0.0 {
            loggradlimit = gradlimit.ln();
        } else {
            gradlimit = 1.0;
            loggradlimit = 0.0;
        }
        for ic in 1..=c.numcontainerforces {
            for ia in m.naf[im]..=m.nal[im] {
                if m.itype[ia] == c.pairrefstore[(2, ic)] {
                    let icref = m.npairref[(m.itype[ia], c.pairrefstore[(1, ic)])];
                    let irefp = p.irefpcof[icref];
                    let expecfuse = self.inp.expecf * p.rmaxex[irefp].abs().sqrt();
                    let mindist = self.nearest_central(im, ia, ic);
                    let dtest = mindist.sqrt() * m.dmax;
                    if !p.harmonic {
                        let rtest = p.rmaxex[irefp];
                        let ddiff = (dtest - rtest) / expecfuse;
                        if gradlimit > 0.0 && ddiff >= loggradlimit {
                            sume += p.pcoefexmax[irefp].abs() * gradlimit * (1.0 + ddiff - loggradlimit);
                        } else {
                            sume += p.pcoefexmax[irefp].abs() * ddiff.exp();
                        }
                    } else {
                        let rtest = p.rmaxex[irefp] - expecfuse;
                        let ddiff = (dtest - rtest) / expecfuse;
                        if ddiff > 0.0 {
                            if gradlimit > 0.0 && ddiff > gradlimit {
                                sume += 0.5 * p.pcoefexmax[irefp].abs() * gradlimit * (2.0 * ddiff - gradlimit);
                            } else {
                                sume += 0.5 * p.pcoefexmax[irefp].abs() * ddiff * ddiff;
                            }
                        }
                    }
                }
            }
        }
        sume
    }

    /// `subroutine dist_containerforce(im, ..., potexmax, dmax, addwav)`: count
    /// contained atoms of molecule `im` that are beyond `rmaxex`.
    pub fn dist_containerforce(&mut self, im: i32, addwav: f32) {
        if self.cont.numcontainerforces == 0 {
            return;
        }
        for ic in 1..=self.cont.numcontainerforces {
            let addpcoef = addwav;
            for ia in self.mol.naf[im]..=self.mol.nal[im] {
                if self.mol.itype[ia] == self.cont.pairrefstore[(2, ic)] {
                    let icref = self.mol.npairref[(self.mol.itype[ia], self.cont.pairrefstore[(1, ic)])];
                    let irefp = self.pot.irefpcof[icref];
                    let mindist = self.nearest_central(im, ia, ic);
                    let dtest = mindist.sqrt() * self.mol.dmax;
                    if dtest > self.pot.rmaxex[irefp] {
                        self.pot.potexmax[icref] += addpcoef;
                    }
                }
            }
        }
    }
}
