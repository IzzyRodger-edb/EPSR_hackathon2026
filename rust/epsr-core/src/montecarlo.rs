//! Monte Carlo moves.
//!
//! Translation of `molfit`, `movei`, `movea`, `movem`, `rotati`, `rotatm`
//! (`molfit.f`).
//!
//! One call to [`State::molfit`] is one Monte Carlo cycle: every molecule is
//! visited once (cell by cell) and given a whole-molecule rotation and
//! translation, accepted with the Metropolis criterion on the change of
//! intermolecular energy.  Every `intra` cycles (`nisum >= intra`) the cycle
//! instead "shakes" the molecules: random displacements of all atoms of a
//! molecule (intramolecular move) followed by rotations of the rotational
//! groups.  Every `rotfreq` cycles the rotational groups are also rotated.
//! Step sizes are adapted per molecule type so that the acceptance ratio
//! tends to `rejrate`.
//!
//! The order of random numbers drawn is identical to the Fortran, so given
//! the same `.ato` file (which stores the generator state) the trajectory is
//! reproduced exactly.

use crate::fortran::fmath::sign;
use crate::state::State;

/// Metropolis acceptance threshold used throughout `molfit`.
#[inline(always)]
fn etest_of(ediff: f32) -> f32 {
    if ediff <= 0.0 {
        10.0
    } else if ediff < 20.0 {
        (-ediff).exp()
    } else {
        0.0
    }
}

impl State {
    /// `subroutine movea(ia, stepmithis)`: random displacement of one atom.
    fn movea(&mut self, ia: i32, stepmithis: f32) {
        let stepsz = stepmithis;
        let r = self.rng.ran1();
        self.mol.za[ia] += stepsz * (2.0 * r - 1.0);
        let r = self.rng.ran1();
        self.mol.ya[ia] += stepsz * (2.0 * r - 1.0);
        let r = self.rng.ran1();
        self.mol.xa[ia] += stepsz * (2.0 * r - 1.0);
    }

    /// `subroutine movei(im, stepmithis)`: displace every atom of molecule `im`.
    fn movei(&mut self, im: i32, stepmithis: f32) {
        for ia in self.mol.naf[im]..=self.mol.nal[im] {
            self.movea(ia, stepmithis);
        }
    }

    /// `subroutine movem(im, stepsz)`: translate molecule `im`.  As in the
    /// Fortran, `stepsz` (passed by reference) is doubled until the move
    /// actually changes the position.
    fn movem(&mut self, im: i32, stepsz: &mut f32) {
        let m = &mut self.mol;
        m.xs = m.xadd[im];
        m.ys = m.yadd[im];
        m.zs = m.zadd[im];
        let r = self.rng.ran1();
        m.zadd[im] += *stepsz * (2.0 * r - 1.0);
        let r = self.rng.ran1();
        m.yadd[im] += *stepsz * (2.0 * r - 1.0);
        let r = self.rng.ran1();
        m.xadd[im] += *stepsz * (2.0 * r - 1.0);
        let mut nochange = m.xs == m.xadd[im] && m.ys == m.yadd[im] && m.zs == m.zadd[im];
        while nochange {
            *stepsz *= 2.0;
            let r = self.rng.ran1();
            m.zadd[im] += *stepsz * (2.0 * r - 1.0);
            let r = self.rng.ran1();
            m.yadd[im] += *stepsz * (2.0 * r - 1.0);
            let r = self.rng.ran1();
            m.xadd[im] += *stepsz * (2.0 * r - 1.0);
            nochange = m.xs == m.xadd[im] && m.ys == m.yadd[im] && m.zs == m.zadd[im];
        }
    }

    /// `subroutine rotati(im, ig, steprithis, actualstep)`: rotate the atoms of
    /// rotational group `ig` about its axis by a random angle.
    fn rotati(&mut self, ig: i32, steprithis: f32) {
        let pi: f32 = std::f32::consts::PI;
        let r = self.rng.ran1();
        let phi = pi * steprithis * (2.0 * r - 1.0);
        let cphi = phi.cos();
        let sphi = sign((1.0 - cphi * cphi).sqrt(), phi);
        let m = &mut self.mol;
        let a1 = m.labaxi[(ig, 1)];
        let a2 = m.labaxi[(ig, 2)];
        let mut pvector = [m.xa[a2] - m.xa[a1], m.ya[a2] - m.ya[a1], m.za[a2] - m.za[a1]];
        let rab2 = crate::ato::scapro(&pvector, &pvector);
        if rab2 > 0.0 {
            let rab = rab2.sqrt();
            for v in pvector.iter_mut() {
                *v /= rab;
            }
            for ia in m.rotf[ig]..=m.rotl[ig] {
                let nat = m.labrot[ia];
                let nax = m.labaxi[(ig, 1)];
                let xvector = [m.xa[nat] - m.xa[nax], m.ya[nat] - m.ya[nax], m.za[nat] - m.za[nax]];
                let frac = crate::ato::scapro(&pvector, &xvector);
                let xpvector = [frac * pvector[0], frac * pvector[1], frac * pvector[2]];
                let xrvector = [xvector[0] - xpvector[0], xvector[1] - xpvector[1], xvector[2] - xpvector[2]];
                let xrpvector = crate::ato::vecpro(&pvector, &xrvector);
                m.xa[nat] = m.xa[nax] + xpvector[0];
                m.ya[nat] = m.ya[nax] + xpvector[1];
                m.za[nat] = m.za[nax] + xpvector[2];
                m.xa[nat] = m.xa[nat] + xrvector[0] * cphi + xrpvector[0] * sphi;
                m.ya[nat] = m.ya[nat] + xrvector[1] * cphi + xrpvector[1] * sphi;
                m.za[nat] = m.za[nat] + xrvector[2] * cphi + xrpvector[2] * sphi;
            }
        }
    }

    /// `subroutine rotatm(im, stepsz)`: rotate molecule `im` about x, y and z.
    fn rotatm(&mut self, im: i32, stepsz: f32) {
        let pi: f32 = std::f32::consts::PI;
        let iref1 = self.mol.naf[im];
        let iref2 = self.mol.nal[im];
        if iref2 > iref1 {
            let r = self.rng.ran1();
            let phi = pi * stepsz * (2.0 * r - 1.0);
            let cphi = phi.cos();
            let sphi = sign((1.0 - cphi * cphi).sqrt(), phi);
            let m = &mut self.mol;
            for ia in iref1..=iref2 {
                let xasave = m.xa[ia] * cphi - m.ya[ia] * sphi;
                m.ya[ia] = m.xa[ia] * sphi + m.ya[ia] * cphi;
                m.xa[ia] = xasave;
            }
            let r = self.rng.ran1();
            let phi = pi * stepsz * (2.0 * r - 1.0);
            let cphi = phi.cos();
            let sphi = sign((1.0 - cphi * cphi).sqrt(), phi);
            let m = &mut self.mol;
            for ia in iref1..=iref2 {
                let zasave = m.za[ia] * cphi - m.xa[ia] * sphi;
                m.xa[ia] = m.za[ia] * sphi + m.xa[ia] * cphi;
                m.za[ia] = zasave;
            }
            let r = self.rng.ran1();
            let phi = pi * stepsz * (2.0 * r - 1.0);
            let cphi = phi.cos();
            let sphi = sign((1.0 - cphi * cphi).sqrt(), phi);
            let m = &mut self.mol;
            for ia in iref1..=iref2 {
                let yasave = m.ya[ia] * cphi - m.za[ia] * sphi;
                m.za[ia] = m.ya[ia] * sphi + m.za[ia] * cphi;
                m.ya[ia] = yasave;
            }
        }
    }

    fn save_coords(&self, im: i32, xr: &mut Vec<f32>, yr: &mut Vec<f32>, zr: &mut Vec<f32>) {
        let m = &self.mol;
        xr.clear();
        yr.clear();
        zr.clear();
        for nat in m.naf[im]..=m.nal[im] {
            xr.push(m.xa[nat]);
            yr.push(m.ya[nat]);
            zr.push(m.za[nat]);
        }
    }
    fn restore_coords(&mut self, im: i32, xr: &[f32], yr: &[f32], zr: &[f32]) {
        let m = &mut self.mol;
        let naf = m.naf[im];
        for (k, nat) in (naf..=m.nal[im]).enumerate() {
            m.xa[nat] = xr[k];
            m.ya[nat] = yr[k];
            m.za[nat] = zr[k];
        }
    }

    /// Internal rotations of the rotational groups of molecule `im` with a
    /// Metropolis test on the intra + intermolecular energy change.
    fn internal_rotation_move(&mut self, im: i32, temp: f32, fadd: &mut f32, faddb: &mut f32, bufs: &mut (Vec<f32>, Vec<f32>, Vec<f32>)) {
        let mt = self.mol.mtype[im];
        let rotstepmin: f32 = 1.0e-5;
        if self.mol.internrot[mt] == 0.0 {
            self.mol.internrot[mt] = rotstepmin;
        }
        self.mol.ninter[(1, mt)] += 1;
        self.save_coords(im, &mut bufs.0, &mut bufs.1, &mut bufs.2);
        self.mol.xs = self.mol.xadd[im];
        self.mol.ys = self.mol.yadd[im];
        self.mol.zs = self.mol.zadd[im];
        for ig in self.mol.nrf[im]..=self.mol.nrl[im] {
            let step = self.mol.internrot[mt];
            self.rotati(ig, step);
        }
        if !self.check_tethered(im, self.mol.tetheredatom[im]) {
            // tether violated: reject without an energy evaluation
            self.mol.nstepreject[im] += 1;
            let (xr, yr, zr) = (std::mem::take(&mut bufs.0), std::mem::take(&mut bufs.1), std::mem::take(&mut bufs.2));
            self.restore_coords(im, &xr, &yr, &zr);
            *bufs = (xr, yr, zr);
            self.mol.xadd[im] = self.mol.xs;
            self.mol.yadd[im] = self.mol.ys;
            self.mol.zadd[im] = self.mol.zs;
            self.mol.ninter[(2, mt)] += 1;
            return;
        }
        let (mut sv, mut sa, mut sd) = (0f32, 0f32, 0f32);
        let fdiffb = self.en_change_intra(im, &mut sv, &mut sa, &mut sd);
        let fdiff = self.en_change(im);
        let ediff = (fdiff + fdiffb) * temp;
        let etest = etest_of(ediff);
        if self.rng.ran1() < etest {
            self.mol.nstepreject[im] = 0;
            *fadd += fdiff;
            *faddb += fdiffb;
            self.erg.save_energies();
            self.save_threeb_energies(im);
            self.cc.save_coord_energies(im);
        } else {
            self.mol.nstepreject[im] += 1;
            let (xr, yr, zr) = (std::mem::take(&mut bufs.0), std::mem::take(&mut bufs.1), std::mem::take(&mut bufs.2));
            self.restore_coords(im, &xr, &yr, &zr);
            *bufs = (xr, yr, zr);
            self.mol.ninter[(2, mt)] += 1;
        }
        self.erg.save_erg_pairs(im);
    }

    /// `subroutine molfit(fadds, faddbs, nisum, intra, rotfreq, ireset)`:
    /// one Monte Carlo cycle.  `fadds`/`faddbs` accumulate the inter- and
    /// intramolecular energy changes per molecule; `nisum` is the cycle
    /// counter (common `/epsrinput_par/`).
    pub fn molfit(&mut self, fadds: &mut f32, faddbs: &mut f32, intra: i32, rotfreq: i32, ireset: i32) {
        let temp = 1.0 / self.mol.simt;
        self.mol.tol = self.mol.tolout / self.mol.dmax;
        for j in 1..=2 {
            for imoltype in 1..=self.mol.nmoltype {
                self.mol.nmovei[(j, imoltype)] = 0;
                self.mol.ninter[(j, imoltype)] = 0;
                self.mol.nwhole[(j, imoltype)] = 0;
                self.mol.nmovem[(j, imoltype)] = 0;
            }
        }
        let mut trstepmax = 1.0 / self.mol.natomr as f32;
        trstepmax = trstepmax.powf(1.0f32 / 3.0f32);
        let trstepmin = (0.0001 * trstepmax).max(1.0e-6);
        let rotstepmax: f32 = 1.0;
        let rotstepmin: f32 = 1.0e-5;
        let mut fadd: f32 = 0.0;
        let mut faddb: f32 = 0.0;
        self.mol.tempang = self.mol.tempangwt;
        self.mol.tempdih = self.mol.tempdihwt * self.mol.tempb;
        let mut bufs: (Vec<f32>, Vec<f32>, Vec<f32>) = (Vec::new(), Vec::new(), Vec::new());
        let nisum = self.inp.nisum;
        if nisum >= intra {
            // ---- molecule "shake": intramolecular moves of every molecule
            for im in 1..=self.mol.nmol {
                self.erg.make_neighbour_list_from_mols(im);
                let mt = self.mol.mtype[im];
                if self.mol.nato[im] > 1 && self.mol.stepmi > 0.0 {
                    let mut stepmithis = self.mol.interntrans[mt];
                    if self.mol.stepmi > 0.0 && stepmithis == 0.0 {
                        stepmithis = trstepmin;
                    }
                    self.mol.interntrans[mt] = stepmithis;
                    self.mol.xs = self.mol.xadd[im];
                    self.mol.ys = self.mol.yadd[im];
                    self.mol.zs = self.mol.zadd[im];
                    self.save_coords(im, &mut bufs.0, &mut bufs.1, &mut bufs.2);
                    self.mol.nmovei[(1, mt)] += 1;
                    self.movei(im, stepmithis);
                    let (mut sv, mut sa, mut sd) = (0f32, 0f32, 0f32);
                    let fdiffb = self.en_change_intra(im, &mut sv, &mut sa, &mut sd);
                    let ediff = fdiffb * temp;
                    let etest = etest_of(ediff);
                    let rtest = self.rng.ran1();
                    if rtest < etest {
                        faddb += fdiffb;
                        self.erg.save_energies();
                        self.centreofmassi(im);
                    } else {
                        self.mol.nmovei[(2, mt)] += 1;
                        let (xr, yr, zr) = (std::mem::take(&mut bufs.0), std::mem::take(&mut bufs.1), std::mem::take(&mut bufs.2));
                        self.restore_coords(im, &xr, &yr, &zr);
                        bufs = (xr, yr, zr);
                    }
                    self.erg.save_erg_pairs(im);
                }
                if self.mol.stepri > 0.0 && self.mol.ngrot[im] > 0 {
                    self.internal_rotation_move(im, temp, &mut fadd, &mut faddb, &mut bufs);
                } else if self.mol.ngrot[im] == 0 {
                    self.mol.internrot[mt] = 0.0;
                }
            }
        } else {
            // ---- normal cycle: whole-molecule rotation + translation, cell by cell
            for ic in 1..=self.cells.ncell {
                if self.cells.nmolsincell[ic] > 0 {
                    for imref in 1..=self.cells.nmolsincell[ic] {
                        let im = self.cells.molsincell[(imref, ic)];
                        let mt = self.mol.mtype[im];
                        self.erg.make_neighbour_list_from_mols(im);
                        if self.mol.stepri > 0.0 && self.mol.ngrot[im] > 0 && nisum == rotfreq * (nisum / rotfreq) {
                            // NOT exercised by the reference test case (see docs/VALIDATION.md)
                            self.internal_rotation_move(im, temp, &mut fadd, &mut faddb, &mut bufs);
                        } else if self.mol.ngrot[im] == 0 {
                            self.mol.internrot[mt] = 0.0;
                        }
                        self.save_coords(im, &mut bufs.0, &mut bufs.1, &mut bufs.2);
                        self.mol.xs = self.mol.xadd[im];
                        self.mol.ys = self.mol.yadd[im];
                        self.mol.zs = self.mol.zadd[im];
                        let (xs, ys, zs) = (self.mol.xs, self.mol.ys, self.mol.zs);
                        if self.mol.stepr > 0.0 && self.mol.nato[im] > 1 {
                            if self.mol.wholerot[mt] <= 0.0 {
                                self.mol.wholerot[mt] = rotstepmin;
                            }
                            self.mol.nwhole[(1, mt)] += 1;
                            let step = self.mol.wholerot[mt];
                            self.rotatm(im, step);
                        } else {
                            self.mol.wholerot[mt] = 0.0;
                        }
                        if self.mol.stepm > 0.0 {
                            if self.mol.wholetrans[mt] <= 0.0 {
                                self.mol.wholetrans[mt] = trstepmin;
                            }
                            self.mol.nmovem[(1, mt)] += 1;
                            let mut step = self.mol.wholetrans[mt];
                            self.movem(im, &mut step);
                            self.mol.wholetrans[mt] = step;
                        } else {
                            self.mol.wholetrans[mt] = 0.0;
                        }
                        if !self.check_tethered(im, self.mol.tetheredatom[im]) {
                            // tether violated: reject without an energy evaluation
                            self.mol.nstepreject[im] += 1;
                            let (xr, yr, zr) = (std::mem::take(&mut bufs.0), std::mem::take(&mut bufs.1), std::mem::take(&mut bufs.2));
                            self.restore_coords(im, &xr, &yr, &zr);
                            bufs = (xr, yr, zr);
                            self.mol.xadd[im] = xs;
                            self.mol.yadd[im] = ys;
                            self.mol.zadd[im] = zs;
                            if self.mol.stepr > 0.0 {
                                self.mol.nwhole[(2, mt)] += 1;
                            }
                            if self.mol.stepm > 0.0 {
                                self.mol.nmovem[(2, mt)] += 1;
                            }
                            continue;
                        }
                        let fdiffb: f32 = 0.0;
                        let fdiff = self.en_change(im);
                        let ediff = (fdiff + fdiffb) * temp;
                        let etest = etest_of(ediff);
                        if self.rng.ran1() < etest {
                            self.mol.nstepreject[im] = 0;
                            fadd += fdiff;
                            faddb += fdiffb;
                            self.erg.save_energies();
                            self.save_threeb_energies(im);
                        } else {
                            self.mol.nstepreject[im] += 1;
                            let (xr, yr, zr) = (std::mem::take(&mut bufs.0), std::mem::take(&mut bufs.1), std::mem::take(&mut bufs.2));
                            self.restore_coords(im, &xr, &yr, &zr);
                            bufs = (xr, yr, zr);
                            self.mol.xadd[im] = xs;
                            self.mol.yadd[im] = ys;
                            self.mol.zadd[im] = zs;
                            if self.mol.stepr > 0.0 {
                                self.mol.nwhole[(2, mt)] += 1;
                            }
                            if self.mol.stepm > 0.0 {
                                self.mol.nmovem[(2, mt)] += 1;
                            }
                        }
                        self.erg.save_erg_pairs(im);
                    }
                }
            }
        }
        if nisum >= intra {
            self.inp.nisum = 1;
        } else {
            self.inp.nisum = nisum + 1;
        }
        if ireset != 1 {
            let rejrate = self.dist.rejrate;
            let m = &mut self.mol;
            let adapt = |n1: i32, n2: i32| -> f32 {
                if n2 > 0 {
                    let rat = rejrate * n1 as f32 / n2 as f32;
                    rat.clamp(0.9, 1.1)
                } else if n1 > 0 {
                    1.1
                } else {
                    1.0
                }
            };
            for imt in 1..=m.nmoltype {
                let rat = adapt(m.nmovei[(1, imt)], m.nmovei[(2, imt)]);
                m.interntrans[imt] *= rat;
                if m.interntrans[imt] > trstepmax {
                    m.interntrans[imt] = trstepmax;
                }
                if m.interntrans[imt] < trstepmin && m.interntrans[imt] != 0.0 {
                    m.interntrans[imt] = trstepmin;
                }
                let rat = adapt(m.ninter[(1, imt)], m.ninter[(2, imt)]);
                m.internrot[imt] *= rat;
                if m.internrot[imt] > rotstepmax {
                    m.internrot[imt] = rotstepmax;
                }
                if m.internrot[imt] < rotstepmin && m.internrot[imt] != 0.0 {
                    m.internrot[imt] = rotstepmin;
                }
                let rat = adapt(m.nwhole[(1, imt)], m.nwhole[(2, imt)]);
                m.wholerot[imt] *= rat;
                if m.wholerot[imt] > rotstepmax {
                    m.wholerot[imt] = rotstepmax;
                }
                if m.wholerot[imt] < rotstepmin && m.wholerot[imt] != 0.0 {
                    m.wholerot[imt] = rotstepmin;
                }
                let rat = adapt(m.nmovem[(1, imt)], m.nmovem[(2, imt)]);
                m.wholetrans[imt] *= rat;
                if m.wholetrans[imt] > trstepmax {
                    m.wholetrans[imt] = trstepmax;
                }
                if m.wholetrans[imt] < trstepmin && m.wholetrans[imt] != 0.0 {
                    m.wholetrans[imt] = trstepmin;
                }
            }
        }
        let nmol = self.mol.nmol as f32;
        *fadds += fadd / nmol;
        *faddbs += faddb / nmol;
    }
}
