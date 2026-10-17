//! Three-body (triplet) potentials — translation of module `threebody_mod`
//! (`threebody_mod.f90`) and of the three-body parts of `read_pcof` /
//! `save_files_others` (`read_and_save_files.f`).
//!
//! A `.pcof` line after the pair potentials
//! ```text
//! tb  A1 A2 A3  amp  awidth  rwidth  rmax12 rmax32  s12 s23 s13  ang1 amp1 [ang2 amp2 ...]
//! ```
//! defines a group of triplets A1–A2–A3 (A2 = centre) with |r1−r2| < rmax12,
//! |r3−r2| < rmax32 and the "same molecule" flags `s12`, `s23`, `s13`
//! (1 = same molecule, 0 = different molecules).  The energy of a triplet is
//! `amp * ampN * 0.5*((cos θ − cosN)/awidth)²` for the requested angle N
//! whose cosine is closest to cos θ.  Triplets are listed per centre molecule
//! by `find_triplets` (called by `fmolec1` and `press`).
//!
//! Fortran quirks reproduced (see the `NB:` comments):
//! * `read_pcof` stores the *cosines* of the requested angles and sets
//!   `awidth = 1.0`, but `save_files_others` (called at the end of every
//!   iteration) re-reads the line from the symbol table and stores the angles
//!   in *radians* and `awidth` in radians.  From the second iteration of a run
//!   on, `en_mole_triplets` therefore compares cos θ with angles in radians.
//! * If the first and third atom types are equal, only triplets with
//!   `ia3 > ia1` are kept.

use crate::error::{EpsrError, Result};
use crate::fortran::fmath::PI4ATAN;
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::{read_int, read_real};
use crate::fortran::output::OutFile;
use crate::fortran::text::{concatenate, fstr, parse, Words};
use crate::fortran::{F1, F2};
use crate::neighbour::minimagec;
use crate::state::State;
use crate::symbols::{LEN_SYMBOL, LEN_VALUE};

/// Module variables of `threebody_mod` (plus the `runepsr17` locals
/// `energthreeb`, `pressthreeb` written to the `.out` file).
#[derive(Clone, Debug, Default)]
pub struct ThreeBody {
    pub nthreeb: i32,
    pub nthreebangsmax: i32,
    pub nthreebneighbmolmax: i32,
    /// `threebgrp(1..3, id)`: the three atom types.
    pub threebgrp: F2<String>,
    pub threebamp: F1<f32>,
    /// `threebrmax(1|2, id)`: maximum 1–2 and 3–2 distances (box units between
    /// `normalise_coordinates` and `unnormalise_coordinates`).
    pub threebrmax: F2<f32>,
    pub threebrwidth: F1<f32>,
    pub threebawidth: F1<f32>,
    pub nthreebangs: F1<i32>,
    pub threebangs: F2<f32>,
    pub threebamps: F2<f32>,
    pub nthreebneighbmol: F1<i32>,
    pub threebneighbmol: F2<i32>,
    pub threebsamemol: F2<i32>,
    /// `ntripletsinmol(im, id)`, `firsttripletinmol(im, id)`.
    pub ntripletsinmol: F2<i32>,
    pub firsttripletinmol: F2<i32>,
    /// `tripletlist(1..3, it, id)`: stored as `tripletlist[id-1][it-1] = [ia1, ia2, ia3]`.
    pub tripletlist: Vec<Vec<[i32; 3]>>,
    pub ergold_triplet: F1<f32>,
    pub ergnew_triplet: F1<f32>,
    pub sume_triplet_dist: F1<f32>,
    pub en_mole_triplet_dist: F1<f32>,
    pub consthreeb: F1<f32>,
    /// `energthreeb` / `pressthreeb` of `runepsr17` (set by `press` only when `nthreeb > 0`).
    pub energthreeb: f32,
    pub pressthreeb: f32,
}

fn input_err(what: &str, s: &str) -> EpsrError {
    EpsrError::Input(format!("read_pcof: cannot read {what} from \"{s}\" in three-body line"))
}

impl State {
    /// The `tb` branch of the three-body / coordination-constraint loop of
    /// `read_pcof`: one three-body group (`nwords > 13` already checked).
    pub(crate) fn read_pcof_threebody(&mut self, sometext: &str, w: &Words) -> Result<()> {
        let pi = PI4ATAN;
        let piconv = pi / 180.0;
        let nwords = w.nwords();
        let tb = &mut self.tb;
        tb.nthreeb += 1;
        let n = tb.nthreeb;
        let nu = n as usize;
        tb.threebgrp.grow(3, nu);
        tb.threebamp.grow(nu);
        tb.threebrmax.grow(2, nu);
        tb.threebawidth.grow(nu);
        tb.threebrwidth.grow(nu);
        tb.nthreebangs.grow(nu);
        tb.threebsamemol.grow(3, nu);
        let word = |k: i32| w.word(sometext, k);
        let rr = |k: i32| read_real(word(k)).map_err(|_| input_err("a real", word(k)));
        let mut iref = 1;
        for i in 1..=3 {
            iref += 1;
            tb.threebgrp[(i, n)] = fstr(word(iref), 3);
            let m = &self.mol;
            if !(1..=m.ntype).any(|it| m.type_[it] == tb.threebgrp[(i, n)]) {
                // find_triplets would print a message and stop processing groups
                return crate::error::unsupported("three-body group whose atom types are not all in the .ato file");
            }
        }
        iref += 1;
        tb.threebamp[n] = rr(iref)?;
        iref += 1;
        tb.threebawidth[n] = rr(iref)?;
        tb.threebawidth[n] = 1.0;
        iref += 1;
        tb.threebrwidth[n] = rr(iref)?;
        for i in 1..=2 {
            iref += 1;
            tb.threebrmax[(i, n)] = rr(iref)?;
        }
        for i in 1..=3 {
            iref += 1;
            tb.threebsamemol[(i, n)] = read_int(word(iref)).map_err(|_| input_err("an integer", word(iref)))?;
        }
        tb.nthreebangs[n] = (nwords - iref) / 2;
        let msg = {
            let mut lw = ListWriter::new();
            lw.chars("read_pcof> No. of tb angles: ").int(n).int(tb.nthreebangs[n]);
            lw.finish()
        };
        if tb.nthreebangs[n] > 9 {
            tb.nthreebangs[n] = 9;
        }
        if tb.nthreebangs[n] > tb.nthreebangsmax {
            tb.nthreebangsmax = tb.nthreebangs[n];
        }
        tb.threebangs.grow(tb.nthreebangsmax as usize, nu);
        tb.threebamps.grow(tb.nthreebangsmax as usize, nu);
        let mut sumamp: f32 = 0.0;
        for i in 1..=tb.nthreebangs[n] {
            iref += 1;
            tb.threebangs[(i, n)] = rr(iref)?;
            tb.threebangs[(i, n)] = (piconv * tb.threebangs[(i, n)]).cos();
            iref += 1;
            tb.threebamps[(i, n)] = rr(iref)?;
            tb.threebamps[(i, n)] = tb.threebamps[(i, n)].abs();
            sumamp += tb.threebamps[(i, n)];
        }
        if sumamp > 0.0 {
            // (normalisation of the amplitudes is commented out in the Fortran)
        } else {
            tb.threebamp[n] = 0.0;
        }
        self.log(&msg);
        // the new symbol
        self.sym.nsymbol += 1;
        let ns = self.sym.nsymbol;
        let irefs = self.sym.ireff(ns, self.sym.nfiles);
        self.sym.symbol[ns] = fstr(w.word(sometext, 1), LEN_SYMBOL);
        self.sym.iflag[ns] = 0;
        self.sym.symbol_default[ns] = String::new();
        self.sym.symbol_type[ns] = 'c';
        let (val, _) = concatenate(2, nwords, w, sometext, LEN_VALUE);
        self.sym.symbol_value[irefs] = val;
        self.sym.symbol_text[ns] = "Three atom labels, amplitude, width, angle width, two dist., and angles".into();
        Ok(())
    }

    /// The `tb` part of the `.pcof` section of `save_files_others`: re-read
    /// each group from its symbol value and write it.
    pub(crate) fn write_pcof_threebody(&mut self, out: &mut OutFile) -> Result<()> {
        if self.tb.nthreeb <= 0 {
            return Ok(());
        }
        let piconv = PI4ATAN / 180.0;
        let mut icount = self.sym.nsymboldata + self.pot.npairspcof;
        for ithreeb in 1..=self.tb.nthreeb {
            icount += 1;
            let iref = self.sym.ireff(icount, self.sym.nfiles);
            let sometext = self.sym.symbol_value[iref].clone();
            let w = parse(&sometext, 11);
            let nwords = w.nwords();
            let tb = &mut self.tb;
            if nwords > 12 {
                let word = |k: i32| w.word(&sometext, k);
                let rr = |k: i32| read_real(word(k)).map_err(|_| input_err("a real", word(k)));
                for i in 1..=3 {
                    tb.threebgrp[(i, ithreeb)] = fstr(word(i), 3);
                }
                let mut i = 4;
                tb.threebamp[ithreeb] = rr(i)?;
                i += 1;
                // NB: here the width is stored in radians (read_pcof stores 1.0)
                tb.threebawidth[ithreeb] = rr(i)?;
                tb.threebawidth[ithreeb] = piconv * tb.threebawidth[ithreeb];
                i += 1;
                tb.threebrwidth[ithreeb] = rr(i)?;
                for j in 1..=2 {
                    i += 1;
                    tb.threebrmax[(j, ithreeb)] = rr(i)?;
                }
                for j in 1..=3 {
                    i += 1;
                    tb.threebsamemol[(j, ithreeb)] = read_int(word(i)).map_err(|_| input_err("an integer", word(i)))?;
                }
                tb.nthreebangs[ithreeb] = (nwords - i) / 2;
                if tb.nthreebangs[ithreeb] > 9 {
                    tb.nthreebangs[ithreeb] = 9;
                }
                if tb.nthreebangs[ithreeb] > tb.nthreebangsmax {
                    // NB: the Fortran sets nthreebangsmax = nthreebangs(nthreeb) here
                    // (and then reallocates); never reached since read_pcof already
                    // took the maximum of the same values.
                    return crate::error::unsupported("three-body group with more angles at save time than at read time");
                }
                for j in 1..=tb.nthreebangs[ithreeb] {
                    i += 1;
                    // NB: angles stored in radians here, not as cosines as in read_pcof.
                    tb.threebangs[(j, ithreeb)] = rr(i)?;
                    tb.threebangs[(j, ithreeb)] = piconv * tb.threebangs[(j, ithreeb)];
                    i += 1;
                    tb.threebamps[(j, ithreeb)] = rr(i)?;
                }
            }
            let tb = &self.tb;
            let n2 = 2 * tb.nthreebangs[ithreeb];
            // NB: the Fortran writes the repeat count with i1.1; for 2*n >= 10 this
            // gives '*', i.e. an unlimited-repeat group, which formats the
            // same items identically.
            let formatline = format!("(a,3(1x,a3),1x,f10.2,2(1x,f7.2),2(1x,f10.5),3(1x,i1),{}(1x,f7.2))", n2);
            let mut v: Vec<FV> = vec![FV::S("tb")];
            for j in 1..=3 {
                v.push(FV::Ch(&tb.threebgrp[(j, ithreeb)], 3));
            }
            v.push(FV::R(tb.threebamp[ithreeb]));
            v.push(FV::R(tb.threebawidth[ithreeb] / piconv));
            v.push(FV::R(tb.threebrwidth[ithreeb]));
            for j in 1..=2 {
                v.push(FV::R(tb.threebrmax[(j, ithreeb)]));
            }
            for j in 1..=3 {
                v.push(FV::I(tb.threebsamemol[(j, ithreeb)] as i64));
            }
            for j in 1..=tb.nthreebangs[ithreeb] {
                v.push(FV::R(tb.threebangs[(j, ithreeb)] / piconv));
                v.push(FV::R(tb.threebamps[(j, ithreeb)]));
            }
            out.line(&fmt(&formatline, &v));
        }
        Ok(())
    }

    /// `subroutine find_triplets(nmol,natom,xadd,yadd,zadd,xa,ya,za,mref,itype,ntype,type)`.
    pub fn find_triplets(&mut self) {
        let m = &self.mol;
        let nmol = m.nmol;
        let tb = &mut self.tb;
        // (ergold_triplet / ergnew_triplet are reallocated: contents undefined until set)
        tb.ergold_triplet.grow(nmol.max(1) as usize);
        tb.ergnew_triplet.grow(nmol.max(1) as usize);
        if tb.nthreeb <= 0 {
            return;
        }
        let nthreeb = tb.nthreeb as usize;
        tb.ntripletsinmol.grow(nmol as usize, nthreeb);
        tb.firsttripletinmol.grow(nmol as usize, nthreeb);
        tb.sume_triplet_dist.grow(nthreeb);
        tb.en_mole_triplet_dist.grow(nthreeb);
        tb.tripletlist.resize(nthreeb, Vec::new());
        tb.nthreebneighbmolmax = 4;
        tb.nthreebneighbmol.grow(nmol as usize);
        tb.threebneighbmol.grow(tb.nthreebneighbmolmax as usize, nmol as usize);
        for im in 1..=nmol {
            tb.nthreebneighbmol[im] = 0;
        }
        let mut msgs: Vec<String> = Vec::new();
        for id in 1..=tb.nthreeb {
            for im1 in 1..=nmol {
                tb.ntripletsinmol[(im1, id)] = 0;
                tb.firsttripletinmol[(im1, id)] = 0;
            }
            let list = &mut tb.tripletlist[(id - 1) as usize];
            list.clear();
            let mut ntriplets: i32 = 0;
            let (mut ifound1, mut ifound2, mut ifound3) = (0, 0, 0);
            let mut it = 0;
            while it < m.ntype && (ifound1 == 0 || ifound2 == 0 || ifound3 == 0) {
                it += 1;
                if m.type_[it] == tb.threebgrp[(1, id)] {
                    ifound1 = it;
                }
                if m.type_[it] == tb.threebgrp[(2, id)] {
                    ifound2 = it;
                }
                if m.type_[it] == tb.threebgrp[(3, id)] {
                    ifound3 = it;
                }
            }
            if ifound1 == 0 || ifound2 == 0 || ifound3 == 0 {
                // The Fortran prints a message and returns, leaving the remaining
                // groups with the triplet lists of the previous call.  Not reached:
                // read_pcof_threebody rejects such groups (unsupported).
                self.log(&format!("\nfind_triplets> atoms not matched to those in .ato file for triplet{id:>3}"));
                return;
            }
            let (ltyp1, ltyp2, ltyp3) = (ifound1, ifound2, ifound3);
            let mut atom1list: Vec<i32> = Vec::new();
            let mut atom2list: Vec<i32> = Vec::new();
            let mut atom3list: Vec<i32> = Vec::new();
            for ia2 in 1..=m.natom {
                if m.itype[ia2] == ltyp1 {
                    atom1list.push(ia2);
                }
                if m.itype[ia2] == ltyp2 {
                    atom2list.push(ia2);
                }
                if m.itype[ia2] == ltyp3 {
                    atom3list.push(ia2);
                }
            }
            msgs.push(format!(" {} {} {} {}", id, atom1list.len(), atom2list.len(), atom3list.len()));
            let natom1list = atom1list.len() as i32;
            let s1 = tb.threebsamemol[(1, id)];
            let s2 = tb.threebsamemol[(2, id)];
            let s3 = tb.threebsamemol[(3, id)];
            let rmax1 = tb.threebrmax[(1, id)];
            let rmax2 = tb.threebrmax[(2, id)];
            for &ia2 in atom2list.iter() {
                let im2 = m.mref[ia2];
                let x2 = m.xa[ia2] + m.xadd[im2];
                let y2 = m.ya[ia2] + m.yadd[im2];
                let z2 = m.za[ia2] + m.zadd[im2];
                let mut iref1: i32 = 0;
                while iref1 < natom1list {
                    let mut foundia1 = false;
                    let mut ia1 = 0;
                    let mut im1 = 0;
                    while iref1 < natom1list && !foundia1 {
                        iref1 += 1;
                        ia1 = atom1list[(iref1 - 1) as usize];
                        im1 = m.mref[ia1];
                        if ia1 != ia2 && ((s1 == 1 && im1 == im2) || (s1 == 0 && im1 != im2)) {
                            let x1 = m.xa[ia1] + m.xadd[im1];
                            let y1 = m.ya[ia1] + m.yadd[im1];
                            let z1 = m.za[ia1] + m.zadd[im1];
                            let mut xdif = x1 - x2;
                            let mut ydif = y1 - y2;
                            let mut zdif = z1 - z2;
                            minimagec(&mut xdif, &mut ydif, &mut zdif);
                            let dsq = xdif * xdif + ydif * ydif + zdif * zdif;
                            let d12 = dsq.sqrt();
                            foundia1 = d12 < rmax1;
                        }
                    }
                    if foundia1 {
                        for &ia3 in atom3list.iter() {
                            let im3 = m.mref[ia3];
                            if ((s2 == 1 && im2 == im3) || (s2 == 0 && im2 != im3))
                                && ((s3 == 1 && im1 == im3) || (s3 == 0 && im1 != im3))
                                && ((ltyp3 == ltyp1 && ia3 > ia1) || (ltyp3 != ltyp1 && ia3 != ia1))
                            {
                                let x3 = m.xa[ia3] + m.xadd[im3];
                                let y3 = m.ya[ia3] + m.yadd[im3];
                                let z3 = m.za[ia3] + m.zadd[im3];
                                let mut xdif = x3 - x2;
                                let mut ydif = y3 - y2;
                                let mut zdif = z3 - z2;
                                minimagec(&mut xdif, &mut ydif, &mut zdif);
                                let dsq = xdif * xdif + ydif * ydif + zdif * zdif;
                                let d32 = dsq.sqrt();
                                let foundia3 = d32 < rmax2;
                                if foundia3 {
                                    ntriplets += 1;
                                    list.push([ia1, ia2, ia3]);
                                    tb.ntripletsinmol[(im2, id)] += 1;
                                    if tb.ntripletsinmol[(im2, id)] == 1 {
                                        tb.firsttripletinmol[(im2, id)] = ntriplets;
                                    }
                                    for imn in [im1, im3] {
                                        if imn != im2 {
                                            let mut it = 0;
                                            let mut found = false;
                                            while it < tb.nthreebneighbmol[im2] && !found {
                                                it += 1;
                                                found = imn == tb.threebneighbmol[(it, im2)];
                                            }
                                            if !found {
                                                tb.nthreebneighbmol[im2] += 1;
                                                if tb.nthreebneighbmol[im2] > tb.nthreebneighbmolmax {
                                                    tb.nthreebneighbmolmax += 4;
                                                    tb.threebneighbmol.grow(tb.nthreebneighbmolmax as usize, nmol as usize);
                                                }
                                                tb.threebneighbmol[(tb.nthreebneighbmol[im2], im2)] = imn;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            msgs.push(format!(" find_triplets> group {} : {} triplets", id, ntriplets));
        }
        for s in msgs {
            self.log(&s);
        }
        self.get_triplet_pot_cons();
    }

    /// `subroutine get_triplet_pot_cons()`: (only printed by the Fortran).
    pub fn get_triplet_pot_cons(&mut self) {
        let pi = PI4ATAN;
        let tb = &mut self.tb;
        if tb.nthreeb <= 0 {
            return;
        }
        tb.consthreeb.grow(tb.nthreeb as usize);
        for id in 1..=tb.nthreeb {
            let mut sumangle: f32 = 0.0;
            let mut sumlorentz: f32 = 0.0;
            let nangle: i32 = 1001;
            for ia in 1..=nangle {
                let angle = (ia - 1) as f32 * pi / (nangle - 1) as f32;
                let sineangle = angle.sin();
                sumangle += sineangle;
                for ib in 1..=tb.nthreebangs[id] {
                    let mut arg = (angle - tb.threebangs[(ib, id)]) / tb.threebawidth[id];
                    arg = 0.5 * arg * arg;
                    if arg < 30.0 {
                        sumlorentz += tb.threebamps[(ib, id)] * (-arg).exp();
                    }
                }
            }
            tb.consthreeb[id] = sumlorentz / sumangle;
        }
        let mut lw = ListWriter::new();
        lw.chars("get_triplet_pot_cons> Potential constants: ");
        for id in 1..=self.tb.nthreeb {
            lw.real(self.tb.consthreeb[id]);
        }
        let s = lw.finish();
        self.log(&s);
    }

    /// `function en_mole_triplets(imol, ...)`: three-body energy of the
    /// triplets centred on molecule `imol`; also sets `en_mole_triplet_dist`.
    pub fn en_mole_triplets(&mut self, imol: i32) -> f32 {
        let m = &self.mol;
        let tb = &mut self.tb;
        if tb.nthreeb <= 0 {
            return 0.0;
        }
        let mut ensum: f32 = 0.0;
        for id in 1..=tb.nthreeb {
            tb.en_mole_triplet_dist[id] = 0.0;
            if tb.threebamp[id].abs() > 0.0 {
                let mut entrisum: f32 = 0.0;
                let itf = tb.firsttripletinmol[(imol, id)];
                let itl = itf + tb.ntripletsinmol[(imol, id)] - 1;
                let list = &tb.tripletlist[(id - 1) as usize];
                let nangs = tb.nthreebangs[id];
                let awidth = tb.threebawidth[id];
                for it in itf..=itl {
                    let [ia1, ia2, ia3] = list[(it - 1) as usize];
                    let im1 = m.mref[ia1];
                    let im2 = m.mref[ia2];
                    let im3 = m.mref[ia3];
                    let x1 = m.xa[ia1] + m.xadd[im1];
                    let x2 = m.xa[ia2] + m.xadd[im2];
                    let x3 = m.xa[ia3] + m.xadd[im3];
                    let y1 = m.ya[ia1] + m.yadd[im1];
                    let y2 = m.ya[ia2] + m.yadd[im2];
                    let y3 = m.ya[ia3] + m.yadd[im3];
                    let z1 = m.za[ia1] + m.zadd[im1];
                    let z2 = m.za[ia2] + m.zadd[im2];
                    let z3 = m.za[ia3] + m.zadd[im3];
                    let mut xdif = x1 - x2;
                    let mut ydif = y1 - y2;
                    let mut zdif = z1 - z2;
                    minimagec(&mut xdif, &mut ydif, &mut zdif);
                    let d12sq = xdif * xdif + ydif * ydif + zdif * zdif;
                    let d12 = d12sq.sqrt();
                    let mut xdif = x3 - x2;
                    let mut ydif = y3 - y2;
                    let mut zdif = z3 - z2;
                    minimagec(&mut xdif, &mut ydif, &mut zdif);
                    let d32sq = xdif * xdif + ydif * ydif + zdif * zdif;
                    let d32 = d32sq.sqrt();
                    let mut xdif = x1 - x3;
                    let mut ydif = y1 - y3;
                    let mut zdif = z1 - z3;
                    minimagec(&mut xdif, &mut ydif, &mut zdif);
                    let d13sq = xdif * xdif + ydif * ydif + zdif * zdif;
                    let mut cost = (d12sq + d32sq - d13sq) / (2.0 * d12 * d32);
                    cost = fmax(-1.0, fmin(1.0, cost));
                    let mut entrimin: f32 = 0.0;
                    let mut ib = 1;
                    let mut ibmin = ib;
                    let mut arg = cost - tb.threebangs[(ib, id)];
                    while ib < nangs {
                        ib += 1;
                        let argtest = cost - tb.threebangs[(ib, id)];
                        if argtest.abs() < arg.abs() {
                            arg = argtest;
                            ibmin = ib;
                        }
                    }
                    arg /= awidth;
                    let term = 0.5 * arg * arg;
                    entrimin += tb.threebamps[(ibmin, id)] * term;
                    let truncatefactor: f32 = 1.0;
                    entrisum += truncatefactor * entrimin;
                    tb.en_mole_triplet_dist[id] += truncatefactor * entrimin;
                }
                if tb.ntripletsinmol[(imol, id)] > 0 {
                    tb.en_mole_triplet_dist[id] = tb.threebamp[id] * tb.en_mole_triplet_dist[id];
                    ensum += tb.threebamp[id] * entrisum;
                }
            }
        }
        ensum
    }

    /// `function pr_mole_triplets(imol, ...)`: three-body virial of molecule `imol`.
    pub fn pr_mole_triplets(&self, imol: i32) -> f64 {
        let m = &self.mol;
        let tb = &self.tb;
        if tb.nthreeb <= 0 {
            return 0.0;
        }
        let mut prsum: f64 = 0.0;
        for id in 1..=tb.nthreeb {
            if tb.ntripletsinmol[(imol, id)] > 0 && tb.threebamp[id].abs() > 0.0 {
                let itf = tb.firsttripletinmol[(imol, id)];
                let itl = itf + tb.ntripletsinmol[(imol, id)] - 1;
                let list = &tb.tripletlist[(id - 1) as usize];
                for it in itf..=itl {
                    let [ia1, ia2, ia3] = list[(it - 1) as usize];
                    let im1 = m.mref[ia1];
                    let im2 = m.mref[ia2];
                    let im3 = m.mref[ia3];
                    let x1 = m.xa[ia1] + m.xadd[im1];
                    let x2 = m.xa[ia2] + m.xadd[im2];
                    let x3 = m.xa[ia3] + m.xadd[im3];
                    let y1 = m.ya[ia1] + m.yadd[im1];
                    let y2 = m.ya[ia2] + m.yadd[im2];
                    let y3 = m.ya[ia3] + m.yadd[im3];
                    let z1 = m.za[ia1] + m.zadd[im1];
                    let z2 = m.za[ia2] + m.zadd[im2];
                    let z3 = m.za[ia3] + m.zadd[im3];
                    let mut xdif = x1 - x2;
                    let mut ydif = y1 - y2;
                    let mut zdif = z1 - z2;
                    minimagec(&mut xdif, &mut ydif, &mut zdif);
                    let d12sq = xdif * xdif + ydif * ydif + zdif * zdif;
                    let d12 = d12sq.sqrt();
                    let mut xdif = x3 - x2;
                    let mut ydif = y3 - y2;
                    let mut zdif = z3 - z2;
                    minimagec(&mut xdif, &mut ydif, &mut zdif);
                    let d32sq = xdif * xdif + ydif * ydif + zdif * zdif;
                    let d32 = d32sq.sqrt();
                    if d12 < tb.threebrmax[(1, id)] && d32 < tb.threebrmax[(2, id)] {
                        let mut xdif = x1 - x3;
                        let mut ydif = y1 - y3;
                        let mut zdif = z1 - z3;
                        minimagec(&mut xdif, &mut ydif, &mut zdif);
                        let d13sq = xdif * xdif + ydif * ydif + zdif * zdif;
                        let cost = (d12sq + d32sq - d13sq) / (2.0 * d12 * d32);
                        let mut xcdif = m.xadd[im1] - m.xadd[im3];
                        let mut ycdif = m.yadd[im1] - m.yadd[im3];
                        let mut zcdif = m.zadd[im1] - m.zadd[im3];
                        minimagec(&mut xcdif, &mut ycdif, &mut zcdif);
                        let prfac: f32 = if d12 > 0.0 && d32 > 0.0 {
                            -(xcdif * xdif + ycdif * ydif + zcdif * zdif) / (d12 * d32)
                        } else {
                            0.0
                        };
                        let mut ib = 1;
                        let mut ibmin = ib;
                        let mut arg = cost - tb.threebangs[(ib, id)];
                        while ib < tb.nthreebangs[id] {
                            ib += 1;
                            let argtest = cost - tb.threebangs[(ib, id)];
                            if argtest.abs() < arg.abs() {
                                arg = argtest;
                                ibmin = ib;
                            }
                        }
                        let pot: f32 = prfac * tb.threebamps[(ibmin, id)] * arg / (tb.threebawidth[id] * tb.threebawidth[id]);
                        prsum += (tb.threebamp[id] * pot) as f64;
                    }
                }
            }
        }
        prsum
    }

    /// `subroutine save_threeb_energies(im)`.
    pub fn save_threeb_energies(&mut self, im: i32) {
        let tb = &mut self.tb;
        if tb.nthreeb > 0 {
            tb.ergold_triplet[im] = tb.ergnew_triplet[im];
            if tb.nthreebneighbmol[im] > 0 {
                for it in 1..=tb.nthreebneighbmol[im] {
                    let imref = tb.threebneighbmol[(it, im)];
                    tb.ergold_triplet[imref] = tb.ergnew_triplet[imref];
                }
            }
        }
    }

    /// The three-body part of `en_change(im)`: adds the change of the triplet
    /// energy of `im` and of its triplet neighbours to `fdiff` (fills `ergnew_triplet`).
    pub fn en_change_triplets(&mut self, im: i32, fdiff: &mut f32) {
        let e = self.en_mole_triplets(im);
        self.tb.ergnew_triplet[im] = e;
        *fdiff = *fdiff + self.tb.ergnew_triplet[im] - self.tb.ergold_triplet[im];
        if self.tb.nthreebneighbmol[im] > 0 {
            for it in 1..=self.tb.nthreebneighbmol[im] {
                let imref = self.tb.threebneighbmol[(it, im)];
                let e = self.en_mole_triplets(imref);
                self.tb.ergnew_triplet[imref] = e;
                *fdiff = *fdiff + self.tb.ergnew_triplet[imref] - self.tb.ergold_triplet[imref];
            }
        }
    }
}

/// Fortran `max(a, b)` for reals as compiled by gfortran (`a > b ? a : b` order kept).
#[inline]
fn fmax(a: f32, b: f32) -> f32 {
    if b > a {
        b
    } else {
        a
    }
}

#[inline]
fn fmin(a: f32, b: f32) -> f32 {
    if b < a {
        b
    } else {
        a
    }
}
