//! Coordination-number constraints — translation of modules
//! `coordconstraint_mod` and `coordconstraints_mod` (`coordconstraint_mod.f90`),
//! plus `en_atom_inter_scale_coord`, `en_gauss`, `en_harmonic` and
//! `en_atom_inter_scale_coordextra` (`molfit.f`).
//!
//! A `.pcof` line
//! ```text
//! coordcon atc nc cctemp ccdistort ccgaussweight ccgausswidth  ats1 ds1 dsmax1 [ats2 ds2 dsmax2 ...]
//!          [angles tbweight tbdistort tbgaussweight tbgausswidth ang1 [ang2 ...]]
//!          [anglesintra wt distort gaussweight gausswidth ats ang [ats ang ...]]
//! ```
//! asks that every atom of type `atc` has `nc` neighbours of the types `ats`
//! (on other molecules).  `find_constraints` (called by `fmolec1` and `press`)
//! assigns explicit "bonds" between atoms; each bond then contributes a
//! Morse-like energy (`en_harmonic`, or `en_gauss` for the derived three-body
//! bonds with negative `tbondcc`) to the molecule energy, and the pair is
//! evaluated with the reference potential without the extra Gaussians and
//! repulsive core (`epot`) instead of `ecoul`.
//!
//! The arrays of each constraint keep their contents between calls (as the
//! Fortran allocatables do); `realloc1`/`realloc2` reproduce
//! `reallocate1d_*`/`reallocate2d_*` exactly (shrinking discards, new
//! elements are zero here — uninitialised in the Fortran).

use crate::error::{EpsrError, Result};
use crate::fortran::fmath::nint;
use crate::fortran::format::{fmt, FV};
use crate::fortran::listread::{read_int, read_real};
use crate::fortran::output::OutFile;
use crate::fortran::text::{fstr, parse_n};
use crate::fortran::{F1, F2};
use crate::neighbour::minimagec;
use crate::state::{Mol, State};

/// One `coordconstraint` object.
#[derive(Clone, Debug, Default)]
pub struct CoordConstraint {
    pub atc: String,
    pub ats: F1<String>,
    pub nc: i32,
    pub alternatenc: i32,
    pub cctemp: f32,
    pub ccgaussweight: f32,
    pub ccgausswidth: f32,
    pub ccdistort: f32,
    pub ns: i32,
    pub atclistsize: i32,
    pub atclist: F1<i32>,
    pub nbondcc: F1<i32>,
    pub nearnbondcc: F1<i32>,
    pub expecnbondcc: F1<i32>,
    pub iabondcc: F2<i32>,
    pub wtcc: F2<i32>,
    pub ds: F1<f32>,
    pub dsmax: F1<f32>,
    pub tbondcc: F2<f32>,
    pub adistcc: F2<f32>,
    pub hastb: bool,
    pub hastbintra: bool,
    pub tbweight: f32,
    pub tbgaussweight: f32,
    pub tbintraweight: f32,
    pub tbintragaussweight: f32,
    pub tbdistort: f32,
    pub tbgausswidth: f32,
    pub tbintradistort: f32,
    pub tbintragausswidth: f32,
    pub ntbangles: i32,
    pub ntbintraangles: i32,
    pub tbangles: F1<f32>,
    pub tbcosines: F1<f32>,
    pub tbintraangles: F1<f32>,
    pub tbintracosines: F1<f32>,
    pub tbintraats: F1<String>,
}

/// Module `coordconstraints_mod` variables (+ `fcoordconstraints` of `fmolec1`).
#[derive(Clone, Debug, Default)]
pub struct CoordConstraints {
    /// `coordconstraints(1..ncoordconstraints)` (stored 0-based).
    pub coordconstraints: Vec<CoordConstraint>,
    pub ncoordconstraints: i32,
    pub ergold_coord: F1<f32>,
    pub ergnew_coord: F1<f32>,
    /// Average coordination-constraint energy per molecule from the last `fmolec1`.
    pub fcoordconstraints: f32,
}

/// Reading memory the Fortran never wrote (an uninitialised allocatable
/// element): the result would be undefined, so stop loudly.
fn ub(what: &str) -> ! {
    panic!("not supported by the Rust port: coordination constraint reads uninitialised data ({what})")
}

/// `reallocate1d_*(array, n)`.
fn realloc1<T: Clone + Default>(a: &mut F1<T>, n: usize) {
    if a.data.len() == n {
        return;
    }
    a.data.resize(n, T::default());
}

/// `reallocate2d_*(array, n1, n2)`: keep `(j,i)` for `j <= min(n1)`, `i <= min(n2)`.
fn realloc2<T: Clone + Default>(a: &mut F2<T>, n1: usize, n2: usize) {
    if a.n1 == n1 && a.n2 == n2 {
        return;
    }
    let mut nd = vec![T::default(); n1 * n2];
    for i in 0..a.n2.min(n2) {
        for j in 0..a.n1.min(n1) {
            nd[j + i * n1] = a.data[j + i * a.n1].clone();
        }
    }
    a.n1 = n1;
    a.n2 = n2;
    a.data = nd;
}

/// `function find_type(alabel, ntype, atype)`.
fn find_type(alabel: &str, m: &Mol) -> i32 {
    let mut found = false;
    let mut i = 0;
    while i < m.ntype && !found {
        i += 1;
        found = alabel == m.type_[i];
    }
    if found {
        i
    } else {
        0
    }
}

fn rd_i(s: &str) -> Result<i32> {
    read_int(s).map_err(|_| EpsrError::Input(format!("coordcon: cannot read integer from '{s}'")))
}
fn rd_r(s: &str) -> Result<f32> {
    read_real(s).map_err(|_| EpsrError::Input(format!("coordcon: cannot read real from '{s}'")))
}

/// `function en_gauss(ddif,wtm,wtms,distort,gausswt,gausswidth)` (`molfit.f`).
pub fn en_gauss(ddif: f32, wtm: f32, wtms: f32, distort: f32, gausswt: f32, gausswidth: f32) -> f32 {
    let arg = wtm * ddif;
    let mut args = arg / distort;
    let mut en: f32 = if args > 30.0 {
        (args / 30.0) * 30.0f32.exp() - 1.0
    } else if args < -30.0 {
        -1.0
    } else {
        args.exp() - 1.0
    };
    en = wtms * (en * en - 1.0);
    if gausswt != 0.0 && gausswidth > 0.0 {
        args = arg / gausswidth;
        args = 0.5 * args * args;
        if args < 30.0 {
            en -= gausswt * (-args).exp();
        }
    }
    en
}

/// `function en_harmonic(ddif,wtm,wtms,distort,gausswt,gausswidth)` (`molfit.f`).
pub fn en_harmonic(ddif: f32, wtm: f32, wtms: f32, distort: f32, gausswt: f32, gausswidth: f32) -> f32 {
    let arg = wtm * ddif;
    let mut args = arg / distort;
    let mut en: f32 = if args > 30.0 {
        (args / 30.0) * 30.0f32.exp() - 1.0
    } else if args < -30.0 {
        -1.0
    } else {
        args.exp() - 1.0
    };
    en = wtms * (en * en - 1.0);
    if gausswt != 0.0 && gausswidth > 0.0 {
        args = arg / gausswidth;
        if args.abs() < 1.0 {
            en += gausswt * args * args;
        } else {
            en += gausswt * args.abs();
        }
    }
    en
}

/// `function en_atom_inter_scale_coord(dtest,tb,wtm,wtms,distort,gausswt,gausswidth)`.
pub fn en_atom_inter_scale_coord(dtest: f32, tb: f32, wtm: f32, wtms: f32, distort: f32, gausswt: f32, gausswidth: f32) -> f32 {
    let ddif = tb.abs() - dtest;
    if tb < 0.0 {
        en_gauss(ddif, wtm.sqrt(), wtms, distort, gausswt, gausswidth)
    } else {
        en_harmonic(ddif, wtm.sqrt(), wtms, distort, gausswt, gausswidth)
    }
}

/// `function en_atom_inter_scale_coordextra(x1,...,scalefac)` (`molfit.f`):
/// as `en_atom_inter_scale` but with the reference potential `epot`
/// (no additional Gaussians / repulsive core) for coordination-bonded pairs.
#[inline]
pub fn en_atom_inter_scale_coordextra(
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
        sume = sume + p.xavs[(iz, icref)] + p.xavs1[(iz, icref)] * dtdif;
    }
    dtest *= scalefac;
    let rat = dtest / stepcl;
    let iz = rat as i32 + 2;
    if iz <= npt2 {
        let dtdif = dtest - xstar[iz];
        sume = sume + p.epot[(iz, icref)] + p.epot1[(iz, icref)] * dtdif;
    }
    sume
}

impl CoordConstraint {
    /// `subroutine init_constraint(this, sometext)`; `sometext` is the
    /// `character*256` text after the word `coordcon`.
    pub fn init_constraint(sometext: &str, log: &dyn Fn(&str)) -> Result<CoordConstraint> {
        let mut this = CoordConstraint::default();
        let w = parse_n(sometext, 10, 256);
        let nwords = w.nwords();
        if nwords < 6 {
            return Err(EpsrError::Input("coordcon line with fewer than 6 values".into()));
        }
        let wd = |i: i32| -> Result<&str> {
            if i >= 1 && i <= nwords {
                Ok(w.word(sometext, i))
            } else {
                Err(EpsrError::Input("coordcon line ends prematurely".into()))
            }
        };
        log(&format!(" init_constraint> {} {}", nwords, w.span(sometext, 1, nwords)));
        this.atc = fstr(wd(1)?, 3);
        this.nc = rd_i(wd(2)?)?;
        this.cctemp = rd_r(wd(3)?)?;
        this.ccdistort = rd_r(wd(4)?)?;
        this.ccgaussweight = rd_r(wd(5)?)?;
        this.ccgausswidth = rd_r(wd(6)?)?;
        let mut iread = 6;
        this.ns = 0;
        this.ntbangles = 0;
        this.ntbintraangles = 0;
        this.tbweight = 0.0;
        this.tbdistort = 1.0;
        this.tbgaussweight = 0.0;
        this.tbgausswidth = 1.0;
        this.tbintraweight = 0.0;
        this.tbintradistort = 1.0;
        this.tbintragaussweight = 0.0;
        this.tbintragausswidth = 1.0;
        this.hastb = false;
        this.hastbintra = false;
        while iread < nwords {
            iread += 1;
            if !this.hastb && !this.hastbintra {
                this.hastb = wd(iread)? == "angles";
                this.hastbintra = wd(iread)? == "anglesintra";
                if !this.hastb && !this.hastbintra {
                    this.ns += 1;
                    let ns = this.ns as usize;
                    realloc1(&mut this.ats, ns);
                    realloc1(&mut this.ds, ns);
                    realloc1(&mut this.dsmax, ns);
                    this.ats[ns] = fstr(wd(iread)?, 3);
                    iread += 1;
                    this.ds[ns] = rd_r(wd(iread)?)?;
                    iread += 1;
                    this.dsmax[ns] = rd_r(wd(iread)?)?;
                }
            } else if this.hastb && !this.hastbintra {
                this.tbweight = rd_r(wd(iread)?)?;
                iread += 1;
                this.tbdistort = rd_r(wd(iread)?)?;
                iread += 1;
                this.tbgaussweight = rd_r(wd(iread)?)?;
                iread += 1;
                this.tbgausswidth = rd_r(wd(iread)?)?;
                while iread < nwords && !this.hastbintra {
                    iread += 1;
                    this.hastbintra = wd(iread)? == "anglesintra";
                    if !this.hastbintra {
                        this.ntbangles += 1;
                        realloc1(&mut this.tbangles, this.ntbangles as usize);
                        this.tbangles[this.ntbangles] = rd_r(wd(iread)?)?;
                    }
                }
            } else if this.hastbintra {
                this.tbintraweight = rd_r(wd(iread)?)?;
                iread += 1;
                this.tbintradistort = rd_r(wd(iread)?)?;
                iread += 1;
                this.tbintragaussweight = rd_r(wd(iread)?)?;
                iread += 1;
                this.tbintragausswidth = rd_r(wd(iread)?)?;
                while iread < nwords - 1 {
                    iread += 1;
                    this.ntbintraangles += 1;
                    let n = this.ntbintraangles as usize;
                    realloc1(&mut this.tbintraats, n);
                    realloc1(&mut this.tbintraangles, n);
                    this.tbintraats[n] = fstr(wd(iread)?, 3);
                    iread += 1;
                    this.tbintraangles[n] = rd_r(wd(iread)?)?;
                }
            }
        }
        log(&format!(" Coordination constraint {} with coordination number {} created.", this.atc, this.nc));
        let piconv: f32 = 4.0 * 1.0f32.atan() / 180.0;
        if this.hastb {
            realloc1(&mut this.tbcosines, this.ntbangles as usize);
            for i in 1..=this.ntbangles {
                this.tbcosines[i] = (piconv * this.tbangles[i]).cos();
            }
        }
        if this.hastbintra {
            realloc1(&mut this.tbintracosines, this.ntbintraangles as usize);
            for i in 1..=this.ntbintraangles {
                this.tbintracosines[i] = (piconv * this.tbintraangles[i]).cos();
            }
        }
        Ok(this)
    }

    /// Grow the four bond arrays to `(nbondccsize, natom)`.
    fn realloc_bonds(&mut self, nbondccsize: i32, natom: i32) {
        let (a, b) = (nbondccsize as usize, natom as usize);
        realloc2(&mut self.iabondcc, a, b);
        realloc2(&mut self.wtcc, a, b);
        realloc2(&mut self.tbondcc, a, b);
        realloc2(&mut self.adistcc, a, b);
    }

    /// Copy bond `src` of atom `ia` to slot `dst`.
    fn copy_bond(&mut self, dst: i32, src: i32, ia: i32) {
        self.iabondcc[(dst, ia)] = self.iabondcc[(src, ia)];
        self.wtcc[(dst, ia)] = self.wtcc[(src, ia)];
        self.tbondcc[(dst, ia)] = self.tbondcc[(src, ia)];
        self.adistcc[(dst, ia)] = self.adistcc[(src, ia)];
    }

    fn set_bond(&mut self, k: i32, ia: i32, ja: i32, wt: i32, tb: f32, ad: f32) {
        self.iabondcc[(k, ia)] = ja;
        self.wtcc[(k, ia)] = wt;
        self.tbondcc[(k, ia)] = tb;
        self.adistcc[(k, ia)] = ad;
    }

    /// `iabondcc(k, ia)` for a read the Fortran may do on a never-written slot.
    fn bond_atom(&self, k: i32, ia: i32) -> i32 {
        let v = self.iabondcc[(k, ia)];
        if v < 1 {
            ub("iabondcc");
        }
        v
    }

    /// The "Neighbouring atoms in same molecule" diagnostic (log only).
    fn check_same(&self, m: &Mol, pass: i32, log: &dyn Fn(&str)) {
        let mut sumsame = 0;
        for ia in 1..=m.natom {
            let nb = self.nbondcc[ia];
            if nb > 0 {
                for j in 1..=nb {
                    if m.mref[ia] == m.mref[self.iabondcc[(j, ia)]] {
                        sumsame += 1;
                    }
                }
                for j in 1..nb {
                    let jm = m.mref[self.iabondcc[(j, ia)]];
                    for k in j + 1..=nb {
                        if m.mref[self.iabondcc[(k, ia)]] == jm {
                            sumsame += 1;
                        }
                    }
                }
            }
        }
        log(&format!(
            " find_constraints - pass {pass}> {} Molecules listed as neighbours to themselves (should be zero) {sumsame}",
            self.atc
        ));
    }

    /// `subroutine find_constraints(this, nmol, natom, naf, nal, xadd, ..., cfrac, dmax)`.
    pub fn find_constraints(&mut self, m: &Mol, log: &dyn Fn(&str)) {
        let natom = m.natom;
        let nmol = m.nmol;
        let atype = |ia: i32| -> &str { &m.type_[m.itype[ia]] };
        let ns = self.ns;
        let nc = self.nc;
        // atomic fraction of the central and the secondary atom types
        let mut cacfrac: f32 = 0.0;
        let it = find_type(&self.atc, m);
        if it > 0 {
            cacfrac = m.cfrac[it];
        }
        let mut sacfrac: f32 = 0.0;
        for ia in 1..=ns {
            let it = find_type(&self.ats[ia], m);
            if it > 0 {
                sacfrac += m.cfrac[it];
            }
        }
        if sacfrac > 0.0 {
            self.alternatenc = nint(nc as f32 * cacfrac / sacfrac);
        }
        let mut nbondccsize = nc.max(self.alternatenc);
        realloc1(&mut self.nbondcc, natom as usize);
        realloc1(&mut self.nearnbondcc, natom as usize);
        realloc1(&mut self.expecnbondcc, natom as usize);
        self.realloc_bonds(nbondccsize, natom);
        self.nbondcc.fill(0);
        self.nearnbondcc.fill(0);
        self.expecnbondcc.fill(0);
        // list of atoms of type atc
        let mut ic = 0;
        for ia in 1..=natom {
            if atype(ia) == self.atc {
                ic += 1;
            }
        }
        self.atclistsize = ic;
        realloc1(&mut self.atclist, self.atclistsize as usize);
        ic = 0;
        for ia in 1..=natom {
            if atype(ia) == self.atc {
                ic += 1;
                self.atclist[ic] = ia;
                self.expecnbondcc[ia] = nc;
            } else {
                let mut it1 = 1;
                let mut found = atype(ia) == self.ats[it1];
                while it1 < ns && !found {
                    it1 += 1;
                    found = atype(ia) == self.ats[it1];
                }
                if found {
                    self.expecnbondcc[ia] = self.alternatenc;
                }
            }
        }
        // ---- pass 1: the nc nearest neighbours (on distinct molecules) of each atc atom
        // NB: ishortest/dshortest are not reset between central atoms (Fortran locals).
        let mut ishortest: F1<i32> = F1::new(nc.max(0) as usize);
        let mut dshortest: F1<f32> = F1::new(nc.max(0) as usize);
        let mut sumoutside = 0;
        for ic in 1..=self.atclistsize {
            let ia = self.atclist[ic];
            let im = m.mref[ia];
            let mut nshortest = 0;
            for ja in 1..=natom {
                let jm = m.mref[ja];
                let mut found = im == jm;
                let mut it = 0;
                while it < nshortest && !found {
                    it += 1;
                    found = jm == m.mref[ishortest[it]];
                }
                if !found {
                    it = 0;
                    while it < ns && !found {
                        it += 1;
                        found = self.ats[it] == atype(ja);
                    }
                    if found {
                        let mut xdif = (m.xadd[im] + m.xa[ia]) - (m.xadd[jm] + m.xa[ja]);
                        let mut ydif = (m.yadd[im] + m.ya[ia]) - (m.yadd[jm] + m.ya[ja]);
                        let mut zdif = (m.zadd[im] + m.za[ia]) - (m.zadd[jm] + m.za[ja]);
                        minimagec(&mut xdif, &mut ydif, &mut zdif);
                        let dsq = xdif * xdif + ydif * ydif + zdif * zdif;
                        let d12 = dsq.sqrt();
                        it = 0;
                        found = false;
                        while it < nshortest && !found {
                            it += 1;
                            found = d12 < dshortest[it];
                        }
                        if found {
                            if nshortest < nc {
                                nshortest += 1;
                            }
                            let mut it2 = nshortest;
                            while it2 > it {
                                let it1 = it2 - 1;
                                ishortest[it2] = ishortest[it1];
                                dshortest[it2] = dshortest[it1];
                                it2 = it1;
                            }
                            ishortest[it] = ja;
                            dshortest[it] = d12;
                        } else if nshortest < nc {
                            nshortest += 1;
                            ishortest[nshortest] = ja;
                            dshortest[nshortest] = d12;
                        }
                    }
                }
            }
            for it in 1..=nc {
                let ja = ishortest[it];
                if ja < 1 {
                    ub("ishortest");
                }
                let mut it1 = 0;
                let mut found = false;
                while it1 < ns && !found {
                    it1 += 1;
                    found = atype(ja) == self.ats[it1];
                }
                if found {
                    if dshortest[it] > self.dsmax[it1] {
                        sumoutside += 1;
                    }
                    let dsfound = self.ds[it1];
                    self.nbondcc[ia] += 1;
                    if self.nbondcc[ia] > nbondccsize {
                        nbondccsize = self.nbondcc[ia];
                        self.realloc_bonds(nbondccsize, natom);
                    }
                    let k = self.nbondcc[ia];
                    self.set_bond(k, ia, ja, 1, dsfound, dshortest[it]);
                    if m.itype[ia] != m.itype[ja] {
                        // insert ia in the list of ja, ordered by distance
                        let mut is = 0;
                        let mut found = false;
                        while is < self.nbondcc[ja] && !found {
                            is += 1;
                            found = dshortest[it] < self.adistcc[(is, ja)];
                        }
                        self.nbondcc[ja] += 1;
                        if self.nbondcc[ja] > nbondccsize {
                            nbondccsize = self.nbondcc[ja];
                            self.realloc_bonds(nbondccsize, natom);
                        }
                        if found {
                            let mut is2 = self.nbondcc[ja];
                            while is2 > is {
                                let is1 = is2 - 1;
                                self.copy_bond(is2, is1, ja);
                                is2 = is1;
                            }
                            self.set_bond(is, ja, ia, 1, dsfound, dshortest[it]);
                        } else {
                            let k = self.nbondcc[ja];
                            self.set_bond(k, ja, ia, 1, dsfound, dshortest[it]);
                        }
                    }
                }
            }
        }
        self.check_same(m, 1, log);
        // ---- pass 2: remove the bonds of over-bonded atoms
        for ia in 1..=natom {
            let mut it = self.nbondcc[ia];
            while it > self.expecnbondcc[ia] {
                let ja = self.iabondcc[(it, ia)];
                let mut it1 = self.nbondcc[ja];
                if it1 < 1 {
                    ub("iabondcc(0,ja)");
                }
                let mut found = ia == self.iabondcc[(it1, ja)];
                // NB: the search stops at it1 = 2 (slot 1 is only tested when it is the last one).
                while it1 > 2 && !found {
                    it1 -= 1;
                    found = ia == self.iabondcc[(it1, ja)];
                }
                if found {
                    for it2 in it1..self.nbondcc[ja] {
                        self.copy_bond(it2, it2 + 1, ja);
                    }
                    self.nbondcc[ja] -= 1;
                }
                it -= 1;
            }
            self.nbondcc[ia] = it;
        }
        self.check_same(m, 2, log);
        // ---- pass 3: complete under-bonded atc atoms with the nearest under-bonded partner
        for ic in 1..=self.atclistsize {
            let ia = self.atclist[ic];
            let mut it = self.nbondcc[ia];
            while it < self.expecnbondcc[ia] {
                it += 1;
                let mut ja = 0;
                let mut found = false;
                let mut it1 = 0;
                // shared test: ja under-bonded, on another molecule, of a secondary type,
                // and its molecule not yet a neighbour of ia (NB: slot 1 of ia is read
                // even when ia has no bonds — the stale value of an earlier call)
                let candidate = |this: &CoordConstraint, ja: i32, it1: &mut i32| -> bool {
                    let mut found = this.nbondcc[ja] < this.expecnbondcc[ja] && m.mref[ja] != m.mref[ia];
                    if found {
                        *it1 = 1;
                        found = atype(ja) == this.ats[*it1];
                        while *it1 < ns && !found {
                            *it1 += 1;
                            found = atype(ja) == this.ats[*it1];
                        }
                        if found {
                            let mut it2 = 1;
                            found = m.mref[ja] != m.mref[this.bond_atom(it2, ia)];
                            while it2 < this.nbondcc[ia] && found {
                                it2 += 1;
                                found = m.mref[ja] != m.mref[this.iabondcc[(it2, ia)]];
                            }
                        }
                    }
                    found
                };
                while ja < natom && !found {
                    ja += 1;
                    found = candidate(self, ja, &mut it1);
                }
                let mut it1save = it1;
                if found {
                    let dist = |ja: i32| -> f32 {
                        let (mi, mj) = (m.mref[ia], m.mref[ja]);
                        let mut xdif = m.xadd[mi] + m.xa[ia] - m.xadd[mj] - m.xa[ja];
                        let mut ydif = m.yadd[mi] + m.ya[ia] - m.yadd[mj] - m.ya[ja];
                        let mut zdif = m.zadd[mi] + m.za[ia] - m.zadd[mj] - m.za[ja];
                        minimagec(&mut xdif, &mut ydif, &mut zdif);
                        xdif * xdif + ydif * ydif + zdif * zdif
                    };
                    let mut dsqmin = dist(ja);
                    let mut jamin = ja;
                    while ja < natom {
                        ja += 1;
                        let found = candidate(self, ja, &mut it1);
                        if found {
                            let dsq = dist(ja);
                            if dsq < dsqmin {
                                jamin = ja;
                                dsqmin = dsq;
                                it1save = it1;
                            }
                        }
                    }
                    let dsv = self.ds[it1save];
                    self.nbondcc[jamin] += 1;
                    let k = self.nbondcc[jamin];
                    self.set_bond(k, jamin, ia, 1, dsv, dsqmin.sqrt());
                    let ad = self.adistcc[(k, jamin)];
                    self.nbondcc[ia] += 1;
                    let k2 = self.nbondcc[ia];
                    self.set_bond(k2, ia, jamin, 1, dsv, ad);
                }
            }
        }
        self.check_same(m, 3, log);
        self.nearnbondcc = self.nbondcc.clone();
        let mut nt = 0;
        let mut sumratio: f32 = 0.0;
        for ia in 1..=natom {
            nt += 1;
            if self.expecnbondcc[ia] == self.nbondcc[ia] {
                sumratio += 1.0;
            }
        }
        if nt > 0 {
            sumratio = 100.0 * sumratio / nt as f32;
        }
        log(&format!(
            " find_constraints: {} Fraction outside specified distance(s) = {}",
            self.atc,
            sumoutside as f32 / (nc * self.atclistsize) as f32
        ));
        log(&format!(" find_constraints: {} Percentage constraints correct = {}", self.atc, sumratio));
        // ---- three-body ("angles") bonds between pairs of neighbours of each atc atom
        if self.hastb {
            for ic in 1..=self.atclistsize {
                let ia = self.atclist[ic];
                let mut i = 0;
                while i < self.nearnbondcc[ia] - 1 {
                    i += 1;
                    let ja1 = self.iabondcc[(i, ia)];
                    let jm1 = m.mref[ja1];
                    let tb1 = self.adistcc[(i, ia)];
                    let mut j = i;
                    while j < self.nearnbondcc[ia] {
                        j += 1;
                        let ja2 = self.iabondcc[(j, ia)];
                        let jm2 = m.mref[ja2];
                        // NB: the Fortran uses adistcc(i,ia) (not j) for the second distance.
                        let tb2 = self.adistcc[(i, ia)];
                        let mut xdif = (m.xadd[jm1] + m.xa[ja1]) - (m.xadd[jm2] + m.xa[ja2]);
                        let mut ydif = (m.yadd[jm1] + m.ya[ja1]) - (m.yadd[jm2] + m.ya[ja2]);
                        let mut zdif = (m.zadd[jm1] + m.za[ja1]) - (m.zadd[jm2] + m.za[ja2]);
                        minimagec(&mut xdif, &mut ydif, &mut zdif);
                        let dist12 = xdif * xdif + ydif * ydif + zdif * zdif;
                        let mut dsq = tb1 * tb1 + tb2 * tb2 - 2.0 * tb1 * tb2 * self.tbcosines[1];
                        let mut tb = dsq.sqrt();
                        let mut dsave = (dist12 - dsq).abs();
                        let mut itb = 1;
                        while itb < self.ntbangles {
                            itb += 1;
                            dsq = tb1 * tb1 + tb2 * tb2 - 2.0 * tb1 * tb2 * self.tbcosines[itb];
                            let dtest = (dist12 - dsq).abs();
                            if dtest < dsave {
                                dsave = dtest;
                                tb = dsq.sqrt();
                            }
                        }
                        let mut itb = self.nearnbondcc[ja1];
                        let mut found = false;
                        while itb < self.nbondcc[ja1] && !found {
                            itb += 1;
                            found = ja2 == self.iabondcc[(itb, ja1)];
                        }
                        if found {
                            let a = self.tbondcc[(itb, ja1)].abs();
                            tb = if tb < a { tb } else { a };
                            self.tbondcc[(itb, ja1)] = -tb;
                            let mut itb = self.nearnbondcc[ja2];
                            let mut found = false;
                            while itb < self.nbondcc[ja2] && !found {
                                itb += 1;
                                found = ja1 == self.iabondcc[(itb, ja2)];
                            }
                            if found {
                                self.tbondcc[(itb, ja2)] = -tb;
                            }
                        } else {
                            let mut itb = 1;
                            found = m.mref[self.bond_atom(itb, ja1)] != m.mref[ja2];
                            while itb < self.nbondcc[ja1] && found {
                                itb += 1;
                                found = m.mref[self.iabondcc[(itb, ja1)]] != m.mref[ja2];
                            }
                            itb = 1;
                            found = found && m.mref[self.bond_atom(itb, ja2)] != m.mref[ja1];
                            while itb < self.nbondcc[ja2] && found {
                                itb += 1;
                                found = found && m.mref[self.iabondcc[(itb, ja2)]] != m.mref[ja1];
                            }
                            if found {
                                self.nbondcc[ja1] += 1;
                                self.nbondcc[ja2] += 1;
                                nbondccsize = nbondccsize.max(self.nbondcc[ja1]).max(self.nbondcc[ja2]);
                                self.realloc_bonds(nbondccsize, natom);
                                let k1 = self.nbondcc[ja1];
                                let k2 = self.nbondcc[ja2];
                                self.set_bond(k1, ja1, ja2, 2, -tb, dist12.sqrt());
                                let ad = self.adistcc[(k1, ja1)];
                                self.set_bond(k2, ja2, ja1, 2, -tb, ad);
                            }
                        }
                    }
                }
            }
            self.check_same(m, 4, log);
        }
        // ---- intramolecular three-body ("anglesintra") bonds
        if self.hastbintra {
            let mut ndonors: F1<i32> = F1::new(nmol as usize);
            let mut nacceptors: F1<i32> = F1::new(nmol as usize);
            let mut ndonorsmax: F1<i32> = F1::new(nmol as usize);
            for itb in 1..=self.ntbintraangles {
                ndonors.fill(0);
                nacceptors.fill(0);
                // donors(k, im) / acceptors(k, im), k = 1..ndonors(im) / nacceptors(im)
                let mut donors: Vec<Vec<i32>> = vec![Vec::new(); nmol as usize + 1];
                let mut acceptors: Vec<Vec<i32>> = vec![Vec::new(); nmol as usize + 1];
                let it1t = find_type(&self.tbintraats[itb], m);
                for im in 1..=nmol {
                    ndonorsmax[im] = 0;
                    for ia in m.naf[im]..=m.nal[im] {
                        if m.itype[ia] == it1t {
                            ndonorsmax[im] += 1;
                        }
                    }
                }
                // the two passes are identical apart from the form of the entry test
                for _pass in 1..=2 {
                    for ic in 1..=self.atclistsize {
                        let ia = self.atclist[ic];
                        let im = m.mref[ia];
                        for ja1 in m.naf[im]..=m.nal[im] {
                            let jm1 = m.mref[ja1];
                            let found = nacceptors[jm1] < ndonorsmax[jm1] && atype(ja1) == self.tbintraats[itb] && ja1 != ia;
                            if !found {
                                continue;
                            }
                            let mut xdif = m.xa[ja1] - m.xa[ia];
                            let mut ydif = m.ya[ja1] - m.ya[ia];
                            let mut zdif = m.za[ja1] - m.za[ia];
                            minimagec(&mut xdif, &mut ydif, &mut zdif);
                            let dsq = xdif * xdif + ydif * ydif + zdif * zdif;
                            let tb1 = dsq.sqrt();
                            let mut inearbond = 0;
                            let mut dsave: f32 = 0.0;
                            let mut jamin = 0;
                            let mut dsqmin: f32 = 0.0;
                            let mut dsqtbmin: f32 = 0.0;
                            while inearbond < self.nearnbondcc[ia] {
                                inearbond += 1;
                                let ja2 = self.iabondcc[(inearbond, ia)];
                                let jm2 = m.mref[ja2];
                                let mut found = ndonors[jm2] < ndonorsmax[jm2];
                                let mut it1 = 0;
                                while it1 < self.nbondcc[ja1] && found {
                                    it1 += 1;
                                    found = ja2 != self.iabondcc[(it1, ja1)];
                                }
                                it1 = 0;
                                while it1 < self.nbondcc[ja2] && found {
                                    it1 += 1;
                                    found = ja1 != self.iabondcc[(it1, ja2)];
                                }
                                let mut ic1 = 0;
                                while ic1 < ndonors[jm2] && found {
                                    ic1 += 1;
                                    found = jm1 != donors[jm2 as usize][(ic1 - 1) as usize];
                                }
                                let mut ic2 = 0;
                                while ic2 < nacceptors[jm2] && found {
                                    ic2 += 1;
                                    found = jm1 != acceptors[jm2 as usize][(ic2 - 1) as usize];
                                }
                                if found {
                                    let tb2 = self.adistcc[(inearbond, ia)];
                                    let dsqtb = tb1 * tb1 + tb2 * tb2 - 2.0 * tb1 * tb2 * self.tbintracosines[itb];
                                    let mut xdif = m.xadd[jm1] + m.xa[ja1] - (m.xadd[jm2] + m.xa[ja2]);
                                    let mut ydif = m.yadd[jm1] + m.ya[ja1] - (m.yadd[jm2] + m.ya[ja2]);
                                    let mut zdif = m.zadd[jm1] + m.za[ja1] - (m.zadd[jm2] + m.za[ja2]);
                                    minimagec(&mut xdif, &mut ydif, &mut zdif);
                                    let dsq = xdif * xdif + ydif * ydif + zdif * zdif;
                                    let dtest = (dsq - dsqtb).abs();
                                    if dtest < dsave || dsave == 0.0 {
                                        dsave = dtest;
                                        jamin = ja2;
                                        dsqmin = dsq;
                                        dsqtbmin = dsqtb;
                                    }
                                }
                            }
                            if jamin != 0 {
                                let ja2 = jamin;
                                let jm2 = m.mref[ja2];
                                let tb = dsqtbmin.sqrt();
                                self.nbondcc[ja1] += 1;
                                self.nbondcc[ja2] += 1;
                                nbondccsize = nbondccsize.max(self.nbondcc[ja1]).max(self.nbondcc[ja2]);
                                self.realloc_bonds(nbondccsize, natom);
                                let k1 = self.nbondcc[ja1];
                                let k2 = self.nbondcc[ja2];
                                self.set_bond(k1, ja1, ja2, 3, -tb, dsqmin.sqrt());
                                let ad = self.adistcc[(k1, ja1)];
                                self.set_bond(k2, ja2, ja1, 3, -tb, ad);
                                nacceptors[jm1] += 1;
                                ndonors[jm2] += 1;
                                acceptors[jm1 as usize].push(jm2);
                                donors[jm2 as usize].push(jm1);
                            }
                        }
                    }
                }
                let mut sumsame = 0;
                let (mut nacceptorstotal, mut ndonorstotal) = (0, 0);
                for im in 1..=nmol {
                    nacceptorstotal += nacceptors[im];
                    ndonorstotal += ndonors[im];
                    for &jm in &acceptors[im as usize] {
                        for &d in &donors[im as usize] {
                            if d == jm {
                                sumsame += 1;
                            }
                        }
                    }
                }
                log(&format!(
                    " find_constraints - pass 5> {} {itb} Number of donor-donor pairs (should be zero) {sumsame}",
                    self.atc
                ));
                log(&format!(" find_constraints - pass 5> Total numbers of acceptors and donors: {nacceptorstotal} {ndonorstotal}"));
            }
        }
    }

    /// `function en_coordcon(this, ia, ...)`: constraint energy of atom `ia`.
    pub fn en_coordcon(&self, ia: i32, m: &Mol) -> f32 {
        let mut en_coordcon: f32 = 0.0;
        let nbond = self.nbondcc[ia];
        if nbond > 0 {
            let im = m.mref[ia];
            let x1 = m.xadd[im] + m.xa[ia];
            let y1 = m.yadd[im] + m.ya[ia];
            let z1 = m.zadd[im] + m.za[ia];
            for j in 1..=nbond {
                let ja = self.iabondcc[(j, ia)];
                let jm = m.mref[ja];
                let bothtethered = m.tethered[im] && m.tethered[jm] && m.qtype[ia] == 1 && m.qtype[ja] == 1;
                if !bothtethered {
                    let x2 = m.xadd[jm] + m.xa[ja];
                    let y2 = m.yadd[jm] + m.ya[ja];
                    let z2 = m.zadd[jm] + m.za[ja];
                    let mut xdif = x2 - x1;
                    let mut ydif = y2 - y1;
                    let mut zdif = z2 - z1;
                    minimagec(&mut xdif, &mut ydif, &mut zdif);
                    let dsq = xdif * xdif + ydif * ydif + zdif * zdif;
                    let dtest = dsq.sqrt();
                    let tb = self.tbondcc[(j, ia)];
                    let (mut wtms, mut distort, mut gweight, mut gwidth) = (self.cctemp, self.ccdistort, self.ccgaussweight, self.ccgausswidth);
                    if self.wtcc[(j, ia)] == 2 {
                        wtms = self.tbweight;
                        distort = self.tbdistort;
                        gweight = self.tbgaussweight;
                        gwidth = self.tbgausswidth;
                    } else if self.wtcc[(j, ia)] == 3 {
                        wtms = self.tbintraweight;
                        distort = self.tbintradistort;
                        gweight = self.tbintragaussweight;
                        gwidth = self.tbintragausswidth;
                    }
                    let wtm = m.amassr[m.npairref[(m.itype[ia], m.itype[ja])]] / tb.abs();
                    let en = en_atom_inter_scale_coord(dtest, tb, wtm, wtms, distort, gweight, gwidth);
                    en_coordcon += en;
                }
            }
        }
        en_coordcon
    }

    /// The `coordcon` line of the `.pcof` file (`save_files_others`): a
    /// `character*256` record written with `'(a)'`, i.e. blank-padded to 256.
    pub fn pcof_line(&self) -> String {
        let formatline = format!("(a,1x,a,1x,i3,2(1x,f12.3,1x,f8.3),{:02}(1x,a,2(1x,f8.3)))", self.ns);
        let mut v = vec![
            FV::S("coordcon"),
            FV::Ch(&self.atc, 3),
            FV::I(self.nc as i64),
            FV::R(self.cctemp),
            FV::R(self.ccdistort),
            FV::R(self.ccgaussweight),
            FV::R(self.ccgausswidth),
        ];
        for i in 1..=self.ns {
            v.push(FV::Ch(&self.ats[i], 3));
            v.push(FV::R(self.ds[i]));
            v.push(FV::R(self.dsmax[i]));
        }
        let shortline = fstr(&fmt(&formatline, &v), 256);
        let mut fullline = shortline;
        if self.hastb {
            // NB: the Fortran format strings of the angle parts end with one ')' too
            // many; gfortran stops at the matching parenthesis.
            let mut v = vec![FV::S("angles"), FV::R(self.tbweight), FV::R(self.tbdistort), FV::R(self.tbgaussweight), FV::R(self.tbgausswidth)];
            for i in 1..=self.ntbangles {
                v.push(FV::R(self.tbangles[i]));
            }
            let fl = format!("(1x,a,2(1x,f12.3,1x,f8.3),{:02}(1x,f8.3))", self.ntbangles);
            let shortline = fstr(&fmt(&fl, &v), 256);
            fullline = fstr(&format!("{}{}", fullline, crate::fortran::format::pad_to(&shortline, 256)), 256);
        }
        if self.hastbintra {
            let fl = format!("(1x,a,2(1x,f12.3,1x,f8.3),{:02}(1x,a3,1x,f8.3))", self.ntbintraangles);
            let mut v = vec![
                FV::S("anglesintra"),
                FV::R(self.tbintraweight),
                FV::R(self.tbintradistort),
                FV::R(self.tbintragaussweight),
                FV::R(self.tbintragausswidth),
            ];
            for i in 1..=self.ntbintraangles {
                v.push(FV::Ch(&self.tbintraats[i], 3));
                v.push(FV::R(self.tbintraangles[i]));
            }
            let shortline = fstr(&fmt(&fl, &v), 256);
            fullline = fstr(&format!("{}{}", fullline, crate::fortran::format::pad_to(&shortline, 256)), 256);
        }
        crate::fortran::format::pad_to(&fullline, 256)
    }
}

impl CoordConstraints {
    /// `function check_coord_bonded(ia, ja, mref, itype)`.
    #[inline]
    pub fn check_coord_bonded(&self, ia: i32, ja: i32, m: &Mol) -> bool {
        let mut check = false;
        let mut ic = 0;
        while ic < self.ncoordconstraints && !check {
            ic += 1;
            let c = &self.coordconstraints[(ic - 1) as usize];
            let mut jb = 0;
            let nbond = c.nbondcc[ia];
            while jb < nbond && !check {
                jb += 1;
                check = ja == c.iabondcc[(jb, ia)];
            }
            if m.itype[ia] != m.itype[ja] && !check {
                jb = 0;
                let nbond = c.nbondcc[ja];
                while jb < nbond && !check {
                    jb += 1;
                    check = ia == c.iabondcc[(jb, ja)];
                }
            }
        }
        check
    }

    /// `function en_mole_coord(im, ...)`: constraint energy of molecule `im`.
    pub fn en_mole_coord(&self, im: i32, m: &Mol) -> f32 {
        let mut e: f32 = 0.0;
        for c in self.coordconstraints.iter().take(self.ncoordconstraints as usize) {
            for ia in m.naf[im]..=m.nal[im] {
                if c.nbondcc[ia] > 0 {
                    let en = c.en_coordcon(ia, m);
                    e += en;
                }
            }
        }
        e
    }

    /// `subroutine save_coord_energies(im)` (called by `molfit` only after an
    /// accepted internal-rotation move; NB: not after whole-molecule moves).
    pub fn save_coord_energies(&mut self, im: i32) {
        if self.ncoordconstraints > 0 {
            self.ergold_coord[im] = self.ergnew_coord[im];
        }
    }

    /// `subroutine initialise_coords(n)`.
    pub fn initialise_coords(&mut self, n: i32) {
        realloc1(&mut self.ergold_coord, n as usize);
        realloc1(&mut self.ergnew_coord, n as usize);
    }
}

impl State {
    /// `read_pcof`: a `coordcon` line (`sometext`, `nwords > 8`).
    pub(crate) fn read_coordcon_line(&mut self, sometext: &str, w: &crate::fortran::text::Words) -> Result<()> {
        self.cc.ncoordconstraints += 1;
        let text = fstr(&w.span(sometext, 2, w.nwords()), 256);
        self.log(&format!(" read_pcof> No. coord. constraints: {} {}", self.cc.ncoordconstraints, text));
        let log = self.logger();
        let c = CoordConstraint::init_constraint(&text, &log)?;
        let n = self.cc.ncoordconstraints as usize;
        // reallocate1d_cc keeps the existing elements
        self.cc.coordconstraints.truncate(n - 1);
        self.cc.coordconstraints.push(c);
        Ok(())
    }

    /// `fmolec1`: `initialise_coords(natom)`, `find_constraints` for every
    /// constraint and `writecoordcons`.
    pub fn setup_coord_constraints(&mut self) -> Result<()> {
        if self.cc.ncoordconstraints > 0 {
            let natom = self.mol.natom;
            self.cc.initialise_coords(natom);
            self.find_all_constraints();
            self.writecoordcons()?;
        }
        Ok(())
    }

    /// `do id=1,ncoordconstraints; call find_constraints(coordconstraints(id),...)`.
    pub fn find_all_constraints(&mut self) {
        let log = self.logger();
        let n = self.cc.ncoordconstraints as usize;
        for c in self.cc.coordconstraints.iter_mut().take(n) {
            c.find_constraints(&self.mol, &log);
        }
    }

    /// `subroutine writecoordcons(natom, itype, atype)`: `coordcon.txt`.
    pub fn writecoordcons(&self) -> Result<()> {
        let m = &self.mol;
        let mut out = OutFile::create(&self.sym.full("coordcon.txt"));
        for ia in 1..=m.natom {
            for c in self.cc.coordconstraints.iter().take(self.cc.ncoordconstraints as usize) {
                let nbond = c.nbondcc[ia];
                let nearnbond = c.nearnbondcc[ia];
                let lab = &m.type_[m.itype[ia]];
                if nbond > 0 {
                    let formatline = format!("(i6,1x,a,1x,i2,1x,i2,{:02}(i6,1x,a,1x,e10.3))", nbond);
                    let mut v = vec![FV::I(ia as i64), FV::Ch(lab, 3), FV::I(nearnbond as i64), FV::I(nbond as i64)];
                    for j in 1..=nbond {
                        let ja = c.iabondcc[(j, ia)];
                        v.push(FV::I(ja as i64));
                        v.push(FV::Ch(&m.type_[m.itype[ja]], 3));
                        v.push(FV::R(c.tbondcc[(j, ia)]));
                    }
                    out.line(&fmt(&formatline, &v));
                } else {
                    out.line(&fmt("(i6,1x,a,1x,i2)", &[FV::I(ia as i64), FV::Ch(lab, 3), FV::I(nbond as i64)]));
                }
            }
        }
        out.close()?;
        Ok(())
    }

    /// `normalise_coordinates` / `unnormalise_coordinates`: the bond lengths
    /// `ds`, `dsmax` of each constraint (`fractcoord`/`wholecoord` with `dmax`).
    pub(crate) fn scale_coordconstraint_lengths(&mut self, normalise: bool) {
        let dmax = self.mol.dmax;
        let n = self.cc.ncoordconstraints as usize;
        for c in self.cc.coordconstraints.iter_mut().take(n) {
            for j in 1..=c.ns {
                if normalise {
                    c.ds[j] /= dmax;
                    c.dsmax[j] /= dmax;
                } else {
                    c.ds[j] *= dmax;
                    c.dsmax[j] *= dmax;
                }
            }
        }
    }

    /// The `coordcon` lines at the end of the `.pcof` file (`save_files_others`).
    pub(crate) fn write_pcof_coordcons(&self, out: &mut OutFile) {
        for c in self.cc.coordconstraints.iter().take(self.cc.ncoordconstraints as usize) {
            out.line(&c.pcof_line());
        }
    }
}
