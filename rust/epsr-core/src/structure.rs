//! Radial distribution functions, structure factors and Fourier transforms.
//!
//! Translation of (all in `epsr_standalone_rev3.f` unless noted):
//! * `inita`    — r and Q grids, Fourier kernel `foqols`, reads the `.inpa` accumulators
//! * `uniform`  — number of atoms per r-shell at uniform density (`ycoh`)
//! * `j1qr`     — ∫ r sin(Qr)/Q dr  (spherical-shell Fourier kernel)
//! * `interp`   — linear interpolation of data onto the simulation Q grid
//! * `atogen`   — histogram the configuration into partial g(r) (`ato_manipulation_rev4.f`)
//! * `gensq`    — partial S(Q) from the histograms
//! * `sumsq`, `sumgr` — running averages (and spreads) of S(Q) and g(r)
//! * `entrpb`, `entrpc` — "noise" (roughness) estimates of g(r)
//! * `stog`, `lorch` — S(Q) → g(r) with a Lorch window (`ato_manipulation_rev4.f`)
//! * `ftback`   — back-transform used to remove long-wavelength background from data

use crate::error::{EpsrError, Result};
use crate::fortran::fmath::nint;
use crate::fortran::listread::{read_int, read_real, read_word, TextFile};
use crate::fortran::text::parse;
use crate::fortran::{F1, F2, F3};
use crate::neighbour::minimagec;
use crate::state::{State, MCOEFF, MPZ};
use rayon::prelude::*;

/// `function j1qr(q, r)` (real*8 result).
fn j1qr(q: f32, r: f32) -> f64 {
    if q > 0.0 {
        let qr = (q * r) as f64;
        (qr.sin() - qr * qr.cos()) / ((q * q * q) as f64)
    } else {
        (r * r * r / 3.0) as f64
    }
}

/// `function pofr(delr)`: Fourier transform of a uniform sphere (Lorch-type
/// broadening function), single precision.
fn pofr(delr: f32) -> f32 {
    if delr <= 0.0 {
        1.0
    } else {
        let delr3 = 3.0 / (delr * delr * delr);
        delr3 * (delr.sin() - delr * delr.cos())
    }
}

/// `subroutine uniform(ycoh, xst, rho, rmaxgr, ndim)`.
pub fn uniform(ycoh: &mut F1<f32>, xst: f32, rho: f32, rmaxgr: f32, ndim: i32) {
    let pi = std::f32::consts::PI;
    let rhofac = 4.0 * pi * rho / 3.0;
    let mut x1: f32 = 0.0;
    let mut x13: f32 = 0.0;
    let mut xstep = xst * 0.5;
    let mut x2 = x1 + xstep;
    ycoh[1] = 0.0;
    let mut j = 1;
    while x2 < rmaxgr && j < ndim {
        j += 1;
        x2 = x1 + xstep;
        let x23 = x2 * x2 * x2;
        ycoh[j] = rhofac * (x23 - x13);
        xstep = xst;
        x1 = x2;
        x13 = x23;
    }
}

/// `SUBROUTINE INTERP(NQD,NZPTS,QD,QI,DATA,ALAMBI,in,NQMIN,NQMAX,...)`:
/// linear interpolation of (`qd`,`data`) onto `qi`; returns the first and
/// last grid points covered by the data.
pub fn interp(nqd: i32, nzpts: i32, qd: &F1<f32>, qi: &F1<f32>, data: &F1<f32>, alambi: &mut F2<f32>, col: i32) -> (i32, i32) {
    let mut q1 = qd[1];
    let mut q2 = qd[2];
    let mut d1 = data[1];
    let mut d2 = data[2];
    let mut n = 2;
    let mut i = 1;
    let mut nqmin = 1;
    loop {
        // label 5
        let qt = qi[i];
        if q1 > qt {
            // label 15
            i += 1;
            nqmin = i;
            continue;
        }
        // label 10
        loop {
            if q2 < qt {
                // label 20
                n += 1;
                if n > nqd {
                    return (nqmin, i - 1);
                }
                d1 = d2;
                q1 = q2;
                q2 = qd[n];
                d2 = data[n];
                continue;
            }
            break;
        }
        let grad = (d2 - d1) / (q2 - q1);
        alambi[(i, col)] = d1 + (qt - q1) * grad;
        i += 1;
        if i > nzpts {
            return (nqmin, i - 1);
        }
    }
}

/// `function lorch(x, pi)`: Lorch window sin(x)/x for x < π.
#[inline(always)]
pub fn lorch(x: f32, pi: f32) -> f32 {
    if x == 0.0 {
        1.0
    } else if x < pi {
        x.sin() / x
    } else {
        0.0
    }
}

/// `subroutine stog(rho,q,sofq,r,gofr,nqmin,nqmax,nr,npairs,...,xwt)`:
/// Fourier transform S(Q) → g(r) (Lorch-windowed) for `npairs` columns.
#[allow(clippy::too_many_arguments)]
pub fn stog(rho: f32, q: &F1<f32>, sofq: &F2<f32>, r: &F1<f32>, gofr: &mut F2<f32>, nqmin: &F1<i32>, nqmax: &F1<i32>, nr: i32, npairs: i32, xwt: &F1<f32>) {
    let pi = std::f32::consts::PI;
    let qstep = q[4] - q[3];
    let sfac = 0.5 * qstep / pi / pi / rho;
    let nq = sofq.n1;
    // The r-rows are independent (per-row scratch, serial inner Q sums), so
    // they are computed on the rayon pool and written back in row order; the
    // result is identical at any thread count.
    let rows: Vec<Vec<f32>> = (1..nr + 1)
        .into_par_iter()
        .map(|j| {
            let mut sum: Vec<f32> = vec![0.0; nq];
            let rr = r[j];
            let mut row: Vec<f32> = Vec::with_capacity(npairs.max(0) as usize);
            for inn in 1..=npairs {
                let qmax = pi / q[nqmax[inn]];
                for i in nqmin[inn]..=nqmax[inn] {
                    let qq = q[i];
                    let prod = if rr > 0.0 {
                        let qr = qq * rr;
                        qq * qr.sin() * lorch(qmax * qq, pi)
                    } else {
                        qq * qq * lorch(qmax * qq, pi)
                    };
                    sum[(i - 1) as usize] = sofq[(i, inn)] * prod;
                }
                let mut sumi: f32 = 0.0;
                for i in nqmin[inn]..=nqmax[inn] {
                    sumi += sum[(i - 1) as usize];
                }
                row.push(if rr > 0.0 { xwt[inn] + sfac * sumi / rr } else { xwt[inn] + sfac * sumi });
            }
            row
        })
        .collect();
    for j in 1..=nr {
        let row = &rows[(j - 1) as usize];
        for inn in 1..=npairs {
            gofr[(j, inn)] = row[(inn - 1) as usize];
        }
    }
}

impl State {
    /// `subroutine inita(ireset,iinit,nq,q,qstep,nz,z,fwhm,fwhmq,intra)`.
    /// Returns `(q, z, nz)`.
    #[allow(clippy::too_many_arguments)]
    pub fn inita(&mut self, ireset: i32, iinit: i32, nq: i32, qstep: f32, fwhm: f32, fwhmq: f32, intra: i32) -> Result<(F1<f32>, F1<f32>, i32)> {
        let pi = std::f32::consts::PI;
        self.dist.entwt[1] = 0.0;
        // ---- EP parameters from the symbol table
        let sv = |s: &State, k: i32| s.sym.symbol_value[k].clone();
        let mut iref = self.sym.ireff(self.sym.iqminref, 1);
        iref += 1;
        self.inp.roverlap = read_real(&sv(self, iref))?;
        iref += 1;
        self.inp.rminfac = read_real(&sv(self, iref))?;
        iref += 1;
        self.inp.rminpt = read_real(&sv(self, iref))?.max(1.0);
        iref += 1;
        self.inp.rmaxpt = read_real(&sv(self, iref))?.max(2.0);
        iref += 1;
        self.inp.rbroad = read_real(&sv(self, iref))?;
        iref += 1;
        self.inp.expecf = 0.3;
        self.inp.expecfpwr = 0.0;
        {
            let t = sv(self, iref);
            let w = parse(&t, 11);
            if w.nwords() > 0 {
                self.inp.expecf = read_real(w.word(&t, 1))?;
            }
            if w.nwords() > 1 {
                self.inp.expecfpwr = read_real(w.word(&t, 2))?;
            }
        }
        iref += 1;
        self.pot.harmonic = false;
        if self.sym.symbol_value[iref] == "harmonic" {
            self.pot.harmonic = true;
        } else {
            self.sym.symbol_value[iref] = "exponential".into();
        }
        iref += 1;
        self.pot.addpotential = read_word(&sv(self, iref))?;
        iref += 1;
        self.inp.rcharge = read_real(&sv(self, iref))?;
        iref += 1;
        self.mol.power = read_real(&sv(self, iref))?;
        if self.mol.power <= 6.0 {
            self.mol.power = 12.0;
        }
        iref += 1;
        let p = &mut self.pot;
        p.pdstep = 0.1;
        p.pdmin = 1.0 + p.pdstep / self.inp.rminpt;
        p.pdmax = self.inp.rmaxpt;
        let rminpt = self.inp.rminpt;
        let rmaxpt = self.inp.rmaxpt;
        p.ncoeffp = (((rminpt / p.pdstep) as i32) + (((rmaxpt / rminpt).ln() / p.pdmin.ln()) as i32)).min(MCOEFF as i32);
        let msg = format!(" inita> ncoeffp = {:>11}", p.ncoeffp);
        self.log(&msg);
        let _ncoeffpread = read_int(&sv(self, iref))?;
        iref += 1;
        let _pdmaxread = read_real(&sv(self, iref))?;
        iref += 1;
        let _pdstepread = read_real(&sv(self, iref))?;
        iref += 1;
        self.pot.npitss = read_int(&sv(self, iref))?;
        iref += 1;
        self.pot.paccept = read_real(&sv(self, iref))?;
        iref += 1;
        self.pot.gaussian = false;
        {
            let v = sv(self, iref);
            if v == "T" || v == "t" || v == "Gaussian" {
                self.pot.gaussian = true;
                self.sym.symbol_value[iref] = "Gaussian".into();
            } else {
                self.sym.symbol_value[iref] = "Poisson".into();
            }
        }
        iref += 1;
        self.pot.psigma2read = read_real(&sv(self, iref))?;
        iref += 1;
        self.pot.refpotfac = read_real(&sv(self, iref))?;
        self.pot.psigma2 = self.pot.psigma2read;
        if self.pot.psigma2 < 0.004 {
            self.pot.psigma2 = 0.004;
            self.pot.psigma2read = self.pot.psigma2;
        }
        // ---- accumulators from the .inpa file
        if ireset == 0 {
            if iinit == 0 {
                let fname = format!("{}.inpa", self.sym.fnamerun);
                let mut f = TextFile::open(std::path::Path::new(&fname))
                    .map_err(|e| EpsrError::Io(format!("cannot open {fname}: {e}")))?;
                {
                    let i = &mut self.inp;
                    let mut r = f.list();
                    r.get_i(&mut i.nisum)?;
                    r.get_i(&mut i.nesum)?;
                    r.get_i(&mut i.nesumt)?;
                    r.get_r(&mut i.ensum)?;
                    r.get_r(&mut i.ensum2)?;
                    r.get_r(&mut i.prsum)?;
                    r.get_r(&mut i.prsum2)?;
                }
                let nfile = self.dist.nfile;
                for id in 1..=nfile {
                    f.list().finish()?;
                    {
                        let mut r = f.list();
                        r.get_i(&mut self.pot.npcoeff[id])?;
                        r.get_r(&mut self.pot.pstepd[id])?;
                        r.get_r(&mut self.dist.sigma2[id])?;
                    }
                    let mut r = f.list();
                    for ic in 1..=self.pot.npcoeff[id] {
                        r.get_r(&mut self.pot.pcoefd[(ic, id)])?;
                    }
                }
                {
                    let mut r = f.list();
                    r.get_r(&mut self.pot.feedbackold)?;
                    for inn in 1..=nfile {
                        r.get_r(&mut self.pot.efacold[inn])?;
                    }
                    r.get_i(&mut self.wts.nqinv)?;
                }
                let nnew = nfile + self.mol.npairs;
                for iq in 1..=self.wts.nqinv {
                    for ip in 1..=self.mol.npairs {
                        {
                            let mut r = f.list();
                            r.get_i(&mut self.wts.qinv[iq])?;
                        }
                        let mut r = f.list();
                        for id in 1..=nnew {
                            r.get_d(&mut self.wts.cwttq[(iq, ip, id)])?;
                        }
                    }
                }
            } else {
                self.reset_accumulators(intra);
            }
            if self.inp.nsums >= 0 {
                self.inp.nsum = self.inp.nsums;
            } else {
                self.inp.nsum = 0;
            }
            self.pot.npits = self.pot.npitss;
        } else {
            self.reset_accumulators(intra);
            self.pot.npits = self.pot.npitss;
        }
        // ---- r grid
        let d = &mut self.dist;
        let inp = &mut self.inp;
        if inp.rmaxgr > 0.5 * d.cellmax || inp.rmaxgr <= 0.0 {
            inp.rmaxgr = 0.5 * d.cellmax;
        }
        let msg = format!(" epsr-inita> rmaxgr,rmaxcell,cell,cellmax = {} {} {} {}", inp.rmaxgr, d.rmaxcell, d.cell, d.cellmax);
        d.npt = 0;
        let celstt = d.cellst;
        let mut x1dold = -0.5 * celstt;
        let mut cellhf = 0.25 * celstt;
        d.x[1] = 0.0;
        let mut j: i32 = 2;
        d.npt1 = 2;
        d.npt2 = 2;
        d.nptmin = 2;
        d.nptmax = 2;
        let mut x1d = x1dold + celstt;
        while (j as usize) < MPZ && x1d < inp.rmaxgr {
            d.x[j] = x1d - cellhf;
            d.entwt[j] = 1.0 / d.x[j];
            x1dold = x1d;
            cellhf = 0.5 * celstt;
            if x1d < inp.rminpt {
                d.nptmin = j + 1;
            }
            if x1d < inp.rmaxpt {
                d.nptmax = j + 1;
            }
            if x1d < inp.rmaxgr {
                d.npt1 = j;
            }
            j += 1;
            x1d = x1dold + celstt;
        }
        d.x[j] = x1d - cellhf;
        d.npt = j;
        d.npt2 = j - 2;
        if d.npt1 > d.npt2 {
            d.npt1 = d.npt2;
        }
        let msg2 = format!("{:>6}{:>6}{:>6}{:>6}{:>6}", d.nptmin, d.nptmax, d.npt1, d.npt2, d.npt);
        if inp.nesumt < 1 {
            inp.nesumt = 1;
        }
        if inp.nesum > inp.nesumt - 1 {
            inp.nesum = inp.nesumt - 1;
        }
        if inp.nisum > intra {
            inp.nisum = intra;
        }
        self.log(&msg);
        self.log(&msg2);
        let l = &self.lat;
        let volume = l.alatx * l.blaty * l.clatz;
        let rhotemp = self.mol.natomr as f32 / volume;
        let rmaxgr = self.inp.rmaxgr;
        uniform(&mut self.dist.ycoh, self.dist.cellst, rhotemp, rmaxgr, MPZ as i32);
        // ---- Q grid and Fourier kernel
        let mut q: F1<f32> = F1::new(crate::state::MPTS);
        for i in 1..=nq {
            q[i] = qstep * (i - 1) as f32;
        }
        let ln2x2sqrt = (2.0f32 * 2.0f32.ln()).sqrt();
        let sigres = 0.5 * fwhm / ln2x2sqrt;
        let sigrsq = 0.5 * fwhmq / ln2x2sqrt;
        let trfac = 4.0 * pi * rhotemp;
        let lorchlim: f32 = 4.4934 / rmaxgr;
        let (lorchqwidth, lorchqwidthq) = (self.inp.lorchqwidth, self.inp.lorchqwidthq);
        let d = &mut self.dist;
        for i in 1..=nq {
            let qq = q[i];
            let mut r1: f32 = 0.0;
            let mut j1qr1 = j1qr(qq, r1);
            for j in 1..=d.npt1 {
                let rr = d.x[j];
                let rstep = rr - r1;
                let r2 = rr + rstep;
                let j1qr2 = j1qr(qq, r2);
                let mut v: f32 = if d.ycoh[j] > 0.0 {
                    ((trfac as f64 * (j1qr2 - j1qr1)) / d.ycoh[j] as f64) as f32
                } else {
                    0.0
                };
                r1 = r2;
                j1qr1 = j1qr2;
                let qsig = sigres + qq * sigrsq;
                let qrs = qsig * rr;
                let qr2 = 0.5 * qrs * qrs;
                let mut sums: f32 = 0.0;
                if qr2 < 30.0 {
                    sums = (-qr2).exp();
                    if lorchqwidth > 0.0 || lorchqwidthq > 0.0 {
                        let bwidth = rr * lorchlim * (lorchqwidth + lorchqwidthq * qq);
                        sums *= pofr(bwidth);
                    }
                }
                v *= sums;
                d.foqols[(j, i)] = v;
            }
        }
        if self.mol.coarsegraining {
            self.inita_coarsegrain(nq, &q)?;
        }
        let nz = self.dist.npt1;
        let mut z: F1<f32> = F1::new(MPZ);
        for i in 1..=nz {
            z[i] = self.dist.x[i];
        }
        Ok((q, z, nz))
    }

    /// Re-initialisation of the accumulators and per-dataset EP fits used by
    /// `inita` when `ireset /= 0` or `iinit /= 0`.
    fn reset_accumulators(&mut self, intra: i32) {
        let i = &mut self.inp;
        i.nsum = 0;
        i.nesum = 0;
        i.nesumt = 5;
        i.nsums = -1;
        i.nisum = intra;
        for id in 1..=self.dist.nfile {
            self.pot.npcoeff[id] = self.pot.ncoeffp;
            self.dist.sigma2[id] = self.pot.psigma2;
            for ic in 1..=self.pot.npcoeff[id] {
                self.pot.pcoefd[(ic, id)] = 0.0;
            }
            self.pot.pstepd[id] = 0.002;
        }
    }

    /// `subroutine atogen(wav, ngrsamples)`: histogram the current
    /// configuration into intramolecular (`wav(1,..)`) and intermolecular
    /// (`wav(2,..)`) pair-distance distributions, normalised by the number of
    /// pairs, and count close contacts (`potex`).
    pub fn atogen(&mut self, wav: &mut F3<f32>, ngrsamples: &mut i32) {
        let npt1 = self.dist.npt1;
        let npairs = self.mol.npairs;
        for iz in 1..=npt1 {
            for inn in 1..=npairs {
                wav[(1, iz, inn)] = 0.0;
                wav[(2, iz, inn)] = 0.0;
            }
        }
        for ic in 1..=npairs {
            self.pot.potex[ic] = 0.0;
            self.pot.potexmax[ic] = 0.0;
        }
        if self.mol.coarsegraining {
            self.atogen_cg_alloc();
        }
        let m = &self.mol;
        let stepcl = self.dist.cellst / m.dmax;
        // intramolecular distances
        for im in 1..=m.nmol {
            if m.nato[im] > 1 {
                for ia in m.naf[im]..m.nal[im] {
                    let x1 = m.xa[ia];
                    let y1 = m.ya[ia];
                    let z1 = m.za[ia];
                    for ja in ia + 1..=m.nal[im] {
                        let it = m.itype[ia];
                        let jt = m.itype[ja];
                        let icref = m.npairref[(it, jt)];
                        let mut xx = (x1 - m.xa[ja]) * m.sizefactor;
                        let mut yy = (y1 - m.ya[ja]) * m.sizefactor;
                        let mut zz = (z1 - m.za[ja]) * m.sizefactor;
                        minimagec(&mut xx, &mut yy, &mut zz);
                        let dtest = (xx * xx + yy * yy + zz * zz).sqrt();
                        let iz1 = nint(dtest / stepcl) + 2;
                        if iz1 > 1 && iz1 <= npt1 {
                            wav[(1, iz1, icref)] += 2.0;
                        }
                    }
                }
            }
        }
        // intermolecular distances, from a sample of origin molecules
        let ntype = m.ntype;
        let nmol = m.nmol;
        let mut csampled: F1<f32> = F1::new(ntype as usize);
        let mut molsampled: F1<i32> = F1::new(m.nmoltype as usize);
        let mut molfraction: F1<f32> = F1::new(m.nmoltype as usize);
        if *ngrsamples <= 0 {
            *ngrsamples = 1000;
        }
        if *ngrsamples > nmol {
            *ngrsamples = nmol;
        }
        let probaccept = *ngrsamples as f32 / nmol as f32;
        self.log(&format!(" atogen> Acceptance probability = {probaccept}"));
        for im in 1..=nmol {
            molfraction[m.mtype[im]] += 1.0;
        }
        for it in 1..=m.nmoltype {
            molfraction[it] /= nmol as f32;
        }
        let mut im1 = 1;
        let mut im2 = nmol;
        let mut addwav: f32 = 1.0;
        if probaccept >= 1.0 {
            im1 = 1;
            im2 = nmol - 1;
            addwav = 2.0;
            for ia in m.naf[nmol]..=m.nal[nmol] {
                csampled[m.itype[ia]] += 1.0;
            }
        }
        let sizefactor_unused = m.sizefactor;
        let _ = sizefactor_unused;
        let dmax = m.dmax;
        // The random draws and sampling decisions consume the shared RNG
        // stream and depend on earlier acceptances, so they stay serial; the
        // pair histograms of the accepted molecules are then accumulated in
        // parallel.  Every accumulation is an exact count (+1.0/+addwav in
        // f32), so merging per-thread histograms in any order reproduces the
        // serial result at any thread count.
        let mut accepted: Vec<i32> = Vec::with_capacity(nmol as usize);
        for im in im1..=im2 {
            let r = self.rng.ran1();
            let m = &self.mol;
            if r < probaccept || molsampled[m.mtype[im]] < 1 || molfraction[m.mtype[im]] < probaccept {
                molsampled[m.mtype[im]] += 1;
                for ia in m.naf[im]..=m.nal[im] {
                    csampled[m.itype[ia]] += 1.0;
                }
                accepted.push(im);
            }
        }
        // Container-force counting is independent of the pair loop (it only
        // bumps potexmax by exact counts), and in the serial code ran for
        // every im whether or not it was sampled.
        if self.cont.numcontainerforces > 0 {
            for im in im1..=im2 {
                self.dist_containerforce(im, 1.0);
            }
        }
        struct Acc {
            wav2: F3<f32>,
            potex: F1<f32>,
            wavinter: Vec<f32>,
            iz1values: F1<i32>,
            icrefvalues: F1<i32>,
        }
        let (wn1, wn2, wn3) = (wav.n1, wav.n2, wav.n3);
        let st: &State = self;
        let new_acc = || {
            // Workers inherit the default (unit cubic) lattice, so make this
            // state's lattice the one minimagec uses on this thread.
            st.sync_lattice();
            Acc {
                wav2: F3::new(wn1, wn2, wn3),
                potex: F1::new(st.pot.potex.len()),
                wavinter: vec![0.0; st.cg.wavinter.len()],
                iz1values: F1::new(st.mol.natom as usize),
                icrefvalues: F1::new(st.mol.natom as usize),
            }
        };
        let acc = accepted
            .par_iter()
            .fold(new_acc, |mut acc, &im| {
                let m = &st.mol;
                let xc1 = m.xadd[im];
                let yc1 = m.yadd[im];
                let zc1 = m.zadd[im];
                let (mut jm1, mut jm2) = (1, nmol);
                if probaccept >= 1.0 {
                    jm1 = im + 1;
                    jm2 = nmol;
                }
                for ia in m.naf[im]..=m.nal[im] {
                    let x1 = m.xa[ia] + xc1;
                    let y1 = m.ya[ia] + yc1;
                    let z1 = m.za[ia] + zc1;
                    let it = m.itype[ia];
                    for jm in jm1..=jm2 {
                        if jm != im {
                            let xc2 = m.xadd[jm];
                            let yc2 = m.yadd[jm];
                            let zc2 = m.zadd[jm];
                            for ja in m.naf[jm]..=m.nal[jm] {
                                let jt = m.itype[ja];
                                let icref = m.npairref[(it, jt)];
                                let x2 = m.xa[ja] + xc2;
                                let y2 = m.ya[ja] + yc2;
                                let z2 = m.za[ja] + zc2;
                                let mut xx = x1 - x2;
                                let mut yy = y1 - y2;
                                let mut zz = z1 - z2;
                                minimagec(&mut xx, &mut yy, &mut zz);
                                let dtest = (xx * xx + yy * yy + zz * zz).sqrt();
                                let iz1 = nint(dtest / stepcl) + 2;
                                if iz1 > 1 && iz1 <= npt1 {
                                    acc.iz1values[ja] = iz1;
                                    acc.icrefvalues[ja] = icref;
                                } else {
                                    acc.iz1values[ja] = 0;
                                }
                                let dtesttest = dtest * dmax;
                                let mut scalefac: f32 = 1.0;
                                if m.coarsegraining {
                                    scalefac = 0.5 * (m.scalefacmol[im] + m.scalefacmol[jm]);
                                }
                                let irefp = st.pot.irefpcof[icref];
                                if dtesttest < scalefac * st.pot.rminex[irefp].abs() {
                                    acc.potex[icref] += 1.0;
                                }
                            }
                        } else {
                            for ja in m.naf[jm]..=m.nal[jm] {
                                acc.iz1values[ja] = 0;
                            }
                        }
                    }
                    for ja in m.naf[jm1]..=m.nal[jm2] {
                        let iz1 = acc.iz1values[ja];
                        if iz1 > 0 {
                            let icref = acc.icrefvalues[ja];
                            acc.wav2[(2, iz1, icref)] += addwav;
                        }
                    }
                    // Repeat for coarse graining (CoarseGrain::atogen_accumulate
                    // against the per-thread wavinter buffer)
                    if m.coarsegraining {
                        let ibinia = st.cg.assignatomscalefac[ia];
                        for ja in m.naf[jm1]..=m.nal[jm2] {
                            let iz1 = acc.iz1values[ja];
                            if iz1 > 0 {
                                let jt = m.itype[ja];
                                let ibinja = st.cg.assignatomscalefac[ja];
                                let icref = acc.icrefvalues[ja];
                                // Pairs are referred in the order it <= jt
                                let k = if it <= jt { st.cg.wi(iz1, icref, ibinia, ibinja) } else { st.cg.wi(iz1, icref, ibinja, ibinia) };
                                acc.wavinter[k] += addwav;
                            }
                        }
                    }
                }
                acc
            })
            .reduce(new_acc, |mut a, b| {
                for (x, y) in a.wav2.data.iter_mut().zip(b.wav2.data.iter()) {
                    *x += *y;
                }
                for (x, y) in a.potex.data.iter_mut().zip(b.potex.data.iter()) {
                    *x += *y;
                }
                for (x, y) in a.wavinter.iter_mut().zip(b.wavinter.iter()) {
                    *x += *y;
                }
                a
            });
        for (x, y) in wav.data.iter_mut().zip(acc.wav2.data.iter()) {
            *x += *y;
        }
        for ic in 1..=npairs {
            self.pot.potex[ic] += acc.potex[ic];
        }
        for (x, y) in self.cg.wavinter.iter_mut().zip(acc.wavinter.iter()) {
            *x += *y;
        }
        // normalise
        let m = &self.mol;
        for ic in 1..=npairs {
            let it = m.n1extr[ic];
            let jt = m.n2extr[ic];
            let mut cprod = m.cfrac[it] * m.cfrac[jt] * m.natomr as f32;
            let mut cprod2 = csampled[it] * m.cfrac[jt];
            if it != jt {
                cprod *= 2.0;
                cprod2 += csampled[jt] * m.cfrac[it];
            }
            if cprod > 0.0 {
                for iz in 1..=npt1 {
                    wav[(1, iz, ic)] /= cprod;
                }
            } else {
                for iz in 1..=npt1 {
                    wav[(1, iz, ic)] = 0.0;
                }
            }
            if cprod2 > 0.0 {
                for iz in 1..=npt1 {
                    wav[(2, iz, ic)] /= cprod2;
                }
                self.pot.potex[ic] = self.pot.potex[ic] * m.cfrac[jt] / cprod2;
                if m.coarsegraining {
                    self.cg.atogen_normalise(&self.dist.ycoh, npt1, ic, it, jt, cprod2);
                }
            } else {
                for iz in 1..=npt1 {
                    wav[(2, iz, ic)] = self.dist.ycoh[iz];
                }
                if m.coarsegraining {
                    self.cg.atogen_normalise(&self.dist.ycoh, npt1, ic, it, jt, cprod2);
                }
            }
        }
    }

    /// `subroutine gensq(nq, q, sofq, wav)`: partial structure factors of the
    /// current configuration.  `sofq(1,..)` intramolecular, `sofq(2,..)` intermolecular.
    pub fn gensq(&self, nq: i32, sofq: &mut F3<f32>, wav: &F3<f32>) {
        let d = &self.dist;
        let m = &self.mol;
        // Pairs are independent; within a pair the Q and r summation order is
        // unchanged, so the result is identical at any thread count.
        let cols: Vec<Vec<(f32, f32)>> = (1..=m.npairs)
            .into_par_iter()
            .map(|ic| {
                let jt = m.n2extr[ic];
                let mut col = Vec::with_capacity(nq as usize);
                for i in 1..=nq {
                    let mut sum1: f32 = 0.0;
                    let mut sum2: f32 = 0.0;
                    let mut sum3: f32 = 0.0;
                    for j in 2..=d.npt1 {
                        let add1 = wav[(1, j, ic)];
                        if add1.abs() > 0.0 {
                            sum1 += add1 * d.foqols[(j, i)];
                        }
                    }
                    for j in 2..=d.npt1 {
                        let add2 = (wav[(2, j, ic)] - d.ycoh[j]) * m.natomsphere[jt];
                        sum2 += add2 * d.foqols[(j, i)];
                        sum3 += wav[(2, j, ic)];
                    }
                    if sum3 <= 0.0 {
                        sum2 = 0.0;
                    }
                    col.push((sum1, sum2));
                }
                col
            })
            .collect();
        for ic in 1..=m.npairs {
            let col = &cols[(ic - 1) as usize];
            for i in 1..=nq {
                let (sum1, sum2) = col[(i - 1) as usize];
                sofq[(1, i, ic)] = sum1;
                sofq[(2, i, ic)] = sum2;
            }
        }
    }

    /// `subroutine entrpc(entrop, wav)`: roughness of the current histograms
    /// (note: also sets the two points beyond `npt2` of `wav(2,..)` to `ycoh`).
    pub fn entrpc(&self, wav: &mut F3<f32>) -> f32 {
        let d = &self.dist;
        let n1 = 2;
        let n2 = d.npt2;
        let mut entrop: f32 = 0.0;
        for id in 1..=self.mol.npairs {
            wav[(2, d.npt2 + 1, id)] = d.ycoh[d.npt2 + 1];
            wav[(2, d.npt2 + 2, id)] = d.ycoh[d.npt2 + 2];
            let mut yav1: f32 = 0.0;
            let mut yav2 = wav[(2, n1, id)];
            for j in n1..=n2 {
                let yav3 = wav[(2, j + 1, id)];
                let rj = 0.5 * (yav3 - yav1).abs();
                let pj = 0.25 * (2.0 * yav2 + yav1 + yav3);
                let dif = yav2 - pj;
                let difa = dif.abs();
                let sent = if difa < 0.5 * rj { dif * dif / rj } else { difa - 0.25 * rj };
                entrop += sent * d.entwt[j];
                yav1 = yav2;
                yav2 = yav3;
            }
        }
        entrop
    }

    /// `subroutine entrpb(entrop, nsum)`: roughness of the averaged g(r)
    /// (the Fortran loops over the first `nfile` pairs).
    pub fn entrpb(&self, nsum: i32) -> f32 {
        let d = &self.dist;
        if nsum <= 0 {
            return 0.0;
        }
        let n1 = 2;
        let n2 = d.npt2;
        let mut entrop: f32 = 0.0;
        for id in 1..=d.nfile {
            let mut yav1: f32 = 0.0;
            let mut yav2 = d.gorinter[(n1, id)] * d.ycoh[n1];
            for j in n1..=n2 {
                let yav3 = d.gorinter[(j + 1, id)] * d.ycoh[j + 1];
                let rj = 0.5 * (yav3 - yav1).abs();
                let pj = 0.25 * (2.0 * yav2 + yav1 + yav3);
                let dif = yav2 - pj;
                let difa = dif.abs();
                let sent = if difa < 0.5 * rj { dif * dif / rj } else { difa - 0.25 * rj };
                entrop += sent * d.entwt[j];
                yav1 = yav2;
                yav2 = yav3;
            }
        }
        let ysuma: f32 = 1.0;
        entrop / nsum as f32 / ysuma
    }

    /// `subroutine sumsq(nq, qwind)`: update the running averages of the
    /// partial S(Q) with the current configuration (weights `fc1`, `fc2`).
    pub fn sumsq(&mut self, nq: i32, qwind: &F3<f32>, fc1: f64, fc2: f64) {
        let d = &mut self.dist;
        for id in 1..=self.mol.npairs {
            for j in 1..=nq {
                let add1 = [qwind[(1, j, id)] as f64, qwind[(2, j, id)] as f64];
                let mut sumtot1 = [d.soqintra[(j, id)] as f64, d.soqinter[(j, id)] as f64];
                let mut sumtot2 = [d.sq2intra[(j, id)] as f64, d.sq2inter[(j, id)] as f64];
                let mut errort = [0f64; 2];
                for k in 0..2 {
                    let add2 = add1[k] * add1[k];
                    sumtot2[k] = sumtot2[k] * sumtot2[k] + sumtot1[k] * sumtot1[k];
                    sumtot1[k] = fc2 * sumtot1[k] + fc1 * add1[k];
                    sumtot2[k] = fc2 * sumtot2[k] + fc1 * add2;
                    errort[k] = sumtot2[k] - sumtot1[k] * sumtot1[k];
                    errort[k] = if errort[k] > 0.0 { errort[k].sqrt() } else { 0.0 };
                }
                d.soqintra[(j, id)] = sumtot1[0] as f32;
                d.soqinter[(j, id)] = sumtot1[1] as f32;
                d.sq2intra[(j, id)] = errort[0] as f32;
                d.sq2inter[(j, id)] = errort[1] as f32;
            }
        }
    }

    /// `subroutine sumgr(nz, wav)`: update the running averages of the partial
    /// g(r) and the running coordination integral `entsig`.
    pub fn sumgr(&mut self, nz: i32, wav: &F3<f32>, fc1: f64, fc2: f64) {
        let d = &mut self.dist;
        for id in 1..=self.mol.npairs {
            d.entsig[(1, id)] = 0.0;
            let mut zerolimit: f64 = 0.0;
            for j in 2..=nz {
                let add1: [f64; 2] = if d.ycoh[j] > 0.0 {
                    [(wav[(1, j, id)] / d.ycoh[j]) as f64, (wav[(2, j, id)] / d.ycoh[j]) as f64]
                } else {
                    [0.0, 0.0]
                };
                let mut sumtot1 = [d.gorintra[(j, id)] as f64, d.gorinter[(j, id)] as f64];
                let mut sumtot2 = [d.gr2intra[(j, id)] as f64, d.gr2inter[(j, id)] as f64];
                let mut errort = [0f64; 2];
                for k in 0..2 {
                    let add2 = add1[k] * add1[k];
                    sumtot2[k] = sumtot2[k] * sumtot2[k] + sumtot1[k] * sumtot1[k];
                    sumtot1[k] = fc2 * sumtot1[k] + fc1 * add1[k];
                    sumtot2[k] = fc2 * sumtot2[k] + fc1 * add2;
                    errort[k] = sumtot2[k] - sumtot1[k] * sumtot1[k];
                    errort[k] = if errort[k] > 0.0 { errort[k].sqrt() } else { 0.0 };
                }
                d.gorintra[(j, id)] = sumtot1[0] as f32;
                d.gorinter[(j, id)] = sumtot1[1] as f32;
                d.gr2intra[(j, id)] = errort[0] as f32;
                d.gr2inter[(j, id)] = errort[1] as f32;
                zerolimit += (d.ycoh[j] * (d.gorinter[(j, id)] - 1.0)) as f64;
                d.entsig[(j, id)] = zerolimit as f32;
            }
        }
    }

    /// `subroutine ftback(in, rmin, q, sofq, grt)`: Fourier transform the
    /// difference `sofq` to r-space below `rmin` and back again — the part of
    /// the data/simulation difference that can only come from background
    /// (unphysical r < rmin structure).  `sofq` is replaced by that back-transform.
    pub fn ftback(&self, inn: i32, rmin: f32, q: &F1<f32>, sofq: &mut F1<f32>, grt: &mut F1<f32>) {
        let d = &self.dist;
        let (nqmin, nqmax) = (d.nqmin[inn], d.nqmax[inn]);
        for i in nqmin..=nqmax {
            sofq[i] *= q[i];
        }
        let pi: f32 = 3.141592653;
        let qstep = q[4] - q[3];
        let sfac = 0.5 * qstep / pi / pi / d.rho;
        let mut jsave = 2;
        for j in 2..=d.npt1 {
            let rr = d.x[j];
            if rr >= rmin {
                break;
            }
            jsave = j;
            let mut sum: f32 = 0.0;
            for i in nqmin..=nqmax {
                let qq = q[i];
                let qr = qq * rr;
                sum += sofq[i] * qr.sin();
            }
            grt[j] = sfac * sum;
        }
        let gfac = 4.0 * pi * d.rho;
        for i in nqmin..=nqmax {
            let qq = q[i];
            if qq > 0.0 {
                let mut r1: f32 = 0.0;
                let mut sum: f32 = 0.0;
                for j in 2..=jsave {
                    let rr = d.x[j];
                    let r2 = 2.0 * rr - r1;
                    let rstep = r2 - r1;
                    let qr = qq * rr;
                    sum += grt[j] * qr.sin() * rstep;
                    r1 = r2;
                }
                sofq[i] = gfac * sum / qq;
            } else {
                let mut r1: f32 = 0.0;
                let mut sum: f32 = 0.0;
                for j in 2..=jsave {
                    let rr = d.x[j];
                    let r2 = 2.0 * rr - r1;
                    let rstep = r2 - r1;
                    sum += grt[j] * rr * rstep;
                    r1 = r2;
                }
                sofq[i] = gfac * sum;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn interp_linear() {
        let mut qd = F1::new(10);
        let mut dat = F1::new(10);
        for i in 1..=4i32 {
            qd[i] = i as f32;
            dat[i] = 2.0 * i as f32;
        }
        let mut qi = F1::new(10);
        for i in 1..=6i32 {
            qi[i] = 0.5 * i as f32;
        }
        let mut out: F2<f32> = F2::new(10, 1);
        let (lo, hi) = interp(4, 6, &qd, &qi, &dat, &mut out, 1);
        assert_eq!((lo, hi), (2, 6));
        assert_eq!(out[(2, 1)], 2.0);
        assert_eq!(out[(5, 1)], 5.0);
    }
    #[test]
    fn pow3_is_left_assoc() {
        use crate::fortran::fmath::pow3;
        assert_eq!(pow3(1.1), 1.1f32 * 1.1 * 1.1);
    }
}
