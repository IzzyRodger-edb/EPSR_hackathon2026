//! Neutron weights of the partial structure factors.
//!
//! Translation of `get_wts`, `wtype` (and the trivial `wxinit`) from
//! `epsr_standalone_rev3.f`.
//!
//! A weights file (`.wts`, produced by the `nwts` utility) contains the
//! density-normalisation factor `rhofac` followed by, for each atom pair,
//! a line ` A1  A2` and a line with the intramolecular and intermolecular
//! weights.  From these and the atom fractions the matrices `ciwt` and `cwt`
//! (weight of each partial in each dataset) are formed.  The inverse weights
//! matrix `cwttq` used to distribute the data–simulation differences over the
//! pairs is read from the `.inpa` file by `inita`, or recomputed here
//! (`invert_cwt`) when `feedback` or the dataset factors changed.
//!
//! X-ray weights files (names containing `XWTS.wts`, made by `xwts`) hold
//! form-factor coefficients per atom type; `wxtype` and `wxcreate` turn
//! them into Q-dependent weights `weightx` (and write `X-raywts.dat`,
//! `X-rayf0x.dat`).

use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, FV};
use crate::fortran::listread::{columns, read_real, ListRead, TextFile};
use crate::fortran::output::OutFile;
use crate::fortran::text::parse;
use crate::fortran::F3;
use crate::fortran::{F1, F2};
use crate::state::{State, MPTS, MQINV, MXFILE};

impl State {
    /// `subroutine wtype(a1, a2, wread, in)`: store the weights of pair (a1,a2) for dataset `inn`.
    fn wtype(&mut self, a1: &str, a2: &str, wread: [f32; 2], inn: i32) {
        let m = &self.mol;
        let mut n1 = 0;
        let mut n2 = 0;
        for i in 1..=m.ntype {
            if a1 == m.type_[i] {
                n1 = i;
            }
            if a2 == m.type_[i] {
                n2 = i;
            }
        }
        if n1 > 0 && n2 > 0 {
            let icref = m.npairref[(n1, n2)];
            for j in 1..=2 {
                self.pot.weight[(j, inn, icref)] = wread[(j - 1) as usize];
            }
        } else {
            self.log(&format!(" Mismatched types between atomic coordinates andneutron weights for file {inn:>5}"));
        }
    }

    /// `subroutine get_wts(q, nq, fnamew, ireset)`.
    pub fn get_wts(&mut self, fnamew: &[String], ireset: i32, q: &F1<f32>, nq: i32) -> Result<()> {
        let npairs = self.mol.npairs;
        let nfile = self.dist.nfile;
        for ic in 1..=npairs {
            self.pot.uwt[ic] = 0.0;
        }
        self.wts.nxfile = 0;
        self.xray.weightx = F3::new(MPTS, npairs.max(1) as usize, MXFILE);
        for ix in 1..=MXFILE {
            for it in 1..=self.mol.ntype {
                for iq in 1..=nq {
                    self.xray.f0x[(iq, it, ix)] = 0.0;
                }
            }
        }
        for inn in 1..=nfile {
            for ic in 1..=npairs {
                for j in 1..=2 {
                    self.pot.weight[(j, inn, ic)] = 0.0;
                }
            }
            let fname = &fnamew[(inn - 1) as usize];
            let mut f = TextFile::open(std::path::Path::new(fname)).map_err(|e| EpsrError::Io(format!("cannot open weights file {fname}: {e}")))?;
            if fname.contains("X-ray.wts") || fname.contains("x-ray.wts") || fname.contains("XWTS.wts") {
                self.read_xray_weights(&mut f, inn, q, nq)?;
                continue;
            }
            {
                let mut r = f.list();
                r.get_r(&mut self.dist.rhofac[inn])?;
            }
            self.wts.xfile[inn] = 0;
            let read_pair = |f: &mut TextFile| -> Option<(String, String, [f32; 2])> {
                let l = f.read_line().ok()?;
                let a1 = columns(&l, 1, 3);
                let a2 = columns(&l, 5, 3);
                let mut r = f.list();
                let w1 = r.real().ok()??;
                let w2 = r.real().ok()??;
                Some((a1, a2, [w1, w2]))
            };
            while let Some((a1, a2, w)) = read_pair(&mut f) {
                self.wtype(&a1, &a2, w, inn);
            }
        }
        if self.wts.nxfile > 0 {
            self.wxcreate(q, nq)?;
        }
        let m = &self.mol;
        for inn in 1..=nfile {
            self.pot.wsum[inn] = 0.0;
            self.dist.xwt[inn] = 0.0;
            for ic in 1..=npairs {
                let it = m.n1extr[ic];
                let jt = m.n2extr[ic];
                let mut cprod = m.cfraccg[it] * m.cfraccg[jt];
                if it != jt {
                    cprod *= 2.0;
                }
                self.wts.ciwt[(inn, ic)] = (cprod * self.pot.weight[(1, inn, ic)]) as f64;
                self.wts.cwt[(inn, ic)] = (cprod * self.pot.weight[(2, inn, ic)]) as f64;
                self.pot.wsum[inn] = (self.pot.wsum[inn] as f64 + self.wts.cwt[(inn, ic)].abs()) as f32;
                self.dist.xwt[inn] = (self.dist.xwt[inn] as f64 + self.wts.cwt[(inn, ic)]) as f32;
            }
        }
        let nqinvstep = if self.wts.nxfile > 0 {
            for ix in 1..=self.wts.nxfile {
                for ip in 1..=npairs {
                    let it = m.n1extr[ip];
                    let jt = m.n2extr[ip];
                    let mut cprod = m.cfraccg[it] * m.cfraccg[jt];
                    if it != jt {
                        cprod *= 2.0;
                    }
                    for iq in 1..=nq {
                        self.xray.weightx[(iq, ip, ix)] = cprod as f64 * self.xray.weightx[(iq, ip, ix)];
                    }
                }
            }
            let qinvstep: f32 = 2.0;
            (qinvstep / (q[2] - q[1])) as i32
        } else {
            nq
        };
        // Does the inverted weights matrix need recomputing?
        let mut iflag = 0;
        if ireset == 1 || (self.pot.feedback - self.pot.feedbackold).abs() > 0.001 {
            iflag = 1;
        } else {
            let mut inn = 1;
            while iflag == 0 && inn <= nfile {
                if (self.pot.efacm[inn] - self.pot.efacold[inn]).abs() > 0.001 {
                    iflag = 1;
                }
                inn += 1;
            }
        }
        if iflag != 0 {
            self.recompute_cwttq(nq, nqinvstep)?;
        }
        let m = &self.mol;
        for inn in 1..=nfile {
            for ic in 1..=npairs {
                let mut v = self.wts.cwttq[(1, ic, inn)];
                if m.n1extr[ic] != m.n2extr[ic] {
                    v *= 0.5;
                }
                self.wts.cwtpot[(ic, inn)] = v;
            }
        }
        Ok(())
    }

    /// X-ray branch of `get_wts`: header (`rhofac`, normalisation) and the
    /// form-factor coefficients of each atom type.
    fn read_xray_weights(&mut self, f: &mut TextFile, inn: i32, q: &F1<f32>, nq: i32) -> Result<()> {
        let l = f.read_line()?;
        let w = parse(&l, 5);
        let (normx, iafformaff): (i32, i32);
        if w.nwords() == 3 {
            let t = w.span(&l, 1, 3);
            let mut r = ListRead::from_str(&t);
            self.dist.rhofac[inn] = r.real()?.unwrap_or(0.0);
            normx = r.int()?.unwrap_or(0);
            iafformaff = r.int()?.unwrap_or(0);
        } else {
            let t = w.span(&l, 1, w.nwords());
            self.dist.rhofac[inn] = read_real(&t)?;
            let rf = self.dist.rhofac[inn];
            let nrho = if rf < 0.0 { 2 } else { 1 };
            iafformaff = ((rf.abs() + 1.0) / 2.0) as i32;
            let nx = (rf.abs() - (2 * (iafformaff - 1)) as f32) as i32;
            normx = if nx == 1 { 2 } else { 1 };
            self.dist.rhofac[inn] = nrho as f32;
        }
        if (self.wts.nxfile as usize) < MXFILE {
            self.wts.nxfile += 1;
            let nx = self.wts.nxfile;
            self.xray.datanorm[nx] = self.dist.rhofac[inn] as i32;
            self.dist.rhofac[inn] = 1.0;
            self.xray.normxx[nx] = normx;
            self.xray.iafformaffx[nx] = iafformaff;
            self.wts.xfile[inn] = nx;
            self.xray.datapoint[nx] = inn;
            loop {
                let lab = match f.read_line() {
                    Ok(l) => l,
                    Err(_) => break,
                };
                let a1 = columns(&lab, 1, 3);
                let rec = {
                    let mut r = f.list();
                    (|| -> crate::fortran::IoRes<(f64, [f64; 5], [f64; 5], f32, f32)> {
                        let cxt = r.dble()?.unwrap_or(0.0);
                        let mut axt = [0f64; 5];
                        let mut bxt = [0f64; 5];
                        for i in 0..5 {
                            axt[i] = r.dble()?.unwrap_or(0.0);
                            bxt[i] = r.dble()?.unwrap_or(0.0);
                        }
                        let alphaxt = r.real()?.unwrap_or(0.0);
                        let deltaxt = r.real()?.unwrap_or(0.0);
                        Ok((cxt, axt, bxt, alphaxt, deltaxt))
                    })()
                };
                match rec {
                    Ok((cxt, axt, bxt, alphaxt, deltaxt)) => self.wxtype(&a1, &axt, &bxt, cxt, alphaxt, deltaxt, nx, q, nq),
                    Err(_) => break,
                }
            }
        } else {
            self.log(&format!("\nget_wts> Exceeded limit of {:>2} X-ray files.", MXFILE));
        }
        Ok(())
    }

    /// `subroutine wxtype(a1,axt,bxt,cxt,alphaxt,deltaxt,in,q,nq)`: form factor of atom type `a1`.
    #[allow(clippy::too_many_arguments)]
    fn wxtype(&mut self, a1: &str, axt: &[f64; 5], bxt: &[f64; 5], cxt: f64, alphaxt: f32, deltaxt: f32, inn: i32, q: &F1<f32>, nq: i32) {
        let one: f64 = 1.0;
        let pi4: f64 = 16.0 * one.atan();
        let mut n1 = 0;
        for i in 1..=self.mol.ntype {
            if a1 == self.mol.type_[i] {
                n1 = i;
            }
        }
        if n1 > 0 {
            let x = &mut self.xray;
            for iq in 1..=nq {
                let qx = q[iq] as f64 / pi4;
                let qx2 = qx * qx;
                let mut sum = cxt;
                for i in 0..5 {
                    let term = -bxt[i] * qx2;
                    if term > -40.0 && term < 40.0 {
                        sum += axt[i] * term.exp();
                    }
                }
                x.f0x[(iq, n1, inn)] = sum;
                let fac: f32 = if alphaxt != 0.0 && deltaxt > 0.0 {
                    let expon = q[iq] * q[iq] / (2.0 * deltaxt * deltaxt);
                    if expon < 20.0 {
                        (1.0 - (alphaxt as f64 / x.f0x[(1, n1, inn)]) * ((-expon).exp() as f64)) as f32
                    } else {
                        1.0
                    }
                } else {
                    1.0
                };
                x.fm0x[(iq, n1, inn)] = fac as f64 * sum;
            }
        } else {
            self.log(&format!(" Mismatched types between atomic coordinates andX-ray weights for file {inn:>5}"));
        }
    }

    /// `subroutine wxcreate(q,nq)`: Q-dependent X-ray weights of every pair,
    /// the self-scattering background, and the files `X-raywts.dat`, `X-rayf0x.dat`.
    fn wxcreate(&mut self, q: &F1<f32>, nq: i32) -> Result<()> {
        let m = &self.mol;
        let amolrat = m.natomr as f32 / m.nmol as f32;
        for ix in 1..=self.wts.nxfile {
            let x = &mut self.xray;
            let nrho = x.datanorm[ix];
            let iafformaff = x.iafformaffx[ix];
            let normx = x.normxx[ix];
            for iq in 1..=nq {
                let mut singlem: f32 = 0.0;
                let mut singlei: f32 = 0.0;
                let mut bbar: f64 = 0.0;
                for it in 1..=m.ntype {
                    if m.qtypea[it] == 0 {
                        let cf = m.cfraccg[it] as f64;
                        let fm = x.fm0x[(iq, it, ix)];
                        let f0 = x.f0x[(iq, it, ix)];
                        singlem = (singlem as f64 + cf * (fm * fm)) as f32;
                        singlei = (singlei as f64 + cf * (f0 * f0)) as f32;
                        let fsel = if iafformaff != 2 { f0 } else { fm };
                        if normx == 0 {
                            bbar += cf;
                        } else if normx == 1 {
                            bbar += cf * fsel;
                        } else if normx == 2 {
                            bbar += cf * (fsel * fsel);
                        }
                    }
                }
                if normx == 1 {
                    bbar *= bbar;
                }
                for ip in 1..=m.npairs {
                    let n1 = m.n1extr[ip];
                    let n2 = m.n2extr[ip];
                    x.weightx[(iq, ip, ix)] = x.fm0x[(iq, n1, ix)] * x.fm0x[(iq, n2, ix)] / bbar;
                    if nrho == 2 && normx == 1 {
                        x.weightx[(iq, ip, ix)] /= amolrat as f64;
                    }
                }
                x.backx[(iq, ix)] = if iafformaff == 1 { (singlem - singlei) as f64 / bbar } else { 0.0 };
                if nrho == 2 && normx == 1 {
                    x.backx[(iq, ix)] /= amolrat as f64;
                }
            }
            for ip in 1..=m.npairs {
                let dp = x.datapoint[ix];
                self.pot.weight[(1, dp, ip)] = x.weightx[(1, ip, ix)] as f32;
                self.pot.weight[(2, dp, ip)] = self.pot.weight[(1, dp, ip)];
            }
            let mut out = OutFile::create(&self.sym.full("X-raywts.dat"));
            for iq in 1..=nq {
                let mut v = vec![FV::R(q[iq])];
                for ip in 1..=m.npairs {
                    v.push(FV::D(x.weightx[(iq, ip, ix)]));
                }
                out.line(&fmt("(13(1x,e13.6))", &v));
            }
            out.close()?;
            let mut out = OutFile::create(&self.sym.full("X-rayf0x.dat"));
            for iq in 1..=nq {
                let mut v = vec![FV::R(q[iq])];
                for it in 1..=m.ntype {
                    v.push(FV::D(x.fm0x[(iq, it, ix)]));
                }
                out.line(&fmt("(13(1x,e13.6))", &v));
            }
            out.close()?;
        }
        Ok(())
    }

    /// The `iflag /= 0` block of `get_wts`: rebuild the inverse weights
    /// matrix `cwttq` (feedback or dataset factors changed, or a reset).
    fn recompute_cwttq(&mut self, nq: i32, nqinvstep: i32) -> Result<()> {
        let nfile = self.dist.nfile;
        let npairs = self.mol.npairs;
        let nnew = npairs + nfile;
        let mut cwtt: F2<f64> = F2::new(npairs as usize, nnew as usize);
        for iff in 1..=nnew {
            for ip in 1..=npairs {
                cwtt[(ip, iff)] = if iff - nfile == ip { 1.0 } else { 0.0 };
            }
        }
        self.wts.nqinv = 0;
        let mut iqinv = 1;
        while iqinv <= nq {
            if (self.wts.nqinv as usize) < MQINV {
                self.wts.nqinv += 1;
                let nqinv = self.wts.nqinv;
                self.wts.qinv[nqinv] = iqinv;
                // cwtpot(in, ic) for this Q (neutron data: independent of Q)
                let mut cwtpot: F2<f64> = F2::new(nfile as usize, npairs as usize);
                let m = &self.mol;
                for inn in 1..=nfile {
                    for ic in 1..=npairs {
                        let it = m.n1extr[ic];
                        let jt = m.n2extr[ic];
                        let cprodrat: f32 = if m.cfraccg[it] > 0.0 && m.cfraccg[jt] > 0.0 {
                            m.cfrac[it] * m.cfrac[jt] / (m.cfraccg[it] * m.cfraccg[jt])
                        } else {
                            0.0
                        };
                        let ix = self.wts.xfile[inn];
                        cwtpot[(inn, ic)] = if ix > 0 {
                            cprodrat as f64 * self.xray.weightx[(iqinv, ic, ix)]
                        } else {
                            cprodrat as f64 * self.wts.cwt[(inn, ic)]
                        };
                    }
                }
                let feedback = self.pot.feedback;
                let efacm = self.pot.efacm.clone();
                invert_cwt(&mut cwtpot, &mut cwtt, feedback, &efacm, nfile, npairs, &mut self.rng);
                for inn in 1..=nnew {
                    for ic in 1..=npairs {
                        self.wts.cwttq[(nqinv, ic, inn)] = cwtt[(ic, inn)];
                    }
                }
            }
            iqinv += nqinvstep;
        }
        Ok(())
    }
}

/// `subroutine invert_cwt(cwt,cwtt,carray,carrayold,feedback,efacm,nfile,npairs,...)`:
/// Monte Carlo search for the matrix `cwtt` that best inverts the weights
/// matrix `cwt` (datasets × pairs), regularised by `1-feedback` times the
/// identity on the pair rows, so that differences in the data can be
/// distributed over the individual pair potentials.
fn invert_cwt(cwt: &mut F2<f64>, cwtt: &mut F2<f64>, feedback: f32, efacm: &F1<f32>, nfile: i32, npairs: i32, rng: &mut crate::rng::Ran1) {
    let niter = 50000;
    let feedback1: f32 = 1.0 - feedback;
    let nnew = npairs + nfile;
    for i in 1..=nfile {
        for j in 1..=npairs {
            cwt[(i, j)] = cwt[(i, j)] * feedback as f64 * efacm[i] as f64;
        }
    }
    let totaltest: f32 = 1.0e-15;
    let mut carray: Vec<f64> = vec![0.0; (nnew + 1) as usize];
    let mut carrayold: F2<f64> = F2::new(nnew as usize, nnew as usize);
    let full_total = |cwt: &F2<f64>, cwtt: &F2<f64>, carrayold: &mut F2<f64>| -> f64 {
        let mut total: f64 = 0.0;
        for i1 in 1..=nnew {
            for i2 in 1..=nnew {
                let sum: f64 = if i1 <= nfile {
                    let mut s: f64 = 0.0;
                    for j in 1..=npairs {
                        s += cwt[(i1, j)] * cwtt[(j, i2)];
                    }
                    s
                } else {
                    let j = i1 - nfile;
                    feedback1 as f64 * cwtt[(j, i2)]
                };
                carrayold[(i1, i2)] = sum;
                if i1 == i2 {
                    total += (sum - 1.0) * (sum - 1.0);
                } else {
                    total += sum * sum;
                }
            }
        }
        total
    };
    let mut total = full_total(cwt, cwtt, &mut carrayold);
    let mut totalold = total;
    let mut amp: f32 = (0.01f32 as f64 * total) as f32;
    let mut temp: f64 = (totalold - nfile as f64) / (nnew * nnew) as f64;
    let mut itimes = 1;
    let mut iter = 0;
    let mut nreject = 0;
    while (totalold - nfile as f64) >= totaltest as f64 && itimes <= 1000 {
        for i2 in 1..=nnew {
            for j in 1..=npairs {
                iter += 1;
                let add: f64 = ((2.0 * rng.ran1() - 1.0) * amp) as f64;
                for i1 in 1..=nfile {
                    carray[i1 as usize] = carrayold[(i1, i2)] + cwt[(i1, j)] * add;
                }
                let i1 = j + nfile;
                carray[i1 as usize] = carrayold[(i1, i2)] + feedback1 as f64 * add;
                let mut tdif: f64 = 0.0;
                for i1 in 1..=nfile {
                    let sum = carray[i1 as usize];
                    let sumold = carrayold[(i1, i2)];
                    if i1 == i2 {
                        tdif = tdif + (sum - 1.0) * (sum - 1.0) - (sumold - 1.0) * (sumold - 1.0);
                    } else {
                        tdif = tdif + sum * sum - sumold * sumold;
                    }
                }
                let i1 = j + nfile;
                let sum = carray[i1 as usize];
                let sumold = carrayold[(i1, i2)];
                if i1 == i2 {
                    tdif = tdif + (sum - 1.0) * (sum - 1.0) - (sumold - 1.0) * (sumold - 1.0);
                } else {
                    tdif = tdif + sum * sum - sumold * sumold;
                }
                tdif /= temp;
                let iaccept = if tdif <= 0.0 {
                    true
                } else if tdif > 20.0 {
                    false
                } else {
                    (rng.ran1() as f64) <= (-tdif).exp()
                };
                if iaccept {
                    total += tdif * temp;
                    totalold = total;
                    cwtt[(j, i2)] += add;
                    for i1 in 1..=nfile {
                        carrayold[(i1, i2)] = carray[i1 as usize];
                    }
                    let i1 = j + nfile;
                    carrayold[(i1, i2)] = carray[i1 as usize];
                } else {
                    nreject += 1;
                }
            }
        }
        if iter >= niter {
            let mut fac: f32 = if nreject > 0 { 0.5 * (iter - nreject) as f32 / nreject as f32 } else { 4.0 };
            if fac < 0.95 {
                fac = 0.95;
            }
            if fac > 1.05 {
                fac = 1.05;
            }
            amp *= fac;
            itimes += 1;
            iter = 0;
            nreject = 0;
            total = full_total(cwt, cwtt, &mut carrayold);
            totalold = total;
            temp = (totalold - nfile as f64) / (nnew * nnew) as f64;
        }
    }
    for j1 in 1..=npairs {
        for i in 1..=nfile {
            cwtt[(j1, i)] = cwtt[(j1, i)] * feedback as f64 * efacm[i] as f64;
        }
        for i in nfile + 1..=nnew {
            cwtt[(j1, i)] *= feedback1 as f64;
        }
    }
}
