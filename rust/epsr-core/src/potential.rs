//! Reference and empirical potentials — the heart of EPSR.
//!
//! Translation of (`epsr_standalone_rev3.f` unless noted):
//! * `initpt`            — truncation windows, r* grid, initial potentials
//! * `truncate`, `truncatec`, `chsph` (`ato_manipulation_rev4.f`)
//! * `generate_mp`       — molecular (intramolecular non-bonded) potential (`ato_manipulation_rev4.f`)
//! * `generate_rp`       — reference potential: Lennard-Jones (n-6) + Coulomb + repulsive core
//! * `generate_ep`       — empirical potential from the coefficients `pcoefp`
//! * `calc_pwrexpr_grad` — r-space basis functions (power-exponential "Poisson" functions)
//! * `calc_pwrexpq`      — their Fourier transforms
//! * `dnext`, `expand`, `fit_expansion`, `smoop` (`auxiliary_routines.f`)
//! * `ptgnsq`            — fit the data–simulation difference in Q and add it to the EP
//! * `repulsivept`, `multiply_pot_coeffs`, `lsparabola1`
//! * module `lngamma_mod` (`fctrln`)
//!
//! ### How the EP is refined (`ptgnsq`)
//! For each dataset the difference between the measured structure factor
//! and the simulated one is fitted (by a Monte Carlo search, `fit_expansion`)
//! with a sum of Q-space basis functions `expsq`.  The fitted coefficients
//! are distributed over the atom pairs with the inverted weights matrix
//! `cwtpot` and added to the pair coefficients `pcoefp`, scaled by the
//! feedback factor `potfac`.  `generate_ep` then rebuilds the r-space
//! empirical potential `xavs` (and its gradient `xavs1`) from `pcoefp`,
//! rescaling the coefficients so the potential's amplitude does not exceed
//! the requested absolute energy `ereq`.

use crate::error::Result;
use crate::fortran::fmath::{pow4, powi_var, sign};
use crate::fortran::format::{fmt, FV};
use crate::fortran::listread::{IoRes, TextFile};
use crate::fortran::output::OutFile;
use crate::fortran::{F1, F2};
use crate::state::{State, MFILE};
use rayon::prelude::*;

/// `function truncate(xx, rminpt, rmaxpt)`: smooth cosine window 1 → 0 between rminpt and rmaxpt.
pub fn truncate(xx: f32, rminpt: f32, rmaxpt: f32) -> f32 {
    let pi = std::f32::consts::PI;
    let decay = rmaxpt - rminpt;
    if xx >= rmaxpt {
        0.0
    } else if xx <= rminpt {
        1.0
    } else {
        let rat = (xx - rminpt) * pi / decay;
        0.5 * (1.0 + rat.cos())
    }
}

/// `function truncatec(xx, rmaxpt)`: polynomial window used for the Coulomb term.
pub fn truncatec(xx: f32, rmaxpt: f32) -> f32 {
    if xx >= rmaxpt {
        0.0
    } else {
        let rat = xx / rmaxpt;
        let t = pow4(1.0 - rat);
        t * (1.0 + rat * (1.6 + rat * 0.4))
    }
}

/// `function chsph(sig1, sig2, r)`: Coulomb interaction between two uniformly
/// charged spheres of radii sig1, sig2 at separation r (point charges if radii are 0).
pub fn chsph(sig1: f32, sig2: f32, r: f32) -> Result<f32> {
    let (mut ra, mut rb) = (sig1, sig2);
    if sig2 > sig1 {
        ra = sig2;
        rb = sig1;
    }
    if r >= ra + rb {
        Ok(1.0 / r)
    } else if r <= ra - rb {
        if ra > 0.0 {
            Ok((3.0 * (ra * ra - 0.2 * rb * rb) - r * r) / (2.0 * ra * ra * ra))
        } else {
            Ok(1.0e20)
        }
    } else {
        let aup1 = ra - r;
        let alo1 = -rb;
        let aup2 = rb;
        let alo2 = ra - r;
        Ok(sum1(r, ra, rb, aup1, alo1) + sum2(r, ra, rb, aup2, alo2))
    }
}

/// `function sum1(r,a,b,aup,alo)` — part of the overlapping-sphere Coulomb integral.
fn sum1(r: f32, a: f32, b: f32, aup: f32, alo: f32) -> f32 {
    let a2 = a * a;
    let b2 = b * b;
    let r2 = r * r;
    let mut cof = [0f32; 7];
    cof[6] = 1.0 / 6.0;
    cof[5] = 3.0 * r / 5.0;
    cof[4] = -(3.0 * (a2 - r2) + b2) / 4.0;
    cof[3] = -r * (3.0 * (a2 + b2) - r2) / 3.0;
    cof[2] = 3.0 * b2 * (a2 - r2) / 2.0;
    cof[1] = r * b2 * (3.0 * a2 - r2);
    let fac = 3.0 / (8.0 * a2 * a * b2 * b * r);
    let mut sumu: f32 = 0.0;
    let mut suml: f32 = 0.0;
    for i in 1..=6 {
        let ir = 7 - i;
        sumu = aup * (sumu + cof[ir]);
        suml = alo * (suml + cof[ir]);
    }
    fac * (sumu - suml)
}

/// `function sum2(r,a,b,aup,alo)`.
fn sum2(r: f32, _a: f32, b: f32, aup: f32, alo: f32) -> f32 {
    let b2 = b * b;
    let mut cof = [0f32; 4];
    cof[3] = -1.0 / 3.0;
    cof[2] = 0.0;
    cof[1] = b2;
    let fac = 3.0 / (4.0 * r * b2 * b);
    let mut sumu: f32 = 0.0;
    let mut suml: f32 = 0.0;
    for i in 1..=3 {
        let ir = 4 - i;
        sumu = aup * (sumu + cof[ir]);
        suml = alo * (suml + cof[ir]);
    }
    fac * (sumu - suml)
}

/// `function gauss6pot(r,sig,sigg)`: "modified Morse" additional potential
/// (Gaussian core minus a dispersion term from overlapping spheres).
fn gauss6pot(r: f32, sig: f32, sigg: f32) -> Result<f32> {
    let x = r / sig;
    let coul = chsph(0.5, 0.5, x)?;
    let c3 = coul * coul * coul;
    let london = -(c3 * c3);
    let sigg = if sigg > 0.29 { 0.29 } else { sigg };
    let ampgln = 0.5 / (sigg * sigg);
    let mut arg = x / sigg;
    arg = 0.5 * arg * arg - ampgln;
    let term = if arg < 30.0 { (-arg).exp() } else { 0.0 };
    Ok(term + london)
}

/// `function dnext(d, dfac, dstep, rootd)`: next basis-function position;
/// spacing grows geometrically (factor `dfac`) once that exceeds `dstep`.
fn dnext(d: f32, dfac: f32, dstep: f32, rootd: &mut f32) -> f32 {
    let mut dn = d + dstep;
    let dtest = d * dfac;
    if dtest > dn {
        dn = dtest;
        *rootd *= dfac.sqrt();
    }
    dn
}

/// `function expand(nterms, j, expr, coeff, summax)`: Σ_i expr((j-1)*nterms+i) * coeff(i).
#[inline(always)]
fn expand(nterms: i32, j: i32, expr: &F1<f32>, coeff: &F1<f32>, summax: f32) -> f32 {
    let n = nterms as usize;
    let base = ((j - 1) * nterms) as usize;
    let e = &expr.as_slice()[base..base + n];
    let c = &coeff.as_slice()[..n];
    let mut sum: f32 = 0.0;
    for k in 0..n {
        sum += e[k] * c[k];
        if sum.abs() > summax {
            sum = sign(summax, sum);
        }
    }
    sum
}

/// `SUBROUTINE SMOOp(NQMIN, NQMAX, NSMOO, DIFF)`: 1-2-1 smoothing.
fn smoop(nqmin: i32, nqmax: i32, nsmoo: i32, diff: &mut F1<f32>) {
    let n1 = nqmin;
    let n2 = nqmax;
    if n2 - 2 < n1 || nsmoo == 0 {
        return;
    }
    for _it in 1..=nsmoo {
        let mut sm: f32;
        let mut sn = diff[n1];
        sm = 0.0;
        let mut sum1 = 0.25 * (sm + sn);
        for i in n1..=n2 {
            sm = sn;
            let j = i + 1;
            sn = if j <= n2 { diff[j] } else { 0.0 };
            let sum2 = 0.25 * (sm + sn);
            diff[i] = sum1 + sum2;
            sum1 = sum2;
        }
    }
}

/// `module lngamma_mod`, `fctrln(nterms)`: table of ln((i-1)!) for i = 1..nterms+1.
fn fctrln(nterms: i32) -> Vec<f64> {
    let mgamma = (nterms + 1) as usize;
    let mut gamln = vec![0f64; mgamma];
    let mut g: f64 = 0.0;
    let mut a: f64 = 0.0;
    for v in gamln.iter_mut() {
        *v = g;
        a += 1.0;
        g += a.ln();
    }
    gamln
}

/// `subroutine calc_pwrexpr_grad(...)`: r-space EP basis functions (gradient
/// form): Poisson-like functions r^n exp(-r/σ) with widths growing with
/// position.  Fills `expr((i-1)*nterms+m)`.
#[allow(clippy::too_many_arguments)]
fn calc_pwrexpr_grad(lpt: i32, nterms: i32, dmin: f32, dstep: f32, dmax: f32, sigmab: f32, rho: f32, xin: &F1<f32>, expr: &mut F1<f32>, rmin: f32, fbroad: f32, rbroad: f32) {
    let pi: f32 = 3.141592653;
    let expmax: f32 = 25.0;
    let ntmax = (dmax / sigmab + 3.0) as i32;
    let gamln = fctrln(ntmax);
    let mut d = dstep;
    let mut rootd: f32 = 1.0;
    for m in 1..=nterms {
        let sigma = rootd * fbroad * sigmab;
        let n: i32 = if sigma <= d && sigma > 0.0 && m > 1 { (d / sigma - 1.0) as i32 } else { 0 };
        let nref = n + 3;
        let gfac = 4.0 * pi * rho * sigma * sigma * sigma;
        let gfacln = gfac.ln();
        let mut sigln = (gamln[(nref - 1) as usize] + gfacln as f64) as f32;
        if m == 1 {
            sigln -= rmin / sigma;
        }
        for i in 1..=lpt {
            let jref = (i - 1) * nterms + m;
            let r = xin[i];
            let mut fact = rbroad * r;
            let rsigma = r / sigma;
            let expon: f32;
            if r > 0.0 {
                let rsigmaln = rsigma.ln();
                expon = n as f32 * rsigmaln - rsigma - sigln - fact;
                fact = n as f32 / r - 1.0 / sigma - rbroad;
            } else if n == 0 {
                expon = -rsigma - sigln - fact;
                fact = -(1.0 / sigma + rbroad);
            } else if n == 1 {
                expon = -rsigma - sigln - fact;
                fact = 1.0 / sigma;
            } else {
                expon = -2.0 * expmax;
                fact = 0.0;
            }
            expr[jref] = if expon > -expmax { fact * expon.exp() } else { 0.0 };
        }
        d = dnext(d, dmin, dstep, &mut rootd);
    }
}

/// `subroutine calc_gaussq(lpt,nterms,dstep,dmax,sigma,rho,xin,expr,rmin,ncfirst)`:
/// Q-space transforms of the Gaussian EP basis (`Gaussian T`).  `rmin` is
/// always 0 at the call site, so `ncfirst` is never changed (omitted).
fn calc_gaussq(lpt: i32, nterms: i32, dstep: f32, sigma: f32, xin: &F1<f32>, expr: &mut F1<f32>) {
    let expmax: f32 = 25.0;
    for m in 1..=nterms {
        let d = m as f32 * dstep;
        let sigd = sigma;
        for i in 1..=lpt {
            let jref = (i - 1) * nterms + m;
            let q = xin[i];
            let qsigd = q * sigd;
            let qs2 = 0.5 * qsigd * qsigd;
            if qs2 < expmax {
                let qd = q * d;
                let sinqd = if qd > 0.0 { qd.sin() / qd } else { 1.0 };
                expr[jref] = sinqd * (-qs2).exp();
            } else {
                expr[jref] = 0.0;
            }
        }
    }
}

/// `subroutine calc_gaussr_grad(lpt,nterms,dstep,dmax,sigmab,rho,xin,expr,rmin,ncfirst,fbroad,rbroad)`:
/// r-space gradients of the Gaussian EP basis functions.  `rmin` = 0 at the
/// call site (`ncfirst` unchanged); `dmax` and `rbroad` are unused by the Fortran.
#[allow(clippy::too_many_arguments)]
fn calc_gaussr_grad(lpt: i32, nterms: i32, dstep: f32, sigmab: f32, rho: f32, xin: &F1<f32>, expr: &mut F1<f32>, fbroad: f32) {
    let pi = crate::fortran::fmath::PI4ATAN;
    let piroot = (0.5 * pi).sqrt();
    let gfac = piroot / (4.0 * pi * pi * rho);
    let expmax: f32 = 25.0;
    for m in 1..=nterms {
        let d = m as f32 * dstep;
        let sigd = fbroad * sigmab;
        let gfacs = gfac / sigd;
        let sigd2 = sigd * sigd;
        for i in 1..=lpt {
            let jref = (i - 1) * nterms + m;
            let r = xin[i];
            let rmd = (r - d) / sigd;
            let rpd = (r + d) / sigd;
            let rmd2 = 0.5 * rmd * rmd;
            let rpd2 = 0.5 * rpd * rpd;
            let mut sum: f32 = 0.0;
            let gfact;
            if r.abs() > 0.0 && d.abs() > 0.0 {
                gfact = gfacs / (r * d);
                if rmd2 < expmax {
                    sum -= (1.0 / r + (r - d) / sigd2) * (-rmd2).exp();
                }
                if rpd2 < expmax {
                    sum += (1.0 / r + (r + d) / sigd2) * (-rpd2).exp();
                }
            } else {
                gfact = 0.0;
            }
            expr[jref] = gfact * sum;
        }
    }
}

/// `subroutine calc_pwrexpq(...)`: Q-space transforms of the EP basis functions.
#[allow(clippy::too_many_arguments)]
fn calc_pwrexpq(lpt: i32, nterms: i32, dmin: f32, dstep: f32, sigmab: f32, xin: &F1<f32>, expr: &mut F1<f32>) {
    let mut d = dstep;
    let mut rootd: f32 = 1.0;
    for m in 1..=nterms {
        let sigma = rootd * sigmab;
        let n: i32 = if sigma <= d && sigma > 0.0 && m > 1 { (d / sigma - 1.0) as i32 } else { 0 };
        let an2 = (n + 2) as f32;
        let an2fac = 1.0 / an2;
        for i in 1..=lpt {
            let q = xin[i];
            let qsigma = q * sigma;
            let qs2 = qsigma * qsigma;
            let onemqs2 = 1.0 - qs2;
            let onepqs2 = 1.0 + qs2;
            let rtonepqs2 = onepqs2.sqrt();
            let alpha = qsigma.atan();
            let alphan = n as f32 * alpha;
            let term = if qsigma > 0.0 { 2.0 * alphan.cos() + onemqs2 * alphan.sin() / qsigma } else { an2 };
            let pn = term * an2fac / powi_var(rtonepqs2, 4 + n);
            let jref = (i - 1) * nterms + m;
            expr[jref] = pn;
        }
        d = dnext(d, dmin, dstep, &mut rootd);
    }
}

/// Result of `fit_expansion`.
struct FitResult {
    chisqnew: f32,
}

impl State {
    /// `subroutine fit_expansion(yin,expr,yout,wt,wtold,stepr,acceptr,rejrat,chisqnew,nmin,npt,ncfirst,nterms,nits,nfreq)`:
    /// Monte Carlo fit of coefficients `wt` so that Σ wt·expr ≈ yin.
    #[allow(clippy::too_many_arguments)]
    fn fit_expansion(
        &mut self,
        yin: &F1<f32>,
        expr: &F1<f32>,
        yout: &mut F1<f32>,
        wt: &mut F1<f32>,
        wtold: &mut F1<f32>,
        stepr: &mut f32,
        acceptr: &mut f32,
        nmin: i32,
        npt: i32,
        ncfirst: i32,
        nterms: i32,
        nits: i32,
        nfreq: i32,
    ) -> FitResult {
        let summax = 25.0f32.exp();
        let nstep = 100;
        let ntry = 150;
        *acceptr = 0.05;
        let npoints = npt - nmin + 1;
        let mut wtmax: f32 = 0.0;
        for it in 1..=nterms {
            wt[it] = wtold[it];
            wtmax = wt[it].abs().max(wtmax);
        }
        if *stepr < 1.0e-6 * wtmax {
            *stepr = 2.0e-6 * wtmax;
        }
        let mut chisqnew: f32;
        // Chi-squared over the fit range.  The expand() rows are independent
        // (the running clamp is within one row), so they are evaluated on the
        // rayon pool; the squared residuals are then summed serially in index
        // order, so the result is identical at any thread count.
        let chisq = |wt: &F1<f32>| -> f32 {
            let terms: Vec<f32> = (nmin..npt + 1)
                .into_par_iter()
                .with_min_len(64)
                .map(|i| {
                    let term = yin[i] - expand(nterms, i, expr, wt, summax);
                    term * term
                })
                .collect();
            let mut c: f32 = 0.0;
            for t in &terms {
                c += *t;
            }
            c
        };
        if nits > 0 {
            let mut icount = 0;
            chisqnew = chisq(wt);
            let mut chisqold = chisqnew;
            let mut chisqtest: f32 = -1.0;
            let mut itry = 0;
            let mut chisqfinish: f32 = 0.0;
            let mut rejrat: f32 = 0.0;
            while itry < 10 || (itry < ntry && chisqtest < 1.0e-4) {
                itry += 1;
                let accept = *acceptr * chisqold / npoints as f32;
                rejrat = 0.0;
                let chisqstart = chisqold / npoints as f32;
                for _it in 1..=nstep {
                    for ic in ncfirst..=nterms {
                        if wtold[ic].abs() < 1.0 {
                            wt[ic] = wtold[ic] + *stepr * (2.0 * self.rng.ran1() - 1.0);
                        } else {
                            wt[ic] = wtold[ic] * (1.0 + *stepr * (2.0 * self.rng.ran1() - 1.0));
                        }
                    }
                    chisqnew = chisq(wt);
                    let exponent = (chisqnew - chisqold) / accept;
                    let test = if exponent > 0.0 && exponent <= 20.0 {
                        (-exponent).exp() - self.rng.ran1()
                    } else if exponent > 20.0 {
                        -1.0
                    } else {
                        1.0
                    };
                    if test >= 0.0 {
                        icount += 1;
                        if icount == nfreq {
                            icount = 0;
                            smoop(ncfirst, nterms, 1, wt);
                            chisqnew = chisq(wt);
                        }
                        chisqold = chisqnew;
                        for ic in 1..=nterms {
                            wtold[ic] = wt[ic];
                        }
                    } else {
                        rejrat += 1.0;
                    }
                }
                rejrat /= nstep as f32;
                chisqfinish = chisqold / npoints as f32;
                chisqtest = if chisqstart > 0.0 { (chisqfinish - chisqstart) / chisqstart } else { chisqfinish - chisqstart };
                let rejfac = if rejrat > 0.0 { 0.75 / rejrat } else { 2.0 };
                *stepr *= rejfac;
            }
            self.log(&format!(" fit_expansion> {itry:>11} {chisqfinish} {chisqtest} {rejrat}"));
            for ic in 1..=nterms {
                wt[ic] = wtold[ic];
            }
        }
        chisqnew = 0.0;
        for i in 1..=npt {
            let term = expand(nterms, i, expr, wt, summax);
            yout[i] = term;
            if i >= nmin {
                let dd = yin[i] - term;
                chisqnew += dd * dd;
            }
        }
        FitResult { chisqnew }
    }

    /// `subroutine initpt(ireset,rminpt,rmaxpt,expecf,polar,smoos,rminp,rbroad)`.
    pub fn initpt(&mut self) -> Result<()> {
        let rminpt = self.inp.rminpt;
        let rmaxpt = self.inp.rmaxpt;
        self.log(&format!(" initpt> rminpt,rmaxpt: {rminpt} {rmaxpt}"));
        let d = &mut self.dist;
        let p = &mut self.pot;
        let dmax = self.mol.dmax;
        for j in 1..=d.npt {
            d.xstar[j] = d.x[j] / dmax;
        }
        for j in 1..=d.npt {
            let xx = d.x[j];
            p.wfac[j] = truncate(xx, rminpt, rmaxpt);
            p.wfacc[j] = truncatec(xx, rmaxpt);
        }
        p.wfac[d.npt] = 0.0;
        p.wfacc[d.npt] = 0.0;
        let expecf = self.inp.expecf;
        for ip in 1..=self.mol.npairs {
            let irefp = p.irefpcof[ip];
            for j in 1..=d.npt {
                if d.x[j] <= (p.rminex[irefp] + 4.0 * expecf).abs() && d.x[j] > d.x[d.nptmax] {
                    d.nptmax = j;
                }
                p.xav[(j, ip)] = 0.0;
                p.xavs[(j, ip)] = 0.0;
                p.xav1[(j, ip)] = 0.0;
                p.xavs1[(j, ip)] = 0.0;
            }
        }
        self.generate_rp()?;
        self.generate_ep(false)?;
        Ok(())
    }

    /// `subroutine generate_mp()`: intramolecular non-bonded potential
    /// (WCA-truncated n-6 core with parameters `ecore`, `dcore`).
    pub fn generate_mp(&mut self) {
        let m = &self.mol;
        let d = &self.dist;
        let p = &mut self.pot;
        let power = m.power;
        let rat = 6.0 / power;
        let rat1 = rat.sqrt();
        let rat2 = 1.0 / rat1;
        let power1 = (power + 6.0) / (power - 6.0);
        let rat3 = rat1.powf(power1);
        let mut prefac = (rat2 - rat1) * rat3;
        prefac = 1.0 / prefac;
        let rat4 = 1.0 / (power - 6.0);
        let rmaxwca = (m.dcore / m.sizefactor) * (power / 6.0).powf(rat4);
        for ic in 1..=m.npairs {
            let ect = m.ecore;
            let dc = m.dcore / m.sizefactor;
            for j in 2..=d.npt {
                let x2 = d.x[j];
                let ddif = dc / x2;
                let sumlj = if ddif > 0.0 && ect > 0.0 && x2 < rmaxwca {
                    let fac = prefac * ect;
                    let xx = ddif.ln();
                    let mut pwr1 = 0.5 * (power + 6.0) * xx;
                    if pwr1.abs() > 20.0 {
                        pwr1 = sign(20.0, pwr1);
                    }
                    let mut pwr2 = 0.5 * (power - 6.0) * xx;
                    if pwr2.abs() > 40.0 {
                        pwr2 = sign(40.0, pwr2);
                    }
                    2.0 * fac * pwr2.sinh() * pwr1.exp() + ect
                } else {
                    0.0
                };
                p.xav[(j, ic)] = sumlj * p.wfac[j];
            }
            p.xav[(1, ic)] = p.xav[(2, ic)];
            let mut sumlj1 = p.xav[(1, ic)];
            let mut x1 = d.xstar[1];
            for j in 2..=d.npt {
                let x2 = d.xstar[j];
                let xdif = x2 - x1;
                let sumlj2 = p.xav[(j, ic)];
                p.xav1[(j - 1, ic)] = (sumlj2 - sumlj1) / xdif;
                x1 = x2;
                sumlj1 = sumlj2;
            }
            p.xav1[(d.npt, ic)] = 0.0;
        }
        // get_containerforce_list(nmol,ic,rmaxex(irefp),it,jt,...): in the Fortran
        // `it`/`jt` are never set in generate_mp (uninitialised locals); 0 is used
        // here, so no atoms match and every list is empty.
        for ic in 1..=self.mol.npairs {
            let rmaxex = self.pot.rmaxex[self.pot.irefpcof[ic]];
            self.get_containerforce_list(ic, rmaxex, 0, 0);
        }
    }

    /// `subroutine generate_rp()`: reference potential (LJ n-6 + Coulomb
    /// [+ exponential repulsive core where `pcoefex` ≠ 0]) and gradients.
    pub fn generate_rp(&mut self) -> Result<()> {
        let mut gradlimit = self.inp.expecfpwr;
        let loggradlimit;
        if gradlimit > 0.0 {
            loggradlimit = gradlimit.ln();
        } else {
            gradlimit = 0.0;
            loggradlimit = 0.0;
        }
        self.generate_mp();
        let m = &self.mol;
        let d = &self.dist;
        let expecf = self.inp.expecf;
        let power = m.power;
        let rat = 6.0 / power;
        let rat1 = rat.sqrt();
        let rat2 = 1.0 / rat1;
        let power1 = (power + 6.0) / (power - 6.0);
        let rat3 = rat1.powf(power1);
        let mut prefac = (rat2 - rat1) * rat3;
        prefac = 1.0 / prefac;
        for ic in 1..=m.npairs {
            let p = &mut self.pot;
            let irefp = p.irefpcof[ic];
            let it = m.n1extr[ic];
            let jt = m.n2extr[ic];
            let ra = m.qradius[it];
            let rb = m.qradius[jt];
            let ect = m.ecor[ic];
            let qct = m.qcor[ic] * p.coulfac[irefp];
            let dc = m.dcor[ic];
            for j in 2..=d.npt {
                let x2 = d.x[j];
                let ddif = dc / x2;
                let mut sumlj = if ddif > 0.0 && ect > 0.0 {
                    let fac = prefac * ect;
                    let xx = ddif.ln();
                    let mut pwr1 = 0.5 * (power + 6.0) * xx;
                    if pwr1.abs() > 20.0 {
                        pwr1 = sign(20.0, pwr1);
                    }
                    let mut pwr2 = 0.5 * (power - 6.0) * xx;
                    if pwr2.abs() > 40.0 {
                        pwr2 = sign(40.0, pwr2);
                    }
                    2.0 * fac * pwr2.sinh() * pwr1.exp()
                } else {
                    0.0
                };
                let xc = x2;
                let mut sumc = qct * chsph(ra, rb, xc)?;
                sumlj = p.refpotfac * sumlj * p.wfac[j];
                sumc = p.refpotfac * sumc * p.wfacc[j];
                p.epot[(j, ic)] = sumlj + sumc;
                p.ecoul[(j, ic)] = sumlj + sumc;
                // additional Gaussians (the first may be a modified Morse function)
                let mut ig = 0;
                while ig < p.npgauss[irefp] {
                    ig += 1;
                    let wpg = p.wpgauss[(ig, irefp)].abs();
                    if wpg > 0.0 {
                        let rpg = p.rpgauss[(ig, irefp)].abs();
                        let hpg = p.hpgauss[(ig, irefp)];
                        let mut add: f32 = 0.0;
                        if ig == 1 && p.addpotential == "modmorse" {
                            add = hpg.abs() * gauss6pot(d.x[j], rpg, wpg)?;
                        } else {
                            let rdiff = (d.x[j] - rpg) / wpg;
                            let term = hpg.abs().ln() - 0.5 * rdiff * rdiff;
                            if term.abs() < 25.0 {
                                add = sign(1.0, hpg) * term.exp();
                            }
                        }
                        p.ecoul[(j, ic)] += add * p.wfac[j];
                    }
                }
                if m.coul_intra {
                    p.xav[(j, ic)] = m.coul_intra_fac * p.ecoul[(j, ic)];
                }
            }
            // Extra repulsive core where intensity was found below rminex.
            // Exercised by EPSR-test-ex1 (rminex > 0 on three pairs).
            if p.pcoefex[irefp].abs() > 0.0 && p.rminex[irefp].abs() > 0.0 {
                let expecfuse = expecf * p.rminex[irefp].abs().sqrt();
                if !p.harmonic {
                    let rminextest = p.rminex[irefp].abs();
                    for j in 1..=d.npt2 {
                        let ddiff = (rminextest - d.x[j]) / expecfuse;
                        let add = if gradlimit > 0.0 && ddiff >= loggradlimit {
                            p.pcoefex[irefp] * gradlimit * (1.0 + ddiff - loggradlimit)
                        } else {
                            p.pcoefex[irefp] * ddiff.exp()
                        };
                        p.ecoul[(j, ic)] += add * p.wfac[j];
                    }
                } else {
                    // Gaussian (harmonic) repulsive core
                    let rminextest = p.rminex[irefp].abs() + expecfuse;
                    for j in 1..=d.npt2 {
                        let ddiff = (rminextest - d.x[j]) / expecfuse;
                        if ddiff > 0.0 {
                            let add = if gradlimit > 0.0 && ddiff >= gradlimit {
                                0.5 * p.pcoefex[irefp] * gradlimit * (2.0 * ddiff - gradlimit)
                            } else {
                                0.5 * p.pcoefex[irefp] * ddiff * ddiff
                            };
                            p.ecoul[(j, ic)] += add * p.wfac[j];
                        }
                    }
                }
            }
            p.epot[(1, ic)] = p.epot[(2, ic)];
            p.ecoul[(1, ic)] = p.ecoul[(2, ic)];
            p.xav[(1, ic)] = p.xav[(2, ic)];
            let mut sumlj1 = p.epot[(1, ic)];
            let mut sumc1 = p.ecoul[(1, ic)];
            let mut sumx1 = p.xav[(1, ic)];
            let mut x1 = d.xstar[1];
            for j in 2..=d.npt {
                let x2 = d.xstar[j];
                let xdif = x2 - x1;
                let sumlj2 = p.epot[(j, ic)];
                let sumc2 = p.ecoul[(j, ic)];
                let sumx2 = p.xav[(j, ic)];
                p.epot1[(j - 1, ic)] = (sumlj2 - sumlj1) / xdif;
                p.ecoul1[(j - 1, ic)] = (sumc2 - sumc1) / xdif;
                p.xav1[(j - 1, ic)] = (sumx2 - sumx1) / xdif;
                x1 = x2;
                sumlj1 = sumlj2;
                sumc1 = sumc2;
                sumx1 = sumx2;
            }
            p.epot1[(d.npt, ic)] = 0.0;
            p.ecoul1[(d.npt, ic)] = 0.0;
            p.xav[(d.npt, ic)] = 0.0;
        }
        // get_containerforce_list is called at the end of each pair; it only
        // reads rmaxex and the atom types, so calling it afterwards is equivalent.
        for ic in 1..=self.mol.npairs {
            let rmaxex = self.pot.rmaxex[self.pot.irefpcof[ic]];
            let (it, jt) = (self.mol.n1extr[ic], self.mol.n2extr[ic]);
            self.get_containerforce_list(ic, rmaxex, it, jt);
        }
        Ok(())
    }

    /// `subroutine generate_ep(changepotfac)`: rebuild the empirical potential
    /// `xavs` (and gradient `xavs1`, total `xavtot`) from the coefficients.
    pub fn generate_ep(&mut self, changepotfac: bool) -> Result<()> {
        let summax = 25.0f32.exp();
        let nterms = self.pot.ncoeffp;
        let decayp = self.pot.psigma2;
        let npt = self.dist.npt;
        let npt2 = self.dist.npt2;
        self.pot.alloc_basis(npt.max(1) as usize + 1, nterms as usize);
        {
            let p = &mut self.pot;
            let (pdmin, pdstep, pdmax) = (p.pdmin, p.pdstep, p.pdmax);
            if p.gaussian {
                calc_gaussr_grad(npt2, nterms, pdstep, decayp, self.dist.rho, &self.dist.x, &mut p.expgr, 1.0);
            } else {
                calc_pwrexpr_grad(npt2, nterms, pdmin, pdstep, pdmax, decayp, self.dist.rho, &self.dist.x, &mut p.expgr, 0.0, 1.0, self.inp.rbroad);
            }
        }
        let npairs = self.mol.npairs;
        for ip in 1..=npairs {
            let p = &mut self.pot;
            let irefp = p.irefpcof[ip];
            for ic in 2..=nterms {
                p.pwtold[ic] = p.pcoefp[(ic, irefp)];
            }
            p.pwtold[1] = 0.0;
            for j in 1..=npt2 {
                let term = expand(nterms, j, &p.expgr, &p.pwtold, summax);
                p.xavs1[(j, ip)] = term * p.pcoefp[(1, irefp)];
            }
        }
        let dmax = self.mol.dmax;
        for ip in 1..=npairs {
            let p = &mut self.pot;
            let xstar = &self.dist.xstar;
            for j in 1..=npt2 {
                p.xavs1[(j, ip)] *= dmax;
            }
            p.xavs[(npt, ip)] = 0.0;
            p.xavs1[(npt, ip)] = 0.0;
            let mut x1 = xstar[npt];
            let mut xts1 = p.xavs1[(npt, ip)];
            for j in 1..=npt2 {
                let jref = npt - j;
                let x2 = xstar[jref];
                let xts2 = p.xavs1[(jref, ip)];
                p.xavs[(jref, ip)] = p.xavs[(jref + 1, ip)] - 0.5 * (x1 - x2) * (xts1 + xts2);
                x1 = x2;
                xts1 = xts2;
            }
            for j in 1..=npt2 {
                p.xavs[(j, ip)] *= p.wfac[j];
            }
        }
        let mut energabs = self.abs_energy_ep();
        let mut pressfac: f32 = 1.0;
        let p = &self.pot;
        if p.absolute_energy.abs() == 0.0 {
            pressfac = 0.0;
            energabs = 0.0;
        } else if energabs.abs() > 0.0 {
            if energabs.abs() < p.erequnit.abs() {
                pressfac = p.erequnit.abs() / energabs.abs();
            } else {
                pressfac = p.absolute_energy.abs() / energabs.abs();
                if pressfac > 1.0 && p.ereqstep == 0.0 {
                    pressfac = 1.0;
                }
            }
        }
        self.log(&format!(" generate_ep> {} {} {}", p.absolute_energy, energabs, pressfac));
        self.multiply_pot_coeffs(pressfac);
        for ip in 1..=npairs {
            let p = &mut self.pot;
            let xstar = &self.dist.xstar;
            let _irefp = p.irefpcof[ip];
            for j in 1..=npt2 {
                p.xavs[(j, ip)] *= pressfac;
                p.xavtot[(j, ip)] = p.ecoul[(j, ip)] + p.xavs[(j, ip)];
            }
            let mut sums1 = p.xavs[(1, ip)];
            let mut sumtot1 = p.xavtot[(1, ip)];
            let mut x1 = xstar[1];
            for j in 2..=npt {
                let x2 = xstar[j];
                let xdif = x2 - x1;
                let sums2 = p.xavs[(j, ip)];
                let sumtot2 = p.xavtot[(j, ip)];
                p.xavs1[(j - 1, ip)] = (sums2 - sums1) / xdif;
                p.xavtot1[(j - 1, ip)] = (sumtot2 - sumtot1) / xdif;
                x1 = x2;
                sums1 = sums2;
                sumtot1 = sumtot2;
            }
            p.xavs1[(npt, ip)] = 0.0;
            p.xavtot1[(npt, ip)] = 0.0;
        }
        // NB: generate_ep is only ever called with changepotfac = .false. in the
        // Fortran (initpt and ptgnsq); this branch is translated but unreachable.
        if changepotfac && self.pot.absolute_energy > 0.0 {
            let p = &mut self.pot;
            if p.potfac > 0.000001 && p.potfacamp > 0.0 {
                let mut absediff = (p.absolute_energy - energabs) / (p.potfacamp * p.absolute_energy).abs();
                let msg = format!(" changepotfac> {} {} {} {} {}", p.potfacamp, p.potfac, p.absolute_energy, energabs, absediff);
                if absediff.abs() > 0.1 {
                    absediff = sign(0.1, absediff);
                }
                p.potfac *= 1.0 + absediff;
                p.potfac = p.potfac.max(0.000001);
                p.potfac = p.potfac.min(1.0);
                self.log(&msg);
            }
        }
        Ok(())
    }

    /// `subroutine multiply_pot_coeffs(epfact)`.
    pub fn multiply_pot_coeffs(&mut self, epfact: f32) {
        let p = &mut self.pot;
        let nterms = p.ncoeffp;
        for ip in 1..=self.mol.npairs {
            let irefp = p.irefpcof[ip];
            for ic in 2..=nterms {
                p.pcoefp[(ic, irefp)] *= epfact;
            }
        }
    }

    /// `subroutine repulsivept(potadd)`: adjust the repulsive-core amplitudes
    /// `pcoefex` towards the requested close-contact numbers `coordex`, then
    /// regenerate the reference potential.
    pub fn repulsivept(&mut self, potadd: f32) -> Result<()> {
        if potadd > 0.0 {
            let simt = self.mol.simt;
            for ip in 1..=self.mol.npairs {
                let p = &mut self.pot;
                let irefp = p.irefpcof[ip];
                let msg = format!(" repulsivept potex> {} {}", p.potex[ip], p.coordex[irefp]);
                if p.rminex[irefp] > 0.0 {
                    // Exercised by EPSR-test-ex1 (rminex > 0 on three pairs)
                    let mut pchange = p.potex[ip] - p.coordex[irefp];
                    if pchange == 0.0 && p.pcoefex[irefp].abs() > 0.0001 {
                        p.potex[ip] = -0.1;
                        p.pcoefex[irefp] *= 1.0 + p.potex[ip];
                    } else {
                        if p.coordex[irefp] > 1.0 {
                            pchange /= p.coordex[irefp];
                        }
                        pchange = sign((0.005 * simt).min(0.025 * simt * pchange.abs()), pchange);
                        p.pcoefex[irefp] += pchange;
                    }
                }
                p.pcoefexmax[irefp] = p.potexmax[ip];
                self.log(&msg);
            }
            self.generate_rp()?;
        }
        Ok(())
    }

    /// `subroutine ptgnsq(...)`: refine the empirical potential from the
    /// differences `data` (data minus simulated intramolecular part) and
    /// `qwind` (simulated intermolecular partials).  Writes `delfit.dat`.
    #[allow(clippy::too_many_arguments)]
    pub fn ptgnsq(&mut self, nq: i32, q: &F1<f32>, data: &F2<f32>, qwind: &F2<f32>, szeros: &F1<f32>, rminp: &F1<f32>, potnoi: &mut F1<f32>) -> Result<()> {
        let nfile = self.dist.nfile;
        let qminptemp = if self.dist.qminp <= 0.0 { self.dist.qminm } else { self.dist.qminp };
        let potadd = self.pot.potfac;
        let mut qtemp: F2<f32> = F2::new(nq.max(1) as usize, MFILE);
        let mut fittemp: F2<f32> = F2::new(nq.max(1) as usize, MFILE);
        let mut nmini: F1<i32> = F1::new(MFILE);
        for id in 1..=nfile {
            nmini[id] = 2;
            let ix = self.wts.xfile[id];
            // Set up the gradient for extrapolating to Q=0
            let (grads0, conss0): (f32, f32) = if szeros[id] != 0.0 && self.dist.nqmin[id] > 1 {
                let nqmn = self.dist.nqmin[id];
                ((data[(nqmn, id)] - szeros[id]) / (q[nqmn] - q[1]), szeros[id])
            } else {
                (0.0, 0.0)
            };
            for i in 1..=nq {
                let mut sum: f32 = 0.0;
                for ic in 1..=self.mol.npairs {
                    let w = if ix > 0 { self.xray.weightx[(i, ic, ix)] } else { self.wts.cwt[(id, ic)] };
                    sum = (sum as f64 + w * qwind[(i, ic)] as f64) as f32;
                }
                let nqmn = self.dist.nqmin[id];
                if q[i] < qminptemp.max(q[nqmn]) {
                    qtemp[(i, id)] = if szeros[id] != 0.0 { conss0 + grads0 * (q[i] - q[1]) - sum } else { 0.0 };
                } else if i >= nqmn && i <= self.dist.nqmax[id] {
                    qtemp[(i, id)] = data[(i, id)] - sum;
                } else {
                    qtemp[(i, id)] = 0.0;
                }
            }
        }
        let formatout = fmt("(a,i2.2,a)", &[FV::S("("), FV::I((2 * MFILE + 1) as i64), FV::S("(1x,e12.5))")]);
        if self.pot.feedback > 0.0 && potadd > 0.0 {
            let npt2 = self.dist.npt2;
            for id in 1..=nfile {
                for j in 1..=npt2 {
                    if self.dist.x[j] < rminp[id] {
                        nmini[id] = j + 1;
                    }
                }
                self.log(&format!(" Minimum radius integer of potential: {:>5}", nmini[id]));
                let nfreq = 10;
                let nterms = self.pot.npcoeff[id];
                let mut stepp = self.pot.pstepd[id];
                let decayp = self.dist.sigma2[id];
                let mut acceptr = self.pot.paccept;
                let mut pwt = std::mem::take(&mut self.pot.pwt);
                let mut pwtold = std::mem::take(&mut self.pot.pwtold);
                for ic in 1..=nterms {
                    pwtold[ic] = self.pot.pcoefd[(ic, id)];
                    pwt[ic] = pwtold[ic];
                }
                let nqmx = self.dist.nqmax[id];
                self.pot.alloc_basis(nqmx.max(1) as usize + 1, nterms as usize);
                let mut expsq = std::mem::take(&mut self.pot.expsq);
                if self.pot.gaussian {
                    calc_gaussq(nqmx, nterms, self.pot.pdstep, decayp, q, &mut expsq);
                } else {
                    calc_pwrexpq(nqmx, nterms, self.pot.pdmin, self.pot.pdstep, decayp, q, &mut expsq);
                }
                let mut dt: F1<f32> = F1::new(nqmx.max(1) as usize);
                for i in 1..=nqmx {
                    dt[i] = -qtemp[(i, id)];
                }
                let mut dts: F1<f32> = F1::new(nqmx.max(1) as usize);
                // NB: the Fortran calls fit_expansion with identical arguments for szeros = 0 and /= 0.
                let npits = self.pot.npits;
                let fr = self.fit_expansion(&dt, &expsq, &mut dts, &mut pwt, &mut pwtold, &mut stepp, &mut acceptr, 1, nqmx, 1, nterms, npits, nfreq);
                self.pot.expsq = expsq;
                for i in 1..=nqmx {
                    fittemp[(i, id)] = dts[i];
                }
                self.pot.pstepd[id] = stepp;
                potnoi[id] = fr.chisqnew;
                let pfact = potadd;
                for ic in 2..=nterms {
                    self.pot.pcoefd[(ic, id)] = pwtold[ic];
                    for ip in 1..=self.mol.npairs {
                        let p = &mut self.pot;
                        let irefp = p.irefpcof[ip];
                        let add = ((pfact * p.efacp[irefp]) as f64 * self.wts.cwtpot[(ip, id)]) * pwtold[ic] as f64;
                        p.pcoefp[(ic, irefp)] = (p.pcoefp[(ic, irefp)] as f64 + add) as f32;
                    }
                }
                self.pot.pwt = pwt;
                self.pot.pwtold = pwtold;
                let mut rfactorep: f32 = 0.0;
                let nqmn = if szeros[id] == 0.0 { self.dist.nqmin[id] } else { 1 };
                for i in nqmn..=nq {
                    let qt = qtemp[(i, id)];
                    let ft = if i <= fittemp.n1 as i32 { fittemp[(i, id)] } else { 0.0 };
                    if qt.abs() > 0.0 {
                        rfactorep += (ft + qt).abs() / qt.abs();
                    } else {
                        rfactorep += (ft + qt).abs();
                    }
                }
                rfactorep /= (nq - nqmn + 1) as f32;
                self.log(&format!(" ptgnsq> {id:>11} {potadd} {pfact} {rfactorep}"));
            }
            // delfit.dat: fitted and actual differences for each dataset
            let path = self.sym.full("delfit.dat");
            let mut out = OutFile::create(&path);
            for i in 1..=nq {
                let mut v: Vec<FV> = vec![FV::R(q[i])];
                for id in 1..=nfile {
                    v.push(FV::R(fittemp[(i, id)]));
                    v.push(FV::R(-qtemp[(i, id)]));
                    v.push(FV::R((fittemp[(i, id)] + qtemp[(i, id)]).abs()));
                }
                out.line(&fmt(&formatout, &v));
            }
            out.close()?;
            self.generate_ep(false)?;
        }
        self.pot.npits = self.pot.npitss;
        Ok(())
    }
}

/// `subroutine lsparabola1(nchan,fname,nav,navmax,epcfrac,a,b,c,rfacwrite,aav,bav,lastqf,nextqf)`:
/// read the `.erg` history (one line per iteration), fit a straight line to
/// the quality factor and write it to the `.qdr` file.
pub struct ParabolaFit {
    pub nav: i32,
    pub a: f32,
    pub b: f32,
    pub c: f32,
    pub rfacwrite: f32,
    pub aav: f32,
    pub bav: f32,
    pub lastqf: f32,
    pub nextqf: f32,
}

pub fn lsparabola1(erg: &mut TextFile, qdr_path: &std::path::Path, navmax: i32, epcfrac: f32, log: &dyn Fn(&str)) -> Result<ParabolaFit> {
    let navmax_u = navmax.max(1) as usize;
    let mut xval = vec![0f32; navmax_u];
    let mut yval = vec![0f32; navmax_u];
    let mut aval = vec![0f32; navmax_u];
    let mut bval = vec![0f32; navmax_u];
    let mut xmax: f32 = 0.0;
    let mut nav: i32 = 0;
    let mut nsum = 0;
    let mut lastqf: f32 = 0.0;
    let mut nextqf: f32 = 0.0;
    let mut rfacwrite: f32 = 0.0;
    let mut navref = 0;
    loop {
        let rec: IoRes<[f32; 10]> = {
            let mut r = erg.list();
            (|| {
                let mut v = [0f32; 10];
                for x in v.iter_mut() {
                    *x = r.real()?.unwrap_or(0.0);
                }
                Ok(v)
            })()
        };
        let v = match rec {
            Ok(v) => v,
            Err(_) => break,
        };
        rfacwrite = v[2];
        let x = v[3];
        let mut y = v[5];
        let (aread, bread) = (v[7], v[8]);
        nsum += 1;
        y = epcfrac * y + (1.0 - epcfrac) * rfacwrite;
        let thisqf = y;
        if nsum == 1 {
            nextqf = thisqf;
        }
        if x > xmax {
            xmax = x;
        }
        if nsum > 1 {
            lastqf = nextqf;
            nextqf = thisqf;
            nav += 1;
            navref = nav - navmax * ((nav - 1) / navmax);
            let k = (navref - 1) as usize;
            xval[k] = x;
            yval[k] = y;
            aval[k] = aread;
            bval[k] = bread;
        }
    }
    log(&format!("lsparabola1> {nsum:>5} {nav:>5} {navref:>5}"));
    erg.backspace();
    let nsum = nav.min(navmax);
    let mut sum: f64 = 0.0;
    let mut sumx1: f64 = 0.0;
    let mut sumy: f64 = 0.0;
    let mut sumx2: f64 = 0.0;
    let mut sumyx1: f64 = 0.0;
    let mut aav: f32 = 0.0;
    let mut bav: f32 = 0.0;
    if nsum > 0 {
        for i in 0..nsum as usize {
            let x = xval[i];
            let y = yval[i];
            sum += 1.0;
            sumy += y as f64;
            sumx1 += x as f64;
            sumyx1 += (y * x) as f64;
            let xm = (x * x) as f64;
            sumx2 += xm;
            aav += aval[i];
            bav += bval[i].abs();
        }
        aav /= nsum as f32;
        bav /= nsum as f32;
    }
    // quadratic fitting is disabled in the Fortran (delta = 0); fit a straight line
    let (a, b, c): (f32, f32, f32);
    if nav >= 2 {
        c = 0.0;
        let delta = sumx1 * sumx1 - sum * sumx2;
        if delta.abs() > 1.0e-10 {
            let bplus = sumy * sumx1;
            let bminus = sumyx1 * sum;
            b = if bplus != bminus { ((bplus - bminus) / delta) as f32 } else { 0.0 };
        } else {
            b = 0.0;
        }
        sumy -= b as f64 * sumx1;
        a = if sum.abs() > 0.0 { (sumy / sum) as f32 } else { 0.0 };
    } else {
        c = 0.0;
        b = 0.0;
        let delta = sum;
        a = if delta.abs() > 0.0 { (sumy / delta) as f32 } else { 0.0 };
    }
    log(&format!(" lsparabola1> {nsum} 0 {a} {b} {c}"));
    let mut out = OutFile::create(qdr_path);
    let npt = 100;
    let xstep = xmax / npt as f32;
    let (emin, femin) = (0.0f32, 0.0f32);
    for ix in 1..=npt {
        let x = ix as f32 * xstep;
        out.line(&fmt("(7(1x,e13.6))", &[FV::R(a), FV::R(b), FV::R(c), FV::R(x), FV::R(a + x * (b + x * c)), FV::R(emin), FV::R(femin)]));
    }
    out.close()?;
    Ok(ParabolaFit { nav, a, b, c, rfacwrite, aav, bav, lastqf, nextqf })
}
