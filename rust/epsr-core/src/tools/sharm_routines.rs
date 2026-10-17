//! Routines shared by the spherical-harmonic / spatial-density programs
//! (`sharm`, `sdf`, `sdfcube`, `mapgr`, `mapgrentropy`).
//!
//! Translation of
//! * module `sharm_mod` and the routines of `sharm_routines_rev1.f`
//!   (`get_l_and_n_values`, `form_cgcoeffs`, `getrefint`, `getrefintnon`,
//!   `form_eulerangles`, `get_max_size_dlmn`, `form_dlmn`,
//!   `form_sum_prod_dlmn_cgcoeff`);
//! * module `lngamma_mod` (`fctrln`);
//! * `dlm`, `maxdlm`, `cgcoef`, `maxcgcoef`, `rotateeuler` (`auxiliary_routines.f`);
//! * `outgofzs_sharm` (`read_and_save_files.f`), `r_int_titles` (`auxiliary_routines.f`);
//! * `smoosharm`, `sumoverm2` (`plotsharm_routines.f`);
//! * `find_plot_par` (`plot_routines.f`) used by the `mapgr`/`mapgrentropy`
//!   branch of `check_symbol`;
//! * libgfortran's `cosd`/`sind` intrinsics (`_gfortran_cosd_r4`,
//!   `_gfortran_sind_r4`, from `libgfortran/intrinsics/trigd.inc`), which
//!   gfortran >= 14 uses for `cosd`/`sind` in `plotsharm_routines.f` (where they
//!   are not declared `external`).

use crate::ato::{scapro, vecpro};
use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::{read_int, TextFile};
use crate::fortran::output::OutFile;
use crate::fortran::text::{indexlastdot, parse};
use crate::fortran::{F1, F2, F3};
use crate::output::{r_int, w_int, NWRMAX};
use crate::state::State;

// ---------------------------------------------------------------------------
// lngamma_mod
// ---------------------------------------------------------------------------

/// `module lngamma_mod`: table of ln(n!) — `gamln(i)` = ln((i-1)!).
#[derive(Clone, Debug, Default)]
pub struct LnGamma {
    /// `gamln(1..mgamma)` stored 0-based.
    pub gamln: Vec<f64>,
}

impl LnGamma {
    /// `subroutine fctrln(nterms)`.
    pub fn fctrln(&mut self, nterms: i32) {
        let mgamma = (nterms + 1).max(0) as usize;
        let mut v = Vec::with_capacity(mgamma);
        let mut g: f64 = 0.0;
        let mut a: f64 = 0.0;
        for _ in 0..mgamma {
            v.push(g);
            a += 1.0;
            g += a.ln();
        }
        self.gamln = v;
    }
    #[inline(always)]
    pub fn mgamma(&self) -> i32 {
        self.gamln.len() as i32
    }
    /// `gamln(i)` (1-based).
    #[inline(always)]
    pub fn g(&self, i: i32) -> f64 {
        self.gamln[(i - 1) as usize]
    }
}

/// `function dlm(l,m,n,costh)`: rotation matrix element d^l_{mn}(theta).
pub fn dlm(lg: &LnGamma, l: i32, m: i32, n: i32, costh: f32) -> f32 {
    let l2 = 2 * l;
    let lpm = l + m;
    if lpm > l2 || lpm < 0 {
        return 0.0;
    }
    let lmm = l - m;
    if lmm > l2 || lmm < 0 {
        return 0.0;
    }
    let lpn = l + n;
    if lpn > l2 || lpn < 0 {
        return 0.0;
    }
    let lmn = l - n;
    if lmn > l2 || lmn < 0 {
        return 0.0;
    }
    let mmn = m - n;
    let mg = lg.mgamma();
    if lpm + 1 > mg || lmm + 1 > mg || lpn + 1 > mg || lmn + 1 > mg {
        panic!("dlm> {} {} {} {} {} {} {} (log gamma table too small)", l, m, n, lpm + 1, lmm + 1, lpn + 1, lmn + 1);
    }
    let faco: f64 = 0.5 * (lg.g(lpm + 1) + lg.g(lmm + 1) + lg.g(lpn + 1) + lg.g(lmn + 1));
    let kmax = lpm.min(lmn) + 1;
    let kmin = 0.max(mmn) + 1;
    let (dlcth2, dlsth2): (f64, f64);
    if costh >= 1.0 {
        dlcth2 = 0.0;
        dlsth2 = -1.0e30;
    } else if costh <= -1.0 {
        dlcth2 = -1.0e30;
        dlsth2 = 0.0;
    } else {
        let cth2 = ((1.0f32 + costh) * 0.5f32).sqrt() as f64;
        dlcth2 = cth2.ln();
        let sth2 = ((1.0f32 - costh) * 0.5f32).sqrt() as f64;
        dlsth2 = sth2.ln();
    }
    let l1 = lpm + lmn;
    let mut sums: f64 = 0.0;
    let mut nsign: i32 = if (kmin - 1) % 2 == 0 { 1 } else { -1 };
    for k1 in kmin..=kmax {
        let k = k1 - 1;
        let k2 = 2 * k;
        let lpmmk = lpm - k + 1;
        let lmnmk = lmn - k + 1;
        let kmmmn = k - mmn + 1;
        let fac = faco + dlcth2 * ((l1 - k2) as f64) + dlsth2 * ((k2 - mmn) as f64) - lg.g(lpmmk) - lg.g(lmnmk) - lg.g(kmmmn) - lg.g(k1);
        sums += (nsign as f64) * fac.exp();
        nsign = -nsign;
    }
    sums as f32
}

/// `function maxdlm(l,m,n)`: largest ln-gamma index needed by `dlm`.
pub fn maxdlm(l: i32, m: i32, n: i32) -> i32 {
    let l2 = 2 * l;
    let lpm = l + m;
    if lpm > l2 || lpm < 0 {
        return 0;
    }
    let lmm = l - m;
    if lmm > l2 || lmm < 0 {
        return 0;
    }
    let lpn = l + n;
    if lpn > l2 || lpn < 0 {
        return 0;
    }
    let lmn = l - n;
    if lmn > l2 || lmn < 0 {
        return 0;
    }
    let mmn = m - n;
    let mut r = (lpm + 1).max(lmm + 1).max(lpn + 1).max(lmn + 1);
    let kmax = lpm.min(lmn) + 1;
    let kmin = 0.max(mmn) + 1;
    for k1 in kmin..=kmax {
        let k = k1 - 1;
        let lpmmk = lpm - k + 1;
        let lmnmk = lmn - k + 1;
        let kmmmn = k - mmn + 1;
        r = r.max(lpmmk).max(lmnmk).max(kmmmn).max(k1);
    }
    r
}

/// `function cgcoef(l1,l2,l,m1,m2)`: Clebsch-Gordan coefficient.
pub fn cgcoef(lg: &LnGamma, l1: i32, l2: i32, l: i32, m1: i32, m2: i32) -> f32 {
    let m = m1 + m2;
    let ia = l1 + l2 - l + 1;
    if ia <= 0 {
        return 0.0;
    }
    let ib = l + l1 - l2 + 1;
    if ib <= 0 {
        return 0.0;
    }
    let ic = l + l2 - l1 + 1;
    if ic <= 0 {
        return 0.0;
    }
    if m.abs() > l {
        return 0.0;
    }
    let aid = (2 * l + 1) as f32;
    let ie = l1 + l2 + l + 2;
    let mut e1: f32 = (lg.g(ia) + lg.g(ib) + lg.g(ic) - lg.g(ie)) as f32;
    let ja = l1 + m1 + 1;
    let jb = l1 - m1 + 1;
    let jc = l2 + m2 + 1;
    let jd = l2 - m2 + 1;
    let je = l + m + 1;
    let jf = l - m + 1;
    e1 = (e1 as f64 + lg.g(ja) + lg.g(jb) + lg.g(jc) + lg.g(jd) + lg.g(je) + lg.g(jf)) as f32;
    e1 = 0.5f32 * (e1 + aid.ln());
    let kt1 = l - l2 + m1;
    let kt2 = l - l1 - m2;
    let kmin = 0.max(-kt1).max(-kt2) + 1;
    let kt3 = l1 + l2 - l + 1;
    let kt4 = l1 - m1 + 1;
    let kt5 = l2 + m2 + 1;
    let kmax = kt3.min(kt4).min(kt5);
    let mut phase: i32 = if (kmin - 1) % 2 == 0 { 1 } else { -1 };
    let mut sums: f32 = 0.0;
    for k1 in kmin..=kmax {
        let k = k1 - 1;
        let e2: f32 = (lg.g(k1) + lg.g(kt1 + k1) + lg.g(kt2 + k1) + lg.g(kt3 - k) + lg.g(kt4 - k) + lg.g(kt5 - k)) as f32;
        sums += (phase as f32) * (e1 - e2).exp();
        phase = -phase;
    }
    sums
}

/// `function maxcgcoef(l1,l2,l,m1,m2)`.
pub fn maxcgcoef(l1: i32, l2: i32, l: i32, m1: i32, m2: i32) -> i32 {
    let m = m1 + m2;
    let ia = l1 + l2 - l + 1;
    if ia <= 0 {
        return 0;
    }
    let ib = l + l1 - l2 + 1;
    if ib <= 0 {
        return 0;
    }
    let ic = l + l2 - l1 + 1;
    if ic <= 0 {
        return 0;
    }
    if m.abs() > l {
        return 0;
    }
    let ie = l1 + l2 + l + 2;
    let mut r = 0.max(ia).max(ib).max(ic).max(ie);
    let jb = l1 - m1 + 1;
    let jc = l2 + m2 + 1;
    let jd = l2 - m2 + 1;
    let je = l + m + 1;
    let jf = l - m + 1;
    r = r.max(jb).max(jc).max(jd).max(je).max(jf);
    let kt1 = l - l2 + m1;
    let kt2 = l - l1 - m2;
    let kmin = 0.max(-kt1).max(-kt2) + 1;
    let kt3 = l1 + l2 - l + 1;
    let kt4 = l1 - m1 + 1;
    let kt5 = l2 + m2 + 1;
    let kmax = kt3.min(kt4).min(kt5);
    for k1 in kmin..=kmax {
        let k = k1 - 1;
        r = r.max(k1).max(kt1 + k1).max(kt2 + k1).max(kt3 - k).max(kt4 - k).max(kt5 - k);
    }
    r
}

/// `subroutine rotateeuler(xin,yin,zin,phi,theta,chi,xout,yout,zout)`.
pub fn rotateeuler(xin: f32, yin: f32, zin: f32, phi: f32, theta: f32, chi: f32) -> (f32, f32, f32) {
    let calpha = phi.cos();
    let salpha = phi.sin();
    let cbeta = theta.cos();
    let sbeta = theta.sin();
    let cgamma = chi.cos();
    let sgamma = chi.sin();
    let ax = calpha * cbeta * cgamma - salpha * sgamma;
    let bx = -(calpha * cbeta * sgamma) - salpha * cgamma;
    let cx = calpha * sbeta;
    let ay = salpha * cbeta * cgamma + calpha * sgamma;
    let by = -(salpha * cbeta * sgamma) + calpha * cgamma;
    let cy = salpha * sbeta;
    let az = -(sbeta * cgamma);
    let bz = sbeta * sgamma;
    let cz = cbeta;
    (ax * xin + bx * yin + cx * zin, ay * xin + by * yin + cy * zin, az * xin + bz * yin + cz * zin)
}

// ---------------------------------------------------------------------------
// libgfortran cosd / sind (REAL(4))
// ---------------------------------------------------------------------------

const PIO180H: f32 = 1.74560547e-02;
const PIO180L: f32 = -2.76216747e-06;
const COSD30: f32 = 8.66025388e-01;

#[inline(always)]
fn d2r(x: f32) -> f32 {
    x.mul_add(PIO180H, x * PIO180L)
}

/// `_gfortran_sind_r4`.
pub fn gf_sind(x: f32) -> f32 {
    if !x.is_finite() {
        return x - x;
    }
    let mut s: f32 = 1.0f32.copysign(x);
    let ax = x.abs();
    if ax < 0.03125 {
        return d2r(x);
    }
    let mut x = ax % 360.0;
    let n = x as i32;
    if (n as f32) - x == 0.0 && n % 30 == 0 {
        if n % 180 == 0 {
            x = 0.0;
            if n == 180 {
                s = -s;
            }
        } else if n % 90 == 0 {
            x = if n == 90 { 1.0 } else { -1.0 };
        } else if n % 60 == 0 {
            x = COSD30;
            if n >= 180 {
                x = -x;
            }
        } else {
            x = if n < 180 { 0.5 } else { -0.5 };
        }
    } else {
        let mut fold_cos = false;
        if x <= 180.0 {
            if x <= 90.0 {
                if x > 45.0 {
                    x = 90.0 - x;
                    fold_cos = true;
                }
            } else if x <= 135.0 {
                x -= 90.0;
                fold_cos = true;
            } else {
                x = 180.0 - x;
            }
        } else if x <= 270.0 {
            if x <= 225.0 {
                x -= 180.0;
            } else {
                x = 270.0 - x;
                fold_cos = true;
            }
            s = -s;
        } else {
            if x <= 315.0 {
                x -= 270.0;
                fold_cos = true;
            } else {
                x = 360.0 - x;
            }
            s = -s;
        }
        x = d2r(x);
        x = if fold_cos { x.cos() } else { x.sin() };
    }
    x * s
}

/// `_gfortran_cosd_r4`.
pub fn gf_cosd(x: f32) -> f32 {
    if !x.is_finite() {
        return x - x;
    }
    let ax = x.abs();
    if ax <= 0.0078125 {
        // 1 - tiny (= 1 in single precision)
        return 1.0;
    }
    let mut x = ax % 360.0;
    let n = x as i32;
    if (n as f32) - x == 0.0 && n % 30 == 0 {
        if n % 180 == 0 {
            x = if n == 180 { -1.0 } else { 1.0 };
        } else if n % 90 == 0 {
            x = 0.0f32.copysign(0.0);
        } else if n % 60 == 0 {
            x = 0.5;
            if n != 60 && n != 300 {
                x = -x;
            }
        } else {
            x = COSD30;
            if n != 30 && n != 330 {
                x = -x;
            }
        }
    } else {
        let mut neg = false;
        let mut fold_sin = false;
        if x <= 180.0 {
            if x <= 90.0 {
                if x > 45.0 {
                    x = 90.0 - x;
                    fold_sin = true;
                }
            } else {
                if x <= 135.0 {
                    x -= 90.0;
                    fold_sin = true;
                } else {
                    x = 180.0 - x;
                }
                neg = true;
            }
        } else if x <= 270.0 {
            if x <= 225.0 {
                x -= 180.0;
            } else {
                x = 270.0 - x;
                fold_sin = true;
            }
            neg = true;
        } else if x <= 315.0 {
            x -= 270.0;
            fold_sin = true;
        } else {
            x = 360.0 - x;
        }
        x = d2r(x);
        x = if fold_sin { x.sin() } else { x.cos() };
        if neg {
            x = -x;
        }
    }
    x
}

// ---------------------------------------------------------------------------
// sharm_mod + common blocks of sharm.inc / sharm_sdf.inc
// ---------------------------------------------------------------------------

/// The data of module `sharm_mod` plus the variables of the common blocks
/// `sharm.inc`/`sharm_sdf.inc` used by the routines below.
#[derive(Clone, Debug, Default)]
pub struct Sharm {
    pub lg: LnGamma,
    // specification of l1,l2,l
    pub nll1v: i32,
    pub nll2v: i32,
    pub nllv: i32,
    pub n1step: i32,
    pub n2step: i32,
    pub ll1v: F1<i32>,
    pub ll2v: F1<i32>,
    pub llv: F1<i32>,
    pub nlmax: i32,
    pub ndist: i32,
    pub ldistinct: i32,
    pub ll1: F1<i32>,
    pub ll2: F1<i32>,
    pub ll: F1<i32>,
    pub n1: F1<i32>,
    pub n2: F1<i32>,
    pub alsec: F1<f32>,
    pub mdlmn: i32,
    pub mz: i32,
    pub mdist: i32,
    pub nm: F1<i32>,
    pub il1rf: F1<i32>,
    pub il2rf: F1<i32>,
    pub ilrf: F1<i32>,
    pub phim: F1<f32>,
    pub cthetm: F1<f32>,
    pub chim: F1<f32>,
    pub xcent: F1<f32>,
    pub ycent: F1<f32>,
    pub zcent: F1<f32>,
    pub philab: F1<f32>,
    pub dlmnc: F2<f32>,
    pub dlmnmc: F3<f32>,
    pub chiphi: F3<f32>,
    pub mcg: i32,
    pub termcg: F1<f32>,
    pub phimt1: F1<f32>,
    pub cthetmt1: F1<f32>,
    pub chimt1: F1<f32>,
    pub phimt2: F1<f32>,
    pub cthetmt2: F1<f32>,
    pub chimt2: F1<f32>,
    pub dlmnmct: F3<f32>,
    pub chiphit: F3<f32>,
    pub z: F1<f32>,
    pub gofz: F2<f32>,
    pub gofzs: F2<f32>,
    pub gzt: F2<f32>,
    pub gzt2: F2<f32>,
    pub titledat: Vec<String>,
    pub performtest: bool,
    pub outputaveragemol: bool,
    // common blocks
    pub pi: f32,
    pub piconv: f32,
    pub nz: i32,
    pub zstep: f32,
    pub zmaxd: f32,
    pub vol: f32,
    pub phim1: f32,
    pub thetam1: f32,
    pub chim1: f32,
    pub phim2: f32,
    pub thetam2: f32,
    pub chim2: f32,
    pub thetallab: f32,
    pub phillab: f32,
}

/// Count the n2 range of `get_l_and_n_values` for given il2, in1.
fn n2_range(n2step: i32, il2: i32, in1: i32) -> (i32, i32, i32) {
    if n2step > 0 {
        let nts = il2 / n2step;
        let n2min = if in1 > 0 { -nts * n2step } else { 0 };
        (n2min, nts * n2step, n2step)
    } else {
        (0, 0, 1)
    }
}

/// Fortran DO loop `do i=a,b,s` (s > 0) as an iterator.
fn do_range(a: i32, b: i32, s: i32) -> impl Iterator<Item = i32> {
    let n = if b >= a { (b - a) / s + 1 } else { 0 };
    (0..n).map(move |k| a + k * s)
}

impl Sharm {
    /// `subroutine get_l_and_n_values()` (also calls `form_cgcoeffs`).
    pub fn get_l_and_n_values(&mut self, log: &dyn Fn(&str)) {
        let mut terms: Vec<(i32, i32, i32, i32, i32)> = Vec::new();
        for ill1 in 1..=self.nll1v {
            let il1 = self.ll1v[ill1];
            for ill2 in 1..=self.nll2v {
                let il2 = self.ll2v[ill2];
                if (self.ldistinct == 0 && il2 <= il1) || self.ldistinct == 1 {
                    let jlmax = (il1 + il2).abs();
                    let jlmin = (il1 - il2).abs();
                    for ill in 1..=self.nllv {
                        let il = self.llv[ill];
                        if il >= jlmin && il <= jlmax {
                            let (n1max, n1s) = if self.n1step > 0 { (il1 + 1, self.n1step) } else { (1, 1) };
                            for jn1 in do_range(1, n1max, n1s) {
                                let in1 = jn1 - 1;
                                let (n2min, n2max, n2s) = n2_range(self.n2step, il2, in1);
                                for in2 in do_range(n2min, n2max, n2s) {
                                    terms.push((il1, il2, il, in1, in2));
                                }
                            }
                        }
                    }
                }
            }
        }
        let nterm = terms.len();
        self.ndist = nterm as i32;
        self.nlmax = nterm as i32;
        self.ll1 = F1::new(nterm);
        self.ll2 = F1::new(nterm);
        self.ll = F1::new(nterm);
        self.n1 = F1::new(nterm);
        self.n2 = F1::new(nterm);
        self.il1rf = F1::new(nterm);
        self.il2rf = F1::new(nterm);
        self.ilrf = F1::new(nterm);
        self.alsec = F1::new(nterm);
        for (k, t) in terms.iter().enumerate() {
            let i = k + 1;
            self.ll1[i] = t.0;
            self.ll2[i] = t.1;
            self.ll[i] = t.2;
            self.n1[i] = t.3;
            self.n2[i] = t.4;
        }
        for i in 1..=self.ndist {
            log(&fmt("(6(1x,i5))", &[FV::I(i as i64), FV::I(self.ll1[i] as i64), FV::I(self.ll2[i] as i64), FV::I(self.ll[i] as i64), FV::I(self.n1[i] as i64), FV::I(self.n2[i] as i64)]));
        }
        self.form_cgcoeffs(log);
    }

    /// `subroutine form_cgcoeffs()`.
    pub fn form_cgcoeffs(&mut self, log: &dyn Fn(&str)) {
        for id in 1..=self.ndist {
            self.alsec[id] = 0.0;
        }
        let mut ic = 0;
        let mut mgammat = 0;
        for id in 1..=self.ndist {
            let il1 = self.ll1[id];
            let in1 = self.n1[id];
            let il2 = self.ll2[id];
            let in2 = self.n2[id];
            let il = self.ll[id];
            self.il1rf[id] = 1 + getrefint(il1, in1, &self.ll1v, self.nll1v, self.n1step) + il1;
            self.il2rf[id] = 1 + getrefint(il2, in2, &self.ll2v, self.nll2v, self.n2step) + il2;
            self.ilrf[id] = 1 + getrefintnon(il, &self.llv, self.nllv) + il;
            for m1 in -il1..=il1 {
                mgammat = mgammat.max(maxdlm(il1, m1, in1));
                mgammat = mgammat.max(maxdlm(il1, m1, -in1));
                let m2max = il2.min(il - m1);
                let m2min = (-il2).max(-il - m1);
                for m2 in m2min..=m2max {
                    ic += 1;
                    let m = m1 + m2;
                    mgammat = mgammat.max(maxcgcoef(il1, il2, il, m1, m2));
                    mgammat = mgammat.max(maxdlm(il2, m2, in2));
                    mgammat = mgammat.max(maxdlm(il2, m2, -in2));
                    mgammat = mgammat.max(maxdlm(il, m, 0));
                }
            }
        }
        log(&format!(" form_cgcoeffs> Creating {:>11}  log gamma values", mgammat));
        self.lg.fctrln(mgammat);
        self.mcg = ic;
        self.termcg = F1::new(ic.max(0) as usize);
        let mut ic = 0;
        for id in 1..=self.ndist {
            let il1 = self.ll1[id];
            let il2 = self.ll2[id];
            let il = self.ll[id];
            for m1 in -il1..=il1 {
                let m2max = il2.min(il - m1);
                let m2min = (-il2).max(-il - m1);
                for m2 in m2min..=m2max {
                    ic += 1;
                    self.termcg[ic] = cgcoef(&self.lg, il1, il2, il, m1, m2);
                    self.alsec[id] += self.termcg[ic] * self.termcg[ic];
                }
            }
        }
    }

    /// `subroutine get_max_size_dlmn()`.
    pub fn get_max_size_dlmn(&mut self) {
        let mut ilm = 0;
        for i in 1..=self.nll1v {
            let il = self.ll1v[i];
            let (ntmax, nts) = if self.n1step > 0 { (il, self.n1step) } else { (0, 1) };
            for _ in do_range(0, ntmax, nts) {
                ilm += 2 * il + 1;
            }
        }
        self.mdlmn = ilm;
        ilm = 0;
        for i in 1..=self.nll2v {
            let il = self.ll2v[i];
            let (ntmax, nts) = if self.n2step > 0 { (il, self.n2step) } else { (0, 1) };
            for _ in do_range(0, ntmax, nts) {
                ilm += 2 * il + 1;
            }
        }
        self.mdlmn = self.mdlmn.max(ilm);
        ilm = 0;
        for ilt in 1..=self.nllv {
            ilm += 2 * self.llv[ilt] + 1;
        }
        self.mdlmn = self.mdlmn.max(ilm);
    }

    /// Allocation of the per-molecule arrays done in `sharm`/`sdf`.
    pub fn alloc_molecules(&mut self, nmol: i32) {
        let n = nmol.max(0) as usize;
        let md = self.mdlmn.max(0) as usize;
        self.nm = F1::new(n);
        self.phim = F1::new(n);
        self.cthetm = F1::new(n);
        self.chim = F1::new(n);
        self.xcent = F1::new(n);
        self.ycent = F1::new(n);
        self.zcent = F1::new(n);
        self.dlmnmc = F3::new(2, n, md);
        self.chiphi = F3::new(2, n, md);
        self.dlmnc = F2::new(2, md);
        self.philab = F1::new(md);
        if self.performtest {
            self.phimt1 = F1::new(n);
            self.cthetmt1 = F1::new(n);
            self.chimt1 = F1::new(n);
            self.phimt2 = F1::new(n);
            self.cthetmt2 = F1::new(n);
            self.chimt2 = F1::new(n);
            self.dlmnmct = F3::new(2, n, md);
            self.chiphit = F3::new(2, n, md);
        }
    }

    /// `subroutine form_dlmn(imol1,ierr)`.
    pub fn form_dlmn(&mut self, imol1: i32, log: &dyn Fn(&str)) -> i32 {
        let (lltt, ntstep): (Vec<i32>, i32) = if self.nm[imol1] == 1 {
            ((1..=self.nll1v).map(|i| self.ll1v[i]).collect(), self.n1step)
        } else {
            ((1..=self.nll2v).map(|i| self.ll2v[i]).collect(), self.n2step)
        };
        let mut ilmnlm = 0;
        for &il in lltt.iter() {
            let (ntmax, nts) = if ntstep > 0 { (il, ntstep) } else { (0, 1) };
            for inn in do_range(0, ntmax, nts) {
                for m in -il..=il {
                    ilmnlm += 1;
                    if ilmnlm > self.mdlmn {
                        let ierr = ilmnlm - self.mdlmn;
                        log(&format!(" form_dlmn> Array reference {:>11}  larger than allowed {:>11}", ilmnlm, self.mdlmn));
                        return ierr;
                    }
                    let second = inn != 0 && (self.nm[imol1] == 2 || (self.nm[imol1] == 1 && self.ldistinct == 0));
                    if self.performtest {
                        self.chiphi[(1, imol1, ilmnlm)] = (m as f32) * self.phimt1[imol1] + (inn as f32) * self.chimt1[imol1];
                        self.dlmnmc[(1, imol1, ilmnlm)] = dlm(&self.lg, il, m, inn, self.cthetmt1[imol1]);
                        self.chiphit[(1, imol1, ilmnlm)] = (m as f32) * self.phimt2[imol1] + (inn as f32) * self.chimt2[imol1];
                        self.dlmnmct[(1, imol1, ilmnlm)] = dlm(&self.lg, il, m, inn, self.cthetmt2[imol1]);
                        if second {
                            self.chiphi[(2, imol1, ilmnlm)] = (m as f32) * self.phimt1[imol1] - (inn as f32) * self.chimt1[imol1];
                            self.dlmnmc[(2, imol1, ilmnlm)] = dlm(&self.lg, il, m, -inn, self.cthetmt1[imol1]);
                            self.chiphit[(2, imol1, ilmnlm)] = (m as f32) * self.phimt2[imol1] - (inn as f32) * self.chimt2[imol1];
                            self.dlmnmct[(2, imol1, ilmnlm)] = dlm(&self.lg, il, m, -inn, self.cthetmt2[imol1]);
                        }
                    } else {
                        self.chiphi[(1, imol1, ilmnlm)] = (m as f32) * self.phim[imol1] + (inn as f32) * self.chim[imol1];
                        self.dlmnmc[(1, imol1, ilmnlm)] = dlm(&self.lg, il, m, inn, self.cthetm[imol1]);
                        if second {
                            self.chiphi[(2, imol1, ilmnlm)] = (m as f32) * self.phim[imol1] - (inn as f32) * self.chim[imol1];
                            self.dlmnmc[(2, imol1, ilmnlm)] = dlm(&self.lg, il, m, -inn, self.cthetm[imol1]);
                        }
                    }
                }
            }
        }
        0
    }

    /// `subroutine form_sum_prod_dlmn_cgcoeff(nmolsharm,ierr)`.
    pub fn form_sum_prod_dlmn_cgcoeff(&mut self, nmolsharm: i32, log: &dyn Fn(&str)) {
        let gfac: f32 = 4.0f32 * self.pi / 3.0f32;
        let ld = self.ldistinct;
        let (mut lsum1, mut lsum2);
        if ld == 0 {
            if self.nm[nmolsharm] == 1 || self.nm[nmolsharm] == 2 {
                lsum1 = 1;
                lsum2 = 1;
            } else {
                lsum1 = 0;
                lsum2 = 0;
            }
        } else {
            lsum1 = 0;
            lsum2 = 0;
            if self.nm[nmolsharm] == 1 {
                lsum1 = 1;
            }
            if self.nm[nmolsharm] == 2 {
                lsum2 = 1;
            }
        }
        let alab = [1.0f32, 0.0, 0.0];
        let blab = [0.0f32, 1.0, 0.0];
        let clab = [0.0f32, 0.0, 1.0];
        let pt = self.performtest;
        for imol1 in 1..nmolsharm {
            if self.nm[imol1] == 0 {
                continue;
            }
            if ld == 0 {
                if self.nm[imol1] == 1 || self.nm[imol1] == 2 {
                    lsum1 += 1;
                    lsum2 += 1;
                }
            } else {
                if self.nm[imol1] == 1 {
                    lsum1 += 1;
                }
                if self.nm[imol1] == 2 {
                    lsum2 += 1;
                }
            }
            let xc1 = self.xcent[imol1];
            let yc1 = self.ycent[imol1];
            let zc1 = self.zcent[imol1];
            for imol2 in imol1 + 1..=nmolsharm {
                if imol1 == imol2 || self.nm[imol2] == 0 {
                    continue;
                }
                let mut iflag1 = 0;
                let mut iflag2 = 0;
                let nm1 = self.nm[imol1];
                let nm2 = self.nm[imol2];
                if ld == 0 && nm1 == 1 && nm2 == 1 {
                    iflag1 = 1;
                    iflag2 = 1;
                }
                if ld == 1 && nm1 == 1 && nm2 == 2 {
                    iflag1 = 1;
                    iflag2 = 0;
                }
                if ld == 1 && nm1 == 2 && nm2 == 1 {
                    iflag1 = 0;
                    iflag2 = 1;
                }
                if iflag1 + iflag2 == 0 {
                    continue;
                }
                let mut xdif = self.xcent[imol2] - xc1;
                let mut ydif = self.ycent[imol2] - yc1;
                let mut zdif = self.zcent[imol2] - zc1;
                crate::neighbour::minimagec(&mut xdif, &mut ydif, &mut zdif);
                let mut pdif = [xdif, ydif, zdif];
                let opdif = scapro(&pdif, &pdif).sqrt();
                let d12 = opdif;
                for v in pdif.iter_mut() {
                    *v /= opdif;
                }
                let iref = (d12 / self.zstep) as i32 + 1;
                if iref > self.nz {
                    continue;
                }
                let mut cz = scapro(&clab, &pdif);
                let mut perp1 = vecpro(&clab, &pdif);
                let operp1 = scapro(&perp1, &perp1).sqrt();
                let mut phil: f32;
                if operp1 > 0.0 {
                    for v in perp1.iter_mut() {
                        *v /= operp1;
                    }
                    let sphil = -scapro(&alab, &perp1);
                    let cphil = scapro(&blab, &perp1);
                    phil = sphil.atan2(cphil);
                } else {
                    phil = 0.0;
                }
                if pt {
                    phil = self.phillab;
                    cz = self.thetallab.cos();
                }
                let mut ilm = 0;
                for ilt in 1..=self.nllv {
                    let il = self.llv[ilt];
                    if il <= self.ll[self.ndist] {
                        for m in -il..=il {
                            ilm += 1;
                            self.philab[ilm] = (m as f32) * phil;
                            self.dlmnc[(1, ilm)] = dlm(&self.lg, il, m, 0, cz);
                            let sgn: f32 = if m % 2 == 0 { 1.0 } else { -1.0 };
                            self.dlmnc[(2, ilm)] = sgn * dlm(&self.lg, il, m, 0, -cz);
                        }
                    }
                }
                // (molecule at origin, other molecule, which dlmnc row)
                let mut passes: Vec<(i32, i32, usize)> = Vec::new();
                if iflag1 == 1 {
                    passes.push((imol1, imol2, 1));
                }
                if iflag2 == 1 && !pt {
                    passes.push((imol2, imol1, 2));
                }
                for (ma, mb, row) in passes {
                    let mut ic = 0;
                    for id in 1..=self.ndist {
                        let il1 = self.ll1[id];
                        let il2 = self.ll2[id];
                        let in2 = self.n2[id];
                        let il = self.ll[id];
                        let il1ref = self.il1rf[id];
                        let il2ref = self.il2rf[id];
                        let ilref = self.ilrf[id];
                        let sec = if in2 < 0 { 2 } else { 1 };
                        let mut gc = self.gofz[(iref, id)];
                        let mut gs = self.gofzs[(iref, id)];
                        for m1 in -il1..=il1 {
                            let m2max = il2.min(il - m1);
                            let m2min = (-il2).max(-il - m1);
                            for m2 in m2min..=m2max {
                                let m = m1 + m2;
                                ic += 1;
                                let ilmn1 = il1ref + m1;
                                let ilmn2 = il2ref + m2;
                                let ilmn = ilref + m;
                                let (phisum, prod) = if pt {
                                    (
                                        self.philab[ilmn] - self.chiphi[(1, ma, ilmn1)] - self.chiphit[(sec, mb, ilmn2)],
                                        self.dlmnmc[(1, ma, ilmn1)] * self.dlmnmct[(sec, mb, ilmn2)],
                                    )
                                } else {
                                    (
                                        self.philab[ilmn] - self.chiphi[(1, ma, ilmn1)] - self.chiphi[(sec, mb, ilmn2)],
                                        self.dlmnmc[(1, ma, ilmn1)] * self.dlmnmc[(sec, mb, ilmn2)],
                                    )
                                };
                                let cosphi = phisum.cos();
                                let sinphi = phisum.sin();
                                let factor = prod * self.dlmnc[(row, ilmn)] * self.termcg[ic];
                                gc += factor * cosphi;
                                gs += factor * sinphi;
                            }
                        }
                        self.gofz[(iref, id)] = gc;
                        self.gofzs[(iref, id)] = gs;
                    }
                }
            }
        }
        for id in 1..=self.ndist {
            log(&fmt(
                "(1x,i6,5(1x,i3),3(1x,i5),1x,e12.5)",
                &[
                    FV::I(id as i64),
                    FV::I(self.ll1[id] as i64),
                    FV::I(self.ll2[id] as i64),
                    FV::I(self.ll[id] as i64),
                    FV::I(self.n1[id] as i64),
                    FV::I(self.n2[id] as i64),
                    FV::I(self.il1rf[id] as i64),
                    FV::I(self.il2rf[id] as i64),
                    FV::I(self.ilrf[id] as i64),
                    FV::R(self.alsec[id]),
                ],
            ));
        }
        let rho: f32 = (lsum1 as f32) * (lsum2 as f32) / self.vol;
        log(&fmt("(1x,'rho = ',e13.6)", &[FV::R(rho)]));
        let gfac1 = gfac * rho;
        for id in 1..=self.ndist {
            let mut z0: f32 = 0.0;
            let mut z03: f32 = 0.0;
            let fac = ((2 * self.ll1[id] + 1) * (2 * self.ll2[id] + 1)) as f32;
            for iz in 1..=self.nz {
                let zm = self.z[iz];
                let zdiff = zm - z0;
                let z1 = zm + zdiff;
                let z13 = z1 * z1 * z1;
                let rat = gfac1 * (z13 - z03);
                self.gofz[(iz, id)] = fac * self.gofz[(iz, id)] / rat;
                self.gofzs[(iz, id)] = fac * self.gofzs[(iz, id)] / rat;
                if id == 1 {
                    self.gofz[(iz, id)] -= 1.0;
                }
                z0 = z1;
                z03 = z13;
            }
        }
    }

    /// Titles `(l1,l2,l;n1,n2)` of the coefficient columns (written in
    /// `run_sharm`/`run_sdf` with a run-time format
    /// `(1h(,i1,1h,,i1,1h,,i1,1h;,i1,1h,,i1,1h))`; see the module notes).
    pub fn make_titles(&mut self) {
        self.titledat = (1..=self.ndist)
            .map(|id| format!("({},{},{};{},{})", self.ll1[id], self.ll2[id], self.ll[id], self.n1[id], self.n2[id]))
            .collect();
    }

    /// `subroutine outgofzs_sharm(fnameinp,ext,nz,ndist,nsum,z,gofz,gofzs,gzt,gzt2,titler,titledat,...)`.
    pub fn outgofzs_sharm(&mut self, fnameinp: &str, ext: &str, nsum: &mut i32, titler: &str, log: &dyn Fn(&str)) -> Result<()> {
        let ilen = indexlastdot(fnameinp) as usize;
        let fname = format!("{}{}", &fnameinp[..ilen], ext);
        let nz = self.nz;
        let ndist = self.ndist;
        let mut ierr = 0;
        let mut lptin = nz;
        if *nsum > 0 {
            // r_int also overwrites z with the values read from the file
            match r_int(&fname, ndist, &mut self.z, &mut self.gzt, &mut self.gzt2, log)? {
                Some(l) => lptin = l,
                None => ierr = 1,
            }
        }
        if ierr != 0 {
            *nsum = -1;
        }
        if *nsum > 0 {
            for id in 1..=ndist {
                for j in 1..=nz {
                    if j > lptin {
                        self.gzt[(j, id)] = 0.0;
                        self.gzt2[(j, id)] = 0.0;
                    }
                }
            }
        } else {
            for id in 1..=ndist {
                for j in 1..=nz {
                    self.gzt[(j, id)] = 0.0;
                    self.gzt2[(j, id)] = 0.0;
                }
            }
        }
        let (fc1, fc2): (f64, f64);
        if *nsum >= 0 {
            *nsum += 1;
            fc1 = (1.0f32 / (*nsum as f32)) as f64;
            fc2 = 1.0 - fc1;
        } else {
            fc1 = 1.0;
            fc2 = 0.0;
        }
        for id in 1..=ndist {
            for i in 1..=nz {
                let a2 = fc2 * self.gzt[(i, id)] as f64;
                let a1 = fc1 * self.gofz[(i, id)] as f64;
                let b2 = fc2 * self.gzt2[(i, id)] as f64;
                let b1 = fc1 * self.gofzs[(i, id)] as f64;
                self.gzt[(i, id)] = (a1 + a2) as f32;
                self.gzt2[(i, id)] = (b1 + b2) as f32;
            }
        }
        w_int(nz, &fname, ndist, &self.z, &self.gzt, &self.gzt2, titler, &self.titledat, log)
    }
}

/// `function getrefint(lreq,nreq,llv,nllv,nstep)`.
pub fn getrefint(lreq: i32, nreq: i32, llv: &F1<i32>, nllv: i32, nstep: i32) -> i32 {
    let mut sum = 0;
    let mut id = 1;
    while id < nllv && llv[id] < lreq {
        let l = llv[id];
        let nts = if nstep > 0 { l / nstep + 1 } else { 1 };
        sum += nts * (2 * l + 1);
        id += 1;
    }
    let l = llv[id];
    let ntmax = nreq.abs();
    if nstep > 0 {
        let nts = ntmax / nstep;
        sum += nts * (2 * l + 1);
    }
    sum
}

/// `function getrefintnon(lreq,llv,nllv)`.
pub fn getrefintnon(lreq: i32, llv: &F1<i32>, nllv: i32) -> i32 {
    let mut sum = 0;
    let mut id = 1;
    while id < nllv && llv[id] < lreq {
        sum += 2 * llv[id] + 1;
        id += 1;
    }
    sum
}

/// `subroutine form_eulerangles(imol1,xt,yt,zt,naxis1,naxis2,ncz,ncy,lztyp,lytyp,alab,blab,clab,phim,cthetm,chim)`.
/// `xt,yt,zt` are 1-based coordinate arrays; returns (phim, cthetm, chim).
#[allow(clippy::too_many_arguments)]
pub fn form_eulerangles(xt: &F1<f32>, yt: &F1<f32>, zt: &F1<f32>, naxis1: &str, naxis2: &str, lztyp: &[i32], lytyp: &[i32]) -> (f32, f32, f32) {
    let alab = [1.0f32, 0.0, 0.0];
    let blab = [0.0f32, 1.0, 0.0];
    let clab = [0.0f32, 0.0, 1.0];
    let mut ifindz = 0;
    let mut ifindy = 0;
    let (mut czx, mut czy, mut czz, mut cyx, mut cyy, mut cyz) = (0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32);
    for &ia1 in lztyp {
        ifindz += 1;
        czx += xt[ia1];
        czy += yt[ia1];
        czz += zt[ia1];
    }
    for &ia1 in lytyp {
        ifindy += 1;
        cyx += xt[ia1];
        cyy += yt[ia1];
        cyz += zt[ia1];
    }
    let mut ct: [f32; 3];
    let mut pt: [f32; 3];
    if ifindz > 0 {
        let f = ifindz as f32;
        ct = [czx / f, czy / f, czz / f];
    } else {
        ct = [0.0, 0.0, 1.0];
    }
    if ifindy > 0 {
        let f = ifindy as f32;
        pt = [cyx / f, cyy / f, cyz / f];
    } else {
        pt = [0.0, 0.0, 0.0];
    }
    let mut oct = scapro(&ct, &ct).sqrt();
    if oct > 0.0 {
        for v in ct.iter_mut() {
            *v /= oct;
        }
        oct = 1.0;
    } else {
        oct = 0.0;
    }
    let mut at = [0.0f32; 3];
    let mut bt = [0.0f32; 3];
    let (oat, obt);
    let opt = scapro(&pt, &pt).sqrt();
    if opt > 0.0 {
        for v in pt.iter_mut() {
            *v /= opt;
        }
        at = vecpro(&pt, &ct);
        let o = scapro(&at, &at).sqrt();
        for v in at.iter_mut() {
            *v /= o;
        }
        oat = 1.0f32;
        bt = vecpro(&ct, &at);
        let o = scapro(&bt, &bt).sqrt();
        for v in bt.iter_mut() {
            *v /= o;
        }
        obt = 1.0f32;
    } else {
        oat = 0.0f32;
        obt = 0.0f32;
    }
    let neg = |v: [f32; 3]| [-v[0], -v[1], -v[2]];
    let mut a1 = [0.0f32; 3];
    let mut b1 = [0.0f32; 3];
    let mut c1 = [0.0f32; 3];
    let (mut oa1, mut ob1, mut oc1) = (0.0f32, 0.0f32, 0.0f32);
    let n1 = naxis1.trim_end();
    let n2 = naxis2.trim_end();
    let mut set = |a: [f32; 3], b: [f32; 3], c: [f32; 3], oa: f32, ob: f32, oc: f32| {
        a1 = a;
        b1 = b;
        c1 = c;
        oa1 = oa;
        ob1 = ob;
        oc1 = oc;
    };
    match (n1, n2) {
        ("z", "y") => set(at, bt, ct, oat, obt, oct),
        ("-z", "y") => set(at, bt, neg(ct), oat, obt, oct),
        ("y", "z") => set(neg(at), ct, bt, oat, oct, obt),
        ("y", "-z") => set(neg(at), ct, neg(bt), oat, oct, obt),
        ("y", "x") => set(bt, ct, at, obt, oct, oat),
        ("x", "y") => set(ct, bt, neg(at), oct, obt, oat),
        ("z", "x") => set(bt, neg(at), ct, obt, oat, oct),
        ("-z", "x") => set(bt, neg(at), neg(ct), obt, oat, oct),
        ("x", "z") => set(ct, at, bt, oct, oat, obt),
        ("x", "-z") => set(ct, at, neg(bt), oct, oat, obt),
        _ => {}
    }
    let (phim, cthetm, chim);
    if oc1 > 0.0 {
        cthetm = scapro(&clab, &c1);
        let mut perp1 = vecpro(&clab, &c1);
        let operp1 = scapro(&perp1, &perp1).sqrt();
        if operp1 > 0.0 {
            for v in perp1.iter_mut() {
                *v /= operp1;
            }
            let sphim = -scapro(&perp1, &alab);
            let cphim = scapro(&perp1, &blab);
            phim = sphim.atan2(cphim);
        } else {
            phim = 0.0;
        }
        if oa1 > 0.0 && ob1 > 0.0 {
            let schim = scapro(&perp1, &a1);
            let cchim = scapro(&perp1, &b1);
            chim = schim.atan2(cchim);
        } else {
            chim = 0.0;
        }
    } else {
        cthetm = 0.0;
        phim = 0.0;
        chim = 0.0;
    }
    (phim, cthetm, chim)
}

/// `subroutine smoosharm(nqmin,nqmax,nsmoo,diff)` (diff 1-based).
pub fn smoosharm(nqmin: i32, nqmax: i32, nsmoo: i32, diff: &mut F1<f32>) {
    let n1 = nqmin + 1;
    let n2 = nqmax - 1;
    diff[nqmax] = 0.0;
    if (n2 - n1) < 2 || nsmoo == 0 {
        return;
    }
    for _ in 0..nsmoo.max(0) {
        let mut sum1 = 0.25f32 * (diff[nqmin] + diff[n1]);
        let mut sm = diff[n1];
        for i in n1..=n2 {
            let sn = diff[i + 1];
            let sum2 = 0.25f32 * (sm + sn);
            diff[i] = sum1 + sum2;
            sm = sn;
            sum1 = sum2;
        }
    }
}

/// `subroutine sumoverm2(l1,l2,l,n1,m2f,m2l,n2,phil,ctl,phi,cz,chi,phs,sumc,sums)`.
#[allow(clippy::too_many_arguments)]
pub fn sumoverm2(lg: &LnGamma, l1: i32, l2: i32, l: i32, n1: i32, m2f: i32, m2l: i32, n2: i32, phil: f32, ctl: f32, phi: f32, cz: f32, chi: f32, phs: f32, sumc: &mut f32, sums: &mut f32) {
    // phssine=-real(sign(1,n1)): sign(1,0) = 1
    let phssine: f32 = if n1 >= 0 { -1.0 } else { 1.0 };
    let m1 = n1;
    for m2 in m2f..=m2l {
        let m = m1 + m2;
        let cg = cgcoef(lg, l1, l2, l, m1, m2);
        let dl2m2n2 = dlm(lg, l2, m2, n2, cz);
        let dlm0 = dlm(lg, l, m, 0, ctl);
        let phisum = (m2 as f32) * phi + (n2 as f32) * chi - (m as f32) * phil;
        let prod = phs * cg * dl2m2n2 * dlm0;
        *sumc += prod * gf_cosd(phisum);
        *sums += phssine * prod * gf_sind(phisum);
    }
}

/// `subroutine r_int_titles(fname,nfile,titlex,titled,nfdim,ierr)`: read the
/// column titles of `fname(1:indexlastdot+1)//nn` files.  Returns (nfile, titlex, titled).
pub fn r_int_titles(fname: &str, nfdim: i32, log: &dyn Fn(&str)) -> (i32, String, Vec<String>) {
    let ilen = (indexlastdot(fname) + 1) as usize;
    let base: String = fname.chars().take(ilen).collect();
    let mut titlex = String::new();
    let mut titled: Vec<String> = Vec::new();
    let mut ierr = 0;
    let mut nwr = NWRMAX;
    let mut it = 0;
    while ierr == 0 && nwr >= NWRMAX {
        it += 1;
        let f = format!("{base}{:02}", it);
        match TextFile::open(std::path::Path::new(&f)) {
            Ok(mut tf) => {
                let sometext = tf.read_line().unwrap_or_default();
                if sometext.contains('#') {
                    let w = crate::fortran::text::parse_n(&sometext, 20, crate::fortran::text::NSTSIZE);
                    if w.nwords() > 1 {
                        titlex = crate::fortran::text::fstr(w.word(&sometext, 2), 14);
                        let mut iw = 2;
                        let mut ir = 0;
                        let mut titleo = Vec::new();
                        while iw < w.nwords() && ir < NWRMAX {
                            iw += 1;
                            ir += 1;
                            titleo.push(crate::fortran::text::fstr(w.word(&sometext, iw), 29));
                        }
                        nwr = ir;
                        for t in titleo {
                            if (titled.len() as i32) < nfdim {
                                titled.push(t);
                            }
                        }
                    } else {
                        ierr = 1;
                    }
                } else {
                    ierr = 1;
                }
            }
            Err(_) => ierr = 1,
        }
    }
    log(&format!(" r_int_titles> nfile: {:>11}", titled.len()));
    (titled.len() as i32, titlex, titled)
}

/// Parse a title `(l1,l2,l;n1,n2)` as `run_mapgr` does: returns the five integers.
pub fn parse_title(t: &str) -> Result<[i32; 5]> {
    // temptitle(ifirst:ilast) slices with index(...) of the delimiters
    let tt = t.trim_end();
    let mut out = [0i32; 5];
    let delims = [',', ',', ';', ',', ')'];
    let b = tt.as_bytes();
    let lentemptitle = b.len() as i32;
    let idx = |from: i32, c: char| -> i32 {
        // index(temptitle(from:lentemptitle), c)
        if from > lentemptitle {
            return 0;
        }
        match tt[(from - 1) as usize..].find(c) {
            Some(p) => p as i32 + 1,
            None => 0,
        }
    };
    let mut ifirst = idx(1, '(') + 1;
    for k in 0..5 {
        if k > 0 {
            // ifirst = ilast + 2 handled below
        }
        let ilast = ifirst + idx(ifirst, delims[k]) - 2;
        let s = crate::fortran::text::substr(tt, ifirst, ilast);
        out[k] = read_int(&s).map_err(|_| EpsrError::Input(format!("cannot read coefficient indices from title {t:?}")))?;
        ifirst = ilast + 2;
    }
    Ok(out)
}

impl State {
    /// The `mapgrentropy`/`mapgr` branch of `check_symbol` (`process_symbols.f`)
    /// for the symbol `shcoeffs`.  Returns `ierr`.
    pub fn check_symbol_mapgr(&mut self, icount: i32, iref: i32) -> Result<i32> {
        if self.sym.symbol[icount] != "shcoeffs" {
            return Ok(0);
        }
        let text = self.sym.symbol_value[self.sym.ireff(self.sym.nsymbol, self.sym.nfiles)].clone();
        let w = parse(&text, 5);
        let ext = if w.nwords() > 0 { w.word(&text, 1).to_string() } else { String::new() };
        let v = self.sym.symbol_value[iref].clone();
        if v == "search" || v == "<undefined>" {
            return crate::error::unsupported("searching for a coefficients file (shcoeffs search)");
        }
        let w = parse(&v, 5);
        if w.nwords() <= 0 {
            self.log(&format!("\n{} No coefficients file specified! - routine will exit\n", self.sym.prompt));
            self.sym.symbol_value[iref] = self.sym.symbol_save.clone();
            return Ok(0);
        }
        let mut word = w.word(&v, 1).to_string();
        if !word.contains('.') {
            word.push_str(&ext);
        }
        self.sym.symbol_value[iref] = word;
        self.find_plot_par()
    }

    /// `subroutine find_plot_par(nin,nout,ierr)` (`plot_routines.f`).
    fn find_plot_par(&mut self) -> Result<i32> {
        let ifound = self.sym.search_for_symbol("shcoeffs", 0);
        let iref = self.sym.ireff(ifound, self.sym.nfiles);
        let sv = self.sym.symbol_value[iref].clone();
        let mut ilen = indexlastdot(&sv) - 1;
        if ilen <= 0 {
            ilen = crate::fortran::text::len_to_blank(&sv) as i32;
        }
        let stem = crate::fortran::text::substr(&sv, 1, ilen);
        let fnametemp = format!("{}{}.dat", self.sym.directory_name, stem);
        let mut f = match TextFile::open(std::path::Path::new(&fnametemp)) {
            Ok(f) => f,
            Err(_) => return crate::error::unsupported("mapgr: SHARM/SDF .dat file of the coefficients file not found (get_nterms)"),
        };
        let mut centrec = String::new();
        let mut centres = String::new();
        while centrec.is_empty() || centres.is_empty() {
            let sometext = match f.read_line() {
                Ok(l) => l,
                Err(_) => {
                    // end of file (ierr /= 0); sometext keeps its previous value in the Fortran; the
                    // keywords below have already been processed for it.
                    break;
                }
            };
            let w = parse(&sometext, 11);
            if w.nwords() == 0 {
                continue;
            }
            let k = w.word(&sometext, 1);
            if k == "ncoeffs" {
                self.sym.symbol_value[iref + 1] = if w.nwords() >= 2 { w.word(&sometext, 2).to_string() } else { String::new() };
            } else if k == "atom-c" {
                centrec = w.span(&sometext, 2, w.nwords()).trim_end().to_string();
            } else if k == "atom-s" {
                centres = w.span(&sometext, 2, w.nwords()).trim_end().to_string();
            }
        }
        if centrec.is_empty() || centres.is_empty() {
            return crate::error::unsupported("mapgr: problems with the SHARM/SDF .dat file (get_nterms)");
        }
        self.sym.symbol_value[iref + 2] = if centrec != centres { "1".into() } else { "0".into() };
        Ok(0)
    }
}

/// Write an xyz file of the averaged central/second molecule as `sharm`, `sdf`
/// and `sdfcube` do.  `atoms` = (atsymbol, x, y, z, attype, centraldefinition).
pub fn write_mol_xyz(path: &std::path::Path, fname_xyz: &str, atoms: &[(String, f32, f32, f32, String, bool)]) -> Result<()> {
    let mut out = OutFile::create(path);
    out.line(&fmt("(i6)", &[FV::I(atoms.len() as i64)]));
    out.line(&format!("Central molecule coordinate file: {}", fname_xyz.trim_end()));
    for (s, x, y, z, t, cd) in atoms {
        let mut textline = fmt("(a3,1x,f10.5,1x,f10.5,1x,f10.5,1x,a3)", &[FV::Ch(s, 3), FV::R(*x), FV::R(*y), FV::R(*z), FV::Ch(t, 3)]);
        textline = textline.trim_end().to_string();
        if *cd {
            textline.push_str("  #CENTRALMOLECULE");
        }
        out.line(textline.trim_end());
    }
    out.close().map_err(|e| EpsrError::Io(format!("{}: {e}", path.display())))
}

/// List-directed integer as written by `write(symbol_value(i),*) n`.
pub fn list_int(n: i32) -> String {
    let mut w = ListWriter::new();
    w.int(n);
    w.finish()
}

/// The coefficient data read by `run_mapgr` / `run_mapgrentropy` (titles,
/// coefficients, smoothing).
pub struct Coeffs {
    pub lg: LnGamma,
    pub nexpan: i32,
    pub ident: i32,
    pub identsign: i32,
    pub l1ex: Vec<i32>,
    pub l2ex: Vec<i32>,
    pub lex: Vec<i32>,
    pub n1ex: Vec<i32>,
    pub n2ex: Vec<i32>,
    pub nr: i32,
    pub rd: F1<f32>,
    pub datac: F2<f32>,
    pub datas: F2<f32>,
}

impl Coeffs {
    /// Common first part of `run_mapgr` and `run_mapgrentropy` (symbols 1-4).
    pub fn read(st: &State) -> Result<Coeffs> {
        let log = st.logger();
        let bad = |what: &str| EpsrError::Input(format!("cannot read {what}"));
        let mut lg = LnGamma::default();
        lg.fctrln(250);
        let infile = st.sym.symbol_value[1].clone();
        let fname = st.sym.full(&infile).to_string_lossy().into_owned();
        let nexpan0 = read_int(&st.sym.symbol_value[2]).map_err(|_| bad("ncoeffs"))?;
        let mut ident = read_int(&st.sym.symbol_value[3]).map_err(|_| bad("ident"))?;
        let identsign;
        if ident.abs() > 0 {
            identsign = if ident >= 0 { 1 } else { -1 };
            ident = ident.abs();
        } else {
            identsign = 1;
            ident = 0;
        }
        let (nexpan, _titlex, titled) = r_int_titles(&fname, 3000, &log);
        let mut c = Coeffs { lg, nexpan, ident, identsign, l1ex: vec![], l2ex: vec![], lex: vec![], n1ex: vec![], n2ex: vec![], nr: 0, rd: F1::new(500), datac: F2::new(500, nexpan.max(nexpan0).max(1) as usize), datas: F2::new(500, nexpan.max(nexpan0).max(1) as usize) };
        for t in titled.iter().take(nexpan as usize) {
            let v = parse_title(t)?;
            c.l1ex.push(v[0]);
            c.l2ex.push(v[1]);
            c.lex.push(v[2]);
            c.n1ex.push(v[3]);
            c.n2ex.push(v[4]);
        }
        // r_int with the file name as left by r_int_titles
        match r_int(&fname, nexpan, &mut c.rd, &mut c.datac, &mut c.datas, &log)? {
            Some(n) => c.nr = n,
            None => return Err(EpsrError::Io(format!("cannot read coefficients file {fname}"))),
        }
        // smoothing
        let sometext = st.sym.symbol_value[4].clone();
        let w = parse(&sometext, 20);
        let mut rminsmoo: f32 = 0.0;
        let (mut nsmoord, mut nsmoosp, mut nsmooor) = (0i32, 0i32, 0i32);
        let mut r = crate::fortran::listread::ListRead::from_str(&sometext);
        if w.nwords() >= 4 {
            r.get_r(&mut rminsmoo).map_err(|_| bad("nsmoo"))?;
            r.get_i(&mut nsmoord).map_err(|_| bad("nsmoo"))?;
            r.get_i(&mut nsmoosp).map_err(|_| bad("nsmoo"))?;
            r.get_i(&mut nsmooor).map_err(|_| bad("nsmoo"))?;
        } else if w.nwords() >= 3 {
            r.get_i(&mut nsmoord).map_err(|_| bad("nsmoo"))?;
            r.get_i(&mut nsmoosp).map_err(|_| bad("nsmoo"))?;
            r.get_i(&mut nsmooor).map_err(|_| bad("nsmoo"))?;
            rminsmoo = 0.0;
        } else {
            r.get_i(&mut nsmoord).map_err(|_| bad("nsmoo"))?;
            rminsmoo = 0.0;
            nsmoosp = nsmoord;
            nsmooor = nsmoord;
        }
        let mut nrminsmoo = 1;
        while c.rd[nrminsmoo] < rminsmoo {
            nrminsmoo += 1;
        }
        let nr = c.nr;
        let mut tempc: F1<f32> = F1::new(500);
        let mut temps: F1<f32> = F1::new(500);
        for k in 1..=nexpan {
            let ku = (k - 1) as usize;
            let mut nsmoo = nsmoord;
            if c.l2ex[ku] == 0 && c.l1ex[ku] != 0 {
                nsmoo = nsmoosp;
            } else if c.l2ex[ku] != 0 && c.l1ex[ku] != 0 {
                nsmoo = nsmooor;
            }
            if nsmoo != 0 {
                for i in 1..=nr {
                    let fact = c.rd[i] * c.rd[i];
                    tempc[i] = c.datac[(i, k)] * fact;
                    temps[i] = c.datas[(i, k)] * fact;
                }
                smoosharm(nrminsmoo, nr, nsmoo, &mut tempc);
                smoosharm(nrminsmoo, nr, nsmoo, &mut temps);
                for i in 1..=nr {
                    let fact = c.rd[i] * c.rd[i];
                    c.datac[(i, k)] = tempc[i] / fact;
                    c.datas[(i, k)] = temps[i] / fact;
                }
            }
        }
        Ok(c)
    }

    /// The angular sum over coefficients of `run_mapgr`/`run_mapgrentropy`
    /// for coefficient k (1-based); returns (sumc, sums, spsumc, spsums).
    #[allow(clippy::too_many_arguments)]
    pub fn angular(&self, k: usize, only: [i32; 6], phil: f32, ctl: f32, phi: f32, cz: f32, chi: f32) -> (f32, f32, f32, f32) {
        let [l1only, l2only, lonly, n1only, n2only, m2fix] = only;
        let lg = &self.lg;
        let mut sumc: f32 = 0.0;
        let mut sums: f32 = 0.0;
        let mut spsumc: f32 = 0.0;
        let mut spsums: f32 = 0.0;
        let ll1 = self.l1ex[k - 1];
        let ll2 = self.l2ex[k - 1];
        let ll = self.lex[k - 1];
        let n1 = self.n1ex[k - 1];
        let n2 = self.n2ex[k - 1];
        let mut l1tmp = ll1;
        let mut l2tmp = ll2;
        let mut n1tmp = n1;
        let mut n2tmp = n2;
        let ok = |l1tmp: i32, l2tmp: i32, n1tmp: i32, n2tmp: i32| -> bool {
            let mut test = (l1only == 0 && l1tmp == 0) || l1only != 0;
            test = test && ((l2only == 0 && l2tmp == 0) || l2only != 0);
            test = test && ((lonly == 0 && ll == 0) || lonly != 0);
            test = test && ((n1only == 0 && n1tmp == 0) || n1only != 0);
            test && ((n2only == 0 && n2tmp == 0) || n2only != 0)
        };
        let ipow = |e: i32| -> f32 {
            if e % 2 == 0 {
                1.0
            } else {
                -1.0
            }
        };
        if ok(l1tmp, l2tmp, n1tmp, n2tmp) {
            let mut phs: f32 = 1.0;
            let phsn = ipow(l1tmp + l2tmp + ll + n1tmp + n2tmp);
            if self.identsign >= 0 {
                let (m2f, m2l) = if m2fix > 0 { (-l2tmp, l2tmp) } else { (0, 0) };
                sumoverm2(lg, l1tmp, l2tmp, ll, n1tmp, m2f, m2l, n2tmp, phil, ctl, phi, cz, chi, phs, &mut sumc, &mut sums);
                if n2 + n1 > 0 {
                    n1tmp = -n1;
                    n2tmp = -n2;
                    phs *= phsn;
                    sumoverm2(lg, l1tmp, l2tmp, ll, n1tmp, m2f, m2l, n2tmp, phil, ctl, phi, cz, chi, phs, &mut sumc, &mut sums);
                }
                if l2tmp == 0 {
                    spsumc = sumc;
                    spsums = sums;
                }
            } else {
                let (m2f, m2l) = if m2fix > 0 { (-l1tmp, l1tmp) } else { (0, 0) };
                phs = ipow(ll);
                sumoverm2(lg, l2tmp, l1tmp, ll, n2tmp, m2f, m2l, n1tmp, phil, ctl, phi, cz, chi, phs, &mut sumc, &mut sums);
                if n2 + n1 > 0 {
                    n1tmp = -n1;
                    n2tmp = -n2;
                    phs *= phsn;
                    sumoverm2(lg, l2tmp, l1tmp, ll, n2tmp, m2f, m2l, n1tmp, phil, ctl, phi, cz, chi, phs, &mut sumc, &mut sums);
                }
                if l1tmp == 0 {
                    spsumc = sumc;
                    spsums = sums;
                }
            }
        }
        if self.ident == 0 && ll1 != ll2 {
            let mut phs = ipow(ll1 + ll2);
            l1tmp = ll2;
            l2tmp = ll1;
            n1tmp = n2;
            n2tmp = n1;
            if ok(l1tmp, l2tmp, n1tmp, n2tmp) {
                let phsn = ipow(l1tmp + l2tmp + ll + n1tmp + n2tmp);
                let (m2f, m2l) = if m2fix > 0 { (-l2tmp, l2tmp) } else { (0, 0) };
                sumoverm2(lg, l1tmp, l2tmp, ll, n1tmp, m2f, m2l, n2tmp, phil, ctl, phi, cz, chi, phs, &mut sumc, &mut sums);
                if n2 + n1 > 0 {
                    n1tmp = -n2;
                    n2tmp = -n1;
                    phs *= phsn;
                    sumoverm2(lg, l1tmp, l2tmp, ll, n1tmp, m2f, m2l, n2tmp, phil, ctl, phi, cz, chi, phs, &mut sumc, &mut sums);
                }
            }
        }
        (sumc, sums, spsumc, spsums)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn special_functions() {
        let mut lg = LnGamma::default();
        lg.fctrln(40);
        // <1 1; 1 -1 | 0 0> = 1/sqrt(3)
        assert!((cgcoef(&lg, 1, 1, 0, 1, -1) - 1.0 / 3.0f32.sqrt()).abs() < 1e-6);
        // d^1_{00}(theta) = cos(theta)
        assert!((dlm(&lg, 1, 0, 0, 0.3) - 0.3).abs() < 1e-6);
        // d^1_{10}(theta) = -sin(theta)/sqrt(2)
        let c: f32 = 0.6;
        assert!((dlm(&lg, 1, 1, 0, c) + (1.0 - c * c).sqrt() / 2.0f32.sqrt()).abs() < 1e-6);
        assert_eq!(gf_cosd(60.0), 0.5);
        assert_eq!(gf_cosd(-180.0), -1.0);
        assert_eq!(gf_sind(30.0), 0.5);
        assert_eq!(gf_sind(-270.0), 1.0);
        assert_eq!(maxcgcoef(1, 1, 0, 1, -1), 4);
    }
    #[test]
    fn titles() {
        assert_eq!(parse_title("(2,2,4;2,-2)").unwrap(), [2, 2, 4, 2, -2]);
        assert_eq!(parse_title("(10,0,10;0,0)").unwrap(), [10, 0, 10, 0, 0]);
    }
}
