//! Crystalline samples: Bragg peaks, diffuse scattering and lattice refinement.
//!
//! Translation of
//! * the `qwidthqmax` / `mplicities` / `hklqmin` / `diffuse` / `nlatrefine` /
//!   `latrefine` keyword reading of `runepsr17` (`epsr_standalone_rev3.f`);
//! * the lattice-refinement loop of `runepsr17` around `atogen`/`gensq`,
//!   including the Bragg-peak block (`.inten` file, Debye-Waller factor,
//!   truncation convolution) and the fit to the data;
//! * modules `bragginten_mod`, `reciprocallattice_mod`, `diffuse_mod`
//!   (`initialise_diffuse_values`) and the routines `crystfrac`,
//!   `bragginten` (`.bragg` file), `sharm_diffuse` (`.DIFFUSE.h01` files) and
//!   `change_boxsize` of `ato_manipulation_rev4.f`;
//! * `convertphibthetacphic` (`auxiliary_routines.f`).
//!
//! In `bragginten`, `sharm_diffuse` and `change_boxsize` `cosd`/`sind` are not
//! declared `external`: gfortran uses its intrinsics (`_gfortran_cosd_r4`, ...),
//! and `atan2d(y,x)` is compiled to `atan2f(y,x)*57.2957763671875`.
//!
//! Fortran quirks reproduced:
//! * `sharm_diffuse` indexes `intenhkl(nhkl,npairs)` as `intenhkl(ip,ib)`
//!   (transposed).  The port reads the same storage element; configurations for
//!   which this runs past the end of the array (undefined behaviour in the
//!   Fortran, typically a segmentation fault) are reported as unsupported.
//! * `sharm_diffuse` passes the COMMON variable `nsums` to `outgofzs_sharm`,
//!   which increments it once per data set (and sets it to -1 when the previous
//!   `.DIFFUSE.h01` file cannot be read): this changes the accumulation mode of
//!   the running averages of the whole run.
//! * The run-time FORMAT for the coefficient titles in `sharm_diffuse`
//!   (`(1h(i1,1h,,...`) lacks a comma after `1h(`; gfortran >= 14 stops the
//!   program with "Missing comma between descriptors".  As for the `sharm` and
//!   `sdf` programs the port implements the intended format
//!   `(1h(,i1,1h,,...)` (see `tools/build_fortran_reference.sh`).
//! * In the lattice refinement `clatnew` is rescaled with `blatfac`, and when
//!   the last trial is rejected the box is restored without recomputing the
//!   structure factors.

use crate::ato::{perpd, scapro, vecpro};
use crate::epsr::Epsr;
use crate::error::{EpsrError, Result};
use crate::fortran::fmath::{cosd_intrinsic as cosd, nint, pow3, sind_intrinsic as sind, PI4ATAN};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::ListRead;
use crate::fortran::output::OutFile;
use crate::fortran::text::parse;
use crate::fortran::{F1, F2, F3};
use crate::neighbour::minimage;
use crate::state::State;
use crate::tools::makelattice::convertalphabetagamma;
use crate::tools::sharm::sharm_routines::{dlm, maxdlm, LnGamma, Sharm};

/// Run parameters (locals of `runepsr17`) and the data of the modules
/// `bragginten_mod`, `reciprocallattice_mod` and `diffuse_mod`.
#[derive(Clone, Debug, Default)]
pub struct Bragg {
    pub qbraggwidth: f32,
    pub qbraggmax: f32,
    pub multa: i32,
    pub multb: i32,
    pub multc: i32,
    pub usea: i32,
    pub useb: i32,
    pub usec: i32,
    pub hklqmin: f32,
    pub debyewaller: f32,
    pub nbragg: i32,
    pub nbraggup: i32,
    // diffuse_mod
    pub diffusedefined: bool,
    pub multdiffa: i32,
    pub multdiffb: i32,
    pub multdiffc: i32,
    pub ldiffmax: i32,
    pub ldiffstep: i32,
    pub mdiffstep: i32,
    // lattice refinement
    pub nlatrefine: i32,
    pub alatstep: f32,
    pub blatstep: f32,
    pub clatstep: f32,
    pub alphastep: f32,
    pub betastep: f32,
    pub gamstep: f32,
    pub refinelattice: bool,
    // bragginten_mod
    pub nhkl: i32,
    pub hihkl: Vec<i32>,
    pub kihkl: Vec<i32>,
    pub lihkl: Vec<i32>,
    pub qhkl: Vec<f32>,
    /// `intenhkl(nhkl,npairs)`, column-major storage.
    pub intenhkl: Vec<f32>,
    // reciprocallattice_mod
    pub recipra: [f32; 3],
    pub reciprb: [f32; 3],
    pub reciprc: [f32; 3],
    pub rbragg: [f32; 3],
    pub volrecipr: f32,
    pub ampra: f32,
    pub amprb: f32,
    pub amprc: f32,
}

fn need<T>(v: Option<T>, what: &str) -> Result<T> {
    v.ok_or_else(|| EpsrError::Input(format!("end of record reading {what}")))
}

/// `subroutine divvec(b,a)`.
fn divvec(b: &mut [f32; 3], a: f32) {
    if a != 0.0 {
        for v in b.iter_mut() {
            *v /= a;
        }
    }
}

/// `subroutine mulvec(b,a)`.
fn mulvec(b: &mut [f32; 3], a: f32) {
    for v in b.iter_mut() {
        *v *= a;
    }
}

/// `subroutine convertphibthetacphic(phib,thetac,phic,alpha,beta,gam)`.
pub fn convertphibthetacphic(phib: f32, thetac: f32, phic: f32) -> (f32, f32, f32) {
    let pi = PI4ATAN;
    let radconv = pi / 180.0;
    let degconv = 180.0 / pi;
    let gam = phib;
    let phicr = phic * radconv;
    let thetacr = thetac * radconv;
    let cxy = thetacr.sin();
    let cx = cxy * phicr.cos();
    let (alpha, beta);
    if cx.abs() < 1.0e-4 {
        alpha = 90.0;
        beta = 90.0;
    } else {
        let betar = cx.acos();
        beta = betar * degconv;
        let phibr = phib * radconv;
        let alphar = (cxy * (phibr - phicr).cos()).acos();
        alpha = alphar * degconv;
    }
    (alpha, beta, gam)
}

/// `subroutine smood(nqmin,nqmax,nsmoo,diff)`: 3-point smoothing keeping the area.
fn smood(nqmin: i32, nqmax: i32, nsmoo: i32, diff: &mut F1<f32>) {
    let n1 = nqmin;
    let n2 = nqmax - 1;
    if n2 < n1 || nsmoo == 0 {
        return;
    }
    for _ in 1..=nsmoo {
        let mut sm = diff[n1];
        diff[n2 + 1] = diff[n2];
        let mut sum1 = 0.5f32 * sm;
        for i in n1..=n2 {
            let j = i + 1;
            let sn = diff[j];
            let sum2 = 0.25f32 * (sm + sn);
            diff[i] = sum1 + sum2;
            sm = sn;
            sum1 = sum2;
        }
    }
}

impl Bragg {
    /// Keyword reading of `runepsr17` from `qwidthqmax` to `latrefine`.  On
    /// entry `iref` is the `revlorch` line; on exit the `latrefine` line.
    pub fn read_keywords(st: &mut State, iref: &mut i32, nq: i32, qstep: f32) -> Result<Bragg> {
        let mut b = Bragg::default();
        let sv = |st: &State, k: i32| st.sym.symbol_value[k].clone();
        // qwidthqmax
        *iref += 1;
        {
            let mut r = ListRead::from_str(&sv(st, *iref));
            b.qbraggwidth = need(r.real()?, "qwidthqmax")?;
            b.qbraggmax = need(r.real()?, "qwidthqmax")?;
        }
        // mplicities
        *iref += 1;
        {
            let t = sv(st, *iref);
            let w = parse(&t, 5);
            let mut r = ListRead::from_str(&w.span(&t, 1, 3));
            b.multa = need(r.int()?, "mplicities")?;
            b.multb = need(r.int()?, "mplicities")?;
            b.multc = need(r.int()?, "mplicities")?;
            if w.nwords() > 5 {
                let mut r = ListRead::from_str(&w.span(&t, 4, 6));
                b.usea = need(r.int()?, "mplicities")?;
                b.useb = need(r.int()?, "mplicities")?;
                b.usec = need(r.int()?, "mplicities")?;
                if b.usea != 0 {
                    b.usea = 1;
                }
                if b.useb != 0 {
                    b.useb = 1;
                }
                if b.usec != 0 {
                    b.usec = 1;
                }
            } else {
                b.usea = 1;
                b.useb = 1;
                b.usec = 1;
            }
        }
        // hklqmin
        *iref += 1;
        {
            let mut r = ListRead::from_str(&sv(st, *iref));
            b.hklqmin = need(r.real()?, "hklqmin")?;
            st.inp.rminbragg = need(r.real()?, "hklqmin")?;
            b.debyewaller = need(r.real()?, "hklqmin")?;
        }
        let rminbragg = st.inp.rminbragg;
        if rminbragg <= 0.0 {
            b.qbraggmax = 0.0;
        }
        b.nbragg = 0;
        b.nbraggup = 0;
        if b.qbraggmax > 0.0 {
            b.nbragg = nint(b.qbraggmax / qstep) + 1;
            b.nbraggup = nint((b.qbraggmax + 3.0 / rminbragg) / qstep) + 1;
            if b.nbraggup > nq {
                b.nbraggup = nq;
            }
        }
        // diffuse
        *iref += 1;
        {
            let t = sv(st, *iref);
            let w = parse(&t, 5);
            b.diffusedefined = w.nwords() > 5;
            if b.diffusedefined {
                let mut r = ListRead::from_str(&t);
                b.multdiffa = need(r.int()?, "diffuse")?;
                b.multdiffb = need(r.int()?, "diffuse")?;
                b.multdiffc = need(r.int()?, "diffuse")?;
                b.ldiffmax = need(r.int()?, "diffuse")?;
                b.ldiffstep = need(r.int()?, "diffuse")?;
                b.mdiffstep = need(r.int()?, "diffuse")?;
                if b.ldiffmax > 40 {
                    b.ldiffmax = 40;
                }
            }
        }
        // nlatrefine, latrefine
        *iref += 1;
        {
            let mut r = ListRead::from_str(&sv(st, *iref));
            b.nlatrefine = need(r.int()?, "nlatrefine")?;
        }
        *iref += 1;
        {
            let t = sv(st, *iref);
            let w = parse(&t, 5);
            if w.nwords() > 5 && b.nlatrefine > 0 {
                let mut r = ListRead::from_str(&w.span(&t, 1, 6));
                let mut vals = [0.0f32; 6];
                let mut ok = true;
                for v in vals.iter_mut() {
                    match r.real() {
                        Ok(Some(x)) => *v = x,
                        _ => {
                            ok = false;
                            break;
                        }
                    }
                }
                if ok {
                    b.alatstep = vals[0];
                    b.blatstep = vals[1];
                    b.clatstep = vals[2];
                    b.alphastep = vals[3];
                    b.betastep = vals[4];
                    b.gamstep = vals[5];
                } else {
                    // NB: values read before the error are kept by the Fortran
                    // only until they are reset here.
                    b.nlatrefine = 0;
                }
            }
        }
        b.refinelattice = b.nlatrefine > 0;
        if b.diffusedefined && b.qbraggmax > 0.0 {
            if b.multdiffa == 0 || b.multdiffb == 0 || b.multdiffc == 0 {
                return crate::error::unsupported("diffuse scattering with zero unit-cell multiplicity (integer division by zero)");
            }
            if b.ldiffstep <= 0 || b.mdiffstep <= 0 {
                return crate::error::unsupported("diffuse scattering with zero or negative l/m step (zero DO-loop step)");
            }
        }
        Ok(b)
    }

    /// `subroutine bragginten(qstep,qmax,qwidthin,fwhmq,multa,multb,multc,usea,useb,usec,hklqmin,nbragg,qbragg,intenbragg)`.
    /// Called with the coordinates in Å; writes the `.bragg` file.
    #[allow(clippy::too_many_arguments)]
    pub fn bragginten(&mut self, st: &State, qstep: f32, qmax: f32, qwidthin: f32, fwhmq: f32, nbragg: i32, qbragg: &F1<f32>, intenbragg: &mut F2<f32>) -> Result<()> {
        let (multa, multb, multc) = (self.multa, self.multb, self.multc);
        let (usea, useb, usec) = (self.usea, self.useb, self.usec);
        let hklqmin = self.hklqmin;
        let m = &st.mol;
        let l = &st.lat;
        let npairs = m.npairs;
        let ntype = m.ntype;
        st.log("bragginten> num_threads_to be used:            1");
        let pi = PI4ATAN;
        let twopi = 2.0f32 * pi;
        // reciprocal lattice vectors
        self.recipra = vecpro(&l.blatc, &l.clatc);
        self.reciprb = vecpro(&l.clatc, &l.alatc);
        self.reciprc = vecpro(&l.alatc, &l.blatc);
        let mut vol = scapro(&l.alatc, &self.recipra);
        divvec(&mut self.recipra, vol);
        divvec(&mut self.reciprb, vol);
        divvec(&mut self.reciprc, vol);
        mulvec(&mut self.recipra, multa as f32 * twopi);
        mulvec(&mut self.reciprb, multb as f32 * twopi);
        mulvec(&mut self.reciprc, multc as f32 * twopi);
        let sqrmoda = scapro(&self.recipra, &self.recipra);
        let sqrmodb = scapro(&self.reciprb, &self.reciprb);
        let sqrmodc = scapro(&self.reciprc, &self.reciprc);
        vol /= (multa * multb * multc) as f32;
        self.rbragg = vecpro(&self.reciprb, &self.reciprc);
        self.volrecipr = scapro(&self.recipra, &self.rbragg);
        self.ampra = sqrmoda.sqrt();
        self.amprb = sqrmodb.sqrt();
        self.amprc = sqrmodc.sqrt();
        let (ra, rb, rc) = (self.recipra, self.reciprb, self.reciprc);
        let mut qwidth = qwidthin.abs() + qmax * fwhmq;
        let qwidthfac: f32 = 10.0;
        let qlim = qmax + qwidthfac * qwidth;
        let hmax = usea * (nint(qlim / ra[0].abs()) + 1);
        let kmax0 = useb * (nint(qlim / rb[1].abs()) + 1);
        let lmax0 = usec * (nint(qlim / rc[2].abs()) + 1);
        let nhklmax = (hmax + 1) * (2 * kmax0 + 1) * (2 * lmax0 + 1);
        st.log(&format!(" Maximum number of hkl values {nhklmax:>12}"));
        self.hihkl.clear();
        self.kihkl.clear();
        self.lihkl.clear();
        self.qhkl.clear();
        let mut hi = 0;
        while hi <= hmax {
            let kmax = useb * (nint((qlim - hi as f32 * ra[1]) / rb[1].abs()) + 1);
            let kmin = if kmax > 0 { -kmax } else { 0 };
            let mut ki = kmin;
            while ki <= kmax {
                let lmax = usec * (nint(((qlim - hi as f32 * ra[2]) - ki as f32 * rb[2]) / rc[2].abs()) + 1);
                let lmin = if lmax > 0 { -lmax } else { 0 };
                let mut li = lmin;
                while li <= lmax {
                    for it in 0..3 {
                        self.rbragg[it] = (hi as f32 * ra[it] + ki as f32 * rb[it]) + li as f32 * rc[it];
                    }
                    let qhkltest = scapro(&self.rbragg, &self.rbragg).sqrt();
                    if qhkltest <= qlim {
                        if self.hihkl.len() as i32 >= nhklmax {
                            return crate::error::unsupported("bragginten: more hkl values than allocated (nhklmax) - out-of-bounds writes in the Fortran");
                        }
                        self.hihkl.push(hi);
                        self.kihkl.push(ki);
                        self.lihkl.push(li);
                        self.qhkl.push(qhkltest);
                    }
                    li += 1;
                }
                ki += 1;
            }
            hi += 1;
        }
        let nhkl = self.hihkl.len() as i32;
        self.nhkl = nhkl;
        st.log(&format!(" Total number of hkl values {nhkl:>12}"));
        let nh = nhkl as usize;
        // unique Q values
        let mut qunique: Vec<f32> = Vec::new();
        let mut nqunique: Vec<i32> = Vec::new();
        for ib in 0..nh {
            let qhkltest = self.qhkl[ib];
            let mut found = None;
            for (iq, q) in qunique.iter().enumerate() {
                if (qhkltest - *q).abs() < 1.0e-5 {
                    found = Some(iq);
                    break;
                }
            }
            match found {
                None => {
                    qunique.push(qhkltest);
                    nqunique.push(1);
                }
                Some(iq) => nqunique[iq] += 1,
            }
        }
        let ntotalqunique = qunique.len();
        let mut ibvalues: Vec<Vec<usize>> = vec![Vec::new(); ntotalqunique];
        for ib in 0..nh {
            let qhkltest = self.qhkl[ib];
            for (iq, q) in qunique.iter().enumerate() {
                if (qhkltest - *q).abs() < 1.0e-5 {
                    ibvalues[iq].push(ib);
                    break;
                }
            }
        }
        // fractional coordinates
        let natom = m.natom;
        let mut fa = F1::<f32>::new(natom.max(1) as usize);
        let mut fb = F1::<f32>::new(natom.max(1) as usize);
        let mut fc = F1::<f32>::new(natom.max(1) as usize);
        for ia in 1..=natom {
            let im = m.mref[ia];
            let xuse = m.xa[ia] + m.xadd[im];
            let yuse = m.ya[ia] + m.yadd[im];
            let zuse = m.za[ia] + m.zadd[im];
            let (mut a, mut b, mut c);
            if l.cubic {
                c = zuse / l.clatz;
                b = yuse / l.blaty;
                a = xuse / l.alatx;
            } else {
                c = zuse / l.clatz;
                b = (yuse - c * l.claty) / l.blaty;
                a = (xuse - b * l.blatx - c * l.clatx) / l.alatx;
            }
            minimage(&mut a, &mut b, &mut c);
            fa[ia] = a;
            fb[ia] = b;
            fc[ia] = c;
        }
        // amplitudes per atom type (sum over the atoms of each type in atom order)
        let nt = ntype.max(1) as usize;
        let mut ramp = F2::<f32>::new(nh.max(1), nt);
        let mut iamp = F2::<f32>::new(nh.max(1), nt);
        for ia in 1..=natom {
            let it = m.itype[ia];
            for ib in 0..nh {
                let huse = (self.hihkl[ib] * multa) as f32;
                let kuse = (self.kihkl[ib] * multb) as f32;
                let luse = (self.lihkl[ib] * multc) as f32;
                let arg = twopi * ((fa[ia] * huse + fb[ia] * kuse) + fc[ia] * luse);
                let ramptemp = arg.cos();
                let iamptemp = -arg.sin();
                ramp[(ib + 1, it)] += ramptemp;
                iamp[(ib + 1, it)] += iamptemp;
            }
        }
        let arg = 1.0f32 / m.natomr as f32;
        let mut intenhkl = vec![0.0f32; nh * npairs.max(0) as usize];
        for ip in 1..=npairs {
            let it = m.n1extr[ip];
            let jt = m.n2extr[ip];
            for ib in 1..=nh {
                let k = (ib - 1) + (ip as usize - 1) * nh;
                let mut v = ramp[(ib, it)] * ramp[(ib, jt)] + iamp[(ib, it)] * iamp[(ib, jt)];
                v *= arg;
                if it != jt {
                    v *= 2.0;
                }
                intenhkl[k] = v;
            }
        }
        let mut intenqunique = F2::<f32>::new(npairs.max(1) as usize, ntotalqunique.max(1));
        for ip in 1..=npairs {
            for iq in 0..ntotalqunique {
                for &ib in &ibvalues[iq] {
                    let v = intenhkl[ib + (ip as usize - 1) * nh];
                    if self.hihkl[ib] > 0 {
                        intenqunique[(ip, iq + 1)] += 2.0 * v;
                    } else {
                        intenqunique[(ip, iq + 1)] += v;
                    }
                }
            }
        }
        for ip in 1..=npairs {
            for ib in 1..=nbragg {
                intenbragg[(ip, ib)] = 0.0;
            }
        }
        // broadening
        let p15 = twopi.powf(1.5);
        let gfactor = self.volrecipr / (p15 * 2.0);
        let gfactor0 = self.volrecipr / p15;
        let lfactor = self.volrecipr / (pi * pi);
        for iq in 0..ntotalqunique {
            let qhkltemp = qunique[iq];
            qwidth = qwidthin.abs() + qhkltemp * fwhmq;
            let mut nfirst = if qwidthin < 0.0 { 1 } else { 1.max(nint((qhkltemp - 15.0 * qwidth) / qstep) - 1) };
            let mut nlast = nbragg.min(nint((qhkltemp + 15.0 * qwidth) / qstep) + 1);
            if qhkltemp < hklqmin {
                nfirst = 1;
                nlast = 0;
            }
            let mut it = nfirst - 1;
            while it < nlast {
                it += 1;
                let qb = qbragg[it];
                let mut sum: f32 = 0.0;
                if qwidthin > 0.0 {
                    if qhkltemp > 0.0 && qb > 0.0 {
                        let diff = (qhkltemp - qb) / qwidth;
                        let add = (qhkltemp + qb) / qwidth;
                        if diff.abs() < 8.0 {
                            sum = (-(0.5 * diff * diff)).exp();
                        }
                        if add.abs() < 8.0 {
                            sum -= (-(0.5 * add * add)).exp();
                        }
                        sum = sum * gfactor / (qhkltemp * qb * qwidth);
                    } else {
                        let add = qhkltemp / qwidth;
                        let diff = qb / qwidth;
                        sum = 0.5 * (diff * diff + add * add);
                        if sum < 32.0 {
                            sum = gfactor0 * (-sum).exp() / pow3(qwidth);
                        } else {
                            sum = 0.0;
                        }
                    }
                } else {
                    let diff = (qhkltemp - qb) / qwidth;
                    let add = (qhkltemp + qb) / qwidth;
                    sum = lfactor / ((1.0 + diff * diff) * (1.0 + add * add) * pow3(qwidth));
                }
                for ip in 1..=npairs {
                    intenbragg[(ip, it)] += intenqunique[(ip, iq + 1)] * sum;
                }
            }
        }
        // subtract the single-atom scattering, normalise per atom
        let mut ip = 0;
        for it in 1..=ntype {
            for jt in it..=ntype {
                ip += 1;
                if it == jt {
                    let cprod = m.cfraccg[it] * m.cfraccg[jt];
                    for ib in 1..=nbragg {
                        intenbragg[(ip, ib)] = if cprod > 0.0 { (intenbragg[(ip, ib)] - m.cfraccg[it]) / cprod } else { 0.0 };
                    }
                } else {
                    let cprod = 2.0 * m.cfraccg[it] * m.cfraccg[jt];
                    for ib in 1..=nbragg {
                        intenbragg[(ip, ib)] = if cprod > 0.0 { intenbragg[(ip, ib)] / cprod } else { 0.0 };
                    }
                }
            }
        }
        // .bragg file
        let fname = format!("{}{}.bragg", st.sym.directory_name, st.sym.fname_run);
        let mut out = OutFile::create(std::path::Path::new(&fname));
        let f101 = "(1x,a,4(1x,e12.5))";
        out.line(&ListWriter::new().chars("# Multiplicities along a, b and c: ").int(multa).int(multb).int(multc).finish());
        out.line(&fmt(f101, &[FV::S("# Volume of unit cell: "), FV::R(vol)]));
        out.line(&fmt(f101, &[FV::S("# Volume of reciprocal lattice: "), FV::R(self.volrecipr)]));
        for (t, v) in [("# Lattice vectors a: ", &l.alatc), ("# Lattice vectors b: ", &l.blatc), ("# Lattice vectors c: ", &l.clatc)] {
            out.line(&fmt(f101, &[FV::S(t), FV::R(v[0]), FV::R(v[1]), FV::R(v[2]), FV::R(scapro(v, v).sqrt())]));
        }
        for (t, v, a) in [
            ("# Reciprocal lattice vectors a: ", &self.recipra, self.ampra),
            ("# Reciprocal lattice vectors b: ", &self.reciprb, self.amprb),
            ("# Reciprocal lattice vectors c: ", &self.reciprc, self.amprc),
        ] {
            out.line(&fmt(f101, &[FV::S(t), FV::R(v[0]), FV::R(v[1]), FV::R(v[2]), FV::R(a)]));
        }
        out.line(&fmt("(a,i0)", &[FV::S("# Number of Bragg refections calculated: "), FV::I(nhkl as i64)]));
        out.line(&fmt("(1x,\"# \",a)", &[FV::S("Actual Bragg Q values")]));
        out.line(&fmt("(1x,\"# \",a)", &[FV::S("Unique Bragg Q values and multiplicities")]));
        for iq in 0..ntotalqunique {
            out.line(&fmt("(1x,\"#\",i6,1x,f10.5,1x,i0)", &[FV::I(iq as i64 + 1), FV::R(qunique[iq]), FV::I(nqunique[iq] as i64)]));
        }
        out.close()?;
        if self.diffusedefined {
            self.intenhkl = intenhkl;
        } else {
            self.hihkl.clear();
            self.kihkl.clear();
            self.lihkl.clear();
            self.qhkl.clear();
            self.intenhkl.clear();
        }
        Ok(())
    }

    /// `subroutine sharm_diffuse(qstep,nq,nsums,nfile)`: spherical-harmonic
    /// coefficients of the non-Bragg reflections, accumulated in
    /// `<run>dNN.DIFFUSE.h01`.
    pub fn sharm_diffuse(&mut self, st: &mut State, qstep: f32, nq: i32) -> Result<()> {
        if !self.diffusedefined || self.nhkl <= 0 {
            return Ok(());
        }
        let nhkl = self.nhkl as usize;
        let npairs = st.mol.npairs;
        let ntype = st.mol.ntype;
        let nfile = st.dist.nfile;
        // NB: intenhkl(nhkl,npairs) is indexed intenhkl(ip,ib) in the Fortran.
        let total = self.intenhkl.len();
        let tidx = |ip: i32, ib: usize| (ip as usize - 1) + (ib - 1) * nhkl;
        if npairs > 0 && tidx(npairs, nhkl) >= total {
            return crate::error::unsupported(
                "sharm_diffuse: intenhkl(ip,ib) runs past the end of intenhkl(nhkl,npairs) (nhkl > npairs) - undefined behaviour in the Fortran",
            );
        }
        let mut ip = 0;
        for it in 1..=ntype {
            for jt in it..=ntype {
                ip += 1;
                let mut cprod = st.mol.cfrac[it] * st.mol.cfrac[jt];
                if it != jt {
                    cprod *= 2.0;
                }
                for ib in 1..=nhkl {
                    let k = tidx(ip, ib);
                    self.intenhkl[k] /= cprod;
                }
            }
        }
        // initialise_diffuse_values(nfile)
        let mut maxlngamma = 0;
        let mut l = 0;
        while l <= self.ldiffmax {
            let mdiffmax = if self.mdiffstep == 0 { 0 } else { l };
            for m in 0..=mdiffmax {
                maxlngamma = maxlngamma.max(maxdlm(l, m, 0)).max(maxdlm(l, -m, 0));
            }
            l += self.ldiffstep;
        }
        let mut lg = LnGamma::default();
        lg.fctrln(maxlngamma);
        let mut l1ex: Vec<i32> = Vec::new();
        let mut n1ex: Vec<i32> = Vec::new();
        let mut l = 0;
        while l <= self.ldiffmax {
            let mdiffmax = if self.mdiffstep == 0 { 0 } else { l };
            let mut m = 0;
            while m <= mdiffmax {
                l1ex.push(l);
                n1ex.push(m);
                m += self.mdiffstep;
            }
            l += self.ldiffstep;
        }
        let ndiffuse = l1ex.len();
        let mut qmaxdiff: f32 = 0.0;
        for ib in 0..nhkl {
            if self.qhkl[ib] > qmaxdiff {
                qmaxdiff = self.qhkl[ib];
            }
        }
        let qdiffstep = self.ampra.max(self.amprb).max(self.amprc);
        let nqdiff = nint(qmaxdiff / qdiffstep) + 1;
        let nqd = nqdiff as usize;
        let mut qdiff = F1::<f32>::new(nqd);
        let mut diffcoeffc = F3::<f32>::new(ndiffuse.max(1), nqd, nfile.max(1) as usize);
        let mut diffcoeffs = F3::<f32>::new(ndiffuse.max(1), nqd, nfile.max(1) as usize);
        let mut ndiffcoeff = F1::<i32>::new(nqd);
        for iq in 1..=nqd {
            qdiff[iq] = (iq - 1) as f32 * qdiffstep;
        }
        // contributions of the non-Bragg reflections
        let mut intenfile = F1::<f32>::new(nfile.max(1) as usize);
        for ib in 1..=nhkl {
            let hi = self.hihkl[ib - 1];
            let ki = self.kihkl[ib - 1];
            let li = self.lihkl[ib - 1];
            if (hi / self.multdiffa) * self.multdiffa != hi || (ki / self.multdiffb) * self.multdiffb != ki || (li / self.multdiffc) * self.multdiffc != li {
                let qh = self.qhkl[ib - 1];
                let iqdiff = nint(qh / qdiffstep) + 1;
                for it in 0..3 {
                    self.rbragg[it] = (hi as f32 * self.recipra[it] + ki as f32 * self.reciprb[it]) + li as f32 * self.reciprc[it];
                }
                let xbragg = self.rbragg[0];
                let ybragg = self.rbragg[1];
                let zbragg = self.rbragg[2];
                let ctl = if qh > 0.0 { zbragg / qh } else { 1.0 };
                let phil = ybragg.atan2(xbragg) * 5.729_577_636_718_75e1_f32;
                for iff in 1..=nfile {
                    let ix = st.wts.xfile[iff];
                    intenfile[iff] = 0.0;
                    if ix > 0 {
                        let iq = nint(qh / qstep) + 1;
                        if iq > 0 && iq <= nq {
                            for ip in 1..=npairs {
                                intenfile[iff] = (intenfile[iff] as f64 + st.xray.weightx[(iq, ip, ix)] * self.intenhkl[tidx(ip, ib)] as f64) as f32;
                            }
                        }
                    } else {
                        for ip in 1..=npairs {
                            intenfile[iff] = (intenfile[iff] as f64 + st.wts.cwt[(iff, ip)] * self.intenhkl[tidx(ip, ib)] as f64) as f32;
                        }
                    }
                }
                for id in 1..=ndiffuse {
                    let l = l1ex[id - 1];
                    let m = n1ex[id - 1];
                    let factor = (2 * l + 1) as f32;
                    let mut dlmcoeff = dlm(&lg, l, m, 0, ctl) * factor;
                    let mut dlmcoeffc = cosd(m as f32 * phil) * dlmcoeff;
                    let mut dlmcoeffs = sind(m as f32 * phil) * dlmcoeff;
                    if hi > 0 {
                        dlmcoeff = dlm(&lg, l, m, 0, -ctl) * factor;
                        dlmcoeffc = cosd(m as f32 * (phil + 180.0)) * dlmcoeff + dlmcoeffc;
                        dlmcoeffs = sind(m as f32 * (phil + 180.0)) * dlmcoeff + dlmcoeffs;
                    }
                    for iff in 1..=nfile {
                        diffcoeffc[(id, iqdiff, iff)] += dlmcoeffc * intenfile[iff];
                        diffcoeffs[(id, iqdiff, iff)] += dlmcoeffs * intenfile[iff];
                    }
                }
                if hi.abs() > 0 {
                    ndiffcoeff[iqdiff] += 2;
                } else {
                    ndiffcoeff[iqdiff] += 1;
                }
            }
        }
        // titles: run-time format (1h(,i1,1h,,i1,1h,,i1,1h;,i1,1h,,i1,1h)) (see module notes)
        let titledat: Vec<String> = (0..ndiffuse).map(|id| format!("({},0,{};{},0)", l1ex[id], l1ex[id], n1ex[id])).collect();
        let titler = "   Q";
        let fname_run = st.sym.fname_run.clone();
        let ncharoutput = fname_run.find(".EPSR").map(|p| p as i64).unwrap_or(-1);
        let outputfile = if ncharoutput >= 0 { fname_run[..ncharoutput as usize].to_string() } else { String::new() };
        let log = st.logger();
        for iff in 1..=nfile {
            let diffusename = format!("d{:02}", iff);
            let fname = format!("{}{}{}.DIFFUSE.h01", st.sym.directory_name, outputfile, diffusename);
            for iq in 1..=nqd {
                if ndiffcoeff[iq] > 0 {
                    for id in 1..=ndiffuse {
                        diffcoeffc[(id, iq, iff)] /= ndiffcoeff[iq] as f32;
                        diffcoeffs[(id, iq, iff)] /= ndiffcoeff[iq] as f32;
                    }
                }
            }
            let mut sh = Sharm { nz: nqdiff, ndist: ndiffuse as i32, ..Default::default() };
            sh.z = qdiff.clone();
            sh.gofz = F2::new(nqd, ndiffuse.max(1));
            sh.gofzs = F2::new(nqd, ndiffuse.max(1));
            sh.gzt = F2::new(nqd, ndiffuse.max(1));
            sh.gzt2 = F2::new(nqd, ndiffuse.max(1));
            let mut tempdiff = F1::<f32>::new(nqd);
            for id in 1..=ndiffuse {
                for iq in 1..=nqd {
                    tempdiff[iq] = diffcoeffc[(id, iq, iff)];
                }
                smood(1, nqdiff, 1, &mut tempdiff);
                for iq in 1..=nqd {
                    sh.gofz[(iq, id)] = tempdiff[iq];
                }
                for iq in 1..=nqd {
                    tempdiff[iq] = diffcoeffs[(id, iq, iff)];
                }
                smood(1, nqdiff, 1, &mut tempdiff);
                for iq in 1..=nqd {
                    sh.gofzs[(iq, id)] = tempdiff[iq];
                }
            }
            sh.titledat = titledat.clone();
            // NB: nsums is the COMMON variable of runepsr17 (modified here).
            sh.outgofzs_sharm(&fname, "h", &mut st.inp.nsums, titler, &log)?;
            // r_int overwrites z with the Q values read from the previous file
            qdiff = sh.z.clone();
        }
        self.qhkl.clear();
        self.intenhkl.clear();
        Ok(())
    }
}

/// `subroutine change_boxsize(afac,bfac,cfac,newphib,newthetac,newphic)` of
/// `ato_manipulation_rev4.f` (gfortran intrinsic `cosd`/`sind`).
pub fn change_boxsize(st: &mut State, afac: f32, bfac: f32, cfac: f32, newphib: f32, newthetac: f32, newphic: f32) {
    let frac = |l: &crate::state::Lattice, x: &mut f32, y: &mut f32, z: &mut f32| {
        let mut fc = *z / l.clatz;
        let mut fb = (*y - fc * l.claty) / l.blaty;
        let mut fa = (*x - fb * l.blatx - fc * l.clatx) / l.alatx;
        minimage(&mut fa, &mut fb, &mut fc);
        *z = fc;
        *y = fb;
        *x = fa;
    };
    let cart = |l: &crate::state::Lattice, x: &mut f32, y: &mut f32, z: &mut f32| {
        let (fa, fb, fc) = (*x, *y, *z);
        *z = fc * l.clatz;
        *y = fb * l.blaty + fc * l.claty;
        *x = fa * l.alatx + fb * l.blatx + fc * l.clatx;
    };
    {
        let (m, l) = (&mut st.mol, &mut st.lat);
        for im in 1..=m.nmol {
            let (mut x, mut y, mut z) = (m.xadd[im], m.yadd[im], m.zadd[im]);
            frac(l, &mut x, &mut y, &mut z);
            (m.xadd[im], m.yadd[im], m.zadd[im]) = (x, y, z);
            let (mut x, mut y, mut z) = (m.phix[im], m.phiy[im], m.phiz[im]);
            frac(l, &mut x, &mut y, &mut z);
            (m.phix[im], m.phiy[im], m.phiz[im]) = (x, y, z);
        }
        l.alat *= afac;
        l.blat *= bfac;
        l.clat *= cfac;
        l.phib = newphib;
        l.thetac = newthetac;
        l.phic = newphic;
        l.alatx = l.alat;
        l.blatx = cosd(l.phib) * l.blat;
        l.blaty = sind(l.phib) * l.blat;
        l.clatx = (sind(l.thetac) * l.clat) * cosd(l.phic);
        l.claty = (sind(l.thetac) * l.clat) * sind(l.phic);
        l.clatz = cosd(l.thetac) * l.clat;
        for im in 1..=m.nmol {
            let (mut x, mut y, mut z) = (m.xadd[im], m.yadd[im], m.zadd[im]);
            cart(l, &mut x, &mut y, &mut z);
            (m.xadd[im], m.yadd[im], m.zadd[im]) = (x, y, z);
            let (mut x, mut y, mut z) = (m.phix[im], m.phiy[im], m.phiz[im]);
            cart(l, &mut x, &mut y, &mut z);
            (m.phix[im], m.phiy[im], m.phiz[im]) = (x, y, z);
        }
        l.alatc = [l.alatx, 0.0, 0.0];
        l.blatc = [l.blatx, l.blaty, 0.0];
        l.clatc = [l.clatx, l.claty, l.clatz];
        l.bxax = l.blatx / l.alatx;
        l.cxax = l.clatx / l.alatx;
        l.cyby = l.claty / l.blaty;
        let perpab = perpd(&l.alatc, &l.blatc, &l.clatc);
        let perpbc = perpd(&l.blatc, &l.clatc, &l.alatc);
        let perpca = perpd(&l.clatc, &l.alatc, &l.blatc);
        m.dmax = perpab.min(perpbc).min(perpca);
        st.dist.cell = m.dmax;
        l.cubic = !(l.alat != l.blat || l.blat != l.clat || l.clat != l.alat || l.phib != 90.0 || l.thetac != 0.0 || l.phic != 0.0);
        st.log(&ListWriter::new().real(perpab).real(perpbc).real(perpca).finish());
    }
    st.sync_lattice();
}

impl Epsr {
    /// The `do while (ilatrefine.lt.nlatrefine+1)` loop of `runepsr17`:
    /// `atogen` + `gensq`, the Bragg-peak correction of the intermolecular
    /// structure factors and the lattice-parameter refinement.
    pub fn structure_and_lattice(&mut self, nq: i32) -> Result<()> {
        let nlatrefine = self.bragg.nlatrefine;
        let refinelattice = self.bragg.refinelattice;
        let mut oldfit: f32 = 0.0;
        let mut alatold = self.st.lat.alat;
        let mut blatold = self.st.lat.blat;
        let mut clatold = self.st.lat.clat;
        let (mut alphaold, mut betaold, mut gamold) = convertphibthetacphic(self.st.lat.phib, self.st.lat.thetac, self.st.lat.phic);
        let mut alatnew = alatold;
        let mut blatnew = blatold;
        let mut clatnew = clatold;
        let mut alphanew = alphaold;
        let mut betanew = betaold;
        let mut gamnew = gamold;
        // newphib/newthetac/newphic are locals of runepsr17 that keep their value
        let (mut newphib, mut newthetac, mut newphic) = (0.0f32, 0.0f32, 0.0f32);
        let mut ilatrefine = 0;
        while ilatrefine < nlatrefine + 1 {
            let mut ngrs = self.ctl.ngrsamples;
            self.st.atogen(&mut self.wav, &mut ngrs);
            self.ctl.ngrsamples = ngrs;
            self.st.gensq(nq, &mut self.qwind, &self.wav);
            if self.bragg.nbragg > 0 {
                self.bragg_peaks(nq)?;
            }
            if refinelattice && ilatrefine <= nlatrefine {
                // a) fit to the data
                let nfile = self.st.dist.nfile;
                let npairs = self.st.mol.npairs;
                let mut rfacraw: f32 = 0.0;
                for inn in 1..=nfile {
                    let mut rfacdraw: f32 = 0.0;
                    let ix = self.st.wts.xfile[inn];
                    let ds = &self.data[(inn - 1) as usize];
                    for i in 1..=nq {
                        let mut sumraw: f32 = 0.0;
                        if ix > 0 {
                            let x = &self.st.xray;
                            for ic in 1..=npairs {
                                sumraw = (sumraw as f64 + x.weightx[(i, ic, ix)] * (self.qwind[(1, i, ic)] + self.qwind[(2, i, ic)]) as f64) as f32;
                            }
                            sumraw = (sumraw as f64 + x.backx[(i, ix)]) as f32;
                        } else {
                            let w = &self.st.wts;
                            for ic in 1..=npairs {
                                sumraw = ((sumraw as f64 + w.ciwt[(inn, ic)] * self.qwind[(1, i, ic)] as f64) + w.cwt[(inn, ic)] * self.qwind[(2, i, ic)] as f64) as f32;
                            }
                        }
                        let diffraw = if i >= ds.nqmind && i <= ds.nqmaxd { self.datai[(i, inn)] - sumraw } else { 0.0 };
                        if i >= self.st.dist.nqmin[inn] && i <= self.st.dist.nqmax[inn] {
                            rfacdraw = diffraw * diffraw + rfacdraw;
                        }
                    }
                    if self.st.dist.nqpts[inn] > 0 {
                        rfacdraw /= self.st.dist.nqpts[inn] as f32;
                    }
                    rfacraw += rfacdraw;
                }
                rfacraw /= nfile as f32;
                let newfit = rfacraw;
                let e6 = |v: f32| fmt("(1x,e12.6)", &[FV::R(v)]);
                if ilatrefine == 0 {
                    self.st.log("refinelattice> Starting values:-");
                    self.st.log(&format!(
                        "refinelattice> Initial lattice pars.: {}{}{}{}{}{}",
                        e6(alatold),
                        e6(blatold),
                        e6(clatold),
                        e6(alphaold),
                        e6(betaold),
                        e6(gamold)
                    ));
                    self.st.log(&format!("refinelattice> Initial fit: {}", e6(newfit)));
                } else {
                    self.st.log(&format!("refinelattice> End of iteration: {ilatrefine:>5}"));
                    self.st.log(&format!(
                        "refinelattice> Old lattice pars.: {}{}{}{}{}{}",
                        e6(alatold),
                        e6(blatold),
                        e6(clatold),
                        e6(alphaold),
                        e6(betaold),
                        e6(gamold)
                    ));
                    self.st.log(&format!(
                        "refinelattice> New lattice pars.: {}{}{}{}{}{}",
                        e6(alatnew),
                        e6(blatnew),
                        e6(clatnew),
                        e6(alphanew),
                        e6(betanew),
                        e6(gamnew)
                    ));
                    self.st.log(&format!("refinelattice> New fit and old fit: {}{}", e6(newfit), e6(oldfit)));
                }
                // b) accept or reject the trial lattice
                if ilatrefine == 0 {
                    oldfit = newfit;
                } else if newfit <= oldfit {
                    oldfit = newfit;
                    alatold = alatnew;
                    blatold = blatnew;
                    clatold = clatnew;
                    alphaold = alphanew;
                    betaold = betanew;
                    gamold = gamnew;
                    if ilatrefine == nlatrefine {
                        ilatrefine += 1;
                    }
                } else {
                    let alatfac = alatold / alatnew;
                    let blatfac = blatold / blatnew;
                    let clatfac = clatold / clatnew;
                    convertalphabetagamma(alphaold, betaold, gamold, &mut newphib, &mut newthetac, &mut newphic);
                    self.st.unnormalise_coordinates();
                    change_boxsize(&mut self.st, alatfac, blatfac, clatfac, newphib, newthetac, newphic);
                    self.st.normalise_coordinates();
                }
                // c) new trial lattice
                if ilatrefine < nlatrefine {
                    let b = &self.bragg;
                    let dmax = self.st.mol.dmax;
                    alatnew = ((self.st.rng.ran1() * 2.0 - 1.0) * b.alatstep) / dmax + alatold;
                    let mut alatfac = alatnew / alatold;
                    let (mut blatfac, mut clatfac);
                    if self.st.lat.cubic {
                        blatnew = alatnew;
                        clatnew = alatnew;
                        blatfac = alatfac;
                        clatfac = alatfac;
                        alphanew = 90.0;
                        betanew = 90.0;
                        gamnew = 90.0;
                    } else {
                        blatnew = ((self.st.rng.ran1() * 2.0 - 1.0) * b.blatstep) / dmax + blatold;
                        blatfac = blatnew / blatold;
                        clatnew = ((self.st.rng.ran1() * 2.0 - 1.0) * b.clatstep) / dmax + clatold;
                        clatfac = clatnew / clatold;
                        alphanew = (self.st.rng.ran1() * 2.0 - 1.0) * b.alphastep + alphaold;
                        betanew = (self.st.rng.ran1() * 2.0 - 1.0) * b.betastep + betaold;
                        gamnew = (self.st.rng.ran1() * 2.0 - 1.0) * b.gamstep + gamold;
                    }
                    convertalphabetagamma(alphanew, betanew, gamnew, &mut newphib, &mut newthetac, &mut newphic);
                    self.st.unnormalise_coordinates();
                    // keep the density constant
                    let natomr = self.st.mol.natomr as f32;
                    let l = &self.st.lat;
                    let rhoold = natomr / ((l.alatx * l.blaty) * l.clatz);
                    change_boxsize(&mut self.st, alatfac, blatfac, clatfac, newphib, newthetac, newphic);
                    let l = &self.st.lat;
                    let rhonew = natomr / ((l.alatx * l.blaty) * l.clatz);
                    let vol = rhonew / rhoold;
                    let dmaxfac = vol.powf(1.0f32 / 3.0f32);
                    alatfac = dmaxfac;
                    blatfac = dmaxfac;
                    clatfac = dmaxfac;
                    newphib = self.st.lat.phib;
                    newthetac = self.st.lat.thetac;
                    newphic = self.st.lat.phic;
                    change_boxsize(&mut self.st, alatfac, blatfac, clatfac, newphib, newthetac, newphic);
                    self.st.normalise_coordinates();
                    alatnew *= alatfac;
                    blatnew *= blatfac;
                    // NB: the Fortran rescales clatnew with blatfac (equal to clatfac here).
                    clatnew *= blatfac;
                }
            }
            ilatrefine += 1;
        }
        Ok(())
    }

    /// The `if(nbragg.gt.0)` block of the refinement loop: Bragg intensities,
    /// `.inten` file, and replacement of the intermolecular S(Q) up to `qbraggmax`.
    fn bragg_peaks(&mut self, nq: i32) -> Result<()> {
        let nbragg = self.bragg.nbragg;
        let nbraggup = self.bragg.nbraggup;
        let npairs = self.st.mol.npairs;
        let qstep = self.ctl.qstep;
        let nbu = nbraggup as usize;
        let mut qbragg = F1::<f32>::new(nbu);
        let mut intenbragg = F2::<f32>::new(npairs.max(1) as usize, nbu);
        let mut intenconv = F2::<f32>::new(npairs.max(1) as usize, nbragg.max(1) as usize);
        for ib in 1..=nbu {
            qbragg[ib] = qstep * (ib - 1) as f32;
        }
        let pi = PI4ATAN;
        let twopi = 2.0f32 * pi;
        let fname = format!("{}{}.inten", self.st.sym.directory_name, self.st.sym.fname_run);
        let mut out = OutFile::create(std::path::Path::new(&fname));
        self.st.unnormalise_coordinates();
        let (qbraggmax, qbraggwidth, fwhmq) = (self.bragg.qbraggmax, self.bragg.qbraggwidth, self.ctl.fwhmq);
        self.bragg.bragginten(&self.st, qstep, qbraggmax, qbraggwidth, fwhmq, nbraggup, &qbragg, &mut intenbragg)?;
        self.st.normalise_coordinates();
        out.line(&ListWriter::new().chars("# ").int(npairs).int(nbragg).int(nbraggup).finish());
        let multfactor = (self.bragg.multa * self.bragg.multb * self.bragg.multc) as f32;
        for ib in 1..=nbraggup {
            for ip in 1..=npairs {
                intenbragg[(ip, ib)] = intenbragg[(ip, ib)] / multfactor - self.qwind[(1, ib, ip)];
            }
        }
        let debyewaller = self.bragg.debyewaller;
        if debyewaller > 0.0 {
            for ib in 1..=nbraggup {
                let mut qbraggfactor = qbragg[ib] * debyewaller;
                qbraggfactor = 0.5 * qbraggfactor * qbraggfactor;
                if qbraggfactor < 40.0 {
                    qbraggfactor = (-qbraggfactor).exp();
                } else {
                    qbraggfactor = 0.0;
                }
                for ip in 1..=npairs {
                    intenbragg[(ip, ib)] *= qbraggfactor;
                }
            }
        }
        let rminbragg = self.st.inp.rminbragg;
        let p15 = twopi.powf(1.5);
        let gfactor = rminbragg.abs() / (p15 * 2.0);
        let gfactor0 = pow3(rminbragg.abs()) / p15;
        for ib in 1..=nbraggup {
            let sumfactor = 4.0 * pi * qbragg[ib] * qstep;
            for it in 1..=nbragg {
                let mut sum: f32 = 0.0;
                if qbragg[ib] > 0.0 && qbragg[it] > 0.0 {
                    let diff = (qbragg[ib] - qbragg[it]) * rminbragg;
                    let add = (qbragg[ib] + qbragg[it]) * rminbragg;
                    if diff.abs() < 8.0 {
                        sum = (-(0.5 * diff * diff)).exp();
                    }
                    if add.abs() < 8.0 {
                        sum -= (-(0.5 * add * add)).exp();
                    }
                    sum = sum * gfactor / qbragg[it];
                } else {
                    let add = qbragg[ib] * rminbragg;
                    let diff = qbragg[it] * rminbragg;
                    sum = 0.5 * (diff * diff + add * add);
                    if sum < 32.0 {
                        sum = gfactor0 * (-sum).exp();
                    } else {
                        sum = 0.0;
                    }
                    sum *= qbragg[ib];
                }
                sum *= sumfactor;
                for ip in 1..=npairs {
                    intenconv[(ip, it)] += (self.qwind[(2, ib, ip)] - intenbragg[(ip, ib)]) * sum;
                }
            }
        }
        for ip in 1..=npairs {
            for ib in 1..=nbragg {
                self.qwind[(2, ib, ip)] = intenbragg[(ip, ib)] + intenconv[(ip, ib)];
            }
        }
        for ib in 1..=nbragg {
            let mut w = ListWriter::new();
            w.real(qbragg[ib]);
            for ip in 1..=npairs {
                w.real(self.qwind[(2, ib, ip)]).real(intenbragg[(ip, ib)]).real(intenconv[(ip, ib)]);
            }
            out.line(&w.finish());
        }
        out.close()?;
        let mut b = std::mem::take(&mut self.bragg);
        let r = b.sharm_diffuse(&mut self.st, qstep, nq);
        self.bragg = b;
        r
    }
}
