//! Column-format output and input files.
//!
//! Translation of `w_int`, `w_intr`, `r_int` (`auxiliary_routines.f`) and of
//! `sumpe`, `norma`, `initg` (`epsr_standalone_rev3.f`).
//!
//! EPSR writes every distribution as an ASCII table: a header line starting
//! with `#` followed by rows `x  y1 e1  y2 e2 ...` (format `e14.7`), at most
//! 100 columns per file; file `name.EXTnn` holds columns `100*(nn-1)+1...`.
//! The extensions used by `epsr` are:
//!
//! | file | x | contents |
//! |------|---|----------|
//! | `.f01` | Q | simulated intermolecular partial S(Q) and spread |
//! | `.s01` | Q | simulated intramolecular partial S(Q) and spread |
//! | `.d01` | Q | partial S(Q) differences (estimated − simulated) |
//! | `.q01` | Q | estimated partial S(Q) from the data |
//! | `.t01` | Q | data minus background, per dataset |
//! | `.u01` | Q | simulated total F(Q), per dataset |
//! | `.v01` | Q | data − simulation, per dataset |
//! | `.g01` | r | intermolecular partial g(r) and spread |
//! | `.y01` | r | intramolecular partial g(r) and spread |
//! | `.l01` | r | running coordination integral of g(r)−1 |
//! | `.r01` | r | partial g(r) estimated from the data |
//! | `.x01` | r | simulated total G(r), per dataset |
//! | `.w01` | r | Fourier transform of the data, per dataset |
//! | `.o01` | r | reference potential (with and without additions) |
//! | `.p01` | r | empirical potential and intramolecular potential |
//! | `.z01` | r | running coordination numbers |

use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::TextFile;
use crate::fortran::output::OutFile;
use crate::fortran::text::indexlastdot;
use crate::fortran::{F1, F2};
use crate::state::State;

/// Maximum number of columns per output file (`nwrmax` in `wrcom.inc`).
pub const NWRMAX: i32 = 100;

/// Build the file name `fname(1:indexlastdot+1)//nn`, e.g. `run.EPSR.g` → `run.EPSR.g01`.
fn numbered(fname: &str, it: i32) -> String {
    let ilen = (indexlastdot(fname) + 1) as usize;
    let base: String = fname.chars().take(ilen).collect();
    format!("{base}{:02}", it)
}

/// `subroutine w_intr(lpt1,lpt2,fname,nfile,xd,yd,ed,titlex,titled,...)`:
/// write rows `lpt1..=lpt2` of `nfile` columns (value, error) against `xd`.
/// [`w_int`] is the special case `lpt1 = 1`.
#[allow(clippy::too_many_arguments)]
pub fn w_intr(lpt1: i32, lpt2: i32, fname: &str, nfile: i32, xd: &F1<f32>, yd: &F2<f32>, ed: &F2<f32>, titlex: &str, titled: &[String], log: &dyn Fn(&str)) -> Result<()> {
    let nt = (nfile + (NWRMAX - 1)) / NWRMAX;
    let mut icount = 0;
    let lpt = lpt2 - lpt1 + 1;
    for it in 1..=nt {
        let mut iwr = nfile - (it - 1) * NWRMAX;
        if iwr > NWRMAX {
            iwr = NWRMAX;
        }
        let nblock = 2 * iwr + 1;
        let formatout = fmt("(a1,i3.3,a11)", &[FV::S("("), FV::I(nblock as i64), FV::S("(1x,e14.7))")]);
        let fname2 = numbered(fname, it);
        log(&format!(" Writing {:>4} points to {}", lpt, fname2));
        let mut out = OutFile::create(std::path::Path::new(&fname2));
        let mut hdr: Vec<FV> = vec![FV::Ch(titlex, 14)];
        let first = icount;
        for ir in 1..=iwr {
            hdr.push(FV::Ch(&titled[(first + ir - 1) as usize], 29));
        }
        out.line(&fmt(&format!("(1h#,1x,a14,{:03}(1x,a29))", iwr), &hdr));
        for i in lpt1..=lpt2 {
            let mut v: Vec<FV> = Vec::with_capacity((2 * iwr + 1) as usize);
            v.push(FV::R(xd[i]));
            for ir in 1..=iwr {
                v.push(FV::R(yd[(i, first + ir)]));
                v.push(FV::R(ed[(i, first + ir)]));
            }
            out.line(&fmt(&formatout, &v));
        }
        icount += iwr;
        out.close().map_err(|e| EpsrError::Io(format!("{fname2}: {e}")))?;
    }
    Ok(())
}

/// `subroutine w_int(lpt,fname,nfile,xd,yd,ed,titlex,titled,...)`.
#[allow(clippy::too_many_arguments)]
pub fn w_int(lpt: i32, fname: &str, nfile: i32, xd: &F1<f32>, yd: &F2<f32>, ed: &F2<f32>, titlex: &str, titled: &[String], log: &dyn Fn(&str)) -> Result<()> {
    w_intr(1, lpt, fname, nfile, xd, yd, ed, titlex, titled, log)
}

/// `subroutine r_int(lpt,fname,nfile,xd,yd,ed,...)`: read a file written by
/// `w_int` back into `yd`/`ed`.  Returns the number of rows, or `None` if the
/// file does not exist.
pub fn r_int(fname: &str, nfile: i32, xd: &mut F1<f32>, yd: &mut F2<f32>, ed: &mut F2<f32>, log: &dyn Fn(&str)) -> Result<Option<i32>> {
    let limit: f32 = 1.0e-15;
    let nt = (nfile + (NWRMAX - 1)) / NWRMAX;
    let mut icount = 0;
    let mut lpt = 0;
    for it in 1..=nt {
        let mut iwr = nfile - (it - 1) * NWRMAX;
        if iwr > NWRMAX {
            iwr = NWRMAX;
        }
        let fname2 = numbered(fname, it);
        let mut f = match TextFile::open(std::path::Path::new(&fname2)) {
            Ok(f) => f,
            Err(_) => return Ok(None),
        };
        let first = f.read_line()?;
        if !first.contains('#') {
            f.backspace();
        }
        let mut rows: Vec<(f32, Vec<(f32, f32)>)> = Vec::new();
        loop {
            let rec = {
                let mut r = f.list();
                (|| -> crate::fortran::IoRes<(f32, Vec<(f32, f32)>)> {
                    let x = r.real()?.unwrap_or(0.0);
                    let mut v = Vec::with_capacity(iwr as usize);
                    for _ in 0..iwr {
                        let mut y = r.real()?.unwrap_or(0.0);
                        let mut e = r.real()?.unwrap_or(0.0);
                        if y.abs() < limit {
                            y = 0.0;
                        }
                        if e.abs() < limit {
                            e = 0.0;
                        }
                        v.push((y, e));
                    }
                    Ok((x, v))
                })()
            };
            match rec {
                Ok(r) => rows.push(r),
                Err(crate::fortran::IoErr::Eof) => break,
                Err(e) => return Err(EpsrError::Input(format!("{fname2}: {e:?}"))),
            }
        }
        lpt = rows.len() as i32;
        log(&format!(" Read {:>4} points from {}", lpt, fname2));
        for ir in 1..=iwr {
            icount += 1;
            for (i, (_, v)) in rows.iter().enumerate() {
                yd[(i + 1, icount)] = v[(ir - 1) as usize].0;
                ed[(i + 1, icount)] = v[(ir - 1) as usize].1;
            }
        }
        for (i, (x, _)) in rows.iter().enumerate() {
            xd[i + 1] = *x;
        }
    }
    Ok(Some(lpt))
}

impl State {
    /// `subroutine initg(nq,nz,inpfile,ileninpfile,nsum,ireset)`: read the
    /// accumulated g(r) and S(Q) of previous iterations (`.g01 .y01 .f01 .s01`).
    /// `nprev` is the Fortran argument `nsum` (`nsum`, or `nesum` when `nsumt < 0`).
    pub fn initg(&mut self, nprev: i32) -> Result<()> {
        if nprev > 0 {
            let run = self.sym.fnamerun.clone();
            let np = self.mol.npairs;
            let log = self.logger();
            let mut rarray = F1::new(crate::state::MPZ);
            let mut qarray = F1::new(crate::state::MPTS);
            let d = &mut self.dist;
            let r1 = r_int(&format!("{run}.g"), np, &mut rarray, &mut d.gorinter, &mut d.gr2inter, &log)?;
            let r2 = r_int(&format!("{run}.y"), np, &mut rarray, &mut d.gorintra, &mut d.gr2intra, &log)?;
            let r3 = r_int(&format!("{run}.f"), np, &mut qarray, &mut d.soqinter, &mut d.sq2inter, &log)?;
            let r4 = r_int(&format!("{run}.s"), np, &mut qarray, &mut d.soqintra, &mut d.sq2intra, &log)?;
            if r1.is_none() || r2.is_none() || r3.is_none() || r4.is_none() {
                return Err(EpsrError::Input("nsumt > 0 but the accumulated .g01/.y01/.f01/.s01 files are missing".into()));
            }
        } else {
            let d = &mut self.dist;
            d.gorinter.fill(0.0);
            d.gr2inter.fill(0.0);
            d.gorintra.fill(0.0);
            d.gr2intra.fill(0.0);
            d.soqinter.fill(0.0);
            d.sq2inter.fill(0.0);
            d.soqintra.fill(0.0);
            d.sq2intra.fill(0.0);
        }
        Ok(())
    }

    /// `subroutine sumpe(nz,z,inpfile,ileninpfile)`: energy contributions of
    /// each pair from the reference and empirical potentials → `.terg`.
    pub fn sumpe(&mut self, nz: i32, z: &F1<f32>) -> Result<()> {
        let m = &self.mol;
        let d = &mut self.dist;
        let p = &self.pot;
        let mut toteref: f32 = 0.0;
        let mut toteemp: f32 = 0.0;
        let rmaxpot = 0.5 * m.dmax;
        for id in 1..=m.npairs {
            d.energyref[id] = 0.0;
            d.energyemp[id] = 0.0;
            d.energytot[id] = 0.0;
            let it = m.n1extr[id];
            let jt = m.n2extr[id];
            let mut cprod = m.cfrac[it] * m.cfrac[jt];
            if it != jt {
                cprod *= 2.0;
            }
            let mut j = 2;
            while j <= nz && z[j] <= rmaxpot {
                let wt = cprod * d.ycoh[j] * d.gorinter[(j, id)];
                d.energyref[id] += wt * p.ecoul[(j, id)];
                d.energyemp[id] += wt * p.xavs[(j, id)];
                j += 1;
            }
            d.energytot[id] = d.energyref[id] + d.energytot[id];
            toteref += d.energyref[id];
            toteemp += d.energyemp[id];
        }
        let totetot = toteref + toteemp;
        let fname = format!("{}.terg", self.sym.fnamerun);
        let mut out = OutFile::create(std::path::Path::new(&fname));
        out.line(&fmt("(3(1x,e13.6))", &[FV::R(toteref), FV::R(toteemp), FV::R(totetot)]));
        for id in 1..=m.npairs {
            out.line(&fmt("(2(1x,a3))", &[FV::Ch(&m.type_[m.n1extr[id]], 3), FV::Ch(&m.type_[m.n2extr[id]], 3)]));
            out.line(&fmt("(3(1x,e13.6))", &[FV::R(d.energyref[id]), FV::R(d.energyemp[id]), FV::R(d.energytot[id])]));
        }
        out.close()?;
        Ok(())
    }

    /// `subroutine norma(ireset,itimes,ntimes)`: write the accumulator file
    /// `.inpa` and the uniform-density file `.uni`, and store the EP
    /// parameters back into the symbol table (for the `.pcof` file).
    pub fn norma(&mut self) -> Result<()> {
        let fname = format!("{}.inpa", self.sym.fnamerun);
        let mut out = OutFile::create(std::path::Path::new(&fname));
        let i = &self.inp;
        let mut w = ListWriter::new();
        w.int(i.nisum).int(i.nesum).int(i.nesumt).real(i.ensum).real(i.ensum2).real(i.prsum).real(i.prsum2);
        out.line(&w.finish());
        let p = &self.pot;
        let nfile = self.dist.nfile;
        for id in 1..=nfile {
            let mut w = ListWriter::new();
            w.int(id);
            out.line(&w.finish());
            let mut w = ListWriter::new();
            w.int(p.npcoeff[id]).real(p.pstepd[id]).real(self.dist.sigma2[id]);
            out.line(&w.finish());
            let mut w = ListWriter::new();
            for ic in 1..=p.npcoeff[id] {
                w.real(p.pcoefd[(ic, id)]);
            }
            out.line(&w.finish());
        }
        let nnew = nfile + self.mol.npairs;
        let mut w = ListWriter::new();
        w.real(p.feedback);
        for j in 1..=nfile {
            w.real(p.efacm[j]);
        }
        w.int(self.wts.nqinv);
        out.line(&w.finish());
        for iq in 1..=self.wts.nqinv {
            for ip in 1..=self.mol.npairs {
                let mut w = ListWriter::new();
                w.int(self.wts.qinv[iq]).int(self.mol.n1extr[ip]).int(self.mol.n2extr[ip]);
                out.line(&w.finish());
                let mut w = ListWriter::new();
                for id in 1..=nnew {
                    w.dble(self.wts.cwttq[(iq, ip, id)]);
                }
                out.line(&w.finish());
            }
        }
        out.close()?;
        // EP parameters back to the symbol table
        let ld_r = |v: f32| {
            let mut w = ListWriter::new();
            w.real(v);
            w.finish()
        };
        let ld_i = |v: i32| {
            let mut w = ListWriter::new();
            w.int(v);
            w.finish()
        };
        let mut iref = self.sym.ireff(self.sym.iqminref, 1);
        let i = &self.inp;
        let s = &mut self.sym;
        iref += 1;
        s.symbol_value[iref] = ld_r(i.roverlap);
        iref += 1;
        s.symbol_value[iref] = ld_r(i.rminfac);
        iref += 1;
        s.symbol_value[iref] = ld_r(i.rminpt);
        iref += 1;
        s.symbol_value[iref] = ld_r(i.rmaxpt);
        iref += 1;
        s.symbol_value[iref] = ld_r(i.rbroad);
        iref += 1;
        s.symbol_value[iref] = fmt("(f7.3,2x,f7.3)", &[FV::R(i.expecf), FV::R(i.expecfpwr)]);
        iref += 3;
        s.symbol_value[iref] = ld_r(i.rcharge);
        iref += 1;
        s.symbol_value[iref] = ld_r(self.mol.power);
        iref += 1;
        s.symbol_value[iref] = ld_i(p.ncoeffp);
        iref += 1;
        s.symbol_value[iref] = ld_r(p.pdmax);
        iref += 1;
        s.symbol_value[iref] = ld_r(p.pdstep);
        iref += 1;
        s.symbol_value[iref] = ld_i(p.npitss);
        iref += 1;
        s.symbol_value[iref] = ld_r(p.paccept);
        iref += 2;
        s.symbol_value[iref] = fmt("(f8.5)", &[FV::R(p.psigma2read)]);
        iref += 1;
        s.symbol_value[iref] = fmt("(f8.5)", &[FV::R(p.refpotfac)]);
        // .uni
        let fname = format!("{}.uni", self.sym.fnamerun);
        let mut out = OutFile::create(std::path::Path::new(&fname));
        let d = &self.dist;
        let mut v: Vec<FV> = Vec::new();
        for j in 1..=d.npt1 {
            v.push(FV::R(d.x[j]));
            v.push(FV::R(d.ycoh[j]));
        }
        out.line(&fmt("(2(1x,e13.6))", &v));
        out.close()?;
        Ok(())
    }
}
