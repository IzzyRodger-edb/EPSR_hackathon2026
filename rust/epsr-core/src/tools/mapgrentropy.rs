//! `mapgrentropy` — translation of `fortran_src/mapgrentropy.f` (+ `plotsharm_routines.f`).
//!
//! Estimates the radial, spatial and total (orientational) two-body excess
//! entropy from the spherical-harmonic coefficients written by `sharm`/`sdf`,
//! by Monte Carlo sampling of random relative orientations.
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> mapgrentropy <input file>`
//! Keyword table: `crate::symtab_tools::MAPGRENTROPY`.
//! Outputs: `<run>.mapgrentropy.out`, `<run>.mapgrentropy.s01`,
//! `negativegr.txt`, and `ran1.dat` (random number generator state, read at
//! start if present).

use super::sharm::sharm_routines::Coeffs;
use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::{read_int, read_real, TextFile};
use crate::fortran::output::OutFile;
use crate::fortran::text::{fstr, indexlastdot};
use crate::fortran::{F1, F2};
use crate::output::w_int;
use crate::rng::{Ran1, NTAB};
use crate::state::State;

/// `initran1(10)` on `ran1.dat` (opened with status='unknown').
fn initran1(path: &std::path::Path) -> Ran1 {
    let vals = (|| -> Option<(i32, i32, i32, [i32; NTAB])> {
        let mut f = TextFile::open(path).ok()?;
        let mut r = f.list();
        let mcall = r.int().ok()??;
        let idum = r.int().ok()??;
        let iy = r.int().ok()??;
        let mut iv = [0i32; NTAB];
        for v in iv.iter_mut() {
            *v = r.int().ok()??;
        }
        Some((mcall, idum, iy, iv))
    })();
    // On a read error the Fortran resets only idum, iy, mcall (iv keeps its
    // previous - here zero - contents); a partially read record keeps the
    // values read so far, which cannot be reproduced exactly here.
    if vals.is_none() && path.exists() && std::fs::metadata(path).map(|m| m.len() > 0).unwrap_or(false) {
        // a non-empty but unreadable ran1.dat
        return Ran1::from_record(None);
    }
    Ran1::from_record(vals)
}

/// Run the `mapgrentropy` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    if !super::sharm::open(ctx, &mut st, &crate::symtab_tools::MAPGRENTROPY)? {
        return Ok(());
    }
    let log = st.logger();
    let pi = std::f32::consts::PI;
    let ranpath = ctx.path("ran1.dat");
    // open(10,file='ran1.dat',status='unknown') creates the file if absent
    if !ranpath.exists() {
        OutFile::create(&ranpath).close()?;
    }
    let mut rng = initran1(&ranpath);
    let fnameinp = st.sym.fnameinp.clone();
    let idot = indexlastdot(&fnameinp);
    let fnameout = if idot > 0 { format!("{}out", &fnameinp[..idot as usize]) } else { format!("{}.out", fnameinp.trim_end()) };
    let c = Coeffs::read(&st)?;
    let nr = c.nr;
    let nexpan = c.nexpan;
    let bad = |what: &str| EpsrError::Input(format!("mapgrentropy: cannot read {what}"));
    let sv = |i: usize| st.sym.symbol_value[i].clone();
    let forcemirrorsym = fstr(&sv(5), 2) == "y";
    let nvary = read_int(&sv(6)).map_err(|_| bad("ntrials"))?;
    let nfreq = read_int(&sv(7)).map_err(|_| bad("nfreq"))?;
    let molndens = read_real(&sv(8)).map_err(|_| bad("molndens"))?;
    let molmass = read_real(&sv(9)).map_err(|_| bad("molmass"))?;
    let moltemp = read_real(&sv(10)).map_err(|_| bad("moltemp"))?;
    let molrot = read_real(&sv(11)).map_err(|_| bad("molrot"))?;
    let molvib = read_real(&sv(12)).map_err(|_| bad("molvib"))?;
    let gasconstant: f32 = 8.314472;
    let mollambda: f32 = 66.2606876f32 / (2.0f32 * pi * molmass * 1.660539f32 * 1.3806503f32 * moltemp).sqrt();
    let transterm: f32 = gasconstant * (2.5f32 - (molndens * (mollambda * (mollambda * mollambda))).ln());
    let rotterm: f32 = gasconstant * (1.5f32 + ((pi * (moltemp * moltemp * moltemp) / molrot).sqrt()).ln());
    let vibterm: f32 = molvib;
    let kineticterm: f32 = transterm + rotterm + vibterm;
    let kt = kineticterm as f64;
    let mut entropyrd: f64 = kt;
    let mut entropysp: f64 = kt;
    let mut entropyto: f64 = kt;
    let mut grw8: F2<f64> = F2::new(500, 3);
    let mut grwerr8: F2<f64> = F2::new(500, 3);
    for it in 1..=3 {
        for i in 1..=nr {
            grw8[(i, it)] = kt;
            grwerr8[(i, it)] = kt;
        }
    }
    let mut spnegative: F1<f32> = F1::new(500);
    let mut tonegative: F1<f32> = F1::new(500);
    let rd = &c.rd;
    let fact: f32 = 2.0f32 * pi * gasconstant * molndens * (rd[3] - rd[2]);
    let mut out = OutFile::create(std::path::Path::new(&fnameout));
    out.line(&fmt("(\"Temperature [K]: \",e12.5)", &[FV::R(moltemp)]));
    out.line(&fmt("(\"Number density [mols/A**3]: \",e12.5)", &[FV::R(molndens)]));
    out.line(&fmt("(\"De Broglie wavelength: \",e12.5)", &[FV::R(mollambda)]));
    out.line(&fmt("(\"Perfect gas entropy, trans. term: \",e12.5)", &[FV::R(transterm)]));
    out.line(&fmt("(\"Perfect gas entropy, rot. term: \",e12.5)", &[FV::R(rotterm)]));
    out.line(&fmt("(\"Perfect gas entropy,vib. term: \",e12.5)", &[FV::R(vibterm)]));
    out.line(&fmt("(\"Perfect gas entropy, total: \",e12.5)", &[FV::R(kineticterm)]));
    let (mut minrd, mut minsp, mut minto, mut maxrd, mut maxsp, mut maxto) = (0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32);
    // Locals the Fortran leaves uninitialised unless a minimum/maximum is
    // found; the reference build writes 0 for them in that case.
    let (mut minrdr, mut minspr, mut mintor, mut maxrdr, mut maxspr, mut maxtor) = (0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32);
    let mut spmin = [0.0f32; 5];
    let mut tomin = [0.0f32; 5];
    let mut spmax = [0.0f32; 5];
    let mut tomax = [0.0f32; 5];
    let only = [1, 1, 1, 1, 1, 1];
    let ne = nexpan.max(0) as usize;
    let mut zc = vec![0.0f32; ne];
    let mut zs = vec![0.0f32; ne];
    let mut zspc = vec![0.0f32; ne];
    let mut zsps = vec![0.0f32; ne];
    let mut temprd: F1<f32> = F1::new(500);
    let mut tempsp: F1<f32> = F1::new(500);
    let mut tempto: F1<f32> = F1::new(500);
    let mut it = 0;
    while it < nvary {
        it += 1;
        let ctl = 2.0f32 * rng.ran1() - 1.0f32;
        let phil = rng.ran1() * 360.0f32;
        let phi = rng.ran1() * 360.0f32;
        let cz = 2.0f32 * rng.ran1() - 1.0f32;
        let chi = rng.ran1() * 360.0f32;
        for k in 1..=ne {
            let (sumc, sums, spsumc, spsums) = c.angular(k, only, phil, ctl, phi, cz, chi);
            zc[k - 1] = sumc;
            zspc[k - 1] = spsumc;
            if forcemirrorsym {
                zs[k - 1] = 0.0;
                zsps[k - 1] = 0.0;
            } else {
                zs[k - 1] = sums;
                zsps[k - 1] = spsums;
            }
        }
        for i in 1..=nr {
            let mut tsumrd: f32 = 0.0;
            let mut tsumsp: f32 = 0.0;
            let mut tsumto: f32 = 0.0;
            for k in 1..=nexpan {
                let ku = (k - 1) as usize;
                let term = c.datac[(i, k)] * zc[ku] + c.datas[(i, k)] * zs[ku];
                tsumto += term;
                if c.l2ex[ku] == 0 {
                    let term = c.datac[(i, k)] * zspc[ku] + c.datas[(i, k)] * zsps[ku];
                    tsumsp += term;
                    if c.l1ex[ku] == 0 {
                        tsumrd += term;
                    }
                }
            }
            temprd[i] = tsumrd;
            tempsp[i] = tsumsp;
            tempto[i] = tsumto;
        }
        let mut sumrd: f64 = kt;
        let mut sumsp: f64 = kt;
        let mut sumto: f64 = kt;
        let deg = |x: f32| x.acos() * 180.0f32 / pi;
        for i in 1..=nr {
            let factr = fact * rd[i] * rd[i];
            // radial
            let temp = temprd[i];
            let term = 1.0f32 + temp;
            if term > maxrd {
                maxrd = term;
                maxrdr = rd[i];
            }
            if term < minrd {
                minrd = term;
                minrdr = rd[i];
            }
            let contribution: f32 = if term > 0.0 { factr * (temp - term * term.ln()) } else { -factr };
            sumrd += contribution as f64;
            grw8[(i, 1)] += sumrd;
            // spatial
            let mut temp = tempsp[i];
            let mut term = 1.0f32 + temp;
            if term > maxsp {
                maxsp = term;
                maxspr = rd[i];
                spmax = [deg(ctl), phil, phi, deg(cz), chi];
            }
            if term < minsp {
                minsp = term;
                minspr = rd[i];
                spmin = [deg(ctl), phil, phi, deg(cz), chi];
            }
            if temp > 0.0 && spnegative[i] > 0.0 {
                let wtnegative = temp * spnegative[i] / term;
                temp += wtnegative;
                spnegative[i] -= wtnegative;
                term = 1.0f32 + temp;
            }
            let contribution: f32 = if term > 0.0 { factr * (temp - term * term.ln()) } else { -factr };
            sumsp += contribution as f64;
            grw8[(i, 2)] += sumsp;
            // total
            let mut temp = tempto[i];
            let mut term = 1.0f32 + temp;
            if term > maxto {
                maxto = term;
                maxtor = rd[i];
                tomax = [deg(ctl), phil, phi, deg(cz), chi];
            }
            if term < minto {
                minto = term;
                mintor = rd[i];
                tomin = [deg(ctl), phil, phi, deg(cz), chi];
            }
            if temp > 0.0 && tonegative[i] > 0.0 {
                let wtnegative = temp * tonegative[i] / term;
                temp += wtnegative;
                tonegative[i] -= wtnegative;
                term = 1.0f32 + temp;
            }
            let contribution: f32 = if term > 0.0 { factr * (temp - term * term.ln()) } else { -factr };
            sumto += contribution as f64;
            grw8[(i, 3)] += sumto;
        }
        entropyrd += sumrd;
        entropysp += sumsp;
        entropyto += sumto;
        if it == nfreq * (it / nfreq) || it == nvary {
            let r = it as f32 as f64;
            out.line(&fmt(
                "(i10,6(1x,e13.6))",
                &[FV::I(it as i64), FV::D(sumrd), FV::D(sumsp), FV::D(sumto), FV::D(entropyrd / r), FV::D(entropysp / r), FV::D(entropyto / r)],
            ));
        }
    }
    out.line("");
    out.line(&fmt("(\"Minimum values, radial, spatial total: \",3(1x,e12.5))", &[FV::R(minrd), FV::R(minsp), FV::R(minto)]));
    out.line(&fmt("(\"Radius values of minima, radial, spatial total: \",3(1x,e12.5))", &[FV::R(minrdr), FV::R(minspr), FV::R(mintor)]));
    out.line("");
    out.line(&fmt("(\"thetal, phil, phim, thetam, chim of spatial minimum: \",5(1x,e12.5))", &spmin.map(FV::R)));
    out.line(&fmt("(\"thetal, phil, phim, thetam, chim of total minimum: \",5(1x,e12.5))", &tomin.map(FV::R)));
    out.line("");
    out.line(&fmt("(\"Maximum values, radial, spatial total: \",3(1x,e12.5))", &[FV::R(maxrd), FV::R(maxsp), FV::R(maxto)]));
    out.line(&fmt("(\"Radius values of maxima, radial, spatial total: \",3(1x,e12.5))", &[FV::R(maxrdr), FV::R(maxspr), FV::R(maxtor)]));
    out.line("");
    out.line(&fmt("(\"thetal, phil, phim, thetam, chim of spatial maximum: \",5(1x,e12.5))", &spmax.map(FV::R)));
    out.line(&fmt("(\"thetal, phil, phim, thetam, chim of total maximum: \",5(1x,e12.5))", &tomax.map(FV::R)));
    out.close()?;
    let mut ng = OutFile::create(&ctx.path("negativegr.txt"));
    for i in 1..=nr {
        ng.line(&fmt("(3(1x,e12.5))", &[FV::R(rd[i]), FV::R(spnegative[i]), FV::R(tonegative[i])]));
    }
    ng.close()?;
    let mut grw: F2<f32> = F2::new(500, 5);
    let mut grwerr: F2<f32> = F2::new(500, 5);
    let rn = nvary as f32 as f64;
    for i in 1..=nr {
        for it in 1..=3 {
            grw[(i, it)] = (grw8[(i, it)] / rn) as f32;
            grwerr[(i, it)] = (grwerr8[(i, it)] / rn) as f32;
        }
        grw[(i, 4)] = grw[(i, 2)] - grw[(i, 1)];
        grwerr[(i, 4)] = 0.0;
        grw[(i, 5)] = grw[(i, 3)] - grw[(i, 2)];
        grwerr[(i, 5)] = 0.0;
    }
    let titledat = vec![
        fmt("(10x,a,11x)", &[FV::S("S_radial")]),
        fmt("(10x,a,10x)", &[FV::S("S_spatial")]),
        fmt("(11x,a,12x)", &[FV::S("S_total")]),
        fmt("(7x,a,7x)", &[FV::S("S_extra_spatial")]),
        fmt("(7x,a,8x)", &[FV::S("S_extra_orient")]),
    ];
    let ilen = indexlastdot(&fnameinp) as usize;
    let fname = format!("{}s", &fnameinp[..ilen]);
    w_int(nr, &fname, 5, rd, &grw, &grwerr, "r[A]", &titledat, &log)?;
    let mut w = ListWriter::new();
    w.int(rng.mcall).int(rng.idum).int(rng.iy);
    for v in rng.iv.iter() {
        w.int(*v);
    }
    let mut rf = OutFile::create(&ranpath);
    rf.line(&w.finish());
    rf.close()?;
    Ok(())
}

/// Program-specific branch of `check_symbol` (`process_symbols.f`, `symbol_setup.eq.'mapgrentropy'`).
/// Returns the Fortran `ierr` (0 = accepted).
pub fn check_symbol_extra(st: &mut State, icount: i32, _ifilecount: i32, iref: i32) -> Result<i32> {
    st.check_symbol_mapgr(icount, iref)
}
