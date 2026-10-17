//! `epsrwts` — translation of `fortran_src/epsrwts.f` (subroutine `epsrwts`).
//!
//! STATUS: translated, verified byte-identical (EPSR-test-ex24).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> epsrwts <file>.ato < answers`
//!
//! The interactive predecessor of `nwts`: the isotopic compositions, the
//! substituted components, the mixing fraction, the output base name and the
//! normalisation flags are read from standard input.  Unlike `nwts` the
//! intramolecular H-H weights of non-exchanging hydrogens are summed isotope
//! by isotope.  Writes `<base>tot.wts` (or `tot1/2/3.wts`), `<base>dif.wts`,
//! `<base>xx/xy/yy.wts` and rewrites the `.ato` file.

use super::nwts::{scatt_len, MISOTOPE};
use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::IoErr;
use crate::fortran::output::OutFile;
use crate::fortran::text::{fstr, parse};
use crate::fortran::{F2, F3};
use crate::state::{MCOMP, MPAIRS};

fn eof() -> EpsrError {
    EpsrError::Input("epsrwts: unexpected end of standard input".into())
}

/// `read(5,*) i` (no error handling: an error stops the Fortran program).
fn rd_int(ctx: &mut ToolContext) -> Result<Option<i32>> {
    let mut r = ctx.stdin.list();
    Ok(r.int()?)
}
fn rd_real(ctx: &mut ToolContext) -> Result<Option<f32>> {
    let mut r = ctx.stdin.list();
    Ok(r.real()?)
}

/// `read(5,*,err=label) x`: `Ok(None)` = conversion error (jump to label),
/// end of file is fatal.
fn rd_int_err(ctx: &mut ToolContext, v: &mut i32) -> Result<bool> {
    let mut r = ctx.stdin.list();
    match r.int() {
        Ok(Some(x)) => {
            *v = x;
            Ok(true)
        }
        Ok(None) => Ok(true),
        Err(IoErr::Eof) => Err(eof()),
        Err(IoErr::Bad(_)) => Ok(false),
    }
}
fn rd_real_err(ctx: &mut ToolContext, v: &mut f32) -> Result<bool> {
    let mut r = ctx.stdin.list();
    match r.real() {
        Ok(Some(x)) => {
            *v = x;
            Ok(true)
        }
        Ok(None) => Ok(true),
        Err(IoErr::Eof) => Err(eof()),
        Err(IoErr::Bad(_)) => Ok(false),
    }
}

/// Run the `epsrwts` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    let readfile0 = match ctx.args.first() {
        Some(a) => a.clone(),
        None => return crate::error::unsupported("epsrwts without a .ato file argument (interactive file search)"),
    };
    // strip any extension (first '.')
    let readfile = match readfile0.find('.') {
        Some(p) if p > 0 => readfile0[..p].to_string(),
        _ => readfile0.clone(),
    };
    let fnameato = format!("{}.ato", st.sym.full(&readfile).to_string_lossy());
    st.sym.fnameato = fnameato.clone();
    st.mol.ntype = 0;
    st.mol.nmoltype = 0;
    st.mol.sizefactor = 1.0;
    let ok = st.minitc(&fnameato);
    if ok.is_err() || st.dist.rho == 0.0 {
        if let Err(EpsrError::Unsupported(s)) = ok {
            return Err(EpsrError::Unsupported(s));
        }
        st.log(&format!("\nepsrwts> Problems with .ato file {fnameato}\nRoutine will exit now\n"));
        return Ok(());
    }
    st.log("\nepsrwts> Program to calculate inter- and intra-molecular weightings \n for DCS, 1st- or 2nd-order difference data.\n");
    st.log("\nepsrwts> Is the output to be per atom (1), or per molecule (2)? ");
    let nrho = rd_int(ctx)?.ok_or_else(eof)?;
    let atmolrat: f32 = st.mol.natomr as f32 / st.mol.nmol as f32;
    let ntype = st.mol.ntype;
    let nt = ntype as usize;
    if (1..=nt).any(|it| st.mol.qtypea[it] != 0) {
        // several arrays of Q atoms are used uninitialised by the Fortran
        return crate::error::unsupported("epsrwts: Q-type atoms");
    }
    let mut norder;
    loop {
        st.log("\nepsrwts> How many samples (1,2, or 3)? (0 to quit) ");
        norder = rd_int(ctx)?.ok_or_else(eof)?;
        if norder < 1 {
            return Ok(());
        }
        if norder <= 3 {
            break;
        }
    }
    let mut iflagsubs = vec![0i32; nt + 1];
    // Fortran locals keep their value when a list-directed read gets a null value
    let mut massno: i32 = 0;
    let mut isubs: i32 = 0;
    let mut nisotope: F2<i32> = F2::new(MCOMP, 2);
    let mut batom: F3<f32> = F3::new(MCOMP, MISOTOPE, 2);
    let mut abund: F3<f32> = F3::new(MCOMP, MISOTOPE, 2);
    let too_many = || EpsrError::Input(format!("epsrwts: more than {MISOTOPE} isotopes for one component"));
    for it in 1..=nt {
        // 206
        let (chsymtest, btest) = loop {
            let chsymtest = fstr(&st.mol.atsymbol[it], 3);
            let (btest, ifind) = scatt_len(&chsymtest, 0);
            if ifind != 0 {
                break (chsymtest, btest);
            }
            st.log(&format!("\nepsrwts> Did not recognise chemical symbol {chsymtest:<3} for atom type {:<3}\n\nepsrwts> Try typing it again: ", st.mol.type_[it]));
            let mut r = ctx.stdin.list();
            match r.word() {
                Ok(Some(w)) => st.mol.atsymbol[it] = fstr(&w, 3),
                Ok(None) => {}
                Err(IoErr::Eof) => return Err(eof()),
                Err(_) => {}
            }
        };
        if st.mol.atsymbol[it] == "H" {
            st.mol.iexchange[it] = -1;
            while st.mol.iexchange[it] != 0 && st.mol.iexchange[it] != 1 {
                let mut r = ctx.stdin.list();
                match r.int() {
                    Ok(Some(x)) => st.mol.iexchange[it] = x,
                    Ok(None) => {}
                    Err(IoErr::Eof) => return Err(eof()),
                    Err(_) => {}
                }
            }
            if st.mol.iexchange[it] != 1 {
                st.mol.iexchange[it] = 0;
            }
        } else {
            st.mol.iexchange[it] = 0;
        }
        let mut asum: f32 = 0.0;
        let mut nsum: usize = 1;
        loop {
            // 205
            if nsum > MISOTOPE {
                return Err(too_many());
            }
            if !rd_int_err(ctx, &mut massno)? {
                continue;
            }
            let mut a = abund[(it, nsum, 1)];
            let okr = rd_real_err(ctx, &mut a)?;
            abund[(it, nsum, 1)] = a;
            if !okr {
                continue;
            }
            if massno == 0 {
                batom[(it, nsum, 1)] = btest;
            } else {
                let (b, ifind) = scatt_len(&chsymtest, massno);
                batom[(it, nsum, 1)] = b;
                if ifind == 0 {
                    continue;
                }
            }
            asum += abund[(it, nsum, 1)];
            if asum < 1.0 {
                nsum += 1;
                continue;
            }
            for _isum in 1..=nsum {
                abund[(it, nsum, 1)] /= asum;
            }
            break;
        }
        nisotope[(it, 1)] = nsum as i32;
    }
    // exchangeable atoms, sample 1
    let average_exchange = |k: usize, nisotope: &mut F2<i32>, batom: &mut F3<f32>, abund: &mut F3<f32>, st: &crate::state::State| {
        let mut bscatsum: f32 = 0.0;
        let mut cfracsum: f32 = 0.0;
        for it in 1..=nt {
            if st.mol.iexchange[it] == 1 {
                cfracsum += st.mol.cfrac[it];
                for ii in 1..=nisotope[(it, k)] as usize {
                    let term = st.mol.cfrac[it] * abund[(it, ii, k)] * batom[(it, ii, k)];
                    bscatsum += term;
                }
            }
        }
        if cfracsum > 0.0 {
            bscatsum /= cfracsum;
        } else {
            bscatsum = 0.0;
        }
        for it in 1..=nt {
            if st.mol.iexchange[it] == 1 {
                nisotope[(it, k)] = 1;
                abund[(it, 1, k)] = 1.0;
                batom[(it, 1, k)] = bscatsum;
            }
        }
    };
    average_exchange(1, &mut nisotope, &mut batom, &mut abund, &st);
    if norder > 1 {
        loop {
            // 204
            if !rd_int_err(ctx, &mut isubs)? {
                continue;
            }
            if isubs >= 1 && isubs <= ntype {
                let is = isubs as usize;
                iflagsubs[is] = 1;
                let mut asum: f32 = 0.0;
                let mut nsum: usize = 1;
                loop {
                    // 208
                    if nsum > MISOTOPE {
                        return Err(too_many());
                    }
                    if !rd_int_err(ctx, &mut massno)? {
                        continue;
                    }
                    let mut a = abund[(is, nsum, 2)];
                    let okr = rd_real_err(ctx, &mut a)?;
                    abund[(is, nsum, 2)] = a;
                    if !okr {
                        continue;
                    }
                    let chsymtest = fstr(&st.mol.atsymbol[is], 3);
                    let (b, ifind) = scatt_len(&chsymtest, massno);
                    batom[(is, nsum, 2)] = b;
                    if ifind == 0 {
                        continue;
                    }
                    asum += abund[(is, nsum, 2)];
                    if asum < 1.0 {
                        nsum += 1;
                        continue;
                    }
                    for _isum in 1..=nsum {
                        abund[(is, nsum, 2)] /= asum;
                    }
                    break;
                }
                nisotope[(is, 2)] = nsum as i32;
                continue;
            }
            break;
        }
        for it in 1..=nt {
            if iflagsubs[it] == 0 {
                nisotope[(it, 2)] = nisotope[(it, 1)];
                for ib in 1..=nisotope[(it, 2)] as usize {
                    batom[(it, ib, 2)] = batom[(it, ib, 1)];
                    abund[(it, ib, 2)] = abund[(it, ib, 1)];
                }
            }
        }
        average_exchange(2, &mut nisotope, &mut batom, &mut abund, &st);
    }
    let (frac1, frac2): (f32, f32) = if norder == 3 {
        let fracmix = rd_real(ctx)?.ok_or_else(eof)?;
        (fracmix, 1.0 - fracmix)
    } else {
        (1.0, 0.0)
    };
    let npart = nt * (nt + 1) / 2;
    let ndata = if norder > 2 { 2 } else { norder } as usize;
    let mut binter: F2<f32> = F2::new(MCOMP, 3);
    let mut wintra: F2<f32> = F2::new(MPAIRS, 3);
    let mut winter: F2<f32> = F2::new(MPAIRS, 3);
    let mut cprod = vec![0.0f32; MPAIRS + 1];
    let mut partlabel = vec![String::new(); MPAIRS + 1];
    let mut wtsum = [0.0f32; 4];
    let cfrac = st.mol.cfrac.clone();
    let iexchange = st.mol.iexchange.clone();
    for idata in 1..=ndata {
        for it1 in 1..=nt {
            binter[(it1, idata)] = 0.0;
            for ib in 1..=nisotope[(it1, idata)] as usize {
                let term1 = abund[(it1, ib, idata)] * batom[(it1, ib, idata)];
                binter[(it1, idata)] += term1;
            }
        }
        let mut ipart = 0usize;
        wtsum[idata] = 0.0;
        for it1 in 1..=nt {
            let binter1 = binter[(it1, idata)];
            for it2 in it1..=nt {
                ipart += 1;
                if idata == 1 {
                    partlabel[ipart] = fmt("(1x,a3,1x,a3)", &[FV::Ch(&st.mol.type_[it1], 3), FV::Ch(&st.mol.type_[it2], 3)]);
                    cprod[ipart] = if it1 == it2 { cfrac[it1] * cfrac[it1] } else { 2.0 * cfrac[it1] * cfrac[it2] };
                }
                let binter2 = binter[(it2, idata)];
                winter[(ipart, idata)] = binter1 * binter2;
                if iexchange[it1] == 0
                    && iexchange[it2] == 0
                    && nisotope[(it1, idata)] == nisotope[(it2, idata)]
                    && st.mol.atsymbol[it1] == "H"
                    && st.mol.atsymbol[it2] == "H"
                {
                    let mut bintra: f32 = 0.0;
                    for ib in 1..=nisotope[(it1, idata)] as usize {
                        let term = abund[(it1, ib, idata)] * batom[(it1, ib, idata)] * batom[(it2, ib, idata)];
                        bintra += term;
                    }
                    wintra[(ipart, idata)] = bintra;
                } else {
                    wintra[(ipart, idata)] = winter[(ipart, idata)];
                }
                let term = cprod[ipart] * winter[(ipart, idata)];
                wtsum[idata] += term;
            }
        }
    }
    if norder == 3 {
        let mut ipart = 0usize;
        wtsum[3] = 0.0;
        for it1 in 1..=nt {
            binter[(it1, 3)] = frac1 * binter[(it1, 1)] + frac2 * binter[(it1, 2)];
        }
        for it1 in 1..=nt {
            for it2 in it1..=nt {
                ipart += 1;
                winter[(ipart, 3)] = binter[(it1, 3)] * binter[(it2, 3)];
                if iexchange[it1] == 0 && iexchange[it2] == 0 {
                    wintra[(ipart, 3)] = frac1 * wintra[(ipart, 1)] + frac2 * wintra[(ipart, 2)];
                } else {
                    wintra[(ipart, 3)] = winter[(ipart, 3)];
                }
                let term = cprod[ipart] * winter[(ipart, 3)];
                wtsum[3] += term;
            }
        }
    }
    // ---- output base name
    let base = loop {
        st.log("\nepsrwts> Type basename of file to output weights to: ");
        let sometext = ctx.stdin.read_line().map_err(|_| eof())?;
        let w = parse(&sometext, 5);
        if w.nwords() > 0 {
            break format!("{}{}", st.sym.directory_name, w.word(&sometext, 1));
        }
    };
    let wline = |out: &mut OutFile, a: f32, b: f32, c: f32| {
        let mut w = ListWriter::new();
        w.real(a).real(b).real(c);
        out.line(&w.finish());
    };
    let one = || {
        let mut w = ListWriter::new();
        w.real(1.0);
        w.finish()
    };
    let n1extr = st.mol.n1extr.clone();
    let n2extr = st.mol.n2extr.clone();
    for idata in 1..=norder as usize {
        let nstruct = rd_int(ctx)?.ok_or_else(eof)?;
        let fnameb = if norder > 1 { fmt("('tot',i1.1,'.wts')", &[FV::I(idata as i64)]) } else { "tot.wts".to_string() };
        let fname = format!("{base}{fnameb}");
        st.log(&format!("\nepsrwts> Writing TOTAL weights to file {fname}"));
        let mut out = OutFile::create(std::path::Path::new(&fname));
        let rhofac: f32 = if nrho == 2 {
            if nstruct == 1 {
                1.0 / atmolrat
            } else {
                atmolrat
            }
        } else {
            1.0
        };
        out.line(&one());
        for ipart in 1..=npart {
            let _it1 = n1extr[ipart];
            let _it2 = n2extr[ipart];
            let mut winterwr = winter[(ipart, idata)];
            let mut wintrawr = wintra[(ipart, idata)];
            let mut relativewt = winter[(ipart, idata)] * cprod[ipart];
            if nstruct == 1 {
                winterwr /= wtsum[idata];
                wintrawr /= wtsum[idata];
                relativewt /= wtsum[idata];
            }
            out.line(&fmt("(a8)", &[FV::Ch(&partlabel[ipart], 8)]));
            wline(&mut out, wintrawr * rhofac, winterwr * rhofac, relativewt);
        }
        out.close()?;
    }
    if norder == 2 {
        let nstruct = rd_int(ctx)?.ok_or_else(eof)?;
        let fname = format!("{base}dif.wts");
        st.log(&format!("\nepsrwts> Writing 1st order difference weights to file {fname}"));
        let mut out = OutFile::create(std::path::Path::new(&fname));
        let rhofac: f32 = if nrho == 2 {
            if nstruct == 1 {
                1.0 / atmolrat
            } else {
                atmolrat
            }
        } else {
            1.0
        };
        out.line(&one());
        let wtsumdif = wtsum[1] - wtsum[2];
        for ipart in 1..=npart {
            let it1 = n1extr[ipart] as usize;
            let it2 = n2extr[ipart] as usize;
            if iflagsubs[it1] == 1 || iflagsubs[it2] == 1 {
                let mut winterwr = winter[(ipart, 1)] - winter[(ipart, 2)];
                let mut wintrawr = wintra[(ipart, 1)] - wintra[(ipart, 2)];
                let mut relativewt = cprod[ipart] * winterwr;
                if nstruct == 1 && wtsumdif != 0.0 {
                    winterwr /= wtsumdif;
                    wintrawr /= wtsumdif;
                    relativewt /= wtsumdif;
                }
                out.line(&fmt("(a8)", &[FV::Ch(&partlabel[ipart], 8)]));
                wline(&mut out, wintrawr * rhofac, winterwr * rhofac, relativewt);
            }
        }
        out.close()?;
    }
    if norder == 3 {
        let rhofac: f32 = if nrho == 2 { 1.0 / atmolrat } else { 1.0 };
        let mut o10 = OutFile::create(std::path::Path::new(&format!("{base}xx.wts")));
        o10.line(&one());
        let mut o11 = OutFile::create(std::path::Path::new(&format!("{base}xy.wts")));
        o11.line(&one());
        let mut o12 = OutFile::create(std::path::Path::new(&format!("{base}yy.wts")));
        o12.line(&one());
        let mut cxbx: f32 = 0.0;
        let mut cyby1: f32 = 0.0;
        let mut cyby2: f32 = 0.0;
        for it in 1..=nt {
            if iflagsubs[it] == 1 {
                cyby1 += cfrac[it] * binter[(it, 1)];
                cyby2 += cfrac[it] * binter[(it, 2)];
            } else {
                cxbx += cfrac[it] * binter[(it, 1)];
            }
        }
        let wtsumyy = frac1 * wtsum[1] + frac2 * wtsum[2] - wtsum[3];
        let wtsumxy = 2.0 * cxbx * (cyby1 - cyby2);
        let wtsumxx = cxbx * cxbx;
        let wtsumyydif = cyby1 * cyby1 - cyby2 * cyby2;
        let wtsumyy1 = cyby1 * cyby1;
        let wtsumxy1 = 2.0 * cxbx * cyby1;
        for ipart in 1..=npart {
            let it1 = n1extr[ipart] as usize;
            let it2 = n2extr[ipart] as usize;
            let (s1, s2) = (iflagsubs[it1], iflagsubs[it2]);
            let (e1, e2) = (iexchange[it1], iexchange[it2]);
            let (winteryy, wintrayy);
            if s1 == 1 && s2 == 1 {
                winteryy = (frac1 * winter[(ipart, 1)] + frac2 * winter[(ipart, 2)] - winter[(ipart, 3)]) / wtsumyy;
                wintrayy = if e1 == 1 || e2 == 1 {
                    (frac1 * wintra[(ipart, 1)] + frac2 * wintra[(ipart, 2)] - wintra[(ipart, 3)]) / wtsumyy
                } else {
                    0.0
                };
            } else {
                winteryy = 0.0;
                wintrayy = 0.0;
            }
            let relativeyy = cprod[ipart] * winteryy;
            let (winterxy, wintraxy);
            if s1 != s2 {
                winterxy = (winter[(ipart, 1)] - winter[(ipart, 2)] - wtsumyydif * winteryy) / wtsumxy;
                wintraxy = (wintra[(ipart, 1)] - wintra[(ipart, 2)] - wtsumyydif * wintrayy) / wtsumxy;
            } else if s1 == 1 && s2 == 1 && e1 == 0 && e2 == 0 {
                winterxy = 0.0;
                wintraxy = (wintra[(ipart, 1)] - wintra[(ipart, 2)] - wtsumyydif * wintrayy) / wtsumxy;
            } else {
                winterxy = 0.0;
                wintraxy = 0.0;
            }
            let relativexy = cprod[ipart] * winterxy;
            let (winterxx, wintraxx);
            if s1 == 0 && s2 == 0 {
                winterxx = (winter[(ipart, 1)] - wtsumxy1 * winterxy - wtsumyy1 * winteryy) / wtsumxx;
                wintraxx = (wintra[(ipart, 1)] - wtsumxy1 * wintraxy - wtsumyy1 * wintrayy) / wtsumxx;
            } else if s1 == 1 && s2 == 1 && e1 == 0 && e2 == 0 {
                winterxx = 0.0;
                wintraxx = (wintra[(ipart, 1)] - wtsumxy1 * wintraxy - wtsumyy1 * wintrayy) / wtsumxx;
            } else {
                winterxx = 0.0;
                wintraxx = 0.0;
            }
            let relativexx = cprod[ipart] * winterxx;
            let lab = fmt("(a8)", &[FV::Ch(&partlabel[ipart], 8)]);
            if winteryy != 0.0 || wintrayy != 0.0 {
                o12.line(&lab);
                wline(&mut o12, wintrayy * rhofac, winteryy * rhofac, relativeyy);
            }
            if winterxy != 0.0 || wintraxy != 0.0 {
                o11.line(&lab);
                wline(&mut o11, wintraxy * rhofac, winterxy * rhofac, relativexy);
            }
            if winterxx != 0.0 || wintraxx != 0.0 {
                o10.line(&lab);
                wline(&mut o10, wintraxx * rhofac, winterxx * rhofac, relativexx);
            }
        }
        o10.close()?;
        o11.close()?;
        o12.close()?;
    }
    st.moloutc(&fnameato)?;
    Ok(())
}
