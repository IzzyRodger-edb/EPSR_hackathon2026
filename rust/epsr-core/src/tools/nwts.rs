//! `nwts` — translation of `fortran_src/nwts.f` (subroutine `nwts`) and of the
//! `symbol_setup.eq.'nwts'` branch of `check_symbol` (`process_symbols.f`).
//!
//! STATUS: translated, verified byte-identical (EPSR-test-ex20, ex21, ex22).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> nwts <input file>`
//! Keyword table: `crate::symtab_tools::NWTS`.
//!
//! Calculates the neutron weights (intra- and inter-molecular scattering
//! length products) of every atom-type pair for one sample, a first-order
//! difference (2 samples) or a second-order difference (3 samples), writes
//! them to `<run>tot.wts` / `<run>tot{n}.wts`, `<run>dif.wts`,
//! `<run>xx.wts`/`xy.wts`/`yy.wts`, rewrites the `.ato` file (the chemical
//! symbols / exchange flags may have changed) and the input file.
//!
//! Also holds `scatt_len` (`ato_manipulation_rev4.f`), used by `epsrwts` and
//! `changeato` as well.

use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::{read_int, read_real};
use crate::fortran::output::OutFile;
use crate::fortran::text::{fstr, parse};
use crate::fortran::{F2, F3};
use crate::state::{State, MCOMP, MPAIRS};
use crate::symtab_tools::NWTS;

#[path = "sears91.rs"]
mod sears91;

/// `parameter(misotope=5)`
pub const MISOTOPE: usize = 5;

/// `function scatt_len(chsym,massno,ifind)`: neutron scattering length of the
/// isotope `massno` (0 = natural element) of chemical symbol `chsym` from the
/// table `sears91_gudrun.inc`.  Returns `(scatt_len, ifind)`.
pub fn scatt_len(chsym: &str, massno: i32) -> (f32, i32) {
    let c = fstr(chsym, 3);
    for (s, m, b) in sears91::SEARS91.iter() {
        if c == *s && massno == *m {
            return (*b, 1);
        }
    }
    (0.0, 0)
}

/// List-directed write of a single integer into a character variable
/// (`write(symbol_value(iref),*) i`), trailing blanks dropped.
pub(crate) fn list_int(v: i32) -> String {
    let mut w = ListWriter::new();
    w.int(v);
    w.finish().trim_end().to_string()
}

/// Read one list-directed character item from standard input into a
/// `character*len` variable (`read(5,*) c`).
fn read_stdin_word(ctx: &mut ToolContext, len: usize) -> Result<String> {
    let mut r = ctx.stdin.list();
    let w = r.word()?.unwrap_or_default();
    Ok(fstr(&w, len))
}

/// Parse a list of (mass number, abundance) pairs as `nwts` does; returns
/// the scattering lengths and normalised abundances.
fn read_isotopes(
    ctx: &mut ToolContext,
    st: &State,
    sometext: &str,
    it: i32,
    chsymtest: &str,
    btest: f32,
) -> Result<Vec<(f32, f32)>> {
    let w = parse(sometext, 5);
    let nwords = w.nwords();
    let mut asum: f32 = 0.0;
    let mut nsum = 0;
    let mut v: Vec<(f32, f32)> = Vec::new();
    while 2 * nsum < nwords {
        nsum += 1;
        let j2 = 2 * nsum;
        let j1 = j2 - 1;
        let mut massno = read_int(w.word(sometext, j1))?;
        // an odd last word is read past the end of the word list: the Fortran
        // then reads ncf/ncl of the next (unset) word
        if j2 > nwords {
            return crate::error::unsupported("nwts: odd number of entries in an abundance list");
        }
        let ab = read_real(w.word(sometext, j2))?;
        asum += ab;
        let b;
        if massno == 0 {
            b = btest;
        } else {
            let (mut bb, mut ifind) = scatt_len(chsymtest, massno);
            while ifind == 0 {
                st.log(&fmt(
                    "(/'nwts> Did not find massno ',i3,' for chemical',' symbol ',a3,' and abundance ',f10.5,//'Try typing again: ')",
                    &[FV::I(massno as i64), FV::Ch(&st.mol.atsymbol[it], 3), FV::R(ab)],
                ));
                let mut r = ctx.stdin.list();
                massno = r.int()?.ok_or(EpsrError::Input("nwts: mass number expected".into()))?;
                drop(r);
                let t = scatt_len(chsymtest, massno);
                bb = t.0;
                ifind = t.1;
            }
            b = bb;
        }
        v.push((b, ab));
    }
    if v.len() > MISOTOPE {
        return Err(EpsrError::Input(format!("nwts: more than {MISOTOPE} isotopes for one component")));
    }
    for e in v.iter_mut() {
        e.1 /= asum;
    }
    Ok(v)
}

/// Run the `nwts` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    if !ctx.open_input(&mut st, &NWTS)? {
        // no input file: the Fortran continues with the defaults, finds
        // ncomponent (1) different from ntype (0) and returns.
        return Ok(());
    }
    st.log(&fmt(
        "(/'nwts> Program to calculate inter- and intra-molecular',' weightings ',/1x,'for DCS, 1st- or 2nd-order difference data.'/)",
        &[],
    ));
    let mut iref: i32 = 1;
    let wtsfile = st.sym.fname_run.clone();
    iref += 1;
    let nrho = read_int(&st.sym.symbol_value[iref])?;
    let ntype = st.mol.ntype;
    let atmolrat: f32 = st.mol.natomr as f32 / st.mol.nmol as f32;
    iref += 1;
    let norder = read_int(&st.sym.symbol_value[iref])?;
    iref += 1;
    let fracmix = read_real(&st.sym.symbol_value[iref])?;
    iref += 1;
    let normtot = read_int(&st.sym.symbol_value[iref])?;
    iref += 1;
    let normdif = read_int(&st.sym.symbol_value[iref])?;
    iref += 1;
    let ncomponent = read_int(&st.sym.symbol_value[iref])?;
    if ncomponent != ntype {
        st.log("\nnwts> Number of components in input file different from number of atom types in .ato file");
        return Ok(());
    }
    let nt = ntype as usize;
    let mut iflagsubs = vec![0i32; nt + 1];
    let mut nisotope: F2<i32> = F2::new(MCOMP, 2);
    let mut batom: F3<f32> = F3::new(MCOMP, MISOTOPE, 2);
    let mut abund: F3<f32> = F3::new(MCOMP, MISOTOPE, 2);
    if norder > 1 {
        st.log("\nnwts> For sample 1, get the scattering lengths for all components in the sample\n");
    } else {
        st.log("\nnwts> Get the scattering lengths for all components in the sample\n");
    }
    for it in 1..=ntype {
        let itu = it as usize;
        iref += 2;
        if st.mol.qtypea[itu] == 0 {
            let mut chsymtest = fstr(&st.mol.atsymbol[itu], 3);
            let massno0 = 0;
            let (mut btest, mut ifind) = scatt_len(&chsymtest, massno0);
            while ifind == 0 {
                st.log(&fmt(
                    "(/'nwts> Did not recognise chemical symbol ',a3,' for atom type ',a3//'nwts> Try typing it again: ')",
                    &[FV::Ch(&chsymtest, 3), FV::Ch(&st.mol.type_[itu], 3)],
                ));
                chsymtest = read_stdin_word(ctx, 3)?;
                let t = scatt_len(&chsymtest, massno0);
                btest = t.0;
                ifind = t.1;
            }
            st.mol.atsymbol[itu] = chsymtest.clone();
            if st.mol.atsymbol[itu] == "H" {
                iref += 1;
                st.mol.iexchange[itu] = read_int(&st.sym.symbol_value[iref])?;
            } else {
                iref += 1;
                st.mol.iexchange[itu] = 0;
            }
            // write(symbol_value(iref),'(i1)') iexchange(it)
            st.sym.symbol_value[iref] = fmt("(i1)", &[FV::I(st.mol.iexchange[itu] as i64)]);
            iref += 1;
            let sometext = st.sym.symbol_value[iref].clone();
            let v = read_isotopes(ctx, &st, &sometext, itu as i32, &chsymtest, btest)?;
            for (j, (b, a)) in v.iter().enumerate() {
                batom[(itu, j + 1, 1)] = *b;
                abund[(itu, j + 1, 1)] = *a;
            }
            nisotope[(itu, 1)] = v.len() as i32;
            log_isotopes(&st, itu, 1, &nisotope, &batom, &abund);
            iref += 1;
            iflagsubs[itu] = read_int(&st.sym.symbol_value[iref])?;
            if iflagsubs[itu] == 1 {
                iref += 1;
                let sometext = st.sym.symbol_value[iref].clone();
                let v = read_isotopes(ctx, &st, &sometext, itu as i32, &chsymtest, btest)?;
                for (j, (b, a)) in v.iter().enumerate() {
                    batom[(itu, j + 1, 2)] = *b;
                    abund[(itu, j + 1, 2)] = *a;
                }
                nisotope[(itu, 2)] = v.len() as i32;
            } else {
                iref += 1;
                nisotope[(itu, 2)] = nisotope[(itu, 1)];
                for isum in 1..=nisotope[(itu, 2)] as usize {
                    batom[(itu, isum, 2)] = batom[(itu, isum, 1)];
                    abund[(itu, isum, 2)] = abund[(itu, isum, 1)];
                }
            }
            log_isotopes(&st, itu, 2, &nisotope, &batom, &abund);
        } else {
            iref += 4;
        }
    }
    // exchangeable atoms: average scattering length
    let mut bscatsum = [0.0f32; 3];
    let mut cfracsum: f32 = 0.0;
    for it in 1..=nt {
        if st.mol.iexchange[it] == 1 {
            let cf = st.mol.cfrac[it];
            cfracsum += cf;
            for ii in 1..=nisotope[(it, 1)] as usize {
                let term = cf * abund[(it, ii, 1)] * batom[(it, ii, 1)];
                bscatsum[1] += term;
            }
            for ii in 1..=nisotope[(it, 2)] as usize {
                let term = cf * abund[(it, ii, 2)] * batom[(it, ii, 2)];
                bscatsum[2] += term;
            }
        }
    }
    if cfracsum > 0.0 {
        for ii in 1..=2 {
            bscatsum[ii] /= cfracsum;
        }
    } else {
        for ii in 1..=2 {
            bscatsum[ii] = 0.0;
        }
    }
    for it in 1..=nt {
        if st.mol.iexchange[it] == 1 {
            nisotope[(it, 1)] = 1;
            abund[(it, 1, 1)] = 1.0;
            batom[(it, 1, 1)] = bscatsum[1];
            nisotope[(it, 2)] = 1;
            abund[(it, 1, 2)] = 1.0;
            batom[(it, 1, 2)] = bscatsum[2];
        }
    }
    let (frac1, frac2): (f32, f32) = if norder == 3 { (fracmix, 1.0 - fracmix) } else { (1.0, 0.0) };
    if norder == 3 && (1..=nt).any(|it| st.mol.qtypea[it] != 0) {
        // binter(it,3) of Q atoms is used uninitialised by the Fortran (normtot=2)
        return crate::error::unsupported("nwts: second-order difference with Q-type atoms");
    }
    let npart = (ntype * (ntype + 1) / 2) as usize;
    let ndata = norder.min(2);
    let mut binter: F2<f32> = F2::new(MCOMP, 3);
    let mut wintra: F2<f32> = F2::new(MPAIRS, 3);
    let mut winter: F2<f32> = F2::new(MPAIRS, 3);
    let mut cprod = vec![0.0f32; MPAIRS + 1];
    let mut partlabel = vec![String::new(); MPAIRS + 1];
    let mut wtsum = [0.0f32; 4];
    let cfrac = st.mol.cfrac.clone();
    let qtypea = st.mol.qtypea.clone();
    let iexchange = st.mol.iexchange.clone();
    for idata in 1..=ndata.max(0) as usize {
        for it1 in 1..=nt {
            binter[(it1, idata)] = 0.0;
            if qtypea[it1] == 0 {
                for ib in 1..=nisotope[(it1, idata)] as usize {
                    let term1 = abund[(it1, ib, idata)] * batom[(it1, ib, idata)];
                    binter[(it1, idata)] += term1;
                }
            }
        }
        let mut ipart = 0usize;
        wtsum[idata] = 0.0;
        for it1 in 1..=nt {
            let binter1 = binter[(it1, idata)];
            if normtot == 2 && qtypea[it1] == 0 {
                wtsum[idata] += cfrac[it1] * binter1 * binter1;
            }
            for it2 in it1..=nt {
                ipart += 1;
                if qtypea[it1] == 0 && qtypea[it2] == 0 {
                    if idata == 1 {
                        partlabel[ipart] = fmt("(1x,a3,1x,a3)", &[FV::Ch(&st.mol.type_[it1], 3), FV::Ch(&st.mol.type_[it2], 3)]);
                        cprod[ipart] = if it1 == it2 { cfrac[it1] * cfrac[it1] } else { 2.0 * cfrac[it1] * cfrac[it2] };
                    }
                    let binter2 = binter[(it2, idata)];
                    winter[(ipart, idata)] = binter1 * binter2;
                    wintra[(ipart, idata)] = winter[(ipart, idata)];
                    if normtot < 2 {
                        let term = cprod[ipart] * winter[(ipart, idata)];
                        wtsum[idata] += term;
                    }
                }
            }
        }
    }
    if norder == 3 {
        wtsum[3] = 0.0;
        for it1 in 1..=nt {
            if qtypea[it1] == 0 {
                binter[(it1, 3)] = frac1 * binter[(it1, 1)] + frac2 * binter[(it1, 2)];
            }
        }
        let mut ipart = 0usize;
        for it1 in 1..=nt {
            if normtot == 2 {
                wtsum[3] += cfrac[it1] * binter[(it1, 3)] * binter[(it1, 3)];
            }
            for it2 in it1..=nt {
                ipart += 1;
                if qtypea[it1] == 0 && qtypea[it2] == 0 {
                    winter[(ipart, 3)] = binter[(it1, 3)] * binter[(it2, 3)];
                    if iexchange[it1] == 0 && iexchange[it2] == 0 {
                        wintra[(ipart, 3)] = frac1 * wintra[(ipart, 1)] + frac2 * wintra[(ipart, 2)];
                    } else {
                        wintra[(ipart, 3)] = winter[(ipart, 3)];
                    }
                    if normtot < 2 {
                        let term = cprod[ipart] * winter[(ipart, 3)];
                        wtsum[3] += term;
                    }
                }
            }
        }
    }
    // ---- write the weights files
    let base = format!("{}{}", st.sym.directory_name, wtsfile);
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
    for idata in 1..=norder.max(0) as usize {
        let nstruct = normtot;
        let fnameb = if norder > 1 { fmt("('tot',i1.1,'.wts')", &[FV::I(idata as i64)]) } else { "tot.wts".to_string() };
        let fname = format!("{base}{fnameb}");
        st.log(&format!("\nnwts> Writing TOTAL weights to file {fname}"));
        let mut out = OutFile::create(std::path::Path::new(&fname));
        let rhofac: f32 = if nrho == 2 {
            if nstruct == 0 {
                atmolrat
            } else if nstruct == 1 {
                1.0 / atmolrat
            } else {
                1.0
            }
        } else {
            1.0
        };
        out.line(&one());
        for ipart in 1..=npart {
            let it1 = st.mol.n1extr[ipart] as usize;
            let it2 = st.mol.n2extr[ipart] as usize;
            if qtypea[it1] == 0 && qtypea[it2] == 0 {
                let mut winterwr = winter[(ipart, idata)];
                let mut wintrawr = wintra[(ipart, idata)];
                let mut relativewt = winter[(ipart, idata)] * cprod[ipart];
                if nstruct > 0 {
                    winterwr /= wtsum[idata];
                    wintrawr /= wtsum[idata];
                    relativewt /= wtsum[idata];
                }
                out.line(&fmt("(a8)", &[FV::Ch(&partlabel[ipart], 8)]));
                wline(&mut out, wintrawr * rhofac, winterwr * rhofac, relativewt);
            }
        }
        out.close()?;
    }
    if norder == 2 {
        let nstruct = normdif;
        let fname = format!("{base}dif.wts");
        st.log(&format!("\nnwts> Writing 1st order difference weights to file {fname}"));
        let mut out = OutFile::create(std::path::Path::new(&fname));
        let rhofac: f32 = if nrho == 2 {
            if nstruct > 0 {
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
            let it1 = st.mol.n1extr[ipart] as usize;
            let it2 = st.mol.n2extr[ipart] as usize;
            if qtypea[it1] == 0 && qtypea[it2] == 0 && (iflagsubs[it1] == 1 || iflagsubs[it2] == 1) {
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
        let fxx = format!("{base}xx.wts");
        st.log(&format!("\n Writing XX weights to file {fxx}"));
        let mut o10 = OutFile::create(std::path::Path::new(&fxx));
        let rhofac: f32 = if nrho == 2 { 1.0 / atmolrat } else { 1.0 };
        o10.line(&one());
        let fxy = format!("{base}xy.wts");
        st.log(&format!("\n Writing XY (=XH) weights to file {fxy}"));
        let mut o11 = OutFile::create(std::path::Path::new(&fxy));
        o11.line(&one());
        let fyy = format!("{base}yy.wts");
        st.log(&format!("\n Writing YY (=HH) weights to file {fyy}"));
        let mut o12 = OutFile::create(std::path::Path::new(&fyy));
        o12.line(&one());
        let mut cxbx: f32 = 0.0;
        let mut cyby1: f32 = 0.0;
        let mut cyby2: f32 = 0.0;
        for it in 1..=nt {
            if qtypea[it] == 0 {
                if iflagsubs[it] == 1 {
                    cyby1 += cfrac[it] * binter[(it, 1)];
                    cyby2 += cfrac[it] * binter[(it, 2)];
                } else {
                    cxbx += cfrac[it] * binter[(it, 1)];
                }
            }
        }
        let wtsumyy = frac1 * wtsum[1] + frac2 * wtsum[2] - wtsum[3];
        let wtsumxy = 2.0 * cxbx * (cyby1 - cyby2);
        let wtsumxx = cxbx * cxbx;
        let wtsumyydif = cyby1 * cyby1 - cyby2 * cyby2;
        let wtsumyy1 = cyby1 * cyby1;
        let wtsumxy1 = 2.0 * cxbx * cyby1;
        for ipart in 1..=npart {
            let it1 = st.mol.n1extr[ipart] as usize;
            let it2 = st.mol.n2extr[ipart] as usize;
            if !(qtypea[it1] == 0 && qtypea[it2] == 0) {
                continue;
            }
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
    // write the .ato file in case anything was changed
    let fnameato = st.sym.fnameato.clone();
    st.moloutc(&fnameato)?;
    // finally rewrite the input file
    st.save_inp_file()?;
    Ok(())
}

fn log_isotopes(st: &State, it: usize, k: usize, nisotope: &F2<i32>, batom: &F3<f32>, abund: &F3<f32>) {
    if !st.verbose {
        return;
    }
    let n = nisotope[(it, k)];
    let mut v: Vec<FV> = vec![FV::I(n as i64)];
    for j in 1..=n as usize {
        v.push(FV::Ch(&st.mol.type_[it], 3));
        v.push(FV::R(batom[(it, j, k)]));
        v.push(FV::R(abund[(it, j, k)]));
    }
    st.log(&fmt("(i3,1x,((a3,1x,f10.5,1x,f10.5)))", &v));
}

/// Program-specific branch of `check_symbol` (`process_symbols.f`, `symbol_setup.eq.'nwts'`).
/// Returns the Fortran `ierr` (0 = accepted).
pub fn check_symbol_extra(st: &mut State, icount: i32, _ifilecount: i32, iref: i32) -> Result<i32> {
    let sym = st.sym.symbol[icount].clone();
    let nup = st.sym.nsymboluptodata;
    if sym == "fnameato" && st.sym.symbol_value[iref] != st.sym.symbol_default[iref] {
        let nfilet = st.mol.ntype;
        st.sym.symbol_value[nup] = if nfilet < 1000 { format!("{nfilet}") } else { "***".into() };
        let mut ifiletemp = st.sym.nfiles;
        st.sym.nfiles = nfilet;
        st.log(&{
            let mut w = ListWriter::new();
            w.int(nfilet).int(ifiletemp).int(st.sym.nfiles);
            w.finish()
        });
        while ifiletemp < st.sym.nfiles {
            ifiletemp += 1;
            let mut isymbol = nup;
            while isymbol < st.sym.nsymboldata {
                isymbol += 1;
                let r = st.sym.ireff(isymbol, ifiletemp);
                st.sym.symbol_value[r] = st.sym.symbol_default[isymbol].clone();
            }
        }
        let mut il = nup;
        let mut if1 = 0;
        while if1 < st.sym.nfiles {
            if1 += 1;
            st.sym.symbol_value[il + 1] = st.mol.type_[if1].clone();
            st.sym.symbol_value[il + 2] = st.mol.atsymbol[if1].clone();
            if st.sym.symbol_value[il + 2] != "H" {
                st.sym.symbol_value[il + 3] = "0".into();
            }
            il += st.sym.nsymbolplus;
        }
    } else if sym == "output" {
        let itest = read_int(&st.sym.symbol_value[iref])?.clamp(1, 2);
        st.sym.symbol_value[iref] = format!("{itest}");
    } else if sym == "nsamples" {
        let nsamples = read_int(&st.sym.symbol_value[iref])?.min(3).max(1);
        st.sym.symbol_value[iref] = format!("{nsamples}");
        if nsamples == 3 {
            st.sym.iflag[icount + 1] = 0;
        } else {
            st.sym.iflag[icount + 1] = 1;
            st.sym.symbol_value[iref + 1] = st.sym.symbol_default[icount + 1].clone();
        }
        if nsamples < 2 {
            st.sym.iflag[nup + 5] = 1;
            st.sym.iflag[nup + 6] = 1;
            let mut il = nup;
            for _if1 in 1..=st.sym.nfiles {
                st.sym.symbol_value[il + 5] = st.sym.symbol_default[nup + 5].clone();
                st.sym.symbol_value[il + 6] = st.sym.symbol_default[nup + 6].clone();
                il += st.sym.nsymbolplus;
            }
        } else {
            st.sym.iflag[nup + 5] = 0;
            st.sym.iflag[nup + 6] = 0;
        }
    } else if sym == "normtot" {
        let itest = read_int(&st.sym.symbol_value[iref])?.min(2).max(0);
        st.sym.symbol_value[iref] = list_int(itest);
    } else if sym == "normdif" {
        let itest = read_int(&st.sym.symbol_value[iref])?.min(1).max(0);
        st.sym.symbol_value[iref] = list_int(itest);
    } else if sym == "iexchange" {
        if st.sym.symbol_value[iref - 1] != "H" {
            st.sym.symbol_value[iref] = "0".into();
        }
        let itest = read_int(&st.sym.symbol_value[iref])?.min(1).max(0);
        st.sym.symbol_value[iref] = list_int(itest);
    } else if sym == "abundances" || sym == "abund2" {
        let sometext = st.sym.symbol_value[iref].clone();
        let w = parse(&sometext, 5);
        if w.nwords() != 2 * (w.nwords() / 2) {
            // interactive re-entry from standard input (not available while reading a file)
            return crate::error::unsupported("nwts: odd number of entries in abundances/abund2");
        }
        if w.nwords() > 0 {
            st.sym.symbol_value[iref] = w.span(&sometext, 1, w.nwords());
        }
    }
    Ok(0)
}
