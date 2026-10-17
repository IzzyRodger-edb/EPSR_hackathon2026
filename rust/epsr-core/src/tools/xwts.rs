//! `xwts` — translation of `fortran_src/xwts.f` (subroutines `xwts` and
//! `form_fact`) and of the `symbol_setup.eq.'xwts'` branch of `check_symbol`
//! (`process_symbols.f`).
//!
//! STATUS: translated, verified byte-identical (EPSR-test-ex23).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> xwts <input file>`
//! Keyword table: `crate::symtab_tools::XWTS`.
//!
//! Looks up the Waasmaier-Kirfel X-ray form-factor coefficients of every
//! atom type in `f0_WaasKirf.dat` (read from the *current working
//! directory*, as the Fortran does) and writes them, with the MAFF alpha and
//! delta values, to `<run>.wts`; then rewrites the input file.

use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::{read_int, read_real, ListRead, TextFile};
use crate::fortran::output::OutFile;
use crate::fortran::text::{fstr, parse};
use crate::state::State;
use crate::symtab_tools::XWTS;

/// `subroutine form_fact(chsym,axt,bxt,cxt,ifind,ndim)` with ndim = 5.
/// Returns `Some((a, b, c))` if the symbol was found (`ifind = 1`).
pub fn form_fact(chsym: &str, log: &dyn Fn(&str)) -> Result<Option<([f64; 5], [f64; 5], f64)>> {
    let chsym = fstr(chsym, 4);
    if chsym.is_empty() {
        return Ok(None);
    }
    let mut f = match TextFile::open(std::path::Path::new("f0_WaasKirf.dat")) {
        Ok(f) => f,
        Err(_) => {
            log(&fmt("(/'form_factor> X-ray file ',a15,' does not exist')", &[FV::Ch("f0_WaasKirf.dat", 15)]));
            return Ok(None);
        }
    };
    loop {
        let sometext = match f.read_line() {
            Ok(l) => l,
            Err(_) => break,
        };
        if sometext.contains("#S") {
            let w = parse(&sometext, 5);
            if w.nwords() < 3 {
                return crate::error::unsupported("xwts: '#S' line with fewer than 3 words in f0_WaasKirf.dat");
            }
            if w.word(&sometext, 3) == chsym {
                let r1 = f.skip_record();
                let r2 = if r1.is_ok() { f.skip_record() } else { r1.clone() };
                let r3 = if r2.is_ok() { f.read_line() } else { Err(r2.clone().unwrap_err()) };
                // the Fortran ignores the iostat of the skips; ierr is that of the last read
                let (ierr_ok, line) = match r3 {
                    Ok(l) => (true, l),
                    Err(_) => (false, String::new()),
                };
                let w2 = parse(&line, 5);
                if ierr_ok && w2.nwords() == 11 {
                    let mut r = ListRead::from_str(&line);
                    let mut a = [0.0f64; 5];
                    let mut b = [0.0f64; 5];
                    for x in a.iter_mut() {
                        *x = r.dble()?.ok_or(EpsrError::Input("form_fact: bad coefficient line".into()))?;
                    }
                    let c = r.dble()?.ok_or(EpsrError::Input("form_fact: bad coefficient line".into()))?;
                    for x in b.iter_mut() {
                        *x = r.dble()?.ok_or(EpsrError::Input("form_fact: bad coefficient line".into()))?;
                    }
                    return Ok(Some((a, b, c)));
                }
                if !ierr_ok {
                    break;
                }
            }
        }
    }
    Ok(None)
}

/// Run the `xwts` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    if !ctx.open_input(&mut st, &XWTS)? {
        return Ok(());
    }
    st.log("\nxwts> Program to setup the X-ray form factors forfor the specified .ato file.\n");
    let mut iref: i32 = 1;
    let wtsfile = st.sym.fname_run.clone();
    iref += 1;
    let nrho = read_int(&st.sym.symbol_value[iref])?;
    iref += 1;
    let normx = read_int(&st.sym.symbol_value[iref])?;
    iref += 1;
    let iafformaff = read_int(&st.sym.symbol_value[iref])?;
    iref += 1;
    let ncomponent = read_int(&st.sym.symbol_value[iref])?;
    let ntype = st.mol.ntype;
    if ncomponent != ntype {
        st.log("\nxwts> Number of components in input file different from number of atom types in .ato file");
        return Ok(());
    }
    let nt = ntype as usize;
    let mut ax = vec![[0.0f64; 5]; nt + 1];
    let mut bx = vec![[0.0f64; 5]; nt + 1];
    let mut cx = vec![0.0f64; nt + 1];
    let mut alphax = vec![0.0f64; nt + 1];
    let mut deltax = vec![0.0f64; nt + 1];
    let mut ffsym = vec![String::new(); nt + 1];
    let log = st.logger();
    for it in 1..=nt {
        iref += 1;
        if st.mol.qtypea[it] == 0 {
            iref += 1;
            let mut chsymtest = fstr(&st.sym.symbol_value[iref], 4);
            let res;
            loop {
                if let Some(r) = form_fact(&chsymtest, &log)? {
                    res = r;
                    break;
                }
                log(&fmt(
                    "(/'xwts> Did not recognise form factor symbol ',a4,' for atom type ',a3//'xwts> Try typing it again: ')",
                    &[FV::Ch(&chsymtest, 4), FV::Ch(&st.mol.type_[it], 3)],
                ));
                let mut r = ctx.stdin.list();
                match r.word() {
                    Ok(Some(w)) => chsymtest = fstr(&w, 4),
                    Ok(None) => {}
                    Err(_) => return Err(EpsrError::Input("xwts: end of input while asking for a form factor symbol".into())),
                }
            }
            st.sym.symbol_value[iref] = chsymtest.clone();
            ax[it] = res.0;
            bx[it] = res.1;
            cx[it] = res.2;
            ffsym[it] = chsymtest;
            iref += 1;
            alphax[it] = read_real(&st.sym.symbol_value[iref])? as f64;
            iref += 1;
            deltax[it] = read_real(&st.sym.symbol_value[iref])? as f64;
        } else {
            // cx(it) is left uninitialised by the Fortran for Q atoms
            return crate::error::unsupported("xwts: Q-type atoms");
        }
    }
    let fname = format!("{}{}.wts", st.sym.directory_name, wtsfile);
    let mut out = OutFile::create(std::path::Path::new(&fname));
    out.line(&fmt("(3(1x,i3))", &[FV::I(nrho as i64), FV::I(normx as i64), FV::I(iafformaff as i64)]));
    for it in 1..=nt {
        out.line(&fmt("(1x,a3,1x,a4)", &[FV::Ch(&st.mol.type_[it], 3), FV::Ch(&ffsym[it], 4)]));
        let mut v: Vec<FV> = vec![FV::D(cx[it])];
        for i in 0..5 {
            v.push(FV::D(ax[it][i]));
            v.push(FV::D(bx[it][i]));
        }
        v.push(FV::D(alphax[it]));
        v.push(FV::D(deltax[it]));
        out.line(&fmt("(11(1x,f11.6))", &v));
    }
    out.close()?;
    st.save_inp_file()?;
    Ok(())
}

/// Program-specific branch of `check_symbol` (`process_symbols.f`, `symbol_setup.eq.'xwts'`).
/// Returns the Fortran `ierr` (0 = accepted).
pub fn check_symbol_extra(st: &mut State, icount: i32, _ifilecount: i32, iref: i32) -> Result<i32> {
    let sym = st.sym.symbol[icount].clone();
    let nup = st.sym.nsymboluptodata;
    if sym == "fnameato" && st.sym.symbol_value[iref] != st.sym.symbol_default[iref] {
        let nfilet = st.mol.ntype;
        st.sym.symbol_value[nup] = if nfilet < 1000 { format!("{nfilet}") } else { "***".into() };
        let mut ifiletemp = st.sym.nfiles;
        st.sym.nfiles = nfilet;
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
            if st.sym.symbol_value[il + 2] == st.sym.symbol_default[nup + 2] {
                st.sym.symbol_value[il + 2] = st.mol.atsymbol[if1].clone();
            }
            il += st.sym.nsymbolplus;
        }
    } else if sym == "output" {
        let itest = read_int(&st.sym.symbol_value[iref])?.min(2).max(1);
        st.sym.symbol_value[iref] = format!("{itest}");
    } else if sym == "normtot" {
        let itest = read_int(&st.sym.symbol_value[iref])?;
        st.log(&{
            let mut w = ListWriter::new();
            w.int(itest);
            w.finish()
        });
        let itest = itest.min(2).max(0);
        st.sym.symbol_value[iref] = format!("{itest}");
    }
    Ok(0)
}
