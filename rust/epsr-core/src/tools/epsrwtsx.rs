//! `epsrwtsx` — translation of `fortran_src/epsrwtsx.f` (subroutine `epsrwtsx`;
//! `form_fact` is shared with [`super::xwts`]).
//!
//! STATUS: translated, verified byte-identical (EPSR-test-ex25).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> epsrwtsx <file>.ato < answers`
//!
//! Interactive predecessor of `xwts`: asks for the form-factor symbol and
//! the MAFF q/delta values of every atom type, the normalisation type and
//! the output base name, and writes `<base>X-ray.wts`.  `f0_WaasKirf.dat` is
//! read from the current working directory.

use super::xwts::form_fact;
use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::IoErr;
use crate::fortran::output::OutFile;
use crate::fortran::text::{fstr, parse};

fn eof() -> EpsrError {
    EpsrError::Input("epsrwtsx: unexpected end of standard input".into())
}

/// Run the `epsrwtsx` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    let readfile0 = match ctx.args.first() {
        Some(a) => a.clone(),
        None => return crate::error::unsupported("epsrwtsx without a .ato file argument (interactive file search)"),
    };
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
        st.log(&format!("\nepsrwtsx> Problems with .ato file {fnameato}\nRoutine will exit now\n"));
        return Ok(());
    }
    let mut nrho = {
        let mut r = ctx.stdin.list();
        r.int()?.ok_or_else(eof)?
    };
    if nrho < 1 || nrho > 2 {
        nrho = 2;
    }
    if nrho == 2 {
        nrho = -1;
    }
    let nt = st.mol.ntype as usize;
    let mut ax = vec![[0.0f64; 5]; nt + 1];
    let mut bx = vec![[0.0f64; 5]; nt + 1];
    let mut cx = vec![0.0f64; nt + 1];
    let mut alphax = vec![0.0f64; nt + 1];
    let mut deltax = vec![0.0f64; nt + 1];
    let mut ffsym = vec![String::new(); nt + 1];
    let log = st.logger();
    let mut chsymtest = String::new();
    for it in 1..=nt {
        if st.mol.qtypea[it] != 0 {
            // cx(it) is left uninitialised by the Fortran for Q atoms
            return crate::error::unsupported("epsrwtsx: Q-type atoms");
        }
        // read(5,*,iostat=ierr) chsymtest
        let read_sym = |ctx: &mut ToolContext, cur: &mut String| -> Result<()> {
            let mut r = ctx.stdin.list();
            match r.word() {
                Ok(Some(w)) => *cur = fstr(&w, 4),
                Ok(None) => {}
                Err(IoErr::Eof) => return Err(eof()),
                Err(_) => {}
            }
            Ok(())
        };
        read_sym(ctx, &mut chsymtest)?;
        let res = loop {
            if let Some(r) = form_fact(&chsymtest, &log)? {
                break r;
            }
            log(&format!("\nepsrwtsx> Did not recognise form factor symbol {chsymtest:<4} for atom type {:<3}", st.mol.type_[it]));
            read_sym(ctx, &mut chsymtest)?;
        };
        ax[it] = res.0;
        bx[it] = res.1;
        cx[it] = res.2;
        ffsym[it] = chsymtest.clone();
        // read(5,*,iostat=ierr) alphax(it),deltax(it)
        let mut r = ctx.stdin.list();
        let mut a = alphax[it];
        let mut d = deltax[it];
        let res = r.get_d(&mut a).and_then(|_| r.get_d(&mut d));
        drop(r);
        alphax[it] = a;
        deltax[it] = d;
        match res {
            Ok(()) => {}
            Err(IoErr::Eof) => return Err(eof()),
            Err(_) => return crate::error::unsupported("epsrwtsx: invalid MAFF q/delta answer"),
        }
    }
    let normx = {
        let mut r = ctx.stdin.list();
        r.int()?.ok_or_else(eof)?
    };
    let mut iafformaff: i32 = 0;
    {
        let mut r = ctx.stdin.list();
        match r.get_i(&mut iafformaff) {
            Ok(()) => {}
            Err(IoErr::Eof) => return Err(eof()),
            Err(_) => return crate::error::unsupported("epsrwtsx: invalid IAFF/MAFF answer"),
        }
    }
    let base = loop {
        let sometext = ctx.stdin.read_line().map_err(|_| eof())?;
        let w = parse(&sometext, 5);
        if w.nwords() > 0 {
            break format!("{}{}", st.sym.directory_name, w.word(&sometext, 1));
        }
    };
    let fname = format!("{base}X-ray.wts");
    st.log(&format!("\nepsrwtsx> Writing X-ray weights to file {fname}"));
    let mut out = OutFile::create(std::path::Path::new(&fname));
    {
        let mut w = ListWriter::new();
        w.int(nrho).int(normx).int(iafformaff);
        out.line(&w.finish());
    }
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
    Ok(())
}
