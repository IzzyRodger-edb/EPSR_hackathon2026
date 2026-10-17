//! `readxyz` — translation of `fortran_src/readxyz.f` (`run_readxyz`,
//! `readxyz`).
//!
//! STATUS: translated; verified byte-identical against the Fortran built with
//! an explicit interface for `ctype` (see below), test cases EPSR-test-ex40/41.
//!
//! Invocation: `epsr-rs <folder> readxyz <file>.xyz`
//!
//! Reads an `.xyz` file (line 1: number of atoms; line 2: title, which must
//! contain `<word> a b c alpha beta gamma` to define the cell, otherwise the
//! six values are read from standard input; then `label x y z` lines, the
//! label being the chemical symbol) and writes `<file>.ato` with every atom as
//! a separate one-atom molecule, using a non-cubic box header.
//!
//! NOTE on the Fortran: `readxyz.f` calls `ctype` (which has optional
//! arguments) with 11 arguments and without an explicit interface.  gfortran
//! then passes the hidden character lengths (3, 3) where `ctype` expects the
//! addresses of `radmaxsph`/`pseudofacsph`, `present()` is true and the
//! program crashes (SIGSEGV) at the first `ctype` call.  This translation
//! implements what the call means with an explicit interface (absent
//! arguments: `radmaxsphere = 0.5*sigcor`, `pseudofacmol = 1`), and was
//! verified against the Fortran compiled that way.

use self::groupc_common as gc;
use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::{ListRead, TextFile};
use crate::fortran::text::{fstr, indexlastdot, parse};
#[allow(unused_imports)]
use crate::state::State;

/// Shared helpers of the group-C conversion programs.
#[path = "groupc_common.rs"]
pub mod groupc_common;

/// Run the `readxyz` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    let (sometext, w) = gc::command_text("readxyz", &ctx.args);
    if w.nwords() < 2 {
        return crate::error::unsupported("readxyz without a file name (interactive file search)");
    }
    let mut readfile = w.word(&sometext, 2).to_string();
    let ilendot = indexlastdot(&readfile) - 1;
    if ilendot > 0 {
        readfile = readfile[..ilendot as usize].to_string();
    }
    let mut fnameato = gc::get_full_filename(&ctx.dir, &readfile);
    fnameato.push_str(".xyz");
    st.mol.ntype = 0;
    st.mol.nmoltype = 0;
    st.mol.sizefactor = 1.0;
    let (rho, ierr) = readxyz(ctx, &mut st, &fnameato)?;
    if rho == 0.0 || ierr != 0 {
        st.log(&fmt("(/'readxyz> Problems with .xyz file ',a/'Routine will exit now'/)", &[FV::Ch(&fnameato, 256)]));
        return Ok(());
    }
    for it in 1..=st.mol.ntype {
        let s = st.mol.atsymbol[it].clone();
        if s == "XXX" || s.is_empty() {
            st.log(&fmt("(/'readxyz> Specify the chemical symbol for atom type ',a3,': ',$)", &[FV::Ch(&st.mol.type_[it], 3)]));
            let v = ctx.stdin.list().word()?.unwrap_or_default();
            st.mol.atsymbol[it] = fstr(&v, 3);
        }
    }
    for im in 1..=st.mol.nmol {
        gc::centreofmass(&mut st, im);
    }
    let ilendot = indexlastdot(&fnameato) - 1;
    let fnameato = format!("{}.ato", &fnameato[..ilendot.max(0) as usize]);
    gc::moloutc(&mut st, &fnameato, None)
}

/// `subroutine readxyz(readfile,rho,ierr)`.  Returns `(rho, ierr)`.
fn readxyz(ctx: &mut ToolContext, st: &mut State, readfile: &str) -> Result<(f32, i32)> {
    let mut f = match TextFile::open(std::path::Path::new(readfile)) {
        Ok(f) => f,
        // rho is undefined in the Fortran here; ierr != 0 makes the caller stop
        Err(_) => return Ok((0.0, 1)),
    };
    let mut ierr: i32 = 0;
    let mut iline = 0;
    let mut natomexpect: i32 = 0;
    {
        let m = &mut st.mol;
        let l = &mut st.lat;
        m.sizefactor = 1.0;
        l.cubic = false;
        m.nmol = 0;
        m.simtk = 300.0;
        l.phib = 90.0;
        l.thetac = 0.0;
        l.phic = 0.0;
        m.tol = 0.0;
        m.tolout = 0.0;
        m.stepmi = 0.0;
        m.stepri = 0.0;
        m.stepr = 0.0;
        m.stepm = 0.0;
        m.tempangreq = 1.0;
        m.tempdihreq = 10.0;
        m.natom = 0;
        m.natomr = 0;
        m.ngrrot = 0;
        m.nrottot = 0;
        m.ntype = 0;
    }
    let mut ia: i32 = 0;
    let mut im: i32 = 0;
    const MATOM: i32 = 150000;
    let mmol = crate::state::MMOL as i32;
    while ierr == 0 && ia < MATOM && im < mmol {
        let sometext = match f.read_line() {
            Ok(s) => s,
            Err(_) => {
                ierr = -1;
                continue;
            }
        };
        let w = parse(&sometext, 11);
        iline += 1;
        if iline == 1 {
            natomexpect = ListRead::from_str(&sometext).int()?.ok_or_else(|| EpsrError::Input("readxyz: cannot read the number of atoms".into()))?;
        } else if iline == 2 {
            let (mut alat, mut blat, mut clat) = (st.lat.alat, st.lat.blat, st.lat.clat);
            let (mut alpha, mut beta, mut gam) = (0f32, 0f32, 0f32);
            if w.nwords() < 7 {
                ierr = 1;
            } else {
                let s = w.span(&sometext, 2, 7);
                let mut r = ListRead::from_str(&s);
                let res = (|| -> crate::fortran::IoRes<()> {
                    r.get_r(&mut alat)?;
                    r.get_r(&mut blat)?;
                    r.get_r(&mut clat)?;
                    r.get_r(&mut alpha)?;
                    r.get_r(&mut beta)?;
                    r.get_r(&mut gam)?;
                    Ok(())
                })();
                ierr = crate::fortran::listread::iostat(&res);
            }
            while ierr != 0 {
                st.log(" readxyz> Need a,b,c axes and 3 angles (alpha,beta,gamma):");
                let mut r = ctx.stdin.list();
                let res = (|| -> crate::fortran::IoRes<()> {
                    r.get_r(&mut alat)?;
                    r.get_r(&mut blat)?;
                    r.get_r(&mut clat)?;
                    r.get_r(&mut alpha)?;
                    r.get_r(&mut beta)?;
                    r.get_r(&mut gam)?;
                    Ok(())
                })();
                ierr = crate::fortran::listread::iostat(&res);
                if ierr < 0 {
                    // the Fortran loops for ever at end of input
                    return Err(EpsrError::Input("readxyz: end of input while reading the cell (a,b,c,alpha,beta,gamma)".into()));
                }
            }
            st.lat.alat = alat;
            st.lat.blat = blat;
            st.lat.clat = clat;
            let (mut phib, mut thetac, mut phic) = (st.lat.phib, st.lat.thetac, st.lat.phic);
            let log = st.logger();
            gc::convertalphabetagamma(alpha, beta, gam, &mut phib, &mut thetac, &mut phic, &log);
            st.lat.phib = phib;
            st.lat.thetac = thetac;
            st.lat.phic = phic;
            gc::setup_lattice(st);
        } else if w.nwords() < 4 {
            ierr = 1;
        } else {
            ia += 1;
            im += 1;
            st.mol.ensure_atoms(ia as usize);
            {
                let s = w.span(&sometext, 1, 4);
                let m = &mut st.mol;
                let mut r = ListRead::from_str(&s);
                let res = (|| -> crate::fortran::IoRes<()> {
                    if let Some(v) = r.word()? {
                        m.alabel[ia] = fstr(&v, 3);
                    }
                    r.get_r(&mut m.xadd[im])?;
                    r.get_r(&mut m.yadd[ia])?;
                    r.get_r(&mut m.zadd[im])?;
                    Ok(())
                })();
                ierr = crate::fortran::listread::iostat(&res);
                m.imageatom[ia] = 0;
                m.xa[ia] = 0.0;
                m.ya[ia] = 0.0;
                m.za[ia] = 0.0;
            }
            st.ftype(ia);
            let m = &mut st.mol;
            if m.qtype[ia] != 1 {
                m.natomr += 1;
            }
            m.natom += 1;
            m.nbond[ia] = 0;
            m.naf[im] = ia;
            m.nal[im] = ia;
            m.nato[im] = 1;
            m.ngrot[im] = 0;
            m.ngrrot += m.ngrot[im];
            m.tethered[im] = false;
            m.phix[im] = 0.0;
            m.phiy[im] = 0.0;
            m.phiz[im] = 0.0;
        }
    }
    let m = &mut st.mol;
    m.nmol = im;
    let rho = m.natomr as f32 / m.volume;
    {
        let mut w = ListWriter::new();
        w.chars(" readxyz> rho = ").real(rho).int(m.natomr).int(m.natom).int(natomexpect);
        st.log(&w.finish());
    }
    let m = &mut st.mol;
    if m.natom >= natomexpect {
        ierr = 0;
    }
    let natom = m.natom;
    let mut iatomicno = vec![0i32; natom.max(0) as usize + 1];
    for ia in 1..=natom {
        let mut ir = 0;
        let mut found = false;
        while !found && ir < gc::MELEMENTS {
            ir += 1;
            found = st.mol.alabel[ia] == gc::chsymbol(ir);
            if found {
                iatomicno[ia as usize] = ir;
            }
        }
        if !found {
            st.log(" readxyz> Cannot read this file: atom label not recognised.");
            return Ok((rho, 1));
        }
    }
    {
        let m = &mut st.mol;
        m.ntypefile = 0;
        let mut ia = 0;
        while ia < m.natom {
            ia += 1;
            let mut itypefile = 0;
            let mut ifound = 0;
            while itypefile < m.ntypefile && ifound == 0 {
                itypefile += 1;
                if m.alabel[ia] == m.typefile[itypefile] {
                    ifound = itypefile;
                }
            }
            if ifound == 0 {
                m.ntypefile += 1;
                let n = m.ntypefile;
                m.typefile[n] = m.alabel[ia].clone();
                m.iatomicnotype[n] = iatomicno[ia as usize];
            }
        }
    }
    set_type_parameters(st);
    let m = &mut st.mol;
    let mut im = 0;
    m.nmoltype = 0;
    while im < m.nmol {
        im += 1;
        let mut imoltype = 0;
        let mut ifound = 0;
        while imoltype < m.nmoltype && ifound == 0 {
            imoltype += 1;
            if m.alabel[m.naf[im]] == m.moltype[imoltype] {
                ifound = imoltype;
                m.mtype[im] = ifound;
            }
        }
        if ifound == 0 {
            m.nmoltype += 1;
            let nt = m.nmoltype;
            m.moltype[nt] = m.alabel[m.naf[im]].clone();
            m.moltethered[nt] = m.tethered[im];
            m.moltemplate[nt] = fmt("(\"moltype\",i2.2)", &[FV::I(nt as i64)]);
            m.mtype[im] = nt;
            m.interntrans[nt] = m.stepmi;
            m.internrot[nt] = m.stepri;
            m.wholerot[nt] = m.stepr;
            m.wholetrans[nt] = m.stepm;
        }
    }
    st.log(&fmt("(/'readxyz> Following molecule types defined initially:-'/)", &[]));
    for imoltype in 1..=st.mol.nmoltype {
        let m = &st.mol;
        st.log(&fmt(
            "(1x,i3,1x,a3,1x,a10,2(1x,e12.5))",
            &[
                FV::I(imoltype as i64),
                FV::Ch(&m.moltype[imoltype], 3),
                FV::Ch(&m.moltemplate[imoltype], 10),
                FV::R(m.wholerot[imoltype]),
                FV::R(m.wholetrans[imoltype]),
            ],
        ));
    }
    Ok((rho, ierr))
}

/// The `ctype` calls of `readxyz`/`readpdb`: default potential parameters
/// from the element of each atom type in the file.
pub(crate) fn set_type_parameters(st: &mut State) {
    for it in 1..=st.mol.ntypefile {
        let iat = st.mol.iatomicnotype[it];
        let sym = gc::chsymbol(iat);
        let (epsilon, sigma, amass): (f32, f32, f32) = if sym == "H" {
            (0.0, 0.0, 2.0)
        } else if sym == "O" {
            (0.65, 3.1, gc::chatwt(iat))
        } else if sym == "N" {
            (0.7, 3.2, gc::chatwt(iat))
        } else if sym == "C" {
            (0.8, 3.7, gc::chatwt(iat))
        } else {
            (0.8, 3.2, gc::chatwt(iat))
        };
        let a1 = st.mol.typefile[it].clone();
        // ctype(typefile,chsymbol,0,epsilon,sigma,amass,1.0,0.0,0.5*sigma,1.0,0.0):
        // qt=1, rqt=0, natomsph=0.5*sigma, rhosph=1, radminsph=0; radmaxsph and
        // pseudofacsph absent -> 0.5*sigcor(n1) (= 0.5*sigma) and 1.0.
        st.ctype(&a1, sym, 0, epsilon, sigma, amass, 1.0, 0.0, 0.5 * sigma, 1.0, 0.0, 0.5 * sigma, 1.0);
    }
}
