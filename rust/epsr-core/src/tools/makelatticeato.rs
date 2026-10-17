//! `makelatticeato` — translation of `fortran_src/makelatticeato.f`.
//!
//! STATUS: translated; verified against the Fortran (EPSR-test-ex38, ex39).
//!
//! Invocation (same as the Fortran program):
//! `epsr-rs <folder> makelatticeato [<file>.ato [na [nb [nc [<output> [y|n]]]]]]`
//! Missing values are read from standard input.
//!
//! The molecules of an existing `.ato` file are treated as the contents of a
//! unit cell (the box) and replicated `na x nb x nc` times; optionally
//! ("convert") every atom is then made a separate, tethered one-atom molecule.
//!
//! Quirk reproduced: the Fortran parses an uninitialised variable (`text`)
//! instead of the convert answer, so the answer is compared over the
//! character range of the *first word of the previous parse* (normally the
//! output-file name, or the command name when that was given on the command
//! line): an answer `y`/`Y` starting in column 1 is accepted.

use super::makeato::setup_common::{command_text, eof_err, get_full_filename, read_text, search_for_files, SysCmds};
use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, FV};
use crate::fortran::listread::read_int;
use crate::fortran::text::{findex, parse, substr, Words};
use crate::state::MMOL;

/// Run the `makelatticeato` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    let sys = SysCmds::load();
    let (mut sometext, mut w) = command_text("makelatticeato", ctx);
    // (ncf(1), ncl(1)) of the latest parse that found a word
    let mut first: (i32, i32) = (w.ncf[0], w.ncl[0]);
    let note = |w: &Words, first: &mut (i32, i32)| {
        if w.nwords() > 0 {
            *first = (w.ncf[0], w.ncl[0]);
        }
    };
    let readfile: String;
    if w.nwords() < 2 {
        loop {
            let (t, iaccept) = search_for_files(&st, ctx, &sys, ".ato")?;
            if iaccept != 'y' {
                st.log("\nmakelatticeato> No .ato files found in the current directory - try another directory.\n");
                return Ok(());
            }
            sometext = t;
            w = parse(&sometext, 5);
            note(&w, &mut first);
            if w.nwords() > 0 {
                break;
            }
            st.log("\nmakelaticeato> No file specified! - try again!\n");
        }
        readfile = w.span(&sometext, 1, w.nwords());
    } else {
        readfile = w.word(&sometext, 2).to_string();
    }
    let rd5 = |ctx: &mut ToolContext, n: usize| -> Result<Vec<i32>> {
        let mut v = vec![0i32; n];
        let mut r = ctx.stdin.list();
        for x in v.iter_mut() {
            r.get_i(x).map_err(eof_err)?;
        }
        Ok(v)
    };
    let (na, nb, nc);
    if w.nwords() > 2 {
        na = read_int(w.word(&sometext, 3))?;
        if w.nwords() > 3 {
            nb = read_int(w.word(&sometext, 4))?;
            if w.nwords() > 4 {
                nc = read_int(w.word(&sometext, 5))?;
            } else {
                st.log("\nmakelattice> Give number of unit cells along c-axis: ");
                nc = rd5(ctx, 1)?[0];
            }
        } else {
            st.log("\nmakelattice> Give number of unit cells along b and c-axes: ");
            let v = rd5(ctx, 2)?;
            nb = v[0];
            nc = v[1];
        }
    } else {
        st.log("\nmakelattice> Give number of unit cells along a, b and c-axes: ");
        let v = rd5(ctx, 3)?;
        na = v[0];
        nb = v[1];
        nc = v[2];
    }
    let mut fnameato = get_full_filename(&sys, &ctx.dir, &readfile);
    if findex(&fnameato, ".ato") == 0 {
        fnameato.push_str(".ato");
    }
    // output file name
    let outfile: String;
    if w.nwords() > 5 {
        outfile = w.word(&sometext, 6).to_string();
    } else {
        loop {
            st.log("\nmakelatticeato> Type name of file to put lattice in: ");
            sometext = read_text(ctx)?;
            w = parse(&sometext, 5);
            note(&w, &mut first);
            if w.nwords() > 0 {
                break;
            }
            st.log("\nmakelaticeato> No file specified! - try again!\n");
        }
        outfile = w.span(&sometext, 1, w.nwords());
    }
    // convert molecules to atoms?
    let answer: String = if w.nwords() > 6 {
        w.word(&sometext, 7).to_string()
    } else {
        st.log("makelatticeato> Convert molecules to atoms? (y or n) ");
        read_text(ctx)?
    };
    let cmp = substr(&answer, first.0, first.1);
    let cmp = cmp.trim_end();
    let convert = cmp == "y" || cmp == "Y";
    st.mol.ntype = 0;
    st.mol.nmoltype = 0;
    st.mol.sizefactor = 1.0;
    match st.minitc(&fnameato) {
        Ok(()) => {}
        Err(EpsrError::Io(_)) => {
            st.log(&format!("\nProblems with .ato file {}\nRoutine will exit now\n", fnameato));
            return Ok(());
        }
        Err(e) => return Err(e),
    }
    if st.dist.rho == 0.0 {
        st.log(&format!("\nProblems with .ato file {}\nRoutine will exit now\n", fnameato));
        return Ok(());
    }
    let nunitcell = st.mol.nmol;
    let nu = nunitcell.max(0) as usize;
    let mut xadds = vec![0f32; nu + 1];
    let mut yadds = vec![0f32; nu + 1];
    let mut zadds = vec![0f32; nu + 1];
    let mut phixs = vec![0f32; nu + 1];
    let mut phiys = vec![0f32; nu + 1];
    let mut phizs = vec![0f32; nu + 1];
    {
        let m = &st.mol;
        for im in 1..=nunitcell {
            let k = im as usize;
            if m.tethered[im] {
                xadds[k] = m.phix[im];
                yadds[k] = m.phiy[im];
                zadds[k] = m.phiz[im];
            } else {
                xadds[k] = m.xadd[im];
                yadds[k] = m.yadd[im];
                zadds[k] = m.zadd[im];
            }
            phixs[k] = m.phix[im];
            phiys[k] = m.phiy[im];
            phizs[k] = m.phiz[im];
        }
    }
    let mut fnameout = get_full_filename(&sys, &ctx.dir, &outfile);
    if findex(&fnameout, ".ato") == 0 {
        fnameout.push_str(".ato");
    }
    st.log(&format!("makelatticeato> Output being written to file: {}", fnameout));
    let ncells = na.max(0) as i64 * nb.max(0) as i64 * nc.max(0) as i64;
    let natom_unit = st.mol.natom;
    if (ncells * nunitcell as i64) as usize > MMOL {
        return Err(EpsrError::Input(format!(
            "makelatticeato: too many molecules ({})",
            ncells * nunitcell as i64
        )));
    }
    let m = &mut st.mol;
    m.ensure_atoms((ncells * natom_unit as i64).max(natom_unit as i64).max(1) as usize);
    let ngrot_unit = m.ngrrot;
    let nrot_unit = m.nrottot;
    m.ensure_rotgrps((ncells * ngrot_unit as i64).max(1) as usize);
    m.ensure_labrot((ncells * nrot_unit as i64).max(1) as usize);
    let mut ia: i32 = 0;
    let mut ig: i32 = 0;
    m.ngrrot = 0;
    m.nrottot = 0;
    let mut im: i32 = 0;
    let (mut xav, mut yav, mut zav) = (0f32, 0f32, 0f32);
    let mut nadd: i32 = 0;
    let acent = 0.5 * na as f32;
    let bcent = 0.5 * nb as f32;
    let ccent = 0.5 * nc as f32;
    let l = &st.lat;
    for kc in 1..=nc {
        let fc = (kc - 1) as f32 - ccent;
        let zlat = fc * l.clatz;
        for kb in 1..=nb {
            let fb = (kb - 1) as f32 - bcent;
            let ylat = fb * l.blaty + fc * l.claty;
            for ka in 1..=na {
                let fa = (ka - 1) as f32 - acent;
                let xlat = fa * l.alatx + fb * l.blatx + fc * l.clatx;
                let mut ja: i32 = 0;
                let mut irgref: i32 = 0;
                for jm in 1..=nunitcell {
                    let k = jm as usize;
                    im += 1;
                    m.nato[im] = m.nato[jm];
                    m.xadd[im] = xlat + xadds[k];
                    m.yadd[im] = ylat + yadds[k];
                    m.zadd[im] = zlat + zadds[k];
                    m.tethered[im] = m.tethered[jm];
                    m.mtype[im] = m.mtype[jm];
                    if m.tethered[im] {
                        m.phix[im] = xlat + phixs[k];
                        m.phiy[im] = ylat + phiys[k];
                        m.phiz[im] = zlat + phizs[k];
                    } else {
                        m.phix[im] = phixs[k];
                        m.phiy[im] = phiys[k];
                        m.phiz[im] = phizs[k];
                    }
                    xav += m.xadd[im];
                    yav += m.yadd[im];
                    zav += m.zadd[im];
                    m.naf[im] = ia + 1;
                    for _jb in 1..=m.nato[jm] {
                        ja += 1;
                        ia += 1;
                        m.alabel[ia] = m.alabel[ja].clone();
                        m.imageatom[ia] = if m.imageatom[ja] > 0 { m.imageatom[ja] + nadd } else { 0 };
                        m.nbond[ia] = m.nbond[ja];
                        m.xa[ia] = m.xa[ja];
                        m.ya[ia] = m.ya[ja];
                        m.za[ia] = m.za[ja];
                        m.itype[ia] = m.itype[ja];
                        m.mref[ia] = im;
                        for ib in 1..=m.nbond[ia] {
                            m.iabond[(ia, ib)] = m.iabond[(ja, ib)] + nadd;
                            m.tbond[(ia, ib)] = m.tbond[(ja, ib)];
                        }
                    }
                    m.ngrot[im] = m.ngrot[jm];
                    if m.ngrot[im] > 0 {
                        m.ngrrot += m.ngrot[im];
                        for _ih in 1..=m.ngrot[im] {
                            ig += 1;
                            irgref += 1;
                            m.mrot[ig] = im;
                            m.rotgrp[ig] = m.rotgrp[irgref].clone();
                            if m.rotgrp[ig] == "DIH" {
                                return crate::error::unsupported("makelatticeato: DIH rotational groups");
                            }
                            for j in 1..=2 {
                                m.labaxi[(ig, j)] = m.labaxi[(irgref, j)] + nadd;
                            }
                            m.nrot[ig] = m.nrot[irgref];
                            m.rotf[ig] = m.nrottot + 1;
                            m.nrottot += m.nrot[ig];
                            m.rotl[ig] = m.nrottot;
                            let mut labrefold = m.rotf[irgref] - 1;
                            let mut labrefnew = m.rotf[ig] - 1;
                            for _j in 1..=m.nrot[ig] {
                                labrefold += 1;
                                labrefnew += 1;
                                m.labrot[labrefnew] = m.labrot[labrefold] + nadd;
                            }
                        }
                    }
                    m.nal[im] = ia;
                }
                nadd = ia;
            }
        }
    }
    m.natom = ia;
    m.nmol = im;
    let xbox = xav / m.nmol as f32;
    let ybox = yav / m.nmol as f32;
    let zbox = zav / m.nmol as f32;
    for im in 1..=m.nmol {
        m.xadd[im] -= xbox;
        m.yadd[im] -= ybox;
        m.zadd[im] -= zbox;
        if m.tethered[im] {
            m.phix[im] -= xbox;
            m.phiy[im] -= ybox;
            m.phiz[im] -= zbox;
        }
    }
    m.natom = ia;
    if ctx.verbose {
        println!("makelatticeato> {:>8} molecules and {:>8} atoms created", m.nmol, m.natom);
    }
    let l = &mut st.lat;
    l.alat *= na as f32;
    l.blat *= nb as f32;
    l.clat *= nc as f32;
    l.cubic = l.cubic && (na == nb && na == nc && nb == nc);
    if l.cubic {
        m.dmax *= na as f32;
    }
    if convert {
        for im in 1..=m.nmol {
            for ia in m.naf[im]..=m.nal[im] {
                m.xa[ia] += m.xadd[im];
                m.ya[ia] += m.yadd[im];
                m.za[ia] += m.zadd[im];
            }
        }
        m.nmol = m.natom;
        for im in 1..=m.nmol {
            let ia = im;
            m.xadd[im] = m.xa[ia];
            m.phix[im] = m.xa[ia];
            m.yadd[im] = m.ya[ia];
            m.phiy[im] = m.ya[ia];
            m.zadd[im] = m.za[ia];
            m.phiz[im] = m.za[ia];
            m.tethered[im] = true;
            m.tetheredatom[im] = 0;
            m.nstepreject[im] = 0;
            m.nato[im] = 1;
            m.naf[im] = ia;
            m.nal[im] = ia;
            m.nbond[im] = 0;
            m.ngrot[im] = 0;
            m.xa[ia] = 0.0;
            m.ya[ia] = 0.0;
            m.za[ia] = 0.0;
            m.imageatom[ia] = 0;
        }
        m.tolout = 0.1;
        m.nmoltype = 0;
        let mut im = 0;
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
                m.moltetheredatom[nt] = m.tetheredatom[im];
                m.moltemplate[nt] = fmt("('moltype',i2.2)", &[FV::I(nt as i64)]);
                m.mtype[im] = nt;
                m.wholerot[nt] = m.stepr;
                m.wholetrans[nt] = m.stepm;
            }
        }
        // reallocate1d_r(rfractmol,nmoltype); rfractmol=0.0
        m.rfractmol = crate::fortran::F1::new(m.nmoltype.max(1) as usize);
    }
    st.moloutc(&fnameout)
}
