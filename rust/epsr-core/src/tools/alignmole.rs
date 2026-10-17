//! `alignmole` — translation of `fortran_src/alignmole.f` (`alignmole`,
//! `displacement`).  A standalone program (not an EPSR shell command).
//!
//! STATUS: translated; verified byte-identical against the Fortran for runs
//! in which no trial rotation is rejected (test cases EPSR-test-ex47/48), see
//! the note below.
//!
//! Invocation (same as the Fortran program):
//! `epsr-rs <folder> alignmole <file>.ato [x|y|z [niter [all | atom numbers...]]]`
//! (missing values are read from standard input).  Every molecule is rotated
//! by Monte Carlo (`rotatm`) so as to minimise the mean square distance of
//! its atoms (all, or the listed atom numbers) from the chosen axis; the
//! `.ato` file is rewritten.
//!
//! NOTE on the Fortran: when a trial rotation is rejected the coordinates
//! are "restored" from the work arrays `xr, yr, zr`, which are allocated but
//! never assigned; the Fortran therefore copies uninitialised heap memory
//! into the molecule (in practice zeros and leftovers of earlier
//! allocations, giving `**************` in the `.ato` file).  This cannot
//! be reproduced; the translation uses zero-filled work arrays (what a fresh
//! allocation would contain) and prints a warning when this happens.
//! Also, `niter < 10` makes the Fortran divide by zero (`naccum = niter/10`);
//! here this is an input error.

use super::readxyz::groupc_common as gc;
use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::ListRead;
use crate::fortran::text::{findex, parse};
#[allow(unused_imports)]
use crate::state::State;

/// `function displacement(im,alignaxis,useall,axisatom,naxisatom)`.
fn displacement(st: &State, im: i32, alignaxis: char, useall: bool, axisatom: &[i32], naxisatom: i32) -> f32 {
    let m = &st.mol;
    let mut sum: f32 = 0.0;
    let terms = |ia: i32, sum: &mut f32| {
        let (a, b) = match alignaxis {
            'x' => (m.ya[ia], m.za[ia]),
            'y' => (m.xa[ia], m.za[ia]),
            _ => (m.xa[ia], m.ya[ia]),
        };
        *sum += a * a;
        *sum += b * b;
    };
    if useall || naxisatom <= 0 {
        for ia in m.naf[im]..=m.nal[im] {
            terms(ia, &mut sum);
        }
        sum / m.nato[im] as f32
    } else {
        for ir in 1..=naxisatom {
            let ia = axisatom[(ir - 1) as usize] + m.naf[im] - 1;
            terms(ia, &mut sum);
        }
        sum / naxisatom as f32
    }
}

/// Run the `alignmole` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    let (mut sometext, mut w) = gc::command_text("alignmole", &ctx.args);
    if w.nwords() < 2 {
        return crate::error::unsupported("alignmole without a file name (interactive file search)");
    }
    let mut readfile = w.word(&sometext, 2).to_string();
    let alignaxis: char = if w.nwords() > 2 {
        w.word(&sometext, 3).chars().next().unwrap_or(' ')
    } else {
        st.log("\nalignmole> Type axis to align to (x, y, or z): ");
        ctx.stdin.read_line()?.chars().next().unwrap_or(' ')
    };
    let niter: i32 = if w.nwords() > 3 {
        ListRead::from_str(w.word(&sometext, 4)).int()?.ok_or_else(|| EpsrError::Input("alignmole: bad number of iterations".into()))?
    } else {
        st.log("\nalignmole> Give number of iterations: ");
        ctx.stdin.list().int()?.ok_or_else(|| EpsrError::Input("alignmole: bad number of iterations".into()))?
    };
    let mut naxisatom: i32 = 0;
    let mut nf = 5;
    let mut useall = true;
    if w.nwords() >= nf {
        useall = findex(w.word(&sometext, nf), "all") > 0;
    }
    let mut axisatom: Option<Vec<i32>> = None;
    while axisatom.is_none() && !useall {
        naxisatom = w.nwords() - nf + 1;
        if naxisatom >= 2 {
            let s = w.span(&sometext, nf, w.nwords());
            let mut r = ListRead::from_str(&s);
            let mut v = vec![0i32; naxisatom as usize];
            let res = (|| -> crate::fortran::IoRes<()> {
                for x in v.iter_mut() {
                    r.get_i(x)?;
                }
                Ok(())
            })();
            if res.is_ok() {
                axisatom = Some(v);
            }
        }
        if axisatom.is_none() {
            st.log("\nalignmole> Give list of atom numbers to define molecular axis\n (Type \"all\" to use all atoms in each molecule) : ");
            sometext = ctx.stdin.read_line()?;
            useall = findex(&sometext, "all") > 0;
            w = parse(&sometext, 10);
            nf = 1;
        }
    }
    let ilendot = findex(&readfile, ".") - 1;
    if ilendot > 0 {
        readfile = readfile[..ilendot as usize].to_string();
    }
    let mut fnameato = gc::get_full_filename(&ctx.dir, &readfile);
    fnameato.push_str(".ato");
    st.mol.ntype = 0;
    st.mol.nmoltype = 0;
    st.mol.sizefactor = 1.0;
    st.minitc(&fnameato)?;
    if st.dist.rho == 0.0 {
        st.log(&fmt("(/'Problems with .ato file ',a/'Routine will exit now'/)", &[FV::Ch(&fnameato, 256)]));
        return Ok(());
    }
    {
        let mut wl = ListWriter::new();
        wl.real(st.mol.dmax).real(st.mol.dmaxmol);
        st.log(&wl.finish());
    }
    if let Some(ax) = &axisatom {
        let mut ierr = 0;
        let mut im = 0;
        while im < st.mol.nmol && ierr == 0 {
            im += 1;
            let mut ia = 0;
            while ia < naxisatom && ierr == 0 {
                ia += 1;
                let a = ax[(ia - 1) as usize];
                if a < 1 || a > st.mol.nato[im] {
                    ierr = 1;
                }
            }
        }
        if ierr != 0 {
            axisatom = None;
            naxisatom = 0;
            useall = true;
        }
    }
    if naxisatom > 0 {
        let mut wl = ListWriter::new();
        wl.int(naxisatom);
        if let Some(ax) = &axisatom {
            for a in ax {
                wl.int(*a);
            }
        }
        st.log(&wl.finish());
    }
    let ax: Vec<i32> = axisatom.clone().unwrap_or_default();
    let steprs = st.mol.stepr;
    if st.mol.stepr <= 0.0 {
        st.mol.stepr = 0.1;
    }
    // xr, yr, zr: allocated, never assigned (see the module documentation)
    let mut xr: Vec<f32> = Vec::new();
    let mut yr: Vec<f32> = Vec::new();
    let mut zr: Vec<f32> = Vec::new();
    let mut warned = false;
    for im in 1..=st.mol.nmol {
        let nato = st.mol.nato[im] as usize;
        if xr.len() < nato {
            xr = vec![0.0; nato];
            yr = vec![0.0; nato];
            zr = vec![0.0; nato];
        }
        let mut ntrials = 0;
        let mut nreject = 0;
        let naccum = niter / 10;
        let mut iter = 0;
        let mut energold = displacement(&st, im, alignaxis, useall, &ax, naxisatom);
        st.log(&fmt("(/'alignmole> Molecule: ',i4,' Initial energy: ',e12.5)", &[FV::I(im as i64), FV::R(energold)]));
        while iter < niter {
            iter += 1;
            ntrials += 1;
            let stepr = st.mol.stepr;
            gc::rotatm(&mut st, im, stepr);
            let energnew = displacement(&st, im, alignaxis, useall, &ax, naxisatom);
            if energnew > energold {
                nreject += 1;
                if !warned {
                    eprintln!("alignmole: WARNING: a rotation was rejected; the Fortran restores the molecule from never-assigned work arrays (undefined result); zeros are used here");
                    warned = true;
                }
                let m = &mut st.mol;
                let mut nat = m.naf[im] - 1;
                for ia in 1..=m.nato[im] {
                    nat += 1;
                    m.xa[nat] = xr[(ia - 1) as usize];
                    m.ya[nat] = yr[(ia - 1) as usize];
                    m.za[nat] = zr[(ia - 1) as usize];
                }
            } else {
                energold = energnew;
            }
            if naccum == 0 {
                return Err(EpsrError::Input("alignmole: the number of iterations must be at least 10 (the Fortran divides by niter/10)".into()));
            }
            if naccum * (iter / naccum) == iter {
                let mut rat = nreject as f32 / ntrials as f32;
                st.log(&fmt(
                    "(/'alignmole> Iteration ',i7,' Energy: ',e12.5,' Reject rate: ',f10.5)",
                    &[FV::I(iter as i64), FV::R(energold), FV::R(rat)],
                ));
                rat /= 0.75;
                if rat > 1.1 {
                    rat = 1.1;
                }
                if rat < 0.9 {
                    rat = 0.9;
                }
                st.mol.stepr /= rat;
            }
        }
    }
    st.mol.stepr = steprs;
    gc::moloutc(&mut st, &fnameato, None)
}
