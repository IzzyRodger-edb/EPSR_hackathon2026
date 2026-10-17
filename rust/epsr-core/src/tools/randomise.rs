//! `randomise` — translation of `fortran_src/introtcluster.f` (subroutine `randomise`).
//!
//! STATUS: translated; verified against the Fortran (EPSR-test-ex35).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> randomise [<file>.ato]`
//! (the file name is asked for on standard input if not given).
//!
//! Every (untethered) molecule is given a random position in the box (if
//! `stepm > 0`) and a random whole-molecule rotation (if `stepr > 0`); the
//! random numbers come from the generator state stored in the `.ato` file,
//! which is rewritten in place.

use super::makeato::setup_common::{command_text, get_full_filename, rotatm, search_for_files, strip_first_dot, SysCmds};
use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::text::parse;
use crate::neighbour::minimagec;

/// Run the `randomise` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    let sys = SysCmds::load();
    let (sometext, w) = command_text("randomise", ctx);
    let readfile: String;
    if w.nwords() < 2 {
        let (t, _iaccept) = search_for_files(&st, ctx, &sys, ".ato")?;
        let w2 = parse(&t, 5);
        if w2.nwords() <= 0 {
            st.log("\nrandomise> No .ato file specified! - routine will exit\n");
            return Ok(());
        }
        readfile = w2.word(&t, 1).to_string();
    } else {
        readfile = w.word(&sometext, 2).to_string();
    }
    let readfile = strip_first_dot(&readfile);
    let fnameato = format!("{}.ato", get_full_filename(&sys, &ctx.dir, &readfile));
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
    st.normalise_coordinates();
    let stepris = st.mol.stepri;
    let steprs = st.mol.stepr;
    if stepris > 0.0 {
        st.mol.stepri = 1.0;
    }
    if steprs > 0.0 {
        st.mol.stepr = 1.0;
    }
    for im in 1..=st.mol.nmol {
        if !st.mol.tethered[im] {
            if st.mol.stepm > 0.0 {
                let fa = st.rng.ran1() - 0.5;
                let fb = st.rng.ran1() - 0.5;
                let fc = st.rng.ran1() - 0.5;
                let l = &st.lat;
                let m = &mut st.mol;
                m.zadd[im] = fc * l.clatz;
                m.yadd[im] = fb * l.blaty + fc * l.claty;
                m.xadd[im] = fa * l.alatx + fb * l.blatx + fc * l.clatx;
                let (mut x, mut y, mut z) = (m.xadd[im], m.yadd[im], m.zadd[im]);
                minimagec(&mut x, &mut y, &mut z);
                m.xadd[im] = x;
                m.yadd[im] = y;
                m.zadd[im] = z;
            }
            if st.mol.stepr > 0.0 {
                let s = st.mol.stepr;
                rotatm(&mut st, im, s);
            }
        }
    }
    st.log(&format!(
        "randomise> {:>5} rotational groups in {:>5} molecules ",
        st.mol.ngrrot, st.mol.nmol
    ));
    st.mol.stepri = stepris;
    st.mol.stepr = steprs;
    st.unnormalise_coordinates();
    st.moloutc(&fnameato)
}
