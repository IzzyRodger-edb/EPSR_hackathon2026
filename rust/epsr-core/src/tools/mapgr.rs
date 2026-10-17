//! `mapgr` — translation of `fortran_src/mapgrquad.f` (+ `plotsharm_routines.f`).
//!
//! Reconstructs g(r, ω) along specified directions / relative orientations
//! from the spherical-harmonic coefficients written by `sharm`/`sdf`.
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> mapgr <input file>`
//! Keyword table: `crate::symtab_tools::MAPGR`.  Output `<run>.mapgr.g01`:
//! column 1 the centres g(r), then one column per view.

use super::sharm::sharm_routines::Coeffs;
use super::ToolContext;
use crate::ato::cosd;
use crate::error::Result;
use crate::fortran::format::{fmt, FV};
use crate::fortran::listread::{read_int, ListRead};
use crate::fortran::text::{fstr, indexlastdot};
use crate::fortran::F2;
use crate::output::w_int;
use crate::state::State;

/// Run the `mapgr` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    if !super::sharm::open(ctx, &mut st, &crate::symtab_tools::MAPGR)? {
        return Ok(());
    }
    let log = st.logger();
    let c = Coeffs::read(&st)?;
    let nr = c.nr;
    let nexpan = c.nexpan;
    let mut iref = 5;
    let forcemirrorsym = fstr(&st.sym.symbol_value[iref], 2) == "y";
    iref += 1;
    let nviews = read_int(&st.sym.symbol_value[iref]).map_err(|_| crate::error::EpsrError::Input("cannot read nviews".into()))?;
    let ncol = (nviews + 1).max(1) as usize;
    let mut grw: F2<f32> = F2::new(500, ncol);
    let grwerr: F2<f32> = F2::new(500, ncol);
    let mut titledat: Vec<String> = Vec::new();
    let mut ituse = 1;
    for i in 1..=nr {
        grw[(i, ituse)] = 1.0f32 + c.datac[(i, 1)];
    }
    titledat.push(fmt("(8x,a,9x)", &[FV::S("g_centres(r)")]));
    let mut grtype: i32 = 0;
    let mut zc = vec![0.0f32; nexpan.max(0) as usize];
    let mut zs = vec![0.0f32; nexpan.max(0) as usize];
    for _it in 1..=nviews {
        iref += 1;
        let mut itest;
        match read_int(&st.sym.symbol_value[iref]) {
            Ok(v) => {
                grtype = v;
                itest = 0;
            }
            Err(_) => itest = 1,
        }
        let _ = itest;
        iref += 1;
        let mut a = [0.0f32; 5];
        {
            let mut r = ListRead::from_str(&st.sym.symbol_value[iref]);
            itest = 0;
            for x in a.iter_mut() {
                match r.real() {
                    Ok(Some(v)) => *x = v,
                    Ok(None) => {}
                    Err(_) => {
                        itest = 1;
                        break;
                    }
                }
            }
        }
        let [thetal, phil, phim, thetam, chim] = a;
        let (mut l1only, mut l2only, mut lonly, mut n1only, mut n2only, mut m2fix) = (1, 1, 1, 1, 1, 1);
        if grtype < 3 {
            l2only = 0;
            n2only = 0;
            m2fix = 0;
        }
        if grtype < 2 {
            l1only = 0;
            lonly = 0;
            n1only = 0;
        }
        if itest != 0 {
            continue;
        }
        let ctl = cosd(thetal);
        let phi = phim;
        let cz = cosd(thetam);
        let chi = chim;
        for k in 1..=nexpan as usize {
            let (sumc, sums, _, _) = c.angular(k, [l1only, l2only, lonly, n1only, n2only, m2fix], phil, ctl, phi, cz, chi);
            zc[k - 1] = sumc;
            zs[k - 1] = if forcemirrorsym { 0.0 } else { sums };
        }
        ituse += 1;
        for i in 1..=nr {
            let mut sumt: f32 = 0.0;
            for k in 1..=nexpan {
                let term = c.datac[(i, k)] * zc[(k - 1) as usize] + c.datas[(i, k)] * zs[(k - 1) as usize];
                sumt += term;
            }
            grw[(i, ituse)] = sumt + 1.0f32;
        }
        titledat.push(fmt("(10x,'view_',i3.3)", &[FV::I((ituse - 1) as i64)]));
    }
    let fnameinp = st.sym.fnameinp.clone();
    let idot = indexlastdot(&fnameinp);
    let fnameout = if idot > 0 { format!("{}g", &fnameinp[..idot as usize]) } else { format!("{}.g", fnameinp.trim_end()) };
    w_int(nr, &fnameout, ituse, &c.rd, &grw, &grwerr, "r[A]", &titledat, &log)?;
    Ok(())
}

/// Program-specific branch of `check_symbol` (`process_symbols.f`, `symbol_setup.eq.'mapgr'`).
/// Returns the Fortran `ierr` (0 = accepted).
pub fn check_symbol_extra(st: &mut State, icount: i32, _ifilecount: i32, iref: i32) -> Result<i32> {
    // (the 'nvary' branch of check_symbol does not apply: mapgr has no such symbol)
    st.check_symbol_mapgr(icount, iref)
}
