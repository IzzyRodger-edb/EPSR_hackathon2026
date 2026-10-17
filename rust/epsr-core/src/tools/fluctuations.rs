//! `fluctuations` — translation of `fortran_src/fluctuations.f`
//! (`run_fluctuations`, `find_fluctuations`).
//!
//! Number fluctuations of molecule centres in spheres of radius `rmax` placed
//! at `nsamples` random positions in the box: the distribution of
//! N - <N> (accumulated in `<run>.FLUCTUATIONS.n01`) and the running averages
//! of N and N^2 together with S(0) = <N^2>/<N> - <N>, which are stored back
//! in the input file (`ntotal`, `nsqtotal`, `szero`).  The random-number
//! state is kept in `randomno.dat`.
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> fluctuations <input file>`
//! Keyword table: `crate::symtab_tools::FLUCTUATIONS`.

use super::partials::anal_common::*;
use super::ToolContext;
use crate::error::Result;
use crate::fortran::fmath::nint;
use crate::fortran::format::{fmt, FV};
use crate::fortran::{F1, F2};
use crate::neighbour::minimagec;
#[allow(unused_imports)]
use crate::state::State;

const MDIST: usize = crate::state::MFILE;
const MZ: usize = 5000;

/// Run the `fluctuations` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    if !ctx.open_input(&mut st, &crate::symtab_tools::FLUCTUATIONS)? {
        return Ok(());
    }
    let mut iref = 2;
    let nz = sv_int(&st, iref)?.min(MZ as i32);
    iref += 1;
    let mut nsumts = sv_int(&st, iref)?;
    let mut nsumt = if nsumts < 0 { 0 } else { nsumts };
    iref += 1;
    let npix = sv_int(&st, iref)?;
    let mut ndist = st.sym.nfiles;
    if ndist > MDIST as i32 {
        ndist = MDIST as i32;
    }
    let nd = ndist as usize + 1;
    let mut damax = vec![0.0f32; nd];
    let mut sumnumber = vec![0.0f64; nd];
    let mut sumnumbersqr = vec![0.0f64; nd];
    let mut szero = vec![0.0f64; nd];
    for id in 1..=ndist {
        let i = id as usize;
        let mut iref = st.sym.ireff(st.sym.nsymboluptodata + 1, id);
        damax[i] = sv_real(&st, iref)?;
        iref += 1;
        sumnumber[i] = sv_dble(&st, iref)?;
        iref += 1;
        sumnumbersqr[i] = sv_dble(&st, iref)?;
        iref += 1;
        szero[i] = sv_dble(&st, iref)?;
        if nsumt == 0 {
            sumnumber[i] = 0.0;
            sumnumbersqr[i] = 0.0;
        }
    }
    let rho = st.mol.nmol as f32 / (st.lat.alatx * st.lat.blaty * st.lat.clatz);
    let zstep = 1.0f32;
    let mut z: F1<f32> = F1::new(MZ);
    for j in 1..=nz {
        z[j] = (j - 1 - nz / 2) as f32 * zstep;
    }
    let mut gofz: F2<f32> = F2::new(MZ, MDIST);
    st.normalise_coordinates();
    initran1_file(ctx, &mut st)?;
    let mut localnumber = vec![0.0f64; nd];
    let mut localnumbersqr = vec![0.0f64; nd];
    find_fluctuations(&mut st, nz, ndist, npix, zstep, &damax, rho, &mut localnumber, &mut localnumbersqr, &mut gofz);

    let titler = "Density";
    let mut titledat = Vec::new();
    for id in 1..=ndist as usize {
        let text = fmt("(e12.5)", &[FV::R(damax[id])]);
        let text = text.trim_start();
        titledat.push(crate::fortran::fstr(&format!("r_max:{}", upto_blank(text)), 29));
    }
    let mut gzt: F2<f32> = F2::new(MZ, MDIST);
    let mut gzt2: F2<f32> = F2::new(MZ, MDIST);
    outgofz(&st, "n", nz, ndist, &mut nsumt, &mut z, &gofz, &mut gzt, &mut gzt2, titler, &titledat)?;

    let frac = (1.0f32 / nsumt as f32) as f64;
    let frac1 = 1.0 - frac;
    for id in 1..=ndist as usize {
        sumnumber[id] = sumnumber[id] * frac1 + localnumber[id] * frac;
        sumnumbersqr[id] = sumnumbersqr[id] * frac1 + localnumbersqr[id] * frac;
        szero[id] = if sumnumber[id] > 0.0 { sumnumbersqr[id] / sumnumber[id] - sumnumber[id] } else { 0.0 };
    }

    let icount = st.sym.search_for_symbol("nsize", 0);
    let mut iref = st.sym.ireff(icount, 1);
    set_sv_int(&mut st, iref, nz);
    iref += 1;
    if nsumts >= 0 {
        nsumts = nsumt;
    }
    set_sv_int(&mut st, iref, nsumts);
    iref += 1;
    set_sv_int(&mut st, iref, npix);
    iref += 1;
    set_sv_int(&mut st, iref, ndist);
    for id in 1..=ndist as usize {
        iref += 1;
        set_sv_real(&mut st, iref, damax[id]);
        iref += 1;
        set_sv_dble(&mut st, iref, sumnumber[id]);
        iref += 1;
        set_sv_dble(&mut st, iref, sumnumbersqr[id]);
        iref += 1;
        set_sv_dble(&mut st, iref, szero[id]);
    }
    st.save_inp_file()?;
    saveran1_file(ctx, &st)
}

/// `subroutine find_fluctuations(rho,sumnumber,sumnumbersqr)`.
#[allow(clippy::too_many_arguments)]
fn find_fluctuations(
    st: &mut State,
    nz: i32,
    ndist: i32,
    npix: i32,
    zstep: f32,
    damax: &[f32],
    rho: f32,
    sumnumber: &mut [f64],
    sumnumbersqr: &mut [f64],
    gofz: &mut F2<f32>,
) {
    let pi = std::f32::consts::PI;
    let dmax = st.mol.dmax;
    let nd = ndist as usize;
    let mut voidtest = vec![0.0f32; nd + 1];
    let mut rhotest = vec![0i32; nd + 1];
    let mut ltyp = vec![0i32; nd + 1];
    for id in 1..=nd {
        voidtest[id] = damax[id] / dmax;
        rhotest[id] = nint(4.0 * pi * rho * damax[id] * damax[id] * damax[id] / 3.0);
        sumnumber[id] = 0.0;
        sumnumbersqr[id] = 0.0;
    }
    let l = st.lat.clone();
    let mut i = 0;
    while i < npix {
        let fa = st.rng.ran1() - 0.5;
        let fb = st.rng.ran1() - 0.5;
        let fc = st.rng.ran1() - 0.5;
        let zcent = fc * l.clatz;
        let ycent = fb * l.blaty + fc * l.claty;
        let xcent = fa * l.alatx + fb * l.blatx + fc * l.clatx;
        i += 1;
        for v in ltyp.iter_mut() {
            *v = 0;
        }
        let m = &st.mol;
        for im in 1..=m.nmol {
            let mut xdif = m.xadd[im] - xcent;
            let mut ydif = m.yadd[im] - ycent;
            let mut zdif = m.zadd[im] - zcent;
            minimagec(&mut xdif, &mut ydif, &mut zdif);
            let diftest = (xdif * xdif + ydif * ydif + zdif * zdif).sqrt();
            for id in 1..=nd {
                if diftest <= voidtest[id] {
                    ltyp[id] += 1;
                }
            }
        }
        for id in 1..=nd {
            sumnumber[id] += ltyp[id] as f64;
            sumnumbersqr[id] += (ltyp[id] * ltyp[id]) as f64;
            let mut iref = ltyp[id] - rhotest[id] + 1 + nz / 2;
            if iref < 1 {
                iref = 1;
            }
            if iref > nz {
                iref = nz;
            }
            gofz[(iref, id as i32)] += 1.0;
        }
    }
    for id in 1..=nd {
        sumnumber[id] /= npix as f64;
        sumnumbersqr[id] /= npix as f64;
    }
    for id in 1..=ndist {
        let mut sum = 0.0f32;
        for iz in 1..=nz {
            sum += gofz[(iz, id)];
        }
        sum *= zstep;
        if sum > 0.0 {
            for iz in 1..=nz {
                gofz[(iz, id)] /= sum;
            }
        } else {
            for iz in 1..=nz {
                gofz[(iz, id)] = 0.0;
            }
        }
    }
}
