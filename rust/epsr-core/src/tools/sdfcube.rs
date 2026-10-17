//! `sdfcube` — translation of `fortran_src/sdfcube.f`.
//!
//! Spatial density function of second molecules around a central molecule,
//! accumulated directly on a cubic grid in the frame of the central molecule
//! and written as a Gaussian `.cube` file (`<run>.SDFCUBE.cube`).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> sdfcube <input file>`
//! Keyword table: `crate::symtab_tools::SDFCUBE`.
//!
//! Accumulation over runs (`nsumt > 0`): `outcube` reads the existing cube and
//! forms `fc1*old + fc*new` with `fc = 1/nsumt`; `fc1` is computed as
//! `fc1 = 1.0 - fc1` from an uninitialised variable, which gfortran -O2 folds
//! to `fc1 = 1.0` (checked in the reference object code) — reproduced here.
//! For identical central and second molecules `averagemol` is read from the
//! `axiss1` symbol (missing `iref=iref+2` in the Fortran) — also reproduced.

use super::sharm::sharm_routines::{form_eulerangles, list_int, rotateeuler, write_mol_xyz};
use super::sharm::{average_molecule, find_type, has_types, open, read_axis, words, Mols, MCENTRE};
use super::ToolContext;
use crate::ato::get_atomicnumber;
use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::{read_int, read_real, TextFile};
use crate::fortran::output::OutFile;
use crate::fortran::text::{findex, indexlastdot};
use crate::fortran::{F1, F3};
use crate::state::State;

const BOHR: f32 = 0.52917721092;

/// Module `sdfcube_mod` + common `sdfcube.inc` variables used here.
#[derive(Default)]
struct Cube {
    /// Folder-qualified path of `fname_cube`.
    cube_path: std::path::PathBuf,
    nx: i32,
    ny: i32,
    nz: i32,
    nx0: i32,
    ny0: i32,
    nz0: i32,
    nsumt: i32,
    xstep: f32,
    ystep: f32,
    zstep: f32,
    natosave: i32,
    xasave: Vec<f32>,
    yasave: Vec<f32>,
    zasave: Vec<f32>,
    atsymbolsave: Vec<String>,
    rhocube: F3<f32>,
    fname_cube: String,
}

/// Run the `sdfcube` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    if !open(ctx, &mut st, &crate::symtab_tools::SDFCUBE)? {
        return Ok(());
    }
    let log = st.logger();
    // The reference build has no system_commands.txt, so get_full_filename
    // never prepends the folder: the Fortran's `fnameinp` (used in the text of
    // the cube and xyz files) is the bare input file name.
    let fnameinp = st.sym.fname_inp.clone();
    let lenfname_cube = findex(&fnameinp, ".SDFCUBE.dat");
    if lenfname_cube <= 0 {
        return Ok(());
    }
    let fname_cube = format!("{}.SDFCUBE.cube", &fnameinp[..(lenfname_cube - 1) as usize]);
    let mut cb = Cube { cube_path: ctx.path(&fname_cube), fname_cube, ..Default::default() };
    let bad = |what: &str| EpsrError::Input(format!("sdfcube: cannot read {what}"));
    cb.nz = read_int(&st.sym.symbol_value[2]).map_err(|_| bad("mpixel"))?;
    let zmax = read_real(&st.sym.symbol_value[3]).map_err(|_| bad("rmax"))?;
    cb.nsumt = read_int(&st.sym.symbol_value[4]).map_err(|_| bad("nsumt"))?;
    let (_, ws) = words(&st, 5, 5);
    let mut lctyp1: Vec<i32> = Vec::new();
    for w in &ws {
        let f = find_type(&st, w);
        if f <= 0 {
            log(&format!("\nsharm> Atom type {} does not exist in .ato file {}", w, st.sym.fname_ato));
            return Ok(());
        }
        if (lctyp1.len() as i32) < MCENTRE {
            lctyp1.push(f);
        }
    }
    let (nzt1, lztyp1) = read_axis(&st, 6, true)?;
    let (nyt1, lytyp1) = read_axis(&st, 7, true)?;
    let (_, ws) = words(&st, 8, 5);
    let mut lctyp2: Vec<i32> = Vec::new();
    for w in &ws {
        let f = find_type(&st, w);
        if f == 0 {
            log(&format!("\nsharm> Atom type {} does not exist in .ato file {}", w, st.sym.fname_ato));
            return Ok(());
        }
        if (lctyp2.len() as i32) < MCENTRE {
            lctyp2.push(f);
        }
    }
    let ldistinct = if lctyp2.len() == lctyp1.len() { if lctyp1.iter().zip(&lctyp2).any(|(a, b)| a != b) { 1 } else { 0 } } else { 1 };
    let (nzt2, lztyp2, nyt2, lytyp2);
    if ldistinct != 0 {
        let a = read_axis(&st, 9, true)?;
        let b = read_axis(&st, 10, true)?;
        nzt2 = a.0;
        lztyp2 = a.1;
        nyt2 = b.0;
        lytyp2 = b.1;
    } else {
        lctyp2 = lctyp1.clone();
        nzt2 = nzt1.clone();
        nyt2 = nyt1.clone();
        lztyp2 = lztyp1.clone();
        lytyp2 = lytyp1.clone();
    }
    // For identical molecules the Fortran does not skip the two axis symbols
    // of the second molecule (no `iref=iref+2`), so `averagemol` is then read
    // from symbol 9 (`axiss1`).
    let outputaveragemol = st.sym.symbol_value[if ldistinct != 0 { 11 } else { 9 }] == "y";
    let nmult = 5;
    let nsum = cb.nsumt / nmult;
    if nsum * nmult == cb.nsumt {
        st.normalise_coordinates();
        cb.zstep = zmax / (cb.nz as f32);
        initcube(&mut cb);
        let mol = Mols { lctyp1: &lctyp1, lctyp2: &lctyp2, nzt1: &nzt1, nyt1: &nyt1, lztyp1: &lztyp1, lytyp1: &lytyp1, nzt2: &nzt2, nyt2: &nyt2, lztyp2: &lztyp2, lytyp2: &lytyp2 };
        sdfcube(ctx, &mut st, &mut cb, &mol, ldistinct, outputaveragemol, &log)?;
        outcube(&mut cb)?;
    }
    cb.nsumt += 1;
    st.sym.symbol_value[4] = list_int(cb.nsumt);
    st.save_inp_file()?;
    Ok(())
}

/// `subroutine initcube(nin)`.
fn initcube(cb: &mut Cube) {
    let mut opened = false;
    if cb.nsumt > 0 {
        if let Ok(mut f) = TextFile::open(&cb.cube_path) {
            opened = true;
            let _ = f.list().finish();
            let _ = f.list().finish();
            {
                let mut r = f.list();
                if let Ok(Some(v)) = r.int() {
                    cb.natosave = v;
                }
            }
            {
                let mut r = f.list();
                if let Ok(Some(v)) = r.int() {
                    cb.nz0 = v;
                }
                if let Ok(Some(v)) = r.real() {
                    cb.zstep = v;
                }
            }
            cb.zstep *= BOHR;
            cb.nz = (cb.nz0 - 1) / 2;
            cb.nx = cb.nz;
            cb.ny = cb.nz;
            cb.nz0 = 2 * cb.nz + 1;
            cb.nx0 = cb.nz0;
            cb.ny0 = cb.nz0;
            cb.xstep = cb.zstep;
            cb.ystep = cb.zstep;
        }
    }
    if !opened {
        cb.nsumt = 0;
        cb.nx = cb.nz;
        cb.ny = cb.nz;
        cb.nz0 = 2 * cb.nz + 1;
        cb.nx0 = cb.nz0;
        cb.ny0 = cb.nz0;
        cb.xstep = cb.zstep;
        cb.ystep = cb.zstep;
    }
    let n = cb.nz0.max(0) as usize;
    cb.rhocube = F3::new(n, n, n);
}

/// `subroutine outcube(nout)`.
fn outcube(cb: &mut Cube) -> Result<()> {
    let path = cb.cube_path.clone();
    if cb.nsumt > 0 {
        match TextFile::open(&path) {
            Ok(mut f) => {
                let bad = |_| EpsrError::Input(format!("sdfcube: error reading {}", cb.fname_cube));
                f.list().finish().map_err(bad)?;
                f.list().finish().map_err(bad)?;
                let natosavet = { f.list().int().map_err(bad)?.unwrap_or(0) };
                for _ in 0..3 {
                    f.list().finish().map_err(bad)?;
                }
                for _ in 0..natosavet {
                    f.list().finish().map_err(bad)?;
                }
                let fc: f32 = 1.0f32 / (cb.nsumt as f32);
                let fc1: f32 = 1.0; // = 1.0 - (uninitialised fc1), folded by gfortran
                for ix in 1..=cb.nx0 {
                    for iy in 1..=cb.nz0 {
                        let mut iz = 1;
                        while iz <= cb.nz0 {
                            let itmax = 9.min(cb.nz0 - iz);
                            let mut t = [0.0f32; 10];
                            {
                                let mut r = f.list();
                                for v in t.iter_mut().take((itmax + 1) as usize) {
                                    if let Some(x) = r.real().map_err(bad)? {
                                        *v = x;
                                    }
                                }
                            }
                            for it in 0..=itmax {
                                let izref = iz + it;
                                cb.rhocube[(ix, iy, izref)] = fc1 * t[it as usize] + fc * cb.rhocube[(ix, iy, izref)];
                            }
                            iz += 10;
                        }
                    }
                }
            }
            Err(_) => cb.nsumt = 0,
        }
    }
    let mut out = OutFile::create(&path);
    out.line(&ListWriter::new().chars("# Cube file output by sdfcube").finish());
    out.line(&ListWriter::new().chars(&format!("Outer loop x, middle loop y, inner loop z, filname: {}", cb.fname_cube.trim_end())).finish());
    out.line(
        &ListWriter::new()
            .int(cb.natosave)
            .real(-((cb.nx as f32) * cb.xstep / BOHR))
            .real(-((cb.ny as f32) * cb.ystep / BOHR))
            .real(-((cb.nz as f32) * cb.zstep / BOHR))
            .finish(),
    );
    out.line(&ListWriter::new().int(cb.nx0).real(cb.xstep / BOHR).real(0.0).real(0.0).finish());
    out.line(&ListWriter::new().int(cb.ny0).real(0.0).real(cb.ystep / BOHR).real(0.0).finish());
    out.line(&ListWriter::new().int(cb.nz0).real(0.0).real(0.0).real(cb.zstep / BOHR).finish());
    for is in 0..cb.natosave.max(0) as usize {
        let iatomicno = get_atomicnumber(cb.atsymbolsave[is].trim_end());
        out.line(&ListWriter::new().int(iatomicno).real(iatomicno as f32).real(cb.xasave[is] / BOHR).real(cb.yasave[is] / BOHR).real(cb.zasave[is] / BOHR).finish());
    }
    for ix in 1..=cb.nx0 {
        for iy in 1..=cb.ny0 {
            let mut iz = 1;
            while iz <= cb.nz0 {
                let itmax = 9.min(cb.nz0 - iz);
                let v: Vec<FV> = (0..=itmax).map(|it| FV::R(cb.rhocube[(ix, iy, iz + it)])).collect();
                out.line(&fmt("(10(1x,e11.4))", &v));
                iz += 10;
            }
        }
    }
    out.close()?;
    Ok(())
}

/// `subroutine sdfcube(nout,fname_inp)`.
fn sdfcube(ctx: &ToolContext, st: &mut State, cb: &mut Cube, mol: &Mols, ldistinct: i32, outputaveragemol: bool, log: &dyn Fn(&str)) -> Result<()> {
    let nmol = st.mol.nmol;
    let mut sh = super::sharm::sharm_routines::Sharm { outputaveragemol, ldistinct, ..Default::default() };
    let n = nmol.max(0) as usize;
    sh.nm = F1::new(n);
    sh.phim = F1::new(n);
    sh.cthetm = F1::new(n);
    sh.chim = F1::new(n);
    sh.xcent = F1::new(n);
    sh.ycent = F1::new(n);
    sh.zcent = F1::new(n);
    let dmax = st.mol.dmax;
    let xstepd = cb.xstep / dmax;
    let ystepd = cb.ystep / dmax;
    let zstepd = cb.zstep / dmax;
    let mut nsmmo1 = 0;
    let mut nsmmo2 = 0;
    for imol1 in 1..=nmol {
        let mut xo1: f32 = 0.0;
        let mut yo1: f32 = 0.0;
        let mut zo1: f32 = 0.0;
        let mut amasssum: f32 = 0.0;
        sh.nm[imol1] = 0;
        if has_types(st, imol1, mol.lctyp1) {
            sh.nm[imol1] = 1;
        } else if has_types(st, imol1, mol.lctyp2) {
            sh.nm[imol1] = 2;
        }
        if sh.nm[imol1] == 1 {
            nsmmo1 += 1;
        }
        if sh.nm[imol1] == 2 {
            nsmmo2 += 1;
        }
        let m = &st.mol;
        let lct: Option<&[i32]> = if sh.nm[imol1] == 1 && !mol.lctyp1.is_empty() {
            Some(mol.lctyp1)
        } else if sh.nm[imol1] == 2 && !mol.lctyp2.is_empty() {
            Some(mol.lctyp2)
        } else {
            None
        };
        for ia1 in m.naf[imol1]..=m.nal[imol1] {
            match lct {
                Some(l) => {
                    for &t in l {
                        if m.itype[ia1] == t {
                            let amass = m.amassa[m.itype[ia1]];
                            xo1 += m.xa[ia1] * amass;
                            yo1 += m.ya[ia1] * amass;
                            zo1 += m.za[ia1] * amass;
                            amasssum += amass;
                        }
                    }
                }
                None => {
                    let amass = m.amassa[m.itype[ia1]];
                    xo1 += m.xa[ia1] * amass;
                    yo1 += m.ya[ia1] * amass;
                    zo1 += m.za[ia1] * amass;
                    amasssum += amass;
                }
            }
        }
        if amasssum > 0.0 && sh.nm[imol1] != 0 {
            xo1 /= amasssum;
            yo1 /= amasssum;
            zo1 /= amasssum;
            sh.xcent[imol1] = m.xadd[imol1] + xo1;
            sh.ycent[imol1] = m.yadd[imol1] + yo1;
            sh.zcent[imol1] = m.zadd[imol1] + zo1;
            let na = m.nato[imol1] as usize;
            let mut xat = F1::new(na);
            let mut yat = F1::new(na);
            let mut zat = F1::new(na);
            for ia in 1..=m.nato[imol1] {
                let iaref = ia + m.naf[imol1] - 1;
                xat[ia] = m.xa[iaref] - xo1;
                yat[ia] = m.ya[iaref] - yo1;
                zat[ia] = m.za[iaref] - zo1;
            }
            let (p, c, ch) = if sh.nm[imol1] == 1 {
                form_eulerangles(&xat, &yat, &zat, mol.nzt1, mol.nyt1, mol.lztyp1, mol.lytyp1)
            } else {
                form_eulerangles(&xat, &yat, &zat, mol.nzt2, mol.nyt2, mol.lztyp2, mol.lytyp2)
            };
            sh.phim[imol1] = p;
            sh.cthetm[imol1] = c;
            sh.chim[imol1] = ch;
        } else {
            sh.nm[imol1] = 0;
        }
    }
    log(&format!(" Number of type 1 molecules = {:>10}\n Number of type 2 molecules = {:>10}", nsmmo1, nsmmo2));
    let fname_inp = st.sym.fname_inp.clone();
    let stem = fname_inp[..(indexlastdot(&fname_inp) - 1).max(0) as usize].to_string();
    if let Some(acc) = average_molecule(st, &sh, 1, mol, false) {
        cb.natosave = acc.len() as i32;
        cb.xasave = acc.iter().map(|a| a.1).collect();
        cb.yasave = acc.iter().map(|a| a.2).collect();
        cb.zasave = acc.iter().map(|a| a.3).collect();
        cb.atsymbolsave = acc.iter().map(|a| a.0.clone()).collect();
        let f = format!("{stem}-first.xyz");
        write_mol_xyz(&ctx.path(&f), &f, &acc)?;
        if ldistinct == 0 {
            let f = format!("{stem}-second.xyz");
            write_mol_xyz(&ctx.path(&f), &f, &acc)?;
        }
    }
    if let Some(acc) = average_molecule(st, &sh, 2, mol, false) {
        let f = format!("{stem}-second.xyz");
        write_mol_xyz(&ctx.path(&f), &f, &acc)?;
    }
    // pair loop
    let (mut lsum1, mut lsum2);
    if ldistinct == 0 {
        if sh.nm[nmol] == 1 || sh.nm[nmol] == 2 {
            lsum1 = 1;
            lsum2 = 1;
        } else {
            lsum1 = 0;
            lsum2 = 0;
        }
    } else {
        lsum1 = 0;
        lsum2 = 0;
        if sh.nm[nmol] == 1 {
            lsum1 = 1;
        }
        if sh.nm[nmol] == 2 {
            lsum2 = 1;
        }
    }
    let nint = crate::fortran::fmath::nint;
    let bin = |cb: &mut Cube, xdifr: f32, ydifr: f32, zdifr: f32| {
        let ix = cb.nx + nint(xdifr / xstepd) + 1;
        if ix > 0 && ix <= cb.nx0 {
            let iy = cb.ny + nint(ydifr / ystepd) + 1;
            if iy > 0 && iy <= cb.ny0 {
                let iz = cb.nz + nint(zdifr / zstepd) + 1;
                if iz > 0 && iz <= cb.nz0 {
                    cb.rhocube[(ix, iy, iz)] += 1.0f32;
                }
            }
        }
    };
    for imol1 in 1..nmol {
        if sh.nm[imol1] == 0 {
            continue;
        }
        if ldistinct == 0 {
            if sh.nm[imol1] == 1 || sh.nm[imol1] == 2 {
                lsum1 += 1;
                lsum2 += 1;
            }
        } else {
            if sh.nm[imol1] == 1 {
                lsum1 += 1;
            }
            if sh.nm[imol1] == 2 {
                lsum2 += 1;
            }
        }
        let xc1 = sh.xcent[imol1];
        let yc1 = sh.ycent[imol1];
        let zc1 = sh.zcent[imol1];
        for imol2 in imol1 + 1..=nmol {
            if imol1 == imol2 || sh.nm[imol2] == 0 {
                continue;
            }
            let nm1 = sh.nm[imol1];
            let nm2 = sh.nm[imol2];
            let (mut iflag1, mut iflag2) = (0, 0);
            if ldistinct == 0 && nm1 == 1 && nm2 == 1 {
                iflag1 = 1;
                iflag2 = 1;
            }
            if ldistinct == 1 && nm1 == 1 && nm2 == 2 {
                iflag1 = 1;
                iflag2 = 0;
            }
            if ldistinct == 1 && nm1 == 2 && nm2 == 1 {
                iflag1 = 0;
                iflag2 = 1;
            }
            let mut xdif = sh.xcent[imol2] - xc1;
            let mut ydif = sh.ycent[imol2] - yc1;
            let mut zdif = sh.zcent[imol2] - zc1;
            crate::neighbour::minimagec(&mut xdif, &mut ydif, &mut zdif);
            if iflag1 == 1 {
                let (mut a, mut b, mut c) = rotateeuler(xdif, ydif, zdif, -sh.chim[imol1], -sh.cthetm[imol1].acos(), -sh.phim[imol1]);
                crate::neighbour::minimagec(&mut a, &mut b, &mut c);
                bin(cb, a, b, c);
            }
            if iflag2 == 1 {
                let (mut a, mut b, mut c) = rotateeuler(-xdif, -ydif, -zdif, -sh.chim[imol2], -sh.cthetm[imol2].acos(), -sh.phim[imol2]);
                crate::neighbour::minimagec(&mut a, &mut b, &mut c);
                bin(cb, a, b, c);
            }
        }
    }
    let vol = st.lat.alatx * (st.lat.blaty * st.lat.clatz);
    let rho: f32 = (lsum1 as f32) * (lsum2 as f32) / vol;
    log(&fmt("(1x,'rho = ',e13.6)", &[FV::R(rho)]));
    let gfac = xstepd * ystepd * zstepd;
    let gfac1 = 1.0f32 / (gfac * rho);
    for v in cb.rhocube.data.iter_mut() {
        *v *= gfac1;
    }
    Ok(())
}
