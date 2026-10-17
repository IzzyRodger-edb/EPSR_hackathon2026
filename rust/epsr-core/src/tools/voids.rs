//! `voids` — translation of `fortran_src/voids.f` (`run_voids`, `find_voids`,
//! `occupypix`, `markpixel`, `search_voids`, `get_void_list`, `ranmine`).
//!
//! The box is divided into pixels; pixels within `radius` of an atom are
//! occupied (with an `include` list only the listed atom types are used and
//! occupied/unoccupied are swapped).  Connected unoccupied pixels form voids.
//! For each distribution three columns are written to `<run>.VOIDS.n01`:
//! the histogram of void sizes, the pixel-pixel correlation function of void
//! pixels (sampled from 100 random origin pixels), and the same restricted
//! to pixel pairs in different voids.  The fraction of the box volume in
//! voids is written to `fort.14`; the random-number state is kept in
//! `randomno.dat`.
//!
//! Note (as in the Fortran): after the run `nsumt` and `ndist` are stored in
//! the *`npix`* and *`voidmax`* entries of the input file (the update starts
//! one entry too early), so `nsumt` itself is never advanced.
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> voids <input file>`
//! Keyword table: `crate::symtab_tools::VOIDS`.

use super::partials::anal_common::*;
use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::fmath::nint;
use crate::fortran::format::{fmt, FV};
use crate::fortran::text::{concatenate, parse};
use crate::fortran::{F1, F2};
use crate::neighbour::minimagec;
#[allow(unused_imports)]
use crate::state::State;

const MDIST: usize = crate::state::MFILE;
const MZ: usize = 5000;
const MPIX: i32 = 201;
const MPIX3: i32 = MPIX * MPIX * MPIX;
const MAXVOID: i32 = 32767;

struct Voids {
    nz: i32,
    ndist: i32,
    npix: i32,
    damax: Vec<f32>,
    fract2: Vec<f32>,
    nexclude: Vec<usize>,
    exclude: Vec<Vec<String>>,
}

/// Run the `voids` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    if !ctx.open_input(&mut st, &crate::symtab_tools::VOIDS)? {
        return Ok(());
    }
    if !st.lat.cubic {
        return Ok(());
    }
    let mut iref = 2;
    let nz = sv_int(&st, iref)?.min(MZ as i32);
    iref += 1;
    let mut npix = sv_int(&st, iref)?;
    npix /= 2;
    let npixtotal = 2 * npix + 1;
    if npixtotal > MPIX {
        npix = (MPIX - 1) / 2;
    }
    iref += 1;
    let _voidmax = sv_real(&st, iref)?;
    iref += 1;
    let mut nsumt = sv_int(&st, iref)?;
    let mut ndist = st.sym.nfiles;
    let mut ndist3 = 3 * ndist;
    if ndist3 > MDIST as i32 {
        ndist3 = MDIST as i32;
    }
    ndist = ndist3 / 3;
    let nd = ndist as usize + 1;
    let mut v = Voids { nz, ndist, npix, damax: vec![0.0; nd], fract2: vec![0.0; nd], nexclude: vec![0; nd], exclude: vec![Vec::new(); nd] };
    for id in 1..=ndist {
        let i = id as usize;
        let mut iref = st.sym.ireff(st.sym.nsymboluptodata + 1, id);
        v.damax[i] = sv_real(&st, iref)?;
        iref += 1;
        let sometext = st.sym.symbol_value[iref].clone();
        let w = parse(&sometext, 5);
        for ie in 1..=w.nwords() {
            let (texclude, _) = concatenate(ie, ie, &w, &sometext, 256);
            let t3 = crate::fortran::fstr(&texclude, 3);
            let mut found = false;
            let mut it = 0;
            while it < st.mol.ntype && !found {
                it += 1;
                found = st.mol.type_[it] == t3;
            }
            if found {
                v.nexclude[i] += 1;
                v.exclude[i].push(t3);
            }
        }
    }
    let dmax = st.mol.dmax;
    let zstep = dmax / nz as f32;
    let zhalf = 0.5 * zstep;
    let mut z: F1<f32> = F1::new(MZ);
    for j in 1..=nz {
        z[j] = j as f32 * zstep - zhalf;
    }
    let mut gofz: F2<f32> = F2::new(MZ, MDIST);
    st.normalise_coordinates();
    initran1_file(ctx, &mut st)?;
    find_voids(&mut st, &mut v, &mut gofz)?;

    let titler = "Void_radius";
    let mut titledat = Vec::new();
    for id in 1..=ndist as usize {
        let t = fmt("('   Radius_',f5.3)", &[FV::R(v.damax[id])]);
        for _ in 0..3 {
            titledat.push(t.clone());
        }
    }
    let mut gzt: F2<f32> = F2::new(MZ, MDIST);
    let mut gzt2: F2<f32> = F2::new(MZ, MDIST);
    outgofz(&st, "n", nz, ndist3, &mut nsumt, &mut z, &gofz, &mut gzt, &mut gzt2, titler, &titledat)?;

    let boxv = dmax * dmax * dmax;
    let mut idref3 = 0;
    for id in 1..=ndist as usize {
        let idref1 = idref3 + 1;
        idref3 = idref1 + 2;
        let mut sum1 = 0.0f32;
        for i in 1..=nz {
            let voidl = z[i];
            let voidv = voidl * voidl * voidl;
            sum1 += voidv * gzt[(i, idref1)];
        }
        v.fract2[id] = sum1 / boxv;
    }

    let icount = st.sym.search_for_symbol("nsize", 0);
    let mut iref = st.sym.ireff(icount, 1);
    set_sv_int(&mut st, iref, nz);
    iref += 1;
    set_sv_int(&mut st, iref, nsumt);
    iref += 1;
    set_sv_int(&mut st, iref, ndist);
    st.save_inp_file()?;
    let mut lines = Vec::new();
    for id in 1..=ndist as usize {
        lines.push(fmt("(/'Distribution ',i3,': fraction of pixels used = ',f10.5)", &[FV::I(id as i64), FV::R(v.fract2[id])]));
    }
    write_fort14(ctx, &lines)?;
    saveran1_file(ctx, &st)
}

/// `subroutine find_voids(nexclude,exclude)`.
fn find_voids(st: &mut State, v: &mut Voids, gofz: &mut F2<f32>) -> Result<()> {
    let nz = v.nz;
    let natom = st.mol.natom;
    let dmax = st.mol.dmax;
    let volatom = 1.0f32 / natom as f32;
    let zstep = 1.0f32 / nz as f32;
    let mut voidpix: Vec<i8> = Vec::new();
    let mut voidnum: Vec<i16> = Vec::new();
    let mut idref3 = 0;
    for id in 1..=v.ndist as usize {
        let idref1 = idref3 + 1;
        let idref2 = idref1 + 1;
        idref3 = idref2 + 1;
        let radtest = v.damax[id] / dmax;
        let mut pixstep = 1.0f32 / (2 * v.npix + 1) as f32;
        let mut npixsphere = nint(radtest / pixstep);
        pixstep = radtest / npixsphere as f32;
        let ratio = 1.0f32 / pixstep;
        let mut npixt = (ratio / 2.0) as i32;
        let remain = ratio - (2 * npixt + 1) as f32;
        st.log(&format!(" {ratio} {npixt} {remain}"));
        if remain > 1.0 {
            npixt += 1;
        }
        if 2 * npixt + 1 > MPIX {
            npixt = (MPIX - 1) / 2;
        }
        pixstep = 1.0f32 / (2 * npixt + 1) as f32;
        npixsphere = nint(radtest / pixstep);
        let radactual = npixsphere as f32 * pixstep;
        let npixtot = 2 * npixt + 1;
        let npixall = npixtot * npixtot * npixtot;
        st.log(&fmt(
            "(//'find_voids> Distribution ',i2//'find_voids> Adjusting number of pixels to give optimum'/ '            fit to specified void distance.',//'find_voids> New number of pixels = ', i8,//'find_voids> Specified radius = ',f10.5,//'find_voids> Actual radius used = ',f10.5)",
            &[FV::I(id as i64), FV::I(npixtot as i64), FV::R(radtest * dmax), FV::R(radactual * dmax)],
        ));
        let np = npixall as usize;
        voidpix.clear();
        voidpix.resize(np + 1, 0);
        voidnum.clear();
        voidnum.resize(np + 1, 0);

        let mut npixfound = 0;
        {
            let m = &st.mol;
            for ia in 1..=natom {
                let mut ifound = 0;
                if v.nexclude[id] > 0 {
                    let mut ie = 0;
                    while ie < v.nexclude[id] && ifound == 0 {
                        ie += 1;
                        if m.type_[m.itype[ia]] == v.exclude[id][ie - 1] {
                            ifound = m.itype[ia];
                        }
                    }
                } else {
                    ifound = m.itype[ia];
                }
                if ifound != 0 {
                    let im = m.mref[ia];
                    let mut xatom = m.xa[ia] + m.xadd[im];
                    let mut yatom = m.ya[ia] + m.yadd[im];
                    let mut zatom = m.za[ia] + m.zadd[im];
                    minimagec(&mut xatom, &mut yatom, &mut zatom);
                    let ih = nint(xatom / pixstep);
                    let ik = nint(yatom / pixstep);
                    let il = nint(zatom / pixstep);
                    occupypix(npixt, npixsphere, &mut npixfound, ih, ik, il, &mut voidpix)?;
                }
            }
        }
        if v.nexclude[id] > 0 {
            for ipix in 1..=np {
                voidpix[ipix] = if voidpix[ipix] != 0 { 0 } else { 1 };
            }
        }
        st.log(&fmt("(/'find_voids> No. of pixels in sphere = ',i8)", &[FV::I(npixfound as i64)]));
        let mut nvoidpix = 0;
        for ipix in 1..=np {
            if voidpix[ipix] == 0 {
                nvoidpix += 1;
            }
        }
        v.fract2[id] = nvoidpix as f32 / npixall as f32;
        if nvoidpix == 0 {
            return crate::error::unsupported("voids: no unoccupied pixels (the Fortran loops forever choosing a random void pixel)");
        }

        let mut nvoidsold = nvoidpix;
        let mut nvoidsnew = 0;
        let mut it = 0;
        while nvoidsold > nvoidsnew {
            if it > 0 {
                nvoidsold = nvoidsnew;
            }
            search_voids(&mut nvoidsnew, npixt, npixall, &voidpix, &mut voidnum)?;
            it += 1;
            st.log(&format!(" {it} {nvoidsold} {nvoidsnew}"));
        }

        // tally the size of each void (in order of first appearance)
        let mut slot = vec![0i32; MAXVOID as usize + 1];
        let mut lenvoid: Vec<i32> = Vec::new();
        for ipix in 1..=np {
            let vn = voidnum[ipix];
            if voidpix[ipix] == 0 && vn > 0 {
                let s = &mut slot[vn as usize];
                if *s == 0 {
                    lenvoid.push(1);
                    *s = lenvoid.len() as i32;
                } else {
                    lenvoid[(*s - 1) as usize] += 1;
                }
            }
        }
        let volpix = pixstep * pixstep * pixstep;
        let _volpixstar = volpix / volatom;
        let third = 1.0f32 / 3.0f32;
        for &len in &lenvoid {
            let mut iref = ((len as f32 * volpix).powf(third) / zstep) as i32 + 1;
            if iref > nz {
                iref = nz;
            }
            gofz[(iref, idref1 as i32)] += 1.0;
        }

        // radial distribution function of void pixels
        let npixtot2 = npixtot * npixtot;
        let nvcount = 100;
        let mut ivcount = 0;
        let pixhwidth = 0.5f32 / npixtot as f32;
        let mut iscount = 0i32;
        let i2 = idref2 as i32;
        let i3 = idref3 as i32;
        let nzf = nz as f32;
        while ivcount < nvcount {
            ivcount += 1;
            let mut iv1;
            loop {
                iv1 = (npixall as f32 * st.rng.ran1() + 1.0) as i32;
                if voidpix[iv1 as usize] == 0 {
                    break;
                }
            }
            let xdif = pixhwidth * (2.0 * st.rng.ran1() - 1.0);
            let ydif = pixhwidth * (2.0 * st.rng.ran1() - 1.0);
            let zdif = pixhwidth * (2.0 * st.rng.ran1() - 1.0);
            let dist = (xdif * xdif + ydif * ydif + zdif * zdif).sqrt();
            let mut iref = (dist * nzf) as i32 + 1;
            if iref < 1 {
                iref = 1;
            }
            if iref > nz {
                iref = nz;
            }
            gofz[(iref, i2)] += 1.0;
            let ilv1 = (iv1 - 1) / npixtot2;
            let ikv1 = (iv1 - 1 - ilv1 * npixtot2) / npixtot;
            let ihv1 = iv1 - 1 - npixtot * ikv1 - npixtot2 * ilv1;
            let xv1 = ihv1 as f32 / npixtot as f32;
            let yv1 = ikv1 as f32 / npixtot as f32;
            let zv1 = ilv1 as f32 / npixtot as f32;
            iscount = 0;
            let vn1 = voidnum[iv1 as usize];
            for iv2 in 1..=npixall {
                if voidpix[iv2 as usize] == 0 && iv1 != iv2 {
                    iscount += 1;
                    let ivt = iv2;
                    let ilv2 = (ivt - 1) / npixtot2;
                    let ikv2 = (ivt - 1 - ilv2 * npixtot2) / npixtot;
                    let ihv2 = ivt - 1 - npixtot * ikv2 - npixtot2 * ilv2;
                    let xv2 = ihv2 as f32 / npixtot as f32;
                    let yv2 = ikv2 as f32 / npixtot as f32;
                    let zv2 = ilv2 as f32 / npixtot as f32;
                    let mut xdif = xv2 - xv1;
                    let mut ydif = yv2 - yv1;
                    let mut zdif = zv2 - zv1;
                    xdif += pixhwidth * (2.0 * st.rng.ran1() - 1.0);
                    ydif += pixhwidth * (2.0 * st.rng.ran1() - 1.0);
                    zdif += pixhwidth * (2.0 * st.rng.ran1() - 1.0);
                    minimagec(&mut xdif, &mut ydif, &mut zdif);
                    let dist = (xdif * xdif + ydif * ydif + zdif * zdif).sqrt();
                    if xdif > 0.5 || ydif > 0.5 || zdif > 0.5 {
                        st.log(&format!(" {ihv1} {ikv1} {ilv1} {ihv2} {ikv2} {ilv2} {npixtot}"));
                    } else {
                        let mut iref = (dist * nzf) as i32 + 1;
                        if iref < 1 {
                            iref = 1;
                        }
                        if iref > nz {
                            iref = nz;
                        }
                        gofz[(iref, i2)] += 1.0;
                        if vn1 != voidnum[iv2 as usize] {
                            gofz[(iref, i3)] += 1.0;
                        }
                    }
                }
            }
        }
        let rhovoid = (nvcount * iscount) as f32;
        let rvstep = 1.0f32 / nz as f32;
        let pi = 3.141592653f32;
        let fact = 3.0f32 / (4.0 * pi * rhovoid);
        let mut rv13 = 0.0f32;
        for iz in 1..=nz {
            let rv2 = iz as f32 * rvstep;
            let rv23 = rv2 * (rv2 * rv2);
            let rvfact = fact / (rv23 - rv13);
            gofz[(iz, i2)] *= rvfact;
            gofz[(iz, i3)] *= rvfact;
            rv13 = rv23;
        }
    }
    Ok(())
}

/// `subroutine occupypix(npix,npixsphere,npixfound,ih,ik,il,voidpix,mpix3)`.
fn occupypix(npix: i32, npixsphere: i32, npixfound: &mut i32, ih: i32, ik: i32, il: i32, voidpix: &mut [i8]) -> Result<()> {
    let npixsphere2 = npixsphere * npixsphere;
    *npixfound = 0;
    let mut mark = |a: i32, b: i32, c: i32, n: &mut i32| markpixel(npix, a, b, c, n, voidpix);
    for ils in 0..=npixsphere {
        let ils2 = ils * ils;
        for iks in 0..=npixsphere {
            let iks2 = iks * iks;
            let iks2pils2 = ils2 + iks2;
            for ihs in 0..=npixsphere {
                let ihs2 = ihs * ihs;
                if iks2pils2 + ihs2 <= npixsphere2 {
                    mark(ih + ihs, ik + iks, il + ils, npixfound)?;
                    if ihs > 0 {
                        mark(ih - ihs, ik + iks, il + ils, npixfound)?;
                    }
                    if iks > 0 {
                        mark(ih + ihs, ik - iks, il + ils, npixfound)?;
                        if ihs > 0 {
                            mark(ih - ihs, ik - iks, il + ils, npixfound)?;
                        }
                    }
                    if ils > 0 {
                        mark(ih + ihs, ik + iks, il - ils, npixfound)?;
                        if ihs > 0 {
                            mark(ih - ihs, ik + iks, il - ils, npixfound)?;
                        }
                        if iks > 0 {
                            mark(ih + ihs, ik - iks, il - ils, npixfound)?;
                            if ihs > 0 {
                                mark(ih - ihs, ik - iks, il - ils, npixfound)?;
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

/// `subroutine markpixel(npix,ihabs,ikabs,ilabs,npixfound,voidpix,mpix3)`.
fn markpixel(npix: i32, ihabs: i32, ikabs: i32, ilabs: i32, npixfound: &mut i32, voidpix: &mut [i8]) -> Result<()> {
    let offset = npix;
    let npixtot = 2 * npix + 1;
    let wrap = |a: i32| a - npixtot * nint(a as f32 / npixtot as f32) + offset;
    let ihref = wrap(ihabs);
    let ikref = wrap(ikabs);
    let ilref = wrap(ilabs);
    let npixref = 1 + ihref + npixtot * (ikref + npixtot * ilref);
    if npixref <= MPIX3 && npixref >= 1 && (npixref as usize) < voidpix.len() {
        voidpix[npixref as usize] = 1;
        *npixfound += 1;
        Ok(())
    } else {
        Err(EpsrError::Input(format!(
            "markpixel> Error generating pixel reference numbers: {ihabs} {ikabs} {ilabs} {ihref} {ikref} {ilref} {npixref}"
        )))
    }
}

/// `subroutine search_voids(nvoids,npix,npixall,voidpix,voidnum,mpix3)`.
fn search_voids(nvoids: &mut i32, npix: i32, npixall: i32, voidpix: &[i8], voidnum: &mut [i16]) -> Result<()> {
    let mut nvoid = *nvoids;
    let mut voidlist = [0i32; 27];
    for ipix in 1..=npixall {
        if voidpix[ipix as usize] == 0 {
            let nvoidlist = get_void_list(ipix, npix, &mut voidlist, voidpix);
            let mut ilowest = MAXVOID;
            for &p in &voidlist[..nvoidlist] {
                let itest = voidnum[p as usize] as i32;
                if itest < ilowest && itest > 0 {
                    ilowest = itest;
                }
            }
            if ilowest == MAXVOID && nvoid < MAXVOID {
                nvoid += 1;
                ilowest = nvoid;
            } else if nvoid >= MAXVOID {
                return Err(EpsrError::Input("search_voids> Unable to allocate any more void numbers! Try using a larger void distance?".into()));
            }
            for &p in &voidlist[..nvoidlist] {
                voidnum[p as usize] = ilowest as i16;
            }
        }
    }
    // count the distinct void numbers of the unoccupied pixels
    let mut seen = vec![false; 65536];
    *nvoids = 0;
    for ipix in 1..=npixall as usize {
        if voidpix[ipix] == 0 {
            let k = voidnum[ipix] as u16 as usize;
            if !seen[k] {
                seen[k] = true;
                *nvoids += 1;
            }
        }
    }
    Ok(())
}

/// `subroutine get_void_list(ipix,npix,nvoidlist,voidlist,voidpix,mpix3)`.
fn get_void_list(ipix: i32, npix: i32, voidlist: &mut [i32; 27], voidpix: &[i8]) -> usize {
    let offset = npix;
    let npixtot = 2 * npix + 1;
    let npixtot2 = npixtot * npixtot;
    let mut n = 0;
    let mut ilc = (ipix - 1) / npixtot2;
    let mut ikc = (ipix - 1 - ilc * npixtot2) / npixtot;
    let mut ihc = ipix - 1 - npixtot * ikc - npixtot2 * ilc;
    ilc -= offset;
    ikc -= offset;
    ihc -= offset;
    let wrap = |a: i32| a - npixtot * nint(a as f32 / npixtot as f32);
    for ilt in ilc - 1..=ilc + 1 {
        let il = wrap(ilt);
        let iladd = (il + offset) * npixtot2;
        for ikt in ikc - 1..=ikc + 1 {
            let ik = wrap(ikt);
            let ikadd = iladd + (ik + offset) * npixtot;
            for iht in ihc - 1..=ihc + 1 {
                let ih = wrap(iht);
                let pixref = 1 + ih + offset + ikadd;
                if voidpix[pixref as usize] == 0 {
                    voidlist[n] = pixref;
                    n += 1;
                }
            }
        }
    }
    n
}
