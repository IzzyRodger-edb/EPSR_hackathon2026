//! `sdf` — translation of `fortran_src/sdf_rev2.f` (+ `sharm_routines_rev1.f`).
//!
//! Like `sharm`, but the "molecules" are groups of atoms found by type and
//! distance criteria (first listed atom = centre), so that e.g. a fragment of
//! a molecule, or a group of atoms of different molecules, can be used.
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> sdf <input file>`
//! Keyword table: `crate::symtab_tools::SDF`.
//!
//! As for `sharm`, the column titles are written as `(l1,l2,l;n1,n2)` (the
//! Fortran run-time format is rejected by gfortran >= 14, see `sharm.rs`).
//! The Fortran quirks reproduced: the check of the remaining intra-group
//! distances uses `zdiff = ytemp(i2) - ztemp(i1)`, and when a second molecule
//! sharing atoms with a first molecule is removed from the list only
//! `nm, phim, cthetm, chim, iasave` are shifted (not the centres or the
//! rotation-matrix arrays); with `test` the orientations used are all zero
//! (the angles read by `run_sdf` are local variables, not the common-block
//! ones used by `sdf`).

use super::sharm::sharm_routines::{form_eulerangles, list_int, rotateeuler, write_mol_xyz, Sharm};
use super::sharm::{alloc_gofz, ch2, find_type, open, read_axis, read_l_values, read_test, words, MCENTRE};
use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::listread::{read_int, read_real};
use crate::fortran::output::OutFile;
use crate::fortran::text::{concatenate, indexlastdot, parse, parse_n};
use crate::fortran::{F1, F2};
use crate::state::State;

/// Definition of the "molecules" (`sharm.inc` variables of SDF).
struct Def {
    lctyp1: Vec<i32>,
    lctyp2: Vec<i32>,
    distfirst: Vec<f32>,
    distsecond: Vec<f32>,
    fractfirst: f32,
    fractsecond: f32,
    nzt1: String,
    nyt1: String,
    nzt2: String,
    nyt2: String,
    lztyp1: Vec<i32>,
    lytyp1: Vec<i32>,
    lztyp2: Vec<i32>,
    lytyp2: Vec<i32>,
}

/// Run the `sdf` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    if !open(ctx, &mut st, &crate::symtab_tools::SDF)? {
        return Ok(());
    }
    let log = st.logger();
    let mut sh = Sharm { pi: std::f32::consts::PI, ..Default::default() };
    sh.piconv = sh.pi / 180.0f32;
    let (mut zmax, mut nsumt, ierr) = read_l_values(&st, &mut sh)?;
    let bad = |what: &str| EpsrError::Input(format!("sdf: cannot read {what}"));
    // central "molecule"
    let (_, ws) = words(&st, 11, 5);
    let mut lctyp1: Vec<i32> = Vec::new();
    let mut ic = 0;
    while ic < ws.len() && (lctyp1.len() as i32) < MCENTRE {
        let f = find_type(&st, &ws[ic]);
        ic += 1;
        if f <= 0 {
            log(&format!("\nsdf> Atom type {} does not exist in .ato file {}", ws[ic - 1], st.sym.fname_ato));
            return Ok(());
        }
        lctyp1.push(f);
    }
    if lctyp1.is_empty() {
        log("\nsdf> None of the specified atom types were recognised. \n - returning to main program");
        return Ok(());
    }
    let natofirst = lctyp1.len() as i32;
    let (_, ws) = words(&st, 12, 5);
    let npairfirstexpec = natofirst * (natofirst - 1) / 2;
    if (ws.len() as i32) < npairfirstexpec {
        log(&format!("\nsdf> Not enough pair distances defined for central molecule. Need: {:>5}", npairfirstexpec));
        return Ok(());
    }
    let mut distfirst = Vec::new();
    for i in 0..npairfirstexpec as usize {
        distfirst.push(read_real(&ws[i]).map_err(|_| bad("pairdist-c"))?);
    }
    let fractfirst = read_real(&st.sym.symbol_value[13]).map_err(|_| bad("fracdist-c"))?;
    let (nzt1, lztyp1) = read_axis(&st, 14, true)?;
    let (nyt1, lytyp1) = read_axis(&st, 15, true)?;
    // second "molecule"
    let (_, ws) = words(&st, 16, 5);
    let mut lctyp2: Vec<i32> = Vec::new();
    for w in &ws {
        let mut it = 0;
        let mut ifound = 0;
        while it < st.mol.ntype && ifound == 0 && (lctyp2.len() as i32) < MCENTRE {
            it += 1;
            if st.mol.type_[it].trim_end() == w {
                ifound = it;
            }
        }
        if ifound == 0 {
            log(&format!("\nsdf> Atom type {} does not exist in .ato file {}", w, st.sym.fname_ato));
            return Ok(());
        }
        lctyp2.push(ifound);
    }
    if lctyp2.is_empty() {
        log("\nsdf> None of the specified atom types were recognised. \n - returning to main program");
        return Ok(());
    }
    let natosecond = lctyp2.len() as i32;
    sh.ldistinct = if natosecond == natofirst { if lctyp1.iter().zip(&lctyp2).any(|(a, b)| a != b) { 1 } else { 0 } } else { 1 };
    let def;
    if sh.ldistinct != 0 {
        let (_, ws) = words(&st, 17, 5);
        let npairsecondexpec = natosecond * (natosecond - 1) / 2;
        if (ws.len() as i32) < npairsecondexpec {
            log(&format!("\nsdf> Not enough pair distances defined for second molecule. Need: {:>5}", npairsecondexpec));
            return Ok(());
        }
        let mut distsecond = Vec::new();
        for i in 0..npairsecondexpec as usize {
            distsecond.push(read_real(&ws[i]).map_err(|_| bad("pairdist-s"))?);
        }
        let fractsecond = read_real(&st.sym.symbol_value[18]).map_err(|_| bad("fracdist-s"))?;
        let (nzt2, lztyp2) = read_axis(&st, 19, false)?;
        let (nyt2, lytyp2) = read_axis(&st, 20, false)?;
        def = Def { lctyp1, lctyp2, distfirst, distsecond, fractfirst, fractsecond, nzt1, nyt1, nzt2, nyt2, lztyp1, lytyp1, lztyp2, lytyp2 };
    } else {
        def = Def {
            lctyp2,
            distsecond: distfirst.clone(),
            fractsecond: fractfirst,
            nzt2: nzt1.clone(),
            nyt2: nyt1.clone(),
            lztyp2: lztyp1.clone(),
            lytyp2: lytyp1.clone(),
            lctyp1,
            distfirst,
            fractfirst,
            nzt1,
            nyt1,
            lztyp1,
            lytyp1,
        };
    }
    // run_sdf does not include sharm_sdf.inc: the test angles it reads are
    // local variables, and the ones used by `sdf` (common block) stay zero.
    let (pt, _ang) = read_test(&st, 21, sh.piconv, ierr);
    sh.performtest = pt;
    sh.outputaveragemol = st.sym.symbol_value[22] == "y";

    let nmult = 5;
    let mut nsum = nsumt / nmult;
    let nt = nsum * nmult;
    if nt == nsumt {
        st.normalise_coordinates();
        let dmax = st.mol.dmax;
        sh.zmaxd = zmax / dmax;
        if sh.zmaxd > 0.5 {
            sh.zmaxd = 0.5;
            zmax = sh.zmaxd * dmax;
        }
        let _ = zmax;
        sh.get_l_and_n_values(&log);
        alloc_gofz(&mut sh);
        if !sdf(ctx, &mut st, &mut sh, &def, &log)? {
            // Fortran STOP
            return Ok(());
        }
        for i in 1..=sh.nz {
            sh.z[i] *= dmax;
        }
        sh.make_titles();
        let fnameinp = st.sym.fnameinp.clone();
        sh.outgofzs_sharm(&fnameinp, "h", &mut nsum, "   r", &log)?;
    }
    if !sh.performtest {
        nsumt += 1;
    }
    st.sym.symbol_value[4] = list_int(nsumt);
    st.sym.symbol_value[5] = list_int(sh.ndist);
    st.save_inp_file()?;
    if sh.performtest {
        log(" sdf> Test performed.");
    }
    Ok(())
}

/// `subroutine sdf(nout,fname_inp)`.  Returns false where the Fortran STOPs.
fn sdf(ctx: &ToolContext, st: &mut State, sh: &mut Sharm, d: &Def, log: &dyn Fn(&str)) -> Result<bool> {
    sh.pi = std::f32::consts::PI;
    sh.piconv = sh.pi / 180.0f32;
    sh.zstep = sh.zmaxd / (sh.nz as f32);
    sh.vol = st.lat.alatx * (st.lat.blaty * st.lat.clatz);
    sh.get_max_size_dlmn();
    let nmol = st.mol.nmol;
    sh.alloc_molecules(nmol);
    let natofirst = d.lctyp1.len() as i32;
    let natosecond = d.lctyp2.len() as i32;
    let npairfirst = d.distfirst.len() as i32;
    let npairsecond = d.distsecond.len() as i32;
    let nia = natofirst.max(natosecond) as usize;
    let mut iasave: F2<i32> = F2::new(nia, nmol.max(0) as usize);
    let mut nsmmo1 = 0;
    let mut nsmmo2 = 0;
    let mut imol1: i32 = 0;
    let mut ncentral = 0;
    let mut nsecond = 0;
    let mc = MCENTRE as usize;
    // common-block arrays (static)
    let mut atomno = vec![0i32; mc + 2];
    let mut xtemp = vec![0.0f32; mc + 1];
    let mut ytemp = vec![0.0f32; mc + 1];
    let mut ztemp = vec![0.0f32; mc + 1];
    let mut xfirstc = vec![0.0f32; mc + 1];
    let mut yfirstc = vec![0.0f32; mc + 1];
    let mut zfirstc = vec![0.0f32; mc + 1];
    let mut xsecondc = vec![0.0f32; mc + 1];
    let mut ysecondc = vec![0.0f32; mc + 1];
    let mut zsecondc = vec![0.0f32; mc + 1];
    let mut symfirst = vec![String::new(); mc + 1];
    let mut typfirst = vec![String::new(); mc + 1];
    let mut symsecond = vec![String::new(); mc + 1];
    let mut typsecond = vec![String::new(); mc + 1];
    let dmax = st.mol.dmax;
    let natom = st.mol.natom;
    for ia1 in 1..=natom {
        let mut iflag;
        let mut itmol = 0;
        while itmol <= sh.ldistinct {
            itmol += 1;
            let (lctyp, dist, fractcentre): (&[i32], &[f32], f32) = if itmol == 1 { (&d.lctyp1, &d.distfirst, d.fractfirst) } else { (&d.lctyp2, &d.distsecond, d.fractsecond) };
            let natomol = lctyp.len() as i32;
            let npairmol = if itmol == 1 { npairfirst } else { npairsecond };
            let distcentre: Vec<f32> = dist.iter().map(|x| x / dmax).collect();
            let m = &st.mol;
            if m.itype[ia1] != lctyp[0] {
                continue;
            }
            if imol1 < nmol {
                imol1 += 1;
            }
            sh.nm[imol1] = itmol;
            let mut nfound = 1;
            atomno[1] = ia1;
            let im1 = m.mref[ia1];
            let xo1 = m.xadd[im1] + m.xa[ia1];
            let yo1 = m.yadd[im1] + m.ya[ia1];
            let zo1 = m.zadd[im1] + m.za[ia1];
            xtemp[1] = 0.0;
            ytemp[1] = 0.0;
            ztemp[1] = 0.0;
            sh.xcent[imol1] = xo1;
            sh.ycent[imol1] = yo1;
            sh.zcent[imol1] = zo1;
            if natomol > nfound {
                for i in 2..=natomol as usize {
                    atomno[i] = 0;
                }
                let mut ia2 = 0;
                while ia2 < natom && nfound < natomol {
                    ia2 += 1;
                    let mut itest = 1;
                    while itest <= natomol && atomno[itest as usize] != ia2 {
                        itest += 1;
                    }
                    if itest > natomol {
                        let mut i = 1;
                        let mut ifound = 0;
                        while i < natomol && ifound == 0 {
                            i += 1;
                            if lctyp[(i - 1) as usize] == m.itype[ia2] && atomno[i as usize] == 0 {
                                let im2 = m.mref[ia2];
                                let xo2 = m.xadd[im2] + m.xa[ia2];
                                let yo2 = m.yadd[im2] + m.ya[ia2];
                                let zo2 = m.zadd[im2] + m.za[ia2];
                                let mut xdiff = xo2 - xo1;
                                let mut ydiff = yo2 - yo1;
                                let mut zdiff = zo2 - zo1;
                                crate::neighbour::minimagec(&mut xdiff, &mut ydiff, &mut zdiff);
                                let dmeas = (xdiff * xdiff + ydiff * ydiff + zdiff * zdiff).sqrt();
                                let dc = distcentre[(i - 2) as usize];
                                let dtestlow = dc * (1.0f32 - fractcentre);
                                let dtesthigh = dc * (1.0f32 + fractcentre);
                                if dmeas >= dtestlow && dmeas <= dtesthigh {
                                    ifound = 1;
                                    atomno[i as usize] = ia2;
                                    xtemp[i as usize] = xdiff;
                                    ytemp[i as usize] = ydiff;
                                    ztemp[i as usize] = zdiff;
                                    nfound += 1;
                                }
                            }
                        }
                    }
                }
                if nfound == natomol {
                    iflag = 1;
                    let mut iref = natomol - 1;
                    let mut i1 = 2usize;
                    let mut i2 = 3usize;
                    while iref < npairmol && iflag == 1 {
                        let mut xdiff = xtemp[i2] - xtemp[i1];
                        let mut ydiff = ytemp[i2] - ytemp[i1];
                        let mut zdiff = ytemp[i2] - ztemp[i1];
                        crate::neighbour::minimagec(&mut xdiff, &mut ydiff, &mut zdiff);
                        let dmeas = (xdiff * xdiff + ydiff * ydiff + zdiff * zdiff).sqrt();
                        iref += 1;
                        let dc = distcentre[(iref - 1) as usize];
                        let dtestlow = dc * (1.0f32 - fractcentre);
                        let dtesthigh = dc * (1.0f32 + fractcentre);
                        if dmeas < dtestlow || dmeas > dtesthigh {
                            iflag = 0;
                        }
                        i2 += 1;
                        if i2 as i32 > natomol {
                            i1 += 1;
                            i2 = i1 + 1;
                        }
                    }
                } else {
                    iflag = 0;
                }
            } else {
                iflag = 1;
            }
            if iflag == 1 {
                for i in 1..=natomol {
                    iasave[(i, imol1)] = atomno[i as usize];
                }
                if sh.nm[imol1] == 1 {
                    nsmmo1 += 1;
                }
                if sh.nm[imol1] == 2 {
                    nsmmo2 += 1;
                }
                let n = natomol as usize;
                let mut xat = F1::new(n);
                let mut yat = F1::new(n);
                let mut zat = F1::new(n);
                for ia in 1..=n {
                    xat[ia] = xtemp[ia];
                    yat[ia] = ytemp[ia];
                    zat[ia] = ztemp[ia];
                }
                let (p, c, ch) = if sh.nm[imol1] == 1 {
                    form_eulerangles(&xat, &yat, &zat, &d.nzt1, &d.nyt1, &d.lztyp1, &d.lytyp1)
                } else {
                    form_eulerangles(&xat, &yat, &zat, &d.nzt2, &d.nyt2, &d.lztyp2, &d.lytyp2)
                };
                sh.phim[imol1] = p;
                sh.cthetm[imol1] = c;
                sh.chim[imol1] = ch;
                if sh.performtest {
                    sh.phimt1[imol1] = sh.phim1;
                    sh.cthetmt1[imol1] = sh.thetam1.cos();
                    sh.chimt1[imol1] = sh.chim1;
                    sh.phimt2[imol1] = sh.phim2;
                    sh.cthetmt2[imol1] = sh.thetam2.cos();
                    sh.chimt2[imol1] = sh.chim2;
                }
                if sh.form_dlmn(imol1, log) != 0 {
                    log(&format!(" sdf> Error when calculating dlm coefficients for molecule {:>11}", imol1));
                    return Ok(true);
                }
                let m = &st.mol;
                for (kind, nato_k) in [(1, natofirst), (2, natosecond)] {
                    if sh.nm[imol1] != kind {
                        continue;
                    }
                    let cnt = if kind == 1 { &mut ncentral } else { &mut nsecond };
                    if !(*cnt == 0 || sh.outputaveragemol) {
                        continue;
                    }
                    *cnt += 1;
                    let first = *cnt == 1;
                    for ia in 1..=natomol as usize {
                        let (xt, yt, zt) = (xtemp[ia], ytemp[ia], ztemp[ia]);
                        let (xr, yr, zr) = if nato_k > 1 { rotateeuler(xt, yt, zt, -sh.chim[imol1], -sh.cthetm[imol1].acos(), -sh.phim[imol1]) } else { (xt, yt, zt) };
                        let (xc, yc, zc, sy, ty) = if kind == 1 {
                            (&mut xfirstc, &mut yfirstc, &mut zfirstc, &mut symfirst, &mut typfirst)
                        } else {
                            (&mut xsecondc, &mut ysecondc, &mut zsecondc, &mut symsecond, &mut typsecond)
                        };
                        if first {
                            xc[ia] = xr;
                            yc[ia] = yr;
                            zc[ia] = zr;
                            sy[ia] = m.atsymbol[m.itype[atomno[ia]]].clone();
                            ty[ia] = m.type_[m.itype[atomno[ia]]].clone();
                        } else {
                            xc[ia] += xr;
                            yc[ia] += yr;
                            zc[ia] += zr;
                        }
                    }
                }
            } else if imol1 > 0 {
                imol1 -= 1;
            }
        }
    }
    let mut nmolsharm = imol1;
    if sh.ldistinct > 0 {
        let mut imol1 = 0;
        while imol1 < nmolsharm {
            imol1 += 1;
            if sh.nm[imol1] != 1 {
                continue;
            }
            let mut imol2 = 0;
            while imol2 < nmolsharm {
                imol2 += 1;
                let mut iflag = 0;
                if imol1 != imol2 && sh.nm[imol2] == 2 {
                    let mut ia1 = 0;
                    while ia1 < natofirst && iflag == 0 {
                        ia1 += 1;
                        let mut ia2 = 0;
                        while ia2 < natosecond && iflag == 0 {
                            ia2 += 1;
                            if iasave[(ia1, imol1)] == iasave[(ia2, imol2)] {
                                iflag = 1;
                            }
                        }
                    }
                }
                if iflag == 1 {
                    log(&format!("sdf> Removing molecule: {:>5}", imol2));
                    let mut imol = imol2;
                    while imol < nmolsharm {
                        let imolt = imol + 1;
                        sh.nm[imol] = sh.nm[imolt];
                        sh.phim[imol] = sh.phim[imolt];
                        sh.cthetm[imol] = sh.cthetm[imolt];
                        sh.chim[imol] = sh.chim[imolt];
                        let na = if sh.nm[imol] == 1 { natofirst } else { natosecond };
                        for ia in 1..=na {
                            iasave[(ia, imol)] = iasave[(ia, imolt)];
                        }
                        imol = imolt;
                    }
                    nmolsharm -= 1;
                    nsmmo2 -= 1;
                    if imol1 >= imol2 {
                        imol1 -= 1;
                    }
                    imol2 -= 1;
                }
            }
        }
    }
    log(&format!(" sdf> Number of type 1 molecules = {:>10}\n sdf> Number of type 2 molecules = {:>10}", nsmmo1, nsmmo2));
    let fname_inp = st.sym.fname_inp.clone();
    let stem = fname_inp[..(indexlastdot(&fname_inp) - 1).max(0) as usize].to_string();
    let mk = |n: i32, xs: &[f32], ys: &[f32], zs: &[f32], sy: &[String], ty: &[String], cnt: i32| -> Vec<(String, f32, f32, f32, String, bool)> {
        (1..=n as usize)
            .map(|ia| {
                let c = cnt as f32;
                (sy[ia].clone(), xs[ia] * dmax / c, ys[ia] * dmax / c, zs[ia] * dmax / c, ty[ia].clone(), true)
            })
            .collect()
    };
    let f1 = format!("{stem}-first.xyz");
    if ncentral > 0 {
        write_mol_xyz(&ctx.path(&f1), &f1, &mk(natofirst, &xfirstc, &yfirstc, &zfirstc, &symfirst, &typfirst, ncentral))?;
    } else {
        OutFile::create(&ctx.path(&f1)).close()?;
        log(" sdf> Error writing first molecule file");
        return Ok(false);
    }
    let f2 = format!("{stem}-second.xyz");
    if sh.ldistinct == 0 && ncentral > 0 {
        write_mol_xyz(&ctx.path(&f2), &f2, &mk(natofirst, &xfirstc, &yfirstc, &zfirstc, &symfirst, &typfirst, ncentral))?;
    } else if sh.ldistinct == 1 && nsecond > 0 {
        write_mol_xyz(&ctx.path(&f2), &f2, &mk(natosecond, &xsecondc, &ysecondc, &zsecondc, &symsecond, &typsecond, nsecond))?;
    } else {
        OutFile::create(&ctx.path(&f2)).close()?;
        log(" sdf> Error writing second molecule file");
        return Ok(false);
    }
    sh.form_sum_prod_dlmn_cgcoeff(nmolsharm, log);
    let _ = ch2;
    Ok(true)
}

/// Program-specific branch of `check_symbol` (`process_symbols.f`, `symbol_setup.eq.'sdf'`).
/// Returns the Fortran `ierr` (0 = accepted).
pub fn check_symbol_extra(st: &mut State, _icount: i32, _ifilecount: i32, iref: i32) -> Result<i32> {
    let prompt = st.sym.prompt.clone();
    let mut msgs: Vec<String> = Vec::new();
    {
        let s = &mut st.sym;
        if iref == 11 || iref == 16 {
            let t = s.symbol_value[iref].clone();
            let w = parse(&t, 5);
            s.symbol_value[iref] = concatenate(1, w.nwords(), &w, &t, 256).0;
        } else if iref == 12 {
            let n = parse_n(&s.symbol_value[11], 5, 256).nwords();
            let nbondexpect = n * (n - 1) / 2;
            if parse_n(&s.symbol_value[12], 5, 256).nwords() < nbondexpect {
                msgs.push(format!("\n{} WARNING! Number of bonds specified is smaller than expected. Should be: {:>3}", prompt, nbondexpect));
            }
        } else if iref >= 16 {
            let naxisexpect = parse_n(&s.symbol_value[16], 5, 256).nwords();
            for k in [19, 20] {
                let t = s.symbol_value[k].clone();
                let w = parse_n(&t, 5, 256);
                if naxisexpect == 1 && w.nwords() > 1 {
                    msgs.push(format!("\n{} WARNING! No need for axis atoms in this case", prompt));
                    s.symbol_value[k] = w.word(&t, 1).to_string();
                }
            }
        }
        // (the Fortran's `else if(iref.eq.17)` branch is unreachable: iref >= 16 is caught above)
        if s.symbol_value[11] == s.symbol_value[16] {
            s.symbol_value[17] = s.symbol_value[12].clone();
            s.symbol_value[18] = s.symbol_value[13].clone();
            s.symbol_value[19] = s.symbol_value[14].clone();
            s.symbol_value[20] = s.symbol_value[15].clone();
        }
    }
    for m in msgs {
        st.log(&m);
    }
    let _ = read_int;
    Ok(0)
}
