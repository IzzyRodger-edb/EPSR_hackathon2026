//! `sharm` — translation of `fortran_src/sharm_rev2.f` (+ `sharm_routines_rev1.f`).
//!
//! Spherical-harmonic expansion coefficients of the orientational pair
//! correlation function of two molecule types, from a `.ato` configuration.
//! The molecular centre is the mass-weighted centre of the listed atom types,
//! the molecular axes are defined by (averages of) atom positions.
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> sharm <input file>`
//! Keyword table: `crate::symtab_tools::SHARM`.
//!
//! Outputs: `<run>.SHARM.h01...` (coefficients, accumulated over configurations,
//! computed when `nsumt` is a multiple of 5), `<run>.SHARM-first.xyz`,
//! `<run>.SHARM-second.xyz` (averaged molecules), and the updated input file.
//!
//! Note: the Fortran writes the column titles with the run-time format
//! `(1h(i1,1h,,...)`, which gfortran >= 14 rejects ("Missing comma between
//! descriptors") — the original program then stops after writing the xyz files.
//! The port writes the titles `(l1,l2,l;n1,n2)` that older compilers produce.

#[path = "sharm_routines.rs"]
pub mod sharm_routines;

use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::listread::{read_int, read_real, ListRead};
use crate::fortran::text::{concatenate, indexlastdot, parse, parse_n};
use crate::fortran::{F1, F2};
use crate::state::State;
use sharm_routines::{form_eulerangles, list_int, rotateeuler, write_mol_xyz, Sharm};

pub(crate) const MCENTRE: i32 = 20;

/// Words of a symbol value as parsed by `parse(...,nspc,256,...)`.
pub(crate) fn words(st: &State, iref: i32, nspc: i32) -> (String, Vec<String>) {
    let v = st.sym.symbol_value[iref].clone();
    let w = parse_n(&v, nspc, 256);
    let ws = (1..=w.nwords()).map(|i| w.word(&v, i).to_string()).collect();
    (v, ws)
}

/// Type number of atom type `name` (exact match), 0 if absent.
pub(crate) fn find_type(st: &State, name: &str) -> i32 {
    let mut it = 0;
    while it < st.mol.ntype {
        it += 1;
        if st.mol.type_[it].trim_end() == name {
            return it;
        }
    }
    0
}

/// `character*2` assignment of the first word.
pub(crate) fn ch2(s: &str) -> String {
    crate::fortran::text::fstr(s, 2)
}

/// Read an axis definition: (axis name, atom numbers).
pub(crate) fn read_axis(st: &State, iref: i32, cap: bool) -> Result<(String, Vec<i32>)> {
    let (_, ws) = words(st, iref, 5);
    let name = ch2(ws.first().map(|s| s.as_str()).unwrap_or(""));
    let mut n = ws.len() as i32 - 1;
    if cap && n > MCENTRE {
        n = MCENTRE;
    }
    let mut v = Vec::new();
    for i in 1..=n {
        let x = read_int(&ws[i as usize]).map_err(|_| EpsrError::Input(format!("bad axis atom number in {:?}", st.sym.symbol_value[iref])))?;
        v.push(x);
    }
    Ok((name, v))
}

/// Read the 8 test angles: returns (performtest, angles in radians).
pub(crate) fn read_test(st: &State, iref: i32, piconv: f32, ierr_prev: i32) -> (bool, [f32; 8]) {
    let (v, ws) = words(st, iref, 5);
    let mut performtest = ws.first().map(|s| s.as_str()).unwrap_or("") != "notest" && ws.len() > 7;
    let mut a = [0.0f32; 8];
    let mut ierr = ierr_prev;
    if performtest {
        let w = parse_n(&v, 5, 256);
        let s = w.span(&v, 1, w.nwords());
        let mut r = ListRead::from_str(&s);
        ierr = 0;
        for x in a.iter_mut() {
            match r.real() {
                Ok(Some(y)) => *x = y,
                Ok(None) => {}
                Err(_) => {
                    ierr = 1;
                    break;
                }
            }
        }
        for x in a.iter_mut() {
            *x *= piconv;
        }
    }
    if ierr != 0 {
        performtest = false;
    }
    (performtest, a)
}

/// The common input of `sharm` and `sdf`: nr, rmax, nsumt, ncoeffs, l1/l2/l values, n steps.
pub(crate) fn read_l_values(st: &State, sh: &mut Sharm) -> Result<(f32, i32, i32)> {
    let bad = |what: &str| EpsrError::Input(format!("cannot read {what}"));
    sh.nz = read_int(&st.sym.symbol_value[2]).map_err(|_| bad("nr"))?;
    let zmax = read_real(&st.sym.symbol_value[3]).map_err(|_| bad("rmax"))?;
    let nsumt = read_int(&st.sym.symbol_value[4]).map_err(|_| bad("nsumt"))?;
    sh.ndist = read_int(&st.sym.symbol_value[5]).map_err(|_| bad("ncoeffs"))?;
    let mut ierr = 0;
    let rd = |iref: i32, nspc: i32, ierr: &mut i32| -> F1<i32> {
        let (_, ws) = words(st, iref, nspc);
        let mut v = F1::new(ws.len());
        for (k, w) in ws.iter().enumerate() {
            match read_int(w) {
                Ok(x) => {
                    v[k + 1] = x;
                    *ierr = 0;
                }
                Err(_) => *ierr = 1,
            }
        }
        v
    };
    sh.ll1v = rd(6, 1, &mut ierr);
    sh.nll1v = sh.ll1v.len() as i32;
    sh.ll2v = rd(7, 1, &mut ierr);
    sh.nll2v = sh.ll2v.len() as i32;
    let mut e2 = 0;
    sh.llv = rd(8, 5, &mut e2);
    if e2 != 0 {
        return Err(bad("lvalues"));
    }
    sh.nllv = sh.llv.len() as i32;
    sh.n1step = read_int(&st.sym.symbol_value[9]).map_err(|_| bad("n1step"))?;
    sh.n2step = read_int(&st.sym.symbol_value[10]).map_err(|_| bad("n2step"))?;
    Ok((zmax, nsumt, ierr))
}

/// Set up z, the gofz arrays and the column titles (`run_sharm`/`run_sdf`).
pub(crate) fn alloc_gofz(sh: &mut Sharm) {
    sh.mz = sh.nz;
    sh.mdist = sh.ndist;
    let (mz, md) = (sh.mz.max(0) as usize, sh.mdist.max(0) as usize);
    sh.z = F1::new(mz);
    sh.gofz = F2::new(mz, md);
    sh.gofzs = F2::new(mz, md);
    sh.gzt = F2::new(mz, md);
    sh.gzt2 = F2::new(mz, md);
    sh.zstep = sh.zmaxd / (sh.nz as f32);
    let zhalf = 0.5f32 * sh.zstep;
    for j in 1..=sh.nz {
        sh.z[j] = sh.zstep * (j as f32) - zhalf;
    }
}

/// Open and read the input file; returns false if the Fortran would return.
pub(crate) fn open(ctx: &ToolContext, st: &mut State, table: &crate::symbols::ToolTable) -> Result<bool> {
    st.install_tool_table(table);
    let name = ctx.input_name().ok_or_else(|| EpsrError::Input(format!("{}: no input file given", table.name)))?.to_string();
    let mut f = match st.setup_input_file(&name) {
        Ok(f) => f,
        Err(EpsrError::Io(_)) => return Ok(false),
        Err(e) => return Err(e),
    };
    let ierrout = st.read_input_file(&mut f)?;
    Ok(ierrout == 0)
}

/// Run the `sharm` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    if !open(ctx, &mut st, &crate::symtab_tools::SHARM)? {
        return Ok(());
    }
    let log = st.logger();
    let mut sh = Sharm { pi: std::f32::consts::PI, ..Default::default() };
    sh.piconv = sh.pi / 180.0f32;
    let (mut zmax, mut nsumt, ierr) = read_l_values(&st, &mut sh)?;
    // molecule 1 centre atom types
    let (_, ws) = words(&st, 11, 5);
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
    let (nzt1, lztyp1) = read_axis(&st, 12, true)?;
    let (nyt1, lytyp1) = read_axis(&st, 13, true)?;
    let (_, ws) = words(&st, 14, 5);
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
    sh.ldistinct = if lctyp2.len() == lctyp1.len() { if lctyp1.iter().zip(&lctyp2).any(|(a, b)| a != b) { 1 } else { 0 } } else { 1 };
    let (nzt2, lztyp2, nyt2, lytyp2);
    if sh.ldistinct != 0 {
        let a = read_axis(&st, 15, true)?;
        let b = read_axis(&st, 16, true)?;
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
    let (pt, ang) = read_test(&st, 17, sh.piconv, ierr);
    sh.performtest = pt;
    sh.phim1 = ang[0];
    sh.thetam1 = ang[1];
    sh.chim1 = ang[2];
    sh.phim2 = ang[3];
    sh.thetam2 = ang[4];
    sh.chim2 = ang[5];
    sh.thetallab = ang[6];
    sh.phillab = ang[7];
    sh.outputaveragemol = st.sym.symbol_value[18] == "y";

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
        let _ = zmax; // (not used further)
        sh.get_l_and_n_values(&log);
        alloc_gofz(&mut sh);
        let mol = Mols { lctyp1: &lctyp1, lctyp2: &lctyp2, nzt1: &nzt1, nyt1: &nyt1, lztyp1: &lztyp1, lytyp1: &lytyp1, nzt2: &nzt2, nyt2: &nyt2, lztyp2: &lztyp2, lytyp2: &lytyp2 };
        sharm(ctx, &mut st, &mut sh, &mol, &log)?;
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
        log(" sharm> Test performed.");
    }
    Ok(())
}

/// Molecule definitions passed to `sharm`.
pub(crate) struct Mols<'a> {
    pub lctyp1: &'a [i32],
    pub lctyp2: &'a [i32],
    pub nzt1: &'a str,
    pub nyt1: &'a str,
    pub lztyp1: &'a [i32],
    pub lytyp1: &'a [i32],
    pub nzt2: &'a str,
    pub nyt2: &'a str,
    pub lztyp2: &'a [i32],
    pub lytyp2: &'a [i32],
}

/// Does molecule `im` contain all the atom types `lctyp`?
pub(crate) fn has_types(st: &State, im: i32, lctyp: &[i32]) -> bool {
    let m = &st.mol;
    let mut found = true;
    for (k, &t) in lctyp.iter().enumerate() {
        if !found {
            break;
        }
        let mut fa = false;
        let mut ia1 = m.naf[im] - 1;
        while ia1 < m.nal[im] && !fa {
            ia1 += 1;
            fa = m.itype[ia1] == t;
        }
        if k == 0 {
            found = fa;
        } else {
            found = found && fa;
        }
    }
    found
}


/// The averaged-molecule xyz output of `sharm`/`sdfcube` for molecules with nm == kind.
/// Returns the averaged atoms (Å) if any.
pub(crate) fn average_molecule(st: &State, sh: &Sharm, kind: i32, mol: &Mols, rotate_single: bool) -> Option<Vec<(String, f32, f32, f32, String, bool)>> {
    let m = &st.mol;
    let mut ncentral = 0;
    let mut acc: Vec<(String, f32, f32, f32, String, bool)> = Vec::new();
    let mut im = 0;
    while im < m.nmol {
        im += 1;
        if sh.nm[im] != kind {
            continue;
        }
        if !(ncentral == 0 || sh.outputaveragemol) {
            continue;
        }
        ncentral += 1;
        if ncentral == 1 {
            acc = vec![(String::new(), 0.0, 0.0, 0.0, String::new(), false); m.nato[im] as usize];
        }
        let xshift = m.xadd[im] - sh.xcent[im];
        let yshift = m.yadd[im] - sh.ycent[im];
        let zshift = m.zadd[im] - sh.zcent[im];
        for ia in m.naf[im]..=m.nal[im] {
            let xt = m.xa[ia] + xshift;
            let yt = m.ya[ia] + yshift;
            let zt = m.za[ia] + zshift;
            let (xr, yr, zr) = if m.nato[im] > 1 || !rotate_single {
                rotateeuler(xt, yt, zt, -sh.chim[im], -sh.cthetm[im].acos(), -sh.phim[im])
            } else {
                (xt, yt, zt)
            };
            let iaref = ia - m.naf[im] + 1;
            let a = &mut acc[(iaref - 1) as usize];
            if ncentral == 1 {
                a.1 = xr;
                a.2 = yr;
                a.3 = zr;
                a.0 = m.atsymbol[m.itype[ia]].clone();
                a.4 = m.type_[m.itype[ia]].clone();
                let mut found = false;
                let mut ic = 0;
                while !found && ic < mol.lctyp1.len() {
                    found = m.itype[ia] == mol.lctyp1[ic];
                    ic += 1;
                }
                let mut iz = 0;
                while !found && iz < mol.lztyp1.len() {
                    found = iaref == mol.lztyp1[iz];
                    iz += 1;
                }
                let mut iy = 0;
                while !found && iy < mol.lytyp1.len() {
                    found = iaref == mol.lytyp1[iy];
                    iy += 1;
                }
                a.5 = found;
            } else {
                a.1 += xr;
                a.2 += yr;
                a.3 += zr;
            }
        }
    }
    if ncentral == 0 {
        return None;
    }
    let dmax = m.dmax;
    let nc = ncentral as f32;
    for a in acc.iter_mut() {
        a.1 = a.1 * dmax / nc;
        a.2 = a.2 * dmax / nc;
        a.3 = a.3 * dmax / nc;
    }
    Some(acc)
}

/// `subroutine sharm(nout,fname_inp)`.
fn sharm(ctx: &ToolContext, st: &mut State, sh: &mut Sharm, mol: &Mols, log: &dyn Fn(&str)) -> Result<()> {
    sh.pi = std::f32::consts::PI;
    sh.piconv = sh.pi / 180.0f32;
    sh.zstep = sh.zmaxd / (sh.nz as f32);
    sh.vol = st.lat.alatx * (st.lat.blaty * st.lat.clatz);
    sh.get_max_size_dlmn();
    let nmol = st.mol.nmol;
    sh.alloc_molecules(nmol);
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
        let lct: &[i32] = if sh.nm[imol1] == 1 { mol.lctyp1 } else { mol.lctyp2 };
        if sh.nm[imol1] != 0 && !lct.is_empty() {
            for ia1 in m.naf[imol1]..=m.nal[imol1] {
                for &t in lct {
                    if m.itype[ia1] == t {
                        let amass = m.amassa[m.itype[ia1]];
                        xo1 += m.xa[ia1] * amass;
                        yo1 += m.ya[ia1] * amass;
                        zo1 += m.za[ia1] * amass;
                        amasssum += amass;
                    }
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
            let n = m.nato[imol1] as usize;
            let mut xat = F1::new(n);
            let mut yat = F1::new(n);
            let mut zat = F1::new(n);
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
            if sh.performtest {
                sh.phimt1[imol1] = sh.phim1;
                sh.cthetmt1[imol1] = sh.thetam1.cos();
                sh.chimt1[imol1] = sh.chim1;
                sh.phimt2[imol1] = sh.phim2;
                sh.cthetmt2[imol1] = sh.thetam2.cos();
                sh.chimt2[imol1] = sh.chim2;
            }
            if sh.form_dlmn(imol1, log) != 0 {
                log(&format!(" sharm> Error when calculating dlm coefficients for molecule {:>11}", imol1));
                return Ok(());
            }
        } else {
            sh.nm[imol1] = 0;
        }
    }
    log(&format!(" sharm> Number of type 1 molecules = {:>10}\n sharm> Number of type 2 molecules = {:>10}", nsmmo1, nsmmo2));
    let fname_inp = st.sym.fname_inp.clone();
    let stem = &fname_inp[..(indexlastdot(&fname_inp) - 1).max(0) as usize];
    if let Some(acc) = average_molecule(st, sh, 1, mol, true) {
        let f = format!("{stem}-first.xyz");
        write_mol_xyz(&ctx.path(&f), &f, &acc)?;
        if sh.ldistinct == 0 {
            let f = format!("{stem}-second.xyz");
            write_mol_xyz(&ctx.path(&f), &f, &acc)?;
        }
    }
    if let Some(acc) = average_molecule(st, sh, 2, mol, true) {
        let f = format!("{stem}-second.xyz");
        write_mol_xyz(&ctx.path(&f), &f, &acc)?;
    }
    sh.form_sum_prod_dlmn_cgcoeff(nmol, log);
    Ok(())
}

/// Program-specific branch of `check_symbol` (`process_symbols.f`, `symbol_setup.eq.'sharm'`).
/// Returns the Fortran `ierr` (0 = accepted).
pub fn check_symbol_extra(st: &mut State, _icount: i32, _ifilecount: i32, iref: i32) -> Result<i32> {
    let s = &mut st.sym;
    if iref == 11 || iref == 14 {
        let t = s.symbol_value[iref].clone();
        let w = parse(&t, 5);
        s.symbol_value[iref] = concatenate(1, w.nwords(), &w, &t, 256).0;
    }
    if iref >= 14 && s.symbol_value[11] == s.symbol_value[14] {
        s.symbol_value[15] = s.symbol_value[12].clone();
        s.symbol_value[16] = s.symbol_value[13].clone();
    }
    Ok(0)
}
