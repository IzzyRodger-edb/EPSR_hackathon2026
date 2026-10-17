//! `makelattice` — translation of `fortran_src/makelattice.f` (and
//! `convertalphabetagamma` from `auxiliary_routines.f`).
//!
//! STATUS: translated; verified against the Fortran (EPSR-test-ex37).
//!
//! Invocation (same as the Fortran program):
//! `epsr-rs <folder> makelattice [<unit cell file> [na [nb [nc]]]]`
//! Missing values (file name, numbers of unit cells, output name) are read
//! from standard input.
//!
//! Unit-cell file layout:
//! ```text
//! a b c                      lattice constants (Å)
//! alpha beta gamma           cell angles (degrees)
//! [absolute|fractional]      optional coordinate-type line (default fractional)
//! n                          atoms in the unit cell
//! label x y z                n lines
//! label symbol [iexchange]   for each atom type (in order of first appearance):
//! eps sigma mass charge qradius
//! ecore dcore                optional
//! random-number state        optional (as at the end of a .ato file)
//! ```
//! Every atom becomes a one-atom molecule tethered to its lattice site; the
//! box is `na x nb x nc` unit cells (written as a non-cubic `.ato` file, since
//! the Fortran never sets `cubic`).

use super::makeato::setup_common::{command_text, eof_err, read_text};
use super::ToolContext;
use crate::ato::{cosd_intrinsic as cosd, sind_intrinsic as sind};
use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, FV};
use crate::fortran::listread::{read_int, ListRead, TextFile};
use crate::fortran::text::{findex, fstr, parse};
use crate::rng::{Ran1, NTAB};

/// `subroutine convertalphabetagamma(alpha,beta,gam,phib,thetac,phic)`.
pub fn convertalphabetagamma(alpha: f32, beta: f32, gam: f32, phib: &mut f32, thetac: &mut f32, phic: &mut f32) {
    let pi = 4.0f32 * 1.0f32.atan();
    let radconv = pi / 180.0;
    let degconv = 180.0 / pi;
    *phib = gam;
    let betar = beta * radconv;
    let cx = betar.cos();
    if cx.abs() < 1.0e-4 {
        *thetac = 0.0;
        *phic = 0.0;
    } else {
        let gamr = gam * radconv;
        let bx = gamr.cos();
        let by = gamr.sin();
        let alphar = alpha * radconv;
        let bdotc = alphar.cos();
        let cy = (bdotc - bx * cx) / by;
        let phicr = cy.atan2(cx);
        *phic = phicr * degconv;
        let cr = (cx * cx + cy * cy).sqrt();
        if cr.abs() <= 1.0 {
            let thetacr = cr.asin();
            *thetac = degconv * thetacr;
        }
    }
}

/// `initran1(10)`: read the generator state from the next record(s) of `f`.
pub fn initran1(f: &mut TextFile) -> Ran1 {
    let mut r = f.list();
    let st = (|| -> crate::fortran::IoRes<(i32, i32, i32, [i32; NTAB])> {
        let mcall = r.int()?.unwrap_or(0);
        let idum = r.int()?.unwrap_or(0);
        let iy = r.int()?.unwrap_or(0);
        let mut iv = [0i32; NTAB];
        for v in iv.iter_mut() {
            *v = r.int()?.unwrap_or(0);
        }
        Ok((mcall, idum, iy, iv))
    })();
    Ran1::from_record(st.ok())
}

/// Run the `makelattice` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    let (mut sometext, mut w) = command_text("makelattice", ctx);
    let mut nwordstest = 2;
    let readfile: String;
    if w.nwords() < 2 {
        nwordstest = 1;
        loop {
            st.log("\nmakelattice> Type name of file containing unit cell coordinates (incl. extension): ");
            sometext = read_text(ctx)?;
            w = parse(&sometext, 5);
            if w.nwords() > 0 {
                break;
            }
            st.log("\nmakelatice> No file specified! - try again!\n");
        }
        // readfile=sometext(ncf(1):ncl(1)); ilenreadfile=ncl(nwords)-ncf(1)+1
        let ilen = w.ncl[(w.nwords() - 1) as usize] - w.ncf[0] + 1;
        let first = w.word(&sometext, 1).to_string();
        readfile = crate::fortran::text::substr(&first, 1, ilen).trim_end().to_string();
    } else {
        readfile = w.word(&sometext, 2).to_string();
    }
    let (na, nb, nc);
    let rd5 = |ctx: &mut ToolContext, n: usize| -> Result<Vec<i32>> {
        let mut v = vec![0i32; n];
        let mut r = ctx.stdin.list();
        for x in v.iter_mut() {
            r.get_i(x).map_err(eof_err)?;
        }
        Ok(v)
    };
    if w.nwords() > nwordstest {
        na = read_int(w.word(&sometext, nwordstest + 1))?;
        if w.nwords() > nwordstest + 1 {
            nb = read_int(w.word(&sometext, nwordstest + 2))?;
            if w.nwords() > nwordstest + 2 {
                nc = read_int(w.word(&sometext, nwordstest + 3))?;
            } else {
                st.log("makelattice> Give number of unit cells along c-axis: ");
                nc = rd5(ctx, 1)?[0];
            }
        } else {
            st.log("makelattice> Give number of unit cells along b and c-axes: ");
            let v = rd5(ctx, 2)?;
            nb = v[0];
            nc = v[1];
        }
    } else {
        st.log("makelattice> Give number of unit cells along a, b and c-axes: ");
        let v = rd5(ctx, 3)?;
        na = v[0];
        nb = v[1];
        nc = v[2];
    }
    let fname = format!("{}{}", ctx.dir, readfile);
    let mut f = match TextFile::open(std::path::Path::new(&fname)) {
        Ok(f) => f,
        Err(_) => {
            st.log(&format!("makelattice> Problems with file: {}", fname));
            return Ok(());
        }
    };
    let (alpha, beta, gam);
    {
        let l = &mut st.lat;
        let mut r = f.list();
        r.get_r(&mut l.alat)?;
        r.get_r(&mut l.blat)?;
        r.get_r(&mut l.clat)?;
        drop(r);
        let mut r = f.list();
        let (mut a, mut b, mut g) = (0f32, 0f32, 0f32);
        r.get_r(&mut a)?;
        r.get_r(&mut b)?;
        r.get_r(&mut g)?;
        alpha = a;
        beta = b;
        gam = g;
        convertalphabetagamma(alpha, beta, gam, &mut l.phib, &mut l.thetac, &mut l.phic);
    }
    st.log(&format!(
        "makelattice> Lattice constants, angles and derived angles: {} {} {} / {} {} {} / {} {} {}",
        st.lat.alat, st.lat.blat, st.lat.clat, alpha, beta, gam, st.lat.phib, st.lat.thetac, st.lat.phic
    ));
    let mut absolute = false;
    let t = f.read_line()?;
    let nunitcell: i32 = if findex(&t, "absolute") > 0 {
        absolute = true;
        let mut v = 0;
        f.list().get_i(&mut v)?;
        v
    } else if findex(&t, "fractional") > 0 {
        let mut v = 0;
        f.list().get_i(&mut v)?;
        v
    } else {
        read_int(&t)?
    };
    let n = nunitcell.max(0) as usize;
    let mut atsymbolunit = vec![String::new(); n + 1];
    let mut xunit = vec![0f32; n + 1];
    let mut yunit = vec![0f32; n + 1];
    let mut zunit = vec![0f32; n + 1];
    st.mol.ntype = 0;
    st.mol.ensure_atoms(n.max(1));
    for i in 1..=nunitcell {
        let t = f.read_line()?;
        let w = parse(&t, 5);
        if w.nwords() > 3 {
            let iu = i as usize;
            atsymbolunit[iu] = fstr(w.word(&t, 1), 3);
            st.mol.alabel[i] = atsymbolunit[iu].clone();
            st.ftype(i);
            let s = w.span(&t, 2, 4);
            let mut r = ListRead::from_str(&s);
            r.get_r(&mut xunit[iu])?;
            r.get_r(&mut yunit[iu])?;
            r.get_r(&mut zunit[iu])?;
        } else {
            st.log("makelattice> Insufficient data supplied in input file.");
            return Ok(());
        }
    }
    st.log(&format!("makelattice> Number of atom types found: {:>3}", st.mol.ntype));
    let (mut epsit, mut sigt, mut amass, mut qt, mut rqt) = (0f32, 0f32, 0f32, 0f32, 0f32);
    for _it in 1..=st.mol.ntype {
        let t = f.read_line().unwrap_or_default();
        let w = parse(&t, 5);
        let a1 = if w.nwords() > 0 { fstr(w.word(&t, 1), 3) } else { String::new() };
        let ats = if w.nwords() > 1 {
            fstr(w.word(&t, 2), 3)
        } else {
            "XXX".to_string()
        };
        let iexch = if w.nwords() > 2 { read_int(w.word(&t, 3))? } else { -1 };
        {
            let mut r = f.list();
            let _ = (|| -> crate::fortran::IoRes<()> {
                r.get_r(&mut epsit)?;
                r.get_r(&mut sigt)?;
                r.get_r(&mut amass)?;
                r.get_r(&mut qt)?;
                r.get_r(&mut rqt)?;
                Ok(())
            })();
        }
        // natomsph is not declared in makelattice.f: implicit INTEGER, so ctype's
        // REAL dummy receives the bit pattern of the integer 1 (= 1.4e-45)
        let natomsph = f32::from_bits(1);
        let radminsph = 0.0;
        let radmaxsph = 0.5 * sigt;
        let pseudofacsph = 1.0;
        let rhosph = 0.0;
        st.ctype(
            &a1,
            &ats,
            iexch,
            epsit,
            sigt,
            amass,
            qt,
            rqt,
            natomsph,
            rhosph,
            radminsph,
            radmaxsph,
            pseudofacsph,
        );
    }
    {
        let mut r = f.list();
        let (mut e, mut d) = (st.mol.ecore, st.mol.dcore);
        let ok = r.get_r(&mut e).and_then(|_| r.get_r(&mut d));
        drop(r);
        if ok.is_ok() {
            st.mol.ecore = e;
            st.mol.dcore = d;
        } else {
            st.mol.ecore = 0.0;
            st.mol.dcore = 0.0;
        }
    }
    st.rng = initran1(&mut f);
    // output file name
    let mut text;
    loop {
        st.log("\nmakelattice> Type name of file to put lattice in: ");
        text = read_text(ctx)?;
        if parse(&text, 5).nwords() > 0 {
            break;
        }
        st.log("\nmakelatice> No file specified! - try again!\n");
    }
    let w = parse(&text, 5);
    let mut readfile = w.word(&text, 1).to_string();
    let ilenuptodot = findex(&readfile, ".") - 1;
    if ilenuptodot > 0 {
        readfile = readfile[..ilenuptodot as usize].to_string();
    }
    let fnameato = format!("{}{}.ato", ctx.dir, readfile);
    st.log(&format!("makelattice> Output being written to file: {}", fnameato));
    // ---- build the lattice
    {
        let l = &mut st.lat;
        l.alatx = l.alat;
        l.blatx = l.blat * cosd(l.phib);
        l.blaty = l.blat * sind(l.phib);
        l.clatx = l.clat * sind(l.thetac) * cosd(l.phic);
        l.claty = l.clat * sind(l.thetac) * sind(l.phic);
        l.clatz = l.clat * cosd(l.thetac);
    }
    let m = &mut st.mol;
    m.nmol = 0;
    m.simtk = 300.0;
    m.tolout = 0.1;
    m.stepmi = 0.01;
    m.stepri = 0.0;
    m.stepr = 0.001;
    m.stepm = 0.001;
    m.tempb = 65.0;
    m.tempangreq = 1.0;
    m.tempdihreq = 10.0;
    m.tempangwt = 1.0;
    m.tempdihwt = 10.0;
    let mut ia: i32 = 0;
    let mut im: i32 = 0;
    let (mut xav, mut yav, mut zav) = (0f32, 0f32, 0f32);
    m.natom = 0;
    m.ngrrot = 0;
    let acent = 0.5 * na as f32;
    let bcent = 0.5 * nb as f32;
    let ccent = 0.5 * nc as f32;
    let ntot = (na.max(0) as i64 * nb.max(0) as i64 * nc.max(0) as i64 * nunitcell.max(0) as i64) as usize;
    if ntot > crate::state::MMOL {
        return Err(EpsrError::Input(format!(
            "makelattice: {ntot} molecules exceed the maximum {}",
            crate::state::MMOL
        )));
    }
    m.ensure_atoms(ntot.max(n).max(1));
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
                for jm in 1..=nunitcell {
                    let ju = jm as usize;
                    im += 1;
                    m.nato[im] = 1;
                    if absolute {
                        m.xadd[im] = xlat + xunit[ju];
                        m.yadd[im] = ylat + yunit[ju];
                        m.zadd[im] = zlat + zunit[ju];
                    } else {
                        m.xadd[im] = xlat + xunit[ju] * l.alatx + yunit[ju] * l.blatx + zunit[ju] * l.clatx;
                        m.yadd[im] = ylat + yunit[ju] * l.blaty + zunit[ju] * l.claty;
                        m.zadd[im] = zlat + zunit[ju] * l.clatz;
                    }
                    xav += m.xadd[im];
                    yav += m.yadd[im];
                    zav += m.zadd[im];
                    m.phix[im] = 0.0;
                    m.phiy[im] = 0.0;
                    m.phiz[im] = 0.0;
                    m.tethered[im] = true;
                    m.naf[im] = ia + 1;
                    for _ib in 1..=m.nato[im] {
                        ia += 1;
                        m.alabel[ia] = atsymbolunit[ju].clone();
                        m.itype[ia] = m.itype[jm];
                        m.nbond[ia] = 0;
                        m.xa[ia] = 0.0;
                        m.ya[ia] = 0.0;
                        m.za[ia] = 0.0;
                    }
                    m.ngrot[im] = 0;
                    m.nal[im] = ia;
                }
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
        m.phix[im] = m.xadd[im];
        m.phiy[im] = m.yadd[im];
        m.phiz[im] = m.zadd[im];
    }
    // molecule types
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
            m.moltemplate[nt] = fmt("('moltype',i2.2)", &[FV::I(nt as i64)]);
            m.mtype[im] = nt;
            m.wholerot[nt] = m.stepr;
            m.wholetrans[nt] = m.stepm;
        }
    }
    let l = &mut st.lat;
    l.alat *= na as f32;
    l.blat *= nb as f32;
    l.clat *= nc as f32;
    st.moloutc(&fnameato)
}
