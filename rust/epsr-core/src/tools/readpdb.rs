//! `readpdb` — translation of `fortran_src/readpdb.f` (`run_readpdb`,
//! `readpdb`, `upcase`).
//!
//! STATUS: translated; verified byte-identical against the Fortran built with
//! an explicit interface for `ctype` (see `readxyz.rs`: the original crashes
//! at the first `ctype` call for the same reason), test cases
//! EPSR-test-ex44/45/49.
//!
//! Invocation: `epsr-rs <folder> readpdb <file>.pdb`
//!
//! Reads `CRYST1`, `ATOM` and `CONECT` records (residue number = molecule;
//! `CONECT` records merge molecules, re-ordering atoms when they are not
//! contiguous), assigns elements from columns 13-14 of the atom name (unknown
//! elements become `Q` dummy atoms), then asks for the van der Waals factor
//! used to detect bonds (`vanderWaalsRadii.txt` must be in the current
//! directory) and whether to create `DIH` rotational groups, generates bonds,
//! bond-angle constraints and dihedral groups, and writes `<file>.ato`.
//!
//! Faithfully reproduced quirks of the Fortran: a molecule of a single atom
//! never gets `nal` set (it keeps the value left in the array, normally 0);
//! `dihedral` is called with absolute atom numbers although it expects
//! numbers relative to the molecule (correct only for molecule 1).

use super::readxyz::groupc_common as gc;
use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::{ListRead, TextFile};
use crate::fortran::text::{findex, fstr, indexlastdot, parse, substr};
use crate::fortran::{F1, F2, F3};
use crate::state::{State, NEAR};

/// Run the `readpdb` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    let (sometext, w) = gc::command_text("readpdb", &ctx.args);
    if w.nwords() < 2 {
        return crate::error::unsupported("readpdb without a file name (interactive file search)");
    }
    let mut readfile = w.word(&sometext, 2).to_string();
    let ilendot = indexlastdot(&readfile) - 1;
    if ilendot > 0 {
        readfile = readfile[..ilendot as usize].to_string();
    }
    let mut fnameato = gc::get_full_filename(&ctx.dir, &readfile);
    fnameato.push_str(".pdb");
    st.mol.ntype = 0;
    st.mol.nmoltype = 0;
    st.mol.sizefactor = 1.0;
    let mut dih = gc::DihData::default();
    let (rho, ierr) = readpdb(ctx, &mut st, &fnameato, &mut dih)?;
    if rho == 0.0 || ierr != 0 {
        st.log(&fmt("(/'readpdb> Problems with .pdb file ',a/'Routine will exit now'/)", &[FV::Ch(&fnameato, 256)]));
        return Ok(());
    }
    for it in 1..=st.mol.ntype {
        let s = st.mol.atsymbol[it].clone();
        if s == "XXX" || s.is_empty() {
            st.log(&fmt("(/'readpdb> Specify the chemical symbol for atom type ',a3,': ',$)", &[FV::Ch(&st.mol.type_[it], 3)]));
            let v = ctx.stdin.list().word()?.unwrap_or_default();
            st.mol.atsymbol[it] = fstr(&v, 3);
        }
    }
    for im in 1..=st.mol.nmol {
        gc::centreofmass(&mut st, im);
    }
    let ilendot = indexlastdot(&fnameato) - 1;
    let fnameato = format!("{}.ato", &fnameato[..ilendot.max(0) as usize]);
    gc::moloutc(&mut st, &fnameato, Some(&dih))
}

/// `upcase(string)`.
fn upcase(s: &str) -> String {
    s.to_ascii_uppercase()
}

/// List-directed character read of a field into a `character*len` variable.
fn read_char_field(field: &str, len: usize) -> Option<String> {
    match ListRead::from_str(field).word() {
        Ok(Some(v)) => Some(fstr(&v, len)),
        _ => None,
    }
}

/// `Fw.d` repeated reads `'(nfw.d)'` of the record `s` into `vals`; returns iostat.
fn read_fixed(s: &str, w: usize, d: usize, vals: &mut [f32]) -> i32 {
    for (k, v) in vals.iter_mut().enumerate() {
        let field = substr(s, (k * w + 1) as i32, ((k + 1) * w) as i32);
        match gc::read_fwd(&field, d) {
            Ok(x) => *v = x,
            Err(e) => return e.iostat(),
        }
    }
    0
}

/// `subroutine readpdb(readfile,rho,ierr)`.  Returns `(rho, ierr)`.
fn readpdb(ctx: &mut ToolContext, st: &mut State, readfile: &str, dih: &mut gc::DihData) -> Result<(f32, i32)> {
    let pi: f32 = 4.0f32 * 1.0f32.atan();
    let piconv: f32 = 180.0 / pi;
    let mut f = match TextFile::open(std::path::Path::new(readfile)) {
        Ok(f) => f,
        Err(_) => return Ok((0.0, 1)),
    };
    let mut ierr: i32 = 0;
    {
        let m = &mut st.mol;
        let l = &mut st.lat;
        m.sizefactor = 1.0;
        l.cubic = false;
        m.nmol = 0;
        m.simtk = 300.0;
        m.tempb = 0.65e2;
        m.tempang = 1.0;
        m.tempdih = 10.0;
        l.phib = 90.0;
        l.thetac = 0.0;
        l.phic = 0.0;
        m.tol = 0.0;
        m.tolout = 0.0;
        m.stepmi = 0.0;
        m.stepri = 0.0;
        m.stepr = 0.0;
        m.stepm = 0.0;
        m.tempangreq = 1.0;
        m.tempdihreq = 10.0;
        m.natom = 0;
        m.natomr = 0;
        m.ngrrot = 0;
        m.nrottot = 0;
        m.ntype = 0;
    }
    let mut ia: i32 = 0;
    let mut im: i32 = 0;
    let mut imsave: i32 = 0;
    const MATOM: i32 = 150000;
    let mmol = crate::state::MMOL as i32;
    // chsymboltemp(matom), ialabel(matom), mref: per-atom work arrays
    let mut chsymboltemp: F1<String> = F1::new(1024);
    let mut ialabel: F1<i32> = F1::new(1024);
    let (mut xat, mut yat, mut zat) = (0f32, 0f32, 0f32);
    let mut chsymbolsave = String::new();
    let mut alabelsave = String::new();
    st.mol.ensure_atoms(1024);
    while ierr == 0 && ia < MATOM && im < mmol {
        let sometext = match f.read_line() {
            Ok(s) => s,
            Err(_) => {
                ierr = -1;
                continue;
            }
        };
        let rec6 = substr(&sometext, 1, 6);
        if findex(&rec6, "CRYST1") > 0 {
            let mut abc = [st.lat.alat, st.lat.blat, st.lat.clat];
            let _ = read_fixed(&substr(&sometext, 7, 33), 9, 3, &mut abc);
            let mut ang = [0f32; 3];
            ierr = read_fixed(&substr(&sometext, 34, 54), 7, 2, &mut ang);
            while ierr != 0 {
                st.log(" readpdb> Need a,b,c axes and 3 angles (alpha,beta,gamma):");
                let mut r = ctx.stdin.list();
                let res = (|| -> crate::fortran::IoRes<()> {
                    r.get_r(&mut abc[0])?;
                    r.get_r(&mut abc[1])?;
                    r.get_r(&mut abc[2])?;
                    r.get_r(&mut ang[0])?;
                    r.get_r(&mut ang[1])?;
                    r.get_r(&mut ang[2])?;
                    Ok(())
                })();
                ierr = crate::fortran::listread::iostat(&res);
                if ierr < 0 {
                    return Err(EpsrError::Input("readpdb: end of input while reading the cell (a,b,c,alpha,beta,gamma)".into()));
                }
            }
            st.lat.alat = abc[0];
            st.lat.blat = abc[1];
            st.lat.clat = abc[2];
            let (mut phib, mut thetac, mut phic) = (st.lat.phib, st.lat.thetac, st.lat.phic);
            let log = st.logger();
            gc::convertalphabetagamma(ang[0], ang[1], ang[2], &mut phib, &mut thetac, &mut phic, &log);
            st.lat.phib = phib;
            st.lat.thetac = thetac;
            st.lat.phic = phic;
            gc::setup_lattice(st);
        } else if findex(&rec6, "ATOM") > 0 {
            // each read sets iostat; only the last one decides
            match ListRead::from_str(&substr(&sometext, 7, 11)).int() {
                Ok(Some(v)) => ia = v,
                _ => {}
            }
            if let Some(v) = read_char_field(&substr(&sometext, 13, 14), 2) {
                chsymbolsave = v;
            }
            if let Some(v) = read_char_field(&substr(&sometext, 13, 16), 4) {
                alabelsave = v;
            }
            match ListRead::from_str(&substr(&sometext, 23, 26)).int() {
                Ok(Some(v)) => im = v,
                _ => {}
            }
            let mut xyz = [xat, yat, zat];
            ierr = read_fixed(&substr(&sometext, 31, 54), 8, 3, &mut xyz);
            xat = xyz[0];
            yat = xyz[1];
            zat = xyz[2];
            st.log(&fmt("(i5,a2,i5,a4)", &[FV::I(ia as i64), FV::Ch(&chsymbolsave, 2), FV::I(im as i64), FV::Ch(&alabelsave, 4)]));
            if ierr == 0 {
                let m = &mut st.mol;
                m.natom += 1;
                let natom = m.natom;
                m.ensure_atoms(natom as usize);
                if chsymboltemp.len() < m.xa.len() {
                    chsymboltemp.grow(m.xa.len());
                    ialabel.grow(m.xa.len());
                }
                m.nbond[natom] = 0;
                m.imageatom[natom] = 0;
                ialabel[natom] = ia;
                chsymboltemp[natom] = chsymbolsave.clone();
                m.alabel[natom] = fstr(&alabelsave, 3);
                m.xa[natom] = xat;
                m.ya[natom] = yat;
                m.za[natom] = zat;
                if im != imsave {
                    m.nmol += 1;
                    let nm = m.nmol;
                    m.naf[nm] = natom;
                    m.xadd[nm] = 0.0;
                    m.yadd[nm] = 0.0;
                    m.zadd[nm] = 0.0;
                    m.nato[nm] = 1;
                    m.ngrot[nm] = 0;
                    m.ngrrot += m.ngrot[nm];
                    m.tethered[nm] = false;
                    m.phix[nm] = 0.0;
                    m.phiy[nm] = 0.0;
                    m.phiz[nm] = 0.0;
                } else {
                    let nm = m.nmol;
                    m.nal[nm] = natom;
                    m.nato[nm] += 1;
                }
                m.mref[natom] = m.nmol;
                imsave = im;
            }
        } else if findex(&rec6, "CONECT") > 0 && sometext.trim_end().len() > 15 {
            conect(st, &sometext, &mut ialabel, &mut chsymboltemp)?;
        }
    }
    let m = &mut st.mol;
    if m.natom > 0 {
        ierr = 0;
    }
    for im in 1..=m.nmol {
        let mut xsum: f32 = 0.0;
        let mut ysum: f32 = 0.0;
        let mut zsum: f32 = 0.0;
        for ib in m.naf[im]..=m.nal[im] {
            xsum += m.xa[ib];
            ysum += m.ya[ib];
            zsum += m.za[ib];
        }
        m.xadd[im] = xsum / m.nato[im] as f32;
        m.yadd[im] = ysum / m.nato[im] as f32;
        m.zadd[im] = zsum / m.nato[im] as f32;
        for ib in m.naf[im]..=m.nal[im] {
            m.xa[ib] -= m.xadd[im];
            m.ya[ib] -= m.yadd[im];
            m.za[ib] -= m.zadd[im];
        }
    }
    let natom = m.natom;
    let mut iatomicno: F1<i32> = F1::new(natom.max(1) as usize);
    m.natomr = 0;
    for ia in 1..=natom {
        let mut ir = 0;
        let mut found = false;
        while !found && ir < gc::MELEMENTS {
            ir += 1;
            found = chsymboltemp[ia] == upcase(gc::chsymbol(ir));
            if found {
                iatomicno[ia] = ir;
            }
        }
        if !found {
            iatomicno[ia] = 8;
            let lab = st.mol.alabel[ia].clone();
            let rest: String = substr(&lab, 2, 3);
            st.mol.alabel[ia] = format!("Q{}", rest).trim_end().to_string();
        }
        st.ftype(ia);
        if st.mol.qtype[ia] != 1 {
            st.mol.natomr += 1;
        }
    }
    let rho = st.mol.natomr as f32 / st.mol.volume;
    {
        let mut w = ListWriter::new();
        w.chars(" readpdb> rho = ").real(rho).int(st.mol.natomr).int(st.mol.natom).int(st.mol.nmol);
        st.log(&w.finish());
    }
    {
        let m = &mut st.mol;
        m.ntypefile = 0;
        let mut ia = 0;
        while ia < m.natom {
            ia += 1;
            let mut itypefile = 0;
            let mut ifound = 0;
            while itypefile < m.ntypefile && ifound == 0 {
                itypefile += 1;
                if m.alabel[ia] == m.typefile[itypefile] {
                    ifound = itypefile;
                }
            }
            if ifound == 0 {
                m.ntypefile += 1;
                let n = m.ntypefile;
                m.typefile[n] = m.alabel[ia].clone();
                m.iatomicnotype[n] = iatomicno[ia];
            }
        }
    }
    super::readxyz::set_type_parameters(st);
    {
        let m = &mut st.mol;
        let mut im = 0;
        m.nmoltype = 0;
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
                m.moltemplate[nt] = fmt("(\"moltype\",i2.2)", &[FV::I(nt as i64)]);
                m.mtype[im] = nt;
                m.interntrans[nt] = m.stepmi;
                m.internrot[nt] = m.stepri;
                m.wholerot[nt] = m.stepr;
                m.wholetrans[nt] = m.stepm;
            }
        }
    }
    st.log(&fmt("(/' readpdb> Following molecule types defined initially:-'/)", &[]));
    // ---- van der Waals radii (read from the current directory)
    let mut vdw = [0f32; 104];
    let mut fv = match TextFile::open(std::path::Path::new("vanderWaalsRadii.txt")) {
        Ok(f) => f,
        Err(_) => {
            st.log(" readpdb> Unable to open vanderWaalsRadii.txt");
            return Ok((rho, 1));
        }
    };
    let mut icount = 0;
    while let Ok(sometext) = fv.read_line() {
        let w = parse(&sometext, 11);
        if w.nwords() > 2 {
            let ir = ListRead::from_str(w.word(&sometext, 1)).int()?.unwrap_or(0);
            if ir > 0 && ir <= gc::MELEMENTS {
                if let Some(v) = ListRead::from_str(w.word(&sometext, 3)).real()? {
                    vdw[ir as usize] = v;
                }
                if ir == 1 {
                    vdw[ir as usize] = 0.0;
                }
                icount += 1;
            }
        }
    }
    st.log(&fmt(
        "(/' readpdb> Number of van der Waals radii read in: ',i3//' Give van der Waals factor (suggest 0.7). This is used to test whether two atoms are bonded or not: ')",
        &[FV::I(icount as i64)],
    ));
    let mut vander_waals_factor: f32 = 0.0;
    let ierrv = {
        let mut r = ctx.stdin.list();
        let res = r.get_r(&mut vander_waals_factor);
        crate::fortran::listread::iostat(&res)
    };
    if ierrv != 0 || vander_waals_factor <= 0.0 {
        vander_waals_factor = 0.7;
    }
    st.log("\n readpdb> Create dihedral angles (y or no): ");
    let yesorno = fstr(&ctx.stdin.read_line()?, 3);
    let dodihedrals = findex(&yesorno, "y") > 0;
    st.normalise_coordinates();
    let near = NEAR as i32;
    let natom = st.mol.natom;
    // `dihedral` is called with absolute atom numbers (see the module notes):
    // make the (zero) storage beyond natom addressable as in the Fortran.
    st.mol.ensure_atoms(2 * natom as usize + 2);
    // bonds from van der Waals radii
    let dmax = st.mol.dmax;
    let m = &mut st.mol;
    for im in 1..=m.nmol {
        for ia in m.naf[im]..=m.nal[im] - 1 {
            for ja in ia + 1..=m.nal[im] {
                let xdif = m.xa[ia] - m.xa[ja];
                let ydif = m.ya[ia] - m.ya[ja];
                let zdif = m.za[ia] - m.za[ja];
                let dtest = dmax * (xdif * xdif + ydif * ydif + zdif * zdif).sqrt();
                let dvanderwaals =
                    vander_waals_factor * (vdw[m.iatomicnotype[m.itype[ia]] as usize] + vdw[m.iatomicnotype[m.itype[ja]] as usize]);
                if dtest < dvanderwaals && m.nbond[ia] < near && m.nbond[ja] < near {
                    m.nbond[ia] += 1;
                    let nb = m.nbond[ia];
                    m.ensure_bonds(nb as usize);
                    m.iabond[(ia, nb)] = ja;
                    m.tbond[(ia, nb)] = dtest / dmax;
                    m.nbond[ja] += 1;
                    let nb = m.nbond[ja];
                    m.ensure_bonds(nb as usize);
                    m.iabond[(ja, nb)] = ia;
                    m.tbond[(ja, nb)] = dtest / dmax;
                }
            }
        }
    }
    // bond angles
    let mut nbondp: F1<i32> = F1::new(m.xa.len().max(natom as usize + 1));
    for ia in 1..=natom {
        nbondp[ia] = m.nbond[ia];
    }
    let mut ir: i32 = 0;
    for im in 1..=st.mol.nmol {
        {
            let m = &mut st.mol;
            for ia in m.naf[im]..=m.nal[im] {
                let mut ib1 = 0;
                while ib1 < nbondp[ia] - 1 {
                    ib1 += 1;
                    let iaref1 = m.iabond[(ia, ib1)];
                    let mut ib2 = ib1;
                    while ib2 < nbondp[ia] {
                        ib2 += 1;
                        let iaref2 = m.iabond[(ia, ib2)];
                        let mut ic = 0;
                        let mut found = false;
                        while !found && ic < m.nbond[iaref1] {
                            ic += 1;
                            found = iaref2 == m.iabond[(iaref1, ic)];
                        }
                        if !found {
                            let xdif = m.xa[iaref1] - m.xa[iaref2];
                            let ydif = m.ya[iaref1] - m.ya[iaref2];
                            let zdif = m.za[iaref1] - m.za[iaref2];
                            let dtest = (xdif * xdif + ydif * ydif + zdif * zdif).sqrt();
                            if m.nbond[iaref1] < near && m.nbond[iaref2] < near {
                                m.nbond[iaref1] += 1;
                                m.nbond[iaref2] += 1;
                                let n1 = m.nbond[iaref1];
                                let n2 = m.nbond[iaref2];
                                m.ensure_bonds(n1.max(n2) as usize);
                                m.iabond[(iaref1, n1)] = iaref2;
                                m.iabond[(iaref2, n2)] = iaref1;
                                m.tbond[(iaref1, n1)] = -dtest;
                                m.tbond[(iaref2, n2)] = -dtest;
                            }
                        }
                    }
                }
            }
        }
        if dodihedrals && st.mol.nato[im] > 3 {
            make_dihedrals(st, dih, &nbondp, im, &mut ir, piconv)?;
        }
        let m = &mut st.mol;
        m.ngrot[im] = ir - m.ngrrot;
        m.ngrrot = ir;
    }
    st.unnormalise_coordinates();
    // unnormalise_coordinates: dihedral angles back to degrees
    let piconv2: f32 = 180.0 / pi;
    let m = &st.mol;
    let mut ir = 0;
    while ir < m.ngrrot {
        ir += 1;
        if m.rotgrp[ir] == "DIH" {
            for k in 1..=m.ndihedralang[ir] {
                dih.dihedralang[(k, ir)] *= piconv2;
            }
        }
    }
    Ok((rho, ierr))
}

/// The `CONECT` branch of `readpdb`: merge the molecules of bonded atoms.
fn conect(st: &mut State, sometext: &str, ialabel: &mut F1<i32>, chsymboltemp: &mut F1<String>) -> Result<()> {
    let m = &mut st.mol;
    let natom = m.natom;
    let nsecond = (sometext.trim_end().len() as i32 - 6) / 5 - 1;
    if nsecond > 4 {
        return crate::error::unsupported("readpdb: CONECT record with more than 4 bonded atoms (isecond has dimension 4)");
    }
    let mut ifirst = ListRead::from_str(&substr(sometext, 7, 11)).int()?.ok_or_else(|| EpsrError::Input("readpdb: bad CONECT record".into()))?;
    let lookup = |lab: i32| -> Option<i32> {
        let mut ia = 0;
        while ia < natom {
            ia += 1;
            if ialabel[ia] == lab {
                return Some(ia);
            }
        }
        None
    };
    if let Some(ia) = lookup(ifirst) {
        ifirst = ia;
    }
    let mut isecond = [0i32; 5];
    let mut ir1 = 12;
    for i in 1..=nsecond {
        let ir2 = ir1 + 4;
        let mut v = ListRead::from_str(&substr(sometext, ir1, ir2)).int()?.ok_or_else(|| EpsrError::Input("readpdb: bad CONECT record".into()))?;
        if let Some(ia) = lookup(v) {
            v = ia;
        }
        isecond[i as usize] = v;
        ir1 = ir2 + 1;
    }
    let mut i = 1;
    while i <= nsecond && ifirst <= natom && isecond[i as usize] <= natom {
        let is = isecond[i as usize];
        if m.mref[ifirst] != m.mref[is] {
            let immin = m.mref[ifirst].min(m.mref[is]);
            let immax = m.mref[ifirst].max(m.mref[is]);
            if m.naf[immax] != m.nal[immin] + 1 {
                let n = m.nato[immax];
                let mut nbondmove = Vec::new();
                let mut imageatommove = Vec::new();
                let mut ialabelmove = Vec::new();
                let mut chsymbolmove = Vec::new();
                let mut alabelmove = Vec::new();
                let mut xamove = Vec::new();
                let mut yamove = Vec::new();
                let mut zamove = Vec::new();
                for ia in 1..=n {
                    let iaref = ia + m.naf[immax] - 1;
                    nbondmove.push(m.nbond[iaref]);
                    imageatommove.push(m.imageatom[iaref]);
                    ialabelmove.push(ialabel[iaref]);
                    chsymbolmove.push(chsymboltemp[iaref].clone());
                    alabelmove.push(m.alabel[iaref].clone());
                    xamove.push(m.xa[iaref]);
                    yamove.push(m.ya[iaref]);
                    zamove.push(m.za[iaref]);
                }
                let mut ia = m.naf[immax] - 1;
                let mut mrefsave = m.mref[ia];
                m.nal[mrefsave] = ia + n;
                while ia > m.nal[immin] {
                    let iaref = ia + n;
                    m.nbond[iaref] = m.nbond[ia];
                    m.imageatom[iaref] = m.imageatom[ia];
                    ialabel[iaref] = ialabel[ia];
                    chsymboltemp[iaref] = chsymboltemp[ia].clone();
                    m.alabel[iaref] = m.alabel[ia].clone();
                    m.xa[iaref] = m.xa[ia];
                    m.ya[iaref] = m.ya[ia];
                    m.za[iaref] = m.za[ia];
                    if m.mref[ia] != mrefsave {
                        mrefsave = m.mref[ia];
                        m.nal[mrefsave] = iaref;
                    } else {
                        m.naf[mrefsave] = iaref;
                    }
                    m.mref[iaref] = mrefsave;
                    ia -= 1;
                }
                for ia in 1..=n {
                    let iaref = ia + m.nal[immin];
                    let k = (ia - 1) as usize;
                    m.nbond[iaref] = nbondmove[k];
                    m.imageatom[iaref] = imageatommove[k];
                    ialabel[iaref] = ialabelmove[k];
                    chsymboltemp[iaref] = chsymbolmove[k].clone();
                    m.alabel[iaref] = alabelmove[k].clone();
                    m.xa[iaref] = xamove[k];
                    m.ya[iaref] = yamove[k];
                    m.za[iaref] = zamove[k];
                    m.mref[iaref] = immin;
                }
            }
            m.nato[immin] += m.nato[immax];
            m.nal[immin] += m.nato[immax];
            for ia in m.naf[immin]..=m.nal[immin] {
                m.mref[ia] = immin;
            }
            let mut im = immax;
            while im < m.nmol {
                let im1 = im + 1;
                m.nato[im] = m.nato[im1];
                m.naf[im] = m.naf[im1];
                m.nal[im] = m.nal[im1];
                for ia in m.naf[im]..=m.nal[im] {
                    m.mref[ia] = im;
                }
                im += 1;
            }
            m.nmol -= 1;
        }
        i += 1;
    }
    Ok(())
}

/// The dihedral-group part of the molecule loop of `readpdb`.
fn make_dihedrals(st: &mut State, dih: &mut gc::DihData, nbondp: &F1<i32>, im: i32, ir: &mut i32, piconv: f32) -> Result<()> {
    let near = NEAR as i32;
    let (naf, nal) = (st.mol.naf[im], st.mol.nal[im]);
    for iaref1 in naf..=nal {
        let mut ib1 = 0;
        while ib1 < nbondp[iaref1] {
            ib1 += 1;
            let iaref2 = st.mol.iabond[(iaref1, ib1)];
            let ir1 = iaref1;
            let ir2 = iaref2;
            let mut iref = 0;
            let mut ifounddih = 0;
            while iref < *ir && ifounddih == 0 {
                iref += 1;
                let m = &st.mol;
                if m.rotgrp[iref] == "DIH"
                    && ((ir1 == m.labaxi[(iref, 1)] && ir2 == m.labaxi[(iref, 2)]) || (ir2 == m.labaxi[(iref, 1)] && ir1 == m.labaxi[(iref, 2)]))
                {
                    ifounddih = iref;
                }
            }
            let mut ib2 = 0;
            while ib2 < nbondp[iaref1] {
                ib2 += 1;
                let ibref1 = st.mol.iabond[(iaref1, ib2)];
                let mut ib3 = 0;
                while ib3 < nbondp[iaref2] && ibref1 != iaref2 {
                    ib3 += 1;
                    let ibref2 = st.mol.iabond[(iaref2, ib3)];
                    let mut ic = 0;
                    let mut found = ibref1 == ibref2 || ibref2 == iaref1;
                    while ic < st.mol.nbond[ibref2] && !found {
                        ic += 1;
                        found = ibref1 == st.mol.iabond[(ibref2, ic)];
                    }
                    if ifounddih > 0 && !found {
                        let mut iref = 0;
                        while iref < st.mol.ndihedralang[ifounddih] && !found {
                            iref += 1;
                            found = (ibref1 == dih.labrotdih[(1, iref, ifounddih)] && ibref2 == dih.labrotdih[(2, iref, ifounddih)])
                                || (ibref2 == dih.labrotdih[(1, iref, ifounddih)] && ibref1 == dih.labrotdih[(2, iref, ifounddih)]);
                        }
                    }
                    if !found {
                        let mut dread = gc::dihedral(st, im, iaref1, iaref2, ibref1, ibref2);
                        {
                            let mut w = ListWriter::new();
                            w.int(ifounddih).int(iaref1).int(iaref2).int(ibref1).int(ibref2).real(dread);
                            st.log(&w.finish());
                        }
                        dread /= piconv;
                        if ifounddih == 0 {
                            let labrott1 = gc::get_rotgrp(st, nbondp, ir1, ir2);
                            let labrott2 = gc::get_rotgrp(st, nbondp, ir2, ir1);
                            let nrott1 = labrott1.len() as i32;
                            let nrott2 = labrott2.len() as i32;
                            *ir += 1;
                            let irr = *ir;
                            ensure_dih(dih, irr);
                            st.mol.ensure_rotgrps(irr as usize);
                            let m = &mut st.mol;
                            m.ndihedralang[irr] = 1;
                            dih.dihedralang[(m.ndihedralang[irr], irr)] = dread;
                            m.rotgrp[irr] = "DIH".to_string();
                            ifounddih = irr;
                            let (a1, a2, lab, l1, l2) = if nrott1 < nrott2 {
                                (ir1, ir2, &labrott1, ibref1, ibref2)
                            } else {
                                (ir2, ir1, &labrott2, ibref2, ibref1)
                            };
                            m.labaxi[(irr, 1)] = a1;
                            m.labaxi[(irr, 2)] = a2;
                            m.nrot[irr] = lab.len() as i32;
                            m.rotf[irr] = m.nrottot + 1;
                            m.nrottot += m.nrot[irr];
                            m.rotl[irr] = m.nrottot;
                            let nr = m.nrottot as usize;
                            m.ensure_labrot(nr);
                            for i in 1..=m.nrot[irr] {
                                let labref = m.rotf[irr] + i - 1;
                                m.labrot[labref] = lab[(i - 1) as usize];
                            }
                            let nd = m.ndihedralang[irr];
                            dih.labrotdih[(1, nd, irr)] = l1;
                            dih.labrotdih[(2, nd, irr)] = l2;
                        } else if st.mol.ndihedralang[ifounddih] < near {
                            let m = &mut st.mol;
                            m.ndihedralang[ifounddih] += 1;
                            let nd = m.ndihedralang[ifounddih];
                            dih.dihedralang[(nd, ifounddih)] = dread;
                            if ir1 == m.labaxi[(ifounddih, 1)] && ir2 == m.labaxi[(ifounddih, 2)] {
                                dih.labrotdih[(1, nd, ifounddih)] = ibref1;
                                dih.labrotdih[(2, nd, ifounddih)] = ibref2;
                            } else if ir2 == m.labaxi[(ifounddih, 1)] && ir1 == m.labaxi[(ifounddih, 2)] {
                                dih.labrotdih[(1, nd, ifounddih)] = ibref2;
                                dih.labrotdih[(2, nd, ifounddih)] = ibref1;
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

/// Grow the `DIH` arrays to hold group `ir` (`near` angles per group).
fn ensure_dih(dih: &mut gc::DihData, ir: i32) {
    let near = NEAR;
    let need = ir as usize;
    if dih.dihedralang.n2 < need {
        let n = need.max(2 * dih.dihedralang.n2).max(16);
        let mut d2: F2<f32> = F2::new(near, n);
        let mut l3: F3<i32> = F3::new(2, near, n);
        for g in 1..=dih.dihedralang.n2 {
            for k in 1..=near {
                d2[(k, g)] = dih.dihedralang[(k, g)];
                l3[(1, k, g)] = dih.labrotdih[(1, k, g)];
                l3[(2, k, g)] = dih.labrotdih[(2, k, g)];
            }
        }
        dih.dihedralang = d2;
        dih.labrotdih = l3;
    }
}
