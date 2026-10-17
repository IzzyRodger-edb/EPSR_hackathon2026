//! `makemole` — translation of `fortran_src/makemole.f` → subroutine
//! `makemole` (`ato_manipulation_rev4.f`) with `get_rotgrp` and
//! `get_whole_rotgrp`.
//!
//! STATUS: translated for `.gmol` molecule files (the type used in
//! EPSR-test) and the original column-map `.mol` format (with its `.atm`
//! file), verified byte-identical (EPSR-test-ex28, ex29).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> makemole <file>.mol [< answers]`
//!
//! Reads the molecule template (atoms with coordinates and bonded neighbours,
//! `bond`, `angle`, `extra`, `rot`, `wholerot`, `introt`, `changelabel`,
//! `potential`, `qradius`, `temperature`, `vibtemp`, `angtemp`, `dihtemp`,
//! `density`, `ecoredcore`, `pseudofac` lines), asks for any missing bond
//! lengths, angles, potential parameters, temperature or density on standard
//! input, derives the intramolecular distance constraints and rotational
//! groups, rewrites the `.mol` file and writes a single-molecule `.ato` file.
//!
//! Not translated (stop with "not supported"): `dihedral` lines (DIH
//! groups, which `moloutc` of the port does not write), `unitcell`
//! (non-cubic boxes), file types other than `.gmol`/`.mol`, and
//! column-map rows inside a `.gmol` file (`nspace` undefined in the Fortran).

use super::ToolContext;
use crate::ato::cosd;
use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::{read_int, read_real, IoErr, TextFile};
use crate::fortran::output::OutFile;
use crate::fortran::text::{fstr, indexlastdot, parse, substr};
use crate::fortran::F1;
use crate::state::{State, MCOMP, NEAR};

const MSTARLABEL: usize = 2000;

fn eof() -> EpsrError {
    EpsrError::Input("makemole: unexpected end of standard input".into())
}

/// `subroutine get_rotgrp(ir1,ir2,labrott,nrott)` (`nbondp` = number of direct bonds).
fn get_rotgrp(st: &State, nbondp: &[i32], ir1: i32, ir2: i32) -> Vec<i32> {
    let m = &st.mol;
    let mut labrott: Vec<i32> = Vec::new();
    let mut iselect = ir1;
    let mut ir2test = false;
    let mut nrotnew = 0usize;
    let add = |labrott: &mut Vec<i32>, iaref: i32| {
        if iaref != ir1 && iaref != ir2 && !labrott.contains(&iaref) {
            labrott.push(iaref);
        }
    };
    while nrotnew <= labrott.len() {
        for iref in 1..=nbondp[iselect as usize] {
            let iaref = m.iabond[(iselect, iref)];
            if !ir2test {
                ir2test = iselect != ir1 && iaref == ir2;
            }
            add(&mut labrott, iaref);
        }
        nrotnew += 1;
        if nrotnew <= labrott.len() {
            iselect = labrott[nrotnew - 1];
        }
    }
    if ir2test {
        for iref in 1..=nbondp[ir2 as usize] {
            let iaref = m.iabond[(ir2, iref)];
            add(&mut labrott, iaref);
        }
    }
    labrott
}

/// `subroutine get_whole_rotgrp(im,ir1,ir2,labrott,nrott)`.
fn get_whole_rotgrp(st: &State, im: i32, ir1: i32, ir2: i32) -> Vec<i32> {
    let m = &st.mol;
    (m.naf[im]..=m.nal[im]).filter(|&ia| ia != ir1 && ia != ir2).collect()
}

/// Read one real from standard input (`read(5,*,iostat=ierr) x`): returns
/// `(value, ierr)`; end of input is fatal (the Fortran would loop forever).
fn rd_real_iostat(ctx: &mut ToolContext, old: f32) -> Result<(f32, i32)> {
    let mut r = ctx.stdin.list();
    match r.real() {
        Ok(Some(v)) => Ok((v, 0)),
        Ok(None) => Ok((old, 0)),
        Err(IoErr::Eof) => Err(eof()),
        Err(e) => Ok((old, e.iostat())),
    }
}

/// `text(1:n)` of an internal write record (blank padded).
fn cut(s: &str, n: i32) -> String {
    substr(s, 1, n)
}

/// Run the `makemole` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    if ctx.args.is_empty() {
        return crate::error::unsupported("makemole without a .mol file argument (interactive file search)");
    }
    // readfile=sometext(ncf(2):ncl(nwords)), ilenreadfile=ncl(2)-ncf(2)+1
    let readfile0 = ctx.args.join(" ");
    let mut readfile = readfile0.clone();
    let mut ext = ".mol".to_string();
    let ilenreadfile0 = ctx.args[0].len();
    let ilenuptodot = indexlastdot(&readfile) - 1;
    if ilenuptodot > 0 {
        let u = ilenuptodot as usize;
        ext = substr(&readfile, u as i32 + 1, ilenreadfile0 as i32);
        if ext != ".mol" {
            st.log(&format!("\nmakemole> Unable to read files of type: {ext}"));
            return Ok(());
        }
        readfile = readfile[..u].to_string();
    }
    let fnamemol = format!("{}{}{}", st.sym.directory_name, readfile, ext);
    let mut f = match TextFile::open(std::path::Path::new(&fnamemol)) {
        Ok(f) => f,
        Err(_) => {
            st.log(&format!("\nmakemole> Could not open makemole file : {fnamemol}"));
            return Ok(());
        }
    };
    let sometext = f.read_line().unwrap_or_default();
    let w = parse(&sometext, 200);
    let filetype;
    let iatomlabeladd: i32;
    let mut nspace: i32 = 0;
    if w.nwords() > 0 && w.word(&sometext, 1).find('.') == Some(0) {
        filetype = w.word(&sometext, 1).to_string();
        if w.nwords() > 1 {
            iatomlabeladd = read_int(w.word(&sometext, 2))?;
        } else {
            st.log(&format!("\nmakemole> First line must have two numbers for file type: {filetype}"));
            return Ok(());
        }
    } else {
        // original column-map .mol file
        filetype = ".mol".to_string();
        iatomlabeladd = 0;
        nspace = read_int(w.word(&sometext, 1))?;
    }
    let mut g = Grid::default();
    if filetype == ".mol" {
        g.ruler = f.read_line().unwrap_or_default();
    }
    if filetype != ".gmol" && filetype != ".mol" {
        return crate::error::unsupported(&format!("makemole: file type {filetype}"));
    }
    st.log(&format!("\nmakemole> File type: {filetype} Atom label add: {iatomlabeladd:>2}"));

    // ---- defaults
    let m = &mut st.mol;
    m.ntype = 0;
    let mut lbond1: Vec<String> = Vec::new();
    let mut lbond2: Vec<String> = Vec::new();
    let mut bread: Vec<f32> = Vec::new();
    let mut langle: Vec<[String; 3]> = Vec::new();
    let mut aread: Vec<f32> = Vec::new();
    let mut extra: Vec<(i32, i32, f32)> = Vec::new();
    let mut rots: Vec<(i32, i32)> = Vec::new();
    let mut wrots: Vec<(i32, i32)> = Vec::new();
    let mut irots: Vec<(i32, i32)> = Vec::new();
    let mut nchread: Vec<(String, String)> = Vec::new();
    struct Pot {
        label: String,
        p: f32,
        s: f32,
        mass: f32,
        q: f32,
        c: String,
        nsph: f32,
        rmin: f32,
        rmax: f32,
        cg: f32,
    }
    let mut pots: Vec<Pot> = Vec::new();
    let mut qrs: Vec<(String, f32)> = Vec::new();
    m.simtk = 300.0;
    m.tempb = 65000.0;
    m.tempangreq = 1.0;
    m.tempdihreq = 10.0;
    st.dist.rho = 0.0;
    st.mol.ecore = 0.0;
    st.mol.dcore = 0.0;
    let mut pseudofac: f32 = 1.0;
    let mut alabelsave: Vec<String> = vec![String::new()];
    let mut ia: i32 = 0;
    loop {
        let sometext = match f.read_line() {
            Ok(l) => l,
            Err(_) => break,
        };
        let w = parse(&sometext, 200);
        let nwords = w.nwords();
        if nwords == 0 {
            continue;
        }
        let wd = |i: i32| -> &str { w.word(&sometext, i) };
        let first = wd(1).to_string();
        match first.as_str() {
            "atom" => {
                if nwords > 6 && aread.len() < MSTARLABEL {
                    if ia as usize >= MSTARLABEL {
                        st.log(&format!("\nmakemole> Attempted to create more than {MSTARLABEL:>5} atoms"));
                        return Ok(());
                    }
                    ia += 1;
                    let m = &mut st.mol;
                    m.ensure_atoms(ia as usize + 1);
                    m.ensure_bonds(NEAR);
                    let mut isref = 3;
                    m.alabel[ia] = fstr(wd(isref), 3);
                    alabelsave.push(m.alabel[ia].clone());
                    isref += 1;
                    m.xa[ia] = read_real(wd(isref))?;
                    isref += 1;
                    m.ya[ia] = read_real(wd(isref))?;
                    isref += 1;
                    m.za[ia] = read_real(wd(isref))?;
                    isref += 1;
                    m.nbond[ia] = read_int(wd(isref))?;
                    let mut ib = 0;
                    while ib < m.nbond[ia] && (ib as usize) < NEAR && isref < nwords {
                        ib += 1;
                        isref += 1;
                        m.iabond[(ia, ib)] = read_int(wd(isref))? + iatomlabeladd;
                    }
                    m.nbond[ia] = ib;
                }
            }
            "bond" => {
                if nwords >= 4 && bread.len() < MSTARLABEL {
                    lbond1.push(fstr(wd(2), 3));
                    lbond2.push(fstr(wd(3), 3));
                    bread.push(read_real(wd(4))?);
                }
            }
            "angle" => {
                if nwords >= 5 && aread.len() < MSTARLABEL {
                    langle.push([fstr(wd(2), 3), fstr(wd(3), 3), fstr(wd(4), 3)]);
                    aread.push(read_real(wd(5))?);
                }
            }
            "extra" => {
                if nwords >= 4 && extra.len() < MSTARLABEL {
                    let a = read_int(wd(2))?;
                    let b = read_int(wd(3))?;
                    let e = read_real(wd(4))?;
                    extra.push((a + iatomlabeladd, b + iatomlabeladd, e));
                }
            }
            "rot" | "wholerot" | "introt" => {
                let list = match first.as_str() {
                    "rot" => &mut rots,
                    "wholerot" => &mut wrots,
                    _ => &mut irots,
                };
                if nwords >= 3 && list.len() < MSTARLABEL {
                    let a = read_int(wd(2))?;
                    let b = read_int(wd(3))?;
                    list.push((a + iatomlabeladd, b + iatomlabeladd));
                }
            }
            "dihedral" => {
                if nwords >= 6 {
                    return crate::error::unsupported("makemole: dihedral angles (DIH rotational groups)");
                }
            }
            "changelabel" => {
                if nwords > 2 {
                    nchread.push((fstr(wd(2), 3), fstr(wd(3), 3)));
                }
            }
            "potential" => {
                if nwords >= 7 && pots.len() < MCOMP {
                    let s = read_real(wd(4))?;
                    let (nsph, rmin, rmax) = if nwords > 9 {
                        (read_real(wd(8))?, read_real(wd(9))?, read_real(wd(10))?)
                    } else {
                        (1.0, 0.0, 0.5 * s)
                    };
                    let cg = if nwords > 10 { read_real(wd(11))? } else { 1.0 };
                    pots.push(Pot {
                        label: fstr(wd(2), 3),
                        p: read_real(wd(3))?,
                        s,
                        mass: read_real(wd(5))?,
                        q: read_real(wd(6))?,
                        c: fstr(wd(7), 3),
                        nsph,
                        rmin,
                        rmax,
                        cg,
                    });
                }
            }
            "qradius" => {
                if nwords >= 3 && qrs.len() < MCOMP {
                    qrs.push((fstr(wd(2), 3), read_real(wd(3))?));
                }
            }
            "temperature" | "vibtemp" | "angtemp" | "dihtemp" => {
                if nwords > 1 {
                    if let Ok(v) = read_real(wd(2)) {
                        if v > 0.0 {
                            let m = &mut st.mol;
                            match first.as_str() {
                                "temperature" => m.simtk = v,
                                "vibtemp" => m.tempb = v,
                                "angtemp" => m.tempangreq = v,
                                _ => m.tempdihreq = v,
                            }
                        }
                    }
                }
            }
            "density" => {
                if nwords > 1 {
                    st.dist.rho = read_real(wd(2))?;
                }
            }
            "ecoredcore" => {
                if nwords > 1 {
                    st.mol.ecore = read_real(wd(2))?;
                }
                if nwords > 2 {
                    st.mol.dcore = read_real(wd(3))?;
                }
            }
            "unitcell" => {
                if nwords > 6 {
                    return crate::error::unsupported("makemole: unitcell (non-cubic boxes)");
                }
            }
            "pseudofac" => {
                if nwords < 2 {
                    return crate::error::unsupported("makemole: pseudofac without a value");
                }
                pseudofac = read_real(wd(2))?;
                if pseudofac <= 0.0 {
                    pseudofac = 1.0;
                }
            }
            _ => {
                if nwords > 1 {
                    if filetype != ".mol" {
                        // nspace is undefined for this file type in the Fortran
                        return crate::error::unsupported("makemole: column-map lines in a .gmol file");
                    }
                    if !g.row(&mut st, &sometext, &w, nspace, &mut ia)? {
                        return Ok(());
                    }
                }
            }
        }
    }
    let natom = ia;
    st.mol.natom = natom;
    st.mol.ensure_atoms(natom as usize + 1);
    st.mol.ensure_bonds(NEAR);
    // ---- bond lengths
    for ia in 1..=natom {
        let la1 = st.mol.alabel[ia].clone();
        let mut ib = 0;
        while ib < st.mol.nbond[ia] {
            ib += 1;
            let la2 = st.mol.alabel[st.mol.iabond[(ia, ib)]].clone();
            let mut ifound = 0usize;
            for k in 0..lbond1.len() {
                if (la1 == lbond1[k] && la2 == lbond2[k]) || (la1 == lbond2[k] && la2 == lbond1[k]) {
                    ifound = k + 1;
                    break;
                }
            }
            if ifound == 0 {
                // the Fortran loop `do while(ifound.eq.0.and.ierr.ne.0)` runs once
                st.log(&format!("\nmakemole> Give the bond length for atom pair {la1:<3} {la2:<3}: "));
                lbond1.push(la1.clone());
                lbond2.push(la2.clone());
                bread.push(0.0);
                ifound = bread.len();
                let (v, _e) = rd_real_iostat(ctx, 0.0)?;
                bread[ifound - 1] = v;
            }
            st.mol.tbond[(ia, ib)] = bread[ifound - 1];
        }
    }
    st.mol.naf[1] = 1;
    st.mol.nal[1] = natom;
    let mut nbondp = vec![0i32; natom as usize + 1];
    for ia in 1..=natom {
        nbondp[ia as usize] = st.mol.nbond[ia];
    }
    // ---- bond angles -> extra distances
    for ia2 in 1..=natom {
        let la2 = st.mol.alabel[ia2].clone();
        let mut ib1 = 0;
        while ib1 < nbondp[ia2 as usize] - 1 {
            ib1 += 1;
            let ia1 = st.mol.iabond[(ia2, ib1)];
            let la1 = st.mol.alabel[ia1].clone();
            let mut ib2 = ib1;
            while ib2 < nbondp[ia2 as usize] {
                ib2 += 1;
                let ia3 = st.mol.iabond[(ia2, ib2)];
                let la3 = st.mol.alabel[ia3].clone();
                let m = &st.mol;
                let iexist = (1..=m.nbond[ia1]).any(|is| m.iabond[(ia1, is)] == ia3);
                if !iexist {
                    let bond21 = m.tbond[(ia2, ib1)].abs();
                    let bond23 = m.tbond[(ia2, ib2)].abs();
                    let mut ifound = 0usize;
                    for k in 0..langle.len() {
                        let a = &langle[k];
                        if (la1 == a[0] && la2 == a[1] && la3 == a[2]) || (la1 == a[2] && la2 == a[1] && la3 == a[0]) {
                            ifound = k + 1;
                            break;
                        }
                    }
                    if ifound == 0 {
                        // runs once (see the bond lengths above)
                        st.log(&format!("\nmakemole> Give the angle subtended by the triplet {la1:<3} {la2:<3} {la3:<3}: "));
                        langle.push([la1.clone(), la2.clone(), la3.clone()]);
                        aread.push(0.0);
                        ifound = aread.len();
                        let (v, _e) = rd_real_iostat(ctx, 0.0)?;
                        aread[ifound - 1] = v;
                    }
                    let bangle = aread[ifound - 1];
                    let mut b13 = bond21 * bond21 + bond23 * bond23 - 2.0 * bond21 * bond23 * cosd(bangle);
                    b13 = b13.sqrt();
                    let m = &mut st.mol;
                    if (m.nbond[ia1] as usize) < NEAR {
                        m.nbond[ia1] += 1;
                        let n = m.nbond[ia1];
                        m.iabond[(ia1, n)] = ia3;
                        m.tbond[(ia1, n)] = -b13;
                    }
                    if (m.nbond[ia3] as usize) < NEAR {
                        m.nbond[ia3] += 1;
                        let n = m.nbond[ia3];
                        m.iabond[(ia3, n)] = ia1;
                        m.tbond[(ia3, n)] = -b13;
                    }
                }
            }
        }
    }
    let mut ir: i32 = 0;
    // ---- extra distances
    {
        let m = &mut st.mol;
        let mut ia1 = 0;
        while ia1 < natom - 1 {
            ia1 += 1;
            let mut ia2 = ia1;
            while ia2 < natom {
                ia2 += 1;
                let iexist = (1..=m.nbond[ia1]).any(|is| m.iabond[(ia1, is)] == ia2);
                if !iexist {
                    let found = extra.iter().find(|e| (ia1 == e.0 && ia2 == e.1) || (ia1 == e.1 && ia2 == e.0));
                    if let Some(e) = found {
                        let d = -(e.2.abs());
                        if (m.nbond[ia1] as usize) < NEAR {
                            m.nbond[ia1] += 1;
                            let n = m.nbond[ia1];
                            m.iabond[(ia1, n)] = ia2;
                            m.tbond[(ia1, n)] = d;
                        }
                        if (m.nbond[ia2] as usize) < NEAR {
                            m.nbond[ia2] += 1;
                            let n = m.nbond[ia2];
                            m.iabond[(ia2, n)] = ia1;
                            m.tbond[(ia2, n)] = d;
                        }
                    }
                }
            }
        }
    }
    // ---- rotational groups
    let mut nrottot: i32 = 0;
    let push_group = |st: &mut State, ir: i32, a1: i32, a2: i32, grp: &[i32], nrottot: &mut i32| {
        let m = &mut st.mol;
        m.ensure_rotgrps(ir as usize + 1);
        m.rotgrp[ir] = "ROT".into();
        m.labaxi[(ir, 1)] = a1;
        m.labaxi[(ir, 2)] = a2;
        m.nrot[ir] = grp.len() as i32;
        m.rotf[ir] = *nrottot + 1;
        *nrottot += m.nrot[ir];
        m.rotl[ir] = *nrottot;
        m.ensure_labrot(*nrottot as usize + 1);
        for (i, &a) in grp.iter().enumerate() {
            let labref = m.rotf[ir] + i as i32;
            m.labrot[labref] = a;
        }
        m.ndihedralang[ir] = 0;
    };
    for &(ir1, ir2) in rots.iter() {
        let g1 = get_rotgrp(&st, &nbondp, ir1, ir2);
        let g2 = get_rotgrp(&st, &nbondp, ir2, ir1);
        ir += 1;
        if g1.len() < g2.len() {
            push_group(&mut st, ir, ir1, ir2, &g1, &mut nrottot);
        } else {
            push_group(&mut st, ir, ir2, ir1, &g2, &mut nrottot);
        }
    }
    for &(ir1, ir2) in wrots.iter() {
        let g2 = get_whole_rotgrp(&st, 1, ir1, ir2);
        ir += 1;
        push_group(&mut st, ir, ir1, ir2, &g2, &mut nrottot);
        // The Fortran assigns rotgrp(ir)='ROT' *before* testsize_rotgrps(ir)
        // enlarges the array, so the label is lost and the newly allocated
        // element (zero bytes with gfortran/glibc) is written to the .ato file.
        st.mol.rotgrp[ir] = "\0\0\0".into();
    }
    for &(ir1, ir2) in irots.iter() {
        let mut g1 = get_rotgrp(&st, &nbondp, ir1, ir2);
        let g2 = get_rotgrp(&st, &nbondp, ir2, ir1);
        g1.retain(|a| g2.contains(a));
        ir += 1;
        push_group(&mut st, ir, ir1, ir2, &g1, &mut nrottot);
    }
    st.mol.ngrot[1] = ir;
    st.mol.ngrrot = ir;
    st.mol.nrottot = nrottot;
    // ---- temperature and density
    if st.mol.simtk <= 0.0 {
        let mut ierr = 1;
        while ierr != 0 {
            let (v, e) = rd_real_iostat(ctx, st.mol.simtk)?;
            st.mol.simtk = v;
            ierr = e;
        }
    }
    if st.dist.rho <= 0.0 {
        let mut ierr = 1;
        while ierr != 0 {
            let (v, e) = rd_real_iostat(ctx, st.dist.rho)?;
            st.dist.rho = v;
            ierr = e;
        }
    }
    // ---- the molecule
    {
        let m = &mut st.mol;
        m.nmol = 1;
        m.xadd[1] = 0.0;
        m.yadd[1] = 0.0;
        m.zadd[1] = 0.0;
        m.phix[1] = 0.0;
        m.phiy[1] = 0.0;
        m.phiz[1] = 0.0;
        m.tethered[1] = false;
        m.nato[1] = natom;
        m.stepmi = 0.1;
        m.stepri = if m.ngrot[1] > 0 { 1.0 } else { 0.0 };
        m.stepr = 0.3;
        m.stepm = 1.0;
        m.ntype = 0;
        m.natomr = 0;
    }
    for ia in 1..=natom {
        if let Some((_, new)) = nchread.iter().find(|(old, _)| st.mol.alabel[ia] == *old) {
            st.mol.alabel[ia] = new.clone();
        }
        st.ftype(ia);
        if st.mol.qtype[ia] != 1 {
            st.mol.natomr += 1;
        }
    }
    let nt = st.mol.ntype as usize;
    {
        let m = &mut st.mol;
        for a in [&mut m.natomsphere, &mut m.radminsphere, &mut m.radmaxsphere, &mut m.pseudofacmol, &mut m.rhosphere] {
            if a.len() < nt {
                a.grow(nt);
            }
        }
    }
    for it in 1..=nt {
        let mut ifound = pots.iter().position(|p| p.label == st.mol.type_[it]).map(|p| p + 1).unwrap_or(0);
        if ifound == 0 {
            // read(5,*,iostat=ierr) pread,sread,mread,qread,cread
            let mut vals = (0.0f32, 0.0f32, 0.0f32, 0.0f32, String::new());
            loop {
                st.log(&format!("\nmakemole> Type epsilon, sigma, mass, charge and chemical symbol for component {:<3}\n\nmakemole> ? ", st.mol.type_[it]));
                let mut r = ctx.stdin.list();
                let res: std::result::Result<(), IoErr> = (|| {
                    r.get_r(&mut vals.0)?;
                    r.get_r(&mut vals.1)?;
                    r.get_r(&mut vals.2)?;
                    r.get_r(&mut vals.3)?;
                    if let Some(w) = r.word()? {
                        vals.4 = fstr(&w, 3);
                    }
                    Ok(())
                })();
                match res {
                    Ok(()) => break,
                    Err(IoErr::Eof) => return Err(eof()),
                    Err(_) => continue,
                }
            }
            pots.push(Pot {
                label: st.mol.type_[it].clone(),
                p: vals.0,
                s: vals.1,
                mass: vals.2,
                q: vals.3,
                c: vals.4.clone(),
                nsph: 1.0,
                rmin: 0.0,
                rmax: 0.5 * vals.1,
                cg: 1.0,
            });
            ifound = pots.len();
        }
        let p = &pots[ifound - 1];
        let m = &mut st.mol;
        m.atsymbol[it] = p.c.clone();
        m.iexchange[it] = 0;
        m.epscor[it] = p.p;
        m.sigcor[it] = p.s;
        m.amassa[it] = p.mass;
        m.charge[it] = p.q;
        m.natomsphere[it] = p.nsph;
        m.radminsphere[it] = p.rmin;
        m.radmaxsphere[it] = p.rmax;
        m.pseudofacmol[it] = p.cg;
        m.qradius[it] = qrs.iter().find(|q| q.0 == m.type_[it]).map(|q| q.1).unwrap_or(0.0);
    }
    // ---- box size (cubic)
    {
        let vol = st.mol.natomr as f32 / st.dist.rho;
        let dmax = vol.powf(1.0f32 / 3.0f32);
        let l = &mut st.lat;
        st.mol.dmax = dmax;
        l.alat = dmax;
        l.blat = dmax;
        l.clat = dmax;
        l.phib = 90.0;
        l.thetac = 0.0;
        l.phic = 0.0;
        l.cubic = true;
        st.mol.dmaxmol = dmax;
    }
    // ---- rewrite the .mol file
    let mut out = OutFile::create(std::path::Path::new(&fnamemol));
    if filetype == ".mol" {
        let fnameatm = format!("{}{}.atm", st.sym.directory_name, readfile);
        let mut atm = OutFile::create(std::path::Path::new(&fnameatm));
        g.write(&mut out, &mut atm, nspace);
        atm.close()?;
    } else {
        let mut w = ListWriter::new();
        w.chars(&filetype).int(iatomlabeladd);
        out.line(&w.finish());
    }
    for ia1 in 1..=natom {
        if filetype != ".gmol" {
            break;
        }
        let m = &st.mol;
        let mut w = ListWriter::new();
        w.chars("atom ").int(ia1 - iatomlabeladd).chars(" ").chars(&format!("{:<3}", alabelsave[ia1 as usize])).chars(" ");
        w.real(m.xa[ia1]).real(m.ya[ia1]).real(m.za[ia1]).int(nbondp[ia1 as usize]);
        for j in 1..=nbondp[ia1 as usize] {
            w.int(m.iabond[(ia1, j)] - iatomlabeladd);
        }
        out.line(&w.finish());
    }
    for k in 0..lbond1.len() {
        let s = fmt("(a,2(1x,a3),1x,f10.5)", &[FV::S("bond"), FV::Ch(&lbond1[k], 3), FV::Ch(&lbond2[k], 3), FV::R(bread[k])]);
        out.line(&cut(&s, 4 + 2 * 4 + 1 + 10));
    }
    for k in 0..langle.len() {
        let a = &langle[k];
        let s = fmt("(a,3(1x,a3),1x,f10.5)", &[FV::S("angle"), FV::Ch(&a[0], 3), FV::Ch(&a[1], 3), FV::Ch(&a[2], 3), FV::R(aread[k])]);
        out.line(&cut(&s, 5 + 3 * 4 + 1 + 10));
    }
    for (name, list) in [("rot", &rots), ("wholerot", &wrots), ("introt", &irots)] {
        for &(a, b) in list.iter() {
            let s = fmt("(a,2(1x,i3))", &[FV::S(name), FV::I((a - iatomlabeladd) as i64), FV::I((b - iatomlabeladd) as i64)]);
            out.line(&cut(&s, name.len() as i32 + 2 * 4));
        }
    }
    for &(a, b, e) in extra.iter() {
        let s = fmt("(a,2(1x,i3),1x,f10.5)", &[FV::S("extra"), FV::I((a - iatomlabeladd) as i64), FV::I((b - iatomlabeladd) as i64), FV::R(e)]);
        out.line(&cut(&s, 5 + 2 * 4 + 1 + 10));
    }
    for (o, n) in nchread.iter() {
        let s = fmt("(a,2(1x,a3))", &[FV::S("changelabel"), FV::Ch(o, 3), FV::Ch(n, 3)]);
        out.line(&cut(&s, 11 + 8));
    }
    out.line(fmt("(a,1x,f8.2)", &[FV::S("pseudofac"), FV::R(pseudofac)]).trim_end());
    let f113 = "(a,1x,a3,4(1x,f12.6),1x,a3,1x,e13.6,2(1x,f12.6),1x,e13.6)";
    for p in pots.iter() {
        let s = fmt(
            f113,
            &[
                FV::S("potential"),
                FV::Ch(&p.label, 3),
                FV::R(p.p),
                FV::R(p.s),
                FV::R(p.mass),
                FV::R(p.q),
                FV::Ch(&p.c, 3),
                FV::R(p.nsph),
                FV::R(p.rmin),
                FV::R(p.rmax),
                FV::R(p.cg),
            ],
        );
        out.line(s.trim_end());
    }
    for (l, q) in qrs.iter() {
        let s = fmt(f113, &[FV::S("qradius"), FV::Ch(l, 3), FV::R(*q)]);
        out.line(&cut(&s, 8 + 4 + 2 * 11));
    }
    let f114 = "(a,1x,e13.6)";
    out.line(&cut(&fmt(f114, &[FV::S("temperature"), FV::R(st.mol.simtk)]), 11 + 14));
    out.line(&cut(&fmt(f114, &[FV::S("vibtemp"), FV::R(st.mol.tempb)]), 7 + 14));
    out.line(&cut(&fmt(f114, &[FV::S("angtemp"), FV::R(st.mol.tempangreq)]), 7 + 14));
    out.line(&cut(&fmt(f114, &[FV::S("dihtemp"), FV::R(st.mol.tempdihreq)]), 7 + 14));
    out.line(&cut(&fmt("(a,1x,f10.5,1x,f10.5)", &[FV::S("ecoredcore"), FV::R(st.mol.ecore), FV::R(st.mol.dcore)]), 10 + 22));
    out.line(&cut(&fmt(f114, &[FV::S("density"), FV::R(st.dist.rho)]), 7 + 14));
    out.close()?;
    // ---- write the .ato file
    let fnameato = format!("{}{}.ato", st.sym.directory_name, readfile);
    {
        let m = &mut st.mol;
        m.nmoltype = 1;
        m.moltemplate[1] = readfile.clone();
        m.interntrans[1] = 0.0;
        m.internrot[1] = 0.0;
        m.wholerot[1] = 0.0;
        m.wholetrans[1] = 0.0;
        // moloutc allocates scalefacmol (= 1.0) when it is not allocated
        m.scalefacmol = F1::new(m.nmol as usize);
        m.scalefacmol.fill(1.0);
    }
    st.moloutc(&fnameato)?;
    Ok(())
}

/// A label (atom position) of the column-map `.mol` format.
#[derive(Clone, Default)]
struct GLabel {
    ir: i32,
    ic: i32,
    n: i32,
    label: String,
    offset: String,
    ianumber: [i32; 10],
}

/// State of the column-map `.mol` reader (`nrow`, `ncolumnmax`, stars and labels).
#[derive(Default)]
struct Grid {
    ruler: String,
    nrow: i32,
    ncolumnmax: i32,
    stars: Vec<(i32, i32)>,
    labels: Vec<GLabel>,
    nread: i32,
}

/// `read(c,*,iostat=ierr) n` for a single character `c` (gfortran behaviour):
/// digits are read, `/` and `,` leave `n` unchanged without error, anything
/// else is an error.
fn read_char_int(c: u8, n: &mut i32) -> bool {
    match c {
        b'0'..=b'9' => {
            *n = (c - b'0') as i32;
            true
        }
        b'/' | b',' => true,
        _ => false,
    }
}

impl Grid {
    /// One row of the column map (the `else if(nwords.gt.1)` branch of the reader).
    /// Returns `Ok(false)` when the Fortran gives up (`return` without output).
    fn row(&mut self, st: &mut State, sometext: &str, w: &crate::fortran::Words, nspace: i32, ia: &mut i32) -> Result<bool> {
        self.nrow += 1;
        let nrow = self.nrow;
        let nwords = w.nwords();
        let mut iw = 1;
        while iw < nwords {
            iw += 1;
            let cf = w.ncf[(iw - 1) as usize];
            let cl = w.ncl[(iw - 1) as usize];
            let c0 = sometext.as_bytes()[(cf - 1) as usize];
            let a1 = if read_char_int(c0, &mut self.nread) {
                fstr(&substr(sometext, cf + 1, cl), 20)
            } else {
                self.nread = 1;
                fstr(&substr(sometext, cf, cl), 20)
            };
            let ncolumn = (cf - 4 + nspace - 1) / nspace;
            if ncolumn > self.ncolumnmax {
                self.ncolumnmax = ncolumn;
            }
            if a1 == "*" {
                self.stars.push((nrow, ncolumn));
                continue;
            }
            let mut lab = GLabel { ir: nrow, ic: ncolumn, n: self.nread, ..Default::default() };
            let (irl1, icl1) = (nrow, ncolumn);
            if substr(&a1, 1, 3) == "map" {
                lab.label = "map".into();
                let nchara1 = crate::fortran::text::findex(&crate::fortran::format::pad_to(&a1, 20), " ") - 1;
                lab.offset = fstr(&substr(&a1, 4, nchara1), 20);
                let ncharoffset = nchara1 - 3;
                let off = crate::fortran::format::pad_to(&lab.offset, 20);
                let fx = |p: &str, sub: &str| crate::fortran::text::findex(p, sub);
                let ilenopenbra = fx(&off, "(");
                let ilencomma = fx(&off, ",");
                let ilenclosebra = fx(&off, ")");
                let bad = |st: &State| {
                    st.log(&format!(
                        "\nmakemole> Unrecognised format for mapped atom {} at row {irl1:>3} and column {icl1:>3}\n\n - check the .mol file",
                        substr(&a1, 1, nchara1)
                    ));
                };
                let (irmap, icmap);
                if ilenopenbra + ilencomma + ilenclosebra > 8 {
                    let r1 = read_int(&substr(&off, ilenopenbra + 1, ilencomma - 1));
                    let r2 = read_int(&substr(&off, ilencomma + 1, ilenclosebra - 1));
                    match (r1, r2) {
                        (r1, Ok(c)) if r1.as_ref().map(|v| *v > 0).unwrap_or(false) && c > 0 => {
                            irmap = r1.unwrap();
                            icmap = c;
                        }
                        _ => {
                            bad(st);
                            return Ok(false);
                        }
                    }
                } else {
                    let sub = substr(&off, 2, ncharoffset);
                    let ilenread = fx(&sub, "+") + fx(&sub, "-");
                    if ilenread > 1 {
                        let ilenread = ilenread + 1;
                        let r1 = read_int(&substr(&off, 1, ilenread - 1));
                        let r2 = read_int(&substr(&off, ilenread, ncharoffset));
                        match (r1, r2) {
                            (Ok(a), Ok(b)) => {
                                irmap = irl1 + a;
                                icmap = icl1 + b;
                            }
                            _ => {
                                // only the iostat of the second read is tested
                                return crate::error::unsupported("makemole: invalid row offset of a mapped atom");
                            }
                        }
                    } else {
                        bad(st);
                        return Ok(false);
                    }
                }
                let found = self.labels.iter().position(|l| l.ir == irmap && l.ic == icmap);
                match found {
                    None => {
                        st.log(&format!(
                            "\nmakemole> Label referred to by mapped atom{} not found at row {irmap:>3} and column {icmap:>3}\n\n - check the .mol file",
                            substr(&a1, 1, nchara1)
                        ));
                        return Ok(false);
                    }
                    Some(k) => {
                        lab.n = self.labels[k].n;
                        lab.ianumber = self.labels[k].ianumber;
                    }
                }
            } else {
                lab.label = fstr(&a1, 3);
                if lab.n > 9 || lab.n < 0 {
                    return crate::error::unsupported("makemole: label repeat count");
                }
                for ia1 in 1..=lab.n as usize {
                    *ia += 1;
                    let m = &mut st.mol;
                    m.ensure_atoms(*ia as usize + 1);
                    m.ensure_bonds(NEAR);
                    lab.ianumber[ia1] = *ia;
                    m.alabel[*ia] = lab.label.clone();
                    m.xa[*ia] = 0.0;
                    m.ya[*ia] = 0.0;
                    m.za[*ia] = 0.0;
                    m.nbond[*ia] = 0;
                }
            }
            self.labels.push(lab);
            let il1 = self.labels.len() - 1;
            // bonds through the stars next to this label
            for is in 0..self.stars.len() {
                let (irs, ics) = self.stars[is];
                let irdif1 = self.labels[il1].ir - irs;
                let icdif1 = self.labels[il1].ic - ics;
                if (irdif1 == 1 && icdif1.abs() < 2) || (irdif1 == 0 && icdif1 == 1) {
                    for il2 in 0..il1 {
                        let irdif2 = (irs - self.labels[il2].ir).abs();
                        let icdif2 = (ics - self.labels[il2].ic).abs();
                        if irdif2 < 2 && icdif2 < 2 {
                            for ia1 in 1..=self.labels[il1].n as usize {
                                let ia1ref = self.labels[il1].ianumber[ia1];
                                for ia2 in 1..=self.labels[il2].n as usize {
                                    let ia2ref = self.labels[il2].ianumber[ia2];
                                    let m = &mut st.mol;
                                    let exists = (1..=m.nbond[ia1ref]).any(|ib| m.iabond[(ia1ref, ib)] == ia2ref);
                                    if !exists {
                                        m.nbond[ia1ref] += 1;
                                        m.nbond[ia2ref] += 1;
                                        let (n1, n2) = (m.nbond[ia1ref], m.nbond[ia2ref]);
                                        m.iabond[(ia1ref, n1)] = ia2ref;
                                        m.iabond[(ia2ref, n2)] = ia1ref;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        Ok(true)
    }

    /// Rewrite of the column map into the `.mol` file and of the atom numbers into the `.atm` file.
    fn write(&self, out: &mut OutFile, atm: &mut OutFile, nspace: i32) {
        let mut w = ListWriter::new();
        w.int(nspace);
        let l = w.finish();
        out.line(&l);
        atm.line(&l);
        let mut ncharmax = (self.ncolumnmax + 1) * nspace;
        if ncharmax < 80 {
            ncharmax = 80;
        }
        let ruler = fstr(&self.ruler, 256);
        out.line(&substr(&ruler, 1, ncharmax));
        atm.line(&substr(&ruler, 1, ncharmax));
        for ir in 1..=self.nrow {
            let mut buf: Vec<u8> = vec![b' '; 3020 + 8];
            let put = |buf: &mut Vec<u8>, at: i32, s: &str| {
                for (k, c) in s.bytes().enumerate() {
                    let p = (at - 1) as usize + k;
                    if p < buf.len() {
                        buf[p] = c;
                    }
                }
            };
            put(&mut buf, 1, &fmt("(i3,1x)", &[FV::I(ir as i64)]));
            let mut ncharmax = 0;
            for &(irs, ics) in self.stars.iter() {
                if irs == ir {
                    let icfirst = (ics - 1) * nspace + 2 + 4;
                    put(&mut buf, icfirst, "*");
                    if icfirst > ncharmax {
                        ncharmax = icfirst;
                    }
                }
            }
            for l in self.labels.iter() {
                if l.ir == ir {
                    let icfirst = (l.ic - 1) * nspace + 2 + 4;
                    if l.n > 1 {
                        put(&mut buf, icfirst - 1, &fmt("(i1)", &[FV::I(l.n as i64)]));
                    }
                    put(&mut buf, icfirst, &crate::fortran::format::pad_to(&l.label, 3));
                    if l.label != "map" {
                        if icfirst + 2 > ncharmax {
                            ncharmax = icfirst + 2;
                        }
                    } else {
                        let off = crate::fortran::format::pad_to(&l.offset, 20);
                        let ilenoffset = crate::fortran::text::findex(&off, " ") - 1;
                        put(&mut buf, icfirst + 3, &substr(&off, 1, ilenoffset));
                        if icfirst + 2 + ilenoffset > ncharmax {
                            ncharmax = icfirst + 2 + ilenoffset;
                        }
                    }
                }
            }
            out.line(&String::from_utf8_lossy(&buf[..ncharmax.max(0) as usize]));
            for l in self.labels.iter() {
                if l.ir == ir {
                    let icfirst = (l.ic - 1) * nspace + 2 + 4;
                    let a = l.ianumber[1];
                    if a < 10 {
                        put(&mut buf, icfirst, &format!("{:<3}", fmt("(i1)", &[FV::I(a as i64)])));
                    } else if a < 100 {
                        put(&mut buf, icfirst, &format!("{:<3}", fmt("(i2)", &[FV::I(a as i64)])));
                    } else {
                        put(&mut buf, icfirst - 1, &fmt("(i3)", &[FV::I(a as i64)]));
                    }
                    if icfirst + 2 > ncharmax {
                        ncharmax = icfirst + 2;
                    }
                }
            }
            atm.line(&String::from_utf8_lossy(&buf[..ncharmax.max(0) as usize]));
        }
    }
}
