//! `fmole` — translation of `fortran_src/fmole.f` (`fmole`, `initpt_fmole`,
//! `update_ato`).
//!
//! STATUS: translated; verified against the Fortran (EPSR-test-ex33, ex34).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> fmole [<file>.ato [<ntimes>]]`
//! Missing arguments are read from standard input.
//!
//! fmole "shakes" every molecule `ntimes` times (intramolecular Monte Carlo
//! moves of all atoms plus, for a single molecule, rotations of its
//! rotational groups) using only the intramolecular energy (harmonic bond and
//! angle restraints + the repulsive WCA core potential built from
//! `ecore`/`dcore`), then rewrites the `.ato` file.  Before that the bond
//! lists, rotational groups and atom parameters are refreshed from each
//! molecule type's template `.ato` file (`update_ato`).
//!
//! NOT PORTED: the Fortran first calls `makemole` for every template that has
//! a `.mol` file in the folder (to regenerate the template `.ato`); `makemole`
//! is a different program of the suite.  If such a `.mol` file exists this
//! port stops with an "unsupported" error (when no `.mol` file exists the
//! Fortran `makemole` does nothing, and that case is reproduced).

use super::makeato::setup_common::{command_text, eof_err, get_full_filename, search_for_files, strip_first_dot, SysCmds};
use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::listread::{columns, read_int, read_real, ListRead, TextFile};
use crate::fortran::text::{indexlastdot, parse};
use crate::neighbour::minimagec;
use crate::potential::{truncate, truncatec};
use crate::state::State;

/// Run the `fmole` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    let sys = SysCmds::load();
    let (mut sometext, mut w) = command_text("fmole", ctx);
    let readfile: String;
    if w.nwords() < 2 {
        let (t, _iaccept) = search_for_files(&st, ctx, &sys, ".ato")?;
        sometext = t;
        w = parse(&sometext, 5);
        if w.nwords() <= 0 {
            st.log("\nfcluster> No .ato file specified! - routine will exit\n");
            return Ok(());
        }
        readfile = w.word(&sometext, 1).to_string();
    } else {
        readfile = w.word(&sometext, 2).to_string();
    }
    // number of shakes
    let mut ntimes: i32 = 0;
    let mut ierr = 1;
    if w.nwords() > 2 {
        if let Ok(v) = read_int(w.word(&sometext, 3)) {
            ntimes = v;
            ierr = 0;
        }
    }
    while ierr != 0 {
        st.log("\nfmole> Type the number of times to perform the shake: ");
        let mut r = ctx.stdin.list();
        match r.int() {
            Ok(v) => {
                if let Some(v) = v {
                    ntimes = v;
                }
                ierr = 0;
            }
            Err(crate::fortran::IoErr::Eof) => return Err(eof_err(crate::fortran::IoErr::Eof)),
            Err(_) => ierr = 1,
        }
    }
    let mut nrenew = 0;
    if nrenew < 1 || nrenew > ntimes {
        nrenew = ntimes;
    }
    let readfile = strip_first_dot(&readfile);
    let fnameato = format!("{}.ato", get_full_filename(&sys, &ctx.dir, &readfile));
    st.log(&format!(" {}\n {}", readfile, fnameato));
    st.mol.ntype = 0;
    st.mol.nmoltype = 0;
    st.mol.sizefactor = 1.0;
    match st.minitc(&fnameato) {
        Ok(()) => {}
        Err(EpsrError::Io(_)) => {
            st.log(&format!("\nfmole> Problems with .ato file {}\nRoutine will exit now\n", fnameato));
            return Ok(());
        }
        Err(e) => return Err(e),
    }
    if st.dist.rho == 0.0 {
        st.log(&format!("\nfmole> Problems with .ato file {}\nRoutine will exit now\n", fnameato));
        return Ok(());
    }
    // makemole for each template (only the "no .mol file" outcome is ported)
    let nmoltypesave = st.mol.nmoltype;
    let templates: Vec<String> = (1..=nmoltypesave).map(|it| st.mol.moltemplate[it].clone()).collect();
    for t in &templates {
        let ilen = crate::fortran::text::len_to_blank(t);
        if ilen > 0 {
            makemole_template(&st, ctx, &t[..ilen])?;
        }
    }
    st.mol.sizefactor = 1.0;
    st.mol.ntype = 0;
    st.mol.nmoltype = 0;
    st.minitc(&fnameato)?;
    update_ato(&mut st, ctx, &sys)?;
    st.normalise_coordinates();
    let rminpt: f32 = 4.0;
    let rmaxpt: f32 = 6.0;
    initpt_fmole(&mut st, rminpt, rmaxpt);
    if st.mol.stepmi <= 0.0 && st.mol.stepri <= 0.0 {
        st.log("\nfmole> Step size is set to zero!\n");
        return Ok(());
    }
    let stepms = st.mol.stepm;
    let stepris = st.mol.stepri;
    let steprs = st.mol.stepr;
    let tempdihsave = st.mol.tempdih;
    st.mol.stepm = 0.0;
    st.mol.stepr = 0.0;
    if st.mol.nmol > 1 {
        st.mol.stepri = 0.0;
    }
    let mut ftest: f32 = 0.0;
    let mut ftestb: f32 = 0.0;
    let intra = 1;
    let rotfreq = 1;
    st.mol.refine_intra = false;
    st.mol.coul_intra = false;
    st.mol.refine_intra_fac = 0.0;
    st.mol.coul_intra_fac = 0.0;
    st.dist.rejrate = 0.8;
    st.inp.nisum = intra;
    let ireset = 0;
    let mut inn = nrenew;
    let mut it = 0;
    let sepmax: f32 = 0.0;
    let nmolcellin = 10;
    st.initialise_cellneighbourlist(nmolcellin, sepmax);
    while it < ntimes {
        it += 1;
        if inn < nrenew {
            inn += 1;
        } else {
            inn = 1;
            st.adjust_intra_wts();
            let (a, b) = st.fmolec1();
            ftest = a;
            ftestb = b;
        }
        if it == 1 {
            st.log(&format!("\nfmole> Starting intramolecular energy: {}\n", ftestb));
        }
        for im in 1..=st.mol.nmol {
            let m = &mut st.mol;
            let (mut x, mut y, mut z) = (m.xadd[im], m.yadd[im], m.zadd[im]);
            minimagec(&mut x, &mut y, &mut z);
            m.xadd[im] = x;
            m.yadd[im] = y;
            m.zadd[im] = z;
        }
        st.molfit(&mut ftest, &mut ftestb, intra, rotfreq, ireset);
        if 100 * (it / 100) == it || it == ntimes {
            let m = &mut st.mol;
            for j in 1..=2usize {
                m.nmoveit[j - 1] = 0;
                for imoltype in 1..=m.nmoltype {
                    m.nmoveit[j - 1] = m.nmoveit[j - 1] + m.nmovei[(j, imoltype)] + m.ninter[(j, imoltype)];
                }
            }
            st.log(&format!(
                "\nfmole> Iteration {} Intramolecular energy: {}\n No. of moves tried: {} No. of moves rejected: {}",
                it, ftestb, st.mol.nmoveit[0], st.mol.nmoveit[1]
            ));
            st.adjust_intra_wts();
        }
    }
    st.mol.stepm = stepms;
    st.mol.stepri = stepris;
    st.mol.stepr = steprs;
    st.mol.tempdih = tempdihsave;
    st.unnormalise_coordinates();
    st.moloutc(&fnameato)?;
    st.log(" Done fmole");
    Ok(())
}

/// The part of `makemole()` executed when fmole calls it with
/// `sometext = 'makemole '//template`: decide the file name and try to open
/// the `.mol` file.  If it cannot be opened makemole returns without doing
/// anything; otherwise this port stops (makemole is not part of this module).
fn makemole_template(st: &State, ctx: &ToolContext, template: &str) -> Result<()> {
    let mut readfile = template.to_string();
    let mut ext = ".mol".to_string();
    let ilenuptodot = indexlastdot(&readfile) - 1;
    if ilenuptodot > 0 {
        ext = readfile[ilenuptodot as usize..].to_string();
        if ext != ".mol" {
            st.log(&format!("\nmakemole> Unable to read files of type: {}", ext));
            return Ok(());
        }
        readfile = readfile[..ilenuptodot as usize].to_string();
    }
    let fname = format!("{}{}{}", ctx.dir, readfile, ext);
    if std::fs::File::open(&fname).is_err() {
        st.log(&format!("\nmakemole> Could not open makemole file : {}", fname));
        return Ok(());
    }
    crate::error::unsupported(&format!(
        "fmole: template {fname} exists, so the Fortran runs makemole on it first - makemole is not part of the fmole port"
    ))
}

/// `subroutine initpt_fmole(rminpt, rmaxpt)`: r grid, truncation windows,
/// zero potentials, and the repulsive core potential (`generate_mp`).
fn initpt_fmole(st: &mut State, rminpt: f32, rmaxpt: f32) {
    let npairs = st.mol.npairs;
    st.pot.alloc(npairs.max(1) as usize, 1);
    let d = &mut st.dist;
    let dmax = st.mol.dmax;
    d.npt = 1000;
    d.cellst = 0.5 * dmax / d.npt as f32;
    d.npt2 = d.npt - 1;
    d.nptmin = 1;
    d.nptmax = d.npt2;
    for j in 1..=d.npt {
        d.x[j] = (j - 1) as f32 * d.cellst;
        if d.x[j] <= rminpt {
            d.nptmin = j;
        }
        if d.x[j] <= rmaxpt {
            d.nptmax = j;
        }
        d.xstar[j] = d.x[j] / dmax;
    }
    let p = &mut st.pot;
    for j in 1..=d.npt2 {
        let xx = d.x[j];
        p.wfac[j] = truncate(xx, rminpt, rmaxpt);
        p.wfacc[j] = truncatec(xx, rmaxpt);
    }
    p.wfac[d.npt] = 0.0;
    for ip in 1..=npairs {
        for j in 1..=d.npt {
            p.xav[(j, ip)] = 0.0;
            p.xavs[(j, ip)] = 0.0;
            p.xav1[(j, ip)] = 0.0;
            p.xavs1[(j, ip)] = 0.0;
            p.epot[(j, ip)] = 0.0;
            p.epot1[(j, ip)] = 0.0;
            p.ecoul[(j, ip)] = 0.0;
            p.ecoul1[(j, ip)] = 0.0;
        }
        p.irefpcof[ip] = ip;
        p.rminex[ip] = 0.0;
        p.pcoefex[ip] = 0.0;
        p.rmaxex[ip] = 0.0;
        p.pcoefexmax[ip] = 0.0;
    }
    st.mol.power = 12.0;
    st.generate_mp();
}

/// A molecule template read by `update_ato`.
#[derive(Default)]
struct Template {
    natom: i32,
    alabel: Vec<String>,
    iabond: Vec<Vec<i32>>,
    tbond: Vec<Vec<f32>>,
    rotgrp: Vec<String>,
    labaxi: Vec<[i32; 2]>,
    labrot: Vec<Vec<i32>>,
}

fn rd_err(e: crate::fortran::IoErr) -> EpsrError {
    EpsrError::from(e)
}

/// Read a template `.ato` file as `update_ato` does (bonds, rotational groups,
/// and the potential parameters of the atom types found in it).
fn read_template(st: &mut State, f: &mut TextFile) -> Result<Template> {
    let mut t = Template::default();
    let first = f.read_line().unwrap_or_default();
    let w = parse(&first, 5);
    if w.nwords() != 3 {
        let _ = f.skip_record();
        let _ = f.skip_record();
        let _ = f.skip_record();
    } else {
        let _ = f.skip_record();
    }
    let mut typetemplate: Vec<String> = Vec::new();
    {
        let mut r = f.list();
        r.get_i(&mut t.natom).map_err(rd_err)?;
    }
    for _ia in 1..=t.natom {
        let l = f.read_line().unwrap_or_default();
        let lab = columns(&l, 1, 3);
        if !typetemplate.contains(&lab) {
            typetemplate.push(lab.clone());
        }
        t.alabel.push(lab);
        f.skip_record().map_err(rd_err)?;
        let mut r = f.list();
        let mut nb = 0;
        r.get_i(&mut nb).map_err(rd_err)?;
        let mut ib = Vec::new();
        let mut tb = Vec::new();
        for _ja in 1..=nb {
            let mut a = 0;
            let mut b = 0f32;
            r.get_i(&mut a).map_err(rd_err)?;
            r.get_r(&mut b).map_err(rd_err)?;
            ib.push(a);
            tb.push(b);
        }
        t.iabond.push(ib);
        t.tbond.push(tb);
    }
    let mut nrotgrp = 0;
    {
        let mut r = f.list();
        r.get_i(&mut nrotgrp).map_err(rd_err)?;
    }
    for _ir in 1..=nrotgrp {
        let l = f.read_line().unwrap_or_default();
        let g = columns(&l, 1, 3);
        let mut ax = [0i32; 2];
        {
            let mut r = f.list();
            r.get_i(&mut ax[0]).map_err(rd_err)?;
            r.get_i(&mut ax[1]).map_err(rd_err)?;
        }
        if g == "DIH" {
            return crate::error::unsupported("update_ato: DIH rotational groups in a template file");
        }
        let s = f.read_line().unwrap_or_default();
        let w = parse(&s, 5);
        let nrot = read_int(w.word(&s, 1))?;
        let mut v: Vec<i32> = Vec::new();
        {
            let span = if w.nwords() >= 2 {
                w.span(&s, 2, w.nwords())
            } else {
                String::new()
            };
            let mut r = ListRead::from_str(&span);
            for _jr in 1..=(w.nwords() - 1) {
                let mut x = 0;
                r.get_i(&mut x).map_err(rd_err)?;
                v.push(x);
            }
        }
        if nrot > w.nwords() - 1 {
            let mut r = f.list();
            for _jr in w.nwords()..=nrot {
                let mut x = 0;
                r.get_i(&mut x).map_err(rd_err)?;
                v.push(x);
            }
        }
        v.truncate(nrot.max(0) as usize);
        t.rotgrp.push(g);
        t.labaxi.push(ax);
        t.labrot.push(v);
    }
    // potential parameters of the template's atom types
    for _it in 0..typetemplate.len() {
        let l = match f.read_line() {
            Ok(l) => l,
            Err(_) => String::new(),
        };
        let labeltest = columns(&l, 1, 3);
        let symboltest = columns(&l, 5, 3);
        let m = &mut st.mol;
        let mut ic = 0;
        let mut is = 0;
        while ic < m.ntype && is == 0 {
            ic += 1;
            if labeltest == m.type_[ic] {
                is = ic;
            }
        }
        if is > 0 {
            m.atsymbol[is] = symboltest;
            let s = f.read_line().unwrap_or_default();
            let w = parse(&s, 11);
            {
                let span = w.span(&s, 1, 5.min(w.nwords().max(1)));
                let mut r = ListRead::from_str(&span);
                let _ = (|| -> crate::fortran::IoRes<()> {
                    r.get_r(&mut m.epscor[is])?;
                    r.get_r(&mut m.sigcor[is])?;
                    r.get_r(&mut m.amassa[is])?;
                    r.get_r(&mut m.charge[is])?;
                    r.get_r(&mut m.qradius[is])?;
                    Ok(())
                })();
            }
            if w.nwords() > 5 {
                if let Ok(v) = read_real(w.word(&s, 6)) {
                    m.natomsphere[is] = v;
                }
            } else {
                m.natomsphere[is] = 1.0;
            }
            if w.nwords() > 7 {
                let span = w.span(&s, 7, 8);
                let mut r = ListRead::from_str(&span);
                let _ = r.get_r(&mut m.radminsphere[is]).and_then(|_| r.get_r(&mut m.radmaxsphere[is]));
                if m.radmaxsphere[is] <= 0.0 {
                    m.radmaxsphere[is] = 0.5 * m.sigcor[is];
                }
            } else {
                m.radminsphere[is] = 0.0;
                m.radmaxsphere[is] = 0.5 * m.sigcor[is];
            }
            if w.nwords() > 8 {
                if let Ok(v) = read_real(w.word(&s, 9)) {
                    m.pseudofacmol[is] = v;
                }
            } else {
                m.pseudofacmol[is] = 1.0;
            }
        } else {
            let _ = f.skip_record();
        }
    }
    Ok(t)
}

/// `subroutine update_ato(ierrout)`: refresh bonds, rotational groups and atom
/// parameters of every molecule from its molecule type's template `.ato` file.
fn update_ato(st: &mut State, ctx: &ToolContext, sys: &SysCmds) -> Result<()> {
    let mut icurrenttemplate = 0;
    st.mol.ngrrot = 0;
    st.mol.nrottot = 0;
    let mut ig: i32 = 0;
    let mut tpl = Template::default();
    for im in 1..=st.mol.nmol {
        let imoltype = st.mol.mtype[im];
        if imoltype != icurrenttemplate {
            tpl.natom = 0;
            icurrenttemplate = imoltype;
            let name = st.mol.moltemplate[imoltype].trim_end().to_string();
            if !name.is_empty() {
                let fname_template = format!("{}.ato", name);
                st.log(&format!("\nupdate_ato> {}", fname_template));
                let fnametemplate = get_full_filename(sys, &ctx.dir, &fname_template);
                match TextFile::open(std::path::Path::new(&fnametemplate)) {
                    Ok(mut f) => {
                        tpl = read_template(st, &mut f)?;
                    }
                    Err(_) => {
                        st.log(&format!("\nupdate_ato> Could not open template file :{}", fnametemplate));
                    }
                }
            }
        }
        let m = &mut st.mol;
        if m.nato[im] == tpl.natom {
            let nadd = m.naf[im] - 1;
            for ia in m.naf[im]..=m.nal[im] {
                let itref = (ia - nadd - 1) as usize;
                if m.alabel[ia] == tpl.alabel[itref] {
                    let nb = tpl.iabond[itref].len();
                    m.nbond[ia] = nb as i32;
                    m.ensure_bonds(nb);
                    for ja in 1..=nb {
                        m.iabond[(ia, ja)] = tpl.iabond[itref][ja - 1] + nadd;
                        m.tbond[(ia, ja)] = tpl.tbond[itref][ja - 1];
                    }
                }
            }
            m.ngrot[im] = tpl.rotgrp.len() as i32;
            m.ngrrot += m.ngrot[im];
            let n = m.ngrrot as usize;
            m.ensure_rotgrps(n);
            for ir in 0..tpl.rotgrp.len() {
                ig += 1;
                m.rotgrp[ig] = tpl.rotgrp[ir].clone();
                m.labaxi[(ig, 1)] = tpl.labaxi[ir][0] + nadd;
                m.labaxi[(ig, 2)] = tpl.labaxi[ir][1] + nadd;
                m.nrot[ig] = tpl.labrot[ir].len() as i32;
                m.rotf[ig] = m.nrottot + 1;
                m.nrottot += m.nrot[ig];
                m.rotl[ig] = m.nrottot;
                let n = m.nrottot as usize;
                m.ensure_labrot(n);
                for j in 1..=m.nrot[ig] {
                    let labref = m.rotf[ig] + j - 1;
                    m.labrot[labref] = tpl.labrot[ir][(j - 1) as usize] + nadd;
                }
                m.ndihedralang[ig] = 0;
            }
        } else {
            m.ngrrot += m.ngrot[im];
            ig = m.ngrrot;
        }
    }
    Ok(())
}
