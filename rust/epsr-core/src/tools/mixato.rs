//! `mixato` — translation of `fortran_src/mixato.f`.
//!
//! STATUS: translated; verified against the Fortran (EPSR-test-ex31, ex32).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> mixato [<n files>]`
//! Answers read from standard input: number of `.ato` files (if not given),
//! then for each file its name (via `search_for_files`) and the number of
//! molecules wanted, then the atomic number density and the output name.
//!
//! The first molecule of each single-molecule-type `.ato` file is replicated
//! the requested number of times (all at that molecule's position); the box
//! size follows from the density.  Atom types, molecule templates and the
//! random-number state of the last file read are carried over.

use super::makeato::setup_common::{
    command_text, eof_err, get_full_filename, read_int5, read_real5, read_text, search_for_files, strip_first_dot, SysCmds,
};
use super::ToolContext;
use crate::ato::{cosd_intrinsic as cosd, sind_intrinsic as sind};
use crate::error::{EpsrError, Result};
use crate::fortran::listread::read_int;
use crate::fortran::text::parse;
use crate::state::MMOL;

/// The saved first molecule of one input file.
#[derive(Clone, Debug, Default)]
pub struct SavedMol {
    pub nato: i32,
    pub xadd: f32,
    pub yadd: f32,
    pub zadd: f32,
    pub phix: f32,
    pub phiy: f32,
    pub phiz: f32,
    pub tethered: bool,
    pub alabel: Vec<String>,
    pub itype: Vec<i32>,
    pub xa: Vec<f32>,
    pub ya: Vec<f32>,
    pub za: Vec<f32>,
    pub iabond: Vec<Vec<i32>>,
    pub tbond: Vec<Vec<f32>>,
    pub rotgrp: Vec<String>,
    pub labaxi: Vec<[i32; 2]>,
    pub labrot: Vec<Vec<i32>>,
}

/// Run the `mixato` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    let sys = SysCmds::load();
    let (sometext, w) = command_text("mixato", ctx);
    let mut nmix: i32 = 0;
    let mut ierr = 1;
    if w.nwords() > 1 {
        if let Ok(v) = read_int(w.word(&sometext, 2)) {
            nmix = v;
            ierr = 0;
        }
    }
    while ierr != 0 {
        st.log("\nmixato> How many .ATO files do you want to mix? ");
        let mut r = ctx.stdin.list();
        match r.int() {
            Ok(v) => {
                if let Some(v) = v {
                    nmix = v;
                }
                ierr = 0;
            }
            Err(crate::fortran::IoErr::Eof) => return Err(eof_err(crate::fortran::IoErr::Eof)),
            Err(_) => ierr = 1,
        }
    }
    let mut nmincl = vec![0i32; nmix.max(0) as usize + 1];
    let mut saved: Vec<SavedMol> = vec![SavedMol::default(); nmix.max(0) as usize + 1];
    let mut nm1 = 0;
    st.mol.ntype = 0;
    st.mol.nmoltype = 0;
    let mut ntypefirst = st.mol.ntype;
    let mut nmoltypefirst = st.mol.nmoltype;
    let mut stepmis: f32 = 0.0;
    let mut stepris: f32 = 0.0;
    let mut steprs: f32 = 0.0;
    let mut stepms: f32 = 0.0;
    let mut tempbs: f32 = 0.0;
    let mut tempdihs: f32 = 0.0;
    let mut imix = 0;
    while imix < nmix {
        imix += 1;
        st.log(&format!("\n\nmixato> Search for .ato file {:>2}", imix));
        let mut nwords = 0;
        let mut text = String::new();
        while nwords <= 0 {
            let (t, iaccept) = search_for_files(&st, ctx, &sys, ".ato")?;
            if iaccept != 'y' {
                st.log("\nmixato> No .ato files found in the current directory - try another directory.\n");
                return Ok(());
            }
            nwords = parse(&t, 5).nwords();
            if nwords <= 0 {
                st.log("\nmixato> No .ato file specified - try again!\n");
            }
            text = t;
        }
        let w = parse(&text, 5);
        let readfile = strip_first_dot(w.word(&text, 1));
        let fnameato = format!("{}.ato", get_full_filename(&sys, &ctx.dir, &readfile));
        st.mol.sizefactor = 1.0;
        match st.minitc(&fnameato) {
            Ok(()) => {
                let m = &st.mol;
                let ntypefile = m.ntype - ntypefirst;
                let nmoltypefile = m.nmoltype - nmoltypefirst;
                if nmoltypefile > 1 {
                    st.log(&format!(
                        "\nmixato> {}.ato contains more than 1 type of molecule\n\nmixato> This is not allowed - start again",
                        readfile
                    ));
                    return Ok(());
                }
                ntypefirst = m.ntype;
                nmoltypefirst = m.nmoltype;
                st.log(&format!("\n{:>4} new atom types in file {}\n", ntypefile, fnameato));
                st.log("\nmixato> How many of these molecules do you want in the mixture? ");
                nmincl[imix as usize] = read_int5(ctx)?;
                nm1 += nmincl[imix as usize];
                if nm1 as usize > MMOL {
                    st.log(&format!(" No. of molecules requested .gt. than allowed dimension:{:>6}", MMOL));
                    return Ok(());
                }
                let m = &mut st.mol;
                if m.stepmi > stepmis {
                    stepmis = m.stepmi;
                }
                if m.stepri > stepris {
                    stepris = m.stepri;
                }
                if m.stepr > steprs {
                    steprs = m.stepr;
                }
                if m.stepm > stepms {
                    stepms = m.stepm;
                }
                if m.tempb > tempbs {
                    tempbs = m.tempb;
                }
                if m.tempdih > tempdihs {
                    tempdihs = m.tempdih;
                }
                let mut s = SavedMol {
                    nato: m.nato[1],
                    xadd: m.xadd[1],
                    yadd: m.yadd[1],
                    zadd: m.zadd[1],
                    phix: m.phix[1],
                    phiy: m.phiy[1],
                    phiz: m.phiz[1],
                    ..Default::default()
                };
                for ia in m.naf[1]..=m.nal[1] {
                    s.alabel.push(m.alabel[ia].clone());
                    s.itype.push(m.itype[ia]);
                    s.xa.push(m.xa[ia]);
                    s.ya.push(m.ya[ia]);
                    s.za.push(m.za[ia]);
                    let mut ib_v = Vec::new();
                    let mut tb_v = Vec::new();
                    for ib in 1..=m.nbond[ia] {
                        ib_v.push(m.iabond[(ia, ib)]);
                        tb_v.push(m.tbond[(ia, ib)]);
                    }
                    s.iabond.push(ib_v);
                    s.tbond.push(tb_v);
                }
                for ig in 1..=m.ngrot[1] {
                    s.rotgrp.push(m.rotgrp[ig].clone());
                    s.labaxi.push([m.labaxi[(ig, 1)], m.labaxi[(ig, 2)]]);
                    let mut v = Vec::new();
                    for j in 1..=m.nrot[ig] {
                        let labref = m.rotf[ig] + j - 1;
                        v.push(m.labrot[labref]);
                    }
                    if m.rotgrp[ig] == "DIH" {
                        return crate::error::unsupported("mixato: DIH rotational groups");
                    }
                    s.labrot.push(v);
                }
                saved[imix as usize] = s;
                // does the template file exist?
                let nt = m.nmoltype;
                let tmpl = m.moltemplate[nt].trim_end().to_string();
                if !std::path::Path::new(&format!("{}.ato", tmpl)).is_file() {
                    st.log(&format!(" mixato> Did not find template file: {}.ato", tmpl));
                    st.mol.moltemplate[nt] = readfile.clone();
                }
                st.log(&format!("\n  mixato> moltemplate file: {}.ato", st.mol.moltemplate[nt].trim_end()));
            }
            Err(EpsrError::Io(_)) => {
                st.mol.ntype = ntypefirst;
                st.mol.nmoltype = nmoltypefirst;
                st.log(&format!("\nProblems with .ato file {} Try again.\n", fnameato));
                imix -= 1;
            }
            Err(e) => return Err(e),
        }
    }
    // ---- build the mixture
    let mut ia: i32 = 0;
    let mut im: i32 = 0;
    let mut ig: i32 = 0;
    {
        let m = &mut st.mol;
        m.natom = 0;
        m.natomr = 0;
        m.natomcg = 0.0;
        m.natomrcg = 0.0;
        m.ngrrot = 0;
        m.nrottot = 0;
    }
    for imix in 1..=nmix {
        let s = saved[imix as usize].clone();
        for _it in 1..=nmincl[imix as usize] {
            if (im as usize) < MMOL {
                im += 1;
            } else {
                st.log("\nmixato> Tried to define more molecules than allowed by dimension\n");
                return Ok(());
            }
            st.mol.ensure_atoms((ia + s.nato) as usize);
            {
                let m = &mut st.mol;
                m.xadd[im] = s.xadd;
                m.yadd[im] = s.yadd;
                m.zadd[im] = s.zadd;
                m.nato[im] = s.nato;
                m.naf[im] = ia + 1;
            }
            let nadd = ia;
            for is in 1..=s.nato {
                ia += 1;
                let k = (is - 1) as usize;
                st.mol.alabel[ia] = s.alabel[k].clone();
                st.ftype(ia);
                let m = &mut st.mol;
                let it = m.itype[ia];
                m.natomcg += m.natomsphere[it] * m.pseudofacmol[it];
                if m.qtype[ia] != 1 {
                    m.natomr += 1;
                    m.natomrcg += m.natomsphere[it] * m.pseudofacmol[it];
                }
                m.xa[ia] = s.xa[k];
                m.ya[ia] = s.ya[k];
                m.za[ia] = s.za[k];
                let nb = s.iabond[k].len();
                m.nbond[ia] = nb as i32;
                m.ensure_bonds(nb);
                for ib in 1..=nb {
                    m.iabond[(ia, ib)] = s.iabond[k][ib - 1] + nadd;
                    m.tbond[(ia, ib)] = s.tbond[k][ib - 1];
                }
            }
            st.mol.nal[im] = ia;
            let ngrot = s.rotgrp.len() as i32;
            st.mol.ngrot[im] = ngrot;
            if ngrot > 0 {
                st.mol.ngrrot += ngrot;
                let n = st.mol.ngrrot as usize;
                st.mol.ensure_rotgrps(n);
                for ih in 1..=ngrot {
                    ig += 1;
                    let h = (ih - 1) as usize;
                    let m = &mut st.mol;
                    m.mrot[ig] = im;
                    m.rotgrp[ig] = s.rotgrp[h].clone();
                    m.labaxi[(ig, 1)] = s.labaxi[h][0] + nadd;
                    m.labaxi[(ig, 2)] = s.labaxi[h][1] + nadd;
                    m.nrot[ig] = s.labrot[h].len() as i32;
                    m.rotf[ig] = m.nrottot + 1;
                    m.nrottot += m.nrot[ig];
                    m.rotl[ig] = m.nrottot;
                    let n = m.nrottot as usize;
                    m.ensure_labrot(n);
                    for j in 1..=m.nrot[ig] {
                        let labref = m.rotf[ig] + j - 1;
                        m.labrot[labref] = s.labrot[h][(j - 1) as usize] + nadd;
                    }
                }
            }
        }
    }
    st.mol.nmol = im;
    st.mol.natom = ia;
    st.mol.ngrrot = ig;
    st.log("\nmixato> Give atomic number density (per A**3) of mixture: ");
    let rhoatom = read_real5(ctx)?;
    let m = &mut st.mol;
    m.dmax = m.natomrcg / rhoatom;
    m.dmax = m.dmax.powf(1.0f32 / 3.0f32);
    m.dmaxmol = m.sizefactor * m.dmax;
    let l = &mut st.lat;
    l.cubic = true;
    l.alat = m.dmax;
    l.blat = m.dmax;
    l.clat = m.dmax;
    l.phib = 90.0;
    l.thetac = 0.0;
    l.phic = 0.0;
    l.alatx = l.alat;
    l.blatx = l.blat * cosd(l.phib);
    l.blaty = l.blat * sind(l.phib);
    l.clatx = l.clat * sind(l.thetac) * cosd(l.phic);
    l.claty = l.clat * sind(l.thetac) * sind(l.phic);
    l.clatz = l.clat * cosd(l.thetac);
    m.stepmi = stepmis;
    m.stepri = 0.3;
    m.stepr = 0.3;
    m.stepm = 1.0;
    m.tempb = tempbs;
    m.tempdih = tempdihs;
    let _ = (stepris, steprs, stepms);
    let mut text;
    loop {
        st.log("\nmixato> Type name of file to put mixture in: ");
        text = read_text(ctx)?;
        if parse(&text, 5).nwords() > 0 {
            break;
        }
        st.log("\nmixato> No .ato file specified - try again!\n");
    }
    let w = parse(&text, 5);
    let readfile = strip_first_dot(w.word(&text, 1));
    let fnameato = format!("{}{}.ato", ctx.dir, readfile);
    st.moloutc(&fnameato)
}
