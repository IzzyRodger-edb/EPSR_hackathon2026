//! `addato` — translation of `fortran_src/addato.f`.
//!
//! STATUS: translated; verified against the Fortran (EPSR-test-ex36).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> addato [<n files>]`
//! Answers read from standard input: number of `.ato` files to add (if not
//! given), for each its name (via `search_for_files`) and the number of
//! molecules, then the name of the `.ato` file to add them to, and the output
//! name.
//!
//! The first molecule of each added file is inserted the requested number of
//! times at random positions with random orientations into the target box
//! (whose size is unchanged).
//!
//! NOT PORTED: tethered molecules (the `.ato` reader of the port rejects
//! them), so the overlap test against tethered molecules never runs.

use super::makeato::setup_common::{
    command_text, eof_err, get_full_filename, read_int5, read_text, rotatm, search_for_files, strip_first_dot, strip_last_dot, SysCmds,
};
use super::mixato::SavedMol;
use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::listread::read_int;
use crate::fortran::text::parse;
use crate::neighbour::minimagec;
use crate::state::MMOL;

/// Max of the step sizes seen so far (`stepmis`, ...).
#[derive(Default)]
struct MaxSteps {
    stepmi: f32,
    stepri: f32,
    stepr: f32,
    stepm: f32,
    tempb: f32,
    tempdih: f32,
}

impl MaxSteps {
    fn update(&mut self, m: &crate::state::Mol) {
        if m.stepmi > self.stepmi {
            self.stepmi = m.stepmi;
        }
        if m.stepri > self.stepri {
            self.stepri = m.stepri;
        }
        if m.stepr > self.stepr {
            self.stepr = m.stepr;
        }
        if m.stepm > self.stepm {
            self.stepm = m.stepm;
        }
        if m.tempb > self.tempb {
            self.tempb = m.tempb;
        }
        if m.tempdih > self.tempdih {
            self.tempdih = m.tempdih;
        }
    }
}

/// Ask for a `.ato` file (`search_for_files` loop); `None` if none accepted.
fn ask_ato(st: &crate::state::State, ctx: &mut ToolContext, sys: &SysCmds) -> Result<Option<String>> {
    loop {
        let (t, iaccept) = search_for_files(st, ctx, sys, ".ato")?;
        if iaccept != 'y' {
            st.log("\naddato> No .ato files found in the current directory - try another directory.\n");
            return Ok(None);
        }
        let w = parse(&t, 5);
        if w.nwords() > 0 {
            return Ok(Some(w.word(&t, 1).to_string()));
        }
        st.log("\naddato> No .ato file specified - try again!\n");
    }
}

/// Run the `addato` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    let sys = SysCmds::load();
    let (sometext, w) = command_text("addato", ctx);
    let mut nmix: i32 = 0;
    let mut ierr = 1;
    if w.nwords() > 1 {
        if let Ok(v) = read_int(w.word(&sometext, 2)) {
            nmix = v;
            ierr = 0;
        }
    }
    while ierr != 0 {
        st.log("\naddato> How many .ATO files do you want to add? ");
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
    let mut nmincl: Vec<i32> = vec![0; nmix.max(0) as usize + 1];
    let mut saved: Vec<SavedMol> = vec![SavedMol::default(); nmix.max(0) as usize + 1];
    let mut nm1 = 0;
    st.mol.ntype = 0;
    st.mol.nmoltype = 0;
    let mut ntypefirst = st.mol.ntype;
    let mut nmoltypefirst = st.mol.nmoltype;
    let mut mx = MaxSteps::default();
    let mut imix: i32 = 0;
    while imix < nmix {
        st.log("\n\naddato> Search for .ato file to be ADDED ");
        let name = match ask_ato(&st, ctx, &sys)? {
            Some(n) => n,
            None => return Ok(()),
        };
        let readfile = strip_last_dot(&name);
        let fnameato = format!("{}.ato", get_full_filename(&sys, &ctx.dir, &readfile));
        st.mol.sizefactor = 1.0;
        match st.minitc(&fnameato) {
            Ok(()) => {
                let ntypefile = st.mol.ntype - ntypefirst;
                let nmoltypefile = st.mol.nmoltype - nmoltypefirst;
                if nmoltypefile > 1 {
                    st.log(&format!(
                        "\naddato> {}.ato contains more than 1 type of molecule\n\naddato> This is not allowed - try again",
                        readfile
                    ));
                    continue;
                }
                ntypefirst = st.mol.ntype;
                nmoltypefirst = st.mol.nmoltype;
                st.log(&format!("\n{:>4} new atom types in file {}\n", ntypefile, fnameato));
                let nt = st.mol.nmoltype;
                if st.mol.moltethered[nt] {
                    return crate::error::unsupported("addato: tethered molecule types");
                }
                st.log("\naddato> How many of these molecules do you want in the mixture? ");
                let nmoladd = read_int5(ctx)?;
                let immax = 1;
                let mut im = 0;
                let mut irgref = 0;
                while im < immax {
                    im += 1;
                    imix += 1;
                    if imix as usize >= saved.len() {
                        saved.resize(imix as usize + 10, SavedMol::default());
                        nmincl.resize(imix as usize + 10, 0);
                    }
                    nmincl[imix as usize] = nmoladd;
                    nm1 += nmoladd;
                    let m = &st.mol;
                    let mut s = SavedMol {
                        nato: m.nato[im],
                        xadd: m.xadd[im],
                        yadd: m.yadd[im],
                        zadd: m.zadd[im],
                        phix: m.phix[im],
                        phiy: m.phiy[im],
                        phiz: m.phiz[im],
                        tethered: m.tethered[im],
                        ..Default::default()
                    };
                    let iaadd = m.naf[im] - 1;
                    for ir in 1..=s.nato {
                        let ia = iaadd + ir;
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
                    for ig in 1..=m.ngrot[im] {
                        irgref += 1;
                        s.rotgrp.push(m.rotgrp[irgref].clone());
                        s.labaxi.push([m.labaxi[(irgref, 1)], m.labaxi[(irgref, 2)]]);
                        // the Fortran copies the atom list with index ig (= irgref here)
                        let mut v = Vec::new();
                        for j in 1..=m.nrot[ig] {
                            let labref = m.rotf[ig] + j - 1;
                            v.push(m.labrot[labref]);
                        }
                        v.resize(m.nrot[irgref].max(0) as usize, 0);
                        if m.rotgrp[ig] == "DIH" {
                            return crate::error::unsupported("addato: DIH rotational groups");
                        }
                        s.labrot.push(v);
                    }
                    saved[imix as usize] = s;
                }
                // open(35,file=moltemplate(nmoltype)//'.ato'): the name is the full
                // character*256 variable (blank padded) followed by '.ato'
                let tmpl = st.mol.moltemplate[nt].clone();
                let padded = if tmpl.len() < 256 {
                    format!("{}{}.ato", tmpl, " ".repeat(256 - tmpl.len()))
                } else {
                    format!("{}.ato", &tmpl[..256])
                };
                if !std::path::Path::new(padded.trim_end()).is_file() {
                    st.mol.moltemplate[nt] = readfile.clone();
                }
                st.log(&format!("\n  addato> moltemplate file: {}.ato", st.mol.moltemplate[nt].trim_end()));
                mx.update(&st.mol);
                if nm1 as usize > MMOL {
                    st.log(&format!(" No. of molecules requested .gt. than allowed dimension:{:>6}", MMOL));
                    return Ok(());
                }
            }
            Err(EpsrError::Io(_)) => {
                st.mol.ntype = ntypefirst;
                st.mol.nmoltype = nmoltypefirst;
                st.log(&format!("\naddato> Problems with .ato file {} Try again.\n", fnameato));
            }
            Err(e) => return Err(e),
        }
    }
    // ---- the .ato file the molecules are added to
    while imix < nmix + 1 {
        st.log("\n\naddato> Search for .ato file to be ADDED to");
        let name = match ask_ato(&st, ctx, &sys)? {
            Some(n) => n,
            None => return Ok(()),
        };
        let readfile = strip_last_dot(&name);
        let fnameato = format!("{}.ato", get_full_filename(&sys, &ctx.dir, &readfile));
        st.mol.sizefactor = 1.0;
        match st.minitc(&fnameato) {
            Ok(()) => {
                imix += 1;
                mx.update(&st.mol);
            }
            Err(EpsrError::Io(_)) => {
                st.mol.ntype = ntypefirst;
                st.mol.nmoltype = nmoltypefirst;
                st.log(&format!("\naddato> Problems with .ato file {} Try again.\n", fnameato));
            }
            Err(e) => return Err(e),
        }
    }
    let ntn = (1..=st.mol.nmol).filter(|&im| st.mol.tethered[im]).count();
    st.log(&format!(" addato> Number of tethered molecules at start = {}", ntn));
    st.normalise_coordinates();
    let mut im = st.mol.nmol;
    let mut ia = st.mol.natom;
    let mut ig = st.mol.ngrrot;
    for imix in 1..=nmix {
        let s = saved[imix as usize].clone();
        for _iloop in 1..=nmincl[imix as usize] {
            if (im as usize) < MMOL {
                im += 1;
            } else {
                st.log("\naddato> Tried to define more molecules than allowed by dimension\n");
                return Ok(());
            }
            if s.tethered {
                return crate::error::unsupported("addato: tethered molecules");
            }
            st.mol.tethered[im] = false;
            {
                // one insertion attempt suffices: no tethered molecules to overlap with
                let fa = st.rng.ran1() - 0.5;
                let fb = st.rng.ran1() - 0.5;
                let fc = st.rng.ran1() - 0.5;
                let l = &st.lat;
                let m = &mut st.mol;
                m.zadd[im] = fc * l.clatz;
                m.yadd[im] = fb * l.blaty + fc * l.claty;
                m.xadd[im] = fa * l.alatx + fb * l.blatx + fc * l.clatx;
                m.phix[im] = 0.0;
                m.phiy[im] = 0.0;
                m.phiz[im] = 0.0;
                let (mut x, mut y, mut z) = (m.xadd[im], m.yadd[im], m.zadd[im]);
                minimagec(&mut x, &mut y, &mut z);
                m.xadd[im] = x;
                m.yadd[im] = y;
                m.zadd[im] = z;
            }
            st.mol.ensure_atoms((ia + s.nato) as usize);
            st.mol.nato[im] = s.nato;
            let nadd = ia;
            st.mol.naf[im] = ia + 1;
            let dmaxmol = st.mol.dmaxmol;
            for is in 1..=s.nato {
                ia += 1;
                let k = (is - 1) as usize;
                st.mol.alabel[ia] = s.alabel[k].clone();
                st.ftype(ia);
                let m = &mut st.mol;
                if m.qtype[ia] != 1 {
                    m.natomr += 1;
                }
                m.xa[ia] = s.xa[k] / dmaxmol;
                m.ya[ia] = s.ya[k] / dmaxmol;
                m.za[ia] = s.za[k] / dmaxmol;
                let nb = s.iabond[k].len();
                m.nbond[ia] = nb as i32;
                m.ensure_bonds(nb);
                for ib in 1..=nb {
                    m.iabond[(ia, ib)] = s.iabond[k][ib - 1] + nadd;
                    m.tbond[(ia, ib)] = s.tbond[k][ib - 1] / dmaxmol;
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
            if !st.mol.tethered[im] {
                let steprs = st.mol.stepr;
                if steprs > 0.0 {
                    st.mol.stepr = 1.0;
                }
                let s1 = st.mol.stepr;
                rotatm(&mut st, im, s1);
                st.mol.stepr = steprs;
            }
        }
    }
    st.mol.nmol = im;
    st.mol.natom = ia;
    st.unnormalise_coordinates();
    let m = &mut st.mol;
    m.stepmi = mx.stepmi;
    m.stepri = 0.3;
    m.stepr = 0.3;
    m.stepm = 1.0;
    m.tempb = mx.tempb;
    m.tempdih = mx.tempdih;
    let _ = (mx.stepri, mx.stepr, mx.stepm);
    let mut text;
    loop {
        st.log("\naddato> Type name of file to put mixture in: ");
        text = read_text(ctx)?;
        if parse(&text, 5).nwords() > 0 {
            break;
        }
        st.log("\naddato> No .ato file specified - try again!\n");
    }
    let w = parse(&text, 5);
    let readfile = strip_first_dot(w.word(&text, 1));
    let fnameato = format!("{}{}.ato", ctx.dir, readfile);
    st.moloutc(&fnameato)
}
