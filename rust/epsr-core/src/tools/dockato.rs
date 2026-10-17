//! `dockato` — translation of `fortran_src/dockato.f` (program `standalone`,
//! subroutine `dockato`, module `dockato_arrays` and the routines
//! `rotatem`, `random_rotatem`, `search_unsaturated`, `find_valence`,
//! `test_overlap`, `find_neighbour`, `bond_molecule`, `replace_references`,
//! `add_neighbours`, `remove_atoms`, `clean_atom_types`), plus
//! `get_direct_bonds`, `get_rotgrp` and `centreofmass` from
//! `ato_manipulation_rev4.f`.
//!
//! A specified number of copies of the first molecule of the docking `.ato`
//! file are docked onto sites of the docked `.ato` file: atom `atomd` of the
//! new molecule is placed on atom `atomv`, the molecule is rotated so that its
//! `atomdr` neighbour lies on the `atomvr` neighbour of `atomv` (same bond
//! length required), an optional second (dihedral) rotation is made using
//! `atomdrr`/`atomvrr`, overlaps are tested, and the two molecules are bonded
//! (atoms `atomd` and `atomvr` are removed).
//!
//! Invocation (same as the Fortran program):
//! `epsr-rs <folder> dockato [dockingfile atomd atomdr atomdrr ndock radmax
//! dockedfile atomv atomvr atomvrr dminrr dmaxrr overlaps outputfile]`.
//! Missing trailing parameters are read from standard input (`ctx.stdin`),
//! except the two `.ato` file names: the interactive file search
//! (`search_for_files`, which runs a shell `ls`) is not ported, so both file
//! names must be given on the command line.

use super::ToolContext;
use crate::ato::{scapro, vecpro};
use crate::error::{EpsrError, Result};
use crate::fortran::fmath::sign;
use crate::fortran::listread::read_real;
use crate::fortran::text::{findex, fstr, parse};
use crate::fortran::{F1, F2, F3};
use crate::neighbour::minimagec;
use crate::state::{State, MCOMP, NEAR};

/// Module `dockato_arrays`: the saved first molecule of the docking file.
#[derive(Default)]
struct Da {
    natos: i32,
    ngrots: i32,
    nadd: i32,
    xadds: f32,
    yadds: f32,
    zadds: f32,
    amassts: f32,
    phixs: f32,
    phiys: f32,
    phizs: f32,
    nbondss: F1<i32>,
    ndihedralangs: F1<i32>,
    iabonds: F2<i32>,
    labaxis: F2<i32>,
    xas: F1<f32>,
    yas: F1<f32>,
    zas: F1<f32>,
    tbonds: F2<f32>,
    overlapfactor: f32,
    alabels: F1<String>,
    rotgrps: F1<String>,
}

fn word(ctx: &ToolContext, nt: i32) -> String {
    ctx.args[(nt - 2) as usize].clone()
}

fn eof<T>(what: &str) -> Result<T> {
    Err(EpsrError::Input(format!("dockato: end of standard input while reading {what}")))
}

/// List-directed `read(5,*) x` of one character item (truncated to `len`).
fn read_word5(ctx: &mut ToolContext, len: usize, what: &str) -> Result<String> {
    let mut r = ctx.stdin.list();
    match r.word() {
        Ok(Some(w)) => Ok(fstr(&w, len)),
        Ok(None) => Ok(String::new()),
        Err(_) => eof(what),
    }
}
fn read_int5(ctx: &mut ToolContext, v: &mut i32, what: &str) -> Result<()> {
    let mut r = ctx.stdin.list();
    match r.get_i(v) {
        Ok(()) => Ok(()),
        Err(crate::fortran::IoErr::Eof) => eof(what),
        Err(e) => Err(e.into()),
    }
}
fn read_real5(ctx: &mut ToolContext, v: &mut f32, what: &str) -> Result<()> {
    let mut r = ctx.stdin.list();
    match r.get_r(v) {
        Ok(()) => Ok(()),
        Err(crate::fortran::IoErr::Eof) => eof(what),
        Err(e) => Err(e.into()),
    }
}

/// Make sure the dihedral arrays of the rotational groups hold `n` groups.
fn ensure_dih(st: &mut State, n: usize) {
    st.mol.ensure_rotgrps(n);
    let m = &mut st.mol;
    let n = n.max(m.rotgrp.len());
    m.dihedralang.grow(NEAR, n);
    let l = &mut m.labrotdih;
    if l.n3 < n || l.n1 < 2 || l.n2 < NEAR {
        let mut nl: F3<i32> = F3::new(2, NEAR, n);
        for k in 1..=l.n3 {
            for j in 1..=l.n2.min(NEAR) {
                for i in 1..=l.n1.min(2) {
                    nl[(i, j, k)] = l[(i, j, k)];
                }
            }
        }
        *l = nl;
    }
}

/// `nbond(ia)`, where `nbond(0)` is the preceding COMMON element
/// `nbonds(matom)`, which is never set (= 0).
fn nbond_at(st: &State, ia: i32) -> i32 {
    if ia == 0 {
        0
    } else {
        st.mol.nbond[ia]
    }
}

/// Run the `dockato` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    {
        let m = &mut st.mol;
        m.natomsphere.grow(MCOMP);
        m.rhosphere.grow(MCOMP);
        m.radminsphere.grow(MCOMP);
        m.radmaxsphere.grow(MCOMP);
        m.pseudofacmol.grow(MCOMP);
    }
    let mut da = Da::default();
    let nwords = 1 + ctx.args.len() as i32;

    let mut dockingfile = String::new();
    let mut dockedfile = String::new();
    let mut outputfile = String::new();
    let mut atomd = String::new();
    let mut atomdr = String::new();
    let mut atomdrr = String::new();
    let mut atomv = String::new();
    let mut atomvr = String::new();
    let mut atomvrr = String::new();
    let mut ilendockingfile = 0usize;
    let mut ilendockedfile = 0usize;
    let mut ilenatomd = 0usize;
    let mut ilenatomdr = 0usize;
    let mut ilenatomdrr = 0usize;
    let mut ilenatomv = 0usize;
    let mut ilenatomvr = 0usize;
    let mut ilenatomvrr = 0usize;
    let mut ilendminrr = 0usize;
    let mut ilendmaxrr = 0usize;
    let mut ilenoutputfile = 0usize;
    let mut ndock: i32 = 0;
    let mut radmax: f32 = 0.0;
    let mut dminrr: f32 = 0.0;
    let mut dmaxrr: f32 = 0.0;
    let mut nooverlapstext = String::from("<undefined>");

    let mut ntest = 1;
    if nwords > ntest {
        let w = word(ctx, ntest + 1);
        ilendockingfile = w.len();
        dockingfile = fstr(&w, 80);
    }
    ntest += 1;
    if nwords > ntest {
        let w = word(ctx, ntest + 1);
        ilenatomd = w.len();
        atomd = fstr(&w, 3);
    }
    ntest += 1;
    if nwords > ntest {
        let w = word(ctx, ntest + 1);
        ilenatomdr = w.len();
        atomdr = fstr(&w, 3);
    }
    ntest += 1;
    if nwords > ntest {
        let w = word(ctx, ntest + 1);
        ilenatomdrr = w.len();
        atomdrr = fstr(&w, 3);
    }
    ntest += 1;
    if nwords > ntest {
        ndock = crate::fortran::listread::read_int(&word(ctx, ntest + 1))?;
    }
    ntest += 1;
    if nwords > ntest {
        radmax = read_real(&word(ctx, ntest + 1))?;
    }
    ntest += 1;
    if nwords > ntest {
        let w = word(ctx, ntest + 1);
        ilendockedfile = w.len();
        dockedfile = fstr(&w, 80);
    }
    ntest += 1;
    if nwords > ntest {
        let w = word(ctx, ntest + 1);
        ilenatomv = w.len();
        atomv = fstr(&w, 3);
    }
    ntest += 1;
    if nwords > ntest {
        let w = word(ctx, ntest + 1);
        ilenatomvr = w.len();
        atomvr = fstr(&w, 3);
    }
    ntest += 1;
    if nwords > ntest {
        let w = word(ctx, ntest + 1);
        ilenatomvrr = w.len();
        atomvrr = fstr(&w, 3);
    }
    ntest += 1;
    if nwords > ntest {
        let w = word(ctx, ntest + 1);
        dminrr = read_real(&w)?;
        ilendminrr = w.len();
    }
    ntest += 1;
    if nwords > ntest {
        let w = word(ctx, ntest + 1);
        dmaxrr = read_real(&w)?;
        ilendmaxrr = w.len();
    }
    ntest += 1;
    if nwords > ntest {
        nooverlapstext = fstr(&word(ctx, ntest + 1), 20);
    }
    ntest += 1;
    if nwords > ntest {
        let w = word(ctx, ntest + 1);
        ilenoutputfile = w.len();
        outputfile = fstr(&w, 80);
    }
    let _ = radmax;

    let nmix = 2;
    st.mol.ntype = 0;
    st.mol.nmoltype = 0;
    let mut ntypefirst = 0;
    let mut imix = 0;
    while imix < nmix {
        imix += 1;
        let text = if imix == 1 {
            if ilendockingfile == 0 {
                st.log("\n\ndockato> Search for .ato file containing docking molecule");
                return crate::error::unsupported("dockato: interactive .ato file search (search_for_files) - give the docking file on the command line");
            }
            dockingfile.clone()
        } else {
            if ilendockedfile == 0 {
                st.log("\n\ndockato> Search for .ato file to be docked onto");
                return crate::error::unsupported("dockato: interactive .ato file search (search_for_files) - give the docked file on the command line");
            }
            dockedfile.clone()
        };
        let w = parse(&text, 5);
        let mut readfile = w.span(&text, 1, 1);
        let ilendot = findex(&readfile, ".") - 1;
        if ilendot > 0 {
            readfile = readfile[..ilendot as usize].to_string();
        }
        let fnameato = format!("{}.ato", st.sym.full(&readfile).to_string_lossy());

        st.mol.nmoltype = 0;
        st.mol.sizefactor = 1.0;
        match st.minitc(&fnameato) {
            Ok(()) => {}
            Err(EpsrError::Io(e)) => {
                st.log(&format!("\nProblems with .ato file {fnameato} Try again.\n"));
                return crate::error::unsupported(&format!(
                    "dockato: {e} (retry via interactive file search is not ported)"
                ));
            }
            Err(e) => return Err(e),
        }
        if st.dist.rho <= 0.0 {
            return crate::error::unsupported("dockato: .ato file with zero density (retry via interactive file search is not ported)");
        }
        let ntypefile = st.mol.ntype - ntypefirst;
        st.log(&format!("\n{:>4} new atom types in file {}\n", ntypefile, fnameato));
        ntypefirst = st.mol.ntype;

        if imix == 1 {
            let m = &st.mol;
            da.xadds = m.xadd[1];
            da.yadds = m.yadd[1];
            da.zadds = m.zadd[1];
            da.phixs = m.phix[1];
            da.phiys = m.phiy[1];
            da.phizs = m.phiz[1];
            da.natos = m.nato[1];
            da.amassts = m.amasst[1];
            let n = da.natos.max(0) as usize;
            da.alabels = F1::new(n);
            da.xas = F1::new(n);
            da.yas = F1::new(n);
            da.zas = F1::new(n);
            da.nbondss = F1::new(n);
            let mut nbmax = 1usize;
            for ia in m.naf[1]..=m.nal[1] {
                nbmax = nbmax.max(m.nbond[ia].max(0) as usize);
            }
            da.iabonds = F2::new(n, nbmax);
            da.tbonds = F2::new(n, nbmax);
            for ia in m.naf[1]..=m.nal[1] {
                da.alabels[ia] = m.alabel[ia].clone();
                da.xas[ia] = m.xa[ia];
                da.yas[ia] = m.ya[ia];
                da.zas[ia] = m.za[ia];
                da.nbondss[ia] = m.nbond[ia];
                if m.nbond[ia] > 0 {
                    for ib in 1..=m.nbond[ia] {
                        da.iabonds[(ia, ib)] = m.iabond[(ia, ib)];
                        da.tbonds[(ia, ib)] = m.tbond[(ia, ib)];
                    }
                }
            }
            da.ngrots = m.ngrot[1];
            let ng = da.ngrots.max(0) as usize;
            da.rotgrps = F1::new(ng);
            da.labaxis = F2::new(ng, 2);
            // In the Fortran ndihedralangs is allocated (reallocate1d_i leaves new
            // values uninitialised) but only set for DIH groups, which minitc of
            // the port rejects.  The gfortran -O2 reference then reads garbage
            // and crashes (index error on dihedralangs(k,ig) of size 0) for any
            // docking molecule with rotational groups.  The port takes the
            // evident intent (0 dihedral angles); verified against a Fortran
            // build patched with `ndihedralangs=0`.
            da.ndihedralangs = F1::new(ng);
            if da.ngrots > 0 {
                st.log(" dockato> NOTE: docking molecule has rotational groups; the original Fortran reads uninitialised ndihedralangs here (assumed 0)");
            }
            if da.ngrots > 0 {
                for ig in 1..=da.ngrots {
                    da.rotgrps[ig] = m.rotgrp[ig].clone();
                    da.labaxis[(ig, 1)] = m.labaxi[(ig, 1)];
                    da.labaxis[(ig, 2)] = m.labaxi[(ig, 2)];
                    if m.rotgrp[ig] == "DIH" {
                        return crate::error::unsupported("dockato: DIH groups in the docking molecule");
                    }
                }
            }
            st.log(&format!(" dockato> natos,ngrots: {} {}", da.natos, da.ngrots));
        }

        let mut ifind = 0;
        let mut ifindr = 0;
        let mut ifindrr = 0;
        while ifind == 0 && ifindr == 0 && ifindrr == 0 {
            let mut prompted = false;
            if imix == 1 {
                if ilenatomd == 0 {
                    atomd = read_word5(ctx, 3, "atomd")?;
                    prompted = true;
                }
                if ilenatomdr == 0 {
                    atomdr = read_word5(ctx, 3, "atomdr")?;
                    prompted = true;
                }
                if ilenatomdrr == 0 {
                    atomdrr = read_word5(ctx, 3, "atomdrr")?;
                    prompted = true;
                }
                if ndock == 0 {
                    read_int5(ctx, &mut ndock, "number of docking molecules")?;
                    prompted = true;
                }
                if radmax == 0.0 {
                    read_real5(ctx, &mut radmax, "maximum distance")?;
                    prompted = true;
                }
                let m = &st.mol;
                let mut ia = 0;
                while ia < m.nato[1] {
                    ia += 1;
                    if m.alabel[ia] == atomd && ifind == 0 {
                        ifind = ia;
                    }
                    if m.alabel[ia] == atomdr && ifindr == 0 {
                        ifindr = ia;
                    }
                    if m.alabel[ia] == atomdrr && ifindrr == 0 {
                        ifindrr = ia;
                    }
                }
            } else {
                if ilenatomv == 0 {
                    atomv = read_word5(ctx, 3, "atomv")?;
                    prompted = true;
                }
                if ilenatomvr == 0 {
                    atomvr = read_word5(ctx, 3, "atomvr")?;
                    prompted = true;
                }
                if ilenatomvrr == 0 {
                    atomvrr = read_word5(ctx, 3, "atomvrr")?;
                    prompted = true;
                }
                if ilendminrr == 0 {
                    read_real5(ctx, &mut dminrr, "minimum dihedral angle")?;
                    prompted = true;
                }
                if ilendmaxrr == 0 {
                    read_real5(ctx, &mut dmaxrr, "maximum dihedral angle")?;
                    prompted = true;
                }
                if findex(&nooverlapstext, "<undefined>") > 0 {
                    match ctx.stdin.read_line() {
                        Ok(l) => nooverlapstext = fstr(&l, 20),
                        Err(_) => return eof("overlap factor"),
                    }
                    prompted = true;
                }
                let m = &st.mol;
                let mut ia = 0;
                while ia < m.natom {
                    ia += 1;
                    if m.alabel[ia] == atomv && ifind == 0 {
                        ifind = ia;
                    }
                    if m.alabel[ia] == atomvr && ifindr == 0 {
                        ifindr = ia;
                    }
                    if m.alabel[ia] == atomvrr && ifindrr == 0 {
                        ifindrr = ia;
                    }
                }
            }
            if ifind == 0 && ifindr == 0 && ifindrr == 0 && !prompted {
                // The Fortran loops forever here.
                return Err(EpsrError::Input("dockato: none of the specified atom labels exist in the .ato file".into()));
            }
        }
    }

    // ---- overlap factor
    let mut imoladd = 0;
    if findex(&nooverlapstext, "yes") > 0 {
        nooverlapstext = "0.0".into();
    } else if findex(&nooverlapstext, "nooverlaps") > 0 {
        nooverlapstext = "1.0".into();
    }
    da.overlapfactor = read_real(&nooverlapstext).unwrap_or(0.0);
    let nooverlaps = da.overlapfactor > 0.0;
    st.log(&format!("dockato> Overlap factor: {:10.5}", da.overlapfactor));

    let mut imv_set: Option<i32> = None;
    let ivcount = 0;
    let mut it = 0;
    let mut nvalencev = 0;
    while it < ndock && ivcount < st.mol.natom {
        it += 1;
        let natomsave = st.mol.natom;
        let natomrsave = st.mol.natomr;
        da.nadd = st.mol.natom;
        let nmolsave = st.mol.nmol;
        let ngrrotsave = st.mol.ngrrot;
        let dmaxsave = st.mol.dmax;

        let mut ia = st.mol.nal[nmolsave];
        let mut im = nmolsave;
        let mut ig = ngrrotsave;
        if (im as usize) < crate::state::MMOL {
            im += 1;
        } else {
            st.log("\ndockato> Tried to define more molecules than allowed by dimension\n");
            return Ok(());
        }
        st.mol.ensure_atoms((natomsave + da.natos) as usize);
        st.mol.ensure_bonds(da.iabonds.n2);
        {
            let m = &mut st.mol;
            m.xadd[im] = da.xadds;
            m.yadd[im] = da.yadds;
            m.zadd[im] = da.zadds;
            m.phix[im] = da.phixs;
            m.phiy[im] = da.phiys;
            m.phiz[im] = da.phizs;
            m.nato[im] = da.natos;
            m.amasst[im] = da.amassts;
            m.natom = natomsave + da.natos;
            m.naf[im] = ia + 1;
        }
        for is in 1..=da.natos {
            ia += 1;
            st.mol.mref[ia] = im;
            st.mol.alabel[ia] = da.alabels[is].clone();
            st.ftype(ia);
            let m = &mut st.mol;
            if m.qtype[ia] != 1 {
                m.natomr += 1;
            }
            m.xa[ia] = da.xas[is];
            m.ya[ia] = da.yas[is];
            m.za[ia] = da.zas[is];
            m.nbond[ia] = da.nbondss[is];
            if m.nbond[ia] > 0 {
                for ib in 1..=m.nbond[ia] {
                    m.iabond[(ia, ib)] = da.iabonds[(is, ib)] + da.nadd;
                    m.tbond[(ia, ib)] = da.tbonds[(is, ib)];
                }
            }
        }
        st.mol.nal[im] = ia;
        st.mol.ngrot[im] = da.ngrots;
        if st.mol.ngrot[im] > 0 {
            st.mol.ngrrot += st.mol.ngrot[im];
            let n = st.mol.ngrrot as usize;
            ensure_dih(&mut st, n);
            let m = &mut st.mol;
            for ih in 1..=m.ngrot[im] {
                ig += 1;
                m.mrot[ig] = im;
                m.rotgrp[ig] = da.rotgrps[ih].clone();
                for j in 1..=2 {
                    m.labaxi[(ig, j)] = da.labaxis[(ih, j)] + da.nadd;
                }
                m.ndihedralang[ig] = da.ndihedralangs[ih];
                // ndihedralangs is always 0 here (see above)
            }
        }
        st.mol.natom = ia;
        st.mol.nmol = im;

        st.mol.dmax = dmaxsave * (st.mol.natomr as f32 / natomrsave as f32).powf(1.0f32 / 3.0f32);

        let mut ifind = 0;
        let mut iav = 0;
        let natom = st.mol.natom;
        let mut iad = natom;
        let mut iadr: i32 = natom;
        let mut iavr: i32 = natomsave;
        let mut iadrr: i32 = natom;
        let mut iavrr: i32 = natomsave;
        let mut dihedralangfixed: f32 = 0.0;
        let mut imv: i32;
        let mut imd = 0;
        let (mut xv, mut yv, mut zv) = (0f32, 0f32, 0f32);

        while iav < natomsave && ifind == 0 {
            if iad >= natom {
                iad = natomsave;
                let r = search_unsaturated(&st, &mut iav, natomsave, &atomv, &atomvr);
                nvalencev = r.0;
                xv = r.1;
                yv = r.2;
                zv = r.3;
            }
            if nvalencev > 0 {
                imv = st.mol.mref[iav];
                imv_set = Some(imv);
                let (nvalenced, xd, yd, zd) = search_unsaturated(&st, &mut iad, natom, &atomd, &atomdr);
                if nvalenced > 0 {
                    imd = st.mol.mref[iad];
                    {
                        let m = &mut st.mol;
                        m.xadd[imd] = da.xadds + xv - xd;
                        m.yadd[imd] = da.yadds + yv - yd;
                        m.zadd[imd] = da.zadds + zv - zd;
                    }
                    iadr = 0;
                    let mut tbondrot: f32 = 1000.0;
                    let mut ib = 0;
                    {
                        let m = &st.mol;
                        while ib < m.nbond[iad] {
                            ib += 1;
                            let tbondtest = m.tbond[(iad, ib)].abs();
                            if m.alabel[m.iabond[(iad, ib)]] == atomdr && tbondtest < tbondrot {
                                iadr = m.iabond[(iad, ib)];
                                tbondrot = tbondtest;
                            }
                        }
                        iavr = 0;
                        ib = 0;
                        while ib < m.nbond[iav] && iavr == 0 {
                            ib += 1;
                            let ivref = m.iabond[(iav, ib)];
                            if m.alabel[ivref] == atomvr && m.tbond[(iav, ib)].abs() == tbondrot {
                                iavr = ivref;
                            }
                        }
                    }
                    if iadr > 0 && iavr > 0 {
                        rotatem(&mut st, imd, iad, iadr, imv, iav, iavr);
                        let r = random_rotatem(&mut st, &atomdrr, iad, iadr, &atomvrr, iav, iavr, dminrr, dmaxrr);
                        iadrr = r.0;
                        iavrr = r.1;
                        dihedralangfixed = r.2;
                        let mut ioverlap = 0;
                        if nooverlaps {
                            ioverlap = test_overlap(&st, imd, 0, natomsave, iav, iavr, iavrr, iad, iadr, iadrr, da.overlapfactor);
                        }
                        if ioverlap == 0 {
                            ifind = 1;
                        } else {
                            ifind = 0;
                            let m = &mut st.mol;
                            m.xadd[imd] = da.xadds;
                            m.yadd[imd] = da.yadds;
                            m.zadd[imd] = da.zadds;
                            let mut ia = m.naf[imd] - 1;
                            for is in 1..=m.nato[imd] {
                                ia += 1;
                                m.xa[ia] = da.xas[is];
                                m.ya[ia] = da.yas[is];
                                m.za[ia] = da.zas[is];
                            }
                        }
                    }
                }
            }
        }
        st.log(&format!("\ndockato> ifind = {ifind:>5}"));
        if ifind == 0 {
            let m = &mut st.mol;
            m.natom = natomsave;
            m.natomr = natomrsave;
            m.nmol = nmolsave;
            m.ngrrot = ngrrotsave;
            m.dmax = dmaxsave;
        } else {
            let fixeddihedral = dminrr == dmaxrr;
            bond_molecule(&mut st, &da, iav, iavr, iavrr, iad, iadr, iadrr, fixeddihedral, dihedralangfixed)?;
            let _ = imd;
            imoladd += 1;
            st.mol.dmaxmol = st.mol.sizefactor * st.mol.dmax;
            iav -= 1;
            let _ = iav;
        }
    }

    // ---- rotational groups: lists of atoms to be rotated
    let nbondp = get_direct_bonds(&mut st);
    st.mol.nrottot = 0;
    let ngrrot = st.mol.ngrrot;
    for ig in 1..=ngrrot {
        let ir1 = st.mol.labaxi[(ig, 1)];
        let ir2 = st.mol.labaxi[(ig, 2)];
        let labrott1 = get_rotgrp(&st, &nbondp, ir1, ir2);
        let labrott2 = get_rotgrp(&st, &nbondp, ir2, ir1);
        let m = &mut st.mol;
        let (a1, a2, lab) = if labrott1.len() < labrott2.len() { (ir1, ir2, labrott1) } else { (ir2, ir1, labrott2) };
        m.labaxi[(ig, 1)] = a1;
        m.labaxi[(ig, 2)] = a2;
        m.nrot[ig] = lab.len() as i32;
        m.rotf[ig] = m.nrottot + 1;
        m.nrottot += m.nrot[ig];
        m.rotl[ig] = m.nrottot;
        let n = m.nrottot as usize;
        m.ensure_labrot(n);
        for i in 1..=m.nrot[ig] {
            let labref = m.rotf[ig] + i - 1;
            m.labrot[labref] = lab[(i - 1) as usize];
        }
    }

    st.log(&format!("\ndockato> Added {imoladd:>5} molecules"));
    if imoladd < ndock {
        st.log("\ndockato> (no more valid sites found for docking.)");
    }

    // ---- output file name
    let mut sometext;
    if ilenoutputfile > 0 {
        sometext = outputfile.clone();
    } else {
        loop {
            match ctx.stdin.read_line() {
                Ok(l) => sometext = l,
                Err(_) => return eof("output file name"),
            }
            if parse(&sometext, 5).nwords() > 0 {
                break;
            }
            st.log("\ndockato> No .ato file specified - try again!\n");
        }
    }
    let w = parse(&sometext, 5);
    let mut readfile = w.span(&sometext, 1, 1);
    let ilendot = findex(&readfile, ".") - 1;
    if ilendot > 0 {
        readfile = readfile[..ilendot as usize].to_string();
    }
    let fnameato = format!("{}{}.ato", st.sym.directory_name, readfile);
    let imv = match imv_set {
        Some(v) => v,
        None => return crate::error::unsupported("dockato: no docking site found at all (molecule index `imv` undefined in the Fortran)"),
    };
    let mt = st.mol.mtype[imv];
    if mt < 1 {
        // The Fortran then writes moltemplate(0), i.e. outside the array
        // (it blanks the moltype labels and whatever precedes /clusterc/).
        return crate::error::unsupported(&format!(
            "dockato: last docked-onto molecule {imv} has no molecule type (Fortran writes moltemplate(0))"
        ));
    }
    st.mol.moltemplate[mt] = readfile;
    st.moloutc(&fnameato)
}

/// `subroutine rotatem(imd,iad,iadr,imv,iav,iavr,ierr)`.
fn rotatem(st: &mut State, imd: i32, iad: i32, iadr: i32, imv: i32, iav: i32, iavr: i32) -> i32 {
    let m = &mut st.mol;
    let xadiff = m.xadd[imv] - m.xadd[imd] + m.xa[iav] - m.xa[iad];
    let yadiff = m.yadd[imv] - m.yadd[imd] + m.ya[iav] - m.ya[iad];
    let zadiff = m.zadd[imv] - m.zadd[imd] + m.za[iav] - m.za[iad];
    let dsq1 = xadiff * xadiff + yadiff * yadiff + zadiff * zadiff;
    let qradiussum = 6.0f32 * (m.qradius[m.itype[iad]] + m.qradius[m.itype[iav]]);
    if dsq1.sqrt() > qradiussum {
        return 1;
    } else if iavr == 0 || iadr == 0 {
        return 1;
    }
    let mut vvector = [
        m.xadd[imv] + m.xa[iavr] - m.xadd[imd] - m.xa[iad],
        m.yadd[imv] + m.ya[iavr] - m.yadd[imd] - m.ya[iad],
        m.zadd[imv] + m.za[iavr] - m.zadd[imd] - m.za[iad],
    ];
    let mut rvector = [m.xa[iadr] - m.xa[iad], m.ya[iadr] - m.ya[iad], m.za[iadr] - m.za[iad]];
    let vnorm = scapro(&vvector, &vvector).sqrt();
    let rnorm = scapro(&rvector, &rvector).sqrt();
    for i in 0..3 {
        vvector[i] /= vnorm;
        rvector[i] /= rnorm;
    }
    let mut cphi = scapro(&vvector, &rvector);
    if cphi.abs() > 1.0 {
        cphi = sign(1.0, cphi);
    }
    let sphi = (1.0f32 - cphi * cphi).sqrt();
    if sphi > 0.001 {
        let mut pvector = vecpro(&rvector, &vvector);
        let rab = scapro(&pvector, &pvector).sqrt();
        pvector[0] /= rab;
        pvector[1] /= rab;
        pvector[2] /= rab;
        for ia in m.naf[imd]..=m.nal[imd] {
            if ia != iad {
                let xvector = [m.xa[ia] - m.xa[iad], m.ya[ia] - m.ya[iad], m.za[ia] - m.za[iad]];
                let frac = scapro(&pvector, &xvector);
                let xpvector = [frac * pvector[0], frac * pvector[1], frac * pvector[2]];
                let xrvector = [xvector[0] - xpvector[0], xvector[1] - xpvector[1], xvector[2] - xpvector[2]];
                let xrpvector = vecpro(&pvector, &xrvector);
                m.xa[ia] = m.xa[iad] + xpvector[0];
                m.ya[ia] = m.ya[iad] + xpvector[1];
                m.za[ia] = m.za[iad] + xpvector[2];
                m.xa[ia] = m.xa[ia] + xrvector[0] * cphi + xrpvector[0] * sphi;
                m.ya[ia] = m.ya[ia] + xrvector[1] * cphi + xrpvector[1] * sphi;
                m.za[ia] = m.za[ia] + xrvector[2] * cphi + xrpvector[2] * sphi;
            }
        }
    } else if cphi < 1.0 {
        for ia in m.naf[imd]..=m.nal[imd] {
            m.xa[ia] = m.xa[iad] - (m.xa[ia] - m.xa[iad]);
            m.ya[ia] = m.ya[iad] - (m.ya[ia] - m.ya[iad]);
            m.za[ia] = m.za[iad] - (m.za[ia] - m.za[iad]);
        }
    }
    0
}

/// `subroutine random_rotatem(...)`: returns `(iadrr, iavrr, dihedralangfixed)`.
#[allow(clippy::too_many_arguments)]
fn random_rotatem(
    st: &mut State,
    atomdrr: &str,
    iad: i32,
    iadr: i32,
    atomvrr: &str,
    iav: i32,
    iavr: i32,
    dminrr: f32,
    dmaxrr: f32,
) -> (i32, i32, f32) {
    let pi = std::f32::consts::PI;
    let r1;
    {
        let rr = st.rng.ran1();
        r1 = rr;
    }
    let m = &mut st.mol;
    let mut rvector = [m.xa[iadr] - m.xa[iad], m.ya[iadr] - m.ya[iad], m.za[iadr] - m.za[iad]];
    let rnorm = scapro(&rvector, &rvector).sqrt();
    for v in rvector.iter_mut() {
        *v /= rnorm;
    }
    let arg = 2.0f32 * pi * (dminrr + r1 * (dmaxrr - dminrr)) / 360.0f32;
    let dihedralangfixed = arg * 180.0f32 / pi;
    let sphistep = arg.sin();
    let cphistep = arg.cos();

    let imd = m.mref[iad];
    let mut iat = m.naf[imd] - 1;
    let mut iadrr = 0;
    let mut dsqmin: f32 = 1.0e30;
    let mut rdpvector = [0f32; 3];
    let mut rvpvector = [0f32; 3];
    while iat < m.nal[imd] && atomdrr != "0" {
        iat += 1;
        if iat != iad && iat != iadr && m.alabel[iat] == atomdrr {
            let dvector = [m.xa[iat] - m.xa[iad], m.ya[iat] - m.ya[iad], m.za[iat] - m.za[iad]];
            let dsq = scapro(&dvector, &dvector);
            if dsq < dsqmin {
                rdpvector = vecpro(&rvector, &dvector);
                let rdpnorm = scapro(&rdpvector, &rdpvector).sqrt();
                for v in rdpvector.iter_mut() {
                    *v /= rdpnorm;
                }
                let cphi = scapro(&rdpvector, &rvector);
                if cphi > -1.0 && cphi < 1.0 {
                    dsqmin = dsq;
                    iadrr = iat;
                }
            }
        }
    }
    let imv = m.mref[iav];
    iat = m.naf[imv] - 1;
    let mut iavrr = 0;
    dsqmin = 1.0e30;
    while iat < m.nal[imv] && atomvrr != "0" {
        iat += 1;
        if iat != iav && iat != iavr && m.alabel[iat] == atomvrr {
            let vvector = [m.xa[iat] - m.xa[iav], m.ya[iat] - m.ya[iav], m.za[iat] - m.za[iav]];
            let dsq = scapro(&vvector, &vvector);
            if dsq < dsqmin {
                rvpvector = vecpro(&rvector, &vvector);
                let rvpnorm = scapro(&rvpvector, &rvpvector).sqrt();
                for v in rvpvector.iter_mut() {
                    *v /= rvpnorm;
                }
                let cphi = scapro(&rvpvector, &rvector);
                if cphi > -1.0 && cphi < 1.0 {
                    dsqmin = dsq;
                    iavrr = iat;
                }
            }
        }
    }
    let cphinew: f32;
    let mut sphinew: f32 = 0.0;
    if iadrr > 0 && iavrr > 0 {
        let mut cphiold = scapro(&rdpvector, &rvpvector);
        if cphiold.abs() >= 1.0 {
            cphiold = sign(1.0, cphiold);
        }
        let drot = vecpro(&rdpvector, &rvpvector);
        let sphiold = scapro(&drot, &rvector);
        cphinew = cphiold * cphistep + sphiold * sphistep;
        sphinew = sphiold * cphistep - cphiold * sphistep;
    } else if iadrr > 0 {
        sphinew = sphistep;
        cphinew = cphistep;
    } else {
        cphinew = 2.0;
    }
    if (-1.0..=1.0).contains(&cphinew) {
        let mut iat = m.naf[imd] - 1;
        while iat < m.nal[imd] {
            iat += 1;
            if iat != iad && iat != iadr {
                let dvector = [m.xa[iat] - m.xa[iad], m.ya[iat] - m.ya[iad], m.za[iat] - m.za[iad]];
                let frac = scapro(&dvector, &rvector);
                let rdp = vecpro(&rvector, &dvector);
                let mut drot = [0f32; 3];
                for i in 0..3 {
                    let dpar = frac * rvector[i];
                    let dperp = dvector[i] - dpar;
                    drot[i] = dpar + dperp * cphinew + rdp[i] * sphinew;
                }
                m.xa[iat] = m.xa[iad] + drot[0];
                m.ya[iat] = m.ya[iad] + drot[1];
                m.za[iat] = m.za[iad] + drot[2];
            }
        }
    }
    // (the Fortran then only prints a check of the new angle)
    (iadrr, iavrr, dihedralangfixed)
}

/// `subroutine search_unsaturated(iav,iavlast,nvalence,atomv,atomvr,xv,yv,zv)`.
fn search_unsaturated(st: &State, iav: &mut i32, iavlast: i32, atomv: &str, atomvr: &str) -> (i32, f32, f32, f32) {
    let m = &st.mol;
    let mut nvalence = 0;
    let mut bondsearch: f32 = 1000.0;
    let mut ia = *iav;
    while ia < iavlast {
        ia += 1;
        if m.alabel[ia] == atomv {
            let (nvalencefound, bondfound) = find_valence(st, ia, atomvr);
            if nvalencefound > 0 && bondfound < bondsearch {
                *iav = ia;
                bondsearch = bondfound;
                nvalence = nvalencefound;
            }
        }
    }
    if nvalence == 0 {
        *iav = ia;
    }
    let imv = m.mref[*iav];
    let xv = m.xadd[imv] + m.xa[*iav];
    let yv = m.yadd[imv] + m.ya[*iav];
    let zv = m.zadd[imv] + m.za[*iav];
    (nvalence, xv, yv, zv)
}

/// `subroutine find_valence(iav,nvalence,atomvr,bondsearch)`.
fn find_valence(st: &State, iav: i32, atomvr: &str) -> (i32, f32) {
    let m = &st.mol;
    let mut ib = 0;
    let mut bondsearch: f32 = 1000.0;
    while ib < m.nbond[iav] {
        ib += 1;
        if m.alabel[m.iabond[(iav, ib)]] == atomvr && m.tbond[(iav, ib)].abs() < bondsearch {
            bondsearch = m.tbond[(iav, ib)].abs();
        }
    }
    ib = 0;
    let mut nvalence = 0;
    while ib < m.nbond[iav] {
        ib += 1;
        if m.alabel[m.iabond[(iav, ib)]] == atomvr && m.tbond[(iav, ib)].abs() == bondsearch {
            nvalence += 1;
        }
    }
    (nvalence, bondsearch)
}

/// `subroutine test_overlap(...)`.
#[allow(clippy::too_many_arguments)]
fn test_overlap(
    st: &State,
    imd: i32,
    natomfirst: i32,
    natomsave: i32,
    iav: i32,
    iavr: i32,
    iavrr: i32,
    iad: i32,
    iadr: i32,
    iadrr: i32,
    overlapfactor: f32,
) -> i32 {
    let m = &st.mol;
    let mut ioverlap = 0;
    let mut inew = m.naf[imd] - 1;
    while inew < m.nal[imd] && ioverlap == 0 {
        inew += 1;
        let ityped = m.itype[inew];
        let xt = m.xadd[imd] + m.xa[inew];
        let yt = m.yadd[imd] + m.ya[inew];
        let zt = m.zadd[imd] + m.za[inew];
        let iadneighb = find_neighbour(st, iad, inew);
        let iadrneighb = find_neighbour(st, iadr, inew);
        let _iadrrneighb = find_neighbour(st, iadrr, inew);
        let mut iat = natomfirst;
        while iat < natomsave && ioverlap == 0 {
            iat += 1;
            let iavneighb = find_neighbour(st, iav, iat);
            let iavrneighb = find_neighbour(st, iavr, iat);
            let _iavrrneighb = find_neighbour(st, iavrr, iat);
            let imt = m.mref[iat];
            let itt = m.itype[iat];
            let epsilont = (m.epscor[itt] * m.epscor[ityped]).sqrt();
            if epsilont > 0.0 && ((iavneighb == 0 || iadrneighb == 0) || (iavrneighb == 0 || iadneighb == 0)) {
                let sigmat = 0.5f32 * overlapfactor * (m.sigcor[itt] + m.sigcor[ityped]) / m.dmax;
                let mut xx = (m.xadd[imt] + m.xa[iat] - xt) / m.dmax;
                let mut yy = (m.yadd[imt] + m.ya[iat] - yt) / m.dmax;
                let mut zz = (m.zadd[imt] + m.za[iat] - zt) / m.dmax;
                minimagec(&mut xx, &mut yy, &mut zz);
                let dtest = (xx * xx + yy * yy + zz * zz).sqrt();
                if dtest < sigmat {
                    ioverlap = 1;
                    st.log(&format!("\ntest_overlap> {iat} {imt} {itt} {inew} {imd} {ityped} {dtest} {sigmat}"));
                }
            }
        }
    }
    ioverlap
}

/// `subroutine find_neighbour(iad,inew,iadneighb)`.
fn find_neighbour(st: &State, iad: i32, inew: i32) -> i32 {
    let m = &st.mol;
    let mut iadneighb = 0;
    if iad == inew {
        iadneighb = 1;
    }
    if iad > 0 {
        let mut ib = 0;
        while ib < m.nbond[iad] && iadneighb == 0 {
            ib += 1;
            if m.iabond[(iad, ib)] == inew {
                iadneighb = 1;
            }
        }
    }
    iadneighb
}

/// Copy rotational group `src` to `dst` (rotgrp, labaxi, dihedrals).
fn copy_group(st: &mut State, dst: i32, src: i32) {
    let m = &mut st.mol;
    m.rotgrp[dst] = m.rotgrp[src].clone();
    for j in 1..=2 {
        m.labaxi[(dst, j)] = m.labaxi[(src, j)];
    }
    m.ndihedralang[dst] = m.ndihedralang[src];
    for k in 1..=m.ndihedralang[dst] {
        m.dihedralang[(k, dst)] = m.dihedralang[(k, src)];
        for j in 1..=2 {
            m.labrotdih[(j, k, dst)] = m.labrotdih[(j, k, src)];
        }
    }
}

/// `subroutine bond_molecule(iav,iavr,iavrr,iad,iadr,iadrr,fixeddihedral,dihedralangfixed)`.
#[allow(clippy::too_many_arguments)]
fn bond_molecule(
    st: &mut State,
    da: &Da,
    iav: i32,
    iavr: i32,
    mut iavrr: i32,
    iad: i32,
    iadr: i32,
    mut iadrr: i32,
    fixeddihedral: bool,
    dihedralangfixed: f32,
) -> Result<()> {
    let imd = st.mol.mref[iad];
    let imv = st.mol.mref[iav];
    let mut ig_opt: Option<i32> = None;
    st.log(&format!(" bond_molecule> 2 {iav} {iavr} {iavrr} {iad} {iadr} {iadrr}"));
    st.log(&format!(" bond_molecule> 3 {} {} {} {} {} {}", imv, imd, st.mol.ngrot[imv], da.natos, da.ngrots, st.mol.ngrrot));
    let n = st.mol.ngrrot.max(0) as usize + 1;
    ensure_dih(st, n);
    if st.mol.mref.len() < n {
        st.mol.mref.grow(n);
    }

    if da.ngrots > 0 {
        let mut iglast = 0;
        let mut im = 0;
        while im < imv {
            im += 1;
            iglast += st.mol.ngrot[im];
        }
        let mut ig = st.mol.ngrrot;
        let mut ig1 = ig - da.ngrots;
        while ig1 > iglast {
            // NB: the Fortran copies mref (the atom -> molecule map) here, not mrot.
            st.mol.mref[ig] = st.mol.mref[ig1];
            copy_group(st, ig, ig1);
            ig -= 1;
            ig1 -= 1;
        }
        st.mol.ngrot[imv] += da.ngrots;
        ig = iglast;
        let mut ih = 0;
        let m = &mut st.mol;
        while ih < da.ngrots {
            ig += 1;
            ih += 1;
            m.mref[ig] = imv;
            m.rotgrp[ig] = da.rotgrps[ih].clone();
            for j in 1..=2 {
                m.labaxi[(ig, j)] = da.labaxis[(ih, j)] + da.nadd;
            }
            m.ndihedralang[ig] = da.ndihedralangs[ih];
        }
        ig_opt = Some(ig);
    }

    if st.mol.nato[imv] > 2 {
        st.mol.ngrot[imv] += 1;
        st.mol.ngrrot += 1;
        let n = st.mol.ngrrot as usize;
        ensure_dih(st, n);
        let mut iglast = 0;
        let mut im = 0;
        while im < imv {
            im += 1;
            iglast += st.mol.ngrot[im];
        }
        let mut ig = st.mol.ngrrot;
        while ig > iglast {
            let ig1 = ig - 1;
            copy_group(st, ig, ig1);
            ig = ig1;
        }
        let m = &mut st.mol;
        m.rotgrp[ig] = "ROT".into();
        m.labaxi[(ig, 1)] = iadr;
        m.labaxi[(ig, 2)] = iav;
        m.ndihedralang[ig] = 0;
        if m.stepri <= 0.0 {
            m.stepri = 1.0;
        }
        ig_opt = Some(ig);
    }

    let mut noverlap: usize = 2;
    let mut ikeep = vec![0i32; NEAR + 1];
    let mut ilose = vec![0i32; NEAR + 1];
    ikeep[1] = iav;
    ilose[1] = iad;
    ikeep[2] = iadr;
    ilose[2] = iavr;
    let mut ioverlaprr = 0;
    if iadrr > 0 && iavrr > 0 {
        let m = &mut st.mol;
        let mut iadt = m.naf[imd] - 1;
        while iadt < m.nal[imd] {
            iadt += 1;
            if iadt != iad && iadt != iadr {
                let x1 = m.xadd[imd] + m.xa[iadt];
                let y1 = m.yadd[imd] + m.ya[iadt];
                let z1 = m.zadd[imd] + m.za[iadt];
                let mut iavt = m.naf[imv] - 1;
                while iavt < m.nal[imv] {
                    iavt += 1;
                    if iavt != iav && iavt != iavr {
                        let x2 = m.xadd[imv] + m.xa[iavt];
                        let y2 = m.yadd[imv] + m.ya[iavt];
                        let z2 = m.zadd[imv] + m.za[iavt];
                        let mut xdiff = (x2 - x1) / m.dmax;
                        let mut ydiff = (y2 - y1) / m.dmax;
                        let mut zdiff = (z2 - z1) / m.dmax;
                        minimagec(&mut xdiff, &mut ydiff, &mut zdiff);
                        let dist = (xdiff * xdiff + ydiff * ydiff + zdiff * zdiff).sqrt() * m.dmax;
                        let dmin = m.qradius[m.itype[iadt]] + m.qradius[m.itype[iavt]];
                        if dist <= dmin {
                            noverlap += 1;
                            if noverlap > NEAR {
                                return Err(EpsrError::Input("bond_molecule: too many overlapping atoms".into()));
                            }
                            ikeep[noverlap] = iavt;
                            ilose[noverlap] = iadt;
                            if iadt == iadrr && iavt == iavrr {
                                ioverlaprr = 1;
                            }
                        }
                    }
                }
            }
        }
    }

    {
        let m = &mut st.mol;
        for ia in m.naf[imd]..=m.nal[imd] {
            m.xa[ia] = m.xa[ia] + m.xadd[imd] - m.xadd[imv];
            m.ya[ia] = m.ya[ia] + m.yadd[imd] - m.yadd[imv];
            m.za[ia] = m.za[ia] + m.zadd[imd] - m.zadd[imv];
            m.mref[ia] = imv;
        }
    }
    for io in 1..=noverlap {
        replace_references(st, ilose[io], ikeep[io]);
        if iavrr == ilose[io] {
            iavrr = ikeep[io];
        }
        if iadrr == ilose[io] {
            iadrr = ikeep[io];
        }
    }
    for io in 1..=noverlap {
        add_neighbours(st, ilose[io], ikeep[io]);
    }

    if ioverlaprr == 0 && fixeddihedral {
        let m = &st.mol;
        let mut ib = 0;
        let mut ifound = 0;
        while ib < nbond_at(st, iavrr) && ifound == 0 {
            ib += 1;
            if m.iabond[(iavrr, ib)] == iadrr {
                ifound = iadrr;
            }
        }
        ib = 0;
        while ib < nbond_at(st, iadrr) && ifound == 0 {
            ib += 1;
            if m.iabond[(iadrr, ib)] == iavrr {
                ifound = iavrr;
            }
        }
        if ifound == 0 {
            let ig = match ig_opt {
                Some(v) => v,
                None => return crate::error::unsupported("dockato: fixed dihedral with no rotational group (uninitialised group index in the Fortran)"),
            };
            let m = &mut st.mol;
            m.rotgrp[ig] = "DIH".into();
            m.ndihedralang[ig] = 1;
            m.dihedralang[(1, ig)] = dihedralangfixed;
            m.labrotdih[(1, 1, ig)] = iadrr;
            m.labrotdih[(2, 1, ig)] = iavrr;
        }
    }

    {
        let m = &mut st.mol;
        m.nato[imv] += m.nato[imd];
        m.nal[imv] += m.nato[imd];
        m.amasst[imv] += m.amasst[imd];
    }

    for io in 1..=noverlap {
        let il = ilose[io];
        remove_atoms(st, il, imv);
        let mut jo = 0;
        while jo < noverlap {
            jo += 1;
            if ikeep[jo] > ilose[io] {
                ikeep[jo] -= 1;
            }
            if ilose[jo] > ilose[io] {
                ilose[jo] -= 1;
            }
        }
    }

    // remove duplicated rotation axes
    {
        let m = &mut st.mol;
        let mut igb = 0;
        while igb < m.ngrrot {
            igb += 1;
            let mut ig = igb;
            while ig < m.ngrrot {
                ig += 1;
                if (m.labaxi[(ig, 1)] == m.labaxi[(igb, 1)] && m.labaxi[(ig, 2)] == m.labaxi[(igb, 2)])
                    || (m.labaxi[(ig, 2)] == m.labaxi[(igb, 1)] && m.labaxi[(ig, 1)] == m.labaxi[(igb, 2)])
                {
                    m.ngrrot -= 1;
                    let mut igt = ig - 1;
                    while igt < m.ngrrot {
                        igt += 1;
                        let igt1 = igt + 1;
                        m.rotgrp[igt] = m.rotgrp[igt1].clone();
                        for j in 1..=2 {
                            m.labaxi[(igt, j)] = m.labaxi[(igt1, j)];
                        }
                        m.ndihedralang[igt] = m.ndihedralang[igt1];
                        if m.ndihedralang[igt] > 0 {
                            for k in 1..=m.ndihedralang[igt] {
                                m.dihedralang[(k, igt)] = m.dihedralang[(k, igt1)];
                                for j in 1..=2 {
                                    m.labrotdih[(j, k, igt)] = m.labrotdih[(j, k, igt1)];
                                }
                            }
                        }
                    }
                    ig -= 1;
                }
            }
        }
        for im in 1..=m.nmol {
            m.ngrot[im] = 0;
        }
        let mut ig = 0;
        while ig < m.ngrrot {
            ig += 1;
            let im1 = m.mref[m.labaxi[(ig, 1)]];
            let im2 = m.mref[m.labaxi[(ig, 2)]];
            if im1 == im2 {
                m.ngrot[im1] += 1;
                m.mrot[ig] = im1;
            }
        }
    }

    // remove atom types no longer present
    {
        let m = &mut st.mol;
        let mut itb = 0;
        while itb < m.ntype {
            itb += 1;
            let mut ia = 0;
            let mut itfound = 0;
            while ia < m.natom && itfound == 0 {
                ia += 1;
                if m.itype[ia] == itb {
                    itfound = itb;
                }
            }
            if itfound == 0 {
                m.ntype -= 1;
                let mut it = itb - 1;
                while it < m.ntype {
                    it += 1;
                    let it1 = it + 1;
                    m.type_[it] = m.type_[it1].clone();
                    m.atsymbol[it] = m.atsymbol[it1].clone();
                    m.iexchange[it] = m.iexchange[it1];
                    m.epscor[it] = m.epscor[it1];
                    m.sigcor[it] = m.sigcor[it1];
                    m.amassa[it] = m.amassa[it1];
                    m.charge[it] = m.charge[it1];
                    m.qradius[it] = m.qradius[it1];
                }
                let mut ia = 0;
                while ia < m.natom {
                    ia += 1;
                    if m.itype[ia] > itb {
                        m.itype[ia] -= 1;
                    }
                }
            }
        }
    }

    centreofmass(st, imv);
    st.mol.nmol -= 1;
    st.log(&format!(" bond_molecule> 5 {} {} {} {} {}", st.mol.nmol, st.mol.natom, st.mol.natomr, st.mol.ntype, st.mol.ngrrot));
    Ok(())
}

/// `subroutine replace_references(iad,iav)`.
fn replace_references(st: &mut State, iad: i32, iav: i32) {
    let m = &mut st.mol;
    let imd = m.mref[iad];
    let imv = m.mref[iav];
    m.itype[iad] = m.itype[iav];
    for ia in 1..=m.natom {
        let mut ib = 0;
        while ib < m.nbond[ia] {
            ib += 1;
            if m.iabond[(ia, ib)] == iad {
                m.iabond[(ia, ib)] = iav;
            }
        }
    }
    let mut ig = 0;
    while ig < m.ngrrot {
        ig += 1;
        if m.mrot[ig] == imd {
            m.mrot[ig] = imv;
        }
        for j in 1..=2 {
            if m.labaxi[(ig, j)] == iad {
                m.labaxi[(ig, j)] = iav;
            }
        }
        if m.ndihedralang[ig] > 0 {
            for k in 1..=m.ndihedralang[ig] {
                for j in 1..=2 {
                    if m.labrotdih[(j, k, ig)] == iad {
                        m.labrotdih[(j, k, ig)] = iav;
                    }
                }
            }
        }
    }
}

/// `subroutine add_neighbours(iad,iav)`.
fn add_neighbours(st: &mut State, iad: i32, iav: i32) {
    let mut ib = 0;
    while ib < st.mol.nbond[iad] {
        ib += 1;
        let iaref = st.mol.iabond[(iad, ib)];
        let mut ibc = 0;
        let mut ifound = 0;
        while ibc < st.mol.nbond[iav] && ifound == 0 {
            ibc += 1;
            if st.mol.iabond[(iav, ibc)] == iaref {
                ifound = 1;
            }
        }
        if ifound == 0 {
            let nb = st.mol.nbond[iav] + 1;
            st.mol.ensure_bonds(nb as usize);
            let m = &mut st.mol;
            m.nbond[iav] = nb;
            m.iabond[(iav, nb)] = iaref;
            m.tbond[(iav, nb)] = m.tbond[(iad, ib)];
        }
        ibc = 0;
        ifound = 0;
        while ibc < st.mol.nbond[iaref] && ifound == 0 {
            ibc += 1;
            if st.mol.iabond[(iaref, ibc)] == iav {
                ifound = 1;
            }
        }
        if ifound == 0 {
            let nb = st.mol.nbond[iaref] + 1;
            st.mol.ensure_bonds(nb as usize);
            let m = &mut st.mol;
            m.nbond[iaref] = nb;
            m.iabond[(iaref, nb)] = iav;
            m.tbond[(iaref, nb)] = m.tbond[(iad, ib)];
        }
    }
}

/// `subroutine remove_atoms(iad,imv)`.
fn remove_atoms(st: &mut State, iad: i32, imv: i32) {
    {
        let m = &mut st.mol;
        m.nato[imv] -= 1;
        m.nal[imv] -= 1;
        if m.qtype[iad] != 1 {
            m.natomr -= 1;
        }
        let mut ia = iad;
        while ia < m.natom {
            let ia1 = ia + 1;
            m.xa[ia] = m.xa[ia1];
            m.ya[ia] = m.ya[ia1];
            m.za[ia] = m.za[ia1];
            m.alabel[ia] = m.alabel[ia1].clone();
            m.itype[ia] = m.itype[ia1];
            m.qtype[ia] = m.qtype[ia1];
            m.nbond[ia] = m.nbond[ia1];
            m.mref[ia] = m.mref[ia1];
            let mut ib = 0;
            while ib < m.nbond[ia] {
                ib += 1;
                m.iabond[(ia, ib)] = m.iabond[(ia1, ib)];
                m.tbond[(ia, ib)] = m.tbond[(ia1, ib)];
            }
            ia = ia1;
        }
        for ia in 1..=m.natom {
            let mut ib = 0;
            while ib < m.nbond[ia] {
                ib += 1;
                if m.iabond[(ia, ib)] > iad {
                    m.iabond[(ia, ib)] -= 1;
                }
            }
        }
        let mut ig = 0;
        while ig < m.ngrrot {
            ig += 1;
            for j in 1..=2 {
                if m.labaxi[(ig, j)] > iad {
                    m.labaxi[(ig, j)] -= 1;
                }
            }
            if m.ndihedralang[ig] > 0 {
                for k in 1..=m.ndihedralang[ig] {
                    for j in 1..=2 {
                        if m.labrotdih[(j, k, ig)] > iad {
                            m.labrotdih[(j, k, ig)] -= 1;
                        }
                    }
                }
            }
        }
        m.natom -= 1;
    }
    clean_atom_types(st);
}

/// `subroutine clean_atom_types()`.
fn clean_atom_types(st: &mut State) {
    let m = &mut st.mol;
    let mut iatomtype = 0;
    while iatomtype < m.ntype {
        iatomtype += 1;
        let mut ia = 1;
        let mut ifoundatomtype = m.type_[iatomtype] == m.alabel[ia];
        while ia < m.natom && !ifoundatomtype {
            ia += 1;
            ifoundatomtype = m.type_[iatomtype] == m.alabel[ia];
        }
        if !ifoundatomtype {
            if iatomtype < m.ntype {
                for it in iatomtype..=m.ntype - 1 {
                    let it1 = it + 1;
                    m.type_[it] = m.type_[it1].clone();
                    m.qtypea[it] = m.qtypea[it1];
                    m.iatomicnotype[it] = m.iatomicnotype[it1];
                    m.atsymbol[it] = m.atsymbol[it1].clone();
                    m.iexchange[it] = m.iexchange[it1];
                    m.epscor[it] = m.epscor[it1];
                    m.sigcor[it] = m.sigcor[it1];
                    m.amassa[it] = m.amassa[it1];
                    m.charge[it] = m.charge[it1];
                    m.qradius[it] = m.qradius[it1];
                    m.natomsphere[it] = m.natomsphere[it1];
                    m.rhosphere[it] = m.rhosphere[it1];
                    m.radminsphere[it] = m.radminsphere[it1];
                    m.radmaxsphere[it] = m.radmaxsphere[it1];
                }
                for ia in 1..=m.natom {
                    if m.itype[ia] > iatomtype {
                        m.itype[ia] -= 1;
                    }
                }
                m.ntype -= 1;
            } else {
                m.ntype -= 1;
            }
        }
    }
}

/// `subroutine centreofmass(im)` (ato_manipulation_rev4.f): recentre molecule
/// `im` on its centre of mass and move the molecule centre accordingly.
fn centreofmass(st: &mut State, im: i32) {
    let m = &mut st.mol;
    let (mut x1, mut y1, mut z1) = (0f32, 0f32, 0f32);
    let mut summass: f32 = 0.0;
    for ia in m.naf[im]..=m.nal[im] {
        if m.qtype[ia] == 0 {
            let amass = m.amassa[m.itype[ia]];
            x1 += m.xa[ia] * amass;
            y1 += m.ya[ia] * amass;
            z1 += m.za[ia] * amass;
            summass += amass;
        }
    }
    if summass > 0.0 {
        x1 /= summass;
        y1 /= summass;
        z1 /= summass;
    }
    for ia in m.naf[im]..=m.nal[im] {
        m.xa[ia] -= x1;
        m.ya[ia] -= y1;
        m.za[ia] -= z1;
    }
    m.xadd[im] += x1;
    m.yadd[im] += y1;
    m.zadd[im] += z1;
}

/// `subroutine get_direct_bonds()`: sort each atom's bonds by increasing
/// |length| and count the direct (positive) bonds; returns `nbondp`.
fn get_direct_bonds(st: &mut State) -> F1<i32> {
    let m = &mut st.mol;
    for ia in 1..=m.natom {
        let mut ib1 = 0;
        let mut ib2 = ib1 + 1;
        while ib2 < m.nbond[ia] {
            ib1 += 1;
            ib2 = ib1 + 1;
            if m.tbond[(ia, ib2)].abs() < m.tbond[(ia, ib1)].abs() {
                let iabondsave = m.iabond[(ia, ib2)];
                let tbondsave = m.tbond[(ia, ib2)];
                m.iabond[(ia, ib2)] = m.iabond[(ia, ib1)];
                m.tbond[(ia, ib2)] = m.tbond[(ia, ib1)];
                m.iabond[(ia, ib1)] = iabondsave;
                m.tbond[(ia, ib1)] = tbondsave;
                ib1 = 0;
                ib2 = ib1 + 1;
            }
        }
    }
    // nbondp is a COMMON array (dimension matom): entries beyond natom stay 0.
    let mut nbondp: F1<i32> = F1::new((m.natom.max(0) as usize).max(m.nbond.len()) + 1);
    for ia in 1..=m.natom {
        nbondp[ia] = 0;
        let mut irefb = 0;
        while irefb < m.nbond[ia] {
            irefb += 1;
            if m.tbond[(ia, irefb)] > 0.0 {
                nbondp[ia] += 1;
            }
        }
    }
    nbondp
}

/// `subroutine get_rotgrp(ir1,ir2,labrott,nrott)`: atoms rotated about the
/// axis ir1-ir2 (on the ir2 side... as found by following direct bonds from ir1).
fn get_rotgrp(st: &State, nbondp: &F1<i32>, ir1: i32, ir2: i32) -> Vec<i32> {
    let m = &st.mol;
    let nbp = |i: i32| -> i32 { if i >= 1 && (i as usize) <= nbondp.len() { nbondp[i] } else { 0 } };
    let mut labrott: Vec<i32> = Vec::new();
    let mut iselect = ir1;
    let mut ir2test = false;
    let mut nrotnew = 0usize;
    while nrotnew <= labrott.len() {
        let mut iref = 0;
        while iref < nbp(iselect) {
            iref += 1;
            let iaref = m.iabond[(iselect, iref)];
            if !ir2test {
                ir2test = iselect != ir1 && iaref == ir2;
            }
            if iaref != ir1 && iaref != ir2 && !labrott.contains(&iaref) {
                labrott.push(iaref);
            }
        }
        nrotnew += 1;
        if nrotnew <= labrott.len() {
            iselect = labrott[nrotnew - 1];
        }
    }
    if ir2test {
        let mut iref = 0;
        while iref < nbp(ir2) {
            iref += 1;
            let iaref = m.iabond[(ir2, iref)];
            if iaref != ir1 && iaref != ir2 && !labrott.contains(&iaref) {
                labrott.push(iaref);
            }
        }
    }
    labrott
}
