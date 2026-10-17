//! `changeato` — translation of `fortran_src/changeato.f` and of the routines
//! it uses from `process_symbols.f`: the interactive keyword loop
//! (`process_symbols` / `process_commands` / `check_symbol` with `nin=5`,
//! `ignore=0`), `write_symbols_changeato`, `change_symbol_changeato`
//! (= the `symbol_setup.eq.'changeato'` branch of `check_symbol`),
//! `change_bond`, `change_label`, `add_bond`, `remove_bond`,
//! `remove_molecule_types`, `remove_molecule`, `cleanup_types`,
//! `save_files`/`force_save_files` + `save_files_changeato`, and
//! `change_boxsize` (`ato_manipulation_rev4.f`, cubic boxes).
//!
//! STATUS: translated (the options listed below), verified byte-identical
//! (EPSR-test-ex26, ex27).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> changeato <file>.ato < answers`
//! Keyword table: `crate::symtab_tools::CHANGEATO`.
//!
//! The program reads the `.ato` file, shows its parameters as keywords and
//! lets the user change them from standard input (`<keyword> <value>`, a bare
//! value for the current keyword, `atom <n>` to select an atom type, `d`/`u`
//! to step, `n` to revert, `save`, `q`, `e`...).  The modified configuration
//! is written back with `moloutc`.
//!
//! Supported actions: all numerical parameters (rho rescales the box,
//! temp, tol, step sizes, vibtemp, angtemp, dihtemp, ecore, dcore), the
//! per-atom-type values (atlabel, epsilon, sigma, atmass, charge, qradius,
//! atsymbol, iexchange, natomsph, radminsph, radmaxsph, cgmolfac),
//! `box cubic` for a cubic box, `bond`, `label`, `addbond`, `removebond`,
//! `remove`.  Not translated (they stop with "not supported"): `tether`,
//! `box noncubic/change`, `findmols`, `getmolfs`, `add_w_rot`, `axisshift`,
//! `axisrotate`, `alignmole`, `setimage`, `unsetimage`, `molscalefc`, and the
//! `ls`/`dir`/`cd`/`pwd` shell commands.

use super::nwts::{list_int, scatt_len};
use super::ToolContext;
use crate::ato::perpd;
use crate::error::{EpsrError, Result};
use crate::fortran::format::{pad_to, ListWriter};
use crate::fortran::listread::{read_int, read_real, read_word, IoErr};
use crate::fortran::text::{concatenate, fstr, parse, Words};
use crate::neighbour::minimage;
use crate::state::State;
use crate::symbols::{LEN_COMMAND, LEN_VALUE};
use crate::symtab_tools::CHANGEATO;

fn eof() -> EpsrError {
    EpsrError::Input("changeato: unexpected end of standard input".into())
}

fn list_real(v: f32) -> String {
    let mut w = ListWriter::new();
    w.real(v);
    w.finish().trim_end().to_string()
}
fn list_chars(s: &str) -> String {
    let mut w = ListWriter::new();
    w.chars(s);
    w.finish().trim_end().to_string()
}

/// `readline(line, prompt)`: `read(*,'(a)') line` (end of file stops the Fortran program).
fn readline(ctx: &mut ToolContext) -> Result<String> {
    ctx.stdin.read_line().map_err(|_| eof())
}

/// Run the `changeato` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    st.install_tool_table(&CHANGEATO);
    let name = ctx.input_name().ok_or_else(|| EpsrError::Input("changeato: no .ato file given".into()))?.to_string();
    {
        let _f = st.setup_input_file(&name)?;
    }
    st.mol.ntype = 0;
    st.mol.nmoltype = 0;
    st.mol.sizefactor = 1.0;
    let fnameinp = st.sym.fnameinp.clone();
    st.log(&fnameinp);
    match st.minitc(&fnameinp) {
        Ok(()) => {}
        Err(EpsrError::Unsupported(s)) => return Err(EpsrError::Unsupported(s)),
        Err(_) => {
            st.log(&format!("\nProblems with .ato file {fnameinp}\nRoutine will exit now\n"));
            return Ok(());
        }
    }
    st.log(&format!("\nchangeato> There are {:>2} types of atom in this file\n", st.mol.ntype));
    st.sym.nfiles = st.mol.ntype;
    let nsymbol = st.sym.nsymbol;
    write_symbols_changeato(&mut st, 1, nsymbol);
    process_symbols(ctx, &mut st)
}

/// `subroutine write_symbols_changeato(nfirst,nlast)`.
fn write_symbols_changeato(st: &mut State, nfirst: i32, nlast: i32) {
    let nfirst = nfirst.max(1);
    let mut ifilecount = if nfirst <= st.sym.nsymboldata { 1 } else { st.sym.nfiles };
    let nlast = nlast.min(st.sym.nsymbol);
    let mut icount = nfirst - 1;
    while icount <= nlast {
        if icount >= st.sym.nsymboldata && ifilecount < st.sym.nfiles {
            ifilecount += 1;
            icount = st.sym.nsymboluptodata;
        }
        icount += 1;
        let iref = st.sym.ireff(icount, ifilecount);
        let sym = st.sym.symbol[icount].clone();
        let m = &st.mol;
        let f = ifilecount;
        let v: Option<String> = match sym.as_str() {
            "bond" | "label" | "addbond" | "removebond" => Some(st.sym.symbol_default[iref].clone()),
            "rho" => Some(list_real(st.dist.rho)),
            "temp" => Some(list_real(m.simtk)),
            "tol" => Some(list_real(m.tolout)),
            "stepmi" => Some(list_real(m.stepmi)),
            "stepri" => Some(list_real(m.stepri)),
            "stepr" => Some(list_real(m.stepr)),
            "stepm" => Some(list_real(m.stepm)),
            "vibtemp" => Some(list_real(m.tempb)),
            "angtemp" => Some(list_real(m.tempangreq)),
            "dihtemp" => Some(list_real(m.tempdihreq)),
            "ecore" => Some(list_real(m.ecore)),
            "dcore" => Some(list_real(m.dcore)),
            "box" => Some(list_chars(if st.lat.cubic { "cubic" } else { "noncubic" })),
            "natomtypes" => Some(list_int(m.ntype)),
            "atlabel" => Some(list_chars(&pad_to(&m.type_[f], 3))),
            "epsilon" => Some(list_real(m.epscor[f])),
            "sigma" => Some(list_real(m.sigcor[f])),
            "atmass" => Some(list_real(m.amassa[f])),
            "charge" => Some(list_real(m.charge[f])),
            "qradius" => Some(list_real(m.qradius[f])),
            "atsymbol" => Some(list_chars(&pad_to(&m.atsymbol[f], 3))),
            "iexchange" => Some(list_int(m.iexchange[f])),
            "natomsph" => Some(list_real(get_or0(&m.natomsphere, f))),
            "radminsph" => Some(list_real(get_or0(&m.radminsphere, f))),
            "radmaxsph" => Some(list_real(get_or0(&m.radmaxsphere, f))),
            "cgmolfac" => Some(list_real(get_or0(&m.pseudofacmol, f))),
            _ => None,
        };
        if sym == "natomtypes" {
            st.sym.nfiles = st.mol.ntype;
        }
        if let Some(v) = v {
            st.sym.symbol_value[iref] = v;
        }
    }
}

fn get_or0(a: &crate::fortran::F1<f32>, i: i32) -> f32 {
    if i >= 1 && (i as usize) <= a.len() {
        a[i]
    } else {
        0.0
    }
}

/// `subroutine process_symbols(5,6,ierrout,0)`.
fn process_symbols(ctx: &mut ToolContext, st: &mut State) -> Result<()> {
    let mut command = String::from("next");
    let mut icount: i32 = 0;
    let mut ifilecount: i32 = 1;
    st.sym.isave = 0;
    st.sym.inotrecog = 0;
    while command != "e" {
        let sometext = readline(ctx)?;
        let w = parse(&sometext, 11);
        process_commands(ctx, st, &sometext, &w, &mut command, &mut icount, &mut ifilecount)?;
    }
    if st.sym.isave != 1 {
        // 'Current data have not been saved. Press <Enter> to save, or q to exit without saving'
        let sometext = ctx.stdin.read_line().unwrap_or_default();
        let w = parse(&sometext, 5);
        let command = if w.nwords() > 0 { w.word(&sometext, 1).to_string() } else { String::new() };
        if command != "q" {
            save_files(ctx, st, false)?;
        }
    }
    Ok(())
}

/// `subroutine process_commands(5,6,command,icount,ifilecount,ierrout,0)`.
fn process_commands(
    ctx: &mut ToolContext,
    st: &mut State,
    sometext: &str,
    w: &Words,
    command: &mut String,
    icount: &mut i32,
    ifilecount: &mut i32,
) -> Result<()> {
    let nwords = w.nwords();
    *command = if nwords > 0 { fstr(w.word(sometext, 1), LEN_COMMAND) } else { "d".into() };
    if command == "d" || command == "u" {
        st.sym.symbol_save.clear();
        let s = &st.sym;
        if command == "d" {
            *icount += 1;
            if *icount > s.nsymboldata && *ifilecount < s.nfiles {
                *ifilecount += 1;
                *icount = s.nsymboluptodata + 1;
            } else if *icount > s.nsymbol {
                *icount = 0;
                *ifilecount = 1;
            }
        } else {
            *icount -= 1;
            if *icount < 0 {
                *icount = s.nsymbol;
                *ifilecount = s.nfiles;
            }
            if *icount <= s.nsymboluptodata && *ifilecount > 1 {
                *icount = s.nsymboldata;
                *ifilecount -= 1;
            }
        }
        *command = "next".into();
        return Ok(());
    }
    if command == "bye" {
        *command = "e".into();
    }
    if command == "default" {
        st.sym.inotrecog = 0;
        let s = &mut st.sym;
        let mut is = 0;
        *ifilecount = 1;
        while is < s.nsymboluptodata {
            is += 1;
            s.symbol_value[is] = s.symbol_default[is].clone();
        }
        while *ifilecount < s.nfiles {
            is += 1;
            let iref = is + s.nsymbolplus * (*ifilecount - 1);
            s.symbol_value[iref] = s.symbol_default[is].clone();
            if is >= s.nsymboldata && *ifilecount < s.nfiles {
                is = s.nsymboluptodata;
                *ifilecount += 1;
            }
        }
        *command = "next".into();
    } else if command == "ls" || command == "dir" || command == "cd" || command == "pwd" {
        return crate::error::unsupported("changeato: ls/dir/cd/pwd commands");
    } else if *command == st.sym.dataprompt {
        st.sym.inotrecog = 0;
        let mut ndata = 1;
        if nwords > 1 {
            if let Ok(n) = read_int(w.word(sometext, 2)) {
                ndata = n;
            }
        }
        if ndata < 1 {
            ndata = 1;
        }
        if ndata > st.sym.nfiles {
            ndata = st.sym.nfiles;
        }
        *ifilecount = ndata;
        *icount = st.sym.nsymboluptodata + 1;
        *command = "next".into();
    } else if command == "help" || command == "?" {
        st.sym.inotrecog = 0;
        *command = "next".into();
    } else if command == "save" {
        st.sym.inotrecog = 0;
        save_files(ctx, st, false)?;
        st.sym.isave = 1;
        *command = "next".into();
    } else if command == "forcesave" {
        st.sym.inotrecog = 0;
        save_files(ctx, st, true)?;
        st.sym.isave = 1;
        *command = "next".into();
    } else if command == "q" {
        st.sym.inotrecog = 0;
        st.sym.isave = 1;
        *command = "e".into();
    } else if command != "e" && command != "next" && command != "n" {
        let ifound = st.sym.search_for_symbol(command, *icount + 1);
        if ifound > 0 {
            let iref = st.sym.ireff(ifound, *ifilecount);
            st.sym.symbol_save = st.sym.symbol_value[iref].clone();
            if nwords > 1 {
                let (v, _) = concatenate(2, nwords, w, sometext, LEN_VALUE);
                st.sym.symbol_value[iref] = v;
                *icount = ifound;
                *command = "store".into();
            } else {
                *icount = ifound;
                *command = "next".into();
            }
            if ifound > st.sym.nsymboldata {
                *ifilecount = st.sym.nfiles;
            }
            if ifound < st.sym.nsymboluptodata + 1 {
                *ifilecount = 1;
            }
        } else if *icount > 0 {
            let iref = st.sym.ireff(*icount, *ifilecount);
            st.sym.symbol_save = st.sym.symbol_value[iref].clone();
            let (v, _) = concatenate(1, nwords, w, sometext, LEN_VALUE);
            st.sym.symbol_value[iref] = v;
            *command = "store".into();
        } else {
            st.log(&format!("\nprocess_commands> Command not recognised: {}:  {}", pad_to(command, LEN_COMMAND), crate::fortran::text::substr(sometext, 1, 40)));
            st.sym.inotrecog += 1;
            *command = "next".into();
        }
    }
    if command == "store" {
        st.sym.inotrecog = 0;
        st.sym.isave = 0;
        let ierr = check_symbol(ctx, st, icount, *ifilecount)?;
        *command = if ierr == 0 { "next".into() } else { "n".into() };
    }
    if command == "n" {
        st.sym.inotrecog = 0;
        if !st.sym.symbol_save.is_empty() && *icount > 0 {
            let iref = st.sym.ireff(*icount, *ifilecount);
            st.sym.symbol_value[iref] = st.sym.symbol_save.clone();
            check_symbol(ctx, st, icount, *ifilecount)?;
        } else {
            *icount += 1;
        }
        *command = "next".into();
    }
    Ok(())
}

/// `subroutine check_symbol(nin,nout,icount,ifilecount,ierr,ignore=0)` for changeato.
fn check_symbol(ctx: &mut ToolContext, st: &mut State, icount: &mut i32, ifilecount: i32) -> Result<i32> {
    if *icount < 1 {
        *icount = 1;
    }
    if *icount > st.sym.nsymbol {
        *icount = st.sym.nsymbol;
    }
    let ic = *icount;
    let iref = st.sym.ireff(ic, ifilecount);
    if st.sym.iflag[ic] > 0 && st.sym.symbol_value[iref] != st.sym.symbol_save {
        st.log(&format!("\n{} Cannot change this value!", st.sym.prompt));
        st.sym.symbol_value[iref] = st.sym.symbol_save.clone();
        return Ok(0);
    }
    let ty = st.sym.symbol_type[ic];
    if ty != 'c' {
        let sometext = st.sym.symbol_value[iref].clone();
        let w = parse(&sometext, 11);
        for iw in 1..=w.nwords() {
            let tok = w.word(&sometext, iw);
            let ok = if ty == 'i' { read_int(tok).is_ok() } else { read_real(tok).is_ok() };
            if !ok {
                st.log(&format!(
                    "\n{} Value {} is not {} - try again!",
                    st.sym.prompt,
                    tok,
                    if ty == 'i' { "an integer" } else { "a real number" }
                ));
                st.sym.symbol_value[iref] = st.sym.symbol_save.clone();
                return Ok(0);
            }
        }
    }
    if ic == st.sym.nsymboluptodata {
        if let Ok(nfilet) = read_int(&st.sym.symbol_value[iref]) {
            if nfilet != st.sym.nfiles {
                return crate::error::unsupported("changeato: changing the number of atom types");
            }
        }
    }
    change_symbol_changeato(ctx, st, ic, ifilecount)
}

/// `subroutine change_symbol_changeato(icount,ifilecount,ierr)`.
fn change_symbol_changeato(ctx: &mut ToolContext, st: &mut State, icount: i32, ifilecount: i32) -> Result<i32> {
    let mut ierr = 0;
    let iref = st.sym.ireff(icount, ifilecount);
    let sym = st.sym.symbol[icount].clone();
    let val = st.sym.symbol_value[iref].clone();
    let first_token_y = || {
        let w = parse(&val, 11);
        w.nwords() > 0 && w.word(&val, 1) == "y"
    };
    // read(symbol_value(iref),*,iostat=ierr) rtemp
    let rd_real = |ierr: &mut i32| -> f32 {
        match read_real(&val) {
            Ok(v) => v,
            Err(e) => {
                *ierr = e.iostat();
                0.0
            }
        }
    };
    let f = ifilecount;
    match sym.as_str() {
        "bond" => {
            if val == "y" {
                change_bond(ctx, st)?;
            }
            st.sym.symbol_value[iref] = st.sym.symbol_save.clone();
        }
        "label" => {
            if val == "y" {
                change_label(ctx, st)?;
            }
            st.sym.symbol_value[iref] = st.sym.symbol_save.clone();
        }
        "remove" => {
            if val == "y" {
                remove_molecule_types(ctx, st)?;
            }
            st.sym.symbol_value[iref] = st.sym.symbol_save.clone();
        }
        "addbond" => {
            if val == "y" {
                add_bond(ctx, st)?;
            }
            st.sym.symbol_value[iref] = st.sym.symbol_save.clone();
        }
        "removebond" => {
            if val == "y" {
                remove_bond(ctx, st)?;
            }
            st.sym.symbol_value[iref] = st.sym.symbol_save.clone();
        }
        "findmols" | "getmolfs" => {
            if val == "y" {
                return crate::error::unsupported(&format!("changeato: {sym}"));
            }
            st.sym.symbol_value[iref] = st.sym.symbol_save.clone();
        }
        "tether" | "add_w_rot" | "axisshift" | "axisrotate" | "alignmole" | "setimage" | "unsetimage" | "molscalefc" => {
            if first_token_y() {
                return crate::error::unsupported(&format!("changeato: {sym}"));
            }
            st.sym.symbol_value[iref] = st.sym.symbol_save.clone();
        }
        "rho" => {
            let rtemp = rd_real(&mut ierr);
            if ierr == 0 && rtemp > 0.0 {
                let l = &st.lat;
                let rhoold = st.mol.natomrcg / (l.alatx * l.blaty * l.clatz);
                st.dist.rho = rtemp;
                let vol = rhoold / st.dist.rho;
                let dmaxfac = vol.powf(1.0f32 / 3.0f32);
                let (phib, thetac, phic) = (st.lat.phib, st.lat.thetac, st.lat.phic);
                change_boxsize(st, dmaxfac, dmaxfac, dmaxfac, phib, thetac, phic)?;
            }
        }
        "temp" | "tol" | "stepmi" | "stepri" | "stepr" | "stepm" | "vibtemp" | "angtemp" | "dihtemp" | "ecore" | "dcore" => {
            let rtemp = rd_real(&mut ierr);
            if ierr == 0 {
                let m = &mut st.mol;
                match sym.as_str() {
                    "temp" => m.simtk = rtemp,
                    "tol" => m.tolout = rtemp,
                    "stepmi" => m.stepmi = rtemp,
                    "stepri" => m.stepri = rtemp,
                    "stepr" => m.stepr = rtemp,
                    "stepm" => m.stepm = rtemp,
                    "vibtemp" => m.tempb = rtemp,
                    "angtemp" => m.tempangreq = rtemp,
                    "dihtemp" => m.tempdihreq = rtemp,
                    "ecore" => m.ecore = rtemp,
                    _ => m.dcore = rtemp,
                }
            }
        }
        "box" => {
            let w = parse(&val, 11);
            let t1 = if w.nwords() > 0 { w.word(&val, 1).to_string() } else { String::new() };
            if (st.lat.cubic && t1 == "noncubic") || t1 == "change" || (!st.lat.cubic && t1 == "cubic") {
                return crate::error::unsupported("changeato: changing the box shape (non-cubic boxes)");
            }
            let l = &st.lat;
            st.dist.rho = st.mol.natomr as f32 / (l.alatx * l.blaty * l.clatz);
            let irefrho = st.sym.irefrho;
            st.sym.symbol_value[irefrho] = list_real(st.dist.rho);
            st.sym.symbol_value[iref] = list_chars(if st.lat.cubic { "cubic" } else { "noncubic" });
        }
        "atlabel" => match read_word(&val) {
            Ok(w) => {
                let ctemp = fstr(&w, 3);
                let m = &mut st.mol;
                for ia in 1..=m.natom {
                    if m.itype[ia] == f {
                        m.alabel[ia] = ctemp.clone();
                    }
                }
                m.type_[f] = ctemp;
            }
            Err(e) => ierr = e.iostat(),
        },
        "epsilon" | "sigma" | "atmass" | "charge" | "qradius" | "natomsph" | "radminsph" | "radmaxsph" => {
            let rtemp = rd_real(&mut ierr);
            if ierr == 0 {
                let m = &mut st.mol;
                match sym.as_str() {
                    "epsilon" => m.epscor[f] = rtemp,
                    "sigma" => m.sigcor[f] = rtemp,
                    "atmass" => m.amassa[f] = rtemp,
                    "charge" => m.charge[f] = rtemp,
                    "qradius" => m.qradius[f] = rtemp,
                    "natomsph" => m.natomsphere[f] = rtemp,
                    "radminsph" => m.radminsphere[f] = rtemp,
                    _ => m.radmaxsphere[f] = rtemp,
                }
            }
        }
        "atsymbol" => match read_word(&val) {
            Ok(w) => {
                let ctemp = fstr(&w, 3);
                let (_b, ifind) = scatt_len(&ctemp, 0);
                if ifind == 0 {
                    let irefl = st.sym.ireff(st.sym.nsymboluptodata + 1, f);
                    st.log(&format!(
                        "\nchangeato> Did not recognise chemical symbol {:<3} for atom type {:<3}",
                        ctemp,
                        crate::fortran::text::substr(&st.sym.symbol_value[irefl], 1, 3)
                    ));
                    ierr = 1;
                } else {
                    st.mol.atsymbol[f] = ctemp;
                }
            }
            Err(e) => ierr = e.iostat(),
        },
        "iexchange" => match read_int(&val) {
            Ok(v) => st.mol.iexchange[f] = v,
            Err(e) => ierr = e.iostat(),
        },
        "cgmolfac" => {
            let rtemp = rd_real(&mut ierr);
            if ierr == 0 && rtemp >= 1.0 {
                st.mol.pseudofacmol[f] = rtemp;
            } else {
                st.mol.pseudofacmol[f] = 1.0;
            }
        }
        _ => {}
    }
    if ierr != 0 {
        st.sym.symbol_value[iref] = st.sym.symbol_save.clone();
        st.sym.symbol_save.clear();
    }
    Ok(0)
}

/// `subroutine change_boxsize(afac,bfac,cfac,newphib,newthetac,newphic)` for
/// cubic boxes (`cosd`/`sind` are the gfortran intrinsics, exact at 0 and 90 degrees).
fn change_boxsize(st: &mut State, afac: f32, bfac: f32, cfac: f32, newphib: f32, newthetac: f32, newphic: f32) -> Result<()> {
    if !(newphib == 90.0 && newthetac == 0.0 && newphic == 0.0) {
        return crate::error::unsupported("changeato: change_boxsize for non-orthogonal boxes");
    }
    let frac = |l: &crate::state::Lattice, x: &mut f32, y: &mut f32, z: &mut f32| {
        let mut fc = *z / l.clatz;
        let mut fb = (*y - fc * l.claty) / l.blaty;
        let mut fa = (*x - fb * l.blatx - fc * l.clatx) / l.alatx;
        minimage(&mut fa, &mut fb, &mut fc);
        *z = fc;
        *y = fb;
        *x = fa;
    };
    let cart = |l: &crate::state::Lattice, x: &mut f32, y: &mut f32, z: &mut f32| {
        let (fa, fb, fc) = (*x, *y, *z);
        *z = fc * l.clatz;
        *y = fb * l.blaty + fc * l.claty;
        *x = fa * l.alatx + fb * l.blatx + fc * l.clatx;
    };
    let (m, l) = (&mut st.mol, &mut st.lat);
    for im in 1..=m.nmol {
        let (mut x, mut y, mut z) = (m.xadd[im], m.yadd[im], m.zadd[im]);
        frac(l, &mut x, &mut y, &mut z);
        m.xadd[im] = x;
        m.yadd[im] = y;
        m.zadd[im] = z;
        let (mut x, mut y, mut z) = (m.phix[im], m.phiy[im], m.phiz[im]);
        frac(l, &mut x, &mut y, &mut z);
        m.phix[im] = x;
        m.phiy[im] = y;
        m.phiz[im] = z;
    }
    l.alat *= afac;
    l.blat *= bfac;
    l.clat *= cfac;
    l.phib = newphib;
    l.thetac = newthetac;
    l.phic = newphic;
    // intrinsic cosd/sind: cosd(90)=0, sind(90)=1, cosd(0)=1, sind(0)=0
    l.alatx = l.alat;
    l.blatx = l.blat * 0.0;
    l.blaty = l.blat * 1.0;
    l.clatx = l.clat * 0.0 * 1.0;
    l.claty = l.clat * 0.0 * 0.0;
    l.clatz = l.clat * 1.0;
    for im in 1..=m.nmol {
        let (mut x, mut y, mut z) = (m.xadd[im], m.yadd[im], m.zadd[im]);
        cart(l, &mut x, &mut y, &mut z);
        m.xadd[im] = x;
        m.yadd[im] = y;
        m.zadd[im] = z;
        let (mut x, mut y, mut z) = (m.phix[im], m.phiy[im], m.phiz[im]);
        cart(l, &mut x, &mut y, &mut z);
        m.phix[im] = x;
        m.phiy[im] = y;
        m.phiz[im] = z;
    }
    l.alatc = [l.alatx, 0.0, 0.0];
    l.blatc = [l.blatx, l.blaty, 0.0];
    l.clatc = [l.clatx, l.claty, l.clatz];
    l.bxax = l.blatx / l.alatx;
    l.cxax = l.clatx / l.alatx;
    l.cyby = l.claty / l.blaty;
    let perpab = perpd(&l.alatc, &l.blatc, &l.clatc);
    let perpbc = perpd(&l.blatc, &l.clatc, &l.alatc);
    let perpca = perpd(&l.clatc, &l.alatc, &l.blatc);
    m.dmax = perpab.min(perpbc).min(perpca);
    st.dist.cellmax = m.dmax;
    l.cubic = !(l.alat != l.blat || l.blat != l.clat || l.clat != l.alat || l.phib != 90.0 || l.thetac != 0.0 || l.phic != 0.0);
    if !l.cubic {
        return crate::error::unsupported("changeato: box became non-cubic");
    }
    Ok(())
}

/// One list-directed READ of integers/reals from standard input (`read(5,*) a,b,...`).
fn read_items(ctx: &mut ToolContext, kinds: &str) -> Result<(Vec<i32>, Vec<f32>)> {
    let mut r = ctx.stdin.list();
    let mut iv = Vec::new();
    let mut rv = Vec::new();
    for k in kinds.chars() {
        if k == 'i' {
            iv.push(r.int().map_err(|e| if e == IoErr::Eof { eof() } else { e.into() })?.unwrap_or(0));
        } else {
            rv.push(r.real().map_err(|e| if e == IoErr::Eof { eof() } else { e.into() })?.unwrap_or(0.0));
        }
    }
    Ok((iv, rv))
}

/// `subroutine change_bond()`.
fn change_bond(ctx: &mut ToolContext, st: &mut State) -> Result<()> {
    let (iv, rv) = read_items(ctx, "iiiir")?;
    let (mut ncha1, mut ncha2, mut nat1, mut nat2, mut bonnew) = (iv[0], iv[1], iv[2], iv[3], rv[0]);
    while ncha1 != 0 {
        let m = &mut st.mol;
        for ia in 1..=m.natom {
            let mmm = m.mref[ia];
            let na1 = ia + 1 - m.naf[mmm];
            if m.nbond[ia] > 0 && (na1 == nat1 || na1 == nat2) {
                for ib in 1..=m.nbond[ia] {
                    let ja = m.iabond[(ia, ib)];
                    let na2 = ja + 1 - m.naf[mmm];
                    if m.itype[ia] == ncha1 && m.itype[ja] == ncha2 && na1 == nat1 && na2 == nat2 {
                        m.tbond[(ia, ib)] = bonnew;
                    }
                    if m.itype[ia] == ncha2 && m.itype[ja] == ncha1 && na1 == nat2 && na2 == nat1 {
                        m.tbond[(ia, ib)] = bonnew;
                    }
                }
            }
        }
        let (iv, rv) = read_items(ctx, "iiiir")?;
        ncha1 = iv[0];
        ncha2 = iv[1];
        nat1 = iv[2];
        nat2 = iv[3];
        bonnew = rv[0];
    }
    Ok(())
}

/// `subroutine add_bond()`.
fn add_bond(ctx: &mut ToolContext, st: &mut State) -> Result<()> {
    let (iv, rv) = read_items(ctx, "iiiir")?;
    let (mut ncha1, mut ncha2, mut nat1, mut nat2, mut bonnew) = (iv[0], iv[1], iv[2], iv[3], rv[0]);
    while ncha1 > 0 {
        let m = &mut st.mol;
        for ia in 1..=m.natom {
            let mmm = m.mref[ia];
            let na1 = ia + 1 - m.naf[mmm];
            if m.itype[ia] == ncha1 && na1 == nat1 {
                m.nbond[ia] += 1;
                let ib = m.nbond[ia];
                m.ensure_bonds(ib as usize);
                m.iabond[(ia, ib)] = m.naf[mmm] + nat2 - 1;
                m.tbond[(ia, ib)] = bonnew;
            }
            if m.itype[ia] == ncha2 && na1 == nat2 {
                m.nbond[ia] += 1;
                let ib = m.nbond[ia];
                m.ensure_bonds(ib as usize);
                m.iabond[(ia, ib)] = m.naf[mmm] + nat1 - 1;
                m.tbond[(ia, ib)] = bonnew;
            }
        }
        let (iv, rv) = read_items(ctx, "iiiir")?;
        ncha1 = iv[0];
        ncha2 = iv[1];
        nat1 = iv[2];
        nat2 = iv[3];
        bonnew = rv[0];
    }
    Ok(())
}

/// `subroutine remove_bond()`.
fn remove_bond(ctx: &mut ToolContext, st: &mut State) -> Result<()> {
    let (iv, _) = read_items(ctx, "iiii")?;
    let (mut ncha1, mut ncha2, mut nat1, mut nat2) = (iv[0], iv[1], iv[2], iv[3]);
    while ncha1 > 0 {
        let m = &mut st.mol;
        for ia in 1..=m.natom {
            let mmm = m.mref[ia];
            let na1 = ia - m.naf[mmm] + 1;
            let target = if m.itype[ia] == ncha1 && na1 == nat1 && m.nbond[ia] > 0 {
                Some((ncha2, nat2))
            } else if m.itype[ia] == ncha2 && na1 == nat2 && m.nbond[ia] > 0 {
                Some((ncha1, nat1))
            } else {
                None
            };
            if let Some((tt, tn)) = target {
                let mut ibsav = 0;
                for ib in 1..=m.nbond[ia] {
                    let ibref = m.iabond[(ia, ib)];
                    let na2 = ibref - m.naf[mmm] + 1;
                    if m.itype[ibref] == tt && na2 == tn {
                        ibsav = ib;
                        break;
                    }
                }
                if ibsav > 0 {
                    m.nbond[ia] -= 1;
                    if m.nbond[ia] >= ibsav {
                        for ib in ibsav..=m.nbond[ia] {
                            let ibref = ib + 1;
                            m.iabond[(ia, ib)] = m.iabond[(ia, ibref)];
                            m.tbond[(ia, ib)] = m.tbond[(ia, ibref)];
                        }
                    }
                }
            }
        }
        let (iv, _) = read_items(ctx, "iiii")?;
        ncha1 = iv[0];
        ncha2 = iv[1];
        nat1 = iv[2];
        nat2 = iv[3];
    }
    Ok(())
}

/// `subroutine change_label()`.
fn change_label(ctx: &mut ToolContext, st: &mut State) -> Result<()> {
    let (iv, _) = read_items(ctx, "ii")?;
    let (it, iamol) = (iv[0], iv[1]);
    if it < 1 || it > st.mol.ntype || iamol < 1 || iamol > st.mol.nmoltype {
        return Ok(());
    }
    let mut changeselect: Vec<i32> = Vec::new();
    {
        let m = &st.mol;
        for im in 1..=m.nmol {
            if m.mtype[im] == iamol {
                for ia in m.naf[im]..=m.nal[im] {
                    if m.itype[ia] == it {
                        changeselect.push(ia);
                    }
                }
            }
        }
    }
    let nthistype = changeselect.len() as i32;
    if nthistype <= 0 {
        st.log(&format!("\nchange_label> Warning! No atoms of type {:<3} at atom number {:>4}", st.mol.type_[it], iamol));
        return Ok(());
    }
    // read(5,'(a)') a1  (character*3)
    let line = ctx.stdin.read_line().map_err(|_| eof())?;
    let a1 = fstr(&line, 3);
    let (_, rv) = read_items(ctx, "r")?;
    let fracchange = rv[0];
    let notochange = nthistype.min((fracchange * nthistype as f32) as i32);
    let mut changethis: Vec<i32> = Vec::new();
    let mut ichange = 0;
    let mut icount: usize = 0;
    while ichange < notochange {
        icount += 1;
        if icount as i32 > nthistype {
            icount = 1;
        }
        let r = st.rng.ran1();
        if (r < fracchange || fracchange >= 1.0) && changeselect[icount - 1] > 0 {
            ichange += 1;
            changethis.push(changeselect[icount - 1]);
            changeselect[icount - 1] = 0;
        }
    }
    let nchange = ichange;
    if nchange <= 0 {
        st.log("\nchange_label> Warning: No valid atoms to change labels found!!\n");
        return Ok(());
    }
    for &ia in changethis.iter() {
        let m = &st.mol;
        let epsit = m.epscor[it];
        let sigit = m.sigcor[it];
        let amass = m.amassa[it];
        let qt = m.charge[it];
        let rqt = m.qradius[it];
        let ats = m.atsymbol[it].clone();
        let iexch = m.iexchange[it];
        let natomsph = get_or0(&m.natomsphere, it);
        let rhosph = get_or0(&m.rhosphere, it);
        let radminsph = get_or0(&m.radminsphere, it);
        let radmaxsph = get_or0(&m.radmaxsphere, it);
        let pseudofacsph = get_or0(&m.pseudofacmol, it);
        st.mol.alabel[ia] = a1.clone();
        st.ftype(ia);
        st.ctype(&a1, &ats, iexch, epsit, sigit, amass, qt, rqt, natomsph, rhosph, radminsph, radmaxsph, pseudofacsph);
    }
    cleanup_types(st);
    let (a, b) = (st.sym.nsymboluptodata, st.sym.nsymbol);
    write_symbols_changeato(st, a, b);
    Ok(())
}

/// `subroutine remove_molecule_types()`.
fn remove_molecule_types(ctx: &mut ToolContext, st: &mut State) -> Result<()> {
    let nmt = st.mol.nmoltype as usize;
    let mut findtype = vec![0i32; nmt.max(1) + 1];
    {
        let m = &st.mol;
        for im in 1..=m.nmol {
            let mut imoltype = 0;
            let mut ifound = 0;
            while imoltype < m.nmoltype && ifound == 0 {
                imoltype += 1;
                if m.alabel[m.naf[im]] == m.moltype[imoltype] {
                    ifound = imoltype;
                }
            }
            if ifound > 0 {
                findtype[imoltype as usize] += 1;
            }
        }
    }
    let mut readagain = true;
    while readagain {
        let sometext = match ctx.stdin.read_line() {
            Ok(l) => l,
            Err(_) => return Err(eof()),
        };
        if !sometext.trim_end().is_empty() {
            let w = parse(&sometext, 5);
            readagain = w.word(&sometext, 1) != "q";
            if readagain {
                let mtyperead = read_int(w.word(&sometext, 1))?;
                if mtyperead < 1 || mtyperead > st.mol.nmoltype {
                    return Err(EpsrError::Input(format!("changeato remove: molecule type {mtyperead} out of range")));
                }
                let mut ntoremove = -1;
                if w.nwords() > 1 {
                    if let Ok(n) = read_int(w.word(&sometext, 2)) {
                        ntoremove = n;
                    }
                }
                if ntoremove < 0 {
                    ntoremove = findtype[mtyperead as usize];
                }
                let moltyperead = st.mol.moltype[mtyperead].clone();
                let mut im = 1;
                let mut nremove = 0;
                while im <= st.mol.nmol {
                    if st.mol.alabel[st.mol.naf[im]] == moltyperead && nremove < ntoremove {
                        remove_molecule(st, im);
                        nremove += 1;
                        findtype[mtyperead as usize] -= 1;
                    } else {
                        im += 1;
                    }
                }
            }
        }
    }
    let (a, b) = (st.sym.nsymboluptodata, st.sym.nsymbol);
    write_symbols_changeato(st, a, b);
    Ok(())
}

/// `subroutine remove_molecule(iremove)` (reproduces the Fortran exactly: the
/// bond partners and rotation-group atom numbers of the molecules that move
/// down are not renumbered, nor are the per-atom `qtype`/`imageatom` and the
/// per-molecule `nstepreject`/`scalefacmol`/`tetheredatom`).
fn remove_molecule(st: &mut State, iremove: i32) {
    let m = &mut st.mol;
    if iremove > 0 && iremove <= m.nmol && m.nmol > 1 {
        let mut ig = 0;
        let mut im = 1;
        while im < iremove {
            ig += m.ngrot[im];
            im += 1;
        }
        m.natom -= m.nato[im];
        let mut ia = m.naf[im] - 1;
        while im < m.nmol {
            m.naf[im] = ia + 1;
            m.nato[im] = m.nato[im + 1];
            m.xadd[im] = m.xadd[im + 1];
            m.yadd[im] = m.yadd[im + 1];
            m.zadd[im] = m.zadd[im + 1];
            m.phix[im] = m.phix[im + 1];
            m.phiy[im] = m.phiy[im + 1];
            m.phiz[im] = m.phiz[im + 1];
            m.tethered[im] = m.tethered[im + 1];
            let mut ic = 0;
            let mut iaref = m.naf[im + 1] - 1;
            while ic < m.nato[im] {
                ic += 1;
                ia += 1;
                iaref += 1;
                m.alabel[ia] = m.alabel[iaref].clone();
                m.xa[ia] = m.xa[iaref];
                m.ya[ia] = m.ya[iaref];
                m.za[ia] = m.za[iaref];
                m.nbond[ia] = m.nbond[iaref];
                for ib in 1..=m.nbond[ia] {
                    m.iabond[(ia, ib)] = m.iabond[(iaref, ib)];
                    m.tbond[(ia, ib)] = m.tbond[(iaref, ib)];
                }
            }
            m.nal[im] = ia;
            let mut igref = ig + m.ngrot[im];
            m.ngrot[im] = m.ngrot[im + 1];
            let mut ib = 0;
            while ib < m.ngrot[im] {
                ib += 1;
                ig += 1;
                igref += 1;
                m.rotgrp[ig] = m.rotgrp[igref].clone();
                m.labaxi[(ig, 1)] = m.labaxi[(igref, 1)];
                m.labaxi[(ig, 2)] = m.labaxi[(igref, 2)];
                m.nrot[ig] = m.nrot[igref];
                m.rotf[ig] = m.rotf[igref];
                m.rotl[ig] = m.rotl[igref];
            }
            im += 1;
        }
        m.nmol -= 1;
    }
    cleanup_types(st);
}

/// `subroutine cleanup_types()`.
fn cleanup_types(st: &mut State) {
    let m = &mut st.mol;
    let nt = m.ntype.max(0) as usize;
    let mut findtype = vec![0i32; nt + 2];
    for im in 1..=m.nmol {
        for ia in m.naf[im]..=m.nal[im] {
            let mut ifound = 0;
            let mut it = 0;
            while it < m.ntype && ifound == 0 {
                it += 1;
                if m.alabel[ia] == m.type_[it] {
                    ifound = it;
                }
            }
            if ifound > 0 {
                findtype[ifound as usize] += 1;
            }
        }
    }
    let mut it = m.ntype;
    while it > 0 {
        if findtype[it as usize] <= 0 {
            let mut jt = it;
            while jt < m.ntype {
                m.type_[jt] = m.type_[jt + 1].clone();
                m.epscor[jt] = m.epscor[jt + 1];
                m.sigcor[jt] = m.sigcor[jt + 1];
                m.amassa[jt] = m.amassa[jt + 1];
                m.charge[jt] = m.charge[jt + 1];
                m.qradius[jt] = m.qradius[jt + 1];
                m.atsymbol[jt] = m.atsymbol[jt + 1].clone();
                m.iexchange[jt] = m.iexchange[jt + 1];
                m.rhosphere[jt] = m.rhosphere[jt + 1];
                m.natomsphere[jt] = m.natomsphere[jt + 1];
                m.radminsphere[jt] = m.radminsphere[jt + 1];
                m.radmaxsphere[jt] = m.radmaxsphere[jt + 1];
                jt += 1;
            }
            m.ntype -= 1;
        }
        it -= 1;
    }
    for ia in 1..=m.natom {
        let mut it = 0;
        let mut ifound = 0;
        while it < m.ntype && ifound == 0 {
            it += 1;
            if m.alabel[ia] == m.type_[it] {
                ifound = it;
            }
        }
        m.itype[ia] = ifound;
    }
    let nmt = m.nmoltype.max(0) as usize;
    let mut findtype = vec![0i32; nmt + 2];
    for im in 1..=m.nmol {
        let mut imoltype = 0;
        let mut ifound = 0;
        while imoltype < m.nmoltype && ifound == 0 {
            imoltype += 1;
            if m.alabel[m.naf[im]] == m.moltype[imoltype] {
                ifound = imoltype;
            }
        }
        if ifound > 0 {
            findtype[imoltype as usize] += 1;
        }
    }
    let mut imoltype = m.nmoltype;
    while imoltype > 0 {
        if findtype[imoltype as usize] <= 0 {
            let mut jt = imoltype;
            while jt < m.nmoltype {
                m.moltype[jt] = m.moltype[jt + 1].clone();
                m.moltethered[jt] = m.moltethered[jt + 1];
                m.moltemplate[jt] = m.moltemplate[jt + 1].clone();
                m.wholerot[jt] = m.wholerot[jt + 1];
                m.wholetrans[jt] = m.wholetrans[jt + 1];
                jt += 1;
            }
            m.nmoltype -= 1;
        }
        imoltype -= 1;
    }
    for im in 1..=m.nmol {
        let mut imoltype = 0;
        let mut ifound = 0;
        while imoltype < m.nmoltype && ifound == 0 {
            imoltype += 1;
            if m.alabel[m.naf[im]] == m.moltype[imoltype] {
                ifound = imoltype;
            }
        }
        m.mtype[im] = ifound;
    }
}

/// `subroutine save_files(nin,nout)` / `force_save_files` followed by
/// `save_files_changeato` (= `moloutc(fnameinp)`).
fn save_files(ctx: &mut ToolContext, st: &mut State, force: bool) -> Result<()> {
    let extinp = st.sym.extinp.clone();
    let mut fname_inp = st.sym.fname_inp.clone();
    loop {
        if fname_inp.is_empty() {
            loop {
                let sometext = readline(ctx)?;
                let w = parse(&sometext, 5);
                if w.nwords() > 0 {
                    fname_inp = concatenate(1, w.nwords(), &w, &sometext, 256).0;
                    break;
                }
            }
        }
        // strip any extension (first '.') and add .ato
        let base = match fname_inp.find('.') {
            Some(p) if p > 0 => fname_inp[..p].to_string(),
            _ => fname_inp.trim_end().to_string(),
        };
        fname_inp = format!("{base}{extinp}");
        st.log(&format!("\n{} Writing to input file \"{}\"", st.sym.prompt, fname_inp));
        let full = st.sym.full(&fname_inp).to_string_lossy().into_owned();
        if force || !std::path::Path::new(&full).exists() {
            if !force {
                // open(14,file=fnameinp,status='new') creates the file
                std::fs::write(&full, b"")?;
            }
            st.sym.fname_inp = fname_inp.clone();
            st.sym.fnameinp = full;
            break;
        }
        let ans = loop {
            let sometext = readline(ctx)?;
            let w = parse(&sometext, 5);
            if w.nwords() > 0 {
                break w.word(&sometext, 1).to_string();
            }
        };
        if ans == "y" {
            st.sym.fname_inp = fname_inp.clone();
            st.sym.fnameinp = full;
            break;
        }
        fname_inp.clear();
    }
    let f = st.sym.fnameinp.clone();
    st.moloutc(&f)
}

/// Program-specific branch of `check_symbol` (`process_symbols.f`,
/// `symbol_setup.eq.'changeato'`).  changeato never reads a keyword file, so
/// this branch is reached only through the interactive loop above
/// ([`change_symbol_changeato`]); when called from the generic reader it does nothing.
pub fn check_symbol_extra(_st: &mut State, _icount: i32, _ifilecount: i32, _iref: i32) -> Result<i32> {
    Ok(0)
}
