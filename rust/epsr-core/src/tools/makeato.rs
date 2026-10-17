//! `makeato` — translation of `fortran_src/makeato.f`.
//!
//! STATUS: translated; verified against the Fortran (EPSR-test-ex30).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> makeato [<name>]`
//! The answers (name if not given, epsilon, sigma, mass, charge, chemical
//! symbol, temperature, number density) are read from standard input.
//!
//! The program builds a `.ato` file with one molecule containing a single atom
//! whose label is the first three characters of the file name.
//!
//! NOTE: the Fortran program never allocates the coarse-graining arrays
//! (`natomsphere`, `pseudofacmol`, normally allocated by `ctype`), and
//! `moloutc` reads them before allocating them, so the reference program
//! crashes with SIGSEGV inside `moloutc`, leaving an empty `.ato` file.  The
//! Rust `moloutc` reproduces that outcome (empty file, error return).

#[path = "setup_common.rs"]
pub mod setup_common;

use super::ToolContext;
use crate::error::Result;
use crate::fortran::fstr;
use crate::fortran::text::{findex, parse};
use setup_common::{command_text, eof_err, read_real5, read_text};

/// Run the `makeato` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    let (sometext, w) = command_text("makeato", ctx);
    let readfile: String;
    if w.nwords() < 2 {
        st.log("\nmakeato> Type name of .ato file to create: ");
        let t = read_text(ctx)?;
        let w2 = parse(&t, 5);
        if w2.nwords() <= 0 {
            st.log("\nmakeato> No .ato file specified! - routine will exit\n");
            return Ok(());
        }
        readfile = w2.word(&t, 1).to_string();
    } else {
        readfile = w.word(&sometext, 2).to_string();
    }
    let mut readfile = readfile;
    let ilenuptodot = findex(&readfile, ".") - 1;
    if ilenuptodot > 0 {
        readfile = readfile[..ilenuptodot as usize].to_string();
    }
    let fnameato = format!("{}{}.ato", ctx.dir, readfile);
    st.log("\nmakeato> MAKEATO will create a .ato file containing 1 molecule\n which itself contains just 1 atom\n");
    let m = &mut st.mol;
    let mut natom = 0;
    m.ntype = 1;
    // type(1)=readfile(1:3)
    m.type_[1] = fstr(&readfile, 3);
    let ntype = m.ntype;
    for i in 1..=ntype {
        st.log(&format!("\nmakeato> Atom type {:>2} has label: {}", i, st.mol.type_[i]));
        let ntypes_i = 1;
        st.log(&format!("\nmakeato> For atom type {}:-", st.mol.type_[i]));
        st.log("\n Give its epsilon [kJ/mol]: ");
        st.mol.epscor[i] = read_real5(ctx)?;
        st.log("\n Give its sigma [A]: ");
        st.mol.sigcor[i] = read_real5(ctx)?;
        st.log("\n Give its atomic mass [amu]: ");
        st.mol.amassa[i] = read_real5(ctx)?;
        st.log("\n Give its electronic charge [e] ");
        st.mol.charge[i] = read_real5(ctx)?;
        st.log("\n Give its chemical symbol: ");
        {
            let mut r = ctx.stdin.list();
            if let Some(s) = r.word().map_err(eof_err)? {
                st.mol.atsymbol[i] = fstr(&s, 3);
            }
        }
        let m = &mut st.mol;
        m.qradius[i] = 0.0;
        m.iexchange[i] = 0;
        for _j in 1..=ntypes_i {
            natom += 1;
            m.ensure_atoms(natom as usize);
            m.alabel[natom] = m.type_[i].clone();
            m.xa[natom] = 0.0;
            m.ya[natom] = 0.0;
            m.za[natom] = 0.0;
            m.itype[natom] = i;
            m.nbond[natom] = 0;
        }
    }
    st.mol.natom = natom;
    // (the bond and bond-angle loops are disabled in the Fortran: ia1=ia2=0)
    let m = &mut st.mol;
    let ig = 0;
    m.ngrrot = 0;
    m.nrottot = 0;
    m.ngrot[1] = ig;
    m.nmol = 1;
    m.xadd[1] = 0.0;
    m.yadd[1] = 0.0;
    m.zadd[1] = 0.0;
    m.phix[1] = 0.0;
    m.phiy[1] = 0.0;
    m.phiz[1] = 0.0;
    m.nato[1] = natom;
    m.naf[1] = 1;
    m.nal[1] = natom;
    m.stepmi = 0.1;
    m.stepri = if m.ngrot[1] > 0 { 1.0 } else { 0.0 };
    m.stepr = 1.0;
    m.stepm = 2.0;
    m.tempb = 65.0;
    st.log("\nmakeato> Type the temperature [K]: ");
    st.mol.simtk = read_real5(ctx)?;
    st.log("\nmakeato> and atomic number density [per A**3]: ");
    let rho = {
        let mut v = 0f32;
        ctx.stdin.list().get_r(&mut v).map_err(eof_err)?;
        v
    };
    let m = &mut st.mol;
    m.tempb = 65.0;
    m.ecore = 0.0;
    m.dcore = 0.0;
    let volume = natom as f32 / rho;
    m.volume = volume;
    m.dmax = volume.powf(1.0f32 / 3.0f32);
    st.lat.cubic = true;
    m.dmaxmol = m.dmax;
    m.nmoltype = 1;
    m.moltemplate[1] = readfile.clone();
    m.wholerot[1] = 0.0;
    m.wholetrans[1] = 0.0;
    m.wholetrans[1] = 0.01;
    st.moloutc(&fnameato)
}
