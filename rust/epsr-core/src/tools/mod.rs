//! The auxiliary EPSR programs (setup utilities and analysis routines).
//!
//! Each program of the Fortran EPSR suite is a module here with a
//! `run(&mut ToolContext)` entry point; `run_tool` dispatches on the command
//! name used by the EPSR shell (`epsr-rs <folder> <command> <file>`).
//!
//! Programs driven by a keyword input file (`.PARTIALS.dat`, `.NWTS.dat`, ...)
//! use [`ToolContext::open_input`], which installs the program's keyword table
//! (`symtab_tools.rs`) and reads the file through the same machinery as
//! `epsr` (`symbols.rs`).  Programs that prompt for their parameters read them
//! from `ctx.stdin`.
//!
//! The status of each program (translated / verified) is listed in
//! `docs/TOOLS.md`.

use crate::error::{EpsrError, Result};
use crate::fortran::listread::TextFile;
use crate::state::State;
use std::path::PathBuf;

pub mod partials;
pub mod coord;
pub mod triangles;
pub mod torangles;
pub mod nwts;
pub mod xwts;
pub mod changeato;
pub mod voids;
pub mod clusters;
pub mod chains;
pub mod rings;
pub mod lsi;
pub mod sdf;
pub mod sharm;
pub mod sdfcube;
pub mod mapgrentropy;
pub mod mapgr;
pub mod fluctuations;
pub mod local_density;
pub mod epsrwts;
pub mod epsrwtsx;
pub mod makeato;
pub mod mixato;
pub mod fmole;
pub mod writexyz;
pub mod readxyz;
pub mod readpdb;
pub mod addato;
pub mod dockato;
pub mod makelattice;
pub mod makelatticeato;
pub mod makemole;
pub mod randomise;
pub mod bonds;
pub mod alignmole;

/// Everything a program needs from its invocation.
pub struct ToolContext {
    /// Folder holding the input files, with a trailing `/`.
    pub dir: String,
    /// Arguments after the command name (normally the input file name).
    pub args: Vec<String>,
    /// Standard input (answers to prompts of interactive programs).
    pub stdin: TextFile,
    /// Echo progress messages.
    pub verbose: bool,
}

impl ToolContext {
    pub fn new(dir: &std::path::Path, args: Vec<String>, stdin: TextFile, verbose: bool) -> ToolContext {
        let mut d = dir.to_string_lossy().into_owned();
        if !d.ends_with('/') {
            d.push('/');
        }
        ToolContext { dir: d, args, stdin, verbose }
    }
    /// A fresh simulation state with the directory and verbosity set.
    pub fn new_state(&self) -> State {
        let mut st = State::new();
        st.verbose = self.verbose;
        st.sym.directory_name = self.dir.clone();
        st
    }
    /// Full path of a file in the folder.
    pub fn path(&self, name: &str) -> PathBuf {
        if name.contains('/') {
            PathBuf::from(name)
        } else {
            PathBuf::from(format!("{}{}", self.dir, name))
        }
    }
    /// The first argument (input file name), if any.
    pub fn input_name(&self) -> Option<&str> {
        self.args.first().map(|s| s.as_str())
    }
    /// `setup_symbols_<prog>` + `setup_input_file` + `read_input_file` for a
    /// keyword-driven program.  Returns `Ok(false)` if the input file does not
    /// exist (the Fortran then continues with default values).
    pub fn open_input(&self, st: &mut State, table: &crate::symbols::ToolTable) -> Result<bool> {
        st.install_tool_table(table);
        let name = self.input_name().ok_or_else(|| EpsrError::Input(format!("{}: no input file given", table.name)))?.to_string();
        match st.setup_input_file(&name) {
            Ok(mut f) => {
                st.read_input_file(&mut f)?;
                Ok(true)
            }
            Err(EpsrError::Io(_)) => Ok(false),
            Err(e) => Err(e),
        }
    }
}

/// Command names handled by [`run_tool`] (in addition to `epsr`).
pub const COMMANDS: &[&str] = &["partials", "coord", "triangles", "torangles", "nwts", "xwts", "changeato", "voids", "clusters", "chains", "rings", "lsi", "sdf", "sharm", "sdfcube", "mapgrentropy", "mapgr", "fluctuations", "local_density", "epsrwts", "epsrwtsx", "makeato", "mixato", "fmole", "writexyz", "readxyz", "readpdb", "addato", "dockato", "makelattice", "makelatticeato", "makemole", "randomise", "bonds", "alignmole"];

/// Run the program `command`.
pub fn run_tool(command: &str, ctx: &mut ToolContext) -> Result<()> {
    match command {
        "partials" => partials::run(ctx),
        "coord" => coord::run(ctx),
        "triangles" => triangles::run(ctx),
        "torangles" => torangles::run(ctx),
        "nwts" => nwts::run(ctx),
        "xwts" => xwts::run(ctx),
        "changeato" => changeato::run(ctx),
        "voids" => voids::run(ctx),
        "clusters" => clusters::run(ctx),
        "chains" => chains::run(ctx),
        "rings" => rings::run(ctx),
        "lsi" => lsi::run(ctx),
        "sdf" => sdf::run(ctx),
        "sharm" => sharm::run(ctx),
        "sdfcube" => sdfcube::run(ctx),
        "mapgrentropy" => mapgrentropy::run(ctx),
        "mapgr" => mapgr::run(ctx),
        "fluctuations" => fluctuations::run(ctx),
        "local_density" => local_density::run(ctx),
        "epsrwts" => epsrwts::run(ctx),
        "epsrwtsx" => epsrwtsx::run(ctx),
        "makeato" => makeato::run(ctx),
        "mixato" => mixato::run(ctx),
        "fmole" => fmole::run(ctx),
        "writexyz" => writexyz::run(ctx),
        "readxyz" => readxyz::run(ctx),
        "readpdb" => readpdb::run(ctx),
        "addato" => addato::run(ctx),
        "dockato" => dockato::run(ctx),
        "makelattice" => makelattice::run(ctx),
        "makelatticeato" => makelatticeato::run(ctx),
        "makemole" => makemole::run(ctx),
        "randomise" => randomise::run(ctx),
        "bonds" => bonds::run(ctx),
        "alignmole" => alignmole::run(ctx),
        _ => Err(EpsrError::Input(format!("unknown command {command:?}"))),
    }
}

/// Program-specific branches of `check_symbol`, dispatched on `symbol_setup`.
pub fn check_symbol_extra(st: &mut State, setup: &str, icount: i32, ifilecount: i32, iref: i32) -> Result<i32> {
    match setup {
        "nwts" => nwts::check_symbol_extra(st, icount, ifilecount, iref),
        "xwts" => xwts::check_symbol_extra(st, icount, ifilecount, iref),
        "changeato" => changeato::check_symbol_extra(st, icount, ifilecount, iref),
        "sdf" => sdf::check_symbol_extra(st, icount, ifilecount, iref),
        "sharm" => sharm::check_symbol_extra(st, icount, ifilecount, iref),
        "mapgrentropy" => mapgrentropy::check_symbol_extra(st, icount, ifilecount, iref),
        "mapgr" => mapgr::check_symbol_extra(st, icount, ifilecount, iref),
        _ => Ok(0),
    }
}
