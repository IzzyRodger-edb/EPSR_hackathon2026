//! Reading and writing the `.EPSR.inp` control file.
//!
//! Translation of:
//! * `setup_symbols_epsr`            (`setup_epsr_rev1.f`)
//! * `setup_input_file`              (`file_and_directory_commands.f`)
//! * `read_input_file`               (`read_and_save_files.f`)
//! * `process_commands`, `check_symbol` (`process_symbols.f`) — the branches
//!   used when reading a file (`ignore = 1`, non-interactive)
//! * the `.EPSR.inp` part of `save_files_others` (`read_and_save_files.f`)
//!
//! EPSR represents its input file as a table of *symbols* (keywords) and
//! *values* (strings).  Symbols `1..=nsymboluptodata` are global, symbols
//! `nsymboluptodata+1..=nsymboldata` are replicated for each dataset ("data n"
//! blocks), and symbols beyond `nsymboldata` hold one entry per atom-pair
//! potential read from the `.pcof` file.  [`Symbols::ireff`] maps
//! (symbol, dataset) to a storage slot.  Values are parsed by the consumers
//! (`runepsr17`, `inita`, ...) by position, so the table layout is part of the
//! file format.

use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, FV};
use crate::fortran::listread::{read_int, read_real, ListRead, TextFile};
use crate::fortran::output::OutFile;
use crate::fortran::text::{concatenate, fstr, indexlastdot, parse, substr};
use crate::state::State;
use crate::symtab::*;

/// Declared lengths of the character variables involved.
pub const LEN_SYMBOL: usize = 11;
pub const LEN_VALUE: usize = 256;
pub const LEN_COMMAND: usize = 20;

impl State {
    /// `subroutine setup_symbols_epsr()`: install the keyword table and defaults.
    pub fn setup_symbols_epsr(&mut self) {
        let s = &mut self.sym;
        for (k, (name, ty, flag, default, text)) in EPSR_SYMBOLS.iter().enumerate() {
            let iref = k + 1;
            s.symbol[iref] = fstr(name, LEN_SYMBOL);
            s.symbol_type[iref] = *ty;
            s.symbol_text[iref] = fstr(text, LEN_VALUE);
            s.iflag[iref] = *flag;
            s.symbol_default[iref] = fstr(default, LEN_VALUE);
        }
        s.iereqfirst = MARK_IEREQFIRST;
        s.iereqlast = MARK_IEREQLAST;
        s.iresetref = MARK_IRESETREF;
        s.iinitref = MARK_IINITREF;
        s.iqminref = MARK_IQMINREF;
        s.iroverlap = MARK_IROVERLAP;
        s.irminfac = MARK_IRMINFAC;
        s.incoeffp = MARK_INCOEFFP;
        s.nsymboluptodata = MARK_NSYMBOLUPTODATA;
        s.nsymboldata = MARK_NSYMBOLDATA;
        s.nsymbol = MARK_NSYMBOL;
        s.nsymbolplus = s.nsymboldata - s.nsymboluptodata;
        for is in 1..=s.nsymbol {
            s.symbol_value[is] = s.symbol_default[is].clone();
        }
        s.dataprompt = "data".into();
        s.extinp = ".EPSR.inp".into();
        s.symbol_setup = "epsr".into();
        s.prompt = "Setup epsr>".into();
        s.nfiles = 1;
    }

    /// `subroutine setup_input_file`: work out the input/run file names and
    /// open the `.EPSR.inp` file.  `name` is the file name given on the command
    /// line (e.g. `200DMSO400Water.EPSR.inp`).  Returns the opened file
    /// positioned after its first record (the run title).
    pub fn setup_input_file(&mut self, name: &str) -> Result<TextFile> {
        let mut fname_inp = name.trim().to_string();
        if !fname_inp.contains('.') {
            fname_inp.push_str(&self.sym.extinp);
        }
        self.sym.fname_inp = fname_inp.clone();
        let full = self.sym.full(&fname_inp);
        self.sym.fnameinp = full.to_string_lossy().into_owned();
        self.log(&format!("\nsetup_input_file> Full filename = {}", self.sym.fnameinp));
        let mut f = TextFile::open(&full)
            .map_err(|e| EpsrError::Io(format!("cannot open input file {}: {e}", self.sym.fnameinp)))?;
        self.log(&format!("\nsetup_input_file> Reading input file: \"{}\"", fname_inp));
        let ilenuptodot = (indexlastdot(&fname_inp) - 1).max(0) as usize;
        let title = {
            let mut r = f.list();
            r.word()?.unwrap_or_default()
        };
        let mut fname_run = title;
        if fname_run != fname_inp[..ilenuptodot] {
            self.log("\nsetup_input_file> Run name in input file is different from filename specified");
            fname_run = fname_inp[..ilenuptodot].to_string();
        }
        self.sym.fname_run = fname_run.clone();
        self.sym.fnamerun = self.sym.full(&fname_run).to_string_lossy().into_owned();
        Ok(f)
    }

    /// `subroutine read_input_file`: process every line of the `.EPSR.inp`
    /// file.  Reading `fnameato`/`fnamepcof` lines triggers reading of the
    /// configuration and potential files (via `check_symbol`).
    pub fn read_input_file(&mut self, f: &mut TextFile) -> Result<i32> {
        let mut ierrout = 0;
        self.pot.npairspcof = 0;
        self.mol.ntype = 0;
        let mut command = String::from("next");
        let mut icount: i32 = 1;
        let mut ifilecount: i32 = 1;
        let mut ierr = 0;
        let mut sometext = String::new();
        while command != "e" && ierr == 0 {
            if icount > self.sym.nsymboldata && ifilecount < self.sym.nfiles {
                ifilecount += 1;
                icount = self.sym.nsymboluptodata + 1;
            } else if icount > self.sym.nsymbol {
                icount = 1;
            }
            match f.read_line() {
                Ok(l) => sometext = l,
                Err(e) => ierr = e.iostat(),
            }
            let w = parse(&sometext, 12);
            if w.nwords() > 0 {
                self.process_commands(&sometext, &w, &mut command, &mut icount, &mut ifilecount, &mut ierrout)?;
            }
        }
        Ok(ierrout)
    }

    /// Process additional input lines (same syntax as the `.EPSR.inp` file:
    /// `keyword value`, `data n`, ...) after the file has been read.
    pub fn apply_input_lines(&mut self, lines: &[String]) -> Result<()> {
        let mut command = String::from("next");
        let mut icount: i32 = 1;
        let mut ifilecount: i32 = 1;
        let mut ierrout = 0;
        for l in lines {
            let w = parse(l, 12);
            if w.nwords() > 0 {
                self.process_commands(l, &w, &mut command, &mut icount, &mut ifilecount, &mut ierrout)?;
                if ierrout != 0 {
                    return Err(EpsrError::Input(format!("input line not recognised: {l:?}")));
                }
            }
        }
        Ok(())
    }

    /// `subroutine process_commands(..., ignore=1)` for non-interactive input.
    fn process_commands(
        &mut self,
        sometext: &str,
        w: &crate::fortran::Words,
        command: &mut String,
        icount: &mut i32,
        ifilecount: &mut i32,
        ierrout: &mut i32,
    ) -> Result<()> {
        *ierrout = 0;
        let nwords = w.nwords();
        *command = if nwords > 0 { fstr(w.word(sometext, 1), LEN_COMMAND) } else { "d".into() };
        if command == "d" || command == "u" {
            // Interactive navigation; only 'd' (accept, go forward) can occur in a file.
            self.sym.symbol_save.clear();
            if command == "d" {
                *icount += 1;
                if *icount > self.sym.nsymboldata && *ifilecount < self.sym.nfiles {
                    *ifilecount += 1;
                    *icount = self.sym.nsymboluptodata + 1;
                } else if *icount > self.sym.nsymbol {
                    *icount = 0;
                    *ifilecount = 1;
                }
            } else {
                *icount -= 1;
                if *icount < 0 {
                    *icount = self.sym.nsymbol;
                    *ifilecount = self.sym.nfiles;
                }
                if *icount <= self.sym.nsymboluptodata && *ifilecount > 1 {
                    *icount = self.sym.nsymboldata;
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
            self.sym.inotrecog = 0;
            let mut is = 0;
            *ifilecount = 1;
            while is < self.sym.nsymboluptodata {
                is += 1;
                self.sym.symbol_value[is] = self.sym.symbol_default[is].clone();
            }
            while *ifilecount < self.sym.nfiles {
                is += 1;
                let iref = is + self.sym.nsymbolplus * (*ifilecount - 1);
                self.sym.symbol_value[iref] = self.sym.symbol_default[is].clone();
                if is >= self.sym.nsymboldata && *ifilecount < self.sym.nfiles {
                    is = self.sym.nsymboluptodata;
                    *ifilecount += 1;
                }
            }
            *command = "next".into();
        } else if command == "ls" || command == "dir" || command == "cd" || command == "pwd" {
            self.sym.inotrecog = 0;
        } else if *command == self.sym.dataprompt {
            self.sym.inotrecog = 0;
            let mut ndata = 1;
            if nwords > 1 {
                if let Ok(n) = read_int(w.word(sometext, 2)) {
                    ndata = n;
                }
            }
            if ndata < 1 {
                ndata = 1;
            }
            if ndata > self.sym.nfiles {
                ndata = self.sym.nfiles;
            }
            *ifilecount = ndata;
            *icount = self.sym.nsymboluptodata + 1;
            *command = "next".into();
        } else if command == "help" || command == "?" || command == "save" || command == "forcesave" {
            self.sym.inotrecog = 0;
        } else if command == "q" {
            self.sym.inotrecog = 0;
            self.sym.isave = 1;
            *command = "e".into();
        } else if command != "e" && command != "next" && command != "n" {
            let ifound = self.sym.search_for_symbol(command, *icount + 1);
            if ifound > 0 {
                let iref = self.sym.ireff(ifound, *ifilecount);
                self.sym.symbol_save = self.sym.symbol_value[iref].clone();
                if nwords > 1 {
                    let (v, _) = concatenate(2, nwords, w, sometext, LEN_VALUE);
                    self.sym.symbol_value[iref] = v;
                    *icount = ifound;
                    *command = "store".into();
                } else {
                    *icount = ifound;
                    *command = "next".into();
                }
                if ifound > self.sym.nsymboldata {
                    *ifilecount = self.sym.nfiles;
                }
                if ifound < self.sym.nsymboluptodata + 1 {
                    *ifilecount = 1;
                }
            } else {
                self.log(&format!(
                    "\nprocess_commands> Command not recognised: {}:  {}",
                    crate::fortran::format::pad_to(command, LEN_COMMAND),
                    substr(sometext, 1, 40)
                ));
                self.sym.inotrecog += 1;
                *ierrout = 1;
                *command = "next".into();
            }
        }
        if command == "store" {
            self.sym.inotrecog = 0;
            self.sym.isave = 0;
            let ierr = self.check_symbol(icount, *ifilecount)?;
            *command = if ierr == 0 { "next".into() } else { "n".into() };
        }
        if command == "n" {
            self.sym.inotrecog = 0;
            if !self.sym.symbol_save.is_empty() && *icount > 0 {
                let iref = self.sym.ireff(*icount, *ifilecount);
                self.sym.symbol_value[iref] = self.sym.symbol_save.clone();
            } else {
                *icount += 1;
            }
            *command = "next".into();
        }
        Ok(())
    }

    /// `subroutine check_symbol` (non-interactive branches used by `epsr`).
    /// Returns the Fortran `ierr` (0 = value accepted).
    fn check_symbol(&mut self, icount: &mut i32, ifilecount: i32) -> Result<i32> {
        if *icount < 1 {
            *icount = 1;
        }
        if *icount > self.sym.nsymbol {
            *icount = self.sym.nsymbol;
        }
        let ic = *icount;
        let iref = self.sym.ireff(ic, ifilecount);
        let mut ierr = 0;
        // Type check of the value(s)
        if self.sym.symbol_type[ic] != 'c' {
            let sometext = self.sym.symbol_value[iref].clone();
            let w = parse(&sometext, 11);
            let mut iw = 0;
            while iw < w.nwords() && ierr == 0 {
                iw += 1;
                let tok = w.word(&sometext, iw);
                let ok = if self.sym.symbol_type[ic] == 'i' { read_int(tok).is_ok() } else { read_real(tok).is_ok() };
                if !ok {
                    ierr = 1;
                    self.log(&format!(
                        "\n{} Value {} is not {} - try again!",
                        self.sym.prompt,
                        tok,
                        if self.sym.symbol_type[ic] == 'i' { "an integer" } else { "a real number" }
                    ));
                }
            }
            if ierr != 0 {
                return Ok(ierr);
            }
        }
        // Change of the number of datasets: resequence the per-dataset values.
        if ic == self.sym.nsymboluptodata {
            if let Ok(nfilet) = read_int(&self.sym.symbol_value[iref]) {
                let s = &mut self.sym;
                let nfadd = (nfilet - s.nfiles) * s.nsymbolplus;
                if nfadd > 0 {
                    let mut isymbol = s.nsymbol + 1;
                    let mut ivalue = s.ireff(s.nsymbol, s.nfiles) + 1;
                    while isymbol > s.nsymboldata + 1 {
                        isymbol -= 1;
                        ivalue -= 1;
                        let ivaluenew = ivalue + nfadd;
                        if ivaluenew as usize <= crate::state::MSYMBOL {
                            s.symbol_value[ivaluenew] = s.symbol_value[ivalue].clone();
                        } else {
                            return Err(EpsrError::Input("tried to generate too many symbol values".into()));
                        }
                    }
                } else if nfadd < 0 {
                    let mut isymbol = s.nsymboldata;
                    let mut ivalue = s.ireff(s.nsymboldata, s.nfiles);
                    while isymbol < s.nsymbol {
                        isymbol += 1;
                        ivalue += 1;
                        let ivaluenew = ivalue + nfadd;
                        s.symbol_value[ivaluenew] = s.symbol_value[ivalue].clone();
                    }
                }
                let mut ifiletemp = s.nfiles;
                s.nfiles = nfilet;
                if nfadd > 0 {
                    while ifiletemp < s.nfiles {
                        ifiletemp += 1;
                        let mut isymbol = s.nsymboluptodata;
                        while isymbol < s.nsymboldata {
                            isymbol += 1;
                            let r = s.ireff(isymbol, ifiletemp);
                            s.symbol_value[r] = s.symbol_default[isymbol].clone();
                        }
                    }
                }
            }
        }
        if self.sym.symbol[ic] == "fnameato" {
            self.read_ato(ic, ifilecount)?;
        }
        if self.sym.symbol[ic] == "fnamepcof" {
            self.read_pcof(ic, ifilecount)?;
        }
        if self.sym.symbol[ic].starts_with("atom") {
            self.check_atom_symbol(ic, iref);
        }
        // Program-specific checks (sharm, sdf, nwts, xwts, mapgr, changeato, ...)
        let setup = self.sym.symbol_setup.clone();
        let e = crate::tools::check_symbol_extra(self, &setup, ic, ifilecount, iref)?;
        if e != 0 {
            return Ok(e);
        }
        if setup == "epsr" {
            // a reset forces re-initialisation
            let irefiinit = self.sym.ireff(self.sym.iinitref, ifilecount);
            let irefireset = self.sym.ireff(self.sym.iresetref, ifilecount);
            if self.sym.symbol_value[irefireset] != "0" {
                self.sym.symbol_value[irefiinit] = "1".into();
            }
            // Ensure the number of potential coefficients does not exceed the allowed dimension
            if self.sym.symbol[ic] == "ncoeffp" {
                let irefncoeffp = self.sym.ireff(self.sym.incoeffp, ifilecount);
                match read_int(&self.sym.symbol_value[irefncoeffp]) {
                    Ok(n) => {
                        self.pot.ncoeffp = n;
                        if n as usize > crate::state::MCOEFF {
                            self.log(&format!(
                                "\ncheck_symbol> Too many potential coefficients defined. Will reduce to maximum allowed ({})",
                                fmt("(i3)", &[FV::I(crate::state::MCOEFF as i64)])
                            ));
                            self.pot.ncoeffp = crate::state::MCOEFF as i32;
                            let mut lw = crate::fortran::format::ListWriter::new();
                            lw.int(self.pot.ncoeffp);
                            let txt = lw.finish();
                            let w = parse(&txt, 5);
                            let (v, _) = concatenate(1, w.nwords(), &w, &txt, LEN_VALUE);
                            self.sym.symbol_value[irefncoeffp] = v;
                        }
                    }
                    // NB: the read error code is returned as ierr of check_symbol
                    Err(_) => return Ok(1),
                }
            }
        }
        Ok(ierr)
    }

    /// ATOM CHECK of `check_symbol`: a symbol whose name starts with `atom` holds
    /// a list of atom types (`*` matches any suffix).  If a type is not present
    /// in the current `.ato` file the previous value is restored.
    fn check_atom_symbol(&mut self, _ic: i32, iref: i32) {
        let sometext = self.sym.symbol_value[iref].clone();
        let w = parse(&sometext, 5);
        let mut iw = 0;
        while iw < w.nwords() {
            iw += 1;
            let word = w.word(&sometext, iw);
            let pat: &str = match word.find('*') {
                Some(p) => &word[..p],
                None => word,
            };
            let ncsearch = pat.len();
            let mut ifound = 0;
            let mut it = 0;
            while it < self.mol.ntype && ifound == 0 {
                it += 1;
                if ncsearch > 0 {
                    // type(it)(1:ncsearch) — blank padded to the search length
                    let t = crate::fortran::format::pad_to(&self.mol.type_[it], 3);
                    let t = crate::fortran::format::pad_to(&t, ncsearch);
                    if t == pat {
                        ifound = it;
                    }
                } else {
                    ifound = it;
                }
            }
            if ifound == 0 {
                self.log(&format!(
                    "\n{} WARNING! Atom type {} not found in .ato file {}",
                    self.sym.prompt, word, self.sym.fname_ato
                ));
                self.sym.symbol_value[iref] = self.sym.symbol_save.clone();
            }
        }
    }

    /// `setup_symbols_<tool>` for an auxiliary program: install its keyword
    /// table (generated in `symtab_tools.rs`) and the input-file extension.
    pub fn install_tool_table(&mut self, t: &ToolTable) {
        let s = &mut self.sym;
        for (k, (name, ty, flag, default, text)) in t.symbols.iter().enumerate() {
            let iref = k + 1;
            s.symbol[iref] = fstr(name, LEN_SYMBOL);
            s.symbol_type[iref] = *ty;
            s.symbol_text[iref] = fstr(text, LEN_VALUE);
            s.iflag[iref] = *flag;
            s.symbol_default[iref] = fstr(default, LEN_VALUE);
        }
        let n = t.symbols.len() as i32;
        s.nsymboluptodata = n;
        s.nsymboldata = n;
        s.nsymbol = n;
        for (name, v) in t.marks {
            match *name {
                "nsymboluptodata" => s.nsymboluptodata = *v,
                "nsymboldata" => s.nsymboldata = *v,
                "nsymbol" => s.nsymbol = *v,
                "irefrho" => s.irefrho = *v,
                _ => {}
            }
        }
        s.nsymbolplus = s.nsymboldata - s.nsymboluptodata;
        for is in 1..=s.nsymbol {
            s.symbol_value[is] = s.symbol_default[is].clone();
        }
        // read(symbol_value(nsymboluptodata),*,iostat=ierr) nfiles
        let nf = if s.nsymboluptodata >= 1 && s.nsymboluptodata <= s.nsymbol { read_int(&s.symbol_value[s.nsymboluptodata]).ok() } else { None };
        s.nfiles = match nf {
            Some(v) if s.nsymboluptodata < s.nsymboldata => v,
            _ => 1,
        };
        if t.notes.iter().any(|n| *n == "nfiles=1") {
            s.nfiles = 1;
        }
        if t.notes.iter().any(|n| *n == "nfiles=0") {
            s.nfiles = 0;
        }
        s.prompt = t.prompt.to_string();
        s.dataprompt = t.dataprompt.to_string();
        s.extinp = t.extinp.to_string();
        s.symbol_setup = t.name.to_string();
    }

    /// The input-file part of `save_files_others`: rewrite the `.EPSR.inp`
    /// (or other program's `.dat`) file with the current values (EPSR updates
    /// e.g. `nsumt`, `rho`, `potfac`).
    pub fn save_inp_file(&mut self) -> Result<()> {
        let s = &self.sym;
        let iqmin = if s.symbol_setup == "epsr" { s.ireff(s.iqminref, s.nfiles) } else { s.nsymbol };
        let mut out = OutFile::create(std::path::Path::new(&s.fnameinp));
        let ilen = (indexlastdot(&s.fname_inp) - 1).max(0) as usize;
        let fname_run = s.fname_inp[..ilen].to_string();
        out.line(&fmt("(a,15x,'Title of this file')", &[FV::S(&fname_run)]));
        let mut icount = 0;
        let mut ifilecount = 1;
        let is_epsr = s.symbol_setup == "epsr";
        // ncf(1), ncl(1) left over from the previous `parse` call: the Fortran
        // uses them when a value is blank (parse then finds no words).
        let mut last1: (i32, i32) = (1, 1);
        while icount < s.nsymboldata {
            icount += 1;
            if is_epsr {
                if icount > s.iereqfirst && icount <= s.iereqlast {
                    icount = s.iereqlast + 1;
                } else if icount > iqmin && icount < s.nsymboluptodata {
                    icount = s.nsymboluptodata;
                }
            }
            let iref = s.ireff(icount, ifilecount);
            if icount == s.nsymboluptodata + 1 {
                out.line(&fmt("(/a,1x,i3/)", &[FV::S(&s.dataprompt), FV::I(ifilecount as i64)]));
            }
            let sometext = &s.symbol_value[iref];
            let w = parse(sometext, 5);
            let val = if w.nwords() > 0 {
                last1 = (w.ncf[0], w.ncl[0]);
                w.span(sometext, 1, w.nwords())
            } else {
                substr(sometext, last1.0, last1.1)
            };
            out.line(&fmt(
                "(a11,1x,a,15x,a)",
                &[FV::Ch(&s.symbol[icount], LEN_SYMBOL), FV::S(&val), FV::Ch(&s.symbol_text[icount], LEN_VALUE)],
            ));
            if icount >= s.nsymboldata && ifilecount < s.nfiles {
                ifilecount += 1;
                icount = s.nsymboluptodata;
            }
        }
        out.line("q");
        out.close()?;
        Ok(())
    }

    /// The `.pcof` part of `save_files_others` (header keywords + pair potentials).
    pub fn save_pcof_file(&mut self) -> Result<()> {
        let s = &self.sym;
        let iqmin = s.ireff(s.iqminref, s.nfiles);
        let mut out = OutFile::create(std::path::Path::new(&s.fnamepcof));
        let ifilecount = s.nfiles; // value left by the loop over datasets above
        let mut icount = s.iereqfirst + 1;
        while icount < s.nsymboluptodata {
            let iref = s.ireff(icount, ifilecount);
            let sometext = &s.symbol_value[iref];
            let w = parse(sometext, 5);
            let val = w.span(sometext, 1, w.nwords().max(1));
            out.line(&fmt(
                "(a11,1x,a,15x,a)",
                &[FV::Ch(&s.symbol[icount], LEN_SYMBOL), FV::S(&val), FV::Ch(&s.symbol_text[icount], LEN_VALUE)],
            ));
            if icount == s.iereqlast {
                icount = iqmin;
            }
            icount += 1;
        }
        out.line("q");
        self.write_pcof_pairs(&mut out)?;
        self.write_pcof_threebody(&mut out)?;
        self.write_pcof_coordcons(&mut out);
        out.close()?;
        Ok(())
    }
}

/// A keyword table of an auxiliary EPSR program (see `symtab_tools.rs`).
pub struct ToolTable {
    /// Program name (`symbol_setup`).
    pub name: &'static str,
    /// Input-file extension, e.g. `.PARTIALS.dat`.
    pub extinp: &'static str,
    pub prompt: &'static str,
    pub dataprompt: &'static str,
    pub dataprompttext: &'static str,
    /// (symbol, type, iflag, default, description)
    pub symbols: &'static [(&'static str, char, i32, &'static str, &'static str)],
    /// Markers `name = iref` in the Fortran (number of entries at that point).
    pub marks: &'static [(&'static str, i32)],
    /// Non-table Fortran statements of the routine (handled by the tool code).
    pub notes: &'static [&'static str],
}

/// Helper used by several readers: list-directed read of all words of a value.
pub fn read_reals(s: &str, n: usize) -> Vec<f32> {
    let mut r = ListRead::from_str(s);
    let mut v = Vec::new();
    for _ in 0..n {
        match r.real() {
            Ok(Some(x)) => v.push(x),
            _ => break,
        }
    }
    v
}
