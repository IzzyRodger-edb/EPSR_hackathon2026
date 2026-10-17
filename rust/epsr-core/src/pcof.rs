//! The potential-coefficients file (`.pcof`).
//!
//! Translation of `read_pcof` and of the `.pcof` part of `save_files_others`
//! (`read_and_save_files.f`).
//!
//! ## Layout
//! ```text
//!  keyword  value(s)        description      (ereq ... refpotfac: EP control parameters)
//!  ...
//!  q
//!  npairs
//!  for each pair:   label1 label2  index  coulfac
//!                   rminex  pcoefex  coordex  [rmaxex pcoefexmax] [Gaussians: r w h ...]
//!                   pcoefp(1..ncoeffp)   (pcoefp(1) = refinement factor, then the EP coefficients)
//! ```
//! Each pair becomes an extra symbol `"L1-L2"` in the symbol table whose
//! value is `"rminex coordex ..."`.

use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::{read_int, read_real, TextFile};
use crate::fortran::output::OutFile;
use crate::fortran::text::{concatenate, fstr, indexlastdot, parse};
use crate::state::State;
use crate::symbols::{LEN_SYMBOL, LEN_VALUE};

impl State {
    /// `subroutine read_pcof`: read the `.pcof` file named by symbol `icount`.
    pub fn read_pcof(&mut self, icount: i32, ifilecount: i32) -> Result<()> {
        let ifound = self.sym.ireff(icount, ifilecount);
        let (irefireset, irefiinit) = (self.sym.ireff(self.sym.iresetref, ifilecount), self.sym.ireff(self.sym.iinitref, ifilecount));
        let ireset = read_int(&self.sym.symbol_value[irefireset])?;
        let iinit = read_int(&self.sym.symbol_value[irefiinit])?;
        let v = self.sym.symbol_value[ifound].clone();
        if v == "search" || v == "<undefined>" {
            return Err(EpsrError::Input("no .pcof file specified (fnamepcof)".into()));
        }
        let w = parse(&v, 5);
        let (mut fname_pcof, _) = concatenate(1, w.nwords(), &w, &v, 256);
        let ilenlastdot = indexlastdot(&fname_pcof) - 1;
        if ilenlastdot > 0 {
            fname_pcof = format!("{}.pcof", &fname_pcof[..ilenlastdot as usize]);
        } else {
            fname_pcof.push_str(".pcof");
        }
        self.sym.fname_pcof = fname_pcof.clone();
        self.sym.fnamepcof = self.sym.full(&fname_pcof).to_string_lossy().into_owned();
        let mut ireset = ireset;
        let mut iinit = iinit;
        let mut f = match TextFile::open(std::path::Path::new(&self.sym.fnamepcof)) {
            Ok(f) => f,
            Err(_) => {
                // The .pcof file does not exist: use defaults or existing values.
                self.log(&format!("\n{} Problems with specified .pcof file: {}\n - will use defaults or existing values", self.sym.prompt, fname_pcof));
                let iref = self.sym.ireff(self.sym.incoeffp, ifilecount);
                if let Ok(n) = read_int(&self.sym.symbol_value[iref]) {
                    self.pot.ncoeffp = n;
                }
                if self.pot.ncoeffp as usize > crate::state::MCOEFF {
                    self.pot.ncoeffp = crate::state::MCOEFF as i32;
                }
                let iref = self.sym.ireff(self.sym.iroverlap, ifilecount);
                if let Ok(x) = read_real(&self.sym.symbol_value[iref]) {
                    self.inp.roverlap = x;
                }
                let iref = self.sym.ireff(self.sym.irminfac, ifilecount);
                if let Ok(x) = read_real(&self.sym.symbol_value[iref]) {
                    self.inp.rminfac = x;
                }
                let mut ierrread = 0;
                if self.pot.npairspcof <= 0 {
                    // force initialisation of the .pcof file and the accumulators
                    ireset = 1;
                    iinit = 1;
                    ierrread = 1;
                }
                if ierrread != 0 || ireset == 1 {
                    self.default_pair_potentials();
                }
                self.sym.symbol_value[ifound] = fname_pcof;
                self.sym.symbol_value[irefireset] = ireset.to_string();
                self.sym.symbol_value[irefiinit] = iinit.to_string();
                return Ok(());
            }
        };
        let mut ierrread = 0;
        // Initialise the three-body potential
        self.tb.nthreeb = 0;
        self.tb.nthreebangsmax = 0;
        self.cc.ncoordconstraints = 0;
        // ---- header keywords
        self.sym.nsymbol = self.sym.nsymboldata;
        let mut icount = self.sym.iereqfirst;
        let mut command = String::from(" ");
        let mut ierr = 0;
        while command != "q" && ierr == 0 {
            icount += 1;
            let mut nwords = 0;
            let mut sometext = String::new();
            let mut w = parse("", 11);
            while nwords == 0 && ierr == 0 {
                match f.read_line() {
                    Ok(l) => sometext = l,
                    Err(e) => {
                        ierr = e.iostat();
                        ierrread += ierr.abs();
                    }
                }
                w = parse(&sometext, 11);
                nwords = w.nwords();
            }
            if nwords == 0 {
                break;
            }
            command = fstr(w.word(&sometext, 1), 20);
            if nwords > 1 {
                icount = self.sym.search_for_symbol(&command, icount);
                if icount != 0 && icount <= self.sym.nsymbol {
                    let iref = self.sym.ireff(icount, ifilecount);
                    let (val, _) = concatenate(2, nwords, &w, &sometext, LEN_VALUE);
                    self.sym.symbol_value[iref] = val;
                }
                if command == "ncoeffp" {
                    self.sym.incoeffp = icount;
                }
                if command == "roverlap" {
                    self.sym.iroverlap = icount;
                }
                if command == "rminfac" {
                    self.sym.irminfac = icount;
                }
            }
        }
        let iref = self.sym.ireff(self.sym.incoeffp, ifilecount);
        match read_int(&self.sym.symbol_value[iref]) {
            Ok(n) => self.pot.ncoeffp = n,
            Err(_) => ierrread += 1,
        }
        if self.pot.ncoeffp as usize > crate::state::MCOEFF {
            self.pot.ncoeffp = crate::state::MCOEFF as i32;
        }
        self.sym.symbol_value[iref] = self.pot.ncoeffp.to_string();
        let iref = self.sym.ireff(self.sym.iroverlap, ifilecount);
        match read_real(&self.sym.symbol_value[iref]) {
            Ok(x) => self.inp.roverlap = x,
            Err(_) => ierrread += 1,
        }
        let iref = self.sym.ireff(self.sym.irminfac, ifilecount);
        match read_real(&self.sym.symbol_value[iref]) {
            Ok(x) => self.inp.rminfac = x,
            Err(_) => ierrread += 1,
        }
        // ---- pair potentials
        self.pot.npairspcof = 0;
        {
            let mut r = f.list();
            match r.int() {
                Ok(Some(n)) => self.pot.npairspcof = n,
                Ok(None) => {}
                Err(_) => ierrread += 1,
            }
        }
        self.log(&format!("\n.pcof file \"{}\"", fname_pcof));
        let ncoeffp = self.pot.ncoeffp;
        let mut ierr_last = 0;
        for irefp in 1..=self.pot.npairspcof {
            self.pot.efacp[irefp] = 1.0;
            let sometext = f.read_line()?;
            let w = parse(&sometext, 11);
            self.pot.tlab1[irefp] = fstr(w.word(&sometext, 1), 3);
            self.pot.tlab2[irefp] = fstr(w.word(&sometext, 2), 3);
            self.pot.coulfac[irefp] = 1.0;
            if w.nwords() > 3 && ireset != 1 {
                if let Ok(x) = read_real(w.word(&sometext, 4)) {
                    self.pot.coulfac[irefp] = x;
                }
            }
            self.sym.nsymbol += 1;
            let ns = self.sym.nsymbol;
            let iref = self.sym.ireff(ns, self.sym.nfiles);
            self.sym.symbol[ns] = fstr(&format!("{}-{}", w.word(&sometext, 1), w.word(&sometext, 2)), LEN_SYMBOL);
            self.sym.iflag[ns] = 0;
            self.sym.symbol_default[ns] = String::new();
            self.sym.symbol_type[ns] = 'r';
            let sometext = f.read_line()?;
            let w = parse(&sometext, 11);
            let p = &mut self.pot;
            let mut iword = 1;
            let mut ierr = 0;
            match read_real(w.word(&sometext, iword)) {
                Ok(x) => p.rminex[irefp] = x,
                Err(_) => ierr = 1,
            }
            iword += 1;
            match read_real(w.word(&sometext, iword)) {
                Ok(x) => p.pcoefex[irefp] = x,
                Err(_) => ierr = 1,
            }
            if p.rminex[irefp] == 0.0 || ireset == 2 {
                p.pcoefex[irefp] = 0.0;
            }
            if iword < w.nwords() {
                iword += 1;
                match read_real(w.word(&sometext, iword)) {
                    Ok(x) => p.coordex[irefp] = x,
                    Err(_) => ierr = 1,
                }
            } else {
                p.coordex[irefp] = 0.0;
            }
            p.rmaxex[irefp] = 0.0;
            p.pcoefexmax[irefp] = 0.0;
            let npgausstest = (w.nwords() - iword) / 3;
            let iwordresidual = w.nwords() - iword - npgausstest * 3;
            if iwordresidual == 2 {
                iword += 1;
                match read_real(w.word(&sometext, iword)) {
                    Ok(x) => p.rmaxex[irefp] = x,
                    Err(_) => ierr = 1,
                }
                iword += 1;
                match read_real(w.word(&sometext, iword)) {
                    Ok(x) => p.pcoefexmax[irefp] = x,
                    Err(_) => ierr = 1,
                }
                if p.rmaxex[irefp] == 0.0 {
                    p.pcoefexmax[irefp] = 0.0;
                }
            }
            // additional Gaussian (or modified Morse) potentials: (r, w, h) triples
            let mut igauss: i32 = 0;
            while iword < w.nwords() && ierr == 0 && (igauss as usize) < crate::state::MPGAUSS {
                igauss += 1;
                let (i1, i2, i3) = (iword + 1, iword + 2, iword + 3);
                let mut rd = |k: i32, dst: &mut f32| {
                    if k <= w.nwords() {
                        match read_real(w.word(&sometext, k)) {
                            Ok(x) => *dst = x,
                            Err(_) => ierr = 1,
                        }
                    } else {
                        ierr = -1;
                    }
                };
                rd(i1, &mut p.rpgauss[(igauss, irefp)]);
                rd(i2, &mut p.wpgauss[(igauss, irefp)]);
                rd(i3, &mut p.hpgauss[(igauss, irefp)]);
                if ierr != 0 {
                    igauss -= 1;
                }
                iword = i3;
            }
            p.npgauss[irefp] = igauss;
            if p.pcoefex[irefp].abs() < 1.0e-10 {
                p.pcoefex[irefp] = 0.0;
            }
            ierrread += ierr;
            // NB: the Fortran tests the first characters of the first two words of
            // `sometext`, which at this point holds the *numbers* line (not the
            // labels), so the q/Q exemption never applies.  Reproduced as is.
            let c1 = w.word(&sometext, 1).chars().next().unwrap_or(' ');
            let c2 = if w.nwords() > 1 { w.word(&sometext, 2).chars().next().unwrap_or(' ') } else { ' ' };
            if self.pot.rminex[irefp].abs() < self.inp.roverlap && !matches!(c1, 'q' | 'Q') && !matches!(c2, 'q' | 'Q') {
                self.pot.rminex[irefp] = crate::fortran::fmath::sign(self.inp.roverlap, self.pot.rminex[irefp]);
            }
            let withmax = self.pot.rmaxex[irefp] > 0.0;
            let txt = if self.pot.npgauss[irefp] == 0 {
                let p = &self.pot;
                if withmax {
                    fmt("((1x,f10.5,a,f10.5,1x,f10.5))", &[FV::R(p.rminex[irefp]), FV::S(" "), FV::R(p.coordex[irefp]), FV::R(p.rmaxex[irefp])])
                } else {
                    fmt("(1x,f10.5,a,f10.5)", &[FV::R(p.rminex[irefp]), FV::S(" "), FV::R(p.coordex[irefp])])
                }
            } else {
                let p = &self.pot;
                let formatline = if withmax {
                    format!("((1x,f10.5,a,f10.5,1x,f10.5),{}(1x,f9.4,1x,f9.4,1x,e11.4))", p.npgauss[irefp])
                } else {
                    format!("((1x,f10.5,a,f10.5),{}(1x,f9.4,1x,f9.4,1x,e11.4))", p.npgauss[irefp])
                };
                let mut v = vec![FV::R(p.rminex[irefp]), FV::S(" "), FV::R(p.coordex[irefp])];
                if withmax {
                    v.push(FV::R(p.rmaxex[irefp]));
                }
                for ig in 1..=p.npgauss[irefp] {
                    v.push(FV::R(p.rpgauss[(ig, irefp)]));
                    v.push(FV::R(p.wpgauss[(ig, irefp)]));
                    v.push(FV::R(p.hpgauss[(ig, irefp)]));
                }
                fmt(&formatline, &v)
            };
            let w2 = parse(&txt, 11);
            let (val, _) = concatenate(1, w2.nwords(), &w2, &txt, LEN_VALUE);
            self.sym.symbol_value[iref] = val;
            self.sym.symbol_text[ns] = "Minimum (and maximum) separation and P,W,H of Gaussians for this atom pair".into();
            // coefficients
            {
                let p = &mut self.pot;
                let mut r = f.list();
                ierr_last = 0;
                for ic in 1..=ncoeffp {
                    if r.get_r(&mut p.pcoefp[(ic, irefp)]).is_err() {
                        ierrread += 1;
                        ierr_last = 1;
                        break;
                    }
                }
            }
            let p = &mut self.pot;
            if ireset != 0 {
                if ireset < 2 {
                    p.pcoefp[(1, irefp)] = 1.0;
                }
                for ic in 2..=ncoeffp {
                    p.pcoefp[(ic, irefp)] = 0.0;
                }
            }
            p.efacp[irefp] = p.pcoefp[(1, irefp)];
            // a negative Gaussian position or width holds the EP of this pair fixed
            for ig in 1..=p.npgauss[irefp] {
                if p.rpgauss[(ig, irefp)] < 0.0 || p.wpgauss[(ig, irefp)] < 0.0 {
                    p.efacp[irefp] = 0.0;
                }
            }
            p.pcoefp[(1, irefp)] = p.efacp[irefp];
        }
        // ---- match pairs of the .ato file with pairs of the .pcof file
        for ip in 1..=self.mol.npairs {
            let it = self.mol.n1extr[ip];
            let jt = self.mol.n2extr[ip];
            let mut irefpfound = 0;
            let mut irefp = 0;
            while irefp < self.pot.npairspcof && irefpfound == 0 {
                irefp += 1;
                let (a, b) = (&self.pot.tlab1[irefp], &self.pot.tlab2[irefp]);
                if (*a == self.mol.type_[it] && *b == self.mol.type_[jt]) || (*a == self.mol.type_[jt] && *b == self.mol.type_[it]) {
                    irefpfound = irefp;
                }
            }
            if irefpfound > 0 {
                self.pot.irefpcof[ip] = irefpfound;
                self.log(&format!(
                    "read_pcof> Pair found: {:>5} {:>5} {} {}",
                    ip, irefpfound, self.pot.tlab1[irefpfound], self.pot.tlab2[irefpfound]
                ));
            } else {
                // create a new set of coefficients for this pair
                self.pot.npairspcof += 1;
                let n = self.pot.npairspcof;
                self.pot.irefpcof[ip] = n;
                let (ti, tj) = (self.mol.type_[it].clone(), self.mol.type_[jt].clone());
                self.pot.tlab1[n] = ti.clone();
                self.pot.tlab2[n] = tj.clone();
                let p = &mut self.pot;
                p.pcoefp[(1, n)] = 1.0;
                for ic in 2..=p.ncoeffp {
                    p.pcoefp[(ic, n)] = 0.0;
                }
                let is_q = |t: &str| matches!(t.chars().next(), Some('q') | Some('Q'));
                if !is_q(&ti) && !is_q(&tj) {
                    p.rminex[n] = self.inp.rminfac * self.mol.dcor[ip];
                    if p.rminex[n] < self.inp.roverlap {
                        p.rminex[n] = self.inp.roverlap;
                    }
                } else {
                    p.rminex[n] = 0.0;
                }
                p.pcoefex[n] = 0.0;
                p.rmaxex[n] = 0.0;
                p.pcoefexmax[n] = 0.0;
                p.npgauss[n] = 0;
                // NB: coulfac, coordex and efacp of the new pair are not set by the
                // Fortran (COMMON data, i.e. 0): the pair has no Coulomb term and its
                // EP is not refined in this run.  Reproduced.
                let mut lw = ListWriter::new();
                lw.chars(&ti).chars(" ").chars(&tj);
                let txt = lw.finish();
                let wd = parse(&txt, 11);
                self.sym.nsymbol += 1;
                let ns = self.sym.nsymbol;
                let iref = self.sym.ireff(ns, self.sym.nfiles);
                self.sym.symbol[ns] = fstr(&format!("{}-{}", wd.word(&txt, 1), wd.word(&txt, 2)), LEN_SYMBOL);
                self.sym.iflag[ns] = 0;
                self.sym.symbol_default[ns] = String::new();
                self.sym.symbol_type[ns] = 'r';
                let mut lw = ListWriter::new();
                lw.real(self.pot.rminex[n]);
                let txt = lw.finish();
                let wd = parse(&txt, 11);
                let (val, _) = concatenate(1, 1, &wd, &txt, LEN_VALUE);
                self.sym.symbol_value[iref] = val;
                self.sym.symbol_text[ns] = "Minimum separation and P,W,H of Gaussians for this atom pair".into();
                self.log(&format!(" {} {}", self.sym.symbol[ns], self.sym.symbol_value[iref]));
            }
        }
        // ---- optional three-body / coordination-constraint sections
        // (`do while(ierr.eq.0)`: ierr is the status of the last coefficient read)
        let mut ierr = ierr_last;
        while ierr == 0 {
            match f.read_line() {
                Ok(l) => {
                    let w = parse(&l, 11);
                    if w.nwords() > 0 {
                        let w1 = w.word(&l, 1);
                        if w1 == "tb" && w.nwords() > 13 {
                            self.read_pcof_threebody(&l, &w)?;
                        } else if w1 == "coordcon" && w.nwords() > 8 {
                            self.read_coordcon_line(&l, &w)?;
                        }
                    }
                }
                Err(_) => ierr = 1,
            }
        }
        // if any error occurred during input, ireset (and iinit) are set to 1
        if ierrread != 0 {
            ireset = 1;
            iinit = 1;
        }
        if ierrread != 0 || ireset == 1 {
            self.default_pair_potentials();
        }
        self.sym.symbol_value[ifound] = fname_pcof;
        self.sym.symbol_value[irefireset] = ireset.to_string();
        self.sym.symbol_value[irefiinit] = iinit.to_string();
        Ok(())
    }

    /// Default pair potentials (the `ierrread.ne.0.or.ireset.eq.1` block of
    /// `read_pcof`): one entry per atom pair of the `.ato` file, minimum
    /// distance `rminfac*dcor` (at least `roverlap`), zero EP coefficients.
    pub(crate) fn default_pair_potentials(&mut self) {
        self.sym.nsymbol = self.sym.nsymboldata;
        self.pot.npairspcof = self.mol.npairs;
        for ip in 1..=self.mol.npairs {
            let it = self.mol.n1extr[ip];
            let jt = self.mol.n2extr[ip];
            let is_q = |t: &str| matches!(t.chars().next(), Some('q') | Some('Q'));
            let p = &mut self.pot;
            if !is_q(&self.mol.type_[it]) && !is_q(&self.mol.type_[jt]) {
                p.rminex[ip] = self.inp.rminfac * self.mol.dcor[ip];
                if p.rminex[ip] < self.inp.roverlap {
                    p.rminex[ip] = self.inp.roverlap;
                }
            } else {
                p.rminex[ip] = 0.0;
            }
            p.pcoefex[ip] = 0.0;
            p.coordex[ip] = 0.0;
            p.rmaxex[ip] = 0.0;
            p.pcoefexmax[ip] = 0.0;
            p.irefpcof[ip] = ip;
            p.coulfac[ip] = 1.0;
            p.tlab1[ip] = self.mol.type_[it].clone();
            p.tlab2[ip] = self.mol.type_[jt].clone();
            p.npgauss[ip] = 0;
            self.sym.nsymbol += 1;
            let ns = self.sym.nsymbol;
            let iref = self.sym.ireff(ns, self.sym.nfiles);
            self.sym.symbol[ns] = fstr(&format!("{}-{}", p.tlab1[ip], p.tlab2[ip]), LEN_SYMBOL);
            self.sym.iflag[ns] = 0;
            self.sym.symbol_default[ns] = String::new();
            // write(sometext,*) rminex(ip),' ',coordex(ip); only the first word is kept
            let mut w = ListWriter::new();
            w.real(p.rminex[ip]).chars(" ").real(p.coordex[ip]);
            let txt = w.finish();
            let wd = parse(&txt, 11);
            let (val, _) = concatenate(1, 1, &wd, &txt, LEN_VALUE);
            self.sym.symbol_value[iref] = val;
            self.sym.symbol_text[ns] = "Minimum (coordination number and maximum) separation and P,W,H of Gaussians for this atom pair".into();
            p.pcoefp[(1, ip)] = 1.0;
            for ic in 2..=p.ncoeffp {
                p.pcoefp[(ic, ip)] = 0.0;
            }
        }
    }

    /// Pair section of the `.pcof` file written by `save_files_others`.
    pub(crate) fn write_pcof_pairs(&mut self, out: &mut OutFile) -> Result<()> {
        let mut w = ListWriter::new();
        w.int(self.pot.npairspcof);
        out.line(&w.finish());
        let mut icount = self.sym.nsymboldata;
        for ip in 1..=self.pot.npairspcof {
            out.line(&fmt(
                "(a3,1x,a3,1x,5x,i0,2x,f10.5)",
                &[FV::Ch(&self.pot.tlab1[ip], 3), FV::Ch(&self.pot.tlab2[ip], 3), FV::I(ip as i64), FV::R(self.pot.coulfac[ip])],
            ));
            icount += 1;
            let iref = self.sym.ireff(icount, self.sym.nfiles);
            let sometext = self.sym.symbol_value[iref].clone();
            let wd = parse(&sometext, 11);
            let p = &mut self.pot;
            let mut iword = 1;
            if let Ok(x) = read_real(wd.word(&sometext, iword)) {
                p.rminex[ip] = x;
            }
            if iword < wd.nwords() {
                iword += 1;
                if let Ok(x) = read_real(wd.word(&sometext, iword)) {
                    p.coordex[ip] = x;
                }
            } else {
                p.coordex[ip] = 0.0;
            }
            let mut ierr0 = 0;
            let npgausstest = (wd.nwords() - iword) / 3;
            let iresidual = wd.nwords() - iword - 3 * npgausstest;
            if iresidual == 1 {
                iword += 1;
                if let Ok(x) = read_real(wd.word(&sometext, iword)) {
                    p.rmaxex[ip] = x;
                } else {
                    ierr0 = 1;
                }
            } else {
                p.rmaxex[ip] = 0.0;
                p.pcoefexmax[ip] = 0.0;
            }
            let mut igauss: i32 = 0;
            let mut ierr = ierr0;
            while iword < wd.nwords() && ierr == 0 && (igauss as usize) < crate::state::MPGAUSS {
                igauss += 1;
                for k in 0..3 {
                    iword += 1;
                    let v = if iword <= wd.nwords() { read_real(wd.word(&sometext, iword)).ok() } else { None };
                    match v {
                        Some(x) => match k {
                            0 => p.rpgauss[(igauss, ip)] = x,
                            1 => p.wpgauss[(igauss, ip)] = x,
                            _ => p.hpgauss[(igauss, ip)] = x,
                        },
                        None => ierr = 1,
                    }
                }
                // NB: the Fortran increments (rather than decrements) igauss on a read error here.
                if ierr != 0 {
                    igauss += 1;
                }
            }
            p.npgauss[ip] = igauss;
            let withmax = p.rmaxex[ip] > 0.0;
            if p.npgauss[ip] == 0 {
                let mut w = ListWriter::new();
                w.real(p.rminex[ip]).real(p.pcoefex[ip]).real(p.coordex[ip]);
                if withmax {
                    w.real(p.rmaxex[ip]).real(p.pcoefexmax[ip]);
                }
                out.line(&w.finish());
            } else {
                let formatline = if withmax {
                    format!("(5(1x,f10.5),{}(1x,f9.4,1x,f9.4,1x,e11.4))", p.npgauss[ip])
                } else {
                    format!("(3(1x,f10.5),{}(1x,f9.4,1x,f9.4,1x,e11.4))", p.npgauss[ip])
                };
                let mut v = vec![FV::R(p.rminex[ip]), FV::R(p.pcoefex[ip]), FV::R(p.coordex[ip])];
                if withmax {
                    v.push(FV::R(p.rmaxex[ip]));
                    v.push(FV::R(p.pcoefexmax[ip]));
                }
                for ig in 1..=p.npgauss[ip] {
                    v.push(FV::R(p.rpgauss[(ig, ip)]));
                    v.push(FV::R(p.wpgauss[(ig, ip)]));
                    v.push(FV::R(p.hpgauss[(ig, ip)]));
                }
                out.line(&fmt(&formatline, &v));
            }
            let mut w = ListWriter::new();
            for ic in 1..=p.ncoeffp {
                w.real(p.pcoefp[(ic, ip)]);
            }
            out.line(&w.finish());
        }
        Ok(())
    }
}
