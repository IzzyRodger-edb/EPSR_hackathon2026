//! The EPSR driver — translation of `subroutine runepsr17`
//! (`epsr_standalone_rev3.f`), the body of the `epsr` command.
//!
//! The Fortran routine is one ~2500-line subroutine.  Here it is split into
//! the phases a caller (CLI or Python) needs:
//!
//! 1. [`Epsr::setup`]   — read `.EPSR.inp` (and through it `.ato`, `.pcof`),
//!    the weights, the data and the accumulators; build grids and potentials
//!    (Fortran lines 272–1054).
//! 2. [`Epsr::iterate`] — one EPSR iteration (`do 3000 iter=1,niter`):
//!    `ntimes` Monte Carlo cycles, accumulate g(r)/S(Q), refine the
//!    empirical potential, write results, update and save the input files.
//! 3. [`Epsr::run`]     — `setup` + `niter` × `iterate`, i.e. the whole command.
//!
//! All intermediate state is public so that it can be inspected between
//! iterations (see `docs/PYTHON_WRAPPER_NOTES.md`).

use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, FV};
use crate::fortran::listread::{parse_real, read_int, read_real, ListRead, TextFile};
use crate::fortran::output::OutFile;
use crate::fortran::fmath::sign;
use crate::fortran::text::{len_to_blank, parse};
use crate::fortran::{F1, F2, F3};
use crate::output::{w_int, w_intr};
use crate::potential::lsparabola1;
use crate::state::{State, MFILE, MPTS, MPZ, MQD};
use crate::structure::{interp, stog};
use std::path::{Path, PathBuf};

/// Parameters of one dataset (the `data n` block of the `.EPSR.inp` file).
#[derive(Clone, Debug, Default)]
pub struct Dataset {
    /// `datafile` as given in the input file.
    pub infili: String,
    pub fnamedata: String,
    pub fnamewts: String,
    pub nrtype: i32,
    pub rshmin: f32,
    pub szeros: f32,
    pub tweak: f32,
    pub efilereq: f32,
    pub porodadd: f32,
    pub nqmind: i32,
    pub nqmaxd: i32,
    /// Porod amplitude for this dataset (`porodamp(in)`).
    pub porodamp: f32,
}

/// Run-control values read in `runepsr17` (local variables of the Fortran).
#[derive(Clone, Debug, Default)]
pub struct Control {
    pub nq: i32,
    pub qstep: f32,
    pub ireset: i32,
    pub iinit: i32,
    pub ntimes: i32,
    pub niter: i32,
    pub intra: i32,
    pub rotfreq: i32,
    pub inter: i32,
    pub nisumm: i32,
    pub ngrsamples: i32,
    pub fwhm: f32,
    pub fwhmq: f32,
    pub nsmoop: i32,
    pub pressreq: f32,
    pub ereqmin: f32,
    pub ereqmax: f32,
    pub ereqstepread: f32,
    pub epcthresh: f32,
    pub ereqbias: f32,
    pub navmaxparab: i32,
    pub nparabstartfactor: i32,
    pub epcrmin: f32,
    pub epcfrac: f32,
    pub num_threads: i32,
    pub nmolcellin: i32,
    pub thresholde: f32,
    pub stepmisave: f32,
    pub steprisave: f32,
    pub rmaxgrsave: f32,
    pub irefsizefactor: i32,
    /// Porod scattering: Q value, power and the Q index nearest to `porodq`.
    pub porodq: f32,
    pub porodpwr: f32,
    pub nporodq: i32,
    /// `control-p` settings (pressure control by box size or repulsive cores).
    pub pressure: crate::pressure::PressureControl,
}

/// Results of the most recent iteration (what is written to `.out` / `.erg`).
#[derive(Clone, Debug, Default)]
pub struct Summary {
    pub iter: i32,
    /// Running intermolecular / intramolecular energy per molecule (`foldsv`, `foldsb`).
    pub foldsv: f32,
    pub foldsb: f32,
    /// Intermolecular energy and pressure from `press` (`energd`, `pressd`).
    pub energd: f32,
    pub pressd: f32,
    /// Absolute energy (amplitude) of the empirical potential.
    pub energabs: f32,
    /// Overall R-factor and per-dataset R-factors.
    pub rfactot: f32,
    pub rfacraw: f32,
    pub rfacd: Vec<f32>,
    pub epcup: f32,
    pub epccoeff: f32,
    /// Roughness of the averaged g(r) (`sentb`).
    pub sentb: f32,
    pub chisq: f32,
}

/// One EPSR simulation.
pub struct Epsr {
    pub st: State,
    pub ctl: Control,
    pub data: Vec<Dataset>,
    /// Simulation Q grid (`q`) and r grid (`z`), and `nz` = number of r points.
    pub q: F1<f32>,
    pub z: F1<f32>,
    pub nz: i32,
    /// Data interpolated onto the Q grid (`datai`) and errors (`errs`).
    pub datai: F2<f32>,
    pub errs: F2<f32>,
    /// Current-configuration histograms `wav(1|2, r, pair)` and structure factors `qwind(1|2, Q, pair)`.
    pub wav: F3<f32>,
    pub qwind: F3<f32>,
    pub summary: Summary,
    /// `.erg` file (one line per iteration; read for the trend fit).
    erg_path: PathBuf,
    parab: Option<crate::potential::ParabolaFit>,
    /// Records for `fort.35` (energy trace written by the Fortran to unit 35).
    pub fort35: Vec<String>,
    nqminm: i32,
    fnamew: Vec<String>,
    /// Last raw R-factor written to `.erg`.  In the Fortran, `rfacraw` is only
    /// computed on the final iteration and is otherwise an *uninitialised*
    /// local; the gfortran build writes the value from the last `.erg` line,
    /// which is what is reproduced here.
    rfacraw_last: f32,
    /// Control level and baseline of the automatic ereq control (written to `.out`).
    ctrllevel: f32,
    ctrlbaseline: f32,
    /// `.erg` must be rewritten from the start on the next write (iinit = 1).
    erg_fresh: bool,
    /// Bragg peaks, diffuse scattering and lattice refinement (`bragg.rs`).
    pub bragg: crate::bragg::Bragg,
}

impl Epsr {
    /// Phase 1 of `runepsr17`: read all input and initialise.  `dir` is the
    /// folder holding the files, `inp_name` the `.EPSR.inp` file name.
    pub fn setup(dir: &Path, inp_name: &str, verbose: bool) -> Result<Epsr> {
        Self::setup_with_lines(dir, inp_name, &[], verbose)
    }

    /// As [`Epsr::setup`], but after reading the `.EPSR.inp` file the
    /// `extra_lines` are processed exactly as if they had been appended to it
    /// (before the final `q`).  This allows parameters to be overridden
    /// without editing the file, e.g. `["ntimes 10", "data 2", "tweak 1.05"]`.
    /// The overridden values are written back to the `.EPSR.inp` file at the
    /// end of the iteration, as EPSR always does.
    pub fn setup_with_lines(dir: &Path, inp_name: &str, extra_lines: &[String], verbose: bool) -> Result<Epsr> {
        let mut st = State::new();
        st.verbose = verbose;
        let mut d = dir.to_string_lossy().into_owned();
        if !d.ends_with('/') {
            d.push('/');
        }
        st.sym.directory_name = d;
        st.setup_symbols_epsr();
        let mut f = st.setup_input_file(inp_name)?;
        st.read_input_file(&mut f)?;
        if !extra_lines.is_empty() {
            st.apply_input_lines(extra_lines)?;
        }
        if st.mol.npairs == 0 {
            return Err(EpsrError::Input("no atom pairs: the .ato file was not read".into()));
        }
        let mut ctl = Control::default();
        let sv = |st: &State, k: i32| st.sym.symbol_value[k].clone();
        // ---- datasets
        let mut iref = st.sym.ireff(st.sym.nsymboluptodata, 1);
        let nfile = read_int(&sv(&st, iref))?;
        if nfile < 1 || nfile as usize > MFILE {
            return Err(EpsrError::Input(format!("ndata = {nfile} out of range")));
        }
        st.dist.nfile = nfile;
        let npairs = st.mol.npairs;
        st.dist.alloc_pairs(npairs as usize);
        st.pot.alloc(npairs as usize, nfile as usize);
        st.wts.alloc(npairs as usize, nfile as usize);
        let mut data: Vec<Dataset> = Vec::new();
        let mut efilereqsum: f32 = 0.0;
        let dirn = st.sym.directory_name.clone();
        for _in in 1..=nfile {
            let mut ds = Dataset::default();
            iref += 1;
            ds.infili = sv(&st, iref);
            ds.fnamedata = format!("{}{}", dirn, &ds.infili[..len_to_blank(&ds.infili)]);
            iref += 1;
            let w = sv(&st, iref);
            ds.fnamewts = format!("{}{}", dirn, &w[..len_to_blank(&w)]);
            iref += 1;
            ds.nrtype = read_int(&sv(&st, iref))?;
            iref += 1;
            ds.rshmin = read_real(&sv(&st, iref))?;
            iref += 1;
            ds.szeros = read_real(&sv(&st, iref))?;
            iref += 1;
            ds.tweak = read_real(&sv(&st, iref))?;
            iref += 1;
            ds.efilereq = read_real(&sv(&st, iref))?;
            efilereqsum += ds.efilereq;
            iref += 1;
            ds.porodadd = read_real(&sv(&st, iref))?;
            data.push(ds);
        }
        if efilereqsum <= 0.0 {
            data[0].efilereq = 1.0;
        }
        // ---- global parameters (in the order of the symbol table)
        let mut iref = st.sym.ireff(1, 1);
        st.pot.feedback = read_real(&sv(&st, iref))?;
        iref += 1;
        {
            let t = sv(&st, iref);
            let w = parse(&t, 5);
            st.pot.potfacamp = parse_real(w.word(&t, 1))?;
            st.pot.potfac = if w.nwords() > 1 { parse_real(w.word(&t, 2))? } else { st.pot.potfacamp };
        }
        if st.pot.potfacamp == 0.0 {
            st.pot.potfac = 0.0;
        }
        if st.pot.potfac > 1.0 {
            st.pot.potfac = 1.0;
        }
        let first_word = |st: &State, k: i32| -> Result<f32> {
            let t = st.sym.symbol_value[k].clone();
            let w = parse(&t, 5);
            Ok(parse_real(w.word(&t, 1))?)
        };
        iref += 1;
        ctl.pressreq = first_word(&st, iref)?;
        st.pot.absolute_energy = ctl.pressreq;
        iref += 1;
        ctl.ereqmin = first_word(&st, iref)?.max(0.0);
        iref += 1;
        ctl.ereqmax = first_word(&st, iref)?.max(0.0);
        iref += 1;
        ctl.ereqstepread = first_word(&st, iref)?.max(0.0);
        st.pot.ereqstep = ctl.ereqstepread;
        st.pot.erequnit = st.pot.ereqstep;
        iref += 1;
        ctl.epcthresh = read_real(&sv(&st, iref))?;
        iref += 1;
        ctl.ereqbias = read_real(&sv(&st, iref))?;
        iref += 1;
        {
            let mut r = ListRead::from_str(&sv(&st, iref));
            ctl.navmaxparab = r.int()?.unwrap_or(0);
            ctl.nparabstartfactor = r.int()?.unwrap_or(0);
        }
        if ctl.navmaxparab <= 0 || ctl.nparabstartfactor < 1 {
            ctl.navmaxparab = 50;
            ctl.nparabstartfactor = 1;
        }
        st.pot.erequnit = st.pot.ereqstep;
        iref += 1;
        ctl.epcrmin = read_real(&sv(&st, iref))?;
        iref += 1;
        ctl.epcfrac = read_real(&sv(&st, iref))?;
        iref += 1;
        ctl.num_threads = read_int(&sv(&st, iref))?;
        iref += 1;
        {
            let t = sv(&st, iref);
            let w = parse(&t, 5);
            ctl.porodq = parse_real(w.word(&t, 1))?;
            if w.nwords() < 2 || ctl.porodq <= 0.0 {
                for ds in data.iter_mut() {
                    ds.porodamp = 0.0;
                }
                ctl.porodpwr = 0.0;
            } else {
                ctl.porodpwr = parse_real(w.word(&t, 2))?;
            }
        }
        iref += 1;
        ctl.nmolcellin = read_int(&sv(&st, iref))?;
        iref += 1;
        ctl.pressure = crate::pressure::PressureControl::parse(&sv(&st, iref));
        iref += 1;
        {
            let t = sv(&st, iref);
            let w = parse(&t, 5);
            let m = &mut st.mol;
            m.refine_intra = false;
            m.coul_intra = false;
            m.refine_intra_fac = 0.0;
            m.coul_intra_fac = 0.0;
            if w.nwords() > 0 {
                if let Ok(v) = parse_real(w.word(&t, 1)) {
                    m.refine_intra_fac = v;
                }
            }
            if m.refine_intra_fac > 0.0 {
                m.refine_intra = true;
            }
            if w.nwords() > 1 {
                if let Ok(v) = parse_real(w.word(&t, 2)) {
                    m.coul_intra_fac = v;
                }
            }
            if m.coul_intra_fac > 0.0 {
                m.coul_intra = true;
            }
            let v = fmt("(f8.3,1x,f8.3)", &[FV::R(m.refine_intra_fac), FV::R(m.coul_intra_fac)]);
            st.sym.symbol_value[iref] = v;
        }
        iref += 1;
        ctl.irefsizefactor = iref;
        {
            let t = sv(&st, iref);
            let w = parse(&t, 5);
            st.mol.sizefactor = parse_real(w.word(&t, 1))?;
            if w.nwords() > 1 {
                let mut dr = parse_real(w.word(&t, 2))?;
                if !(0.1..=1.0).contains(&dr) {
                    dr = 0.9;
                }
                st.mol.decline_rate = dr;
            } else {
                st.mol.decline_rate = 0.9;
            }
            ctl.thresholde = if w.nwords() > 2 { parse_real(w.word(&t, 3))? } else { 0.0 };
        }
        ctl.stepmisave = st.mol.stepmi;
        ctl.steprisave = st.mol.stepri;
        if st.mol.sizefactor > 1.0 {
            st.pot.potfac = 0.0;
            st.mol.stepmi = 0.0;
        } else if st.pot.potfacamp > 0.0 && st.pot.potfac <= 0.0 {
            st.pot.potfac = st.pot.potfacamp;
        }
        st.log(&fmt("('Steps: ',4(1x,e12.5))", &[FV::R(st.mol.stepmi), FV::R(st.mol.stepri), FV::R(ctl.stepmisave), FV::R(ctl.steprisave)]));
        if st.pot.feedback > 1.0 {
            st.pot.feedback = 0.9;
        }
        if st.pot.feedback < 0.0 {
            st.pot.feedback = 0.0;
        }
        if st.pot.feedback <= 0.0 {
            st.pot.potfacamp = 0.0;
            st.pot.potfac = 0.0;
        }
        let mut nonzero = 0;
        let mut sumnonzero: f32 = 0.0;
        for ds in data.iter() {
            if ds.efilereq > 0.0 {
                nonzero += 1;
                sumnonzero += ds.efilereq;
            }
        }
        if nonzero > 0 && sumnonzero > 0.0 {
            sumnonzero /= nonzero as f32;
        } else {
            sumnonzero = 1.0;
        }
        for (k, ds) in data.iter().enumerate() {
            st.pot.efacm[k + 1] = ds.efilereq / sumnonzero;
        }
        iref += 1;
        ctl.nq = read_int(&sv(&st, iref))?.min(MPTS as i32);
        iref += 1;
        ctl.qstep = read_real(&sv(&st, iref))?;
        iref += 1;
        {
            let t = sv(&st, iref);
            let w = parse(&t, 5);
            st.mol.qscaleinputfactor = parse_real(w.word(&t, 1))?;
            st.mol.selfbroadpwr = if w.nwords() > 1 { parse_real(w.word(&t, 2))? } else { 0.0 };
            st.mol.coarsegraining = st.mol.qscaleinputfactor > 1.0;
        }
        iref += 1;
        ctl.ireset = read_int(&sv(&st, iref))?;
        if ctl.ireset == 1 {
            ctl.pressreq = 0.0;
            st.pot.ereqstep = 0.0;
            st.pot.erequnit = 0.0;
        } else if ctl.ireset == 2 && st.pot.ereqstep > 0.0 {
            ctl.pressreq = st.pot.erequnit;
        }
        st.pot.absolute_energy = ctl.pressreq;
        iref += 1;
        ctl.iinit = read_int(&sv(&st, iref))?;
        iref += 1;
        ctl.ntimes = read_int(&sv(&st, iref))?;
        iref += 1;
        ctl.niter = read_int(&sv(&st, iref))?;
        iref += 1;
        st.inp.nsums = read_int(&sv(&st, iref))?;
        iref += 1;
        ctl.intra = read_int(&sv(&st, iref))?;
        iref += 1;
        ctl.rotfreq = read_int(&sv(&st, iref))?;
        iref += 1;
        ctl.inter = read_int(&sv(&st, iref))?.max(1);
        ctl.nisumm = (0.3f32 * ctl.intra as f32) as i32;
        iref += 1;
        st.dist.rho = read_real(&sv(&st, iref))?;
        iref += 1;
        st.dist.cellst = read_real(&sv(&st, iref))?;
        iref += 1;
        ctl.rmaxgrsave = read_real(&sv(&st, iref)).unwrap_or(0.0);
        st.inp.rmaxgr = ctl.rmaxgrsave;
        iref += 1;
        ctl.ngrsamples = read_int(&sv(&st, iref))?;
        iref += 1;
        ctl.fwhm = read_real(&sv(&st, iref))?;
        iref += 1;
        ctl.fwhmq = read_real(&sv(&st, iref))?;
        iref += 1;
        ctl.nsmoop = read_int(&sv(&st, iref))?;
        iref += 3;
        {
            let mut r = ListRead::from_str(&sv(&st, iref));
            st.inp.lorchqwidth = 0.0;
            st.inp.lorchqwidthq = 0.0;
            let _ = r.get_r(&mut st.inp.lorchqwidth).and_then(|_| r.get_r(&mut st.inp.lorchqwidthq));
        }
        let bragg = crate::bragg::Bragg::read_keywords(&mut st, &mut iref, ctl.nq, ctl.qstep)?;
        iref += 1;
        st.dist.rejrate = read_real(&sv(&st, iref))?;
        iref += 1;
        {
            let t = sv(&st, iref);
            let w = parse(&t, 5);
            if w.nwords() > 0 {
                if let Ok(v) = parse_real(w.word(&t, 1)) {
                    st.dist.qminm = v;
                }
                st.dist.qminp = if w.nwords() > 1 { parse_real(w.word(&t, 2)).unwrap_or(st.dist.qminp) } else { 0.0 };
            } else {
                st.dist.qminm = 0.0;
                st.dist.qminp = 0.0;
            }
            st.sym.symbol_value[iref] = fmt("(2(1x,f10.5))", &[FV::R(st.dist.qminm), FV::R(st.dist.qminp)]);
        }
        st.dist.cell = st.dist.cellmax;
        // ---- grids, accumulators, weights
        let (q, z, nz) = st.inita(ctl.ireset, ctl.iinit, ctl.nq, ctl.qstep, ctl.fwhm, ctl.fwhmq, ctl.intra)?;
        let mut nqminm = 1;
        for i in 1..=ctl.nq {
            if q[i] <= st.dist.qminm {
                nqminm = i;
            }
        }
        let fnamew: Vec<String> = data.iter().map(|d| d.fnamewts.clone()).collect();
        st.get_wts(&fnamew, ctl.iinit, &q, ctl.nq)?;
        if ctl.iinit != 0 {
            ctl.inter = st.inp.nesumt;
        } else {
            st.inp.nesumt = ctl.inter;
        }
        st.normalise_coordinates();
        st.get_maxmolsize();
        st.log(&format!(" epsr> Maximum molecular dimension (A): {}", st.mol.maxmolsize * st.mol.dmax));
        st.dist.upol *= st.inp.rcharge;
        let nprev = if st.inp.nsums >= 0 { st.inp.nsum } else { st.inp.nesum };
        st.initg(nprev)?;
        st.initpt()?;
        // Determine position of data for Porod scattering, if applicable.
        ctl.nporodq = 0;
        if ctl.porodq > 0.0 {
            let mut foundporodq = false;
            let mut porodqmin = ctl.porodq;
            let mut iq = 0;
            while iq < ctl.nq && !foundporodq {
                iq += 1;
                let porodqdif = (ctl.porodq - q[iq]).abs();
                foundporodq = porodqdif > porodqmin;
                if !foundporodq {
                    porodqmin = porodqdif;
                    ctl.nporodq = iq;
                }
            }
        }
        // ---- read the datasets
        let mut datai: F2<f32> = F2::new(MPTS, nfile as usize);
        let mut errs: F2<f32> = F2::new(MPTS, nfile as usize);
        for (k, ds) in data.iter_mut().enumerate() {
            let inn = (k + 1) as i32;
            let mut f = TextFile::open(Path::new(&ds.fnamedata)).map_err(|e| EpsrError::Io(format!("cannot open data file {}: {e}", ds.fnamedata)))?;
            st.log(&format!(" {:>11} {}", ds.nrtype, ds.fnamedata));
            let mut qd: F1<f32> = F1::new(MQD + 1);
            let mut dat: F1<f32> = F1::new(MQD + 1);
            let mut err: F1<f32> = F1::new(MQD + 1);
            let nqd: i32 = crate::datafile::read_data_points(&mut f, ds.nrtype, &mut qd, &mut dat, &mut err);
            if ds.nrtype == 3 || ds.nrtype == 5 {
                let mut q1 = qd[1];
                for i in 1..=nqd {
                    let q2 = qd[i + 1];
                    qd[i] = 0.5 * (q1 + q2);
                    q1 = q2;
                }
            }
            for i in 1..=nqd {
                dat[i] = dat[i] * st.dist.rhofac[inn] * ds.tweak;
            }
            st.log(&format!(" Read {:>5} points from file {}", nqd, ds.fnamedata));
            let (mn, mx) = interp(nqd, ctl.nq, &qd, &q, &dat, &mut datai, inn);
            let _ = interp(nqd, ctl.nq, &qd, &q, &err, &mut errs, inn);
            ds.nqmind = mn;
            ds.nqmaxd = mx;
            let dd = &mut st.dist;
            dd.nqmin[inn] = mn;
            dd.nqmax[inn] = mx;
            if dd.nqmax[inn] >= ctl.nq {
                dd.nqmax[inn] = ctl.nq - 1;
            }
            if dd.nqmin[inn] < nqminm {
                dd.nqmin[inn] = nqminm;
            }
            dd.nqpts[inn] = dd.nqmax[inn] - dd.nqmin[inn] + 1;
            // Set the Porod amplitude for this dataset
            if ctl.nporodq >= dd.nqmin[inn] && ctl.nporodq <= dd.nqmax[inn] {
                ds.porodamp = datai[(ctl.nporodq, inn)] + ds.porodadd;
                st.log(&fmt("(\"Porod amp. set to \",f10.5,\" for dataset \",i5)", &[FV::R(ds.porodamp), FV::I(inn as i64)]));
            } else {
                ds.porodamp = 0.0;
            }
        }
        // ---- .erg history and trend line
        let erg_path = PathBuf::from(format!("{}.erg", st.sym.fnamerun));
        let mut parab = None;
        if ctl.iinit != 1 {
            let mut ergf = TextFile::open(&erg_path).unwrap_or_default();
            let qdr = PathBuf::from(format!("{}.qdr", st.sym.fnamerun));
            let verbose = st.verbose;
            let log = move |s: &str| {
                if verbose {
                    println!("{s}")
                }
            };
            parab = Some(lsparabola1(&mut ergf, &qdr, ctl.navmaxparab, ctl.epcfrac, &log)?);
        }
        let np = npairs as usize;
        let ctl_iinit = ctl.iinit;
        let rfacraw_last = parab.as_ref().map(|p| p.rfacwrite).unwrap_or(0.0);
        Ok(Epsr {
            st,
            ctl,
            data,
            q,
            z,
            nz,
            datai,
            errs,
            wav: F3::new(2, MPZ, np),
            qwind: F3::new(2, MPTS, np),
            summary: Summary::default(),
            erg_path,
            parab,
            fort35: Vec::new(),
            nqminm,
            fnamew,
            rfacraw_last,
            ctrllevel: 0.0,
            ctrlbaseline: 0.0,
            erg_fresh: ctl_iinit == 1,
            bragg,
        })
    }

    fn log(&self, s: &str) {
        self.st.log(s);
    }

    /// Phase 2: one EPSR iteration (`iter` = 1..=niter).  When `iter == niter`
    /// all result files are written.
    pub fn iterate(&mut self, iter: i32) -> Result<&Summary> {
        let ctl = self.ctl.clone();
        self.st.sync_lattice();
        let nfile = self.st.dist.nfile;
        let npairs = self.st.mol.npairs;
        let nq = ctl.nq;
        let nz = self.nz;
        if ctl.ireset == 1 || iter == 1 {
            let m = &mut self.st.mol;
            for j in 1..=2 {
                for imt in 1..=m.nmoltype {
                    m.nmovei[(j, imt)] = 0;
                    m.ninter[(j, imt)] = 0;
                    m.nwhole[(j, imt)] = 0;
                    m.nmovem[(j, imt)] = 0;
                }
            }
        }
        let t0 = std::time::Instant::now();
        // ---- neighbour lists and initial energies
        let sepmax = self.st.dist.xstar[self.st.dist.nptmax] + self.st.mol.maxmolsize + 3.0 / self.st.mol.dmax;
        self.log(&format!("epsr> Maximum separation for neighbourlist = {}", fmt("(f10.5)", &[FV::R(sepmax * self.st.mol.dmax)])));
        self.st.initialise_cellneighbourlist(ctl.nmolcellin, sepmax);
        self.wrap_molecules();
        let (mut foldsv, mut foldsb) = self.st.fmolec1();
        self.fort35.push(fmt(
            "(a,i4,a,2(1x,e13.6),a,4(1x,e13.6))",
            &[
                FV::S("#Refinement: "),
                FV::I(iter as i64),
                FV::S(" Initial energies: "),
                FV::R(foldsv),
                FV::R(foldsb),
                FV::S(" Sizefactor et al.: "),
                FV::R(self.st.mol.sizefactor),
                FV::R(self.st.mol.stepmi),
                FV::R(self.st.mol.stepri),
                FV::R(ctl.pressreq),
            ],
        ));
        self.log(&format!("timing> neighbour lists + energies: {:.3} s", t0.elapsed().as_secs_f64()));
        let t1 = std::time::Instant::now();
        // ---- Monte Carlo
        for _itimes in 1..=ctl.ntimes {
            self.st.timer.sumtime1 = 0;
            self.st.timer.sumtime2 = 0;
            self.st.molfit(&mut foldsv, &mut foldsb, ctl.intra, ctl.rotfreq, ctl.ireset);
            let nisum = self.st.inp.nisum;
            if nisum == 1 {
                self.st.adjust_intra_wts();
                self.wrap_molecules();
                let (a, b) = self.st.fmolec1();
                foldsv = a;
                foldsb = b;
            }
            if nisum == 1 && self.st.mol.stepmi > 0.0 {
                self.fort35.push("#Intramolecular shake> ".into());
            }
            if nisum != 1 && (nisum - 1) == ctl.rotfreq * ((nisum - 1) / ctl.rotfreq) && self.st.mol.stepri > 0.0 {
                self.fort35.push("#Intramolecular rotation> ".into());
            }
            self.fort35.push(fmt(
                "(i6,2(1x,e13.6),2(1x,i16))",
                &[FV::I(nisum as i64), FV::R(foldsv), FV::R(foldsb), FV::I(self.st.timer.sumtime1), FV::I(self.st.timer.sumtime2)],
            ));
        }
        self.log(&format!("timing> Monte Carlo cycles: {:.3} s", t1.elapsed().as_secs_f64()));
        let t2 = std::time::Instant::now();
        // ---- accumulate structure
        self.structure_and_lattice(nq)?;
        let _oiseb = self.st.entrpc(&mut self.wav);
        self.st.inp.nsum += 1;
        let nsum = self.st.inp.nsum;
        let (fc1, fc2): (f64, f64) = if self.st.inp.nsums >= 0 {
            let fc1 = (1.0f32 / nsum as f32) as f64;
            (fc1, (nsum - 1) as f32 as f64 * fc1)
        } else {
            // running average over the last `nesumt` configurations
            let nesum = self.st.inp.nesum;
            let fc1 = (1.0f32 / (nesum + 1) as f32) as f64;
            (fc1, nesum as f32 as f64 * fc1)
        };
        self.st.sumsq(nq, &self.qwind, fc1, fc2);
        self.st.sumgr(nz, &self.wav, fc1, fc2);
        let sentb = self.st.entrpb(nsum);
        let (mut energd, pressd, energabs) = self.st.press();
        energd += self.st.dist.upol;
        // running averages of energy and pressure
        {
            let i = &mut self.st.inp;
            if i.nesum == 0 {
                i.ensum = 0.0;
                i.ensum2 = 0.0;
                i.prsum = 0.0;
                i.prsum2 = 0.0;
            }
            if i.nisum >= ctl.nisumm {
                i.nesum += 1;
                if i.nesum > i.nesumt {
                    i.nesum = i.nesumt;
                }
                let nesum1 = (i.nesum - 1) as f32;
                let ne = i.nesum as f32;
                i.ensum = (nesum1 * i.ensum + energd) / ne;
                let ensumc = energd - i.ensum;
                i.ensum2 = (nesum1 * i.ensum2 + ensumc * ensumc) / ne;
                i.prsum = (nesum1 * i.prsum + pressd) / ne;
                let prsumc = pressd - i.prsum;
                i.prsum2 = (nesum1 * i.prsum2 + prsumc * prsumc) / ne;
            }
        }
        let i = &self.st.inp;
        let enav = i.ensum;
        let enav2 = if i.ensum2 > 0.0 { i.ensum2.sqrt() } else { 0.0 };
        let prav = i.prsum;
        let prav2 = if i.prsum2 > 0.0 { i.prsum2.sqrt() } else { 0.0 };
        self.log(&format!("timing> structure accumulation + pressure: {:.3} s", t2.elapsed().as_secs_f64()));
        let t3 = std::time::Instant::now();
        // ---- data minus intramolecular part, background removal
        let mut datat: F2<f32> = F2::new(MPTS, nfile as usize);
        let mut szerostemp: F1<f32> = F1::new(MFILE);
        for (k, ds) in self.data.iter().enumerate() {
            let inn = (k + 1) as i32;
            let ix = self.st.wts.xfile[inn];
            for iq in 1..=nq {
                if iq >= ds.nqmind && iq <= ds.nqmaxd {
                    let mut sum: f32 = 0.0;
                    for ip in 1..=npairs {
                        let w = if ix > 0 { self.st.xray.weightx[(iq, ip, ix)] } else { self.st.wts.ciwt[(inn, ip)] };
                        sum = (sum as f64 + w * self.st.dist.soqintra[(iq, ip)] as f64) as f32;
                    }
                    datat[(iq, inn)] = self.datai[(iq, inn)] - sum;
                } else {
                    datat[(iq, inn)] = 0.0;
                }
            }
            // Remember to subtract the intramolecular term from the zero limit as well!
            szerostemp[inn] = 0.0;
            if ds.szeros != 0.0 {
                let mut sum: f32 = 0.0;
                for ip in 1..=npairs {
                    let w = if ix > 0 { self.st.xray.weightx[(1, ip, ix)] } else { self.st.wts.ciwt[(inn, ip)] };
                    sum = (sum as f64 + w * self.st.dist.soqintra[(1, ip)] as f64) as f32;
                }
                szerostemp[inn] = ds.szeros - sum;
            }
        }
        let mut dts: F1<f32> = F1::new(MPTS);
        let mut dti: F1<f32> = F1::new(MPTS);
        for (k, ds) in self.data.iter().enumerate() {
            let id = (k + 1) as i32;
            let (nqmn, nqmx) = (self.st.dist.nqmin[id], self.st.dist.nqmax[id]);
            let ix = self.st.wts.xfile[id];
            // Set up the gradient for extrapolating to Q=0
            let (grads0, conss0): (f32, f32) = if ds.szeros != 0.0 && nqmn > 1 {
                ((datat[(nqmn, id)] - szerostemp[id]) / (self.q[nqmn] - self.q[1]), szerostemp[id])
            } else {
                (0.0, datat[(nqmn, id)])
            };
            for i in 1..=nq {
                let mut sum: f32 = 0.0;
                for ic in 1..=npairs {
                    let w = if ix > 0 { self.st.xray.weightx[(i, ic, ix)] } else { self.st.wts.cwt[(id, ic)] };
                    sum = (sum as f64 + w * self.st.dist.soqinter[(i, ic)] as f64) as f32;
                }
                dts[i] = if i < nqmn {
                    if ds.szeros != 0.0 {
                        conss0 + grads0 * (self.q[i] - self.q[1]) - sum
                    } else {
                        0.0
                    }
                } else if i <= nqmx {
                    let mut v = datat[(i, id)] - sum;
                    // Add in the Porod scattering if required
                    if ds.porodamp > 0.0 {
                        v += ds.porodamp * (ctl.porodq / self.q[i]).powf(ctl.porodpwr);
                    }
                    v
                } else {
                    0.0
                };
            }
            if ctl.nsmoop > 0 {
                self.st.ftback(id, ds.rshmin, &self.q, &mut dts, &mut dti);
                for i in 1..=nqmx {
                    if i >= nqmn {
                        datat[(i, id)] -= dts[i];
                    }
                }
            }
        }
        let mut soqintertemp: F2<f32> = F2::new(MPTS, npairs as usize);
        for ip in 1..=npairs {
            for iq in 1..=nq {
                soqintertemp[(iq, ip)] = self.qwind[(2, iq, ip)];
            }
        }
        // ---- refine the empirical potential
        let mut potnoi: F1<f32> = F1::new(MFILE);
        let rshmin: F1<f32> = {
            let mut v = F1::new(MFILE);
            for (k, ds) in self.data.iter().enumerate() {
                v[k + 1] = ds.rshmin;
            }
            v
        };
        if ctl.pressreq > 0.0 {
            let q = self.q.clone();
            self.st.ptgnsq(nq, &q, &datat, &soqintertemp, &szerostemp, &rshmin, &mut potnoi)?;
        }
        let potfac = self.st.pot.potfac;
        self.st.repulsivept(potfac)?;
        self.log(&format!("timing> potential refinement: {:.3} s", t3.elapsed().as_secs_f64()));
        let t4 = std::time::Instant::now();
        // ---- results (written on the last iteration)
        let mut rfacd: Vec<f32> = vec![0.0; nfile as usize];
        let mut rfactot: f32 = 0.0;
        let mut rfacraw: f32 = self.rfacraw_last;
        let mut epcup: f32 = self.summary.epcup;
        let mut epccoeff: f32 = self.summary.epccoeff;
        let mut chisq: f32 = 0.0;
        let mut rfacwrite = self.parab.as_ref().map(|p| p.rfacwrite).unwrap_or(0.0);
        if iter == ctl.niter {
            let res = self.write_results(&mut datat)?;
            rfacd = res.0;
            rfactot = res.1;
            rfacraw = res.2;
            self.rfacraw_last = rfacraw;
            epcup = res.3;
            epccoeff = res.4;
            chisq = res.5;
            if let Some(p) = self.parab.as_mut() {
                p.lastqf = p.nextqf;
                p.nextqf = epccoeff;
            }
            if self.st.pot.ereqstep > 0.0 && potfac > 0.0 {
                self.automatic_ereq(energabs);
            }
            self.st.norma()?;
            if self.ctl.ireset != 0 {
                self.ctl.ireset = 0;
            }
            if self.ctl.iinit > 0 {
                self.ctl.iinit = 0;
            }
            if self.ctl.iinit < 0 {
                self.ctl.iinit = 1;
            }
            rfacwrite = rfactot;
        }
        // ---- .erg line
        let (aparab, bparab, cparab, nparab, bav) = match &self.parab {
            Some(p) => (p.a, p.b, p.c, p.nav, p.bav),
            None => (0.0, 0.0, 0.0, 0, 0.0),
        };
        {
            // With iinit = 1 the .erg history is not read and the file is
            // rewritten from the start (the Fortran writes to the freshly opened unit).
            let mut out = if self.erg_fresh { OutFile::create(&self.erg_path) } else { OutFile::append(&self.erg_path) };
            self.erg_fresh = false;
            out.line(&fmt(
                "(10(1x,e13.6))",
                &[
                    FV::R(energd),
                    FV::R(pressd),
                    FV::R(rfacraw),
                    FV::R(energabs),
                    FV::R(self.st.dist.rho),
                    FV::R(epcup),
                    FV::R(epccoeff),
                    FV::R(aparab),
                    FV::R(bparab),
                    FV::R(cparab),
                ],
            ));
            out.close()?;
        }
        // ---- .out summary
        self.write_out_file(iter, nsum, rfacwrite, epcup, energabs, epccoeff, foldsv, foldsb, nparab, bav, sentb, energd, pressd, enav, enav2, prav, prav2, &rfacd)?;
        // ---- update the symbol table and save .inp/.pcof/.ato
        self.update_symbols()?;
        self.st.unnormalise_coordinates();
        if self.st.mol.sizefactor > 1.0 && foldsv < ctl.thresholde {
            self.st.mol.sizefactor *= self.st.mol.decline_rate;
            if self.st.mol.sizefactor < 1.0 {
                self.st.mol.sizefactor = 1.0;
            }
        }
        self.st.sym.symbol_value[ctl.irefsizefactor] =
            fmt("(f10.5,1x,f10.5,1x,e12.5)", &[FV::R(self.st.mol.sizefactor), FV::R(self.st.mol.decline_rate), FV::R(ctl.thresholde)]);
        let irmax = self.st.sym.search_for_symbol("rmaxgr", 0);
        let ir = self.st.sym.ireff(irmax, 1);
        self.st.sym.symbol_value[ir] = fmt("(e13.6)", &[FV::R(ctl.rmaxgrsave)]);
        self.st.mol.stepmi = ctl.stepmisave;
        self.st.mol.stepri = ctl.steprisave;
        self.st.control_pressure(&ctl.pressure, pressd)?;
        self.st.save_inp_file()?;
        self.st.save_pcof_file()?;
        let fato = self.st.sym.fnameato.clone();
        if self.st.mol.nmol > 0 {
            self.st.moloutc(&fato)?;
        }
        self.st.mol.dmaxmol = self.st.mol.dmax * self.st.mol.sizefactor;
        self.st.normalise_coordinates();
        self.log(&format!("timing> output files: {:.3} s", t4.elapsed().as_secs_f64()));
        self.summary = Summary {
            iter,
            foldsv,
            foldsb,
            energd,
            pressd,
            energabs,
            rfactot,
            rfacraw,
            rfacd,
            epcup,
            epccoeff,
            sentb,
            chisq,
        };
        Ok(&self.summary)
    }

    /// Automatic control of the EP amplitude `ereq` (Fortran lines 2180–2262):
    /// step `ereq` up or down depending on the trend of the quality factor
    /// in the `.erg` history (fitted by `lsparabola1`).
    fn automatic_ereq(&mut self, energabs: f32) {
        let ereqstep = self.st.pot.ereqstep;
        let erequnit = self.st.pot.erequnit;
        let nparabstart = self.ctl.navmaxparab * self.ctl.nparabstartfactor;
        let (aparab, bparab, cparab, nparab, lastqf, nextqf) = match &self.parab {
            Some(p) => (p.a, p.b, p.c, p.nav, p.lastqf, p.nextqf),
            None => (0.0, 0.0, 0.0, 0, 0.0, 0.0),
        };
        let mut factor: f32;
        if nparab > nparabstart && nparabstart > 0 {
            if cparab != 0.0 {
                let epcthreshuse = self.ctl.epcthresh;
                let mut grad = bparab + 2.0 * cparab * energabs;
                if aparab < epcthreshuse {
                    grad = epcthreshuse - aparab;
                }
                let eminmaxdif = sign(self.ctl.ereqbias, grad);
                if eminmaxdif.abs() != 0.0 {
                    let term = 1.0 + eminmaxdif * (eminmaxdif - 4.0 * self.st.rng.ran1() + 2.0);
                    factor = 2.0 * ereqstep * (1.0 - term.sqrt()) / eminmaxdif;
                } else {
                    factor = 2.0 * ereqstep * (2.0 * self.st.rng.ran1() - 1.0);
                }
            } else if bparab != 0.0 {
                let fcurrent = aparab + energabs * bparab;
                self.ctrllevel = aparab;
                self.ctrlbaseline = fcurrent;
                let grad = self.ctrlbaseline - self.ctrllevel;
                let eminmaxdif = sign(self.ctl.ereqbias, grad);
                if eminmaxdif.abs() != 0.0 {
                    factor = -ereqstep * sign(1.0, grad);
                } else {
                    factor = ereqstep * (2.0 * self.st.rng.ran1() - 1.0);
                }
            } else {
                factor = ereqstep;
            }
            if nextqf > lastqf {
                factor = 0.0;
            }
        } else {
            factor = ereqstep;
        }
        if factor > 0.0 {
            factor = erequnit.max(factor);
        }
        if energabs > 2.0 * erequnit && factor < 0.0 {
            factor = (-erequnit).min(factor);
        }
        let (ereqmin, ereqmax) = (self.ctl.ereqmin, self.ctl.ereqmax);
        let mut pressreq = erequnit * crate::fortran::fmath::nint((energabs + factor) / erequnit) as f32;
        if pressreq <= ereqmin.max(erequnit) && factor < 0.0 {
            pressreq = ereqmin + 10.0 * erequnit;
        }
        if ereqmax > 0.0 {
            if pressreq >= ereqmax {
                pressreq = ereqmin + 0.75 * (ereqmax - ereqmin);
            }
            pressreq = pressreq.min(ereqmax);
        }
        self.log(&format!(" epsr-automatic> {nparab:>11} {nparabstart:>11} {energabs} {pressreq}"));
        self.ctl.pressreq = pressreq;
    }

    /// Put every molecule centre back into the box (`minimagec` loop).
    fn wrap_molecules(&mut self) {
        let m = &mut self.st.mol;
        for im in 1..=m.nmol {
            let (mut x, mut y, mut z) = (m.xadd[im], m.yadd[im], m.zadd[im]);
            crate::neighbour::minimagec(&mut x, &mut y, &mut z);
            m.xadd[im] = x;
            m.yadd[im] = y;
            m.zadd[im] = z;
            if m.tethered[im] {
                let (mut x, mut y, mut z) = (m.phix[im], m.phiy[im], m.phiz[im]);
                crate::neighbour::minimagec(&mut x, &mut y, &mut z);
                m.phix[im] = x;
                m.phiy[im] = y;
                m.phiz[im] = z;
            }
        }
    }

    /// The `if(iter.eq.niter)` output block of `runepsr17`.  Returns
    /// (rfacd, rfactot, rfacraw, epcup, epccoeff, chisq).
    fn write_results(&mut self, datat: &mut F2<f32>) -> Result<(Vec<f32>, f32, f32, f32, f32, f32)> {
        let nq = self.ctl.nq;
        let nz = self.nz;
        let nfile = self.st.dist.nfile;
        let npairs = self.st.mol.npairs;
        let run = self.st.sym.fnamerun.clone();
        let verbose = self.st.verbose;
        let log = move |s: &str| {
            if verbose {
                println!("{s}")
            }
        };
        let z = self.z.clone();
        self.st.sumpe(nz, &z)?;
        // titles
        let label = |s: &str| -> String { s[..len_to_blank(s).min(s.len())].to_string() };
        let mut titlepsf: Vec<String> = Vec::new();
        for ic in 1..=npairs {
            let irefp = self.st.pot.irefpcof[ic];
            let a1 = label(&self.st.pot.tlab1[irefp]);
            let a2 = label(&self.st.pot.tlab2[irefp]);
            titlepsf.push(format!("   {a1}-{a2}"));
        }
        let titledat: Vec<String> = self
            .data
            .iter()
            .map(|d| {
                let ilent2 = d.infili.find('.').map(|p| p as i32).unwrap_or(-1);
                let ilent1 = (ilent2 - 28).max(1);
                crate::fortran::text::substr(&d.infili, ilent1, ilent2)
            })
            .collect();
        let titler = "   r";
        let titleq = "   Q";
        // estimated partial structure factors from the data
        let mut datap: F2<f32> = F2::new(MPTS, npairs as usize);
        let mut grlim: F1<f32> = F1::new(npairs as usize);
        let mut nqminp: F1<i32> = F1::new(npairs as usize);
        let mut nqmaxp: F1<i32> = F1::new(npairs as usize);
        {
            let st = &self.st;
            let nqinv = st.wts.nqinv;
            for ic in 1..=npairs {
                grlim[ic] = 1.0;
                nqminp[ic] = 1;
                nqmaxp[ic] = nq;
                let mut iqinv: i32 = 2;
                for iq in 1..=nq {
                    let mut datt: f32 = 0.0;
                    if st.wts.nxfile > 0 {
                        // interpolate the inverse weights matrix to this Q
                        while st.wts.qinv[iqinv] < iq && iqinv < nqinv {
                            iqinv += 1;
                        }
                        let ratinter: f32 = (iq - st.wts.qinv[iqinv - 1]) as f32 / (st.wts.qinv[iqinv] - st.wts.qinv[iqinv - 1]) as f32;
                        let ratinter1: f32 = 1.0 - ratinter;
                        for inn in 1..=nfile {
                            let cwttemp = (ratinter1 as f64 * st.wts.cwttq[(iqinv - 1, ic, inn)] + ratinter as f64 * st.wts.cwttq[(iqinv, ic, inn)]) as f32;
                            let mut term = if iq < st.dist.nqmin[inn] {
                                datat[(st.dist.nqmin[inn], inn)] * cwttemp
                            } else {
                                datat[(iq, inn)] * cwttemp
                            };
                            if st.pot.wsum[inn] <= 0.0 {
                                term = 0.0;
                            }
                            datt += term;
                        }
                        if st.pot.feedback < 1.0 {
                            for inn in 1..=npairs {
                                let inref = inn + nfile;
                                let cwttemp =
                                    (ratinter1 as f64 * st.wts.cwttq[(iqinv - 1, ic, inref)] + ratinter as f64 * st.wts.cwttq[(iqinv, ic, inref)]) as f32;
                                datt += st.dist.soqinter[(iq, inn)] * cwttemp;
                            }
                        }
                        datap[(iq, ic)] = datt;
                        continue;
                    }
                    let ratinter1: f32 = 1.0;
                    for inn in 1..=nfile {
                        let cwttemp = (ratinter1 as f64 * st.wts.cwttq[(1, ic, inn)]) as f32;
                        let mut term = if iq < st.dist.nqmin[inn] {
                            datat[(st.dist.nqmin[inn], inn)] * cwttemp
                        } else {
                            datat[(iq, inn)] * cwttemp
                        };
                        if st.pot.wsum[inn] <= 0.0 {
                            term = 0.0;
                        }
                        datt += term;
                    }
                    if st.pot.feedback < 1.0 {
                        for inn in 1..=npairs {
                            let inref = inn + nfile;
                            let cwttemp = (ratinter1 as f64 * st.wts.cwttq[(1, ic, inref)]) as f32;
                            datt += st.dist.soqinter[(iq, inn)] * cwttemp;
                        }
                    }
                    datap[(iq, ic)] = datt;
                }
            }
        }
        let mut datagr: F2<f32> = F2::new(MPZ, npairs as usize);
        stog(self.st.mol.rhocg, &self.q, &datap, &z, &mut datagr, &nqminp, &nqmaxp, nz, npairs, &grlim);
        let mut errsav: F2<f32> = F2::new(MPTS, npairs as usize);
        let mut chisq: f32 = 0.0;
        {
            let d = &self.st.dist;
            for ic in 1..=npairs {
                let mut rf: f32 = 0.0;
                for i in 1..=nq {
                    if datap[(i, ic)] != 0.0 {
                        let dif = datap[(i, ic)] - d.soqinter[(i, ic)];
                        errsav[(i, ic)] = dif;
                        rf += dif * dif;
                    } else {
                        errsav[(i, ic)] = 0.0;
                    }
                }
                rf /= nq as f32;
                chisq += rf;
            }
        }
        let d = &self.st.dist;
        w_int(nq, &format!("{run}.f"), npairs, &self.q, &d.soqinter, &d.sq2inter, titleq, &titlepsf, &log)?;
        w_int(nq, &format!("{run}.s"), npairs, &self.q, &d.soqintra, &d.sq2intra, titleq, &titlepsf, &log)?;
        let mut lpt1 = nq;
        let mut lpt2 = 1;
        for ds in self.data.iter() {
            lpt1 = lpt1.min(ds.nqmind);
            lpt2 = lpt2.max(ds.nqmaxd);
        }
        w_intr(lpt1, lpt2, &format!("{run}.d"), npairs, &self.q, &errsav, &d.sq2inter, titleq, &titlepsf, &log)?;
        w_intr(lpt1, lpt2, &format!("{run}.q"), npairs, &self.q, &datap, &d.sq2inter, titleq, &titlepsf, &log)?;
        // add the intramolecular part back to the data
        let w = &self.st.wts;
        for (k, ds) in self.data.iter().enumerate() {
            let inn = (k + 1) as i32;
            let ix = w.xfile[inn];
            for i in 1..=nq {
                let mut sum: f32 = 0.0;
                if i >= ds.nqmind && i <= ds.nqmaxd {
                    for ic in 1..=npairs {
                        let wt = if ix > 0 { self.st.xray.weightx[(i, ic, ix)] } else { w.ciwt[(inn, ic)] };
                        sum = (sum as f64 + wt * d.soqintra[(i, ic)] as f64) as f32;
                    }
                }
                datat[(i, inn)] += sum;
            }
        }
        w_intr(lpt1, lpt2, &format!("{run}.t"), nfile, &self.q, datat, &self.errs, titleq, &titledat, &log)?;
        // R-factors
        let mut soqism: F2<f32> = F2::new(MPTS, nfile as usize);
        let mut backsc: F2<f32> = F2::new(MPTS, nfile as usize);
        let mut ufile: F2<f32> = F2::new(MPTS, nfile as usize);
        let mut rfacd: Vec<f32> = vec![0.0; nfile as usize];
        let mut rfactot: f32 = 0.0;
        let mut rfacraw: f32 = 0.0;
        for (k, ds) in self.data.iter().enumerate() {
            let inn = (k + 1) as i32;
            let mut rfdraw: f32 = 0.0;
            for i in 1..=nq {
                let mut sum: f32 = 0.0;
                let mut sumraw: f32 = 0.0;
                let ix = w.xfile[inn];
                if ix > 0 {
                    let x = &self.st.xray;
                    for ic in 1..=npairs {
                        sum = (sum as f64 + x.weightx[(i, ic, ix)] * (d.soqintra[(i, ic)] + d.soqinter[(i, ic)]) as f64) as f32;
                        sumraw = (sumraw as f64 + x.weightx[(i, ic, ix)] * (self.qwind[(1, i, ic)] + self.qwind[(2, i, ic)]) as f64) as f32;
                    }
                    sum = (sum as f64 + x.backx[(i, ix)]) as f32;
                    sumraw = (sumraw as f64 + x.backx[(i, ix)]) as f32;
                } else {
                    for ic in 1..=npairs {
                        sum = ((sum as f64 + w.ciwt[(inn, ic)] * d.soqintra[(i, ic)] as f64) + w.cwt[(inn, ic)] * d.soqinter[(i, ic)] as f64) as f32;
                        sumraw = ((sumraw as f64 + w.ciwt[(inn, ic)] * self.qwind[(1, i, ic)] as f64) + w.cwt[(inn, ic)] * self.qwind[(2, i, ic)] as f64) as f32;
                    }
                }
                // Add in the Porod scattering if required
                if self.q[i] > 0.0 && ds.porodamp > 0.0 {
                    let porod = ds.porodamp * (self.ctl.porodq / self.q[i]).powf(self.ctl.porodpwr);
                    sum += porod;
                    sumraw += porod;
                }
                soqism[(i, inn)] = sum;
                ufile[(i, inn)] = sumraw;
                let diffraw;
                if i >= ds.nqmind && i <= ds.nqmaxd {
                    backsc[(i, inn)] = self.datai[(i, inn)] - soqism[(i, inn)];
                    diffraw = self.datai[(i, inn)] - sumraw;
                } else {
                    backsc[(i, inn)] = 0.0;
                    diffraw = 0.0;
                }
                if i >= d.nqmin[inn] && i <= d.nqmax[inn] {
                    rfacd[k] += backsc[(i, inn)] * backsc[(i, inn)];
                    rfdraw += diffraw * diffraw;
                }
            }
            if d.nqpts[inn] > 0 {
                rfacd[k] /= d.nqpts[inn] as f32;
                rfdraw /= d.nqpts[inn] as f32;
            } else {
                rfacd[k] = 0.0;
                rfdraw = 0.0;
            }
            rfactot += rfacd[k];
            rfacraw += rfdraw;
        }
        rfactot /= nfile as f32;
        rfacraw /= nfile as f32;
        w_intr(1, nq, &format!("{run}.u"), nfile, &self.q, &soqism, &self.errs, titleq, &titledat, &log)?;
        w_intr(lpt1, lpt2, &format!("{run}.v"), nfile, &self.q, &backsc, &self.errs, titleq, &titledat, &log)?;
        w_int(nz, &format!("{run}.g"), npairs, &z, &d.gorinter, &d.gr2inter, titler, &titlepsf, &log)?;
        w_int(nz, &format!("{run}.l"), npairs, &z, &d.entsig, &d.entsig, titler, &titlepsf, &log)?;
        w_int(nz, &format!("{run}.y"), npairs, &z, &d.gorintra, &d.gr2intra, titler, &titlepsf, &log)?;
        w_int(nz, &format!("{run}.r"), npairs, &z, &datagr, &d.gr2inter, titler, &titlepsf, &log)?;
        let mut gz2ism: F2<f32> = F2::new(MPZ, nfile as usize);
        let grzerolim: F1<f32> = F1::new(MFILE);
        for inn in 1..=nfile {
            for i in 1..=nz {
                let mut sum1: f32 = 0.0;
                for ic in 1..=npairs {
                    let sum = (w.ciwt[(inn, ic)] * d.gorintra[(i, ic)] as f64) as f32;
                    sum1 += sum;
                }
                gz2ism[(i, inn)] = sum1;
            }
        }
        let mut gofzism: F2<f32> = F2::new(MPZ, nfile as usize);
        stog(d.rho, &self.q, &soqism, &z, &mut gofzism, &d.nqmin, &d.nqmax, nz, nfile, &grzerolim);
        w_int(nz, &format!("{run}.x"), nfile, &z, &gofzism, &gz2ism, titler, &titledat, &log)?;
        let mut dataigr: F2<f32> = F2::new(MPZ, nfile as usize);
        stog(d.rho, &self.q, &self.datai, &z, &mut dataigr, &d.nqmin, &d.nqmax, nz, nfile, &grzerolim);
        w_int(nz, &format!("{run}.w"), nfile, &z, &dataigr, &gz2ism, titler, &titledat, &log)?;
        let p = &self.st.pot;
        w_int(nz, &format!("{run}.o"), npairs, &z, &p.ecoul, &p.epot, titler, &titlepsf, &log)?;
        w_int(nz, &format!("{run}.p"), npairs, &z, &p.xavs, &p.xav, titler, &titlepsf, &log)?;
        // R-space quality coefficient
        let mut datagr2: F2<f32> = F2::new(MPZ, npairs.max(nfile) as usize);
        stog(d.rho, &self.q, &ufile, &z, &mut datagr2, &d.nqmin, &d.nqmax, nz, nfile, &grzerolim);
        let mut nugdif = 0;
        let mut epcup: f32 = 0.0;
        for inn in 1..=nfile {
            let mut ugdif: Vec<f32> = vec![0.0; (nz + 1) as usize];
            for i in 1..=nz {
                let mut gortemp: f32 = 0.0;
                let ix = w.xfile[inn];
                for ip in 1..=npairs {
                    let wt = if ix > 0 { self.st.xray.weightx[(1, ip, ix)] } else { w.cwt[(inn, ip)] };
                    gortemp = (gortemp as f64 + wt * self.wav[(2, i, ip)] as f64) as f32;
                }
                gortemp /= d.ycoh[i];
                if gortemp != 0.0 && z[i] > self.ctl.epcrmin {
                    ugdif[i as usize] = datagr2[(i, inn)] - dataigr[(i, inn)];
                }
            }
            for i in 1..=nz {
                if z[i] > self.ctl.epcrmin && ugdif[i as usize] != 0.0 {
                    nugdif += 1;
                    let prod = ugdif[i as usize];
                    epcup += prod * prod;
                }
            }
        }
        epcup /= nugdif as f32;
        let epccoeff = self.ctl.epcfrac * epcup + (1.0 - self.ctl.epcfrac) * rfacraw;
        // running coordination numbers (.z01)
        let mut titlez: Vec<String> = Vec::new();
        let mut zc1: F2<f32> = F2::new(MPZ, npairs as usize);
        let mut zc2: F2<f32> = F2::new(MPZ, npairs as usize);
        let m = &self.st.mol;
        for ip in 1..=npairs {
            let irefp = p.irefpcof[ip];
            let a1 = label(&p.tlab1[irefp]);
            let a2 = label(&p.tlab2[irefp]);
            let ilent1 = a1.len() as i32;
            let ilent2 = a2.len() as i32;
            let icharttitle = ilent1 + ilent2 + 1;
            let iremain = 14 - icharttitle - 3;
            let ibetween = iremain + 4;
            let formattitle = fmt("('(3x,a,1h-,a,',i2.2,'x,a,1h-,a,',i2.2,'x)')", &[FV::I(ibetween as i64), FV::I(iremain as i64)]);
            let t = fmt(&formattitle, &[FV::S(&a1), FV::S(&a2), FV::S(&a2), FV::S(&a1)]);
            titlez.push(crate::fortran::text::fstr(&t, 29));
            zc1[(1, ip)] = 0.0;
            zc2[(1, ip)] = 0.0;
            let jt1 = m.n2extr[ip];
            let jt2 = m.n1extr[ip];
            for j in 2..=nz {
                let cprod1 = m.cfrac[jt1] * d.ycoh[j];
                let cprod2 = m.cfrac[jt2] * d.ycoh[j];
                zc1[(j, ip)] = zc1[(j - 1, ip)] + d.gorinter[(j, ip)] * cprod1;
                zc2[(j, ip)] = zc2[(j - 1, ip)] + d.gorinter[(j, ip)] * cprod2;
            }
        }
        w_int(nz, &format!("{run}.z"), npairs, &z, &zc1, &zc2, titler, &titlez, &log)?;
        Ok((rfacd, rfactot, rfacraw, epcup, epccoeff, chisq))
    }

    #[allow(clippy::too_many_arguments)]
    fn write_out_file(
        &self,
        iter: i32,
        nsum: i32,
        rfacwrite: f32,
        epcup: f32,
        energabs: f32,
        epccoeff: f32,
        foldsv: f32,
        foldsb: f32,
        nparab: i32,
        bav: f32,
        sentb: f32,
        energd: f32,
        pressd: f32,
        enav: f32,
        enav2: f32,
        prav: f32,
        prav2: f32,
        rfacd: &[f32],
    ) -> Result<()> {
        let st = &self.st;
        let fname = format!("{}.out", st.sym.fnamerun);
        let mut out = OutFile::create(Path::new(&fname));
        let (date, time) = now_strings();
        out.line(&fmt(
            "(/2x,'Current date: ',a2,':',a2,':',a2,' on ',a2,'-',a2,'-',a4/2x,'EPSR version ',a,/)",
            &[
                FV::S(&time[0..2]),
                FV::S(&time[2..4]),
                FV::S(&time[4..6]),
                FV::S(&date[6..8]),
                FV::S(&date[4..6]),
                FV::S(&date[0..4]),
                FV::Ch(crate::VERSION_DATE, 20),
            ],
        ));
        let ctrllevel: f32 = self.ctrllevel;
        let ctrlbaseline: f32 = self.ctrlbaseline;
        let f147 = "(1x,' Trial no. ',i5,5x,' Iteration no. ',i5,5x,' No. of configurations in sum ',i5\
,/1x,' Chi-square ',e13.6,' R-space coeff.: ',e13.6,' Abs. empirical energy: ',1x,e13.6\
,/1x,' Quality factor ',e13.6,' Control level, baseline and threshold: ',3(e13.6,1x)\
,/1x,' Pair energies: repulsive ',e13.6,' harmonic ',e13.6,' per molecule'\
,/1x,' Tethered energy: ',e13.6,' per molecule'\
,/1x,' Coordination constraint energy: ',e13.6,' per molecule'\
,/1x,' Three-body energy: ',e13.6,' per molecule'\
,/1x,' Three-body pressure: ',e13.6,' kbar'\
,/1x,' No. of .erg lines: ',i7,' Average gradient of quality factor: ',5(e13.6,1x))";
        out.line(&fmt(
            f147,
            &[
                FV::I(self.ctl.ntimes as i64),
                FV::I(iter as i64),
                FV::I(nsum as i64),
                FV::R(rfacwrite),
                FV::R(epcup),
                FV::R(energabs),
                FV::R(epccoeff),
                FV::R(ctrllevel),
                FV::R(ctrlbaseline),
                FV::R(self.ctl.epcthresh),
                FV::R(foldsv),
                FV::R(foldsb),
                FV::R(self.st.ftethered),
                FV::R(self.st.cc.fcoordconstraints),
                FV::R(self.st.tb.energthreeb),
                FV::R(self.st.tb.pressthreeb),
                FV::I(nparab as i64),
                FV::R(bav),
            ],
        ));
        if st.tb.nthreeb > 0 {
            let mut w = crate::fortran::format::ListWriter::new();
            w.chars(" Energy for each triplet: ");
            for id in 1..=st.tb.nthreeb {
                w.real(st.tb.sume_triplet_dist[id] / st.mol.nmol as f32);
            }
            out.line(&w.finish());
        }
        out.line(&fmt("(1x,' Sizefactor ',f10.5,' Number density ',e12.5)", &[FV::R(st.mol.sizefactor), FV::R(st.dist.rho)]));
        out.line(&fmt("('  Net charge = ',f10.5)", &[FV::R(st.mol.sumcharge)]));
        let mut v: Vec<FV> = Vec::new();
        for (j, r) in rfacd.iter().enumerate() {
            v.push(FV::I((j + 1) as i64));
            v.push(FV::R(*r));
        }
        out.line(&fmt("((1x,' R-factor for data file ',i2,':',1x,e11.4))", &v));
        out.line(&fmt(
            "(1x,' Noise in average g(r) ',e13.6\
/1x,' Intermolecular energy and U_pol ',e13.6,1x,e13.6,' kJ/mole'\
/1x,' Intermolecular pressure ',e13.6,' kbar'\
/1x,' Average dipole moment and std. ',e13.6,1x,e13.6,' D'\
/1x,' Av. intermolecular energy ',e13.6,'+-',e13.6,' kJ/mole'\
/1x,' Av. intermolecular pressure ',e13.6,'+-',e13.6,' kbar')",
            &[
                FV::R(sentb / st.dist.rho),
                FV::R(energd),
                FV::R(st.dist.upol),
                FV::R(pressd),
                FV::R(st.dist.amubar),
                FV::R(st.dist.amubar2),
                FV::R(enav),
                FV::R(enav2),
                FV::R(prav),
                FV::R(prav2),
            ],
        ));
        let m = &st.mol;
        let mut tot = [[0i32; 4]; 2];
        for (j, row) in tot.iter_mut().enumerate() {
            for imt in 1..=m.nmoltype {
                row[0] += m.nmovei[(j + 1, imt)];
                row[1] += m.ninter[(j + 1, imt)];
                row[2] += m.nwhole[(j + 1, imt)];
                row[3] += m.nmovem[(j + 1, imt)];
            }
        }
        let v: Vec<FV> = tot.iter().flat_map(|r| r.iter().map(|x| FV::I(*x as i64))).collect();
        out.line(&fmt("(1x,' Number of rejects in MOLFIT:-'/((1x,4(i7,1x))))", &v));
        out.close()?;
        Ok(())
    }

    /// Store updated run-control values in the symbol table (Fortran lines 2395–2427).
    fn update_symbols(&mut self) -> Result<()> {
        let st = &mut self.st;
        let ctl = &mut self.ctl;
        let ip = st.sym.search_for_symbol("potfac", 0);
        let mut iref = st.sym.ireff(ip, 1);
        st.sym.symbol_value[iref] = fmt("(2(1x,f10.6))", &[FV::R(st.pot.potfacamp), FV::R(st.pot.potfac)]);
        iref += 1;
        if ctl.pressreq > 100000.0 {
            ctl.pressreq = 100000.0;
        }
        st.sym.symbol_value[iref] = fmt("(e12.5)", &[FV::R(ctl.pressreq)]);
        iref += 1;
        st.sym.symbol_value[iref] = fmt("(f10.2)", &[FV::R(ctl.ereqmin)]);
        iref += 1;
        st.sym.symbol_value[iref] = fmt("(f10.2)", &[FV::R(ctl.ereqmax)]);
        iref += 1;
        st.sym.symbol_value[iref] = fmt("(f6.3,1x,f6.3)", &[FV::R(ctl.ereqstepread), FV::R(st.pot.erequnit)]);
        iref += 3;
        st.sym.symbol_value[iref] = fmt("(i5,1x,i5)", &[FV::I(ctl.navmaxparab as i64), FV::I(ctl.nparabstartfactor as i64)]);
        if st.inp.nsums >= 0 {
            st.inp.nsums = st.inp.nsum;
        }
        let k = st.sym.search_for_symbol("nsumt", 0);
        let r = st.sym.ireff(k, 1);
        st.sym.symbol_value[r] = fmt("(i5)", &[FV::I(st.inp.nsums as i64)]);
        let k = st.sym.search_for_symbol("inter", 0);
        let r = st.sym.ireff(k, 1);
        st.sym.symbol_value[r] = fmt("(i5)", &[FV::I(ctl.inter as i64)]);
        let r = st.sym.ireff(st.sym.iresetref, 1);
        st.sym.symbol_value[r] = fmt("(i2)", &[FV::I(ctl.ireset as i64)]);
        let r = st.sym.ireff(st.sym.iinitref, 1);
        st.sym.symbol_value[r] = fmt("(i2)", &[FV::I(ctl.iinit as i64)]);
        Ok(())
    }

    /// Write the `fort.35` energy trace into the input folder.
    pub fn write_fort35(&self) -> Result<()> {
        let mut out = OutFile::create(&self.st.sym.full("fort.35"));
        for l in &self.fort35 {
            out.line(l);
        }
        out.close()?;
        Ok(())
    }

    /// Phase 3: the whole `epsr` command: `niter` iterations, then `fort.35`.
    pub fn run(&mut self) -> Result<()> {
        for iter in 1..=self.ctl.niter {
            self.iterate(iter)?;
        }
        self.write_fort35()
    }

    /// Names of the weights files of each dataset.
    pub fn weights_files(&self) -> &[String] {
        &self.fnamew
    }
    /// Index of the first Q point above `qmin`.
    pub fn nqminm(&self) -> i32 {
        self.nqminm
    }
}

/// Local date (YYYYMMDD) and time (hhmmss) strings, as returned by the
/// Fortran `date_and_time` intrinsic.
fn now_strings() -> (String, String) {
    let secs = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs() as i64).unwrap_or(0);
    let days = secs.div_euclid(86400);
    let sod = secs.rem_euclid(86400);
    // civil-from-days (UTC)
    let z = days + 719468;
    let era = z.div_euclid(146097);
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    (
        format!("{:04}{:02}{:02}", y, m, d),
        format!("{:02}{:02}{:02}", sod / 3600, (sod % 3600) / 60, sod % 60),
    )
}
