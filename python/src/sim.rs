//! `epsr._core.Simulation`: one EPSR invocation (`Epsr::setup` + `niter` iterations).
//!
//! The methods mirror the engine API (`epsr_core::api`, `epsr_core::epsr`) and keep
//! the Fortran names; the pure-Python package `epsr` provides descriptive names.

use crate::convert::{arr1, arr2, to_py};
use epsr_core::Epsr;
use numpy::{PyArray1, PyArray2, PyArray3, PyArrayMethods, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::{PyIndexError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::path::PathBuf;

#[pyclass(module = "epsr._core")]
pub struct Simulation {
    inner: Epsr,
    /// Last iteration completed by `iterate`/`run` (0 = none yet).
    last_iter: i32,
}

impl Simulation {
    fn check_dataset(&self, id: usize) -> PyResult<()> {
        if id >= self.inner.n_datasets() {
            return Err(PyIndexError::new_err(format!("dataset index {id} out of range (n_datasets = {})", self.inner.n_datasets())));
        }
        Ok(())
    }
    fn per_pair<'py, F>(&self, py: Python<'py>, n: usize, f: F) -> PyResult<Bound<'py, PyArray2<f32>>>
    where
        F: Fn(&Epsr, usize) -> &[f32],
    {
        let rows: Vec<&[f32]> = (0..self.inner.n_pairs()).map(|ip| f(&self.inner, ip)).collect();
        arr2(py, rows, n)
    }
    fn niter(&self) -> i32 {
        self.inner.ctl.niter
    }
}

/// Build a dict from `name => expression` pairs.
macro_rules! dict {
    ($py:expr; $($k:literal => $v:expr),* $(,)?) => {{
        let d = PyDict::new($py);
        $( d.set_item($k, $v)?; )*
        d
    }};
}

#[pymethods]
impl Simulation {
    /// Read all input files in `folder` (`Epsr::setup_with_lines`).
    ///
    /// `extra_lines` are processed as if appended to the `.EPSR.inp` file, e.g.
    /// `["ntimes 10", "data 2", "tweak 1.05"]`. Releases the GIL.
    #[new]
    #[pyo3(signature = (folder, inp_name, extra_lines=Vec::new(), verbose=false))]
    fn new(py: Python<'_>, folder: PathBuf, inp_name: String, extra_lines: Vec<String>, verbose: bool) -> PyResult<Self> {
        if !folder.is_dir() {
            return Err(to_py(epsr_core::EpsrError::Io(format!("simulation folder {} does not exist", folder.display()))));
        }
        let inner = py.detach(move || Epsr::setup_with_lines(&folder, &inp_name, &extra_lines, verbose)).map_err(to_py)?;
        Ok(Simulation { inner, last_iter: 0 })
    }

    /// All `niter` iterations, then `fort.35` (`Epsr::run`). Writes the result files
    /// and updates `.ato`, `.pcof`, `.EPSR.inp` in place. Releases the GIL.
    fn run(&mut self, py: Python<'_>) -> PyResult<()> {
        if self.last_iter != 0 {
            return Err(PyRuntimeError::new_err(format!(
                "run() called after iterate(); this simulation has completed {} of {} iterations",
                self.last_iter,
                self.niter()
            )));
        }
        let inner = &mut self.inner;
        py.detach(|| inner.run()).map_err(to_py)?;
        self.last_iter = self.niter();
        Ok(())
    }

    /// One iteration (`Epsr::iterate`). `iter` must be the next iteration number,
    /// `1..=niter`. After the last one the result files and `fort.35` are written,
    /// exactly as `run()` does. Returns the summary dict. Releases the GIL.
    fn iterate<'py>(&mut self, py: Python<'py>, iter: i32) -> PyResult<Bound<'py, PyDict>> {
        let niter = self.niter();
        if iter < 1 || iter > niter {
            return Err(PyValueError::new_err(format!(
                "iteration {iter} out of range 1..={niter}; create a new Simulation for the next EPSR run"
            )));
        }
        if iter != self.last_iter + 1 {
            return Err(PyValueError::new_err(format!("iterations must be run in order: expected {}, got {iter}", self.last_iter + 1)));
        }
        let inner = &mut self.inner;
        py.detach(|| -> epsr_core::Result<()> {
            inner.iterate(iter)?;
            if iter == niter {
                inner.write_fort35()?;
            }
            Ok(())
        })
        .map_err(to_py)?;
        self.last_iter = iter;
        self.summary(py)
    }

    /// Write `fort.35` (energy trace). Called automatically by `run()` and by the last `iterate()`.
    fn write_fort35(&self) -> PyResult<()> {
        self.inner.write_fort35().map_err(to_py)
    }

    /// Number of the last completed iteration (0 before the first).
    #[getter]
    fn last_iter(&self) -> i32 {
        self.last_iter
    }

    /// Results of the most recent iteration, with the Fortran names.
    #[getter]
    fn summary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let s = &self.inner.summary;
        Ok(dict!(py;
            "iter" => s.iter, "foldsv" => s.foldsv, "foldsb" => s.foldsb, "energd" => s.energd,
            "pressd" => s.pressd, "energabs" => s.energabs, "rfactot" => s.rfactot, "rfacraw" => s.rfacraw,
            "rfacd" => s.rfacd.clone(), "epcup" => s.epcup, "epccoeff" => s.epccoeff, "sentb" => s.sentb,
            "chisq" => s.chisq,
        ))
    }

    /// Run-control values read from the input (`Control`).
    #[getter]
    fn control<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let c = &self.inner.ctl;
        let p = &c.pressure;
        let pressure = dict!(py;
            "control_pressure" => p.control_pressure, "pacontrol" => p.pacontrol, "pbcontrol" => p.pbcontrol,
            "pccontrol" => p.pccontrol, "ppcontrol" => p.ppcontrol, "control_p_value" => p.control_p_value,
            "control_p_step" => p.control_p_step, "npcontrols" => p.npcontrols,
        );
        Ok(dict!(py;
            "nq" => c.nq, "qstep" => c.qstep, "ireset" => c.ireset, "iinit" => c.iinit, "ntimes" => c.ntimes,
            "niter" => c.niter, "intra" => c.intra, "rotfreq" => c.rotfreq, "inter" => c.inter, "nisumm" => c.nisumm,
            "ngrsamples" => c.ngrsamples, "fwhm" => c.fwhm, "fwhmq" => c.fwhmq, "nsmoop" => c.nsmoop,
            "pressreq" => c.pressreq, "ereqmin" => c.ereqmin, "ereqmax" => c.ereqmax, "ereqstepread" => c.ereqstepread,
            "epcthresh" => c.epcthresh, "ereqbias" => c.ereqbias, "navmaxparab" => c.navmaxparab,
            "nparabstartfactor" => c.nparabstartfactor, "epcrmin" => c.epcrmin, "epcfrac" => c.epcfrac,
            "num_threads" => c.num_threads, "nmolcellin" => c.nmolcellin, "thresholde" => c.thresholde,
            "stepmisave" => c.stepmisave, "steprisave" => c.steprisave, "rmaxgrsave" => c.rmaxgrsave,
            "irefsizefactor" => c.irefsizefactor, "porodq" => c.porodq, "porodpwr" => c.porodpwr,
            "nporodq" => c.nporodq, "pressure" => pressure,
        ))
    }

    /// Parameters of each dataset (`data n` blocks), as dicts.
    #[getter]
    fn datasets<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let l = PyList::empty(py);
        for d in &self.inner.data {
            l.append(dict!(py;
                "infili" => d.infili.clone(), "fnamedata" => d.fnamedata.clone(), "fnamewts" => d.fnamewts.clone(),
                "nrtype" => d.nrtype, "rshmin" => d.rshmin, "szeros" => d.szeros, "tweak" => d.tweak,
                "efilereq" => d.efilereq, "porodadd" => d.porodadd, "nqmind" => d.nqmind, "nqmaxd" => d.nqmaxd,
                "porodamp" => d.porodamp,
            ))?;
        }
        Ok(l)
    }

    #[getter]
    fn n_datasets(&self) -> usize {
        self.inner.n_datasets()
    }
    #[getter]
    fn n_pairs(&self) -> usize {
        self.inner.n_pairs()
    }
    #[getter]
    fn nq(&self) -> usize {
        self.inner.nq()
    }
    #[getter]
    fn nr(&self) -> usize {
        self.inner.nr()
    }
    #[getter]
    fn n_molecules(&self) -> usize {
        self.inner.n_molecules()
    }
    #[getter]
    fn n_atoms(&self) -> usize {
        self.inner.n_atoms()
    }
    /// Atomic number density (atoms/Å³).
    #[getter]
    fn rho(&self) -> f32 {
        self.inner.st.dist.rho
    }

    fn atom_types(&self) -> Vec<String> {
        self.inner.atom_types()
    }
    fn pair_labels(&self) -> Vec<(String, String)> {
        self.inner.pair_labels()
    }
    fn atom_labels(&self) -> Vec<String> {
        self.inner.atom_labels()
    }
    /// For each molecule, the 0-based half-open range `(start, stop)` of its atoms.
    fn molecule_atoms(&self) -> Vec<(usize, usize)> {
        self.inner.molecule_atoms().into_iter().map(|r| (r.start, r.end)).collect()
    }
    fn weights_files(&self) -> Vec<String> {
        self.inner.weights_files().to_vec()
    }

    fn q_grid<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f32>> {
        arr1(py, self.inner.q_grid())
    }
    fn r_grid<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f32>> {
        arr1(py, self.inner.r_grid())
    }

    /// `(n_pairs, nr)` running-average intermolecular partial g(r).
    fn gr_inter<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.per_pair(py, self.inner.nr(), |e, ip| e.gr_inter(ip))
    }
    /// `(n_pairs, nr)` running-average intramolecular partial g(r).
    fn gr_intra<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.per_pair(py, self.inner.nr(), |e, ip| e.gr_intra(ip))
    }
    /// `(n_pairs, nq)` running-average intermolecular partial S(Q).
    fn sq_inter<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.per_pair(py, self.inner.nq(), |e, ip| e.sq_inter(ip))
    }
    /// `(n_pairs, nq)` running-average intramolecular partial S(Q).
    fn sq_intra<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.per_pair(py, self.inner.nq(), |e, ip| e.sq_intra(ip))
    }

    /// Measured F(Q) of dataset `id`, interpolated onto the Q grid.
    fn data_fq<'py>(&self, py: Python<'py>, id: usize) -> PyResult<Bound<'py, PyArray1<f32>>> {
        self.check_dataset(id)?;
        Ok(arr1(py, self.inner.data_fq(id)))
    }
    /// Errors of the measured F(Q) of dataset `id` on the Q grid (`errs`).
    fn data_errors<'py>(&self, py: Python<'py>, id: usize) -> PyResult<Bound<'py, PyArray1<f32>>> {
        self.check_dataset(id)?;
        let nq = self.inner.nq();
        Ok(arr1(py, &self.inner.errs.col(id + 1)[..nq]))
    }
    /// Simulated total F(Q) of dataset `id`.
    fn simulated_fq<'py>(&self, py: Python<'py>, id: usize) -> PyResult<Bound<'py, PyArray1<f32>>> {
        self.check_dataset(id)?;
        Ok(PyArray1::from_vec(py, self.inner.simulated_fq(id)))
    }
    /// First and last Q-grid index (0-based, inclusive) fitted for dataset `id`.
    fn data_q_range(&self, id: usize) -> PyResult<(usize, usize)> {
        self.check_dataset(id)?;
        Ok(self.inner.data_q_range(id))
    }

    /// `(2, n_pairs, nr)` histograms of the *current* configuration (`wav`):
    /// index 0 = intermolecular, 1 = intramolecular.
    fn current_gr<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray3<f32>>> {
        let (np, n) = (self.inner.n_pairs(), self.inner.nr());
        let w = &self.inner.wav;
        let v: Vec<Vec<Vec<f32>>> = (1..=2).map(|k| (1..=np).map(|ip| (1..=n).map(|i| w[(k, i, ip)]).collect()).collect()).collect();
        if np == 0 {
            return PyArray1::<f32>::from_vec(py, vec![]).reshape([2, 0, n]);
        }
        PyArray3::from_vec3(py, &v).map_err(|e| PyValueError::new_err(e.to_string()))
    }
    /// `(2, n_pairs, nq)` partial structure factors of the *current* configuration (`qwind`).
    fn current_sq<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray3<f32>>> {
        let (np, n) = (self.inner.n_pairs(), self.inner.nq());
        let w = &self.inner.qwind;
        let v: Vec<Vec<Vec<f32>>> = (1..=2).map(|k| (1..=np).map(|ip| (1..=n).map(|i| w[(k, i, ip)]).collect()).collect()).collect();
        if np == 0 {
            return PyArray1::<f32>::from_vec(py, vec![]).reshape([2, 0, n]);
        }
        PyArray3::from_vec3(py, &v).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// `(n_pairs, nr)` empirical potential (kJ/mol, `xavs`).
    fn empirical_potential<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.per_pair(py, self.inner.nr(), |e, ip| e.empirical_potential(ip))
    }
    /// `(n_pairs, nr)` reference potential (Lennard-Jones + Coulomb, kJ/mol, `ecoul`).
    fn reference_potential<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.per_pair(py, self.inner.nr(), |e, ip| e.reference_potential(ip))
    }
    /// `(n_pairs, nr)` total potential, reference + empirical (kJ/mol, `xavtot`).
    fn total_potential<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let nr = self.inner.nr();
        self.per_pair(py, nr, |e, ip| &e.st.pot.xavtot.col(ip + 1)[..nr])
    }

    /// `(n_pairs, ncoeffp)` EP coefficients `pcoefp`; column 0 is the refinement flag.
    fn ep_coefficients<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let n = self.inner.st.pot.ncoeffp.max(0) as usize;
        self.per_pair(py, n, |e, ip| e.ep_coefficients(ip))
    }

    /// For each pair, the `.pcof` coefficient column it uses (`irefpcof`, 1-based).
    /// Pairs with the same value share one set of EP coefficients.
    fn ep_columns(&self) -> Vec<i32> {
        (1..=self.inner.n_pairs()).map(|ip| self.inner.st.pot.irefpcof[ip]).collect()
    }

    /// Replace the EP coefficients (`(n_pairs, ncoeffp)`, column 0 = refinement
    /// flag) and rebuild the empirical potential (`generate_ep(false)`).
    ///
    /// Pairs matched to the same `.pcof` entry share one coefficient column; for
    /// those the row written last (highest pair index) wins.
    fn set_ep_coefficients(&mut self, coef: PyReadonlyArray2<'_, f32>) -> PyResult<()> {
        let np = self.inner.n_pairs();
        let n = self.inner.st.pot.ncoeffp.max(0) as usize;
        if coef.shape() != [np, n] {
            return Err(PyValueError::new_err(format!("expected coefficients of shape ({np}, {n}), got {:?}", coef.shape())));
        }
        let a = coef.as_array();
        let pot = &mut self.inner.st.pot;
        for ip in 0..np {
            let irefp = pot.irefpcof[ip + 1];
            for ic in 0..n {
                pot.pcoefp[(ic + 1, irefp as usize)] = a[[ip, ic]];
            }
        }
        self.inner.st.generate_ep(false).map_err(to_py)
    }

    /// Multiply every EP coefficient by `factor` (`multiply_pot_coeffs`) and rebuild the EP.
    fn scale_ep(&mut self, factor: f32) -> PyResult<()> {
        self.inner.st.multiply_pot_coeffs(factor);
        self.inner.st.generate_ep(false).map_err(to_py)
    }

    /// Rebuild the empirical potential tables from the coefficients (`generate_ep(false)`).
    fn regenerate_ep(&mut self) -> PyResult<()> {
        self.inner.st.generate_ep(false).map_err(to_py)
    }

    /// Snapshot of the potential parameters (Fortran names).
    fn pot_params<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let p = &self.inner.st.pot;
        Ok(dict!(py;
            "absolute_energy" => p.absolute_energy, "potfac" => p.potfac, "potfacamp" => p.potfacamp,
            "refpotfac" => p.refpotfac, "feedback" => p.feedback, "ncoeffp" => p.ncoeffp, "pdmin" => p.pdmin,
            "pdstep" => p.pdstep, "pdmax" => p.pdmax, "psigma2" => p.psigma2, "gaussian" => p.gaussian,
        ))
    }

    /// Per atom type: label, Lennard-Jones epsilon (kJ/mol) and sigma (Å), charge (e),
    /// mass (amu) and atom fraction.
    fn type_params<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let m = &self.inner.st.mol;
        let l = PyList::empty(py);
        for it in 1..=m.ntype {
            l.append(dict!(py;
                "label" => m.type_[it].clone(), "epsilon" => m.epscor[it], "sigma" => m.sigcor[it],
                "charge" => m.charge[it], "mass" => m.amassa[it], "fraction" => m.cfrac[it],
            ))?;
        }
        Ok(l)
    }

    /// Cubic box side (Å).
    fn box_size(&self) -> f32 {
        self.inner.box_size()
    }
    /// `(n_atoms, 3)` absolute atom positions (Å).
    fn atom_positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let pos = self.inner.atom_positions();
        let rows: Vec<&[f32]> = pos.iter().map(|p| &p[..]).collect();
        arr2(py, rows, 3)
    }

    /// Current value of an input keyword (global, or of 0-based `dataset`).
    #[pyo3(signature = (name, dataset=None))]
    fn get_parameter(&self, name: &str, dataset: Option<usize>) -> PyResult<Option<String>> {
        if let Some(d) = dataset {
            self.check_dataset(d)?;
        }
        Ok(self.inner.get_parameter(name, dataset))
    }
}
