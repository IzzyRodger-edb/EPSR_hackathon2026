//! `epsr._core.Simulation`: one EPSR invocation (`Epsr::setup` + `niter` iterations).
//!
//! The methods mirror the engine API (`epsr_core::api`, `epsr_core::epsr`) and keep
//! the Fortran names; the pure-Python package `epsr` provides descriptive names.
//!
//! Engine arrays are 1-based and column-major (`F1`/`F2`/`F3`); the helpers here
//! copy the used part of them into 0-based numpy arrays.

use crate::convert::{array_1d, array_2d, to_py};
use epsr_core::fortran::F3;
use epsr_core::Epsr;
use numpy::{
    PyArray1, PyArray2, PyArray3, PyArrayMethods, PyReadonlyArray2, PyUntypedArrayMethods,
};
use pyo3::exceptions::{PyIndexError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::path::PathBuf;

/// First index of the `wav`/`qwind` arrays: 1 = intermolecular, 2 = intramolecular.
const INTER_AND_INTRA: std::ops::RangeInclusive<usize> = 1..=2;

#[pyclass(module = "epsr._core")]
pub struct Simulation {
    engine: Epsr,
    /// Last iteration completed by `iterate`/`run` (0 = none yet).
    last_iteration: i32,
}

impl Simulation {
    fn check_dataset(&self, dataset: usize) -> PyResult<()> {
        let n_datasets = self.engine.n_datasets();
        if dataset >= n_datasets {
            return Err(PyIndexError::new_err(format!(
                "dataset index {dataset} out of range (n_datasets = {n_datasets})"
            )));
        }
        Ok(())
    }

    /// A `(n_pairs, n_points)` array with one row per atom-type pair.
    fn one_row_per_pair<'py, RowOfPair>(
        &self,
        py: Python<'py>,
        n_points: usize,
        row_of_pair: RowOfPair,
    ) -> PyResult<Bound<'py, PyArray2<f32>>>
    where
        RowOfPair: Fn(&Epsr, usize) -> &[f32],
    {
        let rows: Vec<&[f32]> = (0..self.engine.n_pairs())
            .map(|pair| row_of_pair(&self.engine, pair))
            .collect();
        array_2d(py, rows, n_points)
    }

    /// Copy an engine array `values(1|2, point, pair)` into a `(2, n_pairs, n_points)` array.
    fn inter_intra_by_pair<'py>(
        &self,
        py: Python<'py>,
        values: &F3<f32>,
        n_points: usize,
    ) -> PyResult<Bound<'py, PyArray3<f32>>> {
        let n_pairs = self.engine.n_pairs();
        if n_pairs == 0 {
            return PyArray1::<f32>::from_vec(py, vec![]).reshape([2, 0, n_points]);
        }
        let nested: Vec<Vec<Vec<f32>>> = INTER_AND_INTRA
            .map(|kind| {
                (1..=n_pairs)
                    .map(|pair| {
                        (1..=n_points)
                            .map(|point| values[(kind, point, pair)])
                            .collect()
                    })
                    .collect()
            })
            .collect();
        PyArray3::from_vec3(py, &nested).map_err(|error| PyValueError::new_err(error.to_string()))
    }

    fn niter(&self) -> i32 {
        self.engine.ctl.niter
    }

    /// Number of EP coefficients per pair, including the refinement flag (`ncoeffp`).
    fn n_coefficients(&self) -> usize {
        self.engine.st.pot.ncoeffp.max(0) as usize
    }
}

/// Build a Python dict from `"name" => expression` pairs.
macro_rules! dict {
    ($py:expr; $($key:literal => $value:expr),* $(,)?) => {{
        let dictionary = PyDict::new($py);
        $( dictionary.set_item($key, $value)?; )*
        dictionary
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
    fn new(
        py: Python<'_>,
        folder: PathBuf,
        inp_name: String,
        extra_lines: Vec<String>,
        verbose: bool,
    ) -> PyResult<Self> {
        if !folder.is_dir() {
            return Err(to_py(epsr_core::EpsrError::Io(format!(
                "simulation folder {} does not exist",
                folder.display()
            ))));
        }
        let engine = py
            .detach(move || Epsr::setup_with_lines(&folder, &inp_name, &extra_lines, verbose))
            .map_err(to_py)?;
        Ok(Simulation {
            engine,
            last_iteration: 0,
        })
    }

    /// All `niter` iterations, then `fort.35` (`Epsr::run`). Writes the result files
    /// and updates `.ato`, `.pcof`, `.EPSR.inp` in place. Releases the GIL.
    fn run(&mut self, py: Python<'_>) -> PyResult<()> {
        if self.last_iteration != 0 {
            return Err(PyRuntimeError::new_err(format!(
                "run() called after iterate(); this simulation has completed {} of {} iterations",
                self.last_iteration,
                self.niter()
            )));
        }
        let engine = &mut self.engine;
        py.detach(|| engine.run()).map_err(to_py)?;
        self.last_iteration = self.niter();
        Ok(())
    }

    /// One iteration (`Epsr::iterate`). `iteration` must be the next iteration number,
    /// `1..=niter`. After the last one the result files and `fort.35` are written,
    /// exactly as `run()` does. Returns the summary dict. Releases the GIL.
    fn iterate<'py>(&mut self, py: Python<'py>, iteration: i32) -> PyResult<Bound<'py, PyDict>> {
        let niter = self.niter();
        if iteration < 1 || iteration > niter {
            return Err(PyValueError::new_err(format!(
                "iteration {iteration} out of range 1..={niter}; create a new Simulation for the next EPSR run"
            )));
        }
        let expected_iteration = self.last_iteration + 1;
        if iteration != expected_iteration {
            return Err(PyValueError::new_err(format!(
                "iterations must be run in order: expected {expected_iteration}, got {iteration}"
            )));
        }
        let engine = &mut self.engine;
        py.detach(|| -> epsr_core::Result<()> {
            engine.iterate(iteration)?;
            if iteration == niter {
                engine.write_fort35()?;
            }
            Ok(())
        })
        .map_err(to_py)?;
        self.last_iteration = iteration;
        self.summary(py)
    }

    /// Write `fort.35` (energy trace). Called automatically by `run()` and by the last `iterate()`.
    fn write_fort35(&self) -> PyResult<()> {
        self.engine.write_fort35().map_err(to_py)
    }

    /// Number of the last completed iteration (0 before the first).
    #[getter]
    fn last_iter(&self) -> i32 {
        self.last_iteration
    }

    /// Results of the most recent iteration, with the Fortran names.
    #[getter]
    fn summary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let summary = &self.engine.summary;
        Ok(dict!(py;
            "iter" => summary.iter,
            "foldsv" => summary.foldsv,
            "foldsb" => summary.foldsb,
            "energd" => summary.energd,
            "pressd" => summary.pressd,
            "energabs" => summary.energabs,
            "rfactot" => summary.rfactot,
            "rfacraw" => summary.rfacraw,
            "rfacd" => summary.rfacd.clone(),
            "epcup" => summary.epcup,
            "epccoeff" => summary.epccoeff,
            "sentb" => summary.sentb,
            "chisq" => summary.chisq,
        ))
    }

    /// Run-control values read from the input (`Control`).
    #[getter]
    fn control<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let control = &self.engine.ctl;
        let pressure_control = &control.pressure;
        let pressure = dict!(py;
            "control_pressure" => pressure_control.control_pressure,
            "pacontrol" => pressure_control.pacontrol,
            "pbcontrol" => pressure_control.pbcontrol,
            "pccontrol" => pressure_control.pccontrol,
            "ppcontrol" => pressure_control.ppcontrol,
            "control_p_value" => pressure_control.control_p_value,
            "control_p_step" => pressure_control.control_p_step,
            "npcontrols" => pressure_control.npcontrols,
        );
        Ok(dict!(py;
            "nq" => control.nq,
            "qstep" => control.qstep,
            "ireset" => control.ireset,
            "iinit" => control.iinit,
            "ntimes" => control.ntimes,
            "niter" => control.niter,
            "intra" => control.intra,
            "rotfreq" => control.rotfreq,
            "inter" => control.inter,
            "nisumm" => control.nisumm,
            "ngrsamples" => control.ngrsamples,
            "fwhm" => control.fwhm,
            "fwhmq" => control.fwhmq,
            "nsmoop" => control.nsmoop,
            "pressreq" => control.pressreq,
            "ereqmin" => control.ereqmin,
            "ereqmax" => control.ereqmax,
            "ereqstepread" => control.ereqstepread,
            "epcthresh" => control.epcthresh,
            "ereqbias" => control.ereqbias,
            "navmaxparab" => control.navmaxparab,
            "nparabstartfactor" => control.nparabstartfactor,
            "epcrmin" => control.epcrmin,
            "epcfrac" => control.epcfrac,
            "num_threads" => control.num_threads,
            "nmolcellin" => control.nmolcellin,
            "thresholde" => control.thresholde,
            "stepmisave" => control.stepmisave,
            "steprisave" => control.steprisave,
            "rmaxgrsave" => control.rmaxgrsave,
            "irefsizefactor" => control.irefsizefactor,
            "porodq" => control.porodq,
            "porodpwr" => control.porodpwr,
            "nporodq" => control.nporodq,
            "pressure" => pressure,
        ))
    }

    /// Parameters of each dataset (`data n` blocks), as dicts.
    #[getter]
    fn datasets<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let datasets = PyList::empty(py);
        for dataset in &self.engine.data {
            datasets.append(dict!(py;
                "infili" => dataset.infili.clone(),
                "fnamedata" => dataset.fnamedata.clone(),
                "fnamewts" => dataset.fnamewts.clone(),
                "nrtype" => dataset.nrtype,
                "rshmin" => dataset.rshmin,
                "szeros" => dataset.szeros,
                "tweak" => dataset.tweak,
                "efilereq" => dataset.efilereq,
                "porodadd" => dataset.porodadd,
                "nqmind" => dataset.nqmind,
                "nqmaxd" => dataset.nqmaxd,
                "porodamp" => dataset.porodamp,
            ))?;
        }
        Ok(datasets)
    }

    #[getter]
    fn n_datasets(&self) -> usize {
        self.engine.n_datasets()
    }
    #[getter]
    fn n_pairs(&self) -> usize {
        self.engine.n_pairs()
    }
    #[getter]
    fn nq(&self) -> usize {
        self.engine.nq()
    }
    #[getter]
    fn nr(&self) -> usize {
        self.engine.nr()
    }
    #[getter]
    fn n_molecules(&self) -> usize {
        self.engine.n_molecules()
    }
    #[getter]
    fn n_atoms(&self) -> usize {
        self.engine.n_atoms()
    }
    /// Atomic number density (atoms/Å³).
    #[getter]
    fn rho(&self) -> f32 {
        self.engine.st.dist.rho
    }

    fn atom_types(&self) -> Vec<String> {
        self.engine.atom_types()
    }
    fn pair_labels(&self) -> Vec<(String, String)> {
        self.engine.pair_labels()
    }
    fn atom_labels(&self) -> Vec<String> {
        self.engine.atom_labels()
    }
    /// For each molecule, the 0-based half-open range `(start, stop)` of its atoms.
    fn molecule_atoms(&self) -> Vec<(usize, usize)> {
        self.engine
            .molecule_atoms()
            .into_iter()
            .map(|atom_range| (atom_range.start, atom_range.end))
            .collect()
    }
    fn weights_files(&self) -> Vec<String> {
        self.engine.weights_files().to_vec()
    }

    fn q_grid<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f32>> {
        array_1d(py, self.engine.q_grid())
    }
    fn r_grid<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f32>> {
        array_1d(py, self.engine.r_grid())
    }

    /// `(n_pairs, nr)` running-average intermolecular partial g(r).
    fn gr_inter<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.one_row_per_pair(py, self.engine.nr(), |engine, pair| engine.gr_inter(pair))
    }
    /// `(n_pairs, nr)` running-average intramolecular partial g(r).
    fn gr_intra<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.one_row_per_pair(py, self.engine.nr(), |engine, pair| engine.gr_intra(pair))
    }
    /// `(n_pairs, nq)` running-average intermolecular partial S(Q).
    fn sq_inter<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.one_row_per_pair(py, self.engine.nq(), |engine, pair| engine.sq_inter(pair))
    }
    /// `(n_pairs, nq)` running-average intramolecular partial S(Q).
    fn sq_intra<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.one_row_per_pair(py, self.engine.nq(), |engine, pair| engine.sq_intra(pair))
    }

    /// Measured F(Q) of a dataset, interpolated onto the Q grid.
    fn data_fq<'py>(&self, py: Python<'py>, dataset: usize) -> PyResult<Bound<'py, PyArray1<f32>>> {
        self.check_dataset(dataset)?;
        Ok(array_1d(py, self.engine.data_fq(dataset)))
    }
    /// Errors of the measured F(Q) of a dataset on the Q grid (`errs`).
    fn data_errors<'py>(
        &self,
        py: Python<'py>,
        dataset: usize,
    ) -> PyResult<Bound<'py, PyArray1<f32>>> {
        self.check_dataset(dataset)?;
        let nq = self.engine.nq();
        let fortran_column = dataset + 1;
        Ok(array_1d(py, &self.engine.errs.col(fortran_column)[..nq]))
    }
    /// Simulated total F(Q) of a dataset.
    fn simulated_fq<'py>(
        &self,
        py: Python<'py>,
        dataset: usize,
    ) -> PyResult<Bound<'py, PyArray1<f32>>> {
        self.check_dataset(dataset)?;
        Ok(PyArray1::from_vec(py, self.engine.simulated_fq(dataset)))
    }
    /// First and last Q-grid index (0-based, inclusive) fitted for a dataset.
    fn data_q_range(&self, dataset: usize) -> PyResult<(usize, usize)> {
        self.check_dataset(dataset)?;
        Ok(self.engine.data_q_range(dataset))
    }

    /// `(2, n_pairs, nr)` histograms of the *current* configuration (`wav`):
    /// index 0 = intermolecular, 1 = intramolecular.
    fn current_gr<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray3<f32>>> {
        self.inter_intra_by_pair(py, &self.engine.wav, self.engine.nr())
    }
    /// `(2, n_pairs, nq)` partial structure factors of the *current* configuration (`qwind`).
    fn current_sq<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray3<f32>>> {
        self.inter_intra_by_pair(py, &self.engine.qwind, self.engine.nq())
    }

    /// `(n_pairs, nr)` empirical potential (kJ/mol, `xavs`).
    fn empirical_potential<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.one_row_per_pair(py, self.engine.nr(), |engine, pair| {
            engine.empirical_potential(pair)
        })
    }
    /// `(n_pairs, nr)` reference potential (Lennard-Jones + Coulomb, kJ/mol, `ecoul`).
    fn reference_potential<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.one_row_per_pair(py, self.engine.nr(), |engine, pair| {
            engine.reference_potential(pair)
        })
    }
    /// `(n_pairs, nr)` total potential, reference + empirical (kJ/mol, `xavtot`).
    fn total_potential<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let nr = self.engine.nr();
        self.one_row_per_pair(py, nr, |engine, pair| {
            &engine.st.pot.xavtot.col(pair + 1)[..nr]
        })
    }

    /// `(n_pairs, ncoeffp)` EP coefficients `pcoefp`; column 0 is the refinement flag.
    fn ep_coefficients<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.one_row_per_pair(py, self.n_coefficients(), |engine, pair| {
            engine.ep_coefficients(pair)
        })
    }

    /// For each pair, the `.pcof` coefficient column it uses (`irefpcof`, 1-based).
    /// Pairs with the same value share one set of EP coefficients.
    fn ep_columns(&self) -> Vec<i32> {
        (1..=self.engine.n_pairs())
            .map(|fortran_pair| self.engine.st.pot.irefpcof[fortran_pair])
            .collect()
    }

    /// Replace the EP coefficients (`(n_pairs, ncoeffp)`, column 0 = refinement
    /// flag) and rebuild the empirical potential (`generate_ep(false)`).
    ///
    /// Pairs matched to the same `.pcof` entry share one coefficient column; for
    /// those the row written last (highest pair index) wins.
    fn set_ep_coefficients(&mut self, coefficients: PyReadonlyArray2<'_, f32>) -> PyResult<()> {
        let n_pairs = self.engine.n_pairs();
        let n_coefficients = self.n_coefficients();
        if coefficients.shape() != [n_pairs, n_coefficients] {
            return Err(PyValueError::new_err(format!(
                "expected coefficients of shape ({n_pairs}, {n_coefficients}), got {:?}",
                coefficients.shape()
            )));
        }
        let new_values = coefficients.as_array();
        let potential = &mut self.engine.st.pot;
        for pair in 0..n_pairs {
            let pcof_column = potential.irefpcof[pair + 1] as usize;
            for coefficient in 0..n_coefficients {
                potential.pcoefp[(coefficient + 1, pcof_column)] = new_values[[pair, coefficient]];
            }
        }
        self.engine.st.generate_ep(false).map_err(to_py)
    }

    /// Multiply every EP coefficient by `factor` (`multiply_pot_coeffs`) and rebuild the EP.
    fn scale_ep(&mut self, factor: f32) -> PyResult<()> {
        self.engine.st.multiply_pot_coeffs(factor);
        self.engine.st.generate_ep(false).map_err(to_py)
    }

    /// Rebuild the empirical potential tables from the coefficients (`generate_ep(false)`).
    fn regenerate_ep(&mut self) -> PyResult<()> {
        self.engine.st.generate_ep(false).map_err(to_py)
    }

    /// Snapshot of the potential parameters (Fortran names).
    fn pot_params<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let potential = &self.engine.st.pot;
        Ok(dict!(py;
            "absolute_energy" => potential.absolute_energy,
            "potfac" => potential.potfac,
            "potfacamp" => potential.potfacamp,
            "refpotfac" => potential.refpotfac,
            "feedback" => potential.feedback,
            "ncoeffp" => potential.ncoeffp,
            "pdmin" => potential.pdmin,
            "pdstep" => potential.pdstep,
            "pdmax" => potential.pdmax,
            "psigma2" => potential.psigma2,
            "gaussian" => potential.gaussian,
        ))
    }

    /// Per atom type: label, Lennard-Jones epsilon (kJ/mol) and sigma (Å), charge (e),
    /// mass (amu) and atom fraction.
    fn type_params<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let molecules = &self.engine.st.mol;
        let atom_types = PyList::empty(py);
        for atom_type in 1..=molecules.ntype {
            atom_types.append(dict!(py;
                "label" => molecules.type_[atom_type].clone(),
                "epsilon" => molecules.epscor[atom_type],
                "sigma" => molecules.sigcor[atom_type],
                "charge" => molecules.charge[atom_type],
                "mass" => molecules.amassa[atom_type],
                "fraction" => molecules.cfrac[atom_type],
            ))?;
        }
        Ok(atom_types)
    }

    /// Cubic box side (Å).
    fn box_size(&self) -> f32 {
        self.engine.box_size()
    }
    /// `(n_atoms, 3)` absolute atom positions (Å).
    fn atom_positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let positions = self.engine.atom_positions();
        let rows: Vec<&[f32]> = positions.iter().map(|position| &position[..]).collect();
        array_2d(py, rows, 3)
    }

    /// Current value of an input keyword (global, or of 0-based `dataset`).
    #[pyo3(signature = (name, dataset=None))]
    fn get_parameter(&self, name: &str, dataset: Option<usize>) -> PyResult<Option<String>> {
        if let Some(dataset_index) = dataset {
            self.check_dataset(dataset_index)?;
        }
        Ok(self.engine.get_parameter(name, dataset))
    }
}
