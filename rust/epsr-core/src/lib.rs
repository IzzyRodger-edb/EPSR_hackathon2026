//! # epsr-core
//!
//! A Rust translation of the core of **EPSR** (Empirical Potential Structure
//! Refinement, A. K. Soper), version 26 (`fortran_src/epsr_standalone_rev3.f`
//! and the routines it calls).
//!
//! EPSR builds an atomistic model of a liquid or amorphous material that is
//! consistent with measured (neutron) diffraction data.  A Monte Carlo
//! simulation is run with a *reference potential* (Lennard-Jones +
//! Coulomb); the difference between the measured and simulated structure
//! factors is then used to build an *empirical potential* that drives the
//! simulation towards agreement with the data.
//!
//! The translation is deliberately literal: variable names, array indexing
//! (1-based, column-major) and the order of floating-point operations follow
//! the Fortran so that, for the same input files, the Rust engine produces
//! byte-identical output files (see `docs/VALIDATION.md`).
//!
//! ## Module map
//!
//! | module | Fortran origin | role |
//! |--------|----------------|------|
//! | [`fortran`]    | (gfortran runtime) | arrays, formatted I/O and arithmetic compatible with gfortran |
//! | [`state`]      | COMMON blocks / modules | all simulation state |
//! | [`rng`]        | `ran1` | random-number generator (state stored in `.ato`) |
//! | [`symtab`], [`symbols`] | `setup_epsr_rev1.f`, `process_symbols.f`, `read_and_save_files.f` | `.EPSR.inp` reading/writing |
//! | [`ato`]        | `ato_read_and_save.f`, `ato_manipulation_rev4.f` | `.ato` configuration files |
//! | [`pcof`]       | `read_and_save_files.f` | `.pcof` potential files |
//! | [`neighbour`]  | `cellneighbourlist.f`, `ergpairreference.f`, `minimagec.f` | periodic boundaries, neighbour lists |
//! | [`energy`]     | `molfit.f`, `ato_manipulation_rev4.f` | energies and pressure |
//! | [`montecarlo`] | `molfit.f` | Monte Carlo moves |
//! | [`structure`]  | `epsr_standalone_rev3.f`, `ato_manipulation_rev4.f` | g(r), S(Q), Fourier transforms |
//! | [`potential`]  | `epsr_standalone_rev3.f` | reference & empirical potentials, EP refinement |
//! | [`weights`]    | `epsr_standalone_rev3.f` | neutron weights |
//! | [`output`]     | `auxiliary_routines.f`, `epsr_standalone_rev3.f` | result files |
//! | [`epsr`]       | `runepsr17` | the `epsr` command (driver) |
//!
//! ## Quick start
//! ```no_run
//! use epsr_core::Epsr;
//! let mut sim = Epsr::setup(std::path::Path::new("EPSR-test"), "200DMSO400Water.EPSR.inp", false).unwrap();
//! sim.run().unwrap();
//! println!("R-factor {}", sim.summary.rfactot);
//! ```

// The translation deliberately keeps Fortran idioms where changing them could
// alter floating-point results or obscure the correspondence with the
// original source: literal constants such as `3.141592653`, explicit index
// loops, `a = a op b` forms and seemingly redundant assignment pairs.
#![allow(
    clippy::approx_constant,
    clippy::excessive_precision,
    clippy::almost_swapped,
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::too_many_arguments,
    clippy::should_implement_trait,
    clippy::needless_late_init,
    clippy::doc_overindented_list_items
)]

pub mod api;
pub mod ato;
pub mod coarsegrain;
pub mod bragg;
pub mod containerforces;
pub mod datafile;
pub mod coordconstraint;
pub mod energy;
pub mod epsr;
pub mod error;
pub mod fortran;
pub mod montecarlo;
pub mod neighbour;
pub mod output;
pub mod pcof;
pub mod potential;
pub mod pressure;
pub mod rng;
pub mod state;
pub mod structure;
pub mod symbols;
pub mod symtab;
pub mod symtab_tools;
pub mod threebody;
pub mod tools;
pub mod weights;

pub use epsr::{Epsr, Summary};
pub use error::{EpsrError, Result};
pub use state::State;

/// EPSR version string written to `.out` files (`version.inc`).
pub const VERSION_DATE: &str = "26: 2017-03-10";
