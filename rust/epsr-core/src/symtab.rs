//! The `.EPSR.inp` keyword table — data generated from
//! `fortran_src/setup_epsr_rev1.f`, `subroutine setup_symbols_epsr`.
//!
//! Order matters: the Fortran reads the values back by position
//! (`iref=iref+1`), so entries must not be reordered.

// Generated from fortran_src/setup_epsr_rev1.f (subroutine setup_symbols_epsr).
// Each entry: (symbol, type, iflag, default, description).
pub const EPSR_SYMBOLS: &[(&str, char, i32, &str, &str)] = &[
    ("feedback", 'r', -1, "0.9", "Confidence factor - should be < 1. [0.8]"),
    ("potfac", 'r', 0, "1.0", ">0.0 to enable potential refinement, 0.0 to inhibit"),
    ("ereq", 'r', 0, "0.0", "EP amplitude[0.0]."),
    ("ereqmin", 'r', 0, "0.0", "Minimum value for ereq [0.0]. Set 0 to ignore."),
    ("ereqmax", 'r', 0, "0.0", "Maximum value for ereq [0.0]. Set 0 to ignore."),
    ("ereqstep", 'r', 0, "0.0", "Step size and unit size for ereq control [0] (Try 1 to initiate control)"),
    ("thresh", 'r', 0, "0.2", "Control baseline [0.2]"),
    ("bias", 'r', 0, "2", "Controls the bias on the steps in ereq. 0 means unbiassed steps. [2]"),
    ("sfreq", 'i', 0, "50 1", "Sampling frequency for trend line, and no. of times kT for automatic ereq to begin [50 1]"),
    ("rspcrmin", 'r', 0, "1.00", "Minimum distance for calculating the R-space coefficient [1.00]"),
    ("rspcfrac", 'r', 0, "0.5", "Fraction of R-space coefficient in control level [0.2]"),
    ("num_threds", 'i', 0, "0", "No. of parallel threads to be used, (0 to let program choose) [0]"),
    ("porodpwr", 'r', 0, "0 0", "Q-value and power for Porod scattering [0 0]"),
    ("nmolcell", 'i', 0, "10", "Average number of molecules in a cell [10]"),
    ("control-p", 'c', 0, "none 0 0", "Control pressure [a, b and/or c axes control plus value plus step, none 0 0]"),
    ("ref_intra", 'r', 0, "0 0", "Weighting on EP and Coulomb terms for intra-molecular structure [0 0]"),
    ("sizefactor", 'r', 0, "1.0 0.9 0.0", "Multiplying factor for box dimension, decline rate, and threshold. [1.0 0.9 0.0]"),
    ("nq", 'i', -1, "600", "Number of Q values. [600]"),
    ("qstep", 'r', -1, "0.05", "Size of Q step [1/A]. [0.05]"),
    ("qinscalefc", 'r', 0, "1.0 0.0", "Scale factor on input Q values and broadening power for self terms. [1.0 0.0]"),
    ("ireset", 'i', 0, "1", "1: complete reset; 2: sets the Empirical Potential to zero"),
    ("iinit", 'i', 0, "1", "Sets accumulators to zero. Recalculates r and Q. [1]"),
    ("ntimes", 'i', 0, "5", "Number of MC cycles between potential refinements. [5]"),
    ("niter", 'i', 0, "1", "Number of potential refinements before exitting. [1]"),
    ("nsumt", 'i', 0, "-1", "Number of iterations already accumulated. [-1 with reset]"),
    ("intra", 'i', 0, "100", "Number of iterations between molecule shakes. [100]"),
    ("rotfreq", 'i', 0, "5", "Number of iterations between internal rotation moves. [5]"),
    ("inter", 'i', 0, "5", "Number of iterations in running averages. [5]"),
    ("rho", 'r', 1, "0.1", "Atomic number density - will be derived from .ato file"),
    ("cellst", 'r', -1, "0.03", "Size of r step [A]. [0.03]"),
    ("rmaxgr", 'r', -1, "0.0", "Range of g(r) and F.T. (0.0 will use half the cell box) [0.0]"),
    ("ngrsamples", 'i', 0, "0", "Requested no. of origin molecules to sample g(r). (0 will use 1000 molecules) [0]"),
    ("fwhm", 'r', 0, "0.0", "Resolution width - Q independent term. [0.0]"),
    ("fwhmq", 'r', 0, "0.02", "Resolution width - Q dependent term. [0.02 for SLS]"),
    ("nsmoop", 'i', 0, "1", "1 means background subtraction is ON, 0 means OFF"),
    ("fnameato", 'c', 0, "<undefined>", "Name of .ato file"),
    ("fnamepcof", 'c', 0, "<undefined>", "Name of potential coefficients file."),
    ("revlorch", 'r', 0, "0.0 0.0", "Broadening factor in Q space. [0.0 0.0]"),
    ("qwidthqmax", 'r', 0, "0.01 0", "Broadening and maximum Q for Bragg peak calculation"),
    ("mplicities", 'i', 0, "1 1 1", "No. of unit cells along a, b and c for Bragg peak calculation"),
    ("hklqmin", 'r', 0, "0.0 0.0 0.0", "Minimum value of qhkl to be used, minimum radius for Bragg g(r), and Debye-Waller factor"),
    ("diffuse", 'i', 0, "0", "No. of unit cells along a, b and c for diffuse scattering calculation, maximum l, steps in l and m [0]"),
    ("nlatrefine", 'i', 0, "0", "No. of lattice refinement steps [0]"),
    ("latrefine", 'c', 0, "0 0 0 0 0 0", "Refinement steps for each of a, b, c, alpha, beta, gamma [0 0 0 0 0 0]"),
    ("rejrate", 'r', 0, "0.75", "Rejection rate [0.75]"),
    ("qmin", 'r', 0, "0.05 0.0", "Qmin for Fourier transforms and for potential fits. [0.05 0.0]"),
    ("roverlap", 'r', -3, "0.0", "Minimum allowed intermolecular separation between two atoms"),
    ("rminfac", 'r', -3, "0.0", "Factor to set the minimum separation between pairs. [0.0]"),
    ("rminpt", 'r', -3, "9", "Radius at which potential truncation begins."),
    ("rmaxpt", 'r', -3, "12", "Radius at which potential truncation goes to 0.0."),
    ("rbroad", 'r', 1, "0.0", "Controls potential decay. [0.0]"),
    ("expecf", 'r', 0, "0.3 0", "Potential decay for short distance repulsive term and gradient limit. [0.3 0]"),
    ("reppottype", 'c', 0, "exponential", "Repulsive potential type: exponential or harmonic. [exponential]"),
    ("addpottype", 'c', 0, "Gaussian", "Additional potential type: Gaussian or modmorse. [Gaussian]"),
    ("rcharge", 'r', 0, "0.0", "Calculates energy due to molecular polarisation. [0.0]"),
    ("power", 'r', 0, "12", "Repulsive power in Lennard-Jones potential. [12]"),
    ("ncoeffp", 'i', 1, "120", "Number of coefficients used to define the EP. Set by program."),
    ("pdmax", 'r', 1, "12.0", "Maximum distance of Empirical Potential. Set by program."),
    ("pdstep", 'r', 1, "0.1", "Spacing between coefficients in r. Set by program."),
    ("npitss", 'i', 0, "1000", "Number of steps for refining the potential. [1000]"),
    ("paccept", 'r', 1, "0.0005", "Acceptance factor for potential refinement. [0.0005]"),
    ("Gaussian", 'c', -3, "Poisson", "Select \"Gaussian\" for Gaussian representation of EP. Otherwise Poisson. [Poisson]"),
    ("psigma2", 'r', 0, "0.01", "Width for empirical potential functions. [0.01]"),
    ("refpotfac", 'r', 0, "1.0", "Factor to apply to reference potential. [1.0]"),
    ("ndata", 'i', -1, "1", "Number of data files to be fit by EPSR"),
    ("datafile", 'c', 0, "<undefined>", "Name of data file to be fit"),
    ("wtsfile", 'c', -1, "<undefined>", "Name of weights file for this data set"),
    ("nrtype", 'i', 0, "5", "Data type - see User Manual for more details"),
    ("rshmin", 'r', 0, "0.7", "Minimum radius [A] - used for background subtraction"),
    ("szeros", 'r', 0, "0.0", "Zero limit - 0 means use first data point for Q=0"),
    ("tweak", 'r', 0, "1.0", "Scaling factor for this data set. [1.0]"),
    ("efilereq", 'r', 0, "1.0", "Requested energy amplitude for this data set [1.0]"),
    ("porodadd", 'r', 0, "0.0", "Additional amplitude of Porod scattering for this dataset [0.0]"),
];

/// Positions (number of entries defined when the marker was set) of the
/// `iresetnumbers.inc` references and table boundaries.
pub const MARK_IEREQFIRST: i32 = 2;
pub const MARK_IEREQLAST: i32 = 9;
pub const MARK_IRESETREF: i32 = 21;
pub const MARK_IINITREF: i32 = 22;
pub const MARK_IQMINREF: i32 = 46;
pub const MARK_IROVERLAP: i32 = 47;
pub const MARK_IRMINFAC: i32 = 48;
pub const MARK_INCOEFFP: i32 = 57;
pub const MARK_NSYMBOLUPTODATA: i32 = 65;
pub const MARK_NSYMBOLDATA: i32 = 73;
pub const MARK_NSYMBOL: i32 = 73;
