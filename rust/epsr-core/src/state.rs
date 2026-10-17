//! Global simulation state.
//!
//! The Fortran program keeps its state in COMMON blocks (`*.inc` files) and
//! module variables.  Each block is represented here by a struct whose field
//! names are identical to the Fortran variable names, so that a translated
//! statement such as `xavs(iz,icref)` reads `self.pot.xavs[(iz, icref)]`.
//!
//! | Rust struct      | Fortran source |
//! |------------------|----------------|
//! | [`Lattice`]      | `lattice.inc` (`/lattice_r/`, `/lattice_l/`) |
//! | [`Mol`]          | `molecom.inc`, `clustercom.inc`, module `molecom_mod` |
//! | [`Dist`]         | `epsr16.inc` (`/dista/`, `/xaver/`, `/idista/`, `/xaver8/`) |
//! | [`Pot`]          | `potential_par.inc` |
//! | [`Weights`]      | `weights.inc` |
//! | [`InputPar`]     | `epsrinput_par.inc` |
//! | [`Symbols`]      | `symbols.inc`, `iresetnumbers.inc`, `commands.inc`, `directory_name.inc` |
//! | [`CellList`]     | module `cellneighbourlist` |
//! | [`ErgPairs`]     | module `ergpairreference` |
//! | [`Timer`]        | module `timer` (wall-clock bookkeeping only) |
//!
//! Array dimensions follow `epsr_dimension.inc`, except that per-atom arrays
//! are sized from the configuration actually read (the Fortran uses fixed
//! maxima of 150000 atoms).

use crate::fortran::{F1, F2, F3};
use crate::rng::Ran1;
use std::path::PathBuf;

// --- epsr_dimension.inc -----------------------------------------------------
/// Maximum number of points in an input data file.
pub const MQD: usize = 5000;
/// Maximum number of r-space points.
pub const MPZ: usize = 3000;
/// Maximum number of Q-space points.
pub const MPTS: usize = 2001;
/// Maximum number of data sets.
pub const MFILE: usize = 30;
/// Maximum number of atom types.
pub const MCOMP: usize = 40;
/// Maximum number of atom-type pairs.
pub const MPAIRS: usize = MCOMP * (MCOMP + 1) / 2;
/// Maximum number of bonds per atom ("near").
pub const NEAR: usize = 200;
/// Maximum number of molecules.
pub const MMOL: usize = 50000;
// --- potential_par.inc ------------------------------------------------------
/// Maximum number of empirical-potential coefficients.
pub const MCOEFF: usize = 1000;
/// Maximum number of additional Gaussian potentials per pair.
pub const MPGAUSS: usize = 5;
// --- weights.inc --------------------------------------------------------------
pub const MQINV: usize = 15;
// --- symbols.inc --------------------------------------------------------------
pub const MSYMBOL: usize = MPAIRS + 36 + 10 * (MFILE - 1);

/// `lattice.inc`: simulation box geometry.  After `normalise_coordinates` the
/// cell vectors are in units of `dmax` (box-fraction units).
#[derive(Clone, Debug, Default)]
pub struct Lattice {
    pub alat: f32,
    pub blat: f32,
    pub clat: f32,
    pub phib: f32,
    pub thetac: f32,
    pub phic: f32,
    pub alatx: f32,
    pub blatx: f32,
    pub clatx: f32,
    pub blaty: f32,
    pub claty: f32,
    pub clatz: f32,
    pub bxax: f32,
    pub cxax: f32,
    pub cyby: f32,
    pub alatc: [f32; 3],
    pub blatc: [f32; 3],
    pub clatc: [f32; 3],
    pub cubic: bool,
}

/// `molecom.inc` + `clustercom.inc` + `molecom_mod`: molecules, atoms, atom
/// types, intramolecular structure and Monte Carlo step sizes.
#[derive(Clone, Debug, Default)]
pub struct Mol {
    // ---- molecom_r
    pub tol: f32,
    pub tolout: f32,
    pub tempb: f32,
    pub tempdih: f32,
    pub tempang: f32,
    pub tempdihreq: f32,
    pub tempangreq: f32,
    pub tempangwt: f32,
    pub tempdihwt: f32,
    pub stepr: f32,
    pub stepri: f32,
    pub stepm: f32,
    pub stepmi: f32,
    /// Atom coordinates relative to the molecule centre (box-fraction units after normalisation).
    pub xa: F1<f32>,
    pub ya: F1<f32>,
    pub za: F1<f32>,
    pub amassa: F1<f32>,
    pub amassr: F1<f32>,
    /// Bond lengths `tbond(ia,j)` (negative = bond-angle constraint).
    pub tbond: F2<f32>,
    pub xs: f32,
    pub ys: f32,
    pub zs: f32,
    pub simt: f32,
    pub simtk: f32,
    pub charge: F1<f32>,
    pub qradius: F1<f32>,
    pub qcor: F1<f32>,
    pub dcor: F1<f32>,
    pub sigcor: F1<f32>,
    pub epscor: F1<f32>,
    pub dscor: F1<f32>,
    pub ecor: F1<f32>,
    pub cfrac: F1<f32>,
    pub cfraccg: F1<f32>,
    pub dcore: f32,
    pub ecore: f32,
    pub power: f32,
    pub sizefactor: f32,
    pub decline_rate: f32,
    pub dmaxmol: f32,
    pub refine_intra_fac: f32,
    pub coul_intra_fac: f32,
    pub sumcharge: f32,
    pub maxmolsize: f32,
    // ---- molecom_c
    pub alabel: F1<String>,
    pub type_: F1<String>,
    pub atsymbol: F1<String>,
    pub typefile: F1<String>,
    // ---- molecom_i
    pub ndihedralang: F1<i32>,
    /// `dihedralang(k,ig)`: required dihedral angle (deg) of DIH rotational groups.
    pub dihedralang: F2<f32>,
    /// `labrotdih(j,k,ig)`: the two atoms defining dihedral `k` of DIH group `ig`.
    pub labrotdih: F3<i32>,
    pub nbond: F1<i32>,
    pub iabond: F2<i32>,
    pub ntype: i32,
    pub itype: F1<i32>,
    pub natom: i32,
    pub natomr: i32,
    pub mref: F1<i32>,
    pub qtype: F1<i32>,
    pub qtypea: F1<i32>,
    pub npairs: i32,
    pub npairref: F2<i32>,
    pub n1extr: F1<i32>,
    pub n2extr: F1<i32>,
    pub iexchange: F1<i32>,
    pub ntypefile: i32,
    pub iatomicnotype: F1<i32>,
    pub imageatom: F1<i32>,
    // ---- molecom_l
    pub refine_intra: bool,
    pub coul_intra: bool,
    // ---- clustercom
    pub moltype: F1<String>,
    pub moltemplate: F1<String>,
    pub nmol: i32,
    pub nato: F1<i32>,
    pub ngrot: F1<i32>,
    pub naf: F1<i32>,
    pub nal: F1<i32>,
    pub nrf: F1<i32>,
    pub nrl: F1<i32>,
    pub mtype: F1<i32>,
    pub nmoltype: i32,
    pub nmovei: F2<i32>,
    pub ninter: F2<i32>,
    pub nwhole: F2<i32>,
    pub nmovem: F2<i32>,
    pub nmoveit: [i32; 2],
    pub nintert: [i32; 2],
    pub nwholet: [i32; 2],
    pub nmovemt: [i32; 2],
    pub nstepreject: F1<i32>,
    pub tetheredatom: F1<i32>,
    pub moltetheredatom: F1<i32>,
    pub tethered: F1<bool>,
    pub moltethered: F1<bool>,
    /// Molecule centres (box-fraction units after normalisation).
    pub xadd: F1<f32>,
    pub yadd: F1<f32>,
    pub zadd: F1<f32>,
    pub phix: F1<f32>,
    pub phiy: F1<f32>,
    pub phiz: F1<f32>,
    pub dmax: f32,
    pub volume: f32,
    pub amasst: F1<f32>,
    pub moldiam: F1<f32>,
    pub amasstsqrt: F1<f32>,
    pub interntrans: F1<f32>,
    pub internrot: F1<f32>,
    pub wholerot: F1<f32>,
    pub wholetrans: F1<f32>,
    // ---- molecom_mod
    pub ngrrot: i32,
    pub rotgrp: F1<String>,
    pub nrot: F1<i32>,
    pub labaxi: F2<i32>,
    pub rotf: F1<i32>,
    pub rotl: F1<i32>,
    pub mrot: F1<i32>,
    pub nrottot: i32,
    pub labrot: F1<i32>,
    pub coarsegraining: bool,
    pub qscaleinputfactor: f32,
    pub selfbroadpwr: f32,
    pub scalefacmol: F1<f32>,
    pub rfractmol: F1<f32>,
    pub natomsphere: F1<f32>,
    pub rhosphere: F1<f32>,
    pub radminsphere: F1<f32>,
    pub radmaxsphere: F1<f32>,
    pub pseudofacmol: F1<f32>,
    pub natomcg: f32,
    pub natomrcg: f32,
    pub rhocg: f32,
}

impl Mol {
    pub fn new() -> Mol {
        Mol {
            amassa: F1::new(MCOMP),
            amassr: F1::new(MPAIRS),
            charge: F1::new(MCOMP),
            qradius: F1::new(MCOMP),
            qcor: F1::new(MPAIRS),
            dcor: F1::new(MPAIRS),
            sigcor: F1::new(MCOMP),
            epscor: F1::new(MCOMP),
            dscor: F1::new(MPAIRS),
            ecor: F1::new(MPAIRS),
            cfrac: F1::new(MCOMP),
            cfraccg: F1::new(MCOMP),
            type_: F1::new(MCOMP),
            atsymbol: F1::new(MCOMP),
            typefile: F1::new(MCOMP),
            qtypea: F1::new(MCOMP),
            npairref: F2::new(MCOMP, MCOMP),
            n1extr: F1::new(MPAIRS),
            n2extr: F1::new(MPAIRS),
            iexchange: F1::new(MCOMP),
            iatomicnotype: F1::new(MCOMP),
            moltype: F1::new(MCOMP),
            moltemplate: F1::new(MCOMP),
            nato: F1::new(MMOL),
            ngrot: F1::new(MMOL),
            naf: F1::new(MMOL),
            nal: F1::new(MMOL),
            nrf: F1::new(MMOL),
            nrl: F1::new(MMOL),
            mtype: F1::new(MMOL),
            nmovei: F2::new(2, MCOMP),
            ninter: F2::new(2, MCOMP),
            nwhole: F2::new(2, MCOMP),
            nmovem: F2::new(2, MCOMP),
            nstepreject: F1::new(MMOL),
            tetheredatom: F1::new(MMOL),
            moltetheredatom: F1::new(MCOMP),
            tethered: F1::new(MMOL),
            moltethered: F1::new(MCOMP),
            xadd: F1::new(MMOL),
            yadd: F1::new(MMOL),
            zadd: F1::new(MMOL),
            phix: F1::new(MMOL),
            phiy: F1::new(MMOL),
            phiz: F1::new(MMOL),
            amasst: F1::new(MMOL),
            moldiam: F1::new(MMOL),
            amasstsqrt: F1::new(MMOL),
            interntrans: F1::new(MCOMP),
            internrot: F1::new(MCOMP),
            wholerot: F1::new(MCOMP),
            wholetrans: F1::new(MCOMP),
            tbond: F2::new(0, 0),
            iabond: F2::new(0, 0),
            ..Default::default()
        }
    }

    /// Make sure per-atom arrays can hold `n` atoms (the Fortran uses a fixed
    /// maximum `matom`; here storage grows as atoms are read).
    pub fn ensure_atoms(&mut self, n: usize) {
        if self.xa.len() >= n {
            return;
        }
        let m = n.max(2 * self.xa.len()).max(1024);
        self.xa.grow(m);
        self.ya.grow(m);
        self.za.grow(m);
        self.alabel.grow(m);
        self.nbond.grow(m);
        self.itype.grow(m);
        self.mref.grow(m);
        self.qtype.grow(m);
        self.imageatom.grow(m);
        let nb = self.tbond.n2.max(1);
        self.tbond.grow(m, nb);
        self.iabond.grow(m, nb);
    }
    /// Make sure the bond tables can hold `n` bonds per atom.
    pub fn ensure_bonds(&mut self, n: usize) {
        if n > self.tbond.n2 {
            let m1 = self.tbond.n1;
            self.tbond.grow(m1, n);
            self.iabond.grow(m1, n);
        }
    }
    /// `testsize_rotgrps`: grow the rotational-group arrays.
    pub fn ensure_rotgrps(&mut self, n: usize) {
        if self.rotgrp.len() < n {
            let m = n.max(2 * self.rotgrp.len()).max(64);
            self.rotgrp.grow(m);
            self.labaxi.grow(m, 2);
            self.nrot.grow(m);
            self.rotf.grow(m);
            self.rotl.grow(m);
            self.mrot.grow(m);
            self.ndihedralang.grow(m);
            self.dihedralang.grow(NEAR, m);
            if self.labrotdih.n3 < m || self.labrotdih.n1 < 2 || self.labrotdih.n2 < NEAR {
                let mut nl: F3<i32> = F3::new(2, NEAR, m);
                let l = &self.labrotdih;
                for k in 1..=l.n3 {
                    for j in 1..=l.n2.min(NEAR) {
                        for i in 1..=l.n1.min(2) {
                            nl[(i, j, k)] = l[(i, j, k)];
                        }
                    }
                }
                self.labrotdih = nl;
            }
        }
    }
    /// `testsize_labrot`.
    pub fn ensure_labrot(&mut self, n: usize) {
        if self.labrot.len() < n {
            let m = n.max(2 * self.labrot.len()).max(64);
            self.labrot.grow(m);
        }
    }
}

/// `epsr16.inc`: r and Q grids, accumulated structure factors and radial
/// distribution functions.
#[derive(Clone, Debug, Default)]
pub struct Dist {
    /// r-values of the distribution histograms (Å).
    pub x: F1<f32>,
    /// `x / dmax` (box-fraction units).
    pub xstar: F1<f32>,
    /// Number of atoms in each r-shell for a uniform density (normalisation of g(r)).
    pub ycoh: F1<f32>,
    pub entwt: F1<f32>,
    /// Fourier kernel `foqols(r, Q)` used to turn histograms into S(Q).
    pub foqols: F2<f32>,
    /// Running average of the intermolecular partial structure factors.
    pub soqinter: F2<f32>,
    pub sq2inter: F2<f32>,
    /// Running average of the intramolecular partial structure factors.
    pub soqintra: F2<f32>,
    pub sq2intra: F2<f32>,
    pub cellst: f32,
    pub sigma2: F1<f32>,
    pub rhofac: F1<f32>,
    pub xwt: F1<f32>,
    pub cell: f32,
    pub rmaxcell: f32,
    pub cellmax: f32,
    pub upol: f32,
    pub rho: f32,
    pub amubar: f32,
    pub amubar2: f32,
    pub rejrate: f32,
    pub qminm: f32,
    pub qminp: f32,
    pub rlog: f32,
    /// Running averages of inter/intra-molecular partial g(r).
    pub gorinter: F2<f32>,
    pub gr2inter: F2<f32>,
    pub gorintra: F2<f32>,
    pub gr2intra: F2<f32>,
    pub energyref: F1<f32>,
    pub energyemp: F1<f32>,
    pub energytot: F1<f32>,
    pub npt: i32,
    pub npt1: i32,
    pub npt2: i32,
    pub nfile: i32,
    pub ndist: i32,
    pub nqmin: F1<i32>,
    pub nqmax: F1<i32>,
    pub nqpts: F1<i32>,
    pub nptmin: i32,
    pub nptmax: i32,
    pub entsig: F2<f32>,
}

impl Dist {
    pub fn new() -> Dist {
        Dist {
            x: F1::new(MPZ),
            xstar: F1::new(MPZ),
            ycoh: F1::new(MPZ),
            entwt: F1::new(MPZ),
            foqols: F2::new(MPZ, MPTS),
            sigma2: F1::new(MFILE),
            rhofac: F1::new(MFILE),
            xwt: F1::new(MFILE),
            nqmin: F1::new(MFILE),
            nqmax: F1::new(MFILE),
            nqpts: F1::new(MPTS),
            energyref: F1::new(MPAIRS),
            energyemp: F1::new(MPAIRS),
            energytot: F1::new(MPAIRS),
            ..Default::default()
        }
    }
    /// Allocate the pair-dependent arrays once the number of pairs is known.
    pub fn alloc_pairs(&mut self, npairs: usize) {
        let np = npairs.max(1);
        self.soqinter = F2::new(MPTS, np);
        self.sq2inter = F2::new(MPTS, np);
        self.soqintra = F2::new(MPTS, np);
        self.sq2intra = F2::new(MPTS, np);
        self.gorinter = F2::new(MPZ, np);
        self.gr2inter = F2::new(MPZ, np);
        self.gorintra = F2::new(MPZ, np);
        self.gr2intra = F2::new(MPZ, np);
        self.entsig = F2::new(MPZ, np);
    }
}

/// `potential_par.inc`: reference and empirical potentials and their coefficients.
#[derive(Clone, Debug, Default)]
pub struct Pot {
    pub tlab1: F1<String>,
    pub tlab2: F1<String>,
    pub addpotential: String,
    pub ncoeffp: i32,
    pub npcoeff: F1<i32>,
    pub npits: i32,
    pub npitss: i32,
    pub npairspcof: i32,
    pub irefpcof: F1<i32>,
    pub npgauss: F1<i32>,
    pub pstepd: F1<f32>,
    pub pwt: F1<f32>,
    pub pwtold: F1<f32>,
    /// r-space basis functions of the empirical potential (`expgr((j-1)*nterms+m)`).
    pub expgr: F1<f32>,
    /// Q-space basis functions (`expsq((i-1)*nterms+m)`).
    pub expsq: F1<f32>,
    pub pdmin: f32,
    pub pdstep: f32,
    pub pdmax: f32,
    pub paccept: f32,
    pub psigma2: f32,
    pub psigma2read: f32,
    /// EP coefficients fitted to each dataset.
    pub pcoefd: F2<f32>,
    /// EP coefficients for each atom pair (the refined empirical potential).
    pub pcoefp: F2<f32>,
    pub rminex: F1<f32>,
    pub coordex: F1<f32>,
    pub rmaxex: F1<f32>,
    pub pcoefex: F1<f32>,
    pub pcoefexmax: F1<f32>,
    pub potex: F1<f32>,
    pub potexmax: F1<f32>,
    pub weight: F3<f32>,
    pub uwt: F1<f32>,
    /// Intramolecular potential table and its gradient (`xav`, `xav1`).
    pub xav: F2<f32>,
    pub xav1: F2<f32>,
    pub wsum: F1<f32>,
    /// Empirical potential table and gradient (`xavs`, `xavs1`).
    pub xavs: F2<f32>,
    pub xavs1: F2<f32>,
    pub wfac: F1<f32>,
    pub wfacc: F1<f32>,
    /// Total potential (reference + empirical).
    pub xavtot: F2<f32>,
    pub xavtot1: F2<f32>,
    /// Reference potential (Lennard-Jones + Coulomb) with and without additions.
    pub epot: F2<f32>,
    pub epot1: F2<f32>,
    pub ecoul: F2<f32>,
    pub ecoul1: F2<f32>,
    pub feedback: f32,
    pub feedbackold: f32,
    pub efacm: F1<f32>,
    pub efacold: F1<f32>,
    pub efacp: F1<f32>,
    pub rpgauss: F2<f32>,
    pub wpgauss: F2<f32>,
    pub hpgauss: F2<f32>,
    pub absolute_energy: f32,
    pub potfacamp: f32,
    pub potfac: f32,
    pub refpotfac: f32,
    pub ereqstep: f32,
    pub erequnit: f32,
    pub coulfac: F1<f32>,
    pub gaussian: bool,
    pub harmonic: bool,
}

impl Pot {
    pub fn new() -> Pot {
        Pot {
            tlab1: F1::new(MPAIRS),
            tlab2: F1::new(MPAIRS),
            npcoeff: F1::new(MFILE),
            irefpcof: F1::new(MPAIRS),
            npgauss: F1::new(MPAIRS),
            pstepd: F1::new(MFILE),
            pwt: F1::new(MCOEFF),
            pwtold: F1::new(MCOEFF),
            pcoefd: F2::new(MCOEFF, MFILE),
            pcoefp: F2::new(MCOEFF, MPAIRS),
            rminex: F1::new(MPAIRS),
            coordex: F1::new(MPAIRS),
            rmaxex: F1::new(MPAIRS),
            pcoefex: F1::new(MPAIRS),
            pcoefexmax: F1::new(MPAIRS),
            potex: F1::new(MPAIRS),
            potexmax: F1::new(MPAIRS),
            uwt: F1::new(MPAIRS),
            wsum: F1::new(MFILE),
            wfac: F1::new(MPZ),
            wfacc: F1::new(MPZ),
            efacm: F1::new(MFILE),
            efacold: F1::new(MFILE),
            efacp: F1::new(MPAIRS),
            rpgauss: F2::new(MPGAUSS, MPAIRS),
            wpgauss: F2::new(MPGAUSS, MPAIRS),
            hpgauss: F2::new(MPGAUSS, MPAIRS),
            coulfac: F1::new(MPAIRS),
            ..Default::default()
        }
    }
    /// Allocate arrays whose size depends on the number of pairs, data files
    /// and coefficients.
    pub fn alloc(&mut self, npairs: usize, nfile: usize) {
        let np = npairs.max(1);
        self.weight = F3::new(2, MFILE.max(nfile), np);
        for a in [
            &mut self.xav,
            &mut self.xav1,
            &mut self.xavs,
            &mut self.xavs1,
            &mut self.xavtot,
            &mut self.xavtot1,
            &mut self.epot,
            &mut self.epot1,
            &mut self.ecoul,
            &mut self.ecoul1,
        ] {
            *a = F2::new(MPZ, np);
        }
    }
    /// Basis-function tables (`expgr`, `expsq`) sized for `npt` points × `nterms`.
    pub fn alloc_basis(&mut self, npts: usize, nterms: usize) {
        let n = npts * nterms.max(1);
        if self.expgr.len() < n {
            self.expgr = F1::new(n);
        }
        if self.expsq.len() < n {
            self.expsq = F1::new(n);
        }
    }
}

/// `weights.inc`: neutron weights and the inverted weights matrix.
#[derive(Clone, Debug, Default)]
pub struct Weights {
    /// Weight of each partial in each dataset (intermolecular).
    pub cwt: F2<f64>,
    /// Weight of each partial in each dataset (intramolecular).
    pub ciwt: F2<f64>,
    /// `cwtpot(ipair, idataset)` — how each dataset's EP contributes to each pair.
    pub cwtpot: F2<f64>,
    /// Inverted weights matrix `cwttq(iq, ipair, idataset or nfile+jpair)`.
    pub cwttq: F3<f64>,
    pub nxfile: i32,
    pub xfile: F1<i32>,
    pub nqinv: i32,
    pub qinv: F1<i32>,
}

impl Weights {
    pub fn new() -> Weights {
        Weights { xfile: F1::new(MFILE), qinv: F1::new(MQINV), ..Default::default() }
    }
    pub fn alloc(&mut self, npairs: usize, nfile: usize) {
        let np = npairs.max(1);
        let nf = nfile.max(1);
        self.cwt = F2::new(nf, np);
        self.ciwt = F2::new(nf, np);
        self.cwtpot = F2::new(np, np.max(nf));
        self.cwttq = F3::new(MQINV, np, np + nf);
    }
}

/// Maximum number of X-ray datasets (`mxfile` in `X-rayf0.inc`).
pub const MXFILE: usize = 2;

/// `X-rayf0.inc`: X-ray form factors and Q-dependent weights.
#[derive(Clone, Debug, Default)]
pub struct XRay {
    /// Atomic form factors `f0x(Q, type, xfile)` and modified ones `fm0x`.
    pub f0x: F3<f64>,
    pub fm0x: F3<f64>,
    /// Q-dependent weights `weightx(Q, pair, xfile)`.
    pub weightx: F3<f64>,
    /// Self-scattering background from mismatched form factors `backx(Q, xfile)`.
    pub backx: F2<f64>,
    /// Dataset number of each X-ray file, normalisation choices.
    pub datapoint: F1<i32>,
    pub datanorm: F1<i32>,
    pub normxx: F1<i32>,
    pub iafformaffx: F1<i32>,
}

impl XRay {
    pub fn new() -> XRay {
        XRay {
            f0x: F3::new(MPTS, MCOMP, MXFILE),
            fm0x: F3::new(MPTS, MCOMP, MXFILE),
            weightx: F3::new(1, 1, 1),
            backx: F2::new(MPTS, MXFILE),
            datapoint: F1::new(MXFILE),
            datanorm: F1::new(MXFILE),
            normxx: F1::new(MXFILE),
            iafformaffx: F1::new(MXFILE),
        }
    }
}

/// `epsrinput_par.inc`: run-control parameters and running sums.
#[derive(Clone, Debug, Default)]
pub struct InputPar {
    pub nisum: i32,
    pub nesum: i32,
    pub nesumt: i32,
    pub nsum: i32,
    pub nsums: i32,
    pub roverlap: f32,
    pub rminfac: f32,
    pub rminpt: f32,
    pub rmaxpt: f32,
    pub rbroad: f32,
    pub expecf: f32,
    pub expecfpwr: f32,
    pub rcharge: f32,
    pub ensum: f32,
    pub ensum2: f32,
    pub prsum: f32,
    pub prsum2: f32,
    pub lorchqwidth: f32,
    pub lorchqwidthq: f32,
    pub rminbragg: f32,
    pub rmaxgr: f32,
}

/// `symbols.inc` and friends: the keyword/value table that represents the
/// `.EPSR.inp` and `.pcof` files, plus file names.
#[derive(Clone, Debug, Default)]
pub struct Symbols {
    pub symbol: F1<String>,
    pub symbol_value: F1<String>,
    pub symbol_text: F1<String>,
    pub symbol_default: F1<String>,
    pub symbol_type: F1<char>,
    pub iflag: F1<i32>,
    pub symbol_save: String,
    pub fname_inp: String,
    pub fnameinp: String,
    pub fname_run: String,
    pub fnamerun: String,
    pub fname_ato: String,
    pub fnameato: String,
    pub fname_pcof: String,
    pub fnamepcof: String,
    pub dataprompt: String,
    pub prompt: String,
    pub extinp: String,
    pub symbol_setup: String,
    pub nsymboluptodata: i32,
    pub nsymboldata: i32,
    pub nsymbol: i32,
    pub nsymbolplus: i32,
    pub nfiles: i32,
    pub isave: i32,
    pub inotrecog: i32,
    pub irefrho: i32,
    // iresetnumbers.inc
    pub iresetref: i32,
    pub iinitref: i32,
    pub iqminref: i32,
    pub iroverlap: i32,
    pub irminfac: i32,
    pub incoeffp: i32,
    pub iereqfirst: i32,
    pub iereqlast: i32,
    // directory_name.inc
    /// Folder containing the input files, *with* a trailing separator.
    pub directory_name: String,
}

impl Symbols {
    pub fn new() -> Symbols {
        Symbols {
            symbol: F1::new(MSYMBOL),
            symbol_value: F1::new(MSYMBOL),
            symbol_text: F1::new(MSYMBOL),
            symbol_default: F1::new(MSYMBOL),
            symbol_type: F1 { data: vec![' '; MSYMBOL] },
            iflag: F1::new(MSYMBOL),
            ..Default::default()
        }
    }
    /// `function ireff(icount, ifilecount)`: storage index of symbol `icount`
    /// for dataset `ifilecount` (per-dataset symbols are replicated).
    pub fn ireff(&self, icount: i32, ifilecount: i32) -> i32 {
        if icount > self.nsymboldata {
            icount + (self.nfiles - 1) * self.nsymbolplus
        } else if icount > self.nsymboluptodata {
            icount + (ifilecount - 1) * self.nsymbolplus
        } else {
            icount
        }
    }
    /// `search_for_symbol(text, ifound)`.
    pub fn search_for_symbol(&self, text: &str, ifound: i32) -> i32 {
        let text = crate::fortran::fstr(text, 20);
        if ifound > 0 && ifound <= self.nsymbol && self.symbol[ifound] == text {
            return ifound;
        }
        let mut icount = 0;
        while icount < self.nsymbol {
            icount += 1;
            if self.symbol[icount] == text {
                return icount;
            }
        }
        0
    }
    /// Full path of a file in the input folder (`get_full_filename`).
    pub fn full(&self, short: &str) -> PathBuf {
        if short.contains('/') {
            PathBuf::from(short)
        } else {
            PathBuf::from(format!("{}{}", self.directory_name, short))
        }
    }
}

/// Module `cellneighbourlist`: partition of the box into cells for fast
/// neighbour searching.
#[derive(Clone, Debug, Default)]
pub struct CellList {
    pub ncelldiva: i32,
    pub ncelldivb: i32,
    pub ncelldivc: i32,
    pub ncelldivsq: i32,
    pub ncell: i32,
    pub ncellneighb: i32,
    pub maxcellno: i32,
    pub nmolcell: i32,
    pub nmolsincell: F1<i32>,
    pub molcellno: F1<i32>,
    pub iacellneighb: F1<i32>,
    pub ibcellneighb: F1<i32>,
    pub iccellneighb: F1<i32>,
    pub molsincell: F2<i32>,
}

/// Module `ergpairreference`: stored pair energies between neighbouring
/// molecules, so that a Monte Carlo move only needs the energy change.
#[derive(Clone, Debug, Default)]
pub struct ErgPairs {
    pub ergpair: Vec<f32>,
    pub nergpair: i32,
    pub nergpairmax: i32,
    pub molergpair: F2<i32>,
    pub ergpairref: F2<i32>,
    pub nergpairm: F1<i32>,
    pub nnbrefmax: i32,
    pub ergold: F1<f32>,
    pub ergnew: F1<f32>,
    pub nbref: F1<i32>,
    pub nnbref: i32,
    pub nmolhere: i32,
}

/// Module `timer`: accumulated wall-clock time (milliseconds), written to
/// `fort.35` only.  Not part of the physics.
#[derive(Clone, Debug, Default)]
pub struct Timer {
    pub sumtime1: i64,
    pub sumtime2: i64,
}

/// Everything the EPSR engine knows; one instance per simulation.
#[derive(Clone, Debug)]
pub struct State {
    pub lat: Lattice,
    pub mol: Mol,
    pub dist: Dist,
    pub pot: Pot,
    pub wts: Weights,
    pub inp: InputPar,
    pub sym: Symbols,
    pub cells: CellList,
    pub erg: ErgPairs,
    pub timer: Timer,
    pub rng: Ran1,
    /// X-ray form factors and weights (`X-rayf0.inc`).
    pub xray: XRay,
    /// Tethering energy per molecule from the last `fmolec1` (`ftethered`).
    pub ftethered: f32,
    /// Container forces (`containerforces_mod`).
    pub cont: crate::containerforces::ContainerForces,
    /// Echo the progress messages the Fortran writes to unit 6 (stdout).
    pub verbose: bool,
    /// Coarse-graining arrays (`molecom_mod`, used when `qinscalefc` > 1).
    pub cg: crate::coarsegrain::CoarseGrain,
    /// Three-body potentials (`threebody_mod`).
    pub tb: crate::threebody::ThreeBody,
    /// Coordination constraints (`coordconstraints_mod`).
    pub cc: crate::coordconstraint::CoordConstraints,
}

impl State {
    /// Progress message (the Fortran `write(6,...)` lines).
    pub fn log(&self, s: &str) {
        if self.verbose {
            println!("{s}");
        }
    }
    /// A logging closure that does not borrow the state.
    pub fn logger(&self) -> impl Fn(&str) {
        let verbose = self.verbose;
        move |s: &str| {
            if verbose {
                println!("{s}");
            }
        }
    }
    pub fn new() -> State {
        State {
            verbose: false,
            lat: Lattice::default(),
            mol: Mol::new(),
            dist: Dist::new(),
            pot: Pot::new(),
            wts: Weights::new(),
            inp: InputPar::default(),
            sym: Symbols::new(),
            cells: CellList::default(),
            erg: ErgPairs::default(),
            timer: Timer::default(),
            rng: Ran1::default(),
            xray: XRay::new(),
            ftethered: 0.0,
            cont: Default::default(),
            cg: Default::default(),
            tb: Default::default(),
            cc: Default::default(),
        }
    }
}

impl Default for State {
    fn default() -> Self {
        State::new()
    }
}
