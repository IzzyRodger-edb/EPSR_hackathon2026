//! Convenience accessors for embedding the engine (e.g. in the Python wrapper).
//!
//! The simulation state uses Fortran-style 1-based, column-major arrays with
//! fixed leading dimensions (see [`crate::state`]).  The methods here expose
//! the commonly needed results as plain 0-based slices / vectors of the
//! *used* length, in physical units (Å, Å⁻¹, kJ/mol), so callers never need
//! to know about the internal layout.
//!
//! All slices borrow from the simulation and reflect its current state: call
//! them after [`Epsr::iterate`] to see that iteration's results.

use crate::epsr::Epsr;

impl Epsr {
    /// Number of datasets (`ndata`).
    pub fn n_datasets(&self) -> usize {
        self.st.dist.nfile as usize
    }
    /// Number of distinct atom-type pairs.
    pub fn n_pairs(&self) -> usize {
        self.st.mol.npairs as usize
    }
    /// Atom-type labels, in type order.
    pub fn atom_types(&self) -> Vec<String> {
        (1..=self.st.mol.ntype).map(|it| self.st.mol.type_[it].clone()).collect()
    }
    /// Labels `(type1, type2)` of each pair, in the order used by all
    /// pair-resolved arrays.
    pub fn pair_labels(&self) -> Vec<(String, String)> {
        let m = &self.st.mol;
        (1..=m.npairs).map(|ip| (m.type_[m.n1extr[ip]].clone(), m.type_[m.n2extr[ip]].clone())).collect()
    }
    /// Number of Q points (`nq`).
    pub fn nq(&self) -> usize {
        self.ctl.nq as usize
    }
    /// Number of r points of the distribution functions (`nz`).
    pub fn nr(&self) -> usize {
        self.nz as usize
    }
    /// Q grid (Å⁻¹), length [`Epsr::nq`].
    pub fn q_grid(&self) -> &[f32] {
        &self.q.as_slice()[..self.nq()]
    }
    /// r grid (Å), length [`Epsr::nr`].
    pub fn r_grid(&self) -> &[f32] {
        &self.z.as_slice()[..self.nr()]
    }
    /// Running-average intermolecular partial g(r) of pair `ip` (0-based).
    pub fn gr_inter(&self, ip: usize) -> &[f32] {
        &self.st.dist.gorinter.col(ip + 1)[..self.nr()]
    }
    /// Running-average intramolecular partial g(r) of pair `ip`.
    pub fn gr_intra(&self, ip: usize) -> &[f32] {
        &self.st.dist.gorintra.col(ip + 1)[..self.nr()]
    }
    /// Running-average intermolecular partial S(Q) of pair `ip`.
    pub fn sq_inter(&self, ip: usize) -> &[f32] {
        &self.st.dist.soqinter.col(ip + 1)[..self.nq()]
    }
    /// Running-average intramolecular partial S(Q) of pair `ip`.
    pub fn sq_intra(&self, ip: usize) -> &[f32] {
        &self.st.dist.soqintra.col(ip + 1)[..self.nq()]
    }
    /// Measured data of dataset `id` (0-based) interpolated onto the Q grid
    /// (after `tweak` and density normalisation); zero outside the data range.
    pub fn data_fq(&self, id: usize) -> &[f32] {
        &self.datai.col(id + 1)[..self.nq()]
    }
    /// Simulated total structure factor of dataset `id`, i.e. Σ (intra + inter)
    /// weighted partials, computed from the current running averages.
    pub fn simulated_fq(&self, id: usize) -> Vec<f32> {
        let inn = id as i32 + 1;
        let w = &self.st.wts;
        let d = &self.st.dist;
        (1..=self.ctl.nq)
            .map(|i| {
                let mut sum: f32 = 0.0;
                for ic in 1..=self.st.mol.npairs {
                    sum = ((sum as f64 + w.ciwt[(inn, ic)] * d.soqintra[(i, ic)] as f64) + w.cwt[(inn, ic)] * d.soqinter[(i, ic)] as f64) as f32;
                }
                sum
            })
            .collect()
    }
    /// First and last Q-grid indices (0-based, inclusive) used when fitting dataset `id`.
    pub fn data_q_range(&self, id: usize) -> (usize, usize) {
        let inn = id + 1;
        ((self.st.dist.nqmin[inn] - 1) as usize, (self.st.dist.nqmax[inn] - 1) as usize)
    }
    /// Empirical potential (kJ/mol) of pair `ip` on the r grid.
    pub fn empirical_potential(&self, ip: usize) -> &[f32] {
        &self.st.pot.xavs.col(ip + 1)[..self.nr()]
    }
    /// Reference potential (Lennard-Jones + Coulomb, kJ/mol) of pair `ip`.
    pub fn reference_potential(&self, ip: usize) -> &[f32] {
        &self.st.pot.ecoul.col(ip + 1)[..self.nr()]
    }
    /// EP coefficients of pair `ip` (`pcoefp(2..ncoeffp)`; `pcoefp(1)` is the refinement flag).
    pub fn ep_coefficients(&self, ip: usize) -> &[f32] {
        let irefp = self.st.pot.irefpcof[ip + 1] as usize;
        let n = self.st.pot.ncoeffp as usize;
        &self.st.pot.pcoefp.col(irefp)[..n]
    }
    /// Number of molecules and atoms.
    pub fn n_molecules(&self) -> usize {
        self.st.mol.nmol as usize
    }
    pub fn n_atoms(&self) -> usize {
        self.st.mol.natom as usize
    }
    /// Atom labels (one per atom).
    pub fn atom_labels(&self) -> Vec<String> {
        (1..=self.st.mol.natom).map(|ia| self.st.mol.alabel[ia].clone()).collect()
    }
    /// For each molecule, the 0-based range of its atoms.
    pub fn molecule_atoms(&self) -> Vec<std::ops::Range<usize>> {
        let m = &self.st.mol;
        (1..=m.nmol).map(|im| (m.naf[im] - 1) as usize..m.nal[im] as usize).collect()
    }
    /// Cubic box side (Å).
    pub fn box_size(&self) -> f32 {
        self.st.mol.dmax
    }
    /// Absolute atom positions (Å), molecule centre + atom offset.  The engine
    /// works in box-fraction units between iterations; this converts back.
    pub fn atom_positions(&self) -> Vec<[f32; 3]> {
        let m = &self.st.mol;
        let mut out = Vec::with_capacity(m.natom as usize);
        for im in 1..=m.nmol {
            for ia in m.naf[im]..=m.nal[im] {
                out.push([
                    m.xadd[im] * m.dmax + m.xa[ia] * m.dmaxmol,
                    m.yadd[im] * m.dmax + m.ya[ia] * m.dmaxmol,
                    m.zadd[im] * m.dmax + m.za[ia] * m.dmaxmol,
                ]);
            }
        }
        out
    }
    /// Current value of an `.EPSR.inp` keyword (global or, with `dataset`,
    /// per-dataset), e.g. `get_parameter("ntimes", None)`.
    pub fn get_parameter(&self, name: &str, dataset: Option<usize>) -> Option<String> {
        let s = &self.st.sym;
        let k = s.search_for_symbol(name, 0);
        if k == 0 {
            return None;
        }
        let iref = s.ireff(k, dataset.map(|d| d as i32 + 1).unwrap_or(1));
        Some(s.symbol_value[iref].trim().to_string())
    }
}
