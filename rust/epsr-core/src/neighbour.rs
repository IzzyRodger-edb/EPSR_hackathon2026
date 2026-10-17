//! Periodic boundaries and neighbour lists.
//!
//! Translation of:
//! * `minimagec`, `minimage`        (`minimagec.f`)
//! * module `cellneighbourlist`, `atomcellno`, `cellcoord` (`cellneighbourlist.f`)
//! * module `ergpairreference`      (`ergpairreference.f`)
//!
//! The box is divided into cells (`initialise_cellneighbourlist`); the
//! molecules whose centres lie in a cell and its neighbouring cells form the
//! candidate neighbour list of a molecule.  `fmolec1` then evaluates and
//! stores every non-zero molecule-pair energy in `ErgPairs`, so that a Monte
//! Carlo move of molecule `im` only needs the pair energies with its stored
//! neighbours (`make_neighbour_list_from_mols`).  The neighbour lists are
//! rebuilt once per EPSR iteration.

use crate::fortran::fmath::nint;
use crate::fortran::{F1, F2};
use crate::state::{CellList, ErgPairs, Lattice, State};

/// `subroutine minimage(fa,fb,fc)`: minimum-image convention in box-fraction units.
#[inline(always)]
pub fn minimage(fa: &mut f32, fb: &mut f32, fc: &mut f32) {
    *fc -= nint(*fc) as f32;
    *fb -= nint(*fb) as f32;
    *fa -= nint(*fa) as f32;
}

/// The box geometry used by [`minimagec`].  In the Fortran `minimagec` reads
/// the lattice from COMMON `/lattice_r/`; the port keeps the active lattice
/// per thread, refreshed by [`State::sync_lattice`] whenever the box changes
/// (reading/writing `.ato` files, (un)normalising coordinates).
#[derive(Clone, Copy)]
struct BoxGeom {
    cubic: bool,
    alatx: f32,
    blatx: f32,
    blaty: f32,
    clatx: f32,
    claty: f32,
    clatz: f32,
}

thread_local! {
    static GEOM: std::cell::Cell<BoxGeom> = const { std::cell::Cell::new(BoxGeom {
        cubic: true, alatx: 1.0, blatx: 0.0, blaty: 1.0, clatx: 0.0, claty: 0.0, clatz: 1.0,
    }) };
}

/// `subroutine minimagec(x,y,z)`: minimum image for a cubic or general box
/// (the active lattice of this thread, see [`State::sync_lattice`]).
#[inline(always)]
pub fn minimagec(x: &mut f32, y: &mut f32, z: &mut f32) {
    let g = GEOM.with(|c| c.get());
    if g.cubic {
        minimage(x, y, z);
    } else {
        let mut fc = *z / g.clatz;
        let mut fb = (*y - fc * g.claty) / g.blaty;
        let mut fa = (*x - fb * g.blatx - fc * g.clatx) / g.alatx;
        minimage(&mut fa, &mut fb, &mut fc);
        *z = fc * g.clatz;
        *y = fb * g.blaty + fc * g.claty;
        *x = fa * g.alatx + fb * g.blatx + fc * g.clatx;
    }
}

impl State {
    /// Make this state's lattice the one used by [`minimagec`] on this thread.
    pub fn sync_lattice(&self) {
        let l = &self.lat;
        let g = BoxGeom { cubic: l.cubic, alatx: l.alatx, blatx: l.blatx, blaty: l.blaty, clatx: l.clatx, claty: l.claty, clatz: l.clatz };
        GEOM.with(|c| c.set(g));
    }
}

/// `integer function cellcoord(x, ncelldiv)`.
fn cellcoord(x: f32, ncelldiv: i32) -> i32 {
    let x1 = x - nint(x) as f32;
    let mut ix = ((x1 + 0.5) * ncelldiv as f32) as i32;
    while ix < 0 {
        ix += ncelldiv;
    }
    while ix >= ncelldiv {
        ix -= ncelldiv;
    }
    ix
}

/// `integer function atomcellno(x,y,z)`.
fn atomcellno(c: &CellList, l: &Lattice, x: f32, y: f32, z: f32) -> i32 {
    let fc = z / l.clatz;
    let fb = (y - fc * l.claty) / l.blaty;
    let fa = (x - fb * l.blatx - fc * l.clatx) / l.alatx;
    let ia = cellcoord(fa, c.ncelldiva);
    let ib = cellcoord(fb, c.ncelldivb);
    let ic = cellcoord(fc, c.ncelldivc);
    1 + ia + c.ncelldiva * (ib + c.ncelldivb * ic)
}

impl CellList {
    /// `subroutine findiaibic(cellno, ia, ib, ic)`.
    fn findiaibic(&self, cellno: i32) -> (i32, i32, i32) {
        let cellno1 = cellno - 1;
        let ic = cellno1 / self.ncelldivsq;
        let ib = (cellno1 - ic * self.ncelldivsq) / self.ncelldiva;
        let ia = cellno1 - self.ncelldiva * (ib + self.ncelldivb * ic);
        (ia, ib, ic)
    }
    /// `subroutine getcellreference(ia, ib, ic, icellref)`.
    fn getcellreference(&self, ia: i32, ib: i32, ic: i32) -> i32 {
        let ia = periodic(ia, self.ncelldiva);
        let ib = periodic(ib, self.ncelldivb);
        let ic = periodic(ic, self.ncelldivc);
        1 + ia + self.ncelldiva * (ib + self.ncelldivb * ic)
    }
    /// `subroutine make_neighbour_list_from_cells(ic, nnbref, nbref)`: all
    /// molecules in cell `ic` and its neighbouring cells.
    pub fn make_neighbour_list_from_cells(&self, ic: i32, nbref: &mut F1<i32>) -> i32 {
        let mut nnbref = 0;
        let (iao, ibo, ico) = self.findiaibic(ic);
        for icref in 1..=self.ncellneighb {
            let icellref =
                self.getcellreference(iao + self.iacellneighb[icref], ibo + self.ibcellneighb[icref], ico + self.iccellneighb[icref]);
            for iref in 1..=self.nmolsincell[icellref] {
                nnbref += 1;
                nbref[nnbref] = self.molsincell[(iref, icellref)];
            }
        }
        nnbref
    }
    /// `subroutine get_nnbrefmax(nnbrefmax)`.
    pub fn get_nnbrefmax(&self) -> i32 {
        let mut nnbrefmax = 0;
        for ic in 1..=self.ncell {
            let mut nnbref = 0;
            let (iao, ibo, ico) = self.findiaibic(ic);
            for icref in 1..=self.ncellneighb {
                let icellref =
                    self.getcellreference(iao + self.iacellneighb[icref], ibo + self.ibcellneighb[icref], ico + self.iccellneighb[icref]);
                nnbref += self.nmolsincell[icellref];
            }
            if nnbref > nnbrefmax {
                nnbrefmax = nnbref;
            }
        }
        nnbrefmax
    }
}

/// `subroutine findperiodiccell(ix, ncelldiv)`.
fn periodic(mut ix: i32, ncelldiv: i32) -> i32 {
    while ix < 0 {
        ix += ncelldiv;
    }
    while ix >= ncelldiv {
        ix -= ncelldiv;
    }
    ix
}

impl State {
    /// `subroutine initialise_cellneighbourlist(nmolcellin, xadd, yadd, zadd, nmol, rmax)`.
    pub fn initialise_cellneighbourlist(&mut self, nmolcellin: i32, rmax: f32) {
        let l = &self.lat;
        let nmol = self.mol.nmol;
        let c = &mut self.cells;
        c.nmolcell = nmolcellin;
        if c.nmolcell > nmol {
            c.nmolcell = nmol;
        }
        let mut ncelldiv = nint((nmol as f32 / c.nmolcell as f32).powf(1.0f32 / 3.0f32));
        if ncelldiv < 1 {
            ncelldiv = 1;
        }
        let cellsize = 2.0 * l.alatx.min(l.blaty).min(l.clatz) / ncelldiv as f32;
        c.ncelldiva = 2 * (l.alatx / cellsize) as i32 + 1;
        c.ncelldivb = 2 * (l.blaty / cellsize) as i32 + 1;
        c.ncelldivc = 2 * (l.clatz / cellsize) as i32 + 1;
        c.ncell = c.ncelldiva * c.ncelldivb * c.ncelldivc;
        c.ncelldivsq = c.ncelldiva * c.ncelldivb;
        let msg = format!(
            "\ninitialise_cellneighbourlist> Cell numbers along a,b,c and total {:>6} {:>6} {:>6}",
            c.ncelldiva, c.ncelldivb, c.ncelldivc
        );
        let ncell = c.ncell as usize;
        c.nmolsincell = F1::new(ncell);
        c.molcellno = F1::new(nmol as usize);
        c.iacellneighb = F1::new(ncell);
        c.ibcellneighb = F1::new(ncell);
        c.iccellneighb = F1::new(ncell);
        c.maxcellno = 0;
        for im in 1..=nmol {
            let no = atomcellno(c, l, self.mol.xadd[im], self.mol.yadd[im], self.mol.zadd[im]);
            c.molcellno[im] = no;
            c.nmolsincell[no] += 1;
            if c.nmolsincell[no] > c.maxcellno {
                c.maxcellno = c.nmolsincell[no];
            }
        }
        c.molsincell = F2::new(c.maxcellno.max(1) as usize, ncell);
        for ic in 1..=c.ncell {
            c.nmolsincell[ic] = 0;
        }
        for im in 1..=nmol {
            let no = c.molcellno[im];
            c.nmolsincell[no] += 1;
            let k = c.nmolsincell[no];
            c.molsincell[(k, no)] = im;
        }
        self.log(&msg);
        self.makecellneighbourlist(rmax);
    }

    /// `subroutine makecellneighbourlist(rmax)`: which relative cell offsets can
    /// contain molecules within `rmax` of any point in the origin cell.
    fn makecellneighbourlist(&mut self, rmax: f32) {
        let l = &self.lat;
        let c = &mut self.cells;
        let cella = 1.0 / c.ncelldiva as f32;
        let cellb = 1.0 / c.ncelldivb as f32;
        let cellc = 1.0 / c.ncelldivc as f32;
        let cella2 = 0.5 * cella;
        let cellb2 = 0.5 * cellb;
        let cellc2 = 0.5 * cellc;
        let mut xc = [0f32; 8];
        let mut yc = [0f32; 8];
        let mut zc = [0f32; 8];
        let mut ifound = 0usize;
        for ic in [-1i32, 1] {
            let cx = ic as f32 * cellc2 * l.clatc[0];
            let cy = ic as f32 * cellc2 * l.clatc[1];
            let cz = ic as f32 * cellc2 * l.clatc[2];
            for ib in [-1i32, 1] {
                let bx = ib as f32 * cellb2 * l.blatc[0];
                let by = ib as f32 * cellb2 * l.blatc[1];
                let bz = ib as f32 * cellb2 * l.blatc[2];
                for ia in [-1i32, 1] {
                    let ax = ia as f32 * cella2 * l.alatc[0];
                    let ay = ia as f32 * cella2 * l.alatc[1];
                    let az = ia as f32 * cella2 * l.alatc[2];
                    xc[ifound] = ax + bx + cx;
                    yc[ifound] = ay + by + cy;
                    zc[ifound] = az + bz + cz;
                    ifound += 1;
                }
            }
        }
        let ncorner = ifound;
        let nrangea = (c.ncelldiva - 1) / 2;
        let nrangeb = (c.ncelldivb - 1) / 2;
        let nrangec = (c.ncelldivc - 1) / 2;
        let mut ifound: i32 = 0;
        let rmaxsq = rmax * rmax;
        for ic in -nrangec..=nrangec {
            let cx = ic as f32 * cellc * l.clatc[0];
            let cy = ic as f32 * cellc * l.clatc[1];
            let cz = ic as f32 * cellc * l.clatc[2];
            for ib in -nrangeb..=nrangeb {
                let bx = ib as f32 * cellb * l.blatc[0];
                let by = ib as f32 * cellb * l.blatc[1];
                let bz = ib as f32 * cellb * l.blatc[2];
                for ia in -nrangea..=nrangea {
                    let ax = ia as f32 * cella * l.alatc[0];
                    let ay = ia as f32 * cella * l.alatc[1];
                    let az = ia as f32 * cella * l.alatc[2];
                    let xdiffc = ax + bx + cx;
                    let ydiffc = ay + by + cy;
                    let zdiffc = az + bz + cz;
                    let mut found = false;
                    let mut if1 = 0;
                    while if1 < ncorner && !found {
                        let mut if2 = 0;
                        while if2 < ncorner && !found {
                            let xdiff = xdiffc + xc[if1] - xc[if2];
                            let ydiff = ydiffc + yc[if1] - yc[if2];
                            let zdiff = zdiffc + zc[if1] - zc[if2];
                            let rdiffsq = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff;
                            found = rdiffsq < rmaxsq || (ia == 0 && ib == 0 && ic == 0);
                            if2 += 1;
                        }
                        if1 += 1;
                    }
                    if found {
                        ifound += 1;
                        c.iacellneighb[ifound] = ia;
                        c.ibcellneighb[ifound] = ib;
                        c.iccellneighb[ifound] = ic;
                    }
                }
            }
        }
        c.ncellneighb = ifound;
        let msg = format!("\nmakecellneighbourlist> Found {:>6} cell neighbours out of {:>6} total", ifound, c.ncell);
        self.log(&msg);
    }

    /// `subroutine initialise_ergpair(nmol)`.
    pub fn initialise_ergpair(&mut self) {
        let nmol = self.mol.nmol;
        let e = &mut self.erg;
        e.nmolhere = nmol;
        e.nnbrefmax = self.cells.get_nnbrefmax();
        e.nergpair = 0;
        e.nergpairmax = 100 * nmol;
        e.ergpair = vec![0.0; e.nergpairmax as usize];
        let nb = e.nnbrefmax.max(1) as usize;
        e.molergpair = F2::new(nb, nmol as usize);
        e.ergpairref = F2::new(nb, nmol as usize);
        e.nergpairm = F1::new(nmol as usize);
        e.nnbref = 0;
        e.ergold = F1::new(nb);
        e.ergnew = F1::new(nb);
        e.nbref = F1::new(nb);
        let msg = format!(" initialise_ergpair> nnbrefmax = {:>11} {:>11}", e.nnbrefmax, nb);
        self.log(&msg);
    }
}

impl ErgPairs {
    /// `subroutine set_ergpair(im, jm, energy)`: store a pair energy (im == jm
    /// stores the intramolecular energy as the first entry of molecule im).
    pub fn set_ergpair(&mut self, im: i32, jm: i32, energy: f32) {
        if self.nergpairm[im] < self.nnbrefmax && self.nergpairm[jm] < self.nnbrefmax && energy != 0.0 {
            self.nergpair += 1;
            if self.nergpair > self.nergpairmax {
                self.nergpairmax = ((self.nergpairmax + 10) * 12) / 10;
                self.ergpair.resize(self.nergpairmax as usize, 0.0);
            }
            let n = self.nergpair;
            self.ergpair[(n - 1) as usize] = energy;
            if jm == im {
                self.nergpairm[im] = 1;
                self.molergpair[(1, im)] = im;
                self.ergpairref[(1, im)] = n;
            } else {
                self.nergpairm[im] += 1;
                self.nergpairm[jm] += 1;
                let (ki, kj) = (self.nergpairm[im], self.nergpairm[jm]);
                self.molergpair[(ki, im)] = jm;
                self.molergpair[(kj, jm)] = im;
                self.ergpairref[(ki, im)] = n;
                self.ergpairref[(kj, jm)] = n;
            }
        }
    }
    /// `subroutine get_neighbsum(neighbsum)`: average number of stored pairs per molecule.
    pub fn get_neighbsum(&self) -> i32 {
        let mut s = 0;
        for im in 1..=self.nmolhere {
            s += self.nergpairm[im];
        }
        if self.nmolhere > 0 {
            s /= self.nmolhere;
        }
        s
    }
    /// `subroutine make_neighbour_list_from_mols(im)`: load the stored pair
    /// energies of molecule `im` into `nbref`/`ergold`/`ergnew`.
    pub fn make_neighbour_list_from_mols(&mut self, im: i32) {
        self.nnbref = self.nergpairm[im];
        let n = self.nnbref.max(0) as usize;
        if self.ergold.len() < n {
            self.ergold = F1::new(n);
            self.ergnew = F1::new(n);
            self.nbref = F1::new(n);
        }
        for ibref in 1..=self.nnbref {
            self.nbref[ibref] = self.molergpair[(ibref, im)];
            let e = self.ergpair[(self.ergpairref[(ibref, im)] - 1) as usize];
            self.ergold[ibref] = e;
            self.ergnew[ibref] = e;
        }
    }
    /// `subroutine save_energies()`: accept the trial energies.
    pub fn save_energies(&mut self) {
        for ibref in 1..=self.nnbref {
            self.ergold[ibref] = self.ergnew[ibref];
        }
    }
    /// `subroutine save_erg_pairs(im)`: write the (accepted) energies back to the store.
    pub fn save_erg_pairs(&mut self, im: i32) {
        for ibref in 1..=self.nnbref {
            let r = self.ergpairref[(ibref, im)];
            self.ergpair[(r - 1) as usize] = self.ergold[ibref];
        }
    }
}
