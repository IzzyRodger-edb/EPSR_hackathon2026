//! Pressure control (`control-p` keyword of the `.EPSR.inp` file).
//!
//! Translation of:
//! * the `control-p` parsing block and the "Adjust boxsize to control
//!   pressure" block at the end of the iteration loop of `runepsr17`
//!   (`epsr_standalone_rev3.f`)
//! * `change_boxsize`, `generate_fractional_coordinates`,
//!   `generate_cartesian_coordinates`   (`ato_manipulation_rev4.f`)
//! * `multiply_pcoefex`                 (`epsr_standalone_rev3.f`)
//!
//! `control-p <axes> <pressure> <step>`: after every iteration the box
//! lengths named by `<axes>` (`a`, `b`, `c`, `ab`, `bc`, `abc`/`all`) are
//! multiplied by `(1 + step*clip(pressd - pressure))**(1/naxes)`, or with
//! `p` the minimum distance `rminex` of every pair with a repulsive core is
//! divided by that factor.

use crate::error::Result;
use crate::fortran::format::ListWriter;
use crate::fortran::fmath::{cosd_intrinsic, sign, sind_intrinsic};
use crate::fortran::listread::ListRead;
use crate::fortran::text::parse;
use crate::neighbour::minimage;
use crate::state::{State, MPZ};
use crate::structure::uniform;

/// The pressure-control settings (locals of `runepsr17`).
#[derive(Clone, Debug, Default)]
pub struct PressureControl {
    pub control_pressure: bool,
    pub pacontrol: bool,
    pub pbcontrol: bool,
    pub pccontrol: bool,
    pub ppcontrol: bool,
    pub control_p_value: f32,
    pub control_p_step: f32,
    pub npcontrols: i32,
}

impl PressureControl {
    /// Parse the value of the `control-p` symbol.
    pub fn parse(sometext: &str) -> PressureControl {
        let mut c = PressureControl::default();
        let w = parse(sometext, 5);
        if w.nwords() > 0 {
            let w1 = w.word(sometext, 1);
            c.control_pressure = w.nwords() > 2 && ["a", "b", "c", "ab", "bc", "ca", "ac", "abc", "all", "p"].contains(&w1);
            if c.control_pressure {
                let s = w.span(sometext, 2, 3);
                let mut r = ListRead::from_str(&s);
                let ok = r.get_r(&mut c.control_p_value).and_then(|_| r.get_r(&mut c.control_p_step)).is_ok();
                c.control_pressure = ok && c.control_p_step > 0.0;
            }
            if c.control_pressure {
                // Bug fix (deliberate divergence): in the Fortran the list of
                // accepted words contains 'ca' but the branch tests 'ac', so
                // 'ca' selected no axis (npcontrols = 0) and 'ac' was not
                // accepted at all.  Both now select the c+a axes.
                let (a, b, cc, p) = match w1 {
                    "a" => (true, false, false, false),
                    "b" => (false, true, false, false),
                    "c" => (false, false, true, false),
                    "ab" => (true, true, false, false),
                    "bc" => (false, true, true, false),
                    "ac" | "ca" => (true, false, true, false),
                    "abc" | "all" => (true, true, true, false),
                    "p" => (false, false, false, true),
                    _ => (false, false, false, false),
                };
                c.pacontrol = a;
                c.pbcontrol = b;
                c.pccontrol = cc;
                c.ppcontrol = p;
            }
        }
        c.npcontrols = [c.pacontrol, c.pbcontrol, c.pccontrol, c.ppcontrol].iter().filter(|&&x| x).count() as i32;
        c
    }
}

impl State {
    /// "Adjust boxsize to control pressure if requested" (end of the
    /// iteration loop of `runepsr17`; coordinates are in Å here).
    pub fn control_pressure(&mut self, pc: &PressureControl, pressd: f32) -> Result<()> {
        if !(pc.control_pressure && pc.npcontrols > 0) {
            return Ok(());
        }
        let mut rat = pressd - pc.control_p_value;
        if rat.abs() > 1.0 {
            rat = sign(1.0, rat);
        }
        rat = pc.control_p_step * rat;
        let dmaxfac = (1.0f32 + rat).powf(1.0f32 / pc.npcontrols as f32);
        let (mut afac, mut bfac, mut cfac, mut pfac) = (1.0f32, 1.0f32, 1.0f32, 1.0f32);
        if pc.pacontrol {
            afac = dmaxfac;
        }
        if pc.pbcontrol {
            bfac = dmaxfac;
        }
        if pc.pccontrol {
            cfac = dmaxfac;
        }
        if pc.ppcontrol {
            pfac = 1.0 / dmaxfac;
        }
        self.log(&format!(" {afac} {bfac} {cfac} {pfac}"));
        if pc.ppcontrol {
            if self.pot.potfac > 0.0 {
                self.multiply_pcoefex(pfac)?;
            }
        } else {
            let (newphib, newthetac, newphic) = (self.lat.phib, self.lat.thetac, self.lat.phic);
            self.change_boxsize(afac, bfac, cfac, newphib, newthetac, newphic);
            let rmaxgr = self.inp.rmaxgr;
            uniform(&mut self.dist.ycoh, self.dist.cellst, self.dist.rho, rmaxgr, MPZ as i32);
        }
        Ok(())
    }

    /// `subroutine generate_fractional_coordinates(x,y,z)`.
    fn generate_fractional_coordinates(&self, x: &mut f32, y: &mut f32, z: &mut f32) {
        let l = &self.lat;
        let mut fc = *z / l.clatz;
        let mut fb = (*y - fc * l.claty) / l.blaty;
        let mut fa = (*x - fb * l.blatx - fc * l.clatx) / l.alatx;
        minimage(&mut fa, &mut fb, &mut fc);
        *z = fc;
        *y = fb;
        *x = fa;
    }

    /// `subroutine generate_cartesian_coordinates(x,y,z)`.
    fn generate_cartesian_coordinates(&self, x: &mut f32, y: &mut f32, z: &mut f32) {
        let l = &self.lat;
        let (fa, fb, fc) = (*x, *y, *z);
        *z = fc * l.clatz;
        *y = fb * l.blaty + fc * l.claty;
        *x = fa * l.alatx + fb * l.blatx + fc * l.clatx;
    }

    /// `subroutine change_boxsize(afac,bfac,cfac,newphib,newthetac,newphic)`:
    /// scale the lattice vectors and keep every molecule at the same
    /// fractional position.
    pub fn change_boxsize(&mut self, afac: f32, bfac: f32, cfac: f32, newphib: f32, newthetac: f32, newphic: f32) {
        let nmol = self.mol.nmol;
        // NB: phix/phiy/phiz are transformed for every molecule, tethered or not.
        for im in 1..=nmol {
            let m = &self.mol;
            let (mut x, mut y, mut z) = (m.xadd[im], m.yadd[im], m.zadd[im]);
            let (mut px, mut py, mut pz) = (m.phix[im], m.phiy[im], m.phiz[im]);
            self.generate_fractional_coordinates(&mut x, &mut y, &mut z);
            self.generate_fractional_coordinates(&mut px, &mut py, &mut pz);
            let m = &mut self.mol;
            m.xadd[im] = x;
            m.yadd[im] = y;
            m.zadd[im] = z;
            m.phix[im] = px;
            m.phiy[im] = py;
            m.phiz[im] = pz;
        }
        {
            // NB: change_boxsize does not declare cosd/sind EXTERNAL, so the
            // gfortran intrinsics are used (exact at 0 and 90 degrees).
            let l = &mut self.lat;
            l.alat *= afac;
            l.blat *= bfac;
            l.clat *= cfac;
            l.phib = newphib;
            l.thetac = newthetac;
            l.phic = newphic;
            l.alatx = l.alat;
            l.blatx = l.blat * cosd_intrinsic(l.phib);
            l.blaty = l.blat * sind_intrinsic(l.phib);
            l.clatx = l.clat * sind_intrinsic(l.thetac) * cosd_intrinsic(l.phic);
            l.claty = l.clat * sind_intrinsic(l.thetac) * sind_intrinsic(l.phic);
            l.clatz = l.clat * cosd_intrinsic(l.thetac);
        }
        for im in 1..=nmol {
            let m = &self.mol;
            let (mut x, mut y, mut z) = (m.xadd[im], m.yadd[im], m.zadd[im]);
            let (mut px, mut py, mut pz) = (m.phix[im], m.phiy[im], m.phiz[im]);
            self.generate_cartesian_coordinates(&mut x, &mut y, &mut z);
            self.generate_cartesian_coordinates(&mut px, &mut py, &mut pz);
            let m = &mut self.mol;
            m.xadd[im] = x;
            m.yadd[im] = y;
            m.zadd[im] = z;
            m.phix[im] = px;
            m.phiy[im] = py;
            m.phiz[im] = pz;
        }
        let l = &mut self.lat;
        l.alatc = [l.alatx, 0.0, 0.0];
        l.blatc = [l.blatx, l.blaty, 0.0];
        l.clatc = [l.clatx, l.claty, l.clatz];
        l.bxax = l.blatx / l.alatx;
        l.cxax = l.clatx / l.alatx;
        l.cyby = l.claty / l.blaty;
        let perpab = crate::ato::perpd(&l.alatc, &l.blatc, &l.clatc);
        let perpbc = crate::ato::perpd(&l.blatc, &l.clatc, &l.alatc);
        let perpca = crate::ato::perpd(&l.clatc, &l.alatc, &l.blatc);
        let msg = format!(" {perpab} {perpbc} {perpca}");
        self.mol.dmax = perpab.min(perpbc).min(perpca);
        self.dist.cell = self.mol.dmax;
        let l = &mut self.lat;
        l.cubic = !(l.alat != l.blat || l.blat != l.clat || l.clat != l.alat || l.phib != 90.0 || l.thetac != 0.0 || l.phic != 0.0);
        self.log(&msg);
        self.sync_lattice();
    }

    /// `subroutine multiply_pcoefex(epfact)`: scale the repulsive-core radius
    /// `rminex` of every pair with a non-zero repulsive core, update the pair
    /// symbol and regenerate the reference potential.
    pub fn multiply_pcoefex(&mut self, epfact: f32) -> Result<()> {
        for ip in 1..=self.mol.npairs {
            let irefp = self.pot.irefpcof[ip];
            let p = &mut self.pot;
            if p.pcoefex[irefp] != 0.0 && (p.pcoefex[irefp].abs() > 0.0001 || epfact > 1.0) && p.rminex[irefp] > 0.0 {
                p.rminex[irefp] *= epfact;
                let rminex = p.rminex[irefp];
                self.log(&format!(" multiply_pcoefex> {irefp:>11} {rminex} {epfact}"));
                let s = &mut self.sym;
                let icount = s.nsymboluptodata + s.nfiles * s.nsymbolplus + irefp;
                let sometext = s.symbol_value[icount].clone();
                let w = parse(&sometext, 5);
                let mut lw = ListWriter::new();
                lw.real(rminex);
                if w.nwords() > 1 {
                    let rest = w.span(&sometext, 2, w.nwords());
                    lw.chars(&rest);
                }
                let v = lw.finish();
                // the internal file is symbol_value (character*256)
                s.symbol_value[icount] = v;
            }
        }
        self.generate_rp()
    }
}
