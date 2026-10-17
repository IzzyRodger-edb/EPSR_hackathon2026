//! The molecular configuration: reading and writing `.ato` files and related
//! bookkeeping.
//!
//! Translation of:
//! * `read_ato`                     (`read_and_save_files.f`)
//! * `minitc`, `moloutc`, `normalise_coordinates`, `unnormalise_coordinates`
//!                                  (`ato_read_and_save.f`)
//! * `ftype`, `ctype`, `centreofmassi`, `get_maxmolsize`, `fractcoord`,
//!   `wholecoord`                   (`ato_manipulation_rev4.f`)
//! * `bondnumbers`                  (`bondnumbers.f`)
//! * `perpd`, `vecpro`, `scapro`, `cosd`, `sind`, `get_atomicnumber`
//!                                  (`auxiliary_routines.f`)
//!
//! ## `.ato` file layout (as read by `minitc`)
//! ```text
//!  nmol  dmax  temperature                      (cubic box; 9 numbers for a general cell)
//!  tol stepmi stepri stepr stepm tempb [tempangreq tempdihreq [tempangwt tempdihwt]]
//!  for each molecule:
//!     nato xcentre ycentre zcentre phix phiy phiz  T/F nstepreject scalefac imol
//!     for each atom:  label  index  imageatom
//!                     x y z                       (relative to the molecule centre, Å)
//!                     nbond (partner, bondlength) x nbond   (negative length = angle restraint)
//!     ngrot, then for each rotational group: ROT|DIH / axis atoms / nrot + atoms
//!  for each atom type: label element iexchange / eps sigma mass charge qradius ...
//!  ecore dcore
//!  random-number-generator state (mcall idum iy iv(32))
//!  for each molecule type: index template step sizes ...
//! ```
//! After reading, all lengths are converted to fractions of the box size
//! `dmax` (`normalise_coordinates`), which is the internal unit of the
//! Monte Carlo code.

use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::{columns, read_int, read_real, ListRead, TextFile};
use crate::fortran::output::OutFile;
use crate::fortran::text::{concatenate, fstr, indexlastdot, parse};
use crate::fortran::F1;
use crate::rng::{Ran1, NTAB};
use crate::state::{State, NEAR};
/// gfortran intrinsic SIND/COSD (used by routines that do not declare `external cosd,sind`).
pub use crate::fortran::fmath::{cosd_intrinsic, sind_intrinsic};

/// `degradconv = 3.141592653/180.0` (single precision parameter).
fn degradconv() -> f32 {
    3.141592653f32 / 180.0f32
}
/// `function cosd(deg)`.
pub fn cosd(deg: f32) -> f32 {
    (deg * degradconv()).cos()
}
/// `function sind(deg)`.
pub fn sind(deg: f32) -> f32 {
    (deg * degradconv()).sin()
}
/// `function scapro(a,b)`: scalar product.
pub fn scapro(a: &[f32; 3], b: &[f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// `subroutine vecpro(a,b,c)`: vector product c = a × b.
pub fn vecpro(a: &[f32; 3], b: &[f32; 3]) -> [f32; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
/// `function perpd(a,b,c)`: perpendicular distance between the faces spanned
/// by a,b along c.
pub fn perpd(a: &[f32; 3], b: &[f32; 3], c: &[f32; 3]) -> f32 {
    let mut t = vecpro(a, b);
    let amp = scapro(&t, &t).sqrt();
    for v in t.iter_mut() {
        *v /= amp;
    }
    scapro(c, &t)
}

/// Chemical symbols of the elements (`chemical_symbols.inc`), used for
/// `get_atomicnumber`.
const ELEMENTS: &[&str] = &[
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti",
    "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",
    "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
    "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb",
    "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No",
    "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

/// `function get_atomicnumber(ats)` (0 if the symbol is unknown).
pub fn get_atomicnumber(ats: &str) -> i32 {
    ELEMENTS.iter().position(|e| *e == ats).map(|p| p as i32 + 1).unwrap_or(0)
}

impl State {
    /// `subroutine read_ato`: read the `.ato` file named by symbol `icount`.
    pub fn read_ato(&mut self, icount: i32, ifilecount: i32) -> Result<()> {
        let ifound = self.sym.ireff(icount, ifilecount);
        let isf = self.sym.search_for_symbol("sizefactor", 0);
        if isf > 0 {
            let iref = self.sym.ireff(isf, ifilecount);
            self.mol.sizefactor = read_real(&self.sym.symbol_value[iref])?;
        } else {
            self.mol.sizefactor = 1.0;
        }
        let v = self.sym.symbol_value[ifound].clone();
        if v == "search" || v == "<undefined>" {
            return Err(EpsrError::Input("no .ato file specified (fnameato)".into()));
        }
        let w = parse(&v, 5);
        let (mut fname_ato, _) = concatenate(1, w.nwords(), &w, &v, 256);
        let ilenlastdot = indexlastdot(&fname_ato) - 1;
        if ilenlastdot > 0 {
            fname_ato = format!("{}.ato", &fname_ato[..ilenlastdot as usize]);
        } else {
            fname_ato.push_str(".ato");
        }
        self.sym.fname_ato = fname_ato.clone();
        self.sym.fnameato = self.sym.full(&fname_ato).to_string_lossy().into_owned();
        self.log(&format!("\nAttempting to read file: {}\n", self.sym.fnameato));
        self.mol.ntype = 0;
        self.mol.nmoltype = 0;
        let path = self.sym.fnameato.clone();
        self.minitc(&path)?;
        self.sym.symbol_value[ifound] = fname_ato;
        let irho = self.sym.search_for_symbol("rho", 0);
        if irho != 0 {
            let iref = self.sym.ireff(irho, ifilecount);
            let mut w = ListWriter::new();
            w.real(self.dist.rho);
            self.sym.symbol_value[iref] = w.finish().trim_end().to_string();
        }
        Ok(())
    }

    /// `subroutine ftype(ia)`: assign the atom type of atom `ia` from its label,
    /// creating a new type if needed.
    pub fn ftype(&mut self, ia: i32) {
        let m = &mut self.mol;
        let lab = m.alabel[ia].clone();
        let first = lab.chars().next().unwrap_or(' ');
        m.qtype[ia] = if first == 'Q' || first == 'q' { 1 } else { 0 };
        let mut ifound = 0;
        let mut it = 0;
        while it < m.ntype && ifound == 0 {
            it += 1;
            if lab == m.type_[it] {
                ifound = it;
            }
        }
        if ifound > 0 {
            m.itype[ia] = ifound;
        } else if (m.ntype as usize) < NEAR {
            m.ntype += 1;
            let nt = m.ntype;
            m.type_[nt] = lab;
            m.qtypea[nt] = m.qtype[ia];
            m.itype[ia] = nt;
        } else {
            self.log(" Tried to generate more than 200 types of atom");
        }
    }

    /// `subroutine ctype(...)`: store the potential parameters of atom type `a1`.
    #[allow(clippy::too_many_arguments)]
    pub fn ctype(
        &mut self,
        a1: &str,
        ats: &str,
        iexch: i32,
        eread: f32,
        wread: f32,
        amasst: f32,
        qt: f32,
        rqt: f32,
        natomsph: f32,
        rhosph: f32,
        radminsph: f32,
        radmaxsph: f32,
        pseudofacsph: f32,
    ) {
        let m = &mut self.mol;
        let nt = m.ntype.max(0) as usize;
        m.natomsphere.grow(nt);
        m.rhosphere.grow(nt);
        m.radminsphere.grow(nt);
        m.radmaxsphere.grow(nt);
        m.pseudofacmol.grow(nt);
        let mut n1 = 0;
        let mut it = 0;
        while it < m.ntype && n1 == 0 {
            it += 1;
            if a1 == m.type_[it] {
                n1 = it;
            }
        }
        if n1 > 0 {
            m.iatomicnotype[n1] = get_atomicnumber(ats);
            m.atsymbol[n1] = ats.to_string();
            m.iexchange[n1] = iexch;
            m.epscor[n1] = eread;
            m.sigcor[n1] = wread;
            m.amassa[n1] = amasst;
            m.charge[n1] = qt;
            m.qradius[n1] = rqt;
            m.natomsphere[n1] = natomsph;
            m.rhosphere[n1] = rhosph;
            m.radminsphere[n1] = radminsph;
            m.radmaxsphere[n1] = radmaxsph;
            m.pseudofacmol[n1] = pseudofacsph;
        } else {
            self.log(&format!(" Mismatched types between atomic coordinates and potential values: {a1}"));
        }
    }

    /// `subroutine bondnumbers`: renumber bond partners, rotation axes and
    /// rotating atoms from per-molecule numbering to global atom numbers.
    pub fn bondnumbers(&mut self) -> Result<()> {
        let m = &mut self.mol;
        let mut ig = 0;
        for im in 1..=m.nmol {
            let mut namin = 0;
            for ia in m.naf[im]..=m.nal[im] {
                if m.nbond[ia] > 0 {
                    for ib in 1..=m.nbond[ia] {
                        let nt = m.iabond[(ia, ib)];
                        if namin == 0 {
                            namin = nt;
                        }
                        if nt < namin {
                            namin = nt;
                        }
                    }
                }
            }
            if namin != 0 && namin != m.naf[im] {
                for ia in m.naf[im]..=m.nal[im] {
                    if m.imageatom[ia] > 0 {
                        m.imageatom[ia] = m.imageatom[ia] - namin + m.naf[im];
                    }
                    if m.nbond[ia] > 0 {
                        for ib in 1..=m.nbond[ia] {
                            let nt = m.iabond[(ia, ib)] - namin;
                            m.iabond[(ia, ib)] = m.naf[im] + nt;
                        }
                    }
                }
            }
            if m.ngrot[im] > 0 {
                for _ih in 1..=m.ngrot[im] {
                    ig += 1;
                    if namin != 0 && namin != m.naf[im] {
                        m.labaxi[(ig, 1)] = m.labaxi[(ig, 1)] - namin + m.naf[im];
                        m.labaxi[(ig, 2)] = m.labaxi[(ig, 2)] - namin + m.naf[im];
                        for j in m.rotf[ig]..=m.rotl[ig] {
                            m.labrot[j] = m.labrot[j] - namin + m.naf[im];
                        }
                        if m.rotgrp[ig] == "DIH" {
                            for k in 1..=m.ndihedralang[ig] {
                                for j in 1..=2 {
                                    m.labrotdih[(j, k, ig)] = m.labrotdih[(j, k, ig)] - namin + m.naf[im];
                                }
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// `subroutine minitc(fname1, rmax, cell, upol, rho, amubar, amubar2, ndist, ierr)`:
    /// read a `.ato` file and derive densities, atom fractions, pair tables
    /// and dipole statistics.
    pub fn minitc(&mut self, fname1: &str) -> Result<()> {
        let pi: f32 = std::f32::consts::PI;
        let mut f = TextFile::open(std::path::Path::new(fname1))
            .map_err(|e| EpsrError::Io(format!("Problems with file {fname1}: {e}")))?;
        // ---- header
        self.lat.cubic = true;
        let sometext = f.read_line()?;
        let w = parse(&sometext, 5);
        if w.nwords() != 3 {
            // general cell: nmol, temperature / a b c / phib thetac phic
            self.lat.cubic = false;
            {
                let mut r = ListRead::from_str(&sometext);
                let _ = r.get_i(&mut self.mol.nmol).and_then(|_| r.get_r(&mut self.mol.simtk));
            }
            let l = &mut self.lat;
            {
                let mut r = f.list();
                let _ = r.get_r(&mut l.alat).and_then(|_| r.get_r(&mut l.blat)).and_then(|_| r.get_r(&mut l.clat));
            }
            {
                let mut r = f.list();
                let _ = r.get_r(&mut l.phib).and_then(|_| r.get_r(&mut l.thetac)).and_then(|_| r.get_r(&mut l.phic));
            }
        } else {
            {
                let mut r = ListRead::from_str(&sometext);
                self.mol.nmol = r.int()?.unwrap_or(0);
                self.mol.dmax = r.real()?.unwrap_or(0.0);
                self.mol.simtk = r.real()?.unwrap_or(0.0);
            }
            let dmax = self.mol.dmax;
            self.lat.alat = dmax;
            self.lat.blat = dmax;
            self.lat.clat = dmax;
            self.lat.phib = 90.0;
            self.lat.thetac = 0.0;
            self.lat.phic = 0.0;
        }
        self.log(&format!("minitc> cubic: {:>10}", if self.lat.cubic { "T" } else { "F" }));
        let nmol = self.mol.nmol;
        if nmol as usize > crate::state::MMOL {
            return Err(EpsrError::Input(format!("too many molecules ({nmol})")));
        }
        self.mol.scalefacmol = F1::new(nmol.max(1) as usize);
        self.mol.scalefacmol.fill(1.0);
        // ---- step sizes and temperatures
        let sometext = f.read_line()?;
        let w = parse(&sometext, 5);
        {
            let m = &mut self.mol;
            let mut r = ListRead::from_str(&sometext);
            let _ = (|| -> crate::fortran::IoRes<()> {
                r.get_r(&mut m.tol)?;
                r.get_r(&mut m.stepmi)?;
                r.get_r(&mut m.stepri)?;
                r.get_r(&mut m.stepr)?;
                r.get_r(&mut m.stepm)?;
                r.get_r(&mut m.tempb)?;
                Ok(())
            })();
            m.tempangreq = 1.0;
            m.tempdihreq = 10.0;
            m.tempangwt = 1.0;
            m.tempdihwt = 1.0;
            if w.nwords() == 7 {
                let s = w.span(&sometext, 7, 7);
                let mut r = ListRead::from_str(&s);
                let _ = r.get_r(&mut m.tempdihreq);
            }
            if w.nwords() > 7 {
                let s = w.span(&sometext, 7, 8);
                let mut r = ListRead::from_str(&s);
                let _ = r.get_r(&mut m.tempangreq).and_then(|_| r.get_r(&mut m.tempdihreq));
            }
            if w.nwords() > 9 {
                let s = w.span(&sometext, 9, 10);
                let mut r = ListRead::from_str(&s);
                let _ = r.get_r(&mut m.tempangwt).and_then(|_| r.get_r(&mut m.tempdihwt));
            }
            m.tempang = m.tempangwt;
            m.tempdih = m.tempdihwt * m.tempb;
            m.tolout = m.tol;
        }
        // ---- molecules
        let mut ia: i32 = 0;
        let mut ig: i32 = 0;
        self.mol.natom = 0;
        self.mol.natomr = 0;
        self.mol.ngrrot = 0;
        self.mol.nrottot = 0;
        for im in 1..=nmol {
            let sometext = f.read_line()?;
            let w = parse(&sometext, 10);
            {
                let m = &mut self.mol;
                let mut r = ListRead::from_str(&sometext);
                let _ = (|| -> crate::fortran::IoRes<()> {
                    r.get_i(&mut m.nato[im])?;
                    r.get_r(&mut m.xadd[im])?;
                    r.get_r(&mut m.yadd[im])?;
                    r.get_r(&mut m.zadd[im])?;
                    r.get_r(&mut m.phix[im])?;
                    r.get_r(&mut m.phiy[im])?;
                    r.get_r(&mut m.phiz[im])?;
                    Ok(())
                })();
                m.tetheredatom[im] = 0;
                if w.nwords() > 7 {
                    let t = w.word(&sometext, 8);
                    m.tethered[im] = t.starts_with('T');
                    if t.len() > 1 {
                        m.tetheredatom[im] = read_int(&t[1..]).unwrap_or(0);
                    }

                } else {
                    m.tethered[im] = false;
                }
                m.nstepreject[im] = if w.nwords() > 8 { read_int(w.word(&sometext, 9))? } else { 0 };
                if w.nwords() > 9 {
                    m.scalefacmol[im] = read_real(w.word(&sometext, 10))?;
                } else {
                    m.scalefacmol[im] = 1.0;
                }
                m.natom += m.nato[im];
                m.naf[im] = ia + 1;
            }
            let nato = self.mol.nato[im];
            self.mol.ensure_atoms((ia + nato) as usize);
            for _ib in 1..=nato {
                ia += 1;
                self.mol.mref[ia] = im;
                let sometext = f.read_line()?;
                let w = parse(&sometext, 10);
                self.mol.alabel[ia] = fstr(&columns(w.word(&sometext, 1), 0, 3), 3);
                if w.nwords() > 2 {
                    let mut v = read_int(w.word(&sometext, 3)).unwrap_or(0);
                    if v < 0 || v > nato {
                        v = 0;
                    }
                    self.mol.imageatom[ia] = v;
                } else {
                    self.mol.imageatom[ia] = 0;
                }
                self.ftype(ia);
                if self.mol.qtype[ia] != 1 {
                    self.mol.natomr += 1;
                }
                {
                    let m = &mut self.mol;
                    let mut r = f.list();
                    r.get_r(&mut m.xa[ia])?;
                    r.get_r(&mut m.ya[ia])?;
                    r.get_r(&mut m.za[ia])?;
                }
                {
                    let mut r = f.list();
                    let mut nb = 0;
                    r.get_i(&mut nb)?;
                    self.mol.nbond[ia] = nb;
                    self.mol.ensure_bonds(nb.max(0) as usize);
                    for j in 1..=nb {
                        let m = &mut self.mol;
                        r.get_i(&mut m.iabond[(ia, j)])?;
                        r.get_r(&mut m.tbond[(ia, j)])?;
                    }
                }
            }
            self.mol.nal[im] = ia;
            {
                let mut r = f.list();
                r.get_i(&mut self.mol.ngrot[im])?;
            }
            self.mol.nrf[im] = ig + 1;
            if self.mol.ngrot[im] > 0 {
                self.mol.ngrrot += self.mol.ngrot[im];
                let n = self.mol.ngrrot as usize;
                self.mol.ensure_rotgrps(n);
                for _ih in 1..=self.mol.ngrot[im] {
                    ig += 1;
                    let m = &mut self.mol;
                    m.mrot[ig] = im;
                    let l = f.read_line()?;
                    m.rotgrp[ig] = columns(&l, 1, 3);
                    {
                        let mut r = f.list();
                        r.get_i(&mut m.labaxi[(ig, 1)])?;
                        r.get_i(&mut m.labaxi[(ig, 2)])?;
                    }
                    if m.rotgrp[ig] == "DIH" {
                        let mut r = f.list();
                        r.get_i(&mut m.ndihedralang[ig])?;
                        drop(r);
                        for k in 1..=m.ndihedralang[ig] {
                            let mut r = f.list();
                            r.get_i(&mut m.labrotdih[(1, k, ig)])?;
                            r.get_i(&mut m.labrotdih[(2, k, ig)])?;
                            r.get_r(&mut m.dihedralang[(k, ig)])?;
                        }
                    } else {
                        m.ndihedralang[ig] = 0;
                    }
                    let sometext = f.read_line()?;
                    let w = parse(&sometext, 5);
                    m.nrot[ig] = read_int(w.word(&sometext, 1))?;
                    m.rotf[ig] = m.nrottot + 1;
                    m.nrottot += m.nrot[ig];
                    m.rotl[ig] = m.nrottot;
                    let nr = m.nrottot as usize;
                    m.ensure_labrot(nr);
                    let mut iw = 1;
                    let mut labref = m.rotf[ig];
                    while iw < w.nwords() {
                        iw += 1;
                        labref = m.rotf[ig] + iw - 2;
                        m.labrot[labref] = read_int(w.word(&sometext, iw))?;
                    }
                    if labref < m.rotl[ig] {
                        return Err(EpsrError::Input(format!("rotational group {ig}: fewer atoms listed than nrot")));
                    }
                }
            }
            self.mol.nrl[im] = ig;
        }
        self.bondnumbers()?;
        // ---- distinct atom labels in file order
        {
            let m = &mut self.mol;
            m.ntypefile = 0;
            let mut ia = 0;
            while ia < m.natom {
                ia += 1;
                let mut itypefile = 0;
                let mut ifound = 0;
                while itypefile < m.ntypefile && ifound == 0 {
                    itypefile += 1;
                    if m.alabel[ia] == m.typefile[itypefile] {
                        ifound = itypefile;
                    }
                }
                if ifound == 0 {
                    m.ntypefile += 1;
                    let n = m.ntypefile;
                    m.typefile[n] = m.alabel[ia].clone();
                }
            }
        }
        // ---- molecule types (identified by the label of the first atom)
        let nmoltypefirst = self.mol.nmoltype;
        {
            let m = &mut self.mol;
            let mut im = 0;
            while im < m.nmol {
                im += 1;
                let mut imoltype = 0;
                let mut ifound = 0;
                while imoltype < m.nmoltype && ifound == 0 {
                    imoltype += 1;
                    if m.alabel[m.naf[im]] == m.moltype[imoltype] {
                        ifound = imoltype;
                        m.mtype[im] = ifound;
                    }
                }
                if ifound == 0 {
                    m.nmoltype += 1;
                    let nt = m.nmoltype;
                    m.moltype[nt] = m.alabel[m.naf[im]].clone();
                    m.moltethered[nt] = m.tethered[im];
                    m.moltetheredatom[nt] = m.tetheredatom[im];
                    m.moltemplate[nt] = fmt("('moltype',i2.2)", &[FV::I(nt as i64)]);
                    m.mtype[im] = nt;
                    m.wholerot[nt] = m.stepr;
                    m.wholetrans[nt] = m.stepm;
                }
            }
            m.rfractmol = F1::new(m.nmoltype.max(1) as usize);
        }
        // ---- atom type parameters
        for _it in 1..=self.mol.ntypefile {
            let sometext = f.read_line()?;
            let w = parse(&sometext, 8);
            let a1 = fstr(w.word(&sometext, 1), 3);
            let ats = if w.nwords() > 1 { fstr(w.word(&sometext, 2), 3) } else { "XXX".to_string() };
            let iexch = if w.nwords() > 2 { read_int(w.word(&sometext, 3))? } else { -1 };
            let sometext = f.read_line()?;
            let w = parse(&sometext, 11);
            let (mut epsit, mut sigt, mut amass, mut qt, mut rqt) = (0f32, 0f32, 0f32, 0f32, 0f32);
            {
                let s = w.span(&sometext, 1, 5);
                let mut r = ListRead::from_str(&s);
                let _ = (|| -> crate::fortran::IoRes<()> {
                    r.get_r(&mut epsit)?;
                    r.get_r(&mut sigt)?;
                    r.get_r(&mut amass)?;
                    r.get_r(&mut qt)?;
                    r.get_r(&mut rqt)?;
                    Ok(())
                })();
            }
            let mut natomsph: f32 = 1.0;
            if w.nwords() > 5 {
                if let Ok(v) = read_real(w.word(&sometext, 6)) {
                    natomsph = v;
                }
            } else {
                natomsph = 1.0;
            }
            let (mut radminsph, mut radmaxsph): (f32, f32);
            if w.nwords() > 7 {
                radminsph = 0.0;
                radmaxsph = 0.0;
                let s = w.span(&sometext, 7, 8);
                let mut r = ListRead::from_str(&s);
                let _ = r.get_r(&mut radminsph).and_then(|_| r.get_r(&mut radmaxsph));
                if radmaxsph <= 0.0 {
                    radmaxsph = 0.5 * sigt;
                }
            } else {
                radminsph = 0.0;
                radmaxsph = 0.5 * sigt;
            }
            let mut pseudofacsph: f32 = 1.0;
            if w.nwords() > 8 {
                if let Ok(v) = read_real(w.word(&sometext, 9)) {
                    pseudofacsph = v;
                }
                if pseudofacsph < 1.0 {
                    pseudofacsph = 1.0;
                }
            }
            let volsph = 4.0 * pi * (crate::fortran::fmath::pow3(radmaxsph) - crate::fortran::fmath::pow3(radminsph)) / 3.0;
            let rhosph: f32 = if natomsph.abs() > 0.0 {
                if natomsph < 0.0 {
                    let rhosph = natomsph.abs();
                    natomsph = rhosph * volsph / pseudofacsph;
                    rhosph
                } else if volsph > 0.0 {
                    (natomsph * pseudofacsph) / volsph
                } else {
                    0.0
                }
            } else {
                0.0
            };
            self.ctype(&a1, &ats, iexch, epsit, sigt, amass, qt, rqt, natomsph, rhosph, radminsph, radmaxsph, pseudofacsph);
        }
        // ---- core parameters, RNG state, molecule-type step sizes
        let core = {
            let mut r = f.list();
            let mut ecore = 0.0f32;
            let mut dcore = 0.0f32;
            r.get_r(&mut ecore).and_then(|_| r.get_r(&mut dcore)).map(|_| (ecore, dcore))
        };
        match core {
            Ok((e, d)) => {
                self.mol.ecore = e;
                self.mol.dcore = d;
                // initran1(10)
                let st = {
                    let mut r = f.list();
                    (|| -> crate::fortran::IoRes<(i32, i32, i32, [i32; NTAB])> {
                        let mcall = r.int()?.unwrap_or(0);
                        let idum = r.int()?.unwrap_or(0);
                        let iy = r.int()?.unwrap_or(0);
                        let mut iv = [0i32; NTAB];
                        for v in iv.iter_mut() {
                            *v = r.int()?.unwrap_or(0);
                        }
                        Ok((mcall, idum, iy, iv))
                    })()
                };
                self.rng = Ran1::from_record(st.ok());
                let mut icount = 0;
                while icount < self.mol.nmoltype - nmoltypefirst {
                    icount += 1;
                    let sometext = match f.read_line() {
                        Ok(s) => s,
                        Err(_) => break,
                    };
                    let w = parse(&sometext, 11);
                    if w.nwords() > 1 {
                        let m = &mut self.mol;
                        let imolread = read_int(w.word(&sometext, 1)).unwrap_or(0);
                        let imoltype = nmoltypefirst + imolread;
                        m.moltemplate[imoltype] = w.word(&sometext, 2).to_string();
                        m.interntrans[imoltype] = 0.0;
                        m.internrot[imoltype] = 0.0;
                        m.wholerot[imoltype] = 0.0;
                        m.wholetrans[imoltype] = 0.0;
                        let rd = |k: i32| read_real(w.word(&sometext, k)).ok();
                        if w.nwords() > 2 {
                            if let Some(v) = rd(3) {
                                m.interntrans[imoltype] = v;
                            }
                        }
                        if w.nwords() > 3 {
                            if let Some(v) = rd(4) {
                                m.internrot[imoltype] = v;
                            }
                        }
                        if w.nwords() > 4 {
                            if let Some(v) = rd(5) {
                                m.wholerot[imoltype] = v;
                            }
                        }
                        if w.nwords() > 5 {
                            if let Some(v) = rd(6) {
                                m.wholetrans[imoltype] = v;
                            }
                        }
                        if w.nwords() > 6 {
                            if let Some(v) = rd(7) {
                                m.rfractmol[imoltype] = v;
                            }
                        }
                    }
                }
            }
            Err(_) => {
                // The Fortran minitc returns here with ierr /= 0: read_ato prints "No atoms
                // were read from this file", restores ntype, and epsr carries on with an
                // unset box (dmax = 0, no lattice, no pair tables) and hangs in
                // initialise_cellneighbourlist.  Not reproducible output, so not ported.
                return crate::error::unsupported(
                    ".ato file without ecore/dcore and random number state (the Fortran epsr ignores the configuration and hangs)",
                );
            }
        }
        // ---- cell vectors
        let l = &mut self.lat;
        l.alatx = l.alat;
        l.blatx = l.blat * cosd(l.phib);
        l.blaty = l.blat * sind(l.phib);
        l.clatx = l.clat * sind(l.thetac) * cosd(l.phic);
        l.claty = l.clat * sind(l.thetac) * sind(l.phic);
        l.clatz = l.clat * cosd(l.thetac);
        l.alatc = [l.alatx, 0.0, 0.0];
        l.blatc = [l.blatx, l.blaty, 0.0];
        l.clatc = [l.clatx, l.claty, l.clatz];
        l.bxax = l.blatx / l.alatx;
        l.cxax = l.clatx / l.alatx;
        l.cyby = l.claty / l.blaty;
        let perpab = perpd(&l.alatc, &l.blatc, &l.clatc);
        let perpbc = perpd(&l.blatc, &l.clatc, &l.alatc);
        let perpca = perpd(&l.clatc, &l.alatc, &l.blatc);
        let m = &mut self.mol;
        m.dmax = perpab.min(perpbc).min(perpca);

        // minitc dummy argument `cell` is the common variable cellmax
        self.dist.cellmax = m.dmax;
        m.dmaxmol = m.dmax * m.sizefactor;
        let xmax = l.alatx.abs() + l.blatx.abs() + l.clatx.abs();
        let ymax = l.blaty.abs() + l.claty.abs();
        let zmax = l.clatz.abs();
        // dummy argument `rmax` is the common variable rmaxcell
        self.dist.rmaxcell = 0.5 * (xmax * xmax + ymax * ymax + zmax * zmax).sqrt();
        m.simt = 0.008314 * m.simtk;
        // ---- atom fractions and density
        m.natomcg = 0.0;
        m.natomrcg = 0.0;
        for it in 1..=m.ntype {
            m.cfrac[it] = 0.0;
            m.cfraccg[it] = 0.0;
        }
        for ia in 1..=m.natom {
            let it = m.itype[ia];
            m.natomcg += m.natomsphere[it] * m.pseudofacmol[it];
            if m.qtype[ia] != 1 {
                m.natomrcg += m.natomsphere[it] * m.pseudofacmol[it];
            }
            m.cfraccg[it] += m.natomsphere[it] * m.pseudofacmol[it];
            m.cfrac[it] += 1.0;
        }
        m.volume = l.alatx * l.blaty * l.clatz;
        self.dist.rho = m.natomr as f32 / m.volume;
        m.rhocg = m.natomrcg / m.volume;
        if m.natomrcg > 0.0 {
            for it in 1..=m.ntype {
                m.cfraccg[it] /= m.natomrcg;
                m.cfrac[it] /= m.natomr as f32;
            }
        } else if m.natomcg > 0.0 {
            for it in 1..=m.ntype {
                m.cfraccg[it] /= m.natomcg;
                m.cfrac[it] /= m.natom as f32;
            }
        }
        // NB: with natomrcg = natomcg = 0 (e.g. every atom type has natomsphere = 0)
        // the fractions are left as plain counts, as in the Fortran.
        m.sumcharge = 0.0;
        for it in 1..=m.ntype {
            m.sumcharge += m.cfrac[it] * m.charge[it];
        }
        for im in 1..=m.nmol {
            m.amasst[im] = 0.0;
            for ia in m.naf[im]..=m.nal[im] {
                if m.qtype[ia] == 0 {
                    let amass = m.amassa[m.itype[ia]];
                    m.amasst[im] += amass;
                }
            }
            m.amasstsqrt[im] = m.amasst[im].sqrt();
            m.moldiam[im] = 0.0;
            for ia in m.naf[im]..=m.nal[im] {
                let sum = m.xa[ia] * m.xa[ia] + m.ya[ia] * m.ya[ia] + m.za[ia] * m.za[ia];
                if sum > m.moldiam[im] {
                    m.moldiam[im] = sum;
                }
            }
            m.moldiam[im] = m.moldiam[im].sqrt();
        }
        // ---- pair tables
        let mut ip = 0;
        for it in 1..=m.ntype {
            for jt in it..=m.ntype {
                ip += 1;
                m.n1extr[ip] = it;
                m.n2extr[ip] = jt;
                m.npairref[(it, jt)] = ip;
                m.npairref[(jt, it)] = ip;
                let ect = (m.epscor[it] * m.epscor[jt]).sqrt();
                m.ecor[ip] = ect;
                m.qcor[ip] = m.charge[it] * m.charge[jt] * 1390.0;
                let dct = 0.5 * (m.sigcor[it] + m.sigcor[jt]);
                m.dcor[ip] = dct;
                m.dscor[ip] = dct * dct;
                let am1 = m.amassa[it];
                let am2 = m.amassa[jt];
                m.amassr[ip] = (am1 * am2 / (am1 + am2)).sqrt();
            }
        }
        m.npairs = ip;
        self.dist.ndist = ip;
        // ---- dipole moments
        let mut amu0: f32 = 1.85;
        amu0 /= 4.80273;
        let mut amusum: f32 = 0.0;
        let d = &mut self.dist;
        d.amubar = 0.0;
        let mut nbar = 0;
        for im in 1..=m.nmol {
            if m.nato[im] > 1 {
                let (mut qx, mut qy, mut qz) = (0f32, 0f32, 0f32);
                for ia in m.naf[im]..=m.nal[im] {
                    let qt = m.charge[m.itype[ia]];
                    qx += qt * m.xa[ia];
                    qy += qt * m.ya[ia];
                    qz += qt * m.za[ia];
                }
                let amu2 = qx * qx + qy * qy + qz * qz;
                let amu = amu2.sqrt();
                let amud = amu * 4.80273;
                nbar += 1;
                d.amubar += amud;
                d.amubar2 += amud * amud;
                let amudif = amu - amu0;
                amusum += amudif * amudif;
            }
        }
        if nbar > 0 {
            d.amubar /= nbar as f32;
            d.amubar2 = d.amubar2 / nbar as f32 - d.amubar * d.amubar;
            if d.amubar2 > 0.0 {
                d.amubar2 = d.amubar2.sqrt();
            } else {
                d.amubar2 = 0.0;
            }
            amusum /= m.nmol as f32;
            d.upol = 0.5 * amusum * 1390.0 / 1.445;
        } else {
            d.amubar = 0.0;
            d.amubar2 = 0.0;
            d.upol = 0.0;
        }
        self.sync_lattice();
        Ok(())
    }

    /// `subroutine normalise_coordinates`: convert lengths to box-fraction units.
    pub fn normalise_coordinates(&mut self) {
        let m = &mut self.mol;
        let l = &mut self.lat;
        let dmax = m.dmax;
        let dmaxmol = m.dmaxmol;
        m.stepmi /= dmaxmol;
        m.stepm /= dmax;
        for v in [&mut l.alatx, &mut l.blatx, &mut l.clatx, &mut l.blaty, &mut l.claty, &mut l.clatz, &mut l.alat, &mut l.blat, &mut l.clat] {
            *v /= dmax;
        }
        for i in 0..3 {
            l.alatc[i] /= dmax;
            l.blatc[i] /= dmax;
            l.clatc[i] /= dmax;
        }
        if m.stepm > 0.25 {
            m.stepm = 0.25;
        }
        for im in 1..=m.nmol {
            m.xadd[im] /= dmax;
            m.yadd[im] /= dmax;
            m.zadd[im] /= dmax;
            if m.tethered[im] {
                m.phix[im] /= dmax;
                m.phiy[im] /= dmax;
                m.phiz[im] /= dmax;
            }
            m.moldiam[im] /= dmaxmol;
            for ia in m.naf[im]..=m.nal[im] {
                m.xa[ia] /= dmaxmol;
                m.ya[ia] /= dmaxmol;
                m.za[ia] /= dmaxmol;
                for ib in 1..=m.nbond[ia] {
                    m.tbond[(ia, ib)] /= dmaxmol;
                }
            }
        }
        for ip in 1..=m.npairs {
            m.amassr[ip] *= dmaxmol;
        }
        let piconv: f32 = std::f32::consts::PI / 180.0;
        for ir in 1..=m.ngrrot {
            if m.rotgrp[ir] == "DIH" {
                for k in 1..=m.ndihedralang[ir] {
                    m.dihedralang[(k, ir)] *= piconv;
                }
            }
        }
        // three-body maximum distances
        for ir in 1..=self.tb.nthreeb {
            for j in 1..=2 {
                self.tb.threebrmax[(j, ir)] /= dmax;
            }
        }
        self.scale_coordconstraint_lengths(true);
        self.sync_lattice();
    }

    /// `subroutine unnormalise_coordinates`: back to Å.
    pub fn unnormalise_coordinates(&mut self) {
        let m = &mut self.mol;
        let l = &mut self.lat;
        let dmax = m.dmax;
        let dmaxmol = m.dmaxmol;
        m.stepmi *= dmaxmol;
        if m.stepm > 0.25 {
            m.stepm = 0.25;
        }
        m.stepm *= dmax;
        for v in [&mut l.alatx, &mut l.blatx, &mut l.clatx, &mut l.blaty, &mut l.claty, &mut l.clatz, &mut l.alat, &mut l.blat, &mut l.clat] {
            *v *= dmax;
        }
        for i in 0..3 {
            l.alatc[i] *= dmax;
            l.blatc[i] *= dmax;
            l.clatc[i] *= dmax;
        }
        for im in 1..=m.nmol {
            m.xadd[im] *= dmax;
            m.yadd[im] *= dmax;
            m.zadd[im] *= dmax;
            if m.tethered[im] {
                m.phix[im] *= dmax;
                m.phiy[im] *= dmax;
                m.phiz[im] *= dmax;
            }
            m.moldiam[im] *= dmaxmol;
            for ia in m.naf[im]..=m.nal[im] {
                m.xa[ia] *= dmaxmol;
                m.ya[ia] *= dmaxmol;
                m.za[ia] *= dmaxmol;
                for ib in 1..=m.nbond[ia] {
                    m.tbond[(ia, ib)] *= dmaxmol;
                }
            }
        }
        for ip in 1..=m.npairs {
            m.amassr[ip] /= dmaxmol;
        }
        let piconv: f32 = 180.0 / std::f32::consts::PI;
        for ir in 1..=m.ngrrot {
            if m.rotgrp[ir] == "DIH" {
                for k in 1..=m.ndihedralang[ir] {
                    m.dihedralang[(k, ir)] *= piconv;
                }
            }
        }
        // three-body maximum distances
        for ir in 1..=self.tb.nthreeb {
            for j in 1..=2 {
                self.tb.threebrmax[(j, ir)] *= dmax;
            }
        }
        self.scale_coordconstraint_lengths(false);
        self.sync_lattice();
    }

    /// `subroutine get_maxmolsize`.
    pub fn get_maxmolsize(&mut self) {
        let m = &mut self.mol;
        m.maxmolsize = 0.0;
        for im in 1..=m.nmol {
            if m.moldiam[im] > m.maxmolsize {
                m.maxmolsize = m.moldiam[im];
            }
        }
        m.maxmolsize *= 2.0;
    }

    /// `subroutine centreofmassi(im)`: recentre molecule `im` on its centre of mass.
    pub fn centreofmassi(&mut self, im: i32) {
        let m = &mut self.mol;
        let (mut x1, mut y1, mut z1) = (0f32, 0f32, 0f32);
        let mut summass: f32 = 0.0;
        for ia in m.naf[im]..=m.nal[im] {
            if m.qtype[ia] == 0 {
                let amass = m.amassa[m.itype[ia]];
                x1 += m.xa[ia] * amass;
                y1 += m.ya[ia] * amass;
                z1 += m.za[ia] * amass;
                summass += amass;
            }
        }
        if summass > 0.0 {
            x1 /= summass;
            y1 /= summass;
            z1 /= summass;
        }
        for ia in m.naf[im]..=m.nal[im] {
            m.xa[ia] -= x1;
            m.ya[ia] -= y1;
            m.za[ia] -= z1;
        }
    }

    /// `subroutine moloutc(fname1)`: write the configuration to a `.ato` file.
    pub fn moloutc(&mut self, fname1: &str) -> Result<()> {
        let mut out = OutFile::create(std::path::Path::new(fname1));
        if self.mol.ngrrot > 0 && self.mol.stepri <= 0.0 {
            self.log("\nmoloutc> WARNING! ROT groups defined BUT step size is ZERO!");
        }
        // The Fortran reads natomsphere(itype(ia)) in the atom loop and only
        // allocates natomsphere afterwards: with the array unallocated (a program
        // that never called ctype, e.g. makeato) gfortran crashes with SIGSEGV,
        // leaving an empty file.  Reproduce that outcome (empty file + error).
        if self.mol.nmol > 0 && (self.mol.natomsphere.len() < self.mol.ntype as usize || self.mol.natomsphere.is_empty()) {
            drop(out);
            return crate::error::unsupported(
                "moloutc: natomsphere/pseudofacmol not allocated - the Fortran program crashes here (SIGSEGV) and leaves an empty .ato file",
            );
        }
        let (m, l) = (&mut self.mol, &mut self.lat);
        if m.scalefacmol.len() < m.nmol as usize {
            m.scalefacmol.grow(m.nmol as usize);
            m.scalefacmol.fill(1.0);
        }
        if l.cubic {
            out.line(&fmt("(1x,i7,7(1x,e15.8))", &[FV::I(m.nmol as i64), FV::R(m.dmax), FV::R(m.simtk)]));
        } else {
            out.line(&fmt("(1x,i7,7(1x,e15.8))", &[FV::I(m.nmol as i64), FV::R(m.simtk)]));
            out.line(&fmt("(1x,6(1x,e15.8))", &[FV::R(l.alat), FV::R(l.blat), FV::R(l.clat)]));
            out.line(&fmt("(1x,6(1x,e15.8))", &[FV::R(l.phib), FV::R(l.thetac), FV::R(l.phic)]));
        }
        // NB: in moloutc `cosd`/`sind` are not declared EXTERNAL, so gfortran
        // uses its intrinsics here (exact at 0 and 90 degrees), unlike minitc.
        use crate::fortran::fmath::{cosd_intrinsic, sind_intrinsic};
        l.alatx = l.alat;
        l.blatx = l.blat * cosd_intrinsic(l.phib);
        l.blaty = l.blat * sind_intrinsic(l.phib);
        l.clatx = l.clat * sind_intrinsic(l.thetac) * cosd_intrinsic(l.phic);
        l.claty = l.clat * sind_intrinsic(l.thetac) * sind_intrinsic(l.phic);
        l.clatz = l.clat * cosd_intrinsic(l.thetac);
        out.line(&fmt(
            "(1x,10(1x,e12.5))",
            &[
                FV::R(m.tolout),
                FV::R(m.stepmi),
                FV::R(m.stepri),
                FV::R(m.stepr),
                FV::R(m.stepm),
                FV::R(m.tempb),
                FV::R(m.tempangreq),
                FV::R(m.tempdihreq),
                FV::R(m.tempangwt),
                FV::R(m.tempdihwt),
            ],
        ));
        m.natomr = 0;
        m.natom = 0;
        m.natomcg = 0.0;
        m.natomrcg = 0.0;
        for it in 1..=m.ntype {
            m.cfrac[it] = 0.0;
            m.cfraccg[it] = 0.0;
        }
        let mut ig = 0;
        for im in 1..=m.nmol {
            if m.tetheredatom[im] < 1 || m.tetheredatom[im] > m.nato[im] {
                m.tetheredatom[im] = 0;
            }
            let tetheredtext = if m.tethered[im] {
                fmt("(l1,i5.5)", &[FV::L(true), FV::I(m.tetheredatom[im] as i64)])
            } else {
                "F".to_string()
            };
            let tetheredtext = tetheredtext.as_str();
            out.line(&fmt(
                "(1x,i0,6(1x,e15.8),2x,a6,1x,i0,1x,f10.5,1x,i5)",
                &[
                    FV::I(m.nato[im] as i64),
                    FV::R(m.xadd[im]),
                    FV::R(m.yadd[im]),
                    FV::R(m.zadd[im]),
                    FV::R(m.phix[im]),
                    FV::R(m.phiy[im]),
                    FV::R(m.phiz[im]),
                    FV::Ch(tetheredtext, 6),
                    FV::I(m.nstepreject[im] as i64),
                    FV::R(m.scalefacmol[im]),
                    FV::I(im as i64),
                ],
            ));
            m.natom += m.nato[im];
            let nadd = 1 - m.naf[im];
            for ia in m.naf[im]..=m.nal[im] {
                if m.imageatom[ia] > 0 {
                    out.line(&fmt(
                        "(1x,a3,2(2x,i5))",
                        &[FV::Ch(&m.alabel[ia], 3), FV::I((ia + nadd) as i64), FV::I((m.imageatom[ia] + nadd) as i64)],
                    ));
                } else {
                    out.line(&fmt("(1x,a3,2(2x,i5))", &[FV::Ch(&m.alabel[ia], 3), FV::I((ia + nadd) as i64), FV::I(0)]));
                }
                out.line(&fmt("(8(1x,f14.8))", &[FV::R(m.xa[ia]), FV::R(m.ya[ia]), FV::R(m.za[ia])]));
                let mut v: Vec<FV> = vec![FV::I(m.nbond[ia] as i64)];
                for j in 1..=m.nbond[ia] {
                    v.push(FV::I((m.iabond[(ia, j)] + nadd) as i64));
                    v.push(FV::R(m.tbond[(ia, j)]));
                }
                out.line(&fmt("(1x,i4,(5(1x,i4,1x,e9.3)))", &v));
                let it = m.itype[ia];
                m.natomcg += m.natomsphere[it] * m.pseudofacmol[it];
                if m.qtype[ia] == 0 {
                    m.natomr += 1;
                    m.natomrcg += m.natomsphere[it] * m.pseudofacmol[it];
                }
                m.cfrac[it] += 1.0;
                m.cfraccg[it] += m.natomsphere[it] * m.pseudofacmol[it];
            }
            out.line(&fmt("(1x,14(i4,1x))", &[FV::I(m.ngrot[im] as i64)]));
            if m.ngrot[im] > 0 {
                for _ih in 1..=m.ngrot[im] {
                    ig += 1;
                    out.line(&fmt("(1x,a3,2(2x,i5))", &[FV::Ch(&m.rotgrp[ig], 3)]));
                    out.line(&fmt(
                        "(1x,14(i4,1x))",
                        &[FV::I((m.labaxi[(ig, 1)] + nadd) as i64), FV::I((m.labaxi[(ig, 2)] + nadd) as i64)],
                    ));
                    if m.rotgrp[ig] == "DIH" {
                        out.line(&fmt("(1x,14(i4,1x))", &[FV::I(m.ndihedralang[ig] as i64)]));
                        for k in 1..=m.ndihedralang[ig] {
                            out.line(&fmt(
                                "(1x,2(i4,1x),f10.3)",
                                &[
                                    FV::I((m.labrotdih[(1, k, ig)] + nadd) as i64),
                                    FV::I((m.labrotdih[(2, k, ig)] + nadd) as i64),
                                    FV::R(m.dihedralang[(k, ig)]),
                                ],
                            ));
                        }
                    }
                    let mut v: Vec<FV> = vec![FV::I(m.nrot[ig] as i64)];
                    for j in m.rotf[ig]..=m.rotl[ig] {
                        v.push(FV::I((m.labrot[j] + nadd) as i64));
                    }
                    out.line(&fmt("(1x,14(i4,1x))", &v));
                }
            }
        }
        for it in 1..=m.ntype {
            out.line(&fmt(
                "(1x,a3,1x,a3,1x,i2)",
                &[FV::Ch(&m.type_[it], 3), FV::Ch(&m.atsymbol[it], 3), FV::I(m.iexchange[it] as i64)],
            ));
            out.line(&fmt(
                "(5(1x,f10.6),1x,e13.6,2(1x,f10.4),1x,e13.6)",
                &[
                    FV::R(m.epscor[it]),
                    FV::R(m.sigcor[it]),
                    FV::R(m.amassa[it]),
                    FV::R(m.charge[it]),
                    FV::R(m.qradius[it]),
                    FV::R(m.natomsphere[it]),
                    FV::R(m.radminsphere[it]),
                    FV::R(m.radmaxsphere[it]),
                    FV::R(m.pseudofacmol[it]),
                ],
            ));
        }
        out.line(&fmt("(8(1x,f14.8))", &[FV::R(m.ecore), FV::R(m.dcore)]));
        // saveran1(10)
        {
            let r = &self.rng;
            let mut w = ListWriter::new();
            w.int(r.mcall).int(r.idum).int(r.iy);
            for v in r.iv.iter() {
                w.int(*v);
            }
            out.line(&w.finish());
        }
        let m = &mut self.mol;
        if m.rfractmol.len() < m.nmoltype as usize {
            m.rfractmol.grow(m.nmoltype as usize);
        }
        let mut moltypeswritten: Vec<i32> = Vec::new();
        for im in 1..=m.nmol {
            let mut found = false;
            let mut imoltype = 0;
            while imoltype < m.nmoltype && !found {
                imoltype += 1;
                found = m.alabel[m.naf[im]] == m.moltype[imoltype];
            }
            if !moltypeswritten.contains(&imoltype) {
                moltypeswritten.push(imoltype);
                let tmpl = &m.moltemplate[imoltype];
                let ilen = crate::fortran::text::len_to_blank(tmpl);
                out.line(&fmt(
                    "(1x,i3,1x,a,4(1x,e13.6),1x,f10.6)",
                    &[
                        FV::I(moltypeswritten.len() as i64),
                        FV::S(&tmpl[..ilen]),
                        FV::R(m.interntrans[imoltype]),
                        FV::R(m.internrot[imoltype]),
                        FV::R(m.wholerot[imoltype]),
                        FV::R(m.wholetrans[imoltype]),
                        FV::R(m.rfractmol[imoltype]),
                    ],
                ));
            }
        }
        let l = &self.lat;
        m.volume = l.alatx * l.blaty * l.clatz;
        let rho = m.natomr as f32 / m.volume;
        m.rhocg = m.natomrcg / m.volume;
        if m.natomr > 0 {
            for it in 1..=m.ntype {
                m.cfrac[it] /= m.natomr as f32;
                m.cfraccg[it] /= m.natomrcg;
            }
        }
        out.line(&fmt("(' Atomic number density (coarse grain and actual) = ',2(1x,e14.7))", &[FV::R(m.rhocg), FV::R(rho)]));
        for it in 1..=m.ntype {
            if m.qtypea[it] == 0 {
                out.line(&fmt(
                    "(1x,a3,2x,i1,2x,e12.5,10x,e12.5,2x,a3)",
                    &[FV::Ch(&m.atsymbol[it], 3), FV::I(0), FV::R(m.cfraccg[it]), FV::R(m.cfrac[it]), FV::Ch(&m.type_[it], 3)],
                ));
            }
        }
        out.close()?;
        self.sync_lattice();
        Ok(())
    }
}
