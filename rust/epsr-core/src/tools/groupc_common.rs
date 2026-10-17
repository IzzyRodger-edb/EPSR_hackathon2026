//! Routines shared by the file-conversion programs `readxyz`, `writexyz`,
//! `readpdb`, `bonds` and `alignmole` that are not (or not completely) part of
//! the core translation:
//!
//! * `chemical_symbols.inc` (`chsymbol`, `chatwt`, 103 elements);
//! * the command text assembled by `get_command_from_args`, `get_full_filename`
//!   (`file_and_directory_commands.f`);
//! * `centreofmass`, `rotatecoords`, `dihedral`, `get_rotgrp`
//!   (`ato_manipulation_rev4.f`), `convertalphabetagamma` (`auxiliary_routines.f`),
//!   `rotatm` (`molfit.f`);
//! * `moloutc` (`ato_read_and_save.f`) including the parts the core version in
//!   `ato.rs` does not need: non-cubic boxes, `DIH` rotational groups, the
//!   default allocation of `scalefacmol`, `natomsphere`..., `rfractmol`, and
//!   the `.lck` file;
//! * the lattice set-up block that `readxyz` and `readpdb` copy from `minitc`;
//! * Fortran `Fw.d` input editing (used by `readpdb`).

use crate::ato::{cosd, perpd, sind};
use crate::error::Result;
use crate::fortran::fmath::sign;
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::output::OutFile;
use crate::fortran::text::{parse, Words};
use crate::fortran::{F1, F2, F3};
use crate::state::State;

/// `chemical_symbols.inc`: `chsymbol(melements)` (character*3).
pub const CHSYMBOL: [&str; 103] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
];
/// `chemical_symbols.inc`: `chatwt(melements)` (atomic weights, real*4).
pub const CHATWT: [f32; 103] = [
    1.00797, 4.00260, 6.94090, 9.01218, 10.81003, 12.01110, 14.00670, 15.99940, 18.99840, 20.17900, 22.98977, 24.30500, 26.98154,
    28.08550, 30.97376, 32.06600, 35.45300, 39.94800, 39.09830, 40.07800, 44.95591, 47.88000, 50.94150, 51.99600, 54.93800, 55.84700,
    58.93320, 58.67000, 63.54600, 65.39000, 69.72300, 72.59000, 74.92159, 78.96000, 79.90400, 83.80000, 85.46780, 87.62000, 88.90585,
    91.22400, 92.90638, 95.94000, 98.90625, 101.07000, 102.90550, 106.42000, 107.86820, 112.41000, 114.82000, 118.71000, 121.75000,
    127.60000, 126.90447, 131.28999, 132.90543, 137.33000, 138.90550, 140.12000, 140.90765, 144.24001, 145.00000, 150.36000, 151.96001,
    157.25000, 158.92534, 162.50000, 164.93031, 167.25999, 168.93420, 173.03999, 174.96700, 178.49001, 180.94791, 183.95093, 186.20700,
    190.20000, 192.22000, 195.08000, 196.96654, 200.59000, 204.38300, 207.20000, 208.98038, 209.00000, 210.00000, 222.00000, 223.00000,
    226.02541, 227.00000, 232.03806, 231.03587, 238.02890, 237.04817, 244.06274, 243.00000, 247.00000, 247.00000, 251.00000, 254.00000,
    257.00000, 256.00000, 254.00000, 258.00000,
];
/// `melements`.
pub const MELEMENTS: i32 = 103;

/// `chsymbol(ir)` (1-based).
pub fn chsymbol(ir: i32) -> &'static str {
    CHSYMBOL[(ir - 1) as usize]
}
/// `chatwt(ir)` (1-based).
pub fn chatwt(ir: i32) -> f32 {
    CHATWT[(ir - 1) as usize]
}

/// `sometext` as left by `get_command_from_args`: the command name followed
/// by the remaining arguments separated by single blanks, and its `parse`.
pub fn command_text(command: &str, args: &[String]) -> (String, Words) {
    let mut s = command.to_string();
    for a in args {
        let w = parse(a, 5);
        if w.nwords() == 0 {
            // text(ncf(1):ncl(nwords)) of a blank argument: nothing useful; keep a blank word
            s.push(' ');
            continue;
        }
        s.push(' ');
        s.push_str(&w.span(a, 1, w.nwords()));
    }
    let w = parse(&s, 5);
    (s, w)
}

/// `get_full_filename(fnameshort, ...)`: prepend the directory name unless the
/// name already contains a path separator.
pub fn get_full_filename(dir: &str, short: &str) -> String {
    if !dir.is_empty() && !short.contains('/') {
        format!("{dir}{short}")
    } else {
        short.to_string()
    }
}

/// `subroutine convertalphabetagamma(alpha,beta,gam,phib,thetac,phic)`.
/// `thetac`, `phic` are left unchanged on the error branches (intent(out)
/// variables that are not assigned keep their storage in practice).
pub fn convertalphabetagamma(alpha: f32, beta: f32, gam: f32, phib: &mut f32, thetac: &mut f32, phic: &mut f32, log: &dyn Fn(&str)) {
    let pi: f32 = 4.0f32 * 1.0f32.atan();
    let radconv = pi / 180.0;
    let degconv = 180.0 / pi;
    *phib = gam;
    let betar = beta * radconv;
    let cx = betar.cos();
    if cx.abs() < 1.0e-4 {
        *thetac = 0.0;
        *phic = 0.0;
    } else {
        let gamr = gam * radconv;
        let bx = gamr.cos();
        let by = gamr.sin();
        let alphar = alpha * radconv;
        let bdotc = alphar.cos();
        let cy = (bdotc - bx * cx) / by;
        if cx.abs() > 1.0e-4 {
            let phicr = cy.atan2(cx);
            *phic = phicr * degconv;
            let cr = (cx * cx + cy * cy).sqrt();
            if cr.abs() <= 1.0 {
                let thetacr = cr.asin();
                *thetac = degconv * thetacr;
            } else {
                log(&fmt("('convertalphbetagamma> Incompatible axes definitions: ',3(1x,f10.5))", &[FV::R(cr), FV::R(*phib), FV::R(*phic)]));
            }
        } else {
            log(&fmt("('convertalphbetagamma> Incompatible axes definitions: ',3(1x,f10.5))", &[FV::R(alpha), FV::R(beta), FV::R(gam)]));
        }
    }
}

/// The lattice block of `readxyz`/`readpdb` (copied from `minitc`): cell
/// vectors from `alat, blat, clat, phib, thetac, phic`, `dmax`, `dmaxmol`,
/// `rmax`, `simt`, `volume`.
pub fn setup_lattice(st: &mut State) {
    let l = &mut st.lat;
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
    let mut w = ListWriter::new();
    w.real(perpab).real(perpbc).real(perpca);
    let verbose = st.verbose;
    if verbose {
        println!("{}", w.finish());
    }
    let m = &mut st.mol;
    m.dmax = perpab.min(perpbc).min(perpca);
    m.dmaxmol = m.dmax * m.sizefactor;
    let xmax = l.alatx.abs() + l.blatx.abs() + l.clatx.abs();
    let ymax = l.blaty.abs() + l.claty.abs();
    let zmax = l.clatz.abs();
    // `cell` and `rmax` are locals of readxyz/readpdb (only printed)
    let _rmax = 0.5 * (xmax * xmax + ymax * ymax + zmax * zmax).sqrt();
    m.simt = 0.008314 * m.simtk;
    m.volume = l.alatx * l.blaty * l.clatz;
}

/// `subroutine centreofmass(im)`: put the centre of mass of molecule `im` at
/// the molecule centre, moving the centre accordingly.
pub fn centreofmass(st: &mut State, im: i32) {
    let m = &mut st.mol;
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
    m.xadd[im] += x1;
    m.yadd[im] += y1;
    m.zadd[im] += z1;
}

/// `subroutine rotatecoords(czphi,szphi,cztheta,sztheta,czchi,szchi,x,y,z)`.
#[allow(clippy::too_many_arguments)]
pub fn rotatecoords(czphi: f32, szphi: f32, cztheta: f32, sztheta: f32, czchi: f32, szchi: f32, x: &mut f32, y: &mut f32, z: &mut f32) {
    let xdum = *x * czphi + *y * szphi;
    let ydum = *y * czphi - *x * szphi;
    let zdum = *z;
    *x = xdum * cztheta - zdum * sztheta;
    *y = ydum;
    *z = zdum * cztheta + xdum * sztheta;
    let xdum = *x * czchi + *y * szchi;
    let ydum = *y * czchi - *x * szchi;
    let zdum = *z;
    *x = xdum;
    *y = ydum;
    *z = zdum;
}

/// `subroutine rotatm(im,stepsz)` (`molfit.f`): random rotation of molecule
/// `im` about the z, y and x axes.
pub fn rotatm(st: &mut State, im: i32, stepsz: f32) {
    let pi: f32 = std::f32::consts::PI;
    let iref1 = st.mol.naf[im];
    let iref2 = st.mol.nal[im];
    if iref2 > iref1 {
        let r = st.rng.ran1();
        let phi = pi * stepsz * (2.0 * r - 1.0);
        let cphi = phi.cos();
        let sphi = sign((1.0 - cphi * cphi).sqrt(), phi);
        let m = &mut st.mol;
        for ia in iref1..=iref2 {
            let xasave = m.xa[ia] * cphi - m.ya[ia] * sphi;
            m.ya[ia] = m.xa[ia] * sphi + m.ya[ia] * cphi;
            m.xa[ia] = xasave;
        }
        let r = st.rng.ran1();
        let phi = pi * stepsz * (2.0 * r - 1.0);
        let cphi = phi.cos();
        let sphi = sign((1.0 - cphi * cphi).sqrt(), phi);
        let m = &mut st.mol;
        for ia in iref1..=iref2 {
            let zasave = m.za[ia] * cphi - m.xa[ia] * sphi;
            m.xa[ia] = m.za[ia] * sphi + m.xa[ia] * cphi;
            m.za[ia] = zasave;
        }
        let r = st.rng.ran1();
        let phi = pi * stepsz * (2.0 * r - 1.0);
        let cphi = phi.cos();
        let sphi = sign((1.0 - cphi * cphi).sqrt(), phi);
        let m = &mut st.mol;
        for ia in iref1..=iref2 {
            let yasave = m.ya[ia] * cphi - m.za[ia] * sphi;
            m.za[ia] = m.ya[ia] * sphi + m.za[ia] * cphi;
            m.ya[ia] = yasave;
        }
    }
}

/// `scapro` on 3-vectors.
fn scapro(a: &[f32; 3], b: &[f32; 3]) -> f32 {
    crate::ato::scapro(a, b)
}

/// `function dihedral(im1,lb1,lb2,ld1,ld2)`: dihedral angle (degrees) about
/// the axis lb1-lb2 defined by ld1 (bonded to lb1) and ld2 (bonded to lb2);
/// atom numbers relative to the first atom of molecule `im1`.
pub fn dihedral(st: &State, im1: i32, lb1: i32, lb2: i32, ld1: i32, ld2: i32) -> f32 {
    let m = &st.mol;
    let pi: f32 = 4.0f32 * 1.0f32.atan();
    let ia1 = m.naf[im1] + lb1 - 1;
    let ia2 = m.naf[im1] + lb2 - 1;
    let x1 = m.xa[ia1];
    let y1 = m.ya[ia1];
    let z1 = m.za[ia1];
    let x2 = m.xa[ia2];
    let y2 = m.ya[ia2];
    let z2 = m.za[ia2];
    let mut cd = [x2 - x1, y2 - y1, z2 - z1];
    let dmag = scapro(&cd, &cd).sqrt();
    cd[0] /= dmag;
    cd[1] /= dmag;
    cd[2] /= dmag;
    let io1 = m.naf[im1] + ld1 - 1;
    let io2 = m.naf[im1] + ld2 - 1;
    let co1 = [m.xa[io1] - x1, m.ya[io1] - y1, m.za[io1] - z1];
    let co2 = [m.xa[io2] - x2, m.ya[io2] - y2, m.za[io2] - z2];
    let mut cv1 = crate::ato::vecpro(&co1, &cd);
    let mut cv2 = crate::ato::vecpro(&co2, &cd);
    let dmag = scapro(&cv1, &cv1).sqrt();
    cv1[0] /= dmag;
    cv1[1] /= dmag;
    cv1[2] /= dmag;
    let dmag = scapro(&cv2, &cv2).sqrt();
    cv2[0] /= dmag;
    cv2[1] /= dmag;
    cv2[2] /= dmag;
    let mut cosphi = scapro(&cv1, &cv2);
    if cosphi > 1.0 {
        cosphi = 1.0;
    }
    if cosphi < -1.0 {
        cosphi = -1.0;
    }
    let cpar = crate::ato::vecpro(&cv1, &cv2);
    let sinphi = scapro(&cpar, &cd);
    sinphi.atan2(cosphi) * 180.0 / pi
}

/// `subroutine get_rotgrp(ir1,ir2,labrott,nrott)`: atoms that rotate with
/// `ir2`... when rotating about the axis ir1-ir2 (uses the direct bonds
/// `nbondp`).
pub fn get_rotgrp(st: &State, nbondp: &F1<i32>, ir1: i32, ir2: i32) -> Vec<i32> {
    let m = &st.mol;
    let mut labrott: Vec<i32> = Vec::new();
    let mut nrott: i32 = 0;
    let mut iselect = ir1;
    let mut ir2test = false;
    let mut nrotnew = 0;
    let add = |labrott: &mut Vec<i32>, nrott: &mut i32, iaref: i32| {
        let mut itest = 0;
        let mut ifound = 0;
        while itest < *nrott && ifound == 0 {
            itest += 1;
            if iaref == labrott[(itest - 1) as usize] {
                ifound = iaref;
            }
        }
        if ifound == 0 {
            *nrott += 1;
            labrott.push(iaref);
        }
    };
    while nrotnew <= nrott {
        let mut iref = 0;
        while iref < nbondp[iselect] {
            iref += 1;
            let iaref = m.iabond[(iselect, iref)];
            if !ir2test {
                ir2test = iselect != ir1 && iaref == ir2;
            }
            if iaref != ir1 && iaref != ir2 {
                add(&mut labrott, &mut nrott, iaref);
            }
        }
        nrotnew += 1;
        if nrotnew <= nrott {
            iselect = labrott[(nrotnew - 1) as usize];
        }
    }
    if ir2test {
        let mut iref = 0;
        while iref < nbondp[ir2] {
            iref += 1;
            let iaref = m.iabond[(ir2, iref)];
            if iaref != ir1 && iaref != ir2 {
                add(&mut labrott, &mut nrott, iaref);
            }
        }
    }
    let mut w = ListWriter::new();
    w.int(ir1).int(ir2).int(nrott);
    for v in &labrott {
        w.int(*v);
    }
    st.log(&w.finish());
    labrott
}

/// `DIH` rotational-group data (`labrotdih(2,near,matom)`,
/// `dihedralang(near,matom)` in `molecom.inc`), which the core `State` does
/// not hold.
#[derive(Clone, Debug, Default)]
pub struct DihData {
    pub labrotdih: F3<i32>,
    pub dihedralang: F2<f32>,
}

/// Fortran `Fw.d` input editing of the field `s` (exactly `w` characters,
/// blank padded): blanks are ignored (BLANK='NULL'); without a decimal point
/// the last `d` digits are the fraction.  An all-blank field reads as zero.
pub fn read_fwd(s: &str, d: usize) -> crate::fortran::IoRes<f32> {
    use crate::fortran::listread::IoErr;
    let t: String = s.chars().filter(|c| *c != ' ').collect();
    if t.is_empty() {
        return Ok(0.0);
    }
    let b = t.as_bytes();
    let mut i = 0;
    let mut mant = String::new();
    if b[0] == b'+' || b[0] == b'-' {
        mant.push(b[0] as char);
        i = 1;
    }
    let mut digits = String::new();
    let mut seen_dot = false;
    let mut frac_digits = 0usize;
    while i < b.len() && (b[i].is_ascii_digit() || b[i] == b'.') {
        if b[i] == b'.' {
            if seen_dot {
                return Err(IoErr::Bad(format!("bad real {s:?}")));
            }
            seen_dot = true;
        } else {
            digits.push(b[i] as char);
            if seen_dot {
                frac_digits += 1;
            }
        }
        i += 1;
    }
    let mut exp: i32 = 0;
    if i < b.len() {
        let mut j = i;
        if matches!(b[j], b'e' | b'E' | b'd' | b'D' | b'q' | b'Q') {
            j += 1;
        } else if !(b[j] == b'+' || b[j] == b'-') {
            return Err(IoErr::Bad(format!("bad real {s:?}")));
        }
        let es = &t[j..];
        exp = es.parse::<i32>().map_err(|_| IoErr::Bad(format!("bad real {s:?}")))?;
    }
    if digits.is_empty() {
        digits.push('0');
    }
    let scale: i32 = if seen_dot { -(frac_digits as i32) } else { -(d as i32) };
    mant.push_str(&digits);
    let txt = format!("{}e{}", mant, exp + scale);
    txt.parse::<f32>().map_err(|_| IoErr::Bad(format!("bad real {s:?}")))
}

/// `subroutine moloutc(fname1)`: write the configuration to a `.ato` file
/// (complete version: cubic and non-cubic headers, `DIH` groups).
pub fn moloutc(st: &mut State, fname1: &str, dih: Option<&DihData>) -> Result<()> {
    {
        let m = &mut st.mol;
        if m.scalefacmol.len() < m.nmol as usize {
            m.scalefacmol = F1::new(m.nmol.max(0) as usize);
            m.scalefacmol.fill(1.0);
        }
    }
    let mut out = OutFile::create(std::path::Path::new(fname1));
    if st.mol.ngrrot > 0 && st.mol.stepri <= 0.0 {
        st.log("\nmoloutc> WARNING! ROT groups defined BUT step size is ZERO!");
    }
    let (m, l) = (&mut st.mol, &mut st.lat);
    let mut ig = 0;
    if l.cubic {
        out.line(&fmt("(1x,i7,7(1x,e15.8))", &[FV::I(m.nmol as i64), FV::R(m.dmax), FV::R(m.simtk)]));
    } else {
        out.line(&fmt("(1x,i7,7(1x,e15.8))", &[FV::I(m.nmol as i64), FV::R(m.simtk)]));
        out.line(&fmt("(1x,6(1x,e15.8))", &[FV::R(l.alat), FV::R(l.blat), FV::R(l.clat)]));
        out.line(&fmt("(1x,6(1x,e15.8))", &[FV::R(l.phib), FV::R(l.thetac), FV::R(l.phic)]));
    }
    l.alatx = l.alat;
    l.blatx = l.blat * cosd(l.phib);
    l.blaty = l.blat * sind(l.phib);
    l.clatx = l.clat * sind(l.thetac) * cosd(l.phic);
    l.claty = l.clat * sind(l.thetac) * sind(l.phic);
    l.clatz = l.clat * cosd(l.thetac);
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
    // natomsphere... are used inside the molecule loop before the
    // allocation check below; make sure they can be indexed (a fresh
    // allocation is done below exactly as the Fortran does).
    let nt = m.ntype.max(0) as usize;
    let fresh_sphere = m.natomsphere.len() < nt;
    if fresh_sphere {
        // The Fortran indexes the unallocated arrays here; in all uses of
        // this routine by the group-C programs they have been allocated by
        // ctype or minitc, so this branch is only reached with natom = 0.
        m.natomsphere.grow(nt);
        m.pseudofacmol.grow(nt);
    }
    m.natomr = 0;
    m.natom = 0;
    m.natomcg = 0.0;
    m.natomrcg = 0.0;
    for it in 1..=m.ntype {
        m.cfrac[it] = 0.0;
        m.cfraccg[it] = 0.0;
    }
    for im in 1..=m.nmol {
        if m.tetheredatom[im] < 1 || m.tetheredatom[im] > m.nato[im] {
            m.tetheredatom[im] = 0;
        }
        if m.tethered[im] {
            return crate::error::unsupported("writing tethered molecules");
        }
        let tetheredtext = "F";
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
                    let d = match dih {
                        Some(d) => d,
                        None => return crate::error::unsupported("moloutc: DIH group without dihedral data"),
                    };
                    out.line(&fmt("(1x,14(i4,1x))", &[FV::I(m.ndihedralang[ig] as i64)]));
                    for k in 1..=m.ndihedralang[ig] {
                        out.line(&fmt(
                            "(1x,2(i4,1x),f10.3)",
                            &[
                                FV::I((d.labrotdih[(1, k, ig)] + nadd) as i64),
                                FV::I((d.labrotdih[(2, k, ig)] + nadd) as i64),
                                FV::R(d.dihedralang[(k, ig)]),
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
    if fresh_sphere {
        m.natomsphere = F1::new(nt);
        m.rhosphere = F1::new(nt);
        m.radminsphere = F1::new(nt);
        m.radmaxsphere = F1::new(nt);
        m.pseudofacmol = F1::new(nt);
        m.natomsphere.fill(1.0);
        for it in 1..=m.ntype {
            m.radmaxsphere[it] = 0.5 * m.sigcor[it];
        }
        m.pseudofacmol.fill(1.0);
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
    {
        let r = &st.rng;
        let mut w = ListWriter::new();
        w.int(r.mcall).int(r.idum).int(r.iy);
        for v in r.iv.iter() {
            w.int(*v);
        }
        out.line(&w.finish());
    }
    let m = &mut st.mol;
    if m.rfractmol.len() < m.nmoltype as usize {
        m.rfractmol = F1::new(m.nmoltype as usize);
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
    let l = &st.lat;
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
    // Open the lock file if it exists and set it to 'F' (unlocked)
    let lck = format!("{}.lck", fname1.trim_end());
    if std::path::Path::new(&lck).exists() {
        std::fs::write(&lck, "F\n")?;
    }
    Ok(())
}
