//! `clusters` — translation of `fortran_src/clusters.f` (`run_clusters`,
//! `find_clusters`), plus the histogram helpers of `read_and_save_files.f`
//! (`outgofz`) shared with `chains`, `rings` and `lsi`.
//!
//! STATUS: translated, verified byte-identical (EPSR-test-ex60, ex61).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> clusters <input file>`
//! Keyword table: `crate::symtab_tools::CLUSTERS`.
//!
//! Outputs: `<run>.CLUSTERS.n01` (cluster size histograms; row 1 holds the
//! fraction of molecules in percolating clusters), the rewritten
//! `<run>.CLUSTERS.dat`, and `fort.14` (the Fortran writes the "fraction
//! containing 2 or more molecules" lines to unit 14 after
//! `save_files_others` has closed it, so they land in `fort.14`).

use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::output::OutFile;
use crate::fortran::text::indexlastdot;
use crate::fortran::{F1, F2};
use crate::neighbour::minimagec;
use crate::state::State;

/// `mdist` of `clusters.inc` (= `mfile`).
const MDIST: usize = 30;
/// `mz` of `clusters.inc`.
const MZ: i32 = 5000;

/// List-directed read of the first item of a symbol value (`read(v,*) x`).
pub(crate) fn rd_int(s: &str) -> Result<i32> {
    let mut r = crate::fortran::listread::ListRead::from_str(s);
    r.int()?.ok_or_else(|| EpsrError::Input(format!("cannot read an integer from {s:?}")))
}
pub(crate) fn rd_real(s: &str) -> Result<f32> {
    let mut r = crate::fortran::listread::ListRead::from_str(s);
    r.real()?.ok_or_else(|| EpsrError::Input(format!("cannot read a real from {s:?}")))
}

/// `write(symbol_value(iref),*) n`
pub(crate) fn list_int(n: i32) -> String {
    let mut w = ListWriter::new();
    w.int(n);
    w.finish().trim_end().to_string()
}

/// `subroutine outgofz(fnameinp,ext,nz,ndist,nsum,z,gofz,gzt,gzt2,titler,titledat,...)`:
/// accumulate the histograms `gofz` with those of previous runs (read from
/// `<inp>.<ext>01` when `nsum > 0`) and write the running average and spread.
#[allow(clippy::too_many_arguments)]
pub(crate) fn outgofz(
    st: &State,
    fnameinp: &str,
    ext: &str,
    nz: i32,
    ndist: i32,
    nsum: &mut i32,
    z: &mut F1<f32>,
    gofz: &F2<f32>,
    gzt: &mut F2<f32>,
    gzt2: &mut F2<f32>,
    titler: &str,
    titledat: &[String],
) -> Result<()> {
    let ilen = indexlastdot(fnameinp).max(0) as usize;
    let fname = format!("{}{}", &fnameinp[..ilen], ext);
    let log = st.logger();
    let mut ierr = 0;
    let lptin;
    if *nsum > 0 {
        match crate::output::r_int(&fname, ndist, z, gzt, gzt2, &log)? {
            Some(l) => lptin = l,
            None => {
                ierr = 1;
                lptin = 0;
            }
        }
    } else {
        lptin = nz;
    }
    if ierr != 0 {
        *nsum = -1;
    }
    if *nsum > 0 {
        for id in 1..=ndist {
            for j in 1..=nz {
                if gzt2[(j, id)] < 1.0e-20 {
                    gzt2[(j, id)] = 0.0;
                }
                if gzt[(j, id)].abs() < 1.0e-20 {
                    gzt[(j, id)] = 0.0;
                }
                if j <= lptin {
                    gzt2[(j, id)] = gzt2[(j, id)] * gzt2[(j, id)] + gzt[(j, id)] * gzt[(j, id)];
                } else {
                    gzt[(j, id)] = 0.0;
                    gzt2[(j, id)] = 0.0;
                }
            }
        }
    } else {
        for id in 1..=ndist {
            for j in 1..=nz {
                gzt[(j, id)] = 0.0;
                gzt2[(j, id)] = 0.0;
            }
        }
    }
    let fc1: f64;
    let fc2: f64;
    if *nsum >= 0 {
        *nsum += 1;
        fc1 = (1.0f32 / *nsum as f32) as f64;
        fc2 = 1.0 - fc1;
    } else {
        fc1 = 1.0;
        fc2 = 0.0;
    }
    let lpt = nz;
    for id in 1..=ndist {
        for i in 1..=nz {
            let a2 = fc2 * gzt[(i, id)] as f64;
            let a1 = fc1 * gofz[(i, id)] as f64;
            let a = a1 + a2;
            let b2 = fc2 * gzt2[(i, id)] as f64;
            let g = gofz[(i, id)];
            let b1 = fc1 * (g * g) as f64;
            let b = b1 + b2;
            let mut error = b - a * a;
            if error > 0.0 {
                error = error.sqrt();
            } else {
                error = 0.0;
            }
            gzt[(i, id)] = a as f32;
            gzt2[(i, id)] = error as f32;
        }
    }
    crate::output::w_int(lpt, &fname, ndist, z, gzt, gzt2, titler, titledat, &log)
}

/// Install a tool's keyword table and read its input file (= `setup_symbols_<prog>`
/// + `setup_input_file` + `read_input_file`), after `fix` has applied the
/// computed defaults of the setup routine.  `Ok(false)`: no input file
/// (the Fortran returns at once).
pub(crate) fn open_input_with(ctx: &ToolContext, st: &mut State, table: &crate::symbols::ToolTable, fix: &dyn Fn(&mut State)) -> Result<bool> {
    st.install_tool_table(table);
    fix(st);
    let name = ctx.input_name().ok_or_else(|| EpsrError::Input(format!("{}: no input file given", table.name)))?.to_string();
    match st.setup_input_file(&name) {
        Ok(mut f) => {
            let ierrout = st.read_input_file(&mut f)?;
            Ok(ierrout == 0)
        }
        Err(EpsrError::Io(_)) => Ok(false),
        Err(e) => Err(e),
    }
}

/// `type(it).eq.symbol_value(iref)` search used to find the atom types.
fn find_type(st: &State, v: &str) -> i32 {
    let v = v.trim_end();
    let mut it = 0;
    while it < st.mol.ntype {
        it += 1;
        if st.mol.type_[it].trim_end() == v {
            return it;
        }
    }
    0
}

/// Run the `clusters` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    if !open_input_with(ctx, &mut st, &crate::symtab_tools::CLUSTERS, &|_| {})? {
        return Ok(());
    }
    let mut iref: i32 = 2;
    let mut nz = rd_int(&st.sym.symbol_value[iref])?;
    if nz > MZ {
        nz = MZ;
    }
    iref += 1;
    let mut nsumt = rd_int(&st.sym.symbol_value[iref])?;
    let mut ndist = st.sym.nfiles;
    if ndist > MDIST as i32 {
        ndist = MDIST as i32;
    }
    let nd = MDIST + 1;
    let mut ltyp = vec![0i32; nd];
    let mut ltyp1 = vec![0i32; nd];
    let mut atom1 = vec![String::new(); nd];
    let mut atom2 = vec![String::new(); nd];
    let mut damin = vec![0f32; nd];
    let mut damax = vec![0f32; nd];
    for id in 1..=ndist {
        let iref1 = st.sym.ireff(st.sym.nsymboluptodata + 1, id);
        let (iref2, iref3, iref4) = (iref1 + 1, iref1 + 2, iref1 + 3);
        let ifound1 = find_type(&st, &st.sym.symbol_value[iref1]);
        let ifound2 = find_type(&st, &st.sym.symbol_value[iref2]);
        if ifound1 == 0 || ifound2 == 0 {
            return crate::error::unsupported("clusters: atom type not found in the .ato file (Fortran would index type(0))");
        }
        let idu = id as usize;
        ltyp[idu] = ifound1;
        ltyp1[idu] = ifound2;
        atom1[idu] = st.mol.type_[ifound1].trim_end().to_string();
        atom2[idu] = st.mol.type_[ifound2].trim_end().to_string();
        damin[idu] = rd_real(&st.sym.symbol_value[iref3])?;
        damax[idu] = rd_real(&st.sym.symbol_value[iref4])?;
    }
    let mzu = MZ as usize;
    let mut z = F1::<f32>::new(mzu);
    for j in 1..=nz {
        z[j] = (j - 1) as f32;
    }
    let mut gofz = F2::<f32>::new(mzu, MDIST);
    let mut gzt = F2::<f32>::new(mzu, MDIST);
    let mut gzt2 = F2::<f32>::new(mzu, MDIST);

    st.normalise_coordinates();
    find_clusters(&st, nz, ndist, &ltyp, &ltyp1, &damin, &damax, &mut gofz);

    let titler = "Cluster_size";
    let mut titledat = vec![String::new(); MDIST];
    for id in 1..=ndist as usize {
        titledat[id - 1] = format!("   {}-{}", atom1[id], atom2[id]);
    }
    let fnameinp = st.sym.fnameinp.clone();
    outgofz(&st, &fnameinp, "n", nz, ndist, &mut nsumt, &mut z, &gofz, &mut gzt, &mut gzt2, titler, &titledat)?;

    let mut fract2 = vec![0f32; nd];
    for id in 1..=ndist {
        let mut sum1 = 0.0f32;
        let mut sum2 = 0.0f32;
        for i in 1..=nz {
            sum1 += z[i] * gzt[(i, id)];
            if z[i] > 1.0 {
                sum2 += z[i] * gzt[(i, id)];
            }
        }
        fract2[id as usize] = sum2 / sum1;
    }

    let icount = st.sym.search_for_symbol("nsize", 0);
    let mut iref = st.sym.ireff(icount, 1);
    st.sym.symbol_value[iref] = list_int(nz);
    iref += 1;
    st.sym.symbol_value[iref] = list_int(nsumt);
    iref += 1;
    st.sym.symbol_value[iref] = list_int(ndist);
    st.save_inp_file()?;

    // unit 14 was closed by save_files_others: these records go to fort.14
    let mut out = OutFile::create(&ctx.path("fort.14"));
    for id in 1..=ndist {
        out.line(&fmt(
            "(/'Cluster ',i3,': fraction containing 2 or more molecules',f10.5)",
            &[FV::I(id as i64), FV::R(fract2[id as usize])],
        ));
    }
    out.close()?;
    Ok(())
}

/// `subroutine find_clusters()`
#[allow(clippy::too_many_arguments)]
fn find_clusters(st: &State, nz: i32, ndist: i32, ltyp: &[i32], ltyp1: &[i32], damin: &[f32], damax: &[f32], gofz: &mut F2<f32>) {
    let m = &st.mol;
    let nmol = m.nmol;
    let nmu = nmol as usize + 2;
    let mut neighb = F2::<i32>::new(nmu, MDIST);
    let mut mbond = F1::<i32>::new(nmu);
    let mut labcluster = F1::<i32>::new(nmu);
    let mut lencluster = F1::<i32>::new(nmu);
    let mut ncllab = F1::<i32>::new(nmu);
    let dmax = m.dmax;
    for id in 1..=ndist {
        let idu = id as usize;
        for i in 1..=nmol {
            mbond[i] = -1;
        }
        for im1 in 1..=nmol - 1 {
            for ia1 in m.naf[im1]..=m.nal[im1] {
                for im2 in im1 + 1..=nmol {
                    for ia2 in m.naf[im2]..=m.nal[im2] {
                        if (m.itype[ia1] == ltyp[idu] && m.itype[ia2] == ltyp1[idu]) || (m.itype[ia2] == ltyp[idu] && m.itype[ia1] == ltyp1[idu]) {
                            if mbond[im1] < 0 {
                                mbond[im1] = 0;
                            }
                            if mbond[im2] < 0 {
                                mbond[im2] = 0;
                            }
                            let xa1 = m.xa[ia1] + m.xadd[im1];
                            let ya1 = m.ya[ia1] + m.yadd[im1];
                            let za1 = m.za[ia1] + m.zadd[im1];
                            let xa2 = m.xa[ia2] + m.xadd[im2];
                            let ya2 = m.ya[ia2] + m.yadd[im2];
                            let za2 = m.za[ia2] + m.zadd[im2];
                            let mut xdif = xa2 - xa1;
                            let mut ydif = ya2 - ya1;
                            let mut zdif = za2 - za1;
                            minimagec(&mut xdif, &mut ydif, &mut zdif);
                            let dsq12 = xdif * xdif + ydif * ydif + zdif * zdif;
                            let da12 = dsq12.sqrt() * dmax;
                            if da12 >= damin[idu] && da12 <= damax[idu] {
                                mbond[im1] += 1;
                                mbond[im2] += 1;
                                if mbond[im1] as usize > MDIST || mbond[im2] as usize > MDIST {
                                    panic!("clusters: more than {MDIST} bonds for one molecule (overflows neighb in the Fortran)");
                                }
                                neighb[(im1, mbond[im1])] = im2;
                                neighb[(im2, mbond[im2])] = im1;
                            }
                        }
                    }
                }
            }
        }
        let mut ncluster = 0;
        for im1 in 1..=nmol {
            if mbond[im1] >= 0 {
                ncluster += 1;
                ncllab[im1] = ncluster;
                labcluster[ncluster] = ncluster;
            } else {
                ncllab[im1] = nmol + 1;
            }
        }
        let mut nclusterold = 0;
        while nclusterold != ncluster {
            nclusterold = ncluster;
            for im1 in 1..=nmol {
                if mbond[im1] >= 0 {
                    let mut nseek = ncllab[im1];
                    let mut inn = 0;
                    while inn < mbond[im1] {
                        inn += 1;
                        let im2 = neighb[(im1, inn)];
                        if ncllab[im2] < nseek {
                            nseek = ncllab[im2];
                        }
                    }
                    ncllab[im1] = nseek;
                    inn = 0;
                    while inn < mbond[im1] {
                        inn += 1;
                        let im2 = neighb[(im1, inn)];
                        ncllab[im2] = nseek;
                    }
                }
            }
            ncluster = 0;
            for im1 in 1..=nmol {
                if mbond[im1] >= 0 {
                    let mut iflag = 0;
                    let mut ic = 0;
                    while ic < ncluster && iflag == 0 {
                        ic += 1;
                        if ncllab[im1] == labcluster[ic] {
                            iflag = ic;
                        }
                    }
                    if iflag == 0 {
                        ncluster += 1;
                        labcluster[ncluster] = ncllab[im1];
                    }
                }
            }
        }
        for ic in 1..=ncluster {
            let mut sum = 0.0f32;
            for im in 1..=nmol {
                if ncllab[im] == labcluster[ic] {
                    sum += 1.0;
                }
            }
            lencluster[ic] = sum as i32;
        }
        let mut lensum = 0;
        for ic in 1..=ncluster {
            lensum += lencluster[ic];
        }
        let avlen = if ncluster > 0 { lensum as f32 / ncluster as f32 } else { 0.0 };
        st.log(&fmt(
            "(/'For cluster ',i2,' there are ',i5,' clusters of average size = ',e12.5,' molecules')",
            &[FV::I(id as i64), FV::I(ncluster as i64), FV::R(avlen)],
        ));
        let mut npercolate = 0;
        for ic in 1..=ncluster {
            let mut im1 = 0;
            let mut iflag = 0;
            while im1 < nmol && iflag == 0 {
                im1 += 1;
                if ncllab[im1] == labcluster[ic] {
                    iflag = im1;
                }
            }
            let (mut xmin, mut ymin, mut zmin, mut xmax, mut ymax, mut zmax) = (0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32, 0.0f32);
            for im2 in 1..=nmol {
                if ncllab[im2] == labcluster[ic] {
                    let mut xdif = m.xadd[im2] - m.xadd[im1];
                    let mut ydif = m.yadd[im2] - m.yadd[im1];
                    let mut zdif = m.zadd[im2] - m.zadd[im1];
                    minimagec(&mut xdif, &mut ydif, &mut zdif);
                    if xdif > xmax {
                        xmax = xdif;
                    }
                    if xdif < xmin {
                        xmin = xdif;
                    }
                    if ydif > ymax {
                        ymax = ydif;
                    }
                    if ydif < ymin {
                        ymin = ydif;
                    }
                    if zdif > zmax {
                        zmax = zdif;
                    }
                    if zdif < zmin {
                        zmin = zdif;
                    }
                }
            }
            let xdif = xmax - xmin;
            let ydif = ymax - ymin;
            let zdif = zmax - zmin;
            let diftest = 1.0 - damax[idu] / dmax;
            if xdif >= diftest || ydif >= diftest || zdif >= diftest {
                npercolate += lencluster[ic];
            }
        }
        let fracpercolate = npercolate as f32 / lensum as f32;
        st.log(&fmt(
            "(/'Clusters> For cluster ',i2,' fraction of molecules which percolate is ',f10.5)",
            &[FV::I(id as i64), FV::R(fracpercolate)],
        ));
        for ic in 1..=ncluster {
            let iref = lencluster[ic] + 1;
            if iref <= nz {
                gofz[(iref, id)] += 1.0;
            } else {
                gofz[(nz, id)] += 1.0;
            }
        }
        for ic in 2..=nz {
            gofz[(ic, id)] /= ncluster as f32;
        }
        gofz[(1, id)] = fracpercolate;
    }
}
