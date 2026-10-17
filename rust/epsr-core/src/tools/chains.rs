//! `chains` — translation of `fortran_src/chains_sp_criterion.f` (`run_chains`,
//! `find_chains`) and of `chains_routines.f` (`match`, `shortest_path`), the
//! latter shared with `rings` and `lsi`.
//!
//! STATUS: translated, verified byte-identical (EPSR-test-ex62, ex63).
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> chains <input file>`
//! Keyword table: `crate::symtab_tools::CHAINS` (the `nsize` text and default
//! are computed from `mchain` as in `setup_symbols_chains`).
//!
//! Outputs: `<run>.CHAINS.n01` (for each distribution: histogram of the number
//! of bonded neighbours, histogram of chain lengths), the rewritten
//! `<run>.CHAINS.dat`, and `fort.14` (mean and spread of each histogram — the
//! Fortran writes them to unit 14 after `save_files_others` has closed it).

use super::clusters::{list_int, open_input_with, outgofz, rd_int, rd_real};
use super::ToolContext;
use crate::error::Result;
use crate::fortran::format::{fmt, FV};
use crate::fortran::output::OutFile;
use crate::fortran::text::{fstr, parse, substr};
use crate::fortran::{F1, F2};
use crate::neighbour::minimagec;
use crate::state::State;
use crate::symbols::LEN_VALUE;

// chains.inc
const MDIST: usize = 20;
const MZ: usize = 200;
const MCHAIN: i32 = 200;
pub(crate) const MBOND: i32 = 30;
pub(crate) const MRING: i32 = 25;
/// `mmol` of epsr_dimension.inc.
pub(crate) const MMOL: i32 = 50000;

/// `logical function match(string1,string2)` (`chains_routines.f`):
/// `string1` is a `character*3` atom type, `string2` a `character*4` pattern
/// in which `*` matches anything that follows.
pub(crate) fn fmatch(string1: &str, string2: &str) -> bool {
    let s1 = fstr(string1, 3);
    let s2 = fstr(string2, 4);
    match s2.find('*') {
        Some(p) => {
            let nsearch = p as i32; // nstar-1
            if nsearch > 0 {
                substr(&s1, 1, nsearch) == substr(&s2, 1, nsearch)
            } else {
                true
            }
        }
        None => s1 == s2,
    }
}

/// Neighbour tables of `chains.inc` / `lsi.inc` (`neighb(mmol,mbond)`, `nneighb(mmol)`).
pub(crate) struct Neigh {
    pub neighb: F2<i32>,
    pub nneighb: F1<i32>,
}

impl Neigh {
    pub fn new(nmol: i32, mbond: i32) -> Neigh {
        Neigh { neighb: F2::new(nmol as usize + 1, mbond as usize), nneighb: F1::new(nmol as usize + 1) }
    }
}

/// Work arrays of `shortest_path`.
pub(crate) struct SpWork {
    pub rneighb1: F1<i32>,
    pub rneighb2: F1<i32>,
    pub sourceatom1: F1<i32>,
    pub sourceatom2: F1<i32>,
    pub spvaluet: F1<i32>,
    /// `nspvaluet(nmol, nspmax)`
    pub nspvaluet: F2<i32>,
}

impl SpWork {
    pub fn new(nmol: i32, nspmax: i32) -> SpWork {
        let n = nmol as usize + 1;
        SpWork {
            rneighb1: F1::new(n),
            rneighb2: F1::new(n),
            sourceatom1: F1::new(n),
            sourceatom2: F1::new(n),
            spvaluet: F1::new(n),
            nspvaluet: F2::new(nmol.max(1) as usize, nspmax.max(1) as usize),
        }
    }
}

/// `subroutine shortest_path(im1,nz,nmol,nneighb,neighb,rneighb1,rneighb2,spvaluet,nspvaluet,mmol,mbond)`:
/// topological distance of every molecule from `im1` (breadth-first, up to `nz` steps),
/// and in `nspvaluet(im,n)` the number of routes of length `n`.
pub(crate) fn shortest_path(im1: i32, nz: i32, nmol: i32, ng: &Neigh, w: &mut SpWork) {
    for im in 1..=nmol {
        w.spvaluet[im] = nmol + 1;
        for icount in 1..=nz {
            w.nspvaluet[(im, icount)] = 0;
        }
    }
    let mut nneighb1 = 1;
    w.rneighb1[nneighb1] = im1;
    w.sourceatom1[nneighb1] = im1;
    w.spvaluet[im1] = 0;
    let mut icount = 0;
    let mut ineighb1 = 0;
    while icount < nz {
        icount += 1;
        let mut nneighb2 = 0;
        while ineighb1 < nneighb1 {
            ineighb1 += 1;
            let r1 = w.rneighb1[ineighb1];
            let mut inn = 0;
            while inn < ng.nneighb[r1] {
                inn += 1;
                let im2 = ng.neighb[(r1, inn)];
                let im2precede = w.sourceatom1[ineighb1];
                if im2 != im2precede {
                    if icount < w.spvaluet[im2] {
                        w.spvaluet[im2] = icount;
                    }
                    w.nspvaluet[(im2, icount)] += 1;
                    let mut found = false;
                    let mut intest = 0;
                    while !found && intest < nneighb1 {
                        intest += 1;
                        found = im2 == w.rneighb1[intest];
                    }
                    intest = 0;
                    while !found && intest < nneighb2 {
                        intest += 1;
                        found = im2 == w.rneighb2[intest];
                    }
                    if !found && nneighb2 < MMOL {
                        nneighb2 += 1;
                        w.rneighb2[nneighb2] = im2;
                        w.sourceatom2[nneighb2] = r1;
                    }
                }
            }
        }
        let mut inn = 0;
        while inn < nneighb2 && nneighb1 < MMOL {
            inn += 1;
            nneighb1 += 1;
            w.rneighb1[nneighb1] = w.rneighb2[inn];
            w.sourceatom1[nneighb1] = w.sourceatom2[inn];
        }
    }
}

/// Bond criteria of one program run (`atom1..3`, `dff*`, `dfs*`, `dft*`,
/// `ang213min/max` of chains.inc / lsi.inc), indexed `[id][iw]` (1-based).
#[allow(dead_code)]
pub(crate) struct Crit {
    pub md: usize,
    pub nwordsave: Vec<i32>,
    pub atom1: Vec<Vec<String>>,
    pub atom2: Vec<Vec<String>>,
    pub atom3: Vec<Vec<String>>,
    pub dff1: Vec<Vec<f32>>,
    pub dff2: Vec<Vec<f32>>,
    pub dfs1: Vec<Vec<f32>>,
    pub dfs2: Vec<Vec<f32>>,
    pub dft1: Vec<Vec<f32>>,
    pub dft2: Vec<Vec<f32>>,
    pub ang213min: Vec<Vec<f32>>,
    pub ang213max: Vec<Vec<f32>>,
}

/// The per-distribution reading loop of `run_chains` / `run_rings` / `run_lsi`
/// (first 11 data symbols).  Returns `None` where the Fortran writes
/// "Not found enough items ..." and returns.
pub(crate) fn read_criteria(st: &State, ndist: i32, md: usize, prog: &str) -> Result<Option<Crit>> {
    let pi = std::f32::consts::PI;
    let piconv = pi / 180.0;
    let s = || vec![vec![String::new(); md + 1]; md + 1];
    let r = || vec![vec![0f32; md + 1]; md + 1];
    let mut c = Crit {
        md,
        nwordsave: vec![0; md + 1],
        atom1: s(),
        atom2: s(),
        atom3: s(),
        dff1: r(),
        dff2: r(),
        dfs1: r(),
        dfs2: r(),
        dft1: r(),
        dft2: r(),
        ang213min: r(),
        ang213max: r(),
    };
    for id in 1..=ndist {
        let idu = id as usize;
        let mut iref = st.sym.ireff(st.sym.nsymboluptodata + 1, id);
        for k in 0..11 {
            if k > 0 {
                iref += 1;
            }
            let sometext = st.sym.symbol_value[iref].clone();
            let w = parse(&sometext, 5);
            let nwords = w.nwords();
            if k == 0 {
                c.nwordsave[idu] = nwords;
            } else if nwords < c.nwordsave[idu] {
                st.log(&fmt(
                    &format!("('{prog}> Not found enough items in distribution ',i2,', symbol: ',a11)"),
                    &[FV::I(id as i64), FV::Ch(&st.sym.symbol[iref], 11)],
                ));
                return Ok(None);
            }
            if nwords as usize > md {
                return crate::error::unsupported("more bond criteria in one distribution than the Fortran arrays hold");
            }
            for iw in 1..=nwords {
                let word = w.word(&sometext, iw);
                let iwu = iw as usize;
                match k {
                    0 => c.atom1[idu][iwu] = fstr(word, 4),
                    1 => c.atom2[idu][iwu] = fstr(word, 4),
                    2 => c.atom3[idu][iwu] = fstr(word, 4),
                    3 => c.dff1[idu][iwu] = rd_real(word)?,
                    4 => c.dff2[idu][iwu] = rd_real(word)?,
                    5 => c.dfs1[idu][iwu] = rd_real(word)?,
                    6 => c.dfs2[idu][iwu] = rd_real(word)?,
                    7 => c.dft1[idu][iwu] = rd_real(word)?,
                    8 => c.dft2[idu][iwu] = rd_real(word)?,
                    9 => c.ang213max[idu][iwu] = (rd_real(word)? * piconv).cos(),
                    _ => c.ang213min[idu][iwu] = (rd_real(word)? * piconv).cos(),
                }
            }
        }
    }
    Ok(Some(c))
}

/// `atom(1:index(atom,' ')-1)` of a `character*4` variable.
pub(crate) fn upto_blank(a: &str) -> String {
    let p = crate::fortran::format::pad_to(a, 4);
    match p.find(' ') {
        Some(i) => p[..i].to_string(),
        None => String::new(), // index(...)=0: substring (1:-1) is empty
    }
}

/// The bond search of `find_chains` / `find_rings` (identical in both):
/// fills `nneighb`/`neighb` for distribution `id`; returns `nfsum`.
pub(crate) fn find_bonds(st: &State, c: &Crit, id: i32, ng: &mut Neigh) -> i32 {
    let m = &st.mol;
    let nmol = m.nmol;
    let idu = id as usize;
    let dmax = m.dmax;
    let dmaxsqr = dmax * dmax;
    let ty = |ia: i32| -> &str { &m.type_[m.itype[ia]] };
    for i in 1..=nmol {
        ng.nneighb[i] = 0;
    }
    let mut nfsum = 0;
    for im1 in 1..=nmol {
        for im2 in 1..=nmol {
            let mut found = im1 == im2;
            let mut imt = 0;
            while imt < ng.nneighb[im1] && !found {
                imt += 1;
                found = im2 == ng.neighb[(im1, imt)];
            }
            found = found || ng.nneighb[im1] == MBOND || ng.nneighb[im2] == MBOND;
            let mut foundiw = found;
            let mut iw: i32 = 0;
            while !foundiw {
                iw += 1;
                let iwu = iw as usize;
                let mut foundia1 = found;
                let mut ia1 = m.naf[im1] - 1;
                while !foundia1 {
                    ia1 += 1;
                    foundia1 = fmatch(ty(ia1), &c.atom1[idu][iwu]);
                    let mut foundib1 = found;
                    let mut ib1 = m.naf[im1] - 1;
                    while foundia1 && !foundib1 {
                        ib1 += 1;
                        foundib1 = fmatch(ty(ib1), &c.atom2[idu][iwu]);
                        if foundib1 {
                            nfsum += 1;
                        }
                        let xc1 = m.xadd[im1];
                        let yc1 = m.yadd[im1];
                        let zc1 = m.zadd[im1];
                        let xa1 = m.xa[ia1] + xc1;
                        let ya1 = m.ya[ia1] + yc1;
                        let za1 = m.za[ia1] + zc1;
                        let xb1 = m.xa[ib1] + xc1;
                        let yb1 = m.ya[ib1] + yc1;
                        let zb1 = m.za[ib1] + zc1;
                        let mut foundia2 = found;
                        let mut ia2 = m.naf[im2] - 1;
                        while foundib1 && !foundia2 {
                            ia2 += 1;
                            foundia2 = fmatch(ty(ia2), &c.atom3[idu][iwu]);
                            if foundia2 {
                                let xa2 = m.xa[ia2] + m.xadd[im2];
                                let ya2 = m.ya[ia2] + m.yadd[im2];
                                let za2 = m.za[ia2] + m.zadd[im2];
                                let mut xdif = xa2 - xa1;
                                let mut ydif = ya2 - ya1;
                                let mut zdif = za2 - za1;
                                minimagec(&mut xdif, &mut ydif, &mut zdif);
                                let dsq13 = (xdif * xdif + ydif * ydif + zdif * zdif) * dmaxsqr;
                                let da13 = dsq13.sqrt();
                                let foundd13 = da13 >= c.dfs1[idu][iwu] && da13 <= c.dfs2[idu][iwu];
                                if foundd13 {
                                    let mut xdif = xa1 - xb1;
                                    let mut ydif = ya1 - yb1;
                                    let mut zdif = za1 - zb1;
                                    minimagec(&mut xdif, &mut ydif, &mut zdif);
                                    let dsq12 = (xdif * xdif + ydif * ydif + zdif * zdif) * dmaxsqr;
                                    let da12 = dsq12.sqrt();
                                    let foundd12 = da12 >= c.dff1[idu][iwu] && da12 <= c.dff2[idu][iwu];
                                    if foundd12 {
                                        let mut xdif = xa2 - xb1;
                                        let mut ydif = ya2 - yb1;
                                        let mut zdif = za2 - zb1;
                                        minimagec(&mut xdif, &mut ydif, &mut zdif);
                                        let dsq23 = (xdif * xdif + ydif * ydif + zdif * zdif) * dmaxsqr;
                                        let da23 = dsq23.sqrt();
                                        let foundd23 = da23 >= c.dft1[idu][iwu] && da23 <= c.dft2[idu][iwu];
                                        let amin = c.ang213min[idu][iwu];
                                        let amax = c.ang213max[idu][iwu];
                                        if foundd23 && amin <= amax && (amax - amin) < 1.9 {
                                            let anglefound = (dsq12 + dsq13 - dsq23) / (2.0 * da12 * da13);
                                            found = anglefound >= amin && anglefound <= amax;
                                        } else {
                                            found = foundd23;
                                        }
                                    }
                                    if found {
                                        ng.nneighb[im1] += 1;
                                        ng.nneighb[im2] += 1;
                                        let (n1, n2) = (ng.nneighb[im1], ng.nneighb[im2]);
                                        ng.neighb[(im1, n1)] = im2;
                                        ng.neighb[(im2, n2)] = im1;
                                    }
                                }
                            }
                            foundia2 = ia2 >= m.nal[im2] || found;
                        }
                        foundib1 = ib1 >= m.nal[im1] || found;
                    }
                    foundia1 = ia1 >= m.nal[im1] || found;
                }
                foundiw = iw >= c.nwordsave[idu] || found;
            }
        }
    }
    nfsum
}

/// The computed defaults of `setup_symbols_chains` / `_rings` / `_lsi` for `nsize`.
pub(crate) fn set_nsize_default(st: &mut State, text: String, n: i32) {
    let iref = 2;
    st.sym.symbol_text[iref] = fstr(&text, LEN_VALUE);
    let mut w = crate::fortran::format::ListWriter::new();
    w.int(n);
    st.sym.symbol_default[iref] = fstr(&w.finish(), LEN_VALUE);
    st.sym.symbol_value[iref] = st.sym.symbol_default[iref].clone();
}

/// Common driver of `run_chains` and `run_rings`.
pub(crate) fn run_chains_or_rings(ctx: &mut ToolContext, rings: bool) -> Result<()> {
    let mut st = ctx.new_state();
    let (table, nzmax, prog) = if rings {
        (&crate::symtab_tools::RINGS, MRING, "rings")
    } else {
        (&crate::symtab_tools::CHAINS, MCHAIN, "chains")
    };
    let fix = |st: &mut State| {
        if rings {
            set_nsize_default(st, fmt("('Maximum ring length (max ',i2,')')", &[FV::I(MRING as i64)]), MRING);
        } else {
            set_nsize_default(st, fmt("('Maximum chain length (max ',i3,')')", &[FV::I(MCHAIN as i64)]), MCHAIN);
        }
    };
    if !open_input_with(ctx, &mut st, table, &fix)? {
        return Ok(());
    }
    let _ = prog;
    let mut iref: i32 = 2;
    let mut nz = rd_int(&st.sym.symbol_value[iref])?;
    if nz > nzmax {
        nz = nzmax;
    }
    iref += 1;
    let mut nsumt = rd_int(&st.sym.symbol_value[iref])?;
    let nfilesmax = (MDIST / 2) as i32;
    if st.sym.nfiles > nfilesmax {
        st.sym.nfiles = nfilesmax;
    }
    let ndist = st.sym.nfiles;
    let ndist3 = 2 * ndist;
    // the message says 'rings>' in both programs
    let c = match read_criteria(&st, ndist, MDIST, "rings")? {
        Some(c) => c,
        None => return Ok(()),
    };
    let mut z = F1::<f32>::new(MZ);
    for j in 1..=nz {
        z[j] = (j - 1) as f32;
    }
    let mut gofz = F2::<f32>::new(MZ, MDIST);
    let mut gzt = F2::<f32>::new(MZ, MDIST);
    let mut gzt2 = F2::<f32>::new(MZ, MDIST);

    st.normalise_coordinates();
    if rings {
        super::rings::find_rings(&st, &c, nz, ndist, &mut gofz);
    } else {
        find_chains(&st, &c, nz, ndist, &mut gofz);
    }

    let titler = if rings { "Ring_size" } else { "Chain_size" };
    let mut titledat = vec![String::new(); MDIST];
    let mut idref2 = 0usize;
    for id in 1..=ndist as usize {
        let idref1 = idref2 + 1;
        idref2 = idref1 + 1;
        let t = fmt(
            "(4x,a,'-',a,'-',a)",
            &[FV::S(&upto_blank(&c.atom1[id][1])), FV::S(&upto_blank(&c.atom2[id][1])), FV::S(&upto_blank(&c.atom3[id][1]))],
        );
        titledat[idref1 - 1] = fstr(&t, 29);
        titledat[idref2 - 1] = titledat[idref1 - 1].clone();
    }
    let fnameinp = st.sym.fnameinp.clone();
    outgofz(&st, &fnameinp, "n", nz, ndist3, &mut nsumt, &mut z, &gofz, &mut gzt, &mut gzt2, titler, &titledat)?;

    let mut chainlen = vec![[0f32; 2]; MDIST + 1];
    for id in 1..=ndist3 {
        let mut anorm = 0.0f32;
        for i in 1..=nz {
            anorm += gzt[(i, id)];
        }
        if anorm > 0.0 {
            let mut sum1 = 0.0f32;
            let mut sum2 = 0.0f32;
            for i in 1..=nz {
                let term = z[i] * gzt[(i, id)] / anorm;
                sum1 += term;
                sum2 += z[i] * term;
            }
            chainlen[id as usize][0] = sum1;
            let std = sum2 - sum1 * sum1;
            chainlen[id as usize][1] = if std > 0.0 { std.sqrt() } else { 0.0 };
        } else {
            chainlen[id as usize] = [0.0, 0.0];
        }
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
    for id in 1..=ndist as usize {
        let idfirst = 2 * id - 1;
        let idlast = idfirst + 1;
        let mut v = vec![FV::I(id as i64)];
        for idr in idfirst..=idlast {
            v.push(FV::R(chainlen[idr][0]));
            v.push(FV::R(chainlen[idr][1]));
        }
        out.line(&fmt("(1x,i3,6(1x,e12.5))", &v));
    }
    out.close()?;
    Ok(())
}

/// Run the `chains` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    run_chains_or_rings(ctx, false)
}

/// `subroutine find_chains()`
fn find_chains(st: &State, c: &Crit, nz: i32, ndist: i32, gofz: &mut F2<f32>) {
    let nmol = st.mol.nmol;
    let mut ng = Neigh::new(nmol, MBOND);
    let nmu = nmol as usize + 1;
    let mut spvalue = F1::<i32>::new(nmu);
    let mut carryon = F1::<bool>::new(nmu);
    let mut lenchain = F1::<i32>::new(MMOL as usize);
    let mut nrlab = vec![0i32; MZ + 2];
    let mut rneighb = vec![0i32; MZ + 2];
    let mut spseq = vec![0i32; MCHAIN as usize + 2];
    let mut idref2 = 0;
    for id in 1..=ndist {
        let idref1 = idref2 + 1;
        idref2 = idref1 + 1;
        let _nfsum = find_bonds(st, c, id, &mut ng);
        let nspmax = nz;
        let mut w = SpWork::new(nmol, nspmax);
        for im in 1..=nmol {
            carryon[im] = true;
        }
        let mut nchain = 0;
        let mut im1 = 0;
        while im1 < nmol && nchain < MMOL {
            im1 += 1;
            if ng.nneighb[im1] == 1 {
                carryon[im1] = false;
                shortest_path(im1, nspmax, nmol, &ng, &mut w);
                for im in 1..=nmol {
                    spvalue[im] = w.spvaluet[im];
                }
                nrlab[0] = im1;
                nrlab[1] = im1;
                spseq[0] = -1;
                spseq[1] = 0;
                rneighb[1] = 0;
                let mut icount: usize = 1;
                while icount > 0 && nchain < MMOL {
                    if rneighb[icount] >= ng.nneighb[nrlab[icount]] {
                        icount -= 1;
                    } else {
                        rneighb[icount] += 1;
                        let im2 = ng.neighb[(nrlab[icount], rneighb[icount])];
                        if carryon[im2] && im2 != nrlab[icount - 1] && (icount as i32) < nz && nchain < MMOL {
                            let mut ifound = 0;
                            let mut ic = 0;
                            while ifound == 0 && ic < icount {
                                ic += 1;
                                if im2 == nrlab[ic] {
                                    ifound = im2;
                                }
                            }
                            if ifound == 0 {
                                if ng.nneighb[im2] == 1 && nchain < MMOL {
                                    nchain += 1;
                                    icount += 1;
                                    nrlab[icount] = im2;
                                    rneighb[icount] = 1;
                                    spseq[icount] = spvalue[im2];
                                    lenchain[nchain] = icount as i32;
                                } else {
                                    let sptest = spvalue[im2] - spseq[icount];
                                    if sptest == 1 {
                                        icount += 1;
                                        nrlab[icount] = im2;
                                        rneighb[icount] = 0;
                                        spseq[icount] = spvalue[im2];
                                    }
                                }
                            }
                        }
                    }
                }
            } else if ng.nneighb[im1] == 0 && nchain < MMOL {
                nchain += 1;
                lenchain[nchain] = 1;
            }
        }
        for im1 in 1..=nmol {
            let iref = ng.nneighb[im1] + 1;
            if iref <= nz {
                gofz[(iref, idref1)] += 1.0;
            }
        }
        for ir in 1..=nchain {
            let iref = lenchain[ir] + 1;
            if iref <= nz {
                gofz[(iref, idref2)] += 1.0;
            }
        }
    }
}
