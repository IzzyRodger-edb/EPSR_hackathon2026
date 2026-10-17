//! Coarse graining — the parts of `epsr` that are active when the first value
//! of the `qinscalefc` keyword (`qscaleinputfactor`) is > 1
//! (`coarsegraining = qscaleinputfactor > 1.0`, module `molecom_mod`).
//!
//! Translation of:
//! * the `if(coarsegraining)` block at the end of `inita`
//!   (`epsr_standalone_rev3.f`): per-molecule scale factors `scalefacmoluse`,
//!   form-factor bins (`nsffbintype`, `assignatomscalefac`), the spherical
//!   form factors `sphereformfactor`, the intramolecular form factors
//!   `intraformfactor`/`intraformfactorself`, and the files `formfactors.txt`,
//!   `intraformfactors.txt`, `intraformfactorsself.txt`;
//! * `function gaussianqr2(q,r)`;
//! * the `wavinter` histograms of `atogen` (`ato_manipulation_rev4.f`).
//!
//! The other effects of coarse graining are in the shared routines:
//! the pair scale factor `scalefac = 0.5*(scalefacmol(im)+scalefacmol(jm))`
//! in `en_mole_pair` / `en_change` (`molfit.f`) and in the `potex` count of
//! `atogen`.
//!
//! NB: in the Fortran the form factors and `wavinter` are only *computed*:
//! the code in `gensq` that would use them is commented out (gensq only fills
//! the unused temporary `prodssff`), so apart from the three `.txt` files they
//! do not influence any result.  They are translated anyway so that the
//! Rust state mirrors the Fortran one.

use crate::error::Result;
use crate::fortran::format::ListWriter;
use crate::fortran::fmath::{pow3, PI4ATAN};
use crate::fortran::output::OutFile;
use crate::fortran::{F1, F2, F3};
use crate::state::State;

/// Coarse-graining arrays of `molecom_mod` that only the coarse-graining code uses.
#[derive(Clone, Debug, Default)]
pub struct CoarseGrain {
    pub scalefacmoluse: F1<f32>,
    pub scalefactypemin: F1<f32>,
    pub scalefactypemax: F1<f32>,
    pub nsffbintype: F1<i32>,
    pub assignatomscalefac: F1<i32>,
    pub maxnscalefactype: i32,
    /// `sphereformfactor(nq, maxnscalefactype, ntype)`
    pub sphereformfactor: F3<f32>,
    /// `intraformfactor(nq, npairs)`
    pub intraformfactor: F2<f32>,
    /// `intraformfactorself(nq, ntype)`
    pub intraformfactorself: F2<f32>,
    pub ntermssff: F2<i32>,
    pub fracntermssff: F2<f32>,
    pub ntermsintra: F1<i32>,
    pub ntermsintraself: F1<i32>,
    /// `wavinter(npt1, npairs, maxnscalefactype, maxnscalefactype)`, column-major.
    pub wavinter: Vec<f32>,
    pub wavinter_dims: [usize; 4],
}

impl CoarseGrain {
    /// Offset of `wavinter(iz, ic, ib, jb)`.
    #[inline(always)]
    pub fn wi(&self, iz: i32, ic: i32, ib: i32, jb: i32) -> usize {
        let d = &self.wavinter_dims;
        let (iz, ic, ib, jb) = ((iz - 1) as usize, (ic - 1) as usize, (ib - 1) as usize, (jb - 1) as usize);
        assert!(iz < d[0] && ic < d[1] && ib < d[2] && jb < d[3], "wavinter index out of range");
        iz + d[0] * (ic + d[1] * (ib + d[2] * jb))
    }
}

/// `function gaussianqr2(q,r)`: form factor of a Gaussian density across a sphere.
pub fn gaussianqr2(q: f32, r: f32) -> f32 {
    if r <= 0.0 {
        return 1.0;
    }
    let qr = q * r;
    let mut arg = 0.51042225f32 * qr;
    arg = 0.5 * arg * arg;
    if arg < 30.0 {
        (-arg).exp()
    } else {
        0.0
    }
}

impl State {
    /// The `if(coarsegraining) then ... end if` block of `inita`.
    pub fn inita_coarsegrain(&mut self, nq: i32, q: &F1<f32>) -> Result<()> {
        let pi = PI4ATAN;
        let m = &self.mol;
        let nmol = m.nmol;
        let natom = m.natom;
        let ntype = m.ntype;
        let npairs = m.npairs;
        let cg = &mut self.cg;
        // Calculate modified scale factors based on number of atoms in each molecule
        cg.scalefacmoluse.grow(nmol.max(1) as usize);
        for im in 1..=nmol {
            let mut totalatoms: f32 = 0.0;
            for ia in m.naf[im]..=m.nal[im] {
                totalatoms += m.natomsphere[m.itype[ia]] * m.pseudofacmol[m.itype[ia]];
            }
            cg.scalefacmoluse[im] = totalatoms.powf(m.selfbroadpwr);
        }
        let qfactor = 4.0 * pi / 3.0;
        // minimum and maximum scale factors for each atom type
        cg.scalefactypemin = F1::new(ntype as usize);
        cg.scalefactypemax = F1::new(ntype as usize);
        cg.nsffbintype = F1::new(ntype as usize);
        cg.scalefactypemin.fill(1.0);
        cg.scalefactypemax.fill(1.0);
        for ia in 1..=natom {
            let im = m.mref[ia];
            let it = m.itype[ia];
            cg.scalefactypemin[it] = m.scalefacmol[im].min(cg.scalefactypemin[it]);
            cg.scalefactypemax[it] = m.scalefacmol[im].max(cg.scalefactypemax[it]);
        }
        let scalefacratio: f32 = 0.2;
        cg.maxnscalefactype = 0;
        for it in 1..=ntype {
            cg.nsffbintype[it] = ((cg.scalefactypemax[it] - cg.scalefactypemin[it]) / scalefacratio) as i32 + 1;
            if cg.nsffbintype[it] > cg.maxnscalefactype {
                cg.maxnscalefactype = cg.nsffbintype[it];
            }
        }
        cg.assignatomscalefac = F1::new(natom.max(1) as usize);
        for ia in 1..=natom {
            let im = m.mref[ia];
            let it = m.itype[ia];
            let ibin = ((m.scalefacmol[im] - cg.scalefactypemin[it]) / scalefacratio) as i32 + 1;
            cg.assignatomscalefac[ia] = ibin;
        }
        let maxn = cg.maxnscalefactype as usize;
        let nqu = nq as usize;
        cg.sphereformfactor = F3::new(nqu, maxn, ntype as usize);
        cg.intraformfactor = F2::new(nqu, npairs as usize);
        cg.intraformfactorself = F2::new(nqu, ntype as usize);
        cg.ntermssff = F2::new(maxn, ntype as usize);
        cg.fracntermssff = F2::new(maxn, ntype as usize);
        cg.ntermsintra = F1::new(npairs as usize);
        cg.ntermsintraself = F1::new(ntype as usize);
        for ia in 1..=natom {
            let im = m.mref[ia];
            let it = m.itype[ia];
            let qfactorit = qfactor * m.rhosphere[it];
            let radlocalmin = m.radminsphere[it];
            let radlocalmax = m.radmaxsphere[it];
            let ibin = cg.assignatomscalefac[ia];
            cg.ntermssff[(ibin, it)] += 1;
            let scalefac = m.scalefacmol[im];
            let selffac = cg.scalefacmoluse[im];
            let mut radusemin = radlocalmin * scalefac;
            let mut radusemax = radlocalmax * scalefac;
            let radusemin3 = qfactorit * pow3(radusemin);
            let radusemax3 = qfactorit * pow3(radusemax);
            radusemin = radlocalmin * selffac;
            radusemax = radlocalmax * selffac;
            for i in 1..=nq {
                let qq = q[i];
                cg.sphereformfactor[(i, ibin, it)] +=
                    radusemax3 * gaussianqr2(qq, radusemax) - radusemin3 * gaussianqr2(qq, radusemin);
            }
        }
        let mut sumatoms: f32 = 0.0;
        for it in 1..=ntype {
            let mut nterms: i32 = 0;
            let mut sumatomsbin: f32 = 0.0;
            for ib in 1..=cg.nsffbintype[it] {
                nterms += cg.ntermssff[(ib, it)];
                sumatomsbin += cg.sphereformfactor[(1, ib, it)];
                if cg.ntermssff[(ib, it)] > 1 {
                    let d = cg.ntermssff[(ib, it)] as f32;
                    for i in 1..=nq {
                        cg.sphereformfactor[(i, ib, it)] /= d;
                    }
                }
            }
            if nterms > 0 {
                sumatoms += m.cfrac[it] * sumatomsbin / nterms as f32;
                for ib in 1..=cg.nsffbintype[it] {
                    cg.fracntermssff[(ib, it)] = cg.ntermssff[(ib, it)] as f32 / nterms as f32;
                }
            } else {
                for ib in 1..=cg.nsffbintype[it] {
                    cg.fracntermssff[(ib, it)] = 0.0;
                }
            }
        }
        sumatoms = sumatoms * m.natomr as f32 / nmol as f32;
        let msg = format!(" inita> Average number of atoms per molecule: {sumatoms}");
        if sumatoms > 0.0 {
            let sfactor = (1.0 / sumatoms).sqrt();
            for v in cg.sphereformfactor.data.iter_mut() {
                *v *= sfactor;
            }
        }
        {
            let mut out = OutFile::create(&self.sym.full("formfactors.txt"));
            for i in 1..=nq {
                let mut w = ListWriter::new();
                w.real(q[i]);
                for it in 1..=ntype {
                    w.real(cg.sphereformfactor[(i, 1, it)]);
                }
                out.line(&w.finish());
            }
            out.close()?;
        }
        // Now do the intra-molecular terms
        for im in 1..=nmol {
            let scalefac = m.scalefacmol[im];
            let selffac = m.scalefacmol[im];
            for ia in m.naf[im]..=m.nal[im] {
                let it = m.itype[ia];
                let qfactorit = qfactor * m.rhosphere[it];
                let mut raditmin = m.radminsphere[it] * scalefac;
                let mut raditmax = m.radmaxsphere[it] * scalefac;
                let raditmin3 = qfactorit * pow3(raditmin);
                let raditmax3 = qfactorit * pow3(raditmax);
                raditmin *= selffac;
                raditmax *= selffac;
                for ja in ia..=m.nal[im] {
                    let jt = m.itype[ja];
                    let ip = m.npairref[(it, jt)];
                    cg.ntermsintra[ip] += 1;
                    if it == jt {
                        let nitterms = m.natomsphere[it] * m.pseudofacmol[it];
                        if nitterms > 0.0 {
                            let niitermstotal = nitterms * (nitterms - 1.0);
                            let niitermsself = m.pseudofacmol[it] * m.natomsphere[it] * (m.natomsphere[it] - 1.0);
                            let selfcorrself = niitermsself / (nitterms * nitterms);
                            let selfcorr = (niitermstotal - niitermsself) / (nitterms * nitterms);
                            cg.ntermsintraself[it] += 1;
                            for i in 1..=nq {
                                let qq = q[i];
                                let formfactorit = raditmax3 * gaussianqr2(qq, raditmax) - raditmin3 * gaussianqr2(qq, raditmin);
                                cg.intraformfactorself[(i, it)] += selfcorrself * formfactorit * formfactorit;
                                cg.intraformfactor[(i, ip)] += selfcorr * formfactorit * formfactorit;
                            }
                        }
                    } else {
                        let qfactorjt = qfactor * m.rhosphere[jt];
                        let mut radjtmin = m.radminsphere[jt] * scalefac;
                        let mut radjtmax = m.radmaxsphere[jt] * scalefac;
                        let radjtmin3 = qfactorjt * pow3(radjtmin);
                        let radjtmax3 = qfactorjt * pow3(radjtmax);
                        radjtmin *= selffac;
                        radjtmax *= selffac;
                        for i in 1..=nq {
                            let qq = q[i];
                            let formfactorit = raditmax3 * gaussianqr2(qq, raditmax) - raditmin3 * gaussianqr2(qq, raditmin);
                            let formfactorjt = radjtmax3 * gaussianqr2(qq, radjtmax) - radjtmin3 * gaussianqr2(qq, radjtmin);
                            cg.intraformfactor[(i, ip)] += formfactorit * formfactorjt;
                        }
                    }
                }
            }
        }
        // Normalise form factors to number of atoms
        for ip in 1..=npairs {
            if cg.ntermsintra[ip] > 0 && sumatoms > 0.0 {
                let sfactor = 1.0 / (cg.ntermsintra[ip] as f32 * sumatoms);
                for i in 1..=nq {
                    cg.intraformfactor[(i, ip)] *= sfactor;
                }
            }
        }
        for it in 1..=ntype {
            if cg.ntermsintraself[it] > 0 && sumatoms > 0.0 {
                let sfactor = 1.0 / (cg.ntermsintraself[it] as f32 * sumatoms);
                for i in 1..=nq {
                    cg.intraformfactorself[(i, it)] = sfactor * cg.intraformfactorself[(i, it)];
                }
            }
        }
        {
            let mut out = OutFile::create(&self.sym.full("intraformfactors.txt"));
            for i in 1..=nq {
                let mut w = ListWriter::new();
                w.real(q[i]);
                for ip in 1..=npairs {
                    w.real(cg.intraformfactor[(i, ip)]);
                }
                out.line(&w.finish());
            }
            out.close()?;
            let mut out = OutFile::create(&self.sym.full("intraformfactorsself.txt"));
            for i in 1..=nq {
                let mut w = ListWriter::new();
                w.real(q[i]);
                for it in 1..=ntype {
                    w.real(cg.intraformfactorself[(i, it)]);
                }
                out.line(&w.finish());
            }
            out.close()?;
        }
        self.log(&msg);
        Ok(())
    }

    /// `atogen`: `if(allocated(wavinter)) deallocate(wavinter); allocate(...); wavinter=0.0`.
    pub fn atogen_cg_alloc(&mut self) {
        let n = self.cg.maxnscalefactype.max(0) as usize;
        let d = [self.dist.npt1.max(0) as usize, self.mol.npairs.max(0) as usize, n, n];
        self.cg.wavinter_dims = d;
        self.cg.wavinter = vec![0.0; d[0] * d[1] * d[2] * d[3]];
    }
}

impl CoarseGrain {
    /// `atogen`, "Repeat for coarse graining": add the pairs of atom `ia`
    /// (type `it`) recorded in `iz1values`/`icrefvalues` to `wavinter`.
    #[allow(clippy::too_many_arguments)]
    pub fn atogen_accumulate(&mut self, itype: &F1<i32>, ia: i32, it: i32, ja1: i32, ja2: i32, iz1values: &F1<i32>, icrefvalues: &F1<i32>, addwav: f32) {
        let ibinia = self.assignatomscalefac[ia];
        for ja in ja1..=ja2 {
            let iz1 = iz1values[ja];
            if iz1 > 0 {
                let jt = itype[ja];
                let ibinja = self.assignatomscalefac[ja];
                let icref = icrefvalues[ja];
                // Pairs are referred in the order it <= jt
                let k = if it <= jt { self.wi(iz1, icref, ibinia, ibinja) } else { self.wi(iz1, icref, ibinja, ibinia) };
                self.wavinter[k] += addwav;
            }
        }
    }

    /// `atogen` normalisation of `wavinter` for pair `ic` (types `it`,`jt`):
    /// divide by `cprod2` when > 0, else set to `ycoh/(nbins_it*nbins_jt)`.
    pub fn atogen_normalise(&mut self, ycoh: &F1<f32>, npt1: i32, ic: i32, it: i32, jt: i32, cprod2: f32) {
        let nbi = self.nsffbintype[it];
        let nbj = self.nsffbintype[jt];
        if cprod2 > 0.0 {
            for iz in 1..=npt1 {
                for ibinia in 1..=nbi {
                    for ibinja in 1..=nbj {
                        let k = self.wi(iz, ic, ibinia, ibinja);
                        self.wavinter[k] /= cprod2;
                    }
                }
            }
        } else {
            let ycohfactor = 1.0 / (nbi * nbj) as f32;
            for iz in 1..=npt1 {
                let ycohmod = ycoh[iz] * ycohfactor;
                for ibinia in 1..=nbi {
                    for ibinja in 1..=nbj {
                        let k = self.wi(iz, ic, ibinia, ibinja);
                        self.wavinter[k] = ycohmod;
                    }
                }
            }
        }
    }
}
