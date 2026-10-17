//! `writexyz` — translation of `fortran_src/writexyz.f` (`run_writexyz`,
//! `writexyz`).
//!
//! STATUS: translated; verified byte-identical against the Fortran, test
//! cases EPSR-test-ex42/43.
//!
//! Invocation (same as the Fortran program):
//! `epsr-rs <folder> writexyz <file>.ato [append|noappend [nmolcentre xmax ymax zmin zmax phi theta chi [excluded types...]]]`
//! Missing values are read from standard input in the same order as the
//! Fortran prompts (append? y/n, centre molecule, the four distances, the
//! three Euler angles, the line of atom types to exclude).
//!
//! Writes `<file>append.xyz` (blocks appended to the existing file), or
//! `<file>NNNNN.xyz` with the block counter kept in `<file>.nblock`.
//! Only cubic boxes are supported (as by the `.ato` reader of the port).

use super::readxyz::groupc_common as gc;
use super::ToolContext;
use crate::ato::{cosd, sind};
use crate::error::{EpsrError, Result};
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::listread::{IoErr, ListRead, TextFile};
use crate::fortran::text::{findex, fstr, parse, substr};
use crate::neighbour::minimage;
#[allow(unused_imports)]
use crate::state::State;

/// Run the `writexyz` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    let (sometext, w) = gc::command_text("writexyz", &ctx.args);
    let mut appendxyz = true;
    let mut appendset = false;
    if w.nwords() > 2 {
        appendxyz = w.word(&sometext, 3) != "noappend";
        appendset = true;
    }
    let (sometextsave, wsave) = if w.nwords() > 3 {
        let s = w.span(&sometext, 4, w.nwords());
        let ws = parse(&s, 5);
        (s, ws)
    } else {
        (String::new(), parse("", 5))
    };
    let nwordssave = wsave.nwords();
    let mut iread: i32 = 1;
    if w.nwords() < 2 {
        return crate::error::unsupported("writexyz without a file name (interactive file search)");
    }
    let mut readfile = w.word(&sometext, 2).to_string();
    let ilendot = findex(&readfile, ".") - 1;
    if ilendot > 0 {
        readfile = readfile[..ilendot as usize].to_string();
    }
    let mut fnameato = gc::get_full_filename(&ctx.dir, &readfile);
    fnameato.push_str(".ato");
    st.mol.ntype = 0;
    st.mol.nmoltype = 0;
    st.mol.sizefactor = 1.0;
    st.minitc(&fnameato)?;
    if st.dist.rho == 0.0 {
        st.log(&fmt("(/'writexyz> Problems with .ato file ',a/'Routine will exit now'/)", &[FV::Ch(&fnameato, 256)]));
        return Ok(());
    }
    for it in 1..=st.mol.ntype {
        let s = st.mol.atsymbol[it].clone();
        if s == "XXX" || s.is_empty() {
            st.log(&fmt("(/'writexyz> Specify the chemical symbol for atom type ',a3,': ',$)", &[FV::Ch(&st.mol.type_[it], 3)]));
            let v = ctx.stdin.list().word()?.unwrap_or_default();
            st.mol.atsymbol[it] = fstr(&v, 3);
        }
    }
    if !appendset {
        st.log("\nwritexyz> Do you want to append to an existing file? (y or n) ");
        let yesorno = ctx.stdin.read_line()?;
        let c = yesorno.chars().next().unwrap_or(' ');
        appendxyz = !(c == 'n' || c == 'N');
    }
    // centre molecule
    let mut nmolcentre: i32 = 0;
    if nwordssave < iread {
        st.log("\nplotato> Specify number of molecule at centre (0 to plot entire box): ");
        ctx.stdin.list().get_i(&mut nmolcentre)?;
    } else {
        ListRead::from_str(wsave.word(&sometextsave, iread)).get_i(&mut nmolcentre)?;
    }
    iread += 1;
    let (mut xaxismax, mut yaxismax, mut zaxismin, mut zaxismax) = (0f32, 0f32, 0f32, 0f32);
    {
        let s;
        let mut r = if nwordssave < iread + 3 {
            st.log("\nplotato> Give maximum distances along x and y axes,\nminimum and maximum distances along the z-axis\n - use zero distances to show only the central molecule: ");
            ctx.stdin.list()
        } else {
            s = wsave.span(&sometextsave, iread, iread + 3);
            ListRead::from_str(&s)
        };
        r.get_r(&mut xaxismax)?;
        r.get_r(&mut yaxismax)?;
        r.get_r(&mut zaxismin)?;
        r.get_r(&mut zaxismax)?;
    }
    iread += 4;
    let (mut zphi, mut ztheta, mut zchi) = (0f32, 0f32, 0f32);
    {
        let s;
        let mut r = if nwordssave < iread + 2 {
            st.log("\nplotato> Give phi, theta and chi for rotation of coordinates (0 0 0) ");
            ctx.stdin.list()
        } else {
            s = wsave.span(&sometextsave, iread, iread + 2);
            ListRead::from_str(&s)
        };
        r.get_r(&mut zphi)?;
        r.get_r(&mut ztheta)?;
        r.get_r(&mut zchi)?;
    }
    iread += 3;
    let czphi = cosd(zphi);
    let szphi = sind(-zphi);
    let cztheta = cosd(ztheta);
    let sztheta = sind(-ztheta);
    let czchi = cosd(zchi);
    let szchi = sind(-zchi);
    for im in 1..=st.mol.nmol {
        gc::centreofmass(&mut st, im);
    }
    let natom = st.mol.natom;
    let n = natom.max(0) as usize + 1;
    let mut iplot = vec![0i32; n];
    let mut xat = vec![0f32; n];
    let mut yat = vec![0f32; n];
    let mut zat = vec![0f32; n];
    let mut iaplot = 0;
    let mut implot = 0;
    let m = &st.mol;
    let (xaddcentre, yaddcentre, zaddcentre) = if nmolcentre > 0 && nmolcentre <= m.nmol {
        (m.xadd[nmolcentre], m.yadd[nmolcentre], m.zadd[nmolcentre])
    } else {
        (0.0, 0.0, 0.0)
    };
    for im in 1..=m.nmol {
        let mut xdif = m.xadd[im] - xaddcentre;
        let mut ydif = m.yadd[im] - yaddcentre;
        let mut zdif = m.zadd[im] - zaddcentre;
        if !st.lat.cubic {
            return crate::error::unsupported("writexyz for non-cubic boxes");
        } else {
            xdif /= m.dmax;
            ydif /= m.dmax;
            zdif /= m.dmax;
            minimage(&mut xdif, &mut ydif, &mut zdif);
            xdif *= m.dmax;
            ydif *= m.dmax;
            zdif *= m.dmax;
        }
        gc::rotatecoords(czphi, szphi, cztheta, sztheta, czchi, szchi, &mut xdif, &mut ydif, &mut zdif);
        let mut ifound = 0;
        if (xdif.abs() <= xaxismax && ydif.abs() <= yaxismax && zdif >= zaxismin && zdif <= zaxismax) || im == nmolcentre {
            ifound = im;
        }
        if ifound > 0 {
            implot += 1;
            for ia in m.naf[im]..=m.nal[im] {
                iaplot += 1;
                iplot[ia as usize] = 1;
            }
        } else {
            for ia in m.naf[im]..=m.nal[im] {
                iplot[ia as usize] = 0;
            }
        }
        for ia in m.naf[im]..=m.nal[im] {
            let k = ia as usize;
            xat[k] = m.xa[ia];
            yat[k] = m.ya[ia];
            zat[k] = m.za[ia];
            gc::rotatecoords(czphi, szphi, cztheta, sztheta, czchi, szchi, &mut xat[k], &mut yat[k], &mut zat[k]);
        }
        for ia in m.naf[im]..=m.nal[im] {
            let k = ia as usize;
            xat[k] += xdif;
            yat[k] += ydif;
            zat[k] += zdif;
        }
    }
    {
        let mut w = ListWriter::new();
        w.int(implot).int(iaplot);
        st.log(&w.finish());
    }
    let mcomp = crate::state::MCOMP as i32;
    let mut itypplot: Vec<i32> = vec![0; mcomp as usize + 1];
    let mut ntypplot = 0;
    for i in 1..=natom {
        let it = m.itype[i];
        if iplot[i as usize] > 0 {
            let mut iflagt = 0;
            let mut j = 0;
            while j < ntypplot && iflagt == 0 {
                j += 1;
                if itypplot[j as usize] == it {
                    iflagt = j;
                }
            }
            if iflagt == 0 && ntypplot < mcomp {
                ntypplot += 1;
                itypplot[ntypplot as usize] = it;
            }
        }
    }
    let exclude_text = if nwordssave < iread {
        let mut v: Vec<FV> = vec![FV::I(ntypplot as i64)];
        for j in 1..=ntypplot {
            v.push(FV::I(itypplot[j as usize] as i64));
            v.push(FV::Ch(&m.type_[itypplot[j as usize]], 3));
        }
        st.log(&fmt("(/'writexyz> ',i3,' components available for plotting:-',/(1x,i3,1x,a3))", &v));
        st.log("\nwritexyz> Give component numbers to EXCLUDE from the plot -  type 0 to include all: ");
        ctx.stdin.read_line()?
    } else {
        substr(&sometextsave, wsave.ncf[(iread - 1) as usize], wsave.ncl[(nwordssave - 1) as usize])
    };
    let wx = parse(&exclude_text, 5);
    let mut nclplt = ntypplot;
    let mut nplt: Vec<i32> = vec![0; mcomp as usize + 2];
    for i in 1..=ntypplot {
        nplt[i as usize] = itypplot[i as usize];
    }
    let mut iw = 0;
    let mut ierr = 0;
    while iw < wx.nwords() && ierr == 0 {
        iw += 1;
        match ListRead::from_str(wx.word(&exclude_text, iw)).int() {
            Ok(Some(itest)) => {
                let mut i = 0;
                while i < nclplt {
                    i += 1;
                    if nplt[i as usize] == itest {
                        let mut j = i - 1;
                        while j < nclplt - 1 {
                            j += 1;
                            nplt[j as usize] = nplt[(j + 1) as usize];
                        }
                        nclplt = j;
                    }
                }
            }
            Ok(None) => {}
            Err(e) => ierr = e.iostat(),
        }
    }
    {
        let mut v: Vec<FV> = vec![FV::I(nclplt as i64)];
        for j in 1..=nclplt {
            v.push(FV::I(nplt[j as usize] as i64));
            v.push(FV::Ch(&m.type_[nplt[j as usize]], 3));
        }
        st.log(&fmt("(/'writexyz> ',i3,' components to be printed:-',/(1x,i3,1x,a3))", &v));
    }
    for ia in 1..=natom {
        let it = m.itype[ia];
        let k = ia as usize;
        if iplot[k] > 0 {
            let mut iflagt = 0;
            let mut j = 0;
            while j < nclplt && iflagt == 0 {
                j += 1;
                if nplt[j as usize] == it {
                    iflagt = j;
                }
            }
            if iflagt == 0 && m.mref[ia] != nmolcentre {
                iplot[k] = 0;
            }
        } else {
            iplot[k] = 0;
        }
    }
    let p = Plot { iplot, xat, yat, zat, nclplt, nmolcentre, appendxyz };
    writexyz(ctx, &st, &readfile, &p)
}

struct Plot {
    iplot: Vec<i32>,
    xat: Vec<f32>,
    yat: Vec<f32>,
    zat: Vec<f32>,
    nclplt: i32,
    nmolcentre: i32,
    appendxyz: bool,
}

/// The atom records of one block (format 113).
fn block_atoms(st: &State, p: &Plot, out: &mut String) {
    let m = &st.mol;
    for ia in 1..=m.natom {
        let k = ia as usize;
        let mut centralmolecule = p.nmolcentre > 0 && p.nmolcentre <= m.nmol;
        if centralmolecule {
            centralmolecule = ia >= m.naf[p.nmolcentre] && ia <= m.nal[p.nmolcentre];
        }
        let mut j = 0;
        let mut iflagt = 0;
        while j < p.nclplt && iflagt == 0 {
            j += 1;
            if p.iplot[k] > 0 || centralmolecule {
                iflagt = j;
                out.push_str(&fmt(
                    "(a3,4x,e10.3,2x,e10.3,2x,e10.3,2x,a3)",
                    &[FV::Ch(&m.atsymbol[m.itype[ia]], 3), FV::R(p.xat[k]), FV::R(p.yat[k]), FV::R(p.zat[k]), FV::Ch(&m.alabel[ia], 3)],
                ));
                out.push('\n');
            }
        }
    }
}

/// `subroutine writexyz(readfile,ilenreadfile)`.
fn writexyz(ctx: &ToolContext, st: &State, readfile: &str, p: &Plot) -> Result<()> {
    let dir = &ctx.dir;
    let natplot = p.iplot.iter().skip(1).filter(|v| **v > 0).count() as i32;
    let header_count = {
        let mut w = ListWriter::new();
        w.int(natplot);
        w.finish()
    };
    if p.appendxyz {
        let fname = format!("{dir}{readfile}append.xyz");
        let path = std::path::Path::new(&fname);
        // open(status='unknown'): read the blocks already present
        let mut f = if path.exists() { TextFile::open(path)? } else { TextFile::from_string("") };
        let mut ierr = 0;
        let mut nblocks = 0;
        let mut ir: i32 = 0;
        let mut natplot_read: Option<i32> = None;
        while ierr == 0 {
            match f.list().int() {
                Ok(Some(v)) => natplot_read = Some(v),
                Ok(None) => {}
                Err(IoErr::Eof) => ierr = -1,
                Err(IoErr::Bad(_)) => {
                    return crate::error::unsupported("writexyz: appending to a file that is not a sequence of .xyz blocks");
                }
            }
            if ierr == 0 {
                ir = 0;
            }
            let np = natplot_read.unwrap_or(0);
            while ierr == 0 && ir < np + 2 {
                ir += 1;
                if f.skip_record().is_err() {
                    ierr = -1;
                }
            }
            if ierr == 0 {
                nblocks += 1;
            }
        }
        // backspace after end of file: positioned at the end of the data.
        // If the last block was incomplete (or nothing was read: natplot is
        // then undefined in the Fortran), rewind and start a new file.
        let mut text = String::new();
        let complete = matches!(natplot_read, Some(np) if ir == np + 2);
        if complete {
            for l in &f.lines {
                text.push_str(l);
                text.push('\n');
            }
        } else {
            nblocks = 0;
        }
        nblocks += 1;
        let blocknumber = fmt("(i5.5)", &[FV::I(nblocks as i64)]);
        text.push_str(&header_count);
        text.push('\n');
        text.push_str(&format!("#Filename: {dir}{readfile}.ato Block {blocknumber}\n"));
        block_atoms(st, p, &mut text);
        text.push('\n');
        std::fs::write(path, text).map_err(|e| EpsrError::Io(format!("{fname}: {e}")))?;
    } else {
        let fname = format!("{dir}{readfile}.nblock");
        let path = std::path::Path::new(&fname);
        let mut nblocks = 0;
        if path.exists() {
            let mut f = TextFile::open(path)?;
            if let Ok(Some(ir)) = f.list().int() {
                nblocks = ir;
            }
        }
        nblocks += 1;
        let mut w = ListWriter::new();
        w.int(nblocks);
        std::fs::write(path, format!("{}\n", w.finish())).map_err(|e| EpsrError::Io(format!("{fname}: {e}")))?;
        let blocknumber = fmt("(i5.5)", &[FV::I(nblocks as i64)]);
        let fname = format!("{dir}{readfile}{blocknumber}.xyz");
        let mut text = String::new();
        text.push_str(&header_count);
        text.push('\n');
        text.push_str(&format!("#Filename: {dir}{readfile}.ato  Block {blocknumber}\n"));
        block_atoms(st, p, &mut text);
        text.push('\n');
        std::fs::write(&fname, text).map_err(|e| EpsrError::Io(format!("{fname}: {e}")))?;
    }
    Ok(())
}
