//! `bonds` — translation of `fortran_src/bondsclus.f` (`bondsclus`).
//!
//! STATUS: translated; verified byte-identical against the Fortran, test case
//! EPSR-test-ex46.
//!
//! Invocation (same as the Fortran program): `epsr-rs <folder> bonds <file>.ato`
//!
//! For each run of consecutive molecules of the same type (identified by the
//! label of the first atom) writes to `<file>.bonds` the matrix of average
//! intramolecular atom–atom distances (upper triangle) and their standard
//! deviations (lower triangle).

use super::readxyz::groupc_common as gc;
use super::ToolContext;
use crate::error::Result;
use crate::fortran::format::{fmt, ListWriter, FV};
use crate::fortran::output::OutFile;
use crate::fortran::text::findex;
use crate::fortran::F2;
#[allow(unused_imports)]
use crate::state::State;

/// Run the `bonds` program.
pub fn run(ctx: &mut ToolContext) -> Result<()> {
    let mut st = ctx.new_state();
    let (sometext, w) = gc::command_text("bonds", &ctx.args);
    if w.nwords() < 2 {
        return crate::error::unsupported("bonds without a file name (interactive file search)");
    }
    let mut readfile = w.word(&sometext, 2).to_string();
    let ilendot = findex(&readfile, ".") - 1;
    if ilendot > 0 {
        readfile = readfile[..ilendot as usize].to_string();
    }
    let mut fnameato = gc::get_full_filename(&ctx.dir, &readfile);
    fnameato.push_str(".ato");
    let fnamebonds = format!("{}{}.bonds", ctx.dir, readfile);
    st.mol.ntype = 0;
    st.mol.nmoltype = 0;
    st.mol.sizefactor = 1.0;
    st.minitc(&fnameato)?;
    if st.dist.rho == 0.0 {
        st.log(&fmt("(/'Problems with .ato file ',a/'Routine will exit now'/)", &[FV::Ch(&fnameato, 256)]));
        return Ok(());
    }
    let mut out = OutFile::create(std::path::Path::new(&fnamebonds));
    let m = &st.mol;
    let mut ibcnt = 0;
    let mut chtype = String::new();
    let mut im = 0;
    let mut natommol: i32 = 0;
    let mut nsum: i32 = 0;
    let mut alabelt: Vec<String> = Vec::new();
    let mut bond: F2<f32> = F2::new(0, 0);
    let write_block = |out: &mut OutFile, natommol: i32, nsum: i32, alabelt: &[String], bond: &mut F2<f32>| {
        let mut v: Vec<FV> = Vec::new();
        for i in 1..=natommol {
            v.push(FV::I(i as i64));
            v.push(FV::Ch(&alabelt[(i - 1) as usize], 3));
        }
        for l in fmt("(//7x,40(1x,i5,a3,2x))", &v).split('\n') {
            out.line(l);
        }
        out.line("");
        for i in 1..=natommol - 1 {
            for j in i + 1..=natommol {
                // the final block uses float(nsum), the others real(nsum): identical
                bond[(i, j)] /= nsum as f32;
                bond[(j, i)] = bond[(j, i)] / nsum as f32 - bond[(i, j)] * bond[(i, j)];
                if bond[(j, i)] > 0.0 {
                    bond[(j, i)] = bond[(j, i)].sqrt();
                } else {
                    bond[(j, i)] = 0.0;
                }
            }
        }
        for i in 1..=natommol {
            let mut v: Vec<FV> = vec![FV::I(i as i64), FV::Ch(&alabelt[(i - 1) as usize], 3)];
            for j in 1..=natommol {
                v.push(FV::R(bond[(i, j)]));
            }
            for l in fmt("(i3,a3,1x,(40(1x,e10.3)))", &v).split('\n') {
                out.line(l);
            }
        }
        out.line("");
    };
    while im < m.nmol {
        im += 1;
        if m.alabel[m.naf[im]] != chtype {
            if ibcnt > 0 {
                let mut wl = ListWriter::new();
                wl.int(ibcnt).int(natommol).int(nsum);
                st.log(&wl.finish());
                write_block(&mut out, natommol, nsum, &alabelt, &mut bond);
            }
            ibcnt += 1;
            chtype = m.alabel[m.naf[im]].clone();
            natommol = m.nato[im];
            alabelt = (1..=natommol).map(|i| m.alabel[i + m.naf[im] - 1].clone()).collect();
            // allocate(bond(natommol,natommol)): fresh (uninitialised) storage; every
            // element that is written out is assigned below.
            bond = F2::new(natommol as usize, natommol as usize);
            nsum = 1;
            for i in 1..=natommol - 1 {
                bond[(i, i)] = 0.0;
                let iref = i + m.naf[im] - 1;
                let x1 = m.xa[iref];
                let y1 = m.ya[iref];
                let z1 = m.za[iref];
                for j in i + 1..=natommol {
                    bond[(j, j)] = 0.0;
                    let jref = j + m.naf[im] - 1;
                    let xd = m.xa[jref] - x1;
                    let yd = m.ya[jref] - y1;
                    let zd = m.za[jref] - z1;
                    let dist = xd * xd + yd * yd + zd * zd;
                    bond[(i, j)] = dist.sqrt();
                    bond[(j, i)] = dist;
                }
            }
        } else {
            nsum += 1;
            for i in 1..=natommol - 1 {
                let iref = i + m.naf[im] - 1;
                let x1 = m.xa[iref];
                let y1 = m.ya[iref];
                let z1 = m.za[iref];
                for j in i + 1..=natommol {
                    let jref = j + m.naf[im] - 1;
                    let xd = m.xa[jref] - x1;
                    let yd = m.ya[jref] - y1;
                    let zd = m.za[jref] - z1;
                    let dist = xd * xd + yd * yd + zd * zd;
                    bond[(i, j)] += dist.sqrt();
                    bond[(j, i)] += dist;
                }
            }
        }
    }
    let mut wl = ListWriter::new();
    wl.int(ibcnt).int(natommol).int(nsum);
    st.log(&wl.finish());
    write_block(&mut out, natommol, nsum, &alabelt, &mut bond);
    out.close()?;
    Ok(())
}
