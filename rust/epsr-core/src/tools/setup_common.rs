//! Routines shared by the box-building programs (`makeato`, `mixato`, `fmole`,
//! `addato`, `randomise`, `makelattice`, `makelatticeato`).
//!
//! Translation of:
//! * the argument handling of `get_command_from_args` (`file_and_directory_commands.f`):
//!   the command and its arguments are joined into `sometext` and parsed;
//! * `get_system_commands` (`system_commands.f`, only the entries used here:
//!   `system_ls`, `system_quotes`, `system_termination`);
//! * `search_for_files`, `get_full_filename` (`file_and_directory_commands.f`);
//! * `rotatm` (`molfit.f`).
//!
//! Like the Fortran, `system_commands.txt` is looked for in the *current working
//! directory*, and `search_for_files` runs the configured directory-listing
//! command through the shell (`call system(...)`), which creates `dirlist.txt`
//! in the current working directory.

use super::ToolContext;
use crate::error::{EpsrError, Result};
use crate::fortran::fmath::sign;
use crate::fortran::listread::{IoErr, TextFile};
use crate::fortran::text::{concatenate, findex, parse, parse_n, substr, Words};
use crate::state::State;

/// The system commands (`system_commands.inc`) needed by these programs.
#[derive(Clone, Debug, Default)]
pub struct SysCmds {
    pub ls: String,
    pub quotes: String,
    pub termination: String,
}

impl SysCmds {
    /// `subroutine get_system_commands`: read `system_commands.txt` from the
    /// current directory; missing entries stay blank (length 0).
    pub fn load() -> SysCmds {
        let mut s = SysCmds::default();
        if let Ok(mut f) = TextFile::open(std::path::Path::new("system_commands.txt")) {
            while let Ok(line) = f.read_line() {
                let w = parse_n(&line, 11, 256);
                if w.nwords() > 1 {
                    let key = w.word(&line, 1);
                    let val = || concatenate(2, w.nwords(), &w, &line, 256).0;
                    match key {
                        "system_ls" => s.ls = val(),
                        "system_termination" => s.termination = val(),
                        "system_quotes" => s.quotes = val(),
                        _ => {}
                    }
                }
            }
        }
        s
    }
}

/// `sometext` and its words as left by `get_command_from_args`: the command
/// name followed by the remaining arguments, separated by single blanks.
pub fn command_text(command: &str, ctx: &ToolContext) -> (String, Words) {
    let mut s = command.to_string();
    for a in &ctx.args {
        let w = parse(a, 5);
        if w.nwords() > 0 {
            s.push(' ');
            s.push_str(&w.span(a, 1, w.nwords()));
        } else {
            s.push(' ');
        }
    }
    let w = parse(&s, 5);
    (s, w)
}

/// Fortran `index(s, sub)` including the zero-length case (`index(s,'')` = 1).
fn findex0(s: &str, sub: &str) -> i32 {
    if sub.is_empty() {
        1
    } else {
        findex(s, sub)
    }
}

/// `read(5,'(a)') sometext` from standard input (a Fortran runtime error on EOF).
pub fn read_text(ctx: &mut ToolContext) -> Result<String> {
    ctx.stdin
        .read_line()
        .map_err(|_| EpsrError::Input("unexpected end of standard input".into()))
}

/// Map an end-of-file on standard input to an error (the Fortran would stop
/// with a runtime error, or loop for ever on an `iostat=` read).
pub fn eof_err(e: IoErr) -> EpsrError {
    match e {
        IoErr::Eof => EpsrError::Input("unexpected end of standard input".into()),
        IoErr::Bad(s) => EpsrError::Input(s),
    }
}

/// `subroutine get_full_filename(fnameshort,ilenshort,fnamelong,ilenlong)`:
/// prepend the directory unless the name already contains the directory
/// termination string.  (With an empty `system_termination`, `index` returns 1
/// and the name is used as given.)
pub fn get_full_filename(sys: &SysCmds, dir: &str, short: &str) -> String {
    let ilentermination = findex0(short, &sys.termination);
    if !dir.is_empty() && ilentermination == 0 {
        format!("{dir}{short}")
    } else {
        short.to_string()
    }
}

/// `subroutine search_for_files(ext,ilenext,iaccept)`: list the files with the
/// extension (via the shell command `system_ls`), offer each in turn, and if
/// none is accepted ask for a name.  Returns (`sometext`, `iaccept`).
pub fn search_for_files(st: &State, ctx: &mut ToolContext, sys: &SysCmds, ext: &str) -> Result<(String, char)> {
    let cmd = format!("{} {}{}*{}{} > dirlist.txt", sys.ls, sys.quotes, ctx.dir, ext, sys.quotes);
    let ok = std::process::Command::new("sh")
        .arg("-c")
        .arg(&cmd)
        .status()
        .map(|s| s.success())
        .unwrap_or(false);
    let mut ierr = if ok { 0 } else { 1 };
    let mut dl = TextFile::default();
    if ierr == 0 {
        match TextFile::open(std::path::Path::new("dirlist.txt")) {
            Ok(f) => dl = f,
            Err(_) => ierr = 1,
        }
    }
    let mut sometext = String::new();
    let mut iaccept = ' ';
    let mut ifound = 0;
    let mut nline: i32 = 0;
    while ierr == 0 && iaccept != 'y' && iaccept != 'e' {
        sometext.clear();
        match dl.read_line() {
            Ok(l) => sometext = l,
            Err(_) => ierr = -1,
        }
        nline += 1;
        let w = parse(&sometext, 5);
        if w.nwords() > 0 {
            ifound += 1;
            let (fname, ilenfname) = concatenate(1, w.nwords(), &w, &sometext, 256);
            let ilenfname = ilenfname;
            // strip everything up to the last directory termination
            let mut idirpath = 1;
            let mut ifirstchar: i32 = 1;
            let lt = sys.termination.len() as i32;
            while idirpath > 0 && ifirstchar <= ilenfname {
                idirpath = findex0(&substr(&fname, ifirstchar, ilenfname), &sys.termination);
                if idirpath > 0 {
                    let next = ifirstchar + idirpath + lt - 1;
                    if next == ifirstchar {
                        return Err(EpsrError::Input(
                            "search_for_files: system_termination is blank - the Fortran loops for ever here".into(),
                        ));
                    }
                    ifirstchar = next;
                }
            }
            let fname = substr(&fname, ifirstchar, ilenfname);
            st.log(&format!("\nFilename: {}? (Type y to accept, u to go back, e to exit) ", fname));
            sometext = fname;
            let a = read_text(ctx)?;
            iaccept = a.chars().next().unwrap_or(' ');
            if iaccept == 'u' {
                dl.backspace();
                dl.backspace();
                nline -= 2;
                if nline < 0 {
                    dl.pos = 0;
                    nline = 0;
                    let mut e = 0;
                    let mut nw = 0;
                    while e == 0 {
                        let t = match dl.read_line() {
                            Ok(l) => l,
                            Err(_) => {
                                e = -1;
                                String::new()
                            }
                        };
                        nline += 1;
                        nw = parse(&t, 5).nwords();
                    }
                    // backspace over the end-of-file condition (no record moved)
                    nline -= 1;
                    if nw <= 0 {
                        dl.backspace();
                        nline -= 1;
                    }
                    ierr = 0;
                }
            }
        }
        if ifound > 0 && ierr != 0 {
            dl.pos = 0;
            ifound = 0;
            ierr = 0;
            nline = 0;
        }
    }
    if ifound == 0 {
        st.log(&format!(
            "\n No files of extension \"{}\" found in directory \"{}\"\n",
            ext, ctx.dir
        ));
        iaccept = 'n';
    }
    if iaccept != 'y' {
        st.log("\n No files selected...\n\n Type the required filename (without extension): ");
        let t = ctx.stdin.read_line();
        let t = match t {
            Ok(l) => l,
            Err(_) => String::new(),
        };
        let w = parse(&t, 5);
        if w.nwords() > 0 {
            iaccept = 'y';
            sometext = concatenate(1, w.nwords(), &w, &t, 256).0;
        } else {
            iaccept = 'n';
            sometext = t;
        }
    }
    Ok((sometext, iaccept))
}

/// `subroutine rotatm(im, stepsz)` (`molfit.f`): random rotation of molecule
/// `im` about the x, y and z axes.
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

/// Strip an extension: `ilendot=index(readfile,'.')-1` (first dot).
pub fn strip_first_dot(readfile: &str) -> String {
    let ilendot = findex(readfile, ".") - 1;
    if ilendot > 0 {
        readfile[..ilendot as usize].to_string()
    } else {
        readfile.to_string()
    }
}

/// Strip an extension at the last dot (`indexlastdot`).
pub fn strip_last_dot(readfile: &str) -> String {
    let ilendot = crate::fortran::text::indexlastdot(readfile) - 1;
    if ilendot > 0 {
        readfile[..ilendot as usize].to_string()
    } else {
        readfile.to_string()
    }
}

/// `read(5,*) x` (list-directed real from standard input).
pub fn read_real5(ctx: &mut ToolContext) -> Result<f32> {
    let mut v = 0f32;
    ctx.stdin.list().get_r(&mut v).map_err(eof_err)?;
    Ok(v)
}

/// `read(5,*) n` (list-directed integer from standard input).
pub fn read_int5(ctx: &mut ToolContext) -> Result<i32> {
    let mut v = 0i32;
    ctx.stdin.list().get_i(&mut v).map_err(eof_err)?;
    Ok(v)
}
