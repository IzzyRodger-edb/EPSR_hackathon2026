//! Run a validation case through the `epsr-core` library API instead of the
//! `epsr-rs` executable: the same calls a Python (PyO3) wrapper makes.
//!
//! usage: epsr-api-driver CASE_DIR WORK_DIR
//!
//! The input files of CASE_DIR are copied to WORK_DIR and the steps of
//! CASE_DIR/steps.txt (see ../run_validation.sh) are executed there:
//!
//! * `run epsr <inp>`        -> `Epsr::setup(dir, inp, false)?.run()?`
//! * `run <tool> args [< f]` -> `tools::run_tool(tool, &mut ToolContext::new(dir, args, stdin, false))`
//! * `set <file> <kw> <val>` -> rewrite one keyword line (as validation/set_keyword.py)
//! * `delete <files>`        -> remove files
//!
//! The process changes its current directory to WORK_DIR (required by the
//! interactive setup programs, as for the Fortran).
//!
//! After every `epsr` run the convergence summary (`sim.summary`) is printed as
//! one tab-separated line, and at the end the fit to each dataset is written to
//! WORK_DIR/api_fit_<n>.txt from `q_grid`, `data_fq`, `simulated_fq`.
use epsr_core::fortran::TextFile;
use epsr_core::tools::{run_tool, ToolContext};
use epsr_core::Epsr;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

type Res<T> = Result<T, Box<dyn std::error::Error>>;

fn set_keyword(path: &Path, kw: &str, value: &str) -> Res<()> {
    let text = fs::read(path)?;
    let text: String = text.iter().map(|&b| b as char).collect(); // latin-1
    let mut lines: Vec<String> = text.split('\n').map(str::to_string).collect();
    let sep = " ".repeat(15);
    for line in lines.iter_mut() {
        if line.trim() == "q" {
            break;
        }
        if line.get(..12.min(line.len())).map(str::trim) == Some(kw) {
            let t = line.trim_end();
            let desc = if t.contains(&sep) { t.rsplit(&sep).next().unwrap_or("") } else { "" };
            *line = format!("{kw:<12}{value}{sep}{desc}");
            let out: Vec<u8> = lines.join("\n").chars().map(|c| c as u8).collect();
            fs::write(path, out)?;
            return Ok(());
        }
    }
    Err(format!("keyword {kw} not found in {}", path.display()).into())
}

fn main() -> Res<()> {
    let a: Vec<String> = std::env::args().skip(1).collect();
    if a.len() != 2 {
        eprintln!("usage: epsr-api-driver CASE_DIR WORK_DIR");
        std::process::exit(2);
    }
    let case = fs::canonicalize(&a[0])?;
    let _ = fs::remove_dir_all(&a[1]);
    fs::create_dir_all(&a[1])?;
    let work = fs::canonicalize(&a[1])?;
    for e in fs::read_dir(&case)? {
        let p = e?.path();
        let name = p.file_name().unwrap().to_string_lossy().to_string();
        if p.is_file() && name != "steps.txt" && name != "README.md" {
            fs::copy(&p, work.join(&name))?;
        }
    }
    // Like the Fortran, the interactive setup programs (mixato, fmole, ...) read
    // system_commands.txt and write dirlist.txt in the *current directory*, and
    // look up typed file names there: run with the folder as the current
    // directory and pass "./" (file names are built as folder + name).
    std::env::set_current_dir(&work)?;
    let dir = PathBuf::from("./");
    let mut last: Option<Epsr> = None;
    println!("# step\titer\tenergd\tfoldsv\tpressd\tenergabs\trfactot\trfacd...");
    for (n, raw) in fs::read_to_string(case.join("steps.txt"))?.lines().enumerate() {
        let line = raw.split('#').next().unwrap().trim();
        if line.is_empty() {
            continue;
        }
        let (kind, rest) = line.split_once(' ').unwrap_or((line, ""));
        eprintln!("step {}: {line}", n + 1);
        match kind {
            "run" => {
                let (cmd_part, stdin_file) = match rest.split_once('<') {
                    Some((c, f)) => (c.trim(), Some(f.trim())),
                    None => (rest.trim(), None),
                };
                let words: Vec<String> = cmd_part.split_whitespace().map(str::to_string).collect();
                let (cmd, args) = (words[0].as_str(), words[1..].to_vec());
                if cmd == "epsr" {
                    let mut sim = Epsr::setup(&dir, &args[0], false)?;
                    sim.run()?;
                    let s = &sim.summary;
                    let rd: Vec<String> = s.rfacd.iter().map(|x| format!("{x:.5e}")).collect();
                    println!(
                        "{}\t{}\t{:.5}\t{:.5}\t{:.5}\t{:.5}\t{:.5e}\t{}",
                        n + 1, s.iter, s.energd, s.foldsv, s.pressd, s.energabs, s.rfactot, rd.join("\t")
                    );
                    last = Some(sim);
                } else {
                    let stdin = match stdin_file {
                        Some(f) => TextFile::open(&case.join(f))?,
                        None => TextFile::from_string(""),
                    };
                    let mut ctx = ToolContext::new(&dir, args, stdin, false);
                    run_tool(cmd, &mut ctx)?;
                }
            }
            "set" => {
                let w: Vec<&str> = rest.split_whitespace().collect();
                set_keyword(&work.join(w[0]), w[1], &w[2..].join(" "))?;
            }
            "delete" => {
                for f in rest.split_whitespace() {
                    let _ = fs::remove_file(work.join(f));
                }
            }
            _ => return Err(format!("unknown step {line}").into()),
        }
    }
    // Fit of the last run, from the in-memory arrays (what a wrapper hands to numpy).
    if let Some(sim) = last {
        let q = sim.q_grid();
        for id in 0..sim.n_datasets() {
            let (d, f) = (sim.data_fq(id), sim.simulated_fq(id));
            let (lo, hi) = sim.data_q_range(id);
            let mut out = fs::File::create(work.join(format!("api_fit_{}.txt", id + 1)))?;
            writeln!(out, "# Q  data_F(Q)  simulated_F(Q)   (fitted range of dataset {})", id + 1)?;
            for i in lo..=hi.min(q.len() - 1) {
                writeln!(out, "{:10.4} {:14.6e} {:14.6e}", q[i], d[i], f[i])?;
            }
        }
        println!("# pairs: {:?}", sim.pair_labels());
    }
    Ok(())
}
