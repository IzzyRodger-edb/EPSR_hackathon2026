//! Run a validation case through the `epsr-core` library and record the state
//! of the simulation after *every* EPSR iteration, for the Python
//! visualisation tools in `visualisation/epsr_vis`.
//!
//! usage: epsr-recorder CASE_DIR WORK_DIR
//!
//! The steps of CASE_DIR/steps.txt are executed in WORK_DIR exactly as
//! `validation/api_driver` does.  The only difference is in `run epsr`: instead
//! of `sim.run()` the recorder calls `sim.iterate(1..=niter)` followed by
//! `sim.write_fort35()` (which is what `run()` does) and reads the in-memory
//! results between iterations.  The accessors only borrow the simulation, so
//! the files written to WORK_DIR are identical to the `epsr-rs` / Fortran ones.
//!
//! Output, in WORK_DIR/trajectory/:
//!
//! * `meta.json`     grids, labels, sizes, fitted Q ranges, molecule of each atom
//! * `history.tsv`   one line per iteration: Summary, ereq, box, density, MC counters
//! * `*.f32`         little-endian float32 arrays, one block per iteration
//!   (shapes in `meta.json`):
//!     - `fq_data`   [ndata, nq]   measured F(Q) on the simulation Q grid
//!     - `fq_err`    [ndata, nq]   its error
//!     - `fq_sim`    [ndata, nq]   simulated F(Q), running average (what EPSR fits)
//!     - `fq_inst`   [ndata, nq]   simulated F(Q) of the current configuration only
//!     - `gr_inter`, `gr_intra`  [npairs, nr]  running-average partial g(r)
//!     - `sq_inter`, `sq_intra`  [npairs, nq]  running-average partial S(Q)
//!     - `ep`, `ref_pot`         [npairs, nr]  empirical and reference potential (kJ/mol)
//!     - `pos`       [natom, 3]    atom positions (Å), molecules whole, centres in the box
//! * `pos_initial.f32` [natom, 3]  configuration before the first EPSR iteration
use epsr_core::fortran::TextFile;
use epsr_core::tools::{run_tool, ToolContext};
use epsr_core::Epsr;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

type Res<T> = Result<T, Box<dyn std::error::Error>>;

/// Same edit as validation/set_keyword.py (copied from validation/api_driver).
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

/// First value of keyword `kw` in an EPSR keyword file (e.g. `ereq` in the `.pcof`).
fn read_keyword(path: &Path, kw: &str) -> Option<String> {
    let text = fs::read(path).ok()?;
    let text: String = text.iter().map(|&b| b as char).collect();
    for line in text.lines() {
        if line.trim() == "q" {
            break;
        }
        if line.get(..12.min(line.len())).map(str::trim) == Some(kw) {
            return line.get(12..).and_then(|v| v.split_whitespace().next()).map(str::to_string);
        }
    }
    None
}

fn json_str(s: &str) -> String {
    format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\""))
}
fn json_list<T, F: Fn(&T) -> String>(v: &[T], f: F) -> String {
    format!("[{}]", v.iter().map(f).collect::<Vec<_>>().join(", "))
}

/// Append-only float32 array files, one per recorded quantity.
struct Recorder {
    dir: PathBuf,
    files: Vec<(String, BufWriter<fs::File>)>,
    history: BufWriter<fs::File>,
    step: usize,
    meta_written: bool,
    shape: (usize, usize, usize, usize), // nq, nr, natom, npairs (checked on every run)
}

impl Recorder {
    fn new(work: &Path) -> Res<Recorder> {
        let dir = work.join("trajectory");
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir)?;
        let mut history = BufWriter::new(fs::File::create(dir.join("history.tsv"))?);
        writeln!(
            history,
            "step\trun\titer\tereq\tfoldsv\tfoldsb\tenergd\tpressd\tenergabs\trfactot\tchisq\tsentb\tepcup\tepccoeff\t\
             box\trho\tsizefactor\tstepm\tstepr\tstepmi\tstepri\t\
             intra_att\tintra_rej\tirot_att\tirot_rej\trot_att\trot_rej\ttrans_att\ttrans_rej\trfacd"
        )?;
        Ok(Recorder { dir, files: Vec::new(), history, step: 0, meta_written: false, shape: (0, 0, 0, 0) })
    }

    fn put(&mut self, name: &str, v: &[f32]) -> Res<()> {
        if !self.files.iter().any(|(n, _)| n == name) {
            let f = BufWriter::new(fs::File::create(self.dir.join(format!("{name}.f32")))?);
            self.files.push((name.to_string(), f));
        }
        let w = &mut self.files.iter_mut().find(|(n, _)| n == name).unwrap().1;
        for x in v {
            w.write_all(&x.to_le_bytes())?;
        }
        Ok(())
    }

    fn write_positions(path: &Path, pos: &[[f32; 3]]) -> Res<()> {
        let mut w = BufWriter::new(fs::File::create(path)?);
        for p in pos {
            for x in p {
                w.write_all(&x.to_le_bytes())?;
            }
        }
        Ok(())
    }

    /// Called once per run, after `Epsr::setup` and before the first iteration.
    fn begin_run(&mut self, sim: &Epsr, case_name: &str, inp: &str) -> Res<()> {
        let shape = (sim.nq(), sim.nr(), sim.n_atoms(), sim.n_pairs());
        if self.meta_written {
            if shape != self.shape {
                return Err(format!("grid/system changed between runs: {:?} -> {:?}", self.shape, shape).into());
            }
            return Ok(());
        }
        self.shape = shape;
        Self::write_positions(&self.dir.join("pos_initial.f32"), &sim.atom_positions())?;
        let mut mol_of_atom = vec![0usize; sim.n_atoms()];
        for (im, r) in sim.molecule_atoms().into_iter().enumerate() {
            for ia in r {
                mol_of_atom[ia] = im;
            }
        }
        let datasets: Vec<String> = (0..sim.n_datasets()).map(|k| sim.data[k].infili.trim().to_string()).collect();
        let ranges: Vec<(usize, usize)> = (0..sim.n_datasets()).map(|k| sim.data_q_range(k)).collect();
        let f = |x: &f32| format!("{x:e}");
        let mut m = String::from("{\n");
        m += &format!("  \"case\": {},\n  \"inp\": {},\n", json_str(case_name), json_str(inp));
        m += &format!(
            "  \"nq\": {}, \"nr\": {}, \"natom\": {}, \"nmol\": {}, \"npairs\": {}, \"ndata\": {},\n",
            sim.nq(),
            sim.nr(),
            sim.n_atoms(),
            sim.n_molecules(),
            sim.n_pairs(),
            sim.n_datasets()
        );
        m += &format!("  \"q\": {},\n  \"r\": {},\n", json_list(sim.q_grid(), f), json_list(sim.r_grid(), f));
        m += &format!("  \"atom_types\": {},\n", json_list(&sim.atom_types(), |s: &String| json_str(s.trim())));
        m += &format!("  \"pairs\": {},\n", json_list(&sim.pair_labels(), |(a, b)| format!("[{}, {}]", json_str(a), json_str(b))));
        m += &format!("  \"datasets\": {},\n", json_list(&datasets, |s: &String| json_str(s)));
        m += &format!("  \"q_range\": {},\n", json_list(&ranges, |(a, b)| format!("[{a}, {b}]")));
        m += &format!("  \"atom_labels\": {},\n", json_list(&sim.atom_labels(), |s: &String| json_str(s.trim())));
        m += &format!("  \"molecule_of_atom\": {},\n", json_list(&mol_of_atom, |i| i.to_string()));
        m += "  \"arrays\": {\n";
        let (nq, nr, na, np, nd) = (sim.nq(), sim.nr(), sim.n_atoms(), sim.n_pairs(), sim.n_datasets());
        let arrays = [
            ("fq_data", nd, nq),
            ("fq_err", nd, nq),
            ("fq_sim", nd, nq),
            ("fq_inst", nd, nq),
            ("gr_inter", np, nr),
            ("gr_intra", np, nr),
            ("sq_inter", np, nq),
            ("sq_intra", np, nq),
            ("ep", np, nr),
            ("ref_pot", np, nr),
            ("pos", na, 3),
        ];
        m += &arrays.iter().map(|(n, a, b)| format!("    \"{n}\": [{a}, {b}]")).collect::<Vec<_>>().join(",\n");
        m += "\n  }\n}\n";
        fs::write(self.dir.join("meta.json"), m)?;
        self.meta_written = true;
        Ok(())
    }

    /// Called after every `sim.iterate(iter)`.
    fn snapshot(&mut self, sim: &Epsr, run: usize, iter: i32, ereq: &str) -> Res<()> {
        self.step += 1;
        let (nq, nd, np) = (sim.nq(), sim.n_datasets(), sim.n_pairs());
        for k in 0..nd {
            self.put("fq_data", sim.data_fq(k))?;
            self.put("fq_err", &sim.errs.col(k + 1)[..nq])?;
            self.put("fq_sim", &sim.simulated_fq(k))?;
            // current configuration only: same weighted sum as simulated_fq, on qwind
            let inn = k as i32 + 1;
            let w = &sim.st.wts;
            let inst: Vec<f32> = (1..=nq as i32)
                .map(|i| {
                    let mut s: f32 = 0.0;
                    for ic in 1..=np as i32 {
                        s = ((s as f64 + w.ciwt[(inn, ic)] * sim.qwind[(1, i, ic)] as f64) + w.cwt[(inn, ic)] * sim.qwind[(2, i, ic)] as f64) as f32;
                    }
                    s
                })
                .collect();
            self.put("fq_inst", &inst)?;
        }
        for ip in 0..np {
            self.put("gr_inter", sim.gr_inter(ip))?;
            self.put("gr_intra", sim.gr_intra(ip))?;
            self.put("sq_inter", sim.sq_inter(ip))?;
            self.put("sq_intra", sim.sq_intra(ip))?;
            self.put("ep", sim.empirical_potential(ip))?;
            self.put("ref_pot", sim.reference_potential(ip))?;
        }
        let pos: Vec<f32> = sim.atom_positions().into_iter().flatten().collect();
        self.put("pos", &pos)?;

        let s = &sim.summary;
        let m = &sim.st.mol;
        // MOLFIT counters of the last MC cycle, summed over molecule types: row 1 attempts, row 2 rejects
        let mut c = [[0i64; 4]; 2];
        for (j, row) in c.iter_mut().enumerate() {
            for imt in 1..=m.nmoltype {
                row[0] += m.nmovei[(j + 1, imt)] as i64;
                row[1] += m.ninter[(j + 1, imt)] as i64;
                row[2] += m.nwhole[(j + 1, imt)] as i64;
                row[3] += m.nmovem[(j + 1, imt)] as i64;
            }
        }
        let rd: Vec<String> = s.rfacd.iter().map(|x| format!("{x:e}")).collect();
        // between iterations lengths are in box-fraction units: convert step sizes to Å
        writeln!(
            self.history,
            "{}\t{}\t{}\t{}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.step,
            run,
            iter,
            ereq,
            s.foldsv,
            s.foldsb,
            s.energd,
            s.pressd,
            s.energabs,
            s.rfactot,
            s.chisq,
            s.sentb,
            s.epcup,
            s.epccoeff,
            m.dmax,
            sim.st.dist.rho,
            m.sizefactor,
            m.stepm * m.dmax,
            m.stepr,
            m.stepmi * m.dmaxmol,
            m.stepri,
            c[0][0],
            c[1][0],
            c[0][1],
            c[1][1],
            c[0][2],
            c[1][2],
            c[0][3],
            c[1][3],
            rd.join(",")
        )?;
        Ok(())
    }

    fn flush(&mut self) -> Res<()> {
        for (_, f) in self.files.iter_mut() {
            f.flush()?;
        }
        self.history.flush()?;
        Ok(())
    }
}

fn main() -> Res<()> {
    let a: Vec<String> = std::env::args().skip(1).collect();
    if a.len() != 2 {
        eprintln!("usage: epsr-recorder CASE_DIR WORK_DIR");
        std::process::exit(2);
    }
    let case = fs::canonicalize(&a[0])?;
    let case_name = case.file_name().unwrap().to_string_lossy().to_string();
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
    // the interactive setup programs need the folder as current directory (validation/README.md)
    std::env::set_current_dir(&work)?;
    let dir = PathBuf::from("./");
    let mut rec = Recorder::new(&work)?;
    let mut nrun = 0usize;
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
                    nrun += 1;
                    let mut sim = Epsr::setup(&dir, &args[0], false)?;
                    rec.begin_run(&sim, &case_name, &args[0])?;
                    let pcof = sim.st.sym.full(&format!("{}.pcof", args[0].trim_end_matches(".EPSR.inp")));
                    let ereq = read_keyword(&pcof, "ereq").unwrap_or_else(|| "0".into());
                    // == sim.run(), with a snapshot after each iteration
                    for iter in 1..=sim.ctl.niter {
                        sim.iterate(iter)?;
                        rec.snapshot(&sim, nrun, iter, &ereq)?;
                    }
                    sim.write_fort35()?;
                    rec.flush()?;
                    let s = &sim.summary;
                    eprintln!("  run {nrun}: {} iterations, energy {:.4}, pressure {:.4}, R-factor {:.4e}", sim.ctl.niter, s.energd, s.pressd, s.rfactot);
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
    rec.flush()?;
    eprintln!("recorded {} iterations in {}", rec.step, rec.dir.display());
    Ok(())
}
