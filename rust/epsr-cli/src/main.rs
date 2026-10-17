//! `epsr-rs` — command-line equivalent of the Fortran EPSR programs.
//!
//! Usage (same arguments as the Fortran executables):
//! ```text
//! epsr-rs <folder> <command> [<input file> ...] [--quiet]
//! ```
//! e.g. `epsr-rs EPSR-test/ epsr 200DMSO400Water.EPSR.inp` or
//! `epsr-rs EPSR-test/ nwts DMSODWaterD.NWTS.dat`.  Programs that ask
//! questions (makeato, mixato, ...) read the answers from standard input.
//! Files in `<folder>` are created/updated in place, exactly as the Fortran
//! programs do.

use std::path::PathBuf;
use std::process::ExitCode;

fn usage() -> ExitCode {
    eprintln!("usage: epsr-rs <folder> <command> [<input file> ...] [--quiet]");
    eprintln!("commands: epsr {}", epsr_core::tools::COMMANDS.join(" "));
    ExitCode::from(2)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let quiet = args.iter().any(|a| a == "--quiet" || a == "-q");
    let pos: Vec<String> = args.iter().filter(|a| *a != "--quiet" && *a != "-q").cloned().collect();
    if pos.len() < 2 {
        return usage();
    }
    let dir = PathBuf::from(&pos[0]);
    let command = pos[1].to_ascii_lowercase();
    let t0 = std::time::Instant::now();
    if command == "epsr" {
        if pos.len() < 3 {
            return usage();
        }
        let mut sim = match epsr_core::Epsr::setup(&dir, &pos[2], !quiet) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("epsr-rs: {e}");
                return ExitCode::FAILURE;
            }
        };
        if let Err(e) = sim.run() {
            eprintln!("epsr-rs: {e}");
            return ExitCode::FAILURE;
        }
        let s = &sim.summary;
        println!(
            "epsr-rs: iteration {} done in {:.2} s  energy/molecule {:.4}  R-factor {:.4e}",
            s.iter,
            t0.elapsed().as_secs_f64(),
            s.energd,
            s.rfactot
        );
        return ExitCode::SUCCESS;
    }
    if !epsr_core::tools::COMMANDS.contains(&command.as_str()) {
        return usage();
    }
    let stdin = epsr_core::fortran::TextFile::from_reader(Box::new(std::io::BufReader::new(std::io::stdin())));
    let mut ctx = epsr_core::tools::ToolContext::new(&dir, pos[2..].to_vec(), stdin, !quiet);
    match epsr_core::tools::run_tool(&command, &mut ctx) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("epsr-rs: {e}");
            ExitCode::FAILURE
        }
    }
}
