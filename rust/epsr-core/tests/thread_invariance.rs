//! Thread-count invariance: every rayon-parallelised region must produce
//! byte-identical output at any thread count.  This runs the DMSO/water test
//! case once in a 1-thread pool and once in a 4-thread pool and compares
//! every file in the two working folders (masking the `.out` date line and
//! the `fort.35` timing columns, which legitimately differ between runs —
//! same masking as `test_case.rs`).
//!
//! Unlike `test_case.rs` this does not depend on the Fortran reference
//! hashes, so it also holds on platforms whose libm differs from glibc
//! (e.g. Windows), where the reference-hash test fails for the serial code
//! already.
//!
//! Run with `cargo test --release` (two full engine runs, ≈ 10 s).

use std::fs;
use std::path::{Path, PathBuf};

/// Same masking as `test_case.rs`: drop the `.out` date line and the timing
/// columns of `fort.35`.
fn mask(name: &str, data: &[u8]) -> Vec<u8> {
    let text: String = data.iter().map(|&b| b as char).collect();
    let mut lines: Vec<String> = text.split('\n').map(|s| s.to_string()).collect();
    if name.ends_with(".out") && lines.len() > 1 {
        lines[1] = "<date>".into();
    }
    if name == "fort.35" {
        lines = lines
            .into_iter()
            .map(|l| if l.starts_with('#') { l } else { l.split_whitespace().take(3).collect::<Vec<_>>().join(" ") })
            .collect();
    }
    lines.join("\n").chars().map(|c| c as u8).collect()
}

fn copy_dir(src: &Path, dst: &Path) {
    fs::create_dir_all(dst).unwrap();
    for e in fs::read_dir(src).unwrap() {
        let e = e.unwrap();
        if e.file_type().unwrap().is_file() {
            fs::copy(e.path(), dst.join(e.file_name())).unwrap();
        }
    }
}

fn run_case(work: &Path, threads: usize) {
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    pool.install(|| {
        let mut sim = epsr_core::Epsr::setup(work, "200DMSO400Water.EPSR.inp", false).expect("setup");
        sim.run().expect("run");
    });
}

fn file_names(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

#[test]
fn output_is_identical_at_any_thread_count() {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let src = manifest.join("../../EPSR-test");
    if !src.exists() {
        eprintln!("EPSR-test folder not found at {} — skipping", src.display());
        return;
    }
    let work1 = std::env::temp_dir().join(format!("epsr-rs-thr1-{}", std::process::id()));
    let work4 = std::env::temp_dir().join(format!("epsr-rs-thr4-{}", std::process::id()));
    let _ = fs::remove_dir_all(&work1);
    let _ = fs::remove_dir_all(&work4);
    copy_dir(&src, &work1);
    copy_dir(&src, &work4);

    run_case(&work1, 1);
    run_case(&work4, 4);

    let names1 = file_names(&work1);
    let names4 = file_names(&work4);
    assert_eq!(names1, names4, "the two runs produced different file sets");

    let mut failures = Vec::new();
    for name in &names1 {
        let a = mask(name, &fs::read(work1.join(name)).unwrap());
        let b = mask(name, &fs::read(work4.join(name)).unwrap());
        if a != b {
            failures.push(name.clone());
        }
    }
    let _ = fs::remove_dir_all(&work1);
    let _ = fs::remove_dir_all(&work4);
    assert!(failures.is_empty(), "{} of {} files differ between 1-thread and 4-thread runs: {:?}", failures.len(), names1.len(), failures);
}
