//! Regression test: run the engine on a copy of `../EPSR-test` (the
//! DMSO/water test case supplied with the repository) and check that every
//! output file is identical to the output of the original Fortran program.
//!
//! The Fortran outputs are represented by hashes in `reference_hashes.txt`
//! (see `tools/make_reference_hashes.py`); the date line of the `.out` file
//! and the timing columns of `fort.35` are masked before hashing.
//!
//! Run with `cargo test --release` (≈ 3 s).

use std::fs;
use std::path::{Path, PathBuf};

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

fn fnv1a64(b: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &c in b {
        h ^= c as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
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

#[test]
fn dmso_water_test_case_matches_fortran() {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let src = manifest.join("../../EPSR-test");
    if !src.exists() {
        eprintln!("EPSR-test folder not found at {} — skipping", src.display());
        return;
    }
    let work = std::env::temp_dir().join(format!("epsr-rs-test-{}", std::process::id()));
    let _ = fs::remove_dir_all(&work);
    copy_dir(&src, &work);

    let mut sim = epsr_core::Epsr::setup(&work, "200DMSO400Water.EPSR.inp", false).expect("setup");
    sim.run().expect("run");

    let reference = fs::read_to_string(manifest.join("tests/reference_hashes.txt")).unwrap();
    let mut failures = Vec::new();
    let mut n = 0;
    for line in reference.lines().filter(|l| !l.starts_with('#') && !l.trim().is_empty()) {
        let (hash, name) = line.split_once(' ').unwrap();
        let data = fs::read(work.join(name)).unwrap_or_else(|_| panic!("missing output {name}"));
        let got = format!("{:016x}", fnv1a64(&mask(name, &data)));
        if got != hash {
            failures.push(name.to_string());
        }
        n += 1;
    }
    // Basic sanity on the API as well.
    assert_eq!(sim.n_pairs(), 21);
    assert_eq!(sim.n_datasets(), 3);
    assert_eq!(sim.q_grid().len(), 600);
    assert!((sim.summary.rfactot - 1.1478e-4).abs() < 1e-7, "rfactot = {}", sim.summary.rfactot);
    let _ = fs::remove_dir_all(&work);
    assert!(failures.is_empty(), "{} of {} files differ from the Fortran reference: {:?}", failures.len(), n, failures);
}

#[test]
fn epsr_is_send() {
    fn assert_send<T: Send>() {}
    assert_send::<epsr_core::Epsr>();
}

#[test]
fn overrides_are_applied() {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let src = manifest.join("../../EPSR-test");
    if !src.exists() {
        return;
    }
    let work = std::env::temp_dir().join(format!("epsr-rs-ovr-{}", std::process::id()));
    let _ = fs::remove_dir_all(&work);
    copy_dir(&src, &work);
    let lines = vec!["ntimes 2".to_string(), "data 2".to_string(), "tweak 1.05".to_string()];
    let sim = epsr_core::Epsr::setup_with_lines(&work, "200DMSO400Water.EPSR.inp", &lines, false).expect("setup");
    assert_eq!(sim.ctl.ntimes, 2);
    assert_eq!(sim.get_parameter("tweak", Some(1)).as_deref(), Some("1.05"));
    assert_eq!(sim.get_parameter("tweak", Some(0)).as_deref(), Some("1.0"));
    assert!((sim.data[1].tweak - 1.05).abs() < 1e-6);
    let _ = fs::remove_dir_all(&work);
}
