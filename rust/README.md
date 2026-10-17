# EPSR in Rust

A Rust translation of the whole **EPSR** suite (Empirical Potential Structure
Refinement, A. K. Soper), version 26: the `epsr` engine implemented by
`fortran_src/epsr_standalone_rev3.f` and the routines it calls, plus all 35
setup and analysis programs.  The original Fortran in `../fortran_src` is
untouched.

**Status:** the Rust reproduces *every* file written by the Fortran
byte-for-byte on the supplied test case (`../EPSR-test`, 200 DMSO + 400 water,
three neutron datasets, 12 chained runs and input variants), on the 115
`EPSR-test-ex*` feature cases (bar 5 deliberate bug fixes, see
[../FORTRAN_BUGS.md](../FORTRAN_BUGS.md)) and at every step of 300-iteration
refinements of water, silica and Cu-BTC against measured data
([../validation/](../validation/README.md)).  See
[docs/VALIDATION.md](docs/VALIDATION.md).

```
rust/
├── Cargo.toml                 workspace
├── epsr-core/                 the engine (library crate)
│   ├── src/                   see docs/ARCHITECTURE.md for a file-by-file guide
│   └── tests/                 end-to-end test against the Fortran reference
├── epsr-cli/                  `epsr-rs` executable (drop-in for `epsr`)
├── tools/                     Fortran reference build + comparison scripts
└── docs/
    ├── ARCHITECTURE.md        what each file/function does, data flow of an iteration@
    ├── FORTRAN_MAP.md         Fortran routine → Rust function; what is not ported
    ├── VALIDATION.md          how the port was verified, and known caveats
    └── PYTHON_WRAPPER_NOTES.md guidance for building the Python wrapper
```

## Build and run

Requires a recent stable Rust toolchain (edition 2021) and a C linker; the only external crate is `rayon`.

```bash
cd rust
cargo build --release
# same arguments as the Fortran program: <folder> epsr <run>.EPSR.inp
./target/release/epsr-rs ../EPSR-test/ epsr 200DMSO400Water.EPSR.inp
```

Like the Fortran, this **updates the files in the folder in place** (`.ato`,
`.pcof`, `.EPSR.inp`, `.EPSR.inpa`, `.erg`, all result files).  Work on a copy
of `EPSR-test` if you want to keep the original.  Add `--quiet` to suppress
progress messages.

Differences from the Fortran command line: `system_commands.txt` is not read
(it only configures shell commands for the interactive EPSR shell), and
`delfit.dat` / `fort.35` are written to the simulation folder rather than to
the current working directory.

From Rust code:

```rust
let mut sim = epsr_core::Epsr::setup(Path::new("run_folder"), "200DMSO400Water.EPSR.inp", false)?;
sim.run()?;                         // or: for i in 1..=sim.ctl.niter { sim.iterate(i)?; }
let gr = sim.gr_inter(0);           // partial g(r) of the first atom pair
```

## Tests

```bash
cargo test --release       # unit tests + full test case vs. Fortran reference (~3 s)
```

To re-validate against the Fortran yourself (needs `gfortran`):

```bash
tools/build_fortran_reference.sh /tmp/epsr-fortran           # builds the original code
tools/chain_test.sh /tmp/epsr-fortran/epsr ../EPSR-test /tmp/chain 10
```

## Scope

The whole EPSR suite is ported:

* **`epsr` engine.** It covers every optional feature: X-ray data, non-cubic boxes, tethered molecules and `DIH`
  groups, container forces, three-body terms, coordination constraints, coarse graining, pressure control, Bragg
  and diffuse scattering, lattice refinement, all data types, and the Gaussian/harmonic potential options. See the
  table in [docs/FORTRAN_MAP.md](docs/FORTRAN_MAP.md).
* **Setup and analysis programs.** All 35 are ported; see [docs/TOOLS.md](docs/TOOLS.md).

Every feature is verified byte-identical against the original Fortran by the `EPSR-test-ex*` cases in
`test_problems/200DMSO400Water/`. [docs/TEST_CASES.md](docs/TEST_CASES.md) indexes them and [docs/VALIDATION.md](docs/VALIDATION.md) explains how.

Inputs on which the Fortran itself crashes, hangs or reads uninitialised memory raise a clear `not supported by
the Rust port: ...` error. These are listed in FORTRAN_MAP.md. The plotting front-ends (PGPLOT/gnuplot) are not
ported.
