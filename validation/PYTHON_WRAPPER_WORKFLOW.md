# EPSR calculations with the Rust port — a worked workflow for the Python wrapper

This document describes, call by call, the two calculations in this folder:
[water](water/) and [silica_bulk](silica_bulk/). They are built from scratch from the
Dissolve examples in `experimental_data/`, then run with the Rust port and with the
original Fortran. A Python wrapper has to drive the same sequence. Everything below
was run and checked to give the same files as the Fortran, byte for byte
(see [README.md](README.md)).

Read `rust/docs/PYTHON_WRAPPER_NOTES.md` first for the engine API. This document
adds what you only find out by setting up a *new* simulation instead of re-running
the supplied test case. §6 lists those pitfalls.

The executable reference for this document is [api_driver/src/main.rs](api_driver/src/main.rs).
It runs a case's `steps.txt` through the `epsr-core` library, with the same calls a
PyO3 wrapper makes. Its output folder is byte-identical to the Fortran run.

---

## 1. Overview: from experimental data to a refined structure

| # | stage | program (Rust call) | reads | writes |
|---|-------|---------------------|-------|--------|
| 1 | molecule template | `makemole` (`run_tool`) | `<mol>.mol` | `<mol>.ato` |
| 2 | build the box | `mixato` (`run_tool`, answers on stdin) | `<mol>.ato` templates | `<box>.ato` |
| 3 | randomise positions/orientations | `randomise` (`run_tool`) | `<box>.ato` | `<box>.ato` |
| 4 | relax molecule geometry (molecules only) | `fmole` (`run_tool`) | `<box>.ato`, template `.ato` | `<box>.ato`, template `.ato` |
| 5 | neutron weights, one per isotopic sample | `nwts` (`run_tool`) | `<sample>.NWTS.dat`, `<box>.ato` | `<sample>.NWTStot.wts` |
| 6 | write the EPSR input | by hand / Python | — | `<box>.EPSR.inp` |
| 7 | first EPSR run (scratch start) | `Epsr::setup` + `run` | `.EPSR.inp`, `.ato`, data, `.wts` | `.pcof`, `.EPSR.inpa`, `.erg`, all results; rewrites `.EPSR.inp`, `.ato` |
| 8 | equilibrate on the reference potential | `Epsr::setup` + `run` (`ereq = 0`) | everything above | same |
| 9 | switch on the empirical potential | edit `ereq` in `.pcof` | | |
| 10 | refine, repeatedly | `Epsr::setup` + `run`, in a loop | | |

Stages 1–5 only need to be repeated when the system changes. Stage 10 is the loop
users spend their time in.

---

## 2. Setup calls (stages 1–5)

All setup programs share one entry point:

```rust
use epsr_core::fortran::TextFile;
use epsr_core::tools::{run_tool, ToolContext};

std::env::set_current_dir(&folder)?;               // REQUIRED, see §6.1
let stdin = TextFile::from_string(answers);         // "" for programs that ask nothing
let mut ctx = ToolContext::new(Path::new("./"), args, stdin, /*verbose*/ false);
run_tool(command, &mut ctx)?;
```

`args` are the words after the command name on the EPSR command line. For example,
`["waterbox", "100"]` for `fmole waterbox 100`.

### 2.1 `makemole` — molecule template

The input is a `.mol` file ("gmol" format). The water template used here
([water/Water.mol](water/Water.mol)) uses the Dissolve geometry and force field:
O–H 0.976 Å, H–O–H 107.134°, OW ε = 0.6503 kJ/mol, σ = 3.165492 Å, q = −0.82 e;
HW q = +0.41 e.

```
 .gmol           0
 atom            1  001    0.00000000       0.00000000       0.00000000               2           2           3
 atom            2  002    0.97600000       0.00000000       0.00000000               1           1
 atom            3  003   -0.28753687       0.93268352       0.00000000               1           1
bond 001 002    0.97600
bond 001 003    0.97600
angle 002 001 003  107.13400
changelabel 001 OW
changelabel 002 HW
changelabel 003 HW
pseudofac     1.00
potential OW      0.650300     3.165492    15.999400    -0.820000 O    0.100000E+01     0.000000     1.550000  0.100000E+01
potential HW      0.000000     0.000000     2.000000     0.410000 H    0.100000E+01     0.000000     0.000000  0.100000E+01
temperature  0.300000E+03
...
```

`atom` lines are `index, label, x, y, z, nneighbours, neighbours...`. `potential`
lines are `label, epsilon [kJ/mol], sigma [Å], mass, charge [e], element`, followed
by four values kept from the EPSR example files.

Single atoms (silica: [Si.mol](silica_bulk/Si.mol), [O.mol](silica_bulk/O.mol)) are
one-atom molecules with zero neighbours:

```
 .gmol           0
 atom            1  001    0.00000000       0.00000000       0.00000000               0
changelabel 001 Si
pseudofac     1.00
potential Si        0.175000     1.030000    28.085501     2.000000 Si   0.100000E+01     0.000000     1.110000  0.100000E+01
temperature  0.298000E+03
...
```

Call: `run_tool("makemole", ToolContext::new("./", ["Water"], "", false))` writes `Water.ato`.

### 2.2 `mixato` — build the box

This program is interactive. Its answers, one per line, are:

| answer | water | silica |
|--------|-------|--------|
| number of molecule types | `1` | `2` |
| per type: template name, number of molecules | `Water`, `1000` | `Si`, `1500`, `O`, `3000` |
| atomic number density [atoms/Å³] | `0.1` | `0.0664` |
| output name (no extension) | `waterbox` | `silica` |

The box is cubic, with its side derived from the density: 31.0723 Å for water and
40.7707 Å for silica. Call: `run_tool("mixato", ctx(args = [], stdin = answers))`.

### 2.3 `randomise <box>`

Randomises molecule positions and orientations. It asks nothing.

### 2.4 `fmole <box> <niter>` (molecular systems only)

Relaxes the intramolecular geometry with the template's bond and angle restraints.
In the water case: `fmole waterbox 100`.
**Remove the template `.mol` first** (§6.2).

### 2.5 `nwts <sample>.NWTS.dat` — weights per isotopic sample

The keyword file lists, for each atom type of the `.ato`, the isotopes and whether
the type is exchangeable. Water, component 2 (HW):

| sample | `abundances` | `iexchange` |
|--------|--------------|-------------|
| H2O | `0 1` (natural) | 1 |
| D2O | `2 1` (²H) | 1 |
| HDO 50:50 | `1 0.5 2 0.5` | 1 |

OW and both silica types use `0 1` (natural) with `iexchange 0`. Full files:
[water/H2O.NWTS.dat](water/H2O.NWTS.dat), [silica_bulk/SiO2.NWTS.dat](silica_bulk/SiO2.NWTS.dat).
`nwts` rewrites the `.dat` file (it pads lines and fills in defaults) and writes
`<sample>.NWTStot.wts`. That file is the `wtsfile` of an EPSR data set.

Sanity check of the HDO file: the HW–HW weight is 0.02147 = (½(−0.374 + 0.667))² × (2/3)².

A wrapper helper can create these `.dat` files from a dict, using the keyword table
`NWTS` in `rust/epsr-core/src/symtab_tools.rs`.

---

## 3. The EPSR input file (stage 6)

The Fortran `epsr` command never creates an `.EPSR.inp`; the interactive EPSR
shell does. A wrapper must write one. The files here were made from the test case's
input (`EPSR-test/200DMSO400Water.EPSR.inp`), changing only:

| keyword | water | silica | why |
|---------|-------|--------|-----|
| title (line 1) | `waterbox.EPSR` | `silica.EPSR` | becomes the run name and the output file prefix |
| `fnameato`, `fnamepcof` | `waterbox.ato/.pcof` | `silica.ato/.pcof` | |
| `ireset`, `iinit` | `1`, `1` | `1`, `1` | scratch start; EPSR sets both back to 0 itself |
| `nsumt` | `-1` | `-1` | no accumulated averages yet |
| `niter` | `1` | `1` | first run is a single refinement |
| `qmin` | `0.05 0.0` | `0.10 0.0` | silica: as Dissolve `XMin 0.1` |
| `ndata` + `data n` blocks | 3 (D2O, HDO, H2O) | 1 | one block per data file |
| `rshmin` | 0.7 | 1.2 | below the first intermolecular peak |

Each `data` block has `datafile`, `wtsfile`, `nrtype 5` (Gudrun `.mint01`),
`rshmin`, `szeros 0.0`, `tweak 1.0`, `efilereq 1.0`, `porodadd 0.0`.
Every other keyword keeps the test-case value: `nq 600`, `qstep 0.05`, `ntimes 5`,
`feedback 0.9`, `fwhmq 0.02`, `cellst 0.03`, `nsmoop 1`, and so on.

The keyword line layout is `<keyword padded to 12><value>` + 15 spaces + description.
The value is read list-directed, so the exact spacing is not critical. EPSR rewrites
the file after every run.

---

## 4. EPSR runs (stages 7–10)

One Fortran invocation of `epsr` is one `setup` + `run`:

```rust
let mut sim = epsr_core::Epsr::setup(Path::new("./"), "waterbox.EPSR.inp", false)?;
sim.run()?;            // ctl.niter iterations, then writes every output file + fort.35
let s = &sim.summary;  // energd, foldsv, pressd, energabs, rfactot, rfacd[..], ...
```

The sequence used in both cases:

| step | edit before the run | run | what happens |
|------|---------------------|-----|--------------|
| 7 | — (no `.pcof` exists) | `epsr` (niter 1) | default pair potentials, Coulomb on, EP = 0; creates `.pcof`, `.inpa`, `.erg`; sets `ireset`/`iinit` to 0 in the `.EPSR.inp` |
| 8 | `.EPSR.inp`: `niter 99` | `epsr` | 99 refinements on the reference potential only (`ereq = 0`) |
| 9 | `.pcof`: `ereq 3.0` (water) / `50.0` (silica) | | the `EReq` values of the Dissolve examples |
| 10 | `.EPSR.inp`: `niter 100` | `epsr`, twice | EP refinement; the second call restarts from the files written by the first |

In a wrapper, the edits are file edits (see `set_keyword` in
[api_driver/src/main.rs](api_driver/src/main.rs) or [set_keyword.py](set_keyword.py)).
`.EPSR.inp` keywords can also be given as `Epsr::setup_with_lines(dir, inp,
&["niter 99".into()], false)`. EPSR writes them back into the `.EPSR.inp` at the end
of the run, so the edit persists. `ereq` lives in the `.pcof`, so the `.pcof` file
must be edited.

The convergence numbers after each run are in `sim.summary`. The same numbers are
appended to `<run>.EPSR.erg` (one line per iteration) and `<run>.EPSR.out`.

### 4.1 Reading results

After `run()` the in-memory accessors hold the same arrays as the output files:

| accessor (Rust) | file | content |
|-----------------|------|---------|
| `q_grid()`, `data_fq(id)`, `simulated_fq(id)`, `data_q_range(id)` | `.u01` (fit), `.v01` (difference) | F(Q) of data set `id`: data and simulation |
| `r_grid()`, `gr_inter(ip)`, `gr_intra(ip)` | `.g01` | partial g(r) of pair `ip` |
| `sq_inter(ip)`, `sq_intra(ip)` | `.f01` (fit to partial S(Q)) | partial S(Q) |
| — | `.q01`, `.r01` | estimated partial S(Q) from the data, and its Fourier transform |
| — | `.w01`, `.x01` | direct Fourier transform of the data; combined partial g(r) |
| `empirical_potential(ip)`, `reference_potential(ip)`, `ep_coefficients(ip)` | `.p01`, `.o01` | potentials [kJ/mol] |
| `atom_positions()`, `box_size()` | `.ato` | configuration |

Each `.?01` file is a text table. Its first line starts with `#` and gives the
column labels (e.g. `#    r   OW-OW   OW-HW   HW-HW`): x, then (value, error) pairs
per pair or data set. `numpy.loadtxt(path, comments="#")`
reads the numbers. The API driver writes `api_fit_<n>.txt` (Q, data, simulation)
from the accessors after the last run, as an example.

---

## 5. Suggested Python shape for this workflow

```python
import epsr                                    # the PyO3 module
w = epsr.prepare("validation/water", "work/water")   # copy inputs to a fresh folder
w.tool("makemole", ["Water"])
w.tool("mixato", answers=["1", "Water", "1000", "0.1", "waterbox"])
w.tool("randomise", ["waterbox"])
w.remove("Water.mol")                          # see 6.2
w.tool("fmole", ["waterbox", "100"])
for s in ("H2O", "D2O", "HDO"):
    w.tool("nwts", [f"{s}.NWTS.dat"])
w.write_epsr_inp("waterbox", data=[("SLS18502-D2O.mint01", "D2O.NWTStot.wts"), ...])
w.epsr("waterbox.EPSR.inp")                    # scratch start
w.set("waterbox.EPSR.inp", niter=99); w.epsr("waterbox.EPSR.inp")
w.set("waterbox.pcof", ereq=3.0)
w.set("waterbox.EPSR.inp", niter=100)
for _ in range(2):
    sim = w.epsr("waterbox.EPSR.inp")          # returns the Simulation: summary + arrays
```

`w.tool` and `w.epsr` must set the process's current directory to the folder
(§6.1), so all calls on one folder belong on one thread.

---

## 6. Pitfalls found while setting up these cases

### 6.1 The setup programs depend on the current working directory

`mixato`, `fmole`, `makeato`, `addato`, `randomise`, `makelattice` and the other
interactive programs behave like the Fortran in three ways:

- they read `system_commands.txt` from the **current directory**;
- they list candidate files with a shell command that writes `dirlist.txt` there;
- they open *typed* template names (e.g. the `Water` answer to `mixato`) relative to it.

With the simulation folder as `ToolContext.dir` but a different current directory,
`mixato` reports "No .ato files found in the current directory" and writes no box.
The Fortran `mixato` does exactly the same (checked), so this is faithful behaviour.

The wrapper must `chdir` into the folder and pass `./`, as the EPSR shell does. The
current directory is process-wide, so **two folders cannot be set up in parallel
threads of one process**. Use separate processes. (This corrects the remark in
§7 of `rust/docs/PYTHON_WRAPPER_NOTES.md` that different folders can run in
parallel threads; it holds for `epsr` itself and the keyword-file programs, not for
the interactive setup programs.)

With no `system_commands.txt`, the listing command has an empty `ls`. The shell then
tries to *execute* the first matching file and prints
`sh: 1: ./Water.ato: Permission denied` on stderr. That message is harmless and also
occurs with the Fortran; the program falls back to the typed answer. A
`system_commands.txt` with a real `system_ls` entry would make the program offer
each listed file in turn instead. That path was not used or tested here, and it
would need different answers.

### 6.2 `fmole` with the template `.mol` present is not ported

The Fortran `fmole` first runs `makemole` on every template that has a `.mol` file
in the folder. That `makemole` runs *inside* `fmole` and shares its state, so the
regenerated template `.ato` contains `fmole`'s random-number state. The Rust port
does not implement this. It stops with
`not supported by the Rust port: fmole: template ./Water.mol exists ...`
(documented at the top of `rust/epsr-core/src/tools/fmole.rs`).
Workaround used here: delete the `.mol` after `makemole` and before `fmole`.
The two codes are then identical.

### 6.3 Starting without a `.pcof` versus with an empty one

- **No `.pcof` file:** the run uses default pair potentials with `coulfac = 1`
  (charges on) and EP coefficients (1, 0, 0, …). This is the correct scratch start.
- **A `.pcof` with the keyword header but zero pairs:** this looks like a natural
  template, but every pair is then created with `coulfac = 0` and `efacp = 0`,
  because the Fortran leaves them unset. Charges are switched off and the EP is
  frozen for that run. The port reproduces this. Don't generate such a file.

The default `.pcof` has `ereq = 0`. The empirical potential only starts to grow
after you set `ereq` (§4, step 9).

### 6.4 `ireset`/`iinit` reset themselves

After a run EPSR writes `ireset 0` and `iinit 0` back into the `.EPSR.inp`. A
wrapper that regenerates the `.EPSR.inp` from a template before each run would
reset the simulation every time. Edit the file EPSR wrote instead.

### 6.5 `nrtype 5` treats Q as histogram bin edges

For Gudrun `.mint01` files (`nrtype 5`, also `nrtype 3`), each Q value is replaced by
the midpoint between it and the next one, and the last point is dropped. The silica
file starts with a `0.0 0.0 0.0` line, which becomes a zero point at Q = 0.0125 Å⁻¹.
It lies below `qmin = 0.1`, so it is not fitted.

### 6.6 Threads and reproducibility

The Fortran is run with `OMP_NUM_THREADS=1`. The Rust engine is single-threaded and
reproduces that build bit for bit. Results depend only on the folder contents: the
random-number state is stored in the `.ato`. Running the same folder twice gives
identical files.

### 6.7 Timing

Measured on the 99-iteration equilibration runs on this machine (WSL2), with the
two codes run one after the other. One refinement takes about 1.3 s (Fortran) and
1.6 s (Rust) for water (3000 atoms), and about 2.0 s / 2.8 s for silica (4500 atoms).
Later runs were slower because other jobs shared the CPU.
