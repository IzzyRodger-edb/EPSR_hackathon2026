# Water from scratch

EPSR refinement of liquid water against neutron data (SANDALS: H2O, D2O and HDO 50:50).
The only input files are the three measured data files (`validation/water/*.mint01`).
[`water_from_scratch.py`](water_from_scratch.py) defines everything else in Python and
writes it out:

| input | defined in the script as | written by |
|---|---|---|
| molecule geometry and restraints | `ATOMS`, `OH_BOND`, `HOH_ANGLE`, `BONDS`, `ANGLES`, `SETTINGS` | `mol_file_text()` → `Water.mol` → `makemole` |
| reference potential (ε, σ, mass, charge per atom type) | `ATOM_TYPES` | the `potential` lines of `Water.mol` |
| box (1000 molecules, 0.1 atoms/Å³) | `N_MOLECULES`, `ATOMIC_DENSITY` | `mixato`, `randomise`, `fmole` |
| neutron weights | `DATASETS` (hydrogen isotope abundances) | `tools.write_keyword_file` → `nwts` |
| EPSR control file | `DATASETS`, `DATA_TYPE`, `R_MIN_BACKGROUND` | `epsr_input_text()`: file names and per-dataset lines only; EPSR supplies every other keyword's default |

Run from the `python/` folder (the `plot` extra gives the figures, `viewer` the 3D image):

```bash
uv run python examples/07_water_from_scratch/water_from_scratch.py           # ~7 minutes
uv run python examples/07_water_from_scratch/water_from_scratch.py --quick   # 10 iterations per run
```

The schedule is that of `validation/water/steps.txt`: a scratch start, 99 iterations on
the reference potential, then two runs of 100 iterations with the empirical potential
(EReq 3 kJ/mol).

## Output

The script deletes and rewrites `output/` on every run:

- `epsr_run/`: the EPSR folder (`.ato`, `.pcof`, `.EPSR.inp`, `.erg`, every `.EPSR.*` result).
- `trajectory/`: every iteration, recorded by `epsr.vis` (reload with `vis.Trajectory.load`).
- `figures/`: `01`–`09` and `summary.csv`/`summary.md` from `vis.make_report`; `10` the partial
  g(r) against the experimental partials of Soper (2000); `11` the 3D configuration; `12` the
  reference and empirical potentials; `fits.png` the final fit to each dataset.
- `summary.txt`: the console summary (type and pair parameters, R-factors, peaks).

The full run reproduces the validated water calculation (Fortran, Rust CLI and wrapper,
`validation/water`). The files in `epsr_run/` are byte-identical to that run except two:
the date line of `waterbox.EPSR.out`, and `waterbox.EPSR.qdr`, which every
`epsr.Simulation` rewrites when it reads the folder back.

`mixato` prints `sh: 1: ./Water.ato: Permission denied`. This is harmless: without a
`system_commands.txt`, its file-listing command is empty (the Fortran does the same).
