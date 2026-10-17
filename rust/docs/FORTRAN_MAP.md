# Fortran → Rust map

## Which Fortran is "the" EPSR

`fortran_src` contains many revisions and copies of the same files
(`* (copy).f`, `*_rev1..rev4.f`, `*_old*.f`, `*_save*.f`, …).  The current
`epsr` program (version 26, `version.inc`) is built from these files — this
is the set compiled by `tools/build_fortran_reference.sh`:

```
epsr_standalone_rev3.f    main program, runepsr17 and most EPSR routines
ato_manipulation_rev4.f   fmolec1, atogen, press, stog, ... (modules diffuse_mod, timer, bragginten_mod)
molfit.f                  Monte Carlo moves and energy functions
read_and_save_files.f     read_input_file, read_ato, read_pcof, save_files_others
ato_read_and_save.f       minitc, moloutc, (un)normalise_coordinates
setup_epsr_rev1.f         setup_symbols_epsr  (setup_epsr.f is an older, incompatible revision)
process_symbols.f         process_commands, check_symbol
file_and_directory_commands.f  get_command_from_args, setup_input_file, ireff, search_for_symbol
auxiliary_routines.f      ran1, w_int, r_int, smoop, vector helpers, ...
cellneighbourlist.f, ergpairreference.f, minimagec.f, parse.f, bondnumbers.f
lngamma_mod.f90, reallocation_routines.f90, molecom_mod.f90, threebody_mod.f90,
coordconstraint_mod.f90, containerforces_mod.f90, stringutility.f90, pfroutines_mod.f90
system_commands.f, FCreadline_windows.f90, get_molecules.f, plot_routines.f  (link-time only)
```

Of ~216 000 lines in `fortran_src`, about 30 000 are reachable from `epsr`
and ~5 000 statements are executed by the test case (measured with gcov).

## Routine map

| Fortran routine | file | Rust |
|-----------------|------|------|
| `program epsr_standalone`, `get_command_from_args` | epsr_standalone_rev3.f, file_and_directory_commands.f | `epsr-cli/src/main.rs` |
| `runepsr17` | epsr_standalone_rev3.f | `epsr::Epsr::{setup, iterate, run}` (+ private `write_results`, `write_out_file`, `update_symbols`) |
| `setup_symbols_epsr` | setup_epsr_rev1.f | `symbols::State::setup_symbols_epsr` + `symtab::EPSR_SYMBOLS` |
| `setup_input_file` | file_and_directory_commands.f | `State::setup_input_file` |
| `read_input_file` | read_and_save_files.f | `State::read_input_file` |
| `process_commands`, `check_symbol` | process_symbols.f | `State::process_commands`, `State::check_symbol` (non-interactive branches) |
| `ireff`, `search_for_symbol`, `get_full_filename` | file_and_directory_commands.f | `Symbols::ireff`, `Symbols::search_for_symbol`, `Symbols::full` |
| `parse`, `concatenate`, `indexlastdot` | parse.f, process_symbols.f, auxiliary_routines.f | `fortran::text::{parse, concatenate, indexlastdot}` |
| `read_ato` | read_and_save_files.f | `State::read_ato` |
| `minitc` | ato_read_and_save.f | `State::minitc` |
| `ftype`, `ctype` | ato_manipulation_rev4.f | `State::ftype`, `State::ctype` |
| `bondnumbers` | bondnumbers.f | `State::bondnumbers` |
| `testsize_rotgrps`, `testsize_labrot`, `reallocate*` | molecom_mod.f90, reallocation_routines.f90 | `Mol::ensure_*`, `F1::grow`, `F2::grow` |
| `read_pcof` | read_and_save_files.f | `State::read_pcof` |
| `save_files_others` | read_and_save_files.f | `State::save_inp_file`, `save_pcof_file`, `write_pcof_pairs`, `moloutc` |
| `moloutc` | ato_read_and_save.f | `State::moloutc` |
| `normalise_coordinates`, `unnormalise_coordinates`, `fractcoord`, `wholecoord` | ato_read_and_save.f, ato_manipulation_rev4.f | `State::normalise_coordinates`, `unnormalise_coordinates` |
| `get_maxmolsize`, `centreofmassi` | ato_manipulation_rev4.f | `State::get_maxmolsize`, `centreofmassi` |
| `ran1`, `initran1`, `saveran1` | auxiliary_routines.f | `rng::Ran1` (`ran1`, `from_record`; save in `moloutc`) |
| `inita`, `j1qr`, `uniform` | epsr_standalone_rev3.f | `State::inita`, `structure::j1qr`, `structure::uniform` |
| `initg`, `r_int` | epsr_standalone_rev3.f, auxiliary_routines.f | `State::initg`, `output::r_int` |
| `initpt`, `truncate`, `truncatec` | epsr_standalone_rev3.f, ato_manipulation_rev4.f | `State::initpt`, `potential::{truncate, truncatec}` |
| `generate_mp` | ato_manipulation_rev4.f | `State::generate_mp` |
| `generate_rp`, `chsph` | epsr_standalone_rev3.f, ato_manipulation_rev4.f | `State::generate_rp`, `potential::chsph` |
| `generate_ep`, `multiply_pot_coeffs`, `abs_energy_ep` | epsr_standalone_rev3.f, ato_manipulation_rev4.f | `State::generate_ep`, `multiply_pot_coeffs`, `abs_energy_ep` |
| `calc_pwrexpr_grad`, `calc_pwrexpq`, `dnext`, `fctrln` | epsr_standalone_rev3.f, lngamma_mod.f90 | `potential::{calc_pwrexpr_grad, calc_pwrexpq, dnext, fctrln}` |
| `fit_expansion`, `expand`, `smoop` | epsr_standalone_rev3.f, auxiliary_routines.f | `State::fit_expansion`, `potential::{expand, smoop}` |
| `ptgnsq` | epsr_standalone_rev3.f | `State::ptgnsq` |
| `repulsivept` | epsr_standalone_rev3.f | `State::repulsivept` |
| `get_wts`, `wtype`, `wxinit` | epsr_standalone_rev3.f | `State::get_wts`, `wtype` (`wxinit` only zeroes X-ray arrays) |
| `interp` | epsr_standalone_rev3.f | `structure::interp` |
| `lsparabola1` | epsr_standalone_rev3.f | `potential::lsparabola1` |
| `initialise_cellneighbourlist`, `makecellneighbourlist`, `make_neighbour_list_from_cells`, `get_nnbrefmax`, `findiaibic`, `getcellreference`, `findperiodiccell`, `atomcellno`, `cellcoord` | cellneighbourlist.f | `neighbour.rs` |
| `initialise_ergpair`, `set_ergpair`, `get_neighbsum`, `make_neighbour_list_from_mols`, `save_energies`, `save_erg_pairs` | ergpairreference.f | `neighbour.rs` (`ErgPairs`) |
| `minimagec`, `minimage` | minimagec.f | `neighbour::{minimagec, minimage}` |
| `fmolec1` | ato_manipulation_rev4.f | `State::fmolec1` |
| `en_atom_inter_scale`, `en_mole_pair`, `en_atom_intra`, `en_mole_intra`, `en_mole_dihedral`, `en_change`, `en_change_intra`, `check_tethered` | molfit.f | `energy.rs` (`check_tethered` is always true for untethered molecules) |
| `molfit`, `movea`, `movei`, `movem`, `rotati`, `rotatm` | molfit.f | `montecarlo.rs` |
| `adjust_intra_wts` | ato_manipulation_rev4.f | `State::adjust_intra_wts` |
| `atogen` | ato_manipulation_rev4.f | `State::atogen` |
| `gensq`, `entrpb`, `entrpc`, `sumsq`, `sumgr` | epsr_standalone_rev3.f | `structure.rs` |
| `press` | ato_manipulation_rev4.f | `State::press` |
| `ftback` | epsr_standalone_rev3.f | `State::ftback` |
| `stog`, `lorch` | ato_manipulation_rev4.f | `structure::{stog, lorch}` |
| `sumpe`, `norma` | epsr_standalone_rev3.f | `State::sumpe`, `State::norma` |
| `w_int`, `w_intr` | auxiliary_routines.f | `output::{w_int, w_intr}` |
| `openf`, `get_date_and_time`, `initialise_timer`, `get_system_commands` | various | not needed (file handling via `OutFile`/`TextFile`; timing via `std::time`) |
| `cosd`, `sind`, `perpd`, `vecpro`, `scapro`, `get_atomicnumber` | auxiliary_routines.f | `ato.rs` |

## Optional features and their test cases

All the optional features of `epsr` are ported. Each one is verified byte-identical against the Fortran by at least one
`EPSR-test-ex*` case; [TEST_CASES.md](TEST_CASES.md) has the full list. Line numbers refer to `epsr_standalone_rev3.f`
unless stated.

| feature | trigger in input | Fortran | Rust | test cases |
|---------|------------------|---------|------|-----------|
| Recompute inverse weights matrix | `feedback`/`efilereq` changed, `iinit 1` (`get_wts` is called with `iinit`) | `invert_cwt`, `get_wts` | `weights.rs` | ex4 ex6 |
| Simulation/EP reset | `ireset 1`/`2`, `iinit 1` | `runepsr17`, `inita`, `read_pcof` | `epsr.rs`, `structure.rs`, `pcof.rs` | ex5 ex6 ex7 |
| Automatic `ereq` control | `ereqstep > 0` | 2180–2262 | `epsr.rs` (`automatic_ereq`) | ex2 ex3 |
| Exponential repulsive core, coordination feedback | `rminex pcoefex coordex` in `.pcof` | `generate_rp`, `repulsivept` | `potential.rs` | ex1 ex90 ex106 |
| Molecule size factor | `sizefactor > 1` | `generate_mp`, `molfit` | `potential.rs`, `montecarlo.rs` | ex8 ex102 |
| X-ray datasets | `*XWTS.wts` weights | `get_wts`, `xfile` branches | `weights.rs`, `epsr.rs` | ex9 ex10 |
| Additional Gaussian / modified Morse potentials | extra triples on `.pcof` pair lines | `read_pcof`, `generate_rp` | `pcof.rs`, `potential.rs` | ex11 ex12 |
| Tethered molecules, `DIH` groups | `T` flag / `DIH` in `.ato` | `molfit.f` | `energy.rs`, `ato.rs` | ex13 |
| Non-cubic boxes | 9-number `.ato` header | `minitc`, `minimagec` | `neighbour.rs`, `ato.rs` | ex14 ex105 ex121 |
| Container forces | 5-number `.pcof` pair lines | `containerforces_mod.f90` | `containerforces.rs` | ex15 ex91 ex92 |
| Harmonic repulsive core | `reppottype harmonic` | `generate_rp` 5839–5853 | `potential.rs` | ex90 ex91 |
| Gradient limit of the repulsive terms | `expecf` second value > 0 | `generate_rp`, `en_containerforce` | `potential.rs`, `containerforces.rs` | ex91 ex92 |
| Gaussian EP basis | `Gaussian T` | `calc_gaussq`, `calc_gaussr_grad` | `potential.rs` | ex93 |
| `nsmoop 0`, `refpotfac`, `rbroad` | keywords | `ftback`, `calc_pwrexpr_grad` | `structure.rs`, `potential.rs` | ex94 |
| Data types `nrtype` 1, 3, 4 | `nrtype` | 919–958 | `datafile.rs` | ex95 |
| `szeros` ≠ 0 | `szeros` | 1584–1611, `ptgnsq` | `epsr.rs`, `potential.rs` | ex96 |
| Porod scattering | `porodpwr`, `porodadd` | 450–455, 883–897 | `epsr.rs` | ex97 |
| Q-space Lorch broadening | `revlorch` | `inita` 2922–2925 | `structure.rs` (`pofr`) | ex98 |
| `qmin`, `rshmin`, `tweak`, `efilereq`, `fwhm`, `rspcrmin`/`rspcfrac` | keywords | various | `epsr.rs`, `potential.rs` | ex99 |
| Coarse graining | `qinscalefc > 1` | `inita` 2931–3160, `atogen`, `molfit.f` | `coarsegrain.rs` | ex100–ex103 |
| Pressure control | `control-p` | end of iteration, `change_boxsize` | `pressure.rs` | ex105–ex107 |
| `ref_intra` | `ref_intra` | `en_atom_intra`, `generate_rp` | `energy.rs`, `potential.rs` | ex107 |
| Missing `.pcof` / missing pairs, `ncoeffp` in input | — | `read_pcof`, `check_symbol` | `pcof.rs`, `symbols.rs` | ex108 ex109 |
| Three-body potentials | `tb` lines in `.pcof` | `threebody_mod.f90` | `threebody.rs` | ex110–ex114 |
| Coordination constraints | `coordcon` lines in `.pcof` | `coordconstraint_mod.f90` | `coordconstraint.rs` | ex115–ex119 |
| Bragg peaks | `qwidthqmax`, `hklqmin`, `mplicities` | `bragginten` (ato_manipulation_rev4.f) | `bragg.rs` | ex120 ex121 ex124 |
| Diffuse scattering | `diffuse` | `sharm_diffuse` (ato_manipulation_rev4.f) | `bragg.rs` | ex122 |
| Lattice refinement | `nlatrefine`, `latrefine` | 1301–1414 | `bragg.rs` (`structure_and_lattice`) | ex123 ex124 |

### Still unsupported

`EpsrError::Unsupported` is now raised only for inputs on which the Fortran itself fails, or that it never produces:

* `.ato` file without the ecore/dcore line and random-number state: the Fortran carries on with `dmax = 0` and hangs
  (`ato.rs`).
* Three-body group with atom types missing from the `.ato`, or with more angles at save time than at read time
  (`threebody.rs`).
* Diffuse scattering with zero multiplicity or a zero l/m step (division by zero / zero DO step), and more hkl values
  than the Fortran allocates (out-of-bounds writes) (`bragg.rs`).
* Coordination constraints that would read never-initialised Fortran data. This panics with a "not supported"
  message (`coordconstraint.rs`).
* `generate_ep(changepotfac=.true.)` is translated but unreachable: no Fortran caller passes `.true.`.
* The interactive editing branches of `process_commands`, which the `epsr` shell uses, not `epsr` itself.
