# The auxiliary EPSR programs (`epsr-core/src/tools/`)

Besides the `epsr` engine, the Fortran EPSR suite has 35 setup and analysis programs that the EPSR
shell runs by command name. All of them are translated. Each one is a module `tools/<command>.rs` with
an entry point `run(&mut ToolContext)`. `tools::run_tool(command, ctx)` dispatches on the command name.

```
epsr-rs <folder> <command> [args...] [--quiet]   < answers     # same arguments as the Fortran program
```

`ToolContext` (`tools/mod.rs`) holds:

* `dir`: the working folder;
* `args`: the words after the command;
* `stdin`: a streaming `TextFile` for the answers the Fortran reads from unit 5;
* `verbose`.

Programs driven by a keyword file (`*.NWTS.dat`, `*.PARTIALS.dat`, ...) call `ctx.open_input(&mut st, &TABLE)`. This installs the
program's keyword table from `symtab_tools.rs`, which `tools/gen/gen_symtabs.py` generates from the Fortran
`setup_symbols_*` routines. It then reads the file with the same `symbols.rs` machinery as `epsr`. On exit the
file is rewritten exactly as `save_files_others` does.

Every program is verified **byte-identical** against the original Fortran, built with gfortran -O2
-fno-tree-vectorize. The test cases are listed in the table below; [TEST_CASES.md](TEST_CASES.md) has the
full list.

## Programs

| command | Fortran source | Rust module | input | test cases |
|---|---|---|---|---|
| `nwts` | nwts.f (+ `check_symbol` branch) | nwts.rs | `.NWTS.dat` keywords | ex20 ex21 ex22 |
| `xwts` | xwts.f (`xwts`, `form_fact`) | xwts.rs | `.XWTS.dat` keywords | ex23 |
| `epsrwts` | epsrwts.f | epsrwts.rs | `.ato` + stdin | ex24 |
| `epsrwtsx` | epsrwtsx.f | epsrwtsx.rs | `.ato` + stdin | ex25 |
| `changeato` | changeato.f + process_symbols.f (interactive loop) | changeato.rs | stdin commands | ex26 ex27 |
| `makemole` | makemole.f → `makemole`, `get_rotgrp`, `get_whole_rotgrp` | makemole.rs | `.mol`/`.gmol` | ex28 ex29 |
| `makeato` | makeato.f | makeato.rs | stdin | ex30 |
| `mixato` | mixato.f | mixato.rs | stdin | ex31 ex32 ex33 ex35 ex39 |
| `fmole` | fmole.f (`fmole`, `initpt_fmole`, `update_ato`) | fmole.rs | stdin | ex33 ex34 ex35 ex39 |
| `randomise` | introtcluster.f (`randomise`) | randomise.rs | stdin | ex35 ex39 |
| `addato` | addato.f | addato.rs | stdin | ex36 |
| `makelattice` | makelattice.f (+ `convertalphabetagamma`) | makelattice.rs | stdin | ex37 |
| `makelatticeato` | makelatticeato.f | makelatticeato.rs | stdin | ex38 ex39 |
| `readxyz` | readxyz.f | readxyz.rs | `.xyz` + stdin | ex40 ex41 |
| `writexyz` | writexyz.f | writexyz.rs | `.ato` + stdin | ex42 ex43 |
| `readpdb` | readpdb.f | readpdb.rs | `.pdb` + stdin | ex44 ex45 ex49 |
| `bonds` | bondsclus.f | bonds.rs | `.ato` + stdin | ex46 |
| `alignmole` | alignmole.f (`alignmole`, `displacement`) | alignmole.rs | stdin | ex47 ex48 |
| `partials` | partials.f | partials.rs | `.PARTIALS.dat` | ex50 ex51 |
| `coord` | coord.f | coord.rs | `.COORD.dat` | ex52 |
| `triangles` | triangles.f | triangles.rs | `.TRIANGLES.dat` | ex53 |
| `torangles` | torangles.f | torangles.rs | `.TORANGLES.dat` | ex54 |
| `voids` | voids.f | voids.rs | `.VOIDS.dat` | ex55 ex56 |
| `fluctuations` | fluctuations.f | fluctuations.rs | `.FLUCTUATIONS.dat` | ex57 ex59 |
| `local_density` | local_density.f | local_density.rs | `.LOCAL_DENSITY.dat` | ex58 |
| `clusters` | clusters.f | clusters.rs | `.CLUSTERS.dat` | ex60 ex61 |
| `chains` | chains_sp_criterion.f + chains_routines.f | chains.rs | `.CHAINS.dat` | ex62 ex63 |
| `rings` | rings_sp_criterion.f | rings.rs | `.RINGS.dat` | ex64 ex65 |
| `lsi` | lsi.f | lsi.rs | `.LSI.dat` | ex66 ex67 |
| `sharm` | sharm_rev2.f + sharm_routines_rev1.f | sharm.rs | `.SHARM.dat` | ex70 ex71 ex72 ex77 ex79 |
| `sdf` | sdf_rev2.f | sdf.rs | `.SDF.dat` | ex73 ex74 ex78 |
| `sdfcube` | sdfcube.f | sdfcube.rs | `.SDFCUBE.dat` | ex75 ex76 ex79 |
| `mapgr` | mapgrquad.f + plotsharm_routines.f | mapgr.rs | `.MAPGR.dat` | ex70–ex74 ex77 |
| `mapgrentropy` | mapgrentropy.f + plotsharm_routines.f | mapgrentropy.rs | `.MAPGRENTROPY.dat` | ex70–ex72 ex77 ex78 |
| `dockato` | dockato.f (`dockato`, `rotatem`, `search_unsaturated`, ...) | dockato.rs | args + stdin | ex80–ex87 |

## Shared helper files

| file | contents | used by |
|---|---|---|
| `anal_common.rs` | symbol-value accessors, centre/type matching, `outgofz` histogram output, `fort.14` writer, `ran1` state files, re-exports the gfortran `sind`/`cosd` intrinsics | partials, coord, triangles, torangles, voids, fluctuations, local_density |
| `groupc_common.rs` | `chsymbol`/`chatwt` (element data), `convertalphabetagamma`/`setup_lattice`, centre of mass, coordinate rotations, `dihedral`, `get_rotgrp`, a `moloutc` variant | readxyz, writexyz, readpdb, bonds, alignmole |
| `setup_common.rs` | prompt/answer reading (`command_text`, `read_text`, `read_real5`...), file-name helpers, `rotatm` | makeato, mixato, fmole, addato, randomise, makelattice, makelatticeato, makemole |
| `sharm_routines.rs` | rotation matrices `dlm`, Clebsch-Gordan coefficients, Euler angles, smoothing, coefficient-file I/O, gfortran `sind`/`cosd` | sharm, sdf, sdfcube, mapgr, mapgrentropy |
| `sears91.rs` | neutron scattering lengths (`sears91_gudrun.inc`) | nwts |
| `chains.rs` | `match`, `shortest_path` (chains_routines.f) | chains, rings, lsi |
| `clusters.rs` | `outgofz` (histogram output of read_and_save_files.f), `open_input_with` | clusters, chains, lsi (rings via chains) |

## Fortran defects found (and how the reference was made)

* **`readxyz`, `readpdb`.** The binaries crash with SIGSEGV at the first call of `ctype`. That routine has optional
  arguments, but it is called with 11 of 13 arguments and no explicit interface.
  * Reference: built with an explicit `INTERFACE` block added (`rust/tools/build_patched_readxyz_readpdb.sh`).
  * Port: behaves like the patched program.
* **`sharm`, `sdf`.** The run-time FORMAT `'(1h(i1...'` lacks a comma after the Hollerith `1h(`. Current gfortran
  rejects it; older compilers accepted it.
  * Reference: built with the comma inserted (`build_fortran_reference.sh` keeps the original as `*.orig`).
* **`makeato`.** The program crashes (SIGSEGV) in `moloutc`: it never calls `ctype`, so `natomsphere`/`pseudofacmol` are
  unallocated. The buffered output is lost, so the `.ato` file is empty.
  * Port: reproduces this outcome, an empty file and an error exit (ex30).
* **`alignmole`.** A rejected trial rotation restores the molecule from work arrays that are allocated but never assigned
  (uninitialised heap memory).
  * Port: zero-filled arrays.
  * Tests: ex47/ex48 are designed so this branch does not change the output.
* **`mapgr`.** `grwerr(:,1)` is never set, so the error of the g_centres column contains stack garbage.
  * Port: writes 0.
  * Comparison: those columns are masked in `compare_outputs.py`.
* **`readpdb`.** `dihedral(im, ...)` is called with absolute atom numbers. The DIH angles are therefore only
  correct for molecule 1.
  * Port: reproduces this (ex44).
* **`dockato`.** Several branches are reproduced as the Fortran behaves; see the ex80–ex87 READMEs.
  * Removing the bonding-site type re-creates it in `ftype`.
  * The undefined `imv`, when no site is found, is reported as unsupported.

Paths that need the interactive file search (`search_for_files`), or that would index arrays out of bounds
in the Fortran, return `EpsrError::Unsupported` with a message. List them with `grep -n "unsupported(" epsr-core/src/tools/*.rs`.
