# EPSR-test-ex30 — `makeato`

**Program:** `makeato` (`fortran_src/makeato.f`): builds a `.ato` file holding one
molecule of one atom (label = first 3 characters of the file name) from answers on stdin
(epsilon, sigma, mass, charge, chemical symbol, temperature, number density).

**Commands** (`commands.txt`):
1. `makeato Ow < makeato_Ow.stdin` — name given on the command line (water-oxygen parameters).
2. `makeato < makeato_Ar.stdin` — name typed on stdin, with an extension (`Argon.ato`) that is stripped.

**Expected outputs:** `Ow.ato`, `Argon.ato` — both **empty**.
The reference Fortran program crashes (SIGSEGV) inside `moloutc`: `makeato` never calls
`ctype`, so the coarse-graining arrays `natomsphere`/`pseudofacmol` are unallocated when
`moloutc` reads them; the crash happens before the buffered output is flushed, leaving an
empty file. The Rust port reproduces this outcome (empty file, error exit
"moloutc: natomsphere/pseudofacmol not allocated ...").

**Verification:** `run_ex_test.sh` → ALL IDENTICAL (Fortran exit 139, Rust exit 1).
