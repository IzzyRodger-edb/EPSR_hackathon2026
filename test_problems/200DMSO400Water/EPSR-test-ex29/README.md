# EPSR-test-ex29 — `makemole`, more options

Program: `makemole`. Hand-made molecule files:

| input | what it tests |
|---|---|
| `Chain.mol < chain.stdin` | `.gmol` ethanol-like chain: two missing bond lengths and four missing angles asked on stdin, `extra` distances, `rot`, `wholerot` (the Fortran loses the `ROT` label of whole-molecule groups — it writes the zero bytes of the newly allocated element; reproduced), `introt` (intersection of the two groups), `qradius`, `changelabel` to a label without potential (parameters asked on stdin, twice), invalid temperature (ignored), missing density (asked), `pseudofac 0` (→ 1) |
| `Tip4p.mol` | 4-site water with a massless charge site `QM` (Q atom, excluded from the atom count/density), potential line with sphere parameters, `qradius` |
| `Grid.mol < grid.stdin` | the original column-map `.mol` format: labels and `*` bonds on a grid (`nspace 2`), a repeat count (`2H`), relative (`map-4+0`) and absolute (`map(3,7)`) mapped atoms, a missing angle asked on stdin; writes the `.atm` atom-number map |

Outputs: `Chain.ato`, `Chain.mol`, `Tip4p.ato`, `Tip4p.mol`, `Grid.ato`,
`Grid.mol`, `Grid.atm`.

Verification: **ALL IDENTICAL** (reference in `fortran_reference.tar.gz`).
