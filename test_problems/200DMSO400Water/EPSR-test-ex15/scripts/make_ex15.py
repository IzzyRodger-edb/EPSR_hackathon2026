#!/usr/bin/env python3
"""Create the EPSR-test-ex15 inputs from EPSR-test: a massless Q site at the
S position of every DMSO molecule, the Q atom type, the seven Q pairs in the
.pcof file (Q-O1 with a container force rmaxex/pcoefexmax), and iinit = 1.
usage: make_ex15.py EPSR-test-dir out-dir RMAXEX PCOEFEXMAX"""

import sys

src, out, rmaxex, pmax = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
L = open(f"{src}/200DMSO400Water.ato").read().split("\n")
o = [L[0], L[1]]
i = 2
while True:
    h = L[i]
    if not h[:3].strip().isdigit():
        break
    nat = int(h[:3])
    blk_end = i + 1
    # atoms: label line, coordinate line, neighbour lines (continued while
    # the next line does not start with a label)
    j = i + 1
    atoms = []
    for a in range(nat):
        start = j
        j += 2
        nn = int(L[j].split()[0])
        need = 1 + 2 * nn
        got = len(L[j].split())
        j += 1
        while got < need:
            got += len(L[j].split())
            j += 1
        atoms.append(L[start:j])
    # rotation groups up to next header or type list
    k = j
    ngr = int(L[k].split()[0])
    k += 1
    for g in range(ngr):
        k += 1  # ROT/DIH label
        n = int(L[k].split()[0])
        k += 1
        got = len(L[k].split())
        k += 1
        while got < n:
            got += len(L[k].split())
            k += 1
    if nat == 10:  # DMSO: add Q at the S position, no bonds
        h = " 11" + h[3:]
        s_xyz = atoms[0][1]
        atoms.append([" Q       11      0", s_xyz, "    0"])
    o.append(h)
    for a in atoms:
        o.extend(a)
    o.extend(L[j:k])
    i = k
# atom types: insert Q after H
rest = L[i:]
r = []
for n, l in enumerate(rest):
    r.append(l)
    if l.startswith(" H   H    0"):
        r.append(rest[n + 1])
        skip = n + 1
        r.append(" Q   Q    0")
        r.append(
            "  0.000000   0.000000   0.000000   0.000000   0.000000  0.100000E+01     0.0000     0.0000  0.100000E+01"
        )
        break
r.extend(rest[skip + 1 :])
o.extend(r)
open(f"{out}/200DMSO400Water.ato", "w").write("\n".join(o))
# pcof: add Q pairs after the existing 21
P = open(f"{src}/200DMSO400Water.pcof").read().split("\n")
idx = next(n for n, l in enumerate(P) if l.strip() == "21")
P[idx] = P[idx].replace("21", "28")
coef = P[idx + 3]
zero = "   1.00000000    " + "   0.00000000    " * (len(coef.split()) - 1)
add = []
for n, t in enumerate(["S", "C", "O", "H", "Q", "O1", "H1"]):
    a, b = ("Q", t) if t not in ("O1", "H1") else ("Q", t)
    add.append(f"{a:<3} {b:<3}      {22+n}     1.00000")
    add.append(
        f"   0.00000000       0.00000000       0.00000000       {rmaxex}       {pmax}"
        if t == "O1"
        else "   0.00000000       0.00000000       0.00000000    "
    )
    add.append(zero)
end = idx + 1 + 21 * 3
P[end:end] = add
open(f"{out}/200DMSO400Water.pcof", "w").write("\n".join(P))
I = (
    open(f"{src}/200DMSO400Water.EPSR.inp")
    .read()
    .replace("iinit       0 ", "iinit       1 ", 1)
)
open(f"{out}/200DMSO400Water.EPSR.inp", "w").write(I)
