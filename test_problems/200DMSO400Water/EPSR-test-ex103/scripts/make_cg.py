#!/usr/bin/env python3
"""Modify a copy of the EPSR-test inputs for a coarse-graining test case.

usage: make_cg.py CASE_DIR MODE
Edits CASE_DIR/200DMSO400Water.EPSR.inp (keyword values), the molecule
header lines of CASE_DIR/200DMSO400Water.ato (10th word = scalefacmol) and,
for some modes, the coarse-graining columns of the atom-type lines
(natomsph radminsph radmaxsph cgmolfac).
"""

import sys, os

case, mode = sys.argv[1], sys.argv[2]
inp = os.path.join(case, "200DMSO400Water.EPSR.inp")
ato = os.path.join(case, "200DMSO400Water.ato")

KW = {
    "ex100": {"qinscalefc": "2.0 0.0", "ntimes": "3"},
    "ex101": {"qinscalefc": "1.5 0.33", "ntimes": "3"},
    "ex102": {
        "qinscalefc": "3.0 0.5",
        "ntimes": "3",
        "sizefactor": "1.20000 0.90000 1.0E+06",
    },
    "ex103": {"qinscalefc": "1.2", "ntimes": "2", "ngrsamples": "150", "iinit": "1"},
}[mode]


def set_kw(lines, key, val):
    for i, l in enumerate(lines):
        if l[:11].strip() == key:
            rest = l[12:]
            # keep the description text after the value
            j = rest.find("               ")
            desc = rest[j:] if j >= 0 else ""
            lines[i] = f"{key:<11} {val}{desc}"
            return
    raise SystemExit(f"keyword {key} not found")


lines = open(inp).read().split("\n")
for k, v in KW.items():
    set_kw(lines, k, v)
open(inp, "w").write("\n".join(lines))


# scale factor per molecule
def scale(im, nato):
    if mode == "ex100":
        return None  # keep 1.0
    if mode == "ex101":
        if nato == 10:  # DMSO
            return [1.0, 1.25, 1.5, 1.7][im % 4]
        return [0.8, 1.0, 1.2][im % 3]  # water
    if mode == "ex102":
        if nato == 10:
            return 1.0 + 0.05 * (im % 9)
        return 1.0
    if mode == "ex103":
        return 1.0 if nato == 10 else [0.9, 1.0, 1.1, 1.3][im % 4]


out = []
im = 0
alines = open(ato).read().split("\n")
for l in alines:
    w = l.split()
    if len(w) >= 10 and w[7] in ("T", "F") and w[0].isdigit():
        im += 1
        s = scale(im, int(w[0]))
        if s is not None:
            # replace the 10th word in place
            pos = 0
            for k in range(10):
                pos = l.index(w[k], pos)
                if k < 9:
                    pos += len(w[k])
            l = l[:pos] + f"{s:.5f}".rjust(len(w[9])) + l[pos + len(w[9]) :]
    out.append(l)

# coarse-graining columns of atom type lines (ex102 / ex103)
TYPES = {
    "ex102": {
        "S": "-0.08 0.3 1.9 2.0",
        "C": "1.0 0.5 2.0 2.0",
        "O": "1.0 0.0 1.5 2.0",
        "H": "1.0 0.0 0.0 2.0",
        "O1": "-0.05 0.0 1.7 1.0",
        "H1": "1.0 0.2 0.6 0.5",
    },
    "ex103": {"O1": "3.0 0.0 1.6 1.0", "H1": "2.0 0.1 0.5 1.0"},
}.get(mode, {})
if TYPES:
    for i in range(len(out) - 1):
        w = out[i].split()
        if len(w) == 3 and w[0] in TYPES and len(out[i + 1].split()) == 9:
            v = out[i + 1].split()
            out[i + 1] = "  " + "  ".join(v[:5]) + "  " + "  ".join(TYPES[w[0]].split())
open(ato, "w").write("\n".join(out))
print(mode, "molecules", im)
