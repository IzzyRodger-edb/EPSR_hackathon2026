#!/usr/bin/env python3
"""Compare the files produced by the Rust and Fortran EPSR runs.

usage: compare_outputs.py REFERENCE_DIR CANDIDATE_DIR [--quiet]

For every file in REFERENCE_DIR that differs from the pristine input (i.e.
every file EPSR writes), the same file in CANDIDATE_DIR is compared:

* byte-identical files are reported as IDENTICAL;
* otherwise the first differing line is shown and, for numeric tables, the
  maximum absolute and relative difference over all numbers is reported.

Lines that legitimately differ between runs are masked before comparing:
* `.out` line 2 (current date/time),
* `fort.35` timing columns (the last two integers of each energy line),
* `run.log` (console output) is ignored.

Exit status is 0 if every compared file is identical after masking.
"""

import os
import re
import sys

IGNORE = {"run.log"}
NUM = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")


def mask(name, lines):
    if name.endswith(".out") and len(lines) > 1 and "Current date" in lines[1]:
        lines = lines[:]
        lines[1] = "<date masked>"
    if name == "fort.35":
        out = []
        for l in lines:
            if l.startswith("#"):
                out.append(l)
            else:
                out.append(" ".join(l.split()[:3]))
        lines = out
    return lines


def mask_more(name, lines):
    # mapgrquad.f never sets grwerr(:,1) and writes stack garbage as the error
    # of the g_centres column (columns 31-45 of the data lines); the port writes 0.
    if name.endswith(".mapgr.g01"):
        lines = [
            l if l.startswith("#") else l[:30] + "<masked>" + l[45:] for l in lines
        ]
    return lines


def numbers(line):
    return [float(x.replace("D", "E").replace("d", "e")) for x in NUM.findall(line)]


def compare(name, ref, cand):
    with open(ref, "rb") as f:
        a = f.read()
    with open(cand, "rb") as f:
        b = f.read()
    if a == b:
        return True, "IDENTICAL"
    la = mask_more(name, mask(name, a.decode("latin-1").split("\n")))
    lb = mask_more(name, mask(name, b.decode("latin-1").split("\n")))
    if la == lb:
        return True, "IDENTICAL (after masking date/timing)"
    first = None
    maxabs = 0.0
    maxrel = 0.0
    ndiff = 0
    for i, (x, y) in enumerate(zip(la, lb)):
        if x != y:
            ndiff += 1
            if first is None:
                first = (i + 1, x, y)
            nx, ny = numbers(x), numbers(y)
            if len(nx) == len(ny):
                for u, v in zip(nx, ny):
                    d = abs(u - v)
                    maxabs = max(maxabs, d)
                    s = max(abs(u), abs(v))
                    if s > 0:
                        maxrel = max(maxrel, d / s)
    msg = f"DIFFERENT: {ndiff} lines differ"
    if len(la) != len(lb):
        msg += f"; line counts {len(la)} vs {len(lb)}"
    msg += f"; max |diff| {maxabs:.3e}, max rel {maxrel:.3e}"
    if first:
        msg += f"\n      first at line {first[0]}:\n      ref : {first[1][:150]!r}\n      cand: {first[2][:150]!r}"
    return False, msg


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    ref, cand = sys.argv[1], sys.argv[2]
    ok_all = True
    for name in sorted(os.listdir(ref)):
        if name in IGNORE:
            continue
        p = os.path.join(ref, name)
        q = os.path.join(cand, name)
        if not os.path.isfile(p):
            continue
        if not os.path.exists(q):
            print(f"{name:32s} MISSING in candidate")
            ok_all = False
            continue
        ok, msg = compare(name, p, q)
        ok_all &= ok
        if not ok or "--quiet" not in sys.argv:
            print(f"{name:32s} {msg}")
    # files written only by the candidate are also differences
    for name in sorted(os.listdir(cand)):
        if name in IGNORE or not os.path.isfile(os.path.join(cand, name)):
            continue
        if not os.path.exists(os.path.join(ref, name)):
            print(f"{name:32s} ONLY in candidate")
            ok_all = False
    print("ALL IDENTICAL" if ok_all else "SOME FILES DIFFER")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
