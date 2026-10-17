"""Copy the 200DMSO400Water run in the current folder to a run with another name.

usage: python3 scripts/clone_run.py NEW 'keyword value' ...
Creates NEW.EPSR.inp (title NEW.EPSR, fnameato NEW.ato, fnamepcof NEW.pcof),
NEW.ato, NEW.pcof and the accumulator files NEW.EPSR.{inpa,erg,g01,y01,f01,s01}
(when present), and replaces the value of each given keyword in NEW.EPSR.inp.
The data and weights files are shared."""

import os, re, shutil, sys

old, new = "200DMSO400Water", sys.argv[1]
for ext in [
    ".ato",
    ".pcof",
    ".EPSR.inpa",
    ".EPSR.erg",
    ".EPSR.g01",
    ".EPSR.y01",
    ".EPSR.f01",
    ".EPSR.s01",
]:
    if os.path.exists(old + ext):
        shutil.copy(old + ext, new + ext)
L = open(old + ".EPSR.inp").read().split("\n")
L[0] = L[0].replace(old, new, 1)
subs = dict(a.split(None, 1) for a in sys.argv[2:])
out = []
for l in L:
    w = l.split(None, 1)
    if w and w[0] in ("fnameato", "fnamepcof"):
        l = l.replace(old, new, 1)
    if w and w[0] in subs:
        m = re.match(r"^(\S+\s+)(.*?)(\s{15}.*)$", l)
        l = m.group(1) + subs[w[0]] + m.group(3)
    out.append(l)
open(new + ".EPSR.inp", "w").write("\n".join(out))
