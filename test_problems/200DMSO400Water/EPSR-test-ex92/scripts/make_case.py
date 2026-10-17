#!/usr/bin/env python3
"""Edit an EPSR-test copy in place.
usage: make_case.py DIR [inp:KEY=VALUE ...] [pcof:KEY=VALUE ...] [pcofline:N=TEXT ...]
  inp:KEY=VALUE   replace the value of keyword KEY in 200DMSO400Water.EPSR.inp
  pcof:KEY=VALUE  replace the value of keyword KEY in the .pcof header
  pcofline:N=TEXT replace line N (1-based) of the .pcof file by TEXT
The value field is the text between the keyword and the run of 15 blanks that
precedes the comment, as EPSR writes it."""

import sys, re

d = sys.argv[1]
files = {"inp": f"{d}/200DMSO400Water.EPSR.inp", "pcof": f"{d}/200DMSO400Water.pcof"}
txt = {k: open(v).read().split("\n") for k, v in files.items()}
for a in sys.argv[2:]:
    kind, rest = a.split(":", 1)
    key, val = rest.split("=", 1)
    if kind == "pcofline":
        txt["pcof"][int(key) - 1] = val
        continue
    L = txt[kind]
    for n, l in enumerate(L):
        if l[:12].strip() == key:
            m = re.match(r"(.{12})(.*?)(               .*)$", l)
            L[n] = m.group(1) + val + m.group(3)
            break
    else:
        sys.exit(f"keyword {key} not found in {kind}")
for k, v in files.items():
    open(v, "w").write("\n".join(txt[k]))
