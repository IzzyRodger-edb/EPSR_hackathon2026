#!/usr/bin/env python3
"""setkw.py FILE [ds:]keyword=value ...  -- replace the value field of a keyword in a .EPSR.inp
(ds = dataset number for per-dataset keywords; the comment text is kept)."""

import sys, re

p = sys.argv[1]
L = open(p).read().split("\n")
for a in sys.argv[2:]:
    k, v = a.split("=", 1)
    ds = 0
    if ":" in k:
        ds, k = k.split(":")
        ds = int(ds)
    cur = 0
    done = False
    for i, l in enumerate(L):
        w = l.split()
        if w[:1] == ["data"] and len(w) > 1:
            cur = int(w[1])
        if w[:1] == [k] and cur == ds:
            m = re.match(r"^(\S+\s+)(.*?)(\s{10,}.*)$", l)
            L[i] = m.group(1) + v + m.group(3)
            done = True
    assert done, a
open(p, "w").write("\n".join(L))
