#!/usr/bin/env python3
"""Set keyword values in an .EPSR.inp file: setkw.py FILE key=value [key=value ...]
The value field is replaced, keeping the 'keyword  value               description' layout.
"""

import sys, re

fn = sys.argv[1]
kv = dict(a.split("=", 1) for a in sys.argv[2:])
out = []
for line in open(fn):
    k = line[:12].strip()
    if k in kv and line.startswith(k):
        m = re.match(r"^(\S+\s+)(.*?)(\s{3,}.*)?$", line.rstrip("\n"))
        desc = m.group(3) or ""
        line = m.group(1) + kv.pop(k) + desc + "\n"
    out.append(line)
if kv:
    sys.exit("keywords not found: %s" % list(kv))
open(fn, "w").write("".join(out))
