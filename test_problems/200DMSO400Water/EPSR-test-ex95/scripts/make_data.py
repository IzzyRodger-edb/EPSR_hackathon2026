#!/usr/bin/env python3
"""Create the nrtype 1, 3 and 4 data files of EPSR-test-ex95 from the
GUDRUN (nrtype 5) .mint01 files of the base test case.

  D6DMSO_D2O.mint01 -> D6DMSO_D2O.nrt1   nrtype 1: "Q S(Q)", no header, point Q
  D6DMSO_HDO.mint01 -> D6DMSO_HDO.nrt3   nrtype 3: GENIE format(i7,3g14.7), 4 header
                                          lines, histogram (bin-boundary) Q as in .mint01
  D6DMSO_H2O.mint01 -> D6DMSO_H2O.nrt4   nrtype 4: GENIE format, 4 header lines,
                                          point Q, plus a rubbish last point (ignored)

Point-mode Q values are the bin centres 0.5*(Q(i)+Q(i+1)) of the histogram Q.
usage (in the case folder): python3 scripts/make_data.py
"""

import sys


def read_mint(fn):
    rows = []
    for line in open(fn):
        if "#" in line or not line.strip():
            continue
        q, s, e = (float(x) for x in line.split()[:3])
        rows.append((q, s, e))
    return rows


def e14(x):
    return f"{x:14.6E}"


def centres(rows):
    return [
        (0.5 * (rows[i][0] + rows[i + 1][0]), rows[i][1], rows[i][2])
        for i in range(len(rows) - 1)
    ]


def genie(fn, title, rows):
    with open(fn, "w") as f:
        f.write(f" {title}\n")
        f.write(" GENIE SHOW DATA listing\n")
        f.write(f" {len(rows)} points\n")
        f.write("      N             Q          S(Q)          E(Q)\n")
        for i, (q, s, e) in enumerate(rows, 1):
            f.write(f"{i:7d}{e14(q)}{e14(s)}{e14(e)}\n")


d2o = read_mint("D6DMSO_D2O.mint01")
with open("D6DMSO_D2O.nrt1", "w") as f:
    for q, s, e in centres(d2o):
        f.write(f" {q:.7E} {s:.7E}\n")
hdo = read_mint("D6DMSO_HDO.mint01")
genie("D6DMSO_HDO.nrt3", "D6DMSO HDO histogram data", hdo)
h2o = centres(read_mint("D6DMSO_H2O.mint01"))
h2o.append((h2o[-1][0] + 0.05, 9.99, 0.0))  # rubbish last point, ignored by EPSR
genie("D6DMSO_H2O.nrt4", "D6DMSO H2O point data", h2o)
