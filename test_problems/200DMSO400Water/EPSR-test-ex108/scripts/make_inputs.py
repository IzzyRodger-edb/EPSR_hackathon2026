# Generates the ex108 inputs from the base EPSR-test files (run in the case folder):
#  * .pcof: drop the pairs S-H1 and O1-H1 (new pairs must be generated) and the
#    'ncoeffp' header line (so the value given in the .EPSR.inp is used to read the
#    coefficients);
#  * .EPSR.inp: add 'ncoeffp 60' before fnameato/fnamepcof, ntimes 2, niter 2.
import re

L = open("200DMSO400Water.pcof").read().split("\n")
out = []
i = 0
while i < len(L):
    l = L[i]
    if l.startswith("ncoeffp "):
        i += 1
        continue
    if re.match(r"^(S   H1|O1  H1) ", l):
        i += 3
        continue
    if l.strip() == "21":
        l = l.replace("21", "19")
    out.append(l)
    i += 1
open("200DMSO400Water.pcof", "w").write("\n".join(out))
L = open("200DMSO400Water.EPSR.inp").read().split("\n")
out = []
for l in L:
    if l.startswith("fnameato "):
        out.append(
            "ncoeffp     60               Number of coefficients used to define the EP."
        )
    l = re.sub(r"^ntimes      5 ", "ntimes      2 ", l)
    l = re.sub(r"^niter       1 ", "niter       2 ", l)
    out.append(l)
open("200DMSO400Water.EPSR.inp", "w").write("\n".join(out))
