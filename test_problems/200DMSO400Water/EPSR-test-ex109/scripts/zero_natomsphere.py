# Set the coarse-graining 'natomsphere' value (6th number of each atom-type
# parameter line) to 0 for every atom type of nosph.ato.
L = open("nosph.ato").read().split("\n")
n = 0
for i in range(len(L) - 30, len(L)):
    w = L[i].split()
    if len(w) >= 6 and w[5] == "0.100000E+01":
        L[i] = L[i].replace("0.100000E+01", "0.000000E+00")
        n += 1
print(n, "atom types changed")
open("nosph.ato", "w").write("\n".join(L))
