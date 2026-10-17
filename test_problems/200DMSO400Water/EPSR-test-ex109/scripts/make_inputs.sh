# Run in the case folder after copying the 16 base files of EPSR-test.
rm -f 200DMSO400Water.pcof 200DMSO400Water.EPSR.erg 200DMSO400Water.EPSR.g01 200DMSO400Water.EPSR.y01 \
      200DMSO400Water.EPSR.f01 200DMSO400Water.EPSR.s01 200DMSO400Water.EPSR.inpa
# .ato step-size line with 7 values (7th = tempdihreq)
sed -i '2s/.*/   0.00000E+00  0.10000E+00  0.30000E+00  0.30000E+00  0.10000E+01  0.65000E+02  0.20000E+02/' 200DMSO400Water.ato
sed -i 's/^ntimes      5 /ntimes      2 /; s/^niter       1 /niter       2 /' 200DMSO400Water.EPSR.inp
python3 - <<'P'
L = open('200DMSO400Water.EPSR.inp').read().split('\n')
i = [k for k, l in enumerate(L) if l.startswith('fnameato ')][0]
L.insert(i, 'ncoeffp     2000               Number of coefficients used to define the EP.')
open('200DMSO400Water.EPSR.inp', 'w').write('\n'.join(L))
P
python3 scripts/clone_run.py nosph
python3 scripts/zero_natomsphere.py
