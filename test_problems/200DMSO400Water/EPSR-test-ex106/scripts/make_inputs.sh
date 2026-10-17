# Run in the case folder after copying the 16 base files of EPSR-test.
sed -i 's/^control-p   none 0 0 /control-p   p 0.0 0.01 /; s/^ntimes      5 /ntimes      2 /; s/^niter       1 /niter       2 /' 200DMSO400Water.EPSR.inp
python3 - <<'P'
L = open('200DMSO400Water.pcof').read().split('\n')
L[67] = '   2.60000000       0.50000000       0.00000000    '   # O-O1:  rminex pcoefex coordex
L[82] = '   2.70000000       0.00000000       0.00000000    '   # O1-O1
open('200DMSO400Water.pcof', 'w').write('\n'.join(L))
P
