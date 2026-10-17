# Run in the case folder after copying the 16 base files of EPSR-test.
sed -i 's/^control-p   none 0 0 /control-p   bc 0.0 0.01 /; s/^ref_intra   0.000    0.000 /ref_intra   0.500    0.300 /; s/^ntimes      5 /ntimes      2 /; s/^niter       1 /niter       2 /' 200DMSO400Water.EPSR.inp
