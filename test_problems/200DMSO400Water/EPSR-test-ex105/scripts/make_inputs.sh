# Run in the case folder after copying the 16 base files of EPSR-test.
sed -i 's/^control-p   none 0 0 /control-p   abc 0.0 0.01 /; s/^ntimes      5 /ntimes      2 /; s/^niter       1 /niter       2 /' 200DMSO400Water.EPSR.inp
python3 scripts/clone_run.py boxa "control-p a 0.0 0.01"
