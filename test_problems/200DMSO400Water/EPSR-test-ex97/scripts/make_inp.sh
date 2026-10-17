# run in the case folder, starting from a copy of EPSR-test/200DMSO400Water.EPSR.inp
python3 -I scripts/setkw.py 200DMSO400Water.EPSR.inp "porodpwr=0.5 3.5" 1:porodadd=0.25 2:porodadd=0.10 3:porodadd=-0.5
