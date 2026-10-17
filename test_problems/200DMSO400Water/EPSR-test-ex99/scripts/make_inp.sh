# run in the case folder, starting from a copy of EPSR-test/200DMSO400Water.EPSR.inp
python3 -I scripts/setkw.py 200DMSO400Water.EPSR.inp "qmin=0.60 0.80" "fwhm=0.03" "rspcrmin=1.50" "rspcfrac=0.3" "qinscalefc=0.98 0.5" 1:rshmin=0.5 2:rshmin=0.9 3:rshmin=1.1 1:tweak=1.03 3:tweak=0.97 2:efilereq=0.5 3:efilereq=0.0
