# run in the case folder, starting from copies of the EPSR-test files
python3 -I scripts/make_data.py
python3 -I scripts/setkw.py 200DMSO400Water.EPSR.inp 1:datafile=D6DMSO_D2O.nrt1 1:nrtype=1 2:datafile=D6DMSO_HDO.nrt3 2:nrtype=3 3:datafile=D6DMSO_H2O.nrt4 3:nrtype=4
