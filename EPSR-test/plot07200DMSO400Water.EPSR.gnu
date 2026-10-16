#reset
set title 'C:\EPSR26\EPSR\WorkedExamples\1DMSO2Water\200DMSO400Water'
set xlabel 'Q [1/A]'
set ylabel 'F(Q)'
set autoscale x
set autoscale y
set nologscale x
set nologscale y
unset grid
set label   1 'D6DMSO_D2O                   ' at  0.100E+02, 0.100E+00
set style line   1 lt   1 pt   8 ps    0.30
set label   2 'D6DMSO_HDO                   ' at  0.100E+02, 0.110E+01
set style line   2 lt   2 pt   8 ps    0.30
set label   3 'D6DMSO_H2O                   ' at  0.100E+02, 0.210E+01
set style line   3 lt   3 pt   8 ps    0.30
plot \
'C:\EPSR26\EPSR\WorkedExamples\1DMSO2Water\200DMSO400Water.EPSR.u01' u   1:( 0.10000E+01*(column(  2)+ 0.000E+00)) notitle w lines ls   1, \
'C:\EPSR26\EPSR\WorkedExamples\1DMSO2Water\200DMSO400Water.EPSR.t01' u   1:(abs(column(  2))>=1.0e-5? 0.10000E+01*column(  2)+ 0.000E+00:1/0) notitle w p ls   1, \
'C:\EPSR26\EPSR\WorkedExamples\1DMSO2Water\200DMSO400Water.EPSR.u01' u   1:( 0.10000E+01*(column(  4)+ 0.100E+01)) notitle w lines ls   2, \
'C:\EPSR26\EPSR\WorkedExamples\1DMSO2Water\200DMSO400Water.EPSR.t01' u   1:(abs(column(  4))>=1.0e-5? 0.10000E+01*column(  4)+ 0.100E+01:1/0) notitle w p ls   2, \
'C:\EPSR26\EPSR\WorkedExamples\1DMSO2Water\200DMSO400Water.EPSR.u01' u   1:( 0.10000E+01*(column(  6)+ 0.200E+01)) notitle w lines ls   3, \
'C:\EPSR26\EPSR\WorkedExamples\1DMSO2Water\200DMSO400Water.EPSR.t01' u   1:(abs(column(  6))>=1.0e-5? 0.10000E+01*column(  6)+ 0.200E+01:1/0) notitle w p ls   3
