#reset
set title 'C:\EPSR26\EPSR\WorkedExamples\1DMSO2Water\200DMSO400Water'
set xlabel 'r [A]'
set ylabel 'f(r)'
set autoscale x
set autoscale y
set nologscale x
set nologscale y
unset grid
set style line   1 lt   1 pt   8 ps    0.30
set style line   2 lt   2 pt   8 ps    0.30
set style line   3 lt   3 pt   8 ps    0.30
plot \
'C:\EPSR26\EPSR\WorkedExamples\1DMSO2Water\200DMSO400Water.EPSR.x01' u   1:( 0.10000E+01*(column(  2)+ 0.000E+00)) title "D6DMSO_D2O" w lines ls   1, \
'C:\EPSR26\EPSR\WorkedExamples\1DMSO2Water\200DMSO400Water.EPSR.w01' u   1:(abs(column(  2))>=1.0e-5? 0.10000E+01*column(  2)+ 0.000E+00:1/0) notitle w p ls   1, \
'C:\EPSR26\EPSR\WorkedExamples\1DMSO2Water\200DMSO400Water.EPSR.x01' u   1:( 0.10000E+01*(column(  4)+ 0.200E+00)) title "D6DMSO_HDO" w lines ls   2, \
'C:\EPSR26\EPSR\WorkedExamples\1DMSO2Water\200DMSO400Water.EPSR.w01' u   1:(abs(column(  4))>=1.0e-5? 0.10000E+01*column(  4)+ 0.200E+00:1/0) notitle w p ls   2, \
'C:\EPSR26\EPSR\WorkedExamples\1DMSO2Water\200DMSO400Water.EPSR.x01' u   1:( 0.10000E+01*(column(  6)+ 0.400E+00)) title "D6DMSO_H2O" w lines ls   3, \
'C:\EPSR26\EPSR\WorkedExamples\1DMSO2Water\200DMSO400Water.EPSR.w01' u   1:(abs(column(  6))>=1.0e-5? 0.10000E+01*column(  6)+ 0.400E+00:1/0) notitle w p ls   3
