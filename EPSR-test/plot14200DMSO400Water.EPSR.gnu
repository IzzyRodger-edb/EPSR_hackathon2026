#reset
set title 'C:\EPSR26\EPSR\WorkedExamples\1DMSO2Water\200DMSO400Water'
set xlabel 'Iteration number'
set ylabel 'U [kJ/mole]'
set autoscale x
set autoscale y
set nologscale x
set nologscale y
unset grid
set style line   1 lt   1 pt   8 ps    0.30
plot \
'C:\EPSR26\EPSR\WorkedExamples\1DMSO2Water\200DMSO400Water.EPSR.erg' u   0:( 0.10000E+01*(column(  1)+ 0.000E+00)) title " Block           1" w lines ls   1
