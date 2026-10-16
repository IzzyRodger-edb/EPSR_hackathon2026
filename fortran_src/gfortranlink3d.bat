@echo off

set cflags=-static -O

gfortran %cflags% -o %EPSRbin%\%1.exe %1.o cosdsind.o wr_int.o randomno.o auxiliary_routines.o myplot3d.o mypgxtal24oct2002.o pgcell.o dsqinf.o -Wl,--subsystem,console -LC:\Users\aks45\MinGW\lib %EPSRbin%\PGPLOT\PGPLOTlib\libpgplot.a -lmingw32 -lgcc -lmsvcrt -luser32 -lkernel32 -lgdi32 

