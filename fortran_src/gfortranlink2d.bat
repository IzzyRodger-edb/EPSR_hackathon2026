@echo off

set cflags=-static -O

gfortran %cflags% -o %EPSRbin%\%1.exe %1.o cosdsind.o wr_int.o randomno.o auxiliary_routines.o myplot2db3.o mypgxtal24oct2002.o pgcell.o dsqinf.o -Wl,--subsystem,console -LC:\Users\aks45\MinGW\lib %EPSRbin%\PGPLOT\PGPLOTlib\libpgplot.a -lmingw32 -lgcc -lmsvcrt -luser32 -lkernel32 -lgdi32 

rem gfortran %cflags% -o %EPSRbin%\%1.exe %1.o cosdsind.o wr_int.o randomno.o auxiliary_routines.o myplot2db3.o mypgxtal24oct2002.o pgcell.o dsqinf.o addgif.o -Wl,--subsystem,console -LC:\Users\aks45\MinGW\lib -L%EPSRbin%\PGPLOT\PGPLOTlib -lcpgplot -lpgplot -lGrWin -mwindows
