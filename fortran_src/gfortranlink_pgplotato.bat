@echo off
set cflags=-static -O
gfortran %cflags% -o %EPSRbin%\%1.exe %1.o cosdsind.o wr_int.o randomno.o auxiliary_routines.o myplot3d_pgplotato.o mypgxtal24oct2002.o pgcell.o dsqinf.o addgif.o -Wl,--subsystem,console -LC:\MinGW\lib -lcpgplot -lpgplot -lGrWin -mwindows
