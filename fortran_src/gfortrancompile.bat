@echo off
set cflags=-c -fbounds-check -fno-backslash -Wall -ffixed-line-length-none -Wno-unused-label -O2
gfortran %cflags% %1
