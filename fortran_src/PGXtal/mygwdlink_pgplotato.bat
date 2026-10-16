@echo off

set F90_PROJ=/compile_only /nologo /warn:nofileopt
set LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /incremental:no /machine:I386 

rem df %F90_PROJ% %1 %2 %3 %4 %5 %6 %7 %8 %9 

link %LINK32_FLAGS% %1.obj cosdsind.obj wr_int.obj randomno.obj PGXtal/myplot3d_pgplotato.obj PGXtal/mypgxtal24oct2002.obj PGXtal/pgcell.obj PGXtal/dsqinf.obj %PGPLOT_DIR%\PGPLOT.lib %PGPLOT_DIR%\GrWin.lib 

rem del if exist %1.obj del %1.obj

