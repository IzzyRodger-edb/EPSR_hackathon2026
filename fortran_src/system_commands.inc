c
c system command to list a directory
c
      character*256 system_ls
      integer*4 ilensystem_ls
c
c system command to change a directory
c
      character*256 system_cd
      integer*4 ilensystem_cd
c
c system command to make a new directory
c
      character*256 system_md
      integer*4 ilensystem_md
c
c system command to edit a file
c
      character*256 system_edit
      integer*4 ilensystem_edit
c
c system command to delete a file or files
c
      character*256 system_del
      integer*4 ilensystem_del
c
c character sequence to signal the end of a directory name
c
      character*256 system_termination
      integer*4 ilensystem_termination
c
c system characters to join two commands in one line
c
      character*256 system_join
      integer*4 ilensystem_join
c
c system command to invoke GNUplot
c
      character*256 system_gnu
      integer*4 ilensystem_gnu
c
c system command to invoke Jmol
c
      character*256 system_jmol
      integer*4 ilensystem_jmol

c Folder where binaries are stored

      character*256 system_binaries
      integer*4 ilensystem_binaries

c Extension to be added to binaries to give executable

      character*256 system_bin_ext
      integer*4 ilensystem_bin_ext
c
c character sequence for command to do 2D and 3D plotting
c
      character*256 system_plot2d,system_plot3d,system_plot3djmol
      integer*4 ilensystem_plot2d,ilensystem_plot3d,ilensystem_plot3djmol
c
c character sequence to plot the .ato file
c
      character*256 system_plotato
      integer*4 ilensystem_plotato
c
c type of output for plot2d and plot3d (e.g. /gif, /cps, etc.)
c
      character*256 system_pgout
      integer*4 ilensystem_pgout
c
c to put quotes around filenames if required
c
      character*256 system_quotes
      integer*4 ilensystem_quotes
c
c to run Mopac
c
      character*256 system_mopac
      integer ilensystem_mopac
c
c to Mopac number of orbitals - this is set when Mopac is compiled.
c
      character*256 system_mopac_norbitals
      integer ilensystem_mopac_norbitals
c
c to set character type for data points in GNUplot
c
      integer*4 gnu_pointtype
c
c to set character size for data points in GNUplot
c
      real*4 gnu_pointsize

 
      common/system_commands_c/system_ls,system_cd
     *,system_termination,system_join,system_gnu,system_jmol
     *,system_binaries,system_bin_ext
     *,system_plot2d,system_plot3d,system_plot3djmol,system_pgout
     *,system_md,system_edit,system_del,system_plotato
     *,system_quotes,system_mopac,system_mopac_norbitals

      common/system_commands_i/ilensystem_ls,ilensystem_cd
     *,ilensystem_termination,ilensystem_join,ilensystem_gnu,ilensystem_jmol
     *,ilensystem_binaries,ilensystem_bin_ext
     *,ilensystem_plot2d,ilensystem_plot3d,ilensystem_plot3djmol,ilensystem_pgout
     *,ilensystem_md,ilensystem_edit,ilensystem_del,ilensystem_plotato
     *,ilensystem_quotes,ilensystem_mopac,ilensystem_mopac_norbitals
     *,gnu_pointtype

      common/system_commands_r/gnu_pointsize
