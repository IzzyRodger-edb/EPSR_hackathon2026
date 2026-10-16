      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_Jmol()

      end      

      subroutine run_Jmol()

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'plot_pars.inc'
      include 'system_commands.inc'
      include 'plotato.inc'

      character*256 plot_command
c
c now generate the command that will actually run this file
c
      plot_command=system_cd(1:ilensystem_cd)
     *//' '
     *//system_quotes(1:ilensystem_quotes)
     *//directory_name(1:nchardir)
     *//system_quotes(1:ilensystem_quotes)
     *//' '//system_join(1:ilensystem_join)
     *//' '//system_jmol(1:ilensystem_jmol)
      write(6,*) plot_command
      call system(plot_command,ierr)
      write(6,*) ierr

      return
      end

