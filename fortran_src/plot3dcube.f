      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_plot3dcube()

      stop
      end      

      subroutine run_plot3dcube()
c
c invokes the system command to perform 3d plotting surfaces using Jmol
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      character*256 plot_output
      integer*4 ierrout
c
c setup the plot3d symbols
c
      call setup_symbols_plot3dcube()
c
c now setup the input file
c
      nin=11
      nout=6
      call setup_input_file(nin,nout,ierrout)
      close(nin)
c
c if an error was encountered on setup, exit now
c
      if(ierrout.ne.0) return

      write(6,102) prompt(1:ilenprompt),fname_inp(1:ilenfname_inp)
102   format(/a,1x,'Will use plot3dcube input file: "',a,'"')

		call map3dcube(fname_inp,plot_output)

c
c now generate the .txt file from which plot3d will get its parameters
c
      open(11,file=directory_name(1:nchardir)//'plot.txt'
     *,status='unknown',iostat=ierr)
      if(ierr.eq.0) then
         write(11,103) plot_output(1:len_trim(plot_output))
103   format(a)
         close(11)
      else
         write(6,104)
104   format(/'run_plot3dcube> Not able to open "plot.txt".')
         return
      endif
c
c now generate the command that will actually run this file
c
      plot_output=system_cd(1:ilensystem_cd)
     *//' '
     *//system_quotes(1:ilensystem_quotes)
     *//directory_name(1:nchardir)
     *//system_quotes(1:ilensystem_quotes)
     *//' '//system_join(1:ilensystem_join)
     *//' '//system_binaries(1:ilensystem_binaries)//system_jmol(1:ilensystem_jmol)
     *//' -s plot.txt'
      write(6,*) plot_output
      call system(plot_output,ierr)
      write(6,*) ierr
      return
      end






