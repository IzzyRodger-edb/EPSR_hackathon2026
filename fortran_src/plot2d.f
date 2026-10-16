      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_plot2d()

      end      

      subroutine run_plot2d()
c
c invokes the system command to perform 3d plotting
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      character*256 plot_output
      integer*4 ierrout
c
c setup the plot2d symbols
c
      call setup_symbols_plot2d()
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
102   format(/a,1x,'Will use plot2d input file: "',a,'"')

c
c now generate the .txt file from which plot2d will get its parameters
c
      open(11,file=directory_name(1:nchardir)//'plot.txt'
     *,status='unknown',iostat=ierr)
      if(ierr.eq.0) then
         write(11,103) fname_inp(1:ilenfname_inp)
103   format(a)
         write(11,103) system_pgout(1:ilensystem_pgout)
         close(11)
      else
         write(6,104)
104   format(/'run_plot2d> Not able to open "plot.txt".')
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
     *//' '//system_binaries(1:ilensystem_binaries)//system_plot2d(1:ilensystem_plot2d)
     *//' < plot.txt'
      write(6,*) plot_output
      call system(plot_output,ierr)
      write(6,*) ierr
      return
      end






