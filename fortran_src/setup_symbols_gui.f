      program setup_standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call setup()

      stop
      end      

c     2007-09-30. Modified read_pcof and write_pcof to include Gaussian values
c     as well as minimum allowed distances

c 2011-12-9. This version prepares a file of values to be read by EPSRshell GUI.
      
      subroutine setup()
c
c sets up and modifies the various input files for EPSR and associated programs 
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      character*20 command
c
c first determine the type of setup required
c
      if(nwords.lt.2) then
         write(6,114)
114   format(/'setup_symbols_gui> You need to specify a command to set up the symbols for')
         return
      else
         nfirst=2
      endif

      command=sometext(ncf(nfirst):ncl(nfirst))
      symbol_setup=command
      ilensetup=index(symbol_setup,' ')-1
c
c first setup the symbols
c
      if(command.eq.'epsr') then

         call setup_symbols_epsr()

      else if(command.eq.'rmcsr') then

         call setup_symbols_rmc()

      else if(command.eq.'partials') then
c
c set up the partials symbol tables and default values
c
         call setup_symbols_partials()

      else if(command.eq.'clusters') then
c
c set up the clusters symbol tables and default values
c
         call setup_symbols_clusters()

      else if(command.eq.'coord') then
c
c set up the coord symbol tables and default values
c
         call setup_symbols_coord()

      else if(command.eq.'chains') then
c
c set up the coord symbol tables and default values
c
         call setup_symbols_chains()

      else if(command.eq.'rings') then
c
c set up the coord symbol tables and default values
c
         call setup_symbols_rings()

      else if(command.eq.'sharm') then
c
c set up the SHARM symbol tables and default values
c
         call setup_symbols_SHARM()

      else if(command.eq.'sdf') then
c
c set up the SDF symbol tables and default values
c
         call setup_symbols_SDF()

      else if(command.eq.'plot2d') then
c
c set up the plot2d symbol tables and default values
c
         call setup_symbols_plot2d()

      else if(command.eq.'plot3d') then
c
c set up the plot3d symbol tables and default values
c
         call setup_symbols_plot3d()

      else if(command.eq.'triangles') then
c
c set up the triangles symbol tables and default values
c
         call setup_symbols_triangles()

      else if(command.eq.'torangles') then
c
c set up the torangles symbol tables and default values
c
         call setup_symbols_torangles()

      else if(command.eq.'voids') then
c
c set up the voids symbol tables and default values
c
         call setup_symbols_voids()

      else if(command.eq.'fluctuations') then
c
c set up the voids symbol tables and default values
c
         call setup_symbols_fluctuations()

      else if(command.eq.'nwts'.or.command.eq.'epsrwts') then
c
c set up the voids symbol tables and default values
c
         call setup_symbols_nwts()

      else if(command.eq.'xwts') then
c
c set up the voids symbol tables and default values
c
         call setup_symbols_xwts()

      else if(command.eq.'q') then
c
c q means quit without doing anything.
c
         return

      else

         write(6,1000)
1000  format(/'setup_symbols_generic> Command not recognised')
         return

      endif

      nout=12
      open(nout,file='symbols_for_GUI.txt',status='unknown',iostat=ierr)
c Because gfortran gives an error when a file opened with status unknown does not exist!
      ierr=0
      call save_symbols(nout)  
      close(nout)      

      return
      end

      subroutine save_symbols(nout)
c
c writes the general input file
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      include 'iresetnumbers.inc'
      integer*4 iqmin

c
c find the symbol reference for qmin
c
      if(symbol_setup.eq.'epsr') then

         iqmin=ireff(iqminref,nfiles)

      else
c
c ensure this value will not be used if it has not been defined
c
         iqmin=nsymbol

      endif

c Write the dataprompt and input file extension

      write(nout,'(1x,a,1x,a)') dataprompt(1:ilendataprompt),extinp(1:ilenextinp)

      icount=0
      ifilecount=1
      do while (icount.lt.nsymboldata)
            icount=icount+1
c
c for EPSR files skip over the bit related to the .pcof file - this is written later
c
            if(symbol_setup.eq.'epsr'.and.icount.gt.iqmin
     *.and.icount.lt.nsymboluptodata) then
               icount=nsymboluptodata
            endif

            iref=ireff(icount,ifilecount)
            sometext=symbol_default(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            write(nout,101,iostat=ierr) symbol(icount)
     *,sometext(ncf(1):ncl(max(1,nwords))),symbol_type(icount),symbol_text(icount)
101   format(a11,1x,a,1x,a,15x,a)
            if(icount.ge.nsymboldata.and.ifilecount.lt.nfiles) then
               ifilecount=ifilecount+1
               icount=nsymboluptodata
            endif
      end do
c
c write a 'q' to signify the end of input
c
      write(nout,102,iostat=ierr)
102   format('q')
c
c if symbol_setup .eq. epsr then save the pcof symbols
c
      if(symbol_setup.eq.'epsr') then

         icount=iqmin+1
         do while (icount.lt.nsymboluptodata)
            iref=ireff(icount,ifilecount)
            sometext=symbol_default(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            write(nout,101) symbol(icount),sometext(ncf(1):ncl(max(1,nwords)))
     *,symbol_type(icount),symbol_text(icount)
            icount=icount+1
         end do
         write(nout,102)

      return
      end


