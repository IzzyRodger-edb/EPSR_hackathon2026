      program setup_standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call setup()

      stop
      stop
      end      

c     2007-09-30. Modified read_pcof and write_pcof to include Gaussian values
c     as well as minimum allowed distances
      
      subroutine setup()
c
c sets up and modifies the various input files for EPSR and associated programs 
c
      use String_Utility
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
         write(6,100) 
100   format(/1x,'What do you want to setup? Options are:-'//
     *1x,'epsr  <filename>  - to create/edit EPSR input filess'/
     *1x,'epsr17 <filename> - to read data from EPSR17 input files'/
     *1x,'epsr16 <filename> - to read data from EPSR16 input files'/
     *1x,'partials <filename> - to create/edit partials input files'/
     *1x,'bragginten <filename> - to create/edit bragginten input files'/
     *1x,'clusters <filename> - to create/edit clusters input files'/
     *1x,'voids <filename> - to create/edit voids input files'/
     *1x,'coord <filename> - to create/edit coord input files'/
     *1x,'chains <filename> - to create/edit chains input files'/
     *1x,'rings <filename> - to create/edit rings input files'/
     *1x,'triangles <filename> - to create/edit triangles input files'/
     *1x,'torangles <filename> - to create/edit torangles input files'/
     *1x,'sharm <filename> - to create etc. SHARM input files'/
     *1x,'sdf <filename> - to create etc. SDF input files'/
     *1x,'fluctuations <filename> - to create fluctuations input files'/
     *1x,'wtsn <filename> - to create neutron weights input files'/
     *1x,'wtsx <filename> - to create x-ray weights input files'/
     *1x,'q - to quit setup/
     *1x,'plot2d <filename> - to create/edit input file for 2D'
     *,' plotting'/
     *1x,'plot3d <filename> - to create/edit input file for 3D'
     *,' plotting'/
     *)
         nwords=0
         do while(nwords.eq.0)
            sometext=' '
            write(6,114)
            read(5,101,iostat=ierr) sometext
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         end do
         nfirst=1
      else
         nfirst=2
      endif
101   format(a)
114   format(/'setup> ',$)
      command=sometext(ncf(nfirst):ncl(nfirst)) !Make sure the symbol is always lower case - to be recognised
      command=StrLowCase(command)
      symbol_setup=command
      ilensetup=index(symbol_setup,' ')-1
c
c first setup the symbols
c
      if(command.eq.'epsr') then

         call setup_symbols_epsr()

      else if(command.eq.'partials') then
c
c set up the partials symbol tables and default values
c
         call setup_symbols_partials()

      else if(command.eq.'bragginten') then
c
c set up the clusters symbol tables and default values
c
         call setup_symbols_bragginten()

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

      else if(command.eq.'wtsn') then
c
c set up the voids symbol tables and default values
c
         call setup_symbols_wtsn()

      else if(command.eq.'wtsx') then
c
c set up the voids symbol tables and default values
c
         call setup_symbols_wtsx()

      else if(command.eq.'q') then
c
c q means quit without doing anything.
c
         return

      else

         write(6,102)
102   format(/'setup> Command not recognised')
         return

      endif
c
c now setup the input file
c
      nin=11
      nout=6
      if(command.ne.'epsr16'.and.command.ne.'epsr17') then
         call setup_input_file(nin,nout,ierrout)
      endif

c
c ierrout .lt. 0 signifies the file could not be opened or on reading encountered an 
c error, but an input file has been successfully defined using the default values,
c so these values are now available for viewing and changing.
c

      if(command.eq.'plot2d') then

         call read_plot2d(nin,nout,ierrout)
         if(ierrout.le.0) call process_plot_symbols(5,6,ierrout,0)

      else if(command.eq.'plot3d') then

         call read_plot3d(nin,nout,ierrout)
         if(ierrout.le.0) call process_plot_symbols(5,6,ierrout,0)

      else

         call read_input_file(nin,nout,ierrout)
         if(ierrout.le.0) call process_symbols(5,6,ierrout,0)

      endif


      return
      end


