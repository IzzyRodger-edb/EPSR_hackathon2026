c$$$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c$$$  RH: setup_symbols.f 
c$$$  for each of the analysis rountines sets up arrays/tables which are
c$$$  used to read the input files (the *.dat files).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setup_symbols_generic(auxprogname,ssg_ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c$$$  Routine to parse a .symbols file to allow external programs to use
c$$$  setup.
c$$$  
c$$$  In the .symbols file required parameters are:
c$$$  - extension   (the extension for the .dat file)
c$$$  - prompt      (the prompt used in the setup program)
c$$$  
c$$$  And if using the analysis routines which calculate for more than
c$$$  one set of input parameters (ndist > 0), then
c$$$  - dataprompt   (the generic name of these different parameter sets)
c$$$  - dataprompttext (a description of this for the  .dat file)
c$$$  are required.
c$$$
c$$$  Other than this all other symbols are declared by prefixing a line
c$$$  with the word "symbol" followed by the value for that symbol. Then
c$$$  on following lines need to specify default, type and text
c$$$  description for each symbol (these are indicated by the key words
c$$$  symbol_default, symbol_type, and symbol_text. There are three
c$$$  types a symbol can have: real -- specified by an r character,
c$$$  integer -- specified by an i, and character -- indicated by a c
c$$$  character.
c$$$
c$$$  Blank lines are ignored.
c$$$  
c$$$  Example .symbols file:
c$$$  
c$$$  extension .NDIST_EXAMPLE.dat
c$$$  prompt setup example>
c$$$  dataprompt ndist_example
c$$$  dataprompttext   ndist_example text, nice
c$$$  
c$$$  symbol fnameato
c$$$  symbol_default <undefined>
c$$$  symbol_type  c               
c$$$  symbol_text  A nice description
c$$$  
c$$$  symbol ndist
c$$$  symbol_default 1
c$$$  symbol_type i               
c$$$  symbol_text Number of distributions/variable sets to calculate over
c$$$  
c$$$  symbol atom1
c$$$  symbol_default XXX
c$$$  symbol_type c               
c$$$  symbol_text Atom1 label
c$$$  
c$$$  symbol atom2
c$$$  symbol_default XXX
c$$$  symbol_type c               
c$$$  symbol_text Atom2 label
c$$$  
c$$$  symbol rmin
c$$$  symbol_default 2.0
c$$$  symbol_type r               
c$$$  symbol_text Minimum distance to be considered bound                                         
c$$$  
c$$$  symbol rmax
c$$$  symbol_default 2.0
c$$$  symbol_type r               
c$$$  symbol_text Maximum distance to be considered bound
c$$$  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      include 'epsr_dimension.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'potential_par.inc'
      
      character*256 auxprogname,sometext2
      character*256 symbols_fname,val
      character*8 ext
      character*20 cmd
      integer*4 ilenext, ilen_sf,nin,lnum,iref,il,
     $     parsed_def, parsed_type, parsed_text,ssg_ierr
      
 102  format(/'setup_symbols_generic> ',a)            
 100  format(a)

      call parse(auxprogname,nwords,ncf,ncl,5,nstsize,nworddim)
c$$$  Need to parse up to space... 
      write (6,102) " command passed in: "//auxprogname(ncf(1):ncl(1))
      write (6,102) " command passed in: "//auxprogname
      

c$$$  extention for template file name
      ext=".symbols"
      ilenext=8
      
c$$$  Get file name to read:
      symbols_fname = auxprogname(ncf(1):ncl(1)) // ext
      ilen_sf = ncl(1) - ncf(1) + ilenext +1

c$$$  initialise variables which need to be checked that they have been
c$$$  parsed:
      extinp=' '
      prompt=' '
      dataprompt=' '
      dataprompttext=' '
      nsymbol=0
      nsymboldata=0
      nsymboluptodata=0

c     open this file if possible
      nin=14
      open(nin,file=symbols_fname(1:ilen_sf),status='old',
     $     iostat=ierr)
      if (ierr.eq.0) then
         write(6,102) "Processing "//symbols_fname(1:ilen_sf)
c$$$  process the .symbols file now
         ierr=0
         lnum = 0
         iref = 0
         do while(ierr.eq.0)
            nwords=0
c$$$  Have to use sometext2 rather than sometext as setup is using
c$$$  sometext
            sometext2 = ' '
            read(nin,100,iostat=ierr) sometext2 ! read next line
            call parse(sometext2,nwords,ncf,ncl,12,nstsize,nworddim)
c$$$  If there are words then process them; ignore blank lines:
            cmd = " "
            if(nwords.gt.0) then
c$$$  cmd is first word on each line; val is the rest of the line
               cmd = sometext2(ncf(1):ncl(1))
               val = sometext2(ncf(2):ncl(nwords))
               if (cmd.eq."extension") then
                  write (6,*) "Found extension, extension = ", val
                  extinp=val
                  ilenextinp=index(extinp,' ')-1
                  
               elseif (cmd.eq."prompt") then
                  write (6,*) "Found prompt, prompt = ", val
                  prompt=val
                  ilenprompt = -1
                  ilenprompt=index(prompt,'>')
                  if (ilenprompt.le.0) then
                     write(6,102) "Prompt must end in a > character" 
                     write(6,102) "Exiting" 
                     ssg_ierr = -1
                     return
                  end if
                  
               elseif (cmd.eq."dataprompt") then
                  write (6,*) "Found dataprompt, dataprompt = ", val
                  dataprompt=val
                  ilendataprompt=index(dataprompt,' ')-1
                  
               elseif (cmd.eq."dataprompttext") then
                  write (6,*) "Found dataprompttext, dataprompttext = "
     $                 , val
                  dataprompttext=val
                  
               elseif (cmd.eq."symbol") then
                  write (6,*) "Found new symbol, symbol = ", val
                  iref = iref + 1 ! found another new symbol 
                  symbol(iref) = val
                  iflag(iref)=0
                  
                  if (val.eq."ndist") then
                     nsymboluptodata = iref
                     write (6,102) "Found ndist symbol, therefore all"//
     $                    "symbols beyond this are part of each "
     $                    //"distribution."
                  end if
                  
c$$$  Check that get default, type, and text -- so read next 3 lines:
                  il = 0
c$$$  variables to check that parsed def:
                  parsed_def = -1 
                  parsed_type = -1
                  parsed_text = -1
                  do while((il.lt.3).and.(ierr.eq.0))
                     sometext2 = ' '
                     nwords = 0
                     read(nin,100,iostat=ierr) sometext2 ! read next line
                     call parse(sometext2,nwords,ncf,ncl,12,nstsize,
     $                    nworddim)
                     if (nwords.ge.0) then ! ignore blank lines
                        cmd = sometext2(ncf(1):ncl(1))
                        val = sometext2(ncf(2):ncl(nwords))

                        if (cmd.eq."symbol_default") then
                           parsed_def = 1
                           symbol_default(iref)=val
                        elseif (cmd.eq."symbol_type") then
                           parsed_type = 1
                           symbol_type(iref)=val
                        elseif (cmd.eq."symbol_text") then
                           parsed_text = 1
                           symbol_text(iref)=val
                        else
                           write(6,102) " "
                           write(6,*) "Unrecognised value for symbol ",
     $                          symbol(iref), " : ", cmd, " in line ",
     $                          lnum
                           write(6,102) "Exiting setup_symbol_generic"
                           ssg_ierr = -1
                           return
                        end if

                        il = il + 1
                     endif
                     lnum = lnum + 1
                  end do        ! end of loop thru' def, type and text lines
                  if (parsed_def.le.0) then
                     write(6,102) " "
                     write(6,*) "No default found for symbol: ",
     $                    symbol(iref)
                     write(6,102) "Exiting setup_symbol_generic"
                     ssg_ierr = -1
                     return
                  end if
                  if (parsed_type.le.0) then
                     write(6,102) " "
                     write(6,*) "No type found for symbol: ",
     $                    symbol(iref)
                     write(6,102) "Exiting setup_symbol_generic"
                     ssg_ierr = -1
                     return
                  end if
                  if (parsed_text.le.0) then
                     write(6,102) " "
                     write(6,*) "No text found for symbol: ",
     $                    symbol(iref)
                     write(6,102) "Exiting setup_symbol_generic"
                     ssg_ierr = -1
                     return
                  end if
                  
               endif
            endif
            lnum = lnum +1
         end do
         
         write(6,102) " "
         write(6,*) "Parsed ",lnum," lines in "//
     $        symbols_fname(1:ilen_sf) 
         write(6,*) "Found ", iref, " symbols "
         
         close(nin)


      else
c$$$  If there is a problem we exit setup:
         write(6,102)"Problem reading file: "
     $        //symbols_fname(1:ilen_sf)
         write(6,102) "Exiting setup_symbol_generic"
         ssg_ierr = -1
         return
      endif

c$$$  Check that have parsed the important variables: prompt,
c$$$  extension...
      if (extinp.eq." ") then
         write(6,102) " "
         write(6,*) "No extension specified in file : ",
     $        symbols_fname(1:ilen_sf)
         write(6,102) "Extension is required."
         write(6,102) "Exiting setup_symbol_generic"
         ssg_ierr = -1

         return
      end if
      if (prompt.eq." ") then
         write(6,102) " "
         write(6,*) "No prompt specified in file : ",
     $        symbols_fname(1:ilen_sf)
         write(6,102) "Prompt is required."
         write(6,102) "Exiting setup_symbol_generic"
         ssg_ierr = -1
         return
      end if

      if ((dataprompt.eq." ").and.(nsymboluptodata.ne.0)) then
         write(6,102) " "
         write(6,*) "No dataprompt specified in file : ",
     $        symbols_fname(1:ilen_sf)
         write(6,102) "If have set ndist then dataprompt is required."
         write(6,102) "Exiting setup_symbol_generic"
         ssg_ierr = -1
         return
      end if

      if ((dataprompttext.eq." ").and.(nsymboluptodata.ne.0)) then
         write(6,102) " "
         write(6,*) "No dataprompttext specified in file : ",
     $        symbols_fname(1:ilen_sf)
         write(6,102) "If have set ndist then dataprompttext"
     $        //" is required."
         write(6,102) "Exiting setup_symbol_generic"
         ssg_ierr = -1
         return
      end if



c$$$  nsymboldata specifies the end of the file -- the number of symbols
c$$$  to read from the file
      nsymboldata=iref 

c$$$  If not set any distributions set nsymboluptodata to be beyond the
c$$$  number of symbols nsymbol:
      if (nsymboluptodata.eq.0) then
         nsymboluptodata=iref + 1
         write(6,102) "ndist not found, therefore no
     $distributions specified"
      end if
      nsymbolplus=nsymboldata-nsymboluptodata
      nsymbol=iref

c$$$  Initialise symbol values and print symbols:
      do is=1,nsymbol
         symbol_value(is)=symbol_default(is)
         write(6,107) is,symbol(is),symbol_value(is),symbol_type(is)
 107     format(1x,i3,1x,a11,1x,a30,1x,a30)
      end do
             
c     get the number of files in this input file
      read(symbol_value(nsymboluptodata),*,iostat=ierr) nfiles
      if(nsymboluptodata.ge.nsymboldata.or.ierr.ne.0) nfiles=1
      
      write(6,102) "Complete setup_symbol_generic",
     $     "-- symbols() have been setup"
      ssg_ierr = 0
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        subroutine setup_symbols_changeato()
        include 'epsr_dimension.inc'
        include 'symbols.inc'
        include 'molecom.inc'
        include 'clustercom.inc'
c       
c       sets the symbols and their default values for input to EPSR.
c       
c       
c       set up the symbol table for the EPSR .inp file
c       

        iref=0
        
        iref=iref+1
        symbol(iref)='bond'
        symbol_type(iref)='c'
        symbol_text(iref)='Intra-molecular bond lengths.
     $Type y to change.'
        iflag(iref)=0
        symbol_default(iref)='n'

        iref=iref+1
        symbol(iref)='label'
        symbol_type(iref)='c'
        symbol_text(iref)='Atom labels. Type y to change.'
        iflag(iref)=0
        symbol_default(iref)='n'

        iref=iref+1
        symbol(iref)='density'
        symbol_type(iref)='r'
        symbol_text(iref)='Density of this .ato file.'
        iflag(iref)=0
        symbol_default(iref)='0.1'

        iref=iref+1
        symbol(iref)='temp'
        symbol_type(iref)='r'
        symbol_text(iref)='Temperature of this .ato file.'
        iflag(iref)=0
        symbol_default(iref)='300'

        iref=iref+1
        symbol(iref)='stepmi'
        symbol_type(iref)='r'
        symbol_text(iref)='Intramolecular translation step.'
        iflag(iref)=0
        symbol_default(iref)='0.1'

        iref=iref+1
        symbol(iref)='stepri'
        symbol_type(iref)='r'
        symbol_text(iref)='Intramolecular rotation step.'
        iflag(iref)=0
        symbol_default(iref)='1.0'

        iref=iref+1
        symbol(iref)='stepr'
        symbol_type(iref)='r'
        symbol_text(iref)='Intermolecular rotation step.'
        iflag(iref)=0
        symbol_default(iref)='1.0'

        iref=iref+1
        symbol(iref)='stepm'
        symbol_type(iref)='r'
        symbol_text(iref)='Intermolecular translation step.'
        iflag(iref)=0
        symbol_default(iref)='0.2'

        iref=iref+1
        symbol(iref)='vtemp'
        symbol_type(iref)='r'
        symbol_text(iref)='Vibrational temperature [65].'
        iflag(iref)=0
        symbol_default(iref)='65'

        iref=iref+1
        symbol(iref)='addbond'
        symbol_type(iref)='c'
        symbol_text(iref)='Create a new intra-molecular bond. Type y.'
        iflag(iref)=0
        symbol_default(iref)='n'

        iref=iref+1
        symbol(iref)='removebond'
        symbol_type(iref)='c'
        symbol_text(iref)='Remove an intra-molecular bond. Type y.'
        iflag(iref)=0
        symbol_default(iref)='n'


        iref=iref+1
        symbol(iref)='natomtypes'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of atom types in this file
     $- set by .ato file'
        iflag(iref)=1
        symbol_default(iref)='1'

        nsymboluptodata=iref

        iref=iref+1
        symbol(iref)='atlabel'
        symbol_type(iref)='c'
        symbol_text(iref)='Atom label in .ato file'
        iflag(iref)=0
        symbol_default(iref)=' '

        iref=iref+1
        symbol(iref)='epsilon'
        symbol_type(iref)='r'
        symbol_text(iref)='Lennard-Jones epsilon [kJ/mol].'
        iflag(iref)=0
        symbol_default(iref)='0.65'

        iref=iref+1
        symbol(iref)='sigma'
        symbol_type(iref)='r'
        symbol_text(iref)='Lennard-Jones sigma [A].'
        iflag(iref)=0
        symbol_default(iref)='3.0'

        iref=iref+1
        symbol(iref)='atmass'
        symbol_type(iref)='r'
        symbol_text(iref)='Atomic mass [amu].'
        iflag(iref)=0
        symbol_default(iref)='12'

        iref=iref+1
        symbol(iref)='charge'
        symbol_type(iref)='r'
        symbol_text(iref)='Coulomb charge [e].'
        iflag(iref)=0
        symbol_default(iref)='0.0'

        iref=iref+1
        symbol(iref)='atsymbol'
        symbol_type(iref)='c'
        symbol_text(iref)='Chemical symbol for this component.'
        iflag(iref)=0
        symbol_default(iref)='XXX'

        iref=iref+1
        symbol(iref)='iexchange'
        symbol_type(iref)='i'
        symbol_text(iref)='= 1 if atom exchanges with other atoms,
     $else 0.'
        iflag(iref)=0
        symbol_default(iref)='-1'

        nsymbol=iref
        nsymboldata=iref
        nsymbolplus=nsymboldata-nsymboluptodata

        nfiles=1

c       
c       initialise the symbol_values to their default values
c       
        do icount=1,nsymbol
           symbol_value(icount)=symbol_default(icount)
        end do

c       
c       signal the type of symbols that have been set up
c       
        prompt='changeato>'
        ilenprompt=index(prompt,'>')
c       
c       set data prompt text
c       
        dataprompt='atom'
        ilendataprompt=4
        dataprompttext='jumps to the specified atom type number'

c       
c       setup the extension for the input file
c       
        extinp='.ato'
        ilenextinp=index(extinp,' ')-1
        
        return
        end

        subroutine setup_symbols_partials()
        include 'epsr_dimension.inc'
        include 'symbols.inc'
c       
c       sets the symbols and their default values for input to clusters.
c       
c       
c       set up the symbol table for the PARTIALS input file
c       
c       iflag decides whether this parameter can be changed
c       a parameter can only be changed if ireset >= iflag
c       

        iref=0

        iref=iref+1
        symbol(iref)='fnameato'
        symbol_type(iref)='c'
        symbol_text(iref)='Name of .ato file '
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='nz'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of r values [max: 1000]'
        iflag(iref)=0
        symbol_default(iref)='300'

        iref=iref+1
        symbol(iref)='nsumt'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of configurations already accumulated'
        iflag(iref)=0
        symbol_default(iref)='-1'

        iref=iref+1
        symbol(iref)='ndist'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of site-site distributions'
        iflag(iref)=1
        symbol_default(iref)='1'

        nsymboluptodata=iref
        nsymboldata=iref
        nsymbolplus=nsymboldata-nsymboluptodata

        nsymbol=iref

c       
c       initialise symbol_value to the default values
c       
        do is=1,nsymbol
           symbol_value(is)=symbol_default(is)
           write(6,100) is,symbol(is),symbol_value(is)
 100       format(1x,i3,1x,a11,1x,a30)
        end do

c       
c       get the number of files in this input file
c       
        read(symbol_value(nsymboluptodata),*,iostat=ierr) nfiles
        if(nsymboluptodata.ge.nsymboldata.or.ierr.ne.0) nfiles=1

c       
c       setup the output prompt
c       
        prompt='setup partials>'
        ilenprompt=index(prompt,'>')
c       
c       set up the data prompt
c       
        dataprompt='partial'
        dataprompttext='jumps to the specified cluster '
        ilendataprompt=index(dataprompt,' ')-1
c       
c       setup the extension for the input file
c       
        extinp='.PARTIALS.dat'
        ilenextinp=index(extinp,' ')-1

        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine setup_symbols_coord()
        include 'epsr_dimension.inc'
        include 'symbols.inc'
c       
c       sets the symbols and their default values for input to clusters.
c       
c       
c       set up the symbol table for the COORD input file
c       
c       iflag decides whether this parameter can be changed
c       a parameter can only be changed if ireset >= iflag
c       

        iref=0

        iref=iref+1
        symbol(iref)='fnameato'
        symbol_type(iref)='c'
        symbol_text(iref)='Name of .ato file '
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='nsize'
        symbol_type(iref)='i'
        symbol_text(iref)='Maximum coordination number (max 200)'
        iflag(iref)=0
        symbol_default(iref)='200'

        iref=iref+1
        symbol(iref)='nsumt'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of configurations already accumulated'
        iflag(iref)=0
        symbol_default(iref)='-1'

        iref=iref+1
        symbol(iref)='ndist'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of coordination number distributions'
        iflag(iref)=0
        symbol_default(iref)='1'

        nsymboluptodata=iref

        iref=iref+1
        symbol(iref)='atom1'
        symbol_type(iref)='c'
        symbol_text(iref)='Atom type 1 to define a coordination number'
        iflag(iref)=0
        symbol_default(iref)='XXX'

        iref=iref+1
        symbol(iref)='atom2'
        symbol_type(iref)='c'
        symbol_text(iref)='Atom type 2 to define a coordination number'
        iflag(iref)=0
        symbol_default(iref)='XXX'

        iref=iref+1
        symbol(iref)='rmin'
        symbol_type(iref)='r'
        symbol_text(iref)='Minimum distance for this bond'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='rmax'
        symbol_type(iref)='r'
        symbol_text(iref)='Maximum distance for this bond'
        iflag(iref)=0
        symbol_default(iref)='2'

        nsymboldata=iref
        nsymbolplus=nsymboldata-nsymboluptodata

        nsymbol=iref

c       
c       initialise symbol_value to the default values
c       
        do is=1,nsymbol
           symbol_value(is)=symbol_default(is)
           write(6,100) is,symbol(is),symbol_value(is)
 100       format(1x,i3,1x,a11,1x,a30)
        end do

c       
c       get the number of files in this input file
c       
        read(symbol_value(nsymboluptodata),*,iostat=ierr) nfiles
        if(nsymboluptodata.ge.nsymboldata.or.ierr.ne.0) nfiles=1

c       
c       setup the output prompt
c       
        prompt='setup coord>'
        ilenprompt=index(prompt,'>')
c       
c       set up the data prompt
c       
        dataprompt='distribution'
        dataprompttext='jumps to the specified distribution '
        ilendataprompt=index(dataprompt,' ')-1
c       
c       setup the extension for the input file
c       
        extinp='.COORD.dat'
        ilenextinp=index(extinp,' ')-1

        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine setup_symbols_clusters()
        include 'epsr_dimension.inc'
        include 'symbols.inc'
c       
c       sets the symbols and their default values for input to clusters.
c       
c       
c       set up the symbol table for the CLUSTERS input file
c       
c       iflag decides whether this parameter can be changed
c       a parameter can only be changed if ireset >= iflag
c       

        iref=0

        iref=iref+1
        symbol(iref)='fnameato'
        symbol_type(iref)='c'
        symbol_text(iref)='Name of .ato file '
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='nsize'
        symbol_type(iref)='i'
        symbol_text(iref)='Maximum cluster size (max 1000)'
        iflag(iref)=0
        symbol_default(iref)='1000'

        iref=iref+1
        symbol(iref)='nsumt'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of configurations already accumulated'
        iflag(iref)=0
        symbol_default(iref)='-1'

        iref=iref+1
        symbol(iref)='ndist'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of cluster distributions'
        iflag(iref)=0
        symbol_default(iref)='1'

        nsymboluptodata=iref

        iref=iref+1
        symbol(iref)='atom1'
        symbol_type(iref)='c'
        symbol_text(iref)='Atom type 1 to define a bond'
        iflag(iref)=0
        symbol_default(iref)='XXX'

        iref=iref+1
        symbol(iref)='atom2'
        symbol_type(iref)='c'
        symbol_text(iref)='Atom type 2 to define a bond'
        iflag(iref)=0
        symbol_default(iref)='XXX'

        iref=iref+1
        symbol(iref)='rmin'
        symbol_type(iref)='r'
        symbol_text(iref)='Minimum distance for this bond'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='rmax'
        symbol_type(iref)='r'
        symbol_text(iref)='Maximum distance for this bond'
        iflag(iref)=0
        symbol_default(iref)='2'
        
        nsymboldata=iref
        nsymbolplus=nsymboldata-nsymboluptodata

        nsymbol=iref

c       
c     initialise symbol_value to the default values
c       
        do is=1,nsymbol
           symbol_value(is)=symbol_default(is)
           write(6,100) is,symbol(is),symbol_value(is)
 100       format(1x,i3,1x,a11,1x,a30)
        end do

c       
c     get the number of files in this input file
c       
        read(symbol_value(nsymboluptodata),*,iostat=ierr) nfiles
        if(nsymboluptodata.ge.nsymboldata.or.ierr.ne.0) nfiles=1

c       
c     setup the output prompt
c       
        prompt='setup clusters>'
        ilenprompt=index(prompt,'>')
c     
c     set up the data prompt
c       
        dataprompt='cluster'
        dataprompttext='Specifies a cluster '
        ilendataprompt=index(dataprompt,' ')-1
c     
c     setup the extension for the input file
c       
        extinp='.CLUSTERS.dat'
        ilenextinp=index(extinp,' ')-1
        
        return
        end

c$$$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine setup_symbols_rings()
        include 'epsr_dimension.inc'
        include 'symbols.inc'
        include 'chains.inc'
c       
c       sets the symbols and their default values for input to clusters.
c       
c       
c       set up the symbol table for the COORD input file
c       
c       iflag decides whether this parameter can be changed
c       a parameter can only be changed if ireset >= iflag
c       

        iref=0

        iref=iref+1
        symbol(iref)='fnameato'
        symbol_type(iref)='c'
        symbol_text(iref)='Name of .ato file '
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='nsize'
        symbol_type(iref)='i'
        write(symbol_text(iref),101) mring
 101    format('Maximum ring length (max ',i2,')')
        iflag(iref)=0
        write(symbol_default(iref),*) mring

        iref=iref+1
        symbol(iref)='nsumt'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of configurations already accumulated'
        iflag(iref)=0
        symbol_default(iref)='-1'

        iref=iref+1
        symbol(iref)='ndist'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of ring distributions'
        iflag(iref)=0
        symbol_default(iref)='1'

        nsymboluptodata=iref

        iref=iref+1
        symbol(iref)='atom1'
        symbol_type(iref)='c'
        symbol_text(iref)='Atom type 1 for first atom in first molecule'
        iflag(iref)=0
        symbol_default(iref)='XXX'

        iref=iref+1
        symbol(iref)='atom2'
        symbol_type(iref)='c'
        symbol_text(iref)='Atom type 2 for second atom in first
     $molecule'
        iflag(iref)=0
        symbol_default(iref)='XXX'

        iref=iref+1
        symbol(iref)='atom3'
        symbol_type(iref)='c'
        symbol_text(iref)='Atom type 3 for atom in second molecule'
        iflag(iref)=0
        symbol_default(iref)='XXX'

        iref=iref+1
        symbol(iref)='rmin13'
        symbol_type(iref)='r'
        symbol_text(iref)='Minimum distance for atom1-atom3'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='rmax13'
        symbol_type(iref)='r'
        symbol_text(iref)='Maximum distance for atom1-atom3'
        iflag(iref)=0
        symbol_default(iref)='2'

        iref=iref+1
        symbol(iref)='rmin23'
        symbol_type(iref)='r'
        symbol_text(iref)='Minimum distance for atom2-atom3'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='rmax23'
        symbol_type(iref)='r'
        symbol_text(iref)='Maximum distance for atom2-atom3'
        iflag(iref)=0
        symbol_default(iref)='2'

        nsymboldata=iref
        nsymbolplus=nsymboldata-nsymboluptodata

        nsymbol=iref

c       
c       initialise symbol_value to the default values
c       
        do is=1,nsymbol
           symbol_value(is)=symbol_default(is)
           write(6,100) is,symbol(is),symbol_value(is)
 100       format(1x,i3,1x,a11,1x,a30)
        end do

c       
c       get the number of files in this input file
c       
        read(symbol_value(nsymboluptodata),*,iostat=ierr) nfiles
        if(nsymboluptodata.ge.nsymboldata.or.ierr.ne.0) nfiles=1

c       
c       setup the output prompt
c       
        prompt='setup rings>'
        ilenprompt=index(prompt,'>')
c       
c       set up the data prompt
c       
        dataprompt='ring'
        dataprompttext='jumps to the specified ring '
        ilendataprompt=index(dataprompt,' ')-1
c       
c       setup the extension for the input file
c       
        extinp='.RINGS.dat'
        ilenextinp=index(extinp,' ')-1

        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine setup_symbols_chains()
        include 'epsr_dimension.inc'
        include 'symbols.inc'
        include 'chains.inc'
c       
c       sets the symbols and their default values for input to clusters.
c       
c       
c       set up the symbol table for the COORD input file
c       
c       iflag decides whether this parameter can be changed
c       a parameter can only be changed if ireset >= iflag
c       

        iref=0

        iref=iref+1
        symbol(iref)='fnameato'
        symbol_type(iref)='c'
        symbol_text(iref)='Name of .ato file '
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='nsize'
        symbol_type(iref)='i'
        write(symbol_text(iref),101) mchain
 101    format('Maximum chain length (max ',i3,')')
        iflag(iref)=0
        write(symbol_default(iref),*) mchain

        iref=iref+1
        symbol(iref)='nsumt'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of configurations already accumulated'
        iflag(iref)=0
        symbol_default(iref)='-1'

        iref=iref+1
        symbol(iref)='ndist'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of chain distributions'
        iflag(iref)=0
        symbol_default(iref)='1'

        nsymboluptodata=iref

        iref=iref+1
        symbol(iref)='atom1'
        symbol_type(iref)='c'
        symbol_text(iref)='Atom type 1 for first atom in first molecule'
        iflag(iref)=0
        symbol_default(iref)='XXX'

        iref=iref+1
        symbol(iref)='atom2'
        symbol_type(iref)='c'
        symbol_text(iref)='Atom type 2 for second atom in first
     $molecule'
        iflag(iref)=0
        symbol_default(iref)='XXX'

        iref=iref+1
        symbol(iref)='atom3'
        symbol_type(iref)='c'
        symbol_text(iref)='Atom type 3 for atom in second molecule'
        iflag(iref)=0
        symbol_default(iref)='XXX'

        iref=iref+1
        symbol(iref)='rmin13'
        symbol_type(iref)='r'
        symbol_text(iref)='Minimum distance for atom1-atom3'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='rmax13'
        symbol_type(iref)='r'
        symbol_text(iref)='Maximum distance for atom1-atom3'
        iflag(iref)=0
        symbol_default(iref)='2'

        iref=iref+1
        symbol(iref)='rmin23'
        symbol_type(iref)='r'
        symbol_text(iref)='Minimum distance for atom2-atom3'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='rmax23'
        symbol_type(iref)='r'
        symbol_text(iref)='Maximum distance for atom2-atom3'
        iflag(iref)=0
        symbol_default(iref)='2'

        nsymboldata=iref
        nsymbolplus=nsymboldata-nsymboluptodata

        nsymbol=iref

c       
c       initialise symbol_value to the default values
c       
        do is=1,nsymbol
           symbol_value(is)=symbol_default(is)
           write(6,100) is,symbol(is),symbol_value(is)
 100       format(1x,i3,1x,a11,1x,a30)
        end do

c       
c       get the number of files in this input file
c       
        read(symbol_value(nsymboluptodata),*,iostat=ierr) nfiles
        if(nsymboluptodata.ge.nsymboldata.or.ierr.ne.0) nfiles=1

c       
c       setup the output prompt
c       
        prompt='setup chains>'
        ilenprompt=index(prompt,'>')
c       
c       set up the data prompt
c       
        dataprompt='chain'
        dataprompttext='jumps to the specified chain '
        ilendataprompt=index(dataprompt,' ')-1
c       
c       setup the extension for the input file
c       
        extinp='.CHAINS.dat'
        ilenextinp=index(extinp,' ')-1

        return
        end

        subroutine setup_symbols_sharm()
        include 'epsr_dimension.inc'
        include 'symbols.inc'
c       
c       sets the symbols and their default values for input to SHARM.
c       
c       iflag decides whether this parameter can be changed
c       a parameter can only be changed if iflag=0
c       

        iref=0

        iref=iref+1
        symbol(iref)='fnameato'
        symbol_type(iref)='c'
        symbol_text(iref)='Name of .ato file '
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='nr'
        symbol_type(iref)='c'
        symbol_text(iref)='Number of radius values (max 200)'
        iflag(iref)=0
        symbol_default(iref)='100'

        iref=iref+1
        symbol(iref)='rmax'
        symbol_type(iref)='r'
        symbol_text(iref)='Maximum radius for spherical
     $harmonic coefficients'
        iflag(iref)=0
        symbol_default(iref)='10'

        iref=iref+1
        symbol(iref)='nsumt'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of configurations already accumulated'
        iflag(iref)=0
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='ncoeffs'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of coefficients
     $(program calculates this)'
        iflag(iref)=1
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='l1values'
        symbol_type(iref)='i'
        symbol_text(iref)='L1 values (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='l2values'
        symbol_type(iref)='i'
        symbol_text(iref)='L2 values (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='lvalues'
        symbol_type(iref)='i'
        symbol_text(iref)='L values (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='n1step'
        symbol_type(iref)='i'
        symbol_text(iref)='Step in N1 values'
        iflag(iref)=0
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='n2step'
        symbol_type(iref)='i'
        symbol_text(iref)='Step in N2 values'
        iflag(iref)=0
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='atom-c'
        symbol_type(iref)='c'
        symbol_text(iref)='Central molecule - list of centre atom types'
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='axisc1'
        symbol_type(iref)='c'
        symbol_text(iref)='First axis definition for central molecule'
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='axisc2'
        symbol_type(iref)='c'
        symbol_text(iref)='Second axis definition for central molecule'
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='atom-s'
        symbol_type(iref)='c'
        symbol_text(iref)='Second molecule - list of centre atom types'
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='axiss1'
        symbol_type(iref)='c'
        symbol_text(iref)='First axis definition for second molecule'
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='axiss2'
        symbol_type(iref)='c'
        symbol_text(iref)='Second axis definition for second molecule'
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        nsymbol=iref
        nsymboldata=iref
        nsymboluptodata=iref+1
        nsymbolplus=nsymboldata-nsymboluptodata

c       
c       now write the defaults to symbol_value
c       
        do is=1,nsymbol
           symbol_value(is)=symbol_default(is)
        end do

c       
c       set the number of files - not used in SHARM
c       
        nfiles=1

c       
c       setup the output prompt
c       
        prompt='Setup sharm>'
        ilenprompt=index(prompt,'>')

c       
c       set up the data prompt - not used in SHARM
c       
        dataprompt=' '
        dataprompttext=' '
        ilendataprompt=0

c       
c       setup the extension for the input file
c       
        extinp='.SHARM.dat'
        ilenextinp=index(extinp,' ')-1

        return
        end

        subroutine setup_symbols_sdf()
        include 'epsr_dimension.inc'
        include 'symbols.inc'
c       
c       sets the symbols and their default values for input to SDF.
c       
c       iflag decides whether this parameter can be changed
c       a parameter can only be changed if iflag=0
c       

        iref=0

        iref=iref+1
        symbol(iref)='fnameato'
        symbol_type(iref)='c'
        symbol_text(iref)='Name of .ato file '
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='nr'
        symbol_type(iref)='c'
        symbol_text(iref)='Number of radius values (max 200)'
        iflag(iref)=0
        symbol_default(iref)='130'

        iref=iref+1
        symbol(iref)='rmax'
        symbol_type(iref)='r'
        symbol_text(iref)='Maximum radius for spherical harmonic
     $coefficients'
        iflag(iref)=0
        symbol_default(iref)='13'

        iref=iref+1
        symbol(iref)='nsumt'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of configurations already accumulated'
        iflag(iref)=0
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='ncoeffs'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of coefficients
     $(program calculates this)'
        iflag(iref)=1
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='l1values'
        symbol_type(iref)='i'
        symbol_text(iref)='L1 values (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='l2values'
        symbol_type(iref)='i'
        symbol_text(iref)='L2 values (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='lvalues'
        symbol_type(iref)='i'
        symbol_text(iref)='L values (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='n1step'
        symbol_type(iref)='i'
        symbol_text(iref)='Step in N1 values'
        iflag(iref)=0
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='n2step'
        symbol_type(iref)='i'
        symbol_text(iref)='Step in N2 values'
        iflag(iref)=0
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='atom-c'
        symbol_type(iref)='c'
        symbol_text(iref)='List of centre atom types.
     $First will be origin.'
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='pairdist-c'
        symbol_type(iref)='r'
        symbol_text(iref)='List of pair distances for centre molecule.'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='fracdist-c'
        symbol_type(iref)='r'
        symbol_text(iref)='Acceptance half-width
     $(fraction of pair distance).'
        iflag(iref)=0
        symbol_default(iref)='0.2'

        iref=iref+1
        symbol(iref)='axisc1'
        symbol_type(iref)='c'
        symbol_text(iref)='First axis definition for centre molecule.'
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='axisc2'
        symbol_type(iref)='c'
        symbol_text(iref)='Second axis definition for centre molecule.'
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='atom-s'
        symbol_type(iref)='c'
        symbol_text(iref)='List of second atom types.
     $First will be origin.'
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='pairdist-s'
        symbol_type(iref)='r'
        symbol_text(iref)='List of pair distances for second molecule.'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='fracdist-s'
        symbol_type(iref)='r'
        symbol_text(iref)='Acceptance half-width (fraction of
     $pair distance).'
        iflag(iref)=0
        symbol_default(iref)='0.2'

        iref=iref+1
        symbol(iref)='axiss1'
        symbol_type(iref)='c'
        symbol_text(iref)='First axis definition for second molecule.'
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='axiss2'
        symbol_type(iref)='c'
        symbol_text(iref)='Second axis definition for second molecule.'
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        nsymbol=iref
        nsymboldata=iref
        nsymboluptodata=iref+1
        nsymbolplus=nsymboldata-nsymboluptodata

c       
c       now write the defaults to symbol_value
c       
        do is=1,nsymbol
           symbol_value(is)=symbol_default(is)
        end do

c       
c       set the number of files - not used in SHARM
c       
        nfiles=1

c       
c       setup the output prompt
c       
        prompt='Setup sdf>'
        ilenprompt=index(prompt,'>')

c       
c       set up the data prompt - not used in SDF
c       
        dataprompt=' '
        dataprompttext=' '
        ilendataprompt=0
c       
c       setup the extension for the input file
c       
        extinp='.SDF.dat'
        ilenextinp=index(extinp,' ')-1

        return
        end
      
      subroutine setup_symbols_plot2d()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
c     
c       sets the symbols and their default values for input to clusters.
c       
c       
c       set up the symbol table for the plot2d input file
c       
c       iflag decides whether this parameter can be changed
c       a parameter can only be changed if iflag=0
c       
c       define symbols and default values
c       
        iref=0

        iref=iref+1
        symbol(iref)='bakcolour'
        symbol_type(iref)='i'
        symbol_text(iref)='Background colour index [1]'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='shcoeffs'
        symbol_type(iref)='c'
        symbol_text(iref)='Name of file containing SHARM coefficients'
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='ncoeffs'
        symbol_type(iref)='i'
        symbol_text(iref)='no. of SHARM coefficients - '
     $ //'determined from coefficients file'
        iflag(iref)=1
        symbol_default(iref)=' '

        iref=iref+1
        symbol(iref)='ident'
        symbol_type(iref)='i'
        symbol_text(iref)='= 0 for identical molecules,'
     *       //' else 1 if different'
        iflag(iref)=0
        symbol_default(iref)=' '
        
        iref=iref+1
        symbol(iref)='setone'
        symbol_type(iref)='i'
        symbol_text(iref)='0 sets first coefficient to zero -
     $normally 1'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='nsmoo'
        symbol_type(iref)='i'
        symbol_text(iref)='number of smoothings on coefficients'
        iflag(iref)=0
        symbol_default(iref)='4'

        iref=iref+1
        symbol(iref)='xdim'
        symbol_type(iref)='i'
        symbol_text(iref)='x dimension of plotting square'
        iflag(iref)=0
        symbol_default(iref)='10'

        iref=iref+1
        symbol(iref)='ydim'
        symbol_type(iref)='i'
        symbol_text(iref)='y dimension of plotting square'
        iflag(iref)=0
        symbol_default(iref)='10'

        iref=iref+1
        symbol(iref)='rmin'
        symbol_type(iref)='r'
        symbol_text(iref)='minimum radius of plot'
        iflag(iref)=0
        symbol_default(iref)='2.0'

        iref=iref+1
        symbol(iref)='ncircles'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of circles at centre of plot (max 25)'
        iflag(iref)=0
        symbol_default(iref)='1'

        nsymboluptodata=iref

        iref=iref+1
        symbol(iref)='xyzcircle'
        symbol_type(iref)='r'
        symbol_text(iref)='(x,y,z) coordinates for circle'
     $  //' (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='0.0 0 0'

        iref=iref+1
        symbol(iref)='radcircle'
        symbol_type(iref)='r'
        symbol_text(iref)='Radius of circle'
        iflag(iref)=0
        symbol_default(iref)='1.0'

        nsymboldata=iref
        nsymbolplus=nsymboldata-nsymboluptodata

        iref=iref+1
        symbol(iref)='nplotxy'
        symbol_type(iref)='i'
        symbol_text(iref)='no. of plots along x- and y-axis
     $[set at 1 1]'
        iflag(iref)=1
        symbol_default(iref)='1 1'

        iref=iref+1
        symbol(iref)='aspect'
        symbol_type(iref)='r'
        symbol_text(iref)='aspect ratio of plot [1.0]'
        iflag(iref)=0
        symbol_default(iref)='1.0'

        iref=iref+1
        symbol(iref)='use_l1_l2'
        symbol_type(iref)='i'
        symbol_text(iref)='use l1 and l2 (1 or 0)'
        iflag(iref)=0
        symbol_default(iref)='1 0'

        iref=iref+1
        symbol(iref)='use_n1_n2'
        symbol_type(iref)='i'
        symbol_text(iref)='use n1 and n2 (1 or 0)'
        iflag(iref)=0
        symbol_default(iref)='1 0'

        iref=iref+1
        symbol(iref)='use_m'
        symbol_type(iref)='i'
        symbol_text(iref)='use m2 (1 or 0)'
        iflag(iref)=0
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='nvary'
        symbol_type(iref)='i'
        symbol_text(iref)='vary thetal (1), phil (2), phim(3),
     $thetam(4), chim(5)'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='ph_th_ch'
        symbol_type(iref)='c'
        symbol_text(iref)='phil, phim, thetam, chim (degrees
     $- separated by spaces'
        iflag(iref)=0
        symbol_default(iref)='0 0 0 0'

        iref=iref+1
        symbol(iref)='axespar'
        symbol_type(iref)='r'
        symbol_text(iref)='axes character size, line width and colour'
     $   //' (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='1.5 2 1'

        iref=iref+1
        symbol(iref)='plottitle'
        symbol_type(iref)='c'
        symbol_text(iref)='plot title'
        iflag(iref)=0
        symbol_default(iref)=' '

        iref=iref+1
        symbol(iref)='titlecoord'
        symbol_type(iref)='r'
        symbol_text(iref)='(x,y) coords. of title, and character size'
     *       //' (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='0.1 -1.3 2'

        iref=iref+1
        symbol(iref)='colour'
        symbol_type(iref)='i'
        symbol_text(iref)='Colour scale [2 - 8]'
        iflag(iref)=0
        symbol_default(iref)='3'

        iref=iref+1
        symbol(iref)='zmin_zmax'
        symbol_type(iref)='r'
        symbol_text(iref)='range of intensites to be plotted'
     *       //' (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='1 5'

        iref=iref+1
        symbol(iref)='contrast'
        symbol_type(iref)='r'
        symbol_text(iref)='Contrast factor [0.0 - 1.0]'
        iflag(iref)=0
        symbol_default(iref)='1.0'

        iref=iref+1
        symbol(iref)='contover'
        symbol_type(iref)='c'
        symbol_text(iref)='Contours overlaid? - set to N -
     $cannot be changed'
        iflag(iref)=1
        symbol_default(iref)='N'

        iref=iref+1
        symbol(iref)='extraline'
        symbol_type(iref)='i'
        symbol_text(iref)='extra lines (0) - cannot be set'
        iflag(iref)=1
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='extratext'
        symbol_type(iref)='i'
        symbol_text(iref)='extra text (0) - cannot be set'
        iflag(iref)=1
        symbol_default(iref)='0'

c       Last symbol is the extension of the file containing the spherical harmonic coefficients
c       Used only for setup, not in the plotting program itself.

        iref=iref+1
        symbol(iref)='extcoeff'
        symbol_type(iref)='c'
        symbol_text(iref)='extension of coefficients file'
        iflag(iref)=0
        symbol_default(iref)='.SHARM.h01'

        nsymbol=iref

c       
c       initialise symbol_value to the default values
c       
        do is=1,nsymbol
           symbol_value(is)=symbol_default(is)
           write(6,100) is,symbol(is),symbol_value(is)
 100       format(1x,i3,1x,a11,1x,a30)
        end do

c       
c       get the number of files in this input file
c       
        read(symbol_value(nsymboluptodata),*,iostat=ierr) nfiles
        if(nsymboluptodata.ge.nsymboldata.or.ierr.ne.0) nfiles=1

c       
c       setup the output prompt
c       
        prompt='setup plot2d>'
        ilenprompt=index(prompt,'>')
c       
c       set up the data prompt
c       
        dataprompt='circle'
        dataprompttext='jumps to the specified circle '
        ilendataprompt=index(dataprompt,' ')-1

c       
c       setup the extension for the input file
c       
        extinp='.plot2d.txt'
        ilenextinp=index(extinp,' ')-1

        return
        end

        subroutine setup_symbols_plot3d()
        include 'epsr_dimension.inc'
        include 'symbols.inc'
c       
c       sets the symbols and their default values for input to clusters.
c       
c       
c       set up the symbol table for the plot3d input file
c       
c       iflag decides whether this parameter can be changed
c       a parameter can only be changed if iflag=0
c       
c       define symbols and default values
c       
        iref=0

        iref=iref+1
        symbol(iref)='shcoeffs'
        symbol_type(iref)='c'
        symbol_text(iref)='Name of file containing spherical
     $harmonic coefficients'
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='ncoeffs'
        symbol_type(iref)='i'
        symbol_text(iref)='no. of coefficients - '
     *       //'determined from coefficients file'
        iflag(iref)=1
        symbol_default(iref)=' '

        iref=iref+1
        symbol(iref)='ident'
        symbol_type(iref)='i'
        symbol_text(iref)='= 0 for identical molecules,'
     *       //' else 1 if different'
        iflag(iref)=0
        symbol_default(iref)=' '
        
        iref=iref+1
        symbol(iref)='setone'
        symbol_type(iref)='i'
        symbol_text(iref)='0 sets first coefficient to zero -
     $normally 1'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='nsmoo'
        symbol_type(iref)='i'
        symbol_text(iref)='number of smoothings on coefficients'
        iflag(iref)=0
        symbol_default(iref)='4'

        iref=iref+1
        symbol(iref)='radmax'
        symbol_type(iref)='r'
        symbol_text(iref)='maximum radius of plotting box'
        iflag(iref)=0
        symbol_default(iref)='5'

        iref=iref+1
        symbol(iref)='nplotxy'
        symbol_type(iref)='i'
        symbol_text(iref)='no. of plots along x- and y-axis
     $[set at 1 1]'
        iflag(iref)=0
        symbol_default(iref)='1 1'

        iref=iref+1
        symbol(iref)='aspect'
        symbol_type(iref)='r'
        symbol_text(iref)='aspect ratio of plot [1.0]'
        iflag(iref)=0
        symbol_default(iref)='1.0'

        iref=iref+1
        symbol(iref)='rmin_rmax'
        symbol_type(iref)='r'
        symbol_text(iref)='minimum and maximum radius of plot'
        iflag(iref)=0
        symbol_default(iref)='2.0 5.0'

        iref=iref+1
        symbol(iref)='surfra'
        symbol_type(iref)='r'
        symbol_text(iref)='fractional isosurface level
     $(-ve for absolute)'
        iflag(iref)=0
        symbol_default(iref)='0.15'

        iref=iref+1
        symbol(iref)='use_l1_l2'
        symbol_type(iref)='i'
        symbol_text(iref)='use l1 and l2 (1 or 0)'
        iflag(iref)=0
        symbol_default(iref)='1 0'

        iref=iref+1
        symbol(iref)='use_n1_n2'
        symbol_type(iref)='i'
        symbol_text(iref)='use n1 and n2 (1 or 0)'
        iflag(iref)=0
        symbol_default(iref)='1 0'

        iref=iref+1
        symbol(iref)='use_m'
        symbol_type(iref)='i'
        symbol_text(iref)='use m2 (1 or 0)'
        iflag(iref)=0
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='nvary'
        symbol_type(iref)='i'
        symbol_text(iref)='vary (thetal, phil) (1),
     $(thetam, phim) (2), (thetam, chim) (3)'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='ph_th_ch'
        symbol_type(iref)='c'
        symbol_text(iref)='phim, thetam, chim (degrees -
     $separated by spaces'
        iflag(iref)=0
        symbol_default(iref)='0 0 0'

        iref=iref+1
        symbol(iref)='nsphere'
        symbol_type(iref)='i'
        symbol_text(iref)='number of spheres at centre of plot (max 25)'
        iflag(iref)=0
        symbol_default(iref)='1'

        nsymboluptodata=iref

        iref=iref+1
        symbol(iref)='radsphere'
        symbol_type(iref)='r'
        symbol_text(iref)='radius of sphere'
        iflag(iref)=0
        symbol_default(iref)='1.0'

        iref=iref+1
        symbol(iref)='rthetaphi'
        symbol_type(iref)='r'
        symbol_text(iref)='(r,theta,phi) coordinates for sphere'
     *       //' (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='0.0 0 0'

        iref=iref+1
        symbol(iref)='rgbsph'
        symbol_type(iref)='r'
        symbol_text(iref)='red green blue fractions for sphere'
     *       //' (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='0.7 0.7 0.7'

        nsymboldata=iref
        nsymbolplus=nsymboldata-nsymboluptodata

        iref=iref+1
        symbol(iref)='axespar'
        symbol_type(iref)='r'
        symbol_text(iref)='axes character size, line width and colour'
     *       //' (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='1.5 2 1'

        iref=iref+1
        symbol(iref)='plottitle'
        symbol_type(iref)='c'
        symbol_text(iref)='plot title'
        iflag(iref)=0
        symbol_default(iref)=' '

        iref=iref+1
        symbol(iref)='titlecoord'
        symbol_type(iref)='r'
        symbol_text(iref)='(x,y) coords. of title, and character size'
     *       //' (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='0.1 -1.3 2'

        iref=iref+1
        symbol(iref)='blank'
        symbol_type(iref)='c'
        symbol_text(iref)='A blank line required by plotting program '
        iflag(iref)=1
        symbol_default(iref)=' '

        iref=iref+1
        symbol(iref)='rgbbak'
        symbol_type(iref)='r'
        symbol_text(iref)='red green blue fractions for background'
     *       //' (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='0.8 0.8 1.0'

        iref=iref+1
        symbol(iref)='ishade'
        symbol_type(iref)='i'
        symbol_text(iref)='ishade (1-8): 0 means no shading,'
     *       //' -ve means inverted shading'
        iflag(iref)=0
        symbol_default(iref)='-1'

        iref=iref+1
        symbol(iref)='rgbobj'
        symbol_type(iref)='r'
        symbol_text(iref)='red green blue fractions for object'
     *       //' (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='1 1 0'

        iref=iref+1
        symbol(iref)='lightcoord'
        symbol_type(iref)='i'
        symbol_text(iref)='(x,y,z) coordinates for light source'
     *       //' (separated by spaces)'
        iflag(iref)=0
        symbol_default(iref)='2 2 0'

        iref=iref+1
        symbol(iref)='fadefc'
        symbol_type(iref)='r'
        symbol_text(iref)='fade factor (0 = no fading, 1=full fading)'
        iflag(iref)=0
        symbol_default(iref)='1.0'

        iref=iref+1
        symbol(iref)='itrans'
        symbol_type(iref)='i'
        symbol_text(iref)='transparency of object
     $(0=0%,1=25%,2=50%,3=75%)'
        iflag(iref)=0
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='appearance'
        symbol_type(iref)='r'
        symbol_text(iref)='difuse, shine, polish and contrast'
        iflag(iref)=0
        symbol_default(iref)='1.0 0.0 1.0 1.0'

        iref=iref+1
        symbol(iref)='rotelev'
        symbol_type(iref)='r'
        symbol_text(iref)='rotation and elevation of viewing
     $point (deg.)'
        iflag(iref)=0
        symbol_default(iref)='15 35'

        iref=iref+1
        symbol(iref)='extraline'
        symbol_type(iref)='i'
        symbol_text(iref)='extra lines (0) - cannot be set'
        iflag(iref)=1
        symbol_default(iref)='0'

        iref=iref+1
        symbol(iref)='extratext'
        symbol_type(iref)='i'
        symbol_text(iref)='extra text (0) - cannot be set'
        iflag(iref)=1
        symbol_default(iref)='0'

c       Last symbol is the extension of the file containing the spherical harmonic coefficients
c       Used only for setup, not in the plotting program itself.

        iref=iref+1
        symbol(iref)='extcoeff'
        symbol_type(iref)='c'
        symbol_text(iref)='extension of coefficients file'
        iflag(iref)=0
        symbol_default(iref)='.SHARM.h01'

        nsymbol=iref

c       
c       initialise symbol_value to the default values
c       
        do is=1,nsymbol
           symbol_value(is)=symbol_default(is)
           write(6,100) is,symbol(is),symbol_value(is)
 100       format(1x,i3,1x,a11,1x,a30)
        end do

c       
c       get the number of files in this input file
c       
        read(symbol_value(nsymboluptodata),*,iostat=ierr) nfiles
        if(nsymboluptodata.ge.nsymboldata.or.ierr.ne.0) nfiles=1

c       
c       setup the output prompt
c       
        prompt='setup plot3d>'
        ilenprompt=index(prompt,'>')
c       
c       set up the data prompt
c       
        dataprompt='sphere'
        dataprompttext='jumps to the specified sphere '
        ilendataprompt=index(dataprompt,' ')-1

c       
c       setup the extension for the input file
c       
        extinp='.plot3d.txt'
        ilenextinp=index(extinp,' ')-1

        return
        end

        subroutine setup_symbols_triangles()
        include 'epsr_dimension.inc'
        include 'symbols.inc'
c       
c       sets the symbols and their default values for input to tri-angles.
c       
c       iflag decides whether this parameter can be changed
c       a parameter can only be changed if ireset >= iflag
c       

        iref=0

        iref=iref+1
        symbol(iref)='fnameato'
        symbol_type(iref)='c'
        symbol_text(iref)='Name of .ato file '
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='nsize'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of theta bins [180]'
        iflag(iref)=0
        symbol_default(iref)='180'

        iref=iref+1
        symbol(iref)='nsumt'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of configurations already accumulated'
        iflag(iref)=0
        symbol_default(iref)='-1'

        iref=iref+1
        symbol(iref)='ndist'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of tri-angle distributions'
        iflag(iref)=0
        symbol_default(iref)='1'

        nsymboluptodata=iref

        iref=iref+1
        symbol(iref)='atom1'
        symbol_type(iref)='c'
        symbol_text(iref)='Atom type 1 to define a triangle'
        iflag(iref)=0
        symbol_default(iref)='XXX'

        iref=iref+1
        symbol(iref)='atom2'
        symbol_type(iref)='c'
        symbol_text(iref)='Atom type 2 to define a triangle'
        iflag(iref)=0
        symbol_default(iref)='XXX'

        iref=iref+1
        symbol(iref)='atom3'
        symbol_type(iref)='c'
        symbol_text(iref)='Atom type 3 to define a triangle'
        iflag(iref)=0
        symbol_default(iref)='XXX'

        iref=iref+1
        symbol(iref)='ltype'
        symbol_type(iref)='i'
        symbol_text(iref)='Single atom3 types [1], or multiple
     $atom3 types [2]?'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='rmin12'
        symbol_type(iref)='r'
        symbol_text(iref)='Minimum distance for atoms 1 and 2'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='rmax12'
        symbol_type(iref)='r'
        symbol_text(iref)='Maximum distance for atoms 1 and 2'
        iflag(iref)=0
        symbol_default(iref)='2'

        iref=iref+1
        symbol(iref)='rmin23'
        symbol_type(iref)='r'
        symbol_text(iref)='Minimum distance for atoms 2 and 3'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='rmax23'
        symbol_type(iref)='r'
        symbol_text(iref)='Maximum distance for atoms 2 and 3'
        iflag(iref)=0
        symbol_default(iref)='2'

        nsymboldata=iref
        nsymbolplus=nsymboldata-nsymboluptodata

        nsymbol=iref

c       
c       initialise symbol_value to the default values
c       
        do is=1,nsymbol
           symbol_value(is)=symbol_default(is)
           write(6,100) is,symbol(is),symbol_value(is)
 100       format(1x,i3,1x,a11,1x,a30)
        end do

c       
c       get the number of files in this input file
c       
        read(symbol_value(nsymboluptodata),*,iostat=ierr) nfiles
        if(nsymboluptodata.ge.nsymboldata.or.ierr.ne.0) nfiles=1

c       
c       setup the output prompt
c       
        prompt='setup triangles>'
        ilenprompt=index(prompt,'>')
c       
c       set up the data prompt
c       
        dataprompt='triangles'
        dataprompttext='jumps to the specified triangle '
        ilendataprompt=index(dataprompt,' ')-1
c       
c       setup the extension for the input file
c       
        extinp='.TRI.dat'
        ilenextinp=index(extinp,' ')-1

        return
        end

        subroutine setup_symbols_torangles()
        include 'epsr_dimension.inc'
        include 'symbols.inc'
c       
c       sets the symbols and their default values for input to torangles.
c       
c       iflag decides whether this parameter can be changed
c       a parameter can only be changed if ireset >= iflag
c       

        iref=0

        iref=iref+1
        symbol(iref)='fnameato'
        symbol_type(iref)='c'
        symbol_text(iref)='Name of .ato file '
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='nsize'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of cos(theta) bins [180]'
        iflag(iref)=0
        symbol_default(iref)='180'

        iref=iref+1
        symbol(iref)='nsumt'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of configurations already accumulated'
        iflag(iref)=0
        symbol_default(iref)='-1'

        iref=iref+1
        symbol(iref)='ndist'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of tri-angle distributions'
        iflag(iref)=0
        symbol_default(iref)='1'

        nsymboluptodata=iref

        iref=iref+1
        symbol(iref)='atom1'
        symbol_type(iref)='c'
        symbol_text(iref)='Atom type to define which molecule type'
        iflag(iref)=0
        symbol_default(iref)='XXX'

        iref=iref+1
        symbol(iref)='ltyp1'
        symbol_type(iref)='i'
        symbol_text(iref)='First atom number to define the axis'
        iflag(iref)=0
        symbol_default(iref)='1'

        iref=iref+1
        symbol(iref)='ltyp2'
        symbol_type(iref)='i'
        symbol_text(iref)='Second atom number to define the axis'
        iflag(iref)=0
        symbol_default(iref)='2'

        iref=iref+1
        symbol(iref)='ltyp3'
        symbol_type(iref)='i'
        symbol_text(iref)='Two atom numbers for vector from atom 1'
        iflag(iref)=0
        symbol_default(iref)='3 3'

        iref=iref+1
        symbol(iref)='ltyp4'
        symbol_type(iref)='i'
        symbol_text(iref)='Two atom numbers for vector from atom 2'
        iflag(iref)=0
        symbol_default(iref)='4 4'

        nsymboldata=iref
        nsymbolplus=nsymboldata-nsymboluptodata

        nsymbol=iref

c       
c       initialise symbol_value to the default values
c       
        do is=1,nsymbol
           symbol_value(is)=symbol_default(is)
           write(6,100) is,symbol(is),symbol_value(is)
 100       format(1x,i3,1x,a11,1x,a30)
        end do

c       
c       get the number of files in this input file
c       
        read(symbol_value(nsymboluptodata),*,iostat=ierr) nfiles
        if(nsymboluptodata.ge.nsymboldata.or.ierr.ne.0) nfiles=1

c       
c       setup the output prompt
c       
        prompt='setup torangles>'
        ilenprompt=index(prompt,'>')
c       
c       set up the data prompt
c       
        dataprompt='torangles'
        dataprompttext='jumps to the specified torangle '
        ilendataprompt=index(dataprompt,' ')-1
c       
c       setup the extension for the input file
c       
        extinp='.TOR.dat'
        ilenextinp=index(extinp,' ')-1

        return
        end

        subroutine setup_symbols_voids()
        include 'epsr_dimension.inc'
        include 'symbols.inc'
c       
c       sets the symbols and their default values for input to voids.
c       
c       
c       set up the symbol table for the VOIDS input file
c       
c       iflag decides whether this parameter can be changed
c       a parameter can only be changed if ireset >= iflag
c       

        iref=0

        iref=iref+1
        symbol(iref)='fnameato'
        symbol_type(iref)='c'
        symbol_text(iref)='Name of .ato file '
        iflag(iref)=0
        symbol_default(iref)='<undefined>'

        iref=iref+1
        symbol(iref)='nsize'
        symbol_type(iref)='i'
        symbol_text(iref)='No. of histogram points (max 1000)'
        iflag(iref)=0
        symbol_default(iref)='1000'

        iref=iref+1
        symbol(iref)='npix'
        symbol_type(iref)='i'
        symbol_text(iref)='No. of pixels for void calculation (max 201)'
        iflag(iref)=0
        symbol_default(iref)='100'

        iref=iref+1
        symbol(iref)='voidmax'
        symbol_type(iref)='r'
        symbol_text(iref)='Maximum void volume (in multiples of atomic volume)'
        iflag(iref)=0
        symbol_default(iref)='500.0'

        iref=iref+1
        symbol(iref)='nsumt'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of configurations already accumulated'
        iflag(iref)=0
        symbol_default(iref)='-1'

        iref=iref+1
        symbol(iref)='ndist'
        symbol_type(iref)='i'
        symbol_text(iref)='Number of void distributions'
        iflag(iref)=0
        symbol_default(iref)='1'

        nsymboluptodata=iref

        iref=iref+1
        symbol(iref)='radius'
        symbol_type(iref)='r'
        symbol_text(iref)='Average distance from atoms to define a void'
        iflag(iref)=0
        symbol_default(iref)='2'

        iref=iref+1
        symbol(iref)='exclude'
        symbol_type(iref)='c'
        symbol_text(iref)='List of atom types to exclude (max 10)'
        iflag(iref)=0
        symbol_default(iref)=' '

        nsymboldata=iref
        nsymbolplus=nsymboldata-nsymboluptodata

        nsymbol=iref

c       
c       initialise symbol_value to the default values
c       
        do is=1,nsymbol
           symbol_value(is)=symbol_default(is)
           write(6,100) is,symbol(is),symbol_value(is)
 100       format(1x,i3,1x,a11,1x,a30)
        end do

c       
c       get the number of files in this input file
c       
        read(symbol_value(nsymboluptodata),*,iostat=ierr) nfiles
        if(nsymboluptodata.ge.nsymboldata.or.ierr.ne.0) nfiles=1

c       
c       setup the output prompt
c       
        prompt='setup voids>'
        ilenprompt=index(prompt,'>')
c       
c       set up the data prompt
c       
        dataprompt='void'
        dataprompttext='Specifies a void '
        ilendataprompt=index(dataprompt,' ')-1
c       
c       setup the extension for the input file
c       
        extinp='.VOIDS.dat'
        ilenextinp=index(extinp,' ')-1

        return
        end

      subroutine setup_symbols_fluctuations()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
c
c sets the symbols and their default values for input to fluctuations.
c
c
c set up the symbol table for the FLUCTUATIONS input file
c
c iflag decides whether this parameter can be changed
c a parameter can only be changed if ireset >= iflag
c

      iref=0

      iref=iref+1
      symbol(iref)='fnameato'
      symbol_type(iref)='c'
      symbol_text(iref)='Name of .ato file '
      iflag(iref)=0
      symbol_default(iref)='<undefined>'

      iref=iref+1
      symbol(iref)='nsumt'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of samples already accumulated'
      iflag(iref)=0
      symbol_default(iref)='0'

      iref=iref+1
      symbol(iref)='nsamples'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of samples to take from this configuration'
      iflag(iref)=0
      symbol_default(iref)='5'

      iref=iref+1
      symbol(iref)='ndist'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of volume sizes'
      iflag(iref)=0
      symbol_default(iref)='1'

      nsymboluptodata=iref

      iref=iref+1
      symbol(iref)='rmax'
      symbol_type(iref)='r'
      symbol_text(iref)='Radius of volume to probe'
      iflag(iref)=0
      symbol_default(iref)='10'

      iref=iref+1
      symbol(iref)='ntotal'
      symbol_type(iref)='r'
      symbol_text(iref)='Sum of N'
      iflag(iref)=0
      symbol_default(iref)='0'

      iref=iref+1
      symbol(iref)='nsqtotal'
      symbol_type(iref)='r'
      symbol_text(iref)='Sum of N^2'
      iflag(iref)=0
      symbol_default(iref)='0'

      nsymboldata=iref
      nsymbolplus=nsymboldata-nsymboluptodata

      nsymbol=iref

c
c initialise symbol_value to the default values
c
      do is=1,nsymbol
         symbol_value(is)=symbol_default(is)
         write(6,100) is,symbol(is),symbol_value(is)
100   format(1x,i3,1x,a11,1x,a30)
      end do

c
c get the number of files in this input file
c
      read(symbol_value(nsymboluptodata),*,iostat=ierr) nfiles
      if(nsymboluptodata.ge.nsymboldata.or.ierr.ne.0) nfiles=1

c
c setup the output prompt
c
      prompt='setup fluctuations>'
      ilenprompt=index(prompt,'>')

c
c set up the data prompt
c
      dataprompt='Size'
      dataprompttext='Specifies a size to measure fluctuations '
      ilendataprompt=index(dataprompt,' ')-1
c
c setup the extension for the input file
c
      extinp='.FLUCTUATIONS.dat'
      ilenextinp=index(extinp,' ')-1

      return
      end


