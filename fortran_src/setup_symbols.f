      subroutine setup_symbols_changeato()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
c
c sets the symbols and their default values for input to EPSR.
c
c
c set up the symbol table for the EPSR .inp file
c

      iref=0

      iref=iref+1
      symbol(iref)='bond'
      symbol_type(iref)='c'
      symbol_text(iref)='Intra-molecular bond lengths. Type y to change.'
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
      symbol(iref)='vibtemp'
      symbol_type(iref)='r'
      symbol_text(iref)='Vibrational energy weight [65].'
      iflag(iref)=0
      symbol_default(iref)='65'

      iref=iref+1
      symbol(iref)='dihtemp'
      symbol_type(iref)='r'
      symbol_text(iref)='Dihedral energy weight [1.0e4].'
      iflag(iref)=0
      symbol_default(iref)='1.0e4'

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
      symbol_text(iref)='Number of atom types in this file - set by .ato file'
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
      symbol_text(iref)='= 1 if atom exchanges with other atoms, else 0.'
      iflag(iref)=0
      symbol_default(iref)='-1'

      nsymbol=iref
      nsymboldata=iref
      nsymbolplus=nsymboldata-nsymboluptodata

      nfiles=1

c
c initialise the symbol_values to their default values
c
      do icount=1,nsymbol
         symbol_value(icount)=symbol_default(icount)
      end do

c
c signal the type of symbols that have been set up
c
      prompt='changeato>'
      ilenprompt=index(prompt,'>')
c
c set data prompt text
c
      dataprompt='atom'
      ilendataprompt=4
      dataprompttext='jumps to the specified atom type number'

c
c setup the extension for the input file
c
      extinp='.ato'
      ilenextinp=index(extinp,' ')-1
      
      return
      end

      subroutine setup_symbols_partials()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
c
c sets the symbols and their default values for input to clusters.
c
c
c set up the symbol table for the PARTIALS input file
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
      prompt='setup partials>'
      ilenprompt=index(prompt,'>')
c
c set up the data prompt
c
      dataprompt='partial'
      dataprompttext='jumps to the specified cluster '
      ilendataprompt=index(dataprompt,' ')-1
c
c setup the extension for the input file
c
      extinp='.PARTIALS.dat'
      ilenextinp=index(extinp,' ')-1

      return
      end

      subroutine setup_symbols_coord()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
c
c sets the symbols and their default values for input to clusters.
c
c
c set up the symbol table for the COORD input file
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
      prompt='setup coord>'
      ilenprompt=index(prompt,'>')
c
c set up the data prompt
c
      dataprompt='distribution'
      dataprompttext='jumps to the specified distribution '
      ilendataprompt=index(dataprompt,' ')-1
c
c setup the extension for the input file
c
      extinp='.COORD.dat'
      ilenextinp=index(extinp,' ')-1

      return
      end

      subroutine setup_symbols_clusters()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
c
c sets the symbols and their default values for input to clusters.
c
c
c set up the symbol table for the CLUSTERS input file
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
      prompt='setup clusters>'
      ilenprompt=index(prompt,'>')
c
c set up the data prompt
c
      dataprompt='cluster'
      dataprompttext='Specifies a cluster '
      ilendataprompt=index(dataprompt,' ')-1
c
c setup the extension for the input file
c
      extinp='.CLUSTERS.dat'
      ilenextinp=index(extinp,' ')-1

      return
      end

      subroutine setup_symbols_rings()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
      include 'chains.inc'
c
c sets the symbols and their default values for input to clusters.
c
c
c set up the symbol table for the COORD input file
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
      symbol(iref)='nsize'
      symbol_type(iref)='i'
      write(symbol_text(iref),101) mring
101   format('Maximum ring length (max ',i2,')')
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
      symbol_text(iref)='Atom type 2 for second atom in first molecule'
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
      prompt='setup rings>'
      ilenprompt=index(prompt,'>')
c
c set up the data prompt
c
      dataprompt='ring'
      dataprompttext='jumps to the specified ring '
      ilendataprompt=index(dataprompt,' ')-1
c
c setup the extension for the input file
c
      extinp='.RINGS.dat'
      ilenextinp=index(extinp,' ')-1

      return
      end
      subroutine setup_symbols_chains()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
      include 'chains.inc'
c
c sets the symbols and their default values for input to clusters.
c
c
c set up the symbol table for the COORD input file
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
      symbol(iref)='nsize'
      symbol_type(iref)='i'
      write(symbol_text(iref),101) mchain
101   format('Maximum chain length (max ',i3,')')
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
      symbol_text(iref)='Atom type 2 for second atom in first molecule'
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
      prompt='setup chains>'
      ilenprompt=index(prompt,'>')
c
c set up the data prompt
c
      dataprompt='chain'
      dataprompttext='jumps to the specified chain '
      ilendataprompt=index(dataprompt,' ')-1
c
c setup the extension for the input file
c
      extinp='.CHAINS.dat'
      ilenextinp=index(extinp,' ')-1

      return
      end

      subroutine setup_symbols_sharm()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
c
c sets the symbols and their default values for input to SHARM.
c
c iflag decides whether this parameter can be changed
c a parameter can only be changed if iflag=0
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
      symbol_text(iref)='Maximum radius for spherical harmonic coefficients'
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
      symbol_text(iref)='Number of coefficients (program calculates this)'
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
c now write the defaults to symbol_value
c
      do is=1,nsymbol
         symbol_value(is)=symbol_default(is)
      end do

c
c set the number of files - not used in SHARM
c
      nfiles=1

c
c setup the output prompt
c
      prompt='Setup sharm>'
      ilenprompt=index(prompt,'>')

c
c set up the data prompt - not used in SHARM
c
      dataprompt=' '
      dataprompttext=' '
      ilendataprompt=0

c
c setup the extension for the input file
c
      extinp='.SHARM.dat'
      ilenextinp=index(extinp,' ')-1

      return
      end

      subroutine setup_symbols_sdf()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
c
c sets the symbols and their default values for input to SDF.
c
c iflag decides whether this parameter can be changed
c a parameter can only be changed if iflag=0
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
      symbol_text(iref)='Maximum radius for spherical harmonic coefficients'
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
      symbol_text(iref)='Number of coefficients (program calculates this)'
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
      symbol_text(iref)='List of centre atom types. First will be origin.'
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
      symbol_text(iref)='Acceptance half-width (fraction of pair distance).'
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
      symbol_text(iref)='List of second atom types. First will be origin.'
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
      symbol_text(iref)='Acceptance half-width (fraction of pair distance).'
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
c now write the defaults to symbol_value
c
      do is=1,nsymbol
         symbol_value(is)=symbol_default(is)
      end do

c
c set the number of files - not used in SHARM
c
      nfiles=1

c
c setup the output prompt
c
      prompt='Setup sdf>'
      ilenprompt=index(prompt,'>')

c
c set up the data prompt - not used in SDF
c
      dataprompt=' '
      dataprompttext=' '
      ilendataprompt=0
c
c setup the extension for the input file
c
      extinp='.SDF.dat'
      ilenextinp=index(extinp,' ')-1

      return
      end

      subroutine setup_symbols_plot2d()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
c
c sets the symbols and their default values for input to clusters.
c
c
c set up the symbol table for the plot2d input file
c
c iflag decides whether this parameter can be changed
c a parameter can only be changed if iflag=0
c
c define symbols and default values
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
     *//'determined from coefficients file'
      iflag(iref)=1
      symbol_default(iref)=' '

      iref=iref+1
      symbol(iref)='ident'
      symbol_type(iref)='i'
      symbol_text(iref)='= 0 for identical molecules,'
     *//' else 1 if different'
      iflag(iref)=0
      symbol_default(iref)=' '
      
      iref=iref+1
      symbol(iref)='setone'
      symbol_type(iref)='i'
      symbol_text(iref)='0 sets first coefficient to zero - normally 1'
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
     *//' (separated by spaces)'
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
      symbol_text(iref)='no. of plots along x- and y-axis [set at 1 1]'
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
      symbol_text(iref)='vary thetal (1), phil (2), phim(3), thetam(4), chim(5)'
      iflag(iref)=0
      symbol_default(iref)='1'

      iref=iref+1
      symbol(iref)='ph_th_ch'
      symbol_type(iref)='c'
      symbol_text(iref)='phil, phim, thetam, chim (degrees - separated by spaces'
      iflag(iref)=0
      symbol_default(iref)='0 0 0 0'

      iref=iref+1
      symbol(iref)='axespar'
      symbol_type(iref)='r'
      symbol_text(iref)='axes character size, line width and colour'
     *//' (separated by spaces)'
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
     *//' (separated by spaces)'
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
     *//' (separated by spaces)'
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
      symbol_text(iref)='Contours overlaid? - set to N - cannot be changed'
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

c Last symbol is the extension of the file containing the spherical harmonic coefficients
c Used only for setup, not in the plotting program itself.

      iref=iref+1
      symbol(iref)='extcoeff'
      symbol_type(iref)='c'
      symbol_text(iref)='extension of coefficients file'
      iflag(iref)=0
      symbol_default(iref)='.SHARM.h01'

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
      prompt='setup plot2d>'
      ilenprompt=index(prompt,'>')
c
c set up the data prompt
c
      dataprompt='circle'
      dataprompttext='jumps to the specified circle '
      ilendataprompt=index(dataprompt,' ')-1

c
c setup the extension for the input file
c
      extinp='.plot2d.txt'
      ilenextinp=index(extinp,' ')-1

      return
      end

      subroutine setup_symbols_plot3d()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
c
c sets the symbols and their default values for input to clusters.
c
c
c set up the symbol table for the plot3d input file
c
c iflag decides whether this parameter can be changed
c a parameter can only be changed if iflag=0
c
c define symbols and default values
c
      iref=0

      iref=iref+1
      symbol(iref)='shcoeffs'
      symbol_type(iref)='c'
      symbol_text(iref)='Name of file containing spherical harmonic coefficients'
      iflag(iref)=0
      symbol_default(iref)='<undefined>'

      iref=iref+1
      symbol(iref)='ncoeffs'
      symbol_type(iref)='i'
      symbol_text(iref)='no. of coefficients - '
     *//'determined from coefficients file'
      iflag(iref)=1
      symbol_default(iref)=' '

      iref=iref+1
      symbol(iref)='ident'
      symbol_type(iref)='i'
      symbol_text(iref)='= 0 for identical molecules,'
     *//' else 1 if different'
      iflag(iref)=0
      symbol_default(iref)=' '
      
      iref=iref+1
      symbol(iref)='setone'
      symbol_type(iref)='i'
      symbol_text(iref)='0 sets first coefficient to zero - normally 1'
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
      symbol_text(iref)='no. of plots along x- and y-axis [set at 1 1]'
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
      symbol_text(iref)='fractional isosurface level (-ve for absolute)'
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
      symbol_text(iref)='vary (thetal, phil) (1), (thetam, phim) (2), '
     *//'(thetam, chim) (3)'
      iflag(iref)=0
      symbol_default(iref)='1'

      iref=iref+1
      symbol(iref)='ph_th_ch'
      symbol_type(iref)='c'
      symbol_text(iref)='phim, thetam, chim (degrees - separated by spaces'
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
     *//' (separated by spaces)'
      iflag(iref)=0
      symbol_default(iref)='0.0 0 0'

      iref=iref+1
      symbol(iref)='rgbsph'
      symbol_type(iref)='r'
      symbol_text(iref)='red green blue fractions for sphere'
     *//' (separated by spaces)'
      iflag(iref)=0
      symbol_default(iref)='0.7 0.7 0.7'

      nsymboldata=iref
      nsymbolplus=nsymboldata-nsymboluptodata

      iref=iref+1
      symbol(iref)='axespar'
      symbol_type(iref)='r'
      symbol_text(iref)='axes character size, line width and colour'
     *//' (separated by spaces)'
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
     *//' (separated by spaces)'
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
     *//' (separated by spaces)'
      iflag(iref)=0
      symbol_default(iref)='0.8 0.8 1.0'

      iref=iref+1
      symbol(iref)='ishade'
      symbol_type(iref)='i'
      symbol_text(iref)='ishade (1-8): 0 means no shading,'
     *//' -ve means inverted shading'
      iflag(iref)=0
      symbol_default(iref)='-1'

      iref=iref+1
      symbol(iref)='rgbobj'
      symbol_type(iref)='r'
      symbol_text(iref)='red green blue fractions for object'
     *//' (separated by spaces)'
      iflag(iref)=0
      symbol_default(iref)='1 1 0'

      iref=iref+1
      symbol(iref)='lightcoord'
      symbol_type(iref)='i'
      symbol_text(iref)='(x,y,z) coordinates for light source'
     *//' (separated by spaces)'
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
      symbol_text(iref)='transparency of object (0=0%,1=25%,2=50%,3=75%)'
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
      symbol_text(iref)='rotation and elevation of viewing point (deg.)'
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

c Last symbol is the extension of the file containing the spherical harmonic coefficients
c Used only for setup, not in the plotting program itself.

      iref=iref+1
      symbol(iref)='extcoeff'
      symbol_type(iref)='c'
      symbol_text(iref)='extension of coefficients file'
      iflag(iref)=0
      symbol_default(iref)='.SHARM.h01'

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
      prompt='setup plot3d>'
      ilenprompt=index(prompt,'>')
c
c set up the data prompt
c
      dataprompt='sphere'
      dataprompttext='jumps to the specified sphere '
      ilendataprompt=index(dataprompt,' ')-1

c
c setup the extension for the input file
c
      extinp='.plot3d.txt'
      ilenextinp=index(extinp,' ')-1

      return
      end

      subroutine setup_symbols_triangles()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
c
c sets the symbols and their default values for input to tri-angles.
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
      symbol_text(iref)='Single atom3 types [1], or multiple atom3 types [2]?'
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
      prompt='setup triangles>'
      ilenprompt=index(prompt,'>')
c
c set up the data prompt
c
      dataprompt='triangles'
      dataprompttext='jumps to the specified triangle '
      ilendataprompt=index(dataprompt,' ')-1
c
c setup the extension for the input file
c
      extinp='.TRI.dat'
      ilenextinp=index(extinp,' ')-1

      return
      end

      subroutine setup_symbols_torangles()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
c
c sets the symbols and their default values for input to torangles.
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
      prompt='setup torangles>'
      ilenprompt=index(prompt,'>')
c
c set up the data prompt
c
      dataprompt='torangles'
      dataprompttext='jumps to the specified torangle '
      ilendataprompt=index(dataprompt,' ')-1
c
c setup the extension for the input file
c
      extinp='.TOR.dat'
      ilenextinp=index(extinp,' ')-1

      return
      end

      subroutine setup_symbols_voids()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
c
c sets the symbols and their default values for input to voids.
c
c
c set up the symbol table for the VOIDS input file
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
      prompt='setup voids>'
      ilenprompt=index(prompt,'>')
c
c set up the data prompt
c
      dataprompt='void'
      dataprompttext='Specifies a void '
      ilendataprompt=index(dataprompt,' ')-1
c
c setup the extension for the input file
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
      symbol(iref)='nsize'
      symbol_type(iref)='i'
      symbol_text(iref)='No. of histogram points (max 1000)'
      iflag(iref)=0
      symbol_default(iref)='1000'

      iref=iref+1
      symbol(iref)='nsumt'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of samples already accumulated'
      iflag(iref)=0
      symbol_default(iref)='-1'

      iref=iref+1
      symbol(iref)='nsamples'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of samples to take from this configuration'
      iflag(iref)=0
      symbol_default(iref)='1000'

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
      dataprompt='size'
      dataprompttext='Specifies a size to measure fluctuations '
      ilendataprompt=index(dataprompt,' ')-1
c
c setup the extension for the input file
c
      extinp='.FLUCTUATIONS.dat'
      ilenextinp=index(extinp,' ')-1

      return
      end

      subroutine setup_symbols_wtsn()
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

      symbol(iref)='datafile'
      symbol_type(iref)='c'
      symbol_text(iref)='Name of data file this applies to'
      iflag(iref)=0
      symbol_default(iref)='<undefined>'

      iref=iref+1
      symbol(iref)='output'
      symbol_type(iref)='i'
      symbol_text(iref)='Is the output to be per atom (1) or per molecule (2)'
      iflag(iref)=0
      symbol_default(iref)='1'

      iref=iref+1
      symbol(iref)='nsamples'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of samples (1,2, or 3)'
      iflag(iref)=0
      symbol_default(iref)='1'

      iref=iref+1
      symbol(iref)='ndist'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of components (determined from ato file)'
      iflag(iref)=1
      symbol_default(iref)='1'

      nsymboluptodata=iref

      iref=iref+1
      symbol(iref)='atom'
      symbol_type(iref)='c'
      symbol_text(iref)='Atom label - set from .ato file'
      iflag(iref)=1
      symbol_default(iref)='0'

      iref=iref+1
      symbol(iref)='iexchange'
      symbol_type(iref)='i'
      symbol_text(iref)='1 if this atom exchanges with other atoms, 0 otherwise'
      iflag(iref)=1
      symbol_default(iref)='0'

      iref=iref+1
      symbol(iref)='abundances'
      symbol_type(iref)='c'
      symbol_text(iref)='List of mass numbers and abundances. 0 for natural isotope'
      iflag(iref)=0
      symbol_default(iref)='0'

      iref=iref+1
      symbol(iref)='isubs'
      symbol_type(iref)='i'
      symbol_text(iref)='This atom isotope substituted? (1 = yes, 0 = no)'
      iflag(iref)=1
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
      prompt='setup wtsn>'
      ilenprompt=index(prompt,'>')

c
c set up the data prompt
c
      dataprompt='atom'
      dataprompttext='Specifies a size to measure fluctuations '
      ilendataprompt=index(dataprompt,' ')-1
c
c setup the extension for the input file
c
      extinp='.WTSN.dat'
      ilenextinp=index(extinp,' ')-1

      return
      end

      subroutine setup_symbols_wtsx()
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

      symbol(iref)='datafile'
      symbol_type(iref)='c'
      symbol_text(iref)='Name of data file this applies to'
      iflag(iref)=0
      symbol_default(iref)='<undefined>'

      iref=iref+1
      symbol(iref)='output'
      symbol_type(iref)='i'
      symbol_text(iref)='Is the output to be per atom (1) or per molecule (2)'
      iflag(iref)=0
      symbol_default(iref)='1'

      iref=iref+1
      symbol(iref)='nsamples'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of samples (1,2, or 3)'
      iflag(iref)=0
      symbol_default(iref)='1'

      iref=iref+1
      symbol(iref)='ndist'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of components (determined from ato file)'
      iflag(iref)=1
      symbol_default(iref)='1'

      nsymboluptodata=iref

      iref=iref+1
      symbol(iref)='atom'
      symbol_type(iref)='c'
      symbol_text(iref)='Atom label - set from .ato file'
      iflag(iref)=1
      symbol_default(iref)='0'

      iref=iref+1
      symbol(iref)='iexchange'
      symbol_type(iref)='i'
      symbol_text(iref)='1 if this atom exchanges with other atoms, 0 otherwise'
      iflag(iref)=1
      symbol_default(iref)='0'

      iref=iref+1
      symbol(iref)='abundances'
      symbol_type(iref)='c'
      symbol_text(iref)='List of mass numbers and abundances. 0 for natural isotope'
      iflag(iref)=0
      symbol_default(iref)='0'

      iref=iref+1
      symbol(iref)='isubs'
      symbol_type(iref)='i'
      symbol_text(iref)='This atom isotope substituted? (1 = yes, 0 = no)'
      iflag(iref)=1
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
      prompt='setup wtsn>'
      ilenprompt=index(prompt,'>')

c
c set up the data prompt
c
      dataprompt='atom'
      dataprompttext='Specifies a size to measure fluctuations '
      ilendataprompt=index(dataprompt,' ')-1
c
c setup the extension for the input file
c
      extinp='.WTSX.dat'
      ilenextinp=index(extinp,' ')-1

      return
      end


