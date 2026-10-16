      subroutine read_epsr16(nin,nout,ierrout)
c
c sets up the EPSR input files
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      character*20 ext
      character*256 text
      character*1 iaccept
      integer*4 ilenext
c
c set up the symbol numbers for the fnameato and fnamepcof symbols
c
      ifnameato=0
      ifnamepcof=0
c
c set the number of potential coefficient pairs to zero - this is used to
c signal whether a valid .pcof file has been read in
c
      npairspcof=0
c
c signal whether any atoms have been read in 
c
      ntype=0
c
c generate the name of the input file
c
c search for the word 'EPSR16' in the input stream, and use the next word to define the 
c input filename, if present
c
      iw=1
      do while(iw.le.nwords.and.sometext(ncf(iw):ncl(iw)).ne.'epsr16')
         iw=iw+1
      end do
      if(iw.ge.nwords) then
c
c no filename found so prompt for it from the user
c
         ext='.inp'
         ilenext=4
         nwords=0
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
499   format(/'changeato> No .inp file specified! - routine will exit'/)
               return
            endif
         end do
c100  format(a) 
         iw=0
      endif
      iw=iw+1
      call concatenate(iw,nwords,fname_inp,ilenfname_inp
     *,nwords,ncf,ncl,sometext)
c
c if this filename has a .ext we do nothing and will attempt to read it as is,
c warts and all!
c
      text=fname_inp
      ilenuptodot=index(fname_inp,'.')-1
c     write(6,*) fname_inp(1:ilenfname_inp),ilenfname_inp,ilenuptodot
      if(ilenuptodot.lt.0) then
c
c .ext not supplied - therefore we need to add on the correct one
c
         ilenuptodot=ilenfname_inp
         fname_inp=fname_inp(1:ilenuptodot)//'.inp'
         ilenfname_inp=ilenuptodot+4
      endif
c
c create the full name of the input file by adding on the directory name if not present
c
      call get_full_filename(fname_inp,ilenfname_inp,fnameinp,ilenfnameinp)
c
c open this file if possible
c
      open(nin,file=fnameinp,status='old',iostat=ierrout)
      if(ierrout.eq.0) then
c
c  If O.K. read in the data from the specified input file
c
            write(nout,102) fname_inp(1:ilenfname_inp)
102   format(/'setup_input_file> Reading input file: "',a,'"')
c
c find the number of characters in this filename up to the last dot
c
            text=fname_inp
            ilenuptodot=indexlastdot(text)-1
c
c read the name of the run - this should correspond to fname_inp but without the 
c last .extension. It will be used to construct the names of input
c and output files (except .ato and .pcof files which may be named independently
c from the name of the run)
c
            read(nin,*,iostat=ierrout) fname_run
            ilenfname_run=index(fname_run,' ')-1
c
c remove the .extension if supplied
c
c
c check that this name corresponds to filename (without the .extension) specified. 
c If not signal the error but carry on, using the correct name when writing the file
c
            if(fname_run(1:ilenfname_run).ne.fname_inp(1:ilenuptodot)) then
                  write(nout,104) 
104   format(/'setup_input_file> Run name in input file is different from '
     *,'filename specified')
                  fname_run=fname_inp(1:ilenuptodot)
                  ilenfname_run=ilenuptodot
            endif
c
c generate the long form of fname_run
c
            call get_full_filename(fname_run,ilenfname_run
     *,fnamerun,ilenfnamerun)
c
c signal that if an error occurred the attempt to open the input file
c was unsuccessful, but we can still read commands from the terminal
c
            if(ierrout.ne.0) then
               ierrout=-1
               return
            endif
c NFILE = number of data files to fit (9 maximum)
            call get_a_line(nin,10,nstsize,36,36,1)
            read(symbol_value(nsymboluptodata),*) nfiles
      write(6,*) nfiles
            do in=1,nfiles
c INFILE = name of data file for S(Q) data
               call get_a_line(nin,10,nstsize,37,37,in)
c FNAMEW = name of file containing weighting factors for this dataset
               call get_a_line(nin,10,nstsize,38,38,in)
c type of input file (1,2, 3 or 4, 5)
c     NRTYPE=1, Q, S(Q) ASCII values, 1 set per line
c     NRTYPE=2, Q, S(Q), E(Q) ASCII values, 1 set per line
c     NRTYPE=3,4 N, Q, S(Q), E(Q) ASCII values, 1 set per line, 4 blank lines
c(GENIE SHOW DATA format, NRTYPE=3 corresponds to histogram format, 4 = point
c mode)
c     NRTYPE=5, Q, S(Q), E(Q) ASCII values, 1 set per line, using GUDRUN .DCS format (histogram)
c RSHMIN = number of nodes in potential
c SZEROS is the Q=0 limit for this file (currently not used. In the past it has been used to 
c normalise g(r) to give the correct zero limit). 
c TWEAK is a simple factor to multiply the data by prior to any processing. Normally 1.0.
c EFILEREQ is the requested energy amplitude for this file - increase it to improve fit. 
c However these values will only work if PRESSREQ is .le. zero. If PRESSREQ is .gt. 0 then
c EFILEREQ is set to PRESSREQ. 
c ** However if EFILEREQ is set to zero, then it will not be set by PRESSREQ. Instead the
c simulation will calculate the structure factor for this dataset but will not refine the 
c empirical potential against these data.
c ** Note that it may not be possible to achieve the specified energies for particular 
c datasets, particularly if a dataset contains more than one partial structure factor. This
c is because the other datasets may be also making demands on each partial, so that the
c different datasets may be in conflict.
               call get_a_line(nin,10,nstsize,39,43,in)
            end do
c FEEDBACK = overall feedback factor. Must be =< 1.0. Should only be set to 1.0 if the
c matrix of weight coefficients is well conditioned.
c Can be increased or decreased to emphasize or de-emphasize the data.
c POTFAC decides whether potential will be refined. POTFAC=0.0 means no refinement,
c (but simulation will proceed as usual) otherwise refinement will take place
c PRESSREQ is the requested absolute energy (kJ/mole) of empirical potential.
c This can be used to control degree of fit and energy/pressure for all datasets together.
c SIZEFACTOR (should be >=1) is used to shrink molecules during initial refinement when
c molecular entanglement may occur. Advisable to set potfac to 0 if sizefactor > 1
c
            call get_a_line(nin,10,nstsize,1,1,1)
c NQ = No. of Q values (Maximum = 1000)
c QSTEP = requested step in  Q
            call get_a_line(nin,10,nstsize,5,6,1)
            read(symbol_value(5),*) nq
            if(nq.gt.mpts) nq=mpts
            write(symbol_value(5),*) nq
c IRESET = 1 to set EP to zero, else IRESET = 0.
c IINIT = 1 to zero all accumulators, else IINIT = 0. IINIT is automatically set whenever
c     IRESET = 1
            call get_a_line(nin,10,nstsize,7,7,1)
c will have to reset since EPSR17 requires different files and format to EPSR16
            ireset=1
            write(symbol_value(7),*) ireset
            write(symbol_value(8),*) ireset
c write these values to the appropriate symbol values
c NTIMES = no. of trials per iteration between potential calculations (typically 5)
c NITER = no. of iterations (typically 1)
c NSUMS is the number of previous iterations that have already been accumulated
c     If NSUMS is -1, then a running average is used (determined from INTER). 
            call get_a_line(nin,10,nstsize,9,11,1)
c INTRA is the number of intermolecular moves that are made before an INTRA-
c molecular move (molecular vibration) is made - typically INTRA=50 or 100
c INTER is the number of iterations used to form running averages - after calling 
c subroutine INITA, the value of NESUMT will be set to INTER
            call get_a_line(nin,10,nstsize,12,12,1)
c number of iterations to occur after intramolecular disordering before 
c accumulating energy and pressure averages.
c RHO no. density in scattering units per A**3 (this is determined from 
c the .ATO file - this value is written here for information only)
            call get_a_line(nin,10,nstsize,14,14,1)
c CELLST = step size in linear region of r - controls number of r values
            call get_a_line(nin,10,nstsize,15,15,1)
c resolution function of data (Gaussian) FWHM(Q) = FWHM + FWHMQ*Q
            call get_a_line(nin,10,nstsize,16,17,1)
c NSMOOp = defines whether a background based on Rmin will be subtracted
c from the data prior to fitting the potential
c (0 for no background subtraction)
            call get_a_line(nin,10,nstsize,18,18,1)
c read filename for molecular coordinates
            call get_a_line(nin,10,nstsize,19,19,1)
            call read_ato(nin,nout,19,1,ierrout)
c read filename for potential coefficients - this is not defined for EPSR16
c so invent a new name
            write(symbol_value(20),*) fname_run(1:ilenfname_run)//'.pcof'
            sometext=symbol_value(20)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            call concatenate(1,nwords,fname_pcof,ilenfname_pcof
     *,nwords,ncf,ncl,sometext)
c
c convert this to a full path name if necessary
c
            call get_full_filename(fname_pcof,ilenfname_pcof
     *,fnamepcof,ilenfnamepcof)

c write these values to the appropriate symbol values
c read minimum Q value for molecular fit
            call get_a_line(nin,10,nstsize,21,21,1)
            close(nin)
      else
            fname_run=fname_inp(1:ilenuptodot)
            ilenfname_run=ilenuptodot
            write(nout,103) fname_inp(1:ilenfname_inp)
103   format(/'setup_input_file> Problems with specified input file: "',a
     *,'"'//'   - will use default values')
c
c generate the long form of fname_run
c
            call get_full_filename(fname_run,ilenfname_run
     *,fnamerun,ilenfnamerun)
c
c signify that the input file could not be opened or read, but a run
c name has been successfully defined with the default values
c
            ierrout=-1
            return
      endif
      return
      end

      subroutine read_epsr17(nin,nout,ierrout)
c
c sets up the EPSR input files
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      character*20 ext
      character*1 iaccept
      character*256 text
      integer*4 ilenext
c
c set up the symbol numbers for the fnameato and fnamepcof symbols
c
      ifnameato=0
      ifnamepcof=0
c
c set the number of potential coefficient pairs to zero - this is used to
c signal whether a valid .pcof file has been read in
c
      npairspcof=0
c
c signal whether any atoms have been read in 
c
      ntype=0
c
c generate the name of the input file
c
c search for the word 'EPSR17' in the input stream, and use the next word to define the 
c input filename, if present
c
      iw=1
      do while(iw.le.nwords.and.sometext(ncf(iw):ncl(iw)).ne.'epsr17')
         iw=iw+1
      end do
      if(iw.ge.nwords) then
c
c no filename found so prompt for it from the user
c
         ext='.inp'
         ilenext=4
         nwords=0
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
499   format(/'changeato> No .inp file specified! - routine will exit'/)
               return
            endif
         end do
c100  format(a) 
         iw=0
      endif
      iw=iw+1
      call concatenate(iw,nwords,fname_inp,ilenfname_inp
     *,nwords,ncf,ncl,sometext)
c
c if this filename has a .ext we do nothing and will attempt to read it as is,
c warts and all!
c
      text=fname_inp
      ilenuptodot=index(fname_inp,'.')-1
c     write(6,*) fname_inp(1:ilenfname_inp),ilenfname_inp,ilenuptodot
      if(ilenuptodot.lt.0) then
c
c .ext not supplied - therefore we need to add on the correct one
c
         ilenuptodot=ilenfname_inp
         fname_inp=fname_inp(1:ilenuptodot)//'.inp'
         ilenfname_inp=ilenuptodot+4
      endif
c
c create the full name of the input file by adding on the directory name if not present
c
      call get_full_filename(fname_inp,ilenfname_inp,fnameinp,ilenfnameinp)
c
c open this file if possible
c
      open(nin,file=fnameinp,status='old',iostat=ierrout)
      if(ierrout.eq.0) then
c
c  If O.K. read in the data from the specified input file
c
            write(nout,102) fname_inp(1:ilenfname_inp)
102   format(/'setup_input_file> Reading input file: "',a,'"')
c
c find the number of characters in this filename up to the last dot
c
            text=fname_inp
            ilenuptodot=indexlastdot(text)-1
c
c read the name of the run - this should correspond to fname_inp but without the 
c last .extension. It will be used to construct the names of input
c and output files (except .ato and .pcof files which may be named independently
c from the name of the run)
c
            read(nin,*,iostat=ierrout) fname_run
            ilenfname_run=index(fname_run,' ')-1
c
c remove the .extension if supplied
c
c
c check that this name corresponds to filename (without the .extension) specified. 
c If not signal the error but carry on, using the correct name when writing the file
c
            if(fname_run(1:ilenfname_run).ne.fname_inp(1:ilenuptodot)) then
                  write(nout,104) 
104   format(/'setup_input_file> Run name in input file is different from '
     *,'filename specified')
                  fname_run=fname_inp(1:ilenuptodot)
                  ilenfname_run=ilenuptodot
            endif
c
c generate the long form of fname_run
c
            call get_full_filename(fname_run,ilenfname_run
     *,fnamerun,ilenfnamerun)
c
c signal that if an error occurred the attempt to open the input file
c was unsuccessful, but we can still read commands from the terminal
c
      write(6,*) ierrout
            if(ierrout.ne.0) then
               ierrout=-1
               return
            endif
c NFILE = number of data files to fit (9 maximum)
            call get_a_line(nin,10,nstsize,36,36,1)
            read(symbol_value(nsymboluptodata),*) nfiles
      write(6,*) nfiles
            do in=1,nfiles
c INFILE = name of data file for S(Q) data
               call get_a_line(nin,10,nstsize,37,37,in)
c FNAMEW = name of file containing weighting factors for this dataset
               call get_a_line(nin,10,nstsize,38,38,in)
c type of input file (1,2, 3 or 4, 5)
c     NRTYPE=1, Q, S(Q) ASCII values, 1 set per line
c     NRTYPE=2, Q, S(Q), E(Q) ASCII values, 1 set per line
c     NRTYPE=3,4 N, Q, S(Q), E(Q) ASCII values, 1 set per line, 4 blank lines
c(GENIE SHOW DATA format, NRTYPE=3 corresponds to histogram format, 4 = point
c mode)
c     NRTYPE=5, Q, S(Q), E(Q) ASCII values, 1 set per line, using GUDRUN .DCS format (histogram)
c RSHMIN = number of nodes in potential
c SZEROS is the Q=0 limit for this file (currently not used. In the past it has been used to 
c normalise g(r) to give the correct zero limit). 
c TWEAK is a simple factor to multiply the data by prior to any processing. Normally 1.0.
c EFILEREQ is the requested energy amplitude for this file - increase it to improve fit. 
c However these values will only work if PRESSREQ is .le. zero. If PRESSREQ is .gt. 0 then
c EFILEREQ is set to PRESSREQ. 
c ** However if EFILEREQ is set to zero, then it will not be set by PRESSREQ. Instead the
c simulation will calculate the structure factor for this dataset but will not refine the 
c empirical potential against these data.
c ** Note that it may not be possible to achieve the specified energies for particular 
c datasets, particularly if a dataset contains more than one partial structure factor. This
c is because the other datasets may be also making demands on each partial, so that the
c different datasets may be in conflict.
               call get_a_line(nin,10,nstsize,39,43,in)
            end do
c FEEDBACK = overall feedback factor. Must be =< 1.0. Should only be set to 1.0 if the
c matrix of weight coefficients is well conditioned.
c Can be increased or decreased to emphasize or de-emphasize the data.
c POTFAC decides whether potential will be refined. POTFAC=0.0 means no refinement,
c (but simulation will proceed as usual) otherwise refinement will take place
c PRESSREQ is the requested absolute energy (kJ/mole) of empirical potential.
c This can be used to control degree of fit and energy/pressure for all datasets together.
c SIZEFACTOR (should be >=1) is used to shrink molecules during initial refinement when
c molecular entanglement may occur. Advisable to set potfac to 0 if sizefactor > 1
c
            call get_a_line(nin,10,nstsize,1,4,1)
c NQ = No. of Q values (Maximum = 1000)
c QSTEP = requested step in  Q
            call get_a_line(nin,10,nstsize,5,6,1)
            read(symbol_value(5),*) nq

            if(nq.gt.mpts) nq=mpts
            write(symbol_value(5),*) nq

c IRESET = 1 to set EP to zero, else IRESET = 0.
c IINIT = 1 to zero all accumulators, else IINIT = 0. IINIT is automatically set whenever
c     IRESET = 1
            call get_a_line(nin,10,nstsize,7,8,1)
            ireset=1
            write(symbol_value(7),*) ireset
            write(symbol_value(8),*) ireset
c NTIMES = no. of trials per iteration between potential calculations (typically 5)
c NITER = no. of iterations (typically 1)
c NSUMS is the number of previous iterations that have already been accumulated
c     If NSUMS is -1, then a running average is used (determined from INTER). 
            call get_a_line(nin,10,nstsize,9,11,1)
c INTRA is the number of intermolecular moves that are made before an INTRA-
c molecular move (molecular vibration) is made - typically INTRA=50 or 100
c INTER is the number of iterations used to form running averages - after calling 
c subroutine INITA, the value of NESUMT will be set to INTER
            call get_a_line(nin,10,nstsize,12,13,1)
c number of iterations to occur after intramolecular disordering before 
c accumulating energy and pressure averages.
c RHO no. density in scattering units per A**3 (this is determined from 
c the .ATO file - this value is written here for information only)
            call get_a_line(nin,10,nstsize,14,14,1)
c CELLST = step size in linear region of r - controls number of r values
            call get_a_line(nin,10,nstsize,15,15,1)
c resolution function of data (Gaussian) FWHM(Q) = FWHM + FWHMQ*Q
            call get_a_line(nin,10,nstsize,16,17,1)
c NSMOOp = defines whether a background based on Rmin will be subtracted
c from the data prior to fitting the potential
c (0 for no background subtraction)
            call get_a_line(nin,10,nstsize,18,18,1)
c read filename for molecular coordinates
            call get_a_line(nin,10,nstsize,19,19,1)
            call read_ato(nin,nout,19,1,ierrout)
c read filename for potential coefficients
            call get_a_line(nin,10,nstsize,20,20,1)
      write(nout,*) symbol_value(20)
         sometext=symbol_value(20)
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         call concatenate(1,nwords,fname_pcof,ilenfname_pcof
     *,nwords,ncf,ncl,sometext)
c
c convert this to a full path name if necessary
c
         call get_full_filename(fname_pcof,ilenfname_pcof
     *,fnamepcof,ilenfnamepcof)

c read minimum Q value for molecular fit
            call get_a_line(nin,10,nstsize,21,21,1)
            close(nin)
      else
            fname_run=fname_inp(1:ilenuptodot)
            ilenfname_run=ilenuptodot
            write(nout,103) fname_inp(1:ilenfname_inp)
103   format(/'setup_input_file> Problems with specified input file: "',a
     *,'"'//'   - will use default values')
c
c generate the long form of fname_run
c
            call get_full_filename(fname_run,ilenfname_run
     *,fnamerun,ilenfnamerun)
c
c signify that the input file could not be opened or read, but a run
c name has been successfully defined with the default values
c
            ierrout=-1
            return
      endif
      return
      end

      subroutine setup_symbols_epsr()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
      include 'iresetnumbers.inc'
c
c sets the symbols and their default values for input to EPSR.
c
c
c set up the symbol table for the EPSR .inp file
c
c iflag decides whether this parameter can be changed
c a parameter can only be changed if ireset >= iflag
c iflag=1 means all accumulators are reset
c iflag=2 means potential coefficients are zerod as well
c
      
      iref=0
      
      iref=iref+1
      symbol(iref)='feedback'
      symbol_type(iref)='r'
      symbol_text(iref)='Confidence factor - should be < 1. [0.8]'
      iflag(iref)=-1
      symbol_default(iref)='0.9'
      
      iref=iref+1
      symbol(iref)='potfac'
      symbol_type(iref)='r'
      symbol_text(iref)='1.0 to enable potential refinement, 0.0 to inhibit'
      iflag(iref)=0
      symbol_default(iref)='1.0'
      
      iref=iref+1
      symbol(iref)='ereq'
      symbol_type(iref)='r'
      symbol_text(iref)='Overall requested energy amplitude'
      iflag(iref)=0
      symbol_default(iref)='0.0'
      
      iref=iref+1
      symbol(iref)='num_threds'
      symbol_type(iref)='i'
      symbol_text(iref)='No. of parallel threads to be used, (0 to let program choose) [1]'
      iflag(iref)=0
      symbol_default(iref)='1'
      
      iref=iref+1
      symbol(iref)='porodpwr'
      symbol_type(iref)='r'
      symbol_text(iref)='Q-value and power for Porod scattering [0]'
      iflag(iref)=0
      symbol_default(iref)='0 0'
      
      iref=iref+1
      symbol(iref)='nmolcell'
      symbol_type(iref)='i'
      symbol_text(iref)='Average number of molecules in a cell [10]'
      iflag(iref)=0
      symbol_default(iref)='10'
      
      iref=iref+1
      symbol(iref)='control-p'
      symbol_type(iref)='c'
      symbol_text(iref)='Control pressure [a, b and/or c axes control plus value plus step, none 0 0]'
      iflag(iref)=0
      symbol_default(iref)='none 0 0'
      
      iref=iref+1
      symbol(iref)='ref_intra'
      symbol_type(iref)='r'
      symbol_text(iref)='Weighting on EP and Coulomb terms for intra-molecular structure [0 0]'
      iflag(iref)=0
      symbol_default(iref)='0 0'
      
      iref=iref+1
      symbol(iref)='sizefactor'
      symbol_type(iref)='r'
      symbol_text(iref)='Multiplying factor for box dimension, decline rate, and threshold. [1.0 0.9 0.0]'
      iflag(iref)=0
      symbol_default(iref)='1.0 0.9 0.0'
      
      iref=iref+1
      symbol(iref)='nq'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of Q values. [600]'
      iflag(iref)=-1
      symbol_default(iref)='600'
      
      iref=iref+1
      symbol(iref)='qstep'
      symbol_type(iref)='r'
      symbol_text(iref)='Size of Q step [1/A]. [0.05]'
      iflag(iref)=-1
      symbol_default(iref)='0.05'
      
      iref=iref+1
      symbol(iref)='ireset'
      symbol_type(iref)='i'
      symbol_text(iref)='1: complete reset; 2: sets the Empirical Potential to zero'
      iflag(iref)=0
      symbol_default(iref)='1'

      iresetref=iref
      
      iref=iref+1
      symbol(iref)='iinit'
      symbol_type(iref)='i'
      symbol_text(iref)='Sets accumulators to zero. Recalculates r and Q. [1]'
      iflag(iref)=0
      symbol_default(iref)='1'
      
      iinitref=iref

      iref=iref+1
      symbol(iref)='ntimes'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of MC cycles between potential refinements. [5]'
      iflag(iref)=0
      symbol_default(iref)='5'
      
      iref=iref+1
      symbol(iref)='niter'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of potential refinements before exitting. [1]'
      iflag(iref)=0
      symbol_default(iref)='1'
      
      iref=iref+1
      symbol(iref)='nsumt'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of iterations already accumulated. [-1 with reset]'
      iflag(iref)=0
      symbol_default(iref)='-1'
      
      iref=iref+1
      symbol(iref)='intra'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of iterations between molecule shakes. [100]'
      iflag(iref)=0
      symbol_default(iref)='100'
      
      iref=iref+1
      symbol(iref)='rotfreq'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of iterations between internal rotation moves. [5]'
      iflag(iref)=0
      symbol_default(iref)='5'
      
      iref=iref+1
      symbol(iref)='inter'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of iterations in running averages. [5]'
      iflag(iref)=0
      symbol_default(iref)='5'
      
      iref=iref+1
      symbol(iref)='rho'
      symbol_type(iref)='r'
      symbol_text(iref)='Atomic number density - will be derived from .ato file'
      iflag(iref)=1
      symbol_default(iref)='0.1'
      
      iref=iref+1
      symbol(iref)='cellst'
      symbol_type(iref)='r'
      symbol_text(iref)='Size of r step [A]. [0.03]'
      iflag(iref)=-1
      symbol_default(iref)='0.03'
      
      iref=iref+1
      symbol(iref)='rmaxgr'
      symbol_type(iref)='r'
      symbol_text(iref)='Range of g(r) and F.T. (0.0 will use half the cell box) [0.0]'
      iflag(iref)=-1
      symbol_default(iref)='0.0'
      
      iref=iref+1
      symbol(iref)='ngrsamples'
      symbol_type(iref)='i'
      symbol_text(iref)='Requested no. of origin molecules to sample g(r). (0 will use 1000 molecules) [0]'
      iflag(iref)=0
      symbol_default(iref)='0'
      
      iref=iref+1
      symbol(iref)='fwhm'
      symbol_type(iref)='r'
      symbol_text(iref)='Resolution width - Q independent term. [0.0]'
      iflag(iref)=0
      symbol_default(iref)='0.0'
      
      iref=iref+1
      symbol(iref)='fwhmq'
      symbol_type(iref)='r'
      symbol_text(iref)='Resolution width - Q dependent term. [0.02 for SLS]'
      iflag(iref)=0
      symbol_default(iref)='0.02'
      
      iref=iref+1
      symbol(iref)='nsmoop'
      symbol_type(iref)='i'
      symbol_text(iref)='1 means background subtraction is ON, 0 means OFF'
      iflag(iref)=0
      symbol_default(iref)='1'
      
      iref=iref+1
      symbol(iref)='fnameato'
      symbol_type(iref)='c'
      symbol_text(iref)='Name of .ato file'
      iflag(iref)=0
      symbol_default(iref)='<undefined>'
      
      iref=iref+1
      symbol(iref)='fnamepcof'
      symbol_type(iref)='c'
      symbol_text(iref)='Name of potential coefficients file.'
      iflag(iref)=0
      symbol_default(iref)='<undefined>'
      
      iref=iref+1
      symbol(iref)='revlorch'
      symbol_type(iref)='r'
      symbol_text(iref)='Broadening factor in Q space. [0.0 0.0]'
      iflag(iref)=0
      symbol_default(iref)='0.0 0.0'
      
      iref=iref+1
      symbol(iref)='qwidthqmax'
      symbol_type(iref)='r'
      symbol_text(iref)='Broadening and maximum Q for Bragg peak calculation'
      iflag(iref)=0
      symbol_default(iref)='0.01 0'

      iref=iref+1
      symbol(iref)='mplicities'
      symbol_type(iref)='i'
      symbol_text(iref)='No. of unit cells along a, b and c for Bragg peak calculation'
      iflag(iref)=0
      symbol_default(iref)='1 1 1'

      iref=iref+1
      symbol(iref)='hklqmin'
      symbol_type(iref)='r'
      symbol_text(iref)='Minimum value of qhkl to be used, minimum radius for Bragg g(r), and Debye-Waller factor'
      iflag(iref)=0
      symbol_default(iref)='0.0 0.0 0.0'

      iref=iref+1
      symbol(iref)='diffuse'
      symbol_type(iref)='i'
      symbol_text(iref)='(h,k,l) patterns for diffuse scattering calculation [0]'
      iflag(iref)=0
      symbol_default(iref)='0'

      iref=iref+1
      symbol(iref)='rejrate'
      symbol_type(iref)='r'
      symbol_text(iref)='Rejection rate [0.75]'
      iflag(iref)=0
      symbol_default(iref)='0.75'

      iref=iref+1
      symbol(iref)='qmin'
      symbol_type(iref)='r'
      symbol_text(iref)='Qmin for Fourier transforms and for potential fits. [0.05 0.0]'
      iflag(iref)=0
      symbol_default(iref)='0.05 0.0'
      
      iqminref=iref

      iref=iref+1
      symbol(iref)='roverlap'
      symbol_type(iref)='r'
      symbol_text(iref)='Minimum allowed intermolecular separation between two atoms'
      iflag(iref)=-3
      symbol_default(iref)='0.0'
      iroverlap=iref
      
      iref=iref+1
      symbol(iref)='rminfac'
      symbol_type(iref)='r'
      symbol_text(iref)='Factor to set the minimum separation between pairs. [0.0]'
      iflag(iref)=-3
      symbol_default(iref)='0.0'
      irminfac=iref     

      iref=iref+1
      symbol(iref)='rminpt'
      symbol_type(iref)='r'
      symbol_text(iref)='Radius at which potential truncation begins.'
      iflag(iref)=-3
      symbol_default(iref)='9'
      
      iref=iref+1
      symbol(iref)='rmaxpt'
      symbol_type(iref)='r'
      symbol_text(iref)='Radius at which potential truncation goes to 0.0.'
      iflag(iref)=-3
      symbol_default(iref)='12'
      
      iref=iref+1
      symbol(iref)='rbroad'
      symbol_type(iref)='r'
      symbol_text(iref)='Controls potential decay. [0.0]'
      iflag(iref)=1
      symbol_default(iref)='0.0'
      
      iref=iref+1
      symbol(iref)='expecf'
      symbol_type(iref)='r'
      symbol_text(iref)='Potential decay for short distance repulsive term and gradient limit. [0.3 0]'
      iflag(iref)=0
      symbol_default(iref)='0.3 0'
      
      iref=iref+1
      symbol(iref)='reppottype'
      symbol_type(iref)='c'
      symbol_text(iref)='Repulsive potential type: exponential or harmonic. [exponential]'
      iflag(iref)=0
      symbol_default(iref)='exponential'
      
      iref=iref+1
      symbol(iref)='addpottype'
      symbol_type(iref)='c'
      symbol_text(iref)='Additional potential type: Gaussian or modmorse. [Gaussian]'
      iflag(iref)=0
      symbol_default(iref)='Gaussian'
      
      iref=iref+1
      symbol(iref)='rcharge'
      symbol_type(iref)='r'
      symbol_text(iref)='Calculates energy due to molecular polarisation. [0.0]'
      iflag(iref)=0
      symbol_default(iref)='0.0'
      
      iref=iref+1
      symbol(iref)='power'
      symbol_type(iref)='r'
      symbol_text(iref)='Repulsive power in Lennard-Jones potential. [12]'
      iflag(iref)=0
      symbol_default(iref)='12'
      
      iref=iref+1
      symbol(iref)='ncoeffp'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of coefficients used to define the EP. Set by program.'
      iflag(iref)=1
      symbol_default(iref)='120'
      incoeffp=iref
      
      iref=iref+1
      symbol(iref)='pdmax'
      symbol_type(iref)='r'
      symbol_text(iref)='Maximum distance of Empirical Potential. Set by program.'
      iflag(iref)=1
      symbol_default(iref)='12.0'
      
      iref=iref+1
      symbol(iref)='pdstep'
      symbol_type(iref)='r'
      symbol_text(iref)='Spacing between coefficients in r. Set by program.'
      iflag(iref)=1
      symbol_default(iref)='0.1'
      
      iref=iref+1
      symbol(iref)='npitss'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of steps for refining the potential. [1000]'
      iflag(iref)=0
      symbol_default(iref)='1000'
      
      iref=iref+1
      symbol(iref)='paccept'
      symbol_type(iref)='r'
      symbol_text(iref)='Acceptance factor for potential refinement. [0.0005]'
      iflag(iref)=1
      symbol_default(iref)='0.0005'
      
      iref=iref+1
      symbol(iref)='Gaussian'
      symbol_type(iref)='c'
      symbol_text(iref)='Select T for Gaussian representation of EP. Otherwise Poisson. [F]'
      iflag(iref)=-3
      symbol_default(iref)='F'
      
      iref=iref+1
      symbol(iref)='psigma2'
      symbol_type(iref)='r'
      symbol_text(iref)='Width for potential function. [0.02]'
      iflag(iref)=0
      symbol_default(iref)='0.02'

      iref=iref+1
      symbol(iref)='refpotfac'
      symbol_type(iref)='r'
      symbol_text(iref)='Factor to apply to reference potential. [1.0]'
      iflag(iref)=0
      symbol_default(iref)='1.0'

      iref=iref+1
      symbol(iref)='ndata'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of data files to be fit by EPSR'
      iflag(iref)=-1
      symbol_default(iref)='1'
      
      nsymboluptodata=iref

      iref=iref+1
      symbol(iref)='datafile'
      symbol_type(iref)='c'
      symbol_text(iref)='Name of data file to be fit'
      iflag(iref)=0
      symbol_default(iref)='<undefined>'
      
      iref=iref+1
      symbol(iref)='wtsfile'
      symbol_type(iref)='c'
      symbol_text(iref)='Name of weights file for this data set'
      iflag(iref)=-1
      symbol_default(iref)='<undefined>'
      
      iref=iref+1
      symbol(iref)='nrtype'
      symbol_type(iref)='i'
      symbol_text(iref)='Data type - see User Manual for more details'
      iflag(iref)=0
      symbol_default(iref)='5'
      
      iref=iref+1
      symbol(iref)='rshmin'
      symbol_type(iref)='r'
      symbol_text(iref)='Minimum radius [A] - used for background subtraction'
      iflag(iref)=0
      symbol_default(iref)='0.7'
      
      iref=iref+1
      symbol(iref)='szeros'
      symbol_type(iref)='r'
      symbol_text(iref)='Zero limit - 0 means use first data point for Q=0'
      iflag(iref)=0
      symbol_default(iref)='0.0'
      
      iref=iref+1
      symbol(iref)='tweak'
      symbol_type(iref)='r'
      symbol_text(iref)='Scaling factor for this data set. [1.0]'
      iflag(iref)=0
      symbol_default(iref)='1.0'
      
      iref=iref+1
      symbol(iref)='efilereq'
      symbol_type(iref)='r'
      symbol_text(iref)='Requested energy amplitude for this data set [1.0]'
      iflag(iref)=0
      symbol_default(iref)='1.0'

      iref=iref+1
      symbol(iref)='porodamp'
      symbol_type(iref)='r'
      symbol_text(iref)='Amplitude of Porod scattering for this dataset [0.0]'
      iflag(iref)=0
      symbol_default(iref)='0.0'

      nsymboldata=iref
      nsymbol=iref
      nsymbolplus=nsymboldata-nsymboluptodata

c
c write the defaults to symbol_value
c
      do is=1,nsymbol
         symbol_value(is)=symbol_default(is)
      end do

      prompt='Setup epsr>'
      ilenprompt=index(prompt,'>')

c
c set data prompt text
c
      dataprompt='data'
      ilendataprompt=4
      dataprompttext='jumps to the specified dataset'

c
c setup the extension for the input file
c
      extinp='.EPSR.inp'
      ilenextinp=index(extinp,' ')-1

c
c default number of files
c
      nfiles=1

      return
      end

