      program epsr_standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call runepsr17()

      end

      subroutine runepsr17()

c EPSR17_3.FOR

c Version 17_3 of EPSR (29/01/2004) for fitting in Q space.

c     Standalone version of EPSR using only ASCII input and output files
c     (file format for input files can be compatible with output from
c     GENIE Show Data command, as well as GUDRUN .mdcs files)

c     Based on EPSR15X, which can be used for an
c     arbitrary unit cell shape, with modified input and output
c     routines to force the simulation box to be cubic.

c     Please refer to AKS in case of questions

c     This (PC) version uses a simplified .INP file, and reads all the Empirical Potential information
c     coefficients from the .PCOF file

c     This version uses input parameters for each dataset to generate the Empirical
c     potential using a series of power-exponential functions in r-space
c     which are fit directly in Q-space. A broadened version of this
c     expansion is then used to generate the perturbation to the empirical
c     potential at each iteration.

c     17/12/2002 Attempt to use a smart inversion method to go from Empirical Potential for each
c     dataset to the Empirical Potential for each atom pair.

c     24/03/2003 In addition estimates the site-site partial structure factors for the data and
c     compares with the simulation.

c     24/10/2003 Variation of method: simulation is treated as starting point for
c     refinement of potentials. (Initial simulation should be with feedback=0.0)
c     Then feedback > 0.0 introduces the data as a perturbation to this. However feedback
c     should not be made equal to unity, unless the matrix of weights is well conditioned
c     otherwise the inverse of the weights matrix is indeterminate.
c     feedback is not allowed to be larger 1.0, otherwise it will be set to 0.1.
c     The overall amplitude of the empirical potential is controlled by the value of
c     PRESSREQ. The EP will grow with each iteration until either the required
c     perturbation goes to zero or the abosolute energy reaches its specified value.
c
c     Individual datasets can be given greater or lesser weighting by specifying the
c     factor (EFACM) greater than or less than 1.0 at the end of the line of numbers
c     for each dataset.

c     12/10/2003 Fix maximum of potential to specified value, not automatically set to
c     half the box size.
c
c     Have introduced neighbour lists and energy array to reduce time spent calculating
c     energies

c     Jan 2004. Introduced a size factor for shrinking molecules during the initial
c     simulations when molecules may overlap and could become entangled if they have a
c     ring structure. Molecules are made smaller by giving a size factor > 1. E.g.
c     SIZEFACTOR=10 means molecules are shrunk to 1/10 of their size compared to
c     intermolecular distances.

c     Jan 2004. Each datafile has now has a factor that can be specified in the .INP file
c     on the line of numbers corresponding to that file.

c     29/01/2004. Modified the input routines so that the program generates the empirical
c     potential from site-site potential coefficients read in from the .PCOF file. The
c     same file is written out at the end containing the new coefficients after
c     modification by the data.
c
c     Shortly after input the site-site potential coefficients are converted to the data
c     potential coefficients (as used in the earlier versions of EPSR) and these are used
c     to generate and modify the EP as this is the most efficient way of performing these
c     actions, particularly for systems where the no. of atom pairs is larger than the
c     number of available datasets, which is the most common scenario. As a check on
c     output the .INPA file will continue to show the data potential coefficients for the
c     time being, but these will no longer be used by the program on input.
c
c     The main reason for this change is so that coefficients generated in one simulation
c     can be used as is (and may be modified) to generate the EP in another. Thus the
c     potential generation is kept completely independent from any parameters used
c     within a particular simulation box, (such as size of box for example).
c
c     IMPORTANT: The potential coefficients are set to zero after IRESET=1. Therefore if
c     running more than one simulation from the same .PCOF file, then only the first
c     simulation should have IRESET=1 initially. All the others should have IINIT=1 to
c     reset all the accumulators, and recalculate the inversion matrix, but not reset the
c     the potential. Otherwise each simulation when it starts will reset the potential.

c     27/06/2006. Program now writes the current version number into the .OUT file.

c     29/06/2006. Today's mod is to exclude the repulsive potential from the calculation
c     of the absolute empirical energy. This is because the repulsive potential can often
c     dominate the absolute magnitude, so that if the repulsive potential increases to
c     push atoms futher apart it can then exceed the maximum absolute energy, causing the
c     empirical absolute energy to decrease, and making the fit worse.
c
c     Also corrected a typo when calculating the energy written to the .OUT file

c     27/09/2006 Modified molout so that it lists the atom number in a molecule
c     alongside its symbol. This is not read on input

c     19/10/2006 Modified ptgnsq so that the repulsive potential is generated even when
c     potfac = 0. To avoid this, set the allowed minimum distances to zero.
c     Also moloutc now prints out the molecule number after the list of coordinates. This
c     is not read on input

c     27/10/2006 Modified running coordination number so that both ab and ba pairs are
c     output.

c     2006-11-12 Previous versions modified fmolec1 so that bonded neighbours
c     are always included in the neighbour list. However this introduced a bug
c     by not cforwarding the correct atom types to the routine which calculates c   the intra-molecular energy.

c     Also a minimum stepsize for translational moves has been set to 1/100th of
c     the maximum stepsize.

c     Also value of 'version' is defined by an include file

c     2006-12-04 Hopefully have corrected pressure calculation, which gave the
c     wrong sign on the pressure for the contribution from the potential.

c     2007-06-13 Worries about the repulsive potential interfering with the EP
c     when potfac is set to zero or very small. Solution: make the increments to
c     the repulsive potential proportional to potfac.

c     However this makes it impossible to run the simulation with the repulsive
c     potential on its own. Next solution: we can set the absolute energy of the
c     EP to zero (ereq), but keep potfac > 0 to run with reference and repulsive
c     potentials alone.

c     However currently setting ereq to zero, has effect of
c     transferring weighting to individual weights, efilereq. However this
c     method of accentuating particular datasets is not very effective, so have
c     resorted to using efilereq to control the relative weight of each dataset
c     in calculating the inversion matrix.

c     Also need to ensure in moloutc that if ROT groups are defined but internal
c     rotational step is set zero, it sends a warning message.

c     ALSO fixed the bug that if X-ray form factor coefficients are negative they
c     are summed correctly!!

c     2007-09-30. Putting in Gaussians into EP as well as exponential repulsive
c     potential. Read in from .pcof file, as for other coefficients.


      use timer
      use cellneighbourlist !Stores and processes the cell neighbour list
      use diffuse_mod
      use threebody_mod
      use molecom_mod

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'epsr16.inc'
      include 'epsrinput_par.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      include 'weights.inc'
      include 'X-rayf0.inc'
      include 'iresetnumbers.inc'
      include 'version.inc'
      include 'lattice.inc'

      real z(mpz),Q(mpts),datai(mpts,mfile),errs(mpts,mfile),ufile(mpts,mfile)
     1,ugdif(mpz,mfile),ugdifsmo(mpz,mfile)
     1,qwind(2,mpts,mpairs),rfac(mpairs)
     1,datat(mpts,mfile),backsc(mpts,mfile),wav(2,mpz,mpairs)
     5,errsav(mpts,mpairs),soqism(mpts,mfile)
     1,rshmin(mfile),szeros(mfile),szerostemp(mfile)
     1,rfacd(mfile),potnoi(mfile)
     1,energf(mfile),pressf(mfile),efac(mpairs),tweak(mfile)
     1,datap(mpts,mpairs),datagr(mpz,mpairs),datagr2(mpz,mpairs)
     1,dataigr(mpz,mfile)
     1,grlim(mpairs),grzerolim(mfile),efilereq(mfile)
     1,stepmisave,steprisave,rmaxgrsave,energabs,energthreeb,pressthreeb,ftethered,fcoordconstraints
     1,epccoeff,epcfrac,ereqstepread,ereqbias,epcthresh,epcthreshuse,epcrmin,epcfac,gortemp,xavstemp,datagrtemp,rtempstep
     1,epcup,epcdn,bias,ereqmin,ereqmax,fractemp,lastqf,nextqf

      integer navmaxparab,nparabstart,nparabstartfactor,nsumparab

      real QD(mqd),DATA(mqd),err(mqd),gofzism(mpz,mfile)
     *,gz2ism(mpz,mfile)

c Temporary store of simulated intermolecular partial structure factors

      real soqintertemp(mpts,mpairs),dts(mpts),dti(mpts)

c Factor used to limit the absolute energy of the empirical potential

      real elimitfac,qbraggfactor

c nrtype is the type of dataset being read - see below.
c nqminp and nqmaxp are used in the Fourier transform of the estimated partial structure
c factors.

      integer nrtype(mfile),nqmind(mfile),nqmaxd(mfile),nqminp(mpairs),nqmaxp(mpairs)
     *,num_threads,nmolcellin,ngrsamples,ngrsamplessave,rotfreq,delnz,j1,j2,j,ilensymbol_value

      character*256 infile,fname,infili(mfile),outputfile
     1,fnamew(mfile),fnamedata(mfile)
     *,fnamewts(mfile)

      character*2 ext1

      character*10 extfile,hashsearch

      character*29 titlepsf(mpairs),titledat(mfile)
      character*256 formattitle

      character*14 titler,titleq
      character*4 atom1,atom2

      character*20 text

      real*8 fc1,fc2,add1(2),add2(2),sumtot1(2),sumtot2(2),errort(2)
      common/sumvalues/fc1,fc2,add1,add2,sumtot1,sumtot2,errort

! Variables for Bragg peak calculation if needed

      integer :: nbragg,multa,multb,multc,usea,useb,usec,ncharoutput
      real :: qbraggwidth,qbraggmax,hklqmin,debyewaller,multfactor
      real, allocatable, dimension(:) :: qbragg
      real, allocatable, dimension(:,:) :: intenbragg,intenconv

! Variables for pressure control

      logical :: control_pressure,pacontrol,pbcontrol,pccontrol,ppcontrol
      real :: control_p_value,control_p_step,ratplimit
     1,dmaxfac,afac,bfac,cfac,pfac,newphib,newthetac,newphic
      character(len=5) control_p_text
      integer :: npcontrols

! Variables for fitting lattice parameters

      logical :: refinelattice
      integer :: nlatrefine,ilatrefine
      real :: alatstep,blatstep,clatstep,alphastep,betastep,gamstep
      real :: alatold,blatold,clatold,alphaold,betaold,gamold,oldfit
      real :: alatnew,blatnew,clatnew,alphanew,betanew,gamnew,newfit
      real :: alatfac,blatfac,clatfac

      character(8)  :: currentdate
      character(10) :: currenttime
      character(5)  :: currentzone

! Variables for calculating the Porod scattering

      real :: porod,porodq,porodpwr,porodqdif,porodqmin
      real, dimension(mfile) :: porodamp,porodadd
      integer :: nporodq
      logical :: foundporodq
c
c EPSR begins....
c

      pi=4.0*atan(1.0)
      degconv=pi/180.0

      elimitfac=0.02

c
c set up the EPSR symbol tables
c
      call setup_symbols_epsr()

c
c set up the EPSR input file
c
      nin=11
      nout=6
      call setup_input_file(nin,nout,ierrout)
c
c if an error occurred on setup, exit now
c
      if(ierrout.ne.0) return

c
c read the input file
c
      call read_input_file(nin,nout,ierrout)
c
c if an error was encountered on read, exit now
c
      if(ierrout.ne.0) return

c
c now read the symbol_values into their appropriate values in the
c program.
c
      icount=nsymboluptodata
      ifilecount=1
      iref=ireff(icount,ifilecount)
c NFILE = number of data files to fit (30 maximum)
      read(symbol_value(iref),*) nfile
      efilereqsum=0.0
      do in=1,nfile
c
c INFILE = name of data file for S(Q) data
c
            iref=iref+1
            INFILI(in)=symbol_value(iref)
c
c add on the directory name
c
            ilen=index(infili(in),' ')-1
            fnamedata(in)=directory_name(1:nchardir)//infili(in)(1:ilen)
c
c FNAMEW = name of file containing weighting factors for this dataset
c
            iref=iref+1
            fnamew(in)=symbol_value(iref)
c
c add on the directory name
c
            ilen=index(fnamew(in),' ')-1
            fnamewts(in)=directory_name(1:nchardir)//fnamew(in)(1:ilen)
c
c type of input file (1,2, 3 or 4, 5)
c     NRTYPE=1, Q, S(Q) ASCII values, 1 set per line
c     NRTYPE=2, Q, S(Q), E(Q) ASCII values, 1 set per line
c     NRTYPE=3,4 N, Q, S(Q), E(Q) ASCII values, 1 set per line, 4 blank lines
c(GENIE SHOW DATA format, NRTYPE=3 corresponds to histogram format, 4 = point
c mode)
c     NRTYPE=5, Q, S(Q), E(Q) ASCII values, 1 set per line, using GUDRUN .DCS format (histogram)
c RSHMIN = minimum radius for this dataset - used in removing the background
c SZEROS is the Q=0 limit for this file. If non-zero it is used to extrapolate the data
c from the smallest Q data point to Q=0. If zero then the last data point is used as the
c q=0 value.
c TWEAK is a simple factor to multiply the data by prior to any processing. Normally 1.0.
c EFILEREQ is the relative importance of this dataset in calculating the inversion matrix
c - increasing it should improve the fit to this dataset in principle, at the expense of
c the fit to the other datasets. Set it to zero if you want to calculate the simulated
c pattern for this dataset but not use the data to calculate the empirical potential. EFILEREQ can take any value but it will be normalised internally to
c prevent the inversion matrix becoming indeterminate.
c
            iref=iref+1
            read(symbol_value(iref),*) nrtype(in)
            iref=iref+1
            read(symbol_value(iref),*) RSHMIN(in)
            iref=iref+1
            read(symbol_value(iref),*) szeros(in)
            iref=iref+1
            read(symbol_value(iref),*) tweak(in)
            iref=iref+1
            read(symbol_value(iref),*) efilereq(in)
            efilereqsum=efilereqsum+efilereq(in)
            iref=iref+1
            read(symbol_value(iref),*) porodadd(in)
      end do

c NB efilereq cannot be zero for every dataset, otherwise the inversion matrix will crash if feedback > 0.0

      if(efilereqsum.le.0.0) then

         efilereq(1)=1.0

      endif
      icount=1
      ifilecount=1
      iref=ireff(icount,ifilecount)
c
c FEEDBACK = confidence in data factor. Must be 0.0 <= FEEDBACK <= 1.0. Should only be set
c to 1.0 if the matrix of weight coefficients is well conditioned.
c Can be increased or decreased to emphasize or de-emphasize the data.
c POTFAC decides whether potential will be refined and how much of the difference
c between data and fit are fed back into the EP at each iteration. Typically POTFAC = 0.1
c or 0.2. If oscillations develop in total energy of simulation, probably means POTFAC is
c too large, so reduce its amplitude. However if POTFAC is too small then simulation will
c take a long time to find a fit. POTFAC=0.0 means no EP refinement, (but simulation will
c proceed as usual).
c PRESSREQ is the requested absolute energy (kJ/mole) of empirical potential.
c This can be used to control degree of fit and energy/pressure for all datasets together.
c SIZEFACTOR (should be >=1) is used to shrink molecules during initial refinement when
c molecular entanglement may occur. Advisable to set potfac to 0 if sizefactor > 1
c
      read(symbol_value(iref),*) feedback
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      read(sometext(ncf(1):ncl(1)),*) potfacamp
      if(nwords.gt.1) then
         read(sometext(ncf(2):ncl(2)),*) potfac
      else
         potfac=potfacamp
      end if
      if(potfacamp.eq.0.0) potfac=0.0
      if(potfac.gt.1.0) potfac=1.0
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      read(sometext(ncf(1):ncl(1)),*) pressreq
      absolute_energy=pressreq
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      read(sometext(ncf(1):ncl(1)),*) ereqmin
      if(ereqmin.lt.0.0) ereqmin=0.0
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      read(sometext(ncf(1):ncl(1)),*) ereqmax
      if(ereqmax.lt.0.0) ereqmax=0.0
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      read(sometext(ncf(1):ncl(1)),*) ereqstepread
      ereqstepread=max(ereqstepread,0.0)
!      ereqstepread=min(ereqstepread,1.0)
      ereqstep=ereqstepread
      erequnit=ereqstep
      iref=iref+1
      read(symbol_value(iref),*) epcthresh
      epcthreshuse=epcthresh
      iref=iref+1
      read(symbol_value(iref),*) ereqbias
      iref=iref+1
      read(symbol_value(iref),*) navmaxparab,nparabstartfactor
      if(navmaxparab.le.0.or.nparabstartfactor.lt.1) then !If a problem with either value, set the defaults.
         navmaxparab=50
         nparabstartfactor=1
      end if
      nparabstart=navmaxparab*nparabstartfactor
      erequnit=ereqstep
      iref=iref+1
      read(symbol_value(iref),*) epcrmin
      iref=iref+1
      read(symbol_value(iref),*) epcfrac
      iref=iref+1
      read(symbol_value(iref),*) num_threads
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      read(sometext(ncf(1):ncl(1)),*) porodq
      if(nwords.lt.2.or.porodq.le.0.0) then
         porodamp=0.0
         porodpwr=0.0
      else
         read(sometext(ncf(2):ncl(2)),*) porodpwr
      end if
      iref=iref+1
      read(symbol_value(iref),*) nmolcellin
      !Determine if pressure control is required
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      control_pressure=.false.
      pacontrol=.false.
      pbcontrol=.false.
      pccontrol=.false.
      ppcontrol=.false.
      control_p_value=0.0
      control_p_step=0.0
      if(nwords.gt.0) then
         control_pressure=nwords.gt.2.and.(
     *sometext(ncf(1):ncl(1)).eq.'a'
     *.or.sometext(ncf(1):ncl(1)).eq.'b'
     *.or.sometext(ncf(1):ncl(1)).eq.'c'
     *.or.sometext(ncf(1):ncl(1)).eq.'ab'
     *.or.sometext(ncf(1):ncl(1)).eq.'bc'
     *.or.sometext(ncf(1):ncl(1)).eq.'ca'
     *.or.sometext(ncf(1):ncl(1)).eq.'abc'
     *.or.sometext(ncf(1):ncl(1)).eq.'all'
     *.or.sometext(ncf(1):ncl(1)).eq.'p')
         if(control_pressure) then
            read(sometext(ncf(2):ncl(3)),*,iostat=ierr) control_p_value,control_p_step
            control_pressure=ierr.eq.0.and.control_p_step.gt.0.0
         end if
         if(control_pressure) then
            if(sometext(ncf(1):ncl(1)).eq.'a') then
               pacontrol=.true.
               pbcontrol=.false.
               pccontrol=.false.
               ppcontrol=.false.
            else if(sometext(ncf(1):ncl(1)).eq.'b') then
               pacontrol=.false.
               pbcontrol=.true.
               pccontrol=.false.
               ppcontrol=.false.
            else if(sometext(ncf(1):ncl(1)).eq.'c') then
               pacontrol=.false.
               pbcontrol=.false.
               pccontrol=.true.
               ppcontrol=.false.
            else if(sometext(ncf(1):ncl(1)).eq.'ab') then
               pacontrol=.true.
               pbcontrol=.true.
               pccontrol=.false.
               ppcontrol=.false.
            else if(sometext(ncf(1):ncl(1)).eq.'bc') then
               pacontrol=.false.
               pbcontrol=.true.
               pccontrol=.true.
               ppcontrol=.false.
            else if(sometext(ncf(1):ncl(1)).eq.'ac') then
               pacontrol=.true.
               pbcontrol=.false.
               pccontrol=.true.
               ppcontrol=.false.
            else if(sometext(ncf(1):ncl(1)).eq.'abc'
     *.or.sometext(ncf(1):ncl(1)).eq.'all') then
               pacontrol=.true.
               pbcontrol=.true.
               pccontrol=.true.
               ppcontrol=.false.
            else if(sometext(ncf(1):ncl(1)).eq.'p') then
!This represents pressure control using the repulsive potential
               pacontrol=.false.
               pbcontrol=.false.
               pccontrol=.false.
               ppcontrol=.true.
            end if
         end if
      end if
      npcontrols=0
      if(pacontrol) npcontrols=npcontrols+1
      if(pbcontrol) npcontrols=npcontrols+1
      if(pccontrol) npcontrols=npcontrols+1
      if(ppcontrol) npcontrols=npcontrols+1
c Determine if to refine the intra-molecular structure
      iref=iref+1
      refine_intra=.false.
      coul_intra=.false.
      refine_intra_fac=0.0
      coul_intra_fac=0.0
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.gt.0) read(sometext(ncf(1):ncl(1)),*,iostat=ierr) refine_intra_fac
      if(refine_intra_fac.gt.0.0)  refine_intra=.true.
      if(nwords.gt.1) read(sometext(ncf(2):ncl(2)),*,iostat=ierr) coul_intra_fac
      if(coul_intra_fac.gt.0.0) coul_intra=.true.
      write(symbol_value(iref),'(f8.3,1x,f8.3)',iostat=ierr) refine_intra_fac, coul_intra_fac
      iref=iref+1
      irefsizefactor=iref
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      read(sometext(ncf(1):ncl(1)),*) sizefactor
      if(nwords.gt.1) then
         read(sometext(ncf(2):ncl(2)),*) decline_rate
! Make sure crazy values are ignored
         if(decline_rate.gt.1.0.or.decline_rate.lt.0.1) decline_rate=0.9
      else
         decline_rate=0.9
      end if
! Set the threshold energy below which the size factor is reduced.
      if(nwords.gt.2) then
         read(sometext(ncf(3):ncl(3)),*) thresholde
      else
         thresholde=0.0
      end if
      !Save the values of stepmi and stepri
      stepmisave=stepmi
      steprisave=stepri
      !If sizefactor > 1.0, then do not intramolecular rotations and refinement.
      ! Also do not update repulsive or empirical potentials.
      if(sizefactor.gt.1.0) then
         potfac=0.0
         stepmi=0.0
!         stepri=0.0
      else if(potfacamp.gt.0.0.and.potfac.le.0.0) then
         potfac=potfacamp
      endif
!      if(sizefactor.gt.2.0) then
!         stepri=0.0
!      end if
      write(6,'("Steps: ",4(1x,e12.5))') stepmi,stepri,stepmisave,steprisave
      if(feedback.gt.1.0) feedback=0.9
      if(feedback.lt.0.0) feedback=0.0
c do not attempt to refine potential if there is no weight given to the data
      if(feedback.le.0.0) then
         potfacamp=0.0
         potfac=0.0
      end if
c POLAR = polarizablity (not currently used)
      polar=0.0

c Set the data importance factors to the requested value. These have to be
c normalised to the number of non-zero values to avoid overemphasizing one dataset
c which could make the inversion matrix indeterminate.

      nonzero=0
      sumnonzero=0.0
      do in=1,nfile
         if(efilereq(in).gt.0.0) then
            nonzero=nonzero+1
            sumnonzero=sumnonzero+efilereq(in)
         endif
      end do

c Create the normalised dataset importance factors, EFACM

      if(nonzero.gt.0.and.sumnonzero.gt.0.0) then
            sumnonzero=sumnonzero/nonzero
      else
            sumnonzero=1.0
      endif
      do in=1,nfile
         efacm(in)=efilereq(in)/sumnonzero
      end do

c NQ = No. of Q values (Maximum = 1000)
c QSTEP = requested step in  Q
      iref=iref+1
      read(symbol_value(iref),*) nq
      if(nq.gt.mpts) nq=mpts
      iref=iref+1
      read(symbol_value(iref),*) qstep
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      read(sometext(ncf(1):ncl(1)),*) qscaleinputfactor
      if(nwords.gt.1) then
         read(sometext(ncf(2):ncl(2)),*) selfbroadpwr
      else
         selfbroadpwr=0.0
      end if
      coarsegraining=qscaleinputfactor.gt.1.0
c NEGAT not currently operational - hold at zero
      negat=0
c IRESET = 1 to set EP to zero and reinitialize certain accummulators
c     else IRESET = 0.
      iref=iref+1
      read(symbol_value(iref),*) ireset
      if(ireset.eq.1) then
         pressreq=0.0
         ereqstep=0.0
         erequnit=0.0
      else if(ireset.eq.2) then
!         pressreq=0.0
         if(ereqstep.gt.0.0) then
            pressreq=erequnit
         end if
      end if
      absolute_energy=pressreq
      iref=iref+1
      read(symbol_value(iref),*) iinit
c NTIMES = no. of trials per iteration (typically 1, but could be more)
c NITER = no. of iterations (typically 1)
c NSUMS is the number of previous iterations that have already been accumulated
c     If NSUMS is -1, then a running average is used to
c     calculate g(r) and S(Q), as per average energy and average pressure. (See comments in
c subroutine INITA.)
      iref=iref+1
      read(symbol_value(iref),*) ntimes
      iref=iref+1
      read(symbol_value(iref),*) niter
      iref=iref+1
      read(symbol_value(iref),*) nsums
c INTRA is the number of intermolecular moves that are made before an INTRA-
c molecular move (molecular vibration) is made - typically INTRA=50 or 100.
! rotfreq is the number of iterations between performing internal molecular rotations
c INTER is the number of iterations the potential is averaged over
      iref=iref+1
      read(symbol_value(iref),*) intra
      iref=iref+1
      read(symbol_value(iref),*) rotfreq
      iref=iref+1
      read(symbol_value(iref),*) inter
c INTER cannot be less than 1
      if(inter.lt.1) inter=1
c number of iterations to occur after intramolecular disordering before
c accumulating energy and pressure averages.
      nisumm=0.3*intra
c RHO no. density in scattering units per A**3 (this is determined by
c the .ATO file - this value is written here for information only)
      iref=iref+1
      read(symbol_value(iref),*) rho
c NDIST = number of r-space distributions (automatically set to npairs)
c CELLST = step size in linear region of r - controls number of r values
c
      iref=iref+1
      read(symbol_value(iref),*) CELLST
! Increase in rmax compared to default (0.5*cell)
      iref=iref+1
      read(symbol_value(iref),*,iostat=ierr) rmaxgrsave
      rmaxgr=rmaxgrsave
      iref=iref+1
      read(symbol_value(iref),*) ngrsamples
c The following factors should not be used in this version of the program.
      efacs=0.0
      pfaclj=0.0
      gfacs=0.0
      gfacmn=0.0
      gfacmx=0.0
      if(gfacs.gt.0.0.and.gfacs.lt.gfacmn) gfacs=gfacmn
      if(gfacs.gt.0.0.and.gfacs.gt.gfacmx) gfacs=gfacmx
      gfac=gfacs
c resolution function of data (Gaussian) FWHM(Q) = FWHM + FWHMQ*Q
      iref=iref+1
      read(symbol_value(iref),*) fwhm
      iref=iref+1
      read(symbol_value(iref),*) fwhmq
c WTENTP is not used
c NSMOOp = defines whether a background based on Rmin will be subtracted
c from the data prior to fitting the potential
c (0 for no background subtraction)
      wtentp=0.0
      iref=iref+1
      read(symbol_value(iref),*) nsmoop
c read broadening of Q space data using revised Lorch function
      iref=iref+3
      lorchqwidth=0.0
      lorchqwidthq=0.0
      read(symbol_value(iref),*,iostat=ierr) lorchqwidth,lorchqwidthq
!Read broadening and maximum Q for Bragg peak calculation, and Debye-Waller factor.
!Use 0 for qbraggmax if you do not want to do this calculation
      iref=iref+1
      read(symbol_value(iref),*) qbraggwidth,qbraggmax
!Read multiplicities for a, b and c axes. This so that when calculating Bragg peak intensities
!proper account is taken of the number of unit cells used to build the EPSR box
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      read(sometext(ncf(1):ncl(3)),*) multa,multb,multc !multipliticities along a, b and c
      if(nwords.gt.5) then
         read(sometext(ncf(4):ncl(6)),*) usea,useb,usec !use this axis if non-zero, otherwise no Bragg values along this axis
         if(usea.ne.0) usea=1 !Must be either 0 or 1
         if(useb.ne.0) useb=1 !Must be either 0 or 1
         if(usec.ne.0) usec=1 !Must be either 0 or 1
      else
         usea=1
         useb=1
         usec=1
      end if
c read minimum Q value to calculate Bragg intensities and radius value to change
c from local structure to long range structure
      iref=iref+1
      read(symbol_value(iref),*) hklqmin,rminbragg,debyewaller
      !Ensure the Bragg intensities are not calculated if rminbragg < 0.0
      if(rminbragg.le.0.0) qbraggmax=0.0
!If qbraggmax .gt. 0, ensure rminbragg is large enough to not affect amplitude at qbraggmax
!      if(qbraggmax.gt.0.0) rminbragg=max(rminbragg,30.0/qbraggmax)
!Allocate required arrays
      nbragg=0
      nbraggup=0
      if(qbraggmax.gt.0.0) then
         nbragg=nint(qbraggmax/qstep)+1
!For the convolution we need to include more Bragg peaks so there is no
!discontinuity at qbraggmax
         nbraggup=nint((qbraggmax+3.0/rminbragg)/qstep)+1
         if(nbraggup.gt.nq) nbraggup=nq
      end if
!Specify whether diffuse scattering will be calculated. If so we need the number of unit cells along each axis, and h,k,l pattern
! e.g. 1 2 0 would mean any h or k for which l=0, 1 1 2 wouild mean any h,k,l for which h=k, etc.
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      diffusedefined=nwords.gt.5
      if(diffusedefined) then
         read(symbol_value(iref),*) multdiffa,multdiffb,multdiffc,ldiffmax,ldiffstep,mdiffstep
         if(ldiffmax.gt.40) ldiffmax=40 !Practical limit based on current accuracy of dlms (generalised spherical harmonics)
      else
         ndiffuse=0
      end if
! Variables for lattice refinement if needed
      iref=iref+1
      read(symbol_value(iref),*) nlatrefine
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      alatstep=0
      blatstep=0
      clatstep=0
      alphastep=0
      betastep=0
      gamstep=0
      if(nwords.gt.5.and.nlatrefine.gt.0) then
         read(sometext(ncf(1):ncl(6)),*,iostat=ierr) alatstep,blatstep,clatstep,alphastep,betastep,gamstep
         if(ierr.ne.0) then
            alatstep=0
            blatstep=0
            clatstep=0
            alphastep=0
            betastep=0
            gamstep=0
            nlatrefine=0
         end if
      end if
      refinelattice=nlatrefine.gt.0
c read requested Monte Carlo rejection ratio
      iref=iref+1
      read(symbol_value(iref),*) rejrate
c read minimum Q values for Fourier transforms and potential fit, qminm and qminp. If qminp is zero, the minimum
! Q for fitting the potential is qminm.
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.gt.0) then
         read(sometext(ncf(1):ncl(1)),*,iostat=ierr) qminm
         if(nwords.gt.1) then
            read(sometext(ncf(2):ncl(2)),*,iostat=ierr) qminp
         else
            qminp=0.0
         end if
      else
         qminm=0.0
         qminp=0.0
      end if
      write(symbol_value(iref),'(2(1x,f10.5))') qminm,qminp
      facpress=0.001
      facpress3=0.01
c
c define pressure and energy factors not currently used, but needed.
c
      eljfac=1.0
      pljfac=1.0
      piFAC=1./SQRT(2.*PI)
      cell=cellmax
      datafc=4.*pi*rho
c
c write the repulsive potential coefficients
c
c     do ip=1,npairs
c           write(6,*) n1extr(ip),n2extr(ip),dcor(ip)
c     end do
c
c initialize distributions
c
      write(6,*) rmaxcell,cellmax
      call inita(ireset,iinit,nq,q,qstep,nz,z,fwhm,fwhmq,intra)
      write(6,*) rminpt,rmaxpt

!Find integer corresponding to minimum Q
      nqminm=1
      do i=1,nq
         if(q(i).le.qminm) nqminm=i
      end do
c
c Set up the weights matrix and its inverse
c
      call get_wts(q,nq,fnamewts,iinit)
c
c set the value of nesumt to inter
c
      if(iinit.ne.0) then
            inter=nesumt
      else
            nesumt=inter
      endif

c
c normalise atom and molecular coordinates to box size
c
      call normalise_coordinates()
! Get the maximum dimension of any molecule
      call get_maxmolsize()
      write(6,*) 'epsr> Maximum molecular dimension (A): ',maxmolsize*dmax
c
c Ignore energy from polarisation unless needed
c
      upol=upol*rcharge
c
c initial g(r) sums
c
      write(6,*) nsum,nesum,nsums
      if(nsums.ge.0) then
            call initg(nq,nz,fnamerun,ilenfnamerun,nsum,ireset)
      else
            call initg(nq,nz,fnamerun,ilenfnamerun,nesum,ireset)
      endif

c initialize reference and empirical potentials

      call initpt(ireset,rminpt,rmaxpt,expecf,polar,szeros,rshmin,rbroad)

! Determine position of data for Porod scattering, if applicable. Assumes q increases monotonically with iq.
      nporodq=0
      if(porodq.gt.0.0) then
         foundporodq=.false.
         porodqdif=porodq
         porodqmin=porodqdif
         iq=0
         do while (iq.lt.nq.and..not.foundporodq)
            iq=iq+1
            porodqdif=abs(porodq-q(iq))
            foundporodq=porodqdif.gt.porodqmin
            if(.not.foundporodq) then
               porodqmin=porodqdif
               nporodq=iq
            end if
         end do
      end if


C read in data to be fitted

      sumerr=0.0
      do in=1,nfile
            infile=fnamedata(in)

c get the length of this name including spaces

            ilendot=index(infile,'.')
            ilentot=ilendot+index(infile(ilendot+1:256),' ')-1
            WRITE(6,*) nrtype(in),' ',infile(1:ilentot)
            OPEN(8,FILE=infile,STATUS='OLD')

c Read a line to see if there are any problems

            read(8,*,iostat=ierr)

c for Genie style files skip the next 4 lines

            if(nrtype(in).gt.2.and.nrtype(in).lt.5) then
                  do i=1,4
                        read(8,*,iostat=ierr)
                  end do
            endif

c for GUDRUN files, ignore any lines that have # in them

            if(nrtype(in).gt.4) then
               read(8,*,iostat=ierr) hashsearch
               do while(index(hashsearch,'#').gt.0.and.ierr.eq.0)
                  read(8,*,iostat=ierr) hashsearch
               end do

c If we got here there was no hash in the line read, so treat the line as good.

            endif

c Backspace the data file 1 line to the first good line

            backspace(8)

            ilen=index(infile,'.')

c read in data to be fitted

            nqd=0
            do while(nqd.lt.mqd.and.ierr.eq.0)
               nqd1=nqd+1
               if(nrtype(in).eq.1) then
                  read(8,*,iostat=ierr) qd(nqd1),data(nqd1)
c     write(6,*) nqd,nqd1,qd(nqd1),data(nqd1)
                  err(nqd1)=0.0
               else if(nrtype(in).eq.2.or.nrtype(in).ge.5) then
                  read(8,*,iostat=ierr) qd(nqd1),data(nqd1)
     1,err(nqd1)
               else
                  read(8,2224,iostat=ierr) n,qd(nqd1),data(nqd1)
     1,err(nqd1)
2224  format(i7,3g14.7)
               endif
               if(ierr.eq.0) nqd=nqd1
            end do

c Ignore the last point which may be rubbish

            if(nrtype(in).gt.2) nqd=nqd-1

c for histogram files, calculate the correct q value

            if(nrtype(in).eq.3.or.nrtype(in).eq.5) then
                  q1=qd(1)
                  do i=1,nqd
                        q2=qd(i+1)
                        qd(i)=0.5*(q1+q2)
                        q1=q2
                  end do
            endif

c multiply data by density factor so that all data are per atom, and input Q values by qscale input factor

            do i=1,nqd
!                  qd(i)=qd(i)*qscaleinputfactor
                  data(i)=data(i)*rhofac(in)*tweak(in)
            end do

            write(6,352) nqd,infile(1:ilentot)
352   format(' Read ',i5,' points from file ',a)

c interpolate input data onto same scale as calculation

            call interp(nqd,nq,qd,q,data,datai,in,nqmind(in),nqmaxd(in)
     1,mqd,mpts,mfile)
            call interp(nqd,nq,qd,q,err,errs,in,nqmind(in),nqmaxd(in)
     1,mqd,mpts,mfile)

c prevent nqmax from being .ge. nq (This is for smoothing purposes in PTGNSQ)

            nqmin(in)=nqmind(in)
            nqmax(in)=nqmaxd(in)
            if(nqmax(in).ge.nq) nqmax(in)=nq-1
            if(nqmin(in).lt.nqminm) nqmin(in)=nqminm
            nqpts(in)=nqmax(in)-nqmin(in)+1
            write(6,958) nqmind(in),nqmaxd(in),nqmin(in),nqmax(in),nqpts(in),nq
958   format(1x,6(1x,i5))

! Set the Porod amplitude for this dataset

            if(nporodq.ge.nqmin(in).and.nporodq.le.nqmax(in)) then
               porodamp(in)=datai(nporodq,in)+porodadd(in)
               write(6,'("Porod amp. set to ",f10.5," for dataset ",i5)') porodamp(in),in
            else
               porodamp(in)=0.0
            end if
      end do
! initialise the r-factor to write
      rfacwrite=0
c
c open file to save energy, pressure and R-factor values
c
      fname=fnamerun(1:ilenfnamerun)//'.erg'
      write(6,912) fname(1:ilenfnamerun+4)
912   format((1x,a))
      call openf(9,fname)
!      open(35,file='dump.erg',status='unknown')
c
c skip to the end of the file unless iinit.eq.1
c
      if(iinit.ne.1) then

c calculate the coefficients of the least squares parabola which goes through the R-factor
c and absolute energies

         fname=fnamerun(1:ilenfnamerun)//'.qdr'
         call lsparabola1(9,fname,nparab,navmaxparab,epcfrac,aparab,bparab,cparab,rfacwrite,aav,bav,lastqf,nextqf)
         ierr=0
!Skip to end of file dump file
!         do while(ierr.eq.0)
!            read(35,*,iostat=ierr)
!         end do
!         backspace(35)

      endif
c
c begin outer loop
c
      stdevo=0.
      stde1o=0.
      qalfco=0.
c
c initialize control on size of fluctuations
c
      trat=0.0
      faco=1.0
      senmin=1.0e20

      do 3000 iter=1,niter
      wtensm=0.
      sumrf=0.
      sument=0.
      topsm=0.
      botsm=0.
      boasum=0.
      qalsum=0.
      qalsm2=0.
      nqals=0
c
c reset s(q) accumulators
c
      if(ireset.eq.1.or.iter.eq.1) then
            do j=1,2
                  do imoltype=1,nmoltype
                     nmovei(j,imoltype)=0
                     ninter(j,imoltype)=0
                     nwhole(j,imoltype)=0
                     nmovem(j,imoltype)=0
                  end do
            end do
      endif

! set the maximum separation of molecules. This corrected for the maximum size of any molecule

      sepmax=xstar(nptmax)+maxmolsize+3.0/dmax
      write(6,'(a,f10.5)') 'epsr> Maximum separation for neighbourlist = ',sepmax*dmax

! Initialise cell neighbour list

      call get_date_and_time('Start initialise_cellneighbourlist#',0)
      call initialise_cellneighbourlist(nmolcellin,xadd,yadd,zadd,nmol,sepmax)
      call get_date_and_time('Finish initialise_cellneighbourlist#',1)
! Ensure molecules are within box.
      do im=1,nmol
         call minimagec(xadd(im),yadd(im),zadd(im))
         if(tethered(im)) call minimagec(phix(im),phiy(im),phiz(im))
      end do
! Initialise various energies
      foldsv=0.0
      foldsb=0.0
      ftethered=0.0
      fcoordconstraints=0.0
      call fmolec1(ierr,num_threads,foldsv,foldsb,ftethered,fcoordconstraints)
      write(6,*) ierr,num_threads,foldsv,foldsb,ftethered,fcoordconstraints,(sume_triplet_dist(id)/real(nmol),id=1,nthreeb)
      write(35,'(a,i4,a,2(1x,e13.6),a,4(1x,e13.6))') '#Refinement: ',iter,' Initial energies: ',foldsv,foldsb,' Sizefactor et al.: ',sizefactor,stepmi,stepri,pressreq

      idone=0
c
c set up theoretical R-factor and calculate renormalization on error bars
c
      theorf=1.0

! Prevent intramolecular rotations if sizefactor > 1 and foldsv > abs(foldsb).
! This to stop molecules adopting unwarranted configurations when intermolecular energy is very large
! (e.g. after a change in sizefactor)
!      if(sizefactor.gt.1.0.and.foldsv.gt.abs(foldsb)) then
!         stepri=0.0
!      else
!         stepri=steprisave
!      end if
      write(6,'("Steps: ",4(1x,e12.5))') stepmi,stepri,stepmisave,steprisave
1999  format(2(1x,e13.6),2(1x,i16),1x,e13.6)
!      return

c exit if an error occurred in fmolec1

      if(ierr.ne.0) return
c     write(6,*) nmol,foldsv*nmol,foldsb*nmol

      call get_date_and_time('Start molfit iterations#',0)
      do 2000 itimes=1,ntimes
c     write(6,*) 'Got to molfit'
! Initialise the intramolecular energy weights
         call initialise_timer()
         call molfit(foldsv,foldsb,nisum,intra,rotfreq,ireset)
!         sumetriplets=0.0
!         do im=1,nmol
!            entest=ergold_triplet(im)
!            entest=en_mole_triplets(im,xadd,yadd,zadd,phix,phiy,phiz,xa,ya,za,mref,tethered)
!            if(im.eq.1) write(6,*) im,ergold_triplet(im),ergnew_triplet(im),entest
!            sumetriplets=sumetriplets+entest
!         end do
!         sumetriplets=sumetriplets/real(nmol)
!If nisum.eq.1, the molecular shake has occurred, so recalculate the energies
         if(nisum.eq.1) then
            call adjust_intra_wts()
            write(6,1999) foldsv,foldsb,sumtime1,sumtime2
c Ensure molecules stay within box.
            do im=1,nmol
               call minimagec(xadd(im),yadd(im),zadd(im))
               if(tethered(im)) call minimagec(phix(im),phiy(im),phiz(im))
            end do
            foldsv=0.0
            foldsb=0.0
            ftethered=0.0
            fcoordconstraints=0.0
            call fmolec1(ierr,num_threads,foldsv,foldsb,ftethered,fcoordconstraints)
         end if
         write(6,1999) foldsv,foldsb,sumtime1,sumtime2
c     write(6,*) nmol,foldsv*nmol,foldsb*nmol
         if(nisum.eq.1.and.stepmi.gt.0.0) write(35,'(a)') '#Intramolecular shake> '
         if(nisum.ne.1.and.(nisum-1).eq.rotfreq*int((nisum-1)/rotfreq).and.stepri.gt.0.0) write(35,'(a)') '#Intramolecular rotation> '
         write(35,'(i6,2(1x,e13.6),2(1x,i16))') nisum,foldsv,foldsb,sumtime1,sumtime2
2000  continue
      call get_date_and_time('Finish molfit iterations#',1)
! Ensure molecules are within box.
!      do im=1,nmol
!         call minimagec(xadd(im),yadd(im),zadd(im))
!         if(tethered(im)) call minimagec(phix(im),phiy(im),phiz(im))
!      end do
! Check final energies are consistent
! Initialise various energies
!      foldsv=0.0
!      foldsb=0.0
!      ftethered=0.0
!      call fmolec1(ierr,num_threads,foldsv,foldsb,ftethered)
!      write(35,'(a,i4,a,2(1x,e13.6),a,3(1x,e13.6))') '#Refinement: ',iter,' Final energies: ',foldsv,foldsb

! set the maximum separation of atoms

!      sepmax=xstar(npt1)

! Initialise cell neighbour list

!      call get_date_and_time('Start initialise_cellneighbourlist#',0)
!      call initialise_cellneighbourlist(nmolcell,xadd,yadd,zadd,nmol,sepmax)
!      call get_date_and_time('Finish initialise_cellneighbourlist#',1)

! *** The following section is needed if we are to refine lattice parameters

      oldfit=0.0
      alatold=alat
      blatold=blat
      clatold=clat
      call convertphibthetacphic(phib,thetac,phic,alphaold,betaold,gamold)
      alatnew=alatold
      blatnew=blatold
      clatnew=clatold
      alphanew=alphaold
      betanew=betaold
      gamnew=gamold
      ilatrefine=0
      do while (ilatrefine.lt.nlatrefine+1)
         call get_date_and_time('epsr> Start atogen#',0)
         call atogen(wav,ngrsamples)
         call get_date_and_time('epsr> Finished atogen#',1)
         write(6,*) 'Got to gensq'
         call gensq(nq,q,qwind,wav)
         !Get the Bragg intensities if required
         if(nbragg.gt.0) then
            write(6,*) qstep,qbraggmax,qbraggwidth,multa,multb,multc,nbragg,nbraggup
            if(allocated(qbragg)) deallocate(qbragg)
            if(allocated(intenbragg)) deallocate(intenbragg)
            if(allocated(intenconv)) deallocate(intenconv)
            allocate(qbragg(nbraggup),intenbragg(npairs,nbraggup),intenconv(npairs,nbragg))
            do ib=1,nbraggup
               qbragg(ib)=qstep*real(ib-1)
            end do
            twopi=2.0*pi
            !Output file name
            outputfile=fname_run(1:ilenfname_run)
            ncharoutput=ilenfname_run
            fname=directory_name(1:nchardir)//outputfile(1:ncharoutput)//'.inten'
            open(unit=25,file=fname,status='unknown')
!            if(rminbragg.le.0.0.or.rminbragg.gt.x(npt1)) rminbragg=x(npt1)
            call unnormalise_coordinates()
            call bragginten(qstep,qbraggmax,qbraggwidth,fwhmq,multa,multb,multc
     *,usea,useb,usec,hklqmin,nbraggup,qbragg,intenbragg)
            call normalise_coordinates()
            !Subtract the single molecule scattering and normalise to the number of unit cells
            write(25,*) '# ',npairs,nbragg,nbraggup
            multfactor=real(multa*multb*multc)
            do ib=1,nbraggup
               do ip=1,npairs
                  intenbragg(ip,ib)=intenbragg(ip,ib)/multfactor-qwind(1,ib,ip)
               end do
            end do
            do ib=1,nbragg
               do ip=1,npairs
                  intenconv(ip,ib)=0.0
               end do
            end do
            !Apply Debye-Waller factor to Bragg intensities
            if(debyewaller.gt.0.0) then
               do ib=1,nbraggup
                  qbraggfactor=qbragg(ib)*debyewaller
                  qbraggfactor=0.5*qbraggfactor*qbraggfactor
                  if(qbraggfactor.lt.40.0) then
                     qbraggfactor=exp(-qbraggfactor)
                  else
                     qbraggfactor=0.0
                  end if
                  do ip=1,npairs
                     intenbragg(ip,ib)=qbraggfactor*intenbragg(ip,ib)
                  end do
               end do
            end if
            !Step through Q values and add broadened intensities to convolution function
            gfactor=abs(rminbragg)/(2.0*(twopi**1.5))
            gfactor0=(abs(rminbragg)**3)/((twopi**1.5))
            do ib=1,nbraggup
               sumfactor=4.0*pi*qbragg(ib)*qstep
               do it=1,nbragg
                  sum=0.0
                  !Gaussian broadening function
                  if(qbragg(ib).gt.0.0.and.qbragg(it).gt.0.0) then
                     diff=(qbragg(ib)-qbragg(it))*rminbragg
                     add=(qbragg(ib)+qbragg(it))*rminbragg
                     if(abs(diff).lt.8.0) sum=exp(-(0.5*diff*diff))
                     if(abs(add).lt.8.0) sum=sum-exp(-(0.5*add*add))
                     sum=sum*gfactor/(qbragg(it))
                  else
                     add=qbragg(ib)*rminbragg
                     diff=qbragg(it)*rminbragg
                     sum=0.5*(diff*diff+add*add)
                     if(sum.lt.32.0) then
                        sum=gfactor0*exp(-sum)
                     else
                        sum=0.0
                     end if
                     sum=sum*qbragg(ib)
                  end if
                  sum=sum*sumfactor
                  do ip=1,npairs
                     intenconv(ip,it)=intenconv(ip,it)+(qwind(2,ib,ip)-intenbragg(ip,ib))*sum
                  end do
               end do
            end do
!Subtract the convolution function
            do ip=1,npairs
               do ib=1,nbragg
!               intenbragg(ip,ib)=intenbragg(ip,ib)-intenconv(ip,ib)
!               qwind(2,ib,ip)=qwind(2,ib,ip)+intenbragg(ip,ib)
                  qwind(2,ib,ip)=intenbragg(ip,ib)+intenconv(ip,ib)
!               qwind(2,ib,ip)=intenbragg(ip,ib)-qwind(1,ib,ip)
               end do
            end do
            do ib=1,nbragg
               write(25,*) qbragg(ib),(qwind(2,ib,ip),intenbragg(ip,ib),intenconv(ip,ib),ip=1,npairs)
            end do
            close(25)
! Output the diffuse scattering if requested
            call sharm_diffuse(qstep,nq,nsums,nfile)
         endif
         if(refinelattice.and.ilatrefine.le.nlatrefine) then
! a) Now calculate fit to data if refining lattice parameters
            rfacraw=0.0
            do in=1,nfile
               rfacdraw=0.0
! Determine whether this is an X-ray dataset
               ix=xfile(in)
               do i=1,nq
                  soqism(i,in)=0.0
                  backsc(i,in)=0.0
                  sum=0.0
                  sumraw=0.0
                  if(ix.gt.0) then
                     do ic=1,npairs
                        sumraw=sumraw+weightx(i,ic,ix)*(qwind(1,i,ic)+qwind(2,i,ic))
                     end do
! Add on the X-ray background due to mismatched form factors.
                     sum=sum+backx(i,ix)
                     sumraw=sumraw+backx(i,ix)
                  else
                     do ic=1,npairs
                        sumraw=sumraw+ciwt(in,ic)*qwind(1,i,ic)+cwt(in,ic)*qwind(2,i,ic)
                     end do
                  endif
                  if(i.ge.nqmind(in).and.i.le.nqmaxd(in)) then
                     diffraw=datai(i,in)-sumraw
                  else
                     diffraw=0.0
                  end if
                  if(i.ge.nqmin(in).and.i.le.nqmax(in)) then
                     rfacdraw=rfacdraw+diffraw*diffraw
                  end if
               end do
               if(nqpts(in).gt.0) then
                  rfacdraw=rfacdraw/nqpts(in)
               endif
               rfacraw=rfacraw+rfacdraw
            end do
            rfacraw=rfacraw/nfile
            newfit=rfacraw
            if(ilatrefine.eq.0) then
               write(6,'(a)') 'refinelattice> Starting values:-'
               write(6,'(a,6(1x,e12.6))') 'refinelattice> Initial lattice pars.: ',alatold,blatold,clatold,alphaold,betaold,gamold
               write(6,'(a,2(1x,e12.6))') 'refinelattice> Initial fit: ',newfit
            else
               write(6,'(a,i5)') 'refinelattice> End of iteration: ',ilatrefine
               write(6,'(a,6(1x,e12.6))') 'refinelattice> Old lattice pars.: ',alatold,blatold,clatold,alphaold,betaold,gamold
               write(6,'(a,6(1x,e12.6))') 'refinelattice> New lattice pars.: ',alatnew,blatnew,clatnew,alphanew,betanew,gamnew
               write(6,'(a,2(1x,e12.6))') 'refinelattice> New fit and old fit: ',newfit,oldfit
            end if
! b) If not the first iteration, compare newfit with oldfit and if better, or oldfit = 0, accept new parameters. Otherwise restore old parameters.
            if(ilatrefine.eq.0) then
               oldfit=newfit
            else if(newfit.le.oldfit) then ! Fit improved so save lattice parameters
               oldfit=newfit
               alatold=alatnew
               blatold=blatnew
               clatold=clatnew
               alphaold=alphanew
               betaold=betanew
               gamold=gamnew
               if(ilatrefine.eq.nlatrefine) ilatrefine=ilatrefine+1 !If this is the last iteration, no need to go back and recalculate S(Q).
            else !Restore original lattice parameters
               alatfac=alatold/alatnew
               blatfac=blatold/blatnew
               clatfac=clatold/clatnew
               call convertalphabetagamma(alphaold,betaold,gamold,newphib,newthetac,newphic)
               call unnormalise_coordinates()
               call change_boxsize(alatfac,blatfac,clatfac,newphib,newthetac,newphic)
               call normalise_coordinates()
            end if
! c) Choose new lattice parameters and try a new box shape/size - skip this step if the last iteration
            if(ilatrefine.lt.nlatrefine) then
               alatnew=alatold+alatstep*(2.0*ran1()-1.0)/dmax
               alatfac=alatnew/alatold
               if(cubic) then
                  blatnew=alatnew
                  clatnew=alatnew
                  blatfac=alatfac
                  clatfac=alatfac
                  alphanew=90.0
                  betanew=90.0
                  gamnew=90.0
               else
                  blatnew=blatold+blatstep*(2.0*ran1()-1.0)/dmax
                  blatfac=blatnew/blatold
                  clatnew=clatold+clatstep*(2.0*ran1()-1.0)/dmax
                  clatfac=clatnew/clatold
                  alphanew=alphaold+alphastep*(2.0*ran1()-1.0)
                  betanew=betaold+betastep*(2.0*ran1()-1.0)
                  gamnew=gamold+gamstep*(2.0*ran1()-1.0)
               end if
               call convertalphabetagamma(alphanew,betanew,gamnew,newphib,newthetac,newphic)
! unnormalise the coordinates before changing box size and shape
               call unnormalise_coordinates()
! Need to ensure density does not change
               rhoold=natomr/(alatx*blaty*clatz)
               call change_boxsize(alatfac,blatfac,clatfac,newphib,newthetac,newphic)
               rhonew=natomr/(alatx*blaty*clatz)
               vol=rhonew/rhoold
               dmaxfac=vol**(1.0/3.0)
               alatfac=dmaxfac
               blatfac=dmaxfac
               clatfac=dmaxfac
               newphib=phib
               newthetac=thetac
               newphic=phic
               call change_boxsize(alatfac,blatfac,clatfac,newphib,newthetac,newphic)
               call normalise_coordinates()
               alatnew=alatnew*alatfac
               blatnew=blatnew*blatfac
               clatnew=clatnew*blatfac
            end if
         end if
         ilatrefine=ilatrefine+1
! e) Go back to *** (atogen) and recalculate fit to data
      end do

! **** End of lattice parameter refinement.

      write(6,*) 'Got to entropc'
      call entrpc(oiseb,wav)
      nb=12
      nb2=15
c
c sum distributions
c
      nsum=nsum+1
      if(nsums.ge.0) then
            fc1=1.0/real(nsum)
            fc2=real(nsum-1)*fc1
      else
            fc1=1.0/real(nesum+1)
            fc2=real(nesum)*fc1
      endif

      ! Update uniform atom distribution if needed
!      if(nisum.ge.intra) then
!            ycoh(1)=0.0
!            atry=5.0/27.0
!            call wholecoord(alatx,dmax)
!            call wholecoord(blaty,dmax)
!            call wholecoord(clatz,dmax)
!            call uniform(ycoh,cellst,rho,cell,mpz)
!            call fractcoord(alatx,dmax)
!            call fractcoord(blaty,dmax)
!            call fractcoord(clatz,dmax)
!            nisum=0
!      end if
c
c sum Q-space distributions
c
      write(6,*) 'Got to sumsq'
      call sumsq(nq,qwind)
c
c sum r-space distributions
c
      write(6,*) 'Got to sumgr'
      call sumgr(nz,wav)
c
c calculate entropy of average distribution
c
      write(6,*) 'Got to entrpb'
      call entrpb(sentb,nsum)
c     endif
      stde1o=oiseb/rho
c
c calculate intermolecular energy and pressure
c
!      if(nthreeb.gt.0) write(6,*) ' Standard deviation (rad^2) for each triplet: ',(sume_triplet_dist(id)/real(nmol),id=1,nthreeb)
      call get_date_and_time('epsr> Start press#',0)
      call press(energd,pressd,energf,pressf,energabs,energthreeb,pressthreeb)
!      absolute_energy=energabs
      call get_date_and_time('epsr> Finished press#',1)
!      if(nthreeb.gt.0) write(6,*) ' Standard deviation (rad^2) for each triplet: ',(sume_triplet_dist(id)/real(nmol),id=1,nthreeb)
c
c correct energy for static polarisation
c
      energd=energd+upol
c
c reset accumulators if necessary
c
      if(nesum.eq.0) then
            ensum=0.0
            ensum2=0.0
            prsum=0.0
            prsum2=0.0
      endif
c
c only sum energies if nisum .ge. nisumm  - this is to prevent incorrect values which
c occur after an intramolecular disordering move from biassing the overall
c energy
c
      if(nisum.ge.nisumm) then
            nesum=nesum+1
            if(nesum.gt.nesumt) nesum=nesumt
            nesum1=nesum-1
            ensum=(nesum1*ensum+energd)/nesum
            ensumc=energd-ensum
            ensum2=(nesum1*ensum2+ensumc*ensumc)/nesum
            prsum=(nesum1*prsum+pressd)/nesum
            prsumc=pressd-prsum
            prsum2=(nesum1*prsum2+prsumc*prsumc)/nesum
      endif
      enav=ensum
      enav2=ensum2
      if(enav2.gt.0.0) then
            enav2=sqrt(enav2)
      else
            enav2=0.0
      endif
      prav=prsum
      prav2=prsum2
      if(prav2.gt.0.0) then
            prav2=sqrt(prav2)
      else
            prav2=0.0
      endif
c
c sum the contributions from the empirical potential
c
      epsum=0
      epsum1=0.0
      do in=1,nfile
            epsum=epsum+abs(energf(in))
            epsum1=epsum1+energf(in)
      end do
c
c set up intermolecular data
c
      do in=1,nfile

         ix=xfile(in)
         if(ix.gt.0) then

c If X-ray data, need to apply the Q-dependent weights


            do iq=1,nq
c
c sum the intramolecular terms and subtract them from the diffraction data
c
               if(iq.ge.nqmind(in).and.iq.le.nqmaxd(in)) then
                  sum=0.0
                  do ip=1,npairs
                     sum=sum+weightx(iq,ip,ix)*soqintra(iq,ip)
                  end do
                  datat(iq,in)=datai(iq,in)-sum
               else
                  datat(iq,in)=0.0
               endif

            end do
!Remember to subtract the intramolecular term from the zero limit as well!
            szerostemp(in)=0.0
            if(szeros(in).ne.0.0) then
               sum=0.0
               do ip=1,npairs
                  sum=sum+weightx(1,ip,ix)*soqintra(1,ip)
               end do
               szerostemp(in)=szeros(in)-sum
            endif

         else

            do iq=1,nq

               if(iq.ge.nqmind(in).and.iq.le.nqmaxd(in)) then
c
c sum the intramolecular terms and subtract them from the diffraction data
c
                  sum=0.0
                  do ip=1,npairs
                     sum=sum+ciwt(in,ip)*soqintra(iq,ip)
                  end do
                  datat(iq,in)=datai(iq,in)-sum
               else
                  datat(iq,in)=0.0
               endif

            end do
!Remember to subtract the intramolecular term from the zero limit as well!
            szerostemp(in)=0.0
            if(szeros(in).ne.0.0) then
               sum=0.0
               do ip=1,npairs
                  sum=sum+ciwt(in,ip)*soqintra(1,ip)
               end do
               szerostemp(in)=szeros(in)-sum
            endif
         endif
      end do

c Set up difference data - simulation and subtract a background from data if needed

      do id=1,nfile

c Form the simulated fit for this dataset

c If X-ray data, need to apply the Q-dependent weights

            ix=xfile(id)

c Set up the gradient for extrapolating to Q=0

            if(szeros(id).ne.0.0.and.nqmin(id).gt.1) then

               grads0=(datat(nqmin(id),id)-szerostemp(id))/(q(nqmin(id))-q(1))
               conss0=szerostemp(id)

            else

               grads0=0.0
               conss0=datat(nqmin(id),id)

            endif

            do i=1,nq
c
c Sum the intermolecular partials for this dataset at this q value
c
               sum=0.0
               if(ix.gt.0) then
                  do ic=1,npairs
                     sum=sum+weightx(i,ic,ix)*soqinter(i,ic)
                  end do
               else
                  do ic=1,npairs
                     sum=sum+cwt(id,ic)*soqinter(i,ic)
                  end do
               endif

c Calculate the difference between data and fit. For Q-values below qmin we take the first
c valid data point as the data below Qmin, unless szeros is non-zero.

               if(i.lt.nqmin(id)) then
                  if(szeros(id).ne.0.0) then
                     dts(i)=conss0+grads0*(q(i)-q(1))-sum
                  else
                     dts(i)=0.0
                  end if
               else if(i.ge.nqmin(id).and.i.le.nqmax(id)) then
                  dts(i)=datat(i,id)-sum
!Add in the Porod scattering if required
                  if(porodamp(id).gt.0.0) then
                     dts(i)=dts(i)+porodamp(id)*(porodq/q(i))**porodpwr
                  end if
               else
                  dts(i)=0.0
               endif

            end do

c
c form a smooth background to data if requested and subtract it
c
            if(nsmoop.gt.0) then
                  call ftback(id,rshmin(id),q,dts,dti)
                  do i=1,nqmax(id)
                     if(i.ge.nqmin(id)) datat(i,id)=datat(i,id)-dts(i)
                  end do
            endif
      end do




c Set up the current simulated partial structure factors

      do ip=1,npairs

         do iq=1,nq

            soqintertemp(iq,ip)=qwind(2,iq,ip)
c           soqintertemp(iq,ip)=soqinter(iq,ip)

         end do

      end do
c
c update potential
c
      wtentt=rho
c     subroutine ptgnsq(nsum,q,wav,data,qwind
c     *,rminpt,rmaxpt,szerostemp,overallfac,gfac
c     *,rminp,potnoi,wtentp
c     *,rbroad,expecf,nsmoop)
      write(6,*) 'Got to ptgnsq'
! Only calculate the new potential if the quality factor got worse at the previous iteration, or automatic adjustment of ereq is off, or
! we are initial phase of automatic adjustment of ereq.
!      if(pressreq.gt.0.0.and.(ereqstep.eq.0.0.or.nextqf.gt.lastqf.or.nparab.lt.nparabstart)) call ptgnsq(nesum,nq,q,wav,datat,soqintertemp
      if(pressreq.gt.0.0) call ptgnsq(nesum,nq,q,wav,datat,soqintertemp
     *,rminpt,rmaxpt,szerostemp,overallfac,gfac
     *,rshmin,potnoi,wtentt
     *,rbroad,expecf,nsmoop)
      !Regenerate the repulsive potential if potfac > 0
      call repulsivept(potfac)
      sumefa=0.0
      nsig=0
      sum1=0.0
      sum2=0.0
      do imoltype=1,nmoltype
         sum1=sum1+nwhole(1,imoltype)+nmovem(1,imoltype)
         sum2=sum2+nwhole(2,imoltype)+nmovem(2,imoltype)
      end do
      if(sum1.gt.0.0) then
         rfrac=sum2/sum1
      else
         rfrac=0.0
      endif

      rhoout=rho*pljfac*pljfac*pljfac
c
c for last attempt values are written to enable continuation
c of present attempt
c calculate chisq at same time
c
      eminmax=0.0
      feminmax=0.0
      if(iter.eq.niter)  then
c
c calculate the potential energies based on the mean g(r) and
c current potentials
c
            call sumpe(nz,z,fnamerun,ilenfnamerun)

c
c invert the data according to the inversion matrix
c
            do ic=1,npairs
c
c set the column title for this pair of atoms
c
               atom1=tlab1(irefpcof(ic))
               atom2=tlab2(irefpcof(ic))
               ilent1=index(atom1,' ')-1
               ilent2=index(atom2,' ')-1
               titlepsf(ic)='   '//atom1(1:ilent1)//'-'//atom2(1:ilent2)
c     write(6,*) ilent1,ilent2,tlab1(ic)(1:3),tlab2(ic)(1:3),titlepsf(ic)(1:20)

c set the assymptotic limit for g(r)

               grlim(ic)=1.0

c Minimum and maximum Q values for transform of partial structure factor data

               nqminp(ic)=1
               nqmaxp(ic)=nq

               iqinv=2
               do iq=1,nq
                  datt=0.0
                  if(nxfile.gt.0) then

c Ensure the current q values of the inversion matrix correspond to this Q value

                     do while(qinv(iqinv).lt.iq.and.iqinv.lt.nqinv)
                        iqinv=iqinv+1
                     end do


c Interpolation ratio

                     ratinter=real(iq-qinv(iqinv-1))/real(qinv(iqinv)-qinv(iqinv-1))
                     ratinter1=1.0-ratinter
                     do in=1,nfile

c Interpolate (linearly) inversion matrix to this Q value

                        cwttemp=ratinter1*cwttq(iqinv-1,ic,in)
     1+ratinter*cwttq(iqinv,ic,in)

                        if(iq.lt.nqmin(in)) then
c
c fill in the data below nqmin
c
                           term=datat(nqmin(in),in)*cwttemp
                        else
                           term=datat(iq,in)*cwttemp
                        endif
                        if(wsum(in).le.0.0) then
                           term=0.0
                        endif
                        datt=datt+term
                     end do
c
c add on the extra terms from the simulation
c
                     if(feedback.lt.1.0) then
                        do in=1,npairs
                           inref=in+nfile
                           cwttemp=ratinter1*cwttq(iqinv-1,ic,inref)
     1+ratinter*cwttq(iqinv,ic,inref)
                           datt=datt+soqinter(iq,in)*cwttemp
                        end do
                     endif
                  else
                     ratinter1=1.0
                     do in=1,nfile

c Interpolate (linearly) inversion matrix to this Q value

                        cwttemp=ratinter1*cwttq(1,ic,in)
                        if(iq.lt.nqmin(in)) then
c
c fill in the data below nqmin
c
                           term=datat(nqmin(in),in)*cwttemp
                        else
                           term=datat(iq,in)*cwttemp
                        endif
                        if(wsum(in).le.0.0) then
                           term=0.0
                        endif
                        datt=datt+term
                     end do
c
c add on the extra terms from the simulation
c
                     if(feedback.lt.1.0) then
                        do in=1,npairs
                           inref=in+nfile
                           cwttemp=ratinter1*cwttq(1,ic,inref)
                           datt=datt+soqinter(iq,in)*cwttemp
                        end do
                     endif
                  endif
                  datap(iq,ic)=datt
               end do
            end do
c
c set up titles for data output files
c
            do in=1,nfile
               ilent2=index(infili(in),'.')-1
               ilent1=ilent2-28
               if(ilent1.lt.1) ilent1=1
               titledat(in)=infili(in)(ilent1:ilent2)
            end do
c
c axis titles
c
            titler='   r'
            titleq='   Q'
c
c perform direct Fourier transform on partial structure factors and inverted data ready for output
c
            call stog(rhocg,q,datap,z,datagr,nqminp,nqmaxp,nz,npairs,mpts,mpz,mpairs,grlim)
            chisq=0.0
            do ic=1,npairs
                  rfac(ic)=0.0
                  do i=1,nq
                        a=soqinter(i,ic)
                        b=sq2inter(i,ic)
                        if(datap(i,ic).ne.0.0) then
                              dif=datap(i,ic)-soqinter(i,ic)
                              errsav(i,ic)=dif
                              rfac(ic)=rfac(ic)+dif*dif
                        else
                              errsav(i,ic)=0.0
                        endif
                  end do
                  rfac(ic)=rfac(ic)/nq
                  chisq=chisq+rfac(ic)
            end do
            infile=infili(nfile)
            ilen=index(infile,'.')
c
c write fit to partial structure factors
c
            write(ext1,331) '.f'
331   format(a2)
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
            call w_int(nq,fname,npairs,q,soqinter,sq2inter
     *,titleq,titlepsf,mpts,mpairs)

            write(ext1,331) '.s'
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
            call w_int(nq,fname,npairs,q,soqintra,sq2intra
     *,titleq,titlepsf,mpts,mpairs)
c
c write difference partial structure factors minus fit
c
            write(ext1,331) '.d'
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
c
c Output will be written between Qmin and Qmax of input datafiles
c
            lpt1=nq
            lpt2=1
            do in=1,nfile
               lpt1=min(lpt1,nqmind(in))
               lpt2=max(lpt2,nqmaxd(in))
            end do
            call w_intr(lpt1,lpt2,fname,npairs,q,errsav,sq2inter
     *,titleq,titlepsf,mpts,mpairs)
c
c write estimated partial structure factors
c
            write(ext1,331) '.q'
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
            call w_intr(lpt1,lpt2,fname,npairs,q,datap,sq2inter
     *,titleq,titlepsf,mpts,mpairs)
c
c write out the input data after interpolation by the program - the so-called
c "total" data and after background subtraction
c
c
c set up intermolecular data
c
            do in=1,nfile

c Identify any X-ray data sets

               ix=xfile(in)

               do i=1,nq
c
c sum the intramolecular terms and subtract them from the diffraction data
c
                  sum=0.0
                  if(i.ge.nqmind(in).and.i.le.nqmaxd(in)) then
                     if(ix.gt.0) then
                        do ic=1,npairs
                           sum=sum+weightx(i,ic,ix)*soqintra(i,ic)
                        end do
                     else
                        do ic=1,npairs
                           sum=sum+ciwt(in,ic)*soqintra(i,ic)
                        end do
                     endif
                  endif
                  datat(i,in)=datat(i,in)+sum
               end do
            end do
            write(ext1,331) '.t'
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
            call w_intr(lpt1,lpt2,fname,nfile,q,datat,errs
     *,titleq,titledat,mpts,mfile)
c
c calculate the fit to the input data
c
            rfactot=0.0
            rfacraw=0.0
            do in=1,nfile
               rfacd(in)=0.0
               rfacdraw=0.0

c Determine whether this is an X-ray dataset

               ix=xfile(in)

               do i=1,nq
                  soqism(i,in)=0.0
                  backsc(i,in)=0.0
                  sum=0.0
                  sumraw=0.0

                  if(ix.gt.0) then
                     do ic=1,npairs
                        sum=sum+weightx(i,ic,ix)
     1*(soqintra(i,ic)+soqinter(i,ic))
                        sumraw=sumraw+weightx(i,ic,ix)
     1*(qwind(1,i,ic)+qwind(2,i,ic))
                     end do

c Add on the X-ray background due to mismatched form factors.

                     sum=sum+backx(i,ix)
                     sumraw=sumraw+backx(i,ix)
                  else
                     do ic=1,npairs
                        sum=sum
     1+ciwt(in,ic)*soqintra(i,ic)+cwt(in,ic)*soqinter(i,ic)
                        sumraw=sumraw
     1+ciwt(in,ic)*qwind(1,i,ic)+cwt(in,ic)*qwind(2,i,ic)
                     end do
                  endif
!Add in the Porod scattering if required
                  if(q(i).gt.0.0.and.porodamp(in).gt.0.0) then
                     porod=porodamp(in)*(porodq/q(i))**porodpwr
                     sum=sum+porod
                     sumraw=sumraw+porod
                  end if

                  soqism(i,in)=sum
                  ufile(i,in)=sumraw
                  if(i.ge.nqmind(in).and.i.le.nqmaxd(in)) then
                     backsc(i,in)=datai(i,in)-soqism(i,in)
                     diffraw=datai(i,in)-sumraw
                  else
                     backsc(i,in)=0.0
                     diffraw=0.0
                  end if
                  if(i.ge.nqmin(in).and.i.le.nqmax(in)) then
                     rfacd(in)=rfacd(in)+backsc(i,in)**2
                     rfacdraw=rfacdraw+diffraw*diffraw
                  end if
               end do
               if(nqpts(in).gt.0) then
                  rfacd(in)=rfacd(in)/nqpts(in)
                  rfacdraw=rfacdraw/nqpts(in)
               else
                  rfacd(in)=0.0
                  rfacdraw=0.0
               endif
               rfactot=rfactot+rfacd(in)
               rfacraw=rfacraw+rfacdraw
            end do
            rfactot=rfactot/nfile
            rfacraw=rfacraw/nfile
c
c write out fit to total data
c
            write(ext1,331) '.u'
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
            call w_intr(1,nq,fname,nfile,q,soqism,errs
     *,titleq,titledat,mpts,mfile)
c
c write out difference between fit to total data and total data
c
            write(ext1,331) '.v'
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
            call w_intr(lpt1,lpt2,fname,nfile,q,backsc,errs
     *,titleq,titledat,mpts,mfile)
c
c write out site-site partial g(r)s as calculated by the program
c
            write(ext1,331) '.g'
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
            call w_int(nz,fname,npairs,z,gorinter,gr2inter
     *,titler,titlepsf,mpz,mpairs)
            write(ext1,331) '.l'
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
            call w_int(nz,fname,npairs,z,entsig,entsig
     *,titler,titlepsf,mpz,mpairs)
            write(ext1,331) '.y'
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
            call w_int(nz,fname,npairs,z,gorintra,gr2intra
     *,titler,titlepsf,mpz,mpairs)
c
c write out the Fourier transform of the estimated partial structure factors
c
            write(ext1,331) '.r'
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
            call w_int(nz,fname,npairs,z,datagr,gr2inter
     *,titler,titlepsf,mpz,mpairs)

c combine partial g(r)'s according to weights matrix

            do in=1,nfile
               do i=1,nz
                  sum1=0.0
                  sum2=0.0
                  do ic=1,npairs
                     sum=ciwt(in,ic)*gorintra(i,ic)
                     sum1=sum1+sum
                     sum2=sum2+sum+cwt(in,ic)*gorinter(i,ic)
                  end do
                  gz2ism(i,in)=sum1

c subtract the large r limit so result oscillates about zero.

c                 gofzism(i,in)=sum2-xwt(in)
               end do

c set the large r limit for Fourier transform of data

               grzerolim(in)=0.0
            end do

c Form the Fourier transform of the fits to the data over the same Q range as the
c data.

!            write(6,*) 'epsr> 2013: ',rhocg,q(10),nqmin,nqmax,nz,nfile,grzerolim
            call stog(rho,q,soqism,z,gofzism,nqmin,nqmax,nz,nfile,mpts,mpz,mfile,grzerolim)
c
c write out the combined partial g(r)'s
c
            write(ext1,331) '.x'
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
            call w_int(nz,fname,nfile,z,gofzism,gz2ism
     *,titler,titledat,mpz,mfile)
c
c perform direct Fourier transform on input data ready for output
c
            call stog(rho,q,datai,z,dataigr,nqmin,nqmax,nz,nfile
     *,mpts,mpz,mfile,grzerolim)
c
c write out direct Fourier transform of input data.
c
            write(ext1,331) '.w'
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
            call w_int(nz,fname,nfile,z,dataigr,gz2ism
     *,titler,titledat,mpz,mfile)
c
c write out the reference potentials
c
            write(ext1,331) '.o'
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
            call w_int(nz,fname,npairs,z,ecoul,epot
     *,titler,titlepsf,mpz,mpairs)
c
c write out the empirical potentials
c
            write(ext1,331) '.p'
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
            call w_int(nz,fname,npairs,z,xavs,xav
     *,titler,titlepsf,mpz,mpairs)
! Calculate a correlation coefficient between empirical potential and difference between simulated
! and estimated total g(r)s
            epccoeff=0.0
            epcup=0.0
            epcdn=0.0
            eminmax=0.0
!            if(energabs.gt.0.0) then
!Transform raw simulated data in Q space to r space, using same Q ranges as data
               call stog(rho,q,ufile,z,datagr2,nqmin,nqmax,nz,nfile,mpts,mpz,mfile,grzerolim)
               nugdif=0
               do in=1,nfile
                  ix=xfile(in)
                  do i=1,nz
                     ugdif(i,in)=0.0
! Get the weighted sum of gofrs at this z value. Only calculate the r-space coefficient if the sum is .ne. 0.0
                     gortemp=0.0
                     do ip=1,npairs
                        if(ix.gt.0) then
                           gortemp=gortemp+weightx(1,ip,ix)*wav(2,i,ip)
                        else
                           gortemp=gortemp+cwt(in,ip)*wav(2,i,ip)
                        end if
                     end do
                     gortemp=gortemp/ycoh(i)
                     if(gortemp.ne.0.0.and.z(i).gt.epcrmin) then
                        prod=datagr2(i,in)-dataigr(i,in)
                        ugdif(i,in)=prod
                     end if
                  end do
!                  epcfac=1.5
!                  delnz=nint(epcfac*pi/(q(nqmax(in))*cellst))+1 !Form a smoothed version of the data averaging over half the FT truncation width.
!                  nterms=2*delnz+1
!                  ifind=0
                  do i=1,nz
!                     ugdifsmo(i,in)=0.0
!                     if(ugdif(i,in).ne.0.0.and.ifind.eq.0) ifind=i !Ignore any low r region where ugdif is zero
!                     if(ifind.ne.0) then
!                        sum1=0.0 !Zero the smoothed ugdif
!                        do j=-delnz,delnz !step over the range of z values to be summed, reflecting at each end to keep total area constant
!                           iref=i+j
!                           if(iref.lt.ifind) then
!                              iref=2*ifind-iref
!                           else if(iref.gt.nz) then
!                              iref=2*nz-iref
!                           end if
!                           sum1=sum1+ugdif(iref,in)
!                        end do
!                        ugdifsmo(i,in)=sum1/real(nterms)!Form top-hat smoothed version of data
!                     end if
                     if(z(i).gt.epcrmin.and.ugdif(i,in).ne.0.0) then
                        nugdif=nugdif+1
                        prod=ugdif(i,in)
                        epcup=epcup+prod*prod
!                        prod=ugdifsmo(i,in)
!                        epcdn=epcdn+prod*prod
                     end if
                  end do
               end do
               epcup=epcup/real(nugdif)
!               fractemp=rfacraw/(rfacraw+epcup)
!               epcdn=epcdn/real(nugdif) !Ostensibly this represents unfittable systematic effects - not currently used
!               epccoeff=epcfrac*fractemp*epcup+(1.0-epcfrac)*(1.0-fractemp)*rfacraw
               epccoeff=epcfrac*epcup+(1.0-epcfrac)*rfacraw
               lastqf=nextqf
               nextqf=epccoeff
!               if(epccoeff.gt.0.0) then
!                  bias=(epcup+epcdn)/epccoeff
!               else
!                  bias=-2.0
!               end if
               if(ereqstep.gt.0.0.and.potfac.gt.0.0) then
                  factor=0.0
                  grad=0.0
                  if(nparab.gt.nparabstart.and.nparabstart.gt.0) then
! Use the parabola to define the best value of ereq (pressreq) if the parabola has a definite minimum,
! the minimum is below the current value of energabs, and energabs.gt.1.0.
                     if(cparab.ne.0.0) then
                        fcurrent=aparab+energabs*(bparab+energabs*cparab) !Fitted correlation at current absolute energy
                        epcthreshuse=epcthresh
!                        eminmax=-bparab/(2.0*cparab) !Position of minimum or maximum ugdif correlation.
!                        feminmax=aparab+eminmax*(bparab+eminmax*cparab) !Fitted correlation at min/max energy
                        grad=bparab+2.0*cparab*energabs !Gradient of fitted correlation at current energy
! The following code will tend to reduce ereq if the parabola is convex: this implies we are in a flatish region so should therefore be reducing ereq.
!                     if(cparab.lt.0.0.and.abs(1.0-eminmax/energabs).lt.0.2) then !Only do this if the maximum energy is within 20% of current energy.
!                        epcthreshuse=max(feminmax,epcthresh)
!                     end if
                        if(aparab.lt.epcthreshuse) grad=epcthreshuse-aparab !Force backwards step in this case.
                        eminmaxdif=sign(ereqbias,grad) !Strength of the bias is determined by the size of the difference between current energy and min/max energy
                        if(abs(eminmaxdif).ne.0.0) then
                           term=1.0+eminmaxdif*(eminmaxdif-4.0*ran1()+2.0) !Choose from a distribution of step sizes, biassed by size of eminmaxdif.
                           factor=2.0*ereqstep*(1.0-sqrt(term))/eminmaxdif
                        else
                           factor=2.0*ereqstep*(2.0*ran1()-1.0) !Choose a step at random in the range +-ereqstep
                        end if
                     else if(bparab.ne.0.0) then
                        fcurrent=aparab+energabs*bparab !Fitted correlation at current absolute energy
                        ctrllevel=aparab
                        ctrlbaseline=fcurrent
!                        ctrlbaseline=fcurrent+epcthresh*bav*energabs
!                        gradav1=epcthreshuse*gradav
!                        grad=bparab+gradav1 !Gradient of fitted correlation at current energy
                        grad=ctrlbaseline-ctrllevel
                        eminmaxdif=sign(ereqbias,grad) !Strength of the bias is determined by the size of the difference between current energy and min/max energy
!                        if(fcurrent.lt.epcthreshuse) then
!                           grad=epcthreshuse-fcurrent !Force backwards step in this case.
!                           eminmaxdif=sign(ereqbias,grad) !Strength of the bias is determined by the size of the difference between current energy and min/max energy
!                        else
!                        if(gradav1.gt.0.0) then
!                           eminmaxdif=ereqbias*grad/gradav1
!                        end if
                        if(abs(eminmaxdif).ne.0.0) then
!                           term=1.0+eminmaxdif*(eminmaxdif-4.0*ran1()+2.0) !Choose from a distribution of step sizes, biassed by size of eminmaxdif.
!                           factor=ereqstep*(1.0-sqrt(term))/eminmaxdif
                           factor=-ereqstep*sign(1.0,grad)
                        else
                           factor=ereqstep*(2.0*ran1()-1.0) !Choose a step at random in the range +-ereqstep
                        end if
      write(6,*) fcurrent,ctrllevel,ctrlbaseline,grad,eminmaxdir,factor
                     else
                        factor=ereqstep
                     end if
!                     if(nextqf.gt.lastqf.and.grad.lt.0.0) then !Only change the EP amplitude if the fit got better. Always change the amplitude if the gradient is +ve.
                     if(nextqf.gt.lastqf) then !Only change the EP amplitude if the fit got better.
                        factor=0.0
                     end if
                  else
                     factor=ereqstep !Assign initial ereq factor - this is such as to prevent epccoeff getting smaller than epcthresh
                  end if
!                  factor=max(factor,-0.05)
!                  factor=min(factor,+0.05)
                  if(factor.gt.0.0) factor=max(erequnit,factor) !Ensure the positive step is always at least 1 unit.
                  if(energabs.gt.2.0*erequnit.and.factor.lt.0.0) factor=min(-erequnit,factor) !Ensure the negative step is always at least -1 unit.
!                  if(pressreq.lt.100.0) then
                     pressreq=erequnit*real(nint((energabs+factor)/erequnit)) !Keep ereq integer multiples of erequnit
!                  else
!                     pressreq=erequnit*real(nint((energabs*(1.0+0.01*factor))/erequnit))
!                  end if
!                  if(pressreq.le.erequnit.and.nparab.gt.navmaxparab.and.nsums.lt.0) iinit=-1! Signal that you want to reset the distributions
!                  if(nparab.gt.nparabstart.and.nparabstart.gt.0) then
                     if(pressreq.le.max(ereqmin,erequnit).and.factor.lt.0.0) then
                        pressreq=ereqmin+10*erequnit  !To stop getting stuck at minimum ereq.
                     end if
                     if(ereqmax.gt.0.0) then
                        if(pressreq.ge.ereqmax) pressreq=ereqmin+0.75*(ereqmax-ereqmin) !To stop getting stuck at maximum ereq.
                        pressreq=min(pressreq,ereqmax)
                     end if
!                  end if
                  write(6,*) 'epsr-automatic> ', nparab,nparabstart,energabs,pressreq
!               end if
!               write(ext1,331) '.c'
!               fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
!               call w_int(nz,fname,nfile,z,ugdif,ugdifsmo,titler,titledat,mpz,mfile)
            end if
c Calculate and write out the intermolecule running coordination numbers. These are
c multiplied by the atomic fraction of the second atom of each pair. The first column
c corresponds to b about a, while the second column corresponds to a about b.

c We need to change the column titles appropriately

            do ip=1,npairs
c
c set the column title for this pair of atoms
c
                  atom1=tlab1(irefpcof(ip))
                  atom2=tlab2(irefpcof(ip))
                  ilent1=index(atom1,' ')-1
                  ilent2=index(atom2,' ')-1

c Calculate the number of spaces between titles

                  icharttitle=ilent1+ilent2+1
                  iremain=14-icharttitle-3
                  ibetween=iremain+4
                  write(formattitle,992) ibetween,iremain
992   format('(3x,a,1h-,a,',i2.2,'x,a,1h-,a,',i2.2,'x)')
                  ilenformat=12+2+10+2+2
                  write(titlepsf(ip),formattitle) atom1(1:ilent1)
     *,atom2(1:ilent2),atom2(1:ilent2),atom1(1:ilent1)

                  datagr(1,ip)=0
                  datagr2(1,ip)=0
                  jt1=n2extr(ip)
                  jt2=n1extr(ip)
                  do j=2,nz
                     cprod1=cfrac(jt1)*ycoh(j)
                     cprod2=cfrac(jt2)*ycoh(j)
                     datagr(j,ip)=datagr(j-1,ip)+gorinter(j,ip)*cprod1
                     datagr2(j,ip)=datagr2(j-1,ip)+gorinter(j,ip)*cprod2
                  end do
            end do
            write(ext1,331) '.z'
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
            call w_int(nz,fname,npairs,z,datagr,datagr2
     *,titler,titlepsf,mpz,mpairs)

C WRITE CURRENT DATA TO INPUT FILE IN A FORM SUITABLE TO BE READ AS NEXT
C INPUT FILE

            call norma(ireset,niter,niter)
c     write(6,*) 'Done Norma'
c
c reset the values of ireset and iinit
c
            if(ireset.ne.0) ireset=0
            if(iinit.gt.0) iinit=0
            if(iinit.lt.0) iinit=1
            rfacwrite=rfactot
      else
            rfacraw=rfacwrite
      endif
c
c output summary of results
c
c
c write out current energy and pressure and Rfactor and density
c
      write(9,594) energd,pressd,rfacraw,energabs,rho,epcup,epccoeff,aparab,bparab,cparab
c     write(6,594) energd,pressd,rfactot,rhoout
594   format(10(1x,e13.6))
      extfile='.out'
      fname=fnamerun(1:ilenfnamerun)//extfile(1:4)
      call openf(10,fname)
      wtent=0.0
      call date_and_time(currentdate,currenttime,currentzone)
      write(10,146) currenttime(1:2),currenttime(3:4),currenttime(5:6)
     1,currentdate(7:8),currentdate(5:6),currentdate(1:4),version
146   format(/2x,'Current date: ',a2,':',a2,':',a2
     1,' on ',a2,'-',a2,'-',a4/2x,'EPSR version ',a,/)
      write(10,147) ntimes,iter,nsum,rfacwrite,epcup,energabs
     1,epccoeff,ctrllevel,ctrlbaseline,epcthresh
     1,foldsv,foldsb,ftethered,fcoordconstraints,energthreeb,pressthreeb
     1,nparab,bav
147   format(1x,' Trial no. ',i5,5x,' Iteration no. ',i5
     1,5x,' No. of configurations in sum ',i5
     2,/1x,' Chi-square ',e13.6,' R-space coeff.: ',e13.6,' Abs. empirical energy: ',1x,e13.6
     3,/1x,' Quality factor ',e13.6,' Control level, baseline and threshold: ',3(e13.6,1x)
     4,/1x,' Pair energies: repulsive ',e13.6,' harmonic ',e13.6,' per molecule'
     5,/1x,' Tethered energy: 'e13.6,' per molecule'
     5,/1x,' Coordination constraint energy: 'e13.6,' per molecule'
     5,/1x,' Three-body energy: 'e13.6,' per molecule'
     5,/1x,' Three-body pressure: 'e13.6,' kbar'
     4,/1x,' No. of .erg lines: ',i7,' Average gradient of quality factor: ',5(e13.6,1x)
     6)
      if(nthreeb.gt.0) write(10,*) ' Energy for each triplet: ',(sume_triplet_dist(id)/real(nmol),id=1,nthreeb)
      write(10,145) sizefactor,rho
145   format(1x,' Sizefactor ',f10.5,' Number density ',e12.5)
      write(10,'("  Net charge = ",f10.5)') sumcharge
      write(10,148) (j,rfacd(j),j=1,nfile)
!      write(10,149) (j,rfac(j),j=1,npairs)
148   format((1x,' R-factor for data file ',i2,':',1x,e11.4))
149   format((1x,' R-factor for partial ',i2,':',1x,e11.4))
      write(10,219) sentb/rho,energd,upol,pressd,amubar,amubar2
     1,enav,enav2,prav,prav2
219   format(1x,' Noise in average g(r) ',e13.6
     1/1x,' Intermolecular energy and U_pol ',e13.6,1x,e13.6,' kJ/mole'
     1/1x,' Intermolecular pressure ',e13.6,' kbar'
     1/1x,' Average dipole moment and std. ',e13.6,1x,e13.6,' D'
c     1/1x,' Exponential energy ',e13.6,' kJ/mole'
c     1/1x,' Exponential pressure ',e13.6,' kbar'
     1/1x,' Av. intermolecular energy ',e13.6,'+-',e13.6,' kJ/mole'
     1/1x,' Av. intermolecular pressure ',e13.6,'+-',e13.6,' kbar')
      do j=1,2
         nmoveit(j)=0
         nintert(j)=0
         nwholet(j)=0
         nmovemt(j)=0
         do imoltype=1,nmoltype
            nmoveit(j)=nmoveit(j)+nmovei(j,imoltype)
            nintert(j)=nintert(j)+ninter(j,imoltype)
            nwholet(j)=nwholet(j)+nwhole(j,imoltype)
            nmovemt(j)=nmovemt(j)+nmovem(j,imoltype)
         end do
      end do
      write(10,220) (nmoveit(j),nintert(j),nwholet(j),nmovemt(j)
     *,j=1,2)
220   format(1x,' Number of rejects in MOLFIT:-'/((1x,4(i7,1x))))
278   format(i5/(6(1x,e12.5)))
      close (10)
c
c write input file
c
c
c first save the internal values to their symbol counterparts
c

      text='potfac'
      icount=0
      call search_for_symbol(text,icount)
      iref=ireff(icount,1)
      write(symbol_value(iref),'(2(1x,f10.6))') potfacamp,potfac
      iref=iref+1
      if(pressreq.gt.100000.0) pressreq=100000.0
      write(symbol_value(iref),'(e12.5)') pressreq
      iref=iref+1
      write(symbol_value(iref),'(f10.2)') ereqmin
      iref=iref+1
      write(symbol_value(iref),'(f10.2)') ereqmax
      iref=iref+1
      write(symbol_value(iref),'(f6.3,1x,f6.3)') ereqstepread,erequnit
      iref=iref+3
      write(symbol_value(iref),'(i5,1x,i5)') navmaxparab,nparabstartfactor
      if(nsums.ge.0) then
            nsums=nsum
      endif
      text='nsumt'
      icount=0
      call search_for_symbol(text,icount)
      iref=ireff(icount,1)
      write(symbol_value(iref),'(i5)') nsums
      text='inter'
      icount=0
      call search_for_symbol(text,icount)
      iref=ireff(icount,1)
      write(symbol_value(iref),'(i5)') inter
      iref=ireff(iresetref,1)
      write(symbol_value(iref),'(i2)') ireset
      iref=ireff(iinitref,1)
      write(symbol_value(iref),'(i2)') iinit
c
c unnormalise the coordinates before writing the .ato file
c
      call unnormalise_coordinates()
c     write(6,*) 'Done unnormalise'

c Generate new sizefactor if needed

      sizefactorold=sizefactor
      if(sizefactor.gt.1.0.and.foldsv.lt.thresholde) then
         sizefactor=sizefactor*decline_rate
         if(sizefactor.lt.1.0) sizefactor=1.0
      end if
      !If a larger positive energy excursion occurs, increase the effective box size.
!      if(sizefactor.gt.1.0.and.foldsv.gt.abs(foldsb)) then
!         sizefactor=sizefactor/decline_rate
!      end if
      write(symbol_value(irefsizefactor),'(f10.5,1x,f10.5,1x,e12.5)') sizefactor,decline_rate,thresholde
!      text='qstep'
!      call search_for_symbol(text,icount)
!      iref=ireff(icount,1)
!      write(symbol_value(iref),'(f10.5)') qstep
      text='rmaxgr'
      call search_for_symbol(text,icount)
      iref=ireff(icount,1)
      write(symbol_value(iref),'(e13.6)') rmaxgrsave
      write(6,'("Steps: ",4(1x,e12.5))') stepmi,stepri,stepmisave,steprisave
      !Restore intramolecular step sizes
      stepmi=stepmisave
      stepri=steprisave

!Adjust boxsize to control pressure if requested

      if(control_pressure.and.npcontrols.gt.0) then
!Set up the smallest value that can be added to 1.0 to give a new value
!         smallest=tiny(1.0)
!         oneplusrat=1.0+smallest
!         do while (oneplusrat.eq.1.0)
!            smallest=smallest*10.0
!            oneplusrat=1.0+smallest
!         end do
         rat=pressd-control_p_value
         if(abs(rat).gt.1.0) rat=sign(1.0,rat)
         rat=control_p_step*rat
!         if(abs(rat).lt.smallest) rat=sign(smallest,rat)
!         if(abs(rat).gt.ratplimit) rat=sign(ratplimit,rat)
         dmaxfac=(1.0+rat)**(1.0/real(npcontrols))
         afac=1.0
         bfac=1.0
         cfac=1.0
         pfac=1.0
         if(pacontrol) afac=dmaxfac
         if(pbcontrol) bfac=dmaxfac
         if(pccontrol) cfac=dmaxfac
         if(ppcontrol) pfac=1.0/dmaxfac
         write(6,*) afac,bfac,cfac,pfac
         if(ppcontrol) then
            !Only do this if potfac is greater than zero, i.e. potential is being modified
            if(potfac.gt.0.0) call multiply_pcoefex(pfac)
         else
            newphib=phib
            newthetac=thetac
            newphic=phic
            write(6,*) alat,blat,clat,phib,thetac,phic
            call change_boxsize(afac,bfac,cfac,newphib,newthetac,newphic)
            write(6,*) alat,blat,clat,phib,thetac,phic
            call uniform(ycoh,cellst,rho,rmaxgr,mpz)
         end if
      end if
!Make sure any changes to the repulsive and Gaussian parameters are recorded
!(Need to implement this - currently only new values of pcoefex will be saved.)

      call save_files_others(14)

c     write(6,*) 'Done save_files_others'

c
c Normalise coordinates ready for the next iteration
c
! Apply current sizefactor in case it has changed
      dmaxmol=dmax*sizefactor
      call normalise_coordinates()
! Recalculate the potential in case dmax has changed
!      call initpt(ireset,rminpt,rmaxpt,expecf,polar,szeros,rshmin,rbroad)

c     write(6,*) 'Done normalise'

3000  continue
      close(9)
!      close(35)
c     write(6,*) 'Done close(9)'
      return

      END

      SUBROUTINE INTERP(NQD,NZPTS,QD,QI,DATA,ALAMBI,in,NQMIN,NQMAX
     1,NQDIM,NQOUTD,mfile)
      DIMENSION QD(NQDIM),QI(NQOUTD),DATA(NQDIM),ALAMBI(NQOUTD,mfile)
C
CINTERPOLATES DATA ON GRID QD ONTO GRID QI. OUTPUT IN ALAMBI
C
      Q1=QD(1)
      Q2=QD(2)
      D1=DATA(1)
      D2=DATA(2)
      N=2
      I=1
      NQMIN=1
5     CONTINUE
      QT=QI(I)
      IF(Q1.GT.QT) GO TO 15
10    CONTINUE
      IF(Q2.LT.QT) GO TO 20
      GRAD=(D2-D1)/(Q2-Q1)
      ALAMBI(I,in)=D1+(QT-Q1)*GRAD
      I=I+1
      IF(I.GT.NZPTS) GO TO 30
      GO TO 5
20    CONTINUE
      N=N+1
      IF(N.GT.NQD) GO TO 30
      D1=D2
      Q1=Q2
      Q2=QD(N)
      D2=DATA(N)
      GO TO 10
15    CONTINUE
      I=I+1
      NQMIN=I
      GO TO 5
30    CONTINUE
      NQMAX=I-1
      RETURN
      END
c
c initialize distribution A
c
      subroutine inita(ireset,iinit,nq,q,qstep,nz,z,fwhm,fwhmq,intra)

      use molecom_mod
      use reallocation_routines

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'epsr16.inc'
      include 'epsrinput_par.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      include 'weights.inc'
      include 'iresetnumbers.inc'
      include 'lattice.inc'

      real q(mpts),z(mpz),bwidth,lorchlim,radlocalmin,radlocalmax,radusemin,radusemax,qrmin,qrmax
     *,qritmin,qritmax,qrjtmin,qrjtmqx,formfactorit,formfactorjt,scalefac,selffac,totalatoms
     *,gaussianqr2,lorentzianqr2
      real*8 j1qr,j1qr_mod,j1qr1,j1qr2
      real nitterms,niitermstotal,niitermsself
      integer ncoeffpread,nterms
      character*256 fname
      character*10 extfile
      logical found

      pi=4.0*atan(1.0)
      entwt(1)=0.0
c     write(6,966) nfile,ndist
966   format(1x,5(1x,i5))
c
c Now read the empirical potential parameters from the relevant symbol_values
c
c ROVERLAP = minimum allowed separation between atoms on different molecules
c RMINFAC = factor on sigmas used to set up minimum separation for indivdidual atom pairs
c RMINPT = point at which the truncation of the potential beings, e.g. 8 A
c RMAXPT = limiting radius to which potential is calculated. Typically this might be 10A
c RBROAD = Decay applied to potential perturbation (Normally small or zero)
c EXPECF = hardness of Exponential potential for controlling g(r) below rmin (normally 0.1)
c expecfpwr = power to increase expecf as a function of minimum distance: expecf=expecf*d**(expecfpwr)
c RCHARGE = 1 to calculate energy due to polarisation of molecule. Currently this only works
c correctly for SPC/E water. Normally set to 0.
      iref=ireff(iqminref,ifilecount)
      iref=iref+1
      read(symbol_value(iref),*) ROVERLAP
      iref=iref+1
      read(symbol_value(iref),*) RMINFAC
      iref=iref+1
      read(symbol_value(iref),*) RMINPT
      rminpt=max(1.0,rminpt) !Make sure silly values can't be specified
      iref=iref+1
      read(symbol_value(iref),*) RMAXPT
      rmaxpt=max(2.0,rmaxpt) !Make sure silly values can't be specified
      iref=iref+1
      read(symbol_value(iref),*) RBROAD
      iref=iref+1
      expecf=0.3
      expecfpwr=0.0
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
      if(nwords.gt.0) read(sometext(ncf(1):ncl(1)),*) expecf
      if(nwords.gt.1) read(sometext(ncf(2):ncl(2)),*) expecfpwr
! Decide whether or not to use a harmonic repulsive force
      iref=iref+1
      harmonic=.false.
      if(symbol_value(iref).eq.'harmonic') then
         harmonic=.true.
      else
         symbol_value(iref)='exponential'
      end if
      iref=iref+1
      read(symbol_value(iref),*) addpotential ! This decides whether to use Gaussian or 'modmorse' for additional potentials
      iref=iref+1
      read(symbol_value(iref),*) RCHARGE
      iref=iref+1
c POWER = power for repulsive term: (Normally 12, must be > 6, otherwise defaults to 12)
      read(symbol_value(iref),*) POWER
      if(power.le.6) power=12
c
c Read number of coefficients, no. of iterations to refine potential, separation in r
c between coefficients, maximum radius to calculate coefficients,
c acceptance factor, and width of functions
c
      iref=iref+1
      pdstep=0.1
      pdmin=1.0+pdstep/rminpt
! pdmax is set by programme to the truncation in the potential
      pdmax=rmaxpt
! ncoeffp is set to 10 times rmaxpt or mcoeff, whichever is smaller
!      ncoeffp=min(int(10.0*rmaxpt+0.0001),mcoeff)
      ncoeffp=min((int(rminpt/pdstep)+int(log(rmaxpt/rminpt)/log(pdmin))),mcoeff)
      write(6,*) "inita> ncoeffp = ",ncoeffp
      read(symbol_value(iref),*) ncoeffpread
      iref=iref+1
      read(symbol_value(iref),*) pdmaxread
      iref=iref+1
      read(symbol_value(iref),*) pdstepread
      iref=iref+1
      read(symbol_value(iref),*) npitss
      iref=iref+1
      read(symbol_value(iref),*) paccept
      iref=iref+1
      !Select whether to use Gaussian or Poisson representation of EP.
      gaussian=.false.
      if(symbol_value(iref).eq.'T'.or.symbol_value(iref).eq.'t'.or.symbol_value(iref).eq.'Gaussian') then
         gaussian=.true.
         symbol_value(iref)='Gaussian'
      else
         symbol_value(iref)='Poisson'
      end if
      iref=iref+1
      read(symbol_value(iref),*) psigma2read
      iref=iref+1
      read(symbol_value(iref),*) refpotfac
!      pdstep=pdmax/ncoeffp
      psigma2=psigma2read
      !To prevent overflows, there is a limit on psigma2
      if(psigma2.lt.0.004) then
         psigma2=0.004
         psigma2read=psigma2
      end if

      if(ireset.eq.0) then
c
c only need to open the .inpa file if iinit.eq.0
c
            if(iinit.eq.0) then
               extfile='.inpa'
               fname=fnamerun(1:ilenfnamerun)//extfile(1:5)
               open(10,file=fname,status='old',iostat=ierr)
c
c if there are problems with the inpa file, then set iinit to 1 and do not attempt to read any
c more from this file
c
               if(ierr.ne.0) then
                  iinit=1
               endif
            endif
c
c if iinit is still zero, then carry on reading .inpa file
c
            if(iinit.eq.0) then
c
c OTHER parameters relate to calculating average energy and pressure
c
               read(10,*) NISUM,NESUM,NESUMT,ENSUM,ENSUM2,PRSUM,PRSUM2
c
c read in existing site-site coefficients from .PCOF file
c
               do id=1,nfile
                  read(10,*)
                  read(10,*) npcoeff(id),pstepd(id),sigma2(id)
c
c read the difference coefficients - these represent the difference between the data
c of the present simulation and the fit.
c
                  read(10,*) (pcoefd(ic,id),ic=1,npcoeff(id))
               end do
c
c read the inverted weights matrix for this potential
c
               read(10,*) feedbackold,(efacold(in),in=1,nfile),nqinv
               nnew=nfile+npairs
               do iq=1,nqinv
               do ip=1,npairs
                  read(10,*) qinv(iq)
                  read(10,*) (cwttq(iq,ip,id),id=1,nnew)
               end do
               end do
            else
               nsum=0
               nesum=0
               nesumt=5
               nsums=-1
               nisum=intra
c
c initialise the number of coefficients for each dataset and
c set the potential coefficients to zero
c
               do id=1,nfile
                  npcoeff(id)=ncoeffp
                  sigma2(id)=psigma2
                  do ic=1,npcoeff(id)
                        pcoefd(ic,id)=0.0
                  end do
                  pstepd(id)=0.002
               end do
            endif
            close(10)
c
c initialise certain accumulators
c
            if(nsums.ge.0) then
                  nsum=nsums
            else
                  nsum=0
            endif
            npits=npitss
      else
            nsum=0
            nesum=0
            nesumt=5
            nsums=-1
            nisum=intra
c
c initialise the number of coefficients for each dataset and
c set the potential coefficients to zero
c
            do id=1,nfile
                  npcoeff(id)=ncoeffp
                  sigma2(id)=psigma2
                  do ic=1,npcoeff(id)
                        pcoefd(ic,id)=0.0
                  end do
                  pstepd(id)=0.002
            end do
            npits=npitss
      end if
c
c rmaxgr should be no bigger than half the box size
c
      if(rmaxgr.gt.0.5*cellmax.or.rmaxgr.le.0.0) rmaxgr=0.5*cellmax
      write(6,*) 'epsr-inita> rmaxgr,rmaxcell,cell,cellmax = ',rmaxgr,rmaxcell,cell,cellmax
c
c set up x values (=r values). For x>rlog the step size increases exponentially
c
      NPT=0
      x3dold=0.
      celstt=cellst
      stpfac=1.+celstt/rlog
      celmin=cell
      x1dold=-0.5*celstt
      cellhf=0.25*celstt
      rminsv=x1dold
      x(1)=0.0
      j=2
      npt1=2
      npt2=2
      nptmin=2
      nptmax=2
      x1d=x1dold+celstt
      do while(j.lt.mpz.and.x1d.lt.rmaxgr)
         x3d=x1d**3/3
c        x3d=x1d
         x(j)=x1d-cellhf
         entwt(j)=1/x(j)
         x1dold=x1d
         x3dold=x3d
c        if(x1d.gt.rlog) celstt=celstt*stpfac
         cellhf=0.5*celstt
c
c save the integer corresponding to the range of the potential
c
         if(x1d.lt.rminpt) nptmin=j+1
         if(x1d.lt.rmaxpt) nptmax=j+1
c
c save the index for x values below the specified maximum value for
c valid Fourier transforms
c
         if(x1d.lt.rmaxgr) npt1=j
c     write(6,*) j,x(j),nptmin,nptmax,npt1
         j=j+1
         x1d=x1dold+celstt
      end do
c
c set the last value of x
c
      x(j)=x1d-cellhf
c     write(6,*) j,x(j),nptmin,nptmax,npt1
      npt=j
      npt2=j-2
      if(npt1.gt.npt2) npt1=npt2
      write(6,959) nptmin,nptmax,npt1,npt2,npt
959   format(5(1x,i5))
      decay=x(nptmax)-x(nptmin)
      if(nesumt.lt.1) nesumt=1
      if(nesum.gt.nesumt-1) nesum=nesumt-1
      if(nisum.gt.intra) nisum=intra
c
c calculate uniform distribution for the box
c To do this try to ensure that at least one element in the nearest shell (27
c boxes within 1 step from origin, including one at origin) is
c sampled on average. We repeat this calculation every intra MC moves so that
c it does not unduly bias the summed g(r).
c
      volume=alatx*blaty*clatz
      rhotemp=real(natomr)/volume
      call uniform(ycoh,cellst,rhotemp,rmaxgr,mpz)
!      if(ireset.ne.0.or.nisum.ge.intra) then
c     if(ireset.ne.0) then
!            ycoh(1)=0.0
!            atry=5.0/27.0
!            call uniform(ycoh,cellst,rhotemp,rmaxgr,mpz)
!            nisum=intra
!      else
!            atry=0.1/27.0
c
c Read the uniform distribution from a file called ****.uni
c
!            extfile='.uni'
!            fname=fnamerun(1:ilenfnamerun)//extfile(1:4)
!            call openf(10,fname)
!            read(10,*,iostat=isave) (a,ycoh(j),j=1,npt1)
!            close(10)
c
c if there was an error encountered (most likely because npt2 has changed)
c we regenerate the uniform distribution
c
!            if(isave.ne.0) then
!               ycoh(1)=0.0
!               atry=5.0/27.0
!               call uniform(ycoh,cellst,rhotemp,rmaxgr,mpz)
!               nisum=intra
!            endif
!      endif

c
c define the q scale on which calculation is made
c
      do i=1,nq
         q(i)=qstep*(i-1)
      end do
c
c sigma for instrumental resolution
c
      sigres=0.5*fwhm/sqrt(2.*alog(2.0))
      sigrsq=0.5*fwhmq/sqrt(2.*alog(2.0))
c
c calculate transform matrices
c
      pi2=2*pi
      trfac=4.0*pi*rhotemp
      lorchlim=4.4934/rmaxgr
      do i=1,nq
            qq=q(i)
            r1=0
            j1qr1=j1qr(qq,r1)
            rstep=x(1)-r1
            DO j=1,npt1
                  rr=x(j)
                  rstep=rr-r1
                  r2=rr+rstep
                  j1qr2=j1qr(qq,r2)
                  if(ycoh(j).gt.0.0) then
                     foqols(j,i)=trfac*(j1qr2-j1qr1)/ycoh(j)
                  else
                     foqols(j,i)=0.0
                  end if
                  r1=r2
                  j1qr1=j1qr2
                  qsig=sigres+qq*sigrsq
                  qrs=qsig*rr
                  qr2=0.5*qrs*qrs
                  sums=0.0
                  if(qr2.lt.30.) then
                     sums=exp(-qr2)
                     if(lorchqwidth.gt.0.0.or.lorchqwidthq.gt.0.0) then
                        bwidth=rr*lorchlim*(lorchqwidth+lorchqwidthq*qq)
                        sums=sums*pofr(bwidth)
                     endif
                  endif
                  foqols(j,i)=foqols(j,i)*sums
            end do
      end do
! If coarse graining, calculate the spherical form factors for each atom type
      if(coarsegraining) then
! Calculate modified scale factors based on number of atoms in each molecule
         call reallocate1d_r(scalefacmoluse,nmol)
         do im=1,nmol
! Get the total number of atoms in this molecule
            totalatoms=0.0
            do ia=naf(im),nal(im)
               totalatoms=totalatoms+natomsphere(itype(ia))*pseudofacmol(itype(ia))
            end do
!selffac is used because intra-molecular atomic separations can be significant, especially for larger molecules
            scalefacmoluse(im)=totalatoms**selfbroadpwr !May need to use cube root? Or extra smearing factor?
         end do
         qfactor=4.0*pi/3.0 !Factor to get average number of atoms per spherical shell
!         qfactor=4.0*pi*(qscaleinputfactor**3) !Factor to get average number of atoms per spherical shell
!         qfactor=(qscaleinputfactor**3) !Factor to get average number of atoms per spherical shell
! In this first part we determine how many spherical form factor bins there will be for each atom type (nscalefacmol), and the bin no. for each molecule
! (All the atoms in a molecule will belong to this same bin, since the scale factor applies to the molecule, not the individual atoms.)
! Determine the minimum and maximum scale factors for each molecular type
         if(allocated(scalefactypemin)) deallocate(scalefactypemin,scalefactypemax,nsffbintype)
         allocate(scalefactypemin(ntype),scalefactypemax(ntype),nsffbintype(ntype))
         scalefactypemin=1.0
         scalefactypemax=1.0
         do ia=1,natom
            im=mref(ia)
            it=itype(ia)
            scalefactypemin(it)=min(scalefacmol(im),scalefactypemin(it))
            scalefactypemax(it)=max(scalefacmol(im),scalefactypemax(it))
         end do
         maxnscalefactype=0
         nsffbintype=0
         ! Get the number of form factor bins for each atom type
         do it=1,ntype
            nsffbintype(it)=int((scalefactypemax(it)-scalefactypemin(it))/scalefacratio)+1 !No. of scale factor bins for each atom type
            if(nsffbintype(it).gt.maxnscalefactype) maxnscalefactype=nsffbintype(it)
!            write(6,*) 'inita> Minimum and maximum scale factors for atom type: ',it,scalefactypemin(it),scalefactypemax(it),nsffbintype(it)
            !Now search through the atom types, and assign the
         end do
         ! Assign the sphere form factor bin for each molecule, and create a store for the no. terms in each form factor bin for each atom type
         if(allocated(assignatomscalefac)) deallocate(assignatomscalefac)
         allocate(assignatomscalefac(natom))
         do ia=1,natom
            im=mref(ia)
            it=itype(ia)
            ibin=int((scalefacmol(im)-scalefactypemin(it))/scalefacratio)+1
            assignatomscalefac(ia)=ibin
!            write(6,*) 'inita> ',ia,ibin
         end do
         if(allocated(sphereformfactor)) deallocate(sphereformfactor,intraformfactor,intraformfactorself,ntermssff,ntermsintra,ntermsintraself)
         allocate(sphereformfactor(nq,maxnscalefactype,ntype),intraformfactor(nq,npairs),intraformfactorself(nq,ntype)
     *,ntermssff(maxnscalefactype,ntype),fracntermssff(maxnscalefactype,ntype),ntermsintra(npairs),ntermsintraself(ntype))
! Zero all the type form factors: will form an average over the different radii for each type
         ntermssff=0
         fracntermssff=0
         ntermsintra=0.0
         ntermsintraself=0.0
         sphereformfactor=0.0
         intraformfactor=0.0
         intraformfactorself=0.0
         do ia=1,natom
            im=mref(ia)
            it=itype(ia)
!            qfactorit=qfactor*cfraccg(it)
            qfactorit=qfactor*rhosphere(it)
            radlocalmin=radminsphere(it)
            radlocalmax=radmaxsphere(it)
            radlocalmin3=radlocalmin**3
            radlocalmax3=radlocalmax**3
!            radlocalmin=0.0
!            radlocalmax=0.5*sigcor(it)
            ibin=assignatomscalefac(ia)
            ntermssff(ibin,it)=ntermssff(ibin,it)+1
            scalefac=scalefacmol(im)
            selffac=scalefacmoluse(im)
!            selffac=1.0
            radusemin=radlocalmin*scalefac
            radusemax=radlocalmax*scalefac
            radusemin3=qfactorit*radusemin**3
            radusemax3=qfactorit*radusemax**3
            radusemin=radlocalmin*selffac
            radusemax=radlocalmax*selffac
!            radusemin=radlocalmin
!            radusemax=radlocalmax
            do i=1,nq
               qq=q(i)
! Weight the sphere form factors by the volume of each sphere
               sphereformfactor(i,ibin,it)=sphereformfactor(i,ibin,it)+(radusemax3*gaussianqr2(qq,radusemax)-radusemin3*gaussianqr2(qq,radusemin))
!               if(isnan(sphereformfactor(i,ibin,it))) then
!                  write(6,*) 'inita-1>',i,ibin,it,nq,ntype,q(i),qrmin,qrmax,radusemin3*gaussianqr2(qrmin),radusemax3*gaussianqr2(qrmax)
!                  stop
!               end if
! Weight the sphere form factors by the average volume of this sphere
!               sphereformfactor(i,it)=sphereformfactor(i,it)+qfactor*rhosphere(it)*(radlocalmax3*gaussianqr2(qrmax)-radlocalmin3*gaussianqr2(qrmin))
!               sphereformfactor(i,it)=sphereformfactor(i,it)+qfactor*gaussianqr2(qrmax)
            end do
         end do
         sumatoms=0 !Will need to normalise everything per atom
         do it=1,ntype
            nterms=0
            sumatomsbin=0.0
            do ib=1,nsffbintype(it)
               nterms=nterms+ntermssff(ib,it)
!               write(6,*) 'inita> ',it,ib,ntermssff(ib,it),radminsphere(it),radmaxsphere(it),rhosphere(it),sphereformfactor(1,ib,it),ntermssff(ib,it)
               sumatomsbin=sumatomsbin+sphereformfactor(1,ib,it) !Q = 0 value corresponds to number of atoms in this type of spherical shell
               if(ntermssff(ib,it).gt.1) then
                  do i=1,nq
                     sphereformfactor(i,ib,it)=sphereformfactor(i,ib,it)/real(ntermssff(ib,it)) !Mean sphere form factor for this form factor bin
                  end do
               end if
            end do
!            write(6,*) 'inita> sumatomsbin: ',it,cfrac(it),nsffbintype(it),sumatomsbin,nterms
! Get fractional ntermssff for each atom type
            if(nterms.gt.0) then
               sumatoms=sumatoms+cfrac(it)*sumatomsbin/real(nterms) !Normalise to number of terms. Don't need atomic fractions here because nterms is already proportional to atom fraction
               do ib=1,nsffbintype(it)
                  fracntermssff(ib,it)=real(ntermssff(ib,it))/real(nterms)
               end do
            else
               do ib=1,nsffbintype(it)
                  fracntermssff(ib,it)=0.0
               end do
            end if
!            write(6,*) 'inita> :',it,(fracntermssff(ib,it),ib=1,nsffbintype(it))
         end do
         sumatoms=sumatoms*real(natomr)/real(nmol) ! convert to average number of atoms per molecule
         write(6,*) 'inita> Average number of atoms per molecule: ',sumatoms
         if(sumatoms.gt.0.0) then
            sfactor=sqrt(1.0/sumatoms) !Eventually everything is normalised to number of atoms in box
            sphereformfactor=sfactor*sphereformfactor
!            do it=1,ntype
!               sfactorit=sfactor/cfraccg(it)! When comparing with data we want the total form factor for the molecule, not weighted by atomic fraction products
!               do ib=1,nsffbintype(it)
!                  do i=1,nq
!                     sphereformfactor(i,ib,it)=sfactorit*sphereformfactor(i,ib,it)
!                  end do
!               end do
!            end do
         end if
         open(10,file='formfactors.txt',status='unknown')
         do i=1,nq
            write(10,*) q(i),(sphereformfactor(i,1,it),it=1,ntype)
         end do
         close(10)
! Now do the intra-molecular terms
         do im=1,nmol
            scalefac=scalefacmol(im)
!selffac is used because intra-molecular atomic separations can be significant, especially for larger molecules
!            selffac=scalefacmoluse(im) !May need to use cube root? Or extra smearing factor?
            selffac=scalefacmol(im) !May need to use cube root? Or extra smearing factor?
            do ia=naf(im),nal(im)
               it=itype(ia)
               qfactorit=qfactor*rhosphere(it)
               raditmin=radminsphere(it)*scalefac
               raditmax=radmaxsphere(it)*scalefac
               raditmin3=qfactorit*raditmin**3
               raditmax3=qfactorit*raditmax**3
!               raditminself=raditmin
!               raditmaxself=raditmax
               raditmin=raditmin*selffac
               raditmax=raditmax*selffac
               do ja=ia,nal(im)
                  jt=itype(ja)
                  ip=npairref(it,jt)
                  ntermsintra(ip)=ntermsintra(ip)+1
                  if(it.eq.jt) then
                     nitterms=natomsphere(it)*pseudofacmol(it) !Total number atoms of this type in pseudo molecule
                     if(nitterms.gt.0) then !If there are no self terms for this atom, ignore it
                        niitermstotal=nitterms*(nitterms-1) !Total number of pairs of type it,it in pseudomolecule
                        niitermsself=pseudofacmol(it)*natomsphere(it)*(natomsphere(it)-1) !Number true "self" pairs in the pseudo molecule
                        selfcorrself=niitermsself/(nitterms*nitterms) !Fractional contribution from true self terms
!                        selfcorr=0.5*(niitermstotal-niitermsself)/(nitterms*nitterms) !Fractional contribution from remainder
                        selfcorr=(niitermstotal-niitermsself)/(nitterms*nitterms) !Fractional contribution from remainder
!                     write(6,*) 'inita> ',it,nitterms,niitermstotal,niitermsself,selfcorrself,selfcorr
                        ntermsintraself(it)=ntermsintraself(it)+1
                        do i=1,nq
                           qq=q(i)
                           formfactorit=(raditmax3*gaussianqr2(qq,raditmax)-raditmin3*gaussianqr2(qq,raditmin))
                           intraformfactorself(i,it)=intraformfactorself(i,it)+selfcorrself*formfactorit*formfactorit
                           intraformfactor(i,ip)=intraformfactor(i,ip)+selfcorr*formfactorit*formfactorit
                        end do
                     end if
                  else
                     qfactorjt=qfactor*rhosphere(jt)
                     radjtmin=radminsphere(jt)*scalefac
                     radjtmax=radmaxsphere(jt)*scalefac
                     radjtmin3=qfactorjt*radjtmin**3
                     radjtmax3=qfactorjt*radjtmax**3
                     radjtmin=radjtmin*selffac
                     radjtmax=radjtmax*selffac
                     do i=1,nq
                        qq=q(i)
                        formfactorit=(raditmax3*gaussianqr2(qq,raditmax)-raditmin3*gaussianqr2(qq,raditmin))
                        formfactorjt=(radjtmax3*gaussianqr2(qq,radjtmax)-radjtmin3*gaussianqr2(qq,radjtmin))
                        intraformfactor(i,ip)=intraformfactor(i,ip)+formfactorit*formfactorjt
                     end do
                  end if
               end do
            end do
         end do
!         write(6,*) 'inita> ',(ntermsintra(ip),ip=1,npairs)
!         write(6,*) 'inita> ',(ntermsintraself(ip),ip=1,ntype)
!Normalise form factors to number of atoms
         do ip=1,npairs
            if(ntermsintra(ip).gt.0.and.sumatoms.gt.0.0) then
               sfactor=1.0/(ntermsintra(ip)*sumatoms)
!               sfactor=1.0/(ntermsintra(ip))
               do i=1,nq
                  intraformfactor(i,ip)=sfactor*intraformfactor(i,ip)
               end do
            end if
         end do
         do it=1,ntype
            if(ntermsintraself(it).gt.0.and.sumatoms.gt.0.0) then
               ip=npairref(it,it)
               sfactor=1.0/(ntermsintraself(it)*sumatoms)
!               sfactor=1.0/(ntermsintra(ip)*sumatoms)
               do i=1,nq
                  intraformfactorself(i,it)=sfactor*intraformfactorself(i,it)
               end do
            end if
         end do
         open(10,file='intraformfactors.txt',status='unknown')
         do i=1,nq
            write(10,*) q(i),(intraformfactor(i,ip),ip=1,npairs)
         end do
         open(10,file='intraformfactorsself.txt',status='unknown')
         do i=1,nq
            write(10,*) q(i),(intraformfactorself(i,it),it=1,ntype)
         end do
         close(10)
      end if
c
c set z scale identical to x scale - inherited from earlier version
c
      nz=npt1
      do i=1,nz
         z(i)=x(i)
      end do
!      stop
      return
      end

      function j1qr(q,r)

      real q,r
      real*8 j1qr,qr

      if(q.gt.0.0) then
         qr=q*r
         j1qr = (sin(qr)-qr*cos(qr))/(q*q*q)
      else
         j1qr = r*r*r/3.0
      end if

      return
      end

      function j1qr_mod(q,r)

      real q,r
      real*8 j1qr_mod,qr

      qr=q*r
      j1qr_mod=0.0
      if(qr.gt.0.0) then
         j1qr_mod = (sin(qr)-qr*cos(qr))/(qr*qr*qr)
      else
         j1qr_mod = 1.0/3.0
      end if

      return
      end

      function gaussianqr2(q,r)

!Gaussian form factor for Gaussian density distribution across sphere

      real q,r,qr,arg
      real gaussianqr2

      gaussianqr2=1.0
      if(r.le.0.0) return !If sigma is zero, uniform scattering
      qr=q*r
!      arg=0.64309*qr !ratio of sigma to radius so that r=0 density is same as uniform sphere (sqrt(2/pi)/3)**(1/3)
      arg=0.51042225*qr !Gives volume of Gaussian half that of sphere, and about 73% of mass inside specified sphere radius
      arg=0.5*arg*arg
      if(arg.lt.30.0) then
         gaussianqr2=exp(-arg) !The 3 is because the number of atoms in the sphere is 4*pi (already included) / 3 (not included) * R**3.
      else
         gaussianqr2=0.0
      end if
      return
      end

      function lorentzianqr2(q,r)

!Gaussian form factor for Gaussian density distribution across sphere

      real q,r,qr,arg
      real lorentzianqr2

      lorentzianqr2=0.0
      if(r.le.0.0) return !If sigma is zero, there can be no scattering
      qr=q*r
!      arg=0.64309*qr !ratio of sigma to radius so that r=0 density is same as uniform sphere (sqrt(2/pi)/3)**(1/3)
      arg=0.59647*qr !Gives volume of Gaussian half that of sphere, and about 73% of mass inside specified sphere radius
      if(arg.lt.30.0) lorentzianqr2=exp(-arg)/3 !The 3 is because the number of atoms in the sphere is 4*pi (already included) / 3 (not included) * R**3.

      return
      end

      function pofr(delr)

      if (delr.le.0.0) then
         pofr=1.0
      else
       delr3=3.0/(delr*delr*delr)
       pofr=delr3*(sin(delr)-delr*cos(delr))
      endif

      return
      end
c
c initialize g(r) sums
c
      subroutine initg(nq,nz,inpfile,ileninpfile,nsum,ireset)
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      character*256 fname,inpfile
      character*2 ext
      real*4 qarray(mpts),rarray(mpz)
      integer*4 ierr
c     write(6,*) nsum
      if(nsum.gt.0) then
c
c read in existing fits in r and Q space
c
            write(ext,100)
100   format('.g')
            fname=inpfile(1:ileninpfile)//ext(1:2)
            call r_int(lptin,fname,npairs,rarray,gorinter,gr2inter
     *,mpz,mpairs,ierr)

            write(ext,101)
101   format('.y')
            fname=inpfile(1:ileninpfile)//ext(1:2)
            call r_int(lptin,fname,npairs,rarray,gorintra,gr2intra
     *,mpz,mpairs,ierr)

            write(ext,102)
102   format('.f')
            fname=inpfile(1:ileninpfile)//ext(1:2)
            call r_int(lptin,fname,npairs,qarray,soqinter,sq2inter
     *,mpts,mpairs,ierr)

            write(ext,103)
103   format('.s')
            fname=inpfile(1:ileninpfile)//ext(1:2)
            call r_int(lptin,fname,npairs,qarray,soqintra,sq2intra
     *,mpts,mpairs,ierr)

      else

c
c zero gorinter
c
            nsum=0
            do id=1,npairs
                  do j=1,nz
                        gorintra(j,id)=0.
                        gr2intra(j,id)=0.
                        gorinter(j,id)=0.
                        gr2inter(j,id)=0.
                  end do
            end do
            do in=1,npairs
                  do i=1,nq
                        soqintra(i,in)=0.
                        sq2intra(i,in)=0.
                        soqinter(i,in)=0.
                        sq2inter(i,in)=0.
                  end do
            end do
      endif
      return
      end
c
c calculate entropy of distribution
c
      subroutine entrpa(n1,n2,dt,entrop,mpts)
      real*4 dt(mpts)
      entrop=0.
      yav1=dt(n1-1)
      yav2=dt(n1)
      ysum=0.0
      do j=n1,n2
            yav3=dt(j+1)
            ysum=ysum+1.0
            rj=0.5*abs(yav3-yav1)
            pj=0.25*(2.*yav2+yav1+yav3)
            dif=yav2-pj
            difa=abs(dif)
            if(difa.lt.0.5*rj) then
                  sent=dif*dif/rj
            else
                  sent=difa-0.25*rj
            endif
            entrop=entrop+sent
            yav1=yav2
            yav2=yav3
      end do
      if(ysum.gt.0.0)   entrop=entrop/ysum
      return
      end
c
c calculate entropy of average distribution
c
      subroutine entrpb(entrop,nsum)
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      if(nsum.le.0) then
            entrop=0.
      else
            n1=2
            n2=npt2
            entrop=0.
            do id=1,nfile
            yav1=0.0
            yav2=gorinter(n1,id)*ycoh(n1)
            do j=n1,n2
                  yav3=gorinter(j+1,id)*ycoh(j+1)
                  rj=0.5*abs(yav3-yav1)
                  pj=0.25*(2.*yav2+yav1+yav3)
                  dif=yav2-pj
                  difa=abs(dif)
                  if(difa.lt.0.5*rj) then
                        sent=dif*dif/rj
                  else
                        sent=difa-0.25*rj
                  endif
                  entrop=entrop+sent*entwt(j)
                  yav1=yav2
                  yav2=yav3
            end do
            end do
            ysuma=1.0
            entrop=entrop/nsum/ysuma
      endif
      return
      end
c
c calculate entropy of molecular distribution
c
      subroutine entrpc(entrop,wav)
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      dimension wav(2,mpz,mpairs)
      n1=2
      n2=npt2
      entrop=0.
      do id=1,npairs
      wav(2,npt2+1,id)=ycoh(npt2+1)
      wav(2,npt2+2,id)=ycoh(npt2+2)
      yav1=0.0
      yav2=wav(2,n1,id)
      do j=n1,n2
            yav3=wav(2,j+1,id)
            rj=0.5*abs(yav3-yav1)
            pj=0.25*(2.*yav2+yav1+yav3)
            dif=yav2-pj
            difa=abs(dif)
            if(difa.lt.0.5*rj) then
                  sent=dif*dif/rj
            else
                  sent=difa-0.25*rj
            endif
            entrop=entrop+sent*entwt(j)
            yav1=yav2
            yav2=yav3
      end do
      end do
      return
      end

c
c sum g(r) and error bars
c
      subroutine sumgr(nz,wav)
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'potential_par.inc'
      real wav(2,mpz,mpairs)

      real*8 fc1,fc2,add1(2),add2(2),sumtot1(2),sumtot2(2),errort(2),zerolimit(mpairs)
      common/sumvalues/fc1,fc2,add1,add2,sumtot1,sumtot2,errort
      character(len=4) :: atom1,atom2
      character(len=9) :: titlepsf
      integer :: ilent1,ilent2,ilent

      do id=1,npairs
            atom1=tlab1(irefpcof(id))
            atom2=tlab2(irefpcof(id))
!      write(6,*) 'sumgr> fc1 and fc2: ',id,nz,npt1,npt2,' ',atom1,' ',atom2
            ilent1=index(atom1,' ')-1
            ilent2=index(atom2,' ')-1
            titlepsf=atom1(1:ilent1)//'-'//atom2(1:ilent2)
            ilent=ilent1+ilent2+1
            entsig(1,id)=0.0
            zerolimit(id)=0.0
            do j=2,nz
                  if(ycoh(j).gt.0.0) then
                        add1(1)=wav(1,j,id)/ycoh(j)
                        add1(2)=wav(2,j,id)/ycoh(j)
                  else
                        add1(1)=0.0
                        add1(2)=0.0
                  endif
                  sumtot1(1)=gorintra(j,id)
                  sumtot1(2)=gorinter(j,id)
                  sumtot2(1)=gr2intra(j,id)
                  sumtot2(2)=gr2inter(j,id)
!      write(6,*) 'sumgr> ',id,j,nz,sumtot1(1),sumtot1(2),sumtot2(1),sumtot2(2)
                  do k=1,2
                     add2(k)=add1(k)**2
                     sumtot2(k)=sumtot2(k)**2+sumtot1(k)**2
                     sumtot1(k)=fc2*sumtot1(k)+fc1*add1(k)
                     sumtot2(k)=fc2*sumtot2(k)+fc1*add2(k)
                     errort(k)=sumtot2(k)-sumtot1(k)**2
                     if(errort(k).gt.0.0) then
                        errort(k)=dsqrt(errort(k))
                     else
                        errort(k)=0.0
                     endif
                  end do
                  gorintra(j,id)=sumtot1(1)
                  gorinter(j,id)=sumtot1(2)
                  gr2intra(j,id)=errort(1)
                  gr2inter(j,id)=errort(2)
                  zerolimit(id)=zerolimit(id)+dble(ycoh(j)*(gorinter(j,id)-1.0))
                  entsig(j,id)=zerolimit(id)
            end do
!            write(6,'(a,i4,1x,a,1x,a,e12.5)') 'sumgr> Atom pair ',id,titlepsf(1:ilent),' zero limit: ',zerolimit(id)
      end do
      return
      end

c
c sum S(Q) and error bars
c
      subroutine sumsq(nq,qwind)
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      real*4 qwind(2,mpts,mpairs)

      real*8 fc1,fc2,add1(2),add2(2),sumtot1(2),sumtot2(2),errort(2)
      common/sumvalues/fc1,fc2,add1,add2,sumtot1,sumtot2,errort

      do id=1,npairs
            do j=1,nq
                  add1(1)=qwind(1,j,id)
                  add1(2)=qwind(2,j,id)
                  sumtot1(1)=soqintra(j,id)
                  sumtot1(2)=soqinter(j,id)
                  sumtot2(1)=sq2intra(j,id)
                  sumtot2(2)=sq2inter(j,id)
                  do k=1,2
                     add2(k)=add1(k)**2
                     sumtot2(k)=sumtot2(k)**2+sumtot1(k)**2
                     sumtot1(k)=fc2*sumtot1(k)+fc1*add1(k)
                     sumtot2(k)=fc2*sumtot2(k)+fc1*add2(k)
                     errort(k)=sumtot2(k)-sumtot1(k)**2
                     if(errort(k).gt.0.0) then
                        errort(k)=dsqrt(errort(k))
                     else
                        errort(k)=0.0
                     endif
                  end do
                  soqintra(j,id)=sumtot1(1)
                  soqinter(j,id)=sumtot1(2)
                  sq2intra(j,id)=errort(1)
                  sq2inter(j,id)=errort(2)
!                  write(25,'(a,2(1x,i7),6(1x,e13.6))') 'sumsq-1> ',id,j,add1(1),add1(2),sumtot1(1),sumtot1(2)
            end do
      end do
      return
      end

c
c sum the potential energies and write a file
c
      subroutine sumpe(nz,z,inpfile,ileninpfile)
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      character*256 inpfile,fname
      character*5 ext
      real*4 z(mpz)
      toteref=0.0
      toteemp=0.0
      totetot=0.0
      rmaxpot=0.5*dmax
      do id=1,npairs
            energyref(id)=0.0
            energyemp(id)=0.0
            energytot(id)=0.0
            it=n1extr(id)
            jt=n2extr(id)
            cprod=cfrac(it)*cfrac(jt)
            if(it.ne.jt) cprod=cprod*2
            j=2
c
c only integrate up to half the box dimension since g(r) may not be
c reliable beyond that distance
c
            do while(j.le.nz.and.z(j).le.rmaxpot)
                  wt=cprod*ycoh(j)*gorinter(j,id)
                  energyref(id)=energyref(id)+wt*ecoul(j,id)
                  energyemp(id)=energyemp(id)+wt*xavs(j,id)
                  j=j+1
            end do
            energytot(id)=energyref(id)+energytot(id)
            toteref=toteref+energyref(id)
            toteemp=toteemp+energyemp(id)
      end do
      totetot=toteref+toteemp
c
c write these results to a file
c
      write(ext,100)
100   format('.terg')
      fname=inpfile(1:ileninpfile)//ext(1:5)
      open(20,file=fname,status='unknown')
      write(20,200) toteref,toteemp,totetot
200   format(3(1x,e13.6))
      do id=1,npairs
            write(20,201) type(n1extr(id)),type(n2extr(id))
            write(20,200) energyref(id),energyemp(id),energytot(id)
      end do
201   format(2(1x,a3))
      close(20)
      return
      end

      subroutine initpt(ireset,rminpt,rmaxpt,expecf,polar,smoos,rminp,rbroad)
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      include 'weights.inc'
      real smoos(mfile),rminp(mfile),rminpt,rmaxpt,expecf,widthptc
      write(6,*) 'initpt> rminpt,rmaxpt: ',rminpt,rmaxpt
      pi=3.141592653
      if(power.gt.0.0) ratt=1.0e10**(1.0/power)
c
c maximum value of a sum - this it to prevent overflows
c
      summax=exp(25.0)
c
c generate renormalized x scale and potential damping function at low r
c
      do j=1,npt
            xstar(j)=x(j)/dmax
      end do
      widthptc=max(3.5,0.2*rmaxpt)
      do j=1,npt
         xx=x(j)

c Lennard-Jones and EP truncation function

         wfac(j)=truncate(xx,rminpt,rmaxpt)

c Coulomb truncation function

!         wfacc(j)=truncatecweeks(xx,widthptc)
         wfacc(j)=truncatec(xx,rmaxpt)

      end do
      wfac(npt)=0.0
      wfacc(npt)=0
c
c set the empirical potentials to zero, and test for the maximum distance required for the potential
c
      do ip=1,npairs
            irefp=irefpcof(ip)
            do j=1,npt
                  if(x(j).le.abs(rminex(irefp)+4.0*expecf).and.x(j).gt.x(nptmax)) nptmax=j
                  xav(j,ip)=0.0
                  xavs(j,ip)=0.0
                  xav1(j,ip)=0.0
                  xavs1(j,ip)=0.0
            end do
      end do
c
c Generate reference potential from parameters
c
      call generate_rp()
c
c generate the EP from the site-site coefficients
c
      call generate_ep(.false.)

      return
      end

      subroutine ptgnsq(nsum,nq,q,wav,data,qwind
     *,rminpt,rmaxpt,szeros,overallfac,gfac
     *,rminp,potnoi,wtentp
     *,rbroad,sigma0,nsmoop)
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      include 'weights.inc'
      include 'X-rayf0.inc'

      real*4 q(mpts),wav(2,mpz,mpairs),data(mpts,mfile)
     1,qwind(mpts,mpairs),qtemp(mpts,mfile),fittemp(mpts,mfile)
     1,dt(mpz),dts(mpz),potnoi(mfile)
     1,szeros(mfile),rminp(mfile),rfactorep(mfile)

      real qminptemp

      integer*4 nmin(mfile),nmini(mfile)
      character*14 formatout

      summax=exp(25.0)
      pi=4.0*atan(1.0)
! Set the minimum Q for fitting the potential
      if(qminp.le.0.0) then
         qminptemp=qminm
      else
         qminptemp=qminp
      end if
      stepcl=cellst
      xavmax=0.0
      qstep=q(3)-q(2)
      qrat=2.0*pi/qstep
      expmax=25.0
      potadd=potfac !Limit the amplitude of the perturbation
      do id=1,nfile
            nmin(id)=2
            nmini(id)=2

c Form the simulated fit for this dataset

c If X-ray data, need to apply the Q-dependent weights

            ix=xfile(id)

c Set up the gradient for extrapolating to Q=0

            if(szeros(id).ne.0.0.and.nqmin(id).gt.1) then

               grads0=(data(nqmin(id),id)-szeros(id))/(q(nqmin(id))-q(1))
               conss0=szeros(id)

            else

               grads0=0.0
               conss0=0.0

            endif

            do i=1,nq
c
c Sum the intermolecular partials for this dataset at this q value
c
               sum=0.0
               if(ix.gt.0) then
                  do ic=1,npairs
                     sum=sum+weightx(i,ic,ix)*qwind(i,ic)
                  end do
               else
                  do ic=1,npairs
                     sum=sum+cwt(id,ic)*qwind(i,ic)
                  end do
               endif

c Calculate the difference between data and fit. For Q-values below qmin we take the first
c valid data point as the data below Qmin, unless szeros is non-zero.

               if(q(i).lt.max(qminptemp,q(nqmin(id)))) then
                  if(szeros(id).ne.0.0) then
                     qtemp(i,id)=conss0+grads0*(q(i)-q(1))-sum
                  else
                     qtemp(i,id)=0.0
                  end if
               else if(i.ge.nqmin(id).and.i.le.nqmax(id)) then
                  qtemp(i,id)=data(i,id)-sum
               else
                  qtemp(i,id)=0.0
               endif

            end do

      end do

c Format for diagnostic files

      write(formatout,332) '(',(2*mfile+1),'(1x,e12.5))'
332   format(a,i2.2,a)

c If X-ray data are present we need to reconstruct the differences using the Q=0
c weights. This is to avoid having to deconvolute the X-ray form factors in r-space.
c This is possibly a controversial practice since it involves using the inversion
c matrix and therefore eliminates any aspects of the data that are incompatible between
c datasets. For compatibility the same is done for datasets even when there are no
c X-ray data present.
c
c     if(nxfile.gt.0) then
c
c First sum the existing differences to make differences in the partial structure
c factors.
c
c        iqinv=2
c        open(10,file='deldat.dat',status='unknown')
c        open(11,file='delpsf.dat',status='unknown')
c        do iq=1,nq
c
c Ensure the current q values of the inversion matrix correspond to this Q value
c
c           do while(qinv(iqinv).lt.iq.and.iqinv.lt.nqinv)
c                 iqinv=iqinv+1
c           end do
c     write(6,*) iq,qinv(iqinv),iqinv
c
c Interpolation ratio
c
c           ratinter=real(iq-qinv(iqinv-1))/real(qinv(iqinv)-qinv(iqinv-1))
c           ratinter1=1.0-ratinter
c
c Correct for multiplication by feedback
c
c           if(feedback.gt.0.0) then
c              ratinter=ratinter/feedback
c              ratinter1=ratinter1/feedback
c           endif
c
c Calculate the partial differences at this Q value
c
c           do ic=1,npairs
c              sum=0.0
c              do id=1,nfile
c
c Interpolate (linearly) inversion matrix to this Q value
c
c                 cwttemp=ratinter1*cwttq(iqinv-1,ic,id)
c     1+ratinter*cwttq(iqinv,ic,id)
c                 sum=sum+cwttemp*qtemp(iq,id)
c              end do
c
c              delpsf(ic)=sum
c
c           end do
c
c           write(11,fmt=formatout) q(iq),(delpsf(ic),ic=1,npairs)
c
c Re-sum these partial structure factor differences to give the differences as
c they would be at Q=0
c
c           do id=1,nfile
c              sum=0.0
c              do ic=1,npairs
c                 sum=sum+cwt(id,ic)*delpsf(ic)
c              end do
c              deldat(id)=sum
c           end do
c           write(10,fmt=formatout) q(iq),(qtemp(iq,id),deldat(id),id=1,nfile)
c           do id=1,nfile
c              qtemp(iq,id)=deldat(id)
c           end do
c        end do
c        close(10)
c        close(11)
c     else
c
c First sum the existing differences to make differences in the partial structure
c factors.
c
c        open(10,file='deldat.dat',status='unknown')
c        open(11,file='delpsf.dat',status='unknown')
c        do iq=1,nq
c
c Interpolation ratio
c
c           ratinter1=1.0
c
c Correct for multiplication by feedback
c
c           if(feedback.gt.0.0) then
c              ratinter1=ratinter1/feedback
c           endif
c
c Calculate the partial differences at this Q value
c
c           do ic=1,npairs
c              sum=0.0
c              do id=1,nfile
c
c Interpolate (linearly) inversion matrix to this Q value
c
c                 cwttemp=ratinter1*cwttq(1,ic,id)
c                 sum=sum+cwttemp*qtemp(iq,id)
c              end do
c
c              delpsf(ic)=sum
c
c           end do
c
c           write(11,fmt=formatout) q(iq),(delpsf(ic),ic=1,npairs)
c
c Re-sum these partial structure factor differences to give the differences as
c they would be at Q=0
c
c           do id=1,nfile
c              sum=0.0
c              do ic=1,npairs
c                 sum=sum+cwt(id,ic)*delpsf(ic)
c              end do
c              deldat(id)=sum
c           end do
c           write(10,fmt=formatout) q(iq),(qtemp(iq,id),deldat(id),id=1,nfile)
c           do id=1,nfile
c              qtemp(iq,id)=deldat(id)
c           end do
c        end do
c        close(10)
c        close(11)
c     endif
c

      if(feedback.gt.0.0.and.potadd.gt.0.0) then
            do id=1,nfile
c
c determine the minimum bin for the intermolecular potential
c
                  do j=1,npt2
                     if(x(j).lt.rminp(id)) nmini(id)=j+1
                  end do
                  write(6,59) nmini(id)
59    format(1x,'Minimum radius integer of potential:',(1x,i5))
c
c parameters for generating the potential
c
c nfreq represents the number of times during the fit the coefficients are
c smoothed
c
                  nfreq=10
                  nterms=npcoeff(id)
                  stepp=pstepd(id)
c
c decay=sigma2(id) controls the width of the Gaussians used to fit the data.
c Typical value should be of order 2*pi/qmax, or a bit less. Obviously reducing
c sigma2 will produce more detail in the fitted s(q)/g(r), but may lead to
c greater possibility of truncation artifacts being reproduced.
c
                  decayp=sigma2(id)
                  acceptr=paccept
c
c set up the transform matrix from power-exponentials to Q-space
c
                  do ic=1,nterms
                        pwtold(ic)=pcoefd(ic,id)
                        pwt(ic)=pwtold(ic)
                  end do
                  if(gaussian) then
                     call calc_gaussq(nqmax(id),nterms,pdstep,pdmax
     *,decayp,rho,q,expsq,0.0,ncfirst)
                  else
                     call calc_pwrexpq(nqmax(id),nterms,pdmin,pdstep,pdmax
     *,decayp,rho,q,expsq,0.0,ncfirst)
                  end if
c
c calculate perturbation to potential
c
c form the difference between data and fit.
c
                  do i=1,nqmax(id)
                        dt(i)=-(qtemp(i,id))
                  end do
c
c fit a series of power-exponentials in q-space to this difference
c
                  if(szeros(id).eq.0.0) then
                     call fit_expansion(dt,expsq,dts,pwt,pwtold
     1,stepp,acceptr,rejrat,chisqnew,1,nqmax(id),1,nterms
     2,npits,nfreq)
                  else
                     call fit_expansion(dt,expsq,dts,pwt,pwtold
     1,stepp,acceptr,rejrat,chisqnew,1,nqmax(id),1,nterms
     2,npits,nfreq)
                  end if
c Save fit

                  do i=1,nqmax(id)
                     fittemp(i,id)=dts(i)
                  end do

c                 write(6,*) decayp,chisqnew
c
c save revised step and chisq values
c
                  pstepd(id)=stepp
                  potnoi(id)=chisqnew
c
c save the new fluctuation coefficients and add them to the empirical
c potential coefficients
c
                  rate=simt
                  pfact=potadd
                  do ic=2,nterms
c
c save the new fluctuation coefficients
c
                        pcoefd(ic,id)=pwtold(ic)
c
c increment the site-site potential coefficients
c
                        do ip=1,npairs
                              irefp=irefpcof(ip)
                              pcoefp(ic,irefp)=pcoefp(ic,irefp)
     *+pfact*efacp(irefp)*cwtpot(ip,id)*pwtold(ic)
                        end do
                  end do
                  rfactorep(id)=0.0
                  if(szeros(id).eq.0.0) then
                     do i=nqmin(id),nq
                        if(abs(qtemp(i,id)).gt.0.0) then
                           rfactorep(id)=rfactorep(id)+abs(fittemp(i,id)+qtemp(i,id))/abs(qtemp(i,id))
                        else
                           rfactorep(id)=rfactorep(id)+abs(fittemp(i,id)+qtemp(i,id))
                        end if
                     end do
                     rfactorep(id)=rfactorep(id)/real(nq-nqmin(id)+1)
                  else
                     do i=1,nq
                        if(abs(qtemp(i,id)).gt.0.0) then
                           rfactorep(id)=rfactorep(id)+abs(fittemp(i,id)+qtemp(i,id))/abs(qtemp(i,id))
                        else
                           rfactorep(id)=rfactorep(id)+abs(fittemp(i,id)+qtemp(i,id))
                        end if
                     end do
                     rfactorep(id)=rfactorep(id)/real(nq)
                  end if
                  write(6,*) 'ptgnsq> ',id,potadd,pfact,rfactorep(id)
            end do
            open(10,file='delfit.dat',status='unknown')
            do i=1,nq
               if(abs(qtemp(i,id)).gt.0.0) then
                  write(10,fmt=formatout) q(i),(fittemp(i,id),-qtemp(i,id),abs(fittemp(i,id)+qtemp(i,id))/abs(qtemp(i,id)),id=1,nfile)
               else
                  write(10,fmt=formatout) q(i),(fittemp(i,id),-qtemp(i,id),abs(fittemp(i,id)+qtemp(i,id)),id=1,nfile)
               end if
            end do
            close(10)

c Regenerate EP

!           call generate_ep(.true.)
            call generate_ep(.false.)

      endif

c
c if no. of iterations to calculate potential is large (ireset=1) reset
c this value to input value for all subsequent iterations
c
      npits=npitss
      return
      end

      subroutine repulsivept(potadd)
      !Regenerates the reference potential including new values of the exponential repulsive terms
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'

      real*4 pchange
c
c modify the repulsive terms as required
c
      if(potadd.gt.0.0) then
            pfact=potadd
            do ip=1,npairs
                  irefp=irefpcof(ip)
                  write(6,*) 'repulsivept potex> ',potex(ip),coordex(irefp)
                  if(rminex(irefp).gt.0.0) then
                     pchange=potex(ip)-coordex(irefp)
! Limit the value of pcoefex becoming absurdly small
                     if(pchange.eq.0.0.and.abs(pcoefex(irefp)).gt.0.0001) then
c
c if potex(ic) is zero, then reduce the current pcoefex(ic) by 10% to
c ensure this term does not dominate the simulation
c
                        potex(ip)=-0.1
                        pcoefex(irefp)=pcoefex(irefp)*(1.0+potex(ip))
                     else
                        !Restrict changes to less than 0.5% of simt
                        if(coordex(irefp).gt.1.0) pchange=pchange/coordex(irefp) !For larger coordination numbers, the change is expressed as a fraction of the requested coordination number
                        pchange=sign(min(0.005*simt,0.025*simt*abs(pchange)),pchange)
 !                       pchange=sign(0.001*simt,pchange)
                        pcoefex(irefp)=pcoefex(irefp)+pchange
                     end if
                  end if
! Repeat for pcoefexmax
                  pcoefexmax(irefp)=potexmax(ip)
!                  if(pcoefexmax(irefp).ge.0.0) then
! Limit the value of pcoefexmax becoming absurdly small
!                     if(potexmax(ip).eq.0.0.and.abs(pcoefexmax(irefp)).gt.0.0001) then
c
c if potexmax(ic) is zero, then reduce the current pcoefexmax(ic) by 10% to
c ensure this term does not dominate the simulation
c
!                        potexmax(ip)=-0.1*pfact
!                        pcoefexmax(irefp)=pcoefexmax(irefp)*(1.0+potexmax(ip))
!                     else
!                        pcoefexmax(irefp)=pcoefexmax(irefp)+potexmax(ip)
!                     end if
!                  end if
            end do
c
c regenerate the reference potential
c
            call generate_rp()

      endif

      return
      end

c
c write a new input file
c
      subroutine norma(ireset,itimes,ntimes)
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'epsr16.inc'
      include 'epsrinput_par.inc'
      include 'molecom.inc'
      include 'potential_par.inc'
      include 'weights.inc'
      include 'iresetnumbers.inc'
      character*256 fname
      character*10 extfile
C
C WRITE CURRENT DATA TO INPUT FILE IN A FORM SUITABLE TO BE READ AS NEXT
C INPUT FILE
C
      if(itimes.eq.ntimes) then
            extfile='.inpa'
            fname=fnamerun(1:ilenfnamerun)//extfile(1:5)
            call OPENf(14,fname)
            write(14,*) nisum,nesum,nesumt,ensum,ensum2,prsum,prsum2
c
c Write the data potential parameters
c
            do id=1,nfile
                  write(14,*) id
                  write(14,*) npcoeff(id),pstepd(id),sigma2(id)
                  write(14,*) (pcoefd(ic,id),ic=1,npcoeff(id))
            end do
            nnew=nfile+npairs
            write(14,*) feedback,(efacm(j),j=1,nfile),nqinv
            do iq=1,nqinv
            do ip=1,npairs
                  write(14,*) qinv(iq),n1extr(ip),n2extr(ip)
                  write(14,*) (cwttq(iq,ip,id),id=1,nnew)
            end do
            end do
            close(14)
c
c write the site-site potential parameters
c
c first save the current internal values to their symbol counterparts
c
            iref=ireff(iqminref,1)
            iref=iref+1
            write(symbol_value(iref),*) roverlap
            iref=iref+1
            write(symbol_value(iref),*) rminfac
            iref=iref+1
            write(symbol_value(iref),*) rminpt
            iref=iref+1
            write(symbol_value(iref),*) rmaxpt
            iref=iref+1
            write(symbol_value(iref),*) rbroad
            iref=iref+1
            write(symbol_value(iref),'(f7.3,2x,f7.3)') expecf,expecfpwr
! reppottype has been assigned at input, so skip over it
            iref=iref+1
! addpotential has been assigned at input, so skip over it
            iref=iref+1
            iref=iref+1
            write(symbol_value(iref),*) rcharge
            iref=iref+1
            write(symbol_value(iref),*) power
            iref=iref+1
            write(symbol_value(iref),*) ncoeffp
            iref=iref+1
            write(symbol_value(iref),*) pdmax
            iref=iref+1
            write(symbol_value(iref),*) pdstep
            iref=iref+1
            write(symbol_value(iref),*) npitss
            iref=iref+1
            write(symbol_value(iref),*) paccept
! Gaussian has been assigned at input, so skip over it
            iref=iref+1
            iref=iref+1
            write(symbol_value(iref),'(f8.5)') psigma2read
            iref=iref+1
            write(symbol_value(iref),'(f8.5)') refpotfac
c
c Write the uniform distribution to a file called EPSR.uni
c
            extfile='.uni'
            fname=fnamerun(1:ilenfnamerun)//extfile(1:4)
            call openf(10,fname)
            write(10,'(2(1x,e13.6))') (x(j),ycoh(j),j=1,npt1)
            close(10)
      endif
      return
      end

      subroutine openf(nunit,fname)
      character*256 fname
      character*10 stat
      stat='unknown'
      open(nunit,file=fname,status=stat,iostat=ierr)
      if(ierr.eq.0) return
      stat='new'
      open(nunit,file=fname,status=stat,iostat=ierr)
      if(ierr.eq.0) return
      write(6,100) fname,stat
100   format(1x,' Open error in OPENF',1x,a50,1x,a10)
      stop
      end

      subroutine get_wts(q,nq,fnamew,ireset)

      use molecom_mod
c
c assign neutron weightings to individual atom pairs - these are initially
c all set to zero
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      include 'weights.inc'
      include 'X-rayf0.inc'
      real wread(2),q(mpts)
      real*8 axt(5),bxt(5),cxt
c weights on individual partials, their transpose, and the cross matrix
      character*256 fnamew(mfile)
      character*3 a1,a2
      integer nq,normx,iafformaff
c
c set some accumulators to zero
c
      do ic=1,npairs
            uwt(ic)=0.0
      end do

c Set the X-ray form factors to zero

      call wxinit(nq)

c Number of X-ray files. This should not exceed mxfile, otherwise only the last one input will
c be used.
      nxfile=0

      do in=1,nfile
            do ic=1,npairs
                  do j=1,2
                        weight(j,in,ic)=0.0
                  end do
            end do
            open(10,file=fnamew(in),status='old')

c decide whether this is an X-ray or neutron wts file

            ilenXray1=index(fnamew(in),'X-ray.wts')
            ilenXray2=index(fnamew(in),'x-ray.wts')
            ilenXray3=index(fnamew(in),'XWTS.wts')
            ilenXray=max(ilenXray1,ilenXray2,ilenXray3)

            if(ilenXray.gt.0) then

c
c read factor for dataset to convert it to per atom. For X-rays this factor will decide
c whether data are normalised to the single atom scattering (1) or square of the mean
c scattering length (2).
c
                  read(10,'(a)') sometext
                  call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)

                  if(nwords.eq.3) then

                     read(sometext(ncf(1):ncl(nwords)),*) rhofac(in),normx,iafformaff

                  else

c input from old epsrwtsx

                     read(sometext(ncf(1):ncl(nwords)),*) rhofac(in)

c Deduce the values of rhofac(in), normx and iafformaff from this one value

c nrho + means data are per atom, -ve data are per molecule

                     if(rhofac(in).lt.0) then

                        nrho=2

                     else

                        nrho=1

                     endif

c iafformaff is 1 means use IAFFs for normalisation, otherwise use MAFFs

                     iafformaff=(abs(rhofac(in))+1)/2

c In the old version normx = 1 means normalise to single atom scattering, otherwise normalise to sum of Fs
c all squared. normx = 0 was not an option

                     normx=abs(rhofac(in))-2*(iafformaff-1)

c However in the new version, normx = 1 means normalise to sum of Fs all squared, and normx = 2 means normalise to single atom scattering

                     if(normx.eq.1) then
                        normx=2
                     else
                        normx=1
                     endif
         write(6,*) 'get_wts> ',rhofac(in),nrho,normx,iafformaff

                     rhofac(in)=nrho

                  endif

                  if(nxfile.lt.mxfile) then
                     nxfile=nxfile+1

c Save the normalisation of the X-ray data.

                     datanorm(nxfile)=rhofac(in)
                     rhofac(in)=1.0
                     normxx(nxfile)=normx
                     iafformaffx(nxfile)=iafformaff
      write(6,*) 'get_wts> ',nxfile,datanorm(nxfile)

c Setup pointers to the X-ray form factor arrays

                     xfile(in)=nxfile

c And to the dataset these weights correspond to.

                     datapoint(nxfile)=in

      write(6,*) 'get_wts> ',in,nxfile,xfile(in),datapoint(nxfile)

c Read form factor coefficients for this atom.

               read(10,108,iostat=ierr) a1,a2
108   format(1x,a3,1x,a4)
                     read(10,*,iostat=ierr) cxt,(axt(i),bxt(i),i=1,5),alphaxt,deltaxt
               do while(ierr.eq.0)

c Assign these X-ray coefficients if the specified atom type matches an atom in the
c .ato file.  Any atoms which are not assigned form factor coefficients are
c automatically assumed not to contribute to this diffraction pattern.
c
                        call wxtype(a1,axt,bxt,cxt,alphaxt,deltaxt,nxfile,q,nq)
                  read(10,108,iostat=ierr) a1,a2
                        read(10,*,iostat=ierr) cxt,(axt(i),bxt(i),i=1,5),alphaxt,deltaxt

                     end do

                  else

                     write(6,109) mxfile
109   format(/'get_wts> Exceeded limit of ',i2,' X-ray files.')

                  endif

            else

                  read(10,*) rhofac(in)
                  xfile(in)=0

c Read partial structure factor assignment and neutron weighting for
c intramolecular and inter molecular terms.
c
                  read(10,107,iostat=ierr) a1,a2
                  read(10,*,iostat=ierr) wread(1),wread(2)
107   format(1x,a3,1x,a3)
                  do while(ierr.eq.0)
c
c assign this neutron weighting if it matches a pair of atoms in the molecule
c Pairs which are not assigned a weighting by this input are left as zero.
c
                     call wtype(a1,a2,wread,in)
                     read(10,107,iostat=ierr) a1,a2
                     read(10,*,iostat=ierr) wread(1),wread(2)

                  end do

            endif
            close (10)
      end do

c Prepare the X-ray weighting factors if present.

      if(nxfile.gt.0) then
            call wxcreate(q,nq)
      end if

      do in=1,nfile

c sum up intermolecular weights and modulus of these for this file

            wsum(in)=0.0
            xwt(in)=0.0
            do ic=1,npairs
                  it=n1extr(ic)
                  jt=n2extr(ic)
                  cprod=cfraccg(it)*cfraccg(jt)
                  if(it.ne.jt) cprod=2.0*cprod
                  ciwt(in,ic)=cprod*weight(1,in,ic)
                  cwt(in,ic)=cprod*weight(2,in,ic)
                  wsum(in)=wsum(in)+abs(cwt(in,ic))
                  xwt(in)=xwt(in)+cwt(in,ic)
            end do
      end do

c Repeat this for the X-ray weights if present

      if(nxfile.gt.0) then

         do ix=1,nxfile
            do ip=1,npairs
               it=n1extr(ip)
               jt=n2extr(ip)
               cprod=cfraccg(it)*cfraccg(jt)
               if(it.ne.jt) cprod=2.0*cprod
               do iq=1,nq
                  weightx(iq,ip,ix)=cprod*weightx(iq,ip,ix)
               end do
            end do
         end do

c Set up the Q-values on which to calculate the inversion matrix.

         qinvstep=2.0
         nqinvstep=int(qinvstep/(q(2)-q(1)))

      else

c If no X-ray data, only need to calculate inversion at one Q value.

         nqinvstep=nq

      endif
c
c if ireset.eq.1, then form the inversion of the coefficients cwtpot
c
c Also calculate the inversion if any of the factors has changed
c
      iflag=0
      if(ireset.eq.1.or.abs(feedback-feedbackold).gt.0.001) then
            iflag=1
      else
            in=1
            do while(iflag.eq.0.and.in.le.nfile)
                  if(abs(efacm(in)-efacold(in)).gt.0.001) iflag=1
                  in=in+1
            end do
      endif
      if(iflag.ne.0) then

c Set the initial values of the inversion matrix
         nnew=npairs+nfile
         do iff=1,nnew
            do ip=1,npairs
c
c for the partials from the simulation, assume the initial
c transform matrix is diagonal, as it would be if there were no
c data
c
               ifref=iff-nfile
               if(ifref.eq.ip) then
                  cwtt(ip,iff)=1.0
               else
                  cwtt(ip,iff)=0.0
               endif
            end do
         end do
c
c generate weights to be used for potential inversion at this Q value
c

         nqinv=0
         do iqinv=1,nq,nqinvstep

            if(nqinv.lt.mqinv) then
               nqinv=nqinv+1
               qinv(nqinv)=iqinv

               write(6,*) 'get_wts> ',nqinv,q(qinv(nqinv)),iqinv,q(iqinv)

               do in=1,nfile

                  ix=xfile(in)
                  if(ix.gt.0) then
                     do ic=1,npairs
                        it=n1extr(ic)
                        jt=n2extr(ic)
                        if(cfraccg(it).gt.0.0.and.cfraccg(jt).gt.0.0) then
                           cprodrat=cfrac(it)*cfrac(jt)/(cfraccg(it)*cfraccg(jt))!Correct for coarse graining if present
                        else
                           cprodrat=0.0
                        end if
                        cwtpot(in,ic)=cprodrat*weightx(iqinv,ic,ix)
                     end do
                  else
                     do ic=1,npairs
                        it=n1extr(ic)
                        jt=n2extr(ic)
                        if(cfraccg(it).gt.0.0.and.cfraccg(jt).gt.0.0) then
                           cprodrat=cfrac(it)*cfrac(jt)/(cfraccg(it)*cfraccg(jt))!Correct for coarse graining if present
                        else
                           cprodrat=0.0
                        end if
                        cwtpot(in,ic)=cprodrat*cwt(in,ic)
                     end do
                  endif
               end do

               call invert_cwt(cwtpot,cwtt,carray,carrayold,feedback
     1,efacm,nfile,npairs,mfile,mpairs)

c Save the results into a big array!

               do in=1,nnew
                  do ic=1,npairs
                     cwttq(nqinv,ic,in)=cwtt(ic,in)
                  end do
               end do
            endif
         end do
      endif

c For the potential inversion coefficients, halve the unlike terms to avoid doubling the
c magnitude of the Empirical Potential.

      do in=1,nfile
         do ic=1,npairs
            cwtpot(ic,in)=cwttq(1,ic,in)
            if(n1extr(ic).ne.n2extr(ic)) then
               cwtpot(ic,in)=0.5*cwtpot(ic,in)
            endif
         end do
      end do
      return
      end

      subroutine wtype(a1,a2,wread,in)
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'potential_par.inc'
      character*3 a1,a2
      dimension wread(2)
      n1=0
      n2=0
      do i=1,ntype
            if(a1.eq.type(i)) n1=i
            if(a2.eq.type(i)) n2=i
      end do
      if(n1.gt.0.and.n2.gt.0) then
            icref=npairref(n1,n2)
            do j=1,2
               weight(j,in,icref)=wread(j)
            end do
      else
            write(6,100) in
100   format(1x,'Mismatched types between atomic coordinates and'
     1,'neutron weights for file ',i5)
      endif
      return
      end

      subroutine wxinit(nq)

c Initialise the X-ray form factors for all components to zero

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'potential_par.inc'
      include 'X-rayf0.inc'
      integer nq
      do ix=1,mxfile
      do it=1,ntype
         do iq=1,nq
            f0x(iq,it,ix)=0.0
         end do
      end do
      end do
      return
      end

      subroutine wxtype(a1,axt,bxt,cxt,alphaxt,deltaxt,in,q,nq)

c Set up the X-ray form factor for the specified atom type

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'potential_par.inc'
      include 'X-rayf0.inc'
      character*3 a1
      real*8 axt(5),bxt(5),cxt,sum,term,qx,qx2,one,pi4
      real q(mpts)
      integer*4 nq
      n1=0
      one=1.0
      pi4=16.0*atan(one)
      do i=1,ntype
         if(a1.eq.type(i)) n1=i
      end do
      if(n1.gt.0) then
         do iq=1,nq
            qx=q(iq)/pi4
            qx2=qx*qx
            sum=cxt
            do i=1,5
               term=-bxt(i)*qx2
               if(term.gt.-40.0.and.term.lt.+40.0) then
                  sum=sum+axt(i)*exp(term)
               endif
            end do
            f0x(iq,n1,in)=sum

c Multiply by the modification function according to
c equation (5) of JCP, vol 113, p 9149 (2000)

            if(alphaxt.ne.0.0.and.deltaxt.gt.0.0) then

               expon=q(iq)*q(iq)/(2.0*deltaxt*deltaxt)
               if(expon.lt.20.0) then
                  fac=1.0-(alphaxt/f0x(1,n1,in))*exp(-expon)
               else
                  fac=1.0
               endif
            else
               fac=1.0
            endif
            fm0x(iq,n1,in)=fac*sum
         end do
      else
         write(6,100) in
100   format(1x,'Mismatched types between atomic coordinates and'
     1,'X-ray weights for file ',i5)
      endif
      return
      end

      subroutine wxcreate(q,nq)

c Set up the X-ray form factor for the specified atom type

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      include 'weights.inc'
      include 'X-rayf0.inc'
      real q(mpts),amolrat
      real*8 bbar
      integer nq

c Get the ratio of atoms to molecules in case it is needed

      amolrat=real(natomr)/real(nmol)

      do ix=1,nxfile

c Determine the type of normalisation

c nrho 1 means data are per atom, 2 data are per molecule

         nrho=datanorm(ix)

c iafformaff is 1 means use IAFFs for normalisation, otherwise use MAFFs

         iafformaff=iafformaffx(ix)

c normx = 0 means no normalisation, normx = 1 means normalise to <F>^2,
c normx = 2 means normalise to <F^2>

         normx=normxx(ix)
         write(6,*) 'wxcreate> ',datanorm(ix),nrho,normx,iafformaff

c Form the Q-dependent weights for each atom pair

         do iq=1,nq

c Form the single atom term

            singlem=0.0
            singlei=0.0
            bbar=0.0
            do it=1,ntype

c Set up the single atom scattering for the modified and independent atom form factors

               if(qtypea(it).eq.0) then
                  singlem=singlem+cfraccg(it)*fm0x(iq,it,ix)**2
                  singlei=singlei+cfraccg(it)*f0x(iq,it,ix)**2
                  if(iafformaff.ne.2) then
                     if(normx.eq.0) then
                        bbar=bbar+cfraccg(it)
                     else if(normx.eq.1) then
                        bbar=bbar+cfraccg(it)*f0x(iq,it,ix)
                 else if(normx.eq.2) then
                        bbar=bbar+cfraccg(it)*f0x(iq,it,ix)**2
                     endif
                  else
                     if(normx.eq.0) then
                        bbar=bbar+cfraccg(it)
                     else if(normx.eq.1) then
                        bbar=bbar+cfraccg(it)*fm0x(iq,it,ix)
               else if(normx.eq.2) then
                        bbar=bbar+cfraccg(it)*fm0x(iq,it,ix)**2
                     endif
                  endif
               endif
            end do

            if(normx.eq.1) bbar=bbar*bbar

            do ip=1,npairs

               n1=n1extr(ip)
               n2=n2extr(ip)
               weightx(iq,ip,ix)=fm0x(iq,n1,ix)*fm0x(iq,n2,ix)/bbar

c Correct for per molecule data if needed

               if(nrho.eq.2.and.normx.eq.1)
     *weightx(iq,ip,ix)=weightx(iq,ip,ix)/amolrat

            end do

c Calculate background shift to simulated structure factor

            if(iafformaff.eq.1) then
               backx(iq,ix)=(singlem-singlei)/bbar
            else
               backx(iq,ix)=0.0
            endif
            if(nrho.eq.2.and.normx.eq.1) backx(iq,ix)=backx(iq,ix)/amolrat

         end do

c Set the Q=0 weights

         do ip=1,npairs
            weight(1,datapoint(ix),ip)=weightx(1,ip,ix)
            weight(2,datapoint(ix),ip)=weight(1,datapoint(ix),ip)
         end do

c Write an output file to check that all is well

         open(10,file='X-raywts.dat',status='unknown')
         do iq=1,nq
            write(10,100) q(iq),(weightx(iq,ip,ix),ip=1,npairs)
100   format(13(1x,e13.6))
         end do
         close(10)

c Write an output file to show the X-ray form factors used

         open(10,file='X-rayf0x.dat',status='unknown')
         do iq=1,nq
            write(10,100) q(iq),(fm0x(iq,it,ix),it=1,ntype)
         end do
         close(10)

      end do
      return
      end

      subroutine gensq(nq,q,sofq,wav)
      use molecom_mod
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      real q(mpts),sofq(2,mpts,mpairs),wav(2,mpz,mpairs),sumwavself(mpairs)
!      open(25,file='dump-25.txt',status='unknown')
!      write(6,*) 'gensq> ',maxnscalefactype
      if(coarsegraining) then
         allocate(prodssff(maxnscalefactype,maxnscalefactype))
      end if
! Set the fudge factors on the self terms
      ycohfactor=1.0
      sumwavself=0.0
      do ic=1,npairs
         it=n1extr(ic)
         jt=n2extr(ic)
         do j=2,npt1
            sumwavself(ic)=sumwavself(ic)+wav(1,j,ic)
         end do
         do i=1,nq
            sum1=0.0
            sum2=0.0
            sum3=0.0
            if(coarsegraining) then
! Calculate the product of form factors at this Q value
               do ibinia=1,nsffbintype(it)
                  do ibinja=1,nsffbintype(jt)
                     prodssff(ibinia,ibinja)=sphereformfactor(i,ibinia,it)*sphereformfactor(i,ibinja,jt)
                  end do
               end do
            end if
            do j=2,npt1
! First do intra-molecular terms
               add1=wav(1,j,ic)
               if(abs(add1).gt.0.0) then
                  sum1=sum1+add1*foqols(j,i)
               end if
            end do
!            if(coarsegraining) then
!               sum2=intraformfactor(i,ic) !This term is added to sum2 so that it gets the correct neutron weighting
!               if(it.eq.jt.and.cfraccg(it).gt.0.0) then
!                  sum1=intraformfactorself(i,it)/cfraccg(it) !Include the self term in this case
!!                  sum1=sum1+intraformfactorself(i,it) !Include the self term in this case
!               end if
!            end if
! Now inter-molecular terms
            do j=2,npt1
!               if(coarsegraining) then
!                  add2=0.0
!                  do ibinia=1,nsffbintype(it)
!                     do ibinja=1,nsffbintype(jt)
!                        term=wavinter(j,ic,ibinia,ibinja)-fracntermssff(ibinia,it)*fracntermssff(ibinja,jt)*ycoh(j)
!                        add2=add2+term*prodssff(ibinia,ibinja)
!                     end do
!                  end do
!                  write(6,*) 'gensq-1> ',it,jt,j,add2,fracntermssff(1,it),fracntermssff(1,jt),prodssff(1,1)
!                  sum2=sum2+add2*foqols(j,i)
!               else
                  add2=(wav(2,j,ic)-ycoh(j))*natomsphere(jt)
                  sum2=sum2+add2*foqols(j,i)
!               end if
               sum3=sum3+wav(2,j,ic)
!                  if(i.eq.2) then
!                     write(25,'(a,3(i7),2(1x,e12.5),5(1x,e12.5))') 'gensq-2> ',ic,i,j,ycoh(j),foqols(j,i),add1,add2,sum1,sum2,sum3
!                  end if
            end do
c
c in the event that there are no intermolecular correlations for a particular atom pair
c (e.g. only one molecule of that type in the box)
c we set the simulated structure factor to zero
c
            if(sum3.le.0.0) sum2=0.0
            sofq(1,i,ic)=sum1
            sofq(2,i,ic)=sum2
         end do
      end do
!      write(6,*) 'gensq> ',(intraformfactorself(1,it),it=1,ntype)
!      write(6,*) 'gensq> ',(intraformfactor(1,ic),ic=1,npairs)
!      write(6,'(a)') 'gensq> '
!      write(6,'((1x,i3,1x,a3,1x,a3,1x,f10.5))') (ic,type(n1extr(ic)),type(n2extr(ic)),sumwavself(ic),ic=1,npairs)
!      write(6,*) 'gensq> ',(sofq(1,1,ic),ic=1,npairs)
!      close(25)
      if(coarsegraining) then
         if(allocated(prodssff)) deallocate(prodssff)
      end if
      return
      end

      subroutine fit_expansion(yin,expr,yout,wt,wtold
     1,stepr,acceptr,rejrat,chisqnew,nmin,npt
     2,ncfirst,nterms,nits,nfreq)
! Subroutine arguments
      real :: yin(*),yout(*),wt(*),wtold(*),expr(*)
      real stepr,acceptr,rejrat,chisqnew
      integer, intent(in) :: nmin,npt,ncfirst,nterms,nits,nfreq
! Internal variables
      integer :: nstep,itry,ntry,its,npoints
      real :: wtmax,chisqold,rejfac,chisqstart,chisqfinish,chisqtest
c
c maximum value of a sum - this it to prevent overflows
c
      summax=exp(25.0)
      nstep=100
      ntry=150
      acceptr=0.05
      npoints=npt-nmin+1
c
c assign the initial coefficients
c
      wtmax=0
      do it=1,nterms
         wt(it)=wtold(it)
         wtmax=max(abs(wt(it)),wtmax)
      end do
      if(stepr.lt.1.0e-6*wtmax) stepr=2.0e-6*wtmax !Prevent step size getting too small!
      if(nits.gt.0) then
         chisqold=0.0
         icount=0
         chisqnew=0.0
         do i=nmin,npt
            term=yin(i)-expand(nterms,i,expr,wt,summax)
            chisqnew=chisqnew+term**2
         end do
         chisqold=chisqnew
         chisqtest=-1.0
         its=0
         itry=0
         do while (itry.lt.10.or.(itry.lt.ntry.and.chisqtest.lt.1.0e-4))
!         do while (itry.lt.ntry.and.abs(chisqtest).gt.1.0e-3.and.chisqtest.lt.0)
!         do while (itry.lt.ntry)
            itry=itry+1
            accept=acceptr*chisqold/real(npoints)
            rejrat=0.0
            chisqstart=chisqold/real(npoints)
            do it=1,nstep
               its=its+1
c
c select some new trial values
c
               do ic=ncfirst,nterms
                  if(abs(wtold(ic)).lt.1.0) then
                        wt(ic)=wtold(ic)
     1+stepr*(2.0*ran1()-1.0)
                  else
                        wt(ic)=wtold(ic)*(1.0
     1+stepr*(2.0*ran1()-1.0))
                  endif
               end do
               chisqnew=0.0
               do i=nmin,npt
                  term=yin(i)-expand(nterms,i,expr,wt,summax)
                  chisqnew=chisqnew+term*term
               end do
               exponent=(chisqnew-chisqold)/accept
               if(exponent.gt.0.0.and.exponent.le.20.0) then
                  test=exp(-exponent)-ran1()
               else if(exponent.gt.20.0) then
                  test=-1.0
               else
                  test=1.0
               endif
               if(test.ge.0.0) then
c
c perform 1 smoothing on the coefficients to prevent them becoming noisy
c every nfreq times a set of moves are accepted
c
                  icount=icount+1
                  if(icount.eq.nfreq) then
                        icount=0
                        call smoop(ncfirst,nterms,1,wt)
c
c regenerate chisqnew for the modified coeffs.
c
                        chisqnew=0.0
                        do i=nmin,npt
                            term=yin(i)-expand(nterms,i,expr,wt,summax)
                            chisqnew=chisqnew+term**2
                        end do
                  endif
c
c accept new values
c
                  chisqold=chisqnew
                  do ic=1,nterms
                        wtold(ic)=wt(ic)
                  end do
               else
                  rejrat=rejrat+1.0
               end if
            end do
!            chisqtest=0.0
!            do i=nmin,npt
!               term=yin(i)-expand(nterms,i,expr,wt,summax)
!               if(abs(yin(i)).gt.0.0) chisqtest=chisqtest+term*term/abs(yin(i))
!            end do
!            chisqtest=chisqtest/real(npt-nmin+1)
            rejrat=rejrat/real(nstep)
!          if(chisqnew.gt.0.0) chisqtest=chisqtest/chisqnew
            chisqfinish=chisqold/real(npoints)
            if(chisqstart.gt.0.0) then
               chisqtest=(chisqfinish-chisqstart)/chisqstart
            else
               chisqtest=(chisqfinish-chisqstart)
            endif
            if(rejrat.gt.0.0) then
               rejfac=0.75/rejrat
            else
               rejfac=2.0
            end if
            stepr=stepr*rejfac
!            write(6,*) 'fit_expansion> ',itry,chisqfinish,chisqtest,rejrat
         end do
         write(6,*) 'fit_expansion> ',itry,chisqfinish,chisqtest,rejrat
         do ic=1,nterms
            wt(ic)=wtold(ic)
         end do
      endif
      chisqnew=0.0
      do i=1,npt
         term=expand(nterms,i,expr,wt,summax)
         yout(i)=term
         if(i.ge.nmin) chisqnew=chisqnew+(yin(i)-term)**2
      end do
      return
      end

      function expand(nterms,j,expr,coeff,summax)
      real*4 coeff(*),expr(*)
      sum=0.0
      do i=1,nterms
            jref=(j-1)*nterms+i
            sum=sum+expr(jref)*coeff(i)
            if(abs(sum).gt.summax) sum=sign(summax,sum)
      end do
      expand=sum
      return
      end

      subroutine calc_gaussq(lpt,nterms,dstep,dmax
     *,sigma,rho,xin,expr,rmin,ncfirst)
      real*4 xin(*),expr(*)
      pi=3.141592653
c
c maximum value of exponent
c
      expmax=25.0
c
c set up array of coefficients
c
      do m=1,nterms
c
c Set the value of d
c
            d=m*dstep
c
c save the integer of the first non-zero coefficient
c
            if(d.lt.rmin) ncfirst=m+1
c
c generate value of sigma at this d value
c
            sigd=sigma
c
c calculate Q dependence of mth term
c
            do i=1,lpt
                  jref=(i-1)*nterms+m
                  q=xin(i)
                  qsigd=q*sigd
                  qs2=0.5*qsigd*qsigd
                  if(qs2.lt.expmax) then
                        qd=q*d
                        if(qd.gt.0.0) then
                              sinqd=sin(qd)/qd
                        else
                              sinqd=1.0
                        endif
                        expr(jref)=sinqd*exp(-qs2)
                  else
                        expr(jref)=0.0
                  endif
            end do
      end do
      return
      end

      subroutine calc_gaussr(lpt,nterms,dstep,dmax
     *,sigmab,rho,xin,expr,rmin,ncfirst,fbroad,rbroad)
      real*4 xin(*),expr(*)
      pi=4.0*atan(1.0)
      piroot=sqrt(0.5*pi)
      gfac=piroot/(4*pi*pi*rho)
c
c maximum value of exponent
c
      expmax=25.0
c
c set up array of coefficients
c
      do m=1,nterms
c
c calculate d dependence of mth term
c
            d=m*dstep
c
c save the integer of the first non-zero coefficient
c
            if(d.lt.rmin) ncfirst=m+1
c
c generate value of sigma at this d value
c
            sigd=fbroad*sigmab
            gfacs=gfac/sigd
            do i=1,lpt
                  jref=(i-1)*nterms+m
                  r=xin(i)
                  rmd=(r-d)/sigd
                  rpd=(r+d)/sigd
                  rmd2=0.5*rmd*rmd
                  rpd2=0.5*rpd*rpd
                  sum=0.0
                  if(r.gt.0.0.and.d.gt.0.0) then
                        gfact=gfacs/(r*d)
                        if(rmd2.lt.expmax) sum=sum+exp(-rmd2)
                        if(rpd2.lt.expmax) sum=sum-exp(-rpd2)
                  else
                        gfact=2.0*gfacs/(sigd*sigd)
                        if(rpd2.lt.expmax) sum=sum+exp(-rpd2)
                  endif
                  expr(jref)=gfact*sum
            end do
      end do
      return
      end

      subroutine calc_gaussr_grad(lpt,nterms,dstep,dmax
     *,sigmab,rho,xin,expr,rmin,ncfirst,fbroad,rbroad)
      real*4 xin(*),expr(*)
      pi=4.0*atan(1.0)
      piroot=sqrt(0.5*pi)
      gfac=piroot/(4*pi*pi*rho)
c
c maximum value of exponent
c
      expmax=25.0
c
c set up array of coefficients
c
      do m=1,nterms
c
c calculate d dependence of mth term
c
            d=m*dstep
c
c save the integer of the first non-zero coefficient
c
            if(d.lt.rmin) ncfirst=m+1
c
c generate value of sigma at this d value
c
            sigd=fbroad*sigmab
            gfacs=gfac/sigd
            sigd2=sigd*sigd
            do i=1,lpt
                  jref=(i-1)*nterms+m
                  r=xin(i)
                  rmd=(r-d)/sigd
                  rpd=(r+d)/sigd
                  rmd2=0.5*rmd*rmd
                  rpd2=0.5*rpd*rpd
                  sum=0.0
                  if(abs(r).gt.0.0.and.abs(d).gt.0.0) then
                        gfact=gfacs/(r*d)
                        if(rmd2.lt.expmax) sum=sum-(1.0/r+(r-d)/sigd2)*exp(-rmd2)
                        if(rpd2.lt.expmax) sum=sum+(1.0/r+(r+d)/sigd2)*exp(-rpd2)
                  else
                        gfact=0.0
                  endif
                  expr(jref)=gfact*sum
            end do
      end do
      return
      end

      subroutine calc_pwrexpq(lpt,nterms,dmin,dstep,dmax
     *,sigmab,rho,xin,expr,rmin,ncfirst)
      real xin(*),expr(*),d,dmin,dstep,dmax,dnext,rootd,sigmab,sigma
      pi=3.141592653
c
c maximum value of exponent
c
      expmax=25.0
c
c set up array of coefficients
c
      d=dstep
      rootd=1.0
      do m=1,nterms
c
c calculate Q dependence of mth term
c
!            d=m*dstep
c
c save the integer of the first non-zero coefficient
c
            if(d.lt.rmin) ncfirst=m+1
            sigma=rootd*sigmab
!            sigma=sigmab
            if(sigma.le.d.and.sigma.gt.0.0.and.m.gt.1) then
                  n=d/sigma-1
            else
                  n=0
            endif
            an2=n+2
            an2fac=1.0/an2
            do i=1,lpt
                  q=xin(i)
                  qsigma=q*sigma
                  qs2=qsigma*qsigma
                  onemqs2=1.0-qs2
                  onepqs2=1.0+qs2
                  onepqs2s=onepqs2*onepqs2
                  rtonepqs2=sqrt(onepqs2)
                  alpha=atan(qsigma)
                  alphan=n*alpha
                  if(qsigma.gt.0.0) then
                        term=2.0*cos(alphan)
     *+onemqs2*sin(alphan)/qsigma
                  else
                        term=an2
                  endif
                  pn=term*an2fac/(rtonepqs2**(4+n))
                  jref=(i-1)*nterms+m
                  expr(jref)=pn
            end do
            d=dnext(d,dmin,dstep,rootd)
      end do
      return
      end

      subroutine calc_pwrexpr(lpt,nterms,dmin,dstep,dmax
     *,sigmab,rho,xin,expr,rmin,ncfirst,fbroad,rbroad)
      use lngamma_mod
      real*4 xin(*),expr(*),d,dmin,dstep,dmax,dnext,rootd,sigmab,sigma
      pi=3.141592653
c
c maximum value of exponent
c
      expmax=25.0
      ntmax=dmax/sigmab+3
      call fctrln(ntmax)
c
c set up array of coefficients
c
      d=dstep
      rootd=1.0
      do m=1,nterms
c
c calculate r dependence of mth term
c
!            d=m*dstep
c
c save the integer of the first non-zero coefficient
c
            if(d.lt.rmin) ncfirst=m+1
c
c generate value of sigma at this d value
c
            sigma=rootd*fbroad*sigmab
            if(sigma.le.d.and.sigma.gt.0.0) then
                  n=d/sigma-1
            else
                  n=0
            endif
            nref=n+3
            gfac=4.0*pi*rho*sigma*sigma*sigma
            gfacln=alog(gfac)
            sigln=gamln(nref)+gfacln
c
c for the first term we add a factor exp(rmin/sigma)
c This is to ensure low  r repulsive potential develops quickly
c
            if(m.eq.1) sigln=sigln-rmin/sigma
            do i=1,lpt
                  jref=(i-1)*nterms+m
                  r=xin(i)
                  fact=rbroad*r
                  rsigma=r/sigma
                  if(r.gt.0.0) then
                        rsigmaln=alog(rsigma)
                        expon=n*rsigmaln-rsigma-sigln-fact
                  else
                        if(n.eq.0) then
                              expon=-rsigma-sigln-fact
                        else
c
c for n>0, the function is identically zero at r=0
c
                              expon=-2.0*expmax
                        endif
                  endif
                  if(expon.gt.-expmax) then
                        expr(jref)=exp(expon)
                  else
                        expr(jref)=0.0
                  endif
            end do
            d=dnext(d,dmin,dstep,rootd)
      end do
      return
      end

      subroutine calc_pwrexpr_grad(lpt,nterms,dmin,dstep,dmax
     *,sigmab,rho,xin,expr,rmin,ncfirst,fbroad,rbroad)
      use lngamma_mod
      real xin(*),expr(*),dmin,dstep,dmax,dnext,rootd,sigmab,sigma
      pi=3.141592653
c
c maximum value of exponent
c
      expmax=25.0
      ntmax=dmax/sigmab+3
      call fctrln(ntmax)
c
c set up array of coefficients
c
      d=dstep
      rootd=1.0
      do m=1,nterms
c
c calculate r dependence of mth term
c
!            d=m*dstep
!            write(6,*) "calc_pwrexp_grad> ", m,d,dmin,dstep,rootd
c
c save the integer of the first non-zero coefficient
c
            if(d.lt.rmin) ncfirst=m+1
c
c generate value of sigma at this d value
c
            sigma=rootd*fbroad*sigmab
            if(sigma.le.d.and.sigma.gt.0.0.and.m.gt.1) then
                  n=d/sigma-1
            else
                  n=0
            endif
            nref=n+3
            gfac=4.0*pi*rho*sigma*sigma*sigma
            gfacln=alog(gfac)
            sigln=gamln(nref)+gfacln
c
c for the first term we add a factor exp(rmin/sigma)
c This is to ensure low  r repulsive potential develops quickly
c
            if(m.eq.1) sigln=sigln-rmin/sigma
            do i=1,lpt
                  jref=(i-1)*nterms+m
                  r=xin(i)
                  fact=rbroad*r
                  rsigma=r/sigma
                  if(r.gt.0.0) then
                        rsigmaln=alog(rsigma)
                        expon=n*rsigmaln-rsigma-sigln-fact
c
c add on the term to give the gradient of the potential
c
                        fact=real(n)/r-1.0/sigma-rbroad
                  else
                        if(n.eq.0) then
                              expon=-rsigma-sigln-fact
                              fact=-(1.0/sigma+rbroad)
                        else if(n.eq.1) then
                              expon=-rsigma-sigln-fact
                              fact=1.0/sigma
                        else
c
c for n>1, the gradient is identically zero at r=0
c
                              expon=-2.0*expmax
                              fact=0.0
                        endif
                  endif
                  if(expon.gt.-expmax) then
                        expr(jref)=fact*exp(expon)
                  else
                        expr(jref)=0.0
                  endif
            end do
            d=dnext(d,dmin,dstep,rootd)
      end do
      return
      end

      function dnext(d,dfac,dstep,rootd)

      real dnext,d,dstep,dfac,dtest,rootd

      dnext=d+dstep
      dtest=d*dfac
      if(dtest.gt.dnext) then
            dnext=dtest
            rootd=rootd*sqrt(dfac)
      endif
      return
      end

      subroutine ftback(in,rmin,q,sofq,grt)
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      dimension q(*),sofq(*),grt(*)
c
c perform a direct transform on sofq for r values, and use rmin to
c generate a background to sofq. Result is in sofq.
c
      do i=nqmin(in),nqmax(in)
            sofq(i)=sofq(i)*q(i)
      end do
      pi=3.141592653
      qstep=q(4)-q(3)
      sfac=0.5*qstep/pi/pi/rho
      jsave=2
      do j=2,npt1
            rr=x(j)
            if(rr.ge.rmin) go to 200
            jsave=j
            sum=0.
            do i=nqmin(in),nqmax(in)
                  qq=q(i)
                  qr=qq*rr
                  sum=sum+sofq(i)*sin(qr)
            end do
            grt(j)=sfac*sum
      end do
200   continue
c
c calculate the S(Q) from the truncated g(r)
c
      gfac=4.0*pi*rho
      do i=nqmin(in),nqmax(in)
            qq=q(i)
            if(qq.gt.0.0) then
                  r1=0
                  sum=0.0
                  do j=2,jsave
                        rr=x(j)
                        r2=2.0*rr-r1
                        rstep=r2-r1
                        qr=qq*rr
                        sum=sum+grt(j)*sin(qr)*rstep
                        r1=r2
                  end do
                  sofq(i)=gfac*sum/qq
            else
                  r1=0
                  sum=0.0
                  do j=2,jsave
                        rr=x(j)
                        r2=2.0*rr-r1
                        rstep=r2-r1
                        sum=sum+grt(j)*rr*rstep
                        r1=r2
                  end do
                  sofq(i)=gfac*sum
            endif
      end do
      return
      end

      subroutine uniform(ycoh,xst,rho,rmaxgr,ndim)
      include 'lattice.inc'
      integer*8 ntries
c
c calculates the uniform distribution for a unit cell of arbitrary shape
c
      dimension ycoh(ndim)

c now calculate the region out to the edge of the box analytically

      pi=4.0*atan(1.0)
      rhofac=4.0*pi*rho/3.0
      x1=0.0
      x13=0.0
      xstep=xst*0.5
      x2=x1+xstep
      ycoh(1)=0
      j=1
      do while(x2.lt.rmaxgr.and.j.lt.ndim)
         j=j+1
         x2=x1+xstep
         x23=x2*x2*x2
         ycoh(j)=rhofac*(x23-x13)
         xstep=xst
         x1=x2
         x13=x23
      end do

      return
      end

      subroutine invert_cwt(cwt,cwtt,carray,carrayold,feedback
     1,efacm,nfile,npairs,mfile,mpairs)
c
c generate the smart tranpose of cwt to take account of the different weights
c of particular partials in different datasets
c
      real*8 cwt(mpairs,mpairs),cwtt(mpairs,mpairs+mfile)
     1,carray(mpairs+mfile),carrayold(mpairs+mfile,mpairs+mfile)
      real*8 temp,sum,sumold,total,totalold,tdif,add
      real*4 efacm(mfile)
c
c the chosen method is to select the values of the transpose matrix by Monte Carlo,
c calculate the cross-correlation matrix, Carray, and aim towards diagonal elements
c of this array to be unity, while off diagonal units are zero. To develop a single
c solution when the problem is overdetermined, the sum of the squares of the tranpose
c coefficients is kept as small as possible.
c
      niter=50000
      ntimes=20
c
c FEEDBACK is the feedback factor used in the simulation. FEEDBACK must not be > 1.0
c
      feedback1=1.0-feedback
c
c we will assume there is, in an addition to the data, a distinct
c set of partials available (in this case from the simulation itself)
c to make up for the lack of distinct information in the data
c
c this may affect the values of the inversion matrix, but is needed to
c force a unique solution. However since we use the inversion matrix only
c on the DIFFERENCE between data a fit, the inversion matrix for these
c pseudo datasets will have no effect on the values of the empirical potential
c since the difference for the pseudo data is by definition zero.
c
      nnew=npairs+nfile
      temp=1.0/real(nnew)
c
c write out input coefficients and form average modulo(weight)
c
      iadd=0
      write(6,*)
      do i=1,nfile
         do j=1,npairs
            cwt(i,j)=cwt(i,j)*feedback*efacm(i)
         end do
         write(6,200) i,(cwt(i,j),j=1,npairs)
      end do
      nfile1=nfile+1
c
c keep refining until a sensible solution is found
c
      totaltest=1.0e-15
c
c sum the cross array
c
      total=0.0
      do i1=1,nnew
         do i2=1,nnew
            if(i1.le.nfile) then
               sum=0.0
               do j=1,npairs
                  sum=sum+cwt(i1,j)*cwtt(j,i2)
               end do
            else
               j=i1-nfile
               sum=feedback1*cwtt(j,i2)
            endif
            carrayold(i1,i2)=sum
            if(i1.eq.i2) then
               total=total+(sum-1.0)**2
            else
               total=total+(sum)**2
            endif
         end do
      end do
      totalold=total
c
c set starting step size
c
      amp=0.01*total
c
c new fitting temperature
c
      temp=(totalold-nfile)/(nnew*nnew)
      itimes=1
      iter=0
      nreject=0
      do while((totalold-nfile).ge.totaltest.and.itimes.le.1000)
c        write(6,*) itimes,iter,nreject,totalold,temp
         do i2=1,nnew
c
c generate new values
c
            do j=1,npairs
               iter=iter+1
               add=(2.0*ran1()-1)*amp
c
c calculate the cross array for this set of values
c
               do i1=1,nfile
                  carray(i1)=carrayold(i1,i2)+cwt(i1,j)*add
               end do
c
c include the extra term for the diagonal cwts
c
               i1=j+nfile
               carray(i1)=carrayold(i1,i2)+feedback1*add
c
c calculate the change in total for the changes in these specific elements of carray
c
               tdif=0.0
               do i1=1,nfile
                  sum=carray(i1)
                  sumold=carrayold(i1,i2)
                  if(i1.eq.i2) then
                        tdif=tdif+(sum-1.0)**2-(sumold-1.0)**2
                  else
                        tdif=tdif+sum**2-sumold**2
                  end if
               end do
               i1=j+nfile
               sum=carray(i1)
               sumold=carrayold(i1,i2)
               if(i1.eq.i2) then
                  tdif=tdif+(sum-1.0)**2-(sumold-1.0)**2
               else
                  tdif=tdif+sum**2-sumold**2
               end if
               tdif=tdif/temp
c
c perform Monte Carlo test to accept or reject this move
c
               iaccept=0
               if(tdif.le.0.0) then
                  iaccept=1
               else if(tdif.gt.20.0) then
                  iaccept=0
               else if(ran1().le.dexp(-tdif)) then
                  iaccept=1
               else
                  iaccept=0
               endif
               if(iaccept.eq.1) then
c
c move accepted - save the new coefficients
c
                  total=total+tdif*temp
                  totalold=total
                  cwtt(j,i2)=cwtt(j,i2)+add
                  do i1=1,nfile
                     carrayold(i1,i2)=carray(i1)
                  end do
c
c include the extra term for the diagonal cwts
c
                  i1=j+nfile
                  carrayold(i1,i2)=carray(i1)
               else
                  nreject=nreject+1
c           write(6,*) iter,nreject,totalold,tdif,temp
               end if
            end do
         end do
c        write(6,*) itimes,iter,nreject,totalold,temp
c
c modify step size to give 50% rejection rate
c
         if(iter.ge.niter) then
            if(nreject.gt.0) then
               fac=0.5*real(iter-nreject)/real(nreject)
            else
               fac=4.0
            endif
            if(fac.lt.0.95) fac=0.95
            if(fac.gt.1.05) fac=1.05
            amp=amp*fac
c           write(6,*)
c           write(6,*) itimes,total,totalold,totalold-nfile
c           write(6,*) nreject,iter,fac,amp,temp
            itimes=itimes+1
            iter=0
            nreject=0
c
c sum the cross array
c
            total=0.0
            do i1=1,nnew
               do i2=1,nnew
                  if(i1.le.nfile) then
                     sum=0.0
                     do j=1,npairs
                        sum=sum+cwt(i1,j)*cwtt(j,i2)
                     end do
                  else
                     j=i1-nfile
                     sum=feedback1*cwtt(j,i2)
                  endif
                  carrayold(i1,i2)=sum
                  if(i1.eq.i2) then
                     total=total+(sum-1.0)**2
                  else
                     total=total+(sum)**2
                  endif
               end do
            end do
            totalold=total
c
c new fitting temperature
c
            temp=(totalold-nfile)/(nnew*nnew)
         endif
      end do
200   format(i5/(6(1x,e12.5)))
      write(6,*)
c     do j=1,npairs
c        write(6,200) j,(cwtt(j,i),i=1,nnew)
c     end do
c
c calculate the cross array to demonstrate that we have indeed separated
c all the partials - this array should be the unit matrix
c
      do j1=1,npairs
         do j2=1,npairs
            sum=0.0
c
c first sum over the data files
c
            do i=1,nfile
               sum=sum+cwtt(j1,i)*cwt(i,j2)
            end do
c
c add on the pseudo terms
c
            i=j2+nfile
            carray(j2)=sum+cwtt(j1,i)*feedback1
         end do
         write(6,200) j1,(carray(j2),j2=1,npairs)
c
c modify the cwtt for i.le.nfile
c
         do i=1,nfile
            cwtt(j1,i)=cwtt(j1,i)*feedback*efacm(i)
         end do
c
c modify the cwtt(j,i) for i .gt. nfile
c
         do i=nfile+1,nnew
            cwtt(j1,i)=cwtt(j1,i)*feedback1
         end do
      end do
      return
      end

      subroutine generate_rp()

      use threebody_mod
      use containerforces_mod !These are attractive forces designed to prevent particular atom pairs from becoming too far apart.

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'epsrinput_par.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'

      real width,pos,gradlimit,loggradlimit,ddiff

! Limit on gradient of repulsive force
      gradlimit=expecfpwr
      if(gradlimit.gt.0.0) then
         loggradlimit=log(gradlimit)
      else
         gradlimit=0.0
         loggradlimit=0.0
      end if

! Generate molecular potential (based on ecore and dcore)

      call generate_mp()

      pi=3.141592653
      if(power.gt.0.0) ratt=1.0e10**(1.0/power)
c
c maximum value of a sum - this it to prevent overflows
c
      summax=exp(25.0)
      nterms=ncoeffp
      decayp=psigma2

c
c set up L-J, Coulomb, Morse and Gaussian potentials between pairs.
c
c first define prefactor for the well depth
c
      rat=6.0/power
      rat1=sqrt(rat)
      rat2=1.0/rat1
      power1=(power+6)/(power-6)
      rat3=rat1**power1
      prefac=(rat2-rat1)*rat3
      prefac=1.0/prefac
      do ic=1,npairs
            irefp=irefpcof(ic)
            write(6,*) 'generate_rp> ',ic,irefp,rminex(irefp),rmaxex(irefp)
            it=n1extr(ic)
            jt=n2extr(ic)
c
c Coulomb interactions are treated as between uniform spherical charge clouds
c of diameter equal to qradius(it)
c
            ra=qradius(it)
            rb=qradius(jt)
            ect=ecor(ic)
            qct=qcor(ic)*coulfac(irefp)
            dc=dcor(ic)
            x1=x(1)
            x2=x(2)
            do j=2,npt
                  x2=x(j)
c
c modified L-J potential n-6 where n>6
c
                  ddif=dc/x2
                  if(ddif.gt.0.0.and.ect.gt.0.0) then
                    fac=prefac*ect
                    xx=alog(ddif)
                    pwr1=0.5*(power+6.0)*xx
                    if(abs(pwr1).gt.20.0) pwr1=sign(20.0,pwr1)
                    pwr2=0.5*(power-6.0)*xx
                    if(abs(pwr2).gt.40.0) pwr2=sign(40.0,pwr2)
                    sumlj=2.0*fac*sinh(pwr2)*exp(pwr1)
                  else
                    sumlj=0.0
                  endif
                  xc=x2
                  sumc=qct*chsph(ra,rb,xc)

c Apply termination functions

                  sumlj=refpotfac*sumlj*wfac(j)
                  sumc=refpotfac*sumc*wfacc(j)

c Save results

                  epot(j,ic)=sumlj+sumc
                  ecoul(j,ic)=sumlj+sumc
!                  if (coul_intra) xav(j,ic)=xav(j,ic)+coul_intra_fac*sumc

c Add in Gaussians or modified Morse functions (Only the first can be a modified Morse function)

                  ig=0
                  do while (ig.lt.npgauss(irefp))

                     ig=ig+1

                     wpg=abs(wpgauss(ig,irefp))
                     if(wpg.gt.0.0) then
                        rpg=abs(rpgauss(ig,irefp))
                        hpg=hpgauss(ig,irefp)
                        add=0.0
                        if(ig.eq.1.and.addpotential.eq.'modmorse') then
                           width=wpg
                           pos=rpg
                           add=abs(hpg)*gauss6pot(x(j),pos,width)
                        else
                           rdiff=(x(j)-rpg)/wpg
                           term=log(abs(hpg))-0.5*rdiff*rdiff
                           if(abs(term).lt.25.0) add=sign(1.0,hpg)*exp(term)
                        endif
c     write(6,*) j,x(j),rpg,wpg,hpg,rdiff,term,add
!                        epot(j,ic)=epot(j,ic)+add*wfac(j)
                        ecoul(j,ic)=ecoul(j,ic)+add*wfac(j)
                     end if
                  end do
! Try this! Put full reference potential into intramolecular terms
                  if (coul_intra) xav(j,ic)=coul_intra_fac*ecoul(j,ic)
                  x1=x2
            end do
c
c add on the extra repulsive term for the intermolecular potential
c then calculate the potential gradients.
c
c only do something if intensity below rminex(ip) was found
c
            if(abs(pcoefex(irefp)).gt.0.0.and.abs(rminex(irefp)).gt.0.0) then
               expecfuse=expecf*sqrt(abs(rminex(irefp)))
               if(.not.harmonic) then
! Exponential repulsive core
                  rminextest=abs(rminex(irefp))
!                  call calc_pwrexpr(npt2,1,pdstep,pdmax
!     *,expecf,rho,x,expgr,rminextest,ncfirst,1.0,rbroad)
                  do j=1,npt2
                        ddiff=(rminextest-x(j))/expecfuse
                        if(gradlimit.gt.0.0.and.ddiff.ge.loggradlimit) then !Use a linear term to prevent exponential divergences
                           add=pcoefex(irefp)*gradlimit*(1.0+ddiff-loggradlimit)
                        else
                           add=pcoefex(irefp)*exp(ddiff)
                        end if
!                        epot(j,ic)=epot(j,ic)+add*wfac(j)
                        ecoul(j,ic)=ecoul(j,ic)+add*wfac(j)
                  end do
               else
! Gaussian repulsive core
                  rminextest=abs(rminex(irefp))+expecfuse
                  do j=1,npt2
                     ddiff=(rminextest-x(j))/expecfuse
                     if(ddiff.gt.0.0) then
                        if(gradlimit.gt.0.0.and.ddiff.ge.gradlimit) then !Use a linear term to prevent exponential divergences
                           add=0.5*pcoefex(irefp)*gradlimit*(2*ddiff-gradlimit)
                        else
                           add=0.5*pcoefex(irefp)*ddiff*ddiff
                        end if
!                        epot(j,ic)=epot(j,ic)+add*wfac(j)
                        ecoul(j,ic)=ecoul(j,ic)+add*wfac(j)
                     end if
                  end do
               end if
            endif
            epot(1,ic)=epot(2,ic)
            ecoul(1,ic)=ecoul(2,ic)
            xav(1,ic)=xav(2,ic)
c
c set up gradients of these potentials
c
            sumlj1=epot(1,ic)
            sumc1=ecoul(1,ic)
            sumx1=xav(1,ic)
            x1=xstar(1)
            x2=xstar(2)
            do j=2,npt
                  x2=xstar(j)
                  xdif=x2-x1
                  sumlj2=epot(j,ic)
                  sumc2=ecoul(j,ic)
                  sumx2=xav(j,ic)
                  epot1(j-1,ic)=(sumlj2-sumlj1)/xdif
                  ecoul1(j-1,ic)=(sumc2-sumc1)/xdif
                  xav1(j-1,ic)=(sumx2-sumx1)/xdif
                  x1=x2
                  sumlj1=sumlj2
                  sumc1=sumc2
                  sumx1=sumx2
            end do
            epot1(npt,ic)=0.0
            ecoul1(npt,ic)=0.0
            xav(npt,ic)=0.0
            !Get a list of atom pairs for this type. If rminex(irefp).ge.0.0 then the list is empty.
            call get_containerforce_list(nmol,ic,rmaxex(irefp),it,jt,itype,qtype,natom,qtypea,ntype)
      end do
      if(numcontainerforces.gt.0) then
         do ic=1,numcontainerforces
            write(6,*) 'generate_rp> ',ic,pairrefstore(1,ic),pairrefstore(2,ic)
         end do
      end if
      return
      end

      subroutine generate_ep(changepotfac)
c
c generate EP from coefficients
c
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'epsrinput_par.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'

      logical, intent(in) :: changepotfac
      real energabs,pressfac,abs_energy_ep,absediff

! First evaluate the empirical potential with current coefficients
c
c maximum value of a sum - this it to prevent overflows
c
      summax=exp(25.0)
c
c sum up the EP's over the different datafiles
c
      nterms=ncoeffp
      decayp=psigma2
c     write(6,*) nterms,decayp,rbroad,rho,xstar(10),dmax
c
c set up the transform matrix from power-exponentials to r-space
c
      if(gaussian) then
         call calc_gaussr_grad(npt2,nterms,pdstep,pdmax
     *,decayp,rho,x,expgr,0.0,ncfirst,1.0,rbroad)
      else
         call calc_pwrexpr_grad(npt2,nterms,pdmin,pdstep,pdmax
     *,decayp,rho,x,expgr,0.0,ncfirst,1.0,rbroad)
      end if
c
c step through the partials and calculate the empirical potential
c
      do ip=1,npairs
            irefp=irefpcof(ip)
c
c Set up the potential coefficients, ignoring the first one
c
            do ic=2,nterms
                  pwtold(ic)=pcoefp(ic,irefp)
            end do
            pwtold(1)=0.0

c
c generate the potential corresponding to these coefficients
c
            sum1=0.0
            sums1=0.0
            do j=1,npt2
                  term=expand(nterms,j,expgr,pwtold,summax)
                  xavs1(j,ip)=term*pcoefp(1,irefp)
            end do
      end do

c Calculate the potential from the gradients.
c
      do ip=1,npairs

            do j=1,npt2
c
c gradients need to be dimensionless, and multiply by truncation function
c
                  xavs1(j,ip)=xavs1(j,ip)*dmax
            end do
c
c from gradients, calculate empirical potentials, integrating
c from the last point
c
            xavs(npt,ip)=0.0
            xavs1(npt,ip)=0.0
            x1=xstar(npt)
            xts1=xavs1(npt,ip)
            do j=1,npt2
                  jref=npt-j
                  x2=xstar(jref)
                  xts2=xavs1(jref,ip)
                  xavs(jref,ip)=xavs(jref+1,ip)-0.5*(x1-x2)*(xts1+xts2)
                  x1=x2
                  xts1=xts2
            end do

c Multiply empirical potential by truncation function

            do j=1,npt2
!               xavs(j,ip)=xavs(j,ip)
               xavs(j,ip)=xavs(j,ip)*wfac(j)
            end do
!            if(ip.eq.1) write(6,*) 'generate_ep> ',ip,irefp,(xavs(j,ip),j=1,20)
      end do
!Determine the absolute energy of the EP
      energabs=abs_energy_ep()
      pressfac=1.0
      if(abs(absolute_energy).eq.0.0) then
         pressfac=0.0
         energabs=0.0
      else if(abs(energabs).gt.0.0) then
         if(abs(energabs).lt.abs(erequnit)) then
            pressfac=abs(erequnit)/abs(energabs)
         else
            pressfac=abs(absolute_energy)/abs(energabs)
            if(pressfac.gt.1.0.and.ereqstep.eq.0.0) pressfac=1.0
         end if
      endif
      write(6,*) 'generate_ep> ', absolute_energy,energabs,pressfac
! Apply potential factor
      call multiply_pot_coeffs(pressfac)
! Now re-evaluate the derivatives with these truncated potentials, multiplying by supplied factor
      do ip=1,npairs
            irefp=irefpcof(ip)
! Multiply by potential factor if required and create the total potential
            if(npgauss(irefp).eq.0) then
               do j=1,npt2
                  xavs(j,ip)=xavs(j,ip)*pressfac
                  xavtot(j,ip)=ecoul(j,ip)+xavs(j,ip)
               end do
            else
               do j=1,npt2
                  xavs(j,ip)=xavs(j,ip)*pressfac
                  xavtot(j,ip)=ecoul(j,ip)+xavs(j,ip)
               end do
            end if
            sums1=xavs(1,ip)
            sumtot1=xavtot(1,ip)
            x1=xstar(1)
            x2=xstar(2)
            do j=2,npt
                  x2=xstar(j)
                  xdif=x2-x1
                  sums2=xavs(j,ip)
                  sumtot2=xavtot(j,ip)
                  xavs1(j-1,ip)=(sums2-sums1)/xdif
                  xavtot1(j-1,ip)=(sumtot2-sumtot1)/xdif
                  x1=x2
                  sums1=sums2
                  sumtot1=sumtot2
            end do
            xavs1(npt,ip)=0.0
!            if(ip.eq.1) write(6,*) 'generate_ep> ',ip,irefp,(xavs1(j,ip),j=1,20)
            xavtot1(npt,ip)=0.0
      end do
!Prove that you have done it!
!      energabs=abs_energy_ep()
!      write(6,*) 'generate_ep> ', absolute_energy,energabs
!! Aim to adjust the size of potfac so that changes of EP absolute energy are no larger than potfacamp*energabs
      if(changepotfac.and.absolute_energy.gt.0.0) then !CAUTION! This doesn't work with absolute_energy - need the absolute energy of the previous potential BEFORE it was perturbed.
         if(potfac.gt.0.000001.and.potfacamp.gt.0.0) then
            absediff=(absolute_energy-energabs)/abs(potfacamp*absolute_energy)
            write(6,*) 'changepotfac> ',potfacamp,potfac,absolute_energy,energabs,absediff
            if(abs(absediff).gt.0.1) absediff=sign(0.1,absediff)
            potfac=potfac*(1.0+absediff)
            potfac=max(potfac,0.000001) !Smallest allowed value of potfac
            potfac=min(potfac,1.0) !Maximum value of potfac
         end if
      end if
      return
      end

      subroutine lsparabola(nchan,fname,nsum,navmax,epcfrac,a,b,c,rfacwrite,aav,bav,lastqf,nextqf)

      character*256 fname
      real emin,femin,rfacmin
      real*8 x,y,sum,sumx1,sumx2,sumx3,sumx4,sumy,sumyx1,sumyx2
     *,delta,cplus,cminus,bplus,bminus,pf1,pf2

      integer, intent(in) :: navmax
      integer, intent(out) :: nsum
      real, intent(in) :: epcfrac
      real :: dump,aread,bread,cread,thisqf
      real, intent(out) :: a,b,c,rfacwrite,aav,bav,lastqf,nextqf

      sum=0
      sumx1=0
      sumx2=0
      sumx3=0
      sumx4=0
      sumy=0
      sumyx1=0
      sumyx2=0
      ierr=0
      xmax=0
      nav=0
      nsum=0
      aav=0.0
      bav=0.0
      lastqf=0.0
      nextqf=0.0
      do while (ierr.eq.0)
         read(nchan,*,iostat=ierr) dump,dump,rfacwrite,x,dump,y,dump,aread,bread,cread !This line reads the coefficient used to fit the quadratic.
         if(ierr.eq.0) then
            nsum=nsum+1
            y=epcfrac*y+(1.0-epcfrac)*rfacwrite
            thisqf=y
            if(nsum.eq.1) nextqf=thisqf
            if(x.gt.xmax) xmax=x
            if(nsum.gt.1) then !Ignore first point!
               lastqf=nextqf
               nextqf=thisqf
               nav=nav+1
               if(nav.gt.navmax) nav=navmax

c form running averages so that results are not overly biassed by earlier values

               pf1=1.0/real(nav)
               pf2=real(nav-1)/real(nav)
               sum=pf2*sum+pf1
               sumy=pf2*sumy+pf1*y
               sumx1=pf2*sumx1+pf1*x
               sumyx1=pf2*sumyx1+pf1*y*x
               xm=x*x
               sumx2=pf2*sumx2+pf1*xm
               sumyx2=pf2*sumyx2+pf1*xm*y
               xm=xm*x
               sumx3=pf2*sumx3+pf1*xm
               xm=xm*x
               sumx4=pf2*sumx4+pf1*xm
               aav=pf2*aav+pf1*aread
               bav=pf2*bav+pf1*abs(bread)
            end if
         endif
      end do

c Backspace the file so that it is before the end of file marker

      backspace(nchan)

c Now work out the coefficients the hard way! The method used here assumes all the sums are
c positive definite

c a) determinant of coefficients

      delta=sum*(sumx2*sumx4-sumx3*sumx3)
     *+sumx1*(sumx3*sumx2-sumx1*sumx4)
     *+sumx2*(sumx1*sumx3-sumx2*sumx2)
      delta=0 !Inhibit quadratic fitting

      if(abs(delta).gt.1.0e-10.and.nsum.ge.navmax) then

c b) quadratic coefficient

         cplus=sumy*sumx1*sumx3+sumyx1*sumx2*sumx1+sumyx2*sum*sumx2
         cminus=sumy*sumx2*sumx2+sumyx1*sumx3*sum+sumyx2*sumx1*sumx1

         if(cplus.ne.cminus) then

            c=(cplus-cminus)/delta

         else

            c=0

         endif

c calculate the other coefficients - a straight line fit if c is zero

         sumy=sumy-c*sumx2
         sumyx1=sumyx1-c*sumx3
         delta=(sumx1*sumx1-sum*sumx2)
         if(abs(delta).gt.1.0e-10) then

            bplus=sumy*sumx1
            bminus=sumyx1*sum
            if(bplus.ne.bminus) then

               b=(bplus-bminus)/delta

            else

               b=0

            endif

         else

            b=0

         endif

         sumy=sumy-b*sumx1

         if(abs(sum).gt.0.0) then

            a=sumy/sum

         else

            a=0

         endif

      else if(nsum.ge.2) then

         c=0

c Try fitting a straight line instead

         delta=(sumx1*sumx1-sum*sumx2)
         if(abs(delta).gt.1.0e-10) then

            bplus=sumy*sumx1
            bminus=sumyx1*sum
            if(bplus.ne.bminus) then

               b=(bplus-bminus)/delta

            else

               b=0

            endif

         else

            b=0

         endif

         sumy=sumy-b*sumx1

         if(abs(sum).gt.0.0) then

            a=sumy/sum

         else

            a=0

         endif

      else

         c=0
         b=0
         delta=sum

         if(abs(delta).gt.0.0) then

            a=sumy/delta

         else

            a=0

         endif

      endif
      write(6,*) 'lsparabola> ',nsum,delta,a,b,c

c write a little file to show the fit

      open(8,file=fname,status='unknown')

      npt=100
      xstep=xmax/npt

      !Calculate position and depth/height of minimum/maximum
      if(c.ne.0.0) then
         emin=-b/(2.0*c)
         femin=(a+emin*(b+emin*c))
      else
         emin=0.0
         femin=0.0
      end if
      do ix=1,npt

         x=ix*xstep

c format here has to match that of the .erg file so that it can be plotted with
c the plot command.

         write(8,100) a,b,c,x,(a+x*(b+x*c)),emin,femin

100   format(7(1x,e13.6))

      end do

      close(8)

      return
      end

      subroutine lsparabola1(nchan,fname,nav,navmax,epcfrac,a,b,c,rfacwrite,aav,bav,lastqf,nextqf)
! In this version instead of using running averages we average only over the last navmax points
      real, dimension(:), allocatable :: xval,yval,aval,bval
      character*256 fname
      integer, intent(in) :: navmax
      integer, intent(out) :: nav
      real, intent(in) :: epcfrac
      real, intent(out) :: a,b,c,rfacwrite,aav,bav,lastqf,nextqf
      real :: x,y,emin,femin,rfacmin,dump,aread,bread,cread,fractemp,thisqf
      real*8 :: xm,sum,sumx1,sumx2,sumx3,sumx4,sumy,sumyx1,sumyx2
     *,delta,cplus,cminus,bplus,bminus
      integer nsum,ierr,i

      allocate(xval(navmax),yval(navmax),aval(navmax),bval(navmax))
      ierr=0
      xmax=0
      nav=0
      nsum=0
      lastqf=0.0
      nextqf=0.0
      do while (ierr.eq.0)
         read(nchan,*,iostat=ierr) dump,dump,rfacwrite,x,dump,y,dump,aread,bread,cread !This line reads the coefficient used to fit the quadratic.
!         fractemp=rfacwrite/(rfacwrite+y)
!         y=epcfrac*fractemp*y+(1.0-epcfrac)*(1.0-fractemp)*rfacwrite
         if(ierr.eq.0) then
            nsum=nsum+1
            y=epcfrac*y+(1.0-epcfrac)*rfacwrite
            thisqf=y
            if(nsum.eq.1) nextqf=thisqf
            if(x.gt.xmax) xmax=x
            if(nsum.gt.1) then !Ignore first point!
               lastqf=nextqf
               nextqf=thisqf
               nav=nav+1
               navref=nav-navmax*int((nav-1)/navmax)
               xval(navref)=x
               yval(navref)=y
               aval(navref)=aread
               bval(navref)=bread
            end if
         end if
      end do
      write(6,'(a,5(1x,i5))') 'lsparabola1> ',nsum,nav,navref
! Backspace the file so that it is before the end of file marker
      backspace(nchan)
      nsum=min(nav,navmax)
      sum=0
      sumx1=0
      sumx2=0
      sumx3=0
      sumx4=0
      sumy=0
      sumyx1=0
      sumyx2=0
      aav=0.0
      bav=0.0
! Form averages
      if(nsum.gt.0) then
         do i=1,nsum
            x=xval(i)
            y=yval(i)
            sum=sum+1
            sumy=sumy+y
            sumx1=sumx1+x
            sumyx1=sumyx1+y*x
            xm=x*x
            sumx2=sumx2+xm
            sumyx2=sumyx2+xm*y
            xm=xm*x
            sumx3=sumx3+xm
            xm=xm*x
            sumx4=sumx4+xm
            aav=aav+aval(i)
            bav=bav+abs(bval(i))
         end do
         aav=aav/real(nsum)
         bav=bav/real(nsum)
      end if
      if(allocated(xval)) deallocate(xval,yval,aval,bval)
c Now work out the coefficients the hard way! The method used here assumes all the sums are
c positive definite

c a) determinant of coefficients

      delta=sum*(sumx2*sumx4-sumx3*sumx3)
     *+sumx1*(sumx3*sumx2-sumx1*sumx4)
     *+sumx2*(sumx1*sumx3-sumx2*sumx2)
      delta=0 !Inhibit quadratic fitting

      if(abs(delta).gt.1.0e-10.and.nsum.ge.navmax) then

c b) quadratic coefficient

         cplus=sumy*sumx1*sumx3+sumyx1*sumx2*sumx1+sumyx2*sum*sumx2
         cminus=sumy*sumx2*sumx2+sumyx1*sumx3*sum+sumyx2*sumx1*sumx1

         if(cplus.ne.cminus) then

            c=(cplus-cminus)/delta

         else

            c=0

         endif

c calculate the other coefficients - a straight line fit if c is zero

         sumy=sumy-c*sumx2
         sumyx1=sumyx1-c*sumx3
         delta=(sumx1*sumx1-sum*sumx2)
         if(abs(delta).gt.1.0e-10) then

            bplus=sumy*sumx1
            bminus=sumyx1*sum
            if(bplus.ne.bminus) then

               b=(bplus-bminus)/delta

            else

               b=0

            endif

         else

            b=0

         endif

         sumy=sumy-b*sumx1

         if(abs(sum).gt.0.0) then

            a=sumy/sum

         else

            a=0

         endif

      else if(nav.ge.2) then

         c=0

c Try fitting a straight line instead

         delta=(sumx1*sumx1-sum*sumx2)
         if(abs(delta).gt.1.0e-10) then

            bplus=sumy*sumx1
            bminus=sumyx1*sum
            if(bplus.ne.bminus) then

               b=(bplus-bminus)/delta

            else

               b=0

            endif

         else

            b=0

         endif

         sumy=sumy-b*sumx1

         if(abs(sum).gt.0.0) then

            a=sumy/sum

         else

            a=0

         endif

      else

         c=0
         b=0
         delta=sum

         if(abs(delta).gt.0.0) then

            a=sumy/delta

         else

            a=0

         endif

      endif
      write(6,*) 'lsparabola1> ',nsum,delta,a,b,c

c write a little file to show the fit

      open(8,file=fname,status='unknown')

      npt=100
      xstep=xmax/npt

      !Calculate position and depth/height of minimum/maximum
      if(c.ne.0.0) then
         emin=-b/(2.0*c)
         femin=(a+emin*(b+emin*c))
      else
         emin=0.0
         femin=0.0
      end if
      do ix=1,npt

         x=ix*xstep

c format here has to match that of the .erg file so that it can be plotted with
c the plot command.

         write(8,100) a,b,c,x,(a+x*(b+x*c)),emin,femin

100   format(7(1x,e13.6))

      end do

      close(8)

      return
      end

      subroutine multiply_pot_coeffs(epfact)
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'epsrinput_par.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'

      real epfact

      nterms=ncoeffp
c
c Step through the defined partials
c
      do ip=1,npairs
         irefp=irefpcof(ip)
c The amplitude of the present EP is multiplied by the specified factor, but only
c if the number of Gaussians defined is zero
         do ic=2,nterms
            pcoefp(ic,irefp)=pcoefp(ic,irefp)*epfact
         end do
      end do

      return
      end

      subroutine multiply_pcoefex(epfact)
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'epsrinput_par.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'

      real epfact

c
c Step through the defined partials
c
      do ip=1,npairs
            irefp=irefpcof(ip)
c The amplitude of the present repulsive potential is multiplied by the specified factor
c However the minimum amplitude is limited to 0.0001. This is to prevent the values
c becoming absurdly small. If this still gives too much pressure then reduce the size
c of rminex.
            if(pcoefex(irefp).ne.0.0) then
               if((abs(pcoefex(irefp)).gt.0.0001.or.epfact.gt.1.0).and.rminex(irefp).gt.0.0) then
!                  write(6,*) irefp,rminex(irefp),epfact
                  rminex(irefp)=rminex(irefp)*epfact
                  write(6,*) 'multiply_pcoefex> ',irefp,rminex(irefp),epfact
! Write the new value of this symbol, keeping the other values intact
                  icount=nsymboluptodata+nfiles*nsymbolplus+irefp
                  sometext=symbol_value(icount)
                  call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
                  if(nwords.gt.1) then
                     write(symbol_value(icount),*,iostat=ierr) rminex(irefp),sometext(ncf(2):ncl(nwords))
                  else
                     write(symbol_value(icount),*,iostat=ierr) rminex(irefp)
                  end if
               endif
            end if
      end do
      call generate_rp()

      return
      end



