      program epsr_standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call runepsr17()

      stop
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

      real z(mpz),Q(mpts),datai(mpts,mfile),errs(mpts,mfile)
     1,qwind(2,mpts,mpairs),rfac(mpairs)
     1,datat(mpts,mfile),backsc(mpts,mfile),wav(2,mpz,mpairs)
     5,errsav(mpts,mpairs),soqism(mpts,mfile)
     1,rshmin(mfile),szeros(mfile)
     1,rfacd(mfile),potnoi(mfile)
     1,energf(mfile),pressf(mfile),efac(mfile),tweak(mfile)
     1,datap(mpts,mpairs),datagr(mpz,mpairs),datagr2(mpz,mpairs)
     1,dataigr(mpz,mfile)
     1,grlim(mpairs),grzerolim(mfile),efilereq(mfile)

      real QD(mqd),DATA(mqd),err(mqd),gofzism(mpz,mfile)
     *,gz2ism(mpz,mfile)

c Temporary store of simulated intermolecular partial structure factors

      real soqintertemp(mpts,mpairs),dts(mpts),dti(mpts)

c Factor used to limit the absolute energy of the empirical potential

      real elimitfac

c nrtype is the type of dataset being read - see below.
c nqminp and nqmaxp are used in the Fourier transform of the estimated partial structure
c factors.

      integer nrtype(mfile),nqminp(mpairs),nqmaxp(mpairs)

      character*256 infile,fname,infili(mfile)
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

c
c EPSR begins....
c

      PI=3.141592653

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
      end do
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
      read(symbol_value(iref),*) potfac
      iref=iref+1
      read(symbol_value(iref),*) pressreq
      iref=iref+1
      read(symbol_value(iref),*) sizefactor
      if(feedback.gt.1.0) feedback=0.9
      if(feedback.lt.0.0) feedback=0.0
c do not attempt to refine potential if there is no weight given to the data
      if(feedback.le.0.0) potfac=0.0
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

      if(nonzero.gt.0) sumnonzero=sumnonzero/nonzero
      do in=1,nfile
         efacm(in)=efilereq(in)/sumnonzero
      end do

c NQ = No. of Q values (Maximum = 1000)
c QSTEP = requested step in  Q
      iref=iref+1
      read(symbol_value(iref),*) nq
      iref=iref+1
      read(symbol_value(iref),*) qstep
      if(nq.gt.mpts) nq=mpts
c NEGAT not currently operational - hold at zero
      negat=0
c IRESET = 1 to set EP to zero and reinitialize certain accummulators
c     else IRESET = 0.
      iref=iref+1
      read(symbol_value(iref),*) ireset
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
c molecular move (molecular vibration) is made - typically INTRA=50 or 100
c INTER is the number of iterations the potential is averaged over
      iref=iref+1
      read(symbol_value(iref),*) intra
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
c read minimum Q value for molecular fit
      iref=iref+3
      read(symbol_value(iref),*) qminm
      facpress=0.001
      facpress3=0.01
c
c define pressure and energy factors not currently used, but needed.
c
      eljfac=1.0
      pljfac=1.0
      piFAC=1./SQRT(2.*PI)
c
c normalise atom and molecular coordinates to box size
c
      call normalise_coordinates
      cell=cellmax
      datafc=4.*pi*rho
c
c write the repulsive potential coefficients
c
c     do ip=1,npairs
c           write(6,*) n1extr(ip),n2extr(ip),dcor(ip)
c     end do
c
c define the q scale on which calculation is made 
c
      nqminm=1
      do 5 i=1,nq
      Q(I)=qstep*(i-1)
      if(q(i).le.qminm) nqminm=i
5     continue
c
c initialize distributions
c
      write(6,*) rmaxcell,cellmax
      call inita(ireset,iinit,nq,q,nz,z,fwhm,fwhmq,intra)
      write(6,*) rminpt,rmaxpt
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

      call initpt(ireset,rminpt,rmaxpt,polar,szeros
     1,rshmin,rbroad,expecf)

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

c multiply data by density factor so that all data are per atom

            do i=1,nqd
                  data(i)=data(i)*rhofac(in)*tweak(in)
            end do

            write(6,352) nqd,infile(1:ilentot)
352   format(' Read ',i5,' points from file ',a)

c interpolate input data onto same scale as calculation

            call interp(nqd,nq,qd,q,data,datai,in,nqmin(in),nqmax(in)
     1,mqd,mpts,mfile)  
            call interp(nqd,nq,qd,q,err,errs,in,nqmin(in),nqmax(in)
     1,mqd,mpts,mfile)

c prevent nqmax from being .ge. nq (This is for smoothing purposes in PTGNSQ)

            if(nqmax(in).ge.nq) nqmax(in)=nq-1
            if(nqmin(in).lt.nqminm) nqmin(in)=nqminm
            nqpts(in)=nqmax(in)-nqmin(in)+1
            write(6,958) nqmin(in),nqmax(in),nqpts(in),nq
958   format(1x,5(1x,i5))
      end do
c
c open file to save energy, pressure and R-factor values
c
      fname=fnamerun(1:ilenfnamerun)//'.erg'
      write(6,912) fname(1:ilenfnamerun+4)
912   format((1x,a))
      call openf(9,fname)
c
c skip to the end of the file unless iinit.eq.1
c
      if(iinit.ne.1) then

c calculate the coefficients of the least squares parabola which goes through the R-factor
c and absolute energies

         fname=fnamerun(1:ilenfnamerun)//'.qdr'
         call lsparabola(9,fname,nparab,aparab,bparab,cparab)

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
                  nmovei(j)=0
                  ninter(j)=0
                  do imoltype=1,nmoltype
                     nwhole(j,imoltype)=0
                     nmovem(j,imoltype)=0
                  end do
            end do
      endif
      idone=0
c
c set up theoretical R-factor and calculate renormalization on error bars
c
      theorf=1.0
c
c calculate initial energy of distribution and make neighbour list
c     write(6,*) 'Got to fmolec1'
      call fmolec1(foldsv,foldsb,ierr)
      write(6,1999) foldsv,foldsb
1999  format(1x,2(e13.6,1x))

c exit if an error occurred in fmolec1

      if(ierr.ne.0) return
c     write(6,*) nmol,foldsv*nmol,foldsb*nmol

      do 2000 itimes=1,ntimes
c     write(6,*) 'Got to molfit'

c Ensure molecules stay within box.

         do im=1,nmol
            call minimagec(xadd(im),yadd(im),zadd(im))
         end do

         call molfit(foldsv,foldsb,nisum,intra,ireset)
         write(6,1999) foldsv,foldsb
c     write(6,*) nmol,foldsv*nmol,foldsb*nmol
2000  continue
      write(6,*) 'Got to atogen'
      call atogen(wav)
      write(6,*) 'Got to entropc'
      call entrpc(oiseb,wav)
      write(6,*) 'Got to gensq'
      call gensq(nq,q,qwind,wav)
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
      write(6,*) 'Got to press'
      call press(energd,pressd,energf,pressf,energabs)
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

c Apply an overall factor to the EP

c Energy factor needed to keep absolute energy at or below requested energy

      if(abs(pressreq).eq.0.0) then
         pressfac=0.0
      else if(abs(energabs).gt.abs(pressreq)) then
         pressfac=elimitfac*(abs(pressreq)-abs(energabs))
         pressfac=pressfac*pressfac
         if(pressfac.lt.25.0) then
            pressfac=exp(-pressfac)
         else
            pressfac=0.0
         endif
      else
         pressfac=1.0
      endif

      do in=1,nfile
         efac(in)=pressfac
      end do

      write(6,*) pressfac
      write(6,*) (efac(in),in=1,nfile),(efacm(in),in=1,nfile)

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
               if(iq.ge.nqmin(in).and.iq.le.nqmax(in)) then
                  sum=0.0
                  do ip=1,npairs
                     sum=sum+weightx(iq,ip,ix)*soqintra(iq,ip)
                  end do
                  datat(iq,in)=datai(iq,in)-sum
               else
                  datat(iq,in)=0.0
               endif

            end do

         else

            do iq=1,nq

               if(iq.ge.nqmin(in).and.iq.le.nqmax(in)) then
c
c sum the intramolecular terms and subtract them from the diffraction data
c
                  sum=0.0
                  do ic=1,npairs
                     sum=sum+ciwt(in,ic)*soqintra(iq,ic)
                  end do
                  datat(iq,in)=datai(iq,in)-sum
               else
                  datat(iq,in)=0.0
               endif

            end do
         endif
      end do

c Set up difference data - simulation and subtract a background from data if needed

      do id=1,nfile

c Form the simulated fit for this dataset

c If X-ray data, need to apply the Q-dependent weights

            ix=xfile(id)

c Set up the gradient for extrapolating to Q=0

            if(szeros(id).ne.0.0.and.nqmin(id).gt.1) then

               grads0=(datat(nqmin(id),id)-szeros(id))/(q(nqmin(id))-q(1))
               conss0=szeros(id)

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
                  dts(i)=conss0+grads0*(q(i)-q(1))-sum
               else if(i.ge.nqmin(id).and.i.le.nqmax(id)) then
                  dts(i)=datat(i,id)-sum
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
c     *,potfac,rminpt,rmaxpt,szeros,overallfac,efac,efacc,gfac
c     *,rminp,potnoi,wtentp
c     *,rbroad,expecf,nsmoop)
      write(6,*) 'Got to ptgnsq'
      call ptgnsq(nesum,nq,q,wav,datat,soqintertemp
     *,potfac,rminpt,rmaxpt,szeros,overallfac,efac,efacc,gfac
     *,rshmin,potnoi,wtentt
     *,rbroad,expecf,nsmoop)
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
c perform direct Fourier transform on inverted data ready for output 
c
            call stog(rho,q,datap,z,datagr,nqminp,nqmaxp,nz,npairs
     *,mpts,mpz,mpairs,grlim)
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
            lpt1=nqmin(1)
            lpt2=nqmax(1)
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
c sum the intramolecular terms and add them from the interference diffraction data
c
                  sum=0.0
                  if(i.ge.nqmin(in).and.i.le.nqmax(in)) then
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

                  soqism(i,in)=sum
                  if(i.ge.nqmin(in).and.i.le.nqmax(in)) then
                     backsc(i,in)=datai(i,in)-soqism(i,in)
                     diffraw=datai(i,in)-sumraw
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
c write out current energy and pressure and Rfactor and density
c
            write(9,594) energd,pressd,rfacraw,energabs
c           write(6,594) energd,pressd,rfactot,rhoout
594   format(6(1x,e13.6))
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

            call stog(rho,q,soqism,z,gofzism,nqmin,nqmax,nz,nfile
     *,mpts,mpz,mfile,grzerolim)
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
c write out the empirical potentials
c
            write(ext1,331) '.p'
            fname=fnamerun(1:ilenfnamerun)//ext1(1:2)
            call w_int(nz,fname,npairs,z,xavs,xav
     *,titler,titlepsf,mpz,mpairs)

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
            if(iinit.ne.0) iinit=0
      endif
c
c output summary of results
c
      extfile='.out'
      fname=fnamerun(1:ilenfnamerun)//extfile(1:4)
      call openf(10,fname)
      wtent=0.0
      write(10,146) version
146   format(/2x,'EPSR version ',a,/)
      write(10,147) itimes,iter,nsum,rfactot
     1,foldsv,foldsb
     1,energabs
      write(10,148) (j,rfacd(j),energf(j),pressf(j),j=1,nfile)
      write(10,149) (j,rfac(j),j=1,npairs)
147   format(1x,' Trial no. ',i5,5x,' Iteration no. ',i5
     1,5x,' No. of trials in sum ',i5
     2,/1x,' Chi-square ',e13.6
     4,/1x,' Energies: repulsive ',e13.6,' harmonic ',e13.6
     4,/1x,' Abs. empirical energy: ',1(1x,e13.6))
148   format((1x,' R-factor, abs. energy and pressure for data file ',i2,':',3(1x,e11.4)))
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
         nwholet(j)=0
         nmovemt(j)=0
         do imoltype=1,nmoltype
            nwholet(j)=nwholet(j)+nwhole(j,imoltype)
            nmovemt(j)=nmovemt(j)+nmovem(j,imoltype)
         end do
      end do
      write(10,220) (nmovei(j),ninter(j),nwholet(j),nmovemt(j)
     *,j=1,2)
220   format(1x,' Number of rejects in MOLFIT:-'/((1x,4(i7,1x))))
278   format(i5/(6(1x,e12.5)))
      close (10)

c Estimate the next value of ereq based on the parabolic fit to the R-factor versus the
c absolute energy. Only do this if ereq is negative and the absolute energy is within 
c 5 kJ/mole or 10% of the requested absolute energy.

      test=abs(energabs-abs(pressreq))

      if(pressreq.lt.0.0.and.(test.lt.0.05*abs(pressreq).or.test.lt.1.0)
     *.and.nparab.gt.20) then

c Write down expected value of ereq to give minimum R-factor and minimum R-factor.
c This can only be done if a definite minimum has been identified

         if(bparab.lt.0.0.and.cparab.gt.0.0) then
            ereqmin=-0.5*bparab/cparab
            if(ereqmin.gt.(energabs+50.0)) ereqmin=energabs+50.0
            rfmin=aparab+ereqmin*(bparab+ereqmin*cparab)

c Set a range of values of ereq to give a 20% fluctuation in R-factor about the minimum - 
c this is to ensure a range of values are sampled. Protect yourself in case fit is almost 
c linear and causes the predicted minimum R-factor to be negative.

            if(rfmin.gt.0.0) then
               erange=sqrt(rfmin*0.2/cparab)
            else
               erange=0.0
            endif

c Throw a dice to specify the next value of ereq within this range

            step=ereqmin+(2.0*ran1()-1)*erange-abs(pressreq)
      write(6,*) ereqmin,rfmin,erange,step

c Prevent the change in EREQ from being too large in the event that EREQ is far
c from the minimum.

            if(abs(step).gt.10.0) step=sign(10.0,step)

         else

c If a minimum in R-factor cannot be found then base the direction of the step on the 
c gradient of the quadratic at the current value of the absolute energy.

            grad=bparab+2.0*cparab*energabs
            step=-ran1()*sign(1.0,grad)

         endif

c "Bounce" the value of EREQ at the value of 5 to prevent it becoming too small or negative

         test=abs(pressreq)+step
         if(test.lt.1.0) then

            step=2.0-2.0*abs(pressreq)-step

         endif

         pressreq=pressreq-step           

      end if
c
c write input file
c
c
c first save the internal values to their symbol counterparts
c

      text='ereq'
      icount=0
      call search_for_symbol(text,icount)
      iref=ireff(icount,1)
      write(symbol_value(iref),*) pressreq

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
      call save_files_others(14)

c     write(6,*) 'Done save_files_others'

      call normalise_coordinates()

c     write(6,*) 'Done normalise'


3000  continue
      close(9)
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
      subroutine inita(ireset,iinit,nq,q,nz,z,fwhm,fwhmq,intra)
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'epsr16.inc'
      include 'epsrinput_par.inc'
      include 'molecom.inc'
      include 'potential_par.inc'
      include 'weights.inc'
      include 'iresetnumbers.inc'
      real q(mpts),z(mpz)
      character*256 fname
      character*10 extfile
      pi=3.141592653
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
c RCHARGE = 1 to calculate energy due to polarisation of molecule. Currently this only works
c correctly for SPC/E water. Normally set to 0.
      iref=ireff(iqminref,1)
      iref=iref+1
      read(symbol_value(iref),*) ROVERLAP
      iref=iref+1
      read(symbol_value(iref),*) RMINFAC
      iref=iref+1
      read(symbol_value(iref),*) RMINPT
      iref=iref+1
      read(symbol_value(iref),*) RMAXPT
      iref=iref+1
      read(symbol_value(iref),*) RBROAD
      iref=iref+1
      read(symbol_value(iref),*) EXPECF
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
      read(symbol_value(iref),*) ncoeffp
      iref=iref+1
      read(symbol_value(iref),*) pdmax
      iref=iref+1
      read(symbol_value(iref),*) pdstep
      iref=iref+1
      read(symbol_value(iref),*) npitss
      iref=iref+1
      read(symbol_value(iref),*) paccept
      iref=iref+1
      read(symbol_value(iref),*) psigma2

c pdmax is set by programme to the truncation in the potential

      pdmax=rmaxpt

      pdstep=pdmax/ncoeffp

c Width of Poisson functions is set by program

      psigma2=0.01*pdstep

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
               nisum=0
               nesum=5
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
            nisum=0
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
c set maximum radius for calculation of g(r)
c
      rmaxgr=amin1(2.0*rmaxpt,rmaxcell)
c
c rmaxgr should be at least as big as half the box size
c
      if(rmaxgr.lt.0.5*cell) rmaxgr=0.5*cell
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
      j=2
      npt1=2
      npt2=2
      nptmin=2
      nptmax=2
      x1d=x1dold+celstt
      do while(j.lt.mpz.and.x1dold.lt.rmaxgr)
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
         if(x1d.lt.0.5*cell) npt1=j+1
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
      npt2=j-1
      if(npt1.gt.npt2) npt1=npt2
      write(6,959) nptmin,nptmax,npt1,npt2,npt
959   format(5(1x,i5))
      decay=x(nptmax)-x(nptmin)
      if(nesumt.lt.1) nesumt=1
      if(nesum.gt.nesumt) nesum=nesumt
      if(nisum.gt.intra) nisum=intra
c
c calculate uniform distribution for the box
c To do this try to ensure that at least one element in the nearest shell (27
c boxes within 1 step from origin, including one at origin) is
c sampled on average. We repeat this calculation every intra MC moves so that
c it does not unduly bias the summed g(r).
c
      if(ireset.ne.0.or.nisum.eq.intra) then
c     if(ireset.ne.0) then
            ycoh(1)=0.0
            atry=5.0/27.0
            call uniform(atry,x,ycoh,cellst,rho,npt2,mpz)
      else
            atry=0.1/27.0
c 
c Read the uniform distribution from a file called ****.uni
c
            extfile='.uni'
            fname=fnamerun(1:ilenfnamerun)//extfile(1:4)
            call openf(10,fname)
            read(10,*,iostat=isave) (ycoh(j),j=1,npt)
            close(10)
c
c if there was an error encountered (most likely because npt2 has changed)
c we regenerate the uniform distribution
c
            if(isave.ne.0) then
               ycoh(1)=0.0
               atry=5.0/27.0
               call uniform(atry,x,ycoh,cellst,rho,npt2,mpz)
            endif
      endif
c
c sigma for instrumental resolution
c
      sigres=0.5*fwhm/sqrt(2.*alog(2.0))
      sigrsq=0.5*fwhmq/sqrt(2.*alog(2.0))
c
c calculate transform matrices
c
      pi2=2*pi
      do i=1,nq
            qq=q(i)
            DO J=1,npt
                  rr=x(j)
                  qr=qq*rr
                  sinqr=sin(qr)
                  qsig=sigres+qq*sigrsq
                  qrs=qsig*rr
                  qr2=0.5*qrs*qrs
                  if(qr2.lt.30.) then
                        if(qr.gt.0.) then
                              sums=sinqr*exp(-qr2)/qr
                        else
                              sums=exp(-qr2)
                        endif
                  else
                        sums=0.
                  endif
                  foqols(j,i)=sums
            end do
      end do
c
c set z scale identical to x scale - inherited from earlier version
c
      nz=npt2
      do i=1,nz
            z(i)=x(i)
      end do
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
      real qarray(mpts),rarray(mpz)
      integer ierr
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
c
c set up weightings on potential, based on values of gorinter
c
      do id=1,npairs
            entsig(1,id)=0.0
            do j=2,nz
                  entsig(j,id)=entsig(j-1,id)+ycoh(j)*gorinter(j,id)
            end do
      end do
      ycohsum=0.0
      do j=2,nz
            ycohsum=ycohsum+ycoh(j)
            if(ycohsum.gt.0.0) then
                  do id=1,npairs
                        entsig(j,id)=entsig(j,id)/ycohsum
                  end do
            endif
      end do
      return
      end
c
c calculate entropy of distribution
c
      subroutine entrpa(n1,n2,dt,entrop,mpts)
      real dt(mpts)
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
      real wav(2,mpz,mfile)

      real*8 fc1,fc2,add1(2),add2(2),sumtot1(2),sumtot2(2),errort(2)
      common/sumvalues/fc1,fc2,add1,add2,sumtot1,sumtot2,errort

      do id=1,npairs
            entsig(1,id)=0.0
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
                  entsig(j,id)=entsig(j-1,id)+ycoh(j)*gorinter(j,id)
            end do
      end do
      ycohsum=0.0
      do j=2,nz
            ycohsum=ycohsum+ycoh(j)
            if(ycohsum.gt.0.0) then
                  do id=1,npairs
                        entsig(j,id)=entsig(j,id)/ycohsum
                  end do
            endif
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
      real qwind(2,mpts,mfile)

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
      real z(mpz)
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

      subroutine initpt(ireset,rminpt,rmaxpt,polar,smoos
     *,rminp,rbroad,sigma0)
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      include 'weights.inc'
      real smoos(mfile),rminp(mfile)
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
      do j=1,npt
         xx=x(j)

c Lennard-Jones truncation function

         wfac(j)=truncate(xx,rminpt,rmaxpt)

c Coulomb truncation function

         wfacc(j)=truncatec(xx,rmaxpt)
      end do
      wfac(npt)=0.0
      wfacc(npt)=0
c
c set the empirical potentials to zero
c
      do ip=1,npairs
            do j=1,npt
                  xav(j,ip)=0.0
                  xavs(j,ip)=0.0
                  xav1(j,ip)=0.0
                  xavs1(j,ip)=0.0
            end do
      end do
c
c generate the EP from the site-site coefficients
c
      if(ireset.eq.0) then
            call generate_ep()
      endif
c 
c set up L-J and Coulomb potential between pairs.
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
            it=n1extr(ic)
            jt=n2extr(ic)
c
c Coulomb interactions are treated as between uniform spherical charge clouds
c of diameter equal to qradius(it)
c
            ra=qradius(it)
            rb=qradius(jt)
            ect=ecor(ic)
            qct=qcor(ic)
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

                  sumlj=sumlj*wfac(j)
                  sumc=sumc*wfacc(j)

c Save results

                  epot(j,ic)=sumlj
                  ecoul(j,ic)=sumlj+sumc
                  x1=x2
            end do
            epot(1,ic)=epot(2,ic)
            ecoul(1,ic)=ecoul(2,ic)
c
c set up gradients of these potentials
c
            sumlj1=epot(1,ic)
            sumc1=ecoul(1,ic)
            x1=xstar(1)
            x2=xstar(2)
            do j=2,npt
                  x2=xstar(j) 
                  xdif=x2-x1
                  sumlj2=epot(j,ic)
                  sumc2=ecoul(j,ic)
                  epot1(j-1,ic)=(sumlj2-sumlj1)/xdif
                  ecoul1(j-1,ic)=(sumc2-sumc1)/xdif
                  x1=x2
                  sumlj1=sumlj2
                  sumc1=sumc2
            end do
            epot1(npt,ic)=0.0
            ecoul1(npt,ic)=0.0
c
c Multiply the derivatives by the damping function, bearing in 
c mind that the derivatives have been evaluated at the bin boundaries.
c
c           do j=1,npt
c                 epot1(j,ic)=epot1(j,ic)*wfac(j)
c                 ecoul1(j,ic)=ecoul1(j,ic)*wfac(j)
c           end do
c
c Reintegrate the derivatives to give the potentials.
c Note that compared to xav1 and xavs1, the derivatives for epot and ecoul
c already correspond to their values at the bin boundaries so the algorithm 
c is slightly different
c
c           epot(npt,ic)=0.0
c           ecoul(npt,ic)=0.0
c           x1=xstar(npt)
c           do j=1,npt2
c                 jref=npt-j
c                 x2=xstar(jref)
c                 xt=epot1(jref,ic)
c                 xts=ecoul1(jref,ic)
c                 epot(jref,ic)=epot(jref+1,ic)-(x1-x2)*xt
c                 ecoul(jref,ic)=ecoul(jref+1,ic)-(x1-x2)*xts
c                 x1=x2
c           end do
      end do
      return
      end

      subroutine ptgnsq(nsum,nq,q,wav,data,qwind
     *,potfac,rminpt,rmaxpt,szeros,overallfac,efac,efacc,gfac
     *,rminp,potnoi,wtentp
     *,rbroad,sigma0,nsmoop)
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      include 'weights.inc'
      include 'X-rayf0.inc'

      real q(mpts),wav(2,mpz,mpairs),data(mpts,mfile)
     1,qwind(mpts,mpairs),qtemp(mpts,mfile),fittemp(mpts,mfile)
     1,dt(mpz),dts(mpz),potnoi(mfile)
     1,szeros(mfile),efac(mfile),rminp(mfile)

      integer nmin(mfile),nmini(mfile)
      character*14 formatout

      summax=exp(25.0)
      pi=3.141592653
      stepcl=cellst
      xavmax=0.0
      qstep=q(3)-q(2)
      qrat=2.0*pi/qstep
      expmax=25.0
      potadd=potfac
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
               conss0=data(nqmin(id),id)

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

               if(i.lt.nqmin(id)) then
                  qtemp(i,id)=conss0+grads0*(q(i)-q(1))-sum
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
                  call calc_pwrexpq(nqmax(id),nterms,pdstep,pdmax
     *,decayp,rho,q,expsq,0.0,ncfirst)
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
                  call fit_expansion(dt,expsq,dts,pwt,pwtold
     1,stepp,acceptr,rejrat,chisqnew,1,nqmax(id),2,nterms
     2,npits,nfreq)

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
                  pfact=efac(id)*potadd
                  do ic=2,nterms
c
c save the new fluctuation coefficients
c
c                       pcoefd(ic,id)=pwtold(ic)
                        pcoefd(ic,id)=efac(id)*pwtold(ic)
c
c increment the site-site potential coefficients, only if no Gaussians are present
c
                        do ip=1,npairs
                              irefp=irefpcof(ip)
c                              if(npgauss(irefp).eq.0) 
c     *pcoefp(ic,irefp)=efac(id)*pcoefp(ic,irefp)+pfact*cwtpot(ip,id)*pwtold(ic)
      pcoefp(ic,irefp)=efac(id)*pcoefp(ic,irefp)+pfact*cwtpot(ip,id)*pwtold(ic)
                        end do
                  end do
      write(6,*) id,efac(id),potadd,pfact,(pcoefp(ic,1),ic=1,5)
            end do
            open(10,file='delfit.dat',status='unknown')
            do i=1,nq
                  write(10,fmt=formatout) q(i),(fittemp(i,id),-qtemp(i,id),id=1,nfile)
            end do
            close(10)
      endif

      if(potadd.gt.0.0) then
c
c modify the repulsive terms as required
c
            pfact=potadd
            do ip=1,npairs
                  irefp=irefpcof(ip)
                  if(potex(ip).eq.0.0) then
c
c if potex(ic) is zero, then reduce the current pcoefex(ic) by 1% to
c ensure this term does not dominate the simulation
c
                        potex(ip)=-0.01*pcoefex(irefp)
                  endif
                  pcoefex(irefp)=pcoefex(irefp)+pfact*potex(ip)
            end do
c
c regenerate the empirical potentials if required
c
            call generate_ep()
      endif
c
c if no. of iterations to calculate potential is large (ireset=1) reset
c this value to input value for all subsequent iterations
c
      npits=npitss
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
            write(symbol_value(iref),*) expecf
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
            iref=iref+1
            write(symbol_value(iref),*) psigma2
c 
c Write the uniform distribution to a file called EPSR.uni
c
            extfile='.uni'
            fname=fnamerun(1:ilenfnamerun)//extfile(1:4)
            call openf(10,fname)
            write(10,*) (ycoh(j),j=1,npt)
            close(10)
      endif
      return
      end
      SUBROUTINE SMOOf(NQMIN,NQMAX,NSMOO,DIFF)
      DIMENSION DIFF(*)
C
C 3-POINT SMOOTHING keeping endpoints fixed
C
      n1=nqmin+1
      n2=nqmax-1
      if(n2.lt.n1.or.nsmoo.eq.0) return
      DO 20 IT=1,NSMOO
      SF=diff(nqmin)
      SM=DIFF(N1)
      sum1=0.25*(sf+sm)
      DO 10 I=N1,N2
      J=I+1
      SN=DIFF(J)
      sum2=0.25*(sm+sn)
      DIFF(I)=sum1+sum2
      SM=SN
      sum1=sum2
10    CONTINUE
20    CONTINUE
      RETURN
      END
      SUBROUTINE SMOOd(NQMIN,NQMAX,NSMOO,DIFF)
      DIMENSION DIFF(*)
C
C 3-POINT SMOOTHING keeping integrated area constant
C
      n1=nqmin
      n2=nqmax
      if(n2.lt.n1.or.nsmoo.eq.0) return
      DO 20 IT=1,NSMOO
      SM=DIFF(N1)
      diff(n2+1)=diff(n2)
      sum1=0.5*sm
      DO 10 I=N1,N2
      J=I+1
      SN=DIFF(J)
      sum2=0.25*(sm+sn)
      DIFF(I)=sum1+sum2
      SM=SN
      sum1=sum2
10    CONTINUE
20    CONTINUE
      RETURN
      END
      SUBROUTINE SMOOp(NQMIN,NQMAX,NSMOO,DIFF)
      DIMENSION DIFF(*)
C
C 3-POINT SMOOTHING keeping integrated area constant, zeroth and n+1'th point
c are held to zero to prevent large deviations developing
C
      n1=nqmin
      n2=nqmax
      if(n2.lt.n1.or.nsmoo.eq.0) return
      DO 20 IT=1,NSMOO
      SM=DIFF(N1)
      sum1=0.25*sm
      DO 10 I=N1,N2
      J=I+1
      if(i.lt.n2) then
            SN=DIFF(J)
      else
            sn=0.0
      endif
      sum2=0.25*(sm+sn)
      DIFF(I)=sum1+sum2
      SM=SN
      sum1=sum2
10    CONTINUE
20    CONTINUE
      RETURN
      END
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
c
c assign neutron weightings to individual atom pairs - these are initially
c all set to zero
c
      include 'epsr_dimension.inc'
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
      integer nq
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
c
c read factor for dataset to convert it to per atom. For X-rays this factor will decide
c whether data are normalised to the single atom scattering (1) or square of the mean
c scattering length (2).
c
            read(10,*) rhofac(in)

c decide whether this is an X-ray or neutron wts file

            ilenXray1=index(fnamew(in),'X-ray.wts')
            ilenXray2=index(fnamew(in),'x-ray.wts')
            ilenXray=max(ilenXray1,ilenXray2)

            if(ilenXray.gt.0) then

                  if(nxfile.lt.mxfile) then
                     nxfile=nxfile+1

c Save the normalisation of the X-ray data.

                     datanorm(nxfile)=rhofac(in)
      write(6,*) 'get_wts> ',nxfile,datanorm(nxfile)
                     rhofac(in)=1.0

c Setup pointers to the X-ray form factor arrays

                     xfile(in)=nxfile

c And to the dataset these weights correspond to.

                     datapoint(nxfile)=in

      write(6,*) in,nxfile,xfile(in),datapoint(nxfile)

                  else

                     write(6,109) mxfile
109   format(/'get_wts> Exceeded limit of ',i2,' X-ray files.')
                     goto 301

                  endif

            else
 
                  xfile(in)=0

            endif

300         continue

            if(ilenXray.gt.0) then


c Read form factor coefficients for this atom.

                  read(10,108,end=301) a1,a2
                  read(10,*) cxt,(axt(i),bxt(i),i=1,5),alphaxt,deltaxt
108   format(1x,a3,1x,a4)

c Assign these X-ray coefficients if the specified atom type matches an atom in the 
c .ato file.  Any atoms which are not assigned form factor coefficients are 
c automatically assumed not to contribute to this diffraction pattern.
c
                  call wxtype(a1,axt,bxt,cxt,alphaxt,deltaxt,nxfile,q,nq)
                  go to 300
            else

c Read partial structure factor assignment and neutron weighting for
c intramolecular and inter molecular terms.
c
                  read(10,107,end=301) a1,a2
                  read(10,*) wread(1),wread(2)
107   format(1x,a3,1x,a3)
c
c assign this neutron weighting if it matches a pair of atoms in the molecule
c Pairs which are not assigned a weighting by this input are left as zero.
c
                  call wtype(a1,a2,wread,in)
                  go to 300
            endif
301         continue
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
                  cprod=cfrac(it)*cfrac(jt)
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
               cprod=cfrac(it)*cfrac(jt)
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

               write(6,*) nqinv,q(qinv(nqinv)),iqinv,q(iqinv)

               do in=1,nfile

                  ix=xfile(in)
                  if(ix.gt.0) then
                     do ic=1,npairs
                        cwtpot(in,ic)=weightx(iqinv,ic,ix)
                     end do
                  else
                     do ic=1,npairs
                        cwtpot(in,ic)=cwt(in,ic)
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

      subroutine press(energy,pressc,energf,pressf,energabs)
c calculate energy and pressure from intermolecular interactions
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      include 'weights.inc'
      include 'en_atom.inc'
      real energf(mfile),pressf(mfile)
      real*8 prpot,presst

      sumc=0.
      presst=0.
      do if=1,nfile
            energf(if)=0.0
            pressf(if)=0.0
      end do
      imold=0
      jmold=0
c
c absolute energy
c
      energabs=0.0
      stepcl=cellst/dmax
      do ia=1,nal(nmol-1)
         im=mref(ia)
         if(im.ne.imold) then
            xc1=xadd(im)
            yc1=yadd(im)
            zc1=zadd(im)
            imold=im
         endif
         it=itype(ia)
         x1=xa(ia)+xc1
         y1=ya(ia)+yc1
         z1=za(ia)+zc1
         jafirst=naf(im+1)
         jalast=natom

c make a list of the neighbours of ia

         call make_neighbour_list(ia,jafirst,jalast,natom,nnbref
     *,nbref,ntableref)
c     write(6,*) ia,nnbref

         ija=0
         do while (ija.lt.nnbref)
            ija=ija+1
            ja=nbref(ija)
            jm=mref(ja)
            if(jm.ne.jmold) then
               xc2=xadd(jm)
               yc2=yadd(jm)
               zc2=zadd(jm)
               jmold=jm
            endif
            jt=itype(ja)
            icref=npairref(it,jt)
c
c set up separation of atom pair
c
            x2=xa(ja)+xc2
            y2=ya(ja)+yc2
            z2=za(ja)+zc2
            xdif=x2-x1
            ydif=y2-y1
            zdif=z2-z1
c
c apply minimum image convention
c
            call minimagec(xdif,ydif,zdif)
c
c separation of centres (for pressure calculation)
c
            xcdif=xc2-xc1
            ycdif=yc2-yc1
            zcdif=zc2-zc1
            call minimagec(xcdif,ycdif,zcdif)
            dsq=xdif*xdif+ydif*ydif+zdif*zdif
            dtest=sqrt(dsq)
c find the r value this corresponds to
            rat=dtest/stepcl
            iz=int(rat)+2
            dtdif=dtest-xstar(iz)
            ipz=iz
            if(dtest.gt.0.0) then
               prfac=(xcdif*xdif+ycdif*ydif+zcdif*zdif)/dtest
            else
               prfac=0.0
            endif
c sum reference potential
            if(iz.le.npt2) then
                  enpot=ecoul(iz,icref)+ecoul1(iz,icref)*dtdif
                  prpot=0.5*(ecoul1(ipz,icref)
     1+ecoul1(ipz-1,icref))*prfac
c
c sum intermolecular energy and pressure from empirical potential
c
                  enwt=xavs(iz,icref)+xavs1(iz,icref)*dtdif
c For absolute empirical energy do not include the repulsive term
                  enabswt=xav(iz,icref)+xav1(iz,icref)*dtdif
                  prwt=0.5*(xavs1(ipz,icref)
     1+xavs1(ipz-1,icref))*prfac
                  sumc=sumc+enwt+enpot
                  energabs=energabs+abs(enabswt)
                  presst=presst+(prwt+prpot)
                  do in=1,nfile
                     energf(in)=energf(in)+abs(cwt(in,icref)*enabswt)
                     pressf(in)=pressf(in)+abs(cwt(in,icref))*prwt
                  end do
            endif
         end do
      end do
c energy in kJ/mole

      energy=sumc/nmol
      energabs=energabs/nmol

c pressure in kbar

      pressc=((nmol*simt-presst/3.0)*100.0)/(6.0222*volume)
      do if=1,nfile
            energf(if)=energf(if)/(nmol*wsum(if))
            pressf(if)=((-pressf(if)/3.0)*100.0)
     1/(6.0222*volume*wsum(if))
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
      integer nq
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

c nrho + means data are per atom, -ve data are per molecule

         nrho=sign(1,datanorm(ix))

c iafformaff is 1 means use IAFFs for normalisation, otherwise use MAFFs

         iafformaff=(abs(datanorm(ix))+1)/2

c normx = 1 means normalise to single atom scattering, otherwise normalise to sum of Fs 
c all squared

         normx=abs(datanorm(ix))-2*(iafformaff-1)
         write(6,*) 'wxcreate> ',datanorm(ix),nrho,normx,iafformaff
      
c Form the Q-dependent weights for each atom pair

         do iq=1,nq

c Form the single atom term

            singlem=0.0
            singlei=0.0
            bbar=0.0
            do it=1,ntype

c Set up the single atom scattering for the modified and independent atom form factors

               singlem=singlem+cfrac(it)*fm0x(iq,it,ix)**2
               singlei=singlei+cfrac(it)*f0x(iq,it,ix)**2
               if(iafformaff.eq.1) then
                  if(normx.eq.1) then
                     bbar=bbar+cfrac(it)*f0x(iq,it,ix)**2
                  else
                     bbar=bbar+cfrac(it)*f0x(iq,it,ix)
                  endif
               else
                  if(normx.eq.1) then
                     bbar=bbar+cfrac(it)*fm0x(iq,it,ix)**2
                  else
                     bbar=bbar+cfrac(it)*fm0x(iq,it,ix)
                  endif
               endif

            end do

            if(normx.eq.2) bbar=bbar*bbar

            do ip=1,npairs

               n1=n1extr(ip)
               n2=n2extr(ip)
               weightx(iq,ip,ix)=fm0x(iq,n1,ix)*fm0x(iq,n2,ix)/bbar

c Correct for per molecule data if needed

               if(nrho.lt.0.and.normx.eq.2) 
     *weightx(iq,ip,ix)=weightx(iq,ip,ix)/amolrat

            end do

c Calculate background shift to simulated structure factor

            if(iafformaff.eq.1) then
               backx(iq,ix)=(singlem-singlei)/bbar
            else
               backx(iq,ix)=0.0
            endif
            if(nrho.lt.0.and.normx.eq.2) backx(iq,ix)=backx(iq,ix)/amolrat

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

      subroutine atogen(wav)
      parameter (pi=3.141592653,pi2=2.0*pi)
c
c if mode = 1 this calculates only intermolecular contribution, Any other
c value and it calculates the total contribution
c
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      real wav(2,mpz,mpairs)
      spi2=sqrt(pi2)
c 
c set weights to zero
c
      do iz=1,npt2
         do in=1,npairs
            wav(1,iz,in)=0.0
            wav(2,iz,in)=0.0
         end do
      end do
      if(npairs.gt.0) then
            do ic=1,npairs
                  potex(ic)=0.0
            end do
      endif
      stepcl=cellst/dmax
c
c generate intramolecular spacings
c
      do im=1,nmol
         if(nato(im).gt.1) then
            do ia=naf(im),nal(im)-1
               x1=xa(ia)
               y1=ya(ia)
               z1=za(ia)
               do ja=ia+1,nal(im)
                     it=itype(ia)
                     jt=itype(ja)
                     icref=npairref(it,jt)
                     x2=xa(ja)
                     y2=ya(ja)
                     z2=za(ja)
                     xx=(x1-x2)*sizefactor
                     yy=(y1-y2)*sizefactor
                     zz=(z1-z2)*sizefactor
c
c apply minimum image convention
c
                     call minimagec(xx,yy,zz)
                     dtest=sqrt(xx*xx+yy*yy+zz*zz)
c find the x value this corresponds to
                     iz1=nint(dtest/stepcl)+2
c
c accumulate the partial g(r) for this pair
c
                     if(iz1.gt.1.and.iz1.le.npt2) then
                        wav(1,iz1,icref)=wav(1,iz1,icref)+2.0
                     endif
               end do
            end do
         endif
      end do
c
c generate intermolecular spacings
c
      do im=1,nmol-1
         xc1=xadd(im)
         yc1=yadd(im)
         zc1=zadd(im)
         do jm=im+1,nmol
            xc2=xadd(jm)
            yc2=yadd(jm)
            zc2=zadd(jm)
            xcdif=xc2-xc1
            ycdif=yc2-yc1
            zcdif=zc2-zc1
c
c apply minimum image convention
c
            call minimagec(xcdif,ycdif,zcdif)
c
c try to avoid including molecules that will make zero contribution
c but we must allow the maximum radius of each molecule
c
            sumdiam=moldiam(im)+moldiam(jm)
            dtsq=xcdif*xcdif+ycdif*ycdif+zcdif*zcdif
            nt=int((sqrt(dtsq)-sumdiam)/stepcl)+2
            npt1test=npt2
            if(nt.le.npt1test) then
               do ia=naf(im),nal(im)
                  x1=xa(ia)+xc1
                  y1=ya(ia)+yc1
                  z1=za(ia)+zc1
                  do ja=naf(jm),nal(jm)
                        it=itype(ia)
                        jt=itype(ja)
                        icref=npairref(it,jt)
                        x2=xa(ja)+xc2
                        y2=ya(ja)+yc2
                        z2=za(ja)+zc2
                        xx=x1-x2
                        yy=y1-y2
                        zz=z1-z2
c
c apply minimum image convention
c
                        call minimagec(xx,yy,zz)
                        dtest=sqrt(xx*xx+yy*yy+zz*zz)
c find the x value this corresponds to
                        iz1=nint(dtest/stepcl)+2
919   format(1x,i5,2(1x,e12.5))
c
c accumulate the partial g(r) for this pair
c
                        if(iz1.gt.1.and.iz1.le.npt2) then
                           wav(2,iz1,icref)=wav(2,iz1,icref)+2.0
                        endif
c
c now check whether this pair of atoms is below a specified minimum distance
c
                        if(dtest*dmax.lt.rminex(irefpcof(icref))) then 
                           potex(icref)=potex(icref)+1.0
                        endif
                  end do
               end do
            endif
         end do
      end do
c
c normalize sums
c
      do ic=1,npairs
         it=n1extr(ic)
         jt=n2extr(ic)
         cprod=cfrac(it)*cfrac(jt)*natomr
         if(it.ne.jt) cprod=cprod*2
         if(cprod.gt.0.0) then
            do iz=1,npt2
               wav(1,iz,ic)=wav(1,iz,ic)/cprod
               wav(2,iz,ic)=wav(2,iz,ic)/cprod
            end do
            potex(ic)=potex(ic)/cprod
         else
            do iz=1,npt2
               wav(1,iz,ic)=0.0
               wav(2,iz,ic)=0.0
            end do
         endif
      end do
      return
      end

      subroutine gensq(nq,q,sofq,wav)
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      dimension q(mpts),sofq(2,mpts,mpairs),wav(2,mpz,mpairs)
      do ic=1,npairs
            do i=1,nq
                  sum1=0.0
                  sum2=0.0
                  sum3=0.0
c
c first check that this dataset contains intermolecular correlations
c
                  do j=2,npt1
                        add1=wav(1,j,ic)
                        add2=wav(2,j,ic)-ycoh(j)
                        sum1=sum1+add1*foqols(j,i)
                        sum2=sum2+add2*foqols(j,i)
                        sum3=sum3+wav(2,j,ic)
c                 write(6,999) i,j,add,foqols(j,i),sum
999   format(1x,2(i5,1x),4(e12.5,1x))
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
      return
      end

      subroutine fit_expansion(yin,expr,yout,wt,wtold
     1,stepr,acceptr,rejrat,chisqnew,nmin,npt
     2,ncfirst,nterms,nits,nfreq)
      real yin(*),yout(*),wt(*),wtold(*),expr(*)
c
c maximum value of a sum - this it to prevent overflows
c
      summax=exp(25.0)
c
c assign the initial coefficients
c
      do it=1,nterms
            wt(it)=wtold(it)
      end do
      if(nits.gt.1) then
          chisqold=0.0
          rejrat=0.0
          icount=0
          do it=1,nits
            chisqnew=0.0
            do i=nmin,npt
                  term=expand(nterms,i,expr,wt,summax)
                  chisqnew=chisqnew+(yin(i)-term)**2
            end do
            if(it.eq.1) chisqold=chisqnew
            accept=acceptr*chisqold
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
                            term=expand(nterms,i,expr,wt,summax)
                            chisqnew=chisqnew+(yin(i)-term)**2
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
          end do
          do ic=1,nterms
            wt(ic)=wtold(ic)
          end do
          rejrat=rejrat/nits
          write(6,*) chisqnew,rejrat
          if(rejrat.gt.0.5) then
            stepr=stepr*0.75/rejrat 
          else
            stepr=1.5*stepr
          endif
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
      real coeff(*),expr(*)
      sum=0.0
      do i=1,nterms
            jref=(j-1)*nterms+i
            sum=sum+expr(jref)*coeff(i)
            if(abs(sum).gt.summax) sum=sign(summax,sum)
      end do
      expand=sum
      return
      end

      subroutine calc_gaussq(lpt,nterms,dmin,dmax,decay,xin,expr)
      real xin(*),expr(*)
      pi=3.141592653
c
c maximum value of exponent
c
      expmax=25.0
c
c step in values of d
c
      if(nterms.gt.1) then
            dstep=(dmax-dmin)/(nterms-1)
      else
            dstep=1.0
      endif
c
c set up array of coefficients
c
      do m=1,nterms
c
c calculate Q dependence of mth term
c
            d=dmin+(m-1)*dstep
            if(d.lt.1.0) then
                  sigd=decay
            else
                  sigd=sqrt(d)*decay
            endif
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

      subroutine calc_gaussr(lpt,nterms,dmin,dmax,decay,rho,xin
     1,expr,rmin,sigfac)
      real xin(*),expr(*)
      pi=3.141592653
      piroot=sqrt(0.5*pi)
      gfac=piroot/(4*pi*pi*rho)
c
c maximum value of exponent
c
      expmax=25.0
c
c step in values of d
c
      if(nterms.gt.1) then
            dstep=(dmax-dmin)/(nterms-1)
      else
            dstep=1.0
      endif
c
c set up array of coefficients
c
      do m=1,nterms
c
c calculate Q dependence of mth term
c
            d=dmin+(m-1)*dstep
            if(d.lt.1.0) then
                  sigd=decay
            else
                  sigd=sqrt(d)*decay
            endif
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
c
c now repeat subtracting the broadened Gaussians if desired
c
            if(sigfac.ne.1.0.and.sigfac.gt.0.0) then
                if(d.ge.rmin) then
                  sigd=sigd*sigfac
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
                      expr(jref)=expr(jref)-gfact*sum
                  end do
                else
                  do i=1,lpt
                      jref=(i-1)*nterms+m
                      expr(jref)=0.0
                  end do
                endif
            endif
      end do
      return
      end

      subroutine calc_pwrexpq(lpt,nterms,dstep,dmax
     *,sigma,rho,xin,expr,rmin,ncfirst)
      real xin(*),expr(*)
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
c calculate Q dependence of mth term
c
            d=m*dstep
c
c save the integer of the first non-zero coefficient
c
            if(d.lt.rmin) ncfirst=m+1
            if(sigma.le.d.and.sigma.gt.0.0) then
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
      end do
      return
      end

      subroutine calc_pwrexpr(lpt,nterms,dstep,dmax
     *,sigmab,rho,xin,expr,rmin,ncfirst,fbroad,rbroad)
      include 'lgamma.inc'
      real xin(*),expr(*)
      pi=3.141592653
c
c maximum value of exponent
c
      expmax=25.0
      dstep2=0.5*dstep
      ntmax=dmax/sigmab+3
      call FCTRLN(ntmax)
c
c set up array of coefficients
c
      do m=1,nterms
c
c calculate r dependence of mth term
c
            d=m*dstep
c
c save the integer of the first non-zero coefficient
c
            if(d.lt.rmin) ncfirst=m+1
c
c generate value of sigma at this d value
c
            sigma=fbroad*sigmab
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
      end do
      return
      end

      subroutine calc_pwrexpr_grad(lpt,nterms,dstep,dmax
     *,sigmab,rho,xin,expr,rmin,ncfirst,fbroad,rbroad)
      include 'lgamma.inc'
      real xin(*),expr(*)
      pi=3.141592653
c
c maximum value of exponent
c
      expmax=25.0
      dstep2=0.5*dstep
      ntmax=dmax/sigmab+3
      call FCTRLN(ntmax)
c
c set up array of coefficients
c
      do m=1,nterms
c
c calculate r dependence of mth term
c
            d=m*dstep
c
c save the integer of the first non-zero coefficient
c
            if(d.lt.rmin) ncfirst=m+1
c
c generate value of sigma at this d value
c
            sigma=fbroad*sigmab
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
      end do
      return
      end

      subroutine stog(rho,q,sofq,r,gofr,nqmin,nqmax,nr,npairs
     *,nqdim,nrdim,npdim,xwt)
      real q(nqdim),r(nrdim),sofq(nqdim,npdim),gofr(nrdim,npdim)
     *,xwt(npdim)
      real lorch
      integer nqmin(npdim),nqmax(npdim)

      pi=3.141592653
      qstep=q(4)-q(3)
      sfac=0.5*qstep/pi/pi/rho
c
c perform a direct transform on sofq for r values, using a Lorch function to remove
c excess truncation oscillations
c
      pi=4*atan(1.0)
      qstep=q(4)-q(3)
      sfac=0.5*qstep/pi/pi/rho
      do j=1,nr
         rr=r(j)
         if(rr.gt.0.0) then
            do in=1,npairs
               gofr(j,in)=0.0
               qmax=pi/q(nqmax(in))
               do i=nqmin(in),nqmax(in)
                  qq=q(i)
                  qr=qq*rr
                  prod=qq*sin(qr)*lorch(qmax*qq,pi)
                  gofr(j,in)=gofr(j,in)+sofq(i,in)*prod
               end do
               gofr(j,in)=xwt(in)+sfac*gofr(j,in)/rr
            end do
         else
            do in=1,npairs 
               gofr(j,in)=0.0
               qmax=pi/q(nqmax(in))
               do i=nqmin(in),nqmax(in)
                  qq=q(i)
                  prod=qq*qq*lorch(qmax*qq,pi)
                  gofr(j,in)=gofr(j,in)+sofq(i,in)*prod
               end do
               gofr(j,in)=xwt(in)+sfac*gofr(j,in)
            end do
         endif
      end do
      return
      end

      function lorch(x,pi)

      real lorch,x,pi

      if(x.eq.0) then

         lorch=1.0

      else if(x.lt.pi) then

         lorch=sin(x)/x

      else

         lorch=0.0

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

      subroutine uniform(atry,xx,ycoh,xst,rho,npt,ndim)
      include 'lattice.inc'
      integer*8 ntries
c
c calculates the uniform distribution for a unit cell of arbitrary shape
c
      dimension xx(ndim),ycoh(ndim)
      volume=alatx*blaty*clatz
c
c each volume element is sampled with probability atry on average
c
      ntries=atry*volume/(xst*xst*xst)
      total=volume*rho
c
c check to prevent round off error
c
      add=total/(ycoh(1)+real(ntries))
      addlim=1.0e-4
      if(add.lt.addlim) then
            ycoh(1)=total/addlim-ntries
            if(ycoh(1).lt.0.0) then
                  ntries=ntries+ycoh(1)
                  ycoh(1)=0.0
            endif
      endif

      write(6,*)
      write(6,*) ntries

      ratio=ycoh(1)/total
      do j=2,npt
            ycoh(j)=ycoh(j)*ratio
      end do
      do i=1,ntries
            fa=ran1()-0.5
            fb=ran1()-0.5
            fc=ran1()-0.5     
            x=fa*alatx+fb*blatx+fc*clatx
            y=fb*blaty+fc*claty
            z=fc*clatz
            dtest=sqrt(x*x+y*y+z*z)
            inz=nint(dtest/xst)+2
            if(inz.le.npt) then
                  ycoh(inz)=ycoh(inz)+1.0
            endif
      end do
      ycoh(1)=ycoh(1)+real(ntries)
c
c renormalise result to number of atoms in box
c
      ratio=total/ycoh(1)

c determine the position of the maximum
      ymax=0.0
      jmax=2
      do j=2,npt
            ycoh(j)=ycoh(j)*ratio
c
c eliminate very small numbers to prevent overflows
c
            if(ycoh(j).lt.1.0e-7) ycoh(j)=0.0
            if(ycoh(j).gt.ymax) then
               ymax=ycoh(j)
               jmax=j
            endif
      end do

c Back off 1 step, in case the maximum was actually past the edge of the box

      jmax=jmax-1

c now calculate the region out to the edge of the box analytically

      pi=3.141592653
      rhofac=4.0*pi*rho/3.0
      x1=0.0
      x13=0.0
      xstep=xst*0.5
      do j=2,jmax
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
      real efacm(mfile)
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
         write(6,200) i,(cwt(i,j),j=1,npairs)
         do j=1,npairs
            cwt(i,j)=cwt(i,j)*feedback*efacm(i)
         end do
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

      subroutine generate_ep()
c
c generate EP from coefficients
c
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'epsrinput_par.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
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
      call calc_pwrexpr_grad(npt2,nterms,pdstep,pdmax
     *,decayp,rho,x,expgr,0.0,ncfirst,4.0,rbroad)
c
c step through the partials and calculate the empirical potential
c
      do ip=1,npairs
            irefp=irefpcof(ip)
c
c Set up the potential coefficients, ignoring the first one
c
            do ic=1,nterms
                  pwtold(ic)=pcoefp(ic,irefp)
            end do
c
c generate the potential corresponding to these coefficients
c
            sum1=0.0
            sums1=0.0
            do j=1,npt2
                  term=expand(nterms,j,expgr,pwtold,summax)
                  xav1(j,ip)=term
                  xavs1(j,ip)=term
            end do
      end do
c
c add on the extra repulsive term for the intermolecular potential
c and the Gaussian(s) if defined, then calculate the potential gradients.
c
      do ip=1,npairs

            irefp=irefpcof(ip)
            call calc_pwrexpr_grad(npt2,1,pdstep,pdmax
     *,expecf,rho,x,expgr,rminex(irefp),ncfirst,1.0,rbroad)
c
c only do something if intensity below rminex(ip) was found
c
            if(pcoefex(irefpcof(ip)).ne.0.0) then
                  do j=1,npt2
                        term=pcoefex(irefp)*expgr(j)
                        xavs1(j,ip)=xavs1(j,ip)+term
                  end do
            endif

c Add in Gaussians

            ig=0
            do while (ig.lt.npgauss(irefp))

               ig=ig+1
               do j=1,npt2

                  wpg=wpgauss(ig,irefp)
                  if(wpg.gt.0.0) then
                     rpg=rpgauss(ig,irefp)
                     hpg=hpgauss(ig,irefp)
                     rdiff=(x(j)-rpg)/wpg

c Gaussian definition for second and subsequent terms

                     if(ig.gt.1) then

                        term=log(abs(hpg))-0.5*rdiff*rdiff
                        if(abs(term).lt.25.0) then
                           add=-sign(1.0,hpg)*exp(term)*rdiff/wpg
c     write(6,*) j,x(j),rpg,wpg,hpg,rdiff,term,add
                           xavs1(j,ip)=xavs1(j,ip)+add
                        endif

                     else

c Morse potential definition for the first term

                        term=-rdiff
                        fac=4*abs(hpg)/wpg
                        if(term.gt.25.0) then
                           term=exp(25.0)
                        else if(term.lt.-25.0) then
                           term=0.0
                        else
                           term=exp(term)
                        endif
                        add=fac*term*(1-2*term)
                        xavs1(j,ip)=xavs1(j,ip)+add

                     endif
                  endif

               end do
            end do
                  
            do j=1,npt2
c
c gradients need to be dimensionless, and multiply by truncation function
c
                  xav1(j,ip)=xav1(j,ip)*dmax
                  xavs1(j,ip)=xavs1(j,ip)*dmax
            end do
c
c from gradients, calculate empirical potentials, integrating
c from the last point
c
            xav(npt,ip)=0.0
            xav1(npt,ip)=0.0
            xavs(npt,ip)=0.0
            xavs1(npt,ip)=0.0
            x1=xstar(npt)
            xt1=xav1(npt,ip)
            xts1=xavs1(npt,ip)
            do j=1,npt2
                  jref=npt-j
                  x2=xstar(jref)
                  xt2=xav1(jref,ip)
                  xts2=xavs1(jref,ip)
                  xav(jref,ip)=xav(jref+1,ip)-0.5*(x1-x2)*(xt1+xt2)
                  xavs(jref,ip)=xavs(jref+1,ip)-0.5*(x1-x2)*(xts1+xts2)
                  x1=x2
                  xt1=xt2
                  xts1=xts2
            end do

c Multiply empirical potential by truncation function

            do j=1,npt2

               xav(j,ip)=xav(j,ip)*wfac(j)
               xavs(j,ip)=xavs(j,ip)*wfac(j)

            end do

c Now re-evaluate the derivatives with these truncated potentials

            sum1=xav(1,ip)
            sums1=xavs(1,ip)
            x1=xstar(1)
            x2=xstar(2)
            do j=2,npt
                  x2=xstar(j) 
                  xdif=x2-x1
                  sum2=xav(j,ip)
                  sums2=xavs(j,ip)
                  xav1(j-1,ip)=(sum2-sum1)/xdif
                  xavs1(j-1,ip)=(sums2-sums1)/xdif
                  x1=x2
                  sum1=sum2
                  sums1=sums2
            end do
            xav1(npt,ip)=0.0
            xavs1(npt,ip)=0.0


      end do
      return
      end

      subroutine lsparabola(nchan,fname,nsum,a,b,c)

      character*256 fname
      real*8 x,y,sum,sumx1,sumx2,sumx3,sumx4,sumy,sumyx1,sumyx2
     *,delta,cplus,cminus,bplus,bminus,pf1,pf2

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
      navmax=100
      nsum=0

      do while (ierr.eq.0)
         read(nchan,*,iostat=ierr) a,a,y,x
         if(ierr.eq.0) then
            nsum=nsum+1
            if(x.gt.xmax) xmax=x
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
         endif
      end do


c Now work out the coefficients the hard way! The method used here assumes all the sums are 
c positive definite

c a) determinant of coefficients

      delta=sum*(sumx2*sumx4-sumx3*sumx3)
     *+sumx1*(sumx3*sumx2-sumx1*sumx4)
     *+sumx2*(sumx1*sumx3-sumx2*sumx2)

      if(abs(delta).gt.1.0e-10.and.nsum.ge.3) then

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
      write(6,*) nsum,delta,a,b,c

c write a little file to show the fit

      open(8,file=fname,status='unknown')

      npt=100
      xstep=xmax/npt

      do ix=1,npt

         x=ix*xstep

c format here has to match that of the .erg file so that it can be plotted with
c the plot command.

         write(8,100) x,x,(a+x*(b+x*c)),x

100   format(4(1x,e13.6))

      end do

      close(8)

      return
      end

