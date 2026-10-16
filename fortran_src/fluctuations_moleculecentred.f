      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_fluctuations()

      stop
      end      

      subroutine run_fluctuations()
c
c programme to estimate the density fluctuations within a specified radius

c Voids are estimated by pixellating the box into a grid of small cubes, then
c estimating the distance from each pixel to the nearest atom. If this distance
c is greater than a specified distance - a nominal radius of the atoms -
c then the pixel is marked as being in a void.

c Then look at each pixel and check the surrounding 26 neighbouring pixels. 
c If any of them are occupied, then both pixels are given the same 
c (normally the lowest non-zero) void number of the two of them.
c of molecules.

c Pixels that are not in a void are put in a void of zero volume.

c Finally a histogram is built up of the number pixels in each void. These 
c numbers are divided by the total number of pixels in the box, so that the 
c effect of different pixel sizes can be checked.

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'clusters.inc'
      include 'ran1com.inc'
      integer*4 ierrout
      character*20 text
      real naverage,naverage2,szero
c
c setup the clusters symbols
c
      call setup_symbols_fluctuations()

c
c now setup the input file
c
      nin=11
      nout=6
      call setup_input_file(nin,nout,ierrout)
c
c if an error was encountered on setup, exit now
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
c Get the number of samples already probed
c
      iref=2
      read(symbol_value(iref),*) nsumt

c Set the number of samples to take from the current configuration

      iref=iref+1
      read(symbol_value(iref),*) npix

c Get the number of different volumes to sample

      ndist=nfiles
      if(ndist.gt.mdist) ndist=mdist

      do id=1,ndist
c
c Set up the radius of the volume to test
c
            iref=ireff(nsymboluptodata+1,id)
            read(symbol_value(iref),*) damax(id)

c make sure this value is not larger than half the box size

            if(damax(id).gt.0.5*dmax) damax(id)=0.5*dmax

c Sum of N so far

            iref=iref+1
            read(symbol_value(iref),*) damin(id)

c Sum of N^2 so far

            iref=iref+1
            read(symbol_value(iref),*) fract2(id)

c Zero these sums if nsumt <=0

            if(nsumt.le.0) then
                  damin(id)=0.0
                  fract2(id)=0.0
            endif
      end do
c
c normalise atomic coordinates to box size
c
      call normalise_coordinates()

c Initialise random number generator - this must be done here, after the .ato
c file has been read, otherwise the random generator is always set to the values
c read in from the .ato file. At the end of voids, the .ato file is not saved.

      open(14,file='randomno.dat',status='unknown')
      call initran1(14)
      write(6,*) mcall,idum,iy
      close(14)

c set up the table of molecular bonds

      call find_fluctuations()

c update the current value of nsumt and ndist

      icount=0
      text='nsumt'
      call search_for_symbol(text,icount)
      iref=ireff(icount,1)
      write(symbol_value(iref),*) nsumt
      iref=iref+1
      write(symbol_value(iref),*) npix
      iref=iref+1
      write(symbol_value(iref),*) ndist

      do id=1,ndist
         iref=iref+1
         write(symbol_value(iref),*) damax(id)
         iref=iref+1
         write(symbol_value(iref),*) damin(id)
         iref=iref+1
         write(symbol_value(iref),*) fract2(id)
      end do

c write input file

      call save_files_others(14)
c
c write out parameters and zero limit
c
      do id=1,ndist
            naverage=damin(id)/nsumt
            naverage2=fract2(id)/nsumt
            szero=(naverage2-naverage*naverage)/naverage
            write(14,*) damax(id),naverage,naverage2,szero
      end do
      close(14)

c Save random number parameters.

      open(14,file='randomno.dat',status='unknown')
      call saveran1(14)
      close(14)

      return
      end

      subroutine find_fluctuations()
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'clusters.inc'
      real*4 voidtest(mdist),xcent,ycent,zcent
      integer nentries,entry(1000)
      logical test

c limit number of entries to 1000 or nmol, whichever is smaller

      if(npix.gt.1000) npix=1000
      if(npix.gt.nmol) npix=nmol

c Setup maximum distance for each distribution

      do id=1,ndist

         voidtest(id)=damax(id)/dmax

      end do
      nentries=0
      i=0
      do while(i.lt.npix)

c Choose a molecule at random

            imol = ran1()*nmol+1

c Check that this molecule has not already been used in the list of entries.

            test = .false.
            itest=0
            do while (.not.test.and.itest.lt.nentries)
                  itest=itest+1
                  test=imol.eq.entry(itest)
            end do
            if(.not.test) then
               nentries=nentries+1
               entry(nentries)=imol

c Zero accumulator for this molecule

               do id=1,ndist

                  ltyp(id)=0

               end do

               i=i+1
               nsumt=nsumt+1

c Look at all the other molecules and pick those that are within the specified distancec
c look at each requested atom radius

               do im = 1,nmol

                  if(im.ne.imol) then

c Calculate the separation of these molecules

                        xdif=xadd(im)-xadd(imol)
                        ydif=yadd(im)-yadd(imol)
                        zdif=zadd(im)-zadd(imol)

c Apply minimum image convention

                        call minimagec(xdif,ydif,zdif)

                        d=sqrt(xdif*xdif+ydif*ydif+zdif*zdif)
                        do id = 1,ndist
                           if(d.le.voidtest(id)) then

                              ltyp(id)=ltyp(id)+1

                           endif
                        end do
                  endif
               end do
                  
c Accumulate N and N^2

               do id=1,ndist

                  damin(id)=damin(id)+ltyp(id)
                  fract2(id)=fract2(id)+ltyp(id)*ltyp(id)
               
               end do
            endif
      end do                        

      return
      end


