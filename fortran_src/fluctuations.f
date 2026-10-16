      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_fluctuations()

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
      include 'lattice.inc'
      include 'ran1com.inc'
      integer*4 ierrout
      character*20 text
      real*8 sumnumber(mdist),sumnumbersqr(mdist),localnumber(mdist),localnumbersqr(mdist)
      real*8 frac,frac1,szero(mdist)
      real damaxmin
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
c set up the number of density bins
c
      pi=4.0*atan(1.0)
      iref=2
      read(symbol_value(iref),*) nz
      nz=min(nz,mz)
c
c Get the number of samples already probed
c
      iref=iref+1
      read(symbol_value(iref),*) nsumts
      if(nsumts.lt.0) then
            nsumt=0
      else
            nsumt=nsumts
      endif


c Set the number of samples to take from the current configuration

      iref=iref+1
      read(symbol_value(iref),*) npix

c Get the number of different volumes to sample

      ndist=nfiles
      if(ndist.gt.mdist) ndist=mdist

c Find the largest value of damax

      damaxmin=0.0
      do id=1,ndist
c
c Set up the radius of the volume to test
c
            iref=ireff(nsymboluptodata+1,id)
            read(symbol_value(iref),*) damax(id)
            if(id.eq.1) then
               damaxmin=damax(id)
            else
               damaxmin=min(damaxmin,damax(id))
            end if

c Average of N so far

            iref=iref+1
            read(symbol_value(iref),*) sumnumber(id)

c Average of N^2 so far

            iref=iref+1
            read(symbol_value(iref),*) sumnumbersqr(id)

c S(0) from number fluctuations

            iref=iref+1
            read(symbol_value(iref),*) szero(id)

            if(nsumt.eq.0) then

                  sumnumber(id)=0.0
                  sumnumbersqr(id)=0.0

            endif

      end do

c
c setup step size and z values
c
      rho=nmol/(alatx*blaty*clatz)
      zstep1=2.0*rho/float(nz)
c      zstep2=1.0/(damaxmin*damaxmin*damaxmin)
c      zstep=max(zstep1,zstep2)

c      zstep=1.0/sqrt(4.0*pi*rho*damaxmin*damaxmin*damaxmin/3.0)
      zstep=1
      write(6,*) nz,ndist,zstep

      do j=1,nz
            z(j)=real(j-1-nz/2)*zstep
      end do

      call initgofz(nz,ndist,gofz,mz,mdist)
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

c set up the table of numbers of molecules in each sample volume

      call find_fluctuations(rho,localnumber,localnumbersqr)

c
c axis title
c
      titler='Density'
      do id=1,ndist
c
c set the column title for this distribution
c
            write(text,100) damax(id)
100   format(e12.5)
            text=adjustl(text)
            ilent=index(text,' ')-1
            titledat(id)='r_max:'//text(1:ilent)
      end do

c output the data

      ext='n'

      call outgofz(fnameinp,ext,nz,ndist,nsumt,z,gofz,gzt,gzt2
     *,titler,titledat,mz,mdist)

c Update the sums

      frac=1.0/real(nsumt)
      frac1=1.0-frac

      do id=1,ndist

            sumnumber(id)=sumnumber(id)*frac1+localnumber(id)*frac
            sumnumbersqr(id)=sumnumbersqr(id)*frac1+localnumbersqr(id)*frac
            if(sumnumber(id).gt.0.0) then
               szero(id)=sumnumbersqr(id)/sumnumber(id)-sumnumber(id)
            else
               szero(id)=0.0
            endif

      end do

c update the current value of nsumt and ndist

      icount=0
      text='nsize'
      call search_for_symbol(text,icount)
      iref=ireff(icount,1)
      write(symbol_value(iref),*) nz
      iref=iref+1
      if(nsumts.ge.0) nsumts=nsumt
      write(symbol_value(iref),*) nsumts
      iref=iref+1
      write(symbol_value(iref),*) npix
      iref=iref+1
      write(symbol_value(iref),*) ndist

      do id=1,ndist
         iref=iref+1
         write(symbol_value(iref),*) damax(id)
         iref=iref+1
         write(symbol_value(iref),*) sumnumber(id)
         iref=iref+1
         write(symbol_value(iref),*) sumnumbersqr(id)
         iref=iref+1
         write(symbol_value(iref),*) szero(id)
      end do

c write input file

      call save_files_others(14)

c Save random number parameters.

      open(14,file='randomno.dat',status='unknown')
      call saveran1(14)
      close(14)

      return
      end

      subroutine find_fluctuations(rho,sumnumber,sumnumbersqr)
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'clusters.inc'
      include 'lattice.inc'
      real*4 rho,voidtest(mdist),xcent,ycent,zcent,pi
      real*8 sumnumber(mdist),sumnumbersqr(mdist)
      integer*4 rhotest(mdist)
      logical test

      pi=4.0*atan(1.0)

c Renormalise maximum radius for each distribution

      do id=1,ndist

         voidtest(id)=damax(id)/dmax

c Setup average number of atoms in this volume

         rhotest(id)=nint(4.0*pi*rho*damax(id)*damax(id)*damax(id)/3.0)

c zero local sums

         sumnumber(id)=0.0
         sumnumbersqr(id)=0.0

      end do

      i=0
      do while(i.lt.npix)

c Choose a position at random within the box

            fa=(ran1()-0.5)
            fb=(ran1()-0.5)
            fc=(ran1()-0.5)

c Convert to Cartesian coordinates

c
c calculate corrected coordinates
c
            zcent=fc*clatz
            ycent=fb*blaty+fc*claty
            xcent=fa*alatx+fb*blatx+fc*clatx

            test = .false.
            if(.not.test) then

               i=i+1

c Zero accumulators

               do id=1,ndist
                  ltyp(id)=0
               end do

c Look at all the other molecules and pick those that are within the box of dimenstion damax(id)/dmax
c look at each requested atom radius

               do im=1,nmol

c Calculate the separation of these molecules

                        xdif=xadd(im)-xcent
                        ydif=yadd(im)-ycent
                        zdif=zadd(im)-zcent

c Apply minimum image convention

                        call minimagec(xdif,ydif,zdif)
                        diftest=sqrt(xdif*xdif+ydif*ydif+zdif*zdif)

                        do id = 1,ndist
                           if(diftest.le.voidtest(id)) then

                              ltyp(id)=ltyp(id)+1

                           endif
                        end do
               end do

c Accumulate number sums

               do id=1,ndist

                  sumnumber(id)=sumnumber(id)+ltyp(id)
                  sumnumbersqr(id)=sumnumbersqr(id)+ltyp(id)*ltyp(id)

c density of this box
            
                  size=damax(id)
                  iref=ltyp(id)-rhotest(id)+1+nz/2
c                  iref=int((real(ltyp(id))-rhotest(id))/(zstep*sqrt(rhotest(id))))+1+nz/2
                  if(iref.lt.1) iref=1
                  if(iref.gt.nz) iref=nz
                  gofz(iref,id)=gofz(iref,id)+1

               end do
            endif
      end do

c Normalise sums to number of terms

      do id=1,ndist

            sumnumber(id)=sumnumber(id)/npix
            sumnumbersqr(id)=sumnumbersqr(id)/npix 

      end do

c Normalise the density distributions to probabilities

      do id=1,ndist

            sum=0.0
            do iz=1,nz
                sum=sum+gofz(iz,id)
            end do
            sum=sum*zstep
            if(sum.gt.0.0) then
               do iz=1,nz
                  gofz(iz,id)=gofz(iz,id)/sum                       
               end do
            else
               do iz=1,nz
                  gofz(iz,id)=0.0
               end do
            endif
      end do
      return
      end


