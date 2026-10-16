      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_voids()

      stop
      end      

      subroutine run_voids()

      use reallocation_routines

      interface

      subroutine find_voids(nexclude,exclude)
         include 'epsr_dimension.inc'
         include 'molecom.inc'  
         include 'clustercom.inc'
         include 'clusters.inc'
         parameter (mpix3=mpix*mpix*mpix,maxvoid=2**15-1)
         integer*1 voidpix(mpix3)
         integer*2 voidnum(mpix3)
         integer*4 ih,ik,il,npixtot,npixall,nvoidpix
         integer*4 labvoid(maxvoid),lenvoid(maxvoid)
         integer(kind=4), dimension(:) :: nexclude
         character(len=3), dimension(:,:) :: exclude
      end subroutine find_voids

      end interface

c
c programme to estimate the distribution of void sizes within a structure.

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
      include 'lattice.inc'
      include 'clusters.inc'
      include 'ran1com.inc'
      integer*4 ierrout,nexcludemax
      character*20 text
      integer(kind=4), dimension(mdist) :: nexclude
      character(len=3), dimension(:,:), allocatable :: exclude
      character*256 texclude
      logical found

      nexcludemax=0
c
c setup the clusters symbols
c
      call setup_symbols_voids()

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
      if(.not.cubic) then
         write(6,*) 'voids> Cannot currently process non-cubic boxes!! Apologies!!'
         return
      end if
c
c set up the number of voids size bins
c
      iref=2
      read(symbol_value(iref),*) nz
      nz=min(nz,mz)
c
c set up number of pixels
c
      iref=iref+1
      read(symbol_value(iref),*) npix

c Ensure the total number of pixels across the box is an odd number, so there is
c always a pixel at the centre of the box.

      npix=npix/2
      npixtotal=2*npix+1
      if(npixtotal.gt.mpix) npix=(mpix-1)/2
      write(6,*) npix,npixtotal

c Setup maximum void volume in units of Angstroms^3

      iref=iref+1
      read(symbol_value(iref),*) voidmax
c
c set up current accumulator
c
      iref=iref+1
      read(symbol_value(iref),*) nsumt
c
c check the number of distributions does not exceed the array size
c
      ndist=nfiles
c
c Two sets of distributions are calculated: 
c 1) a histogram of void sizes
c 2) a void-pixel pair correlation function.
c
      ndist3=3*ndist
      if(ndist3.gt.mdist) ndist3=mdist
      ndist=ndist3/3

      do id=1,ndist
            iref=ireff(nsymboluptodata+1,id)
c
c setup the distance used to define whether a pixel is within an atom
c
            read(symbol_value(iref),*) damax(id)

            iref=iref+1
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.gt.nexcludemax) then
               nexcludemax=nwords
               call reallocate2d_c(exclude,3,nexcludemax,ndist)
            end if
            nexclude(id)=0
            ie=0
            do while(ie.lt.nwords)
               ie=ie+1
               call concatenate(ie,ie,texclude,ilentexclude
     1,nwords,ncf,ncl,sometext)
! Check whether this atom type exists
               found=.false.
               it=0
               do while (it.lt.ntype.and..not.found)
                  it=it+1
                  found=type(it).eq.texclude(1:3)
               end do
               if(found) then
                  nexclude(id)=nexclude(id)+1
                  exclude(nexclude(id),id)=texclude(1:3)
               end if
            end do
            if(nexclude(id).gt.0) then
               write(6,*) id,damax(id),(exclude(ie,id),ie=1,nexclude(id))
            else
               write(6,*) id,damax(id)
            endif

      end do
      write(6,*) nz,ndist

c setup step size and z values
c voidmax represents the maximum void volume sampled in units of the mean volume
c per atom, i.e. 1/rho

      zstep=dmax/nz
      zhalf=0.5*zstep
      do j=1,nz
            z(j)=j*zstep-zhalf
      end do

      call initgofz(nz,ndist3,gofz,mz,mdist)

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
c
c set up the table of molecular bonds
c
      call find_voids(nexclude,exclude)
c
c axis title
c
      titler='Void_radius'
      idref3=0
      do id=1,ndist
      idref1=idref3+1
      idref2=idref1+1
      idref3=idref2+1
c
c set the column title for this pair of atoms
c
         write(titledat(idref1),102) damax(id)
         write(titledat(idref2),102) damax(id)
         write(titledat(idref3),102) damax(id)
102   format('   Radius_',f5.3)
      end do

c output the data

      ext='n'

      call outgofz(fnameinp,ext,nz,ndist3,nsumt,z,gofz,gzt,gzt2
     *,titler,titledat,mz,mdist)

c For each distribution calculate the fraction of box volume occupied by voids

      boxv=dmax*dmax*dmax
      idref3=0
      do id=1,ndist

            idref1=idref3+1
            idref2=idref1+1
            idref3=idref2+1
            sum1=0.0
            do i=1,nz
            
c Void linear dimension

                  voidl=z(i)

c Void volume

                  voidv=voidl*voidl*voidl
                  sum1=sum1+voidv*gzt(i,idref1)

            end do
            fract2(id)=sum1/boxv

      end do

c
c write input file
c

c update the current value of nsumt and ndist

      icount=0
      text='nsize'
      call search_for_symbol(text,icount)
      iref=ireff(icount,1)
      write(symbol_value(iref),*) nz
      iref=iref+1
      write(symbol_value(iref),*) nsumt
      iref=iref+1
      write(symbol_value(iref),*) ndist

      call save_files_others(14)
      
      do id=1,ndist
c
c for each void distribution indicate the fraction of unoccupied pixels
c
            write(14,101) id,fract2(id)
      end do
101   format(/'Distribution ',i3,': fraction of pixels used = '
     *,f10.5)

      close(14)

c Save random number parameters.

      open(14,file='randomno.dat',status='unknown')
      call saveran1(14)
      close(14)

      return
      end

      subroutine find_voids(nexclude,exclude)
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'clusters.inc'

      parameter (mpix3=mpix*mpix*mpix,maxvoid=2**15-1)
      integer*1 voidpix(mpix3)
      integer*2 voidnum(mpix3)
      integer*4 ih,ik,il,npixtot,npixall,nvoidpix
      integer*4 labvoid(maxvoid),lenvoid(maxvoid)
      integer(kind=4), dimension(:) :: nexclude
      character(len=3), dimension(:,:) :: exclude

      nfdim=mdist

c Mean volume per atom

      volatom=1.0/real(natom)

      zstep=1.0/nz

c Keep track of the two distributions

      idref3=0

c
c look at each requested atom radius
c
      do id=1,ndist

c Keep track of the three distributions

         idref1=idref3+1
         idref2=idref1+1
         idref3=idref2+1

c set up the test radius

         radtest=damax(id)/dmax

c Parameters to define the pixel positions 

         pixstep=1.0/real(2*npix+1)

c No. of pixels needed to define sphere.

         npixsphere=nint(radtest/pixstep)

c Need to ensure pixel size is adjusted to give closest possible match to 
c specified radius

c Give the ideal pixel size

         pixstep=radtest/npixsphere

c Convert this to the ideal number of pixels. The final number must be the
c the nearest ODD integer.

         ratio=1.0/pixstep
         npixt=int(ratio/2.0)

c Check the remainder to be sure this the nearest ODD integer

         remain=ratio-(2*npixt+1)
      write(6,*) ratio,npixt,remain
         if(remain.gt.1.0) npixt=npixt+1

         if((2*npixt+1).gt.mpix) npixt=(mpix-1)/2

c Get the correct pixel size

         pixstep=1.0/real(2*npixt+1)

c No. of pixels needed to define sphere.

         npixsphere=nint(radtest/pixstep)
         radactual=npixsphere*pixstep

         npixtot=2*npixt+1

c Total number of pixels

         npixall=npixtot*npixtot*npixtot

         write(6,103) id,npixtot,radtest*dmax,radactual*dmax
103   format(//'find_voids> Distribution ',i2
     1//'find_voids> Adjusting number of pixels to give optimum'
     1/ '            fit to specified void distance.',
     1//'find_voids> New number of pixels = ', i8,
     1//'find_voids> Specified radius = ',f10.5,
     1//'find_voids> Actual radius used = ',f10.5)

c Mark all pixels as unoccupied initially

         do ipix=1,npixall

            voidpix(ipix)=0
            voidnum(ipix)=0

         end do

c Step through all the atoms

         ia=0
         do while(ia.lt.natom)
            ia=ia+1
            ifound=0
            if(nexclude(id).gt.0) then ! If excluded atoms have been defined, the occupied pixels are defined ONLY for those atoms
               ie=0  
               do while(ie.lt.nexclude(id).and.ifound.eq.0)
                  ie=ie+1
                  if(type(itype(ia)).eq.exclude(ie,id)) ifound=itype(ia)
               end do
!      write(6,*) ia,ie,type(itype(ia)),exclude(ie,id),ifound
            else ! Otherwise all atoms create occupied pixels.
               ifound=itype(ia)
            end if
!            write(6,*) ia,ifound
            if(ifound.ne.0) then
               im=mref(ia)
               xatom=xa(ia)+xadd(im)
               yatom=ya(ia)+yadd(im)
               zatom=za(ia)+zadd(im)

c Just as a precaution apply minimum image to make sure this atom is definitely 
c in the range +- 0.5 in each direction

               call minimagec(xatom,yatom,zatom)

c Find the closest pixel

               ih=nint(xatom/pixstep)
               ik=nint(yatom/pixstep)
               il=nint(zatom/pixstep)

c Mark this pixel as occupied, plus all the surrounding pixels within the 
c sphere.

c     write(6,*) xatom,yatom,zatom,ih,ik,il,nvoidpix
               call occupypix(npixt,npixsphere,npixfound,ih,ik,il,voidpix,mpix3)

            endif
         end do

! If nexclude(id) > 0, we are interested in the way the 'excluded' atoms are distributed, 
! so instead of calculating the distribution of unoccuped pixels, we calculate the distribution of occupied pixels.
! To use the same code we now need to define the currently occupied pixels as unoccupied, and the currently unoccupied ones as occupied.

         if(nexclude(id).gt.0) then
            do ipix=1,npixall
               if(voidpix(ipix).ne.0) then
                  voidpix(ipix)=0
               else
                  voidpix(ipix)=1
               end if
            end do
         end if

         write(6,102) npixfound
102   format(/'find_voids> No. of pixels in sphere = ',i8)

c Set the initial number of unoccupied pixels

         nvoidpix=0

c Count the number of unoccupied pixels

         do ipix=1,npixall

            if(voidpix(ipix).eq.0) then
               nvoidpix=nvoidpix+1
            endif

         end do

c Write the total number of unoccupied pixels

         fract2(id)=real(nvoidpix)/real(npixall)

         write(6,101) nvoidpix,npixall,fract2(id)
101   format(/'find_voids> ',i10,' pixels used out of '
     1,i10,' pixels in total.'
     1,//'find_voids> Fraction used = ',f10.5)

c Now correlate neighbouring unoccupied pixels into voids.
c This is done in three stages. First we make a list of the unoccupied pixels
c that are immediately adjacent to each unoccupied pixel. Then we assign all
c of these unoccupied, adjacent pixels the same (lowest) void number. Finally we
c count up the total number of voids. The latter two processes are repeated
c until there is no further reduction in the number of voids.

c Initially we assume each pixel is its own void

         nvoidsold=nvoidpix
         nvoidsnew=0
         it=0

         do while (nvoidsold.gt.nvoidsnew)

            if(it.gt.0) nvoidsold=nvoidsnew

            call search_voids(nvoidsnew,npixt,npixall
     1,voidpix,voidnum,mpix3)

            it=it+1
            write(6,*) it,nvoidsold,nvoidsnew

         end do

c Finally the size of each void is tallied and formed into a histogram.

c First zero the tally and label

         do ivoid=1,nvoidsnew

            labvoid(ivoid)=0
            lenvoid(ivoid)=0

         end do

         nvoid=0

         do ipix=1,npixall

c Get an unoccupied void with finite void number

            if(voidpix(ipix).eq.0.and.voidnum(ipix).gt.0) then

c Check that this void number does not already exist

               iflag=0
               ivoid=0
               do while(ivoid.lt.nvoid.and.iflag.eq.0)

                  ivoid=ivoid+1
                  if(labvoid(ivoid).eq.voidnum(ipix)) iflag=ivoid

               end do

               if(iflag.eq.0) then

                  nvoid=nvoid+1
                  labvoid(nvoid)=voidnum(ipix)
                  lenvoid(nvoid)=1

               else

                  lenvoid(iflag)=lenvoid(iflag)+1

               endif

            endif
         end do

c Form a histogram of number of voids versus volume. 

         ivoid=0
         volpix=pixstep*pixstep*pixstep
         volpixstar=volpix/volatom
         third=1.0/3.0
         do while(ivoid.lt.nvoid)

            ivoid=ivoid+1

            if(ivoid.lt.51) write(6,*) ivoid,lenvoid(ivoid)*volpixstar

c The size of a void is expressed as the cube root of its volume

            iref=int(((real(lenvoid(ivoid))*volpix)**third)/zstep)+1
            if(iref.gt.nz) iref=nz
            gofz(iref,idref1)=gofz(iref,idref1)+1

         end do

c Now calculate the radial distribution function of void pixels. We ensure we
c only count each pair once.

         npixtot=2*npixt+1
         npixtot2=npixtot*npixtot
         ivcount=0
         nvcount=100
         ivoidsum=0

c In the event that the step in z is much larger than the pixel size,
c we will broaden the bin that the data are entered into

         nbroad=nint(real(nz)/real(npixtot))
c Choose an amount to broaden each point
         pixhwidth=0.5/real(npixtot)
         if(nbroad.lt.0) nbroad=0

         do while(ivcount.lt.nvcount)

            ivcount=ivcount+1

c Choose the origin pixel at random
            
            itest=0

            do while (itest.eq.0)

               iv1=npixall*ranmine()+1
               if(voidpix(iv1).eq.0) itest=1

            end do

            if(voidpix(iv1).eq.0) then
!Include the correlation of this pixel with itself.
               iscount=iscount+1
c Add the broadening component
               xdif=pixhwidth*(2.0*ranmine()-1.0)
               ydif=pixhwidth*(2.0*ranmine()-1.0)
               zdif=pixhwidth*(2.0*ranmine()-1.0)
               dist=sqrt(xdif*xdif+ydif*ydif+zdif*zdif)
               iref=int(dist*nz)+1
               if(iref.lt.1) iref=1
               if(iref.gt.nz) iref=nz
               gofz(iref,idref2)=gofz(iref,idref2)+1
               ivoidsum=ivoidsum+1


c Get the coordinates of this pixel - convert these into fractional coordinates

               ilv1=(iv1-1)/npixtot2
               ikv1=(iv1-1-ilv1*npixtot2)/npixtot
               ihv1=iv1-1-npixtot*ikv1-npixtot2*ilv1
               xv1=real(ihv1)/real(npixtot)
               yv1=real(ikv1)/real(npixtot)
               zv1=real(ilv1)/real(npixtot)

               iscount=0
               iv2=0
               do while(iv2.lt.npixall)

                  iv2=iv2+1
                  if(voidpix(iv2).eq.0.and.iv1.ne.iv2) then

                     iscount=iscount+1

c Get the coordinates of this pixel - convert these into fractional coordinates

c                    ivt=int(ranmine()*npixall)+1
                     ivt=iv2
                     ilv2=(ivt-1)/npixtot2
                     ikv2=(ivt-1-ilv2*npixtot2)/npixtot
                     ihv2=ivt-1-npixtot*ikv2-npixtot2*ilv2
c     write(6,*) ivt,ilv2,ikv2,ihv2
                     xv2=real(ihv2)/real(npixtot)
                     yv2=real(ikv2)/real(npixtot)
                     zv2=real(ilv2)/real(npixtot)

c Calculate separation of pixels, applying the minimum image convention

                     xdif=xv2-xv1
                     ydif=yv2-yv1
                     zdif=zv2-zv1

c Add the broadening component

                     xdif=xdif+pixhwidth*(2.0*ranmine()-1.0)
                     ydif=ydif+pixhwidth*(2.0*ranmine()-1.0)
                     zdif=zdif+pixhwidth*(2.0*ranmine()-1.0)
c     xdif=ranmine()-0.5
c     ydif=ranmine()-0.5
c     zdif=ranmine()-0.5
                     call minimagec(xdif,ydif,zdif)
                     dist=sqrt(xdif*xdif+ydif*ydif+zdif*zdif)
c                       write(6,*) ihv1,ikv1,ilv1,ihv2,ikv2,ilv2,npixtot
                     if(xdif.gt.0.5.or.ydif.gt.0.5.or.zdif.gt.0.5) then
                        write(6,*) ihv1,ikv1,ilv1,ihv2,ikv2,ilv2,npixtot
                     else
                        iref=int(dist*nz)+1
                        if(iref.lt.1) iref=1
                        if(iref.gt.nz) iref=nz
                        gofz(iref,idref2)=gofz(iref,idref2)+1
                        if(voidnum(iv1).ne.voidnum(iv2)) then
                        gofz(iref,idref3)=gofz(iref,idref3)+1
                        endif
                        ivoidsum=ivoidsum+1
                     endif

                  endif
               end do
      write(6,*) ivcount,iscount
            endif
         end do

c Normalise this g(r) to 4*pi*rho*r**2delta-r

         rhovoid=real(nvcount*iscount)
         rvstep=1.0/real(nz)
         pi=3.141592653
         fact=3.0/(4.0*pi*rhovoid)
      write(6,*) nvcount,rhovoid,rvstep,pi,fact
         fracfac1=1.0-fract2(id)
         fracfac2=fract2(id)
         fracfac3=fracfac1*fracfac1
         rv1=0.0
         rv13=rv1*(rv1*rv1)
         do iz=1,nz

            rv2=iz*rvstep
            rv23=rv2*(rv2*rv2)
            rvfact=fact/(rv23-rv13)
            gofz(iz,idref2)=gofz(iz,idref2)*rvfact
            gofz(iz,idref3)=gofz(iz,idref3)*rvfact

cc Subtract the 'uninteresting' spherical g(r) term
c
c           rstar=real(iz)/(real(nz)*radactual)
c     write(6,*) iz,rstar,spheregr(rstar),fracfac1,fracfac2,fracfac3
c           if(fract2(id).gt.0.0005) then
c              gofz(iz,idref3)=spheregr(rstar)/frac2(id)
c           else
c              gofz(iz,idref3)=0.0
c           endif
            rv1=rv2
            rv13=rv23

         end do

c Repeat for other distributions

      end do
      return
      end

      function ranmine()

c Get a random number

      include 'ran1com.inc'

      ranmine=ran1()

      return

      end


      subroutine occupypix(npix,npixsphere,npixfound,ih,ik,il
     1,voidpix,mpix3)

c Sets the occupation of pixels within radius npixsphere from the point (ih,ik,il)

c ih,ik and il are specified relative to the centre of the box, whereas pixels are
c counted from a corner. If the calculated pixel position is outside the box then
c the image position is used.

c nvoidpix keeps track of the number of unoccupied pixels

      integer*1 voidpix(mpix3)
      integer*4 npix,npixsphere,npixsphere2,npixfound
      integer*4 ih,ik,il,ils,iks,ihs

c Square no. of pixels in sphere to save square roots later

      npixsphere2=npixsphere*npixsphere
      ils=0
      iks=0
      ihs=0

c Count the number of pixels found in the sphere

      npixfound=0

      do while(ils.le.npixsphere)
         ils2=ils*ils
         iks2pils2=ils2
         iks=0
         do while(iks.le.npixsphere)
            iks2=iks*iks
            iks2pils2=ils2+iks2
            ihs=0
            do while(ihs.le.npixsphere)

               ihs2=ihs*ihs
               ihs2piks2pils2=iks2pils2+ihs2

c Test for inside sphere

               if(ihs2piks2pils2.le.npixsphere2) then

c Found an occupied pixel

                  ihabs=ih+ihs
                  ikabs=ik+iks
                  ilabs=il+ils
                  call markpixel(npix,ihabs,ikabs,ilabs
     1,npixfound,voidpix,mpix3)

c Reverse the sign of ihs, if appropriate

                  if(ihs.gt.0) then
                     ihabs=ih-ihs
                     call markpixel(npix,ihabs,ikabs,ilabs
     1,npixfound,voidpix,mpix3)
                  endif

c Repeat for the reversed sign of iks, if appropriate

                  if(iks.gt.0) then

                     ihabs=ih+ihs
                     ikabs=ik-iks
                     ilabs=il+ils
                     call markpixel(npix,ihabs,ikabs,ilabs
     1,npixfound,voidpix,mpix3)

                     if(ihs.gt.0) then
                        ihabs=ih-ihs
                        call markpixel(npix,ihabs,ikabs,ilabs
     1,npixfound,voidpix,mpix3)

                     endif

                  endif

c Repeat for reversed sign of ils, if appropriate.

                  if(ils.gt.0) then

                     ihabs=ih+ihs
                     ikabs=ik+iks
                     ilabs=il-ils
                     call markpixel(npix,ihabs,ikabs,ilabs
     1,npixfound,voidpix,mpix3)

                     if(ihs.gt.0) then
                        ihabs=ih-ihs
                        call markpixel(npix,ihabs,ikabs,ilabs
     1,npixfound,voidpix,mpix3)
                     endif

                     if(iks.gt.0) then

                        ihabs=ih+ihs
                        ikabs=ik-iks
                        call markpixel(npix,ihabs,ikabs,ilabs
     1,npixfound,voidpix,mpix3)

                        if(ihs.gt.0) then
                           ihabs=ih-ihs
                           call markpixel(npix,ihabs,ikabs,ilabs
     1,npixfound,voidpix,mpix3)
                        endif
                     endif
                  endif
               endif
               ihs=ihs+1
            end do
            iks=iks+1
         end do
         ils=ils+1
      end do
      return
      end

      subroutine markpixel(npix,ihabs,ikabs,ilabs
     1,npixfound,voidpix,mpix3)

c Routine to mark a pixel as occupied

      integer*1 voidpix(mpix3)
      integer*4 npix,offset,npixtot,ihabs,ikabs,ilabs,npixref
     1,ihref,ikref,ilref,npixfound

      offset=npix
      npixtot=2*npix+1

c Apply periodic boundary conditions, and shift coordinates of pixel so that
c origin is at one corner of the box instead of at the centre.

      ihref=ihabs-npixtot*nint(real(ihabs)/real(npixtot))+offset
      ikref=ikabs-npixtot*nint(real(ikabs)/real(npixtot))+offset
      ilref=ilabs-npixtot*nint(real(ilabs)/real(npixtot))+offset

c Generate sequence number of pixel

      npixref=1+ihref+npixtot*(ikref+npixtot*ilref)
      if(npixref.le.mpix3) then
c     write(6,*) ihabs,ikabs,ilabs,ihref,ikref,ilref,npixref
         voidpix(npixref)=1

c Increment pixel counter

         npixfound=npixfound+1

      else

         write(6,100) ihabs,ikabs,ilabs,ihref,ikref,ilref,npixref
100   format('markpixel> Error generating pixel reference numbers:'
     1,/6(1x,i4),1x,i10)
         stop

      endif
      return
      end

      subroutine search_voids(nvoids,npix,npixall
     1,voidpix,voidnum,mpix3)

c Set maximum number of voids to maximum word size for voidnum

      parameter (maxvoid=2**15-1)
      integer*1 voidpix(mpix3)
      integer*2 voidnum(mpix3)
      integer*4 voidlist(27)
      integer*4 nvoid,npix,npixall
      integer*4 npixtot,nvoidlist,voidtally(maxvoid)

      npixtot=2*npix+1
      nvoid=nvoids

c Step through the unoccupied pixels

      do ipix=1,npixall

         if(voidpix(ipix).eq.0) then

c Get the list of unoccupied pixel numbers which are adjacent to ipix

            call get_void_list(ipix,npix,nvoidlist
     1,voidlist,voidpix,mpix3)

c Search for the lowest non-zero void number within this list

            ilowest=maxvoid
            do it=1,nvoidlist

               itest=voidnum(voidlist(it))
               if(itest.lt.ilowest.and.itest.gt.0) ilowest=itest

            end do

c If ilowest is maxvoid, it means none of the unoccupied pixels have been
c given a void number so create one.

            if(ilowest.eq.maxvoid.and.nvoid.lt.maxvoid) then

               nvoid=nvoid+1
               ilowest=nvoid

            else if(nvoid.ge.maxvoid) then

               write(6,100)
100   format(//'search_voids> Unable to allocate any more void numbers!'
     1//'search_voids> Try using a larger void distance?')
               stop

            endif

c Set the void number of all the unoccupied pixels to this (lowest)
c void number.

            do it=1,nvoidlist

               voidnum(voidlist(it))=ilowest

            end do
         endif
      end do

c Now count the total number of voids

      nvoids=0

c voidtally will store the list of void numbers

c Step through the unoccupied pixels

      do ipix=1,npixall

         if(voidpix(ipix).eq.0) then

            ifind=0
            itest=0
            do while(itest.lt.nvoids.and.ifind.eq.0)

               itest=itest+1
               if(voidnum(ipix).eq.voidtally(itest)) ifind=voidnum(ipix)

            end do

            if(ifind.eq.0) then

               nvoids=nvoids+1
               voidtally(nvoids)=voidnum(ipix)

            endif
         endif
      end do

      return
      end

      subroutine get_void_list(ipix,npix,nvoidlist
     1,voidlist,voidpix,mpix3)

c Makes a list of the unoccupied pixels adjacent to ipix 

      integer*1 voidpix(mpix3)
c list of local void pixel numbers
      integer*4 voidlist(27)
      integer*4 ipix,npix,nvoidlist
      integer*4 offset,npixtot,npixtot2
      integer*4 ihc,ikc,ilc,ih,ik,il,pixref

      offset=npix
      npixtot=2*npix+1
      npixtot2=npixtot*npixtot
      nvoidlist=0

c Split ipix into its component indices relative to corner of box

      ilc=(ipix-1)/npixtot2
      ikc=(ipix-1-ilc*npixtot2)/npixtot
      ihc=ipix-1-npixtot*ikc-npixtot2*ilc

c Shift these to centre of box

      ilc=ilc-offset
      ikc=ikc-offset
      ihc=ihc-offset
c     write(6,*) ipix,ihc,ikc,ilc

c Step through the 27 pixels around ipix (including ipix) and determine the 
c smallest, non-zero void number

      ilt=ilc-2
      do while(ilt.lt.ilc+1)

         ilt=ilt+1

c Apply periodic boundary condition to il

         il=ilt-npixtot*nint(real(ilt)/real(npixtot))
         iladd=(il+offset)*npixtot2

         ikt=ikc-2
         do while(ikt.lt.ikc+1)

            ikt=ikt+1
            ik=ikt-npixtot*nint(real(ikt)/real(npixtot))
            ikadd=iladd+(ik+offset)*npixtot

            iht=ihc-2
            do while(iht.lt.ihc+1)

               iht=iht+1
               ih=iht-npixtot*nint(real(iht)/real(npixtot))

               pixref=1+ih+offset+ikadd
c     write(6,*) ih,ik,il,pixref
c Test this pixel for occupation and void number

               if(voidpix(pixref).eq.0) then

                  nvoidlist=nvoidlist+1
                  voidlist(nvoidlist)=pixref

               endif
            end do
         end do
      end do
      return
      end

      function spheregr(rstar)
      parameter (cons1=0.75,cons2=0.0625)
      real spheregr

c Gets the form factor (g(r)) for a solid uniform sphere

      spheregr=0.0
      if(rstar.lt.2.0) then

         spheregr=1-rstar*(cons1-rstar*rstar*cons2)

      endif
      return
      end
  

      
