      module sdfcube_mod

      integer                                      :: natosave,lenfname_cube
      real, dimension(:), allocatable              :: xasave,yasave,zasave
      character(len=3), dimension(:), allocatable  :: atsymbolsave
      real, dimension(:,:,:), allocatable          :: rhocube
      character(len=256)                           :: fname_cube
      integer, dimension(:), allocatable           :: nm !First or second molecule
      real, dimension(:), allocatable              :: phim,cthetm,chim,xcent,ycent,zcent
!Logical variables to perform a test calculation and to output average or single molecules for plotting purposes
      logical                                      :: performtest,outputaveragemol

      end module sdfcube_mod

      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_sdfcube()

      stop
      end
      subroutine run_sdfcube()

      use sdfcube_mod
c
c performs calculation of spherical harmonic coefficients from .ato file
c
c uses arithmetic centre of molecule rather than a particular atom to define
c molecular centre. 
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'sdfcube.inc'
      include 'sharm_sdf.inc'
c
c setup the SHARM symbols
c
      call setup_symbols_sdfcube()
c
c now setup the input file
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
c if an error occurred on read, exit now
c
      if(ierrout.ne.0) return
! Set up the name of the file where the CUBE data will be stored
      lenfname_cube=index(fnameinp,'.SDFCUBE.dat')
      if(lenfname_cube.le.0) return ! Apparently the wrong input file specified!
      fname_cube=fnameinp(1:lenfname_cube-1)//'.SDFCUBE.cube'
      lenfname_cube=len_trim(fname_cube)
      iref=1
c
c set up number of pixels
c
      iref=iref+1
      read(symbol_value(iref),*) nz
c
c set up rmax of calculation
c
      iref=iref+1
      read(symbol_value(iref),*) zmax
c
c set up number of configurations already stored
c
      iref=iref+1
      read(symbol_value(iref),*) nsumt
c
c for molecule 1 read number and type of atoms used to define the molecular
c centre
c
      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
      ncent1=0
c
c skip through the atom types to find out the type number of the specified atom
c
      do ic=1,nwords
         it=0
         ifound=0
         do while(it.lt.ntype.and.ifound.eq.0)
            it=it+1
            if(type(it).eq.symbol_value(iref)(ncf(ic):ncl(ic))) ifound=it
         end do
         if(ifound.le.0) then
            write(6,201) symbol_value(iref)(ncf(ic):ncl(ic))
     *,fname_ato(1:ilenfname_ato)
201   format(/'sharm> Atom type ',a,' does not exist in .ato file ',a)
            return
         else
            if(ncent1.lt.mcentre) then
               ncent1=ncent1+1
               lctyp1(ncent1)=ifound
            endif
         endif
      end do
c
c for molecule 1 read primary axis type (1=X, 2=Y, 3=Z), and
c the number and atom numbers of atoms used to define the molecule's nzt1
c axis
c
      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
      nzt1=symbol_value(iref)(ncf(1):ncl(1))
      ncz1=nwords-1
      if(ncz1.gt.mcentre) ncz1=mcentre
      do i=1,ncz1
         iread=i+1
         read(symbol_value(iref)(ncf(iread):ncl(iread)),*) lztyp1(i)
      end do
c
c for molecule 1 read secondary axis type (1=X, 2=Y, 3=Z), and
c the number and atom numbers of atoms used to define the
c molecule's nyt1 axis
c
      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
      nyt1=symbol_value(iref)(ncf(1):ncl(1))
      ncy1=nwords-1
      if(ncy1.gt.mcentre) ncy1=mcentre
      do i=1,ncy1
         iread=i+1
         read(symbol_value(iref)(ncf(iread):ncl(iread)),*) lytyp1(i)
      end do

c for molecule 2 read number and type of atoms used to define the molecular
c centre. Also check whether this molecule is identical to the first

      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
      ncent2=0
c
c skip through the atom types to find out the type number of the specified atom
c
      do ic=1,nwords
         it=0
         ifound=0
         do while(it.lt.ntype.and.ifound.eq.0)
            it=it+1
            if(type(it).eq.symbol_value(iref)(ncf(ic):ncl(ic))) ifound=it
         end do
         if(ifound.eq.0) then
            write(6,201) symbol_value(iref)(ncf(ic):ncl(ic))
     *,fname_ato(1:ilenfname_ato)
            return
         else
            if(ncent2.lt.mcentre) then
               ncent2=ncent2+1
               lctyp2(ncent2)=ifound
            endif
         endif
      end do

c test for identity of central and second molecules

      ldistinct=0
      if(ncent2.eq.ncent1) then
         ldistinct=0
         do it=1,ncent1
            if(lctyp2(it).ne.lctyp1(it)) ldistinct=1
         end do
      else
         ldistinct=1
      endif

      write(6,*) ncent1,ncent2,ldistinct



      if(ldistinct.ne.0) then

c if molecules are different, then repeat this for the second molecule

c for molecule 2 read primary axis type (1=X, 2=Y, 3=Z), and
c the number and atom numbers of atoms used to define the molecule's nzt1
c axis

         iref=iref+1
         call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
         nzt2=symbol_value(iref)(ncf(1):ncl(1))
         ncz2=nwords-1
         if(ncz2.gt.mcentre) ncz2=mcentre
         do i=1,ncz2
            iread=i+1
            read(symbol_value(iref)(ncf(iread):ncl(iread)),*) lztyp2(i)
         end do

c for molecule 2 read secondary axis type (1=X, 2=Y, 3=Z), and
c the number and atom numbers of atoms used to define the
c molecule's nyt1 axis

         iref=iref+1
         call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
         nyt2=symbol_value(iref)(ncf(1):ncl(1))
         ncy2=nwords-1
         if(ncy2.gt.mcentre) ncy2=mcentre
         do i=1,ncy2
            iread=i+1
            read(symbol_value(iref)(ncf(iread):ncl(iread)),*) lytyp2(i)
         end do
      else

c otherwise set all the values to the same as the central molecule

         ncent2=ncent1
         nzt2=nzt1
         nyt2=nyt1
         if(ncent2.gt.0) then
            do ic=1,ncent1
               lctyp2(ic)=lctyp1(ic)
            end do
         end if
         ncz2=ncz1
         do i=1,ncz1
            lztyp2(i)=lztyp1(i)
         end do
         ncy2=ncy1
         do i=1,iabs(ncy1)
            lytyp2(i)=lytyp1(i)
         end do
      endif
      performtest=.false.
      iref=iref+1
      outputaveragemol=symbol_value(iref).eq.'y'
c     write(6,*) ncent1,ncent2,ncz1,ncy1,ncz2,ncy2
c
c only perform the calculation when nsum is zero or a multiple of nmult
c
      nmult=5
      nsum=nsumt/nmult
      nt=nsum*nmult
      if(nt.eq.nsumt) then
c
c normalise atomic coordinates to box size
c
            call normalise_coordinates()
c
c set up step size and z values
c
            zstep=zmax/real(nz)
            call initcube(nin)
            call sdfcube(nin,fnameinp)
            call outcube(nin)
      endif
      nsumt=nsumt+1
c
c write input file
c
c
c make sure all current values are stored
c
      write(symbol_value(4),*) nsumt
      call save_files_others(14)
      close(14)
      return
      end

      subroutine initcube(nin)

      use sdfcube_mod

      include 'sdfcube.inc'
      real, parameter                              :: pi=3.141592653,bohr=0.52917721092
      integer                                      :: nin,ierr

! If nsumt.gt.0 we need to read in the current values from the specified cube file
      if(nsumt.gt.0) then !Make sure our step sizes are the same as in the existing file
         open(nin,file=fname_cube,status='old',iostat=ierr)
         if(ierr.eq.0) then 
! Read two blank lines
            read(nin,*)
            read(nin,*)
! Number of atoms in central molecule
            read(nin,*) natosave
! Number of pixels along each of x,y,z - we will assume this is the same number for each axis
            read(nin,*) nz0,zstep !zstep in Bohr units
            zstep=zstep*bohr !Convert to Angstroms
            close(nin)
            nz=(nz0-1)/2 
            nx=nz
            ny=nz
            nz0=2*nz+1 !Always an odd number of pixels along each axis (nz0=2*nz+1)
            nx0=nz0
            ny0=nz0
            xstep=zstep
            ystep=zstep
         else !If file could not be opened, reset the distributions 
            nsumt=0
            nx=nz
            ny=nz
            nz0=2*nz+1 !Always an odd number of pixels along each axis (nz0=2*nz+1)
            nx0=nz0
            ny0=nz0
            xstep=zstep
            ystep=zstep
         end if
      else
         nsumt=0
         nx=nz
         ny=nz
         nz0=2*nz+1 !Always an odd number of pixels along each axis (nz0=2*nz+1)
         nx0=nz0
         ny0=nz0
         xstep=zstep
         ystep=zstep
      end if
! Allocate the density array
      if(allocated(rhocube)) deallocate(rhocube)
      allocate(rhocube(nz0,nz0,nz0))
! Set all values to zero
      do ix=1,nx0
         do iy=1,ny0
            do iz=1,nz0
               rhocube(ix,iy,iz)=0.0
            end do
         end do
      end do
      return
      end
            
      subroutine outcube(nout)

      use sdfcube_mod

      include 'sdfcube.inc'
      real, parameter                              :: pi=3.141592653,bohr=0.52917721092
      integer                                      :: nout,natosavet,ierr,get_atomicnumber
      real, dimension(10)                          :: rhocubet

! If nsumt.gt.0 we need to read in the current values from the specified cube file
      if(nsumt.gt.0) then !Make sure our step sizes are the same as in the existing file
         open(nout,file=fname_cube,status='old',iostat=ierr)
         if(ierr.eq.0) then 
! Read two blank lines
            read(nout,*)
            read(nout,*)
! Number of atoms in central molecule
            read(nout,*) natosavet
! Number of pixels along each of x,y,z - we will assume this is the same number for each axis
            read(nout,*) 
            read(nout,*) 
            read(nout,*)
! Coordinates of central molecule - these will be rewritten so ignore current values
            do i=1,natosavet
               read(nout,*) 
            end do
! Get the existing density distribution and add in the new values
            fc=1.0/real(nsumt)
            fc1=1.0-fc1
            do ix=1,nx0
               do iy=1,nz0
                  do iz=1,nz0,10
                     itmax=min(9,nz0-iz)
                     read(nout,*) (rhocubet(1+it),it=0,itmax)
                     do it=0,itmax
                        izref=iz+it
                        rhocube(ix,iy,izref)=fc1*rhocubet(1+it)+fc*rhocube(ix,iy,izref)
                     end do
                  end do
               end do
            end do
            close(nout)
         else !If file could not be opened, reset the distributions 
            nsumt=0
         end if
      end if
! Write out the revised values
      open(nout,file=fname_cube,status='unknown',iostat=ierr)
      if(ierr.eq.0) then
         write(nout,*) '# Cube file output by sdfcube'
         write(nout,*) 'Outer loop x, middle loop y, inner loop z, filname: '//fname_cube(1:len_trim(fname_cube))
         write(nout,*) natosave,-nx*xstep/bohr,-ny*ystep/bohr,-nz*zstep/bohr
         write(nout,*) nx0, xstep/bohr, 0.0, 0.0
         write(nout,*) ny0, 0.0, ystep/bohr, 0.0
         write(nout,*) nz0, 0.0, 0.0, zstep/bohr
         do is=1,natosave
            iatomicno=get_atomicnumber(atsymbolsave(is))
            write(nout,*) iatomicno,real(iatomicno),xasave(is)/bohr,yasave(is)/bohr,zasave(is)/bohr
         end do
         do ix=1,nx0
            do iy=1,ny0
               do iz=1,nz0,10
                  itmax=min(9,nz0-iz)
                  write(nout,'(10(1x,e11.4))') (rhocube(ix,iy,iz+it),it=0,itmax)
               end do
            end do
         end do
         close(nout)
      endif
      return
      end

      subroutine sdfcube(nout,fname_inp)

      use reallocation_routines
      use sdfcube_mod

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      include 'sdfcube.inc'
      include 'sharm_sdf.inc'

c     integer*4 rate,max,time1,time2,time3,time4
      real amass,amasssum

      logical                                      :: centralmolecule,found,foundatom
      character(len=*)                             :: fname_inp
      character(len=256)                           :: fname_xyz,textline !Filename for xyz data, line of text
      integer                                      :: nout !Unit number on which to write xyz file
      real                                         :: xt,yt,zt,xr,yr,zr,xdif,ydif,zdif,xdifr,ydifr,zdifr,xstepd,ystepd,zstepd
      integer                                      :: ncentral,natocentral
      real, dimension(:), allocatable              :: xacentral,yacentral,zacentral
      character(len=3), dimension(:), allocatable  :: atsymbolcentral,attypecentral
      logical, dimension(:), allocatable           :: centraldefinition
      real, dimension(:), allocatable              :: xat,yat,zat
 
      pi=3.141592653
! Arrays to store Euler angles
      if(allocated(nm)) deallocate(nm,phim,cthetm,chim,xcent,ycent,zcent)
      allocate(nm(nmol),phim(nmol),cthetm(nmol),chim(nmol),xcent(nmol),ycent(nmol),zcent(nmol))
c
c set up coordinates of Laboratory X, Y and Z axes
c
      alab(1)=1.0
      alab(2)=0.0
      alab(3)=0.0
      blab(1)=0.0
      blab(2)=1.0
      blab(3)=0.0
      clab(1)=0.0
      clab(2)=0.0
      clab(3)=1.0
! Renormalised step sizes
      xstepd=xstep/dmax
      ystepd=ystep/dmax
      zstepd=zstep/dmax
c
c set up coordinates of molecular centres and Euler angles for each molecule 
c in laboratory coordinate frame
c
      nsmmo1=0
      nsmmo2=0
      do imol1=1,nmol
c define origin of this molecule
            xo1=0.0
            yo1=0.0
            zo1=0.0
            amasssum=0.0
            ncref=0
            nm(imol1)=0
c
c test to see if this is molecule 1 or molecule 2. Need to check that all the specified types occur - not just one
c
            found=.true.
            ic=0
            do while(ic.lt.ncent1.and.found)
               ic=ic+1
               foundatom=.false.
               ia1=naf(imol1)-1
               do while(ia1.lt.nal(imol1).and..not.foundatom)
                  ia1=ia1+1
                  foundatom=itype(ia1).eq.lctyp1(ic)
               end do
               if(ic.eq.1) then
                  found=foundatom
               else
                  found=found.and.foundatom
               end if
            end do
            if(found) then
               nm(imol1)=1
            else !molecule 1 not found, try molecule 2
               found=.true.
               ic=0
               do while(ic.lt.ncent2.and.found)
                  ic=ic+1
                  foundatom=.false.
                  ia1=naf(imol1)-1
                  do while(ia1.lt.nal(imol1).and..not.foundatom)
                     ia1=ia1+1
                     foundatom=itype(ia1).eq.lctyp2(ic)
                  end do
                  if(ic.eq.1) then
                     found=foundatom
                  else
                     found=found.and.foundatom
                  end if
               end do
               if(found) then
                  nm(imol1)=2
               endif
            endif
            if(nm(imol1).eq.1) then
              nsmmo1=nsmmo1+1
            endif
            if(nm(imol1).eq.2) then
              nsmmo2=nsmmo2+1
            endif
            if(nm(imol1).eq.1.and.ncent1.gt.0) then
              do ia1=naf(imol1),nal(imol1)
                do ic=1,ncent1
                  if(itype(ia1).eq.lctyp1(ic)) then
                  amass=amassa(itype(ia1))
                  xo1=xo1+xa(ia1)*amass
                  yo1=yo1+ya(ia1)*amass
                  zo1=zo1+za(ia1)*amass
                  amasssum=amasssum+amass
                  endif
                end do
              end do
            else if(nm(imol1).eq.2.and.ncent2.gt.0) then
              do ia1=naf(imol1),nal(imol1)
                do ic=1,ncent2
                  if(itype(ia1).eq.lctyp2(ic)) then
                  amass=amassa(itype(ia1))
                  xo1=xo1+xa(ia1)*amass
                  yo1=yo1+ya(ia1)*amass
                  zo1=zo1+za(ia1)*amass
                  amasssum=amasssum+amass
                  endif
                end do
              end do
            else 
               do ia1=naf(imol1),nal(imol1)
                  amass=amassa(itype(ia1))
                  xo1=xo1+xa(ia1)*amass
                  yo1=yo1+ya(ia1)*amass
                  zo1=zo1+za(ia1)*amass
                  amasssum=amasssum+amass
               end do
            endif
            if(amasssum.gt.0.0.and.nm(imol1).ne.0) then
               xo1=xo1/amasssum
               yo1=yo1/amasssum
               zo1=zo1/amasssum
c     write(6,991) xo1,yo1,zo1
               xcent(imol1)=xadd(imol1)+xo1
               ycent(imol1)=yadd(imol1)+yo1
               zcent(imol1)=zadd(imol1)+zo1
! Form atomic coordinates relative to this centre
               call reallocate1d_r(xat,nato(imol1))
               call reallocate1d_r(yat,nato(imol1))
               call reallocate1d_r(zat,nato(imol1))
               do ia=1,nato(imol1)
                  iaref=ia+naf(imol1)-1
                  xat(ia)=xa(iaref)-xo1
                  yat(ia)=ya(iaref)-yo1
                  zat(ia)=za(iaref)-zo1
               end do
               if(nm(imol1).eq.1) then
                  call form_eulerangles(imol1,xat,yat,zat,nzt1,nyt1,ncz1,ncy1,lztyp1,lytyp1,alab,blab,clab,phim(imol1),cthetm(imol1),chim(imol1))
               else if(nm(imol1).eq.2) then
                  call form_eulerangles(imol1,xat,yat,zat,nzt2,nyt2,ncz2,ncy2,lztyp2,lytyp2,alab,blab,clab,phim(imol1),cthetm(imol1),chim(imol1))
               end if
            else
               nm(imol1)=0
            endif
499   format(1x,i5,4(1x,e12.5))
      end do
      write(6,333) nsmmo1,nsmmo2
333   format(1x,'Number of type 1 molecules = ',i10
     1/1x, 'Number of type 2 molecules = ',i10)
! Write out coordinates of centralmolecule. We will use the average coordinates of all these molecules.
      im=0
      ncentral=0
      centralmolecule=.false.
      do while(im.lt.nmol)
         im=im+1
         centralmolecule=nm(im).eq.1
         if(centralmolecule) then !Write out coordinates of this molecule
            if(ncentral.eq.0.or.outputaveragemol) then
               ncentral=ncentral+1
               if(ncentral.eq.1) then !First time allocate the arrays to store summed coordinates
                  if(allocated(xacentral)) deallocate(xacentral,yacentral,zacentral,atsymbolcentral,attypecentral,centraldefinition)
                  natocentral=nato(im)
                  allocate(xacentral(natocentral),yacentral(natocentral),zacentral(natocentral)
     *,atsymbolcentral(natocentral),attypecentral(natocentral),centraldefinition(natocentral))
               end if
! Shift coordinates to centre used by SHARM/SDF
               xshift=xadd(im)-xcent(im)
               yshift=yadd(im)-ycent(im)
               zshift=zadd(im)-zcent(im)
               do ia=naf(im),nal(im)
! Temporary coordinates
                  xt=xa(ia)+xshift
                  yt=ya(ia)+yshift
                  zt=za(ia)+zshift
! Need to rotate these to align with the laboratory axes - performed in reverse order about the laboratory axes
                  call rotateeuler(xt,yt,zt,-chim(im),-acos(cthetm(im)),-phim(im),xr,yr,zr)
                  iaref=ia-naf(im)+1
                  if(ncentral.eq.1) then
                     xacentral(iaref)=xr
                     yacentral(iaref)=yr
                     zacentral(iaref)=zr
                     atsymbolcentral(iaref)=atsymbol(itype(ia))
                     attypecentral(iaref)=type(itype(ia))
! Check whether this atom is part of the central molecule definition
                     found=.false.
                       ic=0
                     do while(.not.found.and.ic.lt.ncent1)
                        ic=ic+1
                        found=itype(ia).eq.lctyp1(ic)
                     end do
                     iz=0
                     do while (.not.found.and.iz.lt.ncz1)
                        iz=iz+1
                        found=iaref.eq.lztyp1(iz)
                     end do
                     iy=0
                     do while (.not.found.and.iy.lt.ncy1)
                        iy=iy+1
                        found=iaref.eq.lytyp1(iy)
                     end do
                     centraldefinition(iaref)=found
                  else
                     xacentral(iaref)=xacentral(iaref)+xr
                     yacentral(iaref)=yacentral(iaref)+yr
                     zacentral(iaref)=zacentral(iaref)+zr
                  end if
               end do
            end if
         end if
      end do
! Form average coordinates for central molecule
      if(ncentral.gt.0) then
         do ia=1,natocentral
            xacentral(ia)=xacentral(ia)*dmax/real(ncentral)
            yacentral(ia)=yacentral(ia)*dmax/real(ncentral)
            zacentral(ia)=zacentral(ia)*dmax/real(ncentral)
         end do
! Save these coordinates to output to the cube file
         if(allocated(xasave)) deallocate(xasave,yasave,zasave,atsymbolsave)
         natosave=natocentral
         allocate(xasave(natosave),yasave(natosave),zasave(natosave),atsymbolsave(natosave))
         do ia=1,natocentral
            xasave(ia)=xacentral(ia)   
            yasave(ia)=yacentral(ia)   
            zasave(ia)=zacentral(ia)
            atsymbolsave(ia)=atsymbolcentral(ia)
         end do   
         fname_xyz=fname_inp(1:indexlastdot(fname_inp)-1)//'-first.xyz'
         open(nout,file=fname_xyz,status='unknown',iostat=ierr)
         if(ierr.eq.0) then
! Write number of atoms in molecule
            write(nout,'(i6)') natocentral
! Write information line - not used for plotting
            write(nout,'(a)') 'Central molecule coordinate file: '//fname_xyz(1:len_trim(fname_xyz))
            do ia=1,natocentral
               write(textline,'(a3,1x,f10.5,1x,f10.5,1x,f10.5,1x,a3)') atsymbolcentral(ia),xacentral(ia),yacentral(ia),zacentral(ia),attypecentral(ia)
               ! Append a comment if this atom is part of the central molecule definition
               if(centraldefinition(ia)) then
                  textline=textline(1:len_trim(textline))//'  #CENTRALMOLECULE'
               end if
               write(nout,'(a)') textline(1:len_trim(textline))
            end do
            close(nout)
         end if
         if(ldistinct.eq.0) then !Also write this out as the second molecule
            fname_xyz=fname_inp(1:indexlastdot(fname_inp)-1)//'-second.xyz'
            open(nout,file=fname_xyz,status='unknown',iostat=ierr)
            if(ierr.eq.0) then
! Write number of atoms in molecule
               write(nout,'(i6)') natocentral
! Write information line - not used for plotting
               write(nout,'(a)') 'Central molecule coordinate file: '//fname_xyz(1:len_trim(fname_xyz))
               do ia=1,natocentral
                  write(textline,'(a3,1x,f10.5,1x,f10.5,1x,f10.5,1x,a3)') atsymbolcentral(ia),xacentral(ia),yacentral(ia),zacentral(ia),attypecentral(ia)
                  ! Append a comment if this atom is part of the central molecule definition
                  if(centraldefinition(ia)) then
                     textline=textline(1:len_trim(textline))//'  #CENTRALMOLECULE'
                  end if
                  write(nout,'(a)') textline(1:len_trim(textline))
               end do
               close(nout)
            end if
         end if
      end if
! Now repeat this for second molecules if they are different from first molecules
! Write out coordinates of centralmolecule. We will use the average coordinates of all these molecules.
      im=0
      ncentral=0
      centralmolecule=.false.
      do while(im.lt.nmol)
         im=im+1
         centralmolecule=nm(im).eq.2
         if(centralmolecule) then !Write out coordinates of this molecule
            if(ncentral.eq.0.or.outputaveragemol) then
               ncentral=ncentral+1
               if(ncentral.eq.1) then !First time allocate the arrays to store summed coordinates
                  if(allocated(xacentral)) deallocate(xacentral,yacentral,zacentral,atsymbolcentral,attypecentral,centraldefinition)
                  natocentral=nato(im)
                  allocate(xacentral(natocentral),yacentral(natocentral),zacentral(natocentral)
     *,atsymbolcentral(natocentral),attypecentral(natocentral),centraldefinition(natocentral))
               end if
! Shift coordinates to centre used by SHARM/SDF
               xshift=xadd(im)-xcent(im)
               yshift=yadd(im)-ycent(im)
               zshift=zadd(im)-zcent(im)
               do ia=naf(im),nal(im)
! Temporary coordinates
                  xt=xa(ia)+xshift
                  yt=ya(ia)+yshift
                  zt=za(ia)+zshift
! Need to rotate these to align with the laboratory axes - performed in reverse order about the laboratory axes
                  call rotateeuler(xt,yt,zt,-chim(im),-acos(cthetm(im)),-phim(im),xr,yr,zr)
                  iaref=ia-naf(im)+1
                  if(ncentral.eq.1) then
                     xacentral(iaref)=xr
                     yacentral(iaref)=yr
                     zacentral(iaref)=zr
                     atsymbolcentral(iaref)=atsymbol(itype(ia))
                     attypecentral(iaref)=type(itype(ia))
! Check whether this atom is part of the central molecule definition
                     found=.false.
                        ic=0
                     do while(.not.found.and.ic.lt.ncent1)
                        ic=ic+1
                        found=itype(ia).eq.lctyp1(ic)
                     end do
                     iz=0
                     do while (.not.found.and.iz.lt.ncz1)
                        iz=iz+1
                        found=iaref.eq.lztyp1(iz)
                     end do
                     iy=0
                     do while (.not.found.and.iy.lt.ncy1)
                        iy=iy+1
                        found=iaref.eq.lytyp1(iy)
                     end do
                     centraldefinition(iaref)=found
                  else
                     xacentral(iaref)=xacentral(iaref)+xr
                     yacentral(iaref)=yacentral(iaref)+yr
                     zacentral(iaref)=zacentral(iaref)+zr
                  end if
               end do
            end if
         end if
      end do
! Form average coordinates for central molecule
      if(ncentral.gt.0) then
         do ia=1,natocentral
            xacentral(ia)=xacentral(ia)*dmax/real(ncentral)
            yacentral(ia)=yacentral(ia)*dmax/real(ncentral)
            zacentral(ia)=zacentral(ia)*dmax/real(ncentral)
         end do
         fname_xyz=fname_inp(1:indexlastdot(fname_inp)-1)//'-second.xyz'
         open(nout,file=fname_xyz,status='unknown',iostat=ierr)
         if(ierr.eq.0) then
! Write number of atoms in molecule
            write(nout,'(i6)') natocentral
! Write information line - not used for plotting
            write(nout,'(a)') 'Central molecule coordinate file: '//fname_xyz(1:len_trim(fname_xyz))
            do ia=1,natocentral
               write(textline,'(a3,1x,f10.5,1x,f10.5,1x,f10.5,1x,a3)') atsymbolcentral(ia),xacentral(ia),yacentral(ia),zacentral(ia),attypecentral(ia)
               ! Append a comment if this atom is part of the central molecule definition
               if(centraldefinition(ia)) then
                  textline=textline(1:len_trim(textline))//'  #CENTRALMOLECULE'
               end if
               write(nout,'(a)') textline(1:len_trim(textline))
            end do
            close(nout)
         end if
      end if
c
c initialize counters for number of molecules of each kind in the system
c
      if(ldistinct.eq.0) then
            if(nm(nmol).eq.1.or.nm(nmol).eq.2) then
                  lsum1=1
                  lsum2=1
            else
                  lsum1=0
                  lsum2=0
            endif
      else 
            lsum1=0
            lsum2=0
            if(nm(nmol).eq.1) lsum1=1
            if(nm(nmol).eq.2) lsum2=1
      endif
      write(6,*) nmol,nm(nmol),lsum1,lsum2
      do imol1=1,nmol-1
         if(nm(imol1).ne.0) then
            write(6,*) imol1
c     call system_clock(time3,rate,max)
c
c increment counters for number of molecules of types 1 and 2 in the system
c
            if(ldistinct.eq.0) then
                  if(nm(imol1).eq.1.or.nm(imol1).eq.2) then
                        lsum1=lsum1+1
                        lsum2=lsum2+1
                  endif
            else 
                  if(nm(imol1).eq.1) lsum1=lsum1+1
                  if(nm(imol1).eq.2) lsum2=lsum2+1
            endif
            xc1=xcent(imol1)
            yc1=ycent(imol1)
            zc1=zcent(imol1)
c
c calc now look at second molecules
c
            nm2f=imol1+1
            do imol2=nm2f,nmol
               if(imol1.ne.imol2.and.nm(imol2).ne.0) then
                  iflag1=0
                  iflag2=0
                  nm1=nm(imol1)
                  nm2=nm(imol2)
                  if(ldistinct.eq.0.and.nm1.eq.1.and.nm2.eq.1) then
                     iflag1=1
                     iflag2=1
                  endif
                  if(ldistinct.eq.1.and.nm1.eq.1.and.nm2.eq.2) then
                     iflag1=1
                     iflag2=0
                  endif
                  if(ldistinct.eq.1.and.nm1.eq.2.and.nm2.eq.1) then
                     iflag1=0
                     iflag2=1
                  endif
                  xc2=xcent(imol2)
                  yc2=ycent(imol2)
                  zc2=zcent(imol2)
                  xdif=xc2-xc1
                  ydif=yc2-yc1
                  zdif=zc2-zc1
! Apply periodic boundary conditions
                  call minimagec(xdif,ydif,zdif)
                  if(iflag1.eq.1) then
! Need to rotate coordinates of this molecule onto the frame of the first molecule
                     call rotateeuler(xdif,ydif,zdif,-chim(imol1),-acos(cthetm(imol1)),-phim(imol1),xdifr,ydifr,zdifr)
                     call minimagec(xdifr,ydifr,zdifr)
!                     xdifr=xdif
!                     ydifr=ydif
!                     zdifr=zdif
                     ix=nx+nint(xdifr/xstepd)+1
                     if(ix.gt.0.and.ix.le.nx0) then
                        iy=ny+nint(ydifr/ystepd)+1
                        if(iy.gt.0.and.iy.le.ny0) then
                           iz=nz+nint(zdifr/zstepd)+1
                           if(iz.gt.0.and.iz.le.nz0) then 
!                     write(6,*) xstepd,ystepd,zstepd
!                     write(6,*) imol1,imol2,xdifr,ydifr,zdifr,ix,iy,iz
                              rhocube(ix,iy,iz)=rhocube(ix,iy,iz)+1.0
                           end if
                        end if
                     end if
                  end if
! We add on the case when molecule2 is type 1 and at the origin
                  if(iflag2.eq.1) then
                     call rotateeuler(-xdif,-ydif,-zdif,-chim(imol2),-acos(cthetm(imol2)),-phim(imol2),xdifr,ydifr,zdifr)
                     call minimagec(xdifr,ydifr,zdifr)
!                     xdifr=-xdif
!                     ydifr=-ydif
!                     zdifr=-zdif
                     ix=nx+nint(xdifr/xstepd)+1
                     if(ix.gt.0.and.ix.le.nx0) then
                        iy=ny+nint(ydifr/ystepd)+1
                        if(iy.gt.0.and.iy.le.ny0) then
                           iz=nz+nint(zdifr/zstepd)+1
                           if(iz.gt.0.and.iz.le.nz0) then 
!                  write(6,*) xstepd,ystepd,zstepd
!                  write(6,*) imol2,imol1,xdifr,ydifr,zdifr,ix,iy,iz
                              rhocube(ix,iy,iz)=rhocube(ix,iy,iz)+1.0
                           end if
                        end if
                     end if
!                     write(6,*) imol2,imol1,xdifr,ydifr,zdifr,ix,iy,iz
                  endif
               endif
            end do
c     call system_clock(time4,rate,max)
c     write(6,*) rate,max,time4-time3,nsec
         endif
      end do
      write(6,502) lsum1,lsum2
502   format(2(1x,i5))
c renormalize g(r)'s to density of second atoms in sphere
      vol=(alatx*(blaty*clatz))
      rho=real(lsum1)*real(lsum2)/vol
      write(6,193) rho
193   format(1x,'rho = ',e13.6)
c volume of element
      gfac=xstepd*ystepd*zstepd
      gfac1=1.0/(gfac*rho)
      do ix=1,nx0
         do iy=1,ny0
            do iz=1,nz0
               rhocube(ix,iy,iz)=rhocube(ix,iy,iz)*gfac1
            end do
         end do
      end do
      return
      end










