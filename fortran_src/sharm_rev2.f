      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_sharm()

      stop
      end      

      subroutine run_sharm()
c
c performs calculation of spherical harmonic coefficients from .ato file
c
c uses arithmetic centre of molecule rather than a particular atom to define
c molecular centre. 
c
      use sharm_mod
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'sharm.inc'
      include 'sharm_sdf.inc'
      character*50 formatout
      character*2 ifmt1,ifmt2,ifmt3,ifmt4,il1fmt,il2fmt,ilfmt,in1fmt,in2fmt
c
c setup the SHARM symbols
c
      call setup_symbols_sharm()
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

      pi=4.0*atan(1.0)
      piconv=pi/180.0
c
c read the input file
c
      call read_input_file(nin,nout,ierrout)
c
c if an error occurred on read, exit now
c
      if(ierrout.ne.0) return
      iref=1
c
c set up number of radius values
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
c set up number of coefficients
c
      iref=iref+1
      read(symbol_value(iref),*) ndist
c
c set up l1 values
c
      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,1,256,nworddim)
      nll1v=nwords
      if(allocated(ll1v)) deallocate(ll1v)
      allocate(ll1v(nll1v))
      do i=1,nll1v
         read(symbol_value(iref)(ncf(i):ncl(i)),*,iostat=ierr) ll1v(i)
      end do
c
c set up l2 values
c
      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,1,256,nworddim)
      nll2v=nwords
      if(allocated(ll2v)) deallocate(ll2v)
      allocate(ll2v(nll2v))
      do i=1,nll2v
         read(symbol_value(iref)(ncf(i):ncl(i)),*,iostat=ierr) ll2v(i)
      end do
c
c set up l values
c
      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
      nllv=nwords
      if(allocated(llv)) deallocate(llv)
      allocate(llv(nllv))
      do i=1,nllv
         read(symbol_value(iref)(ncf(i):ncl(i)),*) llv(i)
      end do
c
c set up n1step
c
      iref=iref+1
      read(symbol_value(iref),*) n1step
c
c set up n2step
c
      iref=iref+1
      read(symbol_value(iref),*) n2step
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
         ic=i+1
         read(symbol_value(iref)(ncf(ic):ncl(ic)),*) lztyp1(i)
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
         ic=i+1
         read(symbol_value(iref)(ncf(ic):ncl(ic)),*) lytyp1(i)
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
            ic=i+1
            read(symbol_value(iref)(ncf(ic):ncl(ic)),*) lztyp2(i)
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
            ic=i+1
            read(symbol_value(iref)(ncf(ic):ncl(ic)),*) lytyp2(i)
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
         iref=iref+2
      endif
c     write(6,*) ncent1,ncent2,ncz1,ncy1,ncz2,ncy2
      !Decide whether to use simulated phim, thetam, chim
      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
      performtest=symbol_value(iref)(ncf(1):ncl(1)).ne.'notest'.and.nwords.gt.7
      if(performtest) then
         read(symbol_value(iref)(ncf(1):ncl(nwords)),*,iostat=ierr) phim1,thetam1,chim1,phim2,thetam2,chim2,thetallab,phillab
         phim1=phim1*piconv
         thetam1=thetam1*piconv
         chim1=chim1*piconv
         phim2=phim2*piconv
         thetam2=thetam2*piconv
         chim2=chim2*piconv
         thetallab=thetallab*piconv
         phillab=phillab*piconv
      else
         phim1=0.0
         thetam1=0.0
         chim1=0.0
         phim2=0.0
         thetam2=0.0
         chim2=0.0
         thetallab=0.0
         phillab=0.0
      end if
      if(ierr.ne.0) performtest=.false.
      iref=iref+1
      outputaveragemol=symbol_value(iref).eq.'y'

c Only perform the calculation when nsum is zero or a multiple of nmult

      nmult=5
      nsum=nsumt/nmult
      nt=nsum*nmult
      if(nt.eq.nsumt) then
c
c normalise atomic coordinates to box size
c
         call normalise_coordinates()
         zmaxd=zmax/dmax
         if(zmaxd.gt.0.5) then
            zmaxd=0.5
            zmax=zmaxd*dmax
         endif
! Set up l and n values and allocate some needed arrays
         call get_l_and_n_values()
         mz=nz
         mdist=ndist
         if(allocated(z)) deallocate(z,gofz,gofzs,gzt,gzt2,titledat)
         allocate(z(mz),gofz(mz,mdist),gofzs(mz,mdist),gzt(mz,mdist),gzt2(mz,mdist),titledat(mdist))
! set up step size and z values
         zstep=zmaxd/float(nz)
         zhalf=0.5*zstep
         do j=1,nz
            z(j)=zstep*j-zhalf
         end do
c zero gofz
         call initgofzs(nz,ndist,gofz,gofzs,mz,mdist)
         call sharm(nin,fname_inp)
c
c set the z values
c
         do i=1,nz
            z(i)=z(i)*dmax
         end do
c
c axis title
c
         titler='   r'
c
c possible formats for integers in title
c
         ifmt1='i1'
         ifmt2='i2'
         ifmt3='i3'
         ifmt4='i4'
c
c set the column titles
c
         do id=1,ndist
            il1=ll1(id)
            il2=ll2(id)
            il=ll(id)
            in1=n1(id)
            in2=n2(id)
c
c set up the format for the title
c
            if(il1.lt.10) then
               il1fmt=ifmt1
            else if(il1.lt.100) then
               il1fmt=ifmt2
            else if(il1.lt.1000) then
               il1fmt=ifmt3
            endif
            if(il2.lt.10) then
               il2fmt=ifmt1
            else if(il2.lt.100) then
               il2fmt=ifmt2
            else if(il2.lt.1000) then
               il2fmt=ifmt3
            endif
            if(il.lt.10) then
               ilfmt=ifmt1
            else if(il.lt.100) then
               ilfmt=ifmt2
            else if(il.lt.1000) then
               ilfmt=ifmt3
            endif
            if(in1.lt.10) then
               in1fmt=ifmt1
            else if(in1.lt.100) then
               in1fmt=ifmt2
            else if(in1.lt.1000) then
               in1fmt=ifmt3
            endif
            if(abs(in2).lt.10) then
               if(in2.lt.0) then
                  in2fmt=ifmt2
               else
                  in2fmt=ifmt1
               end if
            else if(abs(in2).lt.100) then
               if(in2.lt.0) then
                  in2fmt=ifmt3
               else
                  in2fmt=ifmt2
               end if
            else if(abs(in2).lt.1000) then
               if(in2.lt.0) then
                  in2fmt=ifmt4
               else
                  in2fmt=ifmt3
               end if
            endif
            formatout='(1h('//il1fmt//',1h,,'//il2fmt//',1h,,'//ilfmt//',1h;,'//in1fmt//',1h,,'//in2fmt//',1h))'
c
c write the title for this distribution
c
            write(titledat(id),fmt=formatout) il1,il2,il,in1,in2
         end do
         ext='h'
         call outgofzs_sharm(fnameinp,ext,nz,ndist,nsum,z,gofz,gofzs,gzt,gzt2,titler,titledat,mz,mdist)
      endif
      if(.not.performtest) nsumt=nsumt+1
c
c write input file
c
c
c make sure all current values are stored
c
      write(symbol_value(4),*) nsumt
      write(symbol_value(5),*) ndist
      write(6,*) ndist,' ',symbol_value(5)
      call save_files_others(14)
      close(14)
      if(performtest) write(6,*) 'sharm> Test performed.'
      return
      end

      subroutine sharm(nout,fname_inp)

      use reallocation_routines
      use sharm_mod

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      include 'sharm.inc'
      include 'sharm_sdf.inc'

c     integer*4 rate,max,time1,time2,time3,time4
      real amass,amasssum

      logical                                      :: centralmolecule,found,foundatom
      character(len=*)                             :: fname_inp
      character(len=256)                           :: fname_xyz,textline !Filename for xyz data, line of text
      integer                                      :: nout,ierr !Unit number on which to write xyz file
      real                                         :: xt,yt,zt,xr,yr,zr
      real, dimension(:), allocatable              :: xacentral,yacentral,zacentral
      integer                                      :: ncentral,natocentral
      character(len=3), dimension(:), allocatable  :: atsymbolcentral,attypecentral
      logical, dimension(:), allocatable           :: centraldefinition
      real, dimension(:), allocatable              :: xat,yat,zat

      pi=4.0*atan(1.0)
      piconv=pi/180.0
      degconv=180.0/pi
      zstep=zmaxd/float(nz)
      gfac=4.0*pi/3.0
      vol=(alatx*(blaty*clatz))
c set up counter for second atoms
      nsec=0
! Arrays to store Euler angles
      if(allocated(nm)) deallocate(nm,phim,cthetm,chim,xcent,ycent,zcent,dlmnmc,chiphi,dlmnc,philab)
      call get_max_size_dlmn()
      allocate(nm(nmol),phim(nmol),cthetm(nmol),chim(nmol),xcent(nmol),ycent(nmol),zcent(nmol))
      allocate(dlmnmc(2,nmol,mdlmn),chiphi(2,nmol,mdlmn),dlmnc(2,mdlmn),philab(mdlmn))
      if(performtest) then
! Allocate suitable arrays if performing a test
         if (allocated(phimt1)) deallocate(phimt1,cthetmt1,chimt1,phimt2,cthetmt2,chimt2,dlmnmct,chiphit)
         allocate(phimt1(nmol),cthetmt1(nmol),chimt1(nmol),phimt2(nmol),cthetmt2(nmol),chimt2(nmol))
         allocate(dlmnmct(2,nmol,mdlmn),chiphit(2,nmol,mdlmn))
      endif
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
!If testing orientations, supply a fictitious but known distribution of angles
            if(performtest) then
!               write(6,499) imol1,phim1,thetam1,chim1,phim2,thetam2,chim2
!               call populate_disk(pi,phim1,thetam1,chim1,phinew,thetanew,chinew,thetaspread,phispread)
! Rotate this molecule onto desired orientation
               phimt1(imol1)=phim1
               cthetmt1(imol1)=cos(thetam1)
               chimt1(imol1)=chim1
!               call populate_disk(pi,phim2,thetam2,chim2,phinew,thetanew,chinew,thetaspread,phispread)
               phimt2(imol1)=phim2
               cthetmt2(imol1)=cos(thetam2)
               chimt2(imol1)=chim2
!               write(6,499) imol1,phim(imol1),cthetm(imol1),chim(imol1),phimt(imol1),cthetmt(imol1),chimt(imol1)
499   format(1x,i5,6(1x,e12.5))
            end if
            call form_dlmn(imol1,ierr)
            if(ierr.ne.0) then
               write(6,*) 'sharm> Error when calculating dlm coefficients for molecule ',imol1
               return
            end if
c
c step through all the atoms and ignore those which do not match the required
c ltype's
         else
            nm(imol1)=0
         endif
      end do
      write(6,333) nsmmo1,nsmmo2
333   format(1x,'sharm> Number of type 1 molecules = ',i10
     1/1x, 'sharm> Number of type 2 molecules = ',i10)
!      stop
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
                  if(nato(im).gt.1) then
                     call rotateeuler(xt,yt,zt,-chim(im),-acos(cthetm(im)),-phim(im),xr,yr,zr)
                  else
                     xr=xt
                     yr=yt
                     zr=zt
                  end if
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
                  else if(outputaveragemol) then
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
                  if(nato(im).gt.1) then
                     call rotateeuler(xt,yt,zt,-chim(im),-acos(cthetm(im)),-phim(im),xr,yr,zr)
                  else
                     xr=xt
                     yr=yt
                     zr=zt
                  end if
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
            endif
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
      nmolsharm=nmol
      call form_sum_prod_dlmn_cgcoeff(nmolsharm,ierr)
      return
      end

