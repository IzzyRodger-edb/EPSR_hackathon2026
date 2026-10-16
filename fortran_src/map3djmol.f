      subroutine map3djmol(infile,plotfile,nextfile,plotnumber)
      
      use molecule_coordinates_mod

      real, parameter :: pi=3.141592653,bohr=0.52917721092
      integer, parameter :: mexpan=3000,mval=10,nworddim=201,nstsize=256
      real datac(500,mexpan),datas(500,mexpan),rd(500)
     1,zcoeffc(mexpan),zcoeffs(mexpan),grw(50,50,500),temp(500),rot,elev,rgbobj(3),rgbbak(3)

      real, dimension(:,:,:), allocatable :: dens
      integer l1ex(mexpan),l2ex(mexpan),lex(mexpan),n1ex(mexpan),n2ex(mexpan)
! Used to specify spheres to be used to plot atoms
      integer :: nspher
      real, dimension(:), allocatable :: rdcoor,thcoor,phcoor,radbal
      real, dimension(:,:), allocatable :: rgbbal
      character(len=3), dimension(:), allocatable :: chsymbol
      character(len=256) :: fname,outfile,moleculename,moleculefile
      integer :: lenmoleculename,get_atomicnumber,leninfile,lenmoleculefile
      character(len=256), save :: firstmoleculename ! This will be used to check that subsequent central molecules are the same as the first
      integer, save :: lenfirstmoleculename,lennextfile,lenmolfile
      character(len=3) :: spin
      character(len=*) :: infile,plotfile,nextfile
      character*80 surfacetype,mesh,transparency
      character(len=256), save :: molfile
      character(len=256), dimension(:), allocatable :: molfileline
      character(len=256) :: molfileline1,molfileline2,molfileline3

      character*10 surfacelabel
      integer nref,nvoxel,mpixel,imesh,itrans,plotnumber,kk0
      integer i,j,k,ir,irm,irp,nth,nthm,nthp,nph,nphm,nphp

      logical newstart,rotatecentralmolecule
c used to specify words in sometext
      character(len=nstsize) sometext
      integer nwords,irmax,ncf(nworddim),ncl(nworddim)
      real scalefactor,translucency,suminterp,sumc,sums,sumt

      external cosd,sind,acosd,asind,atan2d

! Decide whether this is a new plot, or adding to an existing plot
      newstart=nextfile(1:len_trim(nextfile)).eq.'<newstart>'
! Generate an isosurface label
      write(surfacelabel,'(a1,i5.5)') 's',plotnumber
c Initialise random number generator
      open(10,file='ran1.dat',status='unknown')
         call initran1(10)
      close(10)
      nxdim=500
      nfdim=mexpan
c
c initialise log(factorial)
c
      call fctrln(250)
!      write(6,*) 'Type name of file to read plotting parameters from'
!      read(5,98) fname
      fname=infile
! Name of output file - remove any text after the last dot
      lenfirstdot=index(fname,'.CUBE.txt') ! Make sure this is a .CUBE file
      if(lenfirstdot.gt.0) then
         outfile=fname(1:lenfirstdot)//'cube'
      else
         lendfirstdot=indexlastdot(fname)
         if(lenfirstdot.gt.0) then
            outfile=fname(1:lenfirstdot)//'cube'
         else
            outfile=fname(1:len_trim(fname))//'.cube'
         end if
      end if
      lenoutfile=len_trim(outfile)
      if(newstart) then
         plotfile='plot'//outfile(1:indexlastdot(outfile))//'txt'
         molfile=outfile(1:index(outfile,'.cube')-1)//'-mol.cube'
         lenmolfile=len_trim(molfile)
      end if
      nin=10
      open(nin,file=fname,status='old')
      write(6,99)
99    format(1x,'Type name of file containing coefficients')
      read(nin,98) sometext
98    format(a)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      call concatenate(1,nwords,infile,leninfile,nwords,ncf,ncl,sometext)
! Need to check that the same set of SHARM coefficients are used for this as for previous contours
      moleculename=infile
      lenmoleculename=len_trim(moleculename)
! Save this name if this the first plot of the series
      if(newstart) then
         firstmoleculename=moleculename
         lenfirstmoleculename=lenmoleculename
      end if
      if(moleculename(1:lenmoleculename).ne.firstmoleculename(1:lenfirstmoleculename)) then
         write(6,'(a,a,a,a)') 'map3djmol> Molecule symbols: ',moleculename(1:lenmoleculename)
     *,':',firstmoleculename(1:lenfirstmoleculename)
         write(6,'(a/a/a)') ' map3djmol> You need to ensure the central molecule for file '
     *,fname(1:len_trim(fname))
     *,' is the same as for the first plot.'
         return
      endif
      write(6,97)
97    format(1x,'Type no. of expansion coefficients')
      read(nin,*) Nexpan
c
c clean out data arrays
c
      do k=1,nexpan
            do i=1,500
                  datac(i,k)=0.0
                  datas(i,k)=0.0
            end do
      end do
      write(6,96)
96    format(1x,'Are molecules identical? (ident=0 if so, else +1 or -1)')

c A value of ident .lt. 0 will be used to signify that we want molecule 2 to be
c at the origin

      read(nin,*) ident

c save the sign of ident and keep separate from the absolute value

      if(abs(ident).gt.0) then
            identsign=isign(1,ident)
            ident=abs(ident)
      else
            identsign=1
            ident=0
      endif

      write(6,91)
91    format(1x,'Use SHARM/SDF construct (1) or all central atoms (0)')
      read(nin,*) kk0
982   format(1x,4(1x,i5),f12.7)
      if(kk0.lt.0) kk0=0
      if(kk0.gt.2) kk0=2
c
C read in spherical harmonics data
c
      call r_int(nr,infile,nexpan,rd,datac,datas,nxdim,nfdim,ierr)
      do k=1,nexpan
            l1ex(k)=nint(datac(1,k))
            l2ex(k)=nint(datac(2,k))
            lex(k)=nint(datac(3,k))
            n1ex(k)=nint(datac(4,k))
            n2ex(k)=nint(datac(5,k))
      end do
!      do k=1,nexpan
!      write(6,92) l1ex(k),l2ex(k),lex(k),n1ex(k),n2ex(k)
!      end do
!92    format(1x,10(1x,i5))
      write(6,*) 'Type number of smoothings on the data'
      read(nin,*) nsmoo
      write(6,101)
101   format(1x,'Type rmax')
      read(nin,*) xmax
      rstep=rd(3)-rd(2)
      write(6,569)
569   format(1x,'Type number of pixels on each axis')
      read(nin,*) mpixel
      write(6,104)
104   format(1x,'Type minimum and maximum radius to include')
      read(nin,*) rmin,rmax
      write(6,491) 
491   format(1x,'Type fractional isosurface level (-ve for absolute)')
      read(nin,*) surfra
      write(6,95)
95    format(1x,'Use only l1=0 (0) and/or l2=0 (0)?')
      read(nin,*) l1only,l2only
      write(6,93)
93    format(1x,'Use only n1=0 (0) and/or n2=0 (0)?')
      read(nin,*) n1only,n2only
      write(6,94)
94    format(1x,'Use only m2=0 (0)?')
      read(nin,*) m2fix
      write(6,105)
105   format(1x,'Vary thetal and phil(1), thetam and phim (2),',
     1' or thetam and chim (3) for second molecule')
      read(nin,*) nvary
      if(nvary.eq.1) then
            write(6,106)
106   format(1x,'Specify phim, thetam and chim of molecule')
            read(nin,*) phim,thetam,chim
      else if(nvary.eq.2) then
            write(6,107)
107   format(1x,'Specify thetal, phil, and chim of molecule')
            read(nin,*) thetal,phil,chim
      else if(nvary.eq.3) then
            write(6,108)
108   format(1x,'Specify thetal, phil, and phim of molecule')
            read(nin,*) thetal,phil,phim
      endif
121   format(a50)
      write(6,851)
851   format(1x,' Type no. of spheres to plot')
      read(nin,*) nspher
      nsize=nspher
      if(allocated(rdcoor)) deallocate(rdcoor,thcoor,phcoor,radbal,rgbbal,chsymbol)
      allocate(rdcoor(nsize),thcoor(nsize),phcoor(nsize),radbal(nsize),rgbbal(3,nsize),chsymbol(nsize))
      ic=0
      write(6,852)
852   format(1x,' For each sphere give its chemical symbol, and position in (r,theta,phi)')
      piconv=3.141592653/180.0
      do is=1,nspher
         read(nin,'(a)') sometext
         if(kk0.gt.1.and.ic.lt.size(chsymbol)) then
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            ic=ic+1
            read(sometext(ncf(1):ncl(1)),'(a2)') chsymbol(ic)
            read(sometext(ncf(2):ncl(4)),*) rdcoor(ic),thc,phc
            thcoor(ic)=thc*piconv
            phcoor(ic)=phc*piconv
         end if
      end do
      if(kk0.lt.2) then !Use atoms read in from xyz file
         ic=0
         if(nvary.gt.1) then
            moleculefile=infile(1:indexlastdot(infile)-1)//'-second.xyz'
            lenmoleculefile=len_trim(moleculefile)
         else
            moleculefile=infile(1:indexlastdot(infile)-1)//'-first.xyz'
            lenmoleculefile=len_trim(moleculefile)
         end if
         call get_coords_from_xyz(moleculefile,lenmoleculefile)
         nsize=nofatoms
         if(allocated(rdcoor)) deallocate(rdcoor,thcoor,phcoor,radbal,rgbbal,chsymbol)
         allocate(rdcoor(nsize),thcoor(nsize),phcoor(nsize),radbal(nsize),rgbbal(3,nsize),chsymbol(nsize))
         is=0
         do while (is.lt.nofatoms.and.ic.lt.size(chsymbol))
            is=is+1
            if(centralmolecule(is).ge.kk0) then
               ic=ic+1
               chsymbol(ic)=chsymbolt(is)
               rdcoor(ic)=rt(is)
               thcoor(ic)=thetat(is)
               phcoor(ic)=phit(is)
            end if
         end do
      end if
      nspher=ic
!      write(6,853)
!853   format(1x,' Type character size (e.g.1.5), line width (e.g.2)'
!     *,' and colour index')
!      read(nin,*) size,iwidth,icoll
!      write(6,793)
!793   format(1x,' Type title for this plot')
!794   format(a80)
!      write(6,795)
!795   format(1x,' Type x and y coordinates of title, and char. size')
!      read(nin,*) xtr,ytr,csize
      write(6,487)
487   format(1x,' Type 3 RGB Colour indices for background')
      read(nin,*) rgbbak(1),rgbbak(2),rgbbak(3)
!      write(6,485)
!485   format(1x,' Type ishade (1-8): 0 means no shading')
!      read(nin,*) ishade
      write(6,489)
489   format(1x,' Type 3 RGB Colour indices for object')
      read(nin,*) rgbobj(1),rgbobj(2),rgbobj(3)
      write(6,*) rgbobj(1),rgbobj(2),rgbobj(3)
!      write(6,488)
!488   format(1x,' Type 3 (x,y,z) coordinates for light source')
!      read(nin,*) xli,yli,zli
!      write(6,*) xli,yli,zli
!      write(6,486)
!486   format(1x,' Type fade factor (0 = no fading, 1=full fading)')
!      read(nin,*) fadefc
      write(6,484)
484   format(1x,' Transparency of surface (0=opaque,1=translucent,2=transparent)')
      read(nin,*) itrans
      if(itrans.eq.0) then
         transparency='fill'
      else if(itrans.eq.1) then
         transparency='translucent'
      else
         transparency='nofill'
      end if
      lentransparency=len_trim(transparency)
      write(6,*) itrans
      write(6,490)
490   format(1x,' Mesh? (0=nothing,1=mesh,2=dots)')
      read(nin,*) imesh
      if(itrans.gt.1.and.imesh.eq.0) imesh=1 !Ensure something is plotted!
      if(imesh.eq.0) then
         mesh='nomesh nodots'
      else if(imesh.eq.1) then
         mesh='mesh'
      else
         mesh='dots'
      end if
      lenmesh=len_trim(mesh)
      surfacetype=mesh(1:lenmesh)//' '//transparency(1:lentransparency)
      write(6,492)
492   format(1x,' Type rotation, tilt, rotatecentralmolecule and translucency')
      read(nin,'(a)') sometext
      scalefactor=1.0
      translucency=0.5
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.gt.1) then
         read(sometext(ncf(1):ncl(2)), *) rot,elev
         if(nwords.gt.2) then
            read(sometext(ncf(3):ncl(3)),*) scalefactor
            rotatecentralmolecule=scalefactor.gt.0.0
            if(.not.rotatecentralmolecule) scalefactor=1.0
            if(nwords.gt.3) then
               read(sometext(ncf(4):ncl(4)),*) translucency
            end if
         else
            rotatecentralmolecule=.false.
            scalefactor=1.0
            translucency=0.5
         end if
      else
         rot=0.0
         elev=0.0
         rotatecentralmolecule=.false.
         scalefactor=1.0
         translucency=0.5
      end if
      read(nin,'(a)') sometext
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      call concatenate(1,1,spin,ilen,nwords,ncf,ncl,sometext)
      if(spin.ne.'x'.and.spin.ne.'y'.and.spin.ne.'z') spin='off' !Ensure a sensible value of spin has been set. Actual behaviour will be determined by last plot of a sequence.
      read(nin,'(a)') sometext
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      call concatenate(1,nwords,nextfile,lennextfile,nwords,ncf,ncl,sometext)
      close(nin)
! Set up pixels and voxels. Currently a cubic plot area is assumed.
      nx=mpixel
      xstep=xmax/real(nx)
      ny=nx
      nz=nx
      ystep=xstep
      zstep=xstep
      nx0=2*nx+1
      ny0=2*ny+1
      nz0=2*nz+1
      if(allocated(dens)) deallocate(dens)
      allocate(dens(nx0,ny0,nz0))
      nvoxel=nx0*ny0*nz0
c calculate g(r,w) at 10 deg intervals
      nt0=19
      tstep=10
      np0=36
      pstep=10
c
c initial maximum and minimum values
c
      dhigh=1.0
      dlow=1.0
      phimmax=0.0
      thetammax=0.0
      chimmax=0.0
      irmax=0 !Reference to distance at which maximum intensity is found
      do it=1,nt0
         t=real(it-1)*tstep
         do ip=1,np0
            p=(ip-1)*pstep
            if(nvary.eq.1) then
                  thetal=t
                  phil=p
                  phi=phim
                  chi=chim
            else if(nvary.eq.2) then
                  phi=p
                  thetam=t
                  chi=chim
            else if(nvary.eq.3) then
                  phi=phim
                  thetam=t
                  chi=p
            endif
            ctl=cosd(thetal)
            cz=cosd(thetam)
            do k=1,nexpan
               sumc=0.0
               sums=0.0
               ll1=l1ex(k)
               ll2=l2ex(k)
               ll=lex(k)
               n1=n1ex(k)
               n2=n2ex(k)
               l1tmp=ll1
               l2tmp=ll2
               n1tmp=n1
               n2tmp=n2
               if(l1only.eq.0.and.l1tmp.ne.0) go to 595
               if(l2only.eq.0.and.l2tmp.ne.0) go to 595
               if(n1only.eq.0.and.n1tmp.ne.0) go to 595
               if(n2only.eq.0.and.n2tmp.ne.0) go to 595
               phs=1.0
               phsn=((-1)**(l1tmp+l2tmp+ll+n1tmp+n2tmp))

c Normal algorithm for molecule 1 at the origin

               if(identsign.ge.0) then
                  if(m2fix.gt.0) then
                     m2f=-l2temp
                     m2l=l2tmp
                  else
                     m2f=0
                     m2l=0
                  endif
                  call sumoverm2(l1tmp,l2tmp,ll,n1tmp,m2f,m2l,n2tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
c now reverse signs of n1,n2 if either is non-zero
                  if(n1.gt.0) then
                     n1tmp=-n1
                     n2tmp=-n2
                     phs=phs*phsn
                     call sumoverm2(l1tmp,l2tmp,ll,n1tmp,m2f,m2l,n2tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
                  endif
               else
                  if(m2fix.gt.0) then
                     m2f=-l1tmp
                     m2l=l1tmp
                  else
                     m2f=0
                     m2l=0
                  endif

c Phase factor to reverse parity of centre to centre vector

                  phs=(-1)**ll
                  call sumoverm2(l2tmp,l1tmp,ll,n2tmp,m2f,m2l,n1tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
c now reverse signs of n1,n2 if either is non-zero
                  if(n1.gt.0) then
                     n1tmp=-n1
                     n2tmp=-n2
                     phs=phsn*phs
                     call sumoverm2(l2tmp,l1tmp,ll,n2tmp,m2f,m2l,n1tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
                  endif
               endif
c now repeat for interchange of l1,n1,with l2,n2 for identical molecules
c bearing in mind the phase convention on the h's - the case when l1=l2 AND
c n1=n2 doesn't count here. 
c
595            continue
               if(ident.eq.0.and.(ll1.ne.ll2)) then !l1 will always be .ge. l2 (from the way they are calculated).
                  phs=(-1)**(ll1+ll2)
                  l1tmp=ll2
                  l2tmp=ll1
                  n1tmp=n2
                  n2tmp=n1
                  if(l1only.eq.0.and.l1tmp.ne.0) go to 596
                  if(l2only.eq.0.and.l2tmp.ne.0) go to 596
                  if(n1only.eq.0.and.n1tmp.ne.0) go to 596
                  if(n2only.eq.0.and.n2tmp.ne.0) go to 596
                  phsn=((-1)**(l1tmp+l2tmp+ll+n1tmp+n2tmp))
                  if(m2fix.gt.0) then
                     m2f=-l2tmp
                     m2l=l2tmp
                  else
                     m2f=0
                     m2l=0
                  endif
                  call sumoverm2(l1tmp,l2tmp,ll,n1tmp,m2f,m2l,n2tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
c now reverse signs of n1,n2 if either is non-zero
                  if(n1.gt.0) then
                     n1tmp=-n2
                     n2tmp=-n1
                     phs=phs*phsn
                     call sumoverm2(l1tmp,l2tmp,ll,n1tmp,m2f,m2l,n2tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
                  endif
               endif
596            continue
               zcoeffc(k)=sumc
               zcoeffs(k)=sums
!               write(6,*) ll1,ll2,ll,n1,n2,phil,ctl,phi,cz,sumc,sums
            end do
            do i=1,nr
                  sumt=0
                  do k=1,nexpan
                        sumt=sumt+datac(i,k)*zcoeffc(k)+datas(i,k)*zcoeffs(k)
                  end do
c
c test for maximum and minimum in selected r range
c
                  r=rd(i)
                  temp(i)=0.0
                  if(r.ge.rmin) then
                     if(r.le.rmax) then
                        temp(i)=sumt+1.0
                     else
                        temp(i)=1.0
                     end if
                  endif
            end do
c
c perform smoothing on coefficients if needed - this is primarily to remove
c hard edges and so avoid fringes in the maps 
c
            if(nsmoo.gt.0) then
                  call smoosharm(1,nr,nsmoo,temp)
            endif
            do i=1,nr
                  grw(it,ip,i)=temp(i)
                  if(grw(it,ip,i).lt.dlow) dlow=grw(it,ip,i)
                  if(grw(it,ip,i).gt.dhigh) then
                     dhigh=grw(it,ip,i)
                     if(rd(i).ge.rmin) then
                        irmax=i
                        phimmax=phi
                        thetammax=thetam
                        chimmax=chi
                     endif
                  endif
            end do
          end do
      end do
      write(6,*) 'Max. and min. values in this function: ',dhigh,dlow
      if(surfra.ge.0.0) then
            dsurf=dhigh+0.9*(dlow-dhigh)
      else
            dsurf=-surfra
      endif
c
c form the sums of r**2 * g(r,w) between rmin and rmax, above DSURF, and in 
c total
c
      fracold=1.0
      frac=0.9
      ntimes=0
200   continue
      ntimes=ntimes+1
      sumt=0.0
      sump=0.0
      delr=4.0*pi*(rd(3)-rd(2))
      do i=1,nr
            r=rd(i)
            if(r.ge.rmin.and.r.le.rmax) then
                  r2=r*r*delr
                  do ip=1,np0
                  do it=1,nt0
                        test=grw(it,ip,i)
                        sumt=sumt+r2*test
                        if(test.ge.dsurf) sump=sump+r2*test
                  end do
                  end do
            endif
      end do
      sumt=sumt/nt0/np0
      sump=sump/nt0/np0
      write(6,*) 'Total integral between ',rmin,' and ',rmax,sumt
      write(6,*) 'Integral for this map above ',dsurf, 'between '
     *,rmin,' and ',rmax,sump
      if(sumt.gt.0.0)   then
            srat=sump/sumt
            write(6,*) 'Ratio = ',100*srat,'%'
c
c if surfra is +ve it represents the fraction of molecules 
c at this distance
c range which are to be plotted in the contour plot. Therefore adjust the 
c value of dsurf until the fraction of molecules plotted is at the specified 
c level
c
            if(surfra.gt.0.0) then
                  fact=srat/surfra
                  if(fact.gt.1.01.or.fact.lt.0.99) then
c
c estimate next value
c
                        if(srat.lt.surfra) then
c
c increase frac
c
                              frac=0.5*(frac+fracold)
                        else
                              fracold=frac
                              frac=0.9*fracold
                        endif
                        dsurf=dhigh+frac*(dlow-dhigh)
c
c dsurf cannot be above dhigh or below dlow
c
                        if(dsurf.gt.dhigh) dsurf=dhigh
                        if(dsurf.lt.dlow) dsurf=dlow
                        if(ntimes.lt.500) go to 200
                  endif
            endif
      endif
      write(6,*) 'Phi, theta, chi at maximum density, and maximum density: ',phimmax,thetammax,chimmax,dhigh
      phirot=phimmax
      thetarot=thetammax
      chirot=chimmax

c For NVARY = 3 we need to rotate the central molecule about the z-axis by chirot degrees
c because that is the way the orientational spatial density function is plotted, namely the
c distribution out of the z-x plane represents the variation with chim.

      if(nvary.eq.3) phirot=phirot+chirot
!      write(6,*) phirot,thetarot,chirot,dhigh
      phirot=phirot*piconv
      thetarot=thetarot*piconv
      chirot=chirot*piconv
! In case it is needed define the translation for the plot
      xtrans=0.0
      ytrans=0.0
      ztrans=0.0
c
c apply the rotation. Do this in reverse order about the lab. axes
c
      if(nvary.gt.1.and.rotatecentralmolecule) then
         ztrans=rd(irmax)*cosd(thetal)
         rtrans=rd(irmax)*sind(thetal)
         xtrans=rtrans*cosd(phil)
         ytrans=rtrans*sind(phil)
         do is=1,nspher
c rotate about lab z axis
            phcoor(is)=phcoor(is)+chirot
c rotate about lab y axis, then about z axis again
            z=cos(thcoor(is))
            x=sin(thcoor(is))*cos(phcoor(is))
            y=sin(thcoor(is))*sin(phcoor(is))
            d=sqrt(x*x+z*z)
            alpha=atan2(x,z)+thetarot
            znew=d*cos(alpha)
            xnew=d*sin(alpha)
            thcoor(is)=acos(znew)
            phcoor(is)=atan2(y,xnew)+phirot
         end do
      end if
c
c now interpolate onto a cubic grid and determine whether pixel is occupied or not
c
      iref=0
      do k=1,nz0
         z=real(k-nz-1)*zstep
         z2=z*z
         do j=1,ny0
            y=real(j-ny-1)*ystep
!            yy(j)=y
            y2=y*y
            do i=1,nx0
               x=real(i-nx-1)*xstep    
!               xx(i)=x
               x2=x*x
               rxy=sqrt(x2+y2)
               r=sqrt(x2+y2+z2)
c determine pair of rd values
               sum=1.0
               if(r.le.rd(nr)) then
                  ir=nint(r/rstep)+1
                  ir=max(2,ir) !This is to prevent irm referencing undefined grw at low r
                  ir=min(nr-1,ir) !This is to prevent irp referencing undefined grw at high r
                  irp=ir+1
                  irm=ir-1
                  rrat=r/rstep-real(ir-1)
!                  if(abs(rrat).gt.0.5) then
!                       write(6,*) ir,irp,irm,rrat,r,rd(ir),rstep
!                  endif
c determine nearest theta value
                  theta=0.0
                  if(r.gt.0.0) then
                     costh=z/r
                     costh=min(1.0,costh)
                     costh=max(-1.0,costh)
                     theta=acosd(costh)
                  endif
                  nth=nint(theta/tstep)+1
                  nth=max(2,nth)
                  nth=min(nt0-1,nth)
                  nthp=nth+1
                  nthm=nth-1
                  thrat=theta/tstep-real(nth-1)
!                    if(abs(thrat).gt.1.0) then
!                         write(6,*) nth,nthp,nthm,thrat,theta,((nth-1)*tstep),tstep
!                  endif
c determine nearest phi value
                  phi=0
                  if(rxy.gt.0.) then
                     cosphi=x/rxy
                     sinphi=y/rxy
                     phi=atan2d(sinphi,cosphi)
                     if(phi.lt.0.0) phi=phi+360.0
                  end if
                  nph=nint(phi/pstep)+1
                  phrat=phi/pstep-real(nph-1)
                  do while (nph.gt.np0)
                     nph=nph-np0
                  end do
                  do while (nph.lt.1)
                     nph=nph+np0
                  end do
! Try averaging over relevant phi values. We will try to keep the solid angle of this averaging constant
! for different theta.
                  thetatrial=nth*tstep
                  sintheta=sind(thetatrial)
                  if(sintheta.gt.0.0) then
                     nphtrial=nint(0.5/sintheta)
                  else
                     nphtrial=np0/2
                  end if
                  if(nptrial.lt.1) nptrial=1
!      write(6,*) theta,phi,nth,nph,thetatrial,nphtrial
                  nphp=nph+nphtrial
                  nphm=nph-nphtrial
                  if(nphtrial.eq.1) then
                     do while(nphp.gt.np0)
                        nphp=nphp-np0
                     end do
                     do while(nphm.lt.1)
                        nphm=nphm+np0
                     end do
! Perform quadratic interpolation
                     sum=quadinterp(grw,rrat,thrat,phrat,ir,irp,irm,nth,nthp,nthm,nph,nphp,nphm,50,50,500)
                  else
! Form mean over phi values
                     sum=suminterp(grw,rrat,thrat,ir,irp,irm,nth,nthp,nthm,np0,nphp,nphm,50,50,500)
                  end if
               end if
               dens(i,j,k) = sum
!               if(i.eq.nx+1.and.j.eq.ny+1) then
!                  write(6,*) i,j,k
!                  write(6,*) irm,ir,irp,nthm,nth,nthp,nphm,nph,nphp
!                  write(6,*) rrat,thrat,phrat,r,rxy,sum,dens(i,j,k),grw(nth,nph,ir)
!               end if
               if(sum.ge.dsurf) iref=iref+1
            end do
         end do
      end do
      nref=iref
      write(6,999) nref,nvoxel
999   format(1x,' No. of occupied voxels = ',i9,' out of ',i9,' total.')
      nout=nin+1
! Translate the centre of the plot if requested
      xcentre=(xtrans-nx*xstep*scalefactor)/bohr
      ycentre=(ytrans-ny*ystep*scalefactor)/bohr
      zcentre=(ztrans-nz*zstep*scalefactor)/bohr
! The molecular coordinates will be stored in a separate file and added to if more than one surface is to be plotted
      open(nout,file=molfile,status='old',iostat=ierr)
      if(ierr.eq.0.and..not.newstart) then !This file exists and this is not a newstart, so read the existing atoms and add the new ones to it 
         read(nout,'(a)') molfileline1
         read(nout,'(a)') molfileline2
         read(nout,*) nspherold,xcentremol,ycentremol,zcentremol
         if(allocated(molfileline)) deallocate(molfileline)
         allocate(molfileline(nspherold+3))
         do i=1,nspherold+3
            read(nout,'(a)') molfileline(i)
         end do
         rewind(nout)
         write(nout,'(a)') molfileline1(1:len_trim(molfileline1))
         write(nout,'(a)') molfileline2(1:len_trim(molfileline2))
         write(nout,*) nspherold+nspher,xcenterold,ycentreold,zcentreold
         do i=1,nspherold+3
            write(nout,'(a)') molfileline(i)(1:len_trim(molfileline(i)))
         end do
! Add in the new lines
         do is=1,nspher
            xbal=rdcoor(is)*sin(thcoor(is))*cos(phcoor(is))+xtrans
            ybal=rdcoor(is)*sin(thcoor(is))*sin(phcoor(is))+ytrans
            zbal=rdcoor(is)*cos(thcoor(is))+ztrans
            iatomicno=get_atomicnumber(chsymbol(is))
            write(nout,*) iatomicno,real(iatomicno),xbal/bohr,ybal/bohr,zbal/bohr
         end do
         close(nout)
         deallocate(molfileline)
      else
         open(nout,file=molfile,status='unknown',iostat=ierr) 
         write(nout,*) '# Cube file output by map3djmol'
         write(nout,*) 'Outer loop x, middle loop y, inner loop z'
         write(nout,*) nspher,xcentre,ycentre,zcentre
         write(nout,*) 0, xstep*scalefactor/bohr, 0.0, 0.0
         write(nout,*) 0, 0.0, ystep*scalefactor/bohr, 0.0
         write(nout,*) 0, 0.0, 0.0, zstep*scalefactor/bohr
         do is=1,nspher
            xbal=rdcoor(is)*sin(thcoor(is))*cos(phcoor(is))+xtrans
            ybal=rdcoor(is)*sin(thcoor(is))*sin(phcoor(is))+ytrans
            zbal=rdcoor(is)*cos(thcoor(is))+ztrans
            iatomicno=get_atomicnumber(chsymbol(is))
            write(nout,*) iatomicno,real(iatomicno),xbal/bohr,ybal/bohr,zbal/bohr
         end do
         close(nout)
      end if
      !Open the Cube output file. Units need to be in Bohr atomic radii.
      open(nout,file=outfile,status='unknown')
      write(nout,*) '# Cube file output by map3djmol'
      write(nout,*) 'Outer loop x, middle loop y, inner loop z'
      write(nout,*) nspher,xcentre,ycentre,zcentre
      write(nout,*) nx0, xstep*scalefactor/bohr, 0.0, 0.0
      write(nout,*) ny0, 0.0, ystep*scalefactor/bohr, 0.0
      write(nout,*) nz0, 0.0, 0.0, zstep*scalefactor/bohr
      do is=1,nspher
         xbal=rdcoor(is)*sin(thcoor(is))*cos(phcoor(is))+xtrans
         ybal=rdcoor(is)*sin(thcoor(is))*sin(phcoor(is))+ytrans
         zbal=rdcoor(is)*cos(thcoor(is))+ztrans
         iatomicno=get_atomicnumber(chsymbol(is))
         write(nout,*) iatomicno,real(iatomicno),xbal/bohr,ybal/bohr,zbal/bohr
      end do
      do ix=1,nx0
         do iy=1,ny0
            do iz=1,nz0,10
               itmax=min(9,nz0-iz)
               write(nout,'(10(1x,e13.6))') (dens(ix,iy,iz+it),it=0,itmax)
            end do
         end do
      end do
      close(nout)
      if(newstart) then
         open(nout,file=plotfile,status='unknown',iostat=ierr)
      else
         open(nout,file=plotfile,status='old',position='append',iostat=ierr)
      end if
      if(ierr.eq.0) then
         if(newstart) then
            write(nout,795) molfile(1:lenmolfile),rgbbak(1),rgbbak(2),rgbbak(3),translucency
795         format('load "',a,'"; zoom 20; background {'
! Set atom radii to be 15% of VdW radii.
! Put the graph so the z-axis is vertical and the y-axis is horizontal on screen.
     1,f4.1,1x,f4.1,1x,f4.1,'}; spacefill -15; rotate z -90; rotate x -90; set frank off; set defaulttranslucent '
     2,f6.3,';center {0 0 0};set perspectiveDepth off ')
            write(nout,807) rot,elev,spin
807         format('rotate y ',f10.5,'; rotate x ',f10.5,'; spin ',a3)
         end if
         write(nout,797) surfacelabel(1:len_trim(surfacelabel)),dsurf,outfile(1:lenoutfile),rgbobj(1),rgbobj(2),rgbobj(3)
     1,surfacetype(1:len_trim(surfacetype))
797      format('isosurface ',a,1x,f10.5,' "',a,'"; color isosurface {'
     1,f4.1,1x,f4.1,1x,f4.1,'}; isosurface ',a)
      end if
      close(nout)

! Save random number generator

      open(10,file='ran1.dat',status='unknown')
         call saveran1(10)
         close(10)
      end               


