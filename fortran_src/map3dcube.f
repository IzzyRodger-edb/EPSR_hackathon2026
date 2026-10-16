      subroutine map3dcube(infile,plotfile)
		
		parameter(mexpan=3000,mval=10,pi=3.141592653)

      real data(500,mexpan),rd(500),err(500,mexpan)
     1,lvl(3,mval),nl(3),ml(3),zcoeff(mexpan),grw(50,50,500)
     1,th(50),ph(50),xl(2),yl(2),temp(500)
!     2,xx(201),yy(201)
     *,rdcoor(50),thcoor(50),phcoor(50),radbal(50),rgbbal(3,50)
     *,rot,elev,rgbobj(3),rgbbak(3),xlight,ylight,zlight

      real, dimension(:), allocatable :: xpix,ypix,zpix
      real, dimension(:,:,:), allocatable :: dens
      integer, dimension(:,:,:), allocatable :: neighb

!      real dens(0:mpixel2,0:mpixel2,0:mpixel2)
      
!      integer neighb(0:mpixel2,0:mpixel2,0:mpixel2)
 
      integer l1ex(mexpan)
     3,l2ex(mexpan),lex(mexpan),n1ex(mexpan),n2ex(mexpan)

      character*256 infile,fname,outfile,plotfile
      character*80 titl,longtitle1, longtitle2, longtitle3,surfacetype,mesh,transparency
      character*40 xcaption2,ycaption2,zcaption2
     1,xcaption3,ycaption3,zcaption3

      character*10 device
      integer ilendevice,nref,nvoxel,mpixel,mpixel2,imesh,itrans

      external cosd,sind,acosd,asind,atan2d

c Initialise random number generator

      open(10,file='randomno.dat',status='unknown')
         call initran1(10)
        close(10)

      nxdim=500
      nfdim=mexpan
c
c initialise log(factorial)
c
      call fctrln
      write(6,*) 'Type name of file to read plotting parameters from'
      read(5,98) fname
! Name of output file - remove any text after the first dot
      lenfirstdot=index(fname,'.')-1
      if(lenfirstdot.gt.0) then
         outfile=fname(1:lenfirstdot)//'.cube'
      else
         outfile=fname(1:len_trim(fname))//'.cube'
      end if
      lenoutfile=len_trim(outfile)
      plotfile='plot'//fname(1:lenfirstdot)//'.txt'
      nin=10
      open(nin,file=fname,status='old')
      write(6,99)
99    format(1x,'Type name of file containing coefficients')
      read(nin,98) infile
98    format(a)
      write(6,97)
97    format(1x,'Type no. of expansion coefficients')
      read(nin,*) Nexpan
c
c clean out data arrays
c
      do k=1,nexpan
            do i=1,500
                  data(i,k)=0.0
                  err(i,k)=0.0
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
91    format(1x,'Type 0 to set the first coefficient to zero, else 1')
      read(nin,*) kk0
982   format(1x,4(1x,i5),f12.7)
c
C read in data to be fitted
c
      call r_int(nr,infile,nexpan,rd,data,err,nxdim,nfdim,ierr)
      do k=1,nexpan
            l1ex(k)=nint(data(1,k))
            l2ex(k)=nint(data(2,k))
            lex(k)=nint(data(3,k))
            n1ex(k)=nint(data(4,k))
            n2ex(k)=nint(data(5,k))
            if(k.eq.1.and.kk0.eq.0) then
                  do i=1,nr
                        data(i,k)=0.0
                        err(i,k)=0.0
                  end do
            endif
      end do
      do k=1,nexpan
      write(6,92) l1ex(k),l2ex(k),lex(k),n1ex(k),n2ex(k)
      end do
92    format(1x,10(1x,i5))
      write(6,*) 'Type number of smoothings on the data'
      read(nin,*) nsmoo
      if(nsmoo.gt.0) then
      do k=1,nexpan
            do i=1,nr
                  temp(i)=data(i,k)
            end do
c           call smood(1,nr,nsmoo,temp)
            do i=1,nr
                  data(i,k)=temp(i)
            end do
      end do
      endif
      write(6,101)
101   format(1x,'Type rmax')
      read(nin,*) xmax
      rstep=rd(3)-rd(2)
      write(6,569)
569   format(1x,'Type number of pixels on each axis')
      read(nin,*) mpixel
      nx=mpixel
      xstep=xmax/float(nx)
      ny=nx
      nz=nx
      ystep=xstep
      zstep=xstep
      mpixel2=2*mpixel
      if(allocated(dens)) deallocate(dens,neighb)
      allocate(dens(0:mpixel2,0:mpixel2,0:mpixel2),neighb(0:mpixel2,0:mpixel2,0:mpixel2))
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
      if(nspher.gt.0) then
            write(6,852)
852   format(1x,' For each sphere type its radius, postion in '
     *,'(r,theta,phi), and atomicno')
            piconv=3.141592653/180.0
            do is=1,nspher
                  read(nin,*) radbal(is),rdcoor(is),thc
     *,phc,rgbbal(1,is)
                  thcoor(is)=thc*piconv
                  phcoor(is)=phc*piconv
            end do
            write(6,*) phimmax,thetammax,chimmax,dhigh
            phirot=phimmax
            thetarot=thetammax
            chirot=chimmax

c For NVARY = 3 we need to rotate the central molecule about the z-axis by chirot degrees
c because that is the way the orientational spatial density function is plotted, namely the
c distribution out of the z-x plane represents the variation with chim.

            if(nvary.eq.3) phirot=phirot+chirot

            write(6,*) phirot,thetarot,chirot,dhigh
            phirot=phirot*piconv
            thetarot=thetarot*piconv
            chirot=chirot*piconv
c
c apply the rotation. Do this in reverse order about the lab. axes
c
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
487	format(1x,' Type 3 RGB Colour indices for background')
	   read(nin,*) rgbbak(1),rgbbak(2),rgbbak(3)
!	   write(6,485)
!485	format(1x,' Type ishade (1-8): 0 means no shading')
!	   read(nin,*) ishade
	   write(6,489)
489	format(1x,' Type 3 RGB Colour indices for object')
	   read(nin,*) rgbobj(1),rgbobj(2),rgbobj(3)
      write(6,*) rgbobj(1),rgbobj(2),rgbobj(3)
!	   write(6,488)
!488	format(1x,' Type 3 (x,y,z) coordinates for light source')
!	   read(nin,*) xli,yli,zli
!      write(6,*) xli,yli,zli
!	   write(6,486)
!486	format(1x,' Type fade factor (0 = no fading, 1=full fading)')
!	   read(nin,*) fadefc
	   write(6,484)
484	format(1x,' Transparency of surface (0=opaque,1=translucent,2=transparent)')
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
490	format(1x,' Mesh? (0=nothing,1=mesh,2=dots)')
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
492   format(1x,' Type rotation and tilt')
      read(nin,*) rot, elev
      close(nin)
      nx0=2*nx
      ny0=2*ny
      nz0=2*nz
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
      do it=1,nt0
          t=(it-1)*tstep
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
                  ll1=l1ex(k)
                  ll2=l2ex(k)
                  ll=lex(k)
                  n1=n1ex(k)
                  n2=n2ex(k)
                  sum=0.
                  l1tmp=ll1
                  l2tmp=ll2
                  n1tmp=n1
                  n2tmp=n2
                  if(l1only.eq.0.and.l1tmp.ne.0) go to 595
                  if(l2only.eq.0.and.l2tmp.ne.0) go to 595
                  if(n1only.eq.0.and.n1tmp.ne.0) go to 595
                  if(n2only.eq.0.and.n2tmp.ne.0) go to 595
                  phsn=((-1)**(l1tmp+l2tmp+ll+n1tmp+n2tmp))

c Normal algorithm for molecule 1 at the origin

                  if(identsign.ge.0) then

                     if(m2fix.gt.0) then
                        n2f=1
                        n2m=l2tmp+1
                        n2l=2*n2m-1
                     else
                        n2f=1
                        n2m=1
                        n2l=1
                     endif

                     m1=n1tmp
                     do m2t=n2f,n2l
                        m2=m2t-n2m
                        prodo=dlm(l2tmp,m2,n2tmp,cz)
                        m=m1+m2
                        termc=
     1cosd(m2*phi+n2tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                        prod=prodo*termc
                        sum=sum+prod
                     end do
c now reverse signs of n1,n2 if either is non-zero
                     if(n1.gt.0.or.n2.gt.0) then
                        n1tmp=-n1tmp
                        n2tmp=-n2tmp
                        m1=n1tmp
                        do m2t=n2f,n2l
                           m2=m2t-n2m
                           prodo=dlm(l2tmp,m2,n2tmp,cz)
                           m=m1+m2
                           termc=
     1cosd(m2*phi+n2tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                           prod=prodo*termc*phsn
                           sum=sum+prod
                        end do
                     endif
c if BOTH n1 AND n2 are greater than zero we must calculate two more terms
                     if(n1.gt.0.and.n2.gt.0) then
                        n1tmp=-n1tmp
                        n2tmp=n2tmp
                        phs=(-1)**(-n2tmp)
                        m1=n1tmp
                        do m2t=n2f,n2l
                           m2=m2t-n2m
                           prodo=dlm(l2tmp,m2,n2tmp,cz)
                           m=m1+m2
                           termc=
     1cosd(m2*phi+n2tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                           prod=prodo*termc*phs
                           sum=sum+prod
                        end do
                        n1tmp=-n1tmp
                        n2tmp=-n2tmp
                        m1=n1tmp
                        do m2t=n2f,n2l
                           m2=m2t-n2m
                           prodo=dlm(l2tmp,m2,n2tmp,cz)
                           m=m1+m2
                           termc=
     1cosd(m2*phi+n2tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                           prod=prodo*termc*phs*phsn
                           sum=sum+prod
                        end do
                     endif
                  else
                     if(m2fix.gt.0) then
                        n2f=1
                        n2m=l1tmp+1
                        n2l=2*n2m-1
                     else
                        n2f=1
                        n2m=1
                        n2l=1
                     endif

c Phase factor to reverse parity of centre to centre vector

                     phss=(-1)**ll

                     m2=n2tmp
                     do m2t=n2f,n2l
                        m1=m2t-n2m
                        prodo=dlm(l1tmp,m1,n1tmp,cz)
                        m=m1+m2
                        termc=
     1cosd(m1*phi+n1tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                        prod=prodo*termc*phss
                        sum=sum+prod
                     end do
c now reverse signs of n1,n2 if either is non-zero
                     if(n1.gt.0.or.n2.gt.0) then
                        n1tmp=-n1tmp
                        n2tmp=-n2tmp
                        m2=n2tmp
                        do m2t=n2f,n2l
                           m1=m2t-n2m
                           prodo=dlm(l1tmp,m1,n1tmp,cz)
                           m=m1+m2
                           termc=
     1cosd(m1*phi+n1tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                           prod=prodo*termc*phss*phsn
                           sum=sum+prod
                        end do
                     endif
c if BOTH n1 AND n2 are greater than zero we must calculate two more terms
                     if(n1.gt.0.and.n2.gt.0) then
                        n1tmp=-n1tmp
                        n2tmp=n2tmp
                        phs=(-1)**(-n2tmp)
                        m2=n2tmp
                        do m2t=n2f,n2l
                           m1=m2t-n2m
                           prodo=dlm(l1tmp,m1,n1tmp,cz)
                           m=m1+m2
                           termc=
     1cosd(m1*phi+n1tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                           prod=prodo*termc*phss*phs
                           sum=sum+prod
                        end do
                        n1tmp=-n1tmp
                        n2tmp=-n2tmp
                        m2=n2tmp
                        do m2t=n2f,n2l
                           m1=m2t-n2m
                           prodo=dlm(l1tmp,m1,n1tmp,cz)
                           m=m1+m2
                           termc=
     1cosd(m1*phi+n1tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                           prod=prodo*termc*phss*phs*phsn
                           sum=sum+prod
                        end do
                     endif
                  endif
c now repeat for interchange of l1,n1,with l2,n2 for identical molecules
c bearing in mind the phase convention on the h's - the case when l1=l2 AND
c n1=n2 doesn't count here. 
c
595               continue
                  if(ident.eq.0.and.((ll1.ne.ll2).or.(n1.ne.n2)))
     1then
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
                        n2f=1
                        n2m=l2tmp+1
                        n2l=2*n2m-1
                  else
                        n2f=1
                        n2m=1
                        n2l=1
                  endif
                  m1=n1tmp
                  do m2t=n2f,n2l
                        m2=m2t-n2m
                        prodo=dlm(l2tmp,m2,n2tmp,cz)
                        m=m1+m2
                        termc=
     1cosd(m2*phi+n2tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                        prod=prodo*termc*phs
                        sum=sum+prod
                  end do
c now reverse signs of n1,n2 if either is non-zero
                  if(n1.gt.0.or.n2.gt.0) then
                  n1tmp=-n1tmp
                  n2tmp=-n2tmp
                  m1=n1tmp
                  do m2t=n2f,n2l
                        m2=m2t-n2m
                        prodo=dlm(l2tmp,m2,n2tmp,cz)
                        m=m1+m2
                        termc=
     1cosd(m2*phi+n2tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                        prod=prodo*termc*phs*phsn
                        sum=sum+prod
                  end do
                  endif
c if BOTH n1 AND n2 are greater than zero we must calculate two more terms
                  if(n1.gt.0.and.n2.gt.0) then
                  n1tmp=-n1tmp
                  n2tmp=n2tmp
                  phs=phs*(-1)**(-n2tmp)
                  m1=n1tmp
                  do m2t=n2f,n2l
                        m2=m2t-n2m
                        prodo=dlm(l2tmp,m2,n2tmp,cz)
                        m=m1+m2
                        termc=
     1cosd(m2*phi+n2tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                        prod=prodo*termc*phs
                        sum=sum+prod
                  end do
                  n1tmp=-n1tmp
                  n2tmp=-n2tmp
                  m1=n1tmp
                  do m2t=n2f,n2l
                        m2=m2t-n2m
                        prodo=dlm(l2tmp,m2,n2tmp,cz)
                        m=m1+m2
                        termc=
     1cosd(m2*phi+n2tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                        prod=prodo*termc*phs*phsn
                        sum=sum+prod
                  end do
                  endif
                  endif
596               continue
                  zcoeff(k)=sum
            end do
            do i=1,nr
                  sum=0
                  do k=1,nexpan
                        sum=sum+data(i,k)*zcoeff(k)
                  end do
c
c test for maximum and minimum in selected r range
c
                  r=rd(i)
                  if(r.ge.rmin.and.r.le.rmax) then
                        temp(i)=sum+1.0
                  else
                        temp(i)=1.0
                  endif
            end do
c
c perform smoothing on coefficients if needed - this is primarily to remove
c hard edges and so avoid fringes in the maps 
c
            if(nsmoo.gt.0) then
                  call smood(1,nr,nsmoo,temp)
            endif
            do i=1,nr
                  grw(it,ip,i)=temp(i)
                  if(grw(it,ip,i).lt.dlow) dlow=grw(it,ip,i)
                  if(grw(it,ip,i).gt.dhigh) then
                     dhigh=grw(it,ip,i)
                     if(rd(i).ge.rmin) then
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
c                 test=data(1,i)+1.0
c                 sumt=sumt+r2*test
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
c
c now interpolate onto a cubic grid and determine whether pixel is occupied or not
c
      iref=0
      do k=0,nz0
         z=(k-nz)*zstep
         z2=z*z
         do j=0,ny0
            y=(j-ny)*ystep
!            yy(j)=y
            y2=y*y
            do i=0,nx0
               x=(i-nx)*xstep    
!               xx(i)=x
               x2=x*x
               rxy=sqrt(x2+y2)
               r=sqrt(x2+y2+z2)
c determine pair of rd values
               if(r.le.rd(nr)) then
                  ir=max(2,(nint(r/rstep)+1)) !This is to prevent irm referencing undefined grw at low r
                  ir=min(nr-1,ir) !This is to prevent irp referencing undefined grw at high r
                  irp=ir+1
                  irm=ir-1
                  rrat=(r-rd(ir))/rstep
!               if(abs(rrat).gt.0.5) then
c     write(6,*) ir,irp,irm,rrat,r,rd(ir),rstep
!               endif
c determine nearest theta value
                  if(r.gt.0.0) then
                     costh=z/r
                     theta=acosd(costh)
                  else
                     theta=0.0
                  endif
                  nth=nint(theta/tstep)+1
                  if(nth.lt.2) nth=2
                  if(nth.ge.nt0) nth=nt0-1
                  nthp=nth+1
                  nthm=nth-1
                  thrat=theta/tstep-float(nth-1)
!               if(abs(thrat).gt.1.0) then
c     write(6,*) nth,nthp,nthm,thrat,theta,((nth-1)*tstep),tstep
!                  endif
c determine nearest phi value
                  if(rxy.gt.0.) then
                     cosphi=x/rxy
                     sinphi=y/rxy
                     phi=atan2d(sinphi,cosphi)
                     if(phi.lt.0.0) phi=phi+360.0
                  else
                     phi=0.0
                  endif
                  nph=nint(phi/pstep)+1
                  nphs=nph
                  if(nph.gt.np0) nph=1
                  nphp=nph+1
                  if(nphp.gt.np0) nphp=1
                  nphm=nph-1
                  if(nphm.lt.1) nphm=np0
                  phrat=phi/pstep-float(nphs-1)
!               if(abs(phrat).gt.0.5) then
c     write(6,*) nth,nphp,nphm,phrat,phi,((nph-1)*pstep),pstep
!                  endif
c
c Perform quadratic interpolation
c
                  sum=quadinterp(grw,rrat,thrat,phrat,ir,irp,irm,nth,nthp,nthm,nph,nphp,nphm,50,50,500)
c                    sum=grw(nth,nph,ir)
!                  if(r.gt.2.0) then
!                     sum=1.0
!                  else
!                     sum=10.0
!                  end if
               else
                  sum=real(kk0)
               endif
               !List a pixel as occupied or not occupied. Pixels on the edge of the box are always unoccupied
               dens(i,j,k) = sum
               if(sum.ge.dsurf) iref=iref+1
               !Initialise the number of neighbours for this cell
               neighb(i,j,k)=0
            end do
         end do
      end do
      nvoxel=(nx0+1)*(ny0+1)*(nz0+1)
      nref=iref
      write(6,999) nref,nvoxel
999   format(1x,' No. of occupied voxels = ',i9,' out of ',i9,' total.')
      nout=nin+1
      !Open the Cube output file. Units need to be in Bohr atomic radii.
      bohr=0.52917721092
      open(nout,file=outfile,status='unknown')
      write(nout,*) '# Cube file output by map3dcube'
      write(nout,*) 'Outer loop z, middle loop x, inner loop z'
      xcentre=-nx*xstep/bohr
      ycentre=-ny*ystep/bohr
      zcentre=-nz*zstep/bohr
      write(nout,*) nspher,xcentre,ycentre,zcentre
      write(nout,*) nx0+1, xstep/bohr, 0.0, 0.0
      write(nout,*) ny0+1, 0.0, ystep/bohr, 0.0
      write(nout,*) nz0+1, 0.0, 0.0, zstep/bohr
      do is=1,nspher
         xbal=rdcoor(is)*sin(thcoor(is))*cos(phcoor(is))
         ybal=rdcoor(is)*sin(thcoor(is))*sin(phcoor(is))
         zbal=rdcoor(is)*cos(thcoor(is))
         write(nout,*) int(rgbbal(1,is)),rgbbal(1,is),xbal/bohr,ybal/bohr,zbal/bohr
      end do
      do ix=0,nx0
         do iy=0,ny0
            do iz=0,nz0,10
               itmax=min(9,nz0-iz)
               write(nout,*) (dens(ix,iy,iz+it),it=0,itmax)
            end do
         end do
      end do
      close(nout)
      open(nout,file=plotfile,status='unknown',iostat=ierr)
      if(ierr.eq.0) then
         write(nout,795) outfile(1:lenoutfile),rgbbak(1),rgbbak(2),rgbbak(3)
795   format(1x,'load "',a,'"; zoom 20; background {'
! Set atom radii to be 60% of VdW radii.
! Put the graph so the z-axis is vertical and the y-axis is horizontal.
     1,f4.1,1x,f4.1,1x,f4.1,'}; cpk -60; rotate z -90; rotate x -90')
         write(nout,797) dsurf,outfile(1:lenoutfile),rgbobj(1),rgbobj(2),rgbobj(3)
     1,surfacetype(1:len_trim(surfacetype)),rot,elev
797   format(1x,'isosurface ',f10.5,' "',a,'"; color isosurface {'
     1,f4.1,1x,f4.1,1x,f4.1,'}; isosurface ',a
     1,/' rotate y ',f10.5,'; rotate x ',f10.5)
      end if
		close(nout)

! Save random number generator

      open(10,file='randomno.dat',status='unknown')
         call saveran1(10)
      	close(10)
      end               

c calculates the rotation matrix Dlmn(theta)
      function dlm(l,m,n,costh)
      double precision sth2,cth2,sum,fac,fac0,dlcth2,dlsth2,dldif,gamln
        COMMON/L/GAMLN(350)
      dlm=0.
      l2=2*l
      lpm=l+m
      if(lpm.gt.l2.or.lpm.lt.0) return
      lmm=l-m
      if(lmm.gt.l2.or.lmm.lt.0) return
      lpn=l+n
      if(lpn.gt.l2.or.lpn.lt.0) return
      lmn=l-n
      if(lmn.gt.l2.or.lmn.lt.0) return
      mmn=m-n
      faco=0.5*(gamln(lpm+1)+gamln(lmm+1)+gamln(lpn+1)+gamln(lmn+1))
      kmax=min0(lpm,lmn)+1
      kmin=max0(0,mmn)+1
      if(costh.ge.1.) then
            cth2=1.0
            dlcth2=0.0
            sth2=0.0
            dlsth2=-1.0d30
      else if(costh.le.-1.) then
            cth2=0.0
            dlcth2=-1.0d30
            sth2=1.0
            dlsth2=0.0
      else
            cth2=sqrt((1.+costh)*0.5)
            dlcth2=dlog(cth2)
            sth2=sqrt((1.-costh)*0.5)
            dlsth2=dlog(sth2)
      endif
      l1=lpm+lmn
      dldif=dlsth2-dlcth2
      sum=0.
      nsign=(-1)**(kmin-1)
      do k1=kmin,kmax
            k=k1-1
            k2=2*k
            lpmmk=lpm-k+1
            lmnmk=lmn-k+1     
            kmmmn=k-mmn+1
            fac=faco+dlcth2*(l1-k2)+dlsth2*(k2-mmn)
     1-gamln(lpmmk)-gamln(lmnmk)-gamln(kmmmn)-gamln(k1)
            sum=sum+nsign*dexp(fac)
            nsign=nsign*(-1)
      end do
      dlm=sum
      return
      end
      function cgcoef(l1,l2,l,m1,m2)
      common/l/gamln(350)
      real*8 gamln,e1,e2,sum
      cgcoef=0.
      m=m1+m2
      ia=l1+l2-l+1
      ib=l+l1-l2+1
      ic=l+l2-l1+1
c check for closure relation on l's
      if(ia) 150,150,100
100   if(ib) 150,150,110
110   if(ic) 150,150,120
120   continue
      if(iabs(m).gt.l) go to 150
      aid=2*l+1
      ie=l1+l2+l+2
      e1=gamln(ia)+gamln(ib)+gamln(ic)-gamln(ie)
      ja=l1+m1+1
      jb=l1-m1+1
      jc=l2+m2+1
      jd=l2-m2+1
      je=l+m+1
      jf=l-m+1
      e1=e1+gamln(ja)+gamln(jb)+gamln(jc)+gamln(jd)+gamln(je)+gamln(jf)
      e1=0.5*(e1+alog(aid))
      kt1=l-l2+m1
      kt2=l-l1-m2
      kmin=max0(0,-kt1,-kt2)+1
      kt3=l1+l2-l+1
      kt4=l1-m1+1
      kt5=l2+m2+1
      kmax=min0(kt3,kt4,kt5)
      phase=(-1)**(kmin-1)
      sum=0
      do k1=kmin,kmax
            k=k1-1
            e2=gamln(k1)+gamln(kt1+k1)+gamln(kt2+k1)
     c+gamln(kt3-k)+gamln(kt4-k)+gamln(kt5-k)
            sum=sum+phase*dexp(e1-e2)
            phase=-phase
      end do
      cgcoef=sum
150   return
      end
       SUBROUTINE FCTRLN                                                                  
       COMMON/L/GAMLN(350)                                                                
      double precision gamln,a,g
       G=0.0                                                                              
      GAMLN(1)=0.0                                                                        
      A=1.0                                                                               
      DO 100 I=2,350                                                                      
      GAMLN(I)=G+dLOG(A)                                                                  
      A=A+1.0                                                                             
      G=GAMLN(I)                                                                          
  100 CONTINUE                                                                            
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
      function quadinterp(grw,rrat,thrat,phrat,nro,nrp,nrm,ntho,nthp
     1,nthm,npho,nphp,nphm,ntdim,npdim,nrdim)
      real*4 grw(ntdim,npdim,nrdim)
c
c quadratic interpolation along r
c
      delta=rrat
      dp=1.0+delta
      dm=1.0-delta
      dpm=dm*dp
      delta2=0.5*delta
      amm=delta2*(dp*grw(nthm,nphm,nrp)-dm*grw(nthm,nphm,nrm))
     1+dpm*grw(nthm,nphm,nro)
      amo=delta2*(dp*grw(nthm,npho,nrp)-dm*grw(nthm,npho,nrm))
     1+dpm*grw(nthm,npho,nro)
      amp=delta2*(dp*grw(nthm,nphp,nrp)-dm*grw(nthm,nphp,nrm))
     1+dpm*grw(nthm,nphp,nro)
      aom=delta2*(dp*grw(ntho,nphm,nrp)-dm*grw(ntho,nphm,nrm))
     1+dpm*grw(ntho,nphm,nro)
      aoo=delta2*(dp*grw(ntho,npho,nrp)-dm*grw(ntho,npho,nrm))
     1+dpm*grw(ntho,npho,nro)
      aop=delta2*(dp*grw(ntho,nphp,nrp)-dm*grw(ntho,nphp,nrm))
     1+dpm*grw(ntho,nphp,nro)
      apm=delta2*(dp*grw(nthp,nphm,nrp)-dm*grw(nthp,nphm,nrm))
     1+dpm*grw(nthp,nphm,nro)
      apo=delta2*(dp*grw(nthp,npho,nrp)-dm*grw(nthp,npho,nrm))
     1+dpm*grw(nthp,npho,nro)
      app=delta2*(dp*grw(nthp,nphp,nrp)-dm*grw(nthp,nphp,nrm))
     1+dpm*grw(nthp,nphp,nro)
c
c now interpolate in theta
c
      delta=thrat
      dp=1.0+delta
      dm=1.0-delta
      dpm=dm*dp
      delta2=0.5*delta
      am=delta2*(dp*apm-dm*amm)+dpm*aom
      ao=delta2*(dp*apo-dm*amo)+dpm*aoo
      ap=delta2*(dp*app-dm*amp)+dpm*aop
c
c now interpolate in phi
c
      delta=phrat
      dp=1.0+delta
      dm=1.0-delta
      dpm=dm*dp
      delta2=0.5*delta
      quadinterp=delta2*(dp*ap-dm*am)+dpm*ao
      return
      end
