      use lngamma_mod
      parameter(mexpan=3000,mval=10,mpixel=80,mpixel2=2*mpixel
     1,mpixel3=mpixel2**3,pi=3.141592653)
      dimension data(500,mexpan),rd(500),err(500,mexpan)
     1,lvl(3,mval),nl(3),ml(3),zcoeff(mexpan),grw(50,50,500)
     1,th(50),ph(50),xl(2),yl(2),temp(500)
     2,out(0:mpixel2,0:mpixel2,0:mpixel2),xx(0:mpixel2),yy(0:mpixel2),l1ex(mexpan)
     3,l2ex(mexpan),lex(mexpan),n1ex(mexpan),n2ex(mexpan)
     *,rdcoor(50),thcoor(50),phcoor(50),radbal(50),rgbbal(3,50)
      character*256 infile,fname
      character*80 titl,longtitle1, longtitle2, longtitle3
      character*40 xcaption2,ycaption2,zcaption2
     1,xcaption3,ycaption3,zcaption3

      character*10 device
      integer ilendevice

      character(len=14)                      :: titlex
      character(len=29), dimension(mexpan)   :: titled
      character(len=29)                      :: temptitle

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
      call fctrln(250)
      write(6,*) 'Type name of file to read plotting parameters from'
      read(5,98) fname

c Read the plotting device

      read(5,'(a)') device
      ilendevice=index(device,' ')-1
      write(6,*) device(1:ilendevice)

      open(10,file=fname,status='old')
      write(6,99)
99    format(1x,'Type name of file containing coefficients')
      nin=10
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
! Get the titles from the input file
      call r_int_titles(infile,nexpan,titlex,titled,nfdim,ierr)
      do k=1,nexpan
         temptitle=titled(k)
         lentemptitle=len(temptitle)
!         write(6,*) k,temptitle
         ifirst=1
         ifirst=index(temptitle(ifirst:lentemptitle),'(')+1
         ilast=ifirst+index(temptitle(ifirst:lentemptitle),',')-2
!         write(6,*) ifirst,ilast
         read(temptitle(ifirst:ilast),*) l1ex(k)
         ifirst=ilast+2
         ilast=ifirst+index(temptitle(ifirst:lentemptitle),',')-2
!         write(6,*) ifirst,ilast
         read(temptitle(ifirst:ilast),*) l2ex(k)
         ifirst=ilast+2
         ilast=ifirst+index(temptitle(ifirst:lentemptitle),';')-2
!         write(6,*) ifirst,ilast
         read(temptitle(ifirst:ilast),*) lex(k)
         ifirst=ilast+2
         ilast=ifirst+index(temptitle(ifirst:lentemptitle),',')-2
!         write(6,*) ifirst,ilast
         read(temptitle(ifirst:ilast),*) n1ex(k)
         ifirst=ilast+2
         ilast=ifirst+index(temptitle(ifirst:lentemptitle),')')-2
!         write(6,*) ifirst,ilast
         read(temptitle(ifirst:ilast),*) n2ex(k)
         write(6,92) k,l1ex(k),l2ex(k),lex(k),n1ex(k),n2ex(k)
92    format(1x,10(1x,i5))
      end do
c
C read in data to be fitted
c
      call r_int(nr,infile,nexpan,rd,data,err,nxdim,nfdim,ierr)
      if(kk0.eq.0) then
         do i=1,nr
            data(i,1)=0.0
            err(i,1)=0.0
         end do
      endif
      do k=1,nexpan
      write(6,92) l1ex(k),l2ex(k),lex(k),n1ex(k),n2ex(k)
      end do
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
      nx=mpixel
      xstep=xmax/float(nx)
      ny=nx
      nz=nx
      ystep=xstep
      zstep=xstep
      write(6,569)
569   format(1x,'Type number of plots in x and y directions on one page')
      read(nin,*) nplotx,nploty
      aspide=float(nploty)/float(nplotx)
      write(6,568)
568   format(1x,'Type aspect ratio (1.4=Portrait, 0.7=Landscape)')
      read(nin,*) aspect
      asprat=aspide/aspect
      CALL PGBEGIN(0,device(1:ilendevice),1,1)
      CALL PGPAPER(0.0,aspect)
      call pgpage
      yvwid=1.0/float(nploty)
      yvwid2=0.5*yvwid
      if(asprat.lt.1.0) yvwid2=yvwid2*asprat
      xvwid=1.0/float(nplotx)
      xvwid2=0.5*xvwid
      if(asprat.gt.1.0) xvwid2=xvwid2/asprat
      do ipy=1,nploty
      ycent=(nploty-ipy+1-0.5)*yvwid
      yv1=ycent-yvwid2
      yv2=ycent+yvwid2
      do ipx=1,nplotx
c
c set up coordinates of centre of viewport
c
      xcent=(ipx-0.5)*xvwid
      xv1=xcent-xvwid2
      xv2=xcent+xvwid2
      write(6,992) xv1,xv2,yv1,yv2
992   format(4(1x,e13.6))
        CALL PGVPORT(xv1,xv2,yv1,yv2)
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
                        temp(i)=0.0
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
c now interpolate onto a cubic grid
c
      iref=0
      do k=0,nz0
          z=(k-nz)*zstep
          z2=z*z
          do j=0,ny0
            y=(j-ny)*ystep
            yy(j)=y
            y2=y*y
            do i=0,nx0
                  x=(i-nx)*xstep    
                  xx(i)=x
                  x2=x*x
                  rxy=sqrt(x2+y2)
                  r=sqrt(x2+y2+z2)
c determine pair of rd values
                  if(r.le.rd(nr).and.r.ge.rd(1)) then
                  ir=nint((r-rd(1))/rstep)+1
			if(ir.lt.2) ir=2
			if(ir.gt.nr-1) ir=nr-1
                  irp=ir+1
                  irm=ir-1
                  rrat=(r-rd(ir))/rstep
                  if(abs(rrat).gt.0.5) then
c     write(6,*) ir,irp,irm,rrat,r,rd(ir),rstep
                  endif
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
                  thrat=theta/tstep-float(nth-1)
                  if(abs(thrat).gt.1.0) then
c     write(6,*) nth,nthp,nthm,thrat,theta,((nth-1)*tstep),tstep
                  endif
c determine nearest phi value
                  phi=0.0
                  if(rxy.gt.0.) then
                        cosphi=x/rxy
                        sinphi=y/rxy
                        phi=atan2d(sinphi,cosphi)
                        if(phi.lt.0.0) phi=phi+360.0
                  endif
                  nph=nint(phi/pstep)+1
                  phrat=phi/pstep-real(nph-1)
! Impose cyclic boundary conditions on the phis
						do while (nph.gt.np0)
							nph=nph-np0
						end do
						do while (nph.lt.1)
							nph=nph+np0
						end do
                  nphp=nph+1
						do while (nphp.gt.np0)
							nphp=nphp-np0
						end do
                  nphm=nph-1
						do while (nphm.lt.1)
							nphm=nphm+np0
						end do
                  if(abs(phrat).gt.0.5) then
c     write(6,*) nth,nphp,nphm,phrat,phi,((nph-1)*pstep),pstep
                  endif
c
c Perform quadratic interpolation
c
                  sum=quadinterp(grw,rrat,thrat,phrat
     1,ir,irp,irm,nth,nthp,nthm,nph,nphp,nphm,50,50,500)
c                 sum=grw(nth,nph,ir)
                  else
                        sum=0.0
                  endif
                  out(i,j,k)=sum
                  iref=iref+1
            end do
          end do
      end do
      nref=iref
      write(6,999) nref
999   format(1x,' No. of points = ',i9)
1001  continue
      write(6,851)
851   format(1x,' Type no. of spheres to plot')
      read(nin,*) nspher
      if(nspher.gt.0) then
            write(6,852)
852   format(1x,' For each sphere type its radius, postion in '
     *,'(r,theta,phi), and three colour indices')
            piconv=3.141592653/180.0
            do is=1,nspher
                  read(nin,*) radbal(is),rdcoor(is),thc
     *,phc,rgbbal(1,is),rgbbal(2,is),rgbbal(3,is)
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
      write(6,853)
853   format(1x,' Type character size (e.g.1.5), line width (e.g.2)'
     *,' and colour index')
      read(nin,*) size,iwidth,icoll
      write(6,793)
793   format(1x,' Type title for this plot')
      read(nin,794) Titl
794   format(a80)
      write(6,795)
795   format(1x,' Type x and y coordinates of title, and char. size')
      read(nin,*) xtr,ytr,csize
      call plot3d(nin,out,nx0,ny0,nz0,size,iwidth,icoll,xmax,nspher,rdcoor
     *,thcoor,phcoor,radbal,rgbbal,xtr,ytr,csize,titl,dsurf)
      call pgwindow(0.0,1.0,0.0,1.0)
      write(6,796)
796   format(1x,'Extra lines? (>0 = number of, 0 = no)')
      read(nin,*) nexlin
      if(nexlin.gt.0) then
            do i=1,nexlin
                  read(nin,*) xl(1),yl(1),xl(2),yl(2) 
                  call pgline(2,xl,yl)
            end do
      endif
      write(6,797)
797   format(1x,'Extra text? (>0 = number of, 0 = no)')
      read(nin,*) nextex
      if(nextex.gt.0) then
            do i=1,nextex
                  read(nin,*) xt1,yt1,csize
                  read(nin,794) titl
                  call pgsch(csize)
                  call pgtext(xt1,yt1,titl)
            end do
      endif
      end do
      end do
      close(10)

c Save random number generator

      open(10,file='randomno.dat',status='unknown')
         call saveran1(10)
        close(10)

      call pgend
      end               

