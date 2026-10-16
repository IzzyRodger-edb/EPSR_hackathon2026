      use lngamma_mod
      parameter(mexpan=1600,mval=10,mxpix=50,mxout=2*mxpix+1,mypix=50,myout=2*mypix+1
     *,nxdim=500,nfdim=mexpan)
      dimension data(nxdim,mexpan),rd(nxdim),err(nxdim,mexpan),temp(nxdim)
     1,lvl(3,mval),nl(3),ml(3),zcoeff(mexpan),grw(200,nxdim),omeg(200)
     2,out(mxout,myout),eout(mxout,myout),xx(mxout),yy(myout),l1ex(mexpan)
     3,l2ex(mexpan),lex(mexpan),n1ex(mexpan),n2ex(mexpan)
     4,xhole(20),yhole(20),rhole(20),grwerr(200,nxdim)
     5,xholet(20),yholet(20),zholet(20),xl(2),yl(2)
      character*256 infile,fname
      character*80 xlbl,ylbl,titl
      character*80 longtitle1, longtitle2, longtitle3
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
      write(6,554)
554   format(1x,'type colour index for background colour')
      read(10,*) iback
      write(6,99)
99    format(1x,'type name of file for SH coefficients')
      read(10,98) infile
98    format(a50)
      write(6,97)
97    format(1x,'type no. of expansion coefficients')
      read(10,*) Nexpan
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
96    format(1x,'Are molecules identical (ident=0 if so, else 1)')
      read(10,*) ident
      write(6,91)
91    format(1x,'Type 0 to set the first coefficient to zero, else 1')
      read(10,*) kk0
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
      write(6,352)nr,infile                     !print message
352   format(' Read ',i5,' points from file ',a)
      write(6,*) 'Type number of smoothings on the data'
      read(10,*) nsmoo
      write(6,101)
101   format(1x,'type x dimension of plot')
      read(10,*) xdim
      write(6,102)
102   format(1x,'type y dimension of plot')
      read(10,*) ydim
      nx=mxpix
      ny=mypix
      xstep=xdim/(2*nx)
      ystep=ydim/(2*ny)
      write(6,104)
104   format(1x,'type minimum radius to include')
      read(10,*) rmin
      write(6,555)
555   format(1x,'type number of holes to draw')
      read(10,*) nhole
      write(6,556)
556   format(1x,'for each hole type its (x,y,z) coordinates and radius')
      read(10,*) (xholet(i),yholet(i),zholet(i),rhole(i),i=1,nhole)
      iterr=0
200   continue
      write(6,569)
569   format(1x,'Type number of plots in x and y directions on one page')
      read(10,*) nplotx,nploty
      aspide=float(nploty)/float(nplotx)
      write(6,568)
568   format(1x,'Type aspect ratio (1.4=Portrait, 0.7=Landscape)')
      read(10,*) aspect
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
      ycent=(nploty-ipy+0.5)*yvwid
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
      xw1=0
      xw2=1
      yw1=0
      yw2=1
      call pgwindow(xw1,xw2,yw1,yw2)
      call pgsci(iback)
      do ic=1,16
            icref=ic-1
            call pgqcr(icref,cr,cg,cb)
            write(6,*) iback,icref,cr,cg,cb
      end do
      call pgrect(xw1,xw2,yw1,yw2)
      write(6,95)
95    format(1x,'Use only l1=0 (0) and/or l2=0 (0)?')
      read(10,*) l1only,l2only
      write(6,93)
93    format(1x,'Use only n1=0 (0) and/or n2=0 (0)?')
      read(10,*) n1only,n2only
      write(6,94)
94    format(1x,'Use only m2=0 (0)?')
      read(10,*) m2fix
      write(6,105)
105   format(1x,'vary thetal (1), phil(2), phim (3), thetam(4)',
     1' or chim (5) for second molecule')
      read(10,*) nvary
c     write(6,102)
c102  format(1x,'type viewing phi')
c     read(10,*) phiv
      phiv=0.0
      if(nvary.eq.1) then
            write(6,106)
106   format(1x,'specify phil, phim, thetam and chim of molecule')
            read(10,*) philin,phim,thetam,chim
      else if(nvary.eq.2) then
            write(6,107)
107   format(1x,'specify thetal, phim, thetam, and chim of molecule')
            read(10,*) thetal,phim,thetam,chim
      else if(nvary.eq.3) then
            write(6,108)
108   format(1x,'specify thetal, phil, thetam, and chim of molecule')
            read(10,*) thetal,philin,thetam,chim
      else if(nvary.eq.4) then
            write(6,109)
109   format(1x,'specify thetal, phil, phim, and chim of molecule')
            read(10,*) thetal,philin,phim,chim
      else 
            write(6,110)
110   format(1x,'specify thetal, phil, phim, and thetam of molecule')
            read(10,*) thetal,philin,phim,thetam
      endif
      if(nvary.eq.1.or.nvary.eq.4) then
            do ih=1,nhole
                  xhole(ih)=xholet(ih)*cosd(philin)+yholet(ih)
     1*sind(philin)
                  yhole(ih)=zholet(ih)
            end do
      else
            do ih=1,nhole
                  xhole(ih)=yholet(ih)
                  yhole(ih)=xholet(ih)
            end do
      endif
c determine minimum and maximum r value needed
      xmax=nx*xstep
      ymax=ny*ystep
      rmax1=sqrt(xmax*xmax+ymax*ymax)
      rmax2=rd(nr)
      rmax=amin1(rmax1,rmax2)
      nrmin=1
      nrmax=1
      do j=1,nr
            if(rd(j).lt.rmin) nrmin=j
            if(rd(j).le.rmax) nrmax=j
      end do
      nx0=2*nx+1
      ny0=2*ny+1
c calculate g(r,w) at 10 deg intervals
      nz0=145
      zstep=2.5
      do j=1,nz0
            z=(j-1)*zstep
            omeg(j)=z
            if(nvary.eq.1) then
                  thetal=z
                  phil=philin
                  phi=phim
                  chi=chim
            else if(nvary.eq.2) then
                  phil=z
                  phi=phim
                  chi=chim
            else if(nvary.eq.3) then
                  phil=philin
                  phi=z
                  chi=chim
            else if(nvary.eq.4) then
                  phil=philin
                  phi=phim
                  thetam=z
                  chi=chim
            else
                  phil=philin
                  phi=phim
                  chi=z
            endif
            if(thetal.gt.180.) then
                  thetal=360-thetal
                  phil=phil+180.
            endif
            ctl=cosd(thetal)
            if(thetam.gt.180.) then
                  thetam=360.-thetam
                  phi=phi+180.
                  chi=chi+180.
            endif
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
     1cosd(m2*(phi-phil)+n2tmp*chi-m1*phil)
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
     1cosd(m2*(phi-phil)+n2tmp*chi-m1*phil)
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
     1cosd(m2*(phi-phil)+n2tmp*chi-m1*phil)
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
     1cosd(m2*(phi-phil)+n2tmp*chi-m1*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                        prod=prodo*termc*phs*phsn
                        sum=sum+prod
                  end do
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
     1cosd(m2*(phi-phil)+n2tmp*chi-m1*phil)
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
     1cosd(m2*(phi-phil)+n2tmp*chi-m1*phil)
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
     1cosd(m2*(phi-phil)+n2tmp*chi-m1*phil)
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
     1cosd(m2*(phi-phil)+n2tmp*chi-m1*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                        prod=prodo*termc*phs*phsn
                        sum=sum+prod
                  end do
                  endif
                  endif
596               continue
                  zcoeff(k)=sum
            end do
c     write(6,837) (zcoeff(k),k=1,nexpan)
837   format(8(1x,f9.2))
            do i=1,nrmin
                  grw(j,i)=0.
            end do
            do i=nrmax,nr
                  grw(j,i)=0.
            end do
            do i=nrmin,nrmax
                  sum=0
                  sume=0.0
                  do k=1,nexpan
                        sum=sum+data(i,k)*zcoeff(k)
                        et=err(k,i)*zcoeff(k)
                        sume=sume+et*et
                  end do
                  temp(i)=sum+1.0
                  grwerr(j,i)=sqrt(sume)
            end do
c
c perform smoothing on coefficients if needed - this is primarily to remove
c hard edges and so avoid fringes in the maps 
c
            if(nsmoo.gt.0) then
                  call smood(1,nr,nsmoo,temp)
            endif
            do i=1,nr
                  grw(j,i)=temp(i)
            end do
      end do
c now interpolate onto a square grid
      iref=0
      do j=1,ny0
            y=(j-ny-1)*ystep
            yy(j)=y
            y2=y*y
            do i=1,nx0
                  x=(i-nx-1)*xstep  
                  xx(i)=x
                  x2=x*x
                  r=sqrt(x2+y2)
                  if(r.lt.rmin) then
                        sum=-1.0
                        go to 500
                  endif
                  xv=x*cosd(phiv)+y*sind(phiv)
                  yv=y*cosd(phiv)-x*sind(phiv)
                  omega=atan2d(xv,yv)
                  if(omega.lt.0.) omega=omega+360.
c                 write(6,257) xv,yv,omega
257   format(3(1x,f10.4))
                  nom1=(omega+0.5*zstep)/zstep+1
c determine pair of rd values
                  do ir=2,nr
                        is=ir
                        if(r.lt.rd(ir)) go to 20
                  end do
20                continue
                  data1=grw(nom1,is)
                  sum=data1
                  error=grwerr(nom1,is)
                  if(iterr.eq.1) then
                        sumt=sum+error
                        if(sum.lt.-1.0) then
                              if(sumt.gt.-1.0) then
                                    sum=-1.0
                              else
                                    sum=sumt
                              endif
                        endif
                  endif
500               continue
c
c determine distance of this point from each hole
c
                  do ih=1,nhole
                        rh=(x-xhole(ih))**2+(y-yhole(ih))**2
                        rh=sqrt(rh)
                        if(rh.le.rhole(ih)) then
                              sum=99.0
                        endif
                  end do
                  iref=iref+1
                  out(i,j)=sum
                  eout(i,j)=error
            end do
      end do
c
c plot result
c
      n1x=nx0
      n2y=ny0
      size=2.0
      iwidth=4
      xlbl='[\A]'
      ylbl='[\A]'
      titl=' '
      write(6,853)
853   format(1x,' Type character size (e.g.1.5), line width (e.g.2)'
     *,' and colour index')
      read(10,*) size,iwidth,icol
      write(6,793)
793   format(1x,' Type title for this plot')
      read(10,794) Titl
794   format(a80)
      write(6,795)
795   format(1x,' Type x and y coordinates of title, and char. size')
      read(10,*) xtr,ytr,csize
      call pgslw(iwidth)
      call plot(xx,mxout,yy,myout,out,mxout,myout,work,size,iwidth,xlbl,ylbl
     1,xtr,ytr,csize,titl,icol,xv1,xv2,yv1,yv2)
      write(6,796)
      call pgwindow(0.0,1.0,0.0,1.0)
796   format(1x,'Extra lines? (>0 = number of, 0 = no)')
      read(10,*) nexlin
      if(nexlin.gt.0) then
            do i=1,nexlin
                  read(10,*) xl(1),yl(1),xl(2),yl(2)  
                  call pgline(2,xl,yl)
            end do
      endif
      write(6,797)
797   format(1x,'Extra text? (>0 = number of, 0 = no)')
      read(10,*) nextex
      if(nextex.gt.0) then
            do i=1,nextex
                  read(10,*) xt1,yt1,csize
                  read(10,794) titl
                  call pgsch(csize)
                  call pgtext(xt1,yt1,titl)
            end do
      endif
      end do
      end do
      close (10)

c Save random number generator

      open(10,file='randomno.dat',status='unknown')
         call saveran1(10)
        close(10)

      call pgend
      end               
