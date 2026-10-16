      function indexlastdot(text)
c
c finds the position of the last dot in text. 
! If there is no dot, the length of text +1 is returned (as if it had a dot)
c
      character(len=*) text

      ilentext=len_trim(text)
      idot=index(text(1:ilentext),'.')
      indexlastdot=idot
      do while (idot.gt.0.and.indexlastdot.lt.ilentext)
         idot=index(text(indexlastdot+1:ilentext),'.')
         indexlastdot=indexlastdot+idot
      end do
      return
      end

      function get_atomicnumber(ats)

      include 'chemical_symbols.inc'

      integer get_atomicnumber,ir
      character(len=3) ats
      logical found
!Search for atomic number of this atom type
      ir=0
      found=.false.
      do while(.not.found.and.ir.lt.melements)
         ir=ir+1
         found=ats.eq.chsymbol(ir)
      end do
      if(found) then
         get_atomicnumber=ir
      else
         get_atomicnumber=0
      end if
      return
      end

      subroutine convertalphabetagamma(alpha,beta,gam,phib,thetac,phic)
! Converts crystallographic alpha, beta, gamma values to spherical polar coordinates used in EPSR
! The a vector is assumed to lie along the x-axis, and the b vector is in the x-y plane. 
      real, intent(in) :: alpha,beta,gam
      real, intent(out) :: phib,thetac,phic
      real :: pi,radconv,degconv,bx,by,cx,cy,bdotc
      real :: alphar,betar,gamr,phicr,thetacr
      pi=4.0*atan(1.0)
      radconv=pi/180.0
      degconv=180.0/pi
      phib=gam
      betar=beta*radconv
      cx=cos(betar)
      if(abs(cx).lt.1.0e-4) then
         thetac=0.0
         phic=0.0
      else
         gamr=gam*radconv
         bx=cos(gamr)
         by=sin(gamr)
         alphar=alpha*radconv
         bdotc=cos(alphar)
         cy=(bdotc-bx*cx)/by
         if(abs(cx).gt.1.0e-4) then
            phicr=atan2(cy,cx)
            phic=phicr*degconv
            cr=sqrt(cx*cx+cy*cy)
            if(abs(cr).le.1.0) then
               thetacr=asin(cr)
               thetac=degconv*thetacr
            else
               write(6,601) cr,phib,phic
601   format('convertalphbetagamma> Incompatible axes definitions: ',3(1x,f10.5))
               return
            end if
         else
            write(6,601) alpha,beta,gam
            return
         endif
      endif
!      write(6,*) alpha,beta,gam
!      write(6,*) phib,thetac,phic
      return
      end

      subroutine convertphibthetacphic(phib,thetac,phic,alpha,beta,gam)
! Converts spherical polar coordinates used in EPSR to crystallographic alpha, beta, gamma values
! The a vector is assumed to lie along the x-axis, and the b vector is in the x-y plane. 
      real, intent(in) :: phib,thetac,phic
      real, intent(out) :: alpha,beta,gam
      real :: pi,radconv,degconv,bx,by,cx,cy,bdotc
      real :: alphar,betar,gamr,phicr,thetacr
      pi=4.0*atan(1.0)
      radconv=pi/180.0
      degconv=180.0/pi
      gam=phib
      phicr=phic*radconv
      thetacr=thetac*radconv
      cxy=sin(thetacr)
      cx=cxy*cos(phicr)
      if(abs(cx).lt.1.0e-4) then
         alpha=90.0
         beta=90.0
      else
         betar=acos(cx)
         beta=betar*degconv
         phibr=phib*radconv
         alphar=acos(cxy*cos(phibr-phicr))
         alpha=alphar*degconv
      end if
!      write(6,*) phib,thetac,phic
!      write(6,*) alpha,beta,gam
      return
      end

      function cosd(deg)
      parameter (degradconv=3.141592653/180.0)
      cosd=cos(deg*degradconv)
      return
      end
      function sind(deg)
      parameter (degradconv=3.141592653/180.0)
      sind=sin(deg*degradconv)
      return
      end
      function acosd(cosphi)
      parameter (degradconv=3.141592653/180.0)
      acosd=acos(cosphi)
      acosd=acosd/degradconv
      return
      end   
      function asind(sinphi)
      parameter (degradconv=3.141592653/180.0)
      asind=asin(sinphi)
      asind=asind/degradconv
      return
      end   
      function atan2d(sinphi,cosphi)
      parameter (degradconv=3.141592653/180.0)
      atan2d=atan2(sinphi,cosphi)
      atan2d=atan2d/degradconv
      return
      end   

      function perpd(a,b,c)
      real*4 a(3),b(3),c(3),t(3)
c
c calculates the perpendicular projection of vector c onto the direction normal
c to the plane defined by vectors a and b
c
c form vector product a X b
c
      call vecpro(a,b,t)
c
c normalise it to a unit vector
c
      amp=sqrt(scapro(t,t))
      do i=1,3
            t(i)=t(i)/amp
      end do
c
c project c onto this normal vector
c
      perpd=scapro(c,t)
      return
      end   

      subroutine vecpro(a,b,c)
      real*4 a(3),b(3),c(3)
c performs the vector product of c1 x c2 and puts result in c3
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end

      function scapro(a,b)
      !Performs dot product a.b
      real*4 a(3),b(3)
      scapro=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      return
      end

      subroutine mulvec(b,a)
      !Multiplies vector b by factor a
      real*4 a,b(3)
      integer i
      do i=1,3
         b(i)=b(i)*a
      end do
      return
      end

      subroutine divvec(b,a)
      !Divide vector b by a if a .ne. 0.0
      real*4 a,b(3)
      integer i
      if(a.ne.0.0) then
         do i=1,3
            b(i)=b(i)/a 
         end do
      endif
      return
      end

      subroutine normvec(a)
      !Creates a unit vector in the direction of a
      integer i
      real a(3),amp
      amp=sqrt(scapro(a,a))
      if(amp.gt.0.0) then
         do i=1,3
            a(i)=a(i)/amp
         end do
      endif
      return
      end

      subroutine initran1(nchan)
      include 'ran1com.inc'
      integer nchan,ierr
c
c open file containing random number parameters, or else initialise
c
            read(nchan,*,iostat=ierr) mcall,idum,iy,(iv(i),i=1,ntab)
c
c if the number of calls to ran1 exceeds 1.0E8, the sequence is reset, as this
c is approximately 1/20th of the period of the ran1 generator
c
c           if(mcall.gt.100000000) then
c                 idum=0
c                 iy=0
c                 mcall=0
c           endif

c if error or end of file occurs, initialise values

            if(ierr.ne.0) then
 
               idum=0
               iy=0
               mcall=0

            endif

            return

      end

      subroutine saveran1(nchan)
      include 'ran1com.inc'
      integer nchan
c
c open file containing random number parameters, or else initialise
c
      write(nchan,*) mcall,idum,iy,(iv(i),i=1,ntab)
      return
      end

      FUNCTION ran1()
      include 'ran1com.inc'
c
c "Minimal" random number generator of Park and Miller with Bays-Durham shue and
c added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
c the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
c alter idum between successive deviates in a sequence. RNMX should approximate the largest
c floating value that is less than 1.
c
c
c Initialize if idum.le.0
c
      if (idum.le.0.or.iy.eq.0) then 
c
c Be sure to prevent idum = 0.
c
            idum=max(-idum,1) 
c
c Load the shue table (after 8 warm-ups).
c
            do j=NTAB+8,1,-1 
                  k=idum/IQ
                  idum=IA*(idum-k*IQ)-IR*k
                  if (idum.lt.0) idum=idum+IM
                  if (j.le.NTAB) iv(j)=idum
            enddo
            iy=iv(1)
      endif
c
c Start here when not initializing.
c
      k=idum/IQ 
c
c Compute idum=mod(IA*idum,IM) without overflows by Schrage's method.
c
      idum=IA*(idum-k*IQ)-IR*k 
      if (idum.lt.0) idum=idum+IM
c
c Will be in the range 1:NTAB.
c
      j=1+iy/NDIV 
      iy=iv(j) 
c
c Output previously stored value and rell the shue table.
c
      iv(j)=idum
c
c Because users don't expect endpoint values.
c
      ran1=min(AM*iy,RNMX) 
      mcall=mcall+1
      return
      END

      subroutine r_int(lpt,fname,nfile,xd,yd,ed,nxdim,nfdim,ierr) !read data
      include 'wrcom.inc'
      real      xd(nxdim),yd(nxdim,nfdim),ed(nxdim,nfdim)
      character*256 fname,sometext
      integer ierr
      real, parameter :: limit=1.0e-15
      nt=(nfile+(nwrmax-1))/nwrmax
      ilen=indexlastdot(fname)+1
      icount=0
      do it=1,nt
c
c number of elements to read for each x value
c
            iwr=nfile-(it-1)*nwrmax
            if(iwr.gt.nwrmax) iwr=nwrmax
            write(ext,331) it
331   format(i2.2)
            fname=fname(1:ilen)//ext(1:2)
            open(8,file=fname,status='old',iostat=ierr)
            if(ierr.ne.0) return
c
c check to see if the first line is a comment line
c
            read(8,100) sometext
100   format(a)
            ilenhash=index(sometext,'#')
            if(ilenhash.le.0) then
c
c hash not found, therefore assume the line is a line of numbers -
c start the file from the beginning
c
               rewind(8)
            endif
c
c READ it
c
            i=1
200         continue
                  read(8,*,end=201) xo(i),(yo(i,ir),eo(i,ir),ir=1,iwr)
! Remove absurdly small values
                  do ir=1,iwr
                     if(abs(yo(i,ir)).lt.limit) yo(i,ir)=0.0
                     if(abs(eo(i,ir)).lt.limit) eo(i,ir)=0.0
                  end do
                  i=i+1
            go to 200
201         lpt=i-1
            write(6,101) lpt,fname(1:ilen+2)
101         format(1x,'Read ',i4,' points from ',a)
c
c save the data for this file
c
            do ir=1,iwr
                  icount=icount+1
                  do i=1,lpt
                        yd(i,icount)=yo(i,ir)
                        ed(i,icount)=eo(i,ir)
                  end do
            end do
            close(8)
      end do
c
c copy the x values
c
      do i=1,lpt
            xd(i)=xo(i)
      end do
      return
      end

      subroutine r_int_titles(fname,nfile,titlex,titled,nfdim,ierr) !read titles in first row
      include 'wrcom.inc'
      character*14 titlex
      character*29 titled(nfdim)
      character*256 fname
      integer ierr,nfile

      integer, parameter                     :: nworddim=201,nstsize=3020
      character(len=nstsize)                 :: sometext
      integer, dimension(nworddim)           :: ncf,ncl
      integer                                :: nwords

      nt=(nfile+(nwrmax-1))/nwrmax
      ilen=indexlastdot(fname)+1
      icount=0
      it=0
      ierr=0
      nwr=nwrmax
      do while(ierr.eq.0.and.nwr.ge.nwrmax) !Keep reading either until an error occurs (file does not exist),
! or the number of titles read for a file is less than nwrmax (A maximum of nwrmax titles can be written for each file).
         it=it+1
         write(ext,'(i2.2)') it
         fname=fname(1:ilen)//ext(1:2)
!         write(6,*) 'r_int_titles> ',fname(1:len_trim(fname))
         open(8,file=fname,status='old',iostat=ierr)
         if(ierr.eq.0) then
c
c check to see if the first line is a comment line
c
            read(8,'(a)') sometext
            if(index(sometext,'#').gt.0) then
               call parse(sometext,nwords,ncf,ncl,20,nstsize,nworddim)
               if(nwords.gt.1) then
                  titlex=sometext(ncf(2):ncl(2))
                  iw=2
                  ir=0
                  do while(iw.lt.nwords.and.ir.lt.nwrmax)
                     iw=iw+1
                     ir=ir+1
                     titleo(ir)=sometext(ncf(iw):ncl(iw))
                  end do
                  nwr=ir
                  ir=0
                  do while(ir.lt.nwr.and.icount.lt.nfdim)
                     ir=ir+1
                     icount=icount+1
                     titled(icount)=titleo(ir)
                  end do
               else
                  ierr=1
               endif
            else
               ierr=1
            end if
            close(8)
         end if
      end do
      nfile=icount
      write(6,*) 'r_int_titles> nfile: ',nfile
      return
      end

      subroutine w_int(lpt,fname,nfile,xd,yd,ed,titlex,titled,nxdim,nfdim)
      include 'wrcom.inc'
      real xd(nxdim),yd(nxdim,nfdim),ed(nxdim,nfdim)
      character fname*256,formatout*30
      character formattitle*30
      character*14 titlex
      character*29 titled(nfdim)
      nt=(nfile+(nwrmax-1))/nwrmax
      ilen=indexlastdot(fname)+1
      icount=0
c
c copy the x values
c
      do i=1,lpt
            xo(i)=xd(i)
      end do
      do it=1,nt
c
c number of elements to write for each x value
c
            iwr=nfile-(it-1)*nwrmax
            if(iwr.gt.nwrmax) iwr=nwrmax
c
c setup output format specification
c
            nblock=2*iwr+1
            write(formattitle,333) '(1h#,1x,a14,',iwr,'(1x,a29))'
333   format(a12,i3.3,a9)
            write(formatout,332) '(',nblock,'(1x,e14.7))'
            lrec=nblock*15+1
332   format(a1,i3.3,a11)
            write(6,*) formatout
            write(ext,331) it
331   format(i2.2)
            fname=fname(1:ilen)//ext(1:2)
            write(6,101) lpt,fname(1:ilen+2)
101         format(1x,'Writing ',i4,' points to ',a)
            open(8,file=fname,status='unknown',form='formatted'
     1,recl=lrec)
c
c get the data for this file
c
            do ir=1,iwr
                  icount=icount+1
                  titleo(ir)=titled(icount)
                  do i=1,lpt
                        yo(i,ir)=yd(i,icount)
                        eo(i,ir)=ed(i,icount)
                  end do
            end do
c
c write it
c
            write(8,fmt=formattitle,iostat=ierr) titlex,(titleo(ir),ir=1,iwr)
            do i=1,lpt
                  write(8,fmt=formatout,iostat=ierr) xo(i),(yo(i,ir),eo(i,ir),ir=1,iwr)
            end do
            close(8)
      end do
100   format(13(1x,e14.7))
      return
      end
 
      subroutine w_intr(lpt1,lpt2,fname,nfile,xd,yd,ed
     *,titlex,titled,nxdim,nfdim)
      include 'wrcom.inc'
      real xd(nxdim),yd(nxdim,nfdim),ed(nxdim,nfdim)
      character fname*256,formatout*30
      character formattitle*30
      character*14 titlex
      character*29 titled(nfdim)
      nt=(nfile+(nwrmax-1))/nwrmax
      ilen=indexlastdot(fname)+1
      icount=0
      lpt=lpt2-lpt1+1
c
c copy the x values
c
      do i=lpt1,lpt2
            xo(i)=xd(i)
      end do
      do it=1,nt
c
c number of elements to write for each x value
c
            iwr=nfile-(it-1)*nwrmax
            if(iwr.gt.nwrmax) iwr=nwrmax
c
c setup output format specification
c
            nblock=2*iwr+1
            write(formattitle,333) '(1h#,1x,a14,',nblock,'(1x,a29))'
333   format(a12,i3.3,a9)
            write(formatout,332) '(',nblock,'(1x,e14.7))'
            lrec=nblock*15+1
332   format(a1,i3.3,a11)
            write(6,*) formatout
            write(ext,331) it
331   format(i2.2)
            fname=fname(1:ilen)//ext(1:2)
            write(6,101) lpt,fname(1:ilen+2)
101         format(1x,'Writing ',i4,' points to ',a)
            open(8,file=fname,status='unknown',form='formatted'
     1,recl=lrec)
c
c get the data for this file
c
            do ir=1,iwr
                  icount=icount+1
                  titleo(ir)=titled(icount)
                  do i=lpt1,lpt2
                        yo(i,ir)=yd(i,icount)
                        eo(i,ir)=ed(i,icount)
                  end do
            end do
c
c write it
c
            write(8,fmt=formattitle) titlex,(titleo(ir),ir=1,iwr)
            do i=lpt1,lpt2
                  write(8,fmt=formatout) xo(i),(yo(i,ir),eo(i,ir),ir=1,iwr)
            end do
            close(8)
      end do
100   format(13(1x,e14.7))
      return
      end

      subroutine rotateeuler(xin,yin,zin,phi,theta,chi,xout,yout,zout)

      real :: xin,yin,zin,phi,chi,xout,yout,zout
      real :: ax,ay,az,bx,by,bz,cx,cy,cz

!      write(6,*) 'rotateeuler> alpha,beta,gamma: ',phi,theta,chi
      calpha=cos(phi)
      salpha=sin(phi)
      cbeta=cos(theta)
      sbeta=sin(theta)
      cgamma=cos(chi)
      sgamma=sin(chi)
! Set up the 3x3 rotation matrix
      ax=calpha*cbeta*cgamma-salpha*sgamma
      bx=-calpha*cbeta*sgamma-salpha*cgamma
      cx=calpha*sbeta
      ay=salpha*cbeta*cgamma+calpha*sgamma
      by=-salpha*cbeta*sgamma+calpha*cgamma
      cy=salpha*sbeta
      az=-sbeta*cgamma
      bz=sbeta*sgamma
      cz=cbeta
      xout=ax*xin+bx*yin+cx*zin
      yout=ay*xin+by*yin+cy*zin
      zout=az*xin+bz*yin+cz*zin

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
      n2=nqmax-1
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
!
! 3-POINT SMOOTHING holding the points nqmin-1, nqmax+1 to zero
!
      n1=nqmin
      n2=nqmax
      if(n2-2.lt.n1.or.nsmoo.eq.0) return
      DO IT=1,NSMOO
         SM=0
         sn=diff(n1)
         sum1=0.25*(sm+sn)
         DO I=N1,N2
            sm=sn
            J=I+1
            if(j.le.n2) then
               SN=DIFF(J)
            else
               sn=0.0
            end if
            sum2=0.25*(sm+sn)
            DIFF(I)=sum1+sum2
            sum1=sum2
         end do
      end do
      RETURN
      END

      function dlm(l,m,n,costh)

      use lngamma_mod

c calculates the rotation matrix dlmn(theta)

      real(kind=8) sth2,cth2,sums,fac,faco,dlcth2,dlsth2,dldif
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
      if(lpm+1.gt.mgamma.or.lmm+1.gt.mgamma.or.lpn+1.gt.mgamma.or.lmn+1.gt.mgamma) then
         write(6,*) 'dlm> ',l,m,n,lpm+1,lmm+1,lpn+1,lmn+1
         stop
      end if
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
      sums=0.
      nsign=(-1)**(kmin-1)
      do k1=kmin,kmax
            k=k1-1
            k2=2*k
            lpmmk=lpm-k+1
            lmnmk=lmn-k+1     
            kmmmn=k-mmn+1
            fac=faco+dlcth2*(l1-k2)+dlsth2*(k2-mmn)-gamln(lpmmk)-gamln(lmnmk)-gamln(kmmmn)-gamln(k1)
            sums=sums+nsign*dexp(fac)
            nsign=nsign*(-1)
      end do
      dlm=sums
      return
      end

      function maxdlm(l,m,n)

c calculates the maximum factorial required for the rotation matrix dlmn

      integer maxdlm,lpm,lmm,lpn,lmn,mmn,kmin,kmax,l1,k,k1,k2,lpmmk,lmnmk,kmmmn
      maxdlm=0
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
      maxdlm=max(lpm+1,lmm+1,lpn+1,lmn+1)
      kmax=min0(lpm,lmn)+1
      kmin=max0(0,mmn)+1
      l1=lpm+lmn
      do k1=kmin,kmax
            k=k1-1
            k2=2*k
            lpmmk=lpm-k+1
            lmnmk=lmn-k+1     
            kmmmn=k-mmn+1
            maxdlm=max(maxdlm,lpmmk,lmnmk,kmmmn,k1)
      end do
      return
      end

      function cgcoef(l1,l2,l,m1,m2)
! Calculates Clebsch-Gordan coefficient for specified indices
      use lngamma_mod

      real cgcoef,aid,e1,e2,sums
      integer l1,l2,l,m1,m2,m,ia,ib,ic,ie,ja,jb,jc,jd,je,jf,kt1,kt2,kt3,kt4,kt5,kmin,kmax,phase,k1,k

      cgcoef=0.
      m=m1+m2
c check for closure relation on l's
      ia=l1+l2-l+1
      if(ia.le.0) return
      ib=l+l1-l2+1
      if(ib.le.0) return
      ic=l+l2-l1+1
      if(ic.le.0) return
      if(abs(m).gt.l) return
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
      e1=0.5*(e1+log(aid))
      kt1=l-l2+m1
      kt2=l-l1-m2
      kmin=max(0,-kt1,-kt2)+1
      kt3=l1+l2-l+1
      kt4=l1-m1+1
      kt5=l2+m2+1
      kmax=min(kt3,kt4,kt5)
      phase=(-1)**(kmin-1)
      sums=0
      do k1=kmin,kmax
            k=k1-1
            e2=gamln(k1)+gamln(kt1+k1)+gamln(kt2+k1)+gamln(kt3-k)+gamln(kt4-k)+gamln(kt5-k)
            sums=sums+real(phase)*exp(e1-e2)
            phase=-phase
      end do
      cgcoef=sums
      return
      end

      function maxcgcoef(l1,l2,l,m1,m2)
! Calculates the maximum factorial needed for specified Clebsch-Gordan coefficient
!      use lngamma_mod

      integer maxcgcoef
      integer l1,l2,l,m1,m2,m,ia,ib,ic,ie,ja,jb,jc,jd,je,jf,kt1,kt2,kt3,kt4,kt5,kmin,kmax,k1,k

      maxcgcoef=0
      m=m1+m2
c check for closure relation on l's
      ia=l1+l2-l+1
      if(ia.le.0) return
      ib=l+l1-l2+1
      if(ib.le.0) return
      ic=l+l2-l1+1
      if(ic.le.0) return
      if(abs(m).gt.l) return
      aid=2*l+1
      ie=l1+l2+l+2
      maxcgcoef=max(maxcgcoef,ia,ib,ic,ie)
!     e1=gamln(ia)+gamln(ib)+gamln(ic)-gamln(ie)
      ja=l1+m1+1
      jb=l1-m1+1
      jc=l2+m2+1
      jd=l2-m2+1
      je=l+m+1
      jf=l-m+1
      maxcgcoef=max(maxcgcoef,jb,jc,jd,je,jf)
!      e1=e1+gamln(ja)+gamln(jb)+gamln(jc)+gamln(jd)+gamln(je)+gamln(jf)
!      e1=0.5*(e1+log(aid))
      kt1=l-l2+m1
      kt2=l-l1-m2
      kmin=max(0,-kt1,-kt2)+1
      kt3=l1+l2-l+1
      kt4=l1-m1+1
      kt5=l2+m2+1
      kmax=min(kt3,kt4,kt5)
      do k1=kmin,kmax
            k=k1-1
            maxcgcoef=max(maxcgcoef,k1,(kt1+k1),(kt2+k1),(kt3-k),(kt4-k),(kt5-k))
!            e2=gamln(k1)+gamln(kt1+k1)+gamln(kt2+k1)+gamln(kt3-k)+gamln(kt4-k)+gamln(kt5-k)
      end do
      return
      end



