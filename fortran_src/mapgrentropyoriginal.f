      use lngamma_mod
c
c This version of MAP3D puts out the g(r) in specifed directions
c
      parameter(mexpan=3000,mval=10,mpixel=80,mpixel2=(2*mpixel+1)
     1,mpixel3=mpixel2**3)
      dimension data(500,mexpan),rd(500),err(500,mexpan)
     1,lvl(3,mval),nl(3),ml(3),zcoeff(mexpan)
     1,th(50),ph(50),xl(2),yl(2),temp(500)
     2,out(0:mpixel3),xx(201),yy(201),l1ex(mexpan)
     3,l2ex(mexpan),lex(mexpan),n1ex(mexpan),n2ex(mexpan)
     *,rdcoor(50),thcoor(50),phcoor(50),radbal(50),rgbbal(3,50)
      character*256 infile,fname,fnames,fnameout
      character*80 titl,longtitle1, longtitle2, longtitle3
      character*40 xcaption2,ycaption2,zcaption2
     1,xcaption3,ycaption3,zcaption3
      character*1 ext
      character*14 titler
      character*29 titledat(mexpan)

      character(len=14)                      :: titlex
      character(len=29), dimension(mexpan)   :: titled
      character(len=29)                      :: temptitle

      external cosd,sind,acosd,asind,atan2d
		real pi
      integer idot,nfreq
      real(kind=8) entropy,sumr

      pi=4.0*atan(1.0)
      nxdim=500
      nfdim=mexpan
c
c initialise log(factorial)
c
      call fctrln(250)
      open(10,file='ran1.dat',status='unknown')
      call initran1(10)
      close(10)
      write(6,*) 'Type name of file to read parameters from'
      read(5,98) fnames
      idot=indexlastdot(fnames)
      if(idot.gt.0) then
   		fnameout=fnames(1:idot)//'out'
      else
         fnameout=fnames(1:len_trim(fnames))//'.out'
      end if
      open(10,file=fnames,status='old')
      write(6,99)
99    format(1x,'Type name of file containing coefficients')
      read(10,98) infile
98    format(a50)
      write(6,97)
97    format(1x,'Type no. of expansion coefficients')
      read(10,*) nexpan
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
c     write(6,*) ident
      write(6,91)
91    format(1x,'Type 0 to set the first coefficient to zero, else 1')
      read(10,*) kk0
c     write(6,*) kk0
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
      read(10,*) nsmoo
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
      write(6,95)
95    format(1x,'Use only l1=0 (0) and/or l2=0 (0) and/or l = 0 (0)?')
      read(10,*) l1only,l2only,lonly
      write(6,93)
93    format(1x,'Use only n1=0 (0) and/or n2=0 (0)?')
      read(10,*) n1only,n2only
      write(6,94)
94    format(1x,'Use only m2=0 (0)?')
      read(10,*) m2fix
      write(6,105)
105   format(1x,'Specify number of trials to calculate entropy for, and frequency to write')
      read(10,*) nvary,nfreq
121   format(a50)
      close(10)
c calculate g(r,w) at specified angular values
      it=0
      entropy=0.0
      fact=4.0*pi*(rd(3)-rd(2))
      open(10,file=fnameout,status='unknown')
      do while(it.lt.nvary)
            it=it+1
            ctl=2.0*ran1()-1.0
            phil=ran1()*360.0
            phi=ran1()*360.0
            cz=2.0*ran1()-1.0
            chi=ran1()*360.0
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
                  if(lonly.eq.0.and.ll.ne.0) go to 596
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
                        if(abs(m).le.ll) then
                              termc=
     1cosd(m2*phi+n2tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                              prod=prodo*termc
                              sum=sum+prod
                        endif
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
                        if(abs(m).le.ll) then
                              termc=
     1cosd(m2*phi+n2tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                              prod=prodo*termc
                              sum=sum+prod
                        endif
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
                        if(abs(m).le.ll) then
                              termc=
     1cosd(m2*phi+n2tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                              prod=prodo*termc
                              sum=sum+prod
                        endif
                  end do
                  n1tmp=-n1tmp
                  n2tmp=-n2tmp
                  m1=n1tmp
                  do m2t=n2f,n2l
                        m2=m2t-n2m
                        prodo=dlm(l2tmp,m2,n2tmp,cz)
                        m=m1+m2
                        if(abs(m).le.ll) then
                              termc=
     1cosd(m2*phi+n2tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                              prod=prodo*termc
                              sum=sum+prod
                        endif
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
                        if(abs(m).le.ll) then
                              termc=
     1cosd(m2*phi+n2tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                              prod=prodo*termc
                              sum=sum+prod
                        endif
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
                        if(abs(m).le.ll) then
                              termc=
     1cosd(m2*phi+n2tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                              prod=prodo*termc
                              sum=sum+prod
                        endif
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
                        if(abs(m).le.ll) then
                              termc=
     1cosd(m2*phi+n2tmp*chi-m*phil)
     1*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                              prod=prodo*termc
                              sum=sum+prod
                        endif
                  end do
                  n1tmp=-n1tmp
                  n2tmp=-n2tmp
                  m1=n1tmp
                  do m2t=n2f,n2l
                        m2=m2t-n2m
                        prodo=dlm(l2tmp,m2,n2tmp,cz)
                        m=m1+m2
                        if(abs(m).le.ll) then
                              termc=cosd(m2*phi+n2tmp*chi-m*phil)*cgcoef(l1tmp,l2tmp,ll,m1,m2)*dlm(ll,m,0,ctl)
                              prod=prodo*termc
                              sum=sum+prod
                        endif
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
                  temp(i)=sum
            end do
            if(nsmoo.gt.0) then
                  call smood(1,nr,nsmoo,temp)
            endif
            sumr=0.0
            do i=1,nr
                  sum=temp(i)
                  term=1.0+sum
                  if(term.gt.0.0) then
                     contribution=sum-term*log(term)
                  else
                     contribution=-1.0
                  end if
                  sumr=sumr+fact*rd(i)*rd(i)*contribution
            end do
            entropy=entropy+sumr
            if(it.eq.nfreq*int(it/nfreq)) then
               write(6,'(5(e13.6,1x))') ctl,phil,phi,cz,chi
               write(6,'(a,i10,a,2(1x,e13.6))') 'mapgrentropy> Entropy after ',it,' trials:',sumr,entropy/real(it)
               write(10,'(a,i10,a,2(1x,e13.6))') 'mapgrentropy> Entropy after ',it,' trials:',sumr,entropy/real(it)
            end if
      end do
      close(10)
c
c write out the results
c

c set the axis title

!      titler='r[A]'
!      ilen=index(fnames,'.')
!      ext='g'
!      fname=fnames(1:ilen)//ext
!      call w_int(nr,fname,it,rd,grw,err,titler,titledat,nxdim,nfdim)
      open(10,file='ran1.dat',status='unknown')
      call saveran1(10)
      close(10)
      end               

