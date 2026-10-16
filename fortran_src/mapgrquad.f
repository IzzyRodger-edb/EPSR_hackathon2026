      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_mapgr()

      stop
      end      

      subroutine run_mapgr()

      use lngamma_mod
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
c
c This version of MAP3D puts out the g(r) in specifed directions
c
      parameter(mexpan=3000,mval=10,mpixel=80,mpixel2=(2*mpixel+1)
     1,mpixel3=mpixel2**3)
      dimension datac(500,mexpan),datas(500,mexpan),rd(500)
     1,lvl(3,mval),nl(3),ml(3),zcoeffc(mexpan),zcoeffs(mexpan),grw(500,mexpan),grwerr(500,mexpan)
     1,th(50),ph(50),xl(2),yl(2),temp(500)
     2,out(0:mpixel3),xx(201),yy(201),l1ex(mexpan)
     3,l2ex(mexpan),lex(mexpan),n1ex(mexpan),n2ex(mexpan)
     *,rdcoor(50),thcoor(50),phcoor(50),radbal(50),rgbbal(3,50)
     *,tempc(500),temps(500)
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
      character(len=2)                       :: text
      logical                                :: forcemirrorsym,test

      external cosd,sind,acosd,asind,atan2d
      integer idot,nfreq,ident,indentsign,leninfile,lenfname
      integer nsmoord,nsmoosp,nsmooor,grtype,nrminsmoo
c
c setup the mapgrentropy symbols
c
      call setup_symbols_mapgr()
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
      pi=4.0*atan(1.0)
      nxdim=500
      nfdim=mexpan
c
c initialise log(factorial)
c
      call fctrln(250)
      do iref=1,nsymbol
        write(6,'(a,1x,a)') symbol(iref),symbol_value(iref)(1:len_trim(symbol_value(iref)))
      end do
      iref=1
      infile=symbol_value(iref) !Name of file containing SH coefficients
! convert to full name
      leninfile=len_trim(infile)
      call get_full_filename(infile,leninfile,fname,lenfname)
!      write(6,97)
97    format(1x,'Type no. of expansion coefficients')
      iref=iref+1
      read(symbol_value(iref),*) nexpan
c
c clean out data arrays
c
      do k=1,nexpan
            do i=1,500
                  datac(i,k)=0.0
                  datas(i,k)=0.0
            end do
      end do
!      write(6,96)
96    format(1x,'Are molecules identical (ident=0 if so, else 1)')
c A value of ident .lt. 0 will be used to signify that we want molecule 2 to be
c at the origin
      iref=iref+1
      read(symbol_value(iref),*) ident
!      read(10,*) ident
c save the sign of ident and keep separate from the absolute value
      if(abs(ident).gt.0) then
            identsign=isign(1,ident)
            ident=abs(ident)
      else
            identsign=1
            ident=0
      endif
      kk0=1
! Get the titles from the input file
      call r_int_titles(fname,nexpan,titlex,titled,nfdim,ierr)
      do k=1,nexpan
         temptitle=titled(k)
         lentemptitle=len_trim(temptitle)
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
      call r_int(nr,fname,nexpan,rd,datac,datas,nxdim,nfdim,ierr)
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,20,nstsize,nworddim)
      rminsmoo=0.0
      nsmoord=0
      nsmoosp=0
      nsmooor=0
      if(nwords.ge.4) then
         read(sometext,*) rminsmoo,nsmoord,nsmoosp,nsmooor
      else if(nwords.ge.3) then
         read(sometext,*) nsmoord,nsmoosp,nsmooor
         rminsmoo=0.0
      else
         read(sometext,*) nsmoord
         rminsmoo=0.0
         nsmoosp=nsmoord
         nsmooor=nsmoord
      end if
! Determine the minimum radius value at which to start smoothing.
      nrminsmoo=1
      do while(rd(nrminsmoo).lt.rminsmoo)
         nrminsmoo=nrminsmoo+1
      end do
! Perform the smoothing on the coefficients as required
      do k=1,nexpan
!         write(6,92) l1ex(k),l2ex(k),lex(k),n1ex(k),n2ex(k)
         nsmoo=nsmoord
         if(l2ex(k).eq.0.and.l1ex(k).ne.0) then
            nsmoo=nsmoosp
         else if(l2ex(k).ne.0.and.l1ex(k).ne.0) then
            nsmoo=nsmooor
         end if
         if(nsmoo.ne.0) then
            do i=1,nr
               fact=rd(i)*rd(i)
               tempc(i)=datac(i,k)*fact
               temps(i)=datas(i,k)*fact
            end do
            call smoosharm(nrminsmoo,nr,nsmoo,tempc)
            call smoosharm(nrminsmoo,nr,nsmoo,temps)
            do i=1,nr
               fact=rd(i)*rd(i)
               datac(i,k)=tempc(i)/fact
               datas(i,k)=temps(i)/fact
            end do
         endif
      end do
!      write(6,'(a)') 'Apply mirror symmetry (y or no)'
!      read(10,'(a)') text
      iref=iref+1
      read(symbol_value(iref),'(a)') text
      forcemirrorsym=text.eq.'y'
      iref=iref+1
      read(symbol_value(iref),*) nviews
! Add an extra view - the centres distribution function
      ituse=1
      do i=1,nr
         grw(i,ituse)=1.0+datac(i,1)
      end do 
      write(titledat(ituse),'(8x,a,9x)') 'g_centres(r)'
c calculate g(r,w) at specified theta and phi values
      do it=1,nviews
		 iref=iref+1
         read(symbol_value(iref),*,iostat=itest) grtype
         iref=iref+1
         read(symbol_value(iref),*,iostat=itest) thetal,phil,phim,thetam,chim
         l1only=1
         l2only=1
         lonly=1
         n1only=1
         n2only=1
         m2fix=1
		 if(grtype.lt.3) then
			l2only=0
			n2only=0
			m2fix=0
		 end if
		 if(grtype.lt.2) then
            l1only=0
            lonly=0
		 	n1only=0
	     end if
         if(itest.eq.0) then
            ctl=cosd(thetal)
            phi=phim
            cz=cosd(thetam)
            chi=chim
!            write(6,*) acos(ctl)*180.0/pi,phil,phi,acos(cz)*180.0/pi,chi
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
               test=(l1only.eq.0.and.l1tmp.eq.0).or.l1only.ne.0
               test=test.and.((l2only.eq.0.and.l2tmp.eq.0).or.l2only.ne.0)
               test=test.and.((lonly.eq.0.and.ll.eq.0).or.lonly.ne.0)
               test=test.and.((n1only.eq.0.and.n1tmp.eq.0).or.n1only.ne.0)
               test=test.and.((n2only.eq.0.and.n2tmp.eq.0).or.n2only.ne.0)
               if(test) then
                  phs=1.0
                  phsn=((-1)**(l1tmp+l2tmp+ll+n1tmp+n2tmp))

c Normal algorithm for molecule 1 at the origin

                  if(identsign.ge.0) then
                     if(m2fix.gt.0) then
                        m2f=-l2tmp
                        m2l=l2tmp
                     else
                        m2f=0
                        m2l=0
                     end if
                     call sumoverm2(l1tmp,l2tmp,ll,n1tmp,m2f,m2l,n2tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
c now reverse signs of n1,n2 if either is non-zero
                     if((n2+n1).gt.0) then
                        n1tmp=-n1
                        n2tmp=-n2
                        phs=phs*phsn
                        call sumoverm2(l1tmp,l2tmp,ll,n1tmp,m2f,m2l,n2tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
                     end if
                  else
                     if(m2fix.gt.0) then
                        m2f=-l1tmp
                        m2l=l1tmp
                     else
                        m2f=0
                        m2l=0
                     end if

c Phase factor to reverse parity of centre to centre vector

                     phs=(-1)**ll
                     call sumoverm2(l2tmp,l1tmp,ll,n2tmp,m2f,m2l,n1tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
c now reverse signs of n1,n2 if either is non-zero
                     if((n2+n1).gt.0) then
                        n1tmp=-n1
                        n2tmp=-n2
                        phs=phsn*phs
                        call sumoverm2(l2tmp,l1tmp,ll,n2tmp,m2f,m2l,n1tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
                     end if
                  end if
               end if
c now repeat for interchange of l1,n1,with l2,n2 for identical molecules
c bearing in mind the phase convention on the h's - the case when l1=l2 AND
c n1=n2 doesn't count here. 
c
               if(ident.eq.0.and.(ll1.ne.ll2)) then !l1 will always be .ge. l2 (from the way they are calculated).
                  phs=(-1)**(ll1+ll2)
                  l1tmp=ll2
                  l2tmp=ll1
                  n1tmp=n2
                  n2tmp=n1
                  test=(l1only.eq.0.and.l1tmp.eq.0).or.l1only.ne.0
                  test=test.and.((l2only.eq.0.and.l2tmp.eq.0).or.l2only.ne.0)
                  test=test.and.((lonly.eq.0.and.ll.eq.0).or.lonly.ne.0)
                  test=test.and.((n1only.eq.0.and.n1tmp.eq.0).or.n1only.ne.0)
                  test=test.and.((n2only.eq.0.and.n2tmp.eq.0).or.n2only.ne.0)
                  if(test) then
                     phsn=((-1)**(l1tmp+l2tmp+ll+n1tmp+n2tmp))
                     if(m2fix.gt.0) then
                        m2f=-l2tmp
                        m2l=l2tmp
                     else
                        m2f=0
                        m2l=0
                     endif
                     call sumoverm2(l1tmp,l2tmp,ll,n1tmp,m2f,m2l,n2tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
c now reverse signs of n1,n2 if either are non-zero, not including when both are zero
                     if((n2+n1).gt.0) then
                        n1tmp=-n2
                        n2tmp=-n1
                        phs=phs*phsn
                        call sumoverm2(l1tmp,l2tmp,ll,n1tmp,m2f,m2l,n2tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
                     end if
                  end if
               end if
               zcoeffc(k)=sumc
               if(forcemirrorsym) then
                  zcoeffs(k)=0.0
               else
                  zcoeffs(k)=sums
               end if
!               write(6,*) ll1,ll2,ll,n1,n2,phil,ctl,phi,cz,sumc,sums
            end do
            do i=1,nr
               sumt=0
               do k=1,nexpan
                  term=datac(i,k)*zcoeffc(k)+datas(i,k)*zcoeffs(k)
                  sumt=sumt+term
               end do
               temp(i)=sumt+1.0
            end do
            ituse=ituse+1
            do i=1,nr
                  grw(i,ituse)=temp(i)
                  grwerr(i,ituse)=0.0
            end do
! Write the title for this distribution
            write(titledat(ituse),359) ituse-1
359   format(10x,'view_',i3.3)
         end if
      end do
c
c write out the results
c

c set the axis title

      titler='r[A]'
! Generate name of output file
      idot=indexlastdot(fnameinp)
      if(idot.gt.0) then
         fnameout=fnameinp(1:idot)//'g'
      else
         fnameout=fnameinp(1:len_trim(fnameinp))//'.g'
      end if
      write(6,'(a,1x,a)') fnameinp(1:len_trim(fnameinp)),fnameout(1:len_trim(fnameout))
      call w_int(nr,fnameout,ituse,rd,grw,grwerr,titler,titledat,nxdim,nfdim)
      end               

