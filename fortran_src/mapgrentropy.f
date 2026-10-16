      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_mapgrentropy()

      stop
      end      

      subroutine run_mapgrentropy()

      use lngamma_mod
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
c
c This version of MAP3D puts out the g(r) in specifed directions
c
      parameter(mexpan=10000,mval=10,mpixel=80,mpixel2=(2*mpixel+1)
     1,mpixel3=mpixel2**3)
      real datac(500,mexpan),datas(500,mexpan),rd(500),grw(500,5),grwerr(500,5)
     1,zcoeffc(mexpan),zcoeffs(mexpan),zcoeffspc(mexpan),zcoeffsps(mexpan)
     *,sumc,sums,spsumc,spsums
      dimension lvl(3,mval),nl(3),ml(3)
     1,th(50),ph(50),xl(2),yl(2),temprd(500),tempsp(500),tempto(500)
     2,out(0:mpixel3),xx(201),yy(201),l1ex(mexpan)
     3,l2ex(mexpan),lex(mexpan),n1ex(mexpan),n2ex(mexpan)
     *,rdcoor(50),thcoor(50),phcoor(50),radbal(50),rgbbal(3,50)
     *,tempc(500),temps(500)
      character*256 infile,fname,fnameout
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
      real pi,term,temp,tsumrd,tsumsp,tsumto
      real molndens,molmass,moltemp,mollambda,molrot,molvib
      real kineticterm,rotterm,vibterm,transterm
      integer idot,nfreq,ident,indentsign,leninfile,lenfname
      integer nsmoord,nsmoosp,nsmooor,nrminsmoo

      real(kind=8) entropyrd,entropysp,entropyto,sumrd,sumsp,sumto,grw8(500,3),grwerr8(500,3)
      real minrd,minsp,minto,minrdr,minspr,mintor
      real maxrd,maxsp,maxto,maxrdr,maxspr,maxtor
      real thetalspmin,philspmin,phimspmin,thetamspmin,chimspmin
      real thetaltomin,philtomin,phimtomin,thetamtomin,chimtomin
      real spnegative(500),tonegative(500),wtnegative !Store for negative values of g(r) as a function of r. 
! This will be used to correct values at same distance which are above 1.0 to compensate.

c
c setup the mapgrentropy symbols
c
      call setup_symbols_mapgrentropy()
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
      open(10,file='ran1.dat',status='unknown')
      call initran1(10)
      close(10)
      idot=indexlastdot(fnameinp)
      if(idot.gt.0) then
         fnameout=fnameinp(1:idot)//'out'
      else
         fnameout=fnameinp(1:len_trim(fnameinp))//'.out'
      end if
      write(6,'(a,1x,a)') fnameinp(1:len_trim(fnameinp)),fnameout(1:len_trim(fnameout))
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
!      read(10,*) nexpan
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
!      write(6,105)
!105   format(1x,'Specify number of trials to calculate entropy for, and frequency to write')
      iref=iref+1
      read(symbol_value(iref),*) nvary
      iref=iref+1
      read(symbol_value(iref),*) nfreq
      iref=iref+1
      read(symbol_value(iref),*) molndens
      iref=iref+1
      read(symbol_value(iref),*) molmass
      iref=iref+1
      read(symbol_value(iref),*) moltemp
      iref=iref+1
      read(symbol_value(iref),*) molrot
      iref=iref+1
      read(symbol_value(iref),*) molvib
      gasconstant=8.314472
      !De Broglie wavelength [A]
      mollambda=66.2606876/sqrt(2*pi*molmass*1.660539*1.3806503*moltemp)
      transterm=gasconstant*(2.5-log(molndens*(mollambda*(mollambda*mollambda))))
      rotterm=gasconstant*(1.5+log(sqrt(pi*(moltemp**3)/molrot)))
      vibterm=molvib
      kineticterm=transterm+rotterm+vibterm
!      read(10,*) nvary,nfreq
!121   format(a50)
c calculate g(r,w) at specified angular values
      entropyrd=kineticterm !Radial entropy
      entropysp=kineticterm !Spatial entropy
      entropyto=kineticterm !total entropy
      do it=1,3
         do i=1,nr
            grw8(i,it)=kineticterm
            grwerr8(i,it)=kineticterm
         end do
      end do
! zero negative value accumulators
      do i = 1,nr
         spnegative(i)=0.0
         tonegative(i)=0.0
      end do
      fact=2.0*pi*gasconstant*molndens*(rd(3)-rd(2))
      open(10,file=fnameout,status='unknown')
      write(10,'("Temperature [K]: ",e12.5)') moltemp
      write(10,'("Number density [mols/A**3]: ",e12.5)') molndens
      write(10,'("De Broglie wavelength: ",e12.5)') mollambda
      write(10,'("Perfect gas entropy, trans. term: ",e12.5)') transterm
      write(10,'("Perfect gas entropy, rot. term: ",e12.5)') rotterm
      write(10,'("Perfect gas entropy,vib. term: ",e12.5)') vibterm
      write(10,'("Perfect gas entropy, total: ",e12.5)') kineticterm
      minrd=0.0
      minsp=0.0
      minto=0.0
      maxrd=0.0
      maxsp=0.0
      maxto=0.0
      l1only=1
      l2only=1
      lonly=1
      n1only=1
      n2only=1
      m2fix=1
      it=0
      do while(it.lt.nvary)
            it=it+1
            ctl=2.0*ran1()-1.0
            phil=ran1()*360.0
            phi=ran1()*360.0
            cz=2.0*ran1()-1.0
            chi=ran1()*360.0
!            write(6,*) acos(ctl)*180.0/pi,phil,phi,acos(cz)*180.0/pi,chi
            do k=1,nexpan
               sumc=0.0
               sums=0.0
               spsumc=0.0
               spsums=0.0
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
                     if (l2tmp.eq.0) then !Store the spatial terms separately
                        spsumc=sumc
                        spsums=sums
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
                     if (l1tmp.eq.0) then !Store the spatial terms separately since reversing l1 and l2 will not lead to a spatial term.
                        spsumc=sumc
                        spsums=sums
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
c now reverse signs of n1,n2 if either is non-zero
                     if((n2+n1).gt.0) then
                        n1tmp=-n2
                        n2tmp=-n1
                        phs=phs*phsn
                        call sumoverm2(l1tmp,l2tmp,ll,n1tmp,m2f,m2l,n2tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
                     end if
                  end if
               end if
               zcoeffc(k)=sumc
               zcoeffspc(k)=spsumc
               if(forcemirrorsym) then
                  zcoeffs(k)=0.0
                  zcoeffsps(k)=0.0
               else
                  zcoeffs(k)=sums
                  zcoeffsps(k)=spsums
               end if
!               write(6,*) ll1,ll2,ll,n1,n2,phil,ctl,phi,cz,sumc,sums
            end do
            do i=1,nr
                  tsumrd=0.0
                  tsumsp=0.0
                  tsumto=0.0
                  do k=1,nexpan
                        term=datac(i,k)*zcoeffc(k)+datas(i,k)*zcoeffs(k)
                        tsumto=tsumto+term
                        if(l2ex(k).eq.0) then
                           term=datac(i,k)*zcoeffspc(k)+datas(i,k)*zcoeffsps(k)
                           tsumsp=tsumsp+term
                           if(l1ex(k).eq.0) then
                              tsumrd=tsumrd+term
                           end if
                        end if
                  end do
                  temprd(i)=tsumrd
                  tempsp(i)=tsumsp
                  tempto(i)=tsumto
            end do
            sumrd=kineticterm
            sumsp=kineticterm
            sumto=kineticterm
            do i=1,nr
                  factr=fact*rd(i)*rd(i)
                  temp=temprd(i)
                  term=1.0+temp
                  if(term.gt.maxrd) then
                     maxrd=term
                     maxrdr=rd(i)
                  end if
                  if(term.lt.minrd) then
                     minrd=term
                     minrdr=rd(i)
                  end if
                  if(term.gt.0.0) then
                     contribution=factr*(temp-term*log(term))
                  else
                     contribution=-factr
                  end if
                  sumrd=sumrd+contribution
                  grw8(i,1)=grw8(i,1)+sumrd
                  temp=tempsp(i)
                  term=1.0+temp
                  if(term.gt.maxsp) then
                     maxsp=term
                     maxspr=rd(i)
                     thetalspmax=acos(ctl)*180.0/pi
                     philspmax=phil
                     phimspmax=phi
                     thetamspmax=acos(cz)*180.0/pi
                     chimspmax=chi
                  end if
                  if(term.lt.minsp) then
                     minsp=term
                     minspr=rd(i)
                     thetalspmin=acos(ctl)*180.0/pi
                     philspmin=phil
                     phimspmin=phi
                     thetamspmin=acos(cz)*180.0/pi
                     chimspmin=chi
                  end if
!                  if(term.lt.0.0) then
!                     spnegative(i)=spnegative(i)-term
!                  end if
!Make correction for negative regions of g(r)
                  if(temp.gt.0.0.and.spnegative(i).gt.0.0) then
                     wtnegative=temp*spnegative(i)/term
                     temp=temp+wtnegative
                     spnegative(i)=spnegative(i)-wtnegative
                     term=1.0+temp
                  end if
                  if(term.gt.0.0) then
                     contribution=factr*(temp-term*log(term))
                  else
                     contribution=-factr
                  end if
                  sumsp=sumsp+contribution
                  grw8(i,2)=grw8(i,2)+sumsp
                  temp=tempto(i)
                  term=1.0+temp
                  if(term.gt.maxto) then
                     maxto=term
                     maxtor=rd(i)
                     thetaltomax=acos(ctl)*180.0/pi
                     philtomax=phil
                     phimtomax=phi
                     thetamtomax=acos(cz)*180.0/pi
                     chimtomax=chi
                  end if
                  if(term.lt.minto) then
                     minto=term
                     mintor=rd(i)
                     thetaltomin=acos(ctl)*180.0/pi
                     philtomin=phil
                     phimtomin=phi
                     thetamtomin=acos(cz)*180.0/pi
                     chimtomin=chi
                  end if
!                  if(term.lt.0.0) then
!                     tonegative(i)=tonegative(i)-term
!                  end if
!Make correction for negative regions of g(r)
                  if(temp.gt.0.0.and.tonegative(i).gt.0.0) then
                     wtnegative=temp*tonegative(i)/term
                     temp=temp+wtnegative
                     tonegative(i)=tonegative(i)-wtnegative
                     term=1.0+temp
                  end if
                  if(term.gt.0.0) then
                     contribution=factr*(temp-term*log(term))
                  else
                     contribution=-factr
                  end if
                  sumto=sumto+contribution
                  grw8(i,3)=grw8(i,3)+sumto
            end do
            entropyrd=entropyrd+sumrd
            entropysp=entropysp+sumsp
            entropyto=entropyto+sumto
            if(it.eq.nfreq*int(it/nfreq).or.it.eq.nvary) then
               write(6,'(5(e13.6,1x))') ctl,phil,phi,cz,chi
               write(6,'(a,i10,a,3(1x,e13.6))') 'mapgrentropy> Entropy after ',it,' trials:',entropyrd/real(it),entropysp/real(it),entropyto/real(it)
               write(10,'(i10,6(1x,e13.6))') it,sumrd,sumsp,sumto,entropyrd/real(it),entropysp/real(it),entropyto/real(it)
            end if
      end do
      write(10,*)
      write(10,'("Minimum values, radial, spatial total: ",3(1x,e12.5))') minrd,minsp,minto
      write(10,'("Radius values of minima, radial, spatial total: ",3(1x,e12.5))') minrdr,minspr,mintor
      write(10,*)
      write(10,'("thetal, phil, phim, thetam, chim of spatial minimum: ",5(1x,e12.5))') thetalspmin,philspmin,phimspmin,thetamspmin,chimspmin
      write(10,'("thetal, phil, phim, thetam, chim of total minimum: ",5(1x,e12.5))') thetaltomin,philtomin,phimtomin,thetamtomin,chimtomin
      write(10,*)
      write(10,'("Maximum values, radial, spatial total: ",3(1x,e12.5))') maxrd,maxsp,maxto
      write(10,'("Radius values of maxima, radial, spatial total: ",3(1x,e12.5))') maxrdr,maxspr,maxtor
      write(10,*)
      write(10,'("thetal, phil, phim, thetam, chim of spatial maximum: ",5(1x,e12.5))') thetalspmax,philspmax,phimspmax,thetamspmax,chimspmax
      write(10,'("thetal, phil, phim, thetam, chim of total maximum: ",5(1x,e12.5))') thetaltomax,philtomax,phimtomax,thetamtomax,chimtomax
      close(10)
c
c write out the results
c
      open(10,file='negativegr.txt',status='unknown')
      do i=1,nr
         write(10,'(3(1x,e12.5))') rd(i),spnegative(i),tonegative(i)
      end do
      close(10)
      do i=1,nr
         do it=1,3
            grw(i,it)=grw8(i,it)/real(nvary)
            grwerr(i,it)=grwerr8(i,it)/real(nvary)
         end do
         grw(i,4)=grw(i,2)-grw(i,1)
         grwerr(i,4)=0.0
         grw(i,5)=grw(i,3)-grw(i,2)
         grwerr(i,5)=0.0
      end do
      write(6,*) grw(nr,1),entropyrd/real(nvary)

c set the axis title

      titler='r[A]'
      write(titledat(1),'(10x,a,11x)') 'S_radial'
      write(titledat(2),'(10x,a,10x)') 'S_spatial'
      write(titledat(3),'(11x,a,12x)') 'S_total'
      write(titledat(4),'(7x,a,7x)') 'S_extra_spatial'
      write(titledat(5),'(7x,a,8x)') 'S_extra_orient'
      ilen=indexlastdot(fnameinp)
      ext='s'
      fname=fnameinp(1:ilen)//ext
      call w_int(nr,fname,5,rd,grw,grwerr,titler,titledat,nxdim,5)
      open(10,file='ran1.dat',status='unknown')
      call saveran1(10)
      close(10)
      end               

