      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_splot2d()

      stop
      end      

      subroutine run_splot2d()
c
c invokes the system command to perform 3d plotting
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      character*256 plot_output,fname_out
      integer*4 ierrout
c
c setup the plot2d symbols
c
      call setup_symbols_splot2d()
c
c now setup the input file
c
      nin=11
      nout=6
      call setup_input_file(nin,nout,ierrout)
      write(6,102) prompt(1:ilenprompt),fname_inp(1:ilenfname_inp)
102   format(/a,1x,'Will use splot2d input file: "',a,'"')
c
c read the input file
c
      if(ierrout.eq.0) call read_plot_input_file(nin,nout,ierrout)
c
c if an error was encountered on setup, exit now
c
      if(ierrout.ne.0) return
      call map2dsplot(fname_inp,fname_out) !Run the analysis and output the splot file
c
c now generate the command that will actually run this file
c
      plot_output=system_gnu(1:ilensystem_gnu)
     *//' '
     *//system_quotes(1:ilensystem_quotes)
     *//fname_out(1:len_trim(fname_out))
     *//system_quotes(1:ilensystem_quotes)
     *//' -'
      write(6,*) plot_output
      call system(plot_output,ierr)
      return
      end

      subroutine map2dsplot(infile,plotfile)
      
      use lngamma_mod

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'plot_pars.inc'

      real, parameter :: pi=3.141592653
      integer, parameter :: mexpan=3000,mval=10
      real datac(500,mexpan),datas(500,mexpan),rd(500)
     1,zcoeffc(mexpan),zcoeffs(mexpan)
      integer l1ex(mexpan),l2ex(mexpan),lex(mexpan),n1ex(mexpan),n2ex(mexpan)
      real, dimension(:,:), allocatable :: dens,grw,rgbsurf
      real, dimension(:), allocatable :: temp,colourfrac
      character(len=256) :: fname,outfile,datafile,somestring,addstring
      character(len=*) :: infile,plotfile
      character(len=14) :: titlex
      character(len=29), dimension(mexpan) :: titled
      character(len=29) :: temptitle
      integer npixel,mpixel
      integer i,j,k,ir,irm,irp,nth,nthm,nthp,nph,nphm,nphp,xinterp,yinterp
      real :: addcons
      character(len=10) wherecont
      integer :: l1only,l2only,lonly,n1only,n2only,m2fix

      logical forcemirrorsym,test,surfonly,setgrid
c used to specify words in sometext
      integer irmax,ncolour
      real sumc,sums,sumt,xdim,ydim,rmin

      external cosd,sind,acosd,asind,atan2d

      nxdim=500
      nfdim=mexpan
      nin=10
c
c initialise log(factorial)
c
      call fctrln(250)
      fname=infile
! Name of output file - remove any text after the last dot
      lenfirstdot=index(fname,'.splot2d.txt') ! Make sure this is a splot2d inputfile
      if(lenfirstdot.gt.0) then
         outfile=fname(1:lenfirstdot)//'splot2d'
      else
         lendfirstdot=indexlastdot(fname)
         if(lenfirstdot.gt.0) then
            outfile=fname(1:lenfirstdot)//'splot2d'
         else
            outfile=fname(1:len_trim(fname))//'.splot2d'
         end if
      end if
      lenoutfile=len_trim(outfile)
      plotfile='plot'//outfile(1:indexlastdot(outfile))//'txt'
      iref=0
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      call concatenate(1,nwords,datafile,leninfile,nwords,ncf,ncl,sometext)
      iref=iref+1
      read(symbol_value(iref),*) nexpan
c A value of ident .lt. 0 will be used to signify that we want molecule 2 to be
c at the origin
      iref=iref+1
      read(symbol_value(iref),*) ident

c save the sign of ident and keep separate from the absolute value

      if(abs(ident).gt.0) then
            identsign=isign(1,ident)
            ident=abs(ident)
      else
            identsign=1
            ident=0
      endif
! Get the titles from the input file
      call r_int_titles(datafile,nexpan,titlex,titled,nfdim,ierr)
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
C read in spherical harmonics data
c
      call r_int(nr,datafile,nexpan,rd,datac,datas,nxdim,nfdim,ierr)
      rstep=rd(3)-rd(2)
      iref=iref+1
      read(symbol_value(iref),*) nsmoo
      iref=iref+1
      read(symbol_value(iref),*) xdim
      iref=iref+1
      read(symbol_value(iref),*) ydim
      iref=iref+1
      read(symbol_value(iref),*) mpixel
      iref=iref+1
      read(symbol_value(iref),*) rmin,rmax
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.gt.2) then
         read(sometext(ncf(1):ncl(nwords)),*) l1only,l2only,lonly
      else 
         read(sometext(ncf(1):ncl(nwords)),*) l1only,l2only
         lonly=1
      end if
      iref=iref+1
      read(symbol_value(iref),*) n1only,n2only
      iref=iref+1
      read(symbol_value(iref),*) m2fix
      iref=iref+1
      forcemirrorsym=symbol_value(iref).eq.'y'
      iref=iref+1
      read(symbol_value(iref),*) nvary
      iref=iref+1
      if(nvary.eq.1) then
         read(symbol_value(iref),*) philin,phim,thetam,chim
      else if(nvary.eq.2) then
         read(symbol_value(iref),*) thetal,phim,thetam,chim
      else if(nvary.eq.3) then
         read(symbol_value(iref),*) thetal,philin,thetam,chim
      else if(nvary.eq.4) then
         read(symbol_value(iref),*) thetal,philin,phim,chim
      else 
         read(symbol_value(iref),*) thetal,philin,phim,thetam
      endif
      iref=iref+1
      read(symbol_value(iref),*) ncolour
      nsize=ncolour
      if(allocated(colourfrac)) deallocate(colourfrac,rgbsurf)
      allocate(colourfrac(nsize),rgbsurf(3,nsize))
      do is=1,ncolour
         iref=iref+1
         read(symbol_value(iref),*) colourfrac(is)
         iref=iref+1
         read(symbol_value(iref),*) (rgbsurf(i,is),i=1,3)
         write(6,*) colourfrac(is),(abs(rgbsurf(i,is)),i=1,3)
      end do
      iref=iref+1
      read(symbol_value(iref),*) ncontour
      iref=iref+1
      read(symbol_value(iref),*) wherecont
      iref=iref+1
      read(symbol_value(iref),*) addcons
      iref=iref+1
      surfonly=symbol_value(iref).eq.'y'
c
c plot title
c
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.gt.0) then
         call concatenate(1,nwords,title,ilentitle
     *,nwords,ncf,ncl,sometext)
      else
         title=' '
         ilentitle=1
      endif
c
c set up a default title in case it is undefined
c
      if(title(1:ilentitle).eq.'<undefined>') then
         ilenuptodot=indexlastdot(outfile)-1
         title=directory_name(1:nchardir)
     *//outfile(1:ilenuptodot)
         ilentitle=nchardir+ilenuptodot
      endif
c
c x-axis label
c
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.gt.0) then
         call concatenate(1,nwords,xlabel,ilenxlabel
     *,nwords,ncf,ncl,sometext)
      else
         xlabel=' '
         ilenxlabel=1
      endif
c
c y-axis label
c
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.gt.0) then
         call concatenate(1,nwords,ylabel,ilenylabel
     *,nwords,ncf,ncl,sometext)
      else
         ylabel=' '
         ilenylabel=1
      endif
      iref=iref+1
      setgrid=symbol_value(iref).eq.'y'
      iref=iref+1
      read(symbol_value(iref),*) xinterp,yinterp
      if(xinterp.gt.99) xinterp=99
      if(yinterp.gt.99) yinterp=99
      iref=iref+1
      read(symbol_value(iref),*) phiv
! Set up pixels and voxels. Currently a cubic plot area is assumed.
      nx=mpixel
      ny=mpixel
      xstep=xdim/real(nx)
      ystep=ydim/real(ny)
      nx0=2*nx+1
      ny0=2*ny+1
      nz0=145
      zstep=2.5
      if(allocated(dens)) deallocate(dens,grw,temp)
      allocate(dens(nx0,ny0),grw(nz0,nr),temp(nr))
      npixel=nx0*ny0
c
c initial maximum and minimum values
c
      dhigh=1.0
      dlow=1.0
      thetammax=0.0
      irmax=1 !Reference to distance at which maximum intensity is found
      do j=1,nz0
         z=(j-1)*zstep
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
         if(thetam.gt.180.) then
            thetam=360-thetam
            phi=phi+180.
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
            test=(l1only.eq.0.and.l1tmp.eq.0).or.l1only.ne.0
            test=test.and.((l2only.eq.0.and.l2tmp.eq.0).or.l2only.ne.0)
            test=test.and.((lonly.eq.0.and.ll.eq.0).or.lonly.ne.0)
            test=test.and.((n1only.eq.0.and.n1tmp.eq.0).or.n1only.ne.0)
            test=test.and.((n2only.eq.0.and.n2tmp.eq.0).or.n2only.ne.0)
            if(test) then
               phs=1.0
               phsn=((-1)**(l1tmp+l2tmp+ll+n1tmp+n2tmp))
! Normal algorithm for molecule 1 at the origin
               if(identsign.ge.0) then
                  if(m2fix.gt.0) then
                     m2f=-l2temp
                     m2l=l2tmp
                  else
                     m2f=0
                     m2l=0
                  endif
                  call sumoverm2(l1tmp,l2tmp,ll,n1tmp,m2f,m2l,n2tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
! Now reverse signs of n1,n2 if n1 is non-zero
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
! Phase factor to reverse parity of centre to centre vector
                  phs=(-1)**ll
                  call sumoverm2(l2tmp,l1tmp,ll,n2tmp,m2f,m2l,n1tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
! Now reverse signs of n1,n2 if either is non-zero
                  if(n1.gt.0) then
                     n1tmp=-n1
                     n2tmp=-n2
                     phs=phsn*phs
                     call sumoverm2(l2tmp,l1tmp,ll,n2tmp,m2f,m2l,n1tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
                  endif
               endif
            end if
! For identical molecules repeat interchanging of l1,n1 with l2,n2
! bearing in mind the phase convention on the h's - the case when l1=l2 AND
! n1=n2 doesn't count here. 
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
                  if(n1.gt.0) then
                     n1tmp=-n2
                     n2tmp=-n1
                     phs=phs*phsn
                     call sumoverm2(l1tmp,l2tmp,ll,n1tmp,m2f,m2l,n2tmp,phil,ctl,phi,cz,chi,phs,sumc,sums)
                  endif
               endif
            end if
            zcoeffc(k)=sumc
            if(forcemirrorsym) then
               zcoeffs(k)=0.0
            else
               zcoeffs(k)=sums
            end if
!            write(6,*) ll1,ll2,ll,n1,n2,phil,ctl,phi,cz,sumc,sums
         end do
         do i=1,nr
            sumt=0.0
            do k=1,nexpan
               sumt=sumt+datac(i,k)*zcoeffc(k)+datas(i,k)*zcoeffs(k)
!      write(6,*) i,k,l1ex(k),l2ex(k),lex(k),n1ex(k),n2ex(k),zcoeffc(k),zcoeffs(k),datac(i,k),datas(i,k),sumt
            end do
c test for maximum and minimum in selected r range
            r=rd(i)
            temp(i)=0.0
            if(r.ge.rmin) then
               if(r.le.rmax) then
                  temp(i)=sumt+addcons
               else
                  temp(i)=addcons
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
            grw(j,i)=temp(i)
            if(grw(j,i).lt.dlow) dlow=grw(j,i)
            if(grw(j,i).gt.dhigh) then
               dhigh=grw(j,i)
               irmax=i
               thetammax=z
            endif
         end do
      end do
      write(6,*) 'Max. and min. values in this function: ',dhigh,dlow
      write(6,*) 'Radius and angle at maximum density: ',rd(irmax),thetammax
c
c now interpolate onto a square grid
c
      do j=1,ny0
         y=real(j-ny-1)*ystep
         y2=y*y
         do i=1,nx0
            x=real(i-nx-1)*xstep    
            x2=x*x
            r=sqrt(x2+y2)
c determine pair of rd values
            sum=0.0
            if(r.le.rd(nr)) then
               ir=nint(r/rstep)+1
               ir=max(2,ir) !This is to prevent irm referencing undefined grw at low r
               ir=min(nr-1,ir) !This is to prevent irp referencing undefined grw at high r
               irp=ir+1
               irm=ir-1
               rrat=r/rstep-real(ir-1)
               xv=x*cosd(phiv)+y*sind(phiv)
               yv=y*cosd(phiv)-x*sind(phiv)
               theta=atan2d(xv,yv)
               if(theta.lt.0.) theta=theta+360.
               nth=nint(theta/zstep)+1
               if(nth.gt.nz0) nth=nth-nz0
               if(nth.lt.1) nth=nth+nz0
               nthp=nth+1
               if(nthp.gt.nz0) nthp=nthp-nz0
               nthm=nth-1
               if(nthm.lt.1) nthm=nthm+nz0
               thrat=theta/zstep-real(nth-1)
! Perform quadratic interpolation
               sum=quadinterp2d(grw,rrat,thrat,ir,irp,irm,nth,nthp,nthm,nz0,nr)
!
! Form mean over phi values
!              sum=suminterp(grw,rrat,thrat,ir,irp,irm,nth,nthp,nthm,np0,nphp,nphm,50,50,500)
            else
               sum=addcons
            end if
            dens(i,j) = max(sum,0.0)
         end do
      end do
      nout=nin
      !Open the splot output file. 
      open(nout,file=outfile,status='unknown')
      do iy=1,ny0
         y=real(iy-ny-1)*ystep
         do ix=1,nx0
            x=real(ix-ny-1)*xstep
            write(nout,*) x,y,dens(ix,iy)
         end do
         write(nout,*) !Write a blank line for splot format
      end do
      close(nout)
! Write plot file
      open(nout,file=plotfile,status='unknown')
c
c set up a default title in case it is undefined
c
      write(nout,111) title(1:ilentitle)
111   format('reset',/,'unset key'/'set title ',1h',a,1h')
      write(nout,112) xlabel(1:ilenxlabel)
112   format('set xlabel ',1h',a,1h')
      write(nout,122)
122   format('set xtics out')
      write(nout,113) ylabel(1:ilenylabel)
113   format('set ylabel ',1h',a,1h')
      write(nout,123)
123   format('set ytics out')
      write(nout,102)
102   format('set linetype 1 lc rgb "black"'/'set linetype cycle 1')
      if(surfonly) write(nout,103)
103   format('set pm3d map')
      write(nout,'("set pm3d interpolate ",i2,",",i2)') xinterp,yinterp
      write(nout,104)
104   format('set palette model RGB')
!Prepare the string to define the colour palette
      icolour=1
      write(somestring,105) abs(colourfrac(icolour)),(abs(rgbsurf(i,icolour)),i=1,3)
105   format('set palette defined (',f5.2,3(1x,f5.2))
!      write(6,*) colourfrac(icolour),(abs(rgbsurf(i,icolour)),i=1,3)
!      write(6,*) somestring(1:len_trim(somestring))
      do while(icolour.lt.ncolour)
         icolour=icolour+1
         write(addstring,106) abs(colourfrac(icolour)),(abs(rgbsurf(i,icolour)),i=1,3)
106   format(',',4(1x,f5.2))
         somestring=somestring(1:len_trim(somestring))//addstring(1:len_trim(addstring))         
      end do
      somestring=somestring(1:len_trim(somestring))//')'
!      write(6,'(a)') somestring(1:len_trim(somestring))
      write(nout,107) somestring(1:len_trim(somestring))
107   format(a)
      if(ncontour.gt.0) then
         write(nout,108) wherecont(1:len_trim(wherecont))
         write(nout,118) ncontour
108   format('set contour ',a)
118   format('set cntrparam levels ',i5)
      else
         write(nout,109)
109   format('unset contour')
      end if
      if(setgrid) then
         write(nout,119)
119   format('set grid')
      else
         write(nout,129)
129   format('unset grid')
      end if
      write(nout,110) outfile(1:len_trim(outfile))
110   format('splot "',a,'" w pm3d')
      close(nout)
      return
      end               







