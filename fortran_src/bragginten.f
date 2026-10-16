      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call bragginten()

      stop
      end      

      subroutine bragginten()
c
c Program bragginten
c
c Calculates the Bragg intensities from a specified .ato file.
c Intensities are separated into pairs of atom types.
c
c Currently set up for 70 components (number of atom types)
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      character(len=256) :: fname,outputfile,formatout
      integer :: multa,multb,multc !Number of unit cells along each direction a, b, c.
      integer, allocatable, dimension(:) :: h,k,l !hkl values for each Bragg reflection
      integer :: nbragg,nhkl,ncharoutput,it,ib,im,ia,hi,ki,li,hmax,kmax,lmax,nfirst,nlast,nhklmax
      real, dimension(3) :: recipra,reciprb,reciprc,rbragg !Components of the reciprocal lattice vector
      real :: pi,vol,volrecipr,xuse,yuse,zuse,huse,kuse,luse,arg,crystarg !Dummy variables
      real, allocatable, dimension(:,:) :: ramp,iamp  !Real and imaginary parts of scattering amplitude for each atom type and Bragg reflection
      real, allocatable, dimension(:) :: qbragg,qhkl,intenhkl,intenbragg       !Q values of Bragg peaks
      real :: sqrmoda,sqrmodb,sqrmodc,ramptot,iamptot,qbraggtemp,intentemp !Square modulus of each reciprocal lattice vector
      real :: qstep,qmax,qwidth,gfactor,gfactor0,sum,diff,add,lfactor,scatlen(mcomp),qhkltest,scatlevel,ampra,amprb,amprc
      logical :: found

c Get parameters from input file
c
c setup the xwts symbols
c
      call setup_symbols_bragginten()
c
c now setup the input file
c
      nin=11
      nout=6
      call setup_input_file(nin,nout,ierrout)
c
c if an error was encountered on setup, exit now
c
      if(ierrout.ne.0) return
c
c read the input file
c
      call read_input_file(nin,nout,ierrout)
c
c if an error was encountered on read, exit now
c
      if(ierrout.ne.0) return
c
c normalise atomic coordinates to box size
c
      call normalise_coordinates()
      write(6,100)
100   format(/'bragginten> Program to calculate Bragg intensities'
     1,' for the specified .ato file.'/)
99    format(a80)
c
c Set up the base name for the output file. This will be the same as the input file.

      iref=1
      outputfile=fname_run(1:ilenfname_run)
      ncharoutput=ilenfname_run
c
c Read underlying multiplicities along the a, b and c axes - this will have been defined when setting up the .ato file
c
      iref=iref+1
      read(symbol_value(iref),*) multa,multb,multc

c Get the maximum number of h, k, l values

      iref=iref+1
      read(symbol_value(iref),*) hmax,kmax,lmax

c Get q step and q max

      iref=iref+1
      read(symbol_value(iref),*) qstep,qmax

c Get width of broadening function

      iref=iref+1
      read(symbol_value(iref),*) qwidth
      
      iref=iref+1
! Get the scattering lengths

      it=0
      do while(it.lt.nfiles)
         it=it+1
         iref=iref+2
         read(symbol_value(iref),*) scatlen(it)
      end do
      pi=4.0*atan(1.0)
      twopi=2.0*pi
      !Set up the assymtotic scattering level
      scatlevel=0.0
      do it=1,ntype
         if(qtypea(it).eq.0) then
            scatlevel=scatlevel+cfrac(it)*scatlen(it)*scatlen(it)
         end if
      end do

      !Prepare the output format statement
      write(formatout,102) ntype
102   format('(1x,"#",3(1x,i3),1x,f10.5,1x,e12.5,',i2.2,'(1x,2(1x,e12.5)))')
      write(6,*) formatout(1:index(formatout,' '))
      
      !Generate the reciprocal lattice vectors based on supplied information
      call vecpro(blatc,clatc,recipra)
      call vecpro(clatc,alatc,reciprb)
      call vecpro(alatc,blatc,reciprc)
      vol=scapro(alatc,recipra)
      call divvec(recipra,vol)
      call divvec(reciprb,vol)
      call divvec(reciprc,vol)
      !Correct for number of unit cells along each direction
      call mulvec(recipra,real(multa)*twopi)
      call mulvec(reciprb,real(multb)*twopi)
      call mulvec(reciprc,real(multc)*twopi)
      sqrmoda=scapro(recipra,recipra)
      sqrmodb=scapro(reciprb,reciprb)
      sqrmodc=scapro(reciprc,reciprc)
      !Volume of unit cell
      vol=vol/real(multa*multb*multc)
      !Get the volume of the reciprocal lattice
      call vecpro(reciprb,reciprc,rbragg)
      volrecipr=scapro(recipra,rbragg)
      !Get the amplitude of each reciprocal lattice vector
      ampra=sqrt(sqrmoda)
      amprb=sqrt(sqrmodb)
      amprc=sqrt(sqrmodc)
      !Get maximum number of steps in each reciprocal lattice direction
      hmax=nint((qmax+10.0*qwidth)/ampra)+1
      kmax=nint((qmax+10.0*qwidth)/amprb)+1
      lmax=nint((qmax+10.0*qwidth)/amprc)+1
      !Maximum number of hkl values
      nhklmax=(2*hmax+1)*(2*kmax+1)*(2*lmax+1)
      !Allocate arrays
      if(allocated(ramp)) deallocate(ramp)
      if(allocated(iamp)) deallocate(iamp)
      if(allocated(qhkl)) deallocate(qhkl)
      if(allocated(intenhkl)) deallocate(intenhkl)
      allocate(ramp(ntype,nhklmax),iamp(ntype,nhklmax),qhkl(nhklmax),intenhkl(nhklmax))

      !Set these amplitudes to zero
      do ib=1,nhklmax
         do it=1,ntype
            ramp(it,ib)=0.0
            iamp(it,ib)=0.0
         end do
      end do
        
c Create a filename for the output

      fname=directory_name(1:nchardir)//outputfile(1:ncharoutput)//'.txt'
      open(10,file=fname,status='unknown')
      
      write(10,101) '# Average scattering level: ',scatlevel
      write(10,101) '# Volume of unit cell: ',vol
      write(10,101) '# Volume of reciprocal lattice: ',volrecipr
      write(10,101) '# Lattice vectors a: ',alatc(1),alatc(2),alatc(3),sqrt(scapro(alatc,alatc))
      write(10,101) '# Lattice vectors b: ',blatc(1),blatc(2),blatc(3),sqrt(scapro(blatc,blatc))
      write(10,101) '# Lattice vectors c: ',clatc(1),clatc(2),clatc(3),sqrt(scapro(clatc,clatc))
      write(10,101) '# Reciprocal lattice vectors a: ',recipra(1),recipra(2),recipra(3),ampra
      write(10,101) '# Reciprocal lattice vectors b: ',reciprb(1),reciprb(2),reciprb(3),amprb
      write(10,101) '# Reciprocal lattice vectors c: ',reciprc(1),reciprc(2),reciprc(3),amprc
101   format(1x,a,4(1x,e12.5))
!Now run over the allowed h,k,l values and calculate scattering amplitudes and intensities
      ib=0
      hi=-hmax-1
      do while (hi.lt.hmax)
         hi=hi+1
         huse=real(hi*multa)
         ki=-kmax-1
         do while (ki.lt.kmax)
            ki=ki+1
            kuse=real(ki*multb)
            li=-lmax-1
            do while (li.lt.lmax)
               li=li+1
               luse=real(li*multc)
               do it=1,3
                  rbragg(it)=real(hi)*recipra(it)+real(ki)*reciprb(it)+real(li)*reciprc(it)
               end do
!         write(10,*) (rbragg(it),it=1,3)
               qhkltest=sqrt(scapro(rbragg,rbragg))
               if(qhkltest.le.(qmax+5.0*qwidth)) then
                  ib=ib+1
            write(6,*) ib,hi,ki,li
                  qhkl(ib)=qhkltest
                  do im=1,nmol
                     do ia=naf(im),nal(im)
                        it=itype(ia)
                        xuse=xa(ia)+xadd(im)
                        yuse=ya(ia)+yadd(im)
                        zuse=za(ia)+zadd(im)
                        arg=crystarg(xuse,yuse,zuse,huse,kuse,luse)
!               write(6,*) huse,kuse,luse,xuse,yuse,zuse,arg
                        ramp(it,ib)=ramp(it,ib)+cos(arg)
                        iamp(it,ib)=iamp(it,ib)-sin(arg)
                     end do
                  end do
! Sum the real and imaginary parts
                  ramptot=0.0
                  iamptot=0.0
                  do it=1,ntype
                     if(qtypea(it).eq.0) then
                        ramptot=ramptot+scatlen(it)*ramp(it,ib)
                        iamptot=iamptot+scatlen(it)*iamp(it,ib)
                     end if
                  end do
! Total Bragg intensity assuming equal scattering lengths
                  intenhkl(ib)=(ramptot*ramptot+iamptot*iamptot)/(real(natomr))
                  write(10,formatout) hi,ki,li,qhkl(ib),intenhkl(ib),(ramp(it,ib),iamp(it,ib),it=1,ntype)
               end if
            end do
         end do
      end do
      nhkl=ib

      !Set up the Bragg peaks, combining hkl which have the same Q value
      if(allocated(qbragg)) deallocate(qbragg)
      if(allocated(intenbragg)) deallocate(intenbragg)
      nbragg=nint(qmax/qstep)+1
      allocate(qbragg(nbragg),intenbragg(nbragg))
      !Initialise structure factor array
      do ib=1,nbragg
         qbragg(ib)=qstep*real(ib-1)
         intenbragg(ib)=0.0
      end do
      !Factor outside Gaussian
      gfactor=volrecipr/(2.0*(twopi**1.5)*qwidth)
      gfactor0=volrecipr/((twopi**1.5)*(qwidth**3))
      lfactor=volrecipr/((pi**2)*(qwidth**3))
      !Step through hkl values and add broadened intensities to structure factor
      do ib=1,nhkl
         nfirst=max(1,nint((qhkl(ib)-25.0*qwidth)/qstep)+1)
         nlast=min(nbragg,nint((qhkl(ib)+25.0*qwidth)/qstep)+1)
         if(qhkl(ib).eq.0.0) then
            nfirst=1
            nlast=0
         end if
         it=nfirst-1
         do while(it.lt.nlast)
            it=it+1
            sum=0.0
!Gaussian peak shape
!            if(qhkl(ib).gt.0.0.and.qbragg(it).gt.0.0) then
!               diff=(qhkl(ib)-qbragg(it))/qwidth
!               add=(qhkl(ib)+qbragg(it))/qwidth
!               if(abs(diff).lt.8.0) sum=exp(-(0.5*diff*diff))
!               if(abs(add).lt.8.0) sum=sum-exp(-(0.5*add*add))
!               sum=sum*gfactor/(qhkl(ib)*qbragg(it))
!            else
!               add=qhkl(ib)/qwidth
!               diff=qbragg(it)/qwidth
!               sum=0.5*(diff*diff+add*add)
!               if(sum.lt.32.0) then
!                  sum=gfactor0*exp(-sum)
!               else
!                  sum=0.0
!               end if
!            end if
!Lorentzian peak shape
            diff=(qhkl(ib)-qbragg(it))/qwidth
            add=(qhkl(ib)+qbragg(it))/qwidth
            sum=lfactor/((1.0+diff*diff)*(1.0+add*add))
            intenbragg(it)=intenbragg(it)+intenhkl(ib)*sum
         end do
      end do

      do ib=1,nbragg
!Subtract the single atom scattering level
         intenbragg(ib)=intenbragg(ib)-scatlevel
         formatout='(1x,f10.5,1x,e12.5)'
         write(10,formatout) qbragg(ib),intenbragg(ib)
      end do

      close(10)
c
c Rewrite the input file in case anything was changed
c
c      call save_files_others(14)

      return
      end


