      module timer

      implicit none
      
      integer(kind=8) :: sumtime1,sumtime2,time1,time2

      contains

      subroutine initialise_timer()

         sumtime1=0
         sumtime2=0

         return

      end subroutine initialise_timer

      end module timer

      module bragginten_mod

      integer, dimension(:),allocatable                  :: hihkl,kihkl,lihkl
      real, dimension(:,:), allocatable                  :: ramp,iamp,intenhkl,intenqunique    !Real and imaginary parts of scattering amplitude for each atom type and Bragg reflection
      real, dimension(:), allocatable                    :: qhkl,qunique      !Q values of Bragg peaks
      integer                                            :: nhkl,ntotalqunique
      integer, dimension(:), allocatable                 :: nqunique,iqvalues
      integer, dimension(:,:), allocatable :: ibvalues

      end module bragginten_mod

      module reciprocallattice_mod

      real, dimension(3) :: recipra,reciprb,reciprc,rbragg !Components of the reciprocal lattice vector
      real :: volrecipr,ampra,amprb,amprc

      end module reciprocallattice_mod

      module diffuse_mod

      integer                             :: multdiffa,multdiffb,multdiffc,ndiffuse,nqdiff,ldiffmax,ldiffstep,mdiffstep !No. of unit cells along each of h,k,l
      integer, dimension(:), allocatable  :: l1ex,n1ex,ndiffcoeff ! Sets the pattern along each axis
      real, dimension(:,:,:), allocatable :: diffcoeffc,diffcoeffs
      real, dimension(:,:), allocatable   :: diffcoefftc,diffcoeffts,diffcoeffoldc,diffcoeffolds
      real, dimension(:), allocatable     :: qdiff
      real                                :: qdiffstep                    
      logical                             :: diffusedefined
      character(len=3)                    :: diffusename

      contains

      subroutine initialise_diffuse_values(nfile)

      use lngamma_mod
      use bragginten_mod
      use reciprocallattice_mod

      integer mdiffmax,maxlngamma

! First count the number of terms and find the maximum ln(gamma) function required
      ndiffuse=0
      maxlngamma=0
      do l=0,ldiffmax,ldiffstep
         mdiffmax=l
         if(mdiffstep.eq.0) then
            mdiffmax=0
         end if
         do m=0,mdiffmax
            ndiffuse=ndiffuse+1
            maxlngamma=max(maxlngamma,maxdlm(l,m,0),maxdlm(l,-m,0)) !Get the maximum value of the ln(gamma) function to be called
         end do
      end do
      call fctrln(maxlngamma)
! Now allocate arrays and save the values
      if(allocated(l1ex)) deallocate(l1ex,n1ex)
      allocate(l1ex(ndiffuse),n1ex(ndiffuse))
      ndiffuse=0
      do l=0,ldiffmax,ldiffstep
         mdiffmax=l
         if(mdiffstep.eq.0) then
            mdiffmax=0
         end if
         do m=0,mdiffmax,mdiffstep
            ndiffuse=ndiffuse+1
            l1ex(ndiffuse)=l
            n1ex(ndiffuse)=m
         end do
      end do
! Get the maximum Q of this scattering
      qmaxdiff=0.0
      do ib=1,nhkl
         if(qhkl(ib).gt.qmaxdiff) qmaxdiff=qhkl(ib)
      end do
      qdiffstep=max(ampra,amprb,amprc)
      nqdiff=nint(qmaxdiff/qdiffstep)+1
! Set up the arrays to store the radial coefficients
      if(allocated(diffcoeffc)) deallocate(qdiff,diffcoeffc,diffcoeffs,ndiffcoeff)
      allocate(qdiff(nqdiff),diffcoeffc(ndiffuse,nqdiff,nfile),diffcoeffs(ndiffuse,nqdiff,nfile),ndiffcoeff(nqdiff))
      do iff=1,nfile
         do iq=1,nqdiff
            qdiff(iq)=real(iq-1)*qdiffstep
            do id=1,ndiffuse
               ndiffcoeff(iq)=0
               diffcoeffc(id,iq,iff)=0.0
               diffcoeffs(id,iq,iff)=0.0
            end do
         end do
      end do
      return
      end subroutine initialise_diffuse_values
            
      end module diffuse_mod

      subroutine ftype(ia)

      use reallocation_routines
      use molecom_mod
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      character*1 qtypet
c
c is this atom of "Q"-type?
c
      qtypet=alabel(ia)(1:1)
      if(qtypet.eq.'Q'.or.qtypet.eq.'q') then
            qtype(ia)=1
      else
            qtype(ia)=0
      endif
      ifound=0
      it=0
      do while(it.lt.ntype.and.ifound.eq.0)
         it=it+1
         if(alabel(ia).eq.type(it)) ifound=it
      end do
      if(ifound.gt.0) then
         itype(ia)=ifound
      else if(ntype.lt.near) then
         ntype=ntype+1
         type(ntype)=alabel(ia)
         qtypea(ntype)=qtype(ia)
         itype(ia)=ntype
! Allocate space for the coarse-graining parameters that may be defined
         if(ntype.gt.size(rhosphere,1)) then
            call reallocate1d_r(rhosphere,ntype)
            call reallocate1d_r(radminsphere,ntype)
            call reallocate1d_r(radmaxsphere,ntype)
         end if
      else
         write(6,100) near
100   format(1x,'Tried to generate more than ',i3,' types of atom')
         stop
      endif
!      write(6,*) ntype,size(rhosphere)
      return
      end

      subroutine ctype(a1,ats,iexch,eread,wread,amasst,qt,rqt,rhosph,radminsph,radmaxsph)

      use molecom_mod
      include 'epsr_dimension.inc'
      include 'molecom.inc'

      character*3 a1,ats
      integer get_atomicnumber,n1,it
      real, intent(in),optional :: rhosph,radminsph,radmaxsph
      n1=0
      it=0
      do while(it.lt.ntype.and.n1.eq.0)
         it=it+1
         if(a1.eq.type(it)) n1=it
      end do
      if(n1.gt.0) then
            iatomicnotype(n1)=get_atomicnumber(ats)
            atsymbol(n1)=ats
            iexchange(n1)=iexch
            epscor(n1)=eread
            sigcor(n1)=wread
            amassa(n1)=amasst
            charge(n1)=qt
            qradius(n1)=rqt
            if(present(rhosph)) then
               rhosphere(n1)=rhosph
            else
               rhosphere(n1)=0.0
               radminsphere(n1)=0.0
               radmaxsphere(n1)=0.5*sigcor(n1)
            end if
            if(present(radminsph)) then
               radminsphere(n1)=radminsph
            else
               radminsphere(n1)=0.0
            end if
            if(present(radmaxsph)) then
               radmaxsphere(n1)=radmaxsph
            else
               radmaxsphere(n1)=0.5*sigcor(n1)
            end if
            if(radmaxsphere(n1).gt.0.5*sigcor(n1)) radmaxsphere(n1)=0.5*sigcor(n1)
      else
            write(6,100) a1
100   format(1x,'Mismatched types between atomic coordinates and '
     1,'potential values: ',a)
      endif
      return
      end         

      function scatt_len(chsym,massno,ifind)
      character*3 chsym,chtest
      integer*4 massno,masstest
      include 'sears91_gudrun.inc'
      nscatt=mscatt
      ifind=0
      do it=1,nscatt
            chtest=chsymbol(it)
            masstest=chmassno(it)
            if(chsym.eq.chtest.and.massno.eq.masstest) then
                  ifind=1
                  scatt_len=chbscat(it)
                  return
            endif
      end do
      scatt_len=0.0
      return
      end

      subroutine fractcoord(a,b)
      a=a/b
      return
      end

      subroutine wholecoord(a,b)
      a=a*b
      return
      end

      subroutine generate_fractional_coordinates(x,y,z)
      include 'lattice.inc'
      real :: fa,fb,fc
c
c subroutine to generate the fractional coordinates to an arbitrary set of lattice vectors 
c a,b,c, for input Cartesian coordinates.
c
c fractional coordinates relative to (a,b,c)
c
      fc=z/clatz
      fb=(y-fc*claty)/blaty
      fa=(x-fb*blatx-fc*clatx)/alatx
c
c minimum image convention
c
      call minimage(fa,fb,fc)
c
c calculate corrected coordinates
c
      z=fc
      y=fb
      x=fa
   
      return
      end

      subroutine generate_cartesian_coordinates(x,y,z)
      include 'lattice.inc'
      real :: fa,fb,fc
c
c subroutine to generate the fractional coordinates to an arbitrary set of lattice vectors 
c a,b,c, for input Cartesian coordinates.
c
c fractional coordinates relative to (a,b,c)
c
      fc=z
      fb=y
      fa=x
c
c calculate Cartesian coordinates
c
      z=fc*clatz
      y=fb*blaty+fc*claty
      x=fa*alatx+fb*blatx+fc*clatx
   
      return
      end

      subroutine centreofmassi(im)

      use molecom_mod

c For internal atom moves we only the adjust the atom coordinates for any
c change in the centre of mass

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
c
c determine coordinates of molecule centre of mass
c
      x1=0.0
      y1=0.0
      z1=0.0
      summass=0.0
      do ia=naf(im),nal(im)
            if(qtype(ia).eq.0) then
               amass=amassa(itype(ia))
               x1=x1+xa(ia)*amass
               y1=y1+ya(ia)*amass
               z1=z1+za(ia)*amass
               summass=summass+amass
            endif
      end do
      if(summass.gt.0.0) then
         x1=x1/summass
         y1=y1/summass
         z1=z1/summass
      end if
c
c adjust coordinates of atoms to put centre of mass at zero
c
      do ia=naf(im),nal(im)
            xa(ia)=xa(ia)-x1
            ya(ia)=ya(ia)-y1
            za(ia)=za(ia)-z1
      end do

      return
      end

      subroutine centreofmass(im)

      use molecom_mod

c For whole molecule moves or moves where part of a molecule is rotated we need
c to adjust the centre of mass coordinates as well as the atom coordinates

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      
      real coordfactor,x1,y1,z1
c
c determine coordinates of molecule centre of mass
c
      x1=0.0
      y1=0.0
      z1=0.0
      summass=0.0
      do ia=naf(im),nal(im)
            if(qtype(ia).eq.0) then
               amass=amassa(itype(ia))
               x1=x1+xa(ia)*amass
               y1=y1+ya(ia)*amass
               z1=z1+za(ia)*amass
               summass=summass+amass
            endif
      end do
      if(summass.gt.0.0) then
         x1=x1/summass
         y1=y1/summass
         z1=z1/summass
      end if
c
c adjust coordinates of atoms to put centre of mass at zero
! Remember to adjust for molecular size factor if relevant!
c
!      coordfactor=dmaxmol/dmax
!      x1=x1*coordfactor
!      y1=y1*coordfactor
!      z1=z1*coordfactor
      do ia=naf(im),nal(im)
            xa(ia)=xa(ia)-x1
            ya(ia)=ya(ia)-y1
            za(ia)=za(ia)-z1
      end do

c Adjust the centre of mass coordinates

      xadd(im)=xadd(im)+x1
      yadd(im)=yadd(im)+y1
      zadd(im)=zadd(im)+z1

      return
      end

      function en_tethered(im,ia)

      use molecom_mod

c With tethered molecules the energy simply derives from a harmonic force law wrt the tethering position

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      real en_tethered,xdif,ydif,zdif
      integer ia,iaref

c Get the position of this molecule relative to the tether position

      xdif=xadd(im)-phix(im)
      ydif=yadd(im)-phiy(im)
      zdif=zadd(im)-phiz(im)
      if(ia.gt.0.and.ia.le.nato(im)) then
         iaref=naf(im)+ia-1
         xdif=xdif+xa(iaref)
         ydif=ydif+ya(iaref)
         zdif=zdif+za(iaref)
      end if

c Calculate the change in energy 

      call minimagec(xdif,ydif,zdif)
      en_tethered=(xdif*xdif+ydif*ydif+zdif*zdif)

      return    
      end

      function gauss6pot(r,sig,sigg)

      real gauss6pot,sig,sigg,r,coul,london,arg,ampgln,x,term

      x=r/sig
      coul=chsph(0.5,0.5,x)
      london=-(coul**6)
      if(sigg.gt.0.29) sigg=0.29
      ampgln=0.5/(sigg*sigg)
      arg=x/sigg
      arg=0.5*arg*arg-ampgln
      term=0.0
      if(arg.lt.30.0) term=exp(-arg)
      gauss6pot=term+london

      return
      end

      function chsph(sig1,sig2,r)
      ra=sig1
      rb=sig2
      if(sig2.gt.sig1) then
            ra=sig2
            rb=sig1
      endif
      if(r.ge.(ra+rb)) then
            chsph=1.0/r
      else if(r.le.(ra-rb)) then
            if(ra.gt.0.0) then
                  chsph=(3.0*(ra*ra-0.2*rb*rb)-r*r)/(2.0*ra*ra*ra)
            else
                  chsph=1.0e20
            endif
      else
            aup1=ra-r
            alo1=-rb
            aup2=+rb
            alo2=ra-r
            chsph=sum1(r,ra,rb,aup1,alo1)+sum2(r,ra,rb,aup2,alo2)
      endif
      end

      function sum1(r,a,b,aup,alo)
      real*4 cof(10)
      a2=a*a
      b2=b*b
      r2=r*r
      cof(6)=1.0/6.0
      cof(5)=3.0*r/(5.0)
      cof(4)=-(3*(a2-r2)+b2)/4.0
      cof(3)=-r*(3.0*(a2+b2)-r2)/3.0
      cof(2)=3.0*b2*(a2-r2)/2.0
      cof(1)=r*b2*(3.0*a2-r2)
      fac=3/(8.0*a2*a*b2*b*r)
      sumu=0
      suml=0
      do i=1,6
            ir=7-i
            sumu=aup*(sumu+cof(ir))
            suml=alo*(suml+cof(ir))
      end do
      sum1=fac*(sumu-suml)
      return
      end

      function sum2(r,a,b,aup,alo)
      real*4 cof(10)
      a2=a*a
      b2=b*b
      r2=r*r
      cof(3)=-1.0/(3.0)
      cof(2)=0.0
      cof(1)=b2
      fac=3/(4*r*b2*b)
      sumu=0
      suml=0
      do i=1,3
            ir=4-i
            sumu=aup*(sumu+cof(ir))
            suml=alo*(suml+cof(ir))
      end do
      sum2=fac*(sumu-suml)
      return
      end

      subroutine resetq(im)

      use molecom_mod
c
c resets position of Q-type atoms within molecule relative to nuclei
c
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
c
c step through atoms of molecule
c
      do ia=naf(im),nal(im)
            nb=nbond(ia)
c
c if atom is type q, calculate revised coordinates
c
            if(qtype(ia).eq.1.and.nb.gt.0) then
                  wtsum=0.
                  xsum=0.0
                  ysum=0.0
                  zsum=0.0          
                  do ib=1,nb
                        ja=iabond(ia,ib)
                        tb=abs(tbond(ia,ib))
                        wtsum=wtsum+tb
                        xsum=xsum+xa(ja)*tb
                        ysum=ysum+ya(ja)*tb
                        zsum=zsum+za(ja)*tb
                  end do
                  if(wtsum.gt.0.0) then
                        xa(ia)=xsum/wtsum
                        ya(ia)=ysum/wtsum
                        za(ia)=zsum/wtsum
                  endif
            endif
      end do
      return
      end

      function truncate(xx,rminpt,rmaxpt)
c
c potential truncation function. Potential is set to zero beyond rmaxpt. rminpt will
c represent either the decay of the truncation function or the point at which
c the truncation falls to 0.5
c
      real truncate
      pi=4.0*atan(1.0)
      decay=(rmaxpt-rminpt)
      if(xx.ge.rmaxpt) then
         truncate=0.0
         return
      else if(xx.le.rminpt) then
         truncate=1.0
         return
      else
         rat=(xx-rminpt)*pi/decay
         truncate=0.5*(1.0+cos(rat))
         return
      endif
      end

      function truncatec(xx,rmaxpt)
c
c Charged clouds truncation function according to Hummer et. al., JPCM, vol 6 p A141 (1994)
c
      real truncatec
      if(xx.ge.rmaxpt) then
         truncatec=0.0
         return
      else
         rat=xx/rmaxpt
         truncatec=(1.0-rat)**4
         truncatec=truncatec*(1+rat*(1.6+rat*0.4))
         return
      endif
      end

      function truncatecweeks(xx,widthptc)
c
c Charged clouds truncation function according to Hummer et. al., JPCM, vol 6 p A141 (1994)
c
      real truncatecweeks,xx,widthptc
      truncatecweeks=0.0
      if(xx.lt.5.0*widthptc) then
         rat=xx/widthptc
         truncatecweeks=erfc(rat)
      endif
      return
      end

      subroutine makemole()

      use molecom_mod

      INTERFACE 
         subroutine get_rotgrp(ir1,ir2,labrott,nrott)
            integer :: ir1,ir2,nrott
            integer, dimension(:), allocatable, intent(out)  :: labrott
         end subroutine get_rotgrp
      END INTERFACE

      INTERFACE 
         subroutine get_whole_rotgrp(im,ir1,ir2,labrott,nrott)
            integer :: im,ir1,ir2,nrott
            integer, dimension(:), allocatable, intent(out)  :: labrott
         end subroutine get_whole_rotgrp
      END INTERFACE

c
c makes an .ato file from a .mol template file
c
c 7-10-2006 Modified to allow dihedral angles to be specified.

! 07-01-2014 Modified to treat specified dihedral angles as rotational groups, but with the specified 
! dihedral angle written into the .ato file, to be used by EPSR (and fmole possibly?).

c 22-10-2006 Modified to allow vibrational temperature to be specified and output

      parameter (mstarlabel=2000)
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      character*256 readfile,first,fnameatm,ruler
      character*20 ext,filetype
      character*1 iaccept

      integer ilenreadfile,ilenext,ilenfiletype
      character*20 a1,offset(mstarlabel)
      character*3 label(mstarlabel),alabelsave(mstarlabel)
      character*3 la1,la2,la3,la4
c labels for bonds defined in input stream
     *,lbond1(mstarlabel),lbond2(mstarlabel)
c labels for angles defined in input stream - langle0 is the centre of a triplet
     *,langle1(mstarlabel),langle2(mstarlabel),langle3(mstarlabel)
c old labels to be changed
     *,labelold(mstarlabel)
c labels for potential parameters and chemical symbol
      character*3 plabel(mcomp),cread(mcomp)
c labels for charge radii
     *,qrlabel(mcomp)
c labels to be changed on writing the .ato file
     *,labelnew(mstarlabel)

c row and column coordinates of stars (bonds) and labels (atoms)
      integer irstar(mstarlabel),icstar(mstarlabel)
     *,irlabel(mstarlabel),iclabel(mstarlabel)
c atom numbers and number of atoms for a given label
     *,ianumber(9,mstarlabel),ilabel(mstarlabel)
c atom numbers for extra distances to be used to constrain particular atom pairs
     *,iaextra1(mstarlabel),iaextra2(mstarlabel)
c atom numbers to define the axis of rotational groups
     *,iarot1(mstarlabel),iarot2(mstarlabel)
     *,iawrot1(mstarlabel),iawrot2(mstarlabel)
     *,iairot1(mstarlabel),iairot2(mstarlabel)
c atom numbers to define dihedral angles
     *,ndread,idread(4,mstarlabel)
     *,ilenfirst,iatomlabeladd
c temporary rotational groups
      integer                                :: nrott1,nrott2
      integer, dimension(:), allocatable     :: labrott1,labrott2
 
      real*4 bread(mstarlabel),aread(mstarlabel),eread(mstarlabel)
c read values for epsilon, sigma, mass, Coulomb charge, and charge radius
     *,pread(mcomp),sread(mcomp),mread(mcomp),qread(mcomp),qrread(mcomp)
c dihedral angles
     *,dread(mstarlabel)

      logical unitcelldefined,found

      external cosd,sind
c
c initialise random number generator
c
c     call initran1(10)
c
c generate the name of the .mol file
c
      ext='.mol'
      ilenext=4
      if(nwords.lt.2) then
         nwords=0
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499) ext(1:ilenext)
499   format(/'makemole> No ',a,' file specified! - routine will exit'/)
               return
            endif
         end do
         readfile=sometext(ncf(1):ncl(nwords))
         ilenreadfile=ncl(1)-ncf(1)+1
      else
         readfile=sometext(ncf(2):ncl(nwords))
         ilenreadfile=ncl(2)-ncf(2)+1
      endif

c Get the extension that was supplied - .mol files are treated differently from .gmol files

      ilenuptodot=indexlastdot(readfile)-1
      if(ilenuptodot.gt.0) then

         ext=readfile(ilenuptodot+1:ilenreadfile)
         ilenext=ilenreadfile-ilenuptodot
         if(ext(1:ilenext).ne.'.mol') then
        
            write(6,445) ext(1:ilenext)
 445        format(/'makemole> Unable to read files of type: ',a)
            return
            
         endif
         readfile=readfile(1:ilenuptodot)
         ilenreadfile=ilenuptodot

      endif

      fnameato=directory_name(1:nchardir)
     *//readfile(1:ilenreadfile)//ext(1:ilenext)
      ilenfnameato=nchardir+ilenreadfile+ilenext

c
c Attempt to open .mol file
c
      open(10,file=fnameato,status='old',iostat=ierr)

c Any kind of problems and the program returns to calling routine
c and does nothing

      if(ierr.ne.0) then

         write(6,444) fnameato(1:ilenfnameato)
444   format(/'makemole> Could not open makemole file : ',a)
         return

      endif

c Need to determine what file type this is

      read(10,100,iostat=ierr) sometext
      call parse(sometext,nwords,ncf,ncl,200,nstsize,nworddim)

      if(index(sometext(ncf(1):ncl(1)),'.').eq.1) then

         filetype=sometext(ncf(1):ncl(1))
         ilenfiletype=ncl(1)-ncf(1)+1
         if(nwords.gt.1) then

            read(sometext(ncf(2):ncl(2)),*) iatomlabeladd

         else

            write(6,447) filetype(1:ilenfiletype)
 447        format(/'makemole> First line must have two numbers '
     +,'for file type: ',a)
            return

         endif

      else

c Assume original .mol file type

         filetype='.mol'
         ilenfiletype=4
         iatomlabeladd=0
c       
c first number represents the number of characters in a column (including the delimiting space 
c between columns)
c
         read(sometext(ncf(1):ncl(1)),*) nspace

      endif

c
c first stage is to read the molecular template file and generate
c the atoms and bonds reference integers
c
      if(filetype(1:ilenfiletype).eq.'.mol') then
c
c read a ruler line - this can be used to show where the columns occur
c
         read(10,100,iostat=ierr) ruler
         nrow=0
         nrowmax=0
         ncolumnmax=0
         nstar=0
         
      endif

      write(6,446) filetype(1:ilenfiletype),iatomlabeladd
 446  format(/'makemole> File type: ',a,' Atom label add: ',i2)
      
      nlabel=0
      ierr=0
c set the number of atom types to zero at the outset
      ntype=0
c no. of bond values read on input
      nbread=0
c no. of bond angle values read on input
      naread=0
c no. of extra distances read on input
      neread=0
c no. of rotational groups read on input
      nrread=0
      nwrread=0
      nirread=0
c no. of potential values read on input
      npread=0
c no. of charge radius values read on input
      nqrread=0
c initial temperature
      simtk=300
c initial vibrational temperature
      tempb=65000.0
      nvread=0
c initial bond angle temperature (relative)
      tempang=1.0
      nangread=0
c initial dihedral temperature
      tempdih=10.0
      ndihread=0
c initial atom number density
      rho=0.0
c no. of labels to change before writing the .ato file
      nchread=0
c no. of dihedral angles on input
      ndread=0
c assumed values of ecore and dcore unless changed by the input. 
c These are used in the molecule generation program fmole.
      ecore=0.0
      dcore=0.0
      unitcelldefined=.false.
      ia=0
      do while (ierr.eq.0)
c
c read a line of text
c
         sometext=' '
         read(10,100,iostat=ierr) sometext
100   format(a)
c
c parse the text into words. Any gaps of more than 200 spaces will 
c signal the end of the line
c
         call parse(sometext,nwords,ncf,ncl,200,nstsize,nworddim)
         if(nwords.gt.0) then
c
c check that the first word is not a control word
c
         first=sometext(ncf(1):ncl(1))
         ilenfirst=ncl(1)-ncf(1)+1

         if(first(1:ilenfirst).eq.'atom') then

            if(nwords.gt.6.and.naread.lt.mstarlabel) then

             if(ia.ge.mstarlabel) then
                
                write(6,448) mstarlabel
 448            format(/'makemole> Attempted to create more than '
     +,i5,' atoms')
                return

             endif
             
             ia=ia+1
c Ignore the atom number from ghemical
             isref=3
             alabel(ia)=sometext(ncf(isref):ncl(isref))

c Store this label in case it is changed

             alabelsave(ia)=alabel(ia)
c
c set the x and y coordinates according to position in the map
c
             isref=isref+1
             read(sometext(ncf(isref):ncl(isref)),*) xa(ia)
             isref=isref+1
             read(sometext(ncf(isref):ncl(isref)),*) ya(ia)
             isref=isref+1
             read(sometext(ncf(isref):ncl(isref)),*) za(ia)
             isref=isref+1
             read(sometext(ncf(isref):ncl(isref)),*) nbond(ia)
             ib=0
             do while (ib.lt.nbond(ia).and.ib.lt.near.and.isref.lt.nwords)
                
                ib=ib+1
                isref=isref+1
                read(sometext(ncf(isref):ncl(isref)),*) iabond(ia,ib)
                iabond(ia,ib)=iabond(ia,ib)+iatomlabeladd
                
             end do
             nbond(ia)=ib

            endif

         else if(first(1:ilenfirst).eq.'bond') then
c
c assume the next two words are atom type labels, followed by the corresponding
c bond length in A.
c
            if(nwords.ge.4.and.nbread.lt.mstarlabel) then
               nbread=nbread+1
               lbond1(nbread)=sometext(ncf(2):ncl(2))
               lbond2(nbread)=sometext(ncf(3):ncl(3))
               read(sometext(ncf(4):ncl(4)),*) bread(nbread)
            endif

         else if(first(1:ilenfirst).eq.'angle') then
c
c assume the next three words are atom type labels, 
c with the middle one being the one about which the angle subtends
c followed by the corresponding angle in degrees.
c
            if(nwords.ge.5.and.naread.lt.mstarlabel) then
               naread=naread+1
               langle1(naread)=sometext(ncf(2):ncl(2))
               langle2(naread)=sometext(ncf(3):ncl(3))
               langle3(naread)=sometext(ncf(4):ncl(4))
               read(sometext(ncf(5):ncl(5)),*) aread(naread)
            endif

         else if(first(1:ilenfirst).eq.'extra') then
c
c Additional distances that may needed to keep the molecule intact
c assume the next two words are atom type labels, followed by the corresponding
c expected separation in A.
c
            if(nwords.ge.4.and.neread.lt.mstarlabel) then
               neread=neread+1
               read(sometext(ncf(2):ncl(2)),*) iaextra1(neread) 
               read(sometext(ncf(3):ncl(3)),*) iaextra2(neread)
               read(sometext(ncf(4):ncl(4)),*) eread(neread)
               iaextra1(neread)=iaextra1(neread)+iatomlabeladd
               iaextra2(neread)=iaextra2(neread)+iatomlabeladd
            endif

         else if(first(1:ilenfirst).eq.'rot') then
c
c Read two atom numbers to define the axis for rotational groups
c
            if(nwords.ge.3.and.nrread.lt.mstarlabel) then
               nrread=nrread+1
               read(sometext(ncf(2):ncl(2)),*) iarot1(nrread) 
               read(sometext(ncf(3):ncl(3)),*) iarot2(nrread)
               iarot1(nrread)=iarot1(nrread)+iatomlabeladd
               iarot2(nrread)=iarot2(nrread)+iatomlabeladd
            endif

         else if(first(1:ilenfirst).eq.'wholerot') then
c
c Read two atom numbers to define the axis for whole rotational groups
c This is where a whole molecule is rotated about an axis defined by
c two atoms within the molecule
c
            if(nwords.ge.3.and.nwrread.lt.mstarlabel) then
               nwrread=nwrread+1
               read(sometext(ncf(2):ncl(2)),*) iawrot1(nwrread) 
               read(sometext(ncf(3):ncl(3)),*) iawrot2(nwrread)
               iawrot1(nwrread)=iawrot1(nwrread)+iatomlabeladd
               iawrot2(nwrread)=iawrot2(nwrread)+iatomlabeladd
            endif

         else if(first(1:ilenfirst).eq.'introt') then
c
c Read two atom numbers to define the axis for rotational groups
c
            if(nwords.ge.3.and.nirread.lt.mstarlabel) then
               nirread=nirread+1
               read(sometext(ncf(2):ncl(2)),*) iairot1(nirread) 
               read(sometext(ncf(3):ncl(3)),*) iairot2(nirread)
               iairot1(nirread)=iairot1(nirread)+iatomlabeladd
               iairot2(nirread)=iairot2(nirread)+iatomlabeladd
            endif

         else if(first(1:ilenfirst).eq.'dihedral') then
c
c Read four atom numbers to define the dihedral angle plus the dihedral angle. 
c Atoms 2 and 3 define the axis about which the dihedral angle is calculated.
c
            if(nwords.ge.6.and.ndread.lt.mstarlabel) then
               ndread=ndread+1
               do j=1,4
                  read(sometext(ncf(j+1):ncl(j+1)),*) idread(j,ndread)
                  idread(j,ndread)=idread(j,ndread)+iatomlabeladd
               end do
               
               read(sometext(ncf(6):ncl(6)),*) dread(ndread)
            endif

         else if(first(1:ilenfirst).eq.'changelabel') then
c
c Change atom labels: read old atom label and new atom label
c these are only applied after the atom distances have been defined
c
            if(nwords.gt.2) then
               nchread=nchread+1
               labelold(nchread)=sometext(ncf(2):ncl(2))
               labelnew(nchread)=sometext(ncf(3):ncl(3))
            endif

         else if(first(1:ilenfirst).eq.'potential') then
c
c Read epsilon, sigma, mass, charge and chemical symbol of specified label
c
            if(nwords.ge.7.and.npread.lt.mcomp) then
               npread=npread+1
               plabel(npread)=sometext(ncf(2):ncl(2))
               read(sometext(ncf(3):ncl(3)),*) pread(npread) 
               read(sometext(ncf(4):ncl(4)),*) sread(npread)
               read(sometext(ncf(5):ncl(5)),*) mread(npread) 
               read(sometext(ncf(6):ncl(6)),*) qread(npread)
               cread(npread)=sometext(ncf(7):ncl(7))  
            endif

c Read charge radius for specified atom type. If not specified this will
c be set to zero

         else if(first(1:ilenfirst).eq.'qradius') then

            if(nwords.ge.3.and.nqrread.lt.mcomp) then

               nqrread=nqrread+1
               qrlabel(nqrread)=sometext(ncf(2):ncl(2))
               read(sometext(ncf(3):ncl(3)),*) qrread(nqrread)

            endif
c
c set temperature of simulation
c
         else if(first(1:ilenfirst).eq.'temperature') then

            if(nwords.gt.1) then
               read(sometext(ncf(2):ncl(2)),*) simtk
            end if

         else if(first(1:ilenfirst).eq.'vibtemp') then

            if(nwords.gt.1) then
               read(sometext(ncf(2):ncl(2)),*) tempb
               nvread=ncl(2)-ncf(2)+1
            end if

         else if(first(1:ilenfirst).eq.'angtemp') then

            if(nwords.gt.1) then
               read(sometext(ncf(2):ncl(2)),*) tempangreq
               nangread=ncl(2)-ncf(2)+1
            end if

         else if(first(1:ilenfirst).eq.'dihtemp') then

            if(nwords.gt.1) then
               read(sometext(ncf(2):ncl(2)),*) tempdihreq
               ndihread=ncl(2)-ncf(2)+1
            end if

         else if(first(1:ilenfirst).eq.'density') then

            if(nwords.gt.1) then
               read(sometext(ncf(2):ncl(2)),*) rho
            endif
c
c read the potential parameters used to find the molecule
c
         else if(first(1:ilenfirst).eq.'ecoredcore') then

            if(nwords.gt.1) read(sometext(ncf(2):ncl(2)),*) ecore
            if(nwords.gt.2) read(sometext(ncf(3):ncl(3)),*) dcore
c
c read the unit cell parameters used to set the size and shape of the box
c
         else if(first(1:ilenfirst).eq.'unitcell') then

            if(nwords.gt.6) then
               unitcelldefined=.true.
               read(sometext(ncf(2):ncl(7)),*) alat,blat,clat,phib,thetac,phic
            else
               unitcelldefined=.false.
            end if
c
c otherwise treat the words as atom labels or bonds. The first word is the row
c number and so will be ignored on input
c
         else if(nwords.gt.1) then
         nrow=nrow+1
         iw=1
         do while (iw.lt.nwords)
            iw=iw+1
c
c for each word read in, first find out if the first character is a number
c
            read(sometext(ncf(iw):ncf(iw)),*,iostat=ierrread) nread
            if(ierrread.eq.0) then
               a1=sometext(ncf(iw)+1:ncl(iw))
            else
               a1=sometext(ncf(iw):ncl(iw))
               nread=1
            endif
c
c now determine whether this character represents a bond or an atom
c
            ncolumn=(ncf(iw)-4+nspace-1)/nspace
            if(ncolumn.gt.ncolumnmax) ncolumnmax=ncolumn
            if(a1.eq.'*') then
               nstar=nstar+1
               irstar(nstar)=nrow
               icstar(nstar)=ncolumn
            else
c
c set up coordinates of label and no of atoms at this label
c
               nlabel=nlabel+1
               irlabel(nlabel)=nrow
               iclabel(nlabel)=ncolumn
               ilabel(nlabel)=nread
c
c column and row number of current label
c
               il1=nlabel
               icl1=iclabel(il1)
               irl1=irlabel(il1)
c
c if the first three characters of the label are 'map', this label
c represents a mapping to a previous label, given by the row and column offsets 
c from the present position
c
               if(a1(1:3).eq.'map') then
c
c read the row and column offsets
c
                  label(nlabel)=a1(1:3)
                  nchara1=index(a1,' ')-1
                  offset(nlabel)=a1(4:nchara1)
                  nmap=nmap+1
                  ncharoffset=nchara1-3
c
c test to see if this is an absolute reference
c
                  ilenopenbra=index(offset(nlabel),'(')
                  ilencomma=index(offset(nlabel),',')
                  ilenclosebra=index(offset(nlabel),')')
                  ilenread1=ilenopenbra+ilencomma+ilenclosebra
c
c to be a valid reference the sum of the positions of these three characters must
c be at least 1 + 3 + 5 = 9, otherwise there is no room for numbers
c
                  if(ilenread1.gt.8) then
c
c absolute reference
c
                     read(offset(nlabel)(ilenopenbra+1:ilencomma-1)
     *,*,iostat=ierr) irmap
                     read(offset(nlabel)(ilencomma+1:ilenclosebra-1)
     *,*,iostat=ierr) icmap
                     if(ierr.ne.0.or.irmap.le.0.or.icmap.le.0) then
                        write(6,111) a1(1:nchara1),irl1,icl1
                        return
                     endif
                  else
c
c alternatively try see if this is relative reference
c
                     ilenplus=index(offset(nlabel)(2:ncharoffset),'+')
                     ilenminus=index(offset(nlabel)(2:ncharoffset),'-')
                     ilenread=ilenplus+ilenminus
                     if(ilenread.gt.1) then
                        ilenread=ilenread+1
c
c read the offsets
c
                        read(offset(nlabel)(1:ilenread-1)
     *,*,iostat=ierr) nrowoffset
                        read(offset(nlabel)(ilenread:ncharoffset)
     *,*,iostat=ierr) ncoloffset
                        if(ierr.ne.0) then
                           write(6,111) a1(1:nchara1),irl1,icl1
                           return
                        endif
c     write(6,*) irl1,icl1,nrowoffset,ncoloffset
c
c calculate the row and column of the label this map refers to
c
                        irmap=irl1+nrowoffset
                        icmap=icl1+ncoloffset
                     else
                        write(6,111) a1(1:nchara1),irl1,icl1
111   format(/'makemole> Unrecognised format for mapped atom ',a
     *,' at row ',i3,' and column ',i3
     *//' - check the .mol file')
                        return
                     endif
                  endif
c
c find the previous label corresponding to this label - if not found exit with an 
c error message.
c
                  il=0
                  ifound=0
                  do while(il.lt.nlabel-1.and.ifound.eq.0)
                     il=il+1
                     if(irlabel(il).eq.irmap.and.iclabel(il).eq.icmap) ifound=il
                  end do                  
                  if(ifound.eq.0) then
                     write(6,112) a1(1:nchara1),irmap,icmap
112   format(/'makemole> Label referred to by mapped atom',a
     *,' not found at row ',i3,' and column ',i3
     *//' - check the .mol file')
                     return
                  else
                     ilabel(nlabel)=ilabel(ifound)
                     do ia1=1,ilabel(nlabel)
                        ianumber(ia1,nlabel)=ianumber(ia1,ifound)
                     end do
                  endif
               else
c
c otherwise generate the atom numbers for this new atom and zero the bond array
c
                  label(nlabel)=a1
                  do ia1=1,ilabel(nlabel)
                     ia=ia+1
                     ianumber(ia1,nlabel)=ia
                     alabel(ia)=label(nlabel)
c
c set the x and y coordinates according to position in the map
c
                     xa(ia)=0.0
                     ya(ia)=0.0
                     za(ia)=0.0
c                    call ftype(ia)
c                    if(qtype(ia).ne.1) natomr=natomr+1
                     nbond(ia)=0
                  end do
               endif
c
c step through the existing stars to find out which occur as neighbours
c in a previous column on the same row, or in the previous row
c
               is=0
               do while(is.lt.nstar)
                  is=is+1
                  irs=irstar(is)
                  ics=icstar(is)
                  irdif1=irlabel(il1)-irs
                  icdif1=iclabel(il1)-ics
                  if((irdif1.eq.1.and.abs(icdif1).lt.2)
     *.or.(irdif1.eq.0.and.icdif1.eq.1)) then
c
c o.k. found a star in the vicinity - look for another label adjecent to this star (not including
c the present) one in the previous two rows and or columns
c
                     il2=1
                     do while(il2.lt.nlabel)
                        irdif2=abs(irs-irlabel(il2))
                        icdif2=abs(ics-iclabel(il2))
                        if(irdif2.lt.2.and.icdif2.lt.2) then
c
c atom found - set up the bond.
c
c step through all the possible pairs of atoms and identify bonds
c
c     write(6,*) il1,il2
                           do ia1=1,ilabel(il1)
                              ia1ref=ianumber(ia1,il1)
                              do ia2=1,ilabel(il2)
                                 ia2ref=ianumber(ia2,il2)
c
c check that ia2ref does not already exist in the list of neighbours of ia1ref
c
                                 ib=0
                                 ifound=0
                                 do while(ib.lt.nbond(ia1ref).and.ifound.eq.0)
                                    ib=ib+1
                                    if(iabond(ia1ref,ib).eq.ia2ref) ifound=ia2ref
                                 end do
                                 if(ifound.eq.0) then
                                    nbond(ia1ref)=nbond(ia1ref)+1
                                    nbond(ia2ref)=nbond(ia2ref)+1
                                    iabond(ia1ref,nbond(ia1ref))=ia2ref
                                    iabond(ia2ref,nbond(ia2ref))=ia1ref
                                 endif
                              end do
                           end do
                        endif
                        il2=il2+1
                     end do
                  endif
               end do
            endif
         end do
         endif
         endif
      end do
      close(10)
      natom=ia
c
c write out findings so far
c 
!      write(6,*) nrow,ncolumnmax,nstar,nlabel,ntype
!      do it=1,ntype
!           write(6,*) it,' ',type(it)
!      end do
!      do ia=1,natom
!         write(6,*) ia
!         write(6,*) alabel(ia),nbond(ia),(iabond(ia,ib),ib=1,nbond(ia))
!      end do
c
c set the bond lengths
c
      do ia=1,natom
         la1=alabel(ia)
         ib=0
         do while(ib.lt.nbond(ia))
            ib=ib+1
            la2=alabel(iabond(ia,ib))
c
c search through the bonds read in (if any) and assign the appropriate bond length
c
            ibread=0
            ifound=0
            do while(ibread.lt.nbread.and.ifound.eq.0)
               ibread=ibread+1
               if(
     *(la1.eq.lbond1(ibread).and.la2.eq.lbond2(ibread))
     *.or.
     *(la1.eq.lbond2(ibread).and.la2.eq.lbond1(ibread))
     *) ifound=ibread
            end do
            ierr=1 
            do while(ifound.eq.0.and.ierr.ne.0)
               write(6,101) la1,la2
101   format(/'makemole> Give the bond'
     +,' length for atom pair ',a3,1x,a3,': ',$)
               nbread=nbread+1
               ifound=nbread
               lbond1(nbread)=la1
               lbond2(nbread)=la2
               read(5,*,iostat=ierr) bread(nbread)
            end do
            tbond(ia,ib)=bread(ifound)
         end do
      end do
      naf(1)=1
      nal(1)=natom
c
c save the number of primary bonds for this atom
c
      do ia=1,natom
         nbondp(ia)=nbond(ia)
      end do
c
c write out findings so far
c 
c     write(6,*) nrow,ncolumnmax,nstar,nlabel
c     do ia=1,natom
c        write(6,*) ia,alabel(ia),nbond(ia),(iabond(ia,ib),tbond(ia,ib),ib=1,nbond(ia))
c     end do
c
c now step through the atoms and use the angles defined to create extra bonding distances
c
      do ia2=1,natom
         la2=alabel(ia2)
c
c must be at least two bonds for an angle to be definable
c
         ib1=0
         do while(ib1.lt.(nbondp(ia2)-1))
            ib1=ib1+1
            ia1=iabond(ia2,ib1)
            la1=alabel(ia1)
            ib2=ib1
            do while(ib2.lt.nbondp(ia2))
               ib2=ib2+1
               ia3=iabond(ia2,ib2)
               la3=alabel(ia3)
c
c before we set up this bond check that a bond between these two atoms 
c does not already exist
c
               isearch=0
               iexist=0
               do while(isearch.lt.nbond(ia1).and.iexist.eq.0)
                  isearch=isearch+1
                  if(ia3.eq.iabond(ia1,isearch)) iexist=ia3
               end do
               if(iexist.eq.0) then
                  bond21=abs(tbond(ia2,ib1))
                  bond23=abs(tbond(ia2,ib2))
c
c search for the angle between the three atoms
c
c
c search through the angles read in (if any) and assign the appropriate angle
c to this triplet
c
                  iaread=0
                  ifound=0
                  do while(iaread.lt.naread.and.ifound.eq.0)
                   iaread=iaread+1
                   if(
     *(la1.eq.langle1(iaread).and.la2.eq.langle2(iaread)
     *.and.la3.eq.langle3(iaread))
     *.or.
     *(la1.eq.langle3(iaread).and.la2.eq.langle2(iaread)
     *.and.la3.eq.langle1(iaread))
     *) ifound=iaread
                  end do
                  ierr=1 
                  do while(ifound.eq.0.and.ierr.ne.0)
                     write(6,102) la1,la2,la3
 102  format(/'makemole> Give the angle subtended by the triplet '
     *,a3,1x,a3,1x,a3,': ',$)
                       naread=naread+1
                   ifound=naread
                   langle1(naread)=la1
                   langle2(naread)=la2
                   langle3(naread)=la3
                   read(5,*,iostat=ierr) aread(naread)
                  end do
c
c Set up distance between ib1 and ib3
c
                  bangle=aread(ifound)
                  b13=bond21**2+bond23**2-2*bond21*bond23*cosd(bangle)
                  b13=sqrt(b13)
                  if(nbond(ia1).lt.near) then

                   nbond(ia1)=nbond(ia1)+1
                   iabond(ia1,nbond(ia1))=ia3
                   tbond(ia1,nbond(ia1))=-b13

                  endif
                  if(nbond(ia3).lt.near) then

                   nbond(ia3)=nbond(ia3)+1
                   iabond(ia3,nbond(ia3))=ia1
                   tbond(ia3,nbond(ia3))=-b13

                  endif
               endif
            end do
         end do
      end do
c
c save the number of secondary bonds
c
      do ia=1,natom
         nbonds(ia)=nbond(ia)
      end do

c Now setup any dihedral angles requested. These will be defined by means
c of distances between atoms

      id=0
      ir=0
      write(6,*) ndread,naread
      do while (id.lt.ndread)

         id=id+1

         id1=idread(1,id)
         id2=idread(2,id)
         id3=idread(3,id)
         id4=idread(4,id)

c First need to check that id4 is not already listed as a neighbour of id1.
c If it is then nothing is done

         ifound=0
         ib=0
         do while (ib.lt.nbond(id1).and.ifound.eq.0)
            ib=ib+1
            if(iabond(id1,ib).eq.id4) ifound=id4
         end do

         if(ifound.eq.0) then

            ir1=id2
            ir2=id3

! Need to check whether a dihedral group for this pair of axis atoms has already been defined.

            iref=0
            ifounddih=0
            do while (iref.lt.ir.and.ifounddih.eq.0)
               iref=iref+1
               if(rotgrp(iref).eq.'DIH'
     1         .and.((ir1.eq.labaxi(iref,1).and.ir2.eq.labaxi(iref,2))
     1         .or.(ir2.eq.labaxi(iref,1).and.ir1.eq.labaxi(iref,2)))) then
                  ifounddih=iref
               endif
            end do

! Start a new rotational group:
! Get the two possible rotational groups into temporary arrays

            if(ifounddih.eq.0) then
               call get_rotgrp(ir1,ir2,labrott1,nrott1)
               call get_rotgrp(ir2,ir1,labrott2,nrott2)
         write(6,*) nrott1,(labrott1(j),j=1,nrott1)
         write(6,*) nrott2,(labrott2(j),j=1,nrott2)

               ir=ir+1
               ndihedralang(ir)=1
               dihedralang(ndihedralang(ir),ir)=dread(id)
               call testsize_rotgrps(ir)
               rotgrp(ir)='DIH'
               !Use the smallest of the two groups
               if(nrott1.lt.nrott2) then
                  labaxi(ir,1)=ir1
                  labaxi(ir,2)=ir2
                  nrot(ir)=nrott1
                  rotf(ir)=nrottot+1
                  nrottot=nrottot+nrot(ir)
                  rotl(ir)=nrottot
                  call testsize_labrot()
                  do i=1,nrot(ir)
                     labref=rotf(ir)+i-1
                     labrot(labref)=labrott1(i)
                  end do
!Remember to preserve the linking to the axis atoms!!
                  labrotdih(1,ndihedralang(ir),ir)=id1
                  labrotdih(2,ndihedralang(ir),ir)=id4
               else
                  labaxi(ir,1)=ir2
                  labaxi(ir,2)=ir1
                  nrot(ir)=nrott2
                  rotf(ir)=nrottot+1
                  nrottot=nrottot+nrot(ir)
                  rotl(ir)=nrottot
                  call testsize_labrot()
                  do i=1,nrot(ir)
                     labref=rotf(ir)+i-1
                     labrot(labref)=labrott2(i)
                  end do
!Remember to preserve the linking to the axis atoms!!
                  labrotdih(1,ndihedralang(ir),ir)=id4
                  labrotdih(2,ndihedralang(ir),ir)=id1
               end if
            else
! Add the new dihedral angle to the existing list
               if(ndihedralang(ifounddih).lt.near) then
                  ndihedralang(ifounddih)=ndihedralang(ifounddih)+1
                  dihedralang(ndihedralang(ifounddih),ifounddih)=dread(id)
                  if(ir1.eq.labaxi(iref,1).and.ir2.eq.labaxi(iref,2)) then
                     labrotdih(1,ndihedralang(ifounddih),ifounddih)=id1
                     labrotdih(2,ndihedralang(ifounddih),ifounddih)=id4
                  else if((ir2.eq.labaxi(iref,1).and.ir1.eq.labaxi(iref,2))) then
                     labrotdih(1,ndihedralang(ifounddih),ifounddih)=id4
                     labrotdih(2,ndihedralang(ifounddih),ifounddih)=id1
                  end if
               end if
            end if
         end if
      end do
      ngrot(1)=ir
! If no dihedral groups are defined, then set ecore and dcore. This will help to stop overlaps in long chain molecules.
      if(ir.eq.0.and.ecore.eq.0.0) then
         ecore=1.0
         dcore=2.0
      end if
         
c         else if((ierr1.eq.0.and.ierr2.ne.0).or.(nrott1.le.nrott2)) then

c            labaxi(ir,1)=ir2
c            labaxi(ir,2)=ir1
c            nrot(ir)=nrott1
c            do i=1,nrot(ir)

c             labrot(ir,i)=labrott1(i)

c            end do

c     write(6,*) id,id1,id2,id3,id4,ifound

c Find the bond lengths 1-2, 2-3, and 3-4.

!         ib=0
!         bond12=0.0
!         do while (ib.lt.nbond(id1).and.bond12.eq.0.0)
!            ib=ib+1
!            if(iabond(id1,ib).eq.id2) bond12=tbond(id1,ib)
!         end do
                        
!         ib=0
!         bond23=0.0
!         do while (ib.lt.nbond(id2).and.bond23.eq.0.0)
!            ib=ib+1
!            if(iabond(id2,ib).eq.id3) bond23=tbond(id2,ib)
!         end do

!         ib=0
!         bond34=0.0
!         do while (ib.lt.nbond(id3).and.bond34.eq.0.0)
!            ib=ib+1
!            if(iabond(id3,ib).eq.id4) bond34=tbond(id3,ib)
!         end do

c     write(6,*) id,id1,id2,id3,id4,bond12,bond23,bond34

!         if(bond12.gt.0.0.and.bond23.gt.0.0.and.bond34.gt.0.0) then

c Search for the angles 1 - 2 - 3 and 2 - 3 - 4

!            angle123=0.0
!            angle234=0.0
!            iaread=0
!            do while(iaread.lt.naread.and.(angle123.eq.0.0.or.angle234.eq.0.0))

!               iaread=iaread+1
c     write(6,*) la1,la2,la3,la4,iaread,naread
c     write(6,*) langle1(iaread),langle2(iaread),langle3(iaread)
!               if(
!     *(la1.eq.langle1(iaread).and.la2.eq.langle2(iaread)
!     *.and.la3.eq.langle3(iaread))
!     *.or.
!     *(la1.eq.langle3(iaread).and.la2.eq.langle2(iaread)
!     *.and.la3.eq.langle1(iaread))) angle123=aread(iaread)
!               if(
!     *(la2.eq.langle1(iaread).and.la3.eq.langle2(iaread)
!     *.and.la4.eq.langle3(iaread))
!     *.or.
!     *(la2.eq.langle3(iaread).and.la3.eq.langle2(iaread)
!     *.and.la4.eq.langle1(iaread))) angle234=aread(iaread)

!            end do

c     write(6,*) id,id1,id2,id3,id4,angle123,angle234

c Calculate the distance atom id1 to atom id4

!            if(angle123.gt.0.0.and.angle234.gt.0.0) then

!               x1=0.0
!               x4=bond23-bond12*cosd(angle123)-bond34*cosd(angle234)
!               y1=bond12*sind(angle123)
!               y4=bond34*sind(angle234)
!               z1=0.0
!               z4=y4*sind(dread(id))
!               y4=y4*cosd(dread(id))

!               xdiff=x4-x1
!               ydiff=y4-y1
!               zdiff=z4-z1
!               dtest=sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff)

!               if(nbond(id1).lt.near) then
                  
!                  nbond(id1)=nbond(id1)+1
!                  iabond(id1,nbond(id1))=id4
!                  tbond(id1,nbond(id1))=dtest

!               endif
!               if(nbond(id4).lt.near) then
                  
!                  nbond(id4)=nbond(id4)+1
!                  iabond(id4,nbond(id4))=id1
!                  tbond(id4,nbond(id4))=dtest
                  
!               endif
!            endif
!         endif
!      end do

c finally setup any extra distances between particular atoms that have been specified

      ia1=0
      do while(ia1.lt.natom-1)
         ia1=ia1+1
         ia2=ia1
         do while(ia2.lt.natom)
            ia2=ia2+1
c
c before we set up this bond check that a bond between these two atoms 
c does not already exist
c
            isearch=0
            iexist=0
            do while(isearch.lt.nbond(ia1).and.iexist.eq.0)
               isearch=isearch+1
               if(ia2.eq.iabond(ia1,isearch)) iexist=ia2
            end do
            if(iexist.eq.0) then
               ifound=0
               ieread=0
               do while(ieread.lt.neread.and.ifound.eq.0)
                  ieread=ieread+1
                  if(
     *(ia1.eq.iaextra1(ieread).and.ia2.eq.iaextra2(ieread))
     *.or.
     *(ia1.eq.iaextra2(ieread).and.ia2.eq.iaextra1(ieread))
     *) ifound=ieread
               end do
               if(ifound.ne.0) then
                  if(nbond(ia1).lt.near) then

                   nbond(ia1)=nbond(ia1)+1
                   iabond(ia1,nbond(ia1))=ia2
                   tbond(ia1,nbond(ia1))=-abs(eread(ifound))

                  endif
                  if(nbond(ia2).lt.near) then

                   nbond(ia2)=nbond(ia2)+1
                   iabond(ia2,nbond(ia2))=ia1
                   tbond(ia2,nbond(ia2))=-abs(eread(ifound))

                  endif
               endif
            endif
         end do
      end do
c
c now set up the rotational groups
c
      iref=0
      write(6,*) 'Getting rotational groups'
      do while(iref.lt.nrread)
         iref=iref+1
         ir1=iarot1(iref)
         ir2=iarot2(iref)

c Get the two possible rotational groups into temporary arrays

         call get_rotgrp(ir1,ir2,labrott1,nrott1)
         call get_rotgrp(ir2,ir1,labrott2,nrott2)
         ir=ir+1
         call testsize_rotgrps(ir)
         rotgrp(ir)='ROT'
         !Use the smallest of the two groups
         if(nrott1.lt.nrott2) then
            labaxi(ir,1)=ir1
            labaxi(ir,2)=ir2
            nrot(ir)=nrott1
            rotf(ir)=nrottot+1
            nrottot=nrottot+nrot(ir)
            rotl(ir)=nrottot
            call testsize_labrot()
            do i=1,nrot(ir)
               labref=rotf(ir)+i-1
               labrot(labref)=labrott1(i)
            end do
         else
            labaxi(ir,1)=ir2
            labaxi(ir,2)=ir1
            nrot(ir)=nrott2
            rotf(ir)=nrottot+1
            nrottot=nrottot+nrot(ir)
            rotl(ir)=nrottot
            call testsize_labrot()
            do i=1,nrot(ir)
               labref=rotf(ir)+i-1
               labrot(labref)=labrott2(i)
            end do
         end if
         ndihedralang(ir)=0
      end do

c Set up the whole molecule rotational groups

      iref=0
      do while(iref.lt.nwrread)
         iref=iref+1
         ir1=iawrot1(iref)
         ir2=iawrot2(iref)
!      write(6,*) ir,iref,ir1,ir2
! Get all the atoms in the molecule except the two atoms that form the axis.
         call get_whole_rotgrp(1,ir1,ir2,labrott2,nrott2)
         ir=ir+1
         rotgrp(ir)='ROT'
         call testsize_rotgrps(ir)
         labaxi(ir,1)=ir1
         labaxi(ir,2)=ir2
         nrot(ir)=nrott2
         rotf(ir)=nrottot+1
         nrottot=nrottot+nrot(ir)
         rotl(ir)=nrottot
         call testsize_labrot()
         do i=1,nrot(ir)
            labref=rotf(ir)+i-1
            labrot(labref)=labrott2(i)
         end do
         ndihedralang(ir)=0
      end do

c Set up the internal rotational groups: this is between two atoms which are not bonded to each other and
c aims to rotate all the atoms between them. Note that if 
c axis atoms ir1 and ir2 have more than one bond each to the group, then the group cannot be rotated, so should
c be discarded.

      iref=0
      do while(iref.lt.nirread)
         iref=iref+1
         ir1=iairot1(iref)
         ir2=iairot2(iref)
! Get all the atoms in the molecule except the two atoms that form the axis.
         call get_rotgrp(ir1,ir2,labrott1,nrott1)
!         write(6,*) ir1,ir2,nrott1
!         write(6,*) (labrott1(j),j=1,nrott1)
         call get_rotgrp(ir2,ir1,labrott2,nrott2)
!         write(6,*) ir2,ir1,nrott2
!         write(6,*) (labrott2(j),j=1,nrott2)
         !Now compare the two groups and keep atoms which are common to both groups.
         found=.false.
         i1=0
         do while (i1.lt.nrott1)
            i1=i1+1
            found=.false.
            i2=0
            do while(i2.lt.nrott2.and..not.found)
               i2=i2+1
               found=labrott1(i1).eq.labrott2(i2)
            end do
            !If this atom was not found in the second group, we remove it from the first group
            if(.not.found) then
               nrott1=nrott1-1
               do i=i1,nrott1
                  labrott1(i)=labrott1(i+1)
               end do
               i1=i1-1
            end if
         end do
         write(6,*) ir1,ir2,nrott1
         write(6,*) (labrott1(j),j=1,nrott1)
         !need to check that only 1 atom of the group is a neighbour of ir1. Similarly for ir2
         !but there is no easy way to do this since the list of bonds includes bond angles and dihedral angles.
         write(6,*) ir1,ir2,nrott1
         write(6,*) (labrott1(j),j=1,nrott1)
         ir=ir+1
         call testsize_rotgrps(ir)
         rotgrp(ir)='ROT'
         labaxi(ir,1)=ir1
         labaxi(ir,2)=ir2
         nrot(ir)=nrott1
         rotf(ir)=nrottot+1
         nrottot=nrottot+nrot(ir)
         rotl(ir)=nrottot
         call testsize_labrot()
         do i=1,nrot(ir)
            labref=rotf(ir)+i-1
            labrot(labref)=labrott1(i)
         end do
         ndihedralang(ir)=0
      end do
      ngrot(1)=ir
      ngrrot=ir
c
c check that the temperature and density have been set
c
      if(simtk.le.0.0) then
         write(6,116,advance='no')
116   format(/'makemole> Give the temperature'
     *,' (K) of the simulation: ')
         ierr=1
         do while(ierr.ne.0)
            read(5,*,iostat=ierr) simtk
         end do
      endif
      if(rho.le.0.0.and..not.unitcelldefined) then
         write(6,117,advance='no')
117   format(/'makemole> Give the atomic number'
     +,' density of the simulation [ats/A**3]: ')
         ierr=1
         do while(ierr.ne.0)
            read(5,*,iostat=ierr) rho
         end do
      endif
c
c write out findings so far
c 
c     write(6,*) nrow,ncolumnmax,nstar,nlabel
c     do ia=1,natom
c        write(6,*) ia,alabel(ia),nbond(ia)
c     +,(iabond(ia,ib),tbond(ia,ib),ib=1,nbond(ia))
c     end do
c     write(6,*) ngrot(1)
c     do ig=1,ngrot(1)
c        write(6,*) labaxi(ig,1),labaxi(ig,2)
c        write(6,*) nrot(ig),(labrot(ig,j),j=1,nrot(ig))
c     end do
c
c now prepare to write the .ato file
c
c
c set some parameters ready for writing the .ato file
c
      nmol=1
      xadd(1)=0.0
      yadd(1)=0.0
      zadd(1)=0.0
      phix(1)=0.0
      phiy(1)=0.0
      phiz(1)=0.0
      tethered(1)=.false.
      nato(1)=natom
      stepmi=0.1
      if(ngrot(1).gt.0) then
            stepri=1.0
      else
            stepri=0.0
      endif
      stepr=0.3
      stepm=1.0

c Initialise the atom types

      ntype=0
      natomr=0

      do ia=naf(1),nal(1)
c
c change the atom labels if this has been requested
c
         ichread=0
         ifound=0
         do while(ichread.lt.nchread.and.ifound.eq.0)
            ichread=ichread+1
            if(alabel(ia).eq.labelold(ichread)) ifound=ichread
         end do
         if(ifound.gt.0) alabel(ia)=labelnew(ifound)
         call ftype(ia)
         if(qtype(ia).ne.1) natomr=natomr+1

      end do
c
c check that potential parameters have been defined for all the atom types
c
      do it=1,ntype
         ipread=0
         ifound=0
         do while(ipread.lt.npread.and.ifound.eq.0)
            ipread=ipread+1
            if(type(it).eq.plabel(ipread)) ifound=ipread
         end do
      write(6,*) type(it),' ',ifound
         if(ifound.eq.0) then
            ierr=1
            do while(ierr.ne.0)
               write(6,115) type(it)
115   format(/'makemole> Type epsilon, sigma, mass'
     *,', charge and chemical symbol '
     *,'for component ',a3//'makemole> ? ',$)
               read(5,*,iostat=ierr) pread(npread+1),sread(npread+1)
     *,mread(npread+1),qread(npread+1),cread(npread+1)
            end do
            npread=npread+1
            plabel(npread)=type(it)
            ifound=npread
         endif
         atsymbol(it)=cread(ifound)
         iexchange(it)=0
         epscor(it)=pread(ifound)
         sigcor(it)=sread(ifound)
         amassa(it)=mread(ifound)
         charge(it)=qread(ifound)

c Now do the same for charge radius, except that if it is not specified
c we set it to zero.

         ifound=0
         iqrread=0
         do while (iqrread.lt.nqrread.and.ifound.eq.0)
            iqrread=iqrread+1
            if(type(it).eq.qrlabel(iqrread)) ifound=iqrread 
         end do
         if(ifound.eq.0) then
            qradius(it)=0.0
         else
            qradius(it)=qrread(ifound)
         endif
      end do
c
c set the box size
c
      write(6,*) natomr,rho,alat,blat,clat,phib,thetac,phic
      if(.not.unitcelldefined) then
         vol=natomr/rho
         dmax=vol**(1.0/3.0)
         alat=dmax
         blat=dmax
         clat=dmax
         phib=90.0
         thetac=0.0
         phic=0.0
         cubic=.true.
      else
         cubic=.false.
         vol=alat*blat*sind(phib)*clat*cosd(thetac)
         dmax=vol**(1.0/3.0)
      endif
      write(6,*) vol
c
c setup molecular size factor
c
      dmaxmol=dmax

c
c write out the revised .mol file
c
      open(10,file=fnameato,status='unknown',iostat=ierr)
      if(filetype(1:ilenfiletype).eq.'.mol') then
c
c create a second .atm file which will show the atom
c numbers for future reference (e.g. to set up rotational groups, etc.)
c
            fnameatm=directory_name(1:nchardir)
     *//readfile(1:ilenreadfile)//'.atm'
c     write(6,*) fnameatm(1:50)

         open(11,file=fnameatm,status='unknown',iostat=ierr)
c
c first number represents the number of characters in a column (including the delimiting space 
c between columns)
c
         write(10,*) nspace
         write(11,*) nspace
c
c maximum number of characters in a line
c
         ncharmax=(ncolumnmax+1)*nspace
         if(ncharmax.lt.80) ncharmax=80
c
c write a blank line - this can be used for a ruler line to show where the columns occur
c
         write(10,100) ruler(1:ncharmax)
         write(11,100) ruler(1:ncharmax)
c
c step through the rows and find the labels and stars that occur here
c
         do ir=1,nrow
c
c the first four characters are reserved for the row number plus a space
c
            write(sometext,119) ir
 119        format(i3,1x)
            ncharmax=0
            is=0
            do while(is.lt.nstar)
             is=is+1
             if(irstar(is).eq.ir) then
                icfirst=(icstar(is)-1)*nspace+2+4
                sometext(icfirst:icfirst)='*'
                if(icfirst.gt.ncharmax) ncharmax=icfirst
             endif
            end do
            il=0
            do while(il.lt.nlabel)
             il=il+1
             if(irlabel(il).eq.ir) then
                icfirst=(iclabel(il)-1)*nspace+2+4
                if(ilabel(il).gt.1) then
                   write(sometext(icfirst-1:icfirst-1),103) ilabel(il)
 103               format(i1)
                endif
                sometext(icfirst:icfirst+2)=label(il)
                if(label(il).ne.'map') then
                   if(icfirst+2.gt.ncharmax) ncharmax=icfirst+2
                else
                   ilenoffset=index(offset(il),' ')-1
                   sometext(icfirst+3:icfirst+2+ilenoffset)=offset(il)(1:ilenoffset)
                   if((icfirst+2+ilenoffset).gt.ncharmax) ncharmax=icfirst+2+ilenoffset
                endif
             endif
            end do
            write(10,100) sometext(1:ncharmax)
c
c write the atom numbers
c
            il=0
            do while(il.lt.nlabel)
             il=il+1
             if(irlabel(il).eq.ir) then
                icfirst=(iclabel(il)-1)*nspace+2+4
                if(ianumber(1,il).lt.10) then
                   write(sometext(icfirst:icfirst+2),108) ianumber(1,il)
 108               format(i1)
                else if(ianumber(1,il).lt.100) then
                   write(sometext(icfirst:icfirst+2),109) ianumber(1,il)
 109               format(i2)
                else
                   write(sometext(icfirst-1:icfirst+1),110) ianumber(1,il)
 110               format(i3)
                endif
                if(icfirst+2.gt.ncharmax) ncharmax=icfirst+2
             endif
            end do
            write(11,100) sometext(1:ncharmax)
         end do
         close(11)

      else if(filetype(1:ilenfiletype).eq.'.gmol') then

         write(10,*) filetype(1:ilenfiletype),iatomlabeladd
      
c       For gmol file, write the atoms

         ia1=0
         do while (ia1.lt.natom)
         
            ia1=ia1+1
            write(10,*) 'atom ',ia1-iatomlabeladd,' ',alabelsave(ia1),' '
     +,xa(ia1),ya(ia1),za(ia1),nbondp(ia1),((iabond(ia1,j)-iatomlabeladd),j=1,nbondp(ia1))

         end do
      endif
c
c now write the bond lengths, angles, extra distances, rotational groups, and potential parameters
c
      ibread=0
      do while(ibread.lt.nbread)
         ibread=ibread+1
         sometext=' '
         write(sometext,104) 'bond',lbond1(ibread),lbond2(ibread),bread(ibread)
104   format(a,2(1x,a3),1x,f10.5)
         ncharmax=4+2*4+1+10
         write(10,100) sometext(1:ncharmax)
      end do
      iaread=0
      do while(iaread.lt.naread)
         iaread=iaread+1
         sometext=' '
         write(sometext,105) 'angle',langle1(iaread),langle2(iaread),langle3(iaread)
     *,aread(iaread)
105   format(a,3(1x,a3),1x,f10.5)
         ncharmax=5+3*4+1+10
         write(10,100) sometext(1:ncharmax)
      end do
      id=0
      do while(id.lt.ndread)
         id=id+1
         sometext=' '
         write(sometext,121) 'dihedral'
     +,(idread(j,id)-iatomlabeladd,j=1,4),dread(id)
121   format(a,4(1x,i3),1x,f10.5)
         ncharmax=8+4*4+1+10
         write(10,100) sometext(1:ncharmax)
      end do
      irread=0
      do while(irread.lt.nrread)
         irread=irread+1
         sometext=' '
         write(sometext,107) 'rot',iarot1(irread)-iatomlabeladd
     +,iarot2(irread)-iatomlabeladd
107   format(a,2(1x,i3))
         ncharmax=3+2*4
         write(10,100) sometext(1:ncharmax)
      end do
      irread=0
      do while(irread.lt.nwrread)
         irread=irread+1
         sometext=' '
         write(sometext,107) 'wholerot',iawrot1(irread)-iatomlabeladd
     +,iawrot2(irread)-iatomlabeladd
         ncharmax=8+2*4
         write(10,100) sometext(1:ncharmax)
      end do
      irread=0
      do while(irread.lt.nirread)
         irread=irread+1
         sometext=' '
         write(sometext,107) 'introt',iairot1(irread)-iatomlabeladd
     +,iairot2(irread)-iatomlabeladd
         ncharmax=6+2*4
         write(10,100) sometext(1:ncharmax)
      end do
      ieread=0
      do while(ieread.lt.neread)
         ieread=ieread+1
         sometext=' '
         write(sometext,106) 'extra',iaextra1(ieread)-iatomlabeladd
     +,iaextra2(ieread)-iatomlabeladd,eread(ieread)
106   format(a,2(1x,i3),1x,f10.5)
         ncharmax=5+2*4+1+10
         write(10,100) sometext(1:ncharmax)
      end do
c
c write the change labels
c
      ichread=0
      do while(ichread.lt.nchread)
         ichread=ichread+1
         sometext=' '
         write(sometext,118) 'changelabel',labelold(ichread),labelnew(ichread)
118   format(a,2(1x,a3))
         ncharmax=11+8
         write(10,100) sometext(1:ncharmax)
      end do
      ipread=0
      do while(ipread.lt.npread)
         ipread=ipread+1
         sometext=' '
         write(sometext,113) 'potential',plabel(ipread)
     *,pread(ipread),sread(ipread)
     *,mread(ipread),qread(ipread),cread(ipread)
113   format(a,1x,a3,4(1x,e12.5),1x,a3)
         ncharmax=9+4+5*11
         write(10,100) sometext(1:ncharmax)
      end do
      iqrread=0
      do while(iqrread.lt.nqrread)
         iqrread=iqrread+1
         sometext=' '
         write(sometext,113) 'qradius',qrlabel(iqrread)
     *,qrread(iqrread)
         ncharmax=8+4+2*11
         write(10,100) sometext(1:ncharmax)
      end do
c
c write temperature and density
c
      write(sometext,114) 'temperature',simtk
114   format(a,1x,e13.6)
      ncharmax=11+14
      write(10,100) sometext(1:ncharmax)
      write(sometext,114) 'vibtemp',tempb
      ncharmax=7+14
      write(10,100) sometext(1:ncharmax)
      write(sometext,114) 'angtemp',tempangreq
      ncharmax=7+14
      write(10,100) sometext(1:ncharmax)
      write(sometext,114) 'dihtemp',tempdihreq
      ncharmax=7+14
      write(10,100) sometext(1:ncharmax)
c
c write ecore and dcore
c
      write(sometext,120) 'ecoredcore',ecore,dcore
120   format(a,1x,f10.5,1x,f10.5)
      ncharmax=10+22
      write(10,100) sometext(1:ncharmax)
      if(unitcelldefined) then
         write(sometext,122) 'unitcell',alat,blat,clat,phib,thetac,phic
122      format(a,6(1x,f10.5))
         ncharmax=8+66
         write(10,100) sometext(1:ncharmax)
      else
         write(sometext,114) 'density',rho
         ncharmax=7+14
         write(10,100) sometext(1:ncharmax)
      endif
      close(10)

c
c write the .ato file
c
      fnameato=directory_name(1:nchardir)
     *//readfile(1:ilenreadfile)//'.ato'
      ilenfnameato=nchardir+ilenreadfile+4

c
c define the number of molecule types and their template files
c
      nmoltype=1
      moltemplate(nmoltype)=readfile(1:ilenreadfile)
            write(6,*) moltemplate(nmoltype)(1:ilenreadfile)
      write(6,*) fnameato(1:ilenfnameato)
      interntrans(1)=0.0
      internrot(1)=0.0
      wholerot(1)=0.0
      wholetrans(1)=0.0

      call moloutc(fnameato)

      return
      end   

      subroutine fmolec1(ierr,num_threads_input,ftest,ftestb,ftethered)
      
      use cellneighbourlist
      use ergpairreference
      use omp_lib
      use molecom_mod
      use threebody_mod
      use containerforces_mod
      
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'epsrinput_par.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'

      logical ifoundb
      integer im,jm,inn,neighbsum,ic,iref,ierr,num_threads,num_threads_input
      real energy,ftest,ftestb,ftethered
      real sumevib,sumeang,sumedih,ratioang,ratiodih,ratiomin,ratiomax,sumetriplets

!Set the number of threads for parallel sections.
!If num_threads is zero, we determine the maximum possible number of threads and use half of this value

      if(num_threads_input.le.0) then
         num_threads=omp_get_max_threads()/2
      else
         num_threads=num_threads_input
      end if
      if(num_threads.gt.omp_get_max_threads()) num_threads=omp_get_max_threads()      
      call OMP_SET_NUM_THREADS(num_threads)
      write(6,*) 'fmolec1> Number of parallel threads to be used: ',num_threads,omp_get_num_threads(),omp_get_max_threads()
      
      ierr=0

      call get_date_and_time('Start initialise_ergpair#',0)
      call initialise_ergpair(nmol)
      call get_date_and_time('Finish initialise_ergpair#',1)

!     Setup the triplets of atoms if three-body terms are requested
      write(6,'(a,i3,1x,i3)') 'fmolec1> No. of three-body terms = ',nthreeb,nthreebangsmax
      call find_triplets(nmol,natom,xadd,yadd,zadd,xa,ya,za,mref,itype,ntype,type)
!      write(6,*) 'nnbrefmax = ',nnbrefmax
c
c step through atom pairs and calculate separation.
c If within the specified distance range (rmaxpt) then add to the neighbour list
c and calculate inter- or intra-molecular energy as appropriate
c
      ftest=0.0
      ftestb=0.0
! Make sure the intra-molecular energy weights are set
      tempangwt=tempangreq
      tempdihwt=tempdihreq
      tempang=tempangwt
      tempdih=tempdihwt*tempb
      call get_date_and_time('fmolec1> Start energy calculation#',0)
!Do the intramolecular terms first so that they go in the first element of each energy pair list
      do im=1,nmol
         !Get the intramolecular energy of this molecule
         energy=en_mole_intra(im,sumevib,sumeang,sumedih)
         ftestb=ftestb+energy
         jm=im
         call set_ergpair(im,jm,energy)
!         write(6,*) im,nergpairm(im)
!         if(tethered(im)) call check_tethered(im,tetheredatom(im)) ! Set each molecule within specified range of its tethering position if needed.
      end do
      sumetriplets=0.0
      if(nthreeb.gt.0) then
         do id=1,nthreeb
            sume_triplet_dist(id)=0.0
         end do
      end if
      do ic=1,ncell
         if(nmolsincell(ic).gt.0) then
            call make_neighbour_list_from_cells(ic,nnbref,nbref)
            do iref=1,nmolsincell(ic)
               im=molsincell(iref,ic)
!$omp parallel do private(inn,jm) shared(ergold)
               do inn=1,nnbref
                  jm=nbref(inn)
                  ergold(inn)=0.0
                  !Avoid double counting pairs of molecules
                  if(jm.gt.im) then
!      write(6,*) 'fmolec1> ',im,naf(im),nato(im),jm,naf(jm),nato(jm)
                     ergold(inn)=en_mole_pair(im,jm)
                  end if
               end do
!$omp end parallel do
               inn=1
               do while(inn.lt.nnbref)
                  inn=inn+1
                  jm=nbref(inn)
                  if(ergold(inn).ne.0.0.and.jm.gt.im) then
                     call set_ergpair(im,jm,ergold(inn))
                     ftest=ftest+ergold(inn)
                  end if
               end do
! Get the triplet energy for this molecule
               if(nthreeb.gt.0) then
                  ergold_triplet(im)=en_mole_triplets(im,xadd,yadd,zadd,phix,phiy,phiz,xa,ya,za,mref,tethered)
!                  write(6,*) im,ergold_triplet(im),(id,en_mole_triplet_dist(id),id=1,nthreeb)
                  sumetriplets=sumetriplets+ergold_triplet(im)
                  do id=1,nthreeb
                     sume_triplet_dist(id)=sume_triplet_dist(id)+en_mole_triplet_dist(id)
                  end do
!                  write(6,*) im,(id,sume_triplet_dist(id),id=1,nthreeb)
                  ftest=ftest+ergold_triplet(im)
               end if
               if(numcontainerforces.gt.0) then
                  ergold_container(im)=en_containerforce(im,naf(im),nal(im),xadd,yadd,zadd,xa,ya,za
     1,itype,mref,npairref,irefpcof,rmaxex,pcoefexmax,expecf,expecfpwr,dmax,harmonic)
!               write(6,*) im,ergold_container(im)
                  ftest=ftest+ergold_container(im)
               end if
            end do
         endif
      end do
      write(6,*) ftest/real(nmol),ftestb/real(nmol)
!For tethered molecules, need to get their tethering energies
      ftethered=0.0
      do im=1,nmol
         if(tethered(im)) then
            ftethered=ftethered+en_tethered(im,tetheredatom(im))
         end if
      end do
      ftethered=ftethered/real(nmol)
      call get_date_and_time('fmolec1> Finish energy calculation#',1)

c calculate the average number of neighbours per molecule

      call get_neighbsum(neighbsum)
      write(6,101) neighbsum
101   format(/'fmolec1> Average no. of neighbours per molecule = ',i4)
      write(6,102) sumetriplets/real(nmol)
102   format('fmolec1> Average triplet energy per molecule = ',e12.5)
      if(nthreeb.gt.0) write(6,*) 'fmolec1> Total energy for each triplet: ',(sume_triplet_dist(id)/real(nmol),id=1,nthreeb)
      write(6,103) ftethered
103   format('fmolec1> Average tethered energy per molecule = ',e12.5/)

c normalise energies per molecule

      ftest=ftest/nmol
      ftestb=ftestb/nmol

      return
      end

      subroutine get_direct_bonds()

      !Gets the direct bonds between atoms (as opposed to those involving an intermediary such as a bond angle).
      !The assumption is that for a triangle of three atoms, the longest bond of the three is not a direct bond
      !Hence this method would not work if the bondangle were less than ~60 degrees (which is unlikely in practice).

      use molecom_mod
      include 'epsr_dimension.inc'
      include 'molecom.inc'

      integer ia,iaref,iaref1,ib1,ib2,ib,ic,iabondsave
      real tbondsave,tab,tac,tbc
      logical found,founddirect

      !First step is to sort bonds for each atom in order of increasing bond length.
      !This is because to use get_rotgrp we require only the number of direct bonds for each atom - their atom numbers
      !are assumed to be at the beginning of the list, as generated with makemole

      do ia=1,natom
         ib1=0
         ib2=ib1+1
         do while(ib2.lt.nbond(ia))
            ib1=ib1+1
            ib2=ib1+1
            if(abs(tbond(ia,ib2)).lt.abs(tbond(ia,ib1))) then
               iabondsave=iabond(ia,ib2)
               tbondsave=tbond(ia,ib2)
               iabond(ia,ib2)=iabond(ia,ib1)
               tbond(ia,ib2)=tbond(ia,ib1)
               iabond(ia,ib1)=iabondsave
               tbond(ia,ib1)=tbondsave
               !Start again to be sure
               ib1=0
               ib2=ib1+1
            end if
         end do
      end do
                 
      do ia=1,natom
         nbondp(ia)=0
         irefb=0
         do while (irefb.lt.nbond(ia))
            irefb=irefb+1
            !Only direct bonds are positive
            if(tbond(ia,irefb).gt.0.0) nbondp(ia)=nbondp(ia)+1
         end do
!         write(6,*) ia,nbond(ia),(iabond(ia,j),tbond(ia,j),j=1,nbond(ia))
!         write(6,*) ia,nbondp(ia),(iabond(ia,j),tbond(ia,j),j=1,nbondp(ia))
      end do

      return
      end

      subroutine get_rotgrp(ir1,ir2,labrott,nrott)

      use reallocation_routines
      use molecom_mod

      include 'epsr_dimension.inc'
      include 'molecom.inc'

      integer ir1,ir2,nrott,nrottmax
            integer, dimension(:), allocatable, intent(out)           :: labrott !Temporary rotational group


c For each of these atoms check how many atoms would be involved in the rotational
c group. Choose the one with the smallest number of rotatable atoms

      nrott=0
            nrottmax=0
      iselect=ir1
      nrotnew=0
      do while(nrotnew.le.nrott)
c
c step through the primary bonded atoms associated with the current atom iselect.
c If it is not one of the axis atoms and does not already exist in the rotational
c group add it to the group.
c
         iref=0
         do while(iref.lt.nbondp(iselect))
            iref=iref+1
            iaref=iabond(iselect,iref)
c
c check that this atom is not one of the axis atoms
c
            if(iaref.ne.ir1.and.iaref.ne.ir2) then
c 
c now check that it does not already exist within the rotational group
c
               itest=0
               ifound=0
               do while(itest.lt.nrott.and.ifound.eq.0)
                  itest=itest+1
                  if(iaref.eq.labrott(itest)) ifound=iaref
               end do
               if(ifound.eq.0) then
c
c add this atom to the group
c
                  nrott=nrott+1
                              if(nrott.gt.nrottmax) then
                                 nrottmax=nrottmax+10
                                 call reallocate1d_i(labrott,nrottmax)
                              end if
                  labrott(nrott)=iaref
               endif
            endif
!       write(6,*) iselect,nbondp(iselect),iref,iaref,nrott
         end do
         nrotnew=nrotnew+1
         if(nrotnew.le.nrott) then
            iselect=labrott(nrotnew)
         endif
      end do

      return
      end

      subroutine get_whole_rotgrp(im,ir1,ir2,labrott,nrott)

      use reallocation_routines
      use molecom_mod

! Makes a list of all the atoms in the molecule, except the two that form
! the axis

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      integer ir1,ir2,nrott,im,nrottmax
      integer, dimension(:), allocatable, intent(out)        :: labrott


!      write(6,*) ir1,ir2
      nrott=0
      nrottmax=0
      do ia=naf(im),nal(im)
         if(ia.ne.ir1.and.ia.ne.ir2) then
c
c add this atom to the group
c
            nrott=nrott+1
            if(nrott.gt.nrottmax) then
               nrottmax=nrottmax+10
               call reallocate1d_i(labrott,nrottmax)
            end if
            labrott(nrott)=ia
         endif
      end do

      return
      end

      subroutine rotatecoords(czphi,szphi,cztheta,sztheta,czchi,szchi,x,y,z)

c Rotates the coordinates x,y,z by the Euler angles zphi,ztheta,zchi. 
c The rotation zphi about the z-axis is performed first, followed by a rotation ztheta about the new y axis, followed
c by a further rotation zchi about the z-axis.

c Rotated coordinates are put back in x,y,z

      real x,y,z,xdum,ydum,zdum
      real czphi,szphi,cztheta,sztheta,czchi,szchi

c Perform the initial rotation about the z-axis

      xdum=x*czphi+y*szphi
      ydum=y*czphi-x*szphi
      zdum=z

c Perform the rotation ztheta about the y-axis 

      x=xdum*cztheta-zdum*sztheta
      y=ydum
      z=zdum*cztheta+xdum*sztheta

c Perform final the rotation zchi about the z-axis

      xdum=x*czchi+y*szchi
      ydum=y*czchi-x*szchi
      zdum=z

      x=xdum
      y=ydum
      z=zdum

      return
      end     

      function dihedral(im1,lb1,lb2,ld1,ld2)

      use molecom_mod
      !Calculates the dihedral angle about the axis defined by lb1 and lb2 (in the same molecule)
      !by atoms ld1 (bonded to lb1) and ld2 (bonded to lb2)). It is assumed atom numbers have already
      !been correctly assigned to be within the same molecule, so no testing is done to prove this.

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
             
      real dihedral,pi,x1,y1,z1,x2,y2,z2,scapro
      real cd(3),co1(3),co2(3),cv1(3),cv2(3),dmag,cpar(3)
      integer im1,lb1,lb2,ld1,ld2
      integer ia1,ia2

      pi=4.0*atan(1.0)
      ia1=naf(im1)+lb1-1
      ia2=naf(im1)+lb2-1
      x1=xa(ia1)
      y1=ya(ia1)
      z1=za(ia1)
      x2=xa(ia2)
      y2=ya(ia2)
      z2=za(ia2)
c
c define the separation vector cd between the two base atoms and normalize it to 
c unit value
c
      cd(1)=(x2-x1)
      cd(2)=(y2-y1)
      cd(3)=(z2-z1)
      dmag=sqrt(scapro(cd,cd))
      cd(1)=cd(1)/dmag
      cd(2)=cd(2)/dmag
      cd(3)=cd(3)/dmag
c
c now determine the orientation vectors co for base atom 1 and base atom 2, using 
c the values of ltyp2 and ltyp3 
c
      io1=naf(im1)+ld1-1
      io2=naf(im1)+ld2-1
      co1(1)=xa(io1)-x1
      co1(2)=ya(io1)-y1
      co1(3)=za(io1)-z1
      co2(1)=xa(io2)-x2
      co2(2)=ya(io2)-y2
      co2(3)=za(io2)-z2
c
c define new vectors normal to both cd and co
c
      call vecpro(co1,cd,cv1)
      call vecpro(co2,cd,cv2)
c
c renormalize these to unit vectors
c
      dmag=sqrt(scapro(cv1,cv1))
      cv1(1)=cv1(1)/dmag
      cv1(2)=cv1(2)/dmag
      cv1(3)=cv1(3)/dmag
      dmag=sqrt(scapro(cv2,cv2))
      cv2(1)=cv2(1)/dmag
      cv2(2)=cv2(2)/dmag
      cv2(3)=cv2(3)/dmag
c calculate the cosine of the angle between cv1 and cv2
      cosphi=scapro(cv1,cv2)
      if(cosphi.gt.1.0) cosphi=1.0
      if(cosphi.lt.-1.0) cosphi=-1.0
c
c calculate the vector product of cv1 and cv2 - this should be either parallel
c (phi positive) or antiparallel to cd (phi negative)
c
      call vecpro(cv1,cv2,cpar)
      sinphi=scapro(cpar,cd)
      dihedral=atan2(sinphi,cosphi)*180/pi

      return
      end
      
      subroutine generate_mp()

      use molecom_mod
      use containerforces_mod
      
      !Generates a WCA repulsive potential to prevent molecular overlaps

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'epsrinput_par.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'

      pi=4.0*atan(1.0)
      if(power.gt.0.0) ratt=1.0e10**(1.0/power)
c
c maximum value of a sum - this it to prevent overflows
c
      summax=exp(25.0)
      nterms=ncoeffp
      decayp=psigma2
c 
c set up WCA between pairs.
c
c first define prefactor for the well depth
c
      rat=6.0/power
      rat1=sqrt(rat)
      rat2=1.0/rat1
      power1=(power+6)/(power-6)
      rat3=rat1**power1
      prefac=(rat2-rat1)*rat3
      prefac=1.0/prefac
      !Determine minimum of LJ potential - this is maximum non-zero value of WCA
      rat4=1.0/(power-6.0)
      rmaxwca=(dcore/sizefactor)*(power/6.0)**rat4
      do ic=1,npairs
            ect=ecore
            dc=dcore/sizefactor
            x1=x(1)
            x2=x(2)
            do j=2,npt
                  x2=x(j)     
c
c modified L-J potential n-6 where n>6
c
c
c modified L-J potential n-6 where n>6
c
                  ddif=dc/x2
                  if(ddif.gt.0.0.and.ect.gt.0.0.and.x2.lt.rmaxwca) then
                    fac=prefac*ect
                    xx=alog(ddif)
                    pwr1=0.5*(power+6.0)*xx
                    if(abs(pwr1).gt.20.0) pwr1=sign(20.0,pwr1)
                    pwr2=0.5*(power-6.0)*xx
                    if(abs(pwr2).gt.40.0) pwr2=sign(40.0,pwr2)
                    sumlj=2.0*fac*sinh(pwr2)*exp(pwr1)+ect
                  else
                    sumlj=0.0
                  endif
!                  ddif=(dc-x2)
!                  if(ddif.gt.25.0) then
!                    sumlj=ect*summax
!                  else if(ddif.lt.-25.0) then
!                    sumlj=0.0
!                  else
!                    sumlj=ect*exp(ddif)
!                  endif
                  xav(j,ic)=sumlj*wfac(j)
                  x1=x2
            end do
            xav(1,ic)=xav(2,ic)
c
c set up gradients of these potentials
c
            sumlj1=xav(1,ic)
            x1=xstar(1)
            x2=xstar(2)
            do j=2,npt
                  x2=xstar(j) 
                  xdif=x2-x1
                  sumlj2=xav(j,ic)
                  xav1(j-1,ic)=(sumlj2-sumlj1)/xdif
                  x1=x2
                  sumlj1=sumlj2
            end do
            xav1(npt,ic)=0.0
            irefp=irefpcof(ic)
            !Get a list of atom pairs for this type. If rminex(irefp).ge.0.0 then the list is empty.
            call get_containerforce_list(nmol,ic,rmaxex(irefp),it,jt,itype,qtype,natom,qtypea,ntype)
      end do
!      open(20,file='junk.pot',status='unknown')
!      do j=1,npt
!            write(20,*) x(j),xav(j,1),xav1(j,1)
!      end do
!      close(20)

      return
      end
      
      subroutine crystfrac(x,y,z,fa,fb,fc)
      include 'lattice.inc'
      real, intent(in) ::  x,y,z
      real, intent(out) :: fa,fb,fc
! Gets the fractional coordinates of x,y,z in terms of crystallographic axes a,b,c
      if(cubic) then
         fc=z/clatz
         fb=y/blaty
         fa=x/alatx
! Apply minimum image convention
         call minimage(fa,fb,fc)
      else
         fc=z/clatz
         fb=(y-fc*claty)/blaty
         fa=(x-fb*blatx-fc*clatx)/alatx
! Apply minimum image convention
         call minimage(fa,fb,fc)
      endif
      return
      end

      function crystarg(fa,fb,fc,h,k,l)
      include 'lattice.inc'
      real, intent(in) ::fa,fb,fc,h,k,l
      real :: crystarg
      crystarg=twopi*(fa*h+fb*k+fc*l)
      return
      end

      subroutine bragginten(qstep,qmax,qwidthin,fwhmq,multa,multb,multc
     *,hklqmin,nbragg,qbragg,intenbragg)

      use omp_lib
      use molecom_mod
      use bragginten_mod
      use reciprocallattice_mod
      use diffuse_mod
      
!Calculates Bragg intensities up to the specified Q max assuming a unit cell which a/multa, b/multb and c/multc      

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'

      integer :: multa,multb,multc,nbragg,nunit !Number of files, of unit cells along each direction a, b, c, and number of Q values.
      real :: qstep,qmax,qwidth,qwidthin,fwhmq,hklqmin,qwidthfac
      real, dimension(nbragg) :: qbragg       !Q scale for output
      real, dimension(npairs,nbragg) :: intenbragg       !Intensity for each file as a function of Q

      character(len=256) :: fname,outputfile,formatout
      integer :: ncharoutput,iref,it,jt,ib,im,ia,hi,ki,li,nline,nlinemax,hmax,kmin,kmax,lmin,lmax
     *,nfirst,nlast,nhklmax,nsums

      real :: pi,vol,xuse,yuse,zuse,huse,kuse,luse,arg,crystarg !Dummy variables
      real :: sqrmoda,sqrmodb,sqrmodc,ramptot,iamptot,qbraggtemp,intentemp !Square modulus of each reciprocal lattice vector
      real :: gfactor,gfactor0,sum,diff,add,lfactor,qhkltest,cprod,xbragg,ybragg,fc1,fc2,ramptemp,iamptemp
      real, dimension(:), allocatable :: fa,fb,fc!,ramptemp,iamptemp
      integer, dimension(:,:), allocatable :: iatyperef
      integer, dimension(:), allocatable :: natyperef
      logical :: found

      write(6,*) 'bragginten> num_threads_to be used: ',omp_get_max_threads()
      pi=4.0*atan(1.0)
      twopi=2.0*pi
      !Output file name
      outputfile=fname_run(1:ilenfname_run)
      ncharoutput=ilenfname_run
!      outputfile='bragginten'
!      ncharoutput=index(outputfile,' ')-1
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
      qwidth=abs(qwidthin)+qmax*fwhmq
      qwidthfac=5.0
      hmax=nint((qmax+qwidthfac*qwidth)/abs(recipra(1)))+1
      kmax=nint((qmax+qwidthfac*qwidth)/abs(reciprb(2)))+1
      lmax=nint((qmax+qwidthfac*qwidth)/abs(reciprc(3)))+1
      !Maximum number of hkl values
      nhklmax=(hmax+1)*(2*kmax+1)*(2*lmax+1)
      !Allocate arrays
      if(allocated(hihkl)) deallocate(hihkl,kihkl,lihkl,qhkl)
      allocate(hihkl(nhklmax),kihkl(nhklmax),lihkl(nhklmax),qhkl(nhklmax))
!Now run over the allowed h,k,l values and calculate scattering amplitudes and intensities
      ib=0
      hi=0 !Only calculate positive values of hi and use inversion symmetry to get negative hi
      call get_date_and_time('Start setting up hkl values and fractional coordinates#',0)
      do while (hi.le.hmax)
         kmin=-nint((qmax+qwidthfac*qwidth+real(hi)*recipra(2))/abs(reciprb(2)))-1
         kmax=nint((qmax+qwidthfac*qwidth-real(hi)*recipra(2))/abs(reciprb(2)))+1
         ki=kmin
         do while (ki.le.kmax)
            lmin=-nint((qmax+qwidthfac*qwidth+real(hi)*recipra(3)+real(ki)*reciprb(3))/abs(reciprc(3)))-1
            lmax=nint((qmax+qwidthfac*qwidth-real(hi)*recipra(3)-real(ki)*reciprb(3))/abs(reciprc(3)))+1
            li=lmin
            do while (li.le.lmax)
               do it=1,3
                  rbragg(it)=real(hi)*recipra(it)+real(ki)*reciprb(it)+real(li)*reciprc(it)
               end do
!         write(10,*) (rbragg(it),it=1,3)
               qhkltest=sqrt(scapro(rbragg,rbragg))
               if(qhkltest.le.(qmax+qwidthfac*qwidth)) then
                  ib=ib+1
                  hihkl(ib)=hi
                  kihkl(ib)=ki
                  lihkl(ib)=li
!            write(6,*) ib,hi,ki,li
                  qhkl(ib)=qhkltest
               end if
               li=li+1
            end do
            ki=ki+1
         end do
         hi=hi+1
      end do
      nhkl=ib
      if(allocated(fa)) deallocate(fa,fb,fc,ramp,iamp,intenhkl,qunique,nqunique
     *,intenqunique,ibvalues,iqvalues,natyperef,iatyperef)
      allocate(fa(natom),fb(natom),fc(natom),ramp(nhkl,ntype),iamp(nhkl,ntype)
     *,intenhkl(nhkl,npairs),iqvalues(nhkl),qunique(nhkl),nqunique(nhkl),natyperef(ntype)
     *,iatyperef(natom,ntype))
! Set up the unique q values
      ntotalqunique=0
      nquniquemax=0
      do ib=1,nhkl
         qhkltest=qhkl(ib)
         iq=0
         found=.false.
         do while(iq.lt.ntotalqunique.and..not.found)
            iq=iq+1
            found=abs(qhkltest-qunique(iq)).lt.1.0e-5
         end do
         if(.not.found) then
            ntotalqunique=ntotalqunique+1
            qunique(ntotalqunique)=qhkltest
            nqunique(ntotalqunique)=1
            if(nqunique(ntotalqunique).gt.nquniquemax) nquniquemax=nqunique(ntotalqunique)
         else
            nqunique(iq)=nqunique(iq)+1
            if(nqunique(iq).gt.nquniquemax) nquniquemax=nqunique(iq)
         end if
      end do
      allocate(ibvalues(nquniquemax,ntotalqunique))
! Now assign the ib values for each unique Q value.
! Set counters to zero
      do iq=1,ntotalqunique
         nqunique(iq)=0
      end do
! Set up the revised Q scale and intensities
! Store the ib reference numbers for each unique Bragg peak            
      do ib=1,nhkl
         qhkltest=qhkl(ib)
         iq=0
         found=.false.
         do while(iq.lt.ntotalqunique.and..not.found)
            iq=iq+1
            found=abs(qhkltest-qunique(iq)).lt.1.0e-5
         end do
         if(found) then
            nqunique(iq)=nqunique(iq)+1
            ibvalues(nqunique(iq),iq)=ib
            iqvalues(ib)=iq
         end if
      end do
! Set up fractional coordinates
      do ia=1,natom
         im=mref(ia)
         xuse=xa(ia)+xadd(im)
         yuse=ya(ia)+yadd(im)
         zuse=za(ia)+zadd(im)
         call crystfrac(xuse,yuse,zuse,fa(ia),fb(ia),fc(ia))
      end do
      call get_date_and_time('Finished setting up hkl values and fractional coordinates#',1)
      call get_date_and_time('Start calculating Bragg intensities#',0)
      ramp=0.0
      iamp=0.0
      natyperef=0
      do ia=1,natom
         it=itype(ia)
         natyperef(it)=natyperef(it)+1
         iatyperef(natyperef(it),it)=ia
      end do
      do it=1,ntype
!$omp parallel do private(ir,ia,ib,hi,ki,li,huse,kuse,luse,arg,ramptemp,iamptemp) shared(ramp,iamp)
         do ir=1,natyperef(it)
            ia=iatyperef(ir,it)
            do ib=1,nhkl
               hi=hihkl(ib)
               ki=kihkl(ib)
               li=lihkl(ib)
               huse=real(hi*multa)
               kuse=real(ki*multb)
               luse=real(li*multc)
               arg=twopi*(fa(ia)*huse+fb(ia)*kuse+fc(ia)*luse)
               ramptemp=cos(arg)
               iamptemp=-sin(arg)
               ramp(ib,it)=ramp(ib,it)+ramptemp
               iamp(ib,it)=iamp(ib,it)+iamptemp
            end do
         end do
!$omp end parallel do
      end do               
! Bragg intensity for each pair of atom types, averaged for the number of unit cells,
! and per atom
      arg=1.0/real(natomr)
! Sum the real and imaginary parts
!$omp parallel do private(it,jt,ip) shared(intenhkl)
      do ip=1,npairs
         it=n1extr(ip)
         jt=n2extr(ip)
         do ib=1,nhkl
            intenhkl(ib,ip)=ramp(ib,it)*ramp(ib,jt)+iamp(ib,it)*iamp(ib,jt)
            intenhkl(ib,ip)=intenhkl(ib,ip)*arg
            if(it.ne.jt) then
               intenhkl(ib,ip)=intenhkl(ib,ip)*2.0
            end if
         end do
      end do
!$omp end parallel do
!Calculate the total Bragg intensity at each unique Q value
      allocate(intenqunique(npairs,ntotalqunique))
      intenqunique=0.0
!$omp parallel do private(ip,iq,iref) shared(intenhkl,intenqunique)
      do ip=1,npairs
         do iq=1,ntotalqunique
            do iref=1,nqunique(iq)
               ib=ibvalues(iref,iq)
               if(hihkl(ib).gt.0) then
                  intenqunique(ip,iq)=intenqunique(ip,iq)+2.0*intenhkl(ib,ip) !Double intensity to capture inversion symmetry
               else
                  intenqunique(ip,iq)=intenqunique(ip,iq)+intenhkl(ib,ip)
               end if
            end do
         end do
      end do
!$omp end parallel do
      call get_date_and_time('Finished calculating Bragg intensities#',1)
!      if(diffusedefined) then
!         nunit=10
!         do id=1,ndiffuse
!            nunit=nunit+1
!            close(nunit)
!         end do
!      end if
!Zero intensities and get scattering levels for each file
      do ip=1,npairs
         do ib=1,nbragg
            intenbragg(ip,ib)=0.0
         end do
      end do
      call get_date_and_time('Start broadening Bragg intensities#',0)
      !Set up the Bragg peaks, combining hkl which have the same Q value
      !Factor outside Gaussian
!      gfactor=volrecipr/(2.0*(twopi**1.5)*abs(qwidth))
!      gfactor0=volrecipr/((twopi**1.5)*abs(qwidth**3))
!      lfactor=volrecipr/((pi**2)*abs(qwidth**3))
      gfactor=volrecipr/(2.0*(twopi**1.5))
      gfactor0=volrecipr/((twopi**1.5))
      lfactor=volrecipr/((pi**2))
      !Step through hkl values and add broadened intensities to structure factor
      do iq=1,ntotalqunique
         qhkltemp=qunique(iq)
!Q dependent peak width
         qwidth=abs(qwidthin)+qhkltemp*fwhmq
         nfirst=max(1,nint((qhkltemp-15.0*qwidth)/qstep)-1)
         nlast=min(nbragg,nint((qhkltemp+15.0*qwidth)/qstep)+1)
         if(qhkltemp.lt.hklqmin) then
            nfirst=1
            nlast=0
         end if
!         write(6,*) ib,qhkl(ib),qwidth,qstep,nfirst,nlast
         it=nfirst-1
         do while(it.lt.nlast)
            it=it+1
            sum=0.0
            if(qwidthin.gt.0.0) then
!Gaussian peak shape
               if(qhkltemp.gt.0.0.and.qbragg(it).gt.0.0) then
                  diff=(qhkltemp-qbragg(it))/qwidth
                  add=(qhkltemp+qbragg(it))/qwidth
                  if(abs(diff).lt.8.0) sum=exp(-(0.5*diff*diff))
                  if(abs(add).lt.8.0) sum=sum-exp(-(0.5*add*add))
                  sum=sum*gfactor/(qhkltemp*qbragg(it)*qwidth)
               else
                  add=qhkltemp/qwidth
                  diff=qbragg(it)/qwidth
                  sum=0.5*(diff*diff+add*add)
                  if(sum.lt.32.0) then
                     sum=gfactor0*exp(-sum)/(qwidth**3)
                  else
                     sum=0.0
                  end if
               end if
            else
!Lorentzian peak shape
               diff=(qhkltemp-qbragg(it))/qwidth
               add=(qhkltemp+qbragg(it))/qwidth
               sum=lfactor/((1.0+diff*diff)*(1.0+add*add)*(qwidth**3))
            endif
            do ip=1,npairs
               intenbragg(ip,it)=intenbragg(ip,it)+intenqunique(ip,iq)*sum
            end do               
         end do
      end do
      call get_date_and_time('Finished broadening Bragg intensities#',1)
!Subtract the single atom scattering and renormalise the scattering to per atom
      ip=0
      do it=1,ntype
         do jt=it,ntype
            ip=ip+1
            if(it.eq.jt) then
               cprod=cfrac(it)*cfrac(jt)
               do ib=1,nbragg
                  intenbragg(ip,ib)=(intenbragg(ip,ib)-cfrac(it))/cprod 
               end do
            else
               cprod=2.0*cfrac(it)*cfrac(jt)
               do ib=1,nbragg
                  intenbragg(ip,ib)=intenbragg(ip,ib)/cprod 
               end do
            endif
         end do
      end do
c Create a filename for the output
      fname=directory_name(1:nchardir)//outputfile(1:ncharoutput)//'.bragg'
      open(10,file=fname,status='unknown')
      write(10,*) '# Multiplicities along a, b and c: ', multa,multb,multc
      write(10,101) '# Volume of unit cell: ',vol
      write(10,101) '# Volume of reciprocal lattice: ',volrecipr
      write(10,101) '# Lattice vectors a: ',alatc(1),alatc(2),alatc(3),sqrt(scapro(alatc,alatc))
      write(10,101) '# Lattice vectors b: ',blatc(1),blatc(2),blatc(3),sqrt(scapro(blatc,blatc))
      write(10,101) '# Lattice vectors c: ',clatc(1),clatc(2),clatc(3),sqrt(scapro(clatc,clatc))
      write(10,101) '# Reciprocal lattice vectors a: ',recipra(1),recipra(2),recipra(3),ampra
      write(10,101) '# Reciprocal lattice vectors b: ',reciprb(1),reciprb(2),reciprb(3),amprb
      write(10,101) '# Reciprocal lattice vectors c: ',reciprc(1),reciprc(2),reciprc(3),amprc
      write(10,'(a,i0)') '# Number of Bragg refections calculated: ',nhkl
101   format(1x,a,4(1x,e12.5))
      !Prepare the output format statement
      write(formatout,102) npairs
102   format('(1x,"# ",i0,3(1x,i3),1x,f10.5,',i3.3,'(1x,e12.5))')
      write(6,'(a)') formatout(1:len_trim(formatout))
      write(10,'(1x,"# ",a)') 'Actual Bragg Q values'
!      do ib=1,nhkl
!         write(10,formatout) hihkl(ib),kihkl(ib),lihkl(ib),qhkl(ib),(intenhkl(ip,ib),ip=1,npairs)
!         write(10,formatout) ib,hihkl(ib),kihkl(ib),lihkl(ib),qhkl(ib)
!      end do
      write(10,'(1x,"# ",a)') 'Unique Bragg Q values and multiplicities' 
      do iq=1,ntotalqunique
         write(10,'(1x,"#",i6,1x,f10.5,1x,i0)') iq,qunique(iq),nqunique(iq)
!         write(10,'(4(1x,i7),1x,e12.5)') (ib,hihkl(ibvalues(ib,iq)),kihkl(ibvalues(ib,iq)),lihkl(ibvalues(ib,iq))
!     *,intenhkl(ibvalues(ib,iq),npairs),ib=1,nqunique(iq))
      end do
      write(formatout,103) npairs
103   format('(1x,f10.5,',i3.3,'(1x,e12.5))')
      write(6,*) formatout(1:index(formatout,' '))
!      do ib=1,nbragg
!         write(10,formatout) qbragg(ib),(intenbragg(ip,ib),ip=1,npairs)
!      end do
      close(10)

      if(allocated(ramp)) deallocate(ramp,iamp)
      if(.not.diffusedefined) then
         if(allocated(hihkl)) deallocate(hihkl,kihkl,lihkl,qhkl,intenhkl)
      end if
      if(allocated(fa)) deallocate(fa,fb,fc,qunique,nqunique,intenqunique,ibvalues,iqvalues,natyperef,iatyperef)

      return
      end

      function radd(rsum,rterm)
      real :: radd,rsum,rterm
      radd=rsum+rterm
      return
      end

      subroutine sharm_diffuse(qstep,nq,nsums,nfile)
! Outputs diffuse scattering if requested, adding to existing diffuse scattering if accumulating.

      use bragginten_mod
      use reciprocallattice_mod
      use diffuse_mod
      use lngamma_mod

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      include 'weights.inc'
      include 'X-rayf0.inc'

      real                                         :: qstep,qmaxdiff,fc1,fc2,cz,pi
      integer                                      :: nq,nsums,ncharoutput,hi,ki,li,nline,nlinemax,ip,iq,im,id
      character(len=256)                           :: fname,outputfile,formatout,formattitle
      real, dimension(:), allocatable              :: intenfile,tempdiff
      logical                                      :: startfromscratch
      character(len=2)                             :: ifmt1,ifmt2,ifmt3,ifmt4,il1fmt,in1fmt
      character(len=29), dimension(:), allocatable :: titledat
      character(len=14)                            :: titler
      character(len=1)                             :: ext

      if(.not.diffusedefined.or.nhkl.le.0) then !Get out of here if no diffuse scattering has been defined or can be calculated
         return
      end if
      pi=4.0*atan(1.0)
! Normalise intensities to atomic fractions
      ip=0
      do it=1,ntype
         do jt=it,ntype
            ip=ip+1
            cprod=cfrac(it)*cfrac(jt)
            if(it.ne.jt) cprod=2.0*cprod
            do ib=1,nhkl
               intenhkl(ip,ib)=intenhkl(ip,ib)/cprod 
            end do
         end do
      end do
! Get the number of terms
      call initialise_diffuse_values(nfile)
      if(allocated(intenfile)) deallocate(intenfile)
      allocate(intenfile(nfile))
! Step through the Bragg peaks, picking out those that match the requested diff patterns
      do ib=1,nhkl
         hi=hihkl(ib)
         ki=kihkl(ib)
         li=lihkl(ib)
! Calculate contribution to sharm coefficients, avoiding Bragg peaks
         if(multdiffa*(hi/multdiffa).ne.hi.or.multdiffb*(ki/multdiffb).ne.ki.or.multdiffc*(li/multdiffc).ne.li) then
            iqdiff=nint(qhkl(ib)/qdiffstep)+1
            do it=1,3
               rbragg(it)=real(hi)*recipra(it)+real(ki)*reciprb(it)+real(li)*reciprc(it)
            end do
            xbragg=rbragg(1)
            ybragg=rbragg(2)
            zbragg=rbragg(3)
            if(qhkl(ib).gt.0.0) then
               ctl=zbragg/qhkl(ib)
            else
               ctl=1.0
            end if
            phil=atan2d(ybragg,xbragg)
! Weight the intensities for each data file
            do iff=1,nfile
               ix=xfile(iff)
               intenfile(iff)=0.0
               if(ix.gt.0) then
! If X-ray data, need to apply the Q-dependent weights
                  iq=nint(qhkl(ib)/qstep)+1
                  if(iq.gt.0.and.iq.le.nq) then
                     do ip=1,npairs
                        intenfile(iff)=intenfile(iff)+weightx(iq,ip,ix)*intenhkl(ip,ib)
                     end do
                  end if
               else
                  do ip=1,npairs
                     intenfile(iff)=intenfile(iff)+cwt(iff,ip)*intenhkl(ip,ib)
                  end do
               end if
            end do
!            write(6,*) ib,hi,ki,li,(intenfile(iff),iff=1,nfile),ctl,phil
            do id=1,ndiffuse
               l=l1ex(id)
               m=n1ex(id)
               factor=real(2*l+1)
               dlmcoeff=factor*dlm(l,m,0,ctl)
               dlmcoeffc=cosd(m*phil)*dlmcoeff
               dlmcoeffs=sind(m*phil)*dlmcoeff
! If h > 0, add on the term for which (h,k,l) --> (-h,-k,-l) (Inversion image, theta --> pi-theta, phi --> phi + pi)
               if(hi.gt.0) then
                  dlmcoeff=factor*dlm(l,m,0,-ctl)
                  dlmcoeffc=dlmcoeffc+cosd(m*(phil+180.0))*dlmcoeff
                  dlmcoeffs=dlmcoeffs+sind(m*(phil+180.0))*dlmcoeff
               end if
!       if(l.eq.36.and.m.eq.0) write(6,*) ib,hi,ki,li,l,m,ctl,phil,dlmcoeff,dlmcoeffc,dlmcoeffs,intenfile(1)
               do iff=1,nfile
                  diffcoeffc(id,iqdiff,iff)=diffcoeffc(id,iqdiff,iff)+dlmcoeffc*intenfile(iff)
                  diffcoeffs(id,iqdiff,iff)=diffcoeffs(id,iqdiff,iff)+dlmcoeffs*intenfile(iff)
               end do
            end do
            if(abs(hi).gt.0) then
               ndiffcoeff(iqdiff)=ndiffcoeff(iqdiff)+2 !Account for negative hi
            else
               ndiffcoeff(iqdiff)=ndiffcoeff(iqdiff)+1
            end if
!       write(6,*) ib,hi,ki,li,qhkl(ib),qdiffstep,iqdiff,diffcoeffc(1,iqdiff,1),diffcoeffs(1,iqdiff,1),ndiffcoeff(iqdiff)
         end if
      end do
      nblock=2*ndiffuse+1
! Set up a title for the first line
      write(formattitle,333) '(1h#,1x,a14,',nblock,'(1x,a29))'
333   format(a12,i3.3,a9)
      ifmt1='i1'
      ifmt2='i2'
      ifmt3='i3'
      ifmt4='i4'
      if(allocated(titledat)) deallocate(titledat)
      allocate(titledat(ndiffuse))
      do id=1,ndiffuse
         il1=l1ex(id)
         in1=n1ex(id)
         if(il1.lt.10) then
            il1fmt=ifmt1
         else if(il1.lt.100) then
            il1fmt=ifmt2
         else if(il1.lt.1000) then
            il1fmt=ifmt3
         else if(il1.lt.10000) then
            il1fmt=ifmt4
         endif
         if(in1.lt.10) then
            in1fmt=ifmt1
         else if(in1.lt.100) then
            in1fmt=ifmt2
         else if(in1.lt.1000) then
            in1fmt=ifmt3
         else if(in1.lt.10000) then
            in1fmt=ifmt4
         endif
         formatout='(1h('//il1fmt//',1h,,'//ifmt1//',1h,,'//il1fmt//',1h;,'//in1fmt//',1h,,'//ifmt1//',1h))'
         write(titledat(id),formatout) il1,0,il1,in1,0
      end do
      titler='   Q' ! axis title
! Set up the files from and to which the diffuse scattering will be read/written
      if(allocated(diffcoefftc)) deallocate(diffcoefftc,diffcoeffts,diffcoeffoldc,diffcoeffolds,tempdiff)
      allocate(diffcoefftc(nqdiff,ndiffuse),diffcoeffts(nqdiff,ndiffuse)
     *,diffcoeffoldc(nqdiff,ndiffuse),diffcoeffolds(nqdiff,ndiffuse),tempdiff(nqdiff))
      outputfile=fname_run(1:ilenfname_run)
      ncharoutput=index(outputfile,'.EPSR')-1
      do iff=1,nfile
         write(diffusename,'("d",i2.2)') iff
         fname=directory_name(1:nchardir)//outputfile(1:ncharoutput)//diffusename//'.DIFFUSE.h01'
         do iq=1,nqdiff
            if(ndiffcoeff(iq).gt.0) then !If possible average coefficients to number of contributions
               do id=1,ndiffuse
                  diffcoeffc(id,iq,iff)=diffcoeffc(id,iq,iff)/real(ndiffcoeff(iq))
                  diffcoeffs(id,iq,iff)=diffcoeffs(id,iq,iff)/real(ndiffcoeff(iq))
               end do
            end if
         end do
! Smooth the coefficients prior to generating maps. Don't do this very much - it distorts image drastically
         do id=1,ndiffuse
            do iq=1,nqdiff
              tempdiff(iq)=diffcoeffc(id,iq,iff)
            end do
            call smood(1,nqdiff,1,tempdiff)
            do iq=1,nqdiff
               diffcoefftc(iq,id)=tempdiff(iq)
            end do
            do iq=1,nqdiff
              tempdiff(iq)=diffcoeffs(id,iq,iff)
            end do
            call smood(1,nqdiff,1,tempdiff)
            do iq=1,nqdiff
               diffcoeffts(iq,id)=tempdiff(iq)
            end do
         end do
         ext='h'
         call outgofzs_sharm(fname,ext,nqdiff,ndiffuse,nsums,qdiff,diffcoefftc,diffcoeffts
     *,diffcoeffoldc,diffcoeffolds,titler,titledat,nqdiff,ndiffuse)
      end do
      if(allocated(intenfile)) deallocate(intenfile,tempdiff,diffcoefftc,diffcoeffts,diffcoeffoldc,diffcoeffolds)
      if(allocated(intenhkl)) deallocate(qhkl,intenhkl)
      if(allocated(diffcoeffc)) deallocate(qdiff,diffcoeffc,diffcoeffs,ndiffcoeff)
      return
      end

      function jhalf(x)
      real :: jhalf, x
      if(abs(x).gt.0.0) then
         jhalf=sin(x)/x
      else
         jhalf=1.0
      end if
      return
      end

      function truncconv(q,q0,rm,qstep,pi)
! Calculates the truncation convolution function for peak at q0 when maximum radius is rm
      real :: truncconv,q,q0,rm,pi,factor,arg,arg1,arg2,jhalf,qstep
      
      truncconv=0.0
      if(abs(q).gt.0.0) then
         if(abs(q0).gt.0.0) then
            factor=q0*rm*qstep/(pi*q)
            arg1=(q0-q)*rm
            arg2=(q0+q)*rm
            truncconv=factor*(jhalf(arg1)-jhalf(arg2))
         else
            factor=2*q0*q0*rm*qstep/(pi*q*q)
            arg=q*rm
            truncconv=factor*(jhalf(arg)-cos(arg))
         end if
      else if(abs(q0).gt.0.0) then
         factor=2*rm*qstep/(pi)
         arg=q0*rm
         truncconv=factor*(jhalf(arg)-cos(arg))
      else
!         truncconv=2*(rm*rm)*rm/(pi*3.0)
         truncconv=0.0
      end if
      return
      end

      subroutine change_boxsize(afac,bfac,cfac,newphib,newthetac,newphic)

      use molecom_mod
c
c Multiplies the lattice vectors by the specified factors and relocates molecules relative
c to new axes.
c
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      real*4 afac,bfac,cfac,newphib,newthetac,newphic
      
c Get the volume of the current simulation box

!      volold=alatx*blaty*clatz

!Need to place each molecule in same relative position in the box as at present
!First generate the current fractional coordinates

      do im=1,nmol
         call generate_fractional_coordinates(xadd(im),yadd(im),zadd(im))
         call generate_fractional_coordinates(phix(im),phiy(im),phiz(im))
      end do
            
! set up lattice parameters for applying the minimum image convention

      alat=alat*afac
      blat=blat*bfac
      clat=clat*cfac
      phib=newphib
      thetac=newthetac
      phic=newphic
      alatx=alat
      blatx=blat*cosd(phib)
      blaty=blat*sind(phib)
      clatx=clat*sind(thetac)*cosd(phic)
      claty=clat*sind(thetac)*sind(phic)
      clatz=clat*cosd(thetac)

! Now convert fractional coordinates back to Cartesian coordinates

      do im=1,nmol
         call generate_cartesian_coordinates(xadd(im),yadd(im),zadd(im))
         call generate_cartesian_coordinates(phix(im),phiy(im),phiz(im))
      end do

c Get the volume of this new box

!            vol=volold/(alatx*blaty*clatz)

c Scale the box dimensions to give the save volume as previously

!            dmaxfac=vol**(1.0/3.0)
!            alat=alat*dmaxfac
!            alatx=alatx*dmaxfac
!            blat=blat*dmaxfac
!            blatx=blatx*dmaxfac
!            blaty=blaty*dmaxfac
!            clat=clat*dmaxfac
!            clatx=clatx*dmaxfac
!            claty=claty*dmaxfac
!            clatz=clatz*dmaxfac
c
c define the Cartesian coordinates of the lattice vectors
c
      alatc(1)=alatx
      alatc(2)=0.0
      alatc(3)=0.0
      blatc(1)=blatx
      blatc(2)=blaty
      blatc(3)=0.0
      clatc(1)=clatx
      clatc(2)=claty
      clatc(3)=clatz
      bxax=blatx/alatx
      cxax=clatx/alatx
      cyby=claty/blaty
!
c calculate the perpendicular distance between the faces of the lattice
c
      perpab=perpd(alatc,blatc,clatc)
      perpbc=perpd(blatc,clatc,alatc)
      perpca=perpd(clatc,alatc,blatc)
      write(6,*) perpab,perpbc,perpca

c set up maximum radius to be used for Fourier transforms and to truncate
c the potential. This is defined as the minimum perpendicular distance between
c planes. 
c
c G(r) and potentials will be calculated to rmax, or maximum allowed distance
c based on number of points and step size, whichever is smaller.
c
      dmax=amin1(perpab,perpbc,perpca)
!      write(6,*) nmol,dmax
      cell=dmax
      if(alat.ne.blat.or.blat.ne.clat.or.clat.ne.alat) then
         cubic=.false.
      else if(phib.ne.90.0.or.thetac.ne.0.0.or.phic.ne.0.0) then
         cubic=.false.
      else
         cubic=.true.
      end if
      
      return
      end
      
      subroutine atogen(wav,ngrsamples)
      
      use cellneighbourlist
      use ergpairreference
      use omp_lib
      use molecom_mod
      use containerforces_mod
      
c
c if mode = 1 this calculates only intermolecular contribution, Any other
c value and it calculates the total contribution
c
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'

      logical found
      real wav(2,mpz,mpairs),probaccept,rhoactual,pi,pi2,spi2,rmaxhere,addwav,scalefac
      real, dimension(:), allocatable     :: csampled
      integer, dimension(:), allocatable  :: molsampled
      real, dimension(:), allocatable  :: molfraction
      integer, dimension(:), allocatable  :: iz1values,icrefvalues
      integer ic,ngrsamples,iaref,jaref,iref
      real, dimension(:), allocatable     :: natomawaysum

      pi=4.0*atan(1.0)
      pi2=2.0*pi
      spi2=sqrt(pi2)
c 
c set weights to zero
c
      do iz=1,npt1
         do in=1,npairs
            wav(1,iz,in)=0.0
            wav(2,iz,in)=0.0
         end do
      end do
      if(npairs.gt.0) then
            do ic=1,npairs
                  potex(ic)=0.0
            end do
      endif
      stepcl=cellst/dmax
c
c generate intramolecular spacings
c
      do im=1,nmol
         if(nato(im).gt.1) then
            do ia=naf(im),nal(im)-1
               x1=xa(ia)
               y1=ya(ia)
               z1=za(ia)
               do ja=ia+1,nal(im)
                     it=itype(ia)
                     jt=itype(ja)
                     icref=npairref(it,jt)
                     x2=xa(ja)
                     y2=ya(ja)
                     z2=za(ja)
                     xx=(x1-x2)*sizefactor
                     yy=(y1-y2)*sizefactor
                     zz=(z1-z2)*sizefactor
c
c apply minimum image convention
c
                     call minimagec(xx,yy,zz)
                     dtest=sqrt(xx*xx+yy*yy+zz*zz)
c find the x value this corresponds to
                     iz1=nint(dtest/stepcl)+2
c
c accumulate the partial g(r) for this pair
c
                     if(iz1.gt.1.and.iz1.le.npt1) then
                        wav(1,iz1,icref)=wav(1,iz1,icref)+2.0
                     endif
               end do
            end do
         endif
      end do
c
c generate intermolecular spacings
c
      if(allocated(csampled)) deallocate(csampled,iz1values,icrefvalues,molsampled,molfraction)
      allocate(csampled(ntype),iz1values(natom),icrefvalues(natom),molsampled(nmoltype),molfraction(nmoltype))
! Define probability for using a molecule
      if(ngrsamples.le.0) ngrsamples=1000
      if(ngrsamples.gt.nmol) ngrsamples=nmol
      probaccept=real(ngrsamples)/real(nmol)
      write(6,*) 'atogen> Acceptance probability = ',probaccept
      do it=1,ntype
         csampled(it)=0.0
      end do
! Get the molecular fraction of each molecule type
      do it=1,nmoltype
         molsampled(it)=0
         molfraction(it)=0.0
      end do
      do im=1,nmol
         molfraction(mtype(im))=molfraction(mtype(im))+1.0
      end do
      do it=1,nmoltype
         molfraction(it)=molfraction(it)/real(nmol)
      end do
      im1=1
      im2=nmol
      addwav=1.0
      if(probaccept.ge.1.0) then
         im1=1
         im2=nmol-1
         addwav=2.0
! If probaccept.ge.1.0 this means we only access nmol-1 first molecules. Hence we need to count the last molecule when calculating the number of 
! first molecules sampled.
         do ia=naf(nmol),nal(nmol)
            it=itype(ia)
            csampled(it)=csampled(it)+1.0
         end do
      end if
      do im=im1,im2
! Always select molecules that are weakly represented
         if(ran1().lt.probaccept.or.molsampled(mtype(im)).lt.1
     *.or.molfraction(mtype(im)).lt.probaccept) then
            molsampled(mtype(im))=molsampled(mtype(im))+1
            xc1=xadd(im)
            yc1=yadd(im)
            zc1=zadd(im)
            jm1=1
            jm2=nmol
            if(probaccept.ge.1.0) then
               jm1=im+1
               jm2=nmol
            end if
            do ia=naf(im),nal(im)
               it=itype(ia)
               csampled(it)=csampled(it)+1
            end do
            do ia=naf(im),nal(im)
               x1=xa(ia)+xc1
               y1=ya(ia)+yc1
               z1=za(ia)+zc1
               it=itype(ia)
!$omp parallel do private(jm,xc2,yc2,zc2,ja,jaref,jt,icref,x2,y2,z2,xx,yy,zz,dtest,dtesttest,iz1,irefp,scalefac) shared(iz1values,icrefvalues,potex)
               do jm=jm1,jm2
                  if(jm.ne.im) then
                     xc2=xadd(jm)
                     yc2=yadd(jm)
                     zc2=zadd(jm)
                     do ja=naf(jm),nal(jm)
                        jt=itype(ja)
                        icref=npairref(it,jt)
                        x2=xa(ja)+xc2
                        y2=ya(ja)+yc2
                        z2=za(ja)+zc2
                        xx=x1-x2
                        yy=y1-y2
                        zz=z1-z2
c
c apply minimum image convention
c
                        call minimagec(xx,yy,zz)
                        dtest=sqrt(xx*xx+yy*yy+zz*zz)
c find the x value this corresponds to
                        iz1=nint(dtest/stepcl)+2
c      write(6,*) ia,ja,xx,yy,zz,dtest,stepcl,iz1
c
c accumulate the partial g(r) for this pair
c
                        if(iz1.gt.1.and.iz1.le.npt1) then
                           iz1values(ja)=iz1
                           icrefvalues(ja)=icref
                        else
                           iz1values(ja)=0
                        endif
c
c now check whether this pair of atoms is below a specified minimum distance or above a specified maximum distance
c
                        dtesttest=dtest*dmax
                        scalefac=1.0
                        if(coarsegraining) scalefac=0.5*(scalefacmol(im)+scalefacmol(jm))
                        irefp=irefpcof(icref)
                        if(pcoefex(irefp).ge.0.0) then
                           if(rminex(irefp).gt.0.0.and.dtesttest.lt.scalefac*rminex(irefp)) then
                              potex(icref)=potex(icref)+addwav
                           endif
                        end if
                     end do
                  else
                     do ja=naf(jm),nal(jm)
                        iz1values(ja)=0
                     end do
                  end if
               end do
!$omp end parallel do
!               write(6,*) 'atogen> ',im,ia,jm1,jm2,naf(jm1),nal(jm2),addwav
               do ja=naf(jm1),nal(jm2)
                  iz1=iz1values(ja)
                  if(iz1.gt.0) then
                     icref=icrefvalues(ja)
                     wav(2,iz1,icref)=wav(2,iz1,icref)+addwav
                  end if
               end do
            end do
         end if
! Check overlaps if container forces are in use
         if(numcontainerforces.gt.0.0) then !Deal with container forces
            call dist_containerforce(im,naf(im),nal(im),xadd,yadd,zadd,xa,ya,za
     1,itype,mref,npairref,irefpcof,rmaxex,potexmax,dmax,addwav)
         endif
      end do
c
c normalize sums
c
      do ic=1,npairs
         it=n1extr(ic)
         jt=n2extr(ic)
         cprod=cfrac(it)*cfrac(jt)*real(natomr)
         cprod2=csampled(it)*cfrac(jt)
         if(it.ne.jt) then
            cprod=cprod*2
            cprod2=cprod2+csampled(jt)*cfrac(it)
         end if
!         write(6,*) 'atogen> ',ic,cprod,cprod2
         if(cprod.gt.0.0) then
            do iz=1,npt1
               wav(1,iz,ic)=wav(1,iz,ic)/cprod
            end do
         else
            do iz=1,npt1
               wav(1,iz,ic)=0.0
            end do
         endif
         if(cprod2.gt.0.0) then
            do iz=1,npt1
               wav(2,iz,ic)=wav(2,iz,ic)/cprod2
            end do
         else
            do iz=1,npt1
               wav(2,iz,ic)=ycoh(iz)
            end do
         endif
!         write(6,'(4(1x,i4),5(1x,e12.5))') ic,it,jt,natomr,csampled(it),cfrac(it),cfrac(jt),cprod,cprod2
      end do

      if(allocated(csampled)) deallocate(csampled,iz1values,icrefvalues,molsampled,molfraction)
      return
      end

      subroutine press(energy,pressc,energf,pressf,energabs,energthreeb,prthreebout)

c calculate energy and pressure from intermolecular interactions

      use molecom_mod
      use ergpairreference
      use threebody_mod

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      include 'weights.inc'
      real energf(mfile),pressf(mfile),term,cfracprod,stepcl,energthreeb,enthreeb,prthreebout
      real*8 prpot,presst,pressthreeb,prthreeb
      integer im,jm,jref

      stepcl=cellst/dmax
      sumc=0.
      presst=0.
c
c absolute energy of EP coefficients
c
      energabs=abs_energy_ep()
      do im=1,nmol-1
         !Get neighbour list
         call make_neighbour_list_from_mols(im)
         xc1=xadd(im)
         yc1=yadd(im)
         zc1=zadd(im)
!Zero the energy and pressure accumulators
         do jref=1,nnbref
            ergold(jref)=0
            ergnew(jref)=0
         end do
!$omp parallel do private(jref,jm,ia,it,x1,y1,z1,ja,xc2,yc2,zc2,jt,icref,x2,y2,z2,xdif,ydif,zdif,xcdif,ycdif,zcdif,dsq,dtest,rat,iz,dtdif,ipz,prfac,enpot,prpot,enwt,enabswt,prwt) shared(ergold,ergnew)
         do jref=1,nnbref
            jm=nbref(jref)
!Avoid double counting pairs
            if(jm.gt.im) then
               do ia=naf(im),nal(im)
                  it=itype(ia)
                  x1=xa(ia)+xc1
                  y1=ya(ia)+yc1
                  z1=za(ia)+zc1
                  do ja=naf(jm),nal(jm)
                     xc2=xadd(jm)
                     yc2=yadd(jm)
                     zc2=zadd(jm)
                     jt=itype(ja)
                     icref=npairref(it,jt)
c
c set up separation of atom pair
c
                     x2=xa(ja)+xc2
                     y2=ya(ja)+yc2
                     z2=za(ja)+zc2
                     xdif=x2-x1
                     ydif=y2-y1
                     zdif=z2-z1
c
c apply minimum image convention
c
                     call minimagec(xdif,ydif,zdif)
c
c separation of centres (for pressure calculation)
c
                     xcdif=xc2-xc1
                     ycdif=yc2-yc1
                     zcdif=zc2-zc1
                     call minimagec(xcdif,ycdif,zcdif)
                     dsq=xdif*xdif+ydif*ydif+zdif*zdif
                     dtest=sqrt(dsq)
c find the r value this corresponds to
                     rat=dtest/stepcl
                     iz=int(rat)+2
                     dtdif=dtest-xstar(iz)
                     ipz=iz
                     if(dtest.gt.0.0) then
                        prfac=(xcdif*xdif+ycdif*ydif+zcdif*zdif)/dtest
                     else
                        prfac=0.0
                     endif
c sum reference potential
                     if(iz.le.npt1) then
                        enpot=ecoul(iz,icref)+ecoul1(iz,icref)*dtdif
                        prpot=0.5*(ecoul1(ipz,icref)+ecoul1(ipz-1,icref))*prfac
c
c sum intermolecular energy and pressure from empirical potential
c
                        enwt=xavs(iz,icref)+xavs1(iz,icref)*dtdif
c For absolute empirical energy do not include the repulsive term (but this is no longer included
c in xavs, so we can use xavs).
                        enabswt=enwt
                        prwt=0.5*(xavs1(ipz,icref)+xavs1(ipz-1,icref))*prfac
                        ergold(jref)=ergold(jref)+enwt+enpot
                        ergnew(jref)=ergnew(jref)+prwt+prpot
                     end if
                  end do
               end do
            end if
         end do
!$omp end parallel do
!Accumulate values
         do jref=1,nnbref
            sumc=sumc+ergold(jref)
            presst=presst+ergnew(jref)
         end do
      end do
!Add on the triple body terms if present
      energthreeb=0.0
      pressthreeb=0.0 
      if(nthreeb.gt.0) then
         call find_triplets(nmol,natom,xadd,yadd,zadd,xa,ya,za,mref,itype,ntype,type)
         do im=1,nmol
            energthreeb=energthreeb+en_mole_triplets(im,xadd,yadd,zadd,phix,phiy,phiz,xa,ya,za,mref,tethered)
            pressthreeb=pressthreeb+pr_mole_triplets(im,xadd,yadd,zadd,xa,ya,za,mref)
         end do
         sumc=sumc+energthreeb
         presst=presst+pressthreeb
         energthreeb=energthreeb/nmol
         prthreebout=-100*(pressthreeb/3.0)/(6.0222*volume)
      end if
c energy in kJ/mole
      energy=sumc/nmol
c pressure in kbar
      pressc=((nmol*simt-presst/3.0)*100.0)/(6.0222*volume)
      do if=1,nfile
            energf(if)=energf(if)/(nmol*wsum(if))
            pressf(if)=((-pressf(if)/3.0)*100.0)
     1/(6.0222*volume*wsum(if))
      end do

      return
      end

      subroutine stog(rho,q,sofq,r,gofr,nqmin,nqmax,nr,npairs
     *,nqdim,nrdim,npdim,xwt)
      real q(nqdim),r(nrdim),sofq(nqdim,npdim),gofr(nrdim,npdim)
     *,xwt(npdim)
      real lorch,sumi
      integer nqmin(npdim),nqmax(npdim)
      real, dimension (:,:), allocatable :: sum

      pi=3.141592653
      qstep=q(4)-q(3)
      sfac=0.5*qstep/pi/pi/rho
c
c perform a direct transform on sofq for r values, using a Lorch function to remove
c excess truncation oscillations
c
      pi=4*atan(1.0)
      qstep=q(4)-q(3)
      sfac=0.5*qstep/pi/pi/rho
      allocate(sum(nqdim,npairs))
      do j=1,nr
         rr=r(j)
         if(rr.gt.0.0) then
!$omp parallel do private(in,qmax,i,qq,qr,prod) shared(sum)
            do in=1,npairs
               qmax=pi/q(nqmax(in))
!######!$omp parallel do private(i,qq,qr,prod) shared(sum)
               do i=nqmin(in),nqmax(in)
                  qq=q(i)
                  qr=qq*rr
                  prod=qq*sin(qr)*lorch(qmax*qq,pi)
                  sum(i,in)=sofq(i,in)*prod
               end do
            end do
!$omp end parallel do
            do in=1,npairs
               sumi=0.0
               do i=nqmin(in),nqmax(in)
                  sumi=sumi+sum(i,in)
               end do
               gofr(j,in)=xwt(in)+sfac*sumi/rr
            end do
         else
!$omp parallel do private(in,qmax,i,qq,prod) shared(sum)
            do in=1,npairs 
               qmax=pi/q(nqmax(in))
!######!$omp parallel do private(i,qq,prod) shared(sum)
               do i=nqmin(in),nqmax(in)
                  qq=q(i)
                  prod=qq*qq*lorch(qmax*qq,pi)
                  sum(i,in)=sofq(i,in)*prod
               end do
            end do
!$omp end parallel do
            do in=1,npairs
               sumi=0.0
               do i=nqmin(in),nqmax(in)
                  sumi=sumi+sum(i,in)
               end do
               gofr(j,in)=xwt(in)+sfac*sumi
            end do
         endif
      end do
      deallocate(sum)
      return
      end

      function lorch(x,pi)

      real*4 lorch,x,pi

      if(x.eq.0) then

         lorch=1.0

      else if(x.lt.pi) then

         lorch=sin(x)/x

      else

         lorch=0.0

      endif

      return
      end

      real function abs_energy_ep()

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'epsrinput_par.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'

      real energabs,cfracprod,term
      integer ip,irefp,i
c
c absolute energy of EP
c
      energabs=0.0
c Sum the absolute values of the potentials directly
c         
      do ip=1,npairs
         irefp=irefpcof(ip)
         cfracprod=cfrac(n1extr(ip))*cfrac(n2extr(ip))
         do i=1,npt1
            term=ycoh(i)*cfracprod*abs(xavs(i,irefp))
!            term=cfracprod*abs(xavs(i,ip))
            energabs=energabs+term
         end do
      end do
      abs_energy_ep=energabs

      return
      end

      subroutine adjust_intra_wts()

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'

      real sumevib,sumang,sumdih,ratioang,ratiodih,ratiomin,ratiomax
     *,energy,en_mole_intra
      integer im,imoltyp
! Adjust the bond angle energy weight
! First we try to set the intra-molecular weights to give the required energy
! ratios between vibrational, angle and dihedral
      sumevib=0.0
      sumeang=0.0
      sumedih=0.0
      tempang=1.0
      tempdih=10.0
      do im=1,nmol
      !Get the intramolecular energy of this molecule
         energy=en_mole_intra(im,sumevib,sumeang,sumedih)
!      write(6,*) im,nergpairm(im)
      end do
      write(6,*) 'adjust_intra_wts> Vibrational, bond angle and dihedral energies: ',sumevib,sumeang,sumedih
      if(sumevib.gt.0.0) then
         ratioang=sumeang/sumevib
         ratiodih=sumedih/sumevib
      else
         ratioang=1.0
         ratiodih=1.0
      end if
      write(6,*) 'adjust_intra_wts> Relative bond angle and dihedral energies: ',ratioang,ratiodih
!      ratiomax=1.01
!      ratiomin=0.99
! Now calculate the relative weights on the bond angle energies
!      if(tempangreq.gt.0.0) then
!         if(tempangwt.gt.0.0) then
!            ratioang=ratioang/tempangreq
!            if(ratioang.gt.ratiomax) ratioang=ratiomax
!            if(ratioang.lt.ratiomin) ratioang=ratiomin
!            tempangwt=tempangwt*ratioang
!         else
!            tempangwt=tempangreq
!         end if
!      else
!         tempangwt=0.0
!      end if
! Now calculate the relative weights on the dihedral energy
!      if(tempdihreq.gt.0.0) then
!         if(tempdihwt.gt.0.0) then
!            ratiodih=ratiodih/tempdihreq
!            if(ratiodih.gt.ratiomax) ratiodih=ratiomax
!            if(ratiodih.lt.ratiomin) ratiodih=ratiomin
!            tempdihwt=tempdihwt*ratiodih
!         else
!            tempdihwt=tempdihreq
!         end if
!      else
!         tempdihwt=0.0
!      end if
      tempangwt=tempangreq
      tempang=tempangwt
      tempdihwt=tempdihreq
      tempdih=tempdihwt*tempb
      write(6,*) 'adjust_intra_wts> Bond angle and dihedral energy weights: ',tempang,tempdih
      do imoltype=1,nmoltype
         write(6,*) 'adjust_intra_wts> Internal step sizes: ',imoltype,interntrans(imoltype),internrot(imoltype)
      end do
      return
      end

      subroutine get_maxmolsize()

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      real test
      integer im,ia

      maxmolsize=0.0
      do im=1,nmol
         if(moldiam(im).gt.maxmolsize) maxmolsize=moldiam(im)
      end do
      maxmolsize=2.0*maxmolsize
      
      end subroutine get_maxmolsize


