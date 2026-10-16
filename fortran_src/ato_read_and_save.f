      subroutine minitc(fname1,rmax,cell,upol,rho,amubar,amubar2,ndist,ierr)

      use reallocation_routines
      use molecom_mod
      use threebody_mod
c
c this is a modified version of minitx from epsr15x to include a purely
c cubic box
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      include 'lattice.inc'
      character*256 fname1
      character*3 a1,ats
      integer*4 iexch,j,k
      external cosd,sind
      real pi,piconv,natomsph,rhosph,radminsph,radmaxsph,volsph,pseudofacsph

      pi=4.0*atan(1.0)
      piconv=pi/180.0

! Create the lock file and write set to true (= locked)

!      write(6,'(a)') fname1(1:len_trim(fname1))//'.lck'
!      open(10,file=fname1(1:len_trim(fname1))//'.lck',status='unknown',iostat=ierr)
!      write(10,'(a)') 'T'
!      close (10)

      open(10,file=fname1,status='old',iostat=ierr)
c
c if there is a problem with this file, print out an error message and exit
c
      if(ierr.ne.0) then
         ilenfname=index(fname1,'.')
         ilenfname=ilenfname+index(fname1(ilenfname:256),' ')-1
         write(6,599) fname1(1:ilenfname)
599   format(/'Problems with file ',a/)
         return
      endif

c Decide whether this is a cubic or a non-cubic simulation box

      cubic=.true.
      read(10,'(a)',iostat=ierr) sometext
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.ne.3) then

         cubic=.false.
c read no. of atoms in molecule, temperature in Kelvin, and lattice vectors
         read(sometext,*,iostat=ierr) nmol,simtk
         read(10,*,iostat=ierr) alat,blat,clat
         read(10,*,iostat=ierr) phib,thetac,phic

      else

         cubic=.true.
         read(sometext,*,iostat=ierr) nmol,dmax,simtk
         alat=dmax
         blat=dmax
         clat=dmax
         phib=90.0
         thetac=0.0
         phic=0.0

      endif
      write(6,"('minitc> cubic: ',L10)") cubic
! Assign initial scale factors for each molecule if not defined on input
      call reallocate1d_r(scalefacmol,nmol)
      scalefacmol=1.0
c
c read fitting tolerence, step sizes, and temperatures for harmonic terms
c
      read(10,'(a)',iostat=ierr) sometext
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      read(sometext,*,iostat=ierr) tol,stepmi,stepri,stepr,stepm,tempb
      tempangreq=1.0
      tempdihreq=10.0
      tempangwt=1.0
      tempdihwt=1.0
      if(nwords.eq.7) then
         read(sometext(ncf(7):ncl(7)),*,iostat=ierr) tempdihreq
      end if
      if(nwords.gt.7) then
         read(sometext(ncf(7):ncl(8)),*,iostat=ierr) tempangreq,tempdihreq
      end if
      if(nwords.gt.9) then
         read(sometext(ncf(9):ncl(10)),*,iostat=ierr) tempangwt,tempdihwt
      end if
      tempang=tempangwt
      tempdih=tempdihwt*tempb
      tolout=tol
      ia=0
      ig=0
      natom=0
      natomr=0
      ngrrot=0
      nrottot=0
! Reset array size counters
      ngrrotmax=0
      nrottotmax=0
      do im=1,nmol
         read(10,'(a)',iostat=ierr) sometext
         call parse(sometext,nwords,ncf,ncl,10,nstsize,nworddim)
         read(sometext,*,iostat=ierr) nato(im),xadd(im),yadd(im),zadd(im),phix(im),phiy(im)
     1,phiz(im)
         tetheredatom(im)=0
         if(nwords.gt.7) then
            tethered(im)=sometext(ncf(8):ncf(8)).eq.'T'
            if(ncl(8).gt.ncf(8)) then
               read(sometext(ncf(8)+1:ncl(8)),*,iostat=ierr) tetheredatom(im)
               if(ierr.ne.0) tetheredatom(im)=0
            end if
         else
            tethered(im)=.false.
         end if
         if(nwords.gt.8) then
            read(sometext(ncf(9):ncl(9)),*) nstepreject(im)
         else
            nstepreject(im)=0
         end if
         if(nwords.gt.9) then
            read(sometext(ncf(10):ncl(10)),*) scalefacmol(im)
         else
            scalefacmol(im)=1.0
         end if
!         write(6,*) im,nwords,nstepreject(im),scalefacmol(im)
         natom=natom+nato(im)
         naf(im)=ia+1
c set a flag to determine if this molecule has any Q-type atoms
         iflag=0
         do ib=1,nato(im)
            ia=ia+1
            mref(ia)=im
c for each atom read its label
            read(10,'(a)',iostat=ierr) sometext
            call parse(sometext,nwords,ncf,ncl,10,nstsize,nworddim)
            read(sometext(ncf(1):ncl(1)),'(a3)',iostat=ierr) alabel(ia)
            !Read the image atom if present.
            if(nwords.gt.2) then
               read(sometext(ncf(3):ncl(3)),*,iostat=ierr) imageatom(ia)
               !As a precaution, don't allow unphysical values
               if(imageatom(ia).lt.0) imageatom(ia)=0
               if(imageatom(ia).gt.nato(im)) imageatom(ia)=0
            else
               imageatom(ia)=0
            end if
c set up 'type' of atom - this is used in assigning neutron weights
            call ftype(ia)
c
c do not count atoms which begin with 'Q' or 'q' - these are assumed 
c invisible to the scattering radiation and experiment (although they
c can be used to calculate other structural quantities if necessary
c
            if(qtype(ia).ne.1) then
                  natomr=natomr+1
            else
                  iflag=1
            endif
c read its starting coordinates and sigmas
            read(10,*,iostat=ierr) xa(ia),ya(ia),za(ia)
c read the number of atoms to which it must be bound, the atom numbers of
c those neighbours and the anticipated bond length 
            read(10,*,iostat=ierr) nbond(ia),(iabond(ia,j),tbond(ia,j)
     1,j=1,nbond(ia))
         end do
         nal(im)=ia
c
c reset positions of Q-type atoms
c
c     if(iflag.eq.1) call resetq(im)
c read in number of internal rotation groups
         read(10,*,iostat=ierr) ngrot(im)
         nrf(im)=ig+1
         if(ngrot(im).gt.0) then
            ngrrot=ngrrot+ngrot(im)
            call testsize_rotgrps()
            do ih=1,ngrot(im)
                  ig=ig+1
c save molecule in which this rotational group occurs
                  mrot(ig)=im
!                  write(6,*) 'minitc> ',ig,ngrrot,ngrrotmax,size(rotgrp)
                  read(10,'(1x,a3)',iostat=ierr) rotgrp(ig)
c for rotational groups read atom numbers of two atoms to form axis of rotation
                  read(10,*,iostat=ierr) labaxi(ig,1),labaxi(ig,2)
! For DIH groups read number of dihedral angles defined for a given axis, and
! two atom numbers to calculate the dihedral angle and the required dihedral angle 
                  if(rotgrp(ig).eq.'DIH') then
                     read(10,*,iostat=ierr) ndihedralang(ig)
                     do k=1,ndihedralang(ig)
                        read(10,*,iostat=ierr) (labrotdih(j,k,ig),j=1,2),dihedralang(k,ig)
                     end do
                  else
                     ndihedralang(ig)=0
                  end if
c Read number of atoms in group and their numbers
! Read the first line
                  read(10,'(a)',iostat=ierr) sometext
                  call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
                  read(sometext(ncf(1):ncl(1)),*) nrot(ig)
                  rotf(ig)=nrottot+1
                  nrottot=nrottot+nrot(ig)
                  rotl(ig)=nrottot
                  call testsize_labrot()
! Read the remaining words on this line. The first has already been read so we start from 2
                  iw=1
                  labref=rotf(ig)
                  do while(iw.lt.nwords)
                     iw=iw+1
                     labref=rotf(ig)+iw-2
                     read(sometext(ncf(iw):ncl(iw)),*) labrot(labref)
                  end do
! Read the remaining groups directly if needed
                  if(labref.lt.rotl(ig)) then
                     labref=labref+1
                     read(10,*,iostat=ierr) (labrot(j),j=labref,rotl(ig))
                  end if
            end do
         endif
         nrl(im)=ig
      end do
c
c check bond atom numbers and rotational group atom numbers are correct
c if not convert to correct values
c
      call bondnumbers
c
c make a list of the atom types read in from the current file
c
      ntypefile=0
      ia=0
      do while(ia.lt.natom)
         ia=ia+1
         itypefile=0
         ifound=0
         do while(itypefile.lt.ntypefile.and.ifound.eq.0)
            itypefile=itypefile+1
            if(alabel(ia).eq.typefile(itypefile)) ifound=itypefile
         end do
         if(ifound.eq.0) then
            ntypefile=ntypefile+1
            typefile(ntypefile)=alabel(ia)
         endif
      end do
c
c make a list of the molecule types read in from the current file
c
      nmoltypefirst=nmoltype
      im=0
      do while(im.lt.nmol)
         im=im+1
         imoltype=0
         ifound=0
         do while(imoltype.lt.nmoltype.and.ifound.eq.0)
            imoltype=imoltype+1
            if(alabel(naf(im)).eq.moltype(imoltype)) then
               ifound=imoltype
               mtype(im)=ifound
            end if
         end do
         if(ifound.eq.0) then
            nmoltype=nmoltype+1
            moltype(nmoltype)=alabel(naf(im))
            moltethered(nmoltype)=tethered(im)
            moltetheredatom(nmoltype)=tetheredatom(im)
c
c initialise the template filename
c
            write(moltemplate(nmoltype),'("moltype",i2.2)') nmoltype
            mtype(im)=nmoltype

c Initialise the rotational and translational step sizes

            wholerot(nmoltype)=stepr
            wholetrans(nmoltype)=stepm

         endif
      end do
! Initialise the relative broadening of each molecule type if coarse graining
      call reallocate1d_r(rfractmol,nmoltype)
      rfractmol=0.0

c Write these out as a check

      write(6,201)
201   format(/'minitc> Following molecule types defined initially:-'/)

      do imoltype=1,nmoltype

         write(6,202) imoltype,moltype(imoltype),moltemplate(imoltype)(1:10),wholerot(imoltype)
     *,wholetrans(imoltype)
202   format(1x,i3,1x,a3,1x,a10,2(1x,e12.5))

      end do
!      do it=1,ntypefile
!            write(6,'(a,i5,1x,a)') 'minitc> ',it,typefile(it)
!      end do

c
c read Lennard-Jones, mass and charge parameters for each atom type. 
! rhosphere,radminsphere and radmaxsphere (used in coarse graining) are allocated in ctype.
c
      do it=1,ntypefile
            read(10,'(a)',iostat=ierr) sometext
c
c parse this text into words
c
            call parse(sometext,nwords,ncf,ncl,8,nstsize,nworddim)
            a1=sometext(ncf(1):ncl(1))
            if(nwords.gt.1) then
               ats=sometext(ncf(2):ncl(2))
            else
               ats='XXX'
            endif
            if(nwords.gt.2) then
               read(sometext(ncf(3):ncl(3)),*) iexch
            else
               iexch=-1
            endif
            read(10,'(a)',iostat=ierr) sometext
            call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
            read(sometext(ncf(1):ncl(5)),*,iostat=ierr) epsit,sigt,amass,qt,rqt
!Get the coarse graining parameters if present
            if(nwords.gt.5) then
               read(sometext(ncf(6):ncl(6)),*,iostat=ierr) natomsph !If this is negative it represents a density, not the number of atoms in the sphere.
            else
               natomsph=1.0
            end if
            if(nwords.gt.7) then
               read(sometext(ncf(7):ncl(8)),*,iostat=ierr) radminsph,radmaxsph
               if(radmaxsph.le.0.0) radmaxsph=0.5*sigt !In case default value is read in, set this to radius of atom
            else
               radminsph=0.0
               radmaxsph=0.5*sigt
            end if
            if(nwords.gt.8) then
               read(sometext(ncf(9):ncl(9)),*,iostat=ierr) pseudofacsph
               if(pseudofacsph.lt.1.0) pseudofacsph=1.0
            else
               pseudofacsph=1.0
            end if
! Get the density of this atom type
            volsph=4.0*pi*(radmaxsph**3-radminsph**3)/3.0
            if(abs(natomsph).gt.0.0) then
               if(natomsph.lt.0.0) then
                  rhosph=abs(natomsph)
                  natomsph=rhosph*volsph/pseudofacsph
               else if(volsph.gt.0.0) then
                  rhosph=(natomsph*pseudofacsph)/volsph
               else
                  rhosph=0.0
               end if
            else
               rhosph=0.0
            end if
c assign type number
            call ctype(a1,ats,iexch,epsit,sigt,amass,qt,rqt,natomsph,rhosph,radminsph,radmaxsph,pseudofacsph)
      end do
c
c read core energy and diameter - used in fmole for non-bonded atoms
c
      read(10,*,iostat=ierr) ecore,dcore
c
c initialise random number generator
c
      if(ierr.eq.0) then
         call initran1(10)
c
c read the names of the template files if present. Normally these should be 
c specified in the order in which the different molecule types occur in the 
c .ato file.
c
         icount=0
         do while(icount.lt.(nmoltype-nmoltypefirst).and.ierr.eq.0)
            icount=icount+1
            read(10,'(a)',iostat=ierr) sometext
c
c parse this text into words
c
            if(ierr.eq.0) then
               call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
               if(nwords.gt.1) then
c
c get the molecular type this line belongs to
c
                  read(sometext(ncf(1):ncl(1)),*,iostat=ierr) imolread
                  imoltype=nmoltypefirst+imolread
                  moltemplate(imoltype)=sometext(ncf(2):ncl(2))
                  interntrans(imoltype)=0.0
                  internrot(imoltype)=0.0
                  wholerot(imoltype)=0.0
                  wholetrans(imoltype)=0.0
c Read the rotational and translational step sizes for each molecule type if specified
                  if(nwords.gt.2) read(sometext(ncf(3):ncl(3)),*,iostat=ierr) interntrans(imoltype)
                  if(nwords.gt.3) read(sometext(ncf(4):ncl(4)),*,iostat=ierr) internrot(imoltype)
                  if(nwords.gt.4) read(sometext(ncf(5):ncl(5)),*,iostat=ierr) wholerot(imoltype)
                  if(nwords.gt.5) read(sometext(ncf(6):ncl(6)),*,iostat=ierr) wholetrans(imoltype)
                  if(nwords.gt.6) read(sometext(ncf(7):ncl(7)),*,iostat=ierr) rfractmol(imoltype) 
               end if
            end if
         end do

c set the error code back to 0 so that if the template file has not been specified
c routines which do not need this template can still run

         ierr=0
         close (10)

      else

         close (10)
         return

      endif

c Write moltypes out as a check

      write(6,201)

      do imoltype=1,nmoltype

         write(6,203) imoltype,moltype(imoltype)
     *,moltemplate(imoltype)(1:10)
     *,interntrans(imoltype),internrot(imoltype)
     *,wholerot(imoltype),wholetrans(imoltype),rfractmol(imoltype)
203   format(1x,i3,1x,a3,1x,a10,5(1x,e12.5))

      end do
c set up lattice parameters for applying the minimum image convention

      alatx=alat
      blatx=blat*cosd(phib)
      blaty=blat*sind(phib)
      clatx=clat*sind(thetac)*cosd(phic)
      claty=clat*sind(thetac)*sind(phic)
      clatz=clat*cosd(thetac)
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
c
c calculate the perpendicular distance between the faces of the lattice
c
      perpab=perpd(alatc,blatc,clatc)
      perpbc=perpd(blatc,clatc,alatc)
      perpca=perpd(clatc,alatc,blatc)
      write(6,*) perpab,perpbc,perpca
c
c set up maximum radius to be used for Fourier transforms and to truncate
c the potential. This is defined as the minimum perpendicular distance between
c planes. 
c
c G(r) and potentials will be calculated to rmax, or maximum allowed distance
c based on number of points and step size, whichever is smaller.
c
      dmax=min(perpab,perpbc,perpca)
      write(6,*) nmol,dmax
      cell=dmax
c
c setup molecular size factor
c
      dmaxmol=dmax*sizefactor
c
c determine maximum distance across box
c
c G(r) and potentials will be calculated to rmax, or maximum allowed distance
c based on number of points and step size, whichever is smaller.
c
      xmax=abs(alatx)+abs(blatx)+abs(clatx)
      ymax=abs(blaty)+abs(claty)
      zmax=abs(clatz)
      rmax=0.5*sqrt(xmax*xmax+ymax*ymax+zmax*zmax)
      write(6,*) rmax,cell
c
c convert temperature in K to kJ/mole
c
      simt=0.008314*simtk
      write(6,917) nmol,natom,natomr,ntype,dmax
917   format(1x,'no. of molecules to read = ',4(1x,i0),1x,e13.6,1x,i0)
c
c calculate atomic fractions
c
      natomcg=0.0
      natomrcg=0.0
      do it=1,ntype
            cfrac(it)=0.0
            cfraccg(it)=0.0
      end do
      do ia=1,natom
            it=itype(ia)
            natomcg=natomcg+natomsphere(it)*pseudofacmol(it)
            if(qtype(ia).ne.1) natomrcg=natomrcg+natomsphere(it)*pseudofacmol(it)
            cfraccg(it)=cfraccg(it)+natomsphere(it)*pseudofacmol(it) !Beware! This gets confusing because the natomr, cfrac, rho
! represents the fraction of spheres, whereas natomrcg, cfraccg, rhocg represents the fraction of atoms.
            cfrac(it)=cfrac(it)+1.0
      end do
c
c calculate volume of box and atomic number density
c
      volume=alatx*blaty*clatz
      rho=real(natomr)/volume
      rhocg=natomrcg/volume
      write(6,*) volume,natomcg,natomrcg,rho,rhocg
      if(natomrcg.gt.0.0) then
         do it=1,ntype
            cfraccg(it)=cfraccg(it)/natomrcg
            cfrac(it)=cfrac(it)/real(natomr)
         end do
      else if(natomcg.gt.0.0) then
         do it=1,ntype
            cfraccg(it)=cfraccg(it)/natomcg
            cfrac(it)=cfrac(it)/real(natom)
         end do
      end if
      sumcharge=0.0
      do it=1,ntype
            write(6,916) it,type(it),cfrac(it),cfraccg(it)
            sumcharge=sumcharge+cfrac(it)*charge(it)
      end do
916   format(1x,'Atomic fraction (actual and coarse grained)',i4,1x,a,' = ',e12.5,1x,e12.5)
c
c calculate mass of each molecule
c
      do im=1,nmol
            amasst(im)=0.0
            do ia=naf(im),nal(im)
                  if(qtype(ia).eq.0) then
                     amass=amassa(itype(ia))
                     amasst(im)=amasst(im)+amass
                  endif
            end do
            amasstsqrt(im)=sqrt(amasst(im))
c
c set maximum distance of each atom from centre of mass
c
            moldiam(im)=0.0
            do ia=naf(im),nal(im)
                  sum=xa(ia)*xa(ia)+ya(ia)*ya(ia)+za(ia)*za(ia)
                  if(sum.gt.moldiam(im)) moldiam(im)=sum
            end do
            moldiam(im)=sqrt(moldiam(im))
c           write(6,*) im,moldiam(im)
      end do
      ip=0
      do it=1,ntype
            do jt=it,ntype
                  ip=ip+1
                  n1extr(ip)=it
                  n2extr(ip)=jt
                  npairref(it,jt)=ip
                  npairref(jt,it)=ip
                  ect=sqrt(epscor(it)*epscor(jt))
                  ecor(ip)=ect
                  qcor(ip)=charge(it)*charge(jt)*1390
                  dct=0.5*(sigcor(it)+sigcor(jt))
                  dcor(ip)=dct
                  dscor(ip)=dct*dct
                  am1=amassa(it)
                  am2=amassa(jt)
!                  amassr(ip)=tempb*simt*sqrt(am1*am2/(am1+am2))
                  amassr(ip)=sqrt(am1*am2/(am1+am2))
c                 write(6,*) it,jt,sigcor(it),sigcor(jt),dcor(ip)
            end do
      end do
      npairs=ip
      ndist=npairs
      write(6,'(5(1x,i0))') nmol,natomr,ntype,npairs,size(rhosphere)
c
c     calculate dipole moments of molecules. Currently it is assumed the 
c the default dipole moment is that of an isolated water molecule
c
      amu0=1.85
      amu0=amu0/4.80273
      amusum=0.0
      amubar=0.0
      nbar=0
      do im=1,nmol
            if(nato(im).gt.1) then
                  qx=0.0
                  qy=0.0
                  qz=0.0
                  do ia=naf(im),nal(im)
                    qt=charge(itype(ia))
                    qx=qx+qt*xa(ia)
                    qy=qy+qt*ya(ia)
                    qz=qz+qt*za(ia)
                  end do
                  amu2=qx*qx+qy*qy+qz*qz
                  amu=sqrt(amu2)
                  amud=amu*4.80273
                  nbar=nbar+1
                  amubar=amubar+amud
                  amubar2=amubar2+amud*amud
c                 write(6,997) im,amu*4.80273
997   format(1x,'Molecule ',i5,' dipole moment',e13.6)
                  amudif=amu-amu0
                  amusum=amusum+amudif*amudif
            end if
      end do
      if(nbar.gt.0) then
            amubar=amubar/nbar
            amubar2=amubar2/nbar-amubar*amubar  
            if(amubar2.gt.0.0) then
                  amubar2=sqrt(amubar2)
            else  
                  amubar2=0.0
            endif
            amusum=amusum/nmol
c
c calculate static energy correction due to induced dipoles (c.f. Berendsen)
c
            upol=0.5*amusum*1390/1.445
      else
            amubar=0.0
            upol=0.0
      endif
      return
      end

      subroutine normalise_coordinates()

      use molecom_mod
      use threebody_mod
      use coordconstraints_mod
c
c sets up atom and molecule coordinates normalised to box dimension
c
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
c
c fractional step sizes
c
      call fractcoord(stepmi,dmaxmol)
      call fractcoord(stepm,dmax)
      call fractcoord(alatx,dmax)
      call fractcoord(blatx,dmax)
      call fractcoord(clatx,dmax)
      call fractcoord(blaty,dmax)
      call fractcoord(claty,dmax)
      call fractcoord(clatz,dmax)
      call fractcoord(alat,dmax)
      call fractcoord(blat,dmax)
      call fractcoord(clat,dmax)
      do i=1,3
         call fractcoord(alatc(i),dmax)
         call fractcoord(blatc(i),dmax)
         call fractcoord(clatc(i),dmax)
      end do
      if(stepm.gt.0.25) stepm=0.25
      do im=1,nmol
c
c set up fractional coordinates
c
            call fractcoord(xadd(im),dmax)
            call fractcoord(yadd(im),dmax)
            call fractcoord(zadd(im),dmax)
            if(tethered(im)) then
               call fractcoord(phix(im),dmax)
               call fractcoord(phiy(im),dmax)
               call fractcoord(phiz(im),dmax)
            end if
c
c set maximum distance of each atom from centre of mass
c
            call fractcoord(moldiam(im),dmaxmol)
            do ia=naf(im),nal(im)
                  call fractcoord(xa(ia),dmaxmol)
                  call fractcoord(ya(ia),dmaxmol)
                  call fractcoord(za(ia),dmaxmol)
                  if(nbond(ia).gt.0) then
                        do ib=1,nbond(ia)
                            call fractcoord(tbond(ia,ib),dmaxmol)
                        end do
                  end if
            end do
      end do
      do ip=1,npairs
         amassr(ip)=amassr(ip)*dmaxmol
      end do
! If any dihedral angles are defined, convert them to radians
      pi=4.0*atan(1.0)
      piconv=pi/180.0
      ir=0
      do while (ir.lt.ngrrot)
         ir=ir+1
         if(rotgrp(ir).eq.'DIH') then
            do k=1,ndihedralang(ir)
               dihedralang(k,ir)=dihedralang(k,ir)*piconv
            end do
         end if
      end do
! If any three-body potentials are defined, normalise the maximum distances
      if(nthreeb.gt.0) then
         do ir=1,nthreeb
            do j=1,2
                call fractcoord(threebrmax(j,ir),dmax)
            end do
         end do
      end if
!If any coordination constraints are defined, normalise specified bond lengths
      if(ncoordconstraints.gt.0) then
         do ir=1,ncoordconstraints
            do j=1,coordconstraints(ir)%ns
               write(6,*) ir,j,coordconstraints(ir)%ds(j),dmax
               call fractcoord(coordconstraints(ir)%ds(j),dmax)
               call fractcoord(coordconstraints(ir)%dsmax(j),dmax)
            end do
!            call fractcoord(coordconstraints(ir)%ccdistort,dmax)
!            if(coordconstraints(ir)%hasthreebody) then
!               call fractcoord(coordconstraints(ir)%threebodydistort,dmax)
!            endif
         end do
      end if
      return
      end

      subroutine unnormalise_coordinates

      use molecom_mod
      use threebody_mod
      use coordconstraints_mod
c
c sets up atom and molecule coordinates in unnormalised units
c
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
c
c fractional step sizes
c
      call wholecoord(stepmi,dmaxmol)
      if(stepm.gt.0.25) stepm=0.25
      call wholecoord(stepm,dmax)
      call wholecoord(alatx,dmax)
      call wholecoord(blatx,dmax)
      call wholecoord(clatx,dmax)
      call wholecoord(blaty,dmax)
      call wholecoord(claty,dmax)
      call wholecoord(clatz,dmax)
      call wholecoord(alat,dmax)
      call wholecoord(blat,dmax)
      call wholecoord(clat,dmax)
      do i=1,3
         call wholecoord(alatc(i),dmax)
         call wholecoord(blatc(i),dmax)
         call wholecoord(clatc(i),dmax)
      end do
      do im=1,nmol
c
c set up fractional coordinates
c
            call wholecoord(xadd(im),dmax)
            call wholecoord(yadd(im),dmax)
            call wholecoord(zadd(im),dmax)
            if(tethered(im)) then
               call wholecoord(phix(im),dmax)
               call wholecoord(phiy(im),dmax)
               call wholecoord(phiz(im),dmax)
            end if
c
c set maximum distance of each atom from centre of mass
c
            call wholecoord(moldiam(im),dmaxmol)
            do ia=naf(im),nal(im)
                  call wholecoord(xa(ia),dmaxmol)
                  call wholecoord(ya(ia),dmaxmol)
                  call wholecoord(za(ia),dmaxmol)
                  if(nbond(ia).gt.0) then
                        do ib=1,nbond(ia)
                            call wholecoord(tbond(ia,ib),dmaxmol)
                        end do
                  end if
            end do
      end do
      do ip=1,npairs
         amassr(ip)=amassr(ip)/dmaxmol
      end do
! If any dihedral angles are defined, convert them to degrees
      pi=4.0*atan(1.0)
      piconv=180.0/pi
      ir=0
      do while (ir.lt.ngrrot)
         ir=ir+1
         if(rotgrp(ir).eq.'DIH') then
            do k=1,ndihedralang(ir)
               dihedralang(k,ir)=dihedralang(k,ir)*piconv
            end do
         end if
      end do
! If any three-body potentials are defined, unnormalise the maximum distances
      if(nthreeb.gt.0) then
         do ir=1,nthreeb
            do j=1,2
                call wholecoord(threebrmax(j,ir),dmax)
            end do
         end do
      end if
!If any coordination constraints are defined, unnormalise specified bond lengths
      if(ncoordconstraints.gt.0) then
         do ir=1,ncoordconstraints
            do j=1,coordconstraints(ir)%ns
               call wholecoord(coordconstraints(ir)%ds(j),dmax)
               call wholecoord(coordconstraints(ir)%dsmax(j),dmax)
            end do
!            call wholecoord(coordconstraints(ir)%ccdistort,dmax)
!            if(coordconstraints(ir)%hasthreebody) then
!               call wholecoord(coordconstraints(ir)%threebodydistort,dmax)
!            endif
         end do
      end if
      return
      end

      subroutine moloutc(fname1)

      use reallocation_routines
      use molecom_mod
      use threebody_mod

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      integer, dimension(:), allocatable :: moltypeswritten
      character*256 fname1
      character*6 tetheredtext
      logical found
! If molecular scale factors have not been assigned, assign them
!      write(6,*) 'moloutc> ',nmol,size(scalefacmol),allocated(scalefacmol)
!      write(6,*) 'moloutc> ',ntype,size(rhosphere),allocated(rhosphere)
      if(.not.allocated(scalefacmol)) then
         call reallocate1d_r(scalefacmol,nmol)
         scalefacmol=1.0
      else if(size(scalefacmol,1).lt.nmol) then
         call reallocate1d_r(scalefacmol,nmol)
         scalefacmol=1.0
      end if
      open(10,file=fname1,status='unknown',iostat=ierr)
c
c Rewind in case a previous routine has left it open
c
      rewind(10)
c
c if there is a problem with this file, print out an error message and exit
c
      if(ierr.ne.0) then
         ilenfname=index(fname1,'.')
         ilenfname=ilenfname+index(fname1(ilenfname:256),' ')-1
         write(6,599) fname1(1:ilenfname)
599   format(/'moloutc> Problems with file ',a/)
         nmol=0
         return
      endif

c     Write a warning message if ROT groups are defined, but internal rotations are
c     switched off.

      if(ngrrot.gt.0.and.stepri.le.0.0) then
         write(6,598)
598   format(/'moloutc> WARNING! ROT groups defined BUT step size is ZERO!')
      endif

      ig=0
c
c write no. of molecules and tolerance on bond distances
c
      if(cubic) then

         write(10,207) nmol,dmax,simtk

      else

         write(10,207) nmol,simtk
         write(10,208) alat,blat,clat
         write(10,208) phib,thetac,phic

      endif
c set up lattice parameters for calculating volume and density
      alatx=alat
      blatx=blat*cosd(phib)
      blaty=blat*sind(phib)
      clatx=clat*sind(thetac)*cosd(phic)
      claty=clat*sind(thetac)*sind(phic)
      clatz=clat*cosd(thetac)
      write(10,206) tolout,stepmi,stepri,stepr,stepm,tempb,tempangreq,tempdihreq,tempangwt,tempdihwt
      natomr=0
      natom=0
      natomcg=0
      natomrcg=0
      do it=1,ntype
         cfrac(it)=0.0
         cfraccg(it)=0.0
      end do
      do im=1,nmol
      if(tetheredatom(im).lt.1.or.tetheredatom(im).gt.nato(im)) tetheredatom(im)=0
      if(tethered(im)) then
         write(tetheredtext,'(l1,i5.5)') tethered(im),tetheredatom(im)
      else
         write(tetheredtext,'(l1)') tethered(im)
      end if
      write(10,200) nato(im),xadd(im),yadd(im)
     1,zadd(im),phix(im),phiy(im),phiz(im),tetheredtext,nstepreject(im),scalefacmol(im),im
200   format(1x,i0,6(1x,e15.8),2x,a6,1x,i0,1x,f10.5,1x,i5)
c
c integer to convert from absolute atom numbers to relative atom numbers
c
      natom=natom+nato(im)
      nadd=1-naf(im)
      
      do ia=naf(im),nal(im)
c for each atom write its label
         if(imageatom(ia).gt.0) then
            write(10,201) alabel(ia), (ia+nadd),(imageatom(ia)+nadd)
         else
            write(10,201) alabel(ia), (ia+nadd),0
         endif

c write its starting coordinates and sigmas
            write(10,202) xa(ia),ya(ia),za(ia)
c
c write the number of atoms to which it must be bound, the atom numbers of
c those neighbours and the anticipated bond length 
c
            write(10,203) nbond(ia),((iabond(ia,j)+nadd),tbond(ia,j),j=1,nbond(ia))
            natomcg=natomcg+natomsphere(itype(ia))*pseudofacmol(itype(ia))
            if(qtype(ia).eq.0) then
               natomr=natomr+1
               natomrcg=natomrcg+natomsphere(itype(ia))*pseudofacmol(itype(ia))
            endif
            cfrac(itype(ia))=cfrac(itype(ia))+1.0
            cfraccg(itype(ia))=cfraccg(itype(ia))+natomsphere(itype(ia))*pseudofacmol(itype(ia))
      end do
c write number of internal rotation groups
      write(10,204) ngrot(im)
      if(ngrot(im).gt.0) then
            do ih=1,ngrot(im)
                  ig=ig+1
                  write(10,201) rotgrp(ig)
c for rotational groups read atom numbers of two atoms to form axis of rotation
                  write(10,204) (labaxi(ig,1)+nadd),(labaxi(ig,2)+nadd)
! For DIH groups read two atom numbers to calculate the dihedral angle and the required dihedral angle 
                  if(rotgrp(ig).eq.'DIH') then
                     write(10,204,iostat=ierr) ndihedralang(ig)
                     do k=1,ndihedralang(ig)
                        write(10,209,iostat=ierr) ((labrotdih(j,k,ig)+nadd),j=1,2),dihedralang(k,ig)
                     end do
                  end if
c write number of atoms in group and their numbers
                  write(10,204) nrot(ig),((labrot(j)+nadd),j=rotf(ig),rotl(ig))
            end do
      endif
      end do
      if(.not.allocated(natomsphere)) then
         call reallocate1d_r(natomsphere,ntype)
         call reallocate1d_r(rhosphere,ntype)
         call reallocate1d_r(radminsphere,ntype)
         call reallocate1d_r(radmaxsphere,ntype)
         call reallocate1d_r(pseudofacmol,ntype)
         natomsphere=1.0
         rhosphere=0.0
         radminsphere=0.0
         do it=1,ntype
            radmaxsphere(it)=0.5*sigcor(it)
         end do
         pseudofacmol=1.0
      else if(size(natomsphere,1).lt.ntype) then
         call reallocate1d_r(natomsphere,ntype)
         call reallocate1d_r(rhosphere,ntype)
         call reallocate1d_r(radminsphere,ntype)
         call reallocate1d_r(radmaxsphere,ntype)
         call reallocate1d_r(pseudofacmol,ntype)
         natomsphere=1.0
         rhosphere=0.0
         radminsphere=0.0
         do it=1,ntype
            radmaxsphere(it)=0.5*sigcor(it)
         end do
         pseudofacmol=1.0
      end if
      do it=1,ntype
         write(10,211) type(it),atsymbol(it),iexchange(it)
         write(10,'(5(1x,f10.6),1x,e13.6,2(1x,f10.4),1x,e13.6)') epscor(it),sigcor(it),amassa(it),charge(it),qradius(it)
     1,natomsphere(it),radminsphere(it),radmaxsphere(it),pseudofacmol(it)
      end do
c write core energy and diameter which represents the closest overlap of ANY 
c non-bonded atom
      write(10,202) ecore,dcore
c
c save random number generator parameters
c
      call saveran1(10)
c
c write out the molecular template files if defined
!In the order in which they occur in the .ato file!
c
      allocate(moltypeswritten(nmoltype))
!      do imoltype=1,nmoltype
!         moltypeswritten(imoltype)=0
!      end do
! In case the relative broadening for each molecule has not been defined, define them
      if(.not.allocated(rfractmol)) then
         call reallocate1d_r(rfractmol,nmoltype)
         rfractmol=0.0
      else if(size(rfractmol,1).lt.nmoltype) then
         call reallocate1d_r(rfractmol,nmoltype)
         rfractmol=0.0
      end if
      nmoltypeswritten=0
      do im=1,nmol
      !Search for a moltype corresponding to this molecule, but if already written, ignore.
         found = .false.
         imoltype=0
!First find the molecule type of this molecule
         do while (imoltype.lt.nmoltype.and..not.found)
            imoltype=imoltype+1
            found=alabel(naf(im)).eq.moltype(imoltype)
         end do
         !Now search the list of written molecules
         found=.false.
         imoltypewritten=0
         do while(imoltypewritten.lt.nmoltypeswritten.and..not.found)
            imoltypewritten=imoltypewritten+1
            found=imoltype.eq.moltypeswritten(imoltypewritten)
         end do
         if(.not.found) then
            nmoltypeswritten=nmoltypeswritten+1
            ilen=index(moltemplate(imoltype),' ')-1
            write(10,212) nmoltypeswritten,moltemplate(imoltype)(1:ilen)
     *,interntrans(imoltype),internrot(imoltype),wholerot(imoltype),wholetrans(imoltype),rfractmol(imoltype)
            moltypeswritten(nmoltypeswritten)=imoltype
         end if
      end do
c
c calculate volume of box and atomic number density
c
      volume=alatx*blaty*clatz
      rho=real(natomr)/volume
      rhocg=natomrcg/volume
      if(natomr.gt.0) then
         do it=1,ntype
            cfrac(it)=cfrac(it)/real(natomr)
            cfraccg(it)=cfraccg(it)/natomrcg
         end do
      end if
!      write(6,915) rho
      write(10,915) rhocg,rho
915   format(' Atomic number density (coarse grain and actual) = ',2(1x,e14.7))
      !Write out the atomic composition in the form required by Gudrun
      do it=1,ntype
!         if(qtypea(it).eq.0) write(6,916) atsymbol(it),0,cfrac(it),type(it)
         if(qtypea(it).eq.0) write(10,916) atsymbol(it),0,cfraccg(it),cfrac(it),type(it)
      end do
916   format(1x,a3,2x,i1,2x,e12.5,10x,e12.5,2x,a3)

      close(10)

! Open the lock file if it exists and set to false (= unlocked)

      open(10,file=fname1(1:len_trim(fname1))//'.lck',status='old',iostat=ierr)
      if(ierr.eq.0) then
         write(10,'(a)') 'F'
         close (10)
      end if

      return
207   format(1x,i7,7(1x,e15.8))
208   format(1x,6(1x,e15.8))
206   format(1x,10(1x,e12.5))
201   format(1x,a3,2(2x,i5))
202   format(8(1x,f14.8))
203   format(1x,i4,(5(1x,i4,1x,e9.3)))
204   format(1x,14(i4,1x))
209   format(1x,2(i4,1x),f10.3)
205   format(1x,i32)
211   format(1x,a3,1x,a3,1x,i2)
212   format(1x,i3,1x,a,4(1x,e13.6),1x,f10.6)
      end


