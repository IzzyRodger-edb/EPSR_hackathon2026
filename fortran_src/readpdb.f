      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_readpdb()

      end      

      subroutine run_readpdb()

c To generate a 3D plot of a .ato file

c 2006-11-15 Output options are (1) GNUplot, (2) PGPlot, (3) ATOMS

c Also atoms to be plotted will be defined by their types, not atomic symbol

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'plot_pars.inc'
      include 'system_commands.inc'
      include 'plotato.inc'
      include 'lattice.inc'

      character*256 readfile
      character*4 ext
      character*1 iaccept
      character*3020 sometextsave
      character*4 yesorno

      integer iread,nwordssave,ncfsave(nworddim),nclsave(nworddim)

c First identify the number of arguments read in

      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)


c generate the name of the .pdb file

      ext='.pdb'
      ilenext=4
      if(nwords.lt.2) then
         nwords=0
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
499   format(/'readpdb> No .pdb file specified! - routine will exit'/)
               return
            endif
         end do
         readfile=sometext(ncf(1):ncl(1))
         ilenreadfile=ncl(1)-ncf(1)+1
      else
         readfile=sometext(ncf(2):ncl(2))
         ilenreadfile=ncl(2)-ncf(2)+1
      endif
c
c strip out any extension that was supplied - extension must be .xyz, otherwise
c nothing is done.
c
      ilendot=indexlastdot(readfile)-1
      if(ilendot.gt.0) then
         readfile=readfile(1:ilendot)
         ilenreadfile=ilendot
      endif
c
c create the full name of the input file by adding on the directory name if not present
c
      call get_full_filename(readfile,ilenreadfile,fnameato,ilenfnameato)
      fnameato=fnameato(1:ilenfnameato)//'.pdb'
      ilenfnameato=ilenfnameato+4
c
c initialise number of atom and molecule types to zero
c
      ntype=0
      nmoltype=0
c
c read .pdb file
c
      sizefactor=1.0
      call readpdb(fnameato,rho,ierr)
      if(rho.eq.0.0.or.ierr.ne.0) then
         write(6,599) fnameato
599   format(/'readpdb> Problems with .pdb file ',a/'Routine will exit now'/)
         return
      endif
c     call normalise_coordinates
c
c step through the atom types and make sure the chemical symbol has been specified
c
      do it=1,ntype
         if(atsymbol(it).eq.'XXX'.or.atsymbol(it).eq.'   ') then
            write(6,115) type(it)
115   format(/'readpdb> Specify the chemical symbol for atom type ',a3,': ',$)
            read(5,*) atsymbol(it)
         endif
      end do

c Ensure the coordinates of each atom in each molecule are relative to the centre of mass

      do im=1,nmol

         call centreofmass(im)

      end do

c Replace the .pdb extension with the .ato extension

      ilendot=indexlastdot(fnameato)-1
      fnameato=fnameato(1:ilendot)//'.ato'
      call moloutc(fnameato)

      return
      end

      subroutine readpdb(readfile,rho,ierr)

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

      INTERFACE
         function upcase(string)
         character(len=*), intent(in) :: string
         character(len=len(string)) :: upcase
         integer :: i,j
         end function
      END INTERFACE

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'plot_pars.inc'
      include 'system_commands.inc'
      include 'potential_par.inc'
      include 'lattice.inc'
      include 'chemical_symbols.inc'

      character*256 readfile
      integer ilenreadfile,ierr,iline,natomexpect
      real vanderWaalsFactor,rho,xdif,ydif,zdif,pi,piconv,alphaaxes,betaaxes,gammaaxes,dtest,dvanderwaals
      real, dimension(melements) :: vanderWaalsRadius
      integer, dimension(:), allocatable :: iatomicno,ialabel
      character*2 chsymbolsave
      character(len=2), dimension(:), allocatable :: chsymboltemp
      character*3 a1,ats,yesorno
      character*4 alabelsave
      character*10 junk
      integer iexch,j,k,ir,ia,ja,ifirst
      integer, dimension(4)                  :: isecond
      external cosd,sind
      logical found,dodihedrals
c temporary rotational groups
      integer                                :: nrott1,nrott2
      integer, dimension(:), allocatable     :: labrott1,labrott2
!Temporary arrays for moving atoms around
      integer, dimension(:), allocatable                 :: nbondmove,imageatommove,ialabelmove
      real, dimension(:), allocatable                    :: xamove,yamove,zamove
      character(len=2), dimension(:), allocatable        :: chsymbolmove
      character(len=3), dimension(:), allocatable        :: alabelmove

      pi=4.0*atan(1.0)
      piconv=180.0/pi
      ilenreadfile=len_trim(readfile)
      open(10,file=readfile,status='old',iostat=ierr)
!      write(6,*) ' readpdb> Open file ',ierr
      if(ierr.ne.0) return !Any errors opening file, return
      iline=0
      natomexpect=0
      sizefactor=1.0
      cubic=.false.
      nmol=0
      simtk=300.0
      tempb=0.65e2
      tempang=1.0
      tempdih=10.0
      phib=90.0
      thetac=0.0
      phic=0.0
      tol=0.0
      tolout=0.0
!Ensure no atom moves are attempted
      stepmi=0.0
      stepri=0.0
      stepr=0.0
      stepm=0.0
      tempangreq=1.0
      tempdihreq=10.0
      ia=0
      im=0
      ig=0
      natom=0
      natomr=0
      ngrrot=0
      nrottot=0
      ntype=0
! Reset array size counters
      ngrrotmax=0
      nrottotmax=0
      imsave=0
      if(allocated(chsymboltemp)) deallocate(chsymboltemp,ialabel)
      allocate(chsymboltemp(matom),ialabel(matom))
      do while (ierr.eq.0.and.ia.lt.matom.and.im.lt.mmol) 
         read(10,'(a)',iostat=ierr) sometext
!         write(6,*) ' readpdb> Read line: ',iline,ierr
         if(ierr.eq.0) then
!            call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
            iline=iline+1
            if(index(sometext(1:6),'CRYST1').gt.0) then
               read(sometext(7:33),'(3f9.3)',iostat=ierr) alat,blat,clat
               read(sometext(34:54),'(3f7.2)',iostat=ierr) alpha,beta,gam
               do while(ierr.ne.0)
                  write(6,*) 'readpdb> Need a,b,c axes and 3 angles (alpha,beta,gamma):'
                  read(5,*,iostat=ierr) alat,blat,clat,alpha,beta,gam
               end do
! Convert angles to internal spherical polar coordinates
               call convertalphabetagamma(alpha,beta,gam,phib,thetac,phic)
c set up lattice parameters for applying the minimum image convention
               alatx=alat
               blatx=blat*cosd(phib)
               blaty=blat*sind(phib)
               clatx=clat*sind(thetac)*cosd(phic)
               claty=clat*sind(thetac)*sind(phic)
               clatz=clat*cosd(thetac)
c define the Cartesian coordinates of the lattice vectors
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
c calculate the perpendicular distance between the faces of the lattice
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
               dmax=min(perpab,perpbc,perpca)
               write(6,*) dmax
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
c
c calculate volume of box
c
               volume=alatx*blaty*clatz
            else if(index(sometext(1:6),'ATOM').gt.0) then
               !Otherwise we read the list of atomic coordinates
               read(sometext(7:11),*,iostat=ierr) ia
               read(sometext(13:14),*,iostat=ierr) chsymbolsave
               read(sometext(13:16),*,iostat=ierr) alabelsave
               read(sometext(23:26),*,iostat=ierr) im
               read(sometext(31:54),'(3f8.3)',iostat=ierr) xat,yat,zat
               write(6,'(i5,a2,i5,a4)') ia,chsymbolsave,im,alabelsave
               if(ierr.eq.0) then
                  natom=natom+1
                  nbond(natom)=0
                  imageatom(natom)=0
                  ialabel(natom)=ia !The value of ia read in is kept as a (temporary) label to identify this atom uniquely
                  chsymboltemp(natom)=chsymbolsave
                  alabel(natom)=alabelsave
                  xa(natom)=xat
                  ya(natom)=yat
                  za(natom)=zat
c
c do not count atoms which begin with 'Q' or 'q' - these are assumed 
c invisible to the scattering radiation and experiment (although they
c can be used to calculate other structural quantities if necessary
c
                  if(im.ne.imsave) then
                     nmol=nmol+1
                     naf(nmol)=natom
                     xadd(nmol)=0.0
                     yadd(nmol)=0.0
                     zadd(nmol)=0.0
                     nato(nmol)=1
                     ngrot(nmol)=0
                     ngrrot=ngrrot+ngrot(nmol)
                     tethered(nmol)=.false.
                     phix(nmol)=0.0
                     phiy(nmol)=0.0
                     phiz(nmol)=0.0
                  else
                     nal(nmol)=natom
                     nato(nmol)=nato(nmol)+1
                  end if
                  mref(natom)=nmol
                  imsave=im
               end if
            else if(index(sometext(1:6),'CONECT').gt.0.and.len_trim(sometext).gt.15) then !It will be assumed any atom numbers read in here have already been defined in a
!corresponding ATOM line. Otherwise this line will be ignored.
               !Determine the number of numbers on this line
               nsecond=int((len_trim(sometext)-6)/5)-1
               read(sometext(7:11),*) ifirst
               !Get the internal atom number associated with this label
               found=.false.
               ia=0
               do while(.not.found.and.ia.lt.natom)
                  ia=ia+1
                  found=ialabel(ia).eq.ifirst
               end do
               if(found) ifirst=ia
               ir1=12
               do i=1,nsecond
                  ir2=ir1+4
                  read(sometext(ir1:ir2),*) isecond(i)
                  !Get the internal atom number associated with this label
                  found=.false.
                  ia=0
                  do while(.not.found.and.ia.lt.natom)
                     ia=ia+1
                     found=ialabel(ia).eq.isecond(i)
                  end do
                  if(found) isecond(i)=ia
                  ir1=ir2+1
               end do
               write(6,*) ' readpdb> CONECT: ',ifirst,mref(ifirst),(isecond(i),mref(isecond(i)),i=1,nsecond)
               i=1
               do while(i.le.nsecond.and.ifirst.le.natom.and.isecond(i).le.natom)
                  if(mref(ifirst).ne.mref(isecond(i))) then 
                     !Need to merge molecules
                     immin=min(mref(ifirst),mref(isecond(i)))
                     immax=max(mref(ifirst),mref(isecond(i)))
                     !Put all the atoms in immax into immin
                     !Check whether atom numbers in immax are contiguous with those in immin
!               write(6,*) ' readpdb> CONECT: ',ifirst,mref(ifirst),isecond(i),mref(isecond(i)),immin,immax,nal(immin),naf(immax)
                     if(naf(immax).ne.(nal(immin)+1)) then
                        !If not contiguous, need to reorder atoms in list between immin and immax
                        !First save the data in immax
                        if(allocated(nbondmove)) deallocate(nbondmove,imageatommove,ialabelmove,chsymbolmove,alabelmove,xamove,yamove,zamove)
                        allocate(nbondmove(nato(immax)),imageatommove(nato(immax)),ialabelmove(nato(immax)),chsymbolmove(nato(immax))
     1,alabelmove(nato(immax)),xamove(nato(immax)),yamove(nato(immax)),zamove(nato(immax)))
                        do ia=1,nato(immax)
                           iaref=ia+naf(immax)-1
                           write(6,*) immax,ia,iaref,ialabel(iaref)
                           nbondmove(ia)=nbond(iaref)
                           imageatommove(ia)=imageatom(iaref)
                           ialabelmove(ia)=ialabel(iaref)
                           chsymbolmove(ia)=chsymboltemp(iaref)
                           alabelmove(ia)=alabel(iaref)
                           xamove(ia)=xa(iaref)
                           yamove(ia)=ya(iaref)
                           zamove(ia)=za(iaref)
                        end do
                        !Now move atoms, starting from the one above naf(immax) and going backwards
                        ia=naf(immax)-1
                        mrefsave=mref(ia)
                        nal(mrefsave)=ia+nato(immax)
                        do while(ia.gt.nal(immin))
                           iaref=ia+nato(immax)
                           write(6,*) ia,iaref,ialabel(ia)
                           nbond(iaref)=nbond(ia)
                           imageatom(iaref)=imageatom(ia)
                           ialabel(iaref)=ialabel(ia)
                           chsymboltemp(iaref)=chsymboltemp(ia)
                           alabel(iaref)=alabel(ia)
                           xa(iaref)=xa(ia)
                           ya(iaref)=ya(ia)
                           za(iaref)=za(ia)
                           if(mref(ia).ne.mrefsave) then !Remember to reset the integers defining the atom numbers for each molecule.
                              mrefsave=mref(ia)
                              nal(mrefsave)=iaref
                           else
                              naf(mrefsave)=iaref
                           end if
                           mref(iaref)=mrefsave !The molecule numbers for this region do not change
                           ia=ia-1
                        end do
                        !Finally insert the molecule in the space created
                        do ia=1,nato(immax)
                           iaref=ia+nal(immin)
                           write(6,*) ia,iaref,ialabelmove(ia)
                           nbond(iaref)=nbondmove(ia)
                           imageatom(iaref)=imageatommove(ia)
                           ialabel(iaref)=ialabelmove(ia)
                           chsymboltemp(iaref)=chsymbolmove(ia)
                           alabel(iaref)=alabelmove(ia)
                           xa(iaref)=xamove(ia)
                           ya(iaref)=yamove(ia)
                           za(iaref)=zamove(ia)
                           mref(iaref)=immin
                        end do                           
                        deallocate(nbondmove,imageatommove,ialabelmove,chsymbolmove,alabelmove,xamove,yamove,zamove)
                     end if
!               write(6,*) ' readpdb> CONECT: ',immin,immax,naf(immin),nal(immin),nato(immin),naf(immax),nal(immax),nato(immax)
                     nato(immin)=nato(immin)+nato(immax)
                     nal(immin)=nal(immin)+nato(immax)
!               write(6,*) ' readpdb> CONECT: ',immin,immax,naf(immin),nal(immin),nato(immin),naf(immax),nal(immax),nato(immax)
                     !Since we haven't defined any rotational or dihedral groups at this stage we don't need to worry about them
                     do ia=naf(immin),nal(immin)
                        mref(ia)=immin
                     end do
                     !Need to remove 1 molecule from system if immax is less than nmol
                     im=immax
                     do while(im.lt.nmol)
                        im1=im+1
                        nato(im)=nato(im1)
                        naf(im)=naf(im1)
                        nal(im)=nal(im1)
                        do ia=naf(im),nal(im)
                           mref(ia)=im
                        end do
                        im=im+1
                        !No need to worry about molecular coordinates, etc., at this stage as none have been defined.
                     end do
                     nmol=nmol-1
                     do im=1,immax
                        write(6,*) im,naf(im),nal(im),mref(naf(im)),mref(nal(im))
                     end do
                  end if
                  i=i+1
               end do
            end if
         end if
      end do
      close(10)
      if(natom.gt.0) ierr=0 !Reset the error indicator so that calling program does not think read failed.
      do im=1,nmol
         xsum=0.0
         ysum=0.0
         zsum=0.0
         do ib=naf(im),nal(im)
            xsum=xsum+xa(ib)
            ysum=ysum+ya(ib)
            zsum=zsum+za(ib)
         end do
         xadd(im)=xsum/nato(im)
         yadd(im)=ysum/nato(im)
         zadd(im)=zsum/nato(im)
         do ib=naf(im),nal(im)
            xa(ib)=xa(ib)-xadd(im)
            ya(ib)=ya(ib)-yadd(im)
            za(ib)=za(ib)-zadd(im)
         end do
!         write(6,*) im,alabel(naf(im)),naf(im),nal(im),nato(im),xadd(im),yadd(im),zadd(im)
      end do
!      do ia=1,natom
!         write(6,*) natom,ia,alabel(ia),xa(ia),ya(ia),za(ia)
!      end do
!Set up the atomic numbers of all the atoms
      if(allocated(iatomicno)) deallocate(iatomicno)
      allocate(iatomicno(natom))
      natomr=0
      do ia=1,natom
         ir=0
         found=.false.
         do while(.not.found.and.ir.lt.melements)
            ir=ir+1
            !It is assumed the atom labels are also the chemical symbols
!            write(6,'(i5,a2,i3,a2)') ia,chsymboltemp(ia),ir,upcase(chsymbol(ir))
            found=chsymboltemp(ia).eq.upcase(chsymbol(ir)) !For PDB files the first 1 or 2 characters of the atom name are the atomic symbol in upper case
            if(found) then
               iatomicno(ia)=ir
            end if
         end do
         if(.not.found) then
! If no valid atomic symbol found. This implies this atom is a dummy, non-scattering atom, so we convert it a Q atom of type oxygen
            iatomicno(ia)=8
            alabel(ia)(1:1)='Q'
         end if
c set up 'type' of atom - this is used in assigning neutron weights
         call ftype(ia)
         if(qtype(ia).ne.1) then
            natomr=natomr+1
         endif
      end do
      rho=natomr/volume
      write(6,*) ' readpdb> rho = ',rho,natomr,natom,nmol
      do it=1,ntype
         write(6,*) it,type(it)
      end do
      if(allocated(chsymboltemp)) deallocate(chsymboltemp)
      
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
            iatomicnotype(ntypefile)=iatomicno(ia)
            write(6,*) ntypefile,typefile(ntypefile),iatomicnotype(ntypefile)
         endif
      end do
!Use these values to set the potential parameters
      do it=1,ntypefile
         if(chsymbol(iatomicnotype(it)).eq.'H ') then
            epsilon=0
            sigma=0
            amass=2
         else if(chsymbol(iatomicnotype(it)).eq.'O ') then
            epsilon=0.65
            sigma=3.1
            amass=chatwt(iatomicnotype(it))
         else if(chsymbol(iatomicnotype(it)).eq.'N ') then
            epsilon=0.7
            sigma=3.2
            amass=chatwt(iatomicnotype(it))
         else if(chsymbol(iatomicnotype(it)).eq.'C ') then
            epsilon=0.8
            sigma=3.7
            amass=chatwt(iatomicnotype(it))
         else
            epsilon=0.8
            sigma=3.2
            amass=chatwt(iatomicnotype(it))
         end if
         call ctype(typefile(it),chsymbol(iatomicnotype(it))
     *,0,epsilon,sigma,amass,1.0,0.0,0.5*sigma,1.0,0.0)
      end do  
c
c make a list of the molecule types read in from the current file
c
      im=0
      nmoltype=0
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
c
c initialise the template filename
c
            write(moltemplate(nmoltype),'("moltype",i2.2)') nmoltype
            mtype(im)=nmoltype

c Initialise the rotational and translational step sizes

            interntrans(nmoltype)=stepmi
            internrot(nmoltype)=stepri
            wholerot(nmoltype)=stepr
            wholetrans(nmoltype)=stepm

         endif
      end do

c Write these out as a check

      write(6,201)
201   format(/' readpdb> Following molecule types defined initially:-'/)

      do imoltype=1,nmoltype

         write(6,202) imoltype,moltype(imoltype),moltemplate(imoltype)(1:10),wholerot(imoltype)
     *,wholetrans(imoltype)
      end do
202   format(1x,i3,1x,a3,1x,a10,2(1x,e12.5))

! Now setup the bonding parameters within each molecule
! Get the list of van der Waals radii
      open(10,file='vanderWaalsRadii.txt',status='old',iostat=ierr)
      if(ierr.ne.0) then
         write(6,*) ' readpdb> Unable to open vanderWaalsRadii.txt'
         return
      end if
      read(10,'(a)',iostat=ierr) sometext
      icount=0
      do while (ierr.eq.0)
         call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
         if(nwords.gt.2) then !Ignore lines with insufficient data
            !The first word is the atomic number.
            !The second is the atomic symbol - ignored in this case.
            !The third is the van der Waals radius. 
            read(sometext(ncf(1):ncl(1)),*) ir
            if(ir.gt.0.and.ir.le.melements) then
               read(sometext(ncf(3):ncl(3)),*) vanderWaalsRadius(ir)
               if(ir.eq.1) vanderWaalsRadius(ir)=0.0 !The van der Waals radius for H is set to zero.
               icount=icount+1
            end if
         end if
         read(10,'(a)',iostat=ierr) sometext
      end do
      close(10)
      if(icount.gt.0) ierr=0
      write(6,48,advance='no') icount
48    format(/' readpdb> Number of van der Waals radii read in: ',i3
     *//' Give van der Waals factor (suggest 0.7). This is used to test whether two atoms are bonded or not: ')
      read(5,*,iostat=ierr) vanderWaalsFactor
      if(ierr.ne.0.or.vanderWaalsFactor.le.0.0) then
         vanderwaalsfactor=0.7
      end if
      write(6,49,advance='no')
49    format(/' readpdb> Create dihedral angles (y or no): ')
      read(5,'(a)') yesorno
      dodihedrals=index(yesorno,'y').gt.0
      call normalise_coordinates()
      do im=1,nmol
         do ia=naf(im),nal(im)-1
            do ja=ia+1,nal(im)
               xdif=xa(ia)-xa(ja)
               ydif=ya(ia)-ya(ja)
               zdif=za(ia)-za(ja)
               dtest=dmax*sqrt(xdif*xdif+ydif*ydif+zdif*zdif)
               dvanderwaals=vanderWaalsFactor*(vanderWaalsRadius(iatomicnotype(itype(ia)))+vanderWaalsRadius(iatomicnotype(itype(ja))))
               found=dtest.lt.dvanderwaals
!               write(6,*) im,ia,ja,dtest,dvanderwaals,found
               if(found) then
                  if(nbond(ia).lt.near.and.nbond(ja).lt.near) then
                     nbond(ia)=nbond(ia)+1
                     iabond(ia,nbond(ia))=ja
                     tbond(ia,nbond(ia))=dtest/dmax
                     nbond(ja)=nbond(ja)+1
                     iabond(ja,nbond(ja))=ia
                     tbond(ja,nbond(ja))=dtest/dmax
                  end if
               end if
            end do
         end do
      end do
      !Now check for bond angles within each molecule
      do ia=1,natom
         nbondp(ia)=nbond(ia)
      end do
      ir=0
      do im=1,nmol
         do ia=naf(im),nal(im)
            ib1=0
            do while(ib1.lt.nbondp(ia)-1) !Don't want to include any new bond angles that are already defined
               ib1=ib1+1
               iaref1=iabond(ia,ib1)
               ib2=ib1
               do while(ib2.lt.nbondp(ia))
                  ib2=ib2+1
                  iaref2=iabond(ia,ib2)
                  !Need to check that iabond(ia,ib2) does not already exist in iaref1's neighbour list
                  ic=0
                  found=.false.
                  do while(.not.found.and.ic.lt.nbond(iaref1))
                     ic=ic+1
                     found=iaref2.eq.iabond(iaref1,ic)
                  end do
!      write(6,*) im,ia,ib1,iaref1,ib2,iaref2,ic,iabond(iaref1,ic),found
                  if(.not.found) then
                     xdif=xa(iaref1)-xa(iaref2)
                     ydif=ya(iaref1)-ya(iaref2)
                     zdif=za(iaref1)-za(iaref2)
                     dtest=sqrt(xdif*xdif+ydif*ydif+zdif*zdif)
                     if(nbond(iaref1).lt.near.and.nbond(iaref2).lt.near) then
                        nbond(iaref1)=nbond(iaref1)+1
                        nbond(iaref2)=nbond(iaref2)+1
                        iabond(iaref1,nbond(iaref1))=iaref2
                        iabond(iaref2,nbond(iaref2))=iaref1
                        tbond(iaref1,nbond(iaref1))=-dtest
                        tbond(iaref2,nbond(iaref2))=-dtest
                     end if
                  end if
               end do
            end do
         end do
         !Dihedral angles do...
         if(dodihedrals.and.nato(im).gt.3) then !Can only do this if molecule contains at least 4 atoms
         !Also dihedrals can only be created between atoms that do not have bonds
            do iaref1=naf(im),nal(im) !ia1ref will be the first axis atom
               ib1=0
               do while(ib1.lt.nbondp(iaref1)) !Look at the first neighbours of iaref1 to determine the second atom to form the rotation axis (iaref2)
                  ib1=ib1+1
                  iaref2=iabond(iaref1,ib1)
                  ! Need to check whether a dihedral group for this pair of axis atoms has already been defined.
                  ir1=iaref1
                  ir2=iaref2
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
                  ib2=0
                  do while(ib2.lt.nbondp(iaref1)) !Look at other first neighbours of iaref1 to determine the second axis atom (iaref2)
                     ib2=ib2+1
                     ibref1=iabond(iaref1,ib2) !(ibref1 must be a direct bond to iaref1, not via a bond angle, hence us nbondp and not nbond)
                     ib3=0
                     do while(ib3.lt.nbondp(iaref2).and.ibref1.ne.iaref2) !Look at the first neighbours of iaref2 to determine second atom in the dihedral (ibref2)
                        ib3=ib3+1
                        ibref2=iabond(iaref2,ib3)
                        ic=0
                        found=ibref1.eq.ibref2.or.ibref2.eq.iaref1
                        do while(ic.lt.nbond(ibref2).and..not.found) !Need to check that ibref2 has no bonds with ibref1
                           ic=ic+1
                           found=ibref1.eq.iabond(ibref2,ic)
                        end do
                        if(ifounddih.gt.0.and..not.found) then
                           !Need to check whether a dihedral group with these second atoms already exists
                           iref=0
                           do while(iref.lt.ndihedralang(ifounddih).and..not.found)
                              iref=iref+1
                              found=(ibref1.eq.labrotdih(1,iref,ifounddih).and.ibref2.eq.labrotdih(2,iref,ifounddih))
     1.or.(ibref2.eq.labrotdih(1,iref,ifounddih).and.ibref1.eq.labrotdih(2,iref,ifounddih))
                           end do
                        endif
                        if(.not.found) then
                           dread=dihedral(im,iaref1,iaref2,ibref1,ibref2)
                           write(6,*) ifounddih,iaref1,iaref2,ibref1,ibref2,dread
                           dread=dread/piconv!Need to convert angles to radians as these are normalised coordinates
                           if(ifounddih.eq.0) then
                              ! Start a new rotational group:
                              ! Get the two possible rotational groups into temporary arrays
                              call get_rotgrp(ir1,ir2,labrott1,nrott1)
                              call get_rotgrp(ir2,ir1,labrott2,nrott2)
                              write(6,'(i5/(10(1x,i6)))') nrott1,(labrott1(j),j=1,nrott1)
                              write(6,'(i5/(10(1x,i6)))') nrott2,(labrott2(j),j=1,nrott2)
                              ir=ir+1
                              ndihedralang(ir)=1
                              dihedralang(ndihedralang(ir),ir)=dread
                              call testsize_rotgrps(ir)
                              rotgrp(ir)='DIH'
                              ifounddih=ir !Need to record the fact that this is the current found dihedral group
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
                                 labrotdih(1,ndihedralang(ir),ir)=ibref1
                                 labrotdih(2,ndihedralang(ir),ir)=ibref2
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
                                 labrotdih(1,ndihedralang(ir),ir)=ibref2
                                 labrotdih(2,ndihedralang(ir),ir)=ibref1
                              end if
                           else
                              ! Add the new dihedral angle to the existing list
                              if(ndihedralang(ifounddih).lt.near) then
                                 ndihedralang(ifounddih)=ndihedralang(ifounddih)+1
                                 dihedralang(ndihedralang(ifounddih),ifounddih)=dread
                                 if(ir1.eq.labaxi(ifounddih,1).and.ir2.eq.labaxi(ifounddih,2)) then
                                    labrotdih(1,ndihedralang(ifounddih),ifounddih)=ibref1
                                    labrotdih(2,ndihedralang(ifounddih),ifounddih)=ibref2
                                 else if((ir2.eq.labaxi(ifounddih,1).and.ir1.eq.labaxi(ifounddih,2))) then
                                    labrotdih(1,ndihedralang(ifounddih),ifounddih)=ibref2
                                    labrotdih(2,ndihedralang(ifounddih),ifounddih)=ibref1
                                 end if
!                                 write(6,*) labaxi(ifounddih,1),labaxi(ifounddih,2)
!                                 do iref=1,ndihedralang(ifounddih)
!                                    write(6,*) labrotdih(1,iref,ifounddih),labrotdih(2,iref,ifounddih),dihedralang(iref,ifounddih)*piconv
!                                 end do
                              end if
                           end if
                        end if
                     end do
                  end do
               end do
            end do
         end if
         ngrot(im)=ir-ngrrot
         ngrrot=ir
      end do
!      do ia=1,natom
!         nbond(ia)=nbondp(ia)
!      end do
      call unnormalise_coordinates()

      end subroutine readpdb

      function upcase(string)
      character(len=*), intent(in) :: string
      character(len=len(string)) :: upcase
      integer :: i,j
      do i = 1,len(string)
          j = iachar(string(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               upcase(i:i) = achar(iachar(string(i:i))-32)
          else
               upcase(i:i) = string(i:i)
          end if
      end do
      end function upcase

   


