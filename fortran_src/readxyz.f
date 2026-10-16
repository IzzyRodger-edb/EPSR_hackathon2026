      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_readxyz()

      end      

      subroutine run_readxyz()

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


c generate the name of the .xyz file

      ext='.xyz'
      ilenext=4
      if(nwords.lt.2) then
         nwords=0
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
499   format(/'readxyz> No .xyz file specified! - routine will exit'/)
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
      fnameato=fnameato(1:ilenfnameato)//'.xyz'
      ilenfnameato=ilenfnameato+4
c
c initialise number of atom and molecule types to zero
c
      ntype=0
      nmoltype=0
c
c read .xyz file
c
      sizefactor=1.0
      call readxyz(fnameato,rho,ierr)
      if(rho.eq.0.0.or.ierr.ne.0) then
         write(6,599) fnameato
599   format(/'readxyz> Problems with .xyz file ',a/'Routine will exit now'/)
         return
      endif
c     call normalise_coordinates
c
c step through the atom types and make sure the chemical symbol has been specified
c
      do it=1,ntype
         if(atsymbol(it).eq.'XXX'.or.atsymbol(it).eq.'   ') then
            write(6,115) type(it)
115   format(/'readxyz> Specify the chemical symbol for atom type ',a3,': ',$)
            read(5,*) atsymbol(it)
         endif
      end do

c Ensure the coordinates of each atom in each molecule are relative to the centre of mass

      do im=1,nmol

         call centreofmass(im)

      end do

c Replace the .xyz extension with the .ato extension

      ilendot=indexlastdot(fnameato)-1
      fnameato=fnameato(1:ilendot)//'.ato'
      call moloutc(fnameato)

      return
      end

      subroutine readxyz(readfile,rho,ierr)

      use molecom_mod
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
      real vanderWaalsFactor,rho,xdif,ydif,zdif,pi,piconv,alphaaxes,betaaxes,gammaaxes
      real, dimension(melements) :: vanderWaalsRadius
      integer, dimension(:), allocatable :: iatomicno
      character*3 a1,ats
      integer iexch,j,k,ir,ia,ja
      external cosd,sind
      logical found

      ilenreadfile=len_trim(readfile)
      open(10,file=readfile,status='old',iostat=ierr)
!      write(6,*) 'readxyz> Open file ',ierr
      if(ierr.ne.0) return !Any errors opening file, return
      iline=0
      natomexpect=0
      sizefactor=1.0
      cubic=.false.
      nmol=0
      simtk=300.0
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
      do while (ierr.eq.0.and.ia.lt.matom.and.im.lt.mmol) 
         read(10,'(a)',iostat=ierr) sometext
!         write(6,*) 'readxyz> Read line: ',iline,ierr
         if(ierr.eq.0) then
            call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
            iline=iline+1
            if(iline.eq.1) then
               read(sometext,*) natomexpect
            else if(iline.eq.2) then
               if(nwords.lt.7) then
                  ierr=1
               else ! Get the crystallographic axes and angles
                  read(sometext(ncf(2):ncl(7)),*,iostat=ierr) alat,blat,clat,alpha,beta,gam
               endif
               do while(ierr.ne.0)
                  write(6,*) 'readxyz> Need a,b,c axes and 3 angles (alpha,beta,gamma):'
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
c
c calculate volume of box
c
               volume=alatx*blaty*clatz
            else !Otherwise we read the list of atomic coordinates, creating single atom molecules
               if(nwords.lt.4) then
                  ierr=1
               else
                  ia=ia+1
                  im=im+1
                  read(sometext(ncf(1):ncl(4)),*,iostat=ierr) alabel(ia),xadd(im),yadd(ia),zadd(im)
                  imageatom(ia)=0
                  xa(ia)=0.0
                  ya(ia)=0.0
                  za(ia)=0.0
c set up 'type' of atom - this is used in assigning neutron weights
                  call ftype(ia)
c
c do not count atoms which begin with 'Q' or 'q' - these are assumed 
c invisible to the scattering radiation and experiment (although they
c can be used to calculate other structural quantities if necessary
c
                  if(qtype(ia).ne.1) then
                        natomr=natomr+1
                  endif
                  natom=natom+1
                  nbond(ia)=0
                  naf(im)=ia
                  nal(im)=ia
                  nato(im)=1
                  ngrot(im)=0
                  ngrrot=ngrrot+ngrot(im)
                  tethered(im)=.false.
                  phix(im)=0.0
                  phiy(im)=0.0
                  phiz(im)=0.0
               end if
            end if
         end if
      end do
      close(10)
      nmol=im
!      do ia=1,natom
!         write(6,*) natomexpect,ia,alabel(ia),xa(ia),ya(ia),za(ia)
!      end do
      rho=natomr/volume
      write(6,*) 'readxyz> rho = ',rho,natomr,natom,natomexpect
      if(natom.ge.natomexpect) ierr=0
!Set up the atomic numbers of all the atoms
      if(allocated(iatomicno)) deallocate(iatomicno)
      allocate(iatomicno(natom))
      do ia=1,natom
         ir=0
         found=.false.
         do while(.not.found.and.ir.lt.melements)
            ir=ir+1
            !It is assumed the atom labels are also the chemical symbols
            found=alabel(ia).eq.chsymbol(ir)
            if(found) then
               iatomicno(ia)=ir
            end if
         end do
         if(.not.found) then !If a valid symbol is not found, then the file cannot be processed
            ierr=1
            write(6,*) 'readxyz> Cannot read this file: atom label not recognised.'
            return
         end if
      end do
! Use changeato to sort these atoms into molecules
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
201   format(/'readxyz> Following molecule types defined initially:-'/)

      do imoltype=1,nmoltype

         write(6,202) imoltype,moltype(imoltype),moltemplate(imoltype)(1:10),wholerot(imoltype)
     *,wholetrans(imoltype)
202   format(1x,i3,1x,a3,1x,a10,2(1x,e12.5))

      end do
      end subroutine readxyz


