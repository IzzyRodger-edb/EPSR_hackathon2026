      module dockato_arrays

      integer natos,ngrots,natomax,ngrotmax,nrotmax,ndihmax,nadd

      integer, dimension(:), allocatable           :: nbondss,ndihedralangs!(near,mmix)!(ma,mmix)
      integer, dimension(:,:), allocatable         :: iabonds,labaxis!(ma,near,mmix)(near,2,mmix)(near,near,mmix)
      integer, dimension(:,:,:), allocatable       :: labrotdihs!(2,near,near,mmix)

      real                                           :: xadds,yadds,zadds,amassts,phixs,phiys,phizs
      real, dimension(:), allocatable              :: xas,yas,zas!(ma,mmix)
      real, dimension(:,:), allocatable            :: tbonds,dihedralangs!(ma,near,mmix)
      real                                           :: overlapfactor ! Used to control the extent to which overlaps are allowed. 0.0 means all overlaps are allowed.

c temporary rotational groups
      integer                                      :: nrott1,nrott2
      integer, dimension(:), allocatable           :: labrott1,labrott2

      character(len=3), dimension(:), allocatable  :: alabels,rotgrps!(ma,mmix)(near,mmix)

      end module dockato_arrays

      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call dockato()

      end      

      subroutine dockato()
      
      use reallocation_routines
      use molecom_mod
      use dockato_arrays

      INTERFACE 
         subroutine get_rotgrp(ir1,ir2,labrott,nrott)
            integer :: ir1,ir2,nrott
            integer, dimension(:), allocatable, intent(out)  :: labrott
         end subroutine get_rotgrp
      END INTERFACE

c A specified number of the first (docking) molecule .ato file are docked onto the 
c specified site on the  second (docked) molecule .ato file.

c 9-Oct-2006. Modified so that when docking occurs, the added molecule becomes physically
c bonded to the docked molecule. This version of the program is quite different from 
c previous versions. Thus far it will overlap and rotate the molecules as required, but it will
c not bond them, so it cannot be used for more than one pair of molecules, since the bonding sites
c have not been eliminated at this stage. The next version will include the bonding operation.

c The bonding sites are designated by Qv atoms, which are removed once docking
c has occurred. The number of neighbouring Qv atoms is used to designate the degree of
c saturation of each site.

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      character*256 readfile
      character*4 ext

      character*20 nooverlapstext
      logical nooverlaps !If true docking molecule will be checked for overlaps with docked molecule.

c names of docking atoms

      character*3 atomv,atomd,atomvr,atomdr,atomvrr,atomdrr
      character*1 iaccept

      integer ndock
      integer iad,iadr

      real dminrr,dmaxrr ! minimum and maximum distance to determine if a second rotation is possible

c Variables to save filenames and atom labels 

      character*80 dockingfile,dockedfile,outputfile

      integer ilendockingfile,ilendockedfile,ilenatomd,ilenatomv
     *,ilenoutputfile,ilenatomdr,ilenatomvr,ilenatomvrr,ilenatomdrr

      logical fixeddihedral     !Use to indicate whether a specific dihedral has been set. If random
!then no specific intramolecular distance is set to control the dihedral angle.

c Keep the compiler happy

      ilendockingfile=0
      ilendockedfile=0
      ilenatomd=0
      ilenatomdr=0
      ilenatomdrr=0
      ilenatomv=0
      ilenatomvr=0
      ilenatomvrr=0
      ilendminrr=0
      ilendmaxrr=0
      ilenoutputfile=0
      ndock=0
      radmax=0.0
      dminrr=0.0
      dmaxrr=0.0
      xadds=0
      yadds=0
      zadds=0
      phixs=0
      phiys=0
      phizs=0
      natos=0
      amassts=0
      nooverlapstext='<undefined>'

      ntest=1
c get name of docking file - the only first molecule in this file will be used for docking
      if(nwords.gt.ntest) then
         nt=ntest+1
         dockingfile=sometext(ncf(nt):ncl(nt))
         ilendockingfile=ncl(nt)-ncf(nt)+1
      endif
c get name of atom to be used for translating the docking molecule
      ntest=ntest+1
      if(nwords.gt.ntest) then
         nt=ntest+1
         atomd=sometext(ncf(nt):ncl(nt))
         ilenatomd=ncl(nt)-ncf(nt)+1
      endif
c get name of atom to be used for rotating the docking molecule
      ntest=ntest+1
      if(nwords.gt.ntest) then
         nt=ntest+1
         atomdr=sometext(ncf(nt):ncl(nt))
         ilenatomdr=ncl(nt)-ncf(nt)+1
      endif
c get name of second atom to be used to generate a second rotation if desired
      ntest=ntest+1
      if(nwords.gt.ntest) then
         nt=ntest+1
         atomdrr=sometext(ncf(nt):ncl(nt))
         ilenatomdrr=ncl(nt)-ncf(nt)+1
      endif
c get number of molecules to be docked
      ntest=ntest+1
      if(nwords.gt.ntest) then
         nt=ntest+1
         read(sometext(ncf(nt):ncl(nt)),*) ndock
      endif 
c get maximum distance of these molecules from first molecule
      ntest=ntest+1
      if(nwords.gt.ntest) then
         nt=ntest+1
         read(sometext(ncf(nt):ncl(nt)),*) radmax
      endif 
c get name of .ato file to be docked onto.
      ntest=ntest+1
      if(nwords.gt.ntest) then
         nt=ntest+1
         dockedfile=sometext(ncf(nt):ncl(nt))
         ilendockedfile=ncl(nt)-ncf(nt)+1
      endif
c get name of atom to be translated onto
      ntest=ntest+1
      if(nwords.gt.ntest) then
         nt=ntest+1
         atomv=sometext(ncf(nt):ncl(nt))
         ilenatomv=ncl(nt)-ncf(nt)+1
      endif
c get name of atom to be rotated onto
      ntest=ntest+1
      if(nwords.gt.ntest) then
         nt=ntest+1
         atomvr=sometext(ncf(nt):ncl(nt))
         ilenatomvr=ncl(nt)-ncf(nt)+1
      endif
c get name of atom on the system being docked onto to be used to generate
c a second rotation if desired.
      ntest=ntest+1
      if(nwords.gt.ntest) then
         nt=ntest+1
         atomvrr=sometext(ncf(nt):ncl(nt))
         ilenatomvrr=ncl(nt)-ncf(nt)+1
      endif
c get minimum dihedral angle
      ntest=ntest+1
      if(nwords.gt.ntest) then
         nt=ntest+1
         read(sometext(ncf(nt):ncl(nt)),*) dminrr
         ilendminrr=ncl(nt)-ncf(nt)+1
      endif
c get maximum dihedral angle
      ntest=ntest+1
      if(nwords.gt.ntest) then
         nt=ntest+1
         read(sometext(ncf(nt):ncl(nt)),*) dmaxrr
         ilendmaxrr=ncl(nt)-ncf(nt)+1
      endif
c Determine if docking molecule will be checked for overlaps with docked molecule
      ntest=ntest+1
      nooverlaps=.false.
      if(nwords.gt.ntest) then
         nt=ntest+1
         nooverlapstext=sometext(ncf(nt):ncl(nt))
      endif 
c get name of output file
      ntest=ntest+1
      if(nwords.gt.ntest) then
         nt=ntest+1
         outputfile=sometext(ncf(nt):ncl(nt))
         ilenoutputfile=ncl(nt)-ncf(nt)+1
      endif 

      ierr=1
      nmix=2
      natomax=0
      ngrotmax=0
      nrotmax=0
      ndihmax=0
      natos=0
      ngrots=0

c set the number of atom and molecule types to zero at the outset

      ntype=0
      nmoltype=0
      ntypefirst=0
      nmoltypefirst=nmoltype
      imix=0
      do while(imix.lt.nmix)
         imix=imix+1
         if(imix.eq.1) then
            if(ilendockingfile.eq.0) then
               nwords=0
               write(6,699)
699   format(//'dockato> Search for .ato file containing docking molecule')
            else
               nwords=1
               ncf(1)=1
               ncl(1)=ilendockingfile
               sometext(1:ilendockingfile)=dockingfile
            endif
         else
            if(ilendockedfile.eq.0) then
               nwords=0
               write(6,698)
698   format(//'dockato> Search for .ato file to be docked onto')
            else
               nwords=1
               ncf(1)=1
               ncl(1)=ilendockedfile
               sometext(1:ilendockedfile)=dockedfile
            endif
         endif

c generate the name of each .ato file in turn

         ext='.ato'
         ilenext=4
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)

c if no .ato files are available in this directory, abandon this attempt at dockato

            if(iaccept.ne.'y') then
               write(6,700)
700   format(/'dockato> No .ato files found in the current directory - try another directory.'/)
               return
            endif
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
499   format(/'dockato> No .ato file specified - try again!'/)
            endif
         end do
         readfile=sometext(ncf(1):ncl(1))
         ilenreadfile=ncl(1)-ncf(1)+1

c strip out any extension that was supplied - extension must be .ato, otherwise
c nothing is done.

         ilendot=index(readfile,'.')-1
         if(ilendot.gt.0) then
            readfile=readfile(1:ilendot)
            ilenreadfile=ilendot
         endif
c
c create the full name of the input file by adding on the directory name if not present
c
         call get_full_filename(readfile,ilenreadfile,fnameato,ilenfnameato)
         fnameato=fnameato(1:ilenfnameato)//'.ato'
         ilenfnameato=ilenfnameato+4

c read .ato file

         nmoltype=0
         nmoltypefirst=nmoltype
         sizefactor=1.0
         call minitc(fnameato,rmaxcell,cellmax,upol,rho,amubar,amubar2
     *,ndist,ierr)
         if(rho.gt.0.0.and.ierr.eq.0) then
            ntypefile=ntype-ntypefirst
            nmoltypefirst=nmoltype
            write(6,36) ntypefile,fnameato(1:ilenfnameato)
36    format(/1x,i3,' new atom types in file ',a/)
            if(ntypefile.gt.0) then
            do it=ntypefirst+1,ntype
                     write(6,33) it,type(it)
                  end do
         end if
33    format(1x,'Atom type ',i2,' has label ',a3)
            ntypefirst=ntype

c
c save the parameters of the first molecule of the docking molecule file
c
            if(imix.eq.1) then
               xadds=xadd(1)
               yadds=yadd(1)
               zadds=zadd(1)
               phixs=phix(1)
               phiys=phiy(1)
               phizs=phiz(1)
               natos=nato(1)
               amassts=amasst(1)

c Ensure there is enough room to store the coordinates of this molecule.

               if(natos.gt.natomax) then
                  natomax=natos
                  call reallocate1d_c(alabels,3,natomax)
                  call reallocate1d_r(xas,natomax)
                  call reallocate1d_r(yas,natomax)
                  call reallocate1d_r(zas,natomax)
                  call reallocate1d_i(nbondss,natomax)
                  call reallocate2d_i(iabonds,natomax,near)
                  call reallocate2d_r(tbonds,natomax,near)
               end if
               do ia=naf(1),nal(1)
                  alabels(ia)=alabel(ia)
                  xas(ia)=xa(ia)
                  yas(ia)=ya(ia)
                  zas(ia)=za(ia)
                  nbondss(ia)=nbond(ia)
                  if(nbond(ia).gt.0) then
                     do ib=1,nbond(ia)
                        iabonds(ia,ib)=iabond(ia,ib)
                        tbonds(ia,ib)=tbond(ia,ib)
                     end do
                  endif
               end do
c
c rotational groups
c
               ngrots=ngrot(1)
               if(ngrots.gt.ngrotmax) then
                  ngrotmax=ngrots
                  call reallocate1d_c(rotgrps,3,ngrotmax)
                  call reallocate2d_i(labaxis,ngrotmax,2)
                  call reallocate1d_i(ndihedralangs,ngrotmax)
                  call reallocate2d_r(dihedralangs,ndihmax,ngrotmax)
                  call reallocate3d_i(labrotdihs,2,ndihmax,ngrotmax)
               end if
               if(ngrots.gt.0) then
                  do ig=1,ngrots
                     rotgrps(ig)=rotgrp(ig)
                     labaxis(ig,1)=labaxi(ig,1)
                     labaxis(ig,2)=labaxi(ig,2)
!Check if there are dihedral angles defined
                     if(rotgrp(ig).eq.'DIH') then
                        ndihedralangs(ig)=ndihedralang(ig)
                        if(ndihedralangs(ig).gt.ndihmax) then
                           ndihmax=ndihedralang(ig)
                           call reallocate2d_r(dihedralangs,ndihmax,ngrotmax)
                           call reallocate3d_i(labrotdihs,2,ndihmax,ngrotmax)
                        end if
                        do k=1,ndihedralang(ig)
                           dihedralangs(k,ig)=dihedralang(k,ig)
                           do j=1,2
                              labrotdihs(j,k,ig)=labrotdih(j,k,ig)
                           end do
                        end do
                     end if
                  end do
               endif
               write(6,*) 'dockato> natos,ngrots: ',natos,ngrots
            endif

            ifind=0
            ifindr=0
            ifindrr=0
            do while (ifind.eq.0.and.ifindr.eq.0.and.ifindrr.eq.0)
            if(imix.eq.1) then

c Specify the site on this molecule to be used for docking

               if(ilenatomd.eq.0) then
                  write(6,696)
696   format(/'dockato> Specify site to be used for docking: ',$)
                  read(5,*) atomd
               endif

c Specify the site on this molecule to be used for rotating

               if(ilenatomdr.eq.0) then
                  write(6,689)
689   format(/'dockato> Specify site to be used for rotating'
     1,' this molecule: ',$)
                  read(5,*) atomdr
               endif

c Specify the site on this molecule to be used for a second rotation

               if(ilenatomdrr.eq.0) then
                  write(6,687)
687   format(/'dockato> Specify site to be used for a second'
     1,' rotation of this molecule : ',$)
                  read(5,*) atomdrr
               endif

c Determine how many of these docking molecules are needed

               if(ndock.eq.0) then
                  write(6,697)
697   format(/'dockato> How many docking molecules are needed? ',$)
                  read(5,*) ndock
               endif

c Determine maximum distance from the first atom in file for docking purposes.

               if(radmax.eq.0.0) then
                  write(6,685)
685   format(/'dockato> Give maximum distance from central molecule to'
     *,' dock molecules: ',$)
                  read(5,*) radmax
               endif

               iad=0
               iadr=0
               iadrr=0

c Check the current atoms that at least one site of the requested type exists

               ia=0
               do while (ia.lt.nato(1))

                  ia=ia+1
                  if(alabel(ia).eq.atomd.and.ifind.eq.0) then

                     ifind=ia

c Save the atom number

                     iad=ia

c Save the name of the molecular template file

                  endif

                  if(alabel(ia).eq.atomdr.and.ifindr.eq.0) then

                     ifindr=ia

c Save the atom number

                     iadr=ia

                  endif

                  if(alabel(ia).eq.atomdrr.and.ifindrr.eq.0) then

                     ifindrr=ia

c Save the atom number

                     iadrr=ia

                  endif
               end do

c Check that one of the specified atoms on the docking molecule is a 'Qv' atom

c              if(atomd.ne.'Qv'.and.atomdr.eq.'Qv') then
c                 write(6,694)
c                 ifind=0
c                 ifindr=0
c                 ifindrr=0
c                 ilenatomd=0
c                 ilenatomdr=0
c                 ilenatomdrr=0
694   format(/'dockato> Atoms on the docking molecule are incorrectly specified'
     *,' - re-specify types') 
c              endif

            else

c Specify the site on this molecule to be used for docking on to.

               if(ilenatomv.eq.0) then
                  write(6,661)
661   format(/'dockato> Specify site to be docked onto: ',$)
                  read(5,*) atomv
               endif

c Specify the site on this molecule to be rotated on to

               if(ilenatomvr.eq.0) then
                  write(6,688)
688   format(/'dockato> Specify site to be rotated on to: ',$)
                  read(5,*) atomvr
               endif

c Specify the site on this molecule to be rotated on to

               if(ilenatomvrr.eq.0) then
                  write(6,686)
686   format(/'dockato> Specify site to be rotated on to, second rotation: ',$)
                  read(5,*) atomvrr
               endif

c Specify minimum dihedral angle for this atom must be for the docking molecule (0 for no second rotation)

               if(ilendminrr.eq.0) then
                  write(6,684,advance='no')
684   format(/'dockato> Specify minimum dihedral angle for second rotation ')
                  read(5,*) dminrr
               endif

c Specify maximum dihedral angle this atom must be for the docking molecule (0 for no second rotation)

               if(ilendmaxrr.eq.0) then
                  write(6,683,advance='no')
683   format(/'dockato> Specify maximum dihedral angle for second rotation: ')
                  read(5,*) dmaxrr
               endif
c Check whether to allow overlaps or not
               if(index(nooverlapstext,'<undefined>').gt.0) then
                    write(6,784,advance='no')
784     format(/'dockato> Type overlaps factor (0 or yes to allow all overlaps) ')
                    read(5,'(a)') nooverlapstext
               end if               

c Check that the current atoms have at least one site of the requested type
c and check that this atom is charged.

               ia=0
               do while (ia.lt.natom)

                  ia=ia+1
                  if(alabel(ia).eq.atomv.and.ifind.eq.0) then

                     ifind=ia

                  endif

                  if(alabel(ia).eq.atomvr.and.ifindr.eq.0) then

                     ifindr=ia

                  endif

                  if(alabel(ia).eq.atomvrr.and.ifindrr.eq.0) then

                     ifindrr=ia

                  endif

               end do


c Check that one of the specified atoms on the docking molecule is a 'Qv' atom

c              if(atomv.eq.'Qv'.and.atomvr.ne.'Qv') then
c                 write(6,693)
c                 ifind=0
c                 ifindr=0
c                 ifindrr=0
c                 ilenatomv=0
c                 ilenatomvr=0
c                 ilenatomvrr=0
693   format(/'dockato> Atoms on the docked molecule are incorrectly specified'
     *,' - re-specify types')             
c              endif
            endif
            end do

         else

            ntype=ntypefirst
            nmoltype=nmoltypefirst
            write(6,599) fnameato
599   format(/'Problems with .ato file ',a,' Try again.'/)
            if(imix.eq.1) ilendockingfile=0
            if(imix.eq.2) ilendockedfile=0
            imix=imix-1
            
         endif

      end do

c Maximum distance of docking molecule from origin

c Now proceed to dock molecules, 1 molecule at a time

      it=0
      imoladd=0

c Counter for the search for unsaturated sites.

      ivcount=0

c Check whether overlaps are allowed or not

! The overlap factor is set to the value contained in nooverlapstext, unless this value is a character string, in which case overlapfactor is set to 0.0
      if(index(nooverlapstext,'yes').gt.0) then
         nooverlapstext='0.0'
      else if(index(nooverlapstext,'nooverlaps').gt.0) then
         nooverlapstext='1.0'
      end if
      read(nooverlapstext,*,iostat=ierr) overlapfactor
      if(ierr.ne.0) overlapfactor = 0.0
      nooverlaps=overlapfactor.gt.0.0
      write(6,'("dockato> Overlap factor: ",f10.5)') overlapfactor


      do while(it.lt.ndock.and.ivcount.lt.natom)

         it=it+1

c Save the current numbers of atoms, molecules and rotational groups.

         natomsave=natom
         natomrsave=natomr
         nadd=natom
         write(6,*) 'dockato> nadd',nadd
         nmolsave=nmol
         ngrrotsave=ngrrot
         nrottotsave=nrottot
         dmaxsave=dmax

c Add the new molecule at the end of the current list of molecules.

         ia=nal(nmolsave)
         im=nmolsave
         ig=ngrrotsave
         if(im.lt.mmol) then
            im=im+1
         else
            write(6,601)
601   format(/'dockato> Tried to define more molecules than allowed by dimension'/)
            return
         endif

c Centre of mass coordinates

         xadd(im)=xadds
         yadd(im)=yadds
         zadd(im)=zadds
         phix(im)=phixs
         phiy(im)=phiys
         phiz(im)=phizs
         nato(im)=natos
         amasst(im)=amassts
         natom=natomsave+natos

c Atomic coordinates relative to centre of mass

         naf(im)=ia+1
         write(6,*) 'dockato> natos,ngrots: ',natos,ngrots
         do is=1,natos
            if(ia.lt.matom) then
               ia=ia+1
            else
               write(6,602)
602   format(/'dockato> Tried to define more atoms than allowed by dimension'/)
               return
            endif
            mref(ia)=im
            alabel(ia)=alabels(is)
            call ftype(ia)
c     write(6,*) is,ia,itype(ia)
c
c do not count atoms which begin with 'Q' or 'q' - these are assumed 
c invisible to the scattering radiation and experiment (although they
c can be used to calculate other structural quantities if necessary
c
            if(qtype(ia).ne.1) then
               natomr=natomr+1
            endif
            xa(ia)=xas(is)
            ya(ia)=yas(is)
            za(ia)=zas(is)
            nbond(ia)=nbondss(is)
            if(nbond(ia).gt.0) then
               do ib=1,nbond(ia)
                  iabond(ia,ib)=iabonds(is,ib)+nadd
                  tbond(ia,ib)=tbonds(is,ib)
               end do
            endif
         end do
         nal(im)=ia
c
c rotational groups
c
         ngrot(im)=ngrots
         if(ngrot(im).gt.0) then
            ngrrot=ngrrot+ngrot(im)
            call testsize_rotgrps()
            do ih=1,ngrot(im)
               ig=ig+1
c save molecule in which this rotational group occurs
               mrot(ig)=im
               write(6,*) 'dockato> ig, nadd',ig,nadd
               rotgrp(ig)=rotgrps(ih)
c for rotational groups read atom numbers of two atoms to form axis of rotation
               do j=1,2
                  labaxi(ig,j)=labaxis(ih,j)+nadd
               end do
               ndihedralang(ig)=ndihedralangs(ih)
               do k=1,ndihedralang(ig)
                  dihedralang(k,ig)=dihedralangs(k,ih)
                  do j=1,2
                     labrotdih(j,k,ig)=labrotdihs(j,k,ih)+nadd
                  end do
               end do
            end do
         endif
         natom=ia
         nmol=im
!         ngrrot=ig
         write(6,*) 'dockato > ',natom,natomr,nmol,ngrrot,(ngrot(i),i=1,nmol)


c New value of dmax

         dmax=dmaxsave*((real(natomr)/real(natomrsave))**(1.0/3.0))

c        do imt=1,nmol
c           write(6,*)
c           write(6,*) imt,xadd(imt),yadd(imt),zadd(imt)
c           write(6,*)
c           do iat=naf(imt),nal(imt)
c              write(6,*) iat,xa(iat),ya(iat),za(iat),nbond(iat),iabond(iat,1)
c           end do
c        end do
c        write(6,*)

c Of the available atoms search for a suitable site to be docked onto and a site suitable for docking.

         ifind=0
         iav=0
         iavr=natomsave
         iavrr=natomsave
         iad=natom
         iadr=natom
         iadrr=natom

c Identify the next unsaturated site of the specified type within the atoms to be docked onto.

         do while (iav.lt.natomsave.and.ifind.eq.0)

            if(iad.ge.natom) then

               iad=natomsave
               call search_unsaturated(iav,natomsave,nvalencev,atomv,atomvr,xv,yv,zv)
      write(6,*) 'dockato> 1 ',iav,nvalencev,xv,yv,zv

            endif


            if(nvalencev.gt.0) then

c O.K. so found an unsaturated site.
c Now identify an unsaturated site on the docking molecule of specified type

c Molecule number for this site

               imv=mref(iav)
               call search_unsaturated(iad,natom,nvalenced,atomd,atomdr,xd,yd,zd)

      write(6,*) 'dockato> 2 ',iad,nvalenced,xd,yd,zd

               if(nvalenced.gt.0) then

c OK so this atom pair apparently meets the requirements for a valid docking. 
c But we do not yet know if we dock a molecule here whether:-
c (a) it can be rotated onto a valid docking position;
c (b) after translation and rotation it will not overlap existing atoms 
c (as defined by the Lennard-Jones parameters)

c The current plan is to simply dock onto the valid site found, then rotate as required, then
c search for overlaps. The docking is only accepted if no overlaps occur after rotation. 
c Since the molecules are listed in the .ato file in the order in which they are entered, NOT in the
c order in which they appear in the docking sequence, this should ensure all the available sites for 
c docking eventually become occupied. 

c So the next thing to do is to perform a translation of this molecule as required, 
c and then check for overlaps.

                  imd=mref(iad)
      write(6,*) iav,iad,xv,yv,zv,xd,yd,zd
                  xadd(imd)=xadds+xv-xd
                  yadd(imd)=yadds+yv-yd
                  zadd(imd)=zadds+zv-zd
c     write(6,*) xadd(imd),yadd(imd),zadd(imd)

c Now perform rotations as required, if possible

c Search for a site on the docking molecule to be used for defining the rotation.

                  iadr=0
                  tbondrot=1000.0
                  ib=0
                  do while (ib.lt.nbond(iad))

c We choose the atom of the correct type which is the shortest distance from iad.
c This is to ensure that if there is more than one atom of type atomdr in the neighbour list 
c of iad on the docking molecule,we pick the nearest one.

                     ib=ib+1
                     tbondtest=abs(tbond(iad,ib))
                     if(alabel(iabond(iad,ib)).eq.atomdr
     *.and.tbondtest.lt.tbondrot) then
                        iadr=iabond(iad,ib)
                        tbondrot=tbondtest
                     endif
                  end do
      write(6,*) iad,iadr,tbondrot

c Repeat this search on the molecule to be docked onto. In this case we are looking for
c an atom of type atomvr which is a neighbour of iav at a distance of tbondrot.

                  iavr=0
                  ib=0
                  do while (ib.lt.nbond(iav).and.iavr.eq.0)

                     ib=ib+1
                     ivref=iabond(iav,ib)
                     if(alabel(ivref).eq.atomvr.and.abs(tbond(iav,ib)).eq.tbondrot) then

                        iavr=ivref

                     endif
                  end do
      write(6,*) iav,iavr,tbondrot
                  if(iadr.gt.0.and.iavr.gt.0) then

c If iadr and iavr are defined we perform the first rotation

                     call rotatem(imd,iad,iadr,imv,iav,iavr,ierr)

c Now attempt to perform 2nd rotation if required and possible

c Now perform a second rotation if possible about the axis defined by iad to iadr

                     call random_rotatem(atomdrr,iad,iadr,iadrr
     *,atomvrr,iav,iavr,iavrr,dminrr,dmaxrr,dihedralangfixed)

      write(6,*) iav,iavr,iavrr,iad,iadr,iadrr,tbondrot

c Test for overlaps if required

                    ioverlap=0
                    if(nooverlaps) call test_overlap(ioverlap,imd,0,natomsave
     *,iav,iavr,iavrr,iad,iadr,iadrr,overlapfactor)

                     if(ioverlap.eq.0) then

                        ifind=1

                     else

                        ifind=0

c If this rotation is not allowed, restore the coordinates of
c the docking molecule

                        xadd(imd)=xadds
                        yadd(imd)=yadds
                        zadd(imd)=zadds

                        ia=naf(imd)-1
                        do is=1,nato(imd)
                           ia=ia+1
                           xa(ia)=xas(is)
                           ya(ia)=yas(is)
                           za(ia)=zas(is)
                        end do

                     endif
                  endif
               endif
            endif
         end do

c If ifind.eq.0 at this point, we have searched all possible pairs of sites, 
c and not one is allowed.
c Therefore we have to abandon this docking attempt.

c If no site was found, ignore this added molecule.

         write(6,664) ifind
664   format(/'dockato> ifind = ',i5)
         if(ifind.eq.0) then

            natom=natomsave
            natomr=natomrsave
            nmol=nmolsave
            ngrrot=ngrrotsave
            dmax=dmaxsave

         else

c Molecule docking successful. Therefore bond the new molecule to the old. This means eliminating
c atoms iad and iadr, and giving there respective neighbours to the atoms they are bonded to.

            fixeddihedral=dminrr.eq.dmaxrr
            call bond_molecule(iav,iavr,iavrr,iad,iadr,iadrr,fixeddihedral,dihedralangfixed)
            write(6,*) 'dockato> ',nmol,natom,natomr,nato(imv),nal(imv),nato(imd)

            imoladd=imoladd+1
            dmaxmol=sizefactor*dmax

c If search was successful reset docked site counter to ensure it rescans the current
c site, in case another molecule can be docked here.

            iav=iav-1

         endif

      end do

! We now need to go through all rotational groups and setup the lists of atoms to be rotated within each group.

! First need to get the directly bonded atoms before we can do this,
! since get_rotgrp works only on the directly bonded atoms

      call get_direct_bonds()

      nrottot=0
      nrottotmax=0
      do ig=1,ngrrot
         ir1=labaxi(ig,1)
         ir2=labaxi(ig,2)
         call get_rotgrp(ir1,ir2,labrott1,nrott1)
         call get_rotgrp(ir2,ir1,labrott2,nrott2)
         !Use the smallest of the two groups
         if(nrott1.lt.nrott2) then
            labaxi(ig,1)=ir1
            labaxi(ig,2)=ir2
            nrot(ig)=nrott1
            rotf(ig)=nrottot+1
            nrottot=nrottot+nrot(ig)
            rotl(ig)=nrottot
            call testsize_labrot()
            do i=1,nrot(ig)
               labref=rotf(ig)+i-1
      write(6,*) 'dockato> get labrot> ',ig,rotf(ig),rotl(ig),nrottot,nrottotmax,labref
               labrot(labref)=labrott1(i)
            end do
         else
            labaxi(ig,1)=ir2
            labaxi(ig,2)=ir1
            nrot(ig)=nrott2
            rotf(ig)=nrottot+1
            nrottot=nrottot+nrot(ig)
            rotl(ig)=nrottot
            call testsize_labrot()
            do i=1,nrot(ig)
               labref=rotf(ig)+i-1
      write(6,*) 'dockato> get labrot> ',ig,rotf(ig),rotl(ig),nrottot,nrottotmax,labref
               labrot(labref)=labrott2(i)
            end do
         end if
      end do

      write(6,662) imoladd
662   format(/'dockato> Added ',i5,' molecules')
      if(imoladd.lt.ndock) write(6,663)
663   format(/'dockato> (no more valid sites found for docking.)')

c
c get the new filename to put the mixture in
c
      if(ilenoutputfile.gt.0) then
         nwords=1
         ncf(1)=1
         ncl(1)=ilenoutputfile
         sometext(1:ilenoutputfile)=outputfile
      else
         nwords=0
      endif
      do while (nwords.le.0)
         write(6,35)
35    format(/'dockato> Type name of file to put mixture in: ',$)
         read(5,101,iostat=ierr) sometext
101   format(a)
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         if(nwords.le.0) then
            write(6,499)
         endif
      end do
      readfile=sometext(ncf(1):ncl(1))
      ilenreadfile=ncl(1)-ncf(1)+1
c
c strip out any extension that was supplied - extension must be .ato, otherwise
c nothing is done.
c
      ilendot=index(readfile,'.')-1
      if(ilendot.gt.0) then
         readfile=readfile(1:ilendot)
         ilenreadfile=ilendot
      endif
      fnameato=directory_name(1:nchardir)//readfile(1:ilenreadfile)//'.ato'
      ilenfnameato=nchardir+ilenreadfile+4

! Reset the molecule template file for this molecule to the current filename, so that
! any changes to this file will be reflected in any .ato files which are derived from
! this file

      moltemplate(mtype(imv))=readfile(1:ilenreadfile)

c
c write the new .ato file
c
      call moloutc(fnameato)
      return
      end

      subroutine rotatem(imd,iad,iadr,imv,iav,iavr,ierr)

c Rotates molecule imd about atom iad so that atom iadr on imd overlaps atom iavr
c on molecule imv. It is assumed the two molecules have been previously shifted so that
c atom iad on imd overlaps atom iav on imv. If not the routine returns without doing
c anything. Minimum image convention is not used for this check so this method is only
c to be used in conjunction with dockato.

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      real*4 dsq(2)

c vvector defines the vector on the existing molecule on to which the docking molecule
c is to be rotated
c rvector defines the vector on the docking molecule which is to be rotated to coincide
c with vvector.
c pvector is the vector about which the rotation occurs, namely the vector normal to the 
c plane defined by vvector and rvector

      real*4 vvector(3),rvector(3),pvector(3)
c xvector is the vector from centre of rotation to an individual atom
c xpvector is the component of this vector along the axis of rotation
c xrvector is the component of this vector perpendicular to the axis of rotation
c xrpvector is the vector perpendicular to both axis of rotation and xrvector

      real*4 xvector(3),xpvector(3),xrvector(3),xrpvector(3)

c First check that atom iad does indeed overlap atom iav

      write(6,*) imv,iav,iavr,imd,iad,iadr
      write(6,*) itype(iav),itype(iavr),itype(iad),itype(iadr)

      ierr=1
      xadiff=xadd(imv)-xadd(imd)+xa(iav)-xa(iad)
      yadiff=yadd(imv)-yadd(imd)+ya(iav)-ya(iad)
      zadiff=zadd(imv)-zadd(imd)+za(iav)-za(iad)
      dsq(1)=xadiff*xadiff+yadiff*yadiff+zadiff*zadiff
      qradiussum=6.0*(qradius(itype(iad))+qradius(itype(iav)))
      write(6,*) xadiff,yadiff,zadiff,sqrt(dsq(1)),qradiussum
      if(sqrt(dsq(1)).gt.qradiussum) then
         write(6,701)
701   format('rotatem> Molecules do not overlap - cannot rotate')
         return
      else if(iavr.eq.0.or.iadr.eq.0) then
         write(6,702)
702   format('rotatem> No rotation is requested')
         return
      endif

c generate unit vectors for vvector and rvector, measured from the atom
c being docked onto

      vvector(1)=xadd(imv)+xa(iavr)-xadd(imd)-xa(iad)
      vvector(2)=yadd(imv)+ya(iavr)-yadd(imd)-ya(iad)
      vvector(3)=zadd(imv)+za(iavr)-zadd(imd)-za(iad)
      rvector(1)=xa(iadr)-xa(iad)
      rvector(2)=ya(iadr)-ya(iad)
      rvector(3)=za(iadr)-za(iad)
      vnorm=sqrt(scapro(vvector,vvector))
      rnorm=sqrt(scapro(rvector,rvector))
c     write(6,*) (vvector(i),rvector(i),i=1,3)
      do i=1,3
         vvector(i)=vvector(i)/vnorm
         rvector(i)=rvector(i)/rnorm
      end do

c angle between these vectors

      cphi=scapro(vvector,rvector)
      if(abs(cphi).gt.1.0) cphi=sign(1.0,cphi)
      sphi=sqrt(1.0-cphi*cphi)

c No need to rotate if the atoms already overlap

      write(6,*) cphi,sphi
      if(sphi.gt.0.001) then

c define the vector normal to the plane of vvector and rvector about which the molecule
c will be rotated

         call vecpro(rvector,vvector,pvector)

c Define the unit vector in this direction

         rab=sqrt(scapro(pvector,pvector))
      write(6,*) rab,(pvector(i),i=1,3)

         pvector(1)=pvector(1)/rab
         pvector(2)=pvector(2)/rab
         pvector(3)=pvector(3)/rab
      write(6,*) (pvector(i),i=1,3)
c step through atoms in group
         do ia=naf(imd),nal(imd)
            if(ia.ne.iad) then
c Vector from the centre of rotation to this atom
               xvector(1)=xa(ia)-xa(iad)
               xvector(2)=ya(ia)-ya(iad)
               xvector(3)=za(ia)-za(iad)
c Component of this vector parallel to axis of rotation
               frac=scapro(pvector,xvector)
               do i=1,3
                  xpvector(i)=frac*pvector(i)
               enddo
      write(6,*) ia,(xvector(i),xpvector(i),i=1,3)
      
c Component of vector perpendicular to axis of rotation
               do i=1,3
                  xrvector(i)=xvector(i)-xpvector(i)
               enddo
c Vector perpendicular to plane defined by axis of rotation and xrvector
               call vecpro(pvector,xrvector,xrpvector)
c Generate rotated coordinates. The vector component parallel to the axis of rotation
c is not rotated
               xa(ia)=xa(iad)+xpvector(1)
               ya(ia)=ya(iad)+xpvector(2)
               za(ia)=za(iad)+xpvector(3)
               xa(ia)=xa(ia)+xrvector(1)*cphi+xrpvector(1)*sphi
               ya(ia)=ya(ia)+xrvector(2)*cphi+xrpvector(2)*sphi
               za(ia)=za(ia)+xrvector(3)*cphi+xrpvector(3)*sphi
            endif
         end do

      else if(cphi.lt.1.0) then

c Only can get here if sphi=0 and cphi=-1. In that case we simply invert the coordinates
c about the centre of rotation. (Until I think of something better to do!)

         do ia=naf(imd),nal(imd)

            xa(ia)=xa(iad)-(xa(ia)-xa(iad))
            ya(ia)=ya(iad)-(ya(ia)-ya(iad))
            za(ia)=za(iad)-(za(ia)-za(iad))

         end do

      endif
      ierr=0
      return
      end

      subroutine random_rotatem(atomdrr,iad,iadr,iadrr
     *,atomvrr,iav,iavr,iavrr,dminrr,dmaxrr,dihedralangfixed)

c Generate a random dihedral rotation of the docking molecule about the axis defined 
c by iad,iadr.
c dminrr and dmaxrr give the minimum and maximum dihedral angles requested.
c The dihedral angle is measured using the existing axis of rotation and, plus the atom
c of type atomvrr on the docked molecule and the atom of type atomdrr on the docking 
c molecule. 

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      character*3 atomdrr,atomvrr

c vectors to be used for rotating molecules

      real*4 rvector(3),vvector(3),rvpvector(3),dvector(3)
     *,rdpvector(3),drot(3),dminrr,dmaxrr,dihedralangfixed

      write(6,100) iad,iadr,atomdrr,iav,iavr,atomvrr
100   format(/'random_rotatem> Beginning...',2(1x,i5,1x,i5,1x,a3)/)

c Coordinates of axis of rotation

      imd=mref(iad)

c Vector defining axis of rotation

      rvector(1)=xa(iadr)-xa(iad)
      rvector(2)=ya(iadr)-ya(iad)
      rvector(3)=za(iadr)-za(iad)

c Convert to a unit vector in the same direction

      rnorm=sqrt(scapro(rvector,rvector))
      do i=1,3
        rvector(i)=rvector(i)/rnorm
      end do
      
c Angle of rotation, selected at random between dminrr and dmaxrr

      pi=4.0*atan(1.0)
      arg=2.0*pi*(dminrr+ran1()*(dmaxrr-dminrr))/360.0
      dihedralangfixed=arg*180.0/pi
      sphistep=sin(arg)
      cphistep=cos(arg)

c Search for an appropriate atom on the docking molecule. We choose the one that is closest
c to iad

      imd=mref(iad)
      iat=naf(imd)-1
      iadrr=0
      dsqmin=1.0e30

      do while(iat.lt.nal(imd).and.atomdrr.ne.'0')

         iat=iat+1

c The atom must be of the correct type and cannot be either the rotation centre (iad)
c or the atom already rotated (iadr).

         if(iat.ne.iad.and.iat.ne.iadr.and.alabel(iat).eq.atomdrr) then

c Vector defining the position of this atom relative to the centre of rotation

            dvector(1)=xa(iat)-xa(iad)
            dvector(2)=ya(iat)-ya(iad)
            dvector(3)=za(iat)-za(iad)

c Square of distance to iad

            dsq=scapro(dvector,dvector)
            if(dsq.lt.dsqmin) then

c Cross product of rvector onto dvector gives a vector perpendicular to both

               call vecpro(rvector,dvector,rdpvector)

c Normalise the perpendicular vector

               rdpnorm=sqrt(scapro(rdpvector,rdpvector))
               do i=1,3
                 rdpvector(i)=rdpvector(i)/rdpnorm
               end do

c Dot product with rvector should not be >1 or <-1

               cphi=scapro(rdpvector,rvector)
               if(cphi.gt.-1.0.and.cphi.lt.1.0) then

                  dsqmin=dsq
                  iadrr=iat

               endif
            endif
         endif
      end do      

c Search for an appropriate atom on the docked molecule. We choose the atom
c which is closest to iavr

      imv=mref(iav)
      iat=naf(imv)-1
      iavrr=0
      dsqmin=1.0e30

      do while(iat.lt.nal(imv).and.atomvrr.ne.'0')

         iat=iat+1

c The atom must be of the correct type and cannot be either the rotation centre (iav)
c or the atom already rotated onto (iavr).

         if(iat.ne.iav.and.iat.ne.iavr.and.alabel(iat).eq.atomvrr) then

c Vector defining the position of this atom relative to the axis of rotation

            vvector(1)=xa(iat)-xa(iav)
            vvector(2)=ya(iat)-ya(iav)
            vvector(3)=za(iat)-za(iav)

c Square of distance to iad

            dsq=scapro(vvector,vvector)
            if(dsq.lt.dsqmin) then

c Cross product of rvector onto vvector gives a vector perpendicular to both

               call vecpro(rvector,vvector,rvpvector)

c Normalise the perpendicular vector

               rvpnorm=sqrt(scapro(rvpvector,rvpvector))
               do i=1,3
                  rvpvector(i)=rvpvector(i)/rvpnorm
               end do

c Dot product with rvector should not be 1 or -1, otherwise vvector is parallel with
c rvector and so cannot be used to rotate.

               cphi=scapro(rvpvector,rvector)
      write(6,100)
      write(6,*) iat,cphi
               if(cphi.gt.-1.0.and.cphi.lt.1.0) then

                  dsqmin=dsq
                  iavrr=iat

               endif
            endif
         endif
      end do

c Now determine angle of rotation

      if(iadrr.gt.0.and.iavrr.gt.0) then

c In this case we attempt to rotate molecule imd so that the angle between rdpvector
c and rvpvector = arg.  

c Get the current angle between rdpvector and rvpvector

c Cosine of the angle is given by the dot product.

            cphiold=scapro(rdpvector,rvpvector)
            if(abs(cphiold).ge.1.0) then
               cphiold=sign(1.0,cphiold)
            endif

c Need to determine the sine. Form a vector perpendicular to both
c rvpvector and rdpvector - this will be parallel or antiparallel to rvector

            call vecpro(rdpvector,rvpvector,drot)

c If we have to rotate backwards, this will be antiparallel to rvector, otherwise it
c will be parallel. Since rvector is a unit vector, the dot product wih rvector will
c give the sine of the angle between rdpvector and rvpvector

            sphiold=scapro(drot,rvector)

c Define the new cosine and sine required

            cphinew=cphiold*cphistep+sphiold*sphistep
            sphinew=sphiold*cphistep-cphiold*sphistep

      else if(iadrr.gt.0) then

c In this case, if iavrr is not found, simply perform a rotation of the docking
c molecule by the amount phinew

            sphinew=sphistep
            cphinew=cphistep

      else

c Otherwise do nothing

            cphinew=2.0

      endif

c Now perform rotation as requested

      if(cphinew.ge.-1.0.and.cphinew.le.1.0) then

         iat=naf(imd)-1
         do while (iat.lt.nal(imd))

            iat=iat+1
            if(iat.ne.iad.and.iat.ne.iadr) then

c Vector defining the position of this atom relative to the centre of rotation

               dvector(1)=xa(iat)-xa(iad)
               dvector(2)=ya(iat)-ya(iad)
               dvector(3)=za(iat)-za(iad)

c Component of dvector along the axis of rotation

               frac=scapro(dvector,rvector)

c Cross product of rvector onto dvector gives a vector perpendicular to both

               call vecpro(rvector,dvector,rdpvector)

               do i=1,3

c Vector parallel to axis of rotation

                  dpar=frac*rvector(i)
                  
c Vector perpendicular to axis of rotation

                  dperp=dvector(i)-dpar

c Create the rotated vector

                  drot(i)=dpar+dperp*cphinew+rdpvector(i)*sphinew
 
               end do

c Generate rotated coordinates

               xa(iat)=xa(iad)+drot(1)
               ya(iat)=ya(iad)+drot(2)
               za(iat)=za(iad)+drot(3)

            endif
         end do
      endif

c As a test, if iavrr.gt.0, check that the new angle between iavrr and iadrr is 
c as requested.


      cphitest=99.0
      sphitest=99.0
      if(iavrr.gt.0) then

c Vector defining the position of this atom relative to the centre of rotation

         dvector(1)=xa(iadrr)-xa(iad)
         dvector(2)=ya(iadrr)-ya(iad)
         dvector(3)=za(iadrr)-za(iad)

c Cross product of rvector onto dvector gives a vector perpendicular to both

         call vecpro(rvector,dvector,rdpvector)

c Normalise the perpendicular vector

         rdpnorm=sqrt(scapro(rdpvector,rdpvector))
         do i=1,3
            rdpvector(i)=rdpvector(i)/rdpnorm
         end do

c Cosine of the angle is given by the dot product.

         cphitest=scapro(rdpvector,rvpvector)
         if(abs(cphitest).ge.1.0) then
            cphitest=sign(1.0,cphitest)
         endif

c Need to determine the the sine. Form a vector perpendicular to both
c rvpvector and rdpvector - this will be parallel or antiparallel to rvector

         call vecpro(rdpvector,rvpvector,drot)

c If we have to rotate backwards, this will be antiparallel to rvector, otherwise it
c will be parallel

         sphitest=scapro(drot,rvector)

      endif

      write(6,101) iad,iadr,iadrr,iav,iavr,iavrr
101   format(/'random_rotatem> Atoms selected: ',6(1x,i5))
      write(6,102) cphinew,sphinew
102   format(/'random_rotatem> Angle rotated (cosine and sine) ',1x,e13.6,1x,e13.6)
      if(iavrr.gt.0) then
         write(6,103) cphiold,sphiold
103   format(/'random_rotatem> Angle old (cosine and sine) ',1x,e13.6,1x,e13.6)
         write(6,104) cphitest,sphitest
104   format(/'random_rotatem> Angle check (cosine and sine) ',1x,e13.6,1x,e13.6/)
      endif

      return
      end

      subroutine search_unsaturated(iav,iavlast,nvalence,atomv,atomvr,xv,yv,zv)

c Finds the first unsaturated site of the specified type amongst the existing atoms

c On entry iav is the number of the atom to search from. The first search should have iav=0.
c On exit iav is the atom number of the atom found suitable for docking.
c ivcount is > 0 for a successful search, =0 for an unsuccessful search
c atomv is the label of the site to be docked onto
c charged is the charge on the docking site
c xad,yad,zad are the coordinates of the specific atom on the docking molecule to be used
c for docking.
c radmax is the maximum distance of docking molecule from origin
c xvbar, etc. are the coordinates of this docking site
c xaddt, etc. are the proposed coordinates of the docking molecule so that the atom specified by
c xad,etc. overlaps with iav. Note this is BEFORE any attempted rotation of the molecule.

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      character*3 atomv,atomvr

      integer*4 iav,iavlast,imv

c We search the neighbours of iav to see whether there any at the existing value of iav.
c This is so that if a previous search found this
c atom, but discarded because of overlaps, we do not search it again.

      nvalence=0
      nvalencefound=0
      bondsearch=1000.0
      ia=iav

      do while (ia.lt.iavlast)

c Search for an atom of the specified type

         ia=ia+1
         imv=mref(ia)

      write(6,*) 'search_unsaturated> ',ia,iavlast,' ',alabel(ia),' ',atomv,' ',atomvr

         if(alabel(ia).eq.atomv) then

c Check the neighbours of ia and sum the number of ligands as specified by atomvr atoms
c in the neighbourhood of iav

            call find_valence(ia,nvalencefound,atomvr,bondfound)

c If a non-zero valence value was identified, select the shortest bond distance found

            if(nvalencefound.gt.0.and.bondfound.lt.bondsearch) then
               iav=ia
               bondsearch=bondfound
               nvalence=nvalencefound
            endif
         endif
      end do

c If zero valence if found, we save the last value of ia, so that it is not searched again

      if(nvalence.eq.0) then

         iav=ia

      endif

      write(6,*) 'search_unsaturated> ',iav,imv,nvalence,' ',alabel(iav),' ',atomvr

c Coordinates of this atom

      imv=mref(iav)
      xv=xadd(imv)+xa(iav)
      yv=yadd(imv)+ya(iav)
      zv=zadd(imv)+za(iav)


      return
      end

      subroutine find_valence(iav,nvalence,atomvr,bondsearch)

c Search through the neighbours of iav and calculate the number of overlaps,
c mean position of the charge, and total charge. Atoms iavfirst to iavlast are used
c for the search.

c searchmax is the maximum separation of a neighbouring atom from iav for it to be counted 
c as overlapping

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      character*3 atomvr

      integer*4 iav,nvalence

c Look at all the bonded neighbours of iav, if any, and find the shortest distance 
c to atoms of type atomvr. This distance is assumed to be the distance of local valence
c atoms.

      ib=0
      ifind=0
      bondsearch=1000.0

      do while(ib.lt.nbond(iav))

         ib=ib+1
         if(alabel(iabond(iav,ib)).eq.atomvr) then

            ifind=ifind+1
            if(abs(tbond(iav,ib)).lt.bondsearch) bondsearch=abs(tbond(iav,ib))

         endif
      end do
      write(6,100)
      write(6,*) bondsearch,atomvr

c Now search for neighbours of type atomvr which are at this minimum distance.

      ib=0
      nvalence=0

      do while(ib.lt.nbond(iav))

         ib=ib+1
c     write(6,100)
c     write(6,*) ib,iabond(iav,ib),alabel(iabond(iav,ib)),tbond(iav,ib),nvalence

         if(alabel(iabond(iav,ib)).eq.atomvr.and.abs(tbond(iav,ib)).eq.bondsearch) then

            nvalence=nvalence+1

         endif

      end do
100   format(/'find_valence> ',$)
      write(6,100)
      write(6,*) nvalence

      return
      end

      subroutine test_overlap(ioverlap,imd,natomfirst,natomsave
     *,iav,iavr,iavrr,iad,iadr,iadrr,overlapfactor)

c Check for overlaps, according to the Lennard-Jones parameters.

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      real overlapfactor
                     
      ioverlap=0

c Step through the atoms of the new molecule

      inew=naf(imd)-1

      do while (inew.lt.nal(imd).and.ioverlap.eq.0)

         inew=inew+1
         ityped=itype(inew)

c Coordinates of this atom

         xt=xadd(imd)+xa(inew)
         yt=yadd(imd)+ya(inew)
         zt=zadd(imd)+za(inew)

c Test whether this atom is a neighbour of iad, iadr or iadrr

         call find_neighbour(iad,inew,iadneighb)
         call find_neighbour(iadr,inew,iadrneighb)
         call find_neighbour(iadrr,inew,iadrrneighb)

c Step through the other atoms and make sure none of them come within the prescribed
c distance

         iat=natomfirst
         do while(iat.lt.natomsave.and.ioverlap.eq.0)
      
            iat=iat+1

c Test whether this atom is a neighbour of iad, iadr or iadrr

            call find_neighbour(iav,iat,iavneighb)
            call find_neighbour(iavr,iat,iavrneighb)
            call find_neighbour(iavrr,iat,iavrrneighb)

            imt=mref(iat)
            itt=itype(iat)
            epsilont=sqrt(epscor(itt)*epscor(ityped))
c     write(6,*) inew,iat,ityped,itt
c     write(6,*) iadneighb,iadrneighb,iavneighb,iavrneighb,epsilont

c Only do the test if epsilont is non-zero, and iat is not a neighbour of iavr, and
c iat is not a neighbour of iav.

            if(epsilont.gt.0.0.and.((iavneighb.eq.0.or.iadrneighb.eq.0)
     *.or.(iavrneighb.eq.0.or.iadneighb.eq.0))) then
                  
               sigmat=0.5*overlapfactor*(sigcor(itt)+sigcor(ityped))/dmax
               xx=(xadd(imt)+xa(iat)-xt)/dmax
               yy=(yadd(imt)+ya(iat)-yt)/dmax
               zz=(zadd(imt)+za(iat)-zt)/dmax
               call minimagec(xx,yy,zz)
               dtest=sqrt(xx*xx+yy*yy+zz*zz)
               if(dtest.lt.sigmat) then
                  ioverlap=1
      write(6,100,advance='no')
      write(6,*) iat,imt,itt,inew,imd,ityped,dtest,sigmat
               endif
            endif
         end do
      end do

100   format(/'test_overlap> ')

      return
      end

      subroutine find_neighbour(iad,inew,iadneighb)

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      iadneighb=0
 
c We include the case when the atom specified (inew) is equal to the base atom specified (iad) as
c a neighbour, since this atom should not be included in overlap calculation.

      if(iad.eq.inew) iadneighb=1
      if(iad.gt.0) then

         ib=0
         do while (ib.lt.nbond(iad).and.iadneighb.eq.0)

            ib=ib+1
            if(iabond(iad,ib).eq.inew) iadneighb=1

         end do

      endif 

      return
      end   

      subroutine bond_molecule(iav,iavr,iavrr,iad,iadr,iadrr,fixeddihedral,dihedralangfixed)

      use reallocation_routines
      use molecom_mod
      use dockato_arrays

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

c List of atoms to be kept and erradicated

      integer noverlap,ikeep(near),ilose(near)
      real dihedralangfixed

c Temporary store to ensure existing rotational groups get inserted at the correct place

      logical found,foundiav,foundiad,fixeddihedral

c Temporary rotational groups used to reassess atoms to be rotated

c        integer nrott1,nrott2,labrott1(near),labrott2(near)


c Adds the list of atoms in the molecule given by iad,iadr, iadrr, to the molecule given by 
c iav, iavr, iavrr.

c First set all the coordinates of the atoms in the new molecule to be consistent with the
c centre of mass coordinates of the old molecule.

      imd=mref(iad)
      imv=mref(iav)
      write(6,*) 'bond_molecule> 1 ',nmol,(nato(j),j=1,nmol),natom
      write(6,*) 'bond_molecule> 2 ',iav,iavr,iavrr,iad,iadr,iadrr
      write(6,*) 'bond_molecule> 3 ',imv,imd,ngrot(imv),natos,ngrots,ngrrot

c Transfer any existing rotational groups on imd to imv

      if(ngrots.gt.0) then

c Find out the last group corresponding to imv

         iglast=0
         im=0
         do while(im.lt.imv)
      
            im=im+1
            iglast=iglast+ngrot(im)

         end do

c Shift all the existing groups above this up by ngrots.

         ig=ngrrot
         ig1=ig-ngrots
c        write(6,*) 'Add ROT groups> ',ngrrot,ngrots,ig,ig1,iglast
         do while (ig1.gt.iglast)

            mref(ig)=mref(ig1)
            rotgrp(ig)=rotgrp(ig1)
            do j=1,2
               labaxi(ig,j)=labaxi(ig1,j)
            end do
            ndihedralang(ig)=ndihedralang(ig1)
            do k=1,ndihedralang(ig)
               dihedralang(k,ig)=dihedralang(k,ig1)
               do j=1,2
                  labrotdih(j,k,ig)=labrotdih(j,k,ig1)
               end do
            end do
            ig=ig-1
            ig1=ig1-1

         end do

c Finally insert the saved rotational groups into the space created

         ngrot(imv)=ngrot(imv)+ngrots
         ig=iglast
         ih=0

         do while (ih.lt.ngrots)

            ig=ig+1
            ih=ih+1
            mref(ig)=imv
            rotgrp(ig)=rotgrps(ih)
            do j=1,2
               labaxi(ig,j)=labaxis(ih,j)+nadd
            end do
            ndihedralang(ig)=ndihedralangs(ih)
            do k=1,ndihedralang(ig)
               dihedralang(k,ig)=dihedralangs(k,ih)
               do j=1,2
                  labrotdih(j,k,ig)=labrotdihs(j,k,ih)+nadd
               end do
            end do
         end do
         write(6,*) 'bond_molecule> 4 ',imv,imd,ngrot(imv),ngrots,ngrrot
      endif


c If iavrr or iadrr are zero it means no dihedral angle was defined. Therefore we need to generate a
c random rotation group of all the atoms in imv about the bond iadr to iav. However this can only be
c done if the number of atoms on the docked molecule is greater than 2.

c      if(nato(imv).gt.2.and.(iavrr.eq.0.or.iadrr.eq.0)) then

c From 5 Sept 2012 we will always generate a random rotational group about a bond unless the dihedral
c angle has been fixed

! From 20/01/2014, a random rotational group is generated if the dihedral angle has not been fixed,
! Otherwise a dihedral group will be generated.

      if(nato(imv).gt.2) then
         ngrot(imv)=ngrot(imv)+1
         ngrrot=ngrrot+1
         call testsize_rotgrps()

c First find out the last group corresponding to imv

         iglast=0
         im=0
         do while(im.lt.imv)
      
            im=im+1
            iglast=iglast+ngrot(im)

         end do

c Now shift all the existing groups above this up by 1

         ig=ngrrot
         do while (ig.gt.iglast)

            ig1=ig-1
            rotgrp(ig)=rotgrp(ig1)
            do j=1,2
               labaxi(ig,j)=labaxi(ig1,j)
            end do
            ndihedralang(ig)=ndihedralang(ig1)
            do k=1,ndihedralang(ig)
               dihedralang(k,ig)=dihedralang(k,ig1)
               do j=1,2
                  labrotdih(j,k,ig)=labrotdih(j,k,ig1)
               end do
            end do
            ig=ig1

         end do

c Now insert the new rotational group

         rotgrp(ig)='ROT'
         labaxi(ig,1)=iadr
         labaxi(ig,2)=iav
         ndihedralang(ig)=0

c Ensure a finite step size is given

         if(stepri.le.0.0) stepri=1.0

      endif

c If both iadrr and iavrr are .gt. 0 then we add in a bond between these two atoms. 
c However after rotation it is possible these and other atoms now overlap, so we need to eliminate
c duplicate atoms in this case. If iadrr and iavrr are not included in the list of overlaps, then
c 20/01/2014 we set up a dihedral group between them. In all cases we always lose iad and iavr
c and keep iav and iadr

      noverlap=2
      ikeep(1)=iav
      ilose(1)=iad
      ikeep(2)=iadr
      ilose(2)=iavr

      ioverlaprr=0
      if(iadrr.gt.0.and.iavrr.gt.0.) then

c First determine if any atoms on imd overlaps any atoms on molecule imv - if so the atom on imd
c will be amalgamated with the corresponding atom on imv

         iadt=naf(imd)-1
         do while (iadt.lt.nal(imd))

            iadt=iadt+1

c Do not include atoms iad and iadr in the test

            if(iadt.ne.iad.and.iadt.ne.iadr) then

               x1=xadd(imd)+xa(iadt)
               y1=yadd(imd)+ya(iadt)
               z1=zadd(imd)+za(iadt)

               iavt=naf(imv)-1
               do while (iavt.lt.nal(imv))

                  iavt=iavt+1
c     write(6,*) imv,imd,iavt,iadt
                  if(iavt.ne.iav.and.iavt.ne.iavr) then


                     x2=xadd(imv)+xa(iavt)
                     y2=yadd(imv)+ya(iavt)
                     z2=zadd(imv)+za(iavt)
                     xdiff=(x2-x1)/dmax
                     ydiff=(y2-y1)/dmax
                     zdiff=(z2-z1)/dmax
                     call minimagec(xdiff,ydiff,zdiff)
                     dist=sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff)*dmax
                     dmin=qradius(itype(iadt))+qradius(itype(iavt))
                     if(dist.le.dmin) then

                        noverlap=noverlap+1
                        ikeep(noverlap)=iavt
                        ilose(noverlap)=iadt

c If this pair includes iadrr and iavrr, then record the fact

                        if(iadt.eq.iadrr.and.iavt.eq.iavrr) ioverlaprr=1

                     endif
                  endif
               end do
            endif
         end do
      endif
      
c Shift all the coordinates on imd so that they are with reference to imv as the centre.

      do ia=naf(imd),nal(imd)

         xa(ia)=xa(ia)+xadd(imd)-xadd(imv)
         ya(ia)=ya(ia)+yadd(imd)-yadd(imv)
         za(ia)=za(ia)+zadd(imd)-zadd(imv)
         mref(ia)=imv

      end do

c Change references to atoms on imd as required by overlaps

      do ioverlap=1,noverlap

c      write(6,*) ioverlap,ikeep(ioverlap),ilose(ioverlap)

         call replace_references(ilose(ioverlap),ikeep(ioverlap))
         if(iavrr.eq.ilose(ioverlap)) iavrr=ikeep(ioverlap)
         if(iadrr.eq.ilose(ioverlap)) iadrr=ikeep(ioverlap)

      end do

c Transfer neighbours from imd to imv

      do ioverlap=1,noverlap

         call add_neighbours(ilose(ioverlap),ikeep(ioverlap))

      end do


c If atoms iavrr and iadrr do not overlap determine the distance between these two atoms,
c and set up a new bond between them, provided this pair does not already exist.

c From 5 Sept 2012 we only create this distance if a particular
c dihedral angle is specified.

c From 20/01/2014 we create a dihedral group instead of a distance.

      if(ioverlaprr.eq.0.and.fixeddihedral) then

c First check that this pair does not already appear in each others neighbour lists

         ib=0
         ifound=0
         do while (ib.lt.nbond(iavrr).and.ifound.eq.0)

            ib=ib+1
            if(iabond(iavrr,ib).eq.iadrr) ifound=iadrr

         end do
         ib=0
         do while (ib.lt.nbond(iadrr).and.ifound.eq.0)

            ib=ib+1
            if(iabond(iadrr,ib).eq.iavrr) ifound=iavrr

         end do
         if(ifound.eq.0) then

            if(fixeddihedral) then
               rotgrp(ig)='DIH'
               ndihedralang(ig)=1
               dihedralang(1,ig)=dihedralangfixed
               labrotdih(1,1,ig)=iadrr
               labrotdih(2,1,ig)=iavrr
            end if

!            x1=xadd(mref(iadrr))+xa(iadrr)
!            y1=yadd(mref(iadrr))+ya(iadrr)
!            z1=zadd(mref(iadrr))+za(iadrr)
!            x2=xadd(mref(iavrr))+xa(iavrr)
!            y2=yadd(mref(iavrr))+ya(iavrr)
!            z2=zadd(mref(iavrr))+za(iavrr)
!            xdiff=(x2-x1)/dmax
!            ydiff=(y2-y1)/dmax
!            zdiff=(z2-z1)/dmax
!            call minimagec(xdiff,ydiff,zdiff)
!            dist=sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff)*dmax

c Atoms do not overlap, so generate a dihedral group.

!            nbond(iavrr)=nbond(iavrr)+1
!            iabond(iavrr,nbond(iavrr))=iadrr
!            tbond(iavrr,nbond(iavrr))=dist
!            nbond(iadrr)=nbond(iadrr)+1
!            iabond(iadrr,nbond(iadrr))=iavrr
!            tbond(iadrr,nbond(iadrr))=dist

         endif
      endif

c Set up the new number of atoms and mass for molecule imv

      nato(imv)=nato(imv)+nato(imd)
      nal(imv)=nal(imv)+nato(imd)
      amasst(imv)=amasst(imv)+amasst(imd)

c This process should have removed all references to the overlapping atoms so we can remove them
c from the list of atoms.

      do ioverlap=1,noverlap
         call remove_atoms(ilose(ioverlap),imv)
c Make sure we update the list references to atoms lost and kept
         io=0
         do while (io.lt.noverlap)
            io=io+1
            if(ikeep(io).gt.ilose(ioverlap)) ikeep(io)=ikeep(io)-1
            if(ilose(io).gt.ilose(ioverlap)) ilose(io)=ilose(io)-1
         end do
      end do

c     do ia=1,natom
c     write(6,*) ia
c     write(6,*) (iabond(ia,j),j=1,nbond(ia))
c     end do

c Need to step through the rotation groups and make sure the rotation axes are not duplicated

      igb=0
      do while (igb.lt.ngrrot)
         igb=igb+1
         ig=igb
         do while (ig.lt.ngrrot)
            ig=ig+1
c     write(6,*) labaxi(ig,1),labaxi(ig,2),labaxi(igb,1),labaxi(igb,2)
            if((labaxi(ig,1).eq.labaxi(igb,1).and.labaxi(ig,2).eq.labaxi(igb,2))
     *.or.(labaxi(ig,2).eq.labaxi(igb,1).and.labaxi(ig,1).eq.labaxi(igb,2))) then
! Found two axes the same - eliminate the second.
               ngrrot=ngrrot-1
               igt=ig-1
               do while (igt.lt.ngrrot)

                  igt=igt+1
                  igt1=igt+1
                  rotgrp(igt)=rotgrp(igt1)
                  do j=1,2
                     labaxi(igt,j)=labaxi(igt1,j)
                  end do
                  ndihedralang(igt)=ndihedralang(igt1)
                  if(ndihedralang(igt).gt.0) then
                     do k=1,ndihedralang(igt)
                        dihedralang(k,igt)=dihedralang(k,igt1)
                        do j=1,2
                           labrotdih(j,k,igt)=labrotdih(j,k,igt1)
                        end do
                     end do
                  end if
               end do
! Decrement ig by 1, so that it looks at the new rotational group in this place
               ig=ig-1
            endif
         end do
      end do

c Will need to recalculate the number of rotational groups for each molecule

      do im=1,nmol
         ngrot(im)=0
      end do
      ig=0
      do while (ig.lt.ngrrot)
         ig=ig+1
         im1=mref(labaxi(ig,1))
         im2=mref(labaxi(ig,2))
         write(6,*) ig,labaxi(ig,1),im1,labaxi(ig,2),im2
         if(im1.eq.im2) then
             ngrot(im1)=ngrot(im1)+1
             mrot(ig)=im1
         endif
      end do

c Check the existing types to see if any of previous types do not now exist in the .ato file.

c     write(6,*) ntype

      itb=0

      do while (itb.lt.ntype)

         itb=itb+1
         ia=0
         itfound=0
         do while(ia.lt.natom.and.itfound.eq.0) 
            ia=ia+1
c     write(6,*) ia,itype(ia)
            if(itype(ia).eq.itb) itfound=itb
         end do

c Remove this type if it is not found

         if(itfound.eq.0) then
            ntype=ntype-1
            it=itb-1
            do while (it.lt.ntype)
               it=it+1
               it1=it+1
               type(it)=type(it1)
               atsymbol(it)=atsymbol(it1)
               iexchange(it)=iexchange(it1)
               epscor(it)=epscor(it1)
               sigcor(it)=sigcor(it1)
               amassa(it)=amassa(it1)
               charge(it)=charge(it1)
               qradius(it)=qradius(it1)
            end do

c Now shift the type numbers of the atoms if necessary

            ia=0
            do while(ia.lt.natom)
               ia=ia+1
               if(itype(ia).gt.itb) itype(ia)=itype(ia)-1
            end do
         endif
      end do

c Shift all internal atomic coordinates to be relative to centre of mass

      call centreofmass(imv)

c Finally drop the number of molecules by 1.

      nmol=nmol-1

      write(6,*) 'bond_molecule> 5 ',nmol,natom,natomr,ntype,ngrrot

      return
      end

      subroutine replace_references(iad,iav)

      use reallocation_routines
      use molecom_mod

c Replaces references to iad with references to iav.

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      imd=mref(iad)
      imv=mref(iav)

c Need to change the type number as appropriate

      itype(iad)=itype(iav)

c Substitute iav for iad in the bond reference numbers

      do ia=1,natom
         ib=0
         do while (ib.lt.nbond(ia))
            ib=ib+1
            if(iabond(ia,ib).eq.iad) iabond(ia,ib)=iav
         end do
      end do

c Now the rotational groups

      ig=0
      do while (ig.lt.ngrrot)
         ig=ig+1
         if(mrot(ig).eq.imd) mrot(ig)=imv
         do j=1,2
            if(labaxi(ig,j).eq.iad) labaxi(ig,j)=iav
         end do
         if(ndihedralang(ig).gt.0) then
            do k=1,ndihedralang(ig)
               do j=1,2
                  if(labrotdih(j,k,ig).eq.iad) labrotdih(j,k,ig)=iav
               end do
            end do
         end if
      end do

      return
      end

      subroutine add_neighbours(iad,iav)

      use reallocation_routines
      use molecom_mod

c Adds the neighbours of iad to iav, but only if this neighbour does not already appear in the list
c of neighbours for this atom.

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

c The neighbours of iad will become the new neighbours of iav

      ib=0
      do while (ib.lt.nbond(iad))
         ib=ib+1
         iaref=iabond(iad,ib)

c First check that this atom is not already listed as a neighbour to iav: only
c add it to the list if it is not already present.

         ibc=0
         ifound=0
         do while (ibc.lt.nbond(iav).and.ifound.eq.0)
            ibc=ibc+1
            if(iabond(iav,ibc).eq.iaref) ifound=1
         end do
         if(ifound.eq.0) then
            nbond(iav)=nbond(iav)+1
            iabond(iav,nbond(iav))=iaref
            tbond(iav,nbond(iav))=tbond(iad,ib)
         endif
         ibc=0
         ifound=0
         do while (ibc.lt.nbond(iaref).and.ifound.eq.0)
            ibc=ibc+1
            if(iabond(iaref,ibc).eq.iav) ifound=1
         end do
         if(ifound.eq.0) then
            nbond(iaref)=nbond(iaref)+1
            iabond(iaref,nbond(iaref))=iav
            tbond(iaref,nbond(iaref))=tbond(iad,ib)
         endif
      end do

      return
      end

      subroutine remove_atoms(iad,imv)

      use reallocation_routines
      use molecom_mod

c Removes atom iad which is no longer in use. imv is the molecule number on which iad occurs

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      integer :: typesave
      logical :: ifoundtype 

c New numbers of atoms in molecule imv

      nato(imv)=nato(imv)-1
      nal(imv)=nal(imv)-1
      if(qtype(iad).ne.1) natomr=natomr-1

c Shift all the atoms number references above iad down by 1

      ia=iad
      do while(ia.lt.natom)
         ia1=ia+1
         xa(ia)=xa(ia1)
         ya(ia)=ya(ia1)
         za(ia)=za(ia1)
         alabel(ia)=alabel(ia1)
         itype(ia)=itype(ia1)
         qtype(ia)=qtype(ia1)
         nbond(ia)=nbond(ia1)
         mref(ia)=mref(ia1)
         ib=0
         do while (ib.lt.nbond(ia))
            ib=ib+1
            iabond(ia,ib)=iabond(ia1,ib)
            tbond(ia,ib)=tbond(ia1,ib)
         end do
         ia=ia1
      end do

c Now shift all bond reference numbers down by 1 if greater than iad

      do ia=1,natom
         ib=0
         do while(ib.lt.nbond(ia))
            ib=ib+1
            if(iabond(ia,ib).gt.iad) then
               iabond(ia,ib)=iabond(ia,ib)-1
            endif

         end do
      end do

c Repeat the last step for the rotational groups. 

      ig=0
      do while (ig.lt.ngrrot)
         ig=ig+1
         do j=1,2
            if(labaxi(ig,j).gt.iad) labaxi(ig,j)=labaxi(ig,j)-1
         end do
         if(ndihedralang(ig).gt.0) then
            do k=1,ndihedralang(ig)
               do j=1,2
                  if(labrotdih(j,k,ig).gt.iad) labrotdih(j,k,ig)=labrotdih(j,k,ig)-1
               end do
            end do
         end if
      end do

c New numbers of atom in the box

      natom=natom-1

      call clean_atom_types() !Clean up the list of atom types to correspond with what is actually in the box

      return
      end

      subroutine clean_atom_types()

      use reallocation_routines
      use molecom_mod

c Removes atom type iatomtype which is no longer in use. It is assumed this has already been checked.

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      integer iatomtype,it,it1
      logical ifoundatomtype


      iatomtype=0
      do while(iatomtype.lt.ntype)
         iatomtype=iatomtype+1
         ia=1
         ifoundatomtype=type(iatomtype).eq.alabel(ia)
         do while (ia.lt.natom.and..not.ifoundatomtype)
            ia=ia+1
            ifoundatomtype=type(iatomtype).eq.alabel(ia)
         end do
         if(.not.ifoundatomtype) then
            write(6,*) 'clean_atom_type> ',iatomtype,ntype,(' ',type(it),it=1,ntype)
            if(iatomtype.lt.ntype) then
               do it=iatomtype,ntype-1
                  it1=it+1
                  type(it)=type(it1)
                  qtypea(it)=qtypea(it1)
                  iatomicnotype(it)=iatomicnotype(it1)
                  atsymbol(it)=atsymbol(it1)
                  iexchange(it)=iexchange(it1)
                  epscor(it)=epscor(it1)
                  sigcor(it)=sigcor(it1)
                  amassa(it)=amassa(it1)
                  charge(it)=charge(it1)
                  qradius(it)=qradius(it1)
                  natomsphere(it)=natomsphere(it1)            
                  rhosphere(it)=rhosphere(it1)
                  radminsphere(it)=radminsphere(it1)
                  radmaxsphere(it)=radmaxsphere(it1)
               end do
               do ia=1,natom
                  if(itype(ia).gt.iatomtype) itype(ia)=itype(ia)-1
               end do
               ntype=ntype-1
            else
               ntype=ntype-1
            end if
         end if
      end do
      return
      end



