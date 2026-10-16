      subroutine get_molecules()

      use reallocation_routines
c
c resets the value of a specified parameter in the .ato file if requested 
c by calling program
c
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      include 'chemical_symbols.inc'

      integer, dimension(:), allocatable              :: atomtypelist,natomslist
      integer, dimension(:,:), allocatable            :: atomslist
      integer :: lenatomtypelist,it,jt,ia,maxatoms
      real vanderWaalsFactor,xi,yi,zi,xj,yj,zj,xdif,ydif,zdif,dtest
      real, dimension(melements) :: vanderWaalsRadius

      logical                                         :: found

! Arrays for new molecule definitions
      integer, dimension(:), allocatable              :: natonew,nafnew,nalnew,nbondnew,mrefnew,imageatomnew,itypenew
      integer, dimension(:,:), allocatable            :: newmolatomlist,iabondnew  
      real, dimension(:,:), allocatable               :: tbondnew  
      real, dimension(:), allocatable                 :: xanew,yanew,zanew  
      character(len=3), dimension(:), allocatable     :: alabelnew
      integer                                         :: nmolnew
      logical, dimension(:), allocatable              :: tetherednew
      integer, dimension(:), allocatable              :: tetheredatomnew
! Arrays to store old molecule types
      real, dimension(:), allocatable                 :: itold,irold,wrold,wtold !Step sizes for current molecules
      character(len=3), dimension(:), allocatable     :: moltypeold
      character(len=256), dimension(:), allocatable   :: moltemplateold
      integer                                         :: nmoltypeold

      call normalise_coordinates()

!Get the list of van der Waals radii
      open(10,file='vanderWaalsRadii.txt',status='old',iostat=ierr)
      if(ierr.ne.0) then
         write(6,*) 'get_molecules> Unable to open vanderWaalsRadii.txt'
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
               icount=icount+1
            end if
         end if
         read(10,'(a)',iostat=ierr) sometext
      end do
      close(10)
      write(6,48,advance='no') icount
48    format(/1x,'get_molecules> Number of van der Waals radii read in: ',i3
     *//' Give van der Waals factor (suggested value 0.8, but check): ')
      read(5,*) vanderWaalsFactor
      ierr=0
! Allocate space to save molecule template information
      call reallocate1d_r(itold,nmoltype)
      call reallocate1d_r(irold,nmoltype)
      call reallocate1d_r(wrold,nmoltype)
      call reallocate1d_r(wtold,nmoltype)
      call reallocate1d_c(moltypeold,3,nmoltype)
      call reallocate1d_c(moltemplateold,256,nmoltype)
! And save it
      nmoltypeold=nmoltype
      do it=1,nmoltype
         itold(it)=interntrans(it)
         irold(it)=internrot(it)
         wrold(it)=wholerot(it)
         wtold(it)=wholetrans(it)
         moltypeold(it)=moltype(it)
         moltemplateold(it)=moltemplate(it)
      end do
!For this exercise we need to set the coordinates of each atom to their absolute values
! (because we are going to alter the molecule reference numbers)
      do im=1,nmol
         do ia=naf(im),nal(im)
            xa(ia)=xa(ia)+xadd(im)
            ya(ia)=ya(ia)+yadd(im)
            za(ia)=za(ia)+zadd(im)
         end do
         xadd(im)=0.0
         yadd(im)=0.0
         zadd(im)=0.0
      end do
      write(6,24) nmol
24    format(/1x,'get_molecules> There are ',i7,' molecules in this file'/)
      do it=1,nmoltypeold
            write(6,27) it,moltypeold(it)
27    format(1x,'Molecule type ',i2,' has label ',a3)
      end do
      write(6,32) natom,ntype
32    format(/1x,'get_molecules> There are ',i7,' atoms and ',i2,' types of atom in this file'/)
      do it=1,ntype
            write(6,33) it,type(it)
33    format(1x,'Atom type ',i2,' has label ',a3)
      end do
      write(6,381,advance='no')
381   format(/1x,'get_molecules> Give pairs of atom types needed to define bonds (e.g. 3 4 4 2), 0 to exit): ')
      read(5,'(a)') sometext
      write(6,*)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(allocated(atomtypelist)) deallocate(atomtypelist)
      if(nwords.lt.2) then !Need a minimum of two numbers to generate bonds
         return
      else
! Only use the numbers in pairs so ignore any extra words at end
         nwords=2*(nwords/2)
         lenatomtypelist=nwords
         allocate(atomtypelist(nwords))
         read(sometext(ncf(1):ncl(nwords)),*) (atomtypelist(it),it=1,nwords)
      end if
      call reallocate1d_i(natomslist,lenatomtypelist)
      do it=1,lenatomtypelist
! Make a list of atoms of each specified type. First count the number, then save the atom numbers
         natomslist(it)=0
      end do
      do ia=1,natom
         do it=1,lenatomtypelist
            if(itype(ia).eq.atomtypelist(it)) natomslist(it)=natomslist(it)+1
         end do
      end do
! Get the largest of these lists
      maxatoms=0
      do it=1,lenatomtypelist
         maxatoms=max(maxatoms,natomslist(it))
         natomslist(it)=0
      end do
      call reallocate2d_i(atomslist,maxatoms,lenatomtypelist)
! Get the values
      do ia=1,natom
         do it=1,lenatomtypelist
            if(itype(ia).eq.atomtypelist(it)) then
               natomslist(it)=natomslist(it)+1
               atomslist(natomslist(it),it)=ia
            end if
         end do
      end do
      do it=1,lenatomtypelist
         write(6,57) atomtypelist(it),natomslist(it)
57    format(' get_molecules> For atom type ',i2,1x,i5,' atoms found')
         write(6,*) (atomslist(j,it),j=1,natomslist(it))
      end do
! Create an array to store the sequence of atom numbers for each new molecule
! Make sure this is large enough to cope with the existing molecules
      maxatoms=0
      do im=1,nmol
         maxatoms=max(maxatoms,nato(im))
      end do
      call reallocate2d_i(newmolatomlist,maxatoms,nmol)
! Populate this list
      write(6,168) 
168   format(/1x,'get_molecules> Initial molecule no,, number of atoms, atom list'/)
      do im=1,nmol
         do ia=1,nato(im)
            iaref=ia+naf(im)-1
            newmolatomlist(ia,im)=iaref
         end do
!         write(6,*) 'get_molecules> ',im,nato(im),(newmolatomlist(ia,im),ia=1,nato(im))
      end do
! Look for bonded atom pairs between the specified atom types, using atom type pairs in list
      it=1
      do while(it.lt.lenatomtypelist)
         jt=it+1
         do il=1,natomslist(it)
            ia=atomslist(il,it)
            xi=xa(ia) !Coordinates of this atom
            yi=ya(ia)
            zi=za(ia)
            do jl=1,natomslist(jt)
               ja=atomslist(jl,jt)
               xj=xa(ja) !Coordinates of this atom
               yj=ya(ja)
               zj=za(ja)
               if(ia.ne.ja.and.nato(mref(ja)).gt.0) then !Don't include atoms from molecules which no longer have any atoms
                  xdif=(xj-xi)
                  ydif=(yj-yi)
                  zdif=(zj-zi)
                  call minimagec(xdif,ydif,zdif)
                  dtest=dmax*sqrt(xdif*xdif+ydif*ydif+zdif*zdif)
                  im=mref(ia)
                  jm=mref(ja)
                  dvanderwaals=vanderWaalsFactor*(vanderWaalsRadius(iatomicnotype(itype(ia)))+vanderWaalsRadius(iatomicnotype(itype(ja))))
                  found=dtest.lt.dvanderwaals!Need to check these atoms EVEN if they are already part of the same molecule
!                   write(6,*) ia,ja,dtest,dvanderwaals
!                   write(6,*) ia,ja,im,jm,itype(ia),itype(ja),iatomicnotype(itype(ia)),iatomicnotype(itype(ja))
!                   write(6,*) ia,ja,im,jm,nato(im),nato(jm),dtest,dvanderwaals
                  if(found) then
                     write(6,*) ia,ja,im,jm,nato(im),nato(jm),dtest,dvanderwaals
                     if(nbond(ia).lt.near.and.nbond(ja).lt.near) then
                        !need to check that ja is not already a neighbour of ia
                        found=.false.
                        ib=0
                        do while (ib.lt.nbond(ia).and..not.found)
                           ib=ib+1
                           found=iabond(ia,ib).eq.ja
                        end do
                        if(.not.found) then
                           nbond(ia)=nbond(ia)+1
                           iabond(ia,nbond(ia))=ja
                           tbond(ia,nbond(ia))=dtest/dmax
                           nbond(ja)=nbond(ja)+1
                           iabond(ja,nbond(ja))=ia
                           tbond(ja,nbond(ja))=dtest/dmax
                           !Only do the next bit if im.ne.jm
                           if(im.ne.jm) then
                              maxatoms=max(maxatoms,(nato(im)+nato(jm)))
                              call reallocate2d_i(newmolatomlist,maxatoms,nmol)
                              !Need to give all the atoms contained in molecule mref(ja) the new molecule number
                              ireflist=nato(im)
                              do jb=1,nato(jm) !No. of atoms in each molecule is kept up to date
                                 jat=newmolatomlist(jb,jm) !Use newmolatomlist to get the atoms in this molecule because naf and nal are not reliable here
                                 ireflist=ireflist+1
                                 mref(jat)=im
                                 newmolatomlist(ireflist,im)=jat
                              end do
                              nato(im)=nato(im)+nato(jm)
                              nato(jm)=0
                           end if
                        endif
                     end if
!                      write(6,*) ia,ja,im,jm,nato(im),nato(jm)
                  end if
               end if
            end do
         end do
         it=jt+1 ! Get next pair
      end do
! Determine the maximum number of atoms in any molecule, and the number of molecules with non-zero
! numbers of atoms
      maxatoms=0
      nmolnew=0
      write(6,237) 
237   format(/1x,'get_molecules> Final molecule no,, number of atoms, atom list'/)
      do im=1,nmol
         maxatoms=max(maxatoms,nato(im))
!         write(6,*) ' get_molecules> ',im,nato(im),(newmolatomlist(ia,im),ia=1,nato(im))
         if(nato(im).gt.0) then
            nmolnew=nmolnew+1
            write(6,*) 'get_molecules> ',im,nato(im),(newmolatomlist(ia,im),ia=1,nato(im))
         end if
      end do
! Create arrays needed for temporary storage during resorting of atoms
      call reallocate1d_i(natonew,nmolnew)
      call reallocate1d_i(nafnew,nmolnew)
      call reallocate1d_i(nalnew,nmolnew)
      call reallocate1d_i(nbondnew,natom)
      call reallocate1d_i(itypenew,natom)
      call reallocate2d_i(iabondnew,natom,near)
      call reallocate2d_r(tbondnew,natom,near)
      call reallocate1d_i(mrefnew,natom)
      call reallocate1d_r(xanew,natom)
      call reallocate1d_r(yanew,natom)
      call reallocate1d_r(zanew,natom)
      call reallocate1d_c(alabelnew,3,natom)
      call reallocate1d_i(imageatomnew,natom)
      if(allocated(tetherednew)) deallocate(tetherednew,tetheredatomnew)
      allocate(tetherednew(nmolnew),tetheredatomnew(nmolnew))
! Create the new molecules
      imnew=0
      ianew=0
      do im=1,nmol
         if(nato(im).gt.0) then
            imnew=imnew+1
            natonew(imnew)=nato(im)
            tetherednew(imnew)=tethered(im)
            tetheredatomnew(imnew)=tetheredatom(im)
            nafnew(imnew)=ianew+1
            do iaref=1,nato(im)
               ia=newmolatomlist(iaref,im)
               ianew=ianew+1
               alabelnew(ianew)=alabel(ia)
               xanew(ianew)=xa(ia)
               yanew(ianew)=ya(ia)
               zanew(ianew)=za(ia)
               imageatomnew(ianew)=imageatom(ia)
               nbondnew(ianew)=nbond(ia)
               itypenew(ianew)=itype(ia)
               do ib=1,nbondnew(ianew)
! Search for the correct bond reference atom
                  found=.false.
                  is=0
                  do while(is.lt.nato(im).and..not.found)
                     is=is+1
                     found=iabond(ia,ib).eq.newmolatomlist(is,im)
                  end do
                  iabondnew(ianew,ib)=is+nafnew(imnew)-1
                  tbondnew(ianew,ib)=tbond(ia,ib)
               end do
               mrefnew(ianew)=imnew
            end do
            nalnew(imnew)=ianew
         end if
      end do
      nmolnew=imnew
! Finally overwrite the existing .ato file arrays
      nmol=nmolnew
      do im=1,nmol
         nato(im)=natonew(im)
         naf(im)=nafnew(im)
         nal(im)=nalnew(im)
         tethered(im)=tetherednew(im)
         tetheredatom(im)=tetheredatomnew(im)
         xadd(im)=0.0
         yadd(im)=0.0
         zadd(im)=0.0
         phix(im)=0.0
         phiy(im)=0.0
         phiz(im)=0.0
         do ia=naf(im),nal(im)
            alabel(ia)=alabelnew(ia)
            xa(ia)=xanew(ia)
            ya(ia)=yanew(ia)
            za(ia)=zanew(ia)
            imageatom(ia)=imageatomnew(ia) !Need to set this value once new box is complete
            nbond(ia)=nbondnew(ia)
            itype(ia)=itypenew(ia)
            do ib=1,nbond(ia)
               iabond(ia,ib)=iabondnew(ia,ib)
               tbond(ia,ib)=tbondnew(ia,ib)
            end do
            mref(ia)=mrefnew(ia)
         end do
! For each atom in a molecule the coordinates should be relative to the centre of mass - not the periodic image
         do ia=naf(im),nal(im)
            x1=xa(ia)-xa(naf(im))
            y1=ya(ia)-ya(naf(im))
            z1=za(ia)-za(naf(im))
            x2=x1
            y2=y1
            z2=z1
            call minimagec(x2,y2,z2)
            !write(6,*) ia,dmax*sqrt(x1*x1+y1*y1+z1*z1),dmax*sqrt(x2*x2+y2*y2+z2*z2)
            if(abs(x2).lt.abs(x1)) xa(ia)=x2+xa(naf(im))
            if(abs(y2).lt.abs(y1)) ya(ia)=y2+ya(naf(im))
            if(abs(z2).lt.abs(z1)) za(ia)=z2+za(naf(im))
         end do
         call centreofmass(im)
         call minimagec(xadd(im),yadd(im),zadd(im)) !Place molecule at its minimum image position
         if(tethered(im)) then ! Set the tethering points if required.
            phix(im)=xadd(im)
            phiy(im)=yadd(im)
            phiz(im)=zadd(im)
         end if
      end do
      call unnormalise_coordinates()
!
! get_molecules will normally be called either after makelattice - generate a crystalline lattice from a set of crystal
! coordinates - or after readxyz - generate a .ato file from a .xyz file. In the former case the atoms may have sensible
! names, but in the latter case they certainly won't since the xyz file consists of the chemical symbol of each atom
! and its (absolute) xyz coordinates. In that case the chemical symbol becomes the atom label. 
! Hence, in either case, for molecular mixtures, there may be at this point more than one molecule type that begins with
! the same atom label, causing confusion to EPSR, which requires the first atom of each molecule type to have
! a distinct label. Equally there may be several atoms with the same label, e.g. C, but which have different bonding
! patterns.
! Make a new list of moltypes. We need to check the old list to be sure that existing moltypes
! are not overwritten, unless they do not exist.
c
c make a list of the molecule types read in from the current file
c
      nmoltype=0
      im=0
      do while(im.lt.nmol)
         im=im+1
! First check whether this molecule type already occurs in the new list
         imoltype=0
         ifoundnew=0
         do while(imoltype.lt.nmoltype.and.ifoundnew.eq.0)
            imoltype=imoltype+1
            if(alabel(naf(im)).eq.moltype(imoltype)) then
               ifoundnew=imoltype
            end if
         end do
         if(ifoundnew.eq.0) then !Does not already occur, so now check the old list
            imoltype=0
            ifoundold=0
            do while(imoltype.lt.nmoltypeold.and.ifoundold.eq.0)
               imoltype=imoltype+1
               if(alabel(naf(im)).eq.moltypeold(imoltype)) then
                  ifoundold=imoltype
               end if
            end do
            if(ifoundold.eq.0) then !Not found in old list either so create a new mol type
               nmoltype=nmoltype+1
               moltype(nmoltype)=alabel(naf(im))
c initialise the template filename
               write(moltemplate(nmoltype),'("moltype",i2.2)') nmoltype
c Initialise the rotational and translational step sizes
               interntrans(nmoltype)=0.0
               internrot(nmoltype)=0.0
               wholerot(nmoltype)=stepr
               wholetrans(nmoltype)=stepm
            else !Found in old list
! getting values from old list
               nmoltype=nmoltype+1
               moltype(nmoltype)=alabel(naf(im))
c Copy the template filename
               moltemplate(nmoltype)=moltemplateold(ifoundold)
c Copy the rotational and translational step sizes
               interntrans(nmoltype)=itold(ifoundold)
               internrot(nmoltype)=irold(ifoundold)
               wholerot(nmoltype)=wrold(ifoundold)
               wholetrans(nmoltype)=wtold(ifoundold)
            endif
         endif
      end do
      write(6,389) nmol,nmoltype
389   format(/1x,'get_molecules> There are now ',i7,' molecules and ',i2,' types of molecule in this file'/)
      do it=1,nmoltype
         write(6,392) it,moltype(it),moltemplate(it)(1:10),interntrans(it),internrot(it),wholerot(it),wholetrans(it)
      end do
392   format(1x,'Molecule type ',i2,' has label ',a3,1x,a,4(1x,f10.5))
      write(6,413)
413   format(/1x,'get_molecules> *** Use the label command to give atoms distinct names if needed ***')
      return
      end
 
