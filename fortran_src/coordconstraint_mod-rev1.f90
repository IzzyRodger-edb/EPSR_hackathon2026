!Used to store the data and methods associated with setting up a coordination constraint in EPSR
!Started on 6th June 2017

module coordconstraint_mod

implicit none

type coordconstraint

   character(len=3)                                :: atc !Atom type at the centre of the constraint
   character(len=3), dimension(:), allocatable     :: ats !List of possible secondary atom types
   integer                                         :: nc,alternatenc !Expected overall coordination number associated with this constraint, and the alternate coordination number for the coordinating species
   real                                            :: cctemp,ccgaussweight,ccgausswidth !temperature factor associated with this coordination constraint
   real                                            :: ccdistort !Distortion on coordination potential. Small (0.01) = symmetric.
   integer                                         :: ns !Number of secondary atoms associated with this constraint
   integer                                         :: atclistsize
   integer, dimension(:), allocatable              :: atclist,nbondcc,nearnbondcc,expecnbondcc !List of atom numbers of type atc
   integer, dimension(:,:), allocatable            :: iabondcc,wtcc !List of secondary (coordinating) atoms for each atc.
   real, dimension(:), allocatable                 :: ds,dsmax !List of expected bond distances associated with atom types ats
   real, dimension(:,:), allocatable               :: tbondcc,adistcc !List of expected bonds and actual distances for each neighbour.
!   real, dimension(:), allocatable                 :: distcc !Temporary distances used in three-body energy calculation
   logical                                         :: hastb,hastbintra !Determines if calculation of tb terms is to be performed
   real                                            :: tbweight,tbgaussweight,tbintraweight,tbintragaussweight
   real                                            :: tbdistort,tbgausswidth,tbintradistort,tbintragausswidth
   integer                                         :: ntbangles,ntbintraangles
   real, dimension(:), allocatable                 :: tbangles,tbcosines,tbintraangles,tbintracosines
   character(len=3), dimension(:), allocatable     :: tbintraats !List of possible secondary atom types for intra angles

!   integer, dimension(:), allocatable              :: nbondcc !No of neighbours for each atom in the simulation
!   integer, dimension(:,:), allocatable            ::  !List of neighbours for each atom

end type coordconstraint

contains

    subroutine init_constraint(this,sometext)

      use reallocation_routines

      type(coordconstraint) :: this
      character(len=*), intent(in) :: sometext
      integer, parameter :: nworddim=100,nwordspc=10
      integer :: nwords,nchartext,i,iread
      integer, dimension(nworddim) :: ncf,ncl
      real :: piconv

      nchartext = len(sometext)
      call parse(sometext,nwords,ncf,ncl,nwordspc,nchartext,nworddim)
      write(6,*) 'init_constraint> ',nwords,' ',sometext(ncf(1):ncl(nwords))
      this%atc=sometext(ncf(1):ncl(1))
      read(sometext(ncf(2):ncl(2)),*) this%nc
      read(sometext(ncf(3):ncl(3)),*) this%cctemp
      read(sometext(ncf(4):ncl(4)),*) this%ccdistort
      read(sometext(ncf(5):ncl(5)),*) this%ccgaussweight
      read(sometext(ncf(6):ncl(6)),*) this%ccgausswidth
      iread=6
      this%ns=0
      this%ntbangles=0
      this%ntbintraangles=0
      this%tbweight=0.0
      this%tbdistort=1.0
      this%tbgaussweight=0.0
      this%tbgausswidth=1.0
      this%tbintraweight=0.0
      this%tbintradistort=1.0
      this%tbintragaussweight=0.0
      this%tbintragausswidth=1.0
      this%hastb=.false.
      this%hastbintra=.false.
      do while (iread.lt.nwords)
         iread=iread+1
         if(.not.this%hastb.and..not.this%hastbintra) then
            this%hastb=sometext(ncf(iread):ncl(iread)).eq.'angles'
            this%hastbintra=sometext(ncf(iread):ncl(iread)).eq.'anglesintra'
            if(.not.this%hastb.and..not.this%hastbintra) then
               this%ns=this%ns+1
               call reallocate1d_c(this%ats,3,this%ns)
               call reallocate1d_r(this%ds,this%ns)
               call reallocate1d_r(this%dsmax,this%ns)
               this%ats(this%ns)=sometext(ncf(iread):ncl(iread))
               iread=iread+1
               read(sometext(ncf(iread):ncl(iread)),*) this%ds(this%ns)
               iread=iread+1
               read(sometext(ncf(iread):ncl(iread)),*) this%dsmax(this%ns)
            end if
         else if(this%hastb.and..not.this%hastbintra) then
            read(sometext(ncf(iread):ncl(iread)),*) this%tbweight
            iread=iread+1
            read(sometext(ncf(iread):ncl(iread)),*) this%tbdistort
            iread=iread+1
            read(sometext(ncf(iread):ncl(iread)),*) this%tbgaussweight
            iread=iread+1
            read(sometext(ncf(iread):ncl(iread)),*) this%tbgausswidth
            do while (iread.lt.nwords.and..not.this%hastbintra)
               iread=iread+1
               this%hastbintra=sometext(ncf(iread):ncl(iread)).eq.'anglesintra'
               if (.not.this%hastbintra) then
                  this%ntbangles=this%ntbangles+1
                  call reallocate1d_r(this%tbangles,this%ntbangles)
                  read(sometext(ncf(iread):ncl(iread)),*) this%tbangles(this%ntbangles)
               end if
            end do
         else if(this%hastbintra) then
            read(sometext(ncf(iread):ncl(iread)),*) this%tbintraweight
            iread=iread+1
            read(sometext(ncf(iread):ncl(iread)),*) this%tbintradistort
            iread=iread+1
            read(sometext(ncf(iread):ncl(iread)),*) this%tbintragaussweight
            iread=iread+1
            read(sometext(ncf(iread):ncl(iread)),*) this%tbintragausswidth
            do while(iread.lt.nwords-1)
               iread=iread+1
               this%ntbintraangles=this%ntbintraangles+1
               call reallocate1d_c(this%tbintraats,3,this%ntbintraangles)
               call reallocate1d_r(this%tbintraangles,this%ntbintraangles)
               this%tbintraats(this%ntbintraangles)=sometext(ncf(iread):ncl(iread))
               iread=iread+1
               read(sometext(ncf(iread):ncl(iread)),*) this%tbintraangles(this%ntbintraangles)
            end do
         end if            
      end do
      write(6,*) 'Coordination constraint ',this%atc,' with coordination number ' &
         ,this%nc,' created. ',(this%ats(i),' ',this%ds(i),i=1,this%ns)
      if(this%hastb) then
         write(6,*) ' 3-body angles: ',this%ntbangles,(this%tbangles(i),i=1,this%ntbangles)
         piconv=4.0*atan(1.0)/180.0
         call reallocate1d_r(this%tbcosines,this%ntbangles)
         do i = 1,this%ntbangles
            this%tbcosines(i)=cos(piconv*this%tbangles(i))
         end do
      end if
      if(this%hastbintra) then
         write(6,*) ' intra 3-body angles: ', this%ntbintraangles,' ', &
            (this%tbintraats(i),' ',this%tbintraangles(i),i=1,this%ntbintraangles)
         piconv=4.0*atan(1.0)/180.0
         call reallocate1d_r(this%tbintracosines,this%ntbintraangles)
         do i = 1,this%ntbintraangles
            this%tbintracosines(i)=cos(piconv*this%tbintraangles(i))
         end do
      end if
   end subroutine init_constraint

   subroutine find_constraints(this,nmol,natom,naf,nal,xadd,yadd,zadd,xa,ya,za,mref, &
         itype,ntype,atype,cfrac)

      use reallocation_routines

      type(coordconstraint), intent(inout) :: this
      integer, intent(in) :: nmol,natom,ntype
      integer, dimension(:), intent(in) :: naf,nal,mref,itype
      real, dimension(:), intent(in) :: xadd,yadd,zadd,xa,ya,za,cfrac
      character(len=3), dimension(:), intent(in) :: atype
      integer :: ia,ja,ic,it,it1,it2,im,jm,i,j,ja1,ja2,jm1,jm2,is,is1,is2,nt,ka,kamin,jamin,it1save,k,km,ja2save,itbstart
      logical :: found
      integer, dimension(:),allocatable :: ishortest
      real, dimension(:),allocatable :: dshortest
      integer  :: nshortest,nbondccsize,itb,sumoutside,nbond,inearbond,sumsame,inearbondsave
      real :: xdif,ydif,zdif,dsq,d12,dsfound,cacfrac,sacfrac,tbondsave,adistsave,dsqmin,dsqtb,dsqtbmin
      real :: piconv,tb,tb1,tb2,dist1,dist2,dist12,dtest,dsave,sumratio
      write(6,*) 'Coordination constraint ',this%atc,' with coordination number ' &
         ,this%nc,' bonded to: ',(this%ats(i),' ',this%ds(i),i=1,this%ns)
! Get the atomic fraction of the central atom type
      cacfrac=0.0
      it=find_type(this%atc,ntype,atype)
      if(it.gt.0) then
         cacfrac=cfrac(it)
      end if
! Sum the atomic fractions of the second atom types
      sacfrac=0.0
      do ia = 1,this%ns
         it=find_type(this%ats(ia),ntype,atype)
         if(it.gt.0) sacfrac=sacfrac+cfrac(it)
      end do
      if(sacfrac.gt.0.0) then
         this%alternatenc=nint(real(this%nc)*cacfrac/sacfrac)
      end if
      !Initialise arrays to store the bonds for each atom
      nbondccsize=max(this%nc,this%alternatenc)
      call reallocate1d_i(this%nbondcc,natom)
      call reallocate1d_i(this%nearnbondcc,natom)
      call reallocate1d_i(this%expecnbondcc,natom)
      call reallocate2d_i(this%iabondcc,nbondccsize,natom)
      call reallocate2d_i(this%wtcc,nbondccsize,natom)
      call reallocate2d_r(this%tbondcc,nbondccsize,natom)
      call reallocate2d_r(this%adistcc,nbondccsize,natom)
!      call reallocate1d_r(this%distcc,nbondccsize)
      this%nbondcc=0 !Set the inital number of bonds for each atom to zero
      this%nearnbondcc=0
      this%expecnbondcc=0 !Expected number of bonds for this coordination constraint
!First job is to prepare a list of atoms of type atc
      ic=0
      do ia = 1,natom
         if(atype(itype(ia)).eq.this%atc) ic=ic+1
      end do
      this%atclistsize=ic
      call reallocate1d_i(this%atclist,this%atclistsize)
!      call reallocate2d_i(this%atslist,this%nc,this%atclistsize)
      ic=0
!   write(6,*) 'line 89'
      do ia = 1,natom
         if(atype(itype(ia)).eq.this%atc) then
            ic=ic+1
            this%atclist(ic)=ia
            this%expecnbondcc(ia)=this%nc
         else !See if this atom is one of the alternate atom types
            it1=1
            found=atype(itype(ia)).eq.this%ats(it1)
            do while (it1.lt.this%ns.and..not.found)
               it1=it1+1
               found=atype(itype(ia)).eq.this%ats(it1)
            end do
            if(found) this%expecnbondcc(ia)=this%alternatenc
         end if
      end do
!   write(6,*) 'line 96'
!Next is to prepare the list of nearest neigbours of atc up to the limit ns.
      call reallocate1d_i(ishortest,this%nc)
      call reallocate1d_r(dshortest,this%nc)
      sumoutside=0
      do ic=1,this%atclistsize
         ia=this%atclist(ic)
         im=mref(ia)
         nshortest=0
         do ja = 1,natom
            jm=mref(ja)
            !First check ja is not in the same molecule as ia, and also is on a different molecule to those on the list of shortest distances
            found=im.eq.jm
            it=0
            do while (it.lt.nshortest.and..not.found)
               it=it+1
               found=jm.eq.mref(ishortest(it))
!               found=ja.eq.ishortest(it)
            end do
            if(.not.found) then
               it=0
               do while (it.lt.this%ns.and..not.found)
                  it=it+1
                  found=this%ats(it).eq.atype(itype(ja))
               end do
               if (found) then
! Get the separation of ia and ja
                  xdif=(xadd(im)+xa(ia))-(xadd(jm)+xa(ja))
                  ydif=(yadd(im)+ya(ia))-(yadd(jm)+ya(ja))
                  zdif=(zadd(im)+za(ia))-(zadd(jm)+za(ja))
                  call minimagec(xdif,ydif,zdif)
                  dsq=xdif*xdif+ydif*ydif+zdif*zdif
                  d12=sqrt(dsq)
                  !Step through the list of shortest distances and if d12 is smaller than one of these signal the fact
                  it=0
                  found=.false.
                  do while (it.lt.nshortest.and..not.found)
                     it=it+1
                     found=d12.lt.dshortest(it)
                  end do
                  if(found) then
                     if(nshortest.lt.this%nc) nshortest=nshortest+1
                     !move all the existing values up 1 step
                     it2=nshortest
!      write(6,*) ia,ja,d12,dshortest(it),it,it2,nshortest
                     do while (it2.gt.it)
                        it1=it2-1
                        ishortest(it2)=ishortest(it1)
                        dshortest(it2)=dshortest(it1)
                        it2=it1
                     end do
                     !Finally insert the new values at it
                     ishortest(it)=ja
                     dshortest(it)=d12
                  else if(nshortest.lt.this%nc) then !Otherwise tag the value on the end of the list, but only if the number in the list is less than the maximum
                     nshortest=nshortest+1
                     ishortest(nshortest)=ja
                     dshortest(nshortest)=d12
                  end if
!        if(ia.eq.1) write(6,*) 'find_constraints: ',ia,' ',atype(itype(ia)),' ',ja,d12 &
!          ,(ishortest(it),dshortest(it),it=1,nshortest)
               end if
            end if
         end do
!   write(6,*) 'line 155'
         !Finally save the list of nearest neighbours and sum the number of atoms outside the specified limit.
         do it=1,this%nc
            ja=ishortest(it)
            !Get the expected bond distance for atom ja
            it1=0
            found=.false.
            do while (it1.lt.this%ns.and..not.found)
               it1=it1+1
               found=atype(itype(ja)).eq.this%ats(it1)
            end do
!   write(6,*) 'line 225 ',ia,ja,atype(itype(ja)),this%ats(it1)
            if(found) then !Just to be safe, only set the atom bond numbers if the atom is correctly recognised.
               if(dshortest(it).gt.this%dsmax(it1)) sumoutside=sumoutside+1!Count the number of atoms outside the specified maximum distance
               dsfound=this%ds(it1)
!   write(6,*) 'line 227 ',ia,ja,atype(itype(ja)),this%ats(it1),dsfound
!               if(this%nbondcc(ia).lt.this%nc.and.dshortest(it).le.this%dsmax(it1)) then !Only set the bonds if they are within range
               this%nbondcc(ia)=this%nbondcc(ia)+1
               if(this%nbondcc(ia).gt.nbondccsize) then
                  nbondccsize=this%nbondcc(ia)
                  call reallocate2d_i(this%iabondcc,nbondccsize,natom)
                  call reallocate2d_i(this%wtcc,nbondccsize,natom)
                  call reallocate2d_r(this%tbondcc,nbondccsize,natom)
                  call reallocate2d_r(this%adistcc,nbondccsize,natom)
               end if
               this%iabondcc(this%nbondcc(ia),ia)=ja
               this%wtcc(this%nbondcc(ia),ia)=1
               this%tbondcc(this%nbondcc(ia),ia)=dsfound
               this%adistcc(this%nbondcc(ia),ia)=dshortest(it)
!               if(itype(ia).ne.itype(ja).and.this%nbondcc(ja).lt.this%alternatenc) then !Only set the alternate bonds if the atoms are of different types. A maximum of this%alternatenc such bonds are allowed
!               if(itype(ia).ne.itype(ja).and.dshortest(it).le.this%dsmax(it1)) then !Only set the alternate bonds if the atoms are of different types. A maximum of this%alternatenc such bonds are allowed
               if(itype(ia).ne.itype(ja)) then !Only set the alternate bonds if the atoms are of different types. A maximum of this%alternatenc such bonds are allowed
                  !Need to insert this bond in the order of shortest - longest
                  !First find the shortest
                  is=0
                  found=.false.
                  do while (is.lt.this%nbondcc(ja).and..not.found)
                     is=is+1
                     found=dshortest(it).lt.this%adistcc(is,ja)
                  end do
                  if(found) then !If a longer value was found, then insert before the longer value. First move all existing values up 1
                     this%nbondcc(ja)=this%nbondcc(ja)+1
                     if(this%nbondcc(ja).gt.nbondccsize) then
                        nbondccsize=this%nbondcc(ja)
                        call reallocate2d_i(this%iabondcc,nbondccsize,natom)
                        call reallocate2d_i(this%wtcc,nbondccsize,natom)
                        call reallocate2d_r(this%tbondcc,nbondccsize,natom)
                        call reallocate2d_r(this%adistcc,nbondccsize,natom)
                     end if
                     is2=this%nbondcc(ja)
                     do while(is2.gt.is)
                        is1=is2-1
                        this%iabondcc(is2,ja)=this%iabondcc(is1,ja)
                        this%wtcc(is2,ja)=this%wtcc(is1,ja)
                        this%tbondcc(is2,ja)=this%tbondcc(is1,ja)
                        this%adistcc(is2,ja)=this%adistcc(is1,ja)
                        is2=is1
                     end do
                     !Then save the new values at is
                     this%iabondcc(is,ja)=ia
                     this%wtcc(is,ja)=1
                     this%tbondcc(is,ja)=dsfound
                     this%adistcc(is,ja)=dshortest(it)
                  else !Otherwise simply add new values at end of list if the list is not full
                     this%nbondcc(ja)=this%nbondcc(ja)+1
                     if(this%nbondcc(ja).gt.nbondccsize) then
                        nbondccsize=this%nbondcc(ja)
                        call reallocate2d_i(this%iabondcc,nbondccsize,natom)
                        call reallocate2d_i(this%wtcc,nbondccsize,natom)
                        call reallocate2d_r(this%tbondcc,nbondccsize,natom)
                        call reallocate2d_r(this%adistcc,nbondccsize,natom)
                     end if
                     this%iabondcc(this%nbondcc(ja),ja)=ia
                     this%wtcc(this%nbondcc(ja),ja)=1
                     this%tbondcc(this%nbondcc(ja),ja)=dsfound
                     this%adistcc(this%nbondcc(ja),ja)=dshortest(it)
                  end if
               end if
            end if
         end do
      end do   
! Check whether the same molecules appear as a neighbour to each molecule
      nt=0
      sumsame=0
      do ia=1,natom
         nt=nt+1
         if(this%nbondcc(ia).gt.0) then
            do j=1,this%nbondcc(ia)
               if(mref(ia).eq.mref(this%iabondcc(j,ia))) then
                  sumsame=sumsame+1
                  write(6,*) 'find_constraints - pass 1> Neighbouring atoms in same molecule: ',ia,this%iabondcc(j,ia)
               end if
            end do
            do j=1,this%nbondcc(ia)-1
               jm=mref(this%iabondcc(j,ia))
               k=j
               do while (k.lt.this%nbondcc(ia))
                  k=k+1
                  km=mref(this%iabondcc(k,ia))
                  if(km.eq.jm) then
                     sumsame=sumsame+1
                     write(6,*) 'find_constraints - pass 1> Neighbouring atoms in same molecule: ' &
         ,ia,this%iabondcc(j,ia),this%iabondcc(k,ia)
                  end if
               end do
            end do
         end if
      end do     
      write(6,*) 'find_constraints - pass 1> ',this%atc,' Molecules listed as neighbours to themselves (should be zero) ' &
            ,sumsame
! If the number of bonds for a given atom is greater than expected, we need to remove the extra atoms from the list and also from the lists of the atoms
! to which they are bonded.
      do ia=1,natom
         !If nbondcc > expecnbondcc it means this atom occurs too often in the bonding lists of the other atoms, so it should be removed
         it=this%nbondcc(ia)
         do while (it.gt.this%expecnbondcc(ia)) 
            ja=this%iabondcc(it,ia)
            !Remove atom ia as neighbour of ja
            it1=this%nbondcc(ja)
            found=ia.eq.this%iabondcc(it1,ja)
            do while(it1.gt.2.and..not.found)
               it1=it1-1
               found=ia.eq.this%iabondcc(it1,ja)
            end do
            if(found) then
               !Save the values for this atom
               tbondsave=this%tbondcc(it1,ja)
               adistsave=this%adistcc(it1,ja)
               !Remove ia from the list of ja
               do it2=it1,this%nbondcc(ja)-1
                  this%iabondcc(it2,ja)=this%iabondcc(it2+1,ja)
                  this%wtcc(it2,ja)=this%wtcc(it2+1,ja)
                  this%tbondcc(it2,ja)=this%tbondcc(it2+1,ja)
                  this%adistcc(it2,ja)=this%adistcc(it2+1,ja)
               end do
               this%nbondcc(ja)=this%nbondcc(ja)-1
            end if
            it=it-1
         end do
         this%nbondcc(ia)=it !Finally reset the reduced number of bonds for ia.
      end do
! Check whether the same molecules appear as a neighbour to each molecule
      nt=0
      sumsame=0
      do ia=1,natom
         nt=nt+1
         if(this%nbondcc(ia).gt.0) then
            do j=1,this%nbondcc(ia)
               if(mref(ia).eq.mref(this%iabondcc(j,ia))) then
                  sumsame=sumsame+1
                  write(6,*) 'find_constraints - pass 2> Neighbouring atoms in same molecule: ',ia,this%iabondcc(j,ia)
               end if
            end do
            do j=1,this%nbondcc(ia)-1
               jm=mref(this%iabondcc(j,ia))
               k=j
               do while (k.lt.this%nbondcc(ia))
                  k=k+1
                  km=mref(this%iabondcc(k,ia))
                  if(km.eq.jm) then
                     sumsame=sumsame+1
                     write(6,*) 'find_constraints - pass 2> Neighbouring atoms in same molecule: ' &
         ,ia,this%iabondcc(j,ia),this%iabondcc(k,ia)
                  end if
               end do
            end do
         end if
      end do     
      write(6,*) 'find_constraints - pass 2> ',this%atc,' Molecules listed as neighbours to themselves (should be zero) ' &
            ,sumsame
! Need to consider the case where an atom is UNDERbonded - can suitable combinations be found??
      do ic=1,this%atclistsize !For simplicity we only look for atoms of type this%atc which are under-bonded - the opposite bonding is set automatically when a bond is found.
         ia=this%atclist(ic)
         im=mref(ia)
         it=this%nbondcc(ia)
         do while (it.lt.this%expecnbondcc(ia))
!            write(6,*) 'Atom ',ia,' type ',atype(itype(ia)),' ',this%nbondcc(ia),this%expecnbondcc(ia) 
            it=it+1
            ! Need to find an atom of one of the types this%ats which is ALSO underbonded.
            ja=0
            found=.false.
            it1=0
            do while(ja.lt.natom.and..not.found)
               ja=ja+1
!               write(6,*) 'Atom ',ja,' ',atype(itype(ja)),' ',this%nbondcc(ja),this%expecnbondcc(ja)
               found=this%nbondcc(ja).lt.this%expecnbondcc(ja).and.mref(ja).ne.mref(ia) !Don't waste time on atoms which are in the same molecule as ia
               if(found) then
                  it1=1
                  found=atype(itype(ja)).eq.this%ats(it1)
                  do while (it1.lt.this%ns.and..not.found)
                     it1=it1+1
                     found=atype(itype(ja)).eq.this%ats(it1)
                  end do
                  if(found) then
                     !Need to check molecule containing ja is not already listed as a neighbour of the molecule containing ia
                     it2=1
                     found=mref(ja).ne.mref(this%iabondcc(it2,ia))
                     do while (it2.lt.this%nbondcc(ia).and.found)
                        it2=it2+1
                        found=mref(ja).ne.mref(this%iabondcc(it2,ia))
                     end do
                  end if
               end if
            end do
            it1save=it1
!            write(6,*) 'Atom ',ja,' type ',atype(itype(ja)),' ',this%ats(it1save),' ',this%nbondcc(ja),this%expecnbondcc(ja) 
            ! If one was found, now look for other possibilities, choosing the one that is nearest to ia
            if(found) then
               !measure the distance from ia to ja
               xdif=xadd(mref(ia))+xa(ia)-xadd(mref(ja))-xa(ja)
               ydif=yadd(mref(ia))+ya(ia)-yadd(mref(ja))-ya(ja)
               zdif=zadd(mref(ia))+za(ia)-zadd(mref(ja))-za(ja)
               call minimagec(xdif,ydif,zdif)
               dsqmin=xdif*xdif+ydif*ydif+zdif*zdif
               jamin=ja
               !Step through the remaining atoms to see if any are closer
               do while (ja.lt.natom)
                  ja=ja+1
                  found=this%nbondcc(ja).lt.this%expecnbondcc(ja).and.mref(ja).ne.mref(ia)
                  if(found) then !Check that this atom is on the list possible second atoms
                     it1=1
                     found=atype(itype(ja)).eq.this%ats(it1)
                     do while (it1.lt.this%ns.and..not.found)
                        it1=it1+1
                        found=atype(itype(ja)).eq.this%ats(it1)
                     end do
                     if(found) then
                        !Need to check molecule containing ja is not already listed as a neighbour of the molecule containing ia
                        it2=1
                        found=mref(ja).ne.mref(this%iabondcc(it2,ia))
                        do while (it2.lt.this%nbondcc(ia).and.found)
                           it2=it2+1
                           found=mref(ja).ne.mref(this%iabondcc(it2,ia))
                        end do
                     end if
                  end if
                  if(found) then
                     xdif=xadd(mref(ia))+xa(ia)-xadd(mref(ja))-xa(ja)
                     ydif=yadd(mref(ia))+ya(ia)-yadd(mref(ja))-ya(ja)
                     zdif=zadd(mref(ia))+za(ia)-zadd(mref(ja))-za(ja)
                     call minimagec(xdif,ydif,zdif)
                     dsq=xdif*xdif+ydif*ydif+zdif*zdif
                     if(dsq.lt.dsqmin) then
                        jamin=ja
                        dsqmin=dsq
                        it1save=it1
                     end if
                  end if
               end do
               this%nbondcc(jamin)=this%nbondcc(jamin)+1
               this%iabondcc(this%nbondcc(jamin),jamin)=ia
               this%wtcc(this%nbondcc(jamin),jamin)=1
               this%tbondcc(this%nbondcc(jamin),jamin)=this%ds(it1save)
               this%adistcc(this%nbondcc(jamin),jamin)=sqrt(dsqmin)
               this%nbondcc(ia)=this%nbondcc(ia)+1
               this%iabondcc(this%nbondcc(ia),ia)=jamin
               this%wtcc(this%nbondcc(ia),ia)=1
               this%tbondcc(this%nbondcc(ia),ia)=this%ds(it1save)
               this%adistcc(this%nbondcc(ia),ia)=this%adistcc(this%nbondcc(jamin),jamin)
               ! Strictly the distance between ia and ja should be inserted in the correct order.
               ! However from this point onwards adistcc is not used, so we won't bother for the time being.
               ! In any case it is likely this atom (ja) is further away than other atoms already in the list.
            end if
         end do
      end do
! Check whether the same molecules appear as a neighbour to each molecule
      nt=0
      sumsame=0
      do ia=1,natom
         nt=nt+1
         if(this%nbondcc(ia).gt.0) then
            do j=1,this%nbondcc(ia)
               if(mref(ia).eq.mref(this%iabondcc(j,ia))) then
                  sumsame=sumsame+1
                  write(6,*) 'find_constraints - pass 3> Neighbouring atoms in same molecule: ',ia,this%iabondcc(j,ia)
               end if
            end do
            do j=1,this%nbondcc(ia)-1
               jm=mref(this%iabondcc(j,ia))
               k=j
               do while (k.lt.this%nbondcc(ia))
                  k=k+1
                  km=mref(this%iabondcc(k,ia))
                  if(km.eq.jm) then
                     sumsame=sumsame+1
                     write(6,*) 'find_constraints - pass 3> Neighbouring atoms in same molecule: ' &
         ,ia,this%iabondcc(j,ia),this%iabondcc(k,ia)
                  end if
               end do
            end do
         end if
      end do     
      write(6,*) 'find_constraints - pass 3> ',this%atc,' Molecules listed as neighbours to themselves (should be zero) ' &
            ,sumsame
! Save the number of nearest neighbours for each atom.
      this%nearnbondcc=this%nbondcc
!      do ia = 1,natom
!         write(6,*) ia,' ',atype(itype(ia)),' ',this%nearnbondcc(ia),this%nbondcc(ia),(this%iabondcc(j,ia),j=1,this%nbondcc(ia))
!      end do
!      stop
      !Assess how well the bonding criterion has been satisfied
      nt=0
      sumratio=0.0
      do ia=1,natom
         nt=nt+1
         if(this%expecnbondcc(ia).eq.this%nbondcc(ia)) then
            sumratio=sumratio+1
         end if
      end do
      if(nt.gt.0) sumratio=100.0*sumratio/real(nt)
      ia=this%atclist(this%atclistsize)
      write(6,*) 'find_constraints: ',ia,' ',atype(itype(ia)),' ', &
            (ishortest(it),' ',atype(itype(ishortest(it))),' ',dshortest(it),it=1,this%nc)
      write(6,*) 'find_constraints: ',this%atc,' Fraction outside specified distance(s) = ', &
            real(sumoutside)/real(this%nc*this%atclistsize)
      write(6,*) 'find_constraints: ',this%atc,' Percentage constraints correct = ', &
            sumratio
      !If these atoms are bonded via tb angles, need to define the corresponding bond distances
      if(this%hastb) then
         piconv=4.0*atan(1.0)/180.0
         do ic=1,this%atclistsize
            ia=this%atclist(ic)
            im=mref(ia)
            i=0
!            write(6,*) ia,this%nbondcc(ia),this%nearnbondcc(ia)
            do while(i.lt.this%nearnbondcc(ia)-1)
               i=i+1
               ja1=this%iabondcc(i,ia)
               jm1=mref(ja1)
! Get the separation of ia and ja1
!               tb1=this%tbondcc(i,ia) !Use the expected separation distance to determine the three-body distance constraint.
               tb1=this%adistcc(i,ia) !Use the actual separation to determine the three-body distance constraint.
!               tb1=max(this%tbondcc(i,ia),this%adistcc(i,ia)) !Use the maximum of dist1 and tb1 to determine the three-body distance constraint. This prevents very small values of dist1 causing problems.
               j=i
               do while (j.lt.this%nearnbondcc(ia))
                  j=j+1
                  ja2=this%iabondcc(j,ia)
                  jm2=mref(ja2)
! Get the separation of ia and ja2
!                  tb2=this%tbondcc(j,ia) !Use the expected separation to determine the three-body distance constraint
                  tb2=this%adistcc(i,ia) !Use the actual separation to determine the three-body distance constraint.
!                  tb2=max(this%tbondcc(j,ia),this%adistcc(j,ia)) !Use the maximum of dist2 and tb2 to determine the three-body distance constraint, as above.
! Get the separation of ja1 and ja2
                  xdif=(xadd(jm1)+xa(ja1))-(xadd(jm2)+xa(ja2))
                  ydif=(yadd(jm1)+ya(ja1))-(yadd(jm2)+ya(ja2))
                  zdif=(zadd(jm1)+za(ja1))-(zadd(jm2)+za(ja2))
                  call minimagec(xdif,ydif,zdif)
                  dist12=xdif*xdif+ydif*ydif+zdif*zdif
!                  dsq=dist1*dist1+dist2*dist2-2.0*dist1*dist2*this%tbcosines(1)!If there are more than one angle we can check which one this particular pair is closer to.
                  dsq=tb1*tb1+tb2*tb2-2.0*tb1*tb2*this%tbcosines(1)!If there are more than one angle we can check which one this particular pair is closer to.
                  tb=sqrt(dsq)
                  dsave=abs(dist12-dsq)
                  itb=1
! If more than one angle specified, find the closest to dist12
                  do while (itb.lt.this%ntbangles)
                     itb=itb+1
!                     dsq=dist1*dist1+dist2*dist2-2.0*dist1*dist2*this%tbcosines(itb)
                     dsq=tb1*tb1+tb2*tb2-2.0*tb1*tb2*this%tbcosines(itb)!If there are more than one angle we can check which one this particular pair is closer to.
                     dtest=abs(dist12-dsq)
                     if(dtest.lt.dsave) then
                        dsave=dtest
                        tb=sqrt(dsq)
                     end if
                  end do
                  !It is possible atom ja2 already occurs in the coordination shell of ja1 (due to direct bonding to different centre atoms).
                  !To avoid this, we replace the distance between ja1 and ja2 if it is shorter, but don't increase the number of bonds
                  itb=this%nearnbondcc(ja1)
                  found=.false.
                  do while (itb.lt.this%nbondcc(ja1).and..not.found)
                     itb=itb+1
                     found=ja2.eq.this%iabondcc(itb,ja1)
                  end do
                  if(found) then
                     tb=min(abs(this%tbondcc(itb,ja1)),tb)
                     this%tbondcc(itb,ja1)=-tb
                     !Need to repeat this for ja2 - first find which atom number ja1 is
                     itb=this%nearnbondcc(ja2)
                     found=.false.
                     do while (itb.lt.this%nbondcc(ja2).and..not.found)
                        itb=itb+1
                        found=ja1.eq.this%iabondcc(itb,ja2)
                     end do
                     if(found) this%tbondcc(itb,ja2)=-tb
                  else
                     !Otherwise we need to check the molecule containing ja2 does not already occur in the coordination shell of ja1
                     itb=1
                     found=mref(this%iabondcc(itb,ja1)).ne.mref(ja2)
                     do while (itb.lt.this%nbondcc(ja1).and.found)
                        itb=itb+1
                        found=mref(this%iabondcc(itb,ja1)).ne.mref(ja2)
                     end do
                     !Repeat for ja2
                     itb=1
                     found=found.and.mref(this%iabondcc(itb,ja2)).ne.mref(ja1)
                     do while (itb.lt.this%nbondcc(ja2).and.found)
                        itb=itb+1
                        found=found.and.mref(this%iabondcc(itb,ja2)).ne.mref(ja1)
                     end do
                     if(found) then
                        this%nbondcc(ja1)=this%nbondcc(ja1)+1
                        this%nbondcc(ja2)=this%nbondcc(ja2)+1  
                        nbondccsize=max(nbondccsize,this%nbondcc(ja1),this%nbondcc(ja2))
                        call reallocate2d_i(this%iabondcc,nbondccsize,natom)
                        call reallocate2d_i(this%wtcc,nbondccsize,natom)
                        call reallocate2d_r(this%tbondcc,nbondccsize,natom)
                        call reallocate2d_r(this%adistcc,nbondccsize,natom)
                        this%iabondcc(this%nbondcc(ja1),ja1) = ja2
                        this%iabondcc(this%nbondcc(ja2),ja2) = ja1
                        this%wtcc(this%nbondcc(ja1),ja1) = 2
                        this%wtcc(this%nbondcc(ja2),ja2) = 2
                        this%tbondcc(this%nbondcc(ja1),ja1) = -tb!Negative will indicate it has a different weighting to direct bonds.
                        this%tbondcc(this%nbondcc(ja2),ja2) = -tb
                        this%adistcc(this%nbondcc(ja1),ja1) = sqrt(dist12)!Negative will indicate it has a different weighting to direct bonds.
                        this%adistcc(this%nbondcc(ja2),ja2) = this%adistcc(this%nbondcc(ja1),ja1)
                     end if
                  end if
               end do
            end do
         end do
! Check whether the same molecules appear as a neighbour to each molecule
         nt=0
         sumsame=0
         do ia=1,natom
            nt=nt+1
            if(this%nbondcc(ia).gt.0) then
               do j=1,this%nbondcc(ia)
                  if(mref(ia).eq.mref(this%iabondcc(j,ia))) then
                     sumsame=sumsame+1
                     write(6,*) 'find_constraints - pass 4> Neighbouring atoms in same molecule: ',ia,this%iabondcc(j,ia)
                  end if
               end do
               do j=1,this%nbondcc(ia)-1
                  jm=mref(this%iabondcc(j,ia))
                  k=j
                  do while (k.lt.this%nbondcc(ia))
                     k=k+1
                     km=mref(this%iabondcc(k,ia))
                     if(km.eq.jm) then
                        sumsame=sumsame+1
                        write(6,*) 'find_constraints - pass 4> Neighbouring atoms in same molecule: ' &
         ,ia,this%iabondcc(j,ia),this%iabondcc(k,ia)
                     end if
                  end do
               end do
            end if
         end do     
         write(6,*) 'find_constraints - pass 4> ',this%atc,' Molecules listed as neighbours to themselves (should be zero) ' &
            ,sumsame
      end if
      write(6,*) 'Got to tbintra',this%hastbintra
      if(this%hastbintra) then !For intra three-body angles, the object is to identify atoms of type tbintra ats bound (same molecule) to ia
! and identify a neighbour of ia with which they are to form the specified angle. If more than one of these atoms is specified (nbtintraangles>1), then 
! a different neighbour must be chosen for each.
!         do ia = 1,natom
!            write(6,*) ia,' ',atype(itype(ia)),' ',this%nearnbondcc(ia),this%nbondcc(ia),(this%iabondcc(j,ia),j=1,this%nbondcc(ia))
!         end do
         piconv=4.0*atan(1.0)/180.0
         do ic=1,this%atclistsize
            ia=this%atclist(ic)
            im=mref(ia)
            ja1=naf(im)-1
            !Search the atoms of molecule im to finds which match tbintraats. We expect one constraint for each listed intra angle, but only up to the number of atoms of listed types in molecule im
            do while (ja1.lt.nal(im))
               ja1=ja1+1
               jm1=im
! Get the separation of ia and ja1
               xdif=xa(ja1)-xa(ia)
               ydif=ya(ja1)-ya(ia)
               zdif=za(ja1)-za(ia)
               call minimagec(xdif,ydif,zdif)
               dsq=xdif*xdif+ydif*ydif+zdif*zdif
               tb1=sqrt(dsq)
               itb=0
               found=ja1.ne.ia
               !Step through the intra angles
               do while (itb.lt.this%ntbintraangles.and.found)
                  itb=itb+1
                  found=atype(itype(ja1)).eq.this%tbintraats(itb)
                  if(found) then
                     !For the first valid ja1 atom, search the near-neighbour bonds of ia and find that which is the closest to the expected separation from ja1
                     inearbond=0
                     dsave=0.0
                     jamin=0
                     dsqmin=0.0
                     dsqtbmin=0.0
                     do while (inearbond.lt.this%nearnbondcc(ia))
                        inearbond=inearbond+1
                        ja2=this%iabondcc(inearbond,ia)
                        !Check that ja2 is not already a neighbour of ja1
                        it1=1
                        found=ja2.ne.this%iabondcc(it1,ja1)
                        do while (it1.lt.this%nbondcc(ja1).and.found)
                           it1=it1+1
                           found=ja2.ne.this%iabondcc(it1,ja1)
                        end do
                        if(found) then               
                           jm2=mref(ja2)
                           ! Get the separation of ia and ja2
                           tb2=this%adistcc(inearbond,ia)
                           ! Get expected separation of ja1 and ja2
                           dsqtb=tb1*tb1+tb2*tb2-2.0*tb1*tb2*this%tbintracosines(itb)
                           ! Get the separation of ja1 and ja2
                           xdif=xadd(jm1)+xa(ja1)-(xadd(jm2)+xa(ja2))
                           ydif=yadd(jm1)+ya(ja1)-(yadd(jm2)+ya(ja2))
                           zdif=zadd(jm1)+za(ja1)-(zadd(jm2)+za(ja2))
                           call minimagec(xdif,ydif,zdif)
                           dsq=xdif*xdif+ydif*ydif+zdif*zdif
                           dtest=abs(dsq-dsqtb)
                           if(dtest.lt.dsave.or.dsave.eq.0.0) then
                              dsave=dtest
                              jamin=ja2
                              dsqmin=dsq
                              dsqtbmin=dsqtb
                           end if
                        end if
                     end do
                     if(jamin.ne.0) then !At least one valid bond was found
                        ja2=jamin
                        tb=sqrt(dsqtbmin)
                        this%nbondcc(ja1)=this%nbondcc(ja1)+1
                        this%nbondcc(ja2)=this%nbondcc(ja2)+1  
                        nbondccsize=max(nbondccsize,this%nbondcc(ja1),this%nbondcc(ja2))
                        call reallocate2d_i(this%iabondcc,nbondccsize,natom)
                        call reallocate2d_i(this%wtcc,nbondccsize,natom)
                        call reallocate2d_r(this%tbondcc,nbondccsize,natom)
                        call reallocate2d_r(this%adistcc,nbondccsize,natom)
                        this%iabondcc(this%nbondcc(ja1),ja1) = ja2
                        this%iabondcc(this%nbondcc(ja2),ja2) = ja1
                        this%wtcc(this%nbondcc(ja1),ja1) = 3
                        this%wtcc(this%nbondcc(ja2),ja2) = 3
                        this%tbondcc(this%nbondcc(ja1),ja1) = -tb
                        this%tbondcc(this%nbondcc(ja2),ja2) = -tb
                        this%adistcc(this%nbondcc(ja1),ja1) = sqrt(dsqmin)!Negative will indicate it has a different weighting to direct bonds.
                        this%adistcc(this%nbondcc(ja2),ja2) = this%adistcc(this%nbondcc(ja1),ja1)
                     end if
                  end if
               end do
            end do
         end do
!         do ia = 1,natom
!            write(6,*) ia,' ',atype(itype(ia)),' ',this%nearnbondcc(ia),this%nbondcc(ia),(this%iabondcc(j,ia),j=1,this%nbondcc(ia))
!         end do
!         stop
! Check whether the same molecules appear as a neighbour to each molecule
         nt=0
         sumsame=0
         do ia=1,natom
            nt=nt+1
            if(this%nbondcc(ia).gt.0) then
               do j=1,this%nbondcc(ia)
                  if(mref(ia).eq.mref(this%iabondcc(j,ia))) then
                     sumsame=sumsame+1
!                     write(6,*) 'find_constraints - pass 5> Neighbouring atoms in same molecule: ',ia,this%iabondcc(j,ia)
                  end if
               end do
               do j=1,this%nbondcc(ia)-1
                  jm=mref(this%iabondcc(j,ia))
                  k=j
                  do while (k.lt.this%nbondcc(ia))
                     k=k+1
                     km=mref(this%iabondcc(k,ia))
                     if(km.eq.jm) then
                        sumsame=sumsame+1
!                        write(6,*) 'find_constraints - pass 5> Neighbouring atoms in same molecule: ' &
!         ,ia,this%iabondcc(j,ia),this%iabondcc(k,ia)
                     end if
                  end do
               end do
            end if
         end do     
         write(6,*) 'find_constraints - pass 5> ',this%atc,' Molecules listed as neighbours to themselves (may not be zero) ' &
            ,sumsame
      end if
   end subroutine find_constraints

   function find_type(alabel,ntype,atype)

      integer :: find_type
      character(len=3), intent(in) :: alabel
      integer, intent(in) :: ntype
      character(len=3), dimension(:), intent(in) :: atype
      logical found
      integer i

      found=.false.
      find_type=0
      i=0
      do while (i < ntype.and..not.found)
         i=i+1
         found=alabel.eq.atype(i)
      end do
      if(found) find_type=i
      return
   end function find_type

   function en_coordcon(this,ia,xadd,yadd,zadd,xa,ya,za,mref,itype,qtype,npairref,amassr,tethered)
      !Calculates the coordination constraint energy of atom ia in coordination constraint this
      real :: en_coordcon
      type(coordconstraint), intent(inout) :: this
      integer, intent(in) :: ia      
      real, dimension(:), intent(in) :: xadd,yadd,zadd,xa,ya,za,amassr
      integer, dimension(:), intent(in) :: mref,itype,qtype
      integer, dimension(:,:), intent(in) :: npairref
      logical, dimension(:), intent(in) :: tethered
      integer :: im,j,ja,jm,nbond    
      real :: en_atom_inter_scale_coord,x1,y1,z1,x2,y2,z2,tb,wtm,wtms,distort &
       ,xdif,ydif,zdif,dsq,dtest,en,gweight,gwidth
      logical :: bothtethered

      en_coordcon = 0.0
      nbond=this%nbondcc(ia)
      if(nbond.gt.0) then
         im=mref(ia)
         x1=xadd(im)+xa(ia)
         y1=yadd(im)+ya(ia)
         z1=zadd(im)+za(ia)
         do j=1,nbond
            ja=this%iabondcc(j,ia)
            jm=mref(ja)
            bothtethered=tethered(im).and.tethered(jm).and.qtype(ia).eq.1.and.qtype(ja).eq.1 !Don't calculate pair energies if both molecules are tethered and both atoms are Q-type
            if(.not.bothtethered) then
               x2=xadd(jm)+xa(ja)
               y2=yadd(jm)+ya(ja)
               z2=zadd(jm)+za(ja)
               xdif=x2-x1
               ydif=y2-y1
               zdif=z2-z1
! apply minimum image convention
               call minimagec(xdif,ydif,zdif)
               dsq=xdif*xdif+ydif*ydif+zdif*zdif
               dtest=sqrt(dsq)
               tb=this%tbondcc(j,ia)
               wtms=this%cctemp
               distort=this%ccdistort
               gweight=this%ccgaussweight
               gwidth=this%ccgausswidth
               if(this%wtcc(j,ia).eq.2) then
                  wtms=this%tbweight
                  distort=this%tbdistort
                  gweight=this%tbgaussweight
                  gwidth=this%tbgausswidth
               else if (this%wtcc(j,ia).eq.3) then
                  wtms=this%tbintraweight
                  distort=this%tbintradistort
                  gweight=this%tbintragaussweight
                  gwidth=this%tbintragausswidth
               end if
               wtm=amassr(npairref(itype(ia),itype(ja)))/abs(tb)
               en=en_atom_inter_scale_coord(dtest,tb,wtm,wtms,distort,gweight,gwidth)
!     write(6,*) 'en_coordcon> ',ia,ja,tb,dtest,wtm,wtms,width,en
               en_coordcon=en_coordcon+en
            end if
         end do
!         if (this%hastb) then !Calculate the three-body energies
!            do j=1,this%nbondcc(ia)-1
!               ja=this%iabondcc(j,ia)
!               jm=mref(ja)
!               x2=xadd(jm)+xa(ja)
!               y2=yadd(jm)+ya(ja)
!               z2=zadd(jm)+za(ja)
!               tb12=this%distcc(j)
!               do k=j+1,this%nbondcc(ia)
!                  ka=this%iabondcc(k,ia)
!                  km=mref(ka)
!                  bothtethered=tethered(jm).and.tethered(km).and.qtype(ja).eq.1.and.qtype(ka).eq.1 !Don't calculate pair energies if both molecules are tethered and both atoms are Q-type
!                  if(.not.bothtethered) then
!                     x3=xadd(km)+xa(ka)
!                     y3=yadd(km)+ya(ka)
!                     z3=zadd(km)+za(ka)
!                     tb13=this%distcc(k)
!                     xdif=x3-x2
!                     ydif=y3-y2
!                     zdif=z3-z2
! apply minimum image convention
!                     call minimagec(xdif,ydif,zdif)
!                     dsq=xdif*xdif+ydif*ydif+zdif*zdif
!                     tb23=sqrt(dsq)
!                     if(tb23.le.abs(tb12-tb13)) then
!                        cost=0.0
!                     else if(tb23.ge.(tb12+tb13)) then
!                        cost=2.0
!                     else
!                        cost=1.0-(tb12*tb12+tb13*tb13-tb23*tb23)/(2.0*tb12*tb13)
!                     end if
! Reverse sign of cost because cosine goes from 1 to -1 with increasing angle
!                     cdist=1.0-this%tbcosines(1)
!                     csave=abs(cost-cdist)
!                     itb=1
!                     do while (itb.lt.this%ntbangles)
!                        itb=itb+1
!                        ctest=abs(cost-(1.0-this%tbcosines(itb)))
!                        if(ctest.lt.csave) then
!                           csave=ctest
!                           cdist=1.0-this%tbcosines(itb)
!                        end if
!                     end do
!                     wtms=this%tbweight
!                     distort=this%tbdistort
!                     wtm=amassr(npairref(itype(ia),itype(ja)))
!                     en=en_atom_inter_scale_coord(cost,cdist,wtm,wtms,distort)
!   if(ia.eq.1) write(6,*) ia,ja,ka,tb12,tb13,tb23,cost,cdist,wtm,wtms,distort,en
!                     en_coordcon=en_coordcon+en
!                  end if
!               end do
!            end do
!         end if
      end if
   end function en_coordcon

   function en_coordcon_total(this,natom,xadd,yadd,zadd,xa,ya,za,mref,itype &
         ,qtype,npairref,amassr,tethered)
      !Calculates the total coordination constraint energy of all atoms in coordination constraint this
      real :: en_coordcon_total
      type(coordconstraint), intent(in) :: this
      integer, intent(in) :: natom
      real, dimension(:), intent(in) :: xadd,yadd,zadd,xa,ya,za,amassr
      integer, dimension(:), intent(in) :: mref,itype,qtype
      integer, dimension(:,:), intent(in) :: npairref
      logical, dimension(:), intent(in) :: tethered
      integer :: ia,j,ja,im,jm      
      real :: en_atom_inter_scale_coord,x1,y1,z1,x2,y2,z2,tb,wtm,wtms &
         ,xdif,ydif,zdif,dtest,dsq,distort
      logical :: bothtethered
      

      en_coordcon_total = 0.0
      do ia=1,natom
         if(this%nbondcc(ia).gt.0) then
            im=mref(ia)
            x1=xadd(im)+xa(ia)
            y1=yadd(im)+ya(ia)
            z1=zadd(im)+za(ia)
            do j=1,this%nbondcc(ia)
               ja=this%iabondcc(j,ia)
               jm=mref(ja)
!               if(ja.gt.ia) then ! Prevent double counting in sum of all coordination constraint energies
                  tb=this%tbondcc(j,ia)
                  wtms=this%cctemp
                  distort=this%ccdistort
                  if(tb.lt.0.0) then
                     tb=abs(tb)
                     wtms=wtms*this%tbweight
                     distort=this%tbdistort
                  end if
                  wtm=amassr(npairref(itype(ia),itype(ja)))/tb
                  distort=sqrt(wtms)*distort
                  x2=xadd(jm)+xa(ja)
                  y2=yadd(jm)+ya(ja)
                  z2=zadd(jm)+za(ja)
                  xdif=x2-x1
                  ydif=y2-y1
                  zdif=z2-z1
! apply minimum image convention
                  call minimagec(xdif,ydif,zdif)
                  dsq=xdif*xdif+ydif*ydif+zdif*zdif
                  dtest=sqrt(dsq)
                  bothtethered=tethered(im).and.tethered(jm).and.qtype(ia).eq.1.and.qtype(ja).eq.1 !Don't calculate pair energies if both molecules are tethered and both atoms are Q-type
                  if(.not.bothtethered) then
                     en_coordcon_total=en_coordcon_total+en_atom_inter_scale_coord(dtest,tb,wtm,wtms,distort)
                  end if
!               end if
            end do
         end if
      end do
   end function en_coordcon_total

end module coordconstraint_mod
            
module coordconstraints_mod

   use coordconstraint_mod

   implicit none
   type (coordconstraint), dimension(:), allocatable :: coordconstraints
   integer :: ncoordconstraints
   real, dimension(:), allocatable                    :: ergold_coord,ergnew_coord!,sume_coord_dist,en_mole_coord_dist

contains

    subroutine reallocate1d_cc(array,newsize)

        !reallocates an array, new values are left uninitialised
        type(coordconstraint), dimension(:), intent(inout), allocatable :: array
        integer             :: oldsize
        integer, intent(in) :: newsize
        type(coordconstraint), dimension(:), allocatable :: arraysave
        integer :: i
        
        if(allocated(array)) then
            oldsize=size(array)
        else
            oldsize=0
        end if
        !Do nothing if array is already the requested size! 
        if(oldsize.eq.newsize) return
        if(oldsize.gt.0) then
            allocate(arraysave(oldsize))
            i=0
            do while (i.lt.oldsize)
                i=i+1
                arraysave(i)=array(i)
            end do
        end if
        if(allocated(array)) deallocate(array)
        if(newsize.gt.0) then
            allocate(array(newsize))
            if(oldsize.gt.0.and.allocated(arraysave)) then
                i=0
                do while (i.lt.min(oldsize,newsize))
                    i=i+1
                    array(i)=arraysave(i)
                end do
            end if
        end if
        if(allocated(arraysave)) deallocate(arraysave)
        return
    end subroutine reallocate1d_cc

   subroutine initialise_coords(nmol)

      use reallocation_routines

      integer, intent(in) :: nmol

! Initialise energy storage arrays
      call reallocate1d_r(ergold_coord,nmol)
      call reallocate1d_r(ergnew_coord,nmol)

   end subroutine initialise_coords

   function en_mole_coord(im,naf,nal,xadd,yadd,zadd,xa,ya,za,mref,itype,qtype,npairref,amassr,tethered)
      !Total energy of molecule im due all coordination constraints
      real                                   :: en_mole_coord,en
      integer, intent(in)                    :: im
      integer, dimension(:), intent(in)      :: naf,nal,mref,itype,qtype
      integer, dimension(:,:), intent(in)    :: npairref
      real, dimension(:), intent(in)         :: xadd,yadd,zadd,xa,ya,za,amassr
      logical, dimension(:), intent(in)      :: tethered
      integer                                :: ia,ic

!      write(6,*) 'en_mole_coord> ',ncoordconstraints,im
      en_mole_coord = 0.0
      do ic=1,ncoordconstraints
!         write(6,*) 'en_mole_coord> ',ic,im,naf(im),nal(im)
         do ia=naf(im),nal(im)
!            write(6,*) 'en_mole_coord> ',ic,im,ia,coordconstraints(ic)%cctemp,coordconstraints(ic)%ccdistort &
!               ,coordconstraints(ic)%nbondcc(ia)
            if(coordconstraints(ic)%nbondcc(ia).gt.0) then
               en=en_coordcon(coordconstraints(ic),ia,xadd,yadd,zadd,xa,ya,za,mref,itype,qtype,npairref,amassr,tethered)
!         write(6,*) 'en_mole_coord>',ic,im,ia,en
               en_mole_coord=en_mole_coord+en
            end if
         end do
      end do
   end function en_mole_coord

   function en_coord_total(natom,xadd,yadd,zadd,xa,ya,za,mref,itype,qtype,npairref,amassr,tethered)
      !Total energy of all atoms due all coordination constraints
      real                                   :: en_coord_total
      integer, intent(in)                    :: natom
      real, dimension(:), intent(in)         :: xadd,yadd,zadd,xa,ya,za,amassr
      integer, dimension(:), intent(in)      :: mref,itype,qtype
      integer, dimension(:,:), intent(in)    :: npairref
      logical, dimension(:), intent(in)      :: tethered
      integer                                :: ic
      en_coord_total = 0.0
      do ic=1,ncoordconstraints
         en_coord_total=en_coord_total+en_coordcon_total(coordconstraints(ic),natom,xadd,yadd,zadd,xa,ya,za,mref,itype &
         ,qtype,npairref,amassr,tethered)
      end do
   end function en_coord_total

   function en_change_mole_coord(im,naf,nal,xadd,yadd,zadd,xa,ya,za,mref,itype &
         ,qtype,npairref,amassr,tethered)
      !Change in total energy of molecule im due all coordination constraints
      real                                   :: en_change_mole_coord
      integer, intent(in)                    :: im
      integer, dimension(:), intent(in)      :: naf,nal,mref,itype,qtype
      integer, dimension(:,:), intent(in)    :: npairref
      real, dimension(:), intent(in)         :: xadd,yadd,zadd,xa,ya,za,amassr
      logical, dimension(:), intent(in)      :: tethered

      en_change_mole_coord=0
      if(ncoordconstraints.gt.0) then
         ergnew_coord(im)=en_mole_coord(im,naf,nal,xadd,yadd,zadd,xa,ya,za,mref,itype,qtype,npairref,amassr,tethered)
         en_change_mole_coord=ergnew_coord(im)-ergold_coord(im)
      end if
   end function en_change_mole_coord

   subroutine save_coord_energies(im)

      integer, intent(in) :: im

      if(ncoordconstraints.gt.0) then
         ergold_coord(im)=ergnew_coord(im)
      end if

   end subroutine save_coord_energies

   function check_coord_bonded(ia,ja,mref,itype)

      logical :: check_coord_bonded
      integer, intent(in) :: ia,ja
      integer, dimension(:), intent(in) :: mref,itype
      integer :: ic,jb,nbond

      check_coord_bonded=.false.
      ic=0
      do while (ic.lt.ncoordconstraints.and..not.check_coord_bonded)
         ic=ic+1
         jb=0
         nbond=coordconstraints(ic)%nbondcc(ia)
!         nbond=min(coordconstraints(ic)%nbondcc(ia),coordconstraints(ic)%expecnbondcc(ia))
         do while (jb.lt.nbond.and..not.check_coord_bonded)
            jb=jb+1
!            check_coord_bonded=ja.eq.coordconstraints(ic)%iabondcc(jb,ia)
            check_coord_bonded=mref(ja).eq.mref(coordconstraints(ic)%iabondcc(jb,ia))
         end do
         !Now check the reverse if atom types are different
         if(itype(ia).ne.itype(ja).and..not.check_coord_bonded) then
            jb=0
            nbond=coordconstraints(ic)%nbondcc(ja)
!            nbond=min(coordconstraints(ic)%nbondcc(ja),coordconstraints(ic)%expecnbondcc(ja))
            do while (jb.lt.nbond.and..not.check_coord_bonded)
               jb=jb+1
!               check_coord_bonded=ia.eq.coordconstraints(ic)%iabondcc(jb,ja)
               check_coord_bonded=mref(ia).eq.mref(coordconstraints(ic)%iabondcc(jb,ja))
            end do
         end if
      end do
   end function check_coord_bonded

   subroutine writecoordcons(natom,itype,atype)

      integer, intent(in) :: natom
      integer, dimension(:), intent(in) :: itype
      character(len=3), dimension(:), intent(in) :: atype
      integer ia,ic,j,nbond,nearnbond
      character(len=256) :: formatline

      open(29,file='coordcon.txt',status='unknown')
!      do ic=1,ncoordconstraints
!         write(6,*) 'writecoordcon> ' &
!            ,size(coordconstraints(ic)%iabondcc,1),size(coordconstraints(ic)%iabondcc,2) &
!            ,size(coordconstraints(ic)%wtcc,1),size(coordconstraints(ic)%wtcc,2) &
!            ,size(coordconstraints(ic)%tbondcc,1),size(coordconstraints(ic)%tbondcc,2) &
!            ,size(coordconstraints(ic)%adistcc,1),size(coordconstraints(ic)%adistcc,2)
!      end do
      do ia = 1,natom
         do ic=1,ncoordconstraints
            nbond=coordconstraints(ic)%nbondcc(ia)
            nearnbond=coordconstraints(ic)%nearnbondcc(ia)
!            write(6,*) 'writecoordcon> ',ia,atype(itype(ia)),ic,nbond
            if(nbond.gt.0) then
               write(formatline,'(a,i2.2,a)') '(i6,1x,a,1x,i2,1x,i2,',nbond,'(i6,1x,a,1x,e10.3))'
!               write(6,'(i5,1x,a)') nbond,formatline
!               write(6,*) ia,ic,nbond
!               do j=1,nbond
!                  write(6,*) j,coordconstraints(ic)%iabondcc(j,ia)
!                  write(6,*) j,itype(coordconstraints(ic)%iabondcc(j,ia))
!                  write(6,*) j,atype(itype(coordconstraints(ic)%iabondcc(j,ia))) 
!                  write(6,*) j,coordconstraints(ic)%wtcc(j,ia)
!                  write(6,*) j,coordconstraints(ic)%tbondcc(j,ia)
!                  write(6,*) j,coordconstraints(ic)%adistcc(j,ia)
!               end do
               write(29,formatline) ia,atype(itype(ia)),nearnbond,nbond &
                     ,(coordconstraints(ic)%iabondcc(j,ia),atype(itype(coordconstraints(ic)%iabondcc(j,ia))) &
                     ,coordconstraints(ic)%tbondcc(j,ia),j=1,nbond)
!                     ,coordconstraints(ic)%adistcc(j,ia),j=1,nbond)
            else
               write(29,'(i6,1x,a,1x,i2)') ia,atype(itype(ia)),nbond 
            end if
         end do
      end do
      close(29)                
   end subroutine writecoordcons

end module coordconstraints_mod




