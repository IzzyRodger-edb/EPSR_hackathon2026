      module threebody_mod

      integer                                            :: nthreeb,nthreebangsmax,nthreebneighbmolmax !Number of three-body groups, maximum number of three-body atoms and angles
      character(len=3), dimension(:,:), allocatable      :: threebgrp !Three atom types to form the group (3 per group).
      real, dimension(:), allocatable                    :: threebamp! Three-body potential amplitude
      real, dimension(:,:), allocatable                  :: threebrmax !Maximum distances for atom pairs (2 per group). Initially the second will be set equal to the first.
      real, dimension(:), allocatable                    :: threebrwidth,threebawidth !Nominal width of three body potential function in radius and angel (deg.)
      integer, dimension(:), allocatable                 :: nthreebangs !Number of angles in group
      real, dimension(:,:), allocatable                  :: threebangs,threebamps !List of requested angles and amplitudes in group.

      integer, dimension(:), allocatable                 :: nthreebneighbmol !Number of neighbouring molecules, for which three-body energies can be defined
      integer, dimension(:,:), allocatable               :: threebatom !List of atom numbers that can form a three-body potential
      integer, dimension(:,:), allocatable               :: threebneighbmol !List of neighbours for each molecule in a triplet
      integer, dimension(:,:,:), allocatable             :: threebatomneighb !List of acceptable neighbours of each atom in threebatom   
      integer, dimension(:,:), allocatable               :: threebsamemol ! Array of 3 1s or 0s to say whether atoms 1 and 2, 2 and 3, and 1 and 3 are on same molecule. 1 means they are on the same molecule.
! The following arrays are used to make lists of atoms of the specified types

      integer                                            :: ltyp1,ltyp2,ltyp3
      integer, dimension(:), allocatable                 :: atom1list,atom2list,atom3list
      integer                                            :: natom1list,natom2list,natom3list

! The following arrays list the number of triplets for each molecule, the list of triplets, and the first entry in this list for each molecule
      integer, dimension(:,:), allocatable               :: ntripletsinmol
      integer, dimension(:,:,:), allocatable             :: tripletlist
      integer, dimension(:,:), allocatable               :: firsttripletinmol
      integer                                            :: ntripletsinmolmax

! The following arrays store the triplet energies of each molecule

      real, dimension(:), allocatable                    :: ergold_triplet,ergnew_triplet,sume_triplet_dist,en_mole_triplet_dist
      real, dimension(:), allocatable                    :: consthreeb !Stores the constant for the potential for each triplet distribution
      
      contains

      subroutine find_triplets(nmol,natom,xadd,yadd,zadd,xa,ya,za,mref,itype,ntype,type)

      use reallocation_routines

      real, dimension(:)  :: xadd,yadd,zadd,xa,ya,za
      integer :: nmol,natom,ntype,im1,im2,im3,ia1,ia2,ia3
      integer, dimension(:) :: mref,itype
      character(len=3), dimension(:) :: type

      logical foundia1,foundia3,found

! Initialise certain arrays

      if(allocated(ergold_triplet)) deallocate(ergold_triplet,ergnew_triplet)
      allocate(ergold_triplet(nmol),ergnew_triplet(nmol))

      if(nthreeb.le.0) return

      call reallocate2d_i(ntripletsinmol,nmol,nthreeb)
      call reallocate2d_i(firsttripletinmol,nmol,nthreeb)
      call reallocate1d_r(sume_triplet_dist,nthreeb)
      call reallocate1d_r(en_mole_triplet_dist,nthreeb)
      ntripletsmax=10000
      call reallocate3d_i(tripletlist,3,ntripletsmax,nthreeb)
! Arrays which will list the triplet neighbours of each molecule
      nthreebneighbmolmax=4
      call reallocate1d_i(nthreebneighbmol,nmol)
      call reallocate2d_i(threebneighbmol,nthreebneighbmolmax,nmol)
! Set the initial values
      do im=1,nmol
         nthreebneighbmol(im)=0
      end do

! First identify the type numbers of the three specified atoms

      do id = 1,nthreeb
         do im1=1,nmol
            ntripletsinmol(im1,id)=0
            firsttripletinmol(im1,id)=0
         end do
         ntriplets=0
         ifound1=0
         ifound2=0
         ifound3=0
         it=0
         do while(it.lt.ntype.and.(ifound1.eq.0.or.ifound2.eq.0.or.ifound3.eq.0))
            it=it+1
            write(6,'(1x,i3,1x,a,3(1x,a),9(1x,e12.5))') it,type(it),threebgrp(1,id),threebgrp(2,id),threebgrp(3,id) &
   ,(threebangs(j,id),j=1,nthreebangs(id))
            if(type(it).eq.threebgrp(1,id)) ifound1=it
            if(type(it).eq.threebgrp(2,id)) ifound2=it
            if(type(it).eq.threebgrp(3,id)) ifound3=it
         end do
         if(ifound1.eq.0.or.ifound2.eq.0.or.ifound3.eq.0) then
            write(6,100) id
100   format(/'find_triplets> atoms not matched to those in .ato file for triplet',i3)
            return
         endif    
         ltyp1=ifound1
         ltyp2=ifound2
         ltyp3=ifound3

!Make a list of all the atoms of type ltyp1, ltype2, and ltyp3

         natom1list=0
         natom2list=0
         natom3list=0
         do ia2=1,natom
            if(itype(ia2).eq.ltyp1) natom1list=natom1list+1
            if(itype(ia2).eq.ltyp2) natom2list=natom2list+1
            if(itype(ia2).eq.ltyp3) natom3list=natom3list+1
         end do
         write(6,'(4(i2,1x),3(i5,1x))') id,ltyp1,ltyp2,ltyp3,natom1list,natom2list,natom3list
         call reallocate1d_i(atom1list,natom1list)
         write(6,'(2(i2,1x),2(i5,1x))') id,ltyp1,natom1list,size(atom1list)
         call reallocate1d_i(atom2list,natom2list)
         write(6,'(2(i2,1x),2(i5,1x))') id,ltyp1,natom2list,size(atom2list)
         call reallocate1d_i(atom3list,natom3list)
         write(6,'(2(i2,1x),2(i5,1x))') id,ltyp1,natom3list,size(atom3list)
!Repeat this saving the atom numbers                  
         iref1=0
         iref2=0
         iref3=0
         do ia2=1,natom
            if(itype(ia2).eq.ltyp1) then
               iref1=iref1+1
               atom1list(iref1)=ia2
            endif
            if(itype(ia2).eq.ltyp2) then
               iref2=iref2+1
               atom2list(iref2)=ia2
            endif
            if(itype(ia2).eq.ltyp3) then
               iref3=iref3+1
               atom3list(iref3)=ia2
            end if
         end do
         write(6,*) id,natom1list,natom2list,natom3list

         ! For each type 2 atom we need to determine which type 1 and type 3 atoms are within the specified distances,
         ! on the understanding that the type 1 and type 3 atoms cannot be the same atom as type 2, nor can
         ! the atom of type 3 be the same as the atom of type 1
         do iref2=1,natom2list
            ia2=atom2list(iref2)
            im2=mref(ia2)
!            write(6,*) 'Molecule: ',im2
            !Get the coordinates of this atom
            x2=xa(ia2)+xadd(im2)
            y2=ya(ia2)+yadd(im2)
            z2=za(ia2)+zadd(im2)
            iref1=0
            !Search for all possible pairs of type 1 and type 3 atoms that satisfy the distance constraints
            do while (iref1.lt.natom1list)
               !Search type 1 atoms
               foundia1=.false.
               do while(iref1.lt.natom1list.and..not.foundia1)
                  iref1=iref1+1
                  ia1=atom1list(iref1)
                  im1=mref(ia1)
                  if(ia1.ne.ia2.and.((threebsamemol(1,id).eq.1.and.im1.eq.im2).or.(threebsamemol(1,id).eq.0.and.im1.ne.im2))) then
                     !Get the coordinates of this atom
                     x1=xa(ia1)+xadd(im1)
                     y1=ya(ia1)+yadd(im1)
                     z1=za(ia1)+zadd(im1)
                     xdif=x1-x2
                     ydif=y1-y2
                     zdif=z1-z2
! apply periodic boundary conditions
                     call minimagec(xdif,ydif,zdif)
                     dsq=xdif*xdif+ydif*ydif+zdif*zdif
                     d12=sqrt(dsq)
!                     add12=threebrmax(1,id)*threebrwidth(id)
                     foundia1=d12.lt.(threebrmax(1,id))
                  end if
               end do
! No point in searching for type 3 atoms if no type 1 was found
               if(foundia1) then
!            write(6,*) 'Molecule: ',im2,im1,ia1
                  iref3=0
                  foundia3=.false.
                  do iref3=1,natom3list
                     ia3=atom3list(iref3)
                     im3=mref(ia3)
!Need to ensure unique combinations of ia3 and ia1 are selected
                     if(((threebsamemol(2,id).eq.1.and.im2.eq.im3).or.(threebsamemol(2,id).eq.0.and.im2.ne.im3)) &
      .and.((threebsamemol(3,id).eq.1.and.im1.eq.im3).or.(threebsamemol(3,id).eq.0.and.im1.ne.im3)) &
      .and.((ltyp3.eq.ltyp1.and.ia3.gt.ia1).or.(ltyp3.ne.ltyp1.and.ia3.ne.ia1))) then
                        !Get the coordinates of this atom
                        x3=xa(ia3)+xadd(im3)
                        y3=ya(ia3)+yadd(im3)
                        z3=za(ia3)+zadd(im3)
                        xdif=x3-x2
                        ydif=y3-y2
                        zdif=z3-z2
! apply periodic boundary conditions
                        call minimagec(xdif,ydif,zdif)
                        dsq=xdif*xdif+ydif*ydif+zdif*zdif
                        d32=sqrt(dsq)
!                        add32=threebrmax(2,id)*threebrwidth(id)
                        foundia3=d32.lt.(threebrmax(2,id))
!.and.d32.gt.(threebrmax(2,id)-add32)
                        if(foundia3) then
                           ntriplets=ntriplets+1
                           if(ntriplets.gt.ntripletsmax) then
                              ntripletsmax=ntriplets+10000
                              call reallocate3d_i(tripletlist,3,ntripletsmax,nthreeb)
                           end if
                           tripletlist(1,ntriplets,id)=ia1
                           tripletlist(2,ntriplets,id)=ia2
                           tripletlist(3,ntriplets,id)=ia3
                           ntripletsinmol(im2,id)=ntripletsinmol(im2,id)+1
                           if(ntripletsinmol(im2,id).eq.1) firsttripletinmol(im2,id)=ntriplets
!                           write(6,*) im2,ia1,ia2,ia3
                           foundia3=.false.
! Need to check that im1 and im3 are listed as triplet neighbours of im2, assuming they are not the same molecule
                           if(im1.ne.im2) then
                              it=0
                              found=.false.
                              do while(it.lt.nthreebneighbmol(im2).and..not.found)
                                 it=it+1
                                 found=im1.eq.threebneighbmol(it,im2)
                              end do
                              if(.not.found) then
                                 nthreebneighbmol(im2)=nthreebneighbmol(im2)+1
!                                 write(6,*) ' find_triplets> ',im1,im2,im3,nthreebneighbmol(im2),nthreebneighbmolmax
                                 if(nthreebneighbmol(im2).gt.nthreebneighbmolmax) then
                                    nthreebneighbmolmax=nthreebneighbmolmax+4
                                    call reallocate2d_i(threebneighbmol,nthreebneighbmolmax,nmol)
                                 end if
!                                 write(6,*) ' find_triplets> ',im1,im2,im3,nthreebneighbmol(im2),nthreebneighbmolmax
                                 threebneighbmol(nthreebneighbmol(im2),im2)=im1
                              end if
                           end if
                           if(im3.ne.im2) then
                              it=0
                              found=.false.
                              do while(it.lt.nthreebneighbmol(im2).and..not.found)
                                 it=it+1
                                 found=im3.eq.threebneighbmol(it,im2)
                              end do
                              if(.not.found) then
                                 nthreebneighbmol(im2)=nthreebneighbmol(im2)+1
!                                 write(6,*) ' find_triplets> ',im1,im2,im3,nthreebneighbmol(im2),nthreebneighbmolmax
                                 if(nthreebneighbmol(im2).gt.nthreebneighbmolmax) then
                                    nthreebneighbmolmax=nthreebneighbmolmax+4
                                    call reallocate2d_i(threebneighbmol,nthreebneighbmolmax,nmol)
                                 end if
!                                 write(6,*) ' find_triplets> ',im1,im2,im3,nthreebneighbmol(im2),nthreebneighbmolmax
                                 threebneighbmol(nthreebneighbmol(im2),im2)=im3
                              end if
                           end if
                        end if
                     end if
                  end do
               end if
            end do
         end do
!         do im2=1,min(20,nmol)
!            write(6,'(a,i3,1x,i7,1x,i9,1x,i9)') 'find_triplets> No. of triplets in mol: ' &
!               ,id,im2,ntripletsinmol(im2,id),firsttripletinmol(im2,id)
!         end do
      end do
      call get_triplet_pot_cons()
!      do im=1,nmol
!         write(6,*) im,nthreebneighbmol(im),(threebneighbmol(it,im),it=1,nthreebneighbmol(im))
!      end do     

      end subroutine find_triplets

      subroutine get_triplet_pot_cons()

      use reallocation_routines

      real sumangle,sumlorentz,pi,angle,sineangle,sumtemp
      integer nangle

      pi=4.0*atan(1.0)
      if(nthreeb.le.0) return
      call reallocate1d_r(consthreeb,nthreeb)
      do id=1,nthreeb
         sumangle=0.0
         sumlorentz=0.0
         nangle=1001
         do ia=1,nangle
            angle=real(ia-1)*pi/real(nangle-1)
            sineangle=sin(angle)
            sumangle=sumangle+sineangle
            do ib=1,nthreebangs(id)
               arg=(angle-threebangs(ib,id))/threebawidth(id)
               arg=0.5*arg*arg
               if(arg.lt.30.) then
                  sumlorentz=sumlorentz+threebamps(ib,id)*exp(-arg)
               end if
!               sumlorentz=sumlorentz+threebamps(ib,id)*sineangle/(1.0+arg*arg)
            end do
         end do
         consthreeb(id)=sumlorentz/sumangle
      end do
      write(6,*) 'get_triplet_pot_cons> Potential constants: ',(consthreeb(id),id=1,nthreeb)
      return
      end subroutine get_triplet_pot_cons


      function en_mole_triplets(imol,xadd,yadd,zadd,phix,phiy,phiz,xa,ya,za,mref,tethered)
! Calculates the three-body energy associated with molecule imol
      use reallocation_routines

      real :: en_mole_triplets,ensum,x1,x2,x3,y1,y2,y3,z1,z2,z3,d12sq,d32sq,d13sq,d12,d32,d12star,d32star
      real :: cost,cosa,pot,sigma,width,test,ratio,entrimin,entrisum,truncatefactor,arg,argtest
      real, dimension(:)  :: xadd,yadd,zadd,xa,ya,za,phix,phiy,phiz
      integer :: imol,im1,im2,im3,ia1,ia2,ia3,itf,itl,i,imin,iia2ia3pairs,nia2ia3pairs,icount,ib,ibmin
      integer, dimension(:) :: mref
      logical  :: found
      logical, dimension(:) :: tethered

! Range of triplets to be calculated

      en_mole_triplets=0.0
      if(nthreeb.le.0) return
      ensum=0.0
      do id=1,nthreeb
         en_mole_triplet_dist(id)=0.0
         if(abs(threebamp(id)).gt.0.0) then
! OK now step through each set of pairs and find the angle value which is closest to the specified angle
            entrisum=0.0
            itf=firsttripletinmol(imol,id)
            itl=itf+ntripletsinmol(imol,id)-1
            do it=itf,itl
! Get coordinates of atoms in triplet
               ia1=tripletlist(1,it,id)
               ia2=tripletlist(2,it,id)
               ia3=tripletlist(3,it,id)
               im1=mref(ia1)
               im2=mref(ia2)
               im3=mref(ia3)
               x1=xa(ia1)+xadd(im1)
               x2=xa(ia2)+xadd(im2)
               x3=xa(ia3)+xadd(im3)
               y1=ya(ia1)+yadd(im1)
               y2=ya(ia2)+yadd(im2)
               y3=ya(ia3)+yadd(im3)
               z1=za(ia1)+zadd(im1)
               z2=za(ia2)+zadd(im2)
               z3=za(ia3)+zadd(im3)
               distmax12=threebrmax(1,id)*(1.0+threebrwidth(id))
               distmax32=threebrmax(2,id)*(1.0+threebrwidth(id))
! Get the distances 1 - 2, 2 - 3
               xdif=x1-x2
               ydif=y1-y2
               zdif=z1-z2
               call minimagec(xdif,ydif,zdif)
               d12sq=xdif*xdif+ydif*ydif+zdif*zdif
               d12=sqrt(d12sq)
               xdif=x3-x2
               ydif=y3-y2
               zdif=z3-z2
               call minimagec(xdif,ydif,zdif)
               d32sq=xdif*xdif+ydif*ydif+zdif*zdif
               d32=sqrt(d32sq)
! Get distance 1 - 3
               xdif=x1-x3
               ydif=y1-y3
               zdif=z1-z3
               call minimagec(xdif,ydif,zdif)
               d13sq=xdif*xdif+ydif*ydif+zdif*zdif
               cost=(d12sq+d32sq-d13sq)/(2.0*d12*d32)
               cost=max(-1.0,min(1.0,cost))
!               cosa=acos(cost)
               entrimin=0.0
               ib=1
               ibmin=ib
               arg=cost-threebangs(ib,id)
               do while (ib.lt.nthreebangs(id))
                  ib=ib+1
                  argtest=cost-threebangs(ib,id)
                  if(abs(argtest).lt.abs(arg)) then
                     arg=argtest
                     ibmin=ib
                  end if
               end do
               arg=arg/threebawidth(id)
               term=(0.5*arg*arg)
               entrimin=entrimin+threebamps(ibmin,id)*term
               truncatefactor=1.0
               entrisum=entrisum+truncatefactor*entrimin
               en_mole_triplet_dist(id)=en_mole_triplet_dist(id)+truncatefactor*entrimin
!               if(im2.lt.21) then 
!                  write(6,*) id,im1,im2,im3,ia1,ia2,ia3,d12,d32,cosa,threebangs(imin,id),entrisum
!                  write(6,*) im1,im2,xadd(im1),yadd(im1),zadd(im1),xadd(im2),yadd(im2),zadd(im2),entrimin,ensum
!                  write(6,*) xa(ia1),ya(ia1),za(ia1),xa(ia2),ya(ia2),za(ia2),xa(ia3),ya(ia3),za(ia3)
!                  write(6,*) id,ia1,ia2,ia3,x2,x3,y2,y3,z2,z3,entrisum
!                  write(6,*),d12,d32,cosa,entrimin
!               end if
            end do
            if(ntripletsinmol(imol,id).gt.0) then
               en_mole_triplet_dist(id)=threebamp(id)*en_mole_triplet_dist(id)
               ensum=ensum+threebamp(id)*entrisum
            end if
         end if
      end do
!      if(imol.lt.21) write(6,*) imol,ensum,(id,en_mole_triplet_dist(id),id=1,nthreeb)
      en_mole_triplets=ensum

      end function en_mole_triplets

      function en_mole_triplets_old(imol,xadd,yadd,zadd,phix,phiy,phiz,xa,ya,za,mref,tethered)
! Calculates the three-body energy associated with molecule imol
      use reallocation_routines

      real :: en_mole_triplets,ensum,x1,x2,x3,y1,y2,y3,z1,z2,z3,d12sq,d32sq,d13sq,d12,d32,d12star,d32star
      real :: cost,cosa,pot,sigma,width,test,ratio
      real, dimension(:)  :: xadd,yadd,zadd,xa,ya,za,phix,phiy,phiz
      integer :: imol,im1,im2,im3,ia1,ia2,ia3,itf,itl,i,imin,iia2ia3pairs,nia2ia3pairs,icount
      integer, dimension(:) :: mref
      integer, dimension(:,:), allocatable :: ia2ia3pairs
      logical  :: found
      logical, dimension(:) :: tethered

! Range of triplets to be calculated

      en_mole_triplets=0.0
      if(nthreeb.le.0) return
      ensum=0.0
      do id=1,nthreeb
         en_mole_triplet_dist(id)=0.0
         if(ntripletsinmol(imol,id).gt.0.and.abs(threebamp(id)).gt.0.0) then
            call reallocate2d_i(ia2ia3pairs,2,ntripletsinmol(imol,id))
            itf=firsttripletinmol(imol,id)
            itl=itf+ntripletsinmol(imol,id)-1
            !Step through listed triplets.
! It will noted that if nato(imol)>1, there may be more than one triplet associated with the atom pair ia2,ia3.
! We will need to identify the one of these pairs that is closest to the specified bond angle.
            nia2ia3pairs=0
            do it=itf,itl
! Get coordinates of atoms in triplet
               ia1=tripletlist(1,it,id)
               ia2=tripletlist(2,it,id)
               ia3=tripletlist(3,it,id)
               found=.false.
               iia2ia3pairs=0
               do while(iia2ia3pairs.lt.nia2ia3pairs.and..not.found)
                  iia2ia3pairs=iia2ia3pairs+1
                  found=ia2.eq.ia2ia3pairs(1,iia2ia3pairs).and.ia3.eq.ia2ia3pairs(2,iia2ia3pairs)
               end do
               if(.not.found) then
                  nia2ia3pairs=nia2ia3pairs+1
                  ia2ia3pairs(1,nia2ia3pairs)=ia2
                  ia2ia3pairs(2,nia2ia3pairs)=ia3
               end if
            end do
! OK now step through each set of pairs and find the angle value which is closest to the specified angle
            do ipairs=1,nia2ia3pairs
!               write(6,*) ipairs,mref(ia2ia3pairs(1,ipairs)),mref(ia2ia3pairs(2,ipairs)),ia2ia3pairs(2,ipairs),ia2ia3pairs(2,ipairs)
               icount=0
               imin=0
               enmin=0
               itmin=0
               d12min=0
               d32min=0
               entrisum=0
               do it=itf,itl
! Get coordinates of atoms in triplet
                  ia1=tripletlist(1,it,id)
                  ia2=tripletlist(2,it,id)
                  ia3=tripletlist(3,it,id)
                  if(ia2.eq.ia2ia3pairs(1,ipairs).and.ia3.eq.ia2ia3pairs(2,ipairs)) then
                     im1=mref(ia1)
                     im2=mref(ia2)
                     im3=mref(ia3)
                     x1=xa(ia1)+xadd(im1)
                     x2=xa(ia2)+xadd(im2)
                     x3=xa(ia3)+xadd(im3)
                     y1=ya(ia1)+yadd(im1)
                     y2=ya(ia2)+yadd(im2)
                     y3=ya(ia3)+yadd(im3)
                     z1=za(ia1)+zadd(im1)
                     z2=za(ia2)+zadd(im2)
                     z3=za(ia3)+zadd(im3)
                     distmax12=threebrmax(1,id)*(1.0+threebrwidth(id))
                     distmax32=threebrmax(2,id)*(1.0+threebrwidth(id))
! Get the distances 1 - 2, 2 - 3
                     xdif=x1-x2
                     ydif=y1-y2
                     zdif=z1-z2
                     call minimagec(xdif,ydif,zdif)
                     d12sq=xdif*xdif+ydif*ydif+zdif*zdif
                     d12=sqrt(d12sq)
                     xdif=x3-x2
                     ydif=y3-y2
                     zdif=z3-z2
                     call minimagec(xdif,ydif,zdif)
                     d32sq=xdif*xdif+ydif*ydif+zdif*zdif
                     d32=sqrt(d32sq)
! Get distance 1 - 3
                     xdif=x1-x3
                     ydif=y1-y3
                     zdif=z1-z3
                     call minimagec(xdif,ydif,zdif)
                     d13sq=xdif*xdif+ydif*ydif+zdif*zdif
                     cost=(d12sq+d32sq-d13sq)/(2.0*d12*d32)
                     cosa=acos(cost)
                     entrimin=0
! Get the closest value from the list of requested angles. Also save the value of ia1 that corresponds to this minimum.
                     do ib=1,nthreebangs(id)
                        icount=icount+1
                        test=abs(cosa-threebangs(ib,id))
                        if(ib.eq.1) then
                           entrimin=test
                        else if (test.lt.entrimin) then
                           entrimin=test
                        end if
!                        if(icount.eq.1) then
!                           enmin=test
!                           itmin=it
!                           d12min=d12
!                           d32min=d32
!                        else if(enmin.gt.test) then
!                           enmin=test
!                           itmin=it
!                           d12min=d12
!                           d32min=d32
!                        end if
                     end do
                     truncatefactor=1.0
!                     if(distmax12.gt.0.0.and.im1.ne.im2) then
! Get relative distances
!                        d12star=d12/distmax12
!                        truncatefactor=truncatefactor*truncate_threeb(d12star)
!                     end if
!                     if(distmax32.gt.0.0.and.im2.ne.im3) then
!                        d32star=d32/distmax32
!                        truncatefactor=truncatefactor*truncate_threeb(d32star)
!                     end if
! Use a harmonic three-body potential to emphasize the angular part.

                     entrisum=entrisum+entrimin*entrimin*truncatefactor
!                     write(6,*) 'en_mole_triplets> ',iia2ia3pairs,icount,imol,im1,im2,im3,ia1,ia2,ia3,cosa,enmin
                  end if
               end do
!               ia1=tripletlist(1,itmin,id)
!               ia2=tripletlist(2,itmin,id)
!               ia3=tripletlist(3,itmin,id)
!               im1=mref(ia1)
!               im2=mref(ia2)
!               im3=mref(ia3)
!               write(6,*) 'en_mole_triplets> ',imol,im1,im2,im3,ia1,ia2,ia3,enmin
! If im3 and im2 are the same, then ia2 and ia3 are on the same molecule. If, in addition, the requested angle is zero
! it means we expect ia3 to lie between ia1 and ia2. In that case, instead of calculating the angle, we attempt to tether ia3 to
! its respective position directly between ia1 and ia2
!               if(im2.eq.im3.and.nthreebangs(id).eq.1.and.threebangs(1,id).eq.0.0) then
!                  ratio=d32min/d12min
!                  ratio=0.5
!                  if(tethered(im1)) then
!                     x1=phix(im1)+xa(ia1)
!                     y1=phiy(im1)+ya(ia1)
!                     z1=phiz(im1)+za(ia1)
!                  else
!                     x1=xa(ia1)+xadd(im1)
!                     y1=ya(ia1)+yadd(im1)
!                     z1=za(ia1)+zadd(im1)
!                  endif
!                  if(tethered(im2)) then
!                     x2=phix(im2)+xa(ia2)
!                     y2=phiy(im2)+ya(ia2)
!                     z2=phiz(im2)+za(ia2)
!                  else
!                     x2=xa(ia2)+xadd(im2)
!                     y2=ya(ia2)+yadd(im2)
!                     z2=za(ia2)+zadd(im2)
!                  endif
!                  x3=xa(ia3)+xadd(im3)
!                  y3=ya(ia3)+yadd(im3)
!                  z3=za(ia3)+zadd(im3)
!                  xdif=x3-x2-ratio*(x1-x2)
!                  ydif=y3-y2-ratio*(y1-y2)
!                  zdif=z3-z2-ratio*(z1-z2)
!                  call minimagec(xdif,ydif,zdif)
!                  ensum=ensum+threebamp(id)*(xdif*xdif+ydif*ydif+zdif*zdif)
!               else
               ensum=ensum+threebamp(id)*entrisum
               en_mole_triplet_dist(id)=en_mole_triplet_dist(id)+threebamp(id)*entrisum
            end do
         end if
      end do
      en_mole_triplets=ensum

      end function en_mole_triplets_old

      function pr_mole_triplets(imol,xadd,yadd,zadd,xa,ya,za,mref)
! Calculates the three-body pressure associated with molecule imol
      use reallocation_routines

      real(kind=8) :: pr_mole_triplets,prsum
      real :: x1,x2,x3,y1,y2,y3,z1,z2,z3,d12sq,d32sq,d13sq,d12,d13,d32,d12star,d32star
      real :: cost,cosa,pot,sigma,width
      real, dimension(:)  :: xadd,yadd,zadd,xa,ya,za
      integer :: imol,im1,im2,im3,ia1,ia2,ia3,itf,itl,i,imin
      integer, dimension(:) :: mref
      real :: xdif,ydif,zdif,xcdif,ycdif,zcdif,prfac

! Range of triplets to be calculated

      pr_mole_triplets=0.0
      if(nthreeb.le.0) return
      prsum=0.0
      do id=1,nthreeb
         if(ntripletsinmol(imol,id).gt.0.and.abs(threebamp(id)).gt.0.0) then
            itf=firsttripletinmol(imol,id)
            itl=itf+ntripletsinmol(imol,id)-1
            !Step through listed triplets
            do it=itf,itl
! Get coordinates of atoms in triplet
               ia1=tripletlist(1,it,id)
               ia2=tripletlist(2,it,id)
               ia3=tripletlist(3,it,id)
               im1=mref(ia1)
               im2=mref(ia2)
               im3=mref(ia3)
               x1=xa(ia1)+xadd(im1)
               x2=xa(ia2)+xadd(im2)
               x3=xa(ia3)+xadd(im3)
               y1=ya(ia1)+yadd(im1)
               y2=ya(ia2)+yadd(im2)
               y3=ya(ia3)+yadd(im3)
               z1=za(ia1)+zadd(im1)
               z2=za(ia2)+zadd(im2)
               z3=za(ia3)+zadd(im3)
! Get the distances 1 - 2, 2 - 3
               xdif=x1-x2
               ydif=y1-y2
               zdif=z1-z2
               call minimagec(xdif,ydif,zdif)
               d12sq=xdif*xdif+ydif*ydif+zdif*zdif
               d12=sqrt(d12sq)
               xdif=x3-x2
               ydif=y3-y2
               zdif=z3-z2
               call minimagec(xdif,ydif,zdif)
               d32sq=xdif*xdif+ydif*ydif+zdif*zdif
               d32=sqrt(d32sq)
               if(d12.lt.threebrmax(1,id).and.d32.lt.threebrmax(2,id)) then
! Get relative distances
                  d12star=d12/threebrmax(1,id)
                  d32star=d32/threebrmax(2,id)
! Get distance 1 - 3
                  xdif=x1-x3
                  ydif=y1-y3
                  zdif=z1-z3
                  call minimagec(xdif,ydif,zdif)
                  d13sq=xdif*xdif+ydif*ydif+zdif*zdif
!                  d13=sqrt(d13sq)
                  cost=(d12sq+d32sq-d13sq)/(2.0*d12*d32)
                  xcdif=xadd(im1)-xadd(im3)
                  ycdif=yadd(im1)-yadd(im3)
                  zcdif=zadd(im1)-zadd(im3)
                  call minimagec(xcdif,ycdif,zcdif)
                  if(d12.gt.0.0.and.d32.gt.0.0) then
                     prfac=-(xcdif*xdif+ycdif*ydif+zcdif*zdif)/(d12*d32)!Don't forget d(cost)/d13
                  else
                     prfac=0.0
                  endif
                  entrimin=0.0
                  ib=1
                  ibmin=ib
                  arg=cost-threebangs(ib,id)
                  do while (ib.lt.nthreebangs(id))
                     ib=ib+1
                     argtest=cost-threebangs(ib,id)
                     if(abs(argtest).lt.abs(arg)) then
                        arg=argtest
                        ibmin=ib
                     end if
                  end do
                  pot=prfac*threebamps(ibmin,id)*arg/(threebawidth(id)*threebawidth(id))
                  prsum=prsum+threebamp(id)*pot
               end if
            end do
         end if
      end do
      pr_mole_triplets=prsum

      end function pr_mole_triplets

      function truncate_threeb(x)

      real :: truncate_threeb,x,arg

      truncate_threeb=0.0
      if(x.ge.1.0) return
      arg=1.0/(x-1.0)+1.0
      if(arg.gt.-30.0) then
         truncate_threeb=exp(arg)
      else
         truncate_threeb=0.0
      end if

      end function truncate_threeb

      function mod_morse(x,sigma,width)

      real :: mod_morse,x,sigma,width,arg

      mod_morse=0.0
      if(width.gt.0.0) then
         arg=(sigma-x)/width
         if(arg.gt.30.0) then
            mod_morse=(arg-1.0)*exp(30.0)
         else if(arg.lt.-30.0) then
            mod_morse=0.0
         else
            mod_morse=(arg-1.0)*exp(arg)
         end if
      end if

      end function mod_morse

      function mod_morse1(x,sigma,width)

      real :: mod_morse1,x,sigma,width,arg

      mod_morse1=0.0
      if(width.gt.0.0) then
         arg=(sigma-x)/width
         arg=0.5*arg*arg
!!         arg=-abs(arg)
         if(arg.gt.30.0) then
            mod_morse1=0.0
         else
!!            mod_morse1=(arg-1.0)*exp(arg)
            mod_morse1=exp(-arg)
         end if
!         mod_morse1=1.0/(1.0+arg*arg)
      end if
!      write(6,*) x,sigma,width,mod_morse1

      end function mod_morse1

      function grad_mod_morse1(x,sigma,width)

      real :: grad_mod_morse1,x,sigma,width,arg

      grad_mod_morse1=0.0
      if(width.gt.0.0) then
         arg=(sigma-x)/width
!         arg=-abs(arg)
         if(arg.lt.-30.0) then
            grad_mod_morse1=0.0
         else
!            if(x.gt.sigma) then
               grad_mod_morse1=-arg*exp(arg) !Gradient is +ve if x > sigma
!            else
!               grad_mod_morse1=arg*exp(arg) !Gradient is -ve if x < sigma
!            end if
            grad_mod_morse1=grad_mod_morse1/width
         end if
      end if

      end function grad_mod_morse1

      subroutine save_threeb_energies(im)

      integer im

      if(nthreeb.gt.0) then
         ergold_triplet(im)=ergnew_triplet(im)
! IMPORTANT! Include the neighbours of im which may have been affected by the movement of im!
         if(nthreebneighbmol(im).gt.0) then
            do it=1,nthreebneighbmol(im)
               imref=threebneighbmol(it,im)
               ergold_triplet(imref)=ergnew_triplet(imref)
            end do
         end if
      end if
      return
      end subroutine save_threeb_energies


end module threebody_mod


