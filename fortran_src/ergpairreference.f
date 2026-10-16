      module ergpairstore
      
      implicit none
      
      real, dimension(:), allocatable :: ergold,ergnew
      integer, dimension(:), allocatable :: nbref
      integer :: nnbref
      
      contains
      
      subroutine initialise_ergpairstore(nnbrefmax)
      
      integer :: nnbrefmax
      
      nnbref=0
      if(allocated(ergold)) deallocate(ergold,ergnew,nbref)
      allocate(ergold(nnbrefmax),ergnew(nnbrefmax),nbref(nnbrefmax))
      
      end subroutine initialise_ergpairstore

      subroutine set_ergnew()

      integer ibref

      do ibref=1,nnbref
         ergnew(ibref)=ergold(ibref)
      end do
      return

      end subroutine set_ergnew
      
      subroutine save_energies()
      
      integer ibref
      
      do ibref=1,nnbref
         ergold(ibref)=ergnew(ibref)
      end do
      return
      
      end subroutine save_energies
      
      end module ergpairstore              

      module ergpairreference

!This version uses a single dimension array to store pair energies

      implicit none

!Stores and manipulates the arrays to store and reference the set of pair energies
      real, dimension(:), allocatable :: ergpair
      integer :: nergpair,nergpairmax
      integer, dimension(:,:), allocatable :: molergpair,ergpairref
      integer, dimension(:), allocatable :: nergpairm
      integer :: nnbrefmax
      real, dimension(:), allocatable :: ergold,ergnew
      integer, dimension(:), allocatable :: nbref
      integer :: nnbref,nmolhere

      contains

      subroutine initialise_ergpair(nmol)

      use cellneighbourlist

      integer im,nmol
      
      nmolhere=nmol
      call get_nnbrefmax(nnbrefmax)
      !Allocate the arrays to store pair energies and reference integers
      if(allocated(ergpair)) deallocate(ergpair)
      if(allocated(molergpair)) deallocate(molergpair)
      if(allocated(ergpairref)) deallocate(ergpairref)
      if(allocated(nergpairm)) deallocate(nergpairm)
      nergpair=0
      nergpairmax=100*nmol
      allocate(ergpair(nergpairmax),molergpair(nnbrefmax,nmol),ergpairref(nnbrefmax,nmol),nergpairm(nmol))
      do im=1,nmol
         nergpairm(im)=0
      end do
      nnbref=0
      if(allocated(ergold)) deallocate(ergold,ergnew,nbref)
      allocate(ergold(nnbrefmax),ergnew(nnbrefmax),nbref(nnbrefmax))
      write(6,*) 'initialise_ergpair> nnbrefmax = ',nnbrefmax,size(ergnew)
      return

      end subroutine initialise_ergpair
      
      subroutine set_ergpair(im,jm,energy)
      
      real :: energy
      integer :: im,jm
      
      if(nergpairm(im).lt.nnbrefmax.and.nergpairm(jm).lt.nnbrefmax.and.energy.ne.0.0) then
!Store the energy
         nergpair=nergpair+1
         if(nergpair.gt.nergpairmax) then
            call reallocate_ergpair()
         end if
         ergpair(nergpair)=energy
!Update the energy reference arrays
         if(jm.eq.im) then
            nergpairm(im)=1
            molergpair(nergpairm(im),im)=im
            ergpairref(nergpairm(im),im)=nergpair
         else
            nergpairm(im)=nergpairm(im)+1
            nergpairm(jm)=nergpairm(jm)+1
            molergpair(nergpairm(im),im)=jm
            molergpair(nergpairm(jm),jm)=im
            ergpairref(nergpairm(im),im)=nergpair
            ergpairref(nergpairm(jm),jm)=nergpair
         end if
      end if
      return

      end subroutine set_ergpair
      
      subroutine reallocate_ergpair()

      real, dimension(:), allocatable :: ergpairsave
      integer :: i,nergpairmaxsave
      
      nergpairmaxsave=nergpairmax
!      write(6,*) 'reallocate_ergpair> ',nergpairmax
      allocate(ergpairsave(nergpairmaxsave))
      do i=1,nergpairmaxsave
         ergpairsave(i)=ergpair(i)
      end do
      deallocate(ergpair)
      nergpairmax=((nergpairmax+10)*12)/10
!      write(6,*) 'reallocate_ergpair> ',nergpairmax
      allocate(ergpair(nergpairmax))
      do i=1,nergpairmaxsave
         ergpair(i)=ergpairsave(i)
      end do
      deallocate(ergpairsave)
      return
      
      end subroutine reallocate_ergpair
            
      subroutine get_neighbsum(neighbsum)
      
      integer im,nmol,neighbsum
      
      neighbsum=0
      nmol=nmolhere
      do im=1,nmol
         neighbsum=neighbsum+nergpairm(im)
      end do
      if(nmol.gt.0) neighbsum=neighbsum/nmol
      return
      end subroutine get_neighbsum 
           
      subroutine make_neighbour_list_from_mols(im)
      
      integer :: im,ibref
      
      nnbref=nergpairm(im)
      if(.not.allocated(ergold)) then
         allocate(ergold(nnbref),ergnew(nnbref),nbref(nnbref))
      else if(size(ergold).lt.nnbref) then
         deallocate(ergold,ergnew,nbref)
         allocate(ergold(nnbref),ergnew(nnbref),nbref(nnbref))
      end if
      do ibref=1,nnbref
         nbref(ibref)=molergpair(ibref,im)
         ergold(ibref)=ergpair(ergpairref(ibref,im))
         ergnew(ibref)=ergold(ibref)
      end do
      return
      
      end subroutine make_neighbour_list_from_mols
      
      subroutine save_energies()
      
      integer ibref
      
      do ibref=1,nnbref
         ergold(ibref)=ergnew(ibref)
      end do
      return
      
      end subroutine save_energies

      subroutine save_erg_pairs(im)
      
      integer im,ibref
      
      do ibref=1,nnbref
         ergpair(ergpairref(ibref,im))=ergold(ibref)
      end do
      return
      
      end subroutine save_erg_pairs
           
      end module ergpairreference


