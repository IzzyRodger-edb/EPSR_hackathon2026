      module menergpairs

      implicit none
      
      real, dimension(:), allocatable :: ergnew
      real, dimension(:,:), allocatable :: ergpair
     integer :: menerg, menergadd, mergnew

      contains

      subroutine initialise_ergpair(nmol)

      integer, intent(in) :: nmol

      write(6,*) menerg,menergadd
      if(menerg.eq.0) then
         menerg=350*natom !Based on ~700 atom neighbours per atom - need half of this for pairs.
      endif
      write(6,*) menerg,menergadd
      if(allocated(ergpair)) deallocate(ergpair)
      allocate(ergpair(menerg,nmol))

      return
      end subroutine initialise_ergpair

      subroutine initialise_menerg(nmol)

      integer, intent(in) :: nmol
      integer :: ilenfname,ilenlong,ierr
      character(len=256) :: fname,fnamelong

! Create the filename to be opened to read menerg

      fname='menerg.dat'
      menerg=0
      ilenfname=len_trim(fname)
      call get_full_filename(fname,ilenfname,fnamelong,ilenlong)
        write(6,'(a,a)') 'initialise_menerg> filename: ',fnamelong(1:ilenlong)
      open(20,file=fnamelong,status='old',iostat=ierr)
      if (ierr.eq.0) then
         read(20,*,iostat=ierr) menerg
         if(ierr.ne.0.or.menerg.lt.500*natom) then
            menerg=500
         endif
         read(20,*,iostat=ierr) menergadd
         if(ierr.ne.0.or.menergadd.lt.100000) then
            menergadd=100
         endif
        write(6,*) 'initialise_menerg> ',menerg,menergadd

      else

         menerg=500
         menergadd=500

      endif
      close(20)
      if(allocated(ergpair)) deallocate(ergpair)
      allocate(ergpair(menerg,mol))
      if(allocated(ergnew)) deallocate(ergnew)
      allocate(ergnew(menerg))

      return
      end subroutine initialise_menerg

      subroutine save_menerg()

      integer :: ilenfname,ilenlong,ierr
      character(len=256) :: fname,fnamelong

! Create the filename to be opened to save menerg

      fname='menerg.dat'
      ilenfname=len_trim(fname)
      call get_full_filename(fname,ilenfname,fnamelong,ilenlong)
      open(20,file=fnamelong,status='unknown',iostat=ierr)
      if (ierr.eq.0) then
         write(20,*,iostat=ierr) menerg
         write(20,*,iostat=ierr) menergadd
      endif
      close(20)

      return
      end subroutine save_menerg

      subroutine reallocate_ergpair(nmol,nergref)

      integer, intent(in) :: nmol,nergref
      real, dimension(:,:), allocatable :: ergpairsave
      real, dimension(:), allocatable :: ergnewsave
      integer :: i,im

      do while (nergref.gt.menerg)

c First save the current values

         if(menerg.gt.0) then
            allocate(ergnewsave(menerg))
            allocate(ergpairsave(menerg,nmol)
            do i=1,menerg
                  ergnewsave(i)=ergnew(i)
            end do
            do im=1,nmol
               do i=1,menerg
                  ergpairsave(i,im)=ergpair(i,im)
               end do
            end do
         end if
c Deallocate the existing values
         if(allocated(ergnew)) deallocate(ergnew)
         if(allocated(nergnew)) deallocate(nergnew)
         if(allocated(ergpair)) deallocate(ergpair)
c Allocate the new array
         allocate(ergnew(menerg+menergadd))
         allocate(ergpair(menerg+menergadd,nmol))
c Save it
         if(menerg.gt.0) then
            do i=1,menerg
                  ergnew(i)=ergnewsave(i)
            end do
            do im=1,nmol
               do i=1,menerg
                  ergpair(i,im)=ergpairsave(i,im)
               end do
            end do
         end if
c Finally deallocate the save array and increment maximum number of elements in ergpair.
         if(allocated(ergnewsave)) deallocate(ergnewsave)
         if(allocated(ergpairsave)) deallocate(ergpairsave)
         write(6,*) 'ergpair> ',menerg,menerg+menergadd
         menerg=menerg+menergadd

      end do
      return
      end subroutine reallocate_ergpair
         
      subroutine save_energies()
      
      include 'epsr_dimension.inc'
      include 'potential_par.inc'
      integer*4 jref
c
c saves the new energy in array ergpair
c

c     ediffsum=0
      do jref=1,nergnewref
         if(nergnew(jref).gt.nergref) then
            write(6,*) 'save_energies> ERROR? ',jref,nergnew(jref),nergref,ergnew(jref)
         endif
c        ediff=ergnew(jref)-ergpair(nergnew(jref))
c        ediffsum=ediffsum+ediff
         if(nergnew(jref).gt.menerg) call reallocate_ergpair(nergnew(jref))
         ergpair(nergnew(jref))=ergnew(jref)
      end do

      return
      end

      subroutine get_ergpair(im)
      
      integer, intent(in) :: im
      
      !Gets the pair energies for specified molecule         
      end module menergpairs


