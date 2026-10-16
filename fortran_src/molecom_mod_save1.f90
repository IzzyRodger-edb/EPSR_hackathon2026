		module molecom_mod

	   use reallocation_routines

      integer                                            :: ngrrot, ngrrotmax! Number of rotational groups and current maximum number
      character(len=3), dimension(:), allocatable        :: rotgrp !Rotational group label (currently only ROT or DIH)
      integer, dimension(:), allocatable                 :: nrot !Number of rotatable atoms in each rotational group
      integer, dimension(:,:), allocatable               :: labaxi !Pair of atoms about which rotations are made
      integer, dimension(:), allocatable                 :: rotf,rotl !First and last in list of rotatable atoms for each group
      integer, dimension(:), allocatable                 :: mrot !Molecule number for each rotational group

      integer                                            :: nrottot,nrottotmax !Total number of rotatable atoms in all groups     
      integer, dimension(:), allocatable                 :: labrot !List of rotatable atoms for all groups

      contains

      subroutine testsize_rotgrps(newsizein)

	      use reallocation_routines
         
         integer, optional          :: newsizein
         integer                    :: newsize
         if(present(newsizein)) then
            newsize=newsizein
         else
            newsize=ngrrot
         end if
         do while(newsize.gt.ngrrotmax)
            ngrrotmax=newsize
            call reallocate1d_c(rotgrp,len(rotgrp),ngrrotmax)
            call reallocate2d_i(labaxi,ngrrotmax,2)
            call reallocate1d_i(nrot,ngrrotmax)
            call reallocate1d_i(rotf,ngrrotmax)
            call reallocate1d_i(rotl,ngrrotmax)
            call reallocate1d_i(mrot,ngrrotmax)
         end do

      end subroutine testsize_rotgrps

      subroutine testsize_labrot(newsizein)

	      use reallocation_routines
         
         integer, optional          :: newsizein
         integer                    :: newsize

         if(present(newsizein)) then
            newsize=newsizein
         else
            newsize=nrottot
         end if
         do while(newsize.gt.nrottotmax)
            nrottotmax=newsize
            call reallocate1d_i(labrot,nrottotmax)
         end do

      end subroutine testsize_labrot

      end module molecom_mod


