      module molecom_mod

      integer                                            :: ngrrot, ngrrotmax! Number of rotational groups and current maximum number
      character(len=3), dimension(:), allocatable        :: rotgrp !Rotational group label (currently only ROT or DIH)
      integer, dimension(:), allocatable                 :: nrot !Number of rotatable atoms in each rotational group
      integer, dimension(:,:), allocatable               :: labaxi !Pair of atoms about which rotations are made
      integer, dimension(:), allocatable                 :: rotf,rotl !First and last in list of rotatable atoms for each group
      integer, dimension(:), allocatable                 :: mrot !Molecule number for each rotational group

      integer                                            :: nrottot,nrottotmax !Total number of rotatable atoms in all groups     
      integer, dimension(:), allocatable                 :: labrot !List of rotatable atoms for all groups
      !Coarse-graining variables
      logical                                            :: coarsegraining !True if coarse graining is in place
      real                                               :: qscaleinputfactor!Scale factor on input Q scale - needed for coarse graining
      real                                               :: selfbroadpwr !Broadening power for self terms - attempts to mimic molecular structure
      real, dimension(:,:,:), allocatable                :: sphereformfactor !Used if coarse graining is on.
      real, dimension(:,:), allocatable                  :: intraformfactor,intraformfactorself !Used if coarse graining is on.
      real, dimension(:), allocatable                    :: scalefacmol,scalefacmoluse !Scale factor assigned to size of each molecule
      real, dimension(:), allocatable                    :: rfractmol !Relative broadening for each molecule type.
      real, dimension(:), allocatable                    :: natomsphere,rhosphere,radminsphere,radmaxsphere !Atomic number density, and minimum and maximum radii for coarse graining 
      real, dimension(:), allocatable                    :: scalefactypemin,scalefactypemax !Maximum and minimum scale factor for each molecule type
      real                                               :: scalefacratio=0.2 !Used to determine the number of distinct spherical form factors for a given molecule type
      integer, dimension(:), allocatable                 :: nsffbintype !Number of spherical form factor bins for a given atom type
      integer, dimension(:), allocatable                 :: assignatomscalefac !Assigned spherical form factor bin for each atom
      integer                                            :: maxnscalefactype !Maximum number of distinct spherical form factors
      integer, dimension(:,:), allocatable               :: ntermssff !No. of terms for each sphere form factor bin and each atom type
      real, dimension(:,:), allocatable                  :: fracntermssff !Fractional no. of terms for each sphere form factor bin and each atom type
      real, dimension(:,:,:,:), allocatable              :: wavinter! Used for distinguishing spherical form factors of different scale factors when coarse graining
      real, dimension(:,:), allocatable                  :: prodssff! Temporary form factor product used by gensq
      integer, dimension(:), allocatable                 :: ntermsintra,ntermsintraself !No. of terms for intra and intra self correlations
      real, dimension(:), allocatable                    :: pseudofacmol !No. of real atoms to be used for each pseudo atom. Normally should be 1.0, but increase for coarse graining.
      real                                               :: natomcg,natomrcg,rhocg!Counts total number of coarse grain atoms

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
         !We only reallocate labrot if the new size is bigger than the old size
         do while(newsize.gt.nrottotmax)
            nrottotmax=newsize
            call reallocate1d_i(labrot,nrottotmax)
         end do

      end subroutine testsize_labrot

      end module molecom_mod


