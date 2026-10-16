   module moltemplates_mod

	   use reallocation_routines

      integer                                         :: nmoltemplates
      integer                                         :: natotmax,nbondtmax,ngrrottmax,nrottmax,ndihtmax

      character(len=3), dimension(:), allocatable     :: typet!(near)
      character(len=256), dimension(:), allocatable   :: moltemplatefile! Filename where the template is stored, in case it is needed to upgrade current values (fmole).
      character(len=3), dimension(:,:), allocatable   :: rotgrpt!(near)

      integer, dimension(:), allocatable              :: natot,natomtypet!No. of atoms and no. of atom types in each template
      integer, dimension(:,:), allocatable            :: nbondt !No. of bonds for each atom in template
      integer, dimension(:,:), allocatable            :: ngrott! No. of rotational groups in each template
      integer, dimension(:,:,:), allocatable          :: iabondt,labrott!(near,near)
      integer, dimension(:,:,:), allocatable          :: labaxit!(near,2)

      real, dimension(:,:,:), allocatable             :: tbondt,dihedralangt!(near,near)

      integer, dimension(:,:), allocatable            :: ndihedralangt!(near,mmix)!(ma,mmix)
      integer, dimension(:,:,:,:), allocatable        :: labrotdiht!(2,near,near,mmix) 

   contains

      subroutine initialise_moltemplates()

         nmoltemplates=0
         natotmax=0
         nbondtmax=0
         ngrrottmax=0
         nrottmax=0
         ndihtmax=0

      end subroutine initialise_moltemplates

      subroutine reallocate_all(newsize)
         !Call this whenever the number of molecule templates increases
         integer :: newsize,nsize
         if(newsize.gt.nmoltemplates) nmoltemplates=newsize
         if(nmoltemplates.gt.0) then
            call reallocate1d_c(typet,len(typet),nmoltemplates)
            call reallocate1d_c(moltemplatefile,len(moltemplatefile),nmoltemplates)
            call reallocate1d_i(natot,nmoltemplates)
            call reallocate1d_i(natomtypet,nmoltemplates)
            call reallocate2d_c(alabelt,len(alabelt),natotmax,nmoltemplates)
            call reallocate2d_i(nbondt,natotmax,nmoltemplates)
            call reallocate3d_i(iabondt,nbondtmax,natotmax,nmoltemplates)
            call reallocate3d_r(tbondt,nbondtmax,natotmax,nmoltemplates)
            call reallocate2d_c(rotgrpt,len(rotgrpt),ngrrottmax,nmoltemplates)
            call reallocate3d_i(labaxit,2,ngrrotmax,nmoltemplates)
            call reallocate2d_i(nrott,ngrrotmax,nmoltemplates)
            call reallocate2d_i(ndihedralangt,ngrrotmax,nmoltemplates)
            call reallocate3d_i(labrott,nrottmax,ngrrotmax,nmoltemplates)
            call reallocate3d_r(dihedralangt,ndihtmax,ngrrotmax,nmoltemplates)
            call reallocate4d_i(labrotdiht,2,ndihtmax,ngrrotmax,nmoltemplates)
         end if
      end subroutine reallocate_all

      subroutine reallocate_atot(newsize)
         !Call this whenever the number of atoms in a molecule changes 
         integer :: newsize
         if(newsize.gt.natotmax) then
            natotmax=newsize
            call reallocate2d_c(alabelt,len(alabelt),natotmax,nmoltemplates)
            call reallocate2d_i(nbondt,natotmax,nmoltemplates)
         end if
      end subroutine reallocate_atot

      subroutine reallocate_bondt(newsize)
         !Call this whenever the number of bonds for an atom changes
         integer  :: newsize
         if(newsize.gt.nbondtmax) then
            nbondtmax=newsize
            call reallocate3d_i(iabondt,nbondtmax,natotmax,nmoltemplates)
            call reallocate3d_r(tbondt,nbondtmax,natotmax,nmoltemplates)
         end if
      end subroutine reallocate_bondt

      subroutine reallocate_rotgrpt(newsize)
         !Call this whenever the number of rotational groups in a molecule changes
         integer  :: newsize
         if(newsize.gt.ngrrottmax) then
            ngrrottmax=newsize
            call reallocate2d_c(rotgrpt,len(rotgrpt),ngrrottmax,nmoltemplates)
            call reallocate3d_i(labaxit,2,ngrrotmax,nmoltemplates)
            call reallocate2d_i(nrott,ngrrotmax,nmoltemplates)
            call reallocate2d_i(ndihedralangt,ngrrotmax,nmoltemplates)
         end if
      end subroutine reallocate_rotgrpt

      subroutine reallocate_labrott(newsize)
         !Call this whenever number of atoms in a rotational group changes
         integer :: newsize
         if(newsize.gt.nrottmax) then
            nrottmax=newsize
            call reallocate3d_i(labrott,nrottmax,ngrrotmax,nmoltemplates)
         end if
      end subroutine reallocate_labrott

      subroutine reallocate_dihedralangt(newsize)
         !Call this whenever the number of dihedral angles changes
         if(newsize.gt.ndihtmax) then
            ndihtmax=newsize
            call reallocate3d_r(dihedralangt,ndihtmax,ngrrotmax,nmoltemplates)
            call reallocate4d_i(labrotdiht,2,ndihtmax,ngrrotmax,nmoltemplates)
         end if
      end subroutines reallocate_dihedralangt

   end module moltemplates_mod


