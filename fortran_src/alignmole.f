      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call alignmole()

      end      

      subroutine alignmole()
c
c Align a specified molecule along a particular axis
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      character(len=256)                  :: readfile
      character(len=4)                    :: ext
      character(len=1)                    :: iaccept,alignaxis
      real                                :: displacement
      logical                             :: useall

      integer                             :: niter,iter,naxisatom,j
      integer, dimension(:), allocatable  :: axisatom       !atom numbers used to define molecular axis
      real, dimension(:), allocatable         :: xr,yr,zr

c
c generate the name of the .ato file
c
      ext='.ato'
      ilenext=4
      if(nwords.lt.2) then
         nwords=0
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
499   format(/'alignmole> No .ato file specified! - routine will exit'/)
               return
            endif
         end do
         readfile=sometext(ncf(1):ncl(1))
         ilenreadfile=ncl(1)-ncf(1)+1
      else
         readfile=sometext(ncf(2):ncl(2))
         ilenreadfile=ncl(2)-ncf(2)+1
      endif
      if(nwords.gt.2) then
         alignaxis=sometext(ncf(3):ncf(3))
      else
         write(6,100,advance='no')
100   format(/'alignmole> Type axis to align to (x, y, or z): ')
         read(5,'(a1)') alignaxis
      end if
      if(nwords.gt.3) then
         read(sometext(ncf(4):ncl(4)),*) niter
      else
         write(6,101,advance='no')
101   format(/'alignmole> Give number of iterations: ')
         read(5,*) niter
      end if

! Get the atom numbers to be used to align the molecules. These numbers have to be common
! to all molecules in the box. Hence best to use this routine for boxes containing only
! one type of molecule.

      naxisatom=0
      nf=5
! See if the word ' all ' appears on its own in sometext - if so we will use all the atoms
! in each molecule to align it
      useall=.true.
      if(nwords.ge.nf) useall=index(sometext(ncf(nf):ncl(nf)),'all').gt.0
      if(allocated(axisatom)) deallocate(axisatom)
      do while (.not.allocated(axisatom).and..not.useall)
         naxisatom=nwords-nf+1
         if(naxisatom.ge.2) then !Must have at least 2 axis atoms
            allocate(axisatom(naxisatom))
            read(sometext(ncf(nf):ncl(nwords)),*,iostat=ierr) (axisatom(j),j=1,naxisatom)
            if(ierr.ne.0) then  !If an error occurred on read, then will have input values again
               deallocate(axisatom)
            endif
         end if
         if(.not.allocated(axisatom)) then
            write(6,104,advance='no')
104   format(/'alignmole> Give list of atom numbers to define molecular axis' /
     *,' (Type "all" to use all atoms in each molecule) : ')
            read(5,'(a)') sometext
            useall=index(sometext,'all').gt.0
            call parse(sometext,nwords,ncf,ncl,10,nstsize,nworddim)
            nf=1
         end if
      end do
c
c strip out any extension that was supplied - extension must be .ato, otherwise
c nothing is done.
c
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
c
c initialise number of atom and molecule types to zero
c
      ntype=0
      nmoltype=0
c
c read .ato file
c
      sizefactor=1.0
      call minitc(fnameato,rmaxcell,cellmax,upol,rho,amubar,amubar2
     *,ndist,ierr)
      if(rho.eq.0.0.or.ierr.ne.0) then
         write(6,599) fnameato
599   format(/'Problems with .ato file ',a/'Routine will exit now'/)
         return
      endif
      write(6,*) dmax,dmaxmol
      
!Now check that the specified axis atoms actually exist for all the molecules.
!If not we use all atoms for all molecules

      if(allocated(axisatom)) then
         im=0
         do while (im.lt.nmol.and.ierr.eq.0)
            im=im+1
            ia=0
            do while (ia.lt.naxisatom.and.ierr.eq.0)
               ia=ia+1
               if(axisatom(ia).lt.1.or.axisatom(ia).gt.nato(im)) ierr=1
            end do
         end do
         if(ierr.ne.0) then
            deallocate(axisatom)
            naxisatom=0
            useall=.true.
         end if
      end if
      if(naxisatom.gt.0) write(6,*) naxisatom,(axisatom(j),j=1,naxisatom)

c Save the rotational step size

      steprs=stepr
      if(stepr.le.0.0) then
         stepr=0.1
      end if

      do im=1,nmol
            if(.not.allocated(xr)) then
                allocate(xr(nato(im)),yr(nato(im)),zr(nato(im)))
            else if (size(xr).lt.nato(im)) then
                deallocate(xr,yr,zr)
                allocate(xr(nato(im)),yr(nato(im)),zr(nato(im)))
            end if
         ntrials=0
         nreject=0
         naccum=niter/10
         iter=0

c Calculate the initial displacement from the specified axis

         energold=displacement(im,alignaxis,useall,axisatom,naxisatom)
         write(6,103) im,energold
103   format(/'alignmole> Molecule: ',i4,' Initial energy: ',e12.5)
         do while (iter.lt.niter)
            iter=iter+1
            ntrials=ntrials+1
c Random rotation about any axis
            call rotatm(im,stepr)
            energnew=displacement(im,alignaxis,useall,axisatom,naxisatom)
            if(energnew.gt.energold) then
c Reject the move
               nreject=nreject+1
c Restore the coordinates
               nat=naf(im)-1
               do ia=1,nato(im)
                  nat=nat+1
                  xa(nat)=xr(ia)
                  ya(nat)=yr(ia)
                  za(nat)=zr(ia)
               end do
            else
c If accepted, save the energy
               energold=energnew
            end if
c Write out a message and reappraise step size when appropriate
            if(naccum*(iter/naccum).eq.iter) then
               rat=real(nreject)/real(ntrials)
               write(6,102) iter,energold,rat
102   format(/'alignmole> Iteration ',i7,' Energy: ',e12.5,' Reject rate: ',f10.5)
               rat=rat/0.75
               if(rat.gt.1.1) rat=1.1
               if(rat.lt.0.9) rat=0.9
               stepr=stepr/rat
            endif
         end do
      end do
c
c set the rotational step sizes back to what they were on input
c
      stepr=steprs
c
c write out the new molecule coordinates
c
      call moloutc(fnameato)
      return
      end
      
      function displacement(im,alignaxis,useall,axisatom,naxisatom)

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      
      real                    :: displacement,sum,diff
      character(len=1)        :: alignaxis
      integer                 :: im,ia,ir,naxisatom
      integer, dimension(*)   :: axisatom
      logical                 :: useall      

      sum=0.0
      displacement=0.0
!      write(6,*) useall,naxisatom,im,nato(im)
      if(useall.or.naxisatom.le.0) then
         if(alignaxis.eq.'x') then
            do ia=naf(im),nal(im)
               diff=ya(ia)
               sum=sum+diff*diff
               diff=za(ia)
               sum=sum+diff*diff
            end do
         else if(alignaxis.eq.'y') then
            do ia=naf(im),nal(im)
               diff=xa(ia)
               sum=sum+diff*diff
               diff=za(ia)
               sum=sum+diff*diff
            end do
         else
            do ia=naf(im),nal(im)
               diff=xa(ia)
               sum=sum+diff*diff
               diff=ya(ia)
               sum=sum+diff*diff
            end do
         end if
         displacement=sum/nato(im)
      else
         if(alignaxis.eq.'x') then
            do ir=1,naxisatom
               ia=axisatom(ir)+naf(im)-1
               diff=ya(ia)
               sum=sum+diff*diff
               diff=za(ia)
               sum=sum+diff*diff
            end do
         else if(alignaxis.eq.'y') then
            do ir=1,naxisatom
               ia=axisatom(ir)+naf(im)-1
               diff=xa(ia)
               sum=sum+diff*diff
               diff=za(ia)
               sum=sum+diff*diff
            end do
         else
            do ir=1,naxisatom
               ia=axisatom(ir)+naf(im)-1
               diff=xa(ia)
               sum=sum+diff*diff
               diff=ya(ia)
               sum=sum+diff*diff
           end do
         end if
         displacement=sum/naxisatom
      endif
!      write(6,*) im,displacement
      return
      end
         
      
      
