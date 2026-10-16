      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call randomise()

      stop
      end      

      subroutine randomise()

      use molecom_mod
c
c perform a random rotation and translation of all molecules in a cluster
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      character*256 readfile
      character*4 ext
      character*1 iaccept
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
499   format(/'randomise> No .ato file specified! - routine will exit'/)
               return
            endif
         end do
         readfile=sometext(ncf(1):ncl(1))
         ilenreadfile=ncl(1)-ncf(1)+1
      else
         readfile=sometext(ncf(2):ncl(2))
         ilenreadfile=ncl(2)-ncf(2)+1
      endif
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
      call normalise_coordinates
c
c save the current rotational step sizes, and substitute 2*pi rotations
c for the present purpose
c
      stepris=stepri
      steprs=stepr
c For rotational moves, only perform the randomisation if the relevant step size is .gt. 0
      if(stepris.gt.0.0) stepri=1.0
      if(steprs.gt.0.0) stepr=1.0
      do im=1,nmol
      
         if(.not.tethered(im)) then

c perform random translation of whole molecule 

            if(stepm.gt.0.0) then

               fa=(ran1()-0.5)
               fb=(ran1()-0.5)
               fc=(ran1()-0.5)

c Convert to Cartesian coordinates

c
c calculate corrected coordinates
c
               zadd(im)=fc*clatz
               yadd(im)=fb*blaty+fc*claty
               xadd(im)=fa*alatx+fb*blatx+fc*clatx

c ensure all molecules are in the central box ( + or - 0.5)

               call minimagec(xadd(im),yadd(im),zadd(im))
c               call minimagec(phix(im),phiy(im),phiz(im))

            endif
c perform random whole molecule rotation if requested
            if(stepr.gt.0.0) then
               call rotatm(im,stepr)
            end if
         end if
      end do
c perform random internal molecule rotation
!      if(ngrrot.gt.0.and.stepri.gt.0.0) then
!            do ig=1,ngrrot
!               if(rotgrp(ig).eq.'ROT') call rotati(ig) !Only rotate free rotation groups at this stage
!            end do
!      endif
      write(6,193) ngrrot,nmol
193   format('randomise> ',i5,' rotational groups in ',i5,' molecules ')
c
c set the rotational step sizes back to what they were on input
c
      stepri=stepris
      stepr=steprs
c
c write out the new molecule coordinates
c
      write(6,*) dmax,dmaxmol
      call unnormalise_coordinates
      call moloutc(fnameato)
      return
      end
