   module molecule_coordinates_mod

      character(len=3), dimension(:), allocatable     :: chsymbolt
      integer, dimension(:), allocatable              :: centralmolecule
      real, dimension(:), allocatable                 :: rt,phit,thetat
      integer                                         :: nofatoms,nread,lenatomsequence
      character(len=256)                              :: xyzfilename,atomsequence

   contains

      subroutine get_coords_from_xyz(fname_run,lenfname_run)

      character(len=*)                                :: fname_run
      integer                                         :: lenfname_run
      real                                            :: xt,yt,zt,cosphi,sinphi   
      integer, parameter                              :: nworddim=201,nstsize=256
      character(len=nstsize)                          :: sometext
      character(len=3)                                :: tempchar
      integer                                         :: nwords,ncf(nworddim),ncl(nworddim)
      real, parameter                                 :: bohr=0.52917721092
      write(6,*) fname_run(1:lenfname_run)
! Attempt to read coordinates of central molecule atoms
      xyzfilename=fname_run(1:indexlastdot(fname_run))//'xyz'
      nread=10
      nofatoms=0
      lenatomsequence=0
      write(6,27) xyzfilename(1:len_trim(xyzfilename))
27    format(/'get_coords_from_xyz> Attempting to read file: ',a)
      open(nread,file=xyzfilename,status='old',iostat=ierr)
      if(ierr.eq.0) then
         read(nread,*,iostat=ierr) nofatoms
         read(nread,*) ! Blank line.
         if(nofatoms.gt.0) then
            if(allocated(chsymbolt)) deallocate(chsymbolt,centralmolecule,rt,phit,thetat)
            allocate(chsymbolt(nofatoms),centralmolecule(nofatoms),rt(nofatoms),phit(nofatoms),thetat(nofatoms))
            do it=1,nofatoms
               read(nread,'(a)') sometext
               call parse(sometext,nwords,ncf,ncl,10,nstsize,nworddim)
               write(6,*) nwords,sometext(1:len_trim(sometext))
               read(sometext(ncf(1):ncl(1)),'(a3)') chsymbolt(it)
               read(sometext(ncf(2):ncl(4)),*) xt,yt,zt
               ierr=1
               if(nwords.gt.4) read(sometext(ncf(5):ncl(5)),*,iostat=ierr) tempchar
               if(ierr.ne.0) tempchar=chsymbolt(it)
               centralmolecule(it)=min(1,index(sometext,'#CENTRALMOLECULE'))
               rt(it)=sqrt(xt*xt+yt*yt+zt*zt)
               thetat(it)=0.0
               phit(it)=0.0
               if(rt(it).gt.0.0) then
                  thetat(it)=acos(zt/rt(it))
                  cosphi=xt/rt(it)
                  sinphi=yt/rt(it)
                  phit(it)=atan2(sinphi,cosphi)
               end if
! Build the atom sequence to identify this molecule
               lentempchar=len_trim(tempchar)
               if(it.eq.1) then
                  atomsequence=tempchar(1:lentempchar)
                  lenatomsequence=lentempchar
               else if(lenatomsequence.lt.len(atomsequence)-lentempchar) then
                  atomsequence=atomsequence(1:lenatomsequence)//tempchar(1:lentempchar)
                  lenatomsequence=lenatomsequence+lentempchar
               endif
           end do
         end if
      end if
      close(nread)
      write(6,*) nofatoms,' atoms read from file ',xyzfilename(1:len_trim(xyzfilename))
      write(6,*) 'Atom sequence ',atomsequence(1:lenatomsequence)

      end subroutine get_coords_from_xyz

   end module molecule_coordinates_mod

