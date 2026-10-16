      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call addato()

      end      

      subroutine addato()
      
      use reallocation_routines
      use molecom_mod
c
c adds up to 20 .ato files to an existing .ato file. If any of these contain more than one molecule, addato will
c use only the first molecule encountered. Only files containing one type
c of molecule can be added.
c The box size of the .ato file to which these molecules cannot be changed, so it is up to the user to add
c the correct number of molecules to get the correct density. The added molecules will be placed at random, except
c that they will not overlap with any tethered molecules - the degree of overlap is determined by dcore which is
c read in from the file being added to.
c
      include 'epsr_dimension.inc'
      parameter(mmix=20,ma=near)
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      character*256 readfile
      character*4 ext
      character*1 iaccept

      integer, dimension(:), allocatable        :: nmincl,natos,ngrots!(mmix)
      integer, dimension(:,:), allocatable      :: nbondss,nrots,itypes,ndihedralangs!(near,mmix)!(ma,mmix)
      integer, dimension(:,:,:), allocatable    :: iabonds,labaxis,labrots!(ma,near,mmix)(near,2,mmix)(near,near,mmix)
      integer, dimension(:,:,:,:), allocatable  :: labrotdihs!(2,near,near,mmix)

c Save coordinates and bonds

      real, dimension(:), allocatable        :: xadds,yadds,zadds,phixs,phiys,phizs!(mmix)
      real, dimension(:,:), allocatable      :: xas,yas,zas!(ma,mmix)
      real, dimension(:,:,:), allocatable    :: tbonds,dihedralangs!(ma,near,mmix)

c Save step sizes and vibrational temperature

      integer natomax,ngrotmax,nrotmax,ndihmax
      real stepmis,stepris,steprs,stepms,tempbs,rminsqr,xdif,ydif,zdif,xdifc,ydifc,zdifc,pi,piconv

      character(len=3), dimension(:,:), allocatable   :: alabels,rotgrps!(ma,mmix)(near,mmix)

      logical                                :: overlap
      
      integer, dimension(:), allocatable     :: tnumber
      integer                                :: it,itimes,ntimes,itn,ntn,im,jm,iat,jat
      
      logical, dimension(:), allocatable     :: tethereds

      ierr=1
      if(nwords.gt.1) then
         read(sometext(ncf(2):ncl(2)),*,iostat=ierr) nmix
      else
         ierr=1
      endif
      do while(ierr.ne.0)
         write(6,580)
580   format(/'addato> How many .ATO files do you want to add? ',$)
         read(5,*,iostat=ierr) nmix
      end do
      call reallocate1d_i(nmincl,nmix)
      call reallocate1d_i(natos,nmix)
      call reallocate1d_i(ngrots,nmix)
      call reallocate1d_r(xadds,nmix)
      call reallocate1d_r(yadds,nmix)
      call reallocate1d_r(zadds,nmix)
      call reallocate1d_r(phixs,nmix)
      call reallocate1d_r(phiys,nmix)
      call reallocate1d_r(phizs,nmix)
      call reallocate1d_l(tethereds,nmix)

      nm1=0
      nmixmax=nmix
      natomax=0
      ngrotmax=0
      nrotmax=0
      ndihmax=0
c
c set the number of atom and molecule types to zero at the outset
c
      ntype=0
      nmoltype=0
      ntypefirst=ntype
      nmoltypefirst=nmoltype

c Set the maximum step sizes to zero

      stepmis=0
      stepris=0
      steprs=0
      stepms=0
      tempbs=0
      tempdihs=0

      imix=0
      do while(imix.lt.nmix)
         write(6,699)
699   format(//'addato> Search for .ato file to be ADDED ')
c
c generate the name of each .ato file in turn
c
         ext='.ato'
         ilenext=4
         nwords=0
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)
c
c if no .ato files are available in this directory, abandon this attempt at addato
c
            if(iaccept.ne.'y') then
               write(6,700)
700   format(/'addato> No .ato files found in the current directory - try another directory.'/)
               return
            endif
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
499   format(/'addato> No .ato file specified - try again!'/)
            endif
         end do
         readfile=sometext(ncf(1):ncl(1))
         ilenreadfile=ncl(1)-ncf(1)+1
c
c strip out any extension that was supplied - extension must be .ato, otherwise
c nothing is done.
c
         ilendot=indexlastdot(readfile)-1
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
c read .ato file
c
         sizefactor=1.0
         call minitc(fnameato,rmaxcell,cellmax,upol,rho,amubar,amubar2
     *,ndist,ierr)
         if(ierr.eq.0) then
            ntypefile=ntype-ntypefirst
            nmoltypefile=nmoltype-nmoltypefirst
c
c because we have to read the entire file to get all the relevant information, a
c .ato file which contains more than 1 molecule type cannot be mixed
c
            if(nmoltypefile.gt.1) then
               write(6,604) readfile(1:ilenreadfile)
604   format(/'addato> ',a,'.ato contains more than 1 type of molecule'
     *//'addato> This is not allowed - try again')
            else
               ntypefirst=ntype
               nmoltypefirst=nmoltype
               write(6,36) ntypefile,fnameato(1:ilenfnameato)
36    format(/1x,i3,' new atom types in file ',a/)
               do it=itype(1),itype(natom)
                  write(6,33) it,type(it)
               end do
33    format(1x,'Atom type ',i2,' has label ',a3)

! If a molecule type is tethered, we use all the molecules, otherwise we use only
! the first molecule and specify the number of times we want this molecule

               if(moltethered(nmoltype)) then 
                  nmix=nmix+nmol-1
                  nmoladd=1
                  immax=nmol
               else
                  write(6,581,advance='no') 
581   format(/'addato> How many of these molecules do you want in the mixture? ')
                  read(5,*) nmoladd
                  immax=1
               endif
               im=0
               irgref=0
               do while (im.lt.immax)
                  im=im+1
                  imix=imix+1
                  if(imix.gt.nmixmax) then
                     nmixmax=nmixmax+10
                     call reallocate1d_i(nmincl,nmixmax)
                     call reallocate1d_i(natos,nmixmax)
                     call reallocate1d_i(ngrots,nmixmax)
                     call reallocate1d_r(xadds,nmixmax)
                     call reallocate1d_r(yadds,nmixmax)
                     call reallocate1d_r(zadds,nmixmax)
                     call reallocate1d_r(phixs,nmixmax)
                     call reallocate1d_r(phiys,nmixmax)
                     call reallocate1d_r(phizs,nmixmax)
                     call reallocate1d_l(tethereds,nmixmax)
                  end if
                  nmincl(imix)=nmoladd
                  nm1=nm1+nmincl(imix)
                  natos(imix)=nato(im)
                  xadds(imix)=xadd(im)
                  yadds(imix)=yadd(im)
                  zadds(imix)=zadd(im)
                  phixs(imix)=phix(im)
                  phiys(imix)=phiy(im)
                  phizs(imix)=phiz(im)
                  tethereds(imix)=tethered(im)
                  iaadd=naf(im)-1
                  if(nato(im).gt.natomax) then
                     natomax=nato(im)
                     call reallocate2d_c(alabels,3,natomax,nmixmax)
                     call reallocate2d_i(itypes,natomax,nmixmax)
                     call reallocate2d_r(xas,natomax,nmixmax)
                     call reallocate2d_r(yas,natomax,nmixmax)
                     call reallocate2d_r(zas,natomax,nmixmax)
                     call reallocate2d_i(nbondss,natomax,nmixmax)
                     call reallocate3d_i(iabonds,natomax,near,nmixmax)
                     call reallocate3d_r(tbonds,natomax,near,nmixmax)
                  end if
                  do ir=1,natos(imix)
                     ia=iaadd+ir
                     alabels(ir,imix)=alabel(ia)
                     itypes(ir,imix)=itype(ia)
                     xas(ir,imix)=xa(ia)
                     yas(ir,imix)=ya(ia)
                     zas(ir,imix)=za(ia)
                     nbondss(ir,imix)=nbond(ia)
                     if(nbond(ia).gt.0) then
                        do ib=1,nbond(ia)
                           iabonds(ir,ib,imix)=iabond(ia,ib)
                           tbonds(ir,ib,imix)=tbond(ia,ib)
                        end do
                     endif
                  end do
c
c rotational groups
c
                  ngrots(imix)=ngrot(im)
                  if(ngrot(im).gt.ngrotmax) then
                     ngrotmax=ngrot(im)
                     call reallocate2d_c(rotgrps,3,ngrotmax,nmixmax)
                     call reallocate2d_i(nrots,ngrotmax,nmixmax)
                     call reallocate3d_i(labaxis,ngrotmax,2,nmixmax)
                     call reallocate2d_i(ndihedralangs,ngrotmax,nmix)
                     call reallocate3d_r(dihedralangs,ndihmax,ngrotmax,nmix)
                     call reallocate4d_i(labrotdihs,2,ndihmax,ngrotmax,nmix)
                     call reallocate3d_i(labrots,ngrotmax,nrotmax,nmix)
                  end if
                  if(ngrot(im).gt.0) then
                     do ig=1,ngrot(im)
                        irgref=irgref+1
                        rotgrps(ig,imix)=rotgrp(irgref)
                        labaxis(ig,1,imix)=labaxi(irgref,1)
                        labaxis(ig,2,imix)=labaxi(irgref,2)
                        nrots(ig,imix)=nrot(irgref)
                        if(nrot(irgref).gt.nrotmax) then
                           nrotmax=nrot(irgref)
                           call reallocate3d_i(labrots,ngrotmax,nrotmax,nmix)
                        end if
                        do j=1,nrot(ig)
                           labref=rotf(ig)+j-1
                           labrots(ig,j,imix)=labrot(labref)
                        end do
!Check if there are dihedral angles defined
                        if(rotgrp(ig).eq.'DIH') then
                           ndihedralangs(ig,imix)=ndihedralang(ig)
                           if(ndihedralang(ig).gt.ndihmax) then
                              ndihmax=ndihedralang(ig)
                              call reallocate3d_r(dihedralangs,ndihmax,ngrotmax,nmix)
                              call reallocate4d_i(labrotdihs,2,ndihmax,ngrotmax,nmix)
                           end if
                           do k=1,ndihedralang(ig)
                              dihedralangs(k,ig,imix)=dihedralang(k,ig)
                              do j=1,2
                                 labrotdihs(j,k,ig,imix)=labrotdih(j,k,ig)
                              end do
                           end do
                        end if
                     end do
                  endif
               end do
            
!Save the template file name for output later
!First find out if the supplied template file actually exists.
               open(35,file=moltemplate(nmoltype)//'.ato',status='old',iostat=ierr)
               close(35)
!               write(6,*) 'addato> Found or not? ',moltemplate(nmoltype)(1:len_trim(moltemplate(nmoltype)))//'.ato',' ',ierr
               if(ierr.ne.0) then
                 !Create a name based on the name of the current file being read from 
                  moltemplate(nmoltype)=readfile(1:ilenreadfile)
               end if
               write(6,*)
               write(6,*) ' addato> moltemplate file: ',moltemplate(nmoltype)(1:len_trim(moltemplate(nmoltype)))//'.ato'
               
c Use the largest of the step sizes

               if(stepmi.gt.stepmis) stepmis=stepmi
               if(stepri.gt.stepris) stepris=stepri
               if(stepr.gt.steprs) steprs=stepr
               if(stepm.gt.stepms) stepms=stepm
               if(tempb.gt.tempbs) tempbs=tempb
               if(tempdih.gt.tempdihs) tempdihs=tempdih

               if(nm1.gt.mmol) then
                  write(6,946) mmol
946   format(1x,'No. of molecules requested .gt. than allowed dimension:',i6)
                  return
               endif
            end if
         else
            ntype=ntypefirst
            nmoltype=nmoltypefirst
            write(6,599) fnameato
599   format(/'addato> Problems with .ato file ',a,' Try again.'/)
         endif
      end do

!Now read the .ato file you want to ADD these molecules to

      do while (imix.lt.nmix+1)
         write(6,698)
698   format(//'addato> Search for .ato file to be ADDED to')
         nwords=0
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)
c
c if no .ato files are available in this directory, abandon this attempt at addato
c
            if(iaccept.ne.'y') then
               write(6,700)
               return
            endif
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
            endif
         end do
         readfile=sometext(ncf(1):ncl(1))
         ilenreadfile=ncl(1)-ncf(1)+1
c
c strip out any extension that was supplied - extension must be .ato, otherwise
c nothing is done.
c
         ilendot=indexlastdot(readfile)-1
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
c read .ato file
c
         sizefactor=1.0
         call minitc(fnameato,rmaxcell,cellmax,upol,rho,amubar,amubar2
     *,ndist,ierr)
         if(ierr.eq.0) then
            imix=imix+1
            ntypefile=ntype-ntypefirst
            nmoltypefile=nmoltype-nmoltypefirst
            write(6,36) ntypefile,fnameato(1:ilenfnameato)
            do it=itype(1),itype(natom)
               write(6,33) it,type(it)
            end do

c Use the largest of the step sizes

            if(stepmi.gt.stepmis) stepmis=stepmi
            if(stepri.gt.stepris) stepris=stepri
            if(stepr.gt.steprs) steprs=stepr
            if(stepm.gt.stepms) stepms=stepm
            if(tempb.gt.tempbs) tempbs=tempb
            if(tempdih.gt.tempdihs) tempdihs=tempdih

         else
            ntype=ntypefirst
            nmoltype=nmoltypefirst
            write(6,599) fnameato
         endif
      end do

c
c O.K. so now generate the new molecular coordinates based on the values for the
c individual molecules read in
c
! First make a list of all tethered molecules: these will be used to test for overlaps

      ntn=0
      do im=1,nmol
         if(tethered(im)) ntn=ntn+1
      end do
      !Gather the tethered molecule numbers
      if(allocated(tnumber)) deallocate(tnumber)
      allocate(tnumber(ntn))
      itn=0
      do im=1,nmol
         if(tethered(im)) then
            itn=itn+1
            tnumber(itn)=im
         end if
      end do
      write(6,*) 'addato> Number of tethered molecules at start = ',ntn
      
      !To be safe normalise the existing coordinates, add normalised coordinates, then unnormalise them
      ! at the end.
      call normalise_coordinates()
      pi=4.0*atan(1.0)
      piconv=pi/180.0
         
      im=nmol
      ia=natom
      ig=ngrrot
      !Box distance squared to test for overlaps
      rminsqr=1.0/(dmax*dmax)
      do imix=1,nmix
         do iloop=1,nmincl(imix)
            if(im.lt.mmol) then
               im=im+1
            else
               write(6,601)
601   format(/'addato> Tried to define more molecules than allowed by dimension'/)
               return
            endif
            if(.not.tethereds(imix)) then
               !Set the tethered flag to false
               tethered(im)=.false.
               overlap=.true.
               itimes=0
               !Allow 1000 attempts before giving up
               ntimes=1000
               do while (overlap.and.itimes.lt.ntimes)
                  itimes=itimes+1
                  !Choose a position at random
                  fa=(ran1()-0.5)
                  fb=(ran1()-0.5)
                  fc=(ran1()-0.5)

c Convert to Cartesian coordinates. Note that by using the a,b,c axes formalism we avoid the need to use
c normalised coordinates for the mininum image convention.

                  zadd(im)=fc*clatz
                  yadd(im)=fb*blaty+fc*claty
                  xadd(im)=fa*alatx+fb*blatx+fc*clatx
                  phix(im)=0.0
                  phiy(im)=0.0
                  phiz(im)=0.0

c ensure all molecules are in the central box ( + or - 0.5)

                  call minimagec(xadd(im),yadd(im),zadd(im))
                  !Step through the tethered atoms test for overlap
                  itn=0
                  overlap=.false.
                  do while(.not.overlap.and.itn.lt.ntn)
                     itn=itn+1
                     jm=tnumber(itn)
                     xdifc=xadd(im)-xadd(jm)
                     ydifc=yadd(im)-yadd(jm)
                     zdifc=zadd(im)-zadd(jm)
                     !Compare the first atoms of each molecule for overlaps
                     iat=0
                     do while(iat.lt.natos(imix).and..not.overlap)
                        iat=iat+1
                        it=itypes(iat,imix)
                        jat=naf(jm)-1
                        do while (jat.lt.nal(jm).and..not.overlap)
                           jat=jat+1
                           jt=itype(jat)
                           xdif=xdifc+xas(iat,imix)/dmaxmol-xa(jat)
                           ydif=ydifc+yas(iat,imix)/dmaxmol-ya(jat)
                           zdif=zdifc+zas(iat,imix)/dmaxmol-za(jat)
                           call minimagec(xdif,ydif,zdif)
                           dtest=xdif*xdif+ydif*ydif+zdif*zdif
                           overlap=ecor(npairref(it,jt)).gt.0.0.and.dtest.lt.rminsqr*dscor(npairref(it,jt))
                           if (overlap) write(6,*) im,jm,iat,jat,it,jt,npairref(it,jt),dtest,dscor(npairref(it,jt)),rminsqr
                        end do
                     end do
                  end do
               end do
            else
               itimes=1
               ntimes=1
               overlap=.false.
               tethered(im)=tethereds(imix)
               xadd(im)=xadds(imix)/dmax
               yadd(im)=yadds(imix)/dmax
               zadd(im)=zadds(imix)/dmax
               phix(im)=phixs(imix)/dmax
               phiy(im)=phiys(imix)/dmax
               phiz(im)=phizs(imix)/dmax
            end if

            if(overlap) then

               write(6,697) ntimes
697   format(/'addato> Could find no insertion point after ',i4,' attempts.')
               im=im-1

            else
                  
               nato(im)=natos(imix)
               nadd=ia
               naf(im)=ia+1
               do is=1,nato(im)
                  if(ia.lt.matom) then
                     ia=ia+1
                  else
                     write(6,602)
602   format(/'addato> Tried to define more atoms than allowed by dimension'/)
                     return
                  endif
                  alabel(ia)=alabels(is,imix)
                  write(6,*) ia,is,imix,alabels(is,imix)
                  call ftype(ia)
c
c do not count atoms which begin with 'Q' or 'q' - these are assumed 
c invisible to the scattering radiation and experiment (although they
c can be used to calculate other structural quantities if necessary
c
                  if(qtype(ia).ne.1) then
                     natomr=natomr+1
                  endif
                  xa(ia)=xas(is,imix)/dmaxmol
                  ya(ia)=yas(is,imix)/dmaxmol
                  za(ia)=zas(is,imix)/dmaxmol
                  nbond(ia)=nbondss(is,imix)
                  if(nbond(ia).gt.0) then
                     do ib=1,nbond(ia)
                        iabond(ia,ib)=iabonds(is,ib,imix)+nadd
                        tbond(ia,ib)=tbonds(is,ib,imix)/dmaxmol
                     end do
                  endif
               end do
               nal(im)=ia
c
c rotational groups
c
               ngrot(im)=ngrots(imix)
               if(ngrot(im).gt.0) then
                  ngrrot=ngrrot+ngrot(im)
                  call testsize_rotgrps()
                  do ih=1,ngrot(im)
                     if(ig.lt.matom) then
                        ig=ig+1
                     else
                        write(6,603)
603   format(/'addato> Tried to define more rotational groups than allowed by dimension'/)
                        return
                     endif
c save molecule in which this rotational group occurs
                     mrot(ig)=im
                     rotgrp(ig)=rotgrps(ih,imix)
c for rotational groups read atom numbers of two atoms to form axis of rotation
                     do j=1,2
                        labaxi(ig,j)=labaxis(ih,j,imix)+nadd
                     end do
c read number of atoms in group and their numbers
                     nrot(ig)=nrots(ih,imix)
                     rotf(ig)=nrottot+1
                     nrottot=nrottot+nrot(ig)
                     rotl(ig)=nrottot
                     call testsize_labrot()
                     do j=1,nrot(ig)
                        labref=rotf(ig)+j-1
                        labrot(labref)=labrots(ih,j,imix)+nadd
                     end do
                     if(rotgrp(ig).eq.'DIH') then
                        ndihedralang(ig)=ndihedralangs(ih,imix)
                        do k=1,ndihedralang(ig)
                           dihedralang(k,ig)=dihedralangs(k,ih,imix)*piconv
                           do j=1,2
                              labrotdih(j,k,ig)=labrotdihs(j,k,ih,imix)+nadd
                           end do
                        end do
                     end if
                  end do
               endif
c perform random whole molecule rotation
               if(.not.tethered(im)) then
                  steprs=stepr
                  if(steprs.gt.0.0) stepr=1.0
                  call rotatm(im,stepr)
                  stepr=steprs
               end if
            end if
         end do
      end do
      nmol=im
      natom=ia

!Unnormalise the coordinates

      call unnormalise_coordinates()

c Leave the step sizes to those found in last file

      stepmi=stepmis
      stepri=0.3
      stepr=0.3
      stepm=1.0
      tempb=tempbs
      tempdih=tempdihs
c
c get the new filename to put the mixture in
c
      nwords=0
      do while (nwords.le.0)
         write(6,35)
35    format(/'addato> Type name of file to put mixture in: ',$)
         read(5,101,iostat=ierr) sometext
101   format(a)
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         if(nwords.le.0) then
            write(6,499)
         endif
      end do
      readfile=sometext(ncf(1):ncl(1))
      ilenreadfile=ncl(1)-ncf(1)+1
c
c strip out any extension that was supplied - extension must be .ato, otherwise
c nothing is done.
c
      ilendot=index(readfile,'.')-1
      if(ilendot.gt.0) then
         readfile=readfile(1:ilendot)
         ilenreadfile=ilendot
      endif
      fnameato=directory_name(1:nchardir)//readfile(1:ilenreadfile)//'.ato'
      ilenfnameato=nchardir+ilenreadfile+4
c
c write the new .ato file
c
      call moloutc(fnameato)
      return
      end   
