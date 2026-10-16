      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call mixato()

      end

      subroutine mixato()

      use reallocation_routines
      use molecom_mod
c
c mixes up to mcomp .ato files. If any of these contain more than one molecule, mixato will fail.
c
      include 'epsr_dimension.inc'
      parameter(mmix=mcomp,ma=near)
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      character*256 readfile
      character*4 ext
      character*1 iaccept

      integer, dimension(:), allocatable        :: nmincl,natos,ngrots,nmoltypes!(mmix)
      integer, dimension(:,:), allocatable      :: nbondss,nrots,ndihedralangs!(near,mmix)!(ma,mmix)
      integer, dimension(:,:,:), allocatable    :: iabonds,labaxis,labrots!(ma,near,mmix)(near,2,mmix)(near,near,mmix)
      integer, dimension(:,:,:,:), allocatable  :: labrotdihs!(2,near,near,mmix)

c Save coordinates and bonds

      real, dimension(:), allocatable        :: xadds,yadds,zadds,phixs,phiys,phizs!(mmix)
      real, dimension(:,:), allocatable      :: xas,yas,zas!(ma,mmix)
      real, dimension(:,:,:), allocatable    :: tbonds,dihedralangs!(ma,near,mmix)

c Save step sizes and vibrational temperature

      integer natomax,ngrotmax,nrotmax,ndihmax
      real stepmis,stepris,steprs,stepms,tempbs,natomrnew

      character(len=3), dimension(:,:), allocatable   :: alabels,rotgrps!(ma,mmix)(near,mmix)

      ierr=1
      if(nwords.gt.1) then
         read(sometext(ncf(2):ncl(2)),*,iostat=ierr) nmix
      else
         ierr=1
      endif
      do while(ierr.ne.0)
         write(6,580)
580   format(/'mixato> How many .ATO files do you want to mix? ',$)
         read(5,*,iostat=ierr) nmix
      end do
      call reallocate1d_i(nmincl,nmix)
      call reallocate1d_i(natos,nmix)
      call reallocate1d_i(ngrots,nmix)
      call reallocate1d_i(nmoltypes,nmix)
      call reallocate1d_r(xadds,nmix)
      call reallocate1d_r(yadds,nmix)
      call reallocate1d_r(zadds,nmix)
      call reallocate1d_r(phixs,nmix)
      call reallocate1d_r(phiys,nmix)
      call reallocate1d_r(phizs,nmix)

      nm1=0
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
         imix=imix+1
         write(6,699) imix
699   format(//'mixato> Search for .ato file ',i2)
c
c generate the name of each .ato file in turn
c
         ext='.ato'
         ilenext=4
         nwords=0
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)
c
c if no .ato files are available in this directory, abandon this attempt at mixato
c
            if(iaccept.ne.'y') then
               write(6,700)
700   format(/'mixato> No .ato files found in the current directory - try another directory.'/)
               return
            endif
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
499   format(/'mixato> No .ato file specified - try again!'/)
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
604   format(/'mixato> ',a,'.ato contains more than 1 type of molecule'
     *//'mixato> This is not allowed - start again')
               return
            endif
            ntypefirst=ntype
            nmoltypefirst=nmoltype
            write(6,36) ntypefile,fnameato(1:ilenfnameato)
36    format(/1x,i3,' new atom types in file ',a/)
            do it=itype(1),itype(natom)
               write(6,33) it,type(it)
            end do
33    format(1x,'Atom type ',i2,' has label ',a3)
            write(6,581)
581   format(/'mixato> How many of these molecules do you want in the mixture? ',$)
            read(5,*) nmincl(imix)
c
c generate the extra molecules based on the last one read in
c
            nm1=nm1+nmincl(imix)
            if(nm1.gt.mmol) then
               write(6,946) mmol
946   format(1x,'No. of molecules requested .gt. than allowed dimension:',i6)
               return
            endif
c
c save the parameters of the first molecule
c
c Use the largest of the step sizes

            if(stepmi.gt.stepmis) stepmis=stepmi
            if(stepri.gt.stepris) stepris=stepri
            if(stepr.gt.steprs) steprs=stepr
            if(stepm.gt.stepms) stepms=stepm
            if(tempb.gt.tempbs) tempbs=tempb
            if(tempdih.gt.tempdihs) tempdihs=tempdih

            natos(imix)=nato(1)
            xadds(imix)=xadd(1)
            yadds(imix)=yadd(1)
            zadds(imix)=zadd(1)
            phixs(imix)=phix(1)
            phiys(imix)=phiy(1)
            phizs(imix)=phiz(1)
            if(nato(1).gt.natomax) then
               natomax=nato(1)
               call reallocate2d_c(alabels,3,natomax,nmix)
               call reallocate2d_r(xas,natomax,nmix)
               call reallocate2d_r(yas,natomax,nmix)
               call reallocate2d_r(zas,natomax,nmix)
               call reallocate2d_i(nbondss,natomax,nmix)
               call reallocate3d_i(iabonds,natomax,near,nmix)
               call reallocate3d_r(tbonds,natomax,near,nmix)
            end if
            do ia=naf(1),nal(1)
               alabels(ia,imix)=alabel(ia)
               xas(ia,imix)=xa(ia)
               yas(ia,imix)=ya(ia)
               zas(ia,imix)=za(ia)
               nbondss(ia,imix)=nbond(ia)
               if(nbond(ia).gt.0) then
                  do ib=1,nbond(ia)
                     iabonds(ia,ib,imix)=iabond(ia,ib)
                     tbonds(ia,ib,imix)=tbond(ia,ib)
                  end do
               endif
            end do
c
c rotational groups
c
            ngrots(imix)=ngrot(1)
            if(ngrot(1).gt.ngrotmax) then
               ngrotmax=ngrot(1)
               call reallocate2d_c(rotgrps,3,ngrotmax,nmix)
               call reallocate2d_i(nrots,ngrotmax,nmix)
               call reallocate3d_i(labaxis,ngrotmax,2,nmix)
               call reallocate2d_i(ndihedralangs,ngrotmax,nmix)
               call reallocate3d_r(dihedralangs,ndihmax,ngrotmax,nmix)
               call reallocate4d_i(labrotdihs,2,ndihmax,ngrotmax,nmix)
               call reallocate3d_i(labrots,ngrotmax,nrotmax,nmix)
            end if
            if(ngrot(1).gt.0) then
               do ig=1,ngrot(1)
                  rotgrps(ig,imix)=rotgrp(ig)
                  labaxis(ig,1,imix)=labaxi(ig,1)
                  labaxis(ig,2,imix)=labaxi(ig,2)
                  nrots(ig,imix)=nrot(ig)
                  if(nrot(ig).gt.nrotmax) then
                     nrotmax=nrot(ig)
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
!Find out if the supplied template file actually exists.
            open(35,file=moltemplate(nmoltype)(1:len_trim(moltemplate(nmoltype)))//'.ato',status='old',iostat=ierr)
            close(35)
            !write(6,*) 'mixato> Found or not? ',moltemplate(nmoltype)(1:len_trim(moltemplate(nmoltype)))//'.ato',' ',ierr
            if(ierr.ne.0) then
                write(6,*) 'mixato> Did not find template file: ',moltemplate(nmoltype)(1:len_trim(moltemplate(nmoltype)))//'.ato'
               !Create a name based on the name of the current file being read from
               moltemplate(nmoltype)=readfile(1:ilenreadfile)
            end if
            write(6,*)
            write(6,*) ' mixato> moltemplate file: ',moltemplate(nmoltype)(1:len_trim(moltemplate(nmoltype)))//'.ato'

         else
            ntype=ntypefirst
            nmoltype=nmoltypefirst
            write(6,599) fnameato
599   format(/'Problems with .ato file ',a,' Try again.'/)
            imix=imix-1
         endif
      end do
c
c O.K. so now generate the new molecular coordinates based on the values for the
c individual molecules read in
c
      ia=0
      im=0
      ig=0
      natom=0
      natomr=0
      natomcg=0
      natomrcg=0
      ngrrot=0
      nrottot=0
      do imix=1,nmix
         do it=1,nmincl(imix)
            if(im.lt.mmol) then
               im=im+1
            else
               write(6,601)
601   format(/'mixato> Tried to define more molecules than allowed by dimension'/)
               return
            endif
            xadd(im)=xadds(imix)
            yadd(im)=yadds(imix)
            zadd(im)=zadds(imix)
            nato(im)=natos(imix)
            nadd=ia
            naf(im)=ia+1
            do is=1,nato(im)
               if(ia.lt.matom) then
                  ia=ia+1
               else
                  write(6,602)
602   format(/'mixato> Tried to define more atoms than allowed by dimension'/)
                  return
               endif
               alabel(ia)=alabels(is,imix)
               call ftype(ia)
               natomcg=natomcg+natomsphere(itype(ia))*pseudofacmol(itype(ia))
c
c do not count atoms which begin with 'Q' or 'q' - these are assumed
c invisible to the scattering radiation and experiment (although they
c can be used to calculate other structural quantities if necessary
c
               if(qtype(ia).ne.1) then
                  natomr=natomr+1
                  natomrcg=natomrcg+natomsphere(itype(ia))*pseudofacmol(itype(ia))
               endif
               xa(ia)=xas(is,imix)
               ya(ia)=yas(is,imix)
               za(ia)=zas(is,imix)
               nbond(ia)=nbondss(is,imix)
               if(nbond(ia).gt.0) then
                  do ib=1,nbond(ia)
                     iabond(ia,ib)=iabonds(is,ib,imix)+nadd
                     tbond(ia,ib)=tbonds(is,ib,imix)
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
                  ig=ig+1
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
                        dihedralang(k,ig)=dihedralangs(k,ih,imix)
                        do j=1,2
                           labrotdih(j,k,ig)=labrotdihs(j,k,ih,imix)+nadd
                        end do
                     end do
                  end if
               end do
            endif
         end do
      end do
      nmol=im
      natom=ia
      ngrrot=ig
      write(6,34)
34    format(/'mixato> Give atomic number density (per A**3) of mixture: ',$)
      read(5,*) rhoatom
      dmax=natomrcg/rhoatom
      dmax=dmax**(1.0/3.0)
c
c normalise all the coordinates to this new value of dmax
c
      dmaxmol=sizefactor*dmax

! Set up the lattice dimensions for calculating the density in moloutc.

      cubic=.true.
      alat=dmax
      blat=dmax
      clat=dmax
      phib=90.0
      thetac=0.0
      phic=0.0
      alatx=alat
      blatx=blat*cosd(phib)
      blaty=blat*sind(phib)
      clatx=clat*sind(thetac)*cosd(phic)
      claty=clat*sind(thetac)*sind(phic)
      clatz=clat*cosd(thetac)


c Set the step sizes to the largest found in all the ato files read in

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
35    format(/'mixato> Type name of file to put mixture in: ',$)
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
