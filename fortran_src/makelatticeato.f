      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call makelatticeato()

      end      

      subroutine makelatticeato()

      use reallocation_routines
      use molecom_mod
c
c Make an EPSR18 ato file from an exising ato file by creating the specified number of unit
c cells along the a,b,c axes.
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'

      character(len=256)                           :: readfile
      integer                                      :: ilenreadfile,ierr,i,nwordstest
      real                                         :: pi,radconv,degconv
c      real, dimension(:), allocatable              :: xunit,yunit,zunit !x,y,z coordinates in unit cell
      integer                                      :: nunitcell,na,nb,nc,ka,kb,kc,jb !Number of atoms in unit cell and number of steps along
      !a, b and c axes
c      character(len=3), dimension(:), allocatable  :: atsymbolunit !Atom symbols in unit cell
      character(len=4)                             :: ext
      character(len=1)                             :: iaccept
      real                                         :: acent,bcent,ccent,fa,fb,fc,xav,yav,zav,xbox,ybox,zbox
      real  xadds(mmol),yadds(mmol),zadds(mmol),phixs(mmol),phiys(mmol),phizs(mmol)
      character(len=256) :: text
      logical :: convert      
      pi=4.0*atan(1.0)
      radconv=pi/180.0
      degconv=180.0/pi
c
c generate the name of the file containing the unit cell coordinates
c
      if(nwords.lt.2) then
         nwords=0
         nwordstest=1

c generate the name of each .ato file in turn

         ext='.ato'
         ilenext=4
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)

c if no .ato files are available in this directory, abandon this attempt at dockato

            if(iaccept.ne.'y') then
               write(6,700)
700   format(/'makelatticeato> No .ato files found in the current directory - try another directory.'/)
               return
            endif
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
499   format(/'makelaticeato> No file specified! - try again!'/)
            endif
         end do
         readfile=sometext(ncf(1):ncl(nwords))
         ilenreadfile=ncl(nwords)-ncf(1)+1
      else
         readfile=sometext(ncf(2):ncl(2))
         ilenreadfile=ncl(2)-ncf(2)+1
      endif
      if(nwords.gt.2) then
         read(sometext(ncf(3):ncl(3)),*) na
         if(nwords.gt.3) then
            read(sometext(ncf(4):ncl(4)),*) nb
            if(nwords.gt.4) then
               read(sometext(ncf(5):ncl(5)),*) nc
               iread=nwordstest+3
            else
               write(6,606,advance='no')
606   format(/'makelattice> Give number of unit cells along c-axis: ')
               read(5,*) nc
            end if
         else
            write(6,607,advance='no')
607   format(/'makelattice> Give number of unit cells along b and c-axes: ')
            read(5,*) nb,nc
         end if
      else
         write(6,608,advance='no')
608   format(/'makelattice> Give number of unit cells along a, b and c-axes: ')
         read(5,*) na,nb,nc
      end if
c
c create the full name of the input file by adding on the directory name if not present
c
      call get_full_filename(readfile,ilenreadfile,fnameato,ilenfnameato)
      if(index(fnameato(1:ilenfnameato),'.ato').eq.0) then
         fnameato=fnameato(1:ilenfnameato)//'.ato'
         ilenfnameato=ilenfnameato+4
      end if
c
c get the filename to put the mixture in
c
      if(nwords.gt.5) then 
         readfile=sometext(ncf(6):ncl(6))
         ilenreadfile=ncl(6)-ncf(6)+1
      else
         nwords=0
         iread=0
         do while (nwords.le.0)
            write(6,35)
35    format(/'makelatticeato> Type name of file to put lattice in: ',$)
            read(5,101,iostat=ierr) sometext
101   format(a)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
            endif
         end do
         readfile=sometext(ncf(1):ncl(nwords))
         ilenreadfile=ncl(nwords)-ncf(1)+1
      end if
! Decide if to convert atoms to molecules
      if(nwords.gt.6) then
         sometext=sometext(ncf(7):ncl(7))
      else
         write(6,613,advance='no')
613   format('makelatticeato> Convert molecules to atoms? (y or n) ')
         read(5,'(a)') sometext
!      write(sometext,'(a)') 'y'
      end if
      call parse(text,nwords,ncf,ncl,5,nstsize,nworddim)
      write(6,*) nwords,' ',sometext(ncf(1):ncl(1))
      convert=sometext(ncf(1):ncl(1)).eq.'y'.or.sometext(ncf(1):ncl(1)).eq.'Y'
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
         write(6,599) fnameato(1:ilenfnameato)
599   format(/'Problems with .ato file ',a/'Routine will exit now'/)
         return
      endif
      write(6,*) dmax,dmaxmol
      write(6,603) alat,blat,clat,phib,thetac,phic
603   format('makelatticeato> Lattice constants, angles and derived angles:'
     13(/3(1x,f10.5)))
   
c Save the coordinates of these molecules

      do im=1,nmol
         if(tethered(im)) then !For tethered we set coordinates to tethering positions to avoid problems with periodic boundary conditions
            xadds(im)=phix(im)
            yadds(im)=phiy(im)
            zadds(im)=phiz(im)
         else
            xadds(im)=xadd(im)
            yadds(im)=yadd(im)
            zadds(im)=zadd(im)
         end if
         phixs(im)=phix(im)
         phiys(im)=phiy(im)
         phizs(im)=phiz(im)
      end do
      !Get the number of molecules in the unit cell
      nunitcell=nmol
c      allocate(atsymbolunit(nunitcell),xunit(nunitcell),yunit(nunitcell),zunit(nunitcell))

      write(6,605) ntype
605   format('makelatticeato> Number of atom types found: ',i3)
      do i=1,ntype
         write(6,'(1x,i3,1x,a3)') i,type(i)
      end do
      
      ! Get full name of output file
      
      call get_full_filename(readfile,ilenreadfile,fnameato,ilenfnameato)
      if(index(fnameato(1:ilenfnameato),'.ato').eq.0) then
         fnameato=fnameato(1:ilenfnameato)//'.ato'
         ilenfnameato=ilenfnameato+4
      end if
      write(6,609) fnameato(1:ilenfnameato)
609   format('makelatticeato> Output being written to file: ',a)
      
      !OK Ready to build lattice with specified number of unit cells
      ia=0
      ig=0
      ngrrot=0
!      ngrrotmax=0
      nrottot=0
!      nrottotmax=0
      im=0
      xav=0.0
      yav=0.0
      zav=0.0
      nadd=0
c
c set up position of centre of box (roughly)
c
      acent=0.5*real(na)
      bcent=0.5*real(nb)
      ccent=0.5*real(nc)
      do kc=1,nc
         fc=real(kc-1)-ccent
         zlat=fc*clatz
         do kb=1,nb
            fb=real(kb-1)-bcent
            ylat=fb*blaty+fc*claty
            do ka=1,na
               fa=real(ka-1)-acent
c
c Cartesian coordinates of this lattice point
c
               xlat=fa*alatx+fb*blatx+fc*clatx
c
c now within this unit cell set up the atom positions
c
               ja=0
               irgref=0
               do jm=1,nunitcell
                  im=im+1
                  nato(im)=nato(jm)
                  xadd(im)=xlat+xadds(jm)
                  yadd(im)=ylat+yadds(jm)
                  zadd(im)=zlat+zadds(jm)
                  tethered(im)=tethered(jm)
                  mtype(im)=mtype(jm)
                  if(tethered(im)) then
                     phix(im)=xlat+phixs(jm)
                     phiy(im)=ylat+phiys(jm)
                     phiz(im)=zlat+phizs(jm)
                  else
                     phix(im)=phixs(jm)
                     phiy(im)=phiys(jm)
                     phiz(im)=phizs(jm)
                  endif                  
c
c form mean coordinate of all the molecules
c
                  xav=xav+xadd(im)
                  yav=yav+yadd(im)
                  zav=zav+zadd(im)
c for each atom set its label
                  naf(im)=ia+1
                  do jb=1,nato(jm)
                     ja=ja+1
                     ia=ia+1
                     alabel(ia)=alabel(ja)
      write(6,*) im,jm,ia,ja,imageatom(ja)
                     if(imageatom(ja).gt.0) then
                        imageatom(ia)=imageatom(ja)+nadd
                     else
                        imageatom(ia)=0
                     end if
                     nbond(ia)=nbond(ja)
                     xa(ia)=xa(ja)
                     ya(ia)=ya(ja)
                     za(ia)=za(ja)
                     itype(ia)=itype(ja)
                     mref(ia)=im
                     if(nbond(ia).gt.0) then
                        do ib=1,nbond(ia)
                           iabond(ia,ib)=iabond(ja,ib)+nadd
                           tbond(ia,ib)=tbond(ja,ib)
                        end do
                     endif
                  end do
c number of internal rotation groups
                  ngrot(im)=ngrot(jm)
                  if(ngrot(im).gt.0) then
                     ngrrot=ngrrot+ngrot(im)
                     call testsize_rotgrps()
                     do ih=1,ngrot(im)
                        ig=ig+1
                        irgref=irgref+1
c save molecule in which this rotational group occurs
                        mrot(ig)=im
                        rotgrp(ig)=rotgrp(irgref)
c for rotational groups read atom numbers of two atoms to form axis of rotation
                        do j=1,2
                           labaxi(ig,j)=labaxi(irgref,j)+nadd
                        end do
c read number of atoms in group and their numbers
                        nrot(ig)=nrot(irgref)
                        rotf(ig)=nrottot+1
                        nrottot=nrottot+nrot(ig)
                        rotl(ig)=nrottot
                        call testsize_labrot()
                        labrefold=rotf(irgref)-1
                        labrefnew=rotf(ig)-1
                        do j=1,nrot(ig)
                           labrefold=labrefold+1
                           labrefnew=labrefnew+1
                           labrot(labrefnew)=labrot(labrefold)+nadd
                        end do
! Do any dihedral groups.
                        if(rotgrp(ig).eq.'DIH') then
                           ndihedralang(ig)=ndihedralang(irgref)
                           do k=1,ndihedralang(ig)
                              dihedralang(k,ig)=dihedralang(k,irgref)
                              do j=1,2
                                 labrotdih(j,k,ig)=labrotdih(j,k,irgref)+nadd
                              end do
                           end do
                        end if
                     end do
                  endif
                  nal(im)=ia
               end do
               nadd=ia
            end do
         end do
      end do
      natom=ia
      nmol=im
c
c coordinates of centre of box
c
      xbox=xav/nmol
      ybox=yav/nmol
      zbox=zav/nmol
c
c shift molecular coordinates so that centre of all molecules is at centre of
c box
c
      do im=1,nmol
            xadd(im)=xadd(im)-xbox
            yadd(im)=yadd(im)-ybox
            zadd(im)=zadd(im)-zbox
            if(tethered(im)) then
               phix(im)=phix(im)-xbox
               phiy(im)=phiy(im)-ybox
               phiz(im)=phiz(im)-zbox
            end if
      end do
      natom=ia
      write(6,610) nmol,natom
610   format('makelatticeato> ',i8,' molecules and ',i8,' atoms created')
      write(6,611) xbox,ybox,zbox
611   format('makelatticeato> Centre of box at: ',3(1x,e12.5))
      alat=alat*na
      blat=blat*nb
      clat=clat*nc
c  Retain the geometry if possible
      cubic=cubic.and.(na.eq.nb.and.na.eq.nc.and.nb.eq.nc)
      if(cubic) dmax=dmax*na
      write(6,612) convert
612   format('makelatticeato> Value of convert = ',L)
      if(convert) then
         !First make all the atom coordinates absolute
         do im=1,nmol
            do ia=naf(im),nal(im)
               xa(ia)=xa(ia)+xadd(im)
               ya(ia)=ya(ia)+yadd(im)
               za(ia)=za(ia)+zadd(im)
            end do
         end do
         !Set the number of molecules equal to number of atoms
         nmol=natom
         !Define the coordinates of each molecule and tether it
         do im=1,nmol
            ia=im
            xadd(im)=xa(ia)
            phix(im)=xa(ia)
            yadd(im)=ya(ia)
            phiy(im)=ya(ia)
            zadd(im)=za(ia)
            phiz(im)=za(ia)
            tethered(im)=.true.
            tetheredatom(im)=0
            nstepreject(im)=0
            nato(im)=1
            naf(im)=ia
            nal(im)=ia
            nbond(im)=0
            ngrot(im)=0
            xa(ia)=0.0
            ya(ia)=0.0
            za(ia)=0.0
            imageatom(ia)=0
         end do
!Set the tolerance low so that atoms cannot move far from their starting positions
         tolout=0.1
         nmoltype=0
c
c make a list of the molecule types read in from the current file
c
         nmoltypefirst=nmoltype
         im=0
         do while(im.lt.nmol)
            im=im+1
            imoltype=0
            ifound=0
            do while(imoltype.lt.nmoltype.and.ifound.eq.0)
               imoltype=imoltype+1
               if(alabel(naf(im)).eq.moltype(imoltype)) then
                  ifound=imoltype
                  mtype(im)=ifound
               end if
            end do
            if(ifound.eq.0) then
               nmoltype=nmoltype+1
               moltype(nmoltype)=alabel(naf(im))
               moltethered(nmoltype)=tethered(im)
               moltetheredatom(nmoltype)=tetheredatom(im)
c
c initialise the template filename
c
               write(moltemplate(nmoltype),'("moltype",i2.2)') nmoltype
               mtype(im)=nmoltype

c Initialise the rotational and translational step sizes

               wholerot(nmoltype)=stepr
               wholetrans(nmoltype)=stepm

            endif
         end do
! Initialise the relative broadening of each molecule type if coarse graining
         call reallocate1d_r(rfractmol,nmoltype)
         rfractmol=0.0

      endif
      !Write out ato file
      call moloutc(fnameato)
      
      return
      end

