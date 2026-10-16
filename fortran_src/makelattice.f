      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call makelattice()

      end      

      subroutine makelattice()
c
c Makes an EPSR18 ato file from a set of supplied lattice coordinates. 
c The unit cell coordinates are supplied in a file alongside
c other potential parameters - this filename will be given to the created .ato file.
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'

      character(len=256)                           :: readfile
      integer                                      :: ilenreadfile,ierr,i,iexch,nwordstest
      real                                         :: pi,radconv,alpha,beta,gamma,alphar,betar
     1,sinphib,cosphib,cosalpha,cosbeta,top,bottom,sum
      real, dimension(:), allocatable              :: xunit,yunit,zunit !x,y,z coordinates in unit cell
      integer                                      :: nunitcell,na,nb,nc,ka,kb,kc !Number of atoms in unit cell and number of steps along
      !a, b and c axes
      character(len=3), dimension(:), allocatable  :: atsymbolunit !Atom symbols in unit cell
      character(len=3)                             :: a1,ats
      real                                         :: epsit,sigt,amass,qt,rqt,acent,bcent,ccent,fa,fb,fc
     1,xav,yav,zav,xbox,ybox,zbox
      
      pi=4.0*atan(1.0)
      radconv=pi/180.0
      degconv=180.0/pi
c
c generate the name of the file containing the unit cell coordinates
c
      nwordstest=2
      if(nwords.lt.2) then
         nwords=0
         nwordstest=1
         do while (nwords.le.0)
            write(6,100,advance='no')
100   format(/'makelattice> Type name of file containing unit cell coordinates (incl. extension): ')
            read(5,101) sometext
101   format(a)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
499   format(/'makelatice> No file specified! - try again!'/)
            endif
         end do
         readfile=sometext(ncf(1):ncl(1))
         ilenreadfile=ncl(nwords)-ncf(1)+1
      else
         readfile=sometext(ncf(2):ncl(2))
         ilenreadfile=ncl(2)-ncf(2)+1
      endif
      write(6,*) nwords,nwordstest
      if(nwords.gt.nwordstest) then
         read(sometext(ncf(nwordstest+1):ncl(nwordstest+1)),*) na
         if(nwords.gt.nwordstest+1) then
            read(sometext(ncf(nwordstest+2):ncl(nwordstest+2)),*) nb
            if(nwords.gt.nwordstest+2) then
               read(sometext(ncf(nwordstest+3):ncl(nwordstest+3)),*) nc
            else
               write(6,606,advance='no')
606   format('makelattice> Give number of unit cells along c-axis: ')
               read(5,*) nc
            end if
         else
            write(6,607,advance='no')
607   format('makelattice> Give number of unit cells along b and c-axes: ')
            read(5,*) nb,nc
         end if
      else
         write(6,608,advance='no')
608   format('makelattice> Give number of unit cells along a, b and c-axes: ')
         read(5,*) na,nb,nc
      end if
c
c Read this file
c
      fnameato=directory_name(1:nchardir)//readfile(1:ilenreadfile)
      ilenfnameato=nchardir+ilenreadfile
      open(10,file=fnameato,status='old',iostat=ierr)
      if(ierr.eq.0) then
         !The lattice a, b and c axes
         read(10,*) alat,blat,clat
         !The lattice internal angles
         read(10,*) alpha,beta,gam
         !Convert these to the internal EPSR18 angles phib,thetac,phic
         call convertalphabetagamma(alpha,beta,gam,phib,thetac,phic)
      else
         write(6,602) fnameato(1:ilenfnameato)
602   format('makelattice> Problems with file: ',a)
         return
      endif
      write(6,603) alat,blat,clat,alpha,beta,gam,phib,thetac,phic
603   format('makelattice> Lattice constants, angles and derived angles:'
     13(/3(1x,f10.5)))
      !If the word 'absolute' is encountered here then absolute Cartesian coordinates are expected
      !If the word 'fractional' is encountered here then fractional unit cell coordinates are expected
      !Default is fractional.
      absolute=.false.
      read(10,101) sometext
      if(index(sometext,'absolute').gt.0) then
         absolute=.true.
         !Get the number of atoms in the unit cell
         read(10,*) nunitcell
      else if(index(sometext,'fractional').gt.0) then
         absolute=.false.
         !Get the number of atoms in the unit cell
         read(10,*) nunitcell
      else
         read(sometext,*) nunitcell
      end if
      allocate(atsymbolunit(nunitcell),xunit(nunitcell),yunit(nunitcell),zunit(nunitcell))
      !Get the unit cell coordinates
      ntype=0
      do i=1,nunitcell
         read(10,101) sometext
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         if(nwords.gt.3) then
            atsymbolunit(i)=sometext(ncf(1):ncl(1))
            !Store this temporarily to set up the list of atom types so that the list of potential
            !parameters can be checked.
            alabel(i)=atsymbolunit(i)
            call ftype(i)
            read(sometext(ncf(2):ncl(4)),*) xunit(i),yunit(i),zunit(i)
         else
            write(6,604)
604   format('makelattice> Insufficient data supplied in input file.')
            return
         endif
      end do
      write(6,605) ntype
605   format('makelattice> Number of atom types found: ',i3)
      do i=1,ntype
         write(6,'(1x,i3,1x,a3)') i,type(i)
      end do
c
c read Lennard-Jones, mass and charge parameters for each atom type
c
      do it=1,ntype
            read(10,101,iostat=ierr) sometext
c
c parse this text into words
c
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            a1=sometext(ncf(1):ncl(1))
            if(nwords.gt.1) then
               ats=sometext(ncf(2):ncl(2))
            else
               ats='XXX'
            endif
            if(nwords.gt.2) then
               read(sometext(ncf(3):ncl(3)),*) iexch
            else
               iexch=-1
            endif
            read(10,*,iostat=ierr) epsit,sigt,amass,qt,rqt
            natomsph=1.0
            radminsph=0.0
            radmaxsph=0.5*sigt
            pseudofacsph=1.0
            rhosph=0.0
c assign type number
            call ctype(a1,ats,iexch,epsit,sigt,amass,qt,rqt,natomsph,rhosph,radminsph,radmaxsph,pseudofacsph)
      end do
c
c read core energy and diameter - used in fmole for non-bonded atoms
c
      read(10,*,iostat=ierr) ecore,dcore
      if(ierr.ne.0) then
         ecore=0.0
         dcore=0.0
      end if
c
c initialise random number generator
c
      call initran1(10)
      close(10)
c
c get the filename to put the mixture in
c
      nwords=0
      do while (nwords.le.0)
         write(6,35)
35    format(/'makelattice> Type name of file to put lattice in: ',$)
         read(5,101,iostat=ierr) sometext
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
      ilenuptodot=index(readfile,'.')-1
      if(ilenuptodot.gt.0) then
         readfile=readfile(1:ilenuptodot)
         ilenreadfile=ilenuptodot
      endif
      fnameato=directory_name(1:nchardir)//readfile(1:ilenreadfile)//'.ato'
      ilenfnameato=nchardir+ilenreadfile+4
      write(6,609) fnameato(1:ilenfnameato)
609   format('makelattice> Output being written to file: ',a)
      
      !OK Ready to build lattice with specified number of unit cells
c
c define the conversion from lattice vectors to Cartesian coordinates
c
      alatx=alat
      blatx=blat*cosd(phib)
      blaty=blat*sind(phib)
      clatx=clat*sind(thetac)*cosd(phic)
      claty=clat*sind(thetac)*sind(phic)
      clatz=clat*cosd(thetac)
      nmol=0
      simtk=300.0
      tolout=0.1
      stepmi=0.01
      stepri=0.0
      stepr=0.001
      stepm=0.001
      tempb=65.0
      tempangreq=1.0
      tempdihreq=10.0
      tempangwt=1.0
      tempdihwt=10.0
      ia=0
      ig=0
      im=0
      xav=0.0
      yav=0.0
      zav=0.0
      natom=0
      ngrrot=0
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
               do jm=1,nunitcell
                  im=im+1
                  nato(im)=1
                  if(absolute) then
                     !Use absolute Cartesian coordinates
                     xadd(im)=xlat+xunit(jm)
                     yadd(im)=ylat+yunit(jm)
                     zadd(im)=zlat+zunit(jm)
                  else
                     !Convert fractional coordinates to Cartesian values
                     xadd(im)=xlat+xunit(jm)*alatx+yunit(jm)*blatx+zunit(jm)*clatx
                     yadd(im)=ylat+yunit(jm)*blaty+zunit(jm)*claty
                     zadd(im)=zlat+zunit(jm)*clatz
                  end if
c
c form mean coordinate of all the molecules
c
                  xav=xav+xadd(im)
                  yav=yav+yadd(im)
                  zav=zav+zadd(im)
                  phix(im)=0.0
                  phiy(im)=0.0
                  phiz(im)=0.0
                  tethered(im)=.true.
c for each atom set its label
                  naf(im)=ia+1
                  do ib=1,nato(im)
                     ia=ia+1
                     alabel(ia)=atsymbolunit(jm)
                     itype(ia)=itype(jm)
                     nbond(ia)=0
                     xa(ia)=0.0
                     ya(ia)=0.0
                     za(ia)=0.0
                  end do
c number of internal rotation groups
                  ngrot(im)=0
                  nal(im)=ia
               end do
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
            phix(im)=xadd(im)
            phiy(im)=yadd(im)
            phiz(im)=zadd(im)
      end do
      write(6,610) nmol,natom
610   format('makelattice> ',1x,i5,' molecules and ',i5,' atoms created')
      write(6,611) xbox,ybox,zbox
611   format('makelattice> Centre of box at: ',3(1x,e12.5))
c
c make a list of the molecule types read in from the current file
c
      nmoltype=0
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

c Write these out as a check

      write(6,201)
201   format(/'makelattice> Following molecule types defined initially:-'/)
      do imoltype=1,nmoltype
         write(6,202) imoltype,moltype(imoltype),moltemplate(imoltype)(1:10),wholerot(imoltype)
     *,wholetrans(imoltype)
202   format(1x,i3,1x,a3,1x,a10,2(1x,e12.5))
      end do
      alat=alat*na
      blat=blat*nb
      clat=clat*nc
      !Write out ato file
      call moloutc(fnameato)
      
      return
      end

