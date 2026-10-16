      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_writexyz()

      end      

      subroutine run_writexyz()

c To generate a 3D plot of a .ato file

c 2006-11-15 Output options are (1) GNUplot, (2) PGPlot, (3) ATOMS

c Also atoms to be plotted will be defined by their types, not atomic symbol

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'plot_pars.inc'
      include 'system_commands.inc'
      include 'plotato.inc'
      include 'lattice.inc'

      character*256 readfile
      character*4 ext
      character*1 iaccept
      character*3020 sometextsave
      character*4 yesorno

      logical appendset


      integer iread,nwordssave,ncfsave(nworddim),nclsave(nworddim)

c First identify the number of arguments read in

      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)

c If number of words is greater than 2 then we need to store the remaining arguments in a text stream

      appendxyz=.true.
      appendset=.false.

c Expect to see the append type as word 3. If not found assume xyz file is to be appended

      if(nwords.gt.2) then

         appendxyz=.not.sometext(ncf(3):ncl(3)).eq.'noappend'
         appendset=.true.

      endif

      if(nwords.gt.3) then

         sometextsave=sometext(ncf(4):ncl(nwords))
         call parse(sometextsave,nwordssave,ncfsave,nclsave,5,nstsize,nworddim)
         write(6,*) nwordssave,' ',sometextsave(ncfsave(1):nclsave(nwordssave))
         iread=1

      else

         sometextsave=''
         nwordssave=0
         write(6,*) nwordssave
         iread=1

      endif

c generate the name of the .ato file

      ext='.ato'
      ilenext=4
      if(nwords.lt.2) then
         nwords=0
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
499   format(/'writexyz> No .ato file specified! - routine will exit'/)
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
599   format(/'writexyz> Problems with .ato file ',a/'Routine will exit now'/)
         return
      endif
c     call normalise_coordinates
c
c step through the atom types and make sure the chemical symbol has been specified
c
      do it=1,ntype
         if(atsymbol(it).eq.'XXX'.or.atsymbol(it).eq.'   ') then
            write(6,115) type(it)
115   format(/'writexyz> Specify the chemical symbol for atom type ',a3,': ',$)
            read(5,*) atsymbol(it)
         endif
      end do

      if(.not.appendset) then

         write(6,121,advance='no')
121   format(/'writexyz> Do you want to append to an existing file? (y or n) ')
         read(5,122) yesorno
122   format(a)
         appendxyz=.not.(yesorno(1:1).eq.'n'.or.yesorno(1:1).eq.'N')
         appendset=.true.

      endif
c
c sort out which atoms will actually be plotted, based on whether all the atoms 
c are to be plotted, a group of them centred on one molecule, or just one molecule
c
      if(nwordssave.lt.iread) then
         write(6,132,advance='no')
132   format(/'plotato> Specify number of molecule at centre (0 to plot entire box): ')
         read(5,*) nmolcentre
      else
         read(sometextsave(ncfsave(iread):nclsave(iread)),*) nmolcentre
      endif
      iread=iread+1

      if(nmolcentre.gt.0) then
         plottype=2
      else
         plottype=1
      endif

c
c sort out which atoms will actually be plotted, based on whether all the atoms 
c are to be plotted, a group of them centred on one molecule, or just one molecule
c
c      if(nmolcentre.gt.0) then
         if(nwordssave.lt.iread+3) then
            write(6,133,advance='no')
133   format(/'plotato> Give maximum distances along x and y axes,'
     */,'minimum and maximum distances along the z-axis'
     */' - use zero distances to show only the central molecule: ')
            read(5,*) xaxismax,yaxismax,zaxismin,zaxismax
         else
            read(sometextsave(ncfsave(iread):nclsave(iread+3)),*) xaxismax,yaxismax,zaxismin,zaxismax
         endif
         iread=iread+4

         if(nwordssave.lt.iread+2) then
            write(6,134,advance='no')
134   format(/'plotato> Give phi, theta and chi for rotation of coordinates (0 0 0) ')
            read(5,*) zphi,ztheta,zchi
         else
            read(sometextsave(ncfsave(iread):nclsave(iread+2)),*) zphi,ztheta,zchi
         endif
         iread=iread+3

c      else

c         xaxismax=dmax
c         yaxismax=dmax
c         zaxismin=-dmax
c         zaxismax=dmax
c         zphi=0.0
c         ztheta=0.0
c         zchi=0.0

c      endif

c Get sines and cosines

      czphi=cosd(zphi)
      szphi=sind(-zphi)
      cztheta=cosd(ztheta)
      sztheta=sind(-ztheta)
      czchi=cosd(zchi)
      szchi=sind(-zchi)
      

c Ensure the coordinates of each atom in each molecule are relative to the centre of mass

      do im=1,nmol

         call centreofmass(im)

      end do
c
c step through molecules and determine which ones are to be plotted
c
      iaplot=0
      implot=0

c Ensure coordinates plotted are relative to molecule at centre

      if(nmolcentre.gt.0.and.nmolcentre.le.nmol) then
         xaddcentre=xadd(nmolcentre)
         yaddcentre=yadd(nmolcentre)
         zaddcentre=zadd(nmolcentre)
      else
         xaddcentre=0.0
         yaddcentre=0.0
         zaddcentre=0.0
      endif
      do im=1,nmol

c Apply minimum image to centre of mass coordinates.

         xdif=xadd(im)-xaddcentre
         ydif=yadd(im)-yaddcentre
         zdif=zadd(im)-zaddcentre
         if(.not.cubic) then

            call minimagec(xdif,ydif,zdif)

         else

            xdif=xdif/dmax
            ydif=ydif/dmax
            zdif=zdif/dmax
            call minimage(xdif,ydif,zdif)
            xdif=xdif*dmax
            ydif=ydif*dmax
            zdif=zdif*dmax
      
         endif
c
c Rotate these coordinates to the new axis
c
         call rotatecoords(czphi,szphi,cztheta,sztheta,czchi,szchi,xdif,ydif,zdif)
c
c if one molecule is to be at centre, check to see if the current molecule 
c is within the correct distance
c
         ifound=0
c         if(plottype.gt.1) then
c            dtest=sqrt(xdif*xdif+ydif*ydif+zdif*zdif)
c            if(dtest.le.radmax.or.im.eq.nmolcentre) then
            if((abs(xdif).le.xaxismax.and.abs(ydif).le.yaxismax
     *.and.zdif.ge.zaxismin.and.zdif.le.zaxismax).or.im.eq.nmolcentre) then
               ifound=im
            endif
c         else
c            ifound=im
c         endif          
         if(ifound.gt.0) then
            implot=implot+1
            do ia=naf(im),nal(im)
               iaplot=iaplot+1
               iplot(ia)=1
            end do
         else
            do ia=naf(im),nal(im)
               iplot(ia)=0
            end do
         endif

c Rotate each molecule onto the new coordinate system

         do ia=naf(im),nal(im)

            xat(ia)=xa(ia)
            yat(ia)=ya(ia)
            zat(ia)=za(ia)
            call rotatecoords(czphi,szphi,cztheta,sztheta,czchi,szchi,xat(ia),yat(ia),zat(ia))

         end do

c
c generate the coordinates of each atom relative to the centre of the plot
c
         do ia=naf(im),nal(im)

            xat(ia)=xat(ia)+xdif
            yat(ia)=yat(ia)+ydif
            zat(ia)=zat(ia)+zdif

         end do
      end do
      write(6,*) implot,iaplot

c Make a list of atom types for those atoms that could be plotted.

      ntypplot=0
      do i=1,natom
         it=itype(i)
         if(iplot(i).gt.0) then
            iflagt=0
            j=0
            do while(j.lt.ntypplot.and.iflagt.eq.0)
               j=j+1
               if(itypplot(j).eq.it) iflagt=j
            end do
            if(iflagt.eq.0.and.ntypplot.lt.mcomp) then
               ntypplot=ntypplot+1
               itypplot(ntypplot)=it
            else
            endif
         endif
      end do

      if(nwordssave.lt.iread) then
         write(6,111) ntypplot,(itypplot(j),type(itypplot(j)),j=1,ntypplot)
111   format(/'writexyz> ',i3,' components available for plotting:-'
     1,/(1x,i3,1x,a3))
         write(6,131,advance='no')
131   format(/'writexyz> Give component numbers to EXCLUDE from the plot - '
     +,' type 0 to include all: ')
         read(5,101) sometext
101   format(a)
      else
         read(sometextsave(ncfsave(iread):nclsave(nwordssave)),101) sometext
      endif
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
c      write(6,*) iread,nwordssave,sometext(ncf(1):ncl(nwords))


c Assume at the outset all atom types are to be plotted

      nclplt=ntypplot
      do i=1,ntypplot

         nplt(i)=itypplot(i)

      end do

      iw=0
      ierr=0
      do while(iw.lt.nwords.and.ierr.eq.0)

         iw=iw+1
         read(sometext(ncf(iw):ncl(iw)),*,iostat=ierr) itest
         if(ierr.eq.0) then
            i=0
            do while(i.lt.nclplt)
               i=i+1
               if(nplt(i).eq.itest) then
                  j=i-1
                  do while (j.lt.nclplt-1)
                     j=j+1
                     nplt(j)=nplt(j+1)
                  end do
                  nclplt=j
               endif
            end do
         endif
      end do
      write(6,151) nclplt,(nplt(j),type(nplt(j)),j=1,nclplt)
151   format(/'writexyz> ',i3,' components to be printed:-'
     1,/(1x,i3,1x,a3))

c Based on this list, revise the values of iplot so that only the specified components
c are plotted.

      do ia=1,natom
         it=itype(ia)
         if(iplot(ia).gt.0) then
            iflagt=0
            j=0
            do while(j.lt.nclplt.and.iflagt.eq.0)
               j=j+1
               if(nplt(j).eq.it) iflagt=j
            end do

c If the type of this atom does not correspond to one of the specified types, we signal
c we do not want to plot this atom

            if(iflagt.eq.0.and.mref(ia).ne.nmolcentre) then
               iplot(ia)=0
            endif
         else
            iplot(ia)=0
         endif

      end do

      call writexyz(readfile,ilenreadfile)

      return
      end


      subroutine writexyz(readfile,ilenreadfile)

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'plot_pars.inc'
      include 'system_commands.inc'
      include 'plotato.inc'

      character*256 readfile
      character*5 blocknumber
      integer ilenreadfile,nblocks,ir
      logical centralmolecule

      ir=0
      if(appendxyz) then

c       
c generate name of appended xyz file.
c
         fnamesplot(0)=directory_name(1:nchardir)//readfile(1:ilenreadfile)//'append.xyz'
         ilenfnamesplot(0)=nchardir+ilenreadfile+4

c Write xyz output file

         open(10,file=fnamesplot(0),status='unknown',iostat=ierr)
         if(ierr.eq.0) then

c Read to the end of file

            nblocks=0
            do while(ierr.eq.0)
               read(10,*,iostat=ierr) natplot
c      write(6,*) natplot,nblocks,ierr
               if(ierr.eq.0) ir=0
               do while(ierr.eq.0.and.ir.lt.natplot+2) 
                  ir=ir+1
                  read(10,*,iostat=ierr) 
c      write(6,*) ir,ierr
               end do
               if(ierr.eq.0) nblocks=nblocks+1
            end do

c Backspace the file to back before the end of file

            backspace(10)

c If an error occurred before the correct number of records were read in, rewind the file to the beginning and start again

c            write(6,*) ir,natplot+2,nblocks
            if(ir.ne.natplot+2) then
               rewind(10)
               nblocks=0
            endif
            nblocks=nblocks+1

c First determine how many atoms to plot
         
            natplot=0
            do ia=1,natom
               if(iplot(ia).gt.0) natplot=natplot+1
            end do
            write(10,*,iostat=ierr) natplot
c            write(6,*) natplot,nblocks,ierr
            write(blocknumber,'(i5.5)') nblocks

c Output some information text

            write(10,100) '#Filename: '//directory_name(1:nchardir)//readfile(1:ilenreadfile)//'.ato Block '//blocknumber
100   format(a)
            do ia=1,natom
               j=0
               iflagt=0
               !Check if a central molecule has been defined.
               centralmolecule=nmolcentre.gt.0.and.nmolcentre.le.nmol
               !Check if the current atom is in the central molecule
               if(centralmolecule) centralmolecule=ia.ge.naf(nmolcentre).and.ia.le.nal(nmolcentre)
               do while(j.lt.nclplt.and.iflagt.eq.0)
                  j=j+1
                  if(iplot(ia).gt.0.or.centralmolecule) then
                     iflagt=j
                     ntag=ntag+1
                     write(10,113) atsymbol(itype(ia)),xat(ia),yat(ia),zat(ia),alabel(ia)
                  endif
113   format(a3,4x,e10.3,2x,e10.3,2x,e10.3,2x,a3)
               end do
            end do

c Write a blank line at end of file

            write(10,*)

         endif
         close(10)

      else

c If separate files are being stored, then we need to get the last block number

         fnamesplot(0)=directory_name(1:nchardir)//readfile(1:ilenreadfile)//'.nblock'

c Write xyz output file

         open(10,file=fnamesplot(0),status='unknown',iostat=ierr)

         nblocks=0
         if(ierr.eq.0) read(10,*,iostat=ierr) ir
         if(ierr.eq.0) nblocks=ir
c      write(6,*) ir,nblocks,ierr
         rewind(10)

         nblocks=nblocks+1
         write(10,*) nblocks 

         close(10)

c Generate the name of the file to write the data to

         write(blocknumber,'(i5.5)') nblocks
         fnamesplot(0)=directory_name(1:nchardir)//readfile(1:ilenreadfile)//blocknumber//'.xyz'
         ilenfnamesplot(0)=nchardir+ilenreadfile+5+4

c Write xyz output file

         open(10,file=fnamesplot(0),status='unknown',iostat=ierr)
      
         if(ierr.eq.0) then

c First determine how many atoms to plot
         
            natplot=0
            do ia=1,natom
               if(iplot(ia).gt.0) natplot=natplot+1
            end do
            write(10,*,iostat=ierr) natplot

c Output some information text

            write(10,100) '#Filename: '//directory_name(1:nchardir)//readfile(1:ilenreadfile)//'.ato  Block '//blocknumber
            do ia=1,natom
               !Check if a central molecule has been defined.
               centralmolecule=nmolcentre.gt.0.and.nmolcentre.le.nmol
               !Check if the current atom is in the central molecule
               if(centralmolecule) centralmolecule=ia.ge.naf(nmolcentre).and.ia.le.nal(nmolcentre)
               j=0
               iflagt=0
               do while(j.lt.nclplt.and.iflagt.eq.0)
                  j=j+1
                  if(iplot(ia).gt.0.or.centralmolecule) then
                     iflagt=j
                     ntag=ntag+1
                     write(10,113) atsymbol(itype(ia)),xat(ia),yat(ia),zat(ia),alabel(ia)
                  endif
               end do
            end do

c Write a blank line at end of file

            write(10,*)

         endif
         close(10)

      endif

      return
      end

