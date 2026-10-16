      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_plotato()

      stop
      end      

      subroutine run_plotato()

c To generate a 3D plot of a .ato file

c 2006-11-15 Output options are (1) GNUplot, (2) PGPlot, (3) ATOMS

c 2011-08-09 Modified to allow a slice of the ato file in a specified plane.

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

      real thestereo(2),phistereo(2)

      integer iread,nwordssave,ncfsave(nworddim),nclsave(nworddim)

      if(nwords.gt.2) then

         sometextsave=sometext(ncf(3):ncl(nwords))
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
      fnameato=fnameato(1:ilenfnameato)//ext
      ilenfnameato=ilenfnameato+ilenext
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
      if(ierr.ne.0) then
         write(6,599) fnameato
599   format(/'Problems with .ato file ',a/'Routine will exit now'/)
         return
      endif
c     call normalise_coordinates
c
c step through the atom types and make sure the chemical symbol has been specified
c
      do it=1,ntype
         if(atsymbol(it).eq.'XXX'.or.atsymbol(it).eq.'   ') then
            write(6,115) type(it)
115   format(/'plotato> Specify the chemical symbol for atom type ',a3,': ',$)
            read(5,*) atsymbol(it)
         endif
      end do

c Decide what kind of plot you want

      if(nwordssave.lt.iread) then
         write(6,450,advance='no')
 450  format(/'plotato> Decide what kind of output you want:-'
     *//' 1 = GNUplot'
     */' 2 = PGplot'
     */' 3 = Jmol'
     *//'plotato> ? ')
         read(5,*) iplottype
      else
         read(sometextsave(ncfsave(iread):nclsave(iread)),*) iplottype
      endif
      iread=iread+1
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

      if(nmolcentre.gt.0) then
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

      else

         xaxismax=dmax
         yaxismax=dmax
         zaxismin=-dmax
         zaxismax=dmax
         zphi=0.0
         ztheta=0.0
         zchi=0.0

      endif

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
     *.and.zdif.ge.zaxismin.and.zdif.le.zaxismax).or.im.eq.nmolcentre.or.nmolcentre.eq.0) then
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
111   format(/'plotato> ',i3,' components available for plotting:-'
     1,/(1x,i3,1x,a3))
         write(6,131,advance='no')
131   format(/'plotato> Give component numbers to EXCLUDE from the plot - '
     +,' type 0 to include all: ')
         read(5,101) sometext
101   format(a)
      else
         read(sometextsave(ncfsave(iread):nclsave(nwordssave)),101) sometext
      endif
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)

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
151   format(/'plotato> ',i3,' components to be plotted:-'
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
c we do not want to plot this atom. However if the atom is one those for the central molecule
c it is always plotted

c            if(iflagt.eq.0.and.mref(ia).ne.nmolcentre) then
            if(iflagt.eq.0) then
               iplot(ia)=0
            endif
         else
            iplot(ia)=0
         endif

      end do

      if(iplottype.ne.3) then

c For each class set up radii and colours for atoms

c First read any preexisting definitions

      call gnuatomsread(unit,posfrac,ptscaler
     *,nclassread,classlab,classtype,classrad
     *,classcolour,classcoloura,classrgb,mclass)

      write(6,127) nclassread
127   format(/'plotato> Number of classes read from gnuatoms.txt = ',i5)

c Setup the maximum and minimum x, y and z values for plot

      xmin=1.0e30
      xmax=-1.0e30
      ymin=1.0e30
      ymax=-1.0e30
      zmin=1.0e30
      zmax=-1.0e30

c Now step through the atoms and set the classes

      nastr=0
      nclass=0

      do i=1,natom

c Signal for which kind of atom to plot. 0 means do not plot

         if(iplot(i).gt.0) then

c Count the number of atoms to be plotted

            nastr=nastr+1

c       So now we need to find an appropriate class from the list read in

            j=0
            ifound=0
            do while (j.lt.nclassread.and.ifound.eq.0)

               j=j+1
               if(classlab(j).eq.atsymbol(itype(i))) ifound=j

            end do
            if(ifound.eq.0.and.nclassread.lt.mclass) then

c       Need to read parameters for this class

               nclassread=nclassread+1
               ifound=nclassread

               write(6,346) atsymbol(itype(i))
 346  format(/'plotato> Give radius and RGB values for class ',a3,': ',$)
               classlab(nclassread)=atsymbol(itype(i))
               read(5,*) classrad(nclassread),(classrgb(j,nclassread),j=1,3)

c       GNUplot character type and line type

               write(6,345) atsymbol(itype(i))
 345  format(/'plotato> Give GNUplot character number and line type for class '
     *,a3,': ',$)
               read(5,*) classtype(nclassread),classcolour(nclassread)

c       Default ATOMS colour

               classcoloura(nclassread)=15

            endif

c Set up the list of classes to be plotted

            iplot(i)=ifound

c Need to scan through the existing class reference integers, if any, to make sure
c this class has not already been listed

            j=0
            iflagt=0
            do while(j.lt.nclass.and.iflagt.eq.0)

               j=j+1
               if(iplotclass(j).eq.ifound) iflagt=j

            end do
            if(iflagt.eq.0.and.nclass.lt.mcomp) then

               nclass=nclass+1
               iplotclass(nclass)=ifound

            endif

c       Check the maximum and minimum values of x, y, and z

            if(xat(i).lt.xmin) xmin=xat(i)
            if(xat(i).gt.xmax) xmax=xat(i)
            if(yat(i).lt.ymin) ymin=yat(i)
            if(yat(i).gt.ymax) ymax=yat(i)
            if(zat(i).lt.zmin) zmin=zat(i)
            if(zat(i).gt.zmax) zmax=zat(i)

         endif
      end do

      write(6,110) (j,classlab(iplotclass(j)),j=1,nclass)
110   format(/'plotato> Following atom classes will be plotted:-'
     1,/(1x,i3,1x,a3))



c Write the list of atom classes for future reference in case it has changed

      call gnuatomswrite(unit,posfrac,ptscaler
     *,nclassread,classlab,classtype,classrad
     *,classcolour,classcoloura,classrgb)

c Max and min of plot. We need to plot in a cube of dimension big enough to fit in all the atoms

      almin=min(xmin,ymin,zmin)
      almax=max(xmax,ymax,zmax)


c       Set these to be integer multiples of unit A

      lmin=int(almin/unit)-1
      lmax=int(almax/unit)+1
      almin=unit*lmin
      almax=unit*lmax

      write(6,126) almin,almax
 126  format(/'plotato> Min and max of plot = ',2(1x,e13.6))

c     write(6,110) ntype,(j,type(j),j=1,ntype)

c Now read any preexisting bond definitions

      call gnubondsread(nclassbond,classlab1,classlab2
     *,classbond1,classbond2,bondrad,bondrgb,boxrgb,boxlwidth
     *,mpairs)

      ib=0

      iread=1
      write(6,121)
121   format(/'plotato> For each new bond: type two atom symbols,'
     */' minimum and maximum lengths of bond, and radius of bond.'
     */' Type 0 0 0 0 0 to end bond input or ste 0 0 0 0 for stereo pairs.')
122   format(/'plotato> ? ',$)

c Signal no stereo plot initially

      nstereo=1

      do while(iread.gt.0.and.ib.lt.mpairs)

         ib1=ib+1
         write(6,122)
         read(5,*) cbstr1(ib1),cbstr2(ib1),dstr1(ib1),dstr2(ib1),radstr(ib1)
         if(cbstr1(ib1).ne.'0'.and.cbstr1(ib1).ne.'ste') then
            ib=ib1
         else if(iplottype.eq.1) then

            if(cbstr1(ib1).eq.'ste') then

c       Signal stereo plot requested for splot.

             nstereo=2
             write(6,325)
 325         format(/'plotato> Type theta and phi values for each stereo plot: ',$)
             read(5,*) thestereo(1),phistereo(1)
     * ,thestereo(2),phistereo(2)

             do j=1,2

c       Ensure thestereo is between 0 and 180

                do while (thestereo(j).gt.180.0)
                   thestereo(j)=thestereo(j)-180.0
                end do
                do while (thestereo(j).lt.0.0)
                   thestereo(j)=thestereo(j)+180.0
                end do  

c       Ensure phistereo is between 0 and 360

                do while (phistereo(j).gt.360.0)
                   phistereo(j)=phistereo(j)-360.0
                end do
                do while (phistereo(j).lt.0.0)
                   phistereo(j)=phistereo(j)+360.0
                end do
             end do

            else

               nstereo=1

            endif 

            iread=0

         else

            iread=0

         endif
      end do

      nbstr=ib

c Check the existing bond list whether any of these already exists. If so then replace existing
c bond with new values. If not then add it to the list.

      ib=0
      do while (ib.lt.nbstr)

         ib=ib+1
         ic=0
         ifound=0
         do while (ic.lt.nclassbond.and.ifound.eq.0)

            ic=ic+1
            if((classlab1(ic).eq.cbstr1(ib).and.classlab2(ic).eq.cbstr2(ib))
     *.or.(classlab1(ic).eq.cbstr2(ib).and.classlab2(ic).eq.cbstr1(ib))) then

               classbond1(ic)=dstr1(ib)
               classbond2(ic)=dstr2(ib)
               bondrad(ic)=radstr(ib)
               ifound=ic
            endif
         end do
         if(ifound.eq.0.and.nclassbond.lt.mpairs) then

            nclassbond=nclassbond+1
            classlab1(nclassbond)=cbstr1(ib)
            classlab2(nclassbond)=cbstr2(ib)
            classbond1(nclassbond)=dstr1(ib)
            classbond2(nclassbond)=dstr2(ib)
            bondrad(nclassbond)=radstr(ib)

         endif
      end do

c Write out the new list of bonds.

      call gnubondswrite(nclassbond,classlab1,classlab2
     *,classbond1,classbond2,bondrad,bondrgb,boxrgb,boxlwidth)

      write(6,124) nclassbond
124   format(/'plotato> Total number of bonds read in = ',i3)

c Define the overall scale factor for atoms

      write(6,349)
 349  format(/'plotato> Give overall scale factor on atom sizes: ',$)
      read(5,*) scaleatoms

c End ignore if iplottype = 3

      endif

c Set viewing angles of plot

      if(iplottype.ne.3.and.nstereo.eq.1) then
         write(6,201)
 201     format(/'plotato> Give plot rotation about x and y (deg): ',$)
         read(5,*) xangle,yangle
         thestereo(1)=xangle
         phistereo(1)=yangle
      endif

c write out the .ato file in case anything has changed
c
c     call unnormalise_coordinates
      call moloutc(fnameato)

      if(iplottype.eq.1) then

         call run_GNUplotato(readfile,ilenreadfile,nstereo,thestereo,phistereo
     *,almin,almax)

        else if(iplottype.eq.2) then

         call run_pgplotato(readfile,ilenreadfile)

      else if(iplottype.eq.3) then

         call run_Jmol(readfile,ilenreadfile)

      endif

      return
      end

      subroutine gnuatomsread(unit,posfrac,ptscaler
     *,nclassread,classlab,classtype,classrad
     *,classcolour,classcoloura,classrgb,mclass)

      character*3 classlab(*)
      real classrad(*),classrgb(3,*)
      integer classtype(*),classcolour(*),classcoloura(*)

c Default settings for the unit scaler, xy plane positioner, and point size factor

      unit=1.0
      posfrac=0.05
      ptscaler=45.0

      open(20,file='gnuatoms.txt',status='old',iostat=ierr)

      read(20,*,iostat=ierr) unitr,posfracr,ptscalerr
      if(ierr.eq.0) then
         unit=unitr
         posfrac=posfracr
         ptscaler=ptscalerr
      endif

      ic=0

      do while (ic.lt.mclass.and.ierr.eq.0)

         ic1=ic+1
         read(20,*,iostat=ierr) classlab(ic1),classtype(ic1)
     *,classrad(ic1),classcolour(ic1),classcoloura(ic1)
     *,(classrgb(j,ic1),j=1,3)

         if(ierr.eq.0) ic=ic1

      end do
      close(20)

      nclassread=ic

      return
      end

      subroutine gnuatomswrite(unit,posfrac,ptscaler
     *,nclassread,classlab,classtype,classrad
     *,classcolour,classcoloura,classrgb)

      character*3 classlab(*)
      real classrad(*),classrgb(3,*)
      integer classtype(*),classcolour(*),classcoloura(*)

      open(20,file='gnuatoms.txt',status='unknown',iostat=ierr)

      write(20,*) unit,posfrac,ptscaler

      ic=0
      do while (ic.lt.nclassread.and.ierr.eq.0)

         ic=ic+1
         write(20,347,iostat=ierr) classlab(ic),classtype(ic),classrad(ic)
     *,classcolour(ic),classcoloura(ic),(classrgb(j,ic),j=1,3)

347   format(1x,a3,1x,i5,1x,f10.5,2(1x,i5),3(1x,f6.3))

      end do
      close(20)

      return
      end

      subroutine gnubondsread(nclassbond,classlab1,classlab2
     *,classbond1,classbond2,bondrad,bondrgb,boxrgb,boxlwidth
     *,mpairs)

      character*3 classlab1(*),classlab2(*)
      real classbond1(*),classbond2(*),bondrad(*)
     *,bondrgb(*),boxrgb(*),boxlwidth

c Initial values of bond and box colours and box line width

      bondrgb(1)=0.2
      bondrgb(2)=0.2
      bondrgb(3)=0.2
      boxrgb(1)=0.1
      boxrgb(2)=0.1
      boxrgb(3)=0.1
      boxlwidth=0.05

      open(20,file='gnubonds.txt',status='old',iostat=ierr)

      read(20,*,iostat=ierr) a,b,c,d,e,f,g

      if(ierr.eq.0) then

         bondrgb(1)=a
         bondrgb(2)=b
         bondrgb(3)=c
         boxrgb(1)=d
         boxrgb(2)=e
         boxrgb(3)=f
         boxlwidth=g

      endif

      ic=0
      do while (ic.lt.mpairs.and.ierr.eq.0)

         ic1=ic+1
         read(20,*,iostat=ierr) classlab1(ic1),classlab2(ic1)
     *,classbond1(ic1),classbond2(ic1),bondrad(ic1)
         if(ierr.eq.0) ic=ic1

      end do
      close(20)

      nclassbond=ic

      return
      end

      subroutine gnubondswrite(nclassbond,classlab1,classlab2
     *,classbond1,classbond2,bondrad,bondrgb,boxrgb,boxlwidth)

      character*3 classlab1(*),classlab2(*)
      real classbond1(*),classbond2(*),bondrad(*)
     *,bondrgb(*),boxrgb(*),boxlwidth

      open(20,file='gnubonds.txt',status='unknown',iostat=ierr)

c Write bond and box colour indices and box line width

      write(20,347) (bondrgb(j),j=1,3),(boxrgb(j),j=1,3),boxlwidth
 347  format(7(1x,f10.5))

      ic=0
      do while (ic.lt.nclassbond.and.ierr.eq.0)

         ic=ic+1
         write(20,348,iostat=ierr) classlab1(ic),classlab2(ic)
     *,classbond1(ic),classbond2(ic),bondrad(ic)
348   format(1x,a3,1x,a3,3(1x,f10.5))

      end do
      close(20)

      return
      end

      subroutine run_GNUplotato(readfile,ilenreadfile
     *,nstereo,thestereo,phistereo,almin,almax)

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'plot_pars.inc'
      include 'system_commands.inc'
      include 'plotato.inc'

      character*256 readfile,plot_command,fnamegnuplot
      character*3 t1,t2
      character*2 extsplot

      integer ilenreadfile,nstereo

      real thestereo(2),phistereo(2)

c Write files suitable for the splot command in gnuplot

c First the list of bonds

      fnamesplot(0)=directory_name(1:nchardir)//readfile(1:ilenreadfile)//'.splb'
      ilenfnamesplot(0)=nchardir+ilenreadfile+5

      open(10,file=fnamesplot(0),status='unknown')

c Write a line of zeros so that splot does not draw cross-lines between bonds

      write(10,222) 0,0.0,0.0,0.0
      write(10,*)
222   format(1x,i5,3(1x,e13.6))

c Step through the atoms and write the coordinates of this atom

      do ia1=1,natom

         if(iplot(ia1).gt.0) then

c Type number for this atom

            t1=atsymbol(itype(ia1))

c Signal that a bond pair has been found

            ipairfound=0

c For this atom we need to see which atoms it is bonded to and list them 

            ia2=ia1
            do while (ia2.lt.natom)

c Reset the signal that a bond pair has been found in case another bond associated with
c ia1 is available

               ifound=0

c Only look at atoms that we want to plot

               ia2=ia2+1
               if(iplot(ia2).gt.0) then

c Type number for this atom

                  t2=atsymbol(itype(ia2))

c Separation between ia1 and ia2

                  xdif=xat(ia1)-xat(ia2)
                  ydif=yat(ia1)-yat(ia2)
                  zdif=zat(ia1)-zat(ia2)
                  dsq=xdif*xdif+ydif*ydif+zdif*zdif
                  dtest=sqrt(dsq)

c Test for bond type numbers and distance between to see whether we should include this
c atom in the list of neighbours

                  ib=0
                  do while(ib.lt.nclassbond.and.ifound.eq.0)
                        
                     ib=ib+1
                     if((t1.eq.classlab1(ib).and.t2.eq.classlab2(ib))
     *.or.
     *(t1.eq.classlab2(ib).and.t2.eq.classlab1(ib))) then

c Test whether the distance constraint for this atom pair is satisfied

                        if(dtest.ge.classbond1(ib).and.dtest.le.classbond2(ib)) then

                           ifound=ia2
                           write(10,222) ia1,xat(ia1),yat(ia1),zat(ia1)
                           write(10,222) ia2,xat(ia2),yat(ia2),zat(ia2)

c Need a blank line after each atom so that splot does not join them up

                           write(10,*)

                        endif
                     endif
                  end do
               endif

c Signal whether any bonds were found for this atom

               if(ipairfound.eq.0) ipairfound=ifound
            end do
         endif
      end do
      close(10)

c Now repeat this for each of the plotting classes specified, so that we can give the
c the different atom types different radii and colours.

      do ic=1,nclass
         write(extsplot,223) ic
223   format(i2.2)
         fnamesplot(ic)=directory_name(1:nchardir)//readfile(1:ilenreadfile)
     *//'.spl'//extsplot
         ilenfnamesplot(ic)=nchardir+ilenreadfile+6

         open(10,file=fnamesplot(ic),status='unknown')

c Step through the atoms and write the coordinates of this atom

         do ia1=1,natom
            if(iplot(ia1).eq.iplotclass(ic)) then

               write(10,222) ia1,xat(ia1),yat(ia1),zat(ia1)
               write(10,*)
            endif

         end do
         close(10)
      end do

c O.K. so now ready to setup the plot file

      fnamegnuplot=directory_name(1:nchardir)//readfile(1:ilenreadfile)
     *//'.gnu'
      ilenfnamegnuplot=nchardir+ilenreadfile+4

      open(11,file=fnamegnuplot,status='unknown')

c write the title and x, y and z labels

      title=directory_name(1:nchardir)//readfile(1:ilenreadfile)
     *//'.ato'
      ilentitle=nchardir+ilenreadfile+4
      write(11,301) 
301   format('reset')

c For stereo plots we set the initial size and origin

      if(nstereo.gt.1) write(11,320)
320   format('set size 1,1'/'set origin 0,0'/'set multiplot')
      x0=0.0
      xsize=0.5-x0
      write(11,*)

      do istereo=1,nstereo

         write(11,*)

         if(nstereo.gt.1) then

            write(11,324) x0
 324        format('set size 0.5,1'/'set origin ',f10.5,',0')
            x0=x0+xsize

         endif
         write(11,321)
 321     format('unset title')
         xlabel='x [\C5]'
         ilenxlabel=5
         write(11,302) xlabel(1:ilenxlabel)
 302     format('set xlabel ',1h',a,1h')

         ylabel='y [\C5]'
         ilenylabel=5
         write(11,303) ylabel(1:ilenylabel)
 303     format('set ylabel ',1h',a,1h')

         zlabel='z [\C5]'
         ilenzlabel=5
         write(11,304) zlabel(1:ilenzlabel)
 304     format('set zlabel ',1h',a,1h')

c       Scale factor on point size

         ptscale=ptscaler/(almax-almin)

c       For the z-axis we need to offset the xy plane slightly so that the z-axis does not clash with the
c       x,y axes. (See GNUplot help under "set ticslevel" for this formula.)

         posxy=posfrac*(almin-almax)+almin

         write(11,305) posfrac
 305     format('set ticslevel ',e13.6)
         write(11,306) posxy,almax
 306     format('set xrange [',e10.3,':',e10.3,']')
         write(11,307) posxy,almax
 307     format('set yrange [',e10.3,':',e10.3,']')
         write(11,308) almin,almax
 308     format('set zrange [',e10.3,':',e10.3,']')
         write(11,322) thestereo(istereo),phistereo(istereo)
 322     format('set view ',e13.6,',',e13.6)

c       Now plot the bonds and atoms

         write(11,*)
         write(11,309)
 309     format('splot \')

         do ifile=0,nclass

            if(ifile.eq.0) then

c       Plot the bonds

             write(11,310) fnamesplot(ifile)(1:ilenfnamesplot(ifile)),bondrad(1)
 310         format(1h',a,1h',' u 2:3:4 notitle w l lt 4 lw ',f10.5,$)

            else

               ic=iplotclass(ifile)

c       Assuming a class was found, which presumeably it must have been defined by now, we plot the
c       atoms belonging to that class.


                write(11,311) fnamesplot(ifile)(1:ilenfnamesplot(ifile))
     *,classtype(ic),classrad(ic)*ptscale,classcolour(ic)
 311            format(1h',a,1h',' u 2:3:4 notitle w p pt ',i5,' ps ',e13.6,' lt ',i5,$)

             endif

c       If this is not the last file, then add the continuation character

            if(ifile.lt.nclass) write(11,312)
 312        format(', \')

         end do
         write(11,*)
         write(11,*)
      end do

c Unset multiplot before finishing if necessary

      if(nstereo.gt.1) write(11,323)
323   format('unset multiplot')

c Finish off the plot file, and generate plot command

      close(11)

      plot_command=system_gnu(1:ilensystem_gnu)
     *//' '
     *//system_quotes(1:ilensystem_quotes)
     *//fnamegnuplot(1:ilenfnamegnuplot)
     *//system_quotes(1:ilensystem_quotes)
     *//' -'

      call system(plot_command,ierr)

      return
      end

      subroutine run_pgplotato(readfile,ilenreadfile)

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'plot_pars.inc'
      include 'system_commands.inc'
      include 'plotato.inc'

      character*256 readfile,plot_command
      character*3 t1,t2
      character*4 boxyesorno

      integer ilenreadfile

c Write files suitable for the MYPGPLOTATO program

c First write the plotting classes specified, so that we can give the
c the different atom types different radii and colours
c
c Also save a reference integer so that we can refer to each atom in the bonds file
c

      natomplot=0

      fnamesplot(1)=directory_name(1:nchardir)//readfile(1:ilenreadfile)
     *//'.pgpa'
      ilenfnamesplot(1)=nchardir+ilenreadfile+5

      open(10,file=fnamesplot(1),status='unknown')

c Step through the atoms and write the coordinates of this atom

      do ia1=1,natom

         if(iplot(ia1).gt.0) then

            ic=iplot(ia1)
            natomplot=natomplot+1
            atomref(ia1)=natomplot
            write(10,223) atsymbol(itype(ia1)),xat(ia1),yat(ia1),zat(ia1)
     * ,classrad(ic)*scaleatoms,(classrgb(j,ic),j=1,3),natomplot

         endif
      end do
      close(10)

 223    format(1x,a3,3(1x,e12.5),4(1x,f6.3),1x,i6)
c

c Now list the bonds

      fnamesplot(0)=directory_name(1:nchardir)//readfile(1:ilenreadfile)//'.pgpb'
      ilenfnamesplot(0)=nchardir+ilenreadfile+5

      open(10,file=fnamesplot(0),status='unknown')

221     format(1x,f10.5,2(1x,i6))

c Step through the atoms and write the coordinates of this atom

      do ia1=1,natom

         if(iplot(ia1).gt.0) then

c       Type name for this atom

            t1=atsymbol(itype(ia1))

c       Signal that a bond pair has been found

            ipairfound=0

c       For this atom we need to see which atoms it is bonded to and list them 

            ia2=ia1
            do while (ia2.lt.natom)

c       Reset the signal that a bond pair has been found in case another bond associated with
c       ia1 is available

             ifound=0

c       Only look at atoms that we want to plot

             ia2=ia2+1
             if(iplot(ia2).gt.0) then

c       Type number for this atom

                t2=atsymbol(itype(ia2))

c       Separation between ia1 and ia2

                xdif=xat(ia1)-xat(ia2)
                ydif=yat(ia1)-yat(ia2)
                zdif=zat(ia1)-zat(ia2)
                dsq=xdif*xdif+ydif*ydif+zdif*zdif
                dtest=sqrt(dsq)

c       Test for bond type numbers and distance between to see whether we should include this
c       atom in the list of neighbours

                ib=0
                do while(ib.lt.nclassbond.and.ifound.eq.0)
                   
                   ib=ib+1
                   if((t1.eq.classlab1(ib).and.t2.eq.classlab2(ib))
     *.or.(t1.eq.classlab2(ib).and.t2.eq.classlab1(ib))) then

c       Test whether the distance constraint for this atom pair is satisfied

                    if(dtest.ge.classbond1(ib).and.dtest.le.classbond2(ib)) then

                       ifound=ia2
                       write(10,221) bondrad(ib),atomref(ia1),atomref(ia2)

                    endif
                   endif
                end do
             endif

c       Signal whether any bonds were found for this atom

             if(ipairfound.eq.0) ipairfound=ifound
            end do
         endif
      end do
      close(10)

      write(6,202)
 202  format(/'plotato> Draw sides of box? (y or n): ',$)
      read(5,*) boxyesorno

c
c now generate the .txt file from which pgplotato will get its parameters
c
      open(11,file=directory_name(1:nchardir)//'mypgplotato.txt'
     *,status='unknown',iostat=ierr)
      if(ierr.eq.0) then
         write(11,103) fnamesplot(1)(1:ilenfnamesplot(1))
         write(11,103) fnamesplot(0)(1:ilenfnamesplot(0))
103   format(a)

         write(11,103) system_pgout(1:ilensystem_pgout)
         write(11,107) xangle,yangle
107   format(7(1x,f10.5))
         write(11,302) boxyesorno
         write(11,107) (bondrgb(j),j=1,3),(boxrgb(j),j=1,3),boxlwidth
302   format(a)
         close(11)
      else
         write(6,104)
104   format(/'plotato> Not able to open "myplotato.txt".')
         return
      endif
c
c write the name of this file to plot.txt
c
        open(11,file=directory_name(1:nchardir)//'plot.txt'
     *,status='unknown',iostat=ierr)
        write(11,103) 'mypgplotato.txt'
      close(11)
c
c now generate the command that will actually run this file
c
      plot_command=system_cd(1:ilensystem_cd)
     *//' '
     *//system_quotes(1:ilensystem_quotes)
     *//directory_name(1:nchardir)
     *//system_quotes(1:ilensystem_quotes)
     *//' '//system_join(1:ilensystem_join)
     *//' '//system_binaries(1:ilensystem_binaries)//system_plotato(1:ilensystem_plotato)
     *//' < plot.txt'
      write(6,*) plot_command
      call system(plot_command,ierr)
      write(6,*) ierr

      return
      end

      subroutine run_Jmol(readfile,ilenreadfile)

      use molecom_mod

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'plot_pars.inc'
      include 'system_commands.inc'
      include 'plotato.inc'
!      include 'lattice.inc'

      character*256 readfile,plot_command
      integer ilenreadfile,nx,ny,nz,nx0,ny0,nz0,npixel,nvoxel,i,j,k,ix,iy,iz,it,itmax
      real :: xtmax,ytmax,ztmax,ltmax,ltstep,bohr,sigma,sum1,xtdif,ytdif,ztdif,rtmax,rtmax2,rtest,denom,gfactor,oversize,scalefc
      real, dimension(:,:,:), allocatable :: dens
c       
c generate name of output file for Jmol.
c
      fnamesplot(0)=directory_name(1:nchardir)//readfile(1:ilenreadfile)//'.cube'
      ilenfnamesplot(0)=len_trim(fnamesplot(0))
! Get some extra information
      write(6,'(a)',advance='no') 'plotato_surface> Give number of pixels in longest direction, sigma for broadening, oversize factor, scale factor : '
      read(5,*) npixel,sigma,oversize,scalefc
! Get maximum dimension in each direction
!      xtmax=0.0
!      ytmax=0.0
!      ztmax=0.0
!      do ia=1,natom
!         if(iplot(ia).gt.0) then
!!            xtmax=max(xtmax,abs(xat(ia)))
!            ytmax=max(ytmax,abs(yat(ia)))
!            ztmax=max(ztmax,abs(zat(ia)))
!         end if
!      end do
!      xtmax=xtmax*oversize
!      ytmax=ytmax*oversize
!      ztmax=ztmax*oversize
      xtmax=dmax*oversize
      ytmax=dmax*oversize
      ztmax=dmax*oversize
      ltmax=max(xtmax,ytmax,ztmax) !Allow a bit extra to avoid truncation at box boundaries
      ltstep=ltmax/real(npixel)
      rtmax=4.0*sigma
      rtmax2=rtmax*rtmax
      denom=0.5/(sigma*sigma)
      pi=4.0*atan(1.0)
      gfactor=1.0/(sigma*sqrt(2.0*pi))
      gfactor=gfactor*gfactor*gfactor !Need cube because this is a three-dimensional density
      nx0=0
      ny0=0
      nz0=0
      if(ltstep.gt.0.0) then
         nx=nint(xtmax/ltstep)
         ny=nint(ytmax/ltstep)
         nz=nint(ztmax/ltstep)
         nx0=2*nx+1
         ny0=2*nx+1
         nz0=2*nx+1
      end if
      nvoxel=nx0*ny0*nz0
      write(6,*) xtmax,ytmax,ztmax,ltmax,ltstep,rtmax
      write(6,*) nx,ny,nz,nx0,ny0,nz0,nvoxel
      bohr=0.52917721092
c Write Gaussian cube output file

      nout=10
      open(nout,file=fnamesplot(0),status='unknown',iostat=ierr)
      if(ierr.eq.0.and.nvoxel.gt.0) then
         if(allocated(dens)) deallocate(dens)
         allocate(dens(nx0,ny0,nz0))
! Write the preamble
         write(nout,*) '# Cube file output by plotato_surface'
         write(nout,*) 'Outer loop x, middle loop y, inner loop z'
         write(nout,*) 0,-xtmax*scalefc/bohr,-ytmax*scalefc/bohr,-ztmax*scalefc/bohr
         write(nout,*) nx0, ltstep*scalefc/bohr, 0.0, 0.0
         write(nout,*) ny0, 0.0, ltstep*scalefc/bohr, 0.0
         write(nout,*) nz0, 0.0, 0.0, ltstep*scalefc/bohr
! Step through the pixels
         do k=1,nz0
            z=real(k-nz-1)*ltstep
            z2=z*z
            do j=1,ny0
               y=real(j-ny-1)*ltstep
               y2=y*y
               do i=1,nx0
                  x=real(i-nx-1)*ltstep    
                  x2=x*x
                  rxy=sqrt(x2+y2)
                  r=sqrt(x2+y2+z2)
                  sum1=0.0
                  do ia=1,natom
                     if(iplot(ia).gt.0) then
                        xtdif=x-xat(ia)
                        ytdif=y-yat(ia)
                        ztdif=z-zat(ia)
                        if(max(xtdif,ytdif,ztdif).lt.rtmax) then
                           rtest=xtdif*xtdif+ytdif*ytdif+ztdif*ztdif
!                           if(rtest.lt.rtmax2) then
                              arg=rtest*denom
                              if(arg.lt.30.0) sum1=sum1+natomsphere(itype(ia))*pseudofacmol(itype(ia))*exp(-arg)
!                           end if
                        end if
                     end if
                  end do
                  dens(i,j,k)=gfactor*sum1
               end do
            end do
         end do
         do ix=1,nx0
            do iy=1,ny0
               do iz=1,nz0,10
                  itmax=min(9,nz0-iz)
                  write(nout,'(10(1x,e13.6))') (dens(ix,iy,iz+it),it=0,itmax)
               end do
            end do
         end do
         close(nout)
      end if

      return
      end

