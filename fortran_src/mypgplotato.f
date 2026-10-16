      include 'epsr_dimension.inc'
      parameter (mbonds=matom*20)

      real*4 xbal(matom),ybal(matom),zbal(matom),radbal(mcomp)
     *,rgbbal(3,mcomp) 
      real*4 radrod(mbonds),rgb(3),bondrgb(3),boxrgb(3),boxlwidth

      integer*4 nbal,nbaltype,itype(matom)

      integer*4 nrod,rod1(mbonds),rod2(mbonds)

      character*256 infile,fname

      character*10 device

      character*3 dummyc,labelbal(mcomp)
		character*4 boxyesorno

      external cosd,sind,acosd,asind,atan2d

c Initialise random number generator

      open(10,file='randomno.dat',status='unknown')
      call initran1(10)
      close(10)

      write(6,98) 'Type name of file to read plotting parameters from'
      read(5,98) fname
      open(10,file=fname,status='old')
      write(6,99)
99    format(1x,'Type name of file containing atom coordinates')
      nin=10
      read(nin,98) infile
98    format(a)

c Read the atom coordinates

c Set the maximum distance found

      dmax=0.0

c Set the smallest radius

      radmin=0.0

c Set the number of atom types found

      nbaltype=0

C Set the number of atoms found

      nbal=0
      
      open(11,file=infile,status='old',iostat=ierr)

         do while (ierr.eq.0.and.nbal.lt.matom)
            ib1=nbal+1
            read(11,*,iostat=ierr) dummyc,xbal(ib1),ybal(ib1),zbal(ib1)
     *,rad,(rgb(j),j=1,3)

            if(ierr.eq.0) then
             nbal=ib1
             if(abs(xbal(nbal)).gt.dmax) dmax=abs(xbal(nbal))
             if(abs(ybal(nbal)).gt.dmax) dmax=abs(ybal(nbal))
             if(abs(zbal(nbal)).gt.dmax) dmax=abs(zbal(nbal))

c Now check the existing types to see which exist

             it=0
             ifound=0
             do while (it.lt.nbaltype.and.ifound.eq.0)
                
                it=it+1
                if(dummyc.eq.labelbal(it)) ifound=it

             end do

             if(ifound.eq.0) then

c Find the smallest radius

                if(nbaltype.eq.1) radmin=rad
                if(rad.lt.radmin) radmin=rad

                nbaltype=nbaltype+1
                labelbal(nbaltype)=dummyc
                radbal(nbaltype)=rad
                do j=1,3
                   rgbbal(j,nbaltype)=rgb(j)
                end do
                ifound=nbaltype

             endif

             itype(nbal)=ifound

            endif
         end do
        
      close(11)

      write(6,97) nbal
97    format(1x,'Number of atoms in file = ',i5)
      
c Round up dmax to the nearest whole number

      dmax=int(dmax)+1

      write(6,101) 2*dmax
101   format(' Size of plot: ',f10.5)

c Return at this point if there is nothing to plot

      if(dmax.eq.0.0) stop

c Read the bonds data

      read(nin,98) infile

      open(11,file=infile,status='old',iostat=ierr)

         nrod=0
         do while (ierr.eq.0)
            ib1=nrod+1
            read(11,*,iostat=ierr) radrod(ib1),rod1(ib1),rod2(ib1)

c Ensure that, if radmin has been set, the bond radius does not exceed the smallest atom radius

            if(radmin.gt.0.0.and.radrod(ib1).gt.radmin) radrod(ib1)=radmin

            if(ierr.eq.0) then
             nrod=ib1
            endif
         end do
        
      close(11)

      write(6,96) nrod
 96   format(1x,'Number of bonds in file = ',i5)

c Read the plotting device

      read(nin,95) device
 95   format(a)
      ilendevice=index(device,' ')-1
      write(6,*) device(1:ilendevice)

c Read the orientation of the plot

      read(nin,*) anglex,angley

c Read whether to plot the box sides

      read(nin,95) boxyesorno

c Read colour indices for bonds and box sides

      read(nin,*) (bondrgb(j),j=1,3),(boxrgb(j),j=1,3),boxlwidth
      
      write(6,569)
569   format(1x,'Type number of plots in x and y directions on one page')
c     read(nin,*) nplotx,nploty
      nplotx=1
      nploty=1
      aspide=float(nploty)/float(nplotx)
      write(6,568)
568   format(1x,'Type aspect ratio (1.4=Portrait, 0.7=Landscape)')
c     read(nin,*) aspect
      aspect=1.0
      asprat=aspide/aspect
        CALL PGBEGIN(0,device(1:ilendevice),1,1)
        CALL PGPAPER(0.0,aspect)
      call pgpage
      yvwid=1.0/float(nploty)
      yvwid2=0.5*yvwid
      if(asprat.lt.1.0) yvwid2=yvwid2*asprat
      xvwid=1.0/float(nplotx)
      xvwid2=0.5*xvwid
      if(asprat.gt.1.0) xvwid2=xvwid2/asprat
      do ipy=1,nploty
         ycent=(nploty-ipy+1-0.5)*yvwid
         yv1=ycent-yvwid2
         yv2=ycent+yvwid2
         do ipx=1,nplotx
c       
c       set up coordinates of centre of viewport
c       
            xcent=(ipx-0.5)*xvwid
            xv1=xcent-xvwid2
            xv2=xcent+xvwid2
            write(6,992) xv1,xv2,yv1,yv2
 992        format(4(1x,e13.6))
            CALL PGVPORT(xv1,xv2,yv1,yv2)
            size=0.1
            iwidth=1
            icoll=1
            call plot3d(nin,size,iwidth,icoll,dmax,nbal
     *,xbal,ybal,zbal,itype,nbaltype,radbal,rgbbal
     *,radrod,rod1,rod2,nrod
     *,anglex,angley,boxyesorno,bondrgb,boxrgb,boxlwidth)
            call pgwindow(0.0,1.0,0.0,1.0)
         end do
      end do
      close(nin)

c Save random number generator

      open(10,file='randomno.dat',status='unknown')
         call saveran1(10)
        close(10)

      call pgend
      end               
