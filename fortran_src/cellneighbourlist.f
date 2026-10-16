      module cellneighbourlist
      
      implicit none

      integer :: ncelldiva,ncelldivb,ncelldivc,ncelldivsq,ncell,ncellneighb,maxcellno,nmolcell
      integer, dimension(:), allocatable :: nmolsincell,molcellno,iacellneighb,ibcellneighb,iccellneighb
      integer, dimension(:,:), allocatable :: molsincell
            
      
      contains
      
!***********************************************************************************************

      subroutine initialise_cellneighbourlist(nmolcellin,xadd,yadd,zadd,nmol,rmax)
      
      implicit none

      include 'lattice.inc'
      
      integer :: nmolcellin,ncelldiv,nmol,ic,im,atomcellno
      real :: cellsize,xadd(*),yadd(*),zadd(*),rmax
      
      !nmolcell is the (approximate) number of molecules in a cell. 
      nmolcell=nmolcellin
      if(nmolcell.gt.nmol) nmolcell=nmol
      !ncelldivin is the number of cells along the shortest axis. We want the cells to be 
      !as close to cubic in each direction so use the Cartesian dimensions to determine cell size
      ncelldiv=nint((real(nmol)/real(nmolcell))**(1.0/3.0))
      if(ncelldiv.lt.1) ncelldiv=1
      cellsize=2*min(alatx,blaty,clatz)/ncelldiv
      ncelldiva=2*int(alatx/cellsize)+1
      ncelldivb=2*int(blaty/cellsize)+1
      ncelldivc=2*int(clatz/cellsize)+1
      ncell=ncelldiva*ncelldivb*ncelldivc
      ncelldivsq=ncelldiva*ncelldivb
      write(6,100) ncelldiva,ncelldivb,ncelldivc,ncell
100   format(/'initialise_cellneighbourlist> Cell numbers along a,b,c and total',3(1x,i6))
      !Assign the array that will store the number of molecules in each cell.
      if(allocated(nmolsincell)) deallocate(nmolsincell)
      if(allocated(molcellno)) deallocate(molcellno)
      if(allocated(iacellneighb)) deallocate(iacellneighb)
      if(allocated(ibcellneighb)) deallocate(ibcellneighb)
      if(allocated(iccellneighb)) deallocate(iccellneighb)
      allocate(nmolsincell(ncell))
      allocate(molcellno(nmol))
      allocate(iacellneighb(ncell))
      allocate(ibcellneighb(ncell))
      allocate(iccellneighb(ncell))
      do ic=1,ncell
         nmolsincell(ic)=0
      end do
      !Assign molecules to cells, storing the number of molecules in each cell
      maxcellno=0
      do im=1,nmol
         molcellno(im)=atomcellno(xadd(im),yadd(im),zadd(im))
         nmolsincell(molcellno(im))=nmolsincell(molcellno(im))+1
         if(nmolsincell(molcellno(im)).gt.maxcellno) maxcellno=nmolsincell(molcellno(im))
      end do
      write(6,101) maxcellno
101   format(/'initialise_cellneighbourlist> Maximum number of molecules per cell: ',i6)
!     */'Cell number and no. of molecules in cell.')
      !Need to make a list of the molecules in each cell. We assume there will never be more than the current maximum
      !in any one cell
      if(allocated(molsincell)) deallocate(molsincell)
      allocate(molsincell(maxcellno,ncell))
      do ic=1,ncell
         nmolsincell(ic)=0
      end do
      do im=1,nmol
         nmolsincell(molcellno(im))=nmolsincell(molcellno(im))+1
         molsincell(nmolsincell(molcellno(im)),molcellno(im))=im
      end do
!      do ic=1,ncell
!         write(6,*) 'Cell ',ic,(molsincell(i,ic),i=1,nmolsincell(ic))
!      end do
      
      !Make the list of neighbouring cells to any one cell at the origin
      call makecellneighbourlist(rmax)
      
      return
      end subroutine initialise_cellneighbourlist 

c************************************************************************************************         

      subroutine makecellneighbourlist(rmax)

      implicit none

c Gets the list of neighbouring cells relative to the one at (0,0,0) which satisfy the required
c distance criterion. Cells which are periodic images of cells already in the list are ignored.

      include 'lattice.inc'
      real :: rmax,rmaxsq,rdiffsq !Maximum separation between any atoms in one cell with those in a neighbouring cell
      real :: xcorner(8),ycorner(8),zcorner(8) !Cartesian coordinates of the corners of each cell relative to the centre of the cell
      real :: cella,cellb,cellc,cella2,cellb2,cellc2 !Half fractions of cell sides
      real :: cx,cy,cz,bx,by,bz,ax,ay,az,xdiff,ydiff,zdiff,xdiffc,ydiffc,zdiffc  !Dummy coordinates
      integer :: ia,ib,ic,i,if1,if2,nrangea,nrangeb,nrangec,ifound,ncorner
      logical :: found

      cella=1.0/real(ncelldiva)
      cellb=1.0/real(ncelldivb)
      cellc=1.0/real(ncelldivc)
      cella2=0.5*cella
      cellb2=0.5*cellb
      cellc2=0.5*cellc
      !Set up coordinates of corner of cell
      ifound=0
      do ic=-1,1,2
         cx=ic*cellc2*clatc(1)
         cy=ic*cellc2*clatc(2)
         cz=ic*cellc2*clatc(3)
         do ib=-1,1,2
            bx=ib*cellb2*blatc(1)
            by=ib*cellb2*blatc(2)
            bz=ib*cellb2*blatc(3)
            do ia=-1,1,2
               ax=ia*cella2*alatc(1)
               ay=ia*cella2*alatc(2)
               az=ia*cella2*alatc(3)
               ifound=ifound+1
               xcorner(ifound)=ax+bx+cx
               ycorner(ifound)=ay+by+cy
               zcorner(ifound)=az+bz+cz
            end do
         end do
      end do
      ncorner=ifound
!      do i=1,ncorner
!         write(6,*) i,xcorner(i),ycorner(i),zcorner(i)
!      end do
      !Now get the list of neighbouring cells, the shortest distance between which is less than rmax
      nrangea=(ncelldiva-1)/2
      nrangeb=(ncelldivb-1)/2
      nrangec=(ncelldivc-1)/2
      ifound=0
      rmaxsq=rmax*rmax
      do ic=-nrangec,nrangec
         cx=ic*cellc*clatc(1)
         cy=ic*cellc*clatc(2)
         cz=ic*cellc*clatc(3)
         do ib=-nrangeb,nrangeb
            bx=ib*cellb*blatc(1)
            by=ib*cellb*blatc(2)
            bz=ib*cellb*blatc(3)
            do ia=-nrangea,nrangea
               ax=ia*cella*alatc(1)
               ay=ia*cella*alatc(2)
               az=ia*cella*alatc(3)
               xdiffc=ax+bx+cx
               ydiffc=ay+by+cy
               zdiffc=az+bz+cz
               !step through all the possible combinations of corners to determine the shortest distance
               found=.false.
               if1=0
               do while(if1.lt.ncorner.and..not.found)
                  if1=if1+1
                  if2=0
                  do while(if2.lt.ncorner.and..not.found)
                     if2=if2+1
                     xdiff=xdiffc+xcorner(if1)-xcorner(if2)
                     ydiff=ydiffc+ycorner(if1)-ycorner(if2)
                     zdiff=zdiffc+zcorner(if1)-zcorner(if2)
                     rdiffsq=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
                     found=rdiffsq.lt.rmaxsq.or.(ia.eq.0.and.ib.eq.0.and.ic.eq.0) !Always include the cell at the origin
                  end do
               end do
               if(found) then
!                  write(6,*) ia,ib,ic
                  ifound=ifound+1
                  iacellneighb(ifound)=ia
                  ibcellneighb(ifound)=ib
                  iccellneighb(ifound)=ic                 
               end if
            end do
         end do
      end do
      ncellneighb=ifound
      write(6,100) ifound,ncell
100   format(/'makecellneighbourlist> Found ',i6,' cell neighbours out of ',i6,' total') 
        
      return
      end subroutine makecellneighbourlist
      
      subroutine make_neighbour_list_from_cells(ic,nnbref,nbref)
     
      integer :: ic,nnbref,iref,icref,icellref
     *,iaorigin,iborigin,icorigin,ianeighb,ibneighb,icneighb
      integer, dimension(*) :: nbref

      nnbref=0
      !Get the coordinates of cell ic
      call findiaibic(ic,iaorigin,iborigin,icorigin)
      !Step through all neighbouring cells, and generate the corresponding cell number
      do icref=1,ncellneighb
         ianeighb=iaorigin+iacellneighb(icref)
         ibneighb=iborigin+ibcellneighb(icref)
         icneighb=icorigin+iccellneighb(icref)
         call getcellreference(ianeighb,ibneighb,icneighb,icellref)
         !Get all the molecules in this cell.
         do iref=1,nmolsincell(icellref)
            nnbref=nnbref+1
            nbref(nnbref)=molsincell(iref,icellref)
         end do
      end do
      return
           
      end subroutine make_neighbour_list_from_cells
      
      subroutine get_nnbrefmax(nnbrefmax)
!To find the maximum number of neighbouring molecules for any cell (including itself).     
      integer :: ic,nnbref,nnbrefmax,icref,icellref
     *,iaorigin,iborigin,icorigin,ianeighb,ibneighb,icneighb

      nnbrefmax=0
      do ic=1,ncell
         nnbref=0
         !Get the coordinates of cell ic
         call findiaibic(ic,iaorigin,iborigin,icorigin)
         !Step through all neighbouring cells, and generate the corresponding cell number
         do icref=1,ncellneighb
            ianeighb=iaorigin+iacellneighb(icref)
            ibneighb=iborigin+ibcellneighb(icref)
            icneighb=icorigin+iccellneighb(icref)
            call getcellreference(ianeighb,ibneighb,icneighb,icellref)
            !Get all the molecules in this cell.
            nnbref=nnbref+nmolsincell(icellref)
         end do
         if(nnbref.gt.nnbrefmax) nnbrefmax=nnbref
      end do
      return
           
      end subroutine get_nnbrefmax
      
      subroutine findiaibic(cellno,ia,ib,ic)

      !Convert cellnumber to the corresponding a,b,c coordinates 

      integer :: cellno,cellno1,ia,ib,ic

      cellno1=cellno-1
      ic=int(cellno1/ncelldivsq)
      ib=int((cellno1-ic*ncelldivsq)/ncelldiva)
      ia=cellno1-ncelldiva*(ib+ncelldivb*ic)

      return
      end subroutine findiaibic

      subroutine getcellreference(ia,ib,ic,icellref)

      integer :: ia,ib,ic,icellref

      !Generates a cell reference based on ia,ib,ic, after correcting the latter for
      !periodic boundary conditions.
      call findperiodiccell(ia,ncelldiva)
      call findperiodiccell(ib,ncelldivb)
      call findperiodiccell(ic,ncelldivc)
      icellref=1+ia+ncelldiva*(ib+ncelldivb*ic)

      return
      end subroutine getcellreference
      
      subroutine findperiodiccell(ix,ncelldiv)

      integer ix,ncelldiv

      !Gets the periodic cell coordinate based number of cell divisions

      do while(ix.lt.0)
         ix=ix+ncelldiv
      end do
      do while(ix.ge.ncelldiv)
         ix=ix-ncelldiv
      end do

      return
      end subroutine findperiodiccell

      end module cellneighbourlist      
      
c************************************************************************************************         

      integer function cellcoord(x,ncelldiv)

c Returns the appropriate cell coordinate, applying periodic boundary conditions

      use cellneighbourlist
      implicit none

      integer ix,ncelldiv
      real x,x1

      x1=x-nint(x)
      ix=int((x1+0.5)*real(ncelldiv))
      do while(ix.lt.0)
         ix=ix+ncelldiv
      end do
      do while(ix.ge.ncelldiv)
         ix=ix-ncelldiv
      end do
      x=x1
      cellcoord=ix

      return
      end function cellcoord

************************************************************************************************         

      integer function atomcellno(x,y,z)

      use cellneighbourlist
      implicit none
      include 'lattice.inc'
      integer :: cellcoord,ia,ib,ic
      real :: x,y,z,fa,fb,fc
      
      !From input Cartesian coordinates, gets cell number

      !First convert coordinates to fractional coordinates along a,b and c axes
      fc=z/clatz
      fb=(y-fc*claty)/blaty
      fa=(x-fb*blatx-fc*clatx)/alatx
      ia=cellcoord(fa,ncelldiva)
      ib=cellcoord(fb,ncelldivb)
      ic=cellcoord(fc,ncelldivc)

      atomcellno=1+ia+ncelldiva*(ib+ncelldivb*ic)
!      write(6,*) ia,ib,ic,atomcellno,ncelldiva,ncelldivb,ncelldivc,ncell

      return
      end function atomcellno



c************************************************************************************************         


      subroutine get_date_and_time(text,start)

      implicit none
      
      character(len=*) text
      character*11 date
      character*12 time
      character*3 month,get_currentmonth

      common /time/startyear,startabsolutetime,startmusecs

      integer values(8),year,daysthisyearsofar,get_daysthisyearsofar
      integer hoursthisyearsofar,minsthisyearsofar,secsthisyearsofar
      integer start,startyear,startabsolutetime,currentyear,sumdays,sumsecs
      integer musecs,startmusecs,i

      call date_and_time(VALUES=values)
!      write(6,'(a,8(1x,i5))') 'get_date_and_time> ',(values(i),i=1,8)

      currentyear=values(1)
      month=get_currentmonth(values(2))

      write(date,701) values(3),month,values(1)
701   format(i2,'-',a3,'-',i4)
      write(time,702) values(5),values(6),values(7),values(8)
702   format(i2.2,':',i2.2,':',i2.2,'.',i3.3)
!      write(6,'(a,1x,a,1x,a,1x,a)') 'get_date_and_time> ',text(1:len_trim(text)),date,time
!      write(6,'(a,1x,a,1x,a,1x,a)') 'get_date_and_time> ',text(1:index(text,'#')-1),date,time
      write(6,703) text(1:index(text,'#')-1),time,date
703   format(a,2x,a,1x,a)

c Generate the total number of days this year so far, not including the current day

      daysthisyearsofar=get_daysthisyearsofar(currentyear,values(2))+values(3)-1
      hoursthisyearsofar=daysthisyearsofar*24+values(5)
      minsthisyearsofar=hoursthisyearsofar*60+values(6)
      secsthisyearsofar=minsthisyearsofar*60+values(7)
      musecs=values(8)
c      write(6,*) get_daysthisyearsofar(year,values(2)),daysthisyearsofar

c Calculate elapsed time

      if(start.eq.0) then

c Initialise start time

         startyear=currentyear
         startabsolutetime=secsthisyearsofar
         startmusecs=values(8)

      else

c Calculate the elapsed time

c First calculate the effect if the year number has changed between start and finish
c We count the number of seconds from the beginning of the start year to the end of the
c year before the current year

         year=startyear
         sumdays=0
         do while(year.lt.currentyear)
            sumdays=sumdays+get_daysthisyearsofar(year,13)
            year=year+1
         end do
c Convert this to seconds
         sumsecs=sumdays*(24*3600)
c Subtract the previous absolute time and add on the number of secs this year
         sumsecs=sumsecs-startabsolutetime+secsthisyearsofar
         if(startmusecs.gt.musecs) then
            sumsecs=sumsecs-1
            musecs=1000+musecs-startmusecs
         else
            musecs=musecs-startmusecs
         endif

         write(6,704) sumsecs,musecs
704   format(1x,'Elapsed time (secs): ',i10,'.',i3.3)

c Store the new values of startyear and startabsolute time
         startyear=currentyear
         startabsolutetime=secsthisyearsofar
         startmusecs=values(8)
      
      endif
      return
      end


      function get_daysthisyearsofar(year,latestmonth)

c Calculates the number of days in the specified year up to but not including latestmonth. 
c To get all months set latestmonth > 12

      implicit none

      integer get_daysthisyearsofar
      integer year,latestmonth,daysthisyearsofar,dayspermonth,i

      if(latestmonth.gt.12) latestmonth = 13
      i=0
      daysthisyearsofar=0
      do while(i.lt.latestmonth-1)

         dayspermonth=0
         i=i+1
         if(i.eq.1) then
            dayspermonth=31
         else if(i.eq.2) then
c Check for leap year
            if((year-4*int(year/4)).eq.0) then
               dayspermonth=29
            else
               dayspermonth=28
            endif
         else if(i.eq.3) then
            dayspermonth=31
         else if(i.eq.4) then
            dayspermonth=30
         else if(i.eq.5) then
            dayspermonth=31
         else if(i.eq.6) then
            dayspermonth=30
         else if(i.eq.7) then
            dayspermonth=31
         else if(i.eq.8) then
            dayspermonth=31
         else if(i.eq.9) then
            dayspermonth=30
         else if(i.eq.10) then
            dayspermonth=31
         else if(i.eq.11) then
            dayspermonth=30
         else if(i.eq.12) then
            dayspermonth=31
         endif
         daysthisyearsofar=daysthisyearsofar+dayspermonth
      end do

      get_daysthisyearsofar=daysthisyearsofar

      return
      end

      function get_currentmonth(values)

      implicit none
      character*3 get_currentmonth,month
      integer values

c Get the current month

      if(values.eq.1) then
         month='JAN'
      else if(values.eq.2) then
         month='FEB'
      else if(values.eq.3) then
         month='MAR'
      else if(values.eq.4) then
         month='APR'
      else if(values.eq.5) then
         month='MAY'
      else if(values.eq.6) then
         month='JUN'
      else if(values.eq.7) then
         month='JUL'
      else if(values.eq.8) then
         month='AUG'
      else if(values.eq.9) then
         month='SEP'
      else if(values.eq.10) then
         month='OCT'
      else if(values.eq.11) then
         month='NOV'
      else if(values.eq.12) then
         month='DEC'
      endif

      get_currentmonth=month

      return
      end
      
