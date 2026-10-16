      module cellneighbourlist
      
      implicit none

      integer :: ncelldiva,ncelldivb,ncelldivc,ncell
      integer, dimension(:), allocatable :: nmolsincell
            
      
      contains
      
!***********************************************************************************************

      subroutine initialise_cellneighbourlist(ncelldivin,xadd,yadd,zadd,nmol,molcellno)
      
      implicit none

      include 'lattice.inc'
      
      integer :: ncelldivin,nmol,molcellno(*),ic,im,maxcellno,atomcellno
      real :: cellsize,xadd(*),yadd(*),zadd(*)
      character(len=256) :: sometext
      
      !ncelldivin is the number of cells along the shortest axis. We want the cells to be 
      !as close to cubic in each direction so use the Cartesian dimensions to determine cell size
      cellsize=min(alatx,blaty,clatz)/ncelldivin
      ncelldiva=nint(alatx/cellsize)
      ncelldivb=nint(blaty/cellsize)
      ncelldivc=nint(clatz/cellsize)
      ncell=ncelldiva*ncelldivb*ncelldivc
      write(6,100) ncelldiva,ncelldivb,ncelldivc,ncell
100   format(/'initialise_cellneighbourlist> Cell numbers along a,b,c and total',3(1x,i6))
      !Assign the array that will store the number of molecules in each cell.
      if(allocated(nmolsincell)) deallocate(nmolsincell)
      allocate(nmolsincell(ncell))
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
101   format(/'initialise_cellneighbourlist> Maximum number of molecules per cell: ',i6
     */'Cell number and no. of molecules in cell.')
      do ic=1,ncell
         write(6,'(2(1x,i6))') ic,nmolsincell(ic)
      end do
      
      return
      end subroutine initialise_cellneighbourlist 


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

c************************************************************************************************         

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
      write(6,*) ia,ib,ic,atomcellno,ncelldiva,ncelldivb,ncelldivc,ncell

      return
      end function atomcellno


c************************************************************************************************         

      subroutine get_date_and_time(text,start)

      implicit none
      
      character*80 text
      character*11 date
      character*12 time
      character*3 month,get_currentmonth

      common /time/startyear,startabsolutetime,startmusecs

      integer values(8),year,daysthisyearsofar,get_daysthisyearsofar
      integer hoursthisyearsofar,minsthisyearsofar,secsthisyearsofar
      integer start,startyear,startabsolutetime,currentyear,sumdays,sumsecs
      integer musecs,startmusecs

      call date_and_time(VALUES=values)

      currentyear=values(1)
      month=get_currentmonth(values(2))

      write(date,701) values(3),month,values(1)
701   format(i2,'-',a3,'-',i4)
      write(time,702) values(5),values(6),values(7),values(8)
702   format(i2.2,':',i2.2,':',i2.2,'.',i3.3)
      write(6,703) text(1:index(text,'#')-1),time,date
703   format(a,2x,a12,1x,a11)

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
      
