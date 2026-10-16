      module cellneighbourlist
      
      integer :: ncelldiva,ncelldivb,ncelldivc,ncell
      integer, dimension(:), allocatable :: nmolsincell
            
      implicit none
      
      contains
      
      subroutine initialise_cellneighbourlist(ncelldivin,xadd,yadd,zadd,nmol,molcellno)
      
      include 'lattice.inc'
      
      integer :: ncelldivin,nmol,molcellno(*),icell,im,maxcellno,atomcellno
      real :: cellsize,xadd(*),yadd(*),zadd(*)
      character(len=256) :: sometext
      
      !ncelldivin is the number of cells along the shortest axis. We want the cells to be 
      !as close to cubic in each direction so use the Cartesian dimensions to determine cell size
      cellsize=min(alatx,blaty,clatz)/ncelldivin
      ncelldiva=nint(alatx/cellsize)
      ncelldivb=nint(blaty/cellsize)
      ncelldivc=nint(clatz/cellsize)
      ncell=ncelldiva*ncelldivb*ncelldivc
      write(6,100) ncelldiva,ncelldivb,ncelldivc,ncelltotal
100   format(/'initialise_cellneighbourlist> Cell numbers along a,b,c and total',3(1x,i))
      !Assign the array that will store the number of molecules in each cell.
      if(allocated(nmolsincell)) deallocate(nmolsincell)
      allocate(nmolsincell(ncell))
      do icell=1,ncell
         nmolsincell(icell)=0
      end do
      !Assign molecules to cells, storing the maximum number of molecules in a cell
      maxcellno=0
      do im=1,nmol
         molcellno(im)=atomcellno(xadd(im),yadd(im),zadd(im))
         if(molcellno(im).gt.maxcellno) maxcellno=molcellno(im)
         nmolsincell(molcellno(im))=nmolsincell(molcellno(im))+1
      end do
      write(6,101) maxcellno
101   format(/'initialise_cellneighbourlist> Maximum number of molecules per cell: ',i
     */'Cell number and no. of molecules in cell.')
      do ic=1,ncell
         write(6,'(2(1x,i))') ic,nmolsincell(ic)
      end do
      
      return
      end subroutine initialise_cellneighbourlist 

c************************************************************************************************         

      function atomcellno(x,y,z)

      implicit none
      include 'lattice.inc'
      integer :: atomcellno,cellcoord,ia,ib,ic
      real :: x,y,z,fa,fb,fc
      
      !From input Cartesian coordinates, gets cell number

      !First convert coordinates to fractional coordinates along a,b and c axes
      fc=z/clatz
      fb=(y-fc*claty)/blaty
      fa=(x-fb*blatx-fc*clatx)/alatx
      ia=getcellcoord(fa,ncelldiva)
      ib=getcellcoord(fb,ncelldivb)
      ic=getcellcoord(fc,ncelldivc)

      atomcellno=1+ia+ncelldivb*(ib+ncelldivc*ic)

      return
      end function atomcellno

c************************************************************************************************         

      function cellcoord(x,ncelldiv)

c Returns the appropriate cell coordinate, apply periodic boundary conditions

      implicit none

      integer ix,ncelldiv,cellcoord
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

      subroutine getcellneighbourlist(ifound,ixcell,iycell,izcell,ncellrange
     1,itype,ntyperange)

      implicit none

c Gets the list of neighbouring cells relative to the one at (0,0,0) which satisfy the required
c distance criterion

      integer ncellrange
      integer ixcell(*),iycell(*),izcell(*),ix,iy,iz
      integer ixsq,iysq,izsq,ncellrangesq
      integer ifound

c itype stores whether a cell is within the inner or outer cell range. Atoms in cells
c outside the inner cell range will treated as one atom at the centre of coordinates of the cell

      integer itype(*),ntyperange,ntyperangesq

c set a maximum separation on the cells. This is based on an atom in the corner of a cubic cell
c being able to see an atom in the opposite corner of a neighbouring cubic cell

      ncellrangesq=ncellrange**2
      ntyperangesq=ntyperange*ntyperange

      ifound=0

c Step through all the possible local cells including the central one,
c eliminating any cells too far away to be useful.

      do iz=-ncellrange,ncellrange
         izsq=0
c Get the shortest possible distance between atoms in these two cells
         if(iz.ne.0) izsq=iz*iz+1-2*abs(iz)
         do iy=-ncellrange,ncellrange
            iysq=izsq
            if(iy.ne.0) iysq=iysq+iy*iy+1-2*abs(iy)
            do ix=-ncellrange,ncellrange
               ixsq=iysq
               if(ix.ne.0) ixsq=ixsq+ix*ix+1-2*abs(ix)
c Eliminate cells that are outside the range
               if(ixsq.le.ncellrangesq) then
                  ifound=ifound+1
                  ixcell(ifound)=ix
                  iycell(ifound)=iy
                  izcell(ifound)=iz
                  if(ixsq.gt.ntyperangesq) then
                     itype(ifound)=1
                  else
                     itype(ifound)=0
                  endif
               endif
            end do
         end do
      end do

      return
      end subroutine getcellneighbourlist

c************************************************************************************************         

      subroutine getcellneighbours(cellnumber,ifound,nfoundcell,ixcell,iycell,izcell
     1,cellneighbours,ncelldiv,xoffset,yoffset,zoffset,itype,itypecell)

      implicit none

c Gets the neighbours of cellnumber, assuming periodic boundary conditions
c The cells are assumed to be numbered first along a, then along b, finally along c

      integer cellnumber,cellnumber1,cellneighbours(*)
      integer nfoundcell,ixcell(*),iycell(*),izcell(*),ncelldiv
      integer ixcentre,iycentre,izcentre,i
      integer xsub,ysub,zsub
      integer ixneighb,iyneighb,izneighb,ncelldivsq
      integer cellnumberneighb
      integer ifound,itest

c itype stores whether a cell is within the inner or outer cell range. itypecell stores
c this for the particular set of current neighbours

      integer itype(*),itypecell(*)

      real xoffset(*),yoffset(*),zoffset(*)

      logical found

      ifound=0
      ncelldivsq=ncelldiv*ncelldiv

c Convert cellnumber to the corresponding x,y,z coordinates 

      cellnumber1=cellnumber-1
      izcentre=int(cellnumber1/ncelldivsq)
      iycentre=int((cellnumber1-izcentre*ncelldivsq)/ncelldiv)
      ixcentre=cellnumber1-ncelldiv*(iycentre+ncelldiv*izcentre)
c      write(6,*) cellnumber1,ixcentre,iycentre,izcentre

c Step through all the cells including the central one, and applying periodic boundary conditions
c and eliminating any cells too far away to be useful.

      i=0
      do while(i.lt.nfoundcell)

         i=i+1
         izneighb=izcell(i)+izcentre
         call getperiodiccell(izneighb,ncelldiv,zsub)
         iyneighb=iycell(i)+iycentre
         call getperiodiccell(iyneighb,ncelldiv,ysub)
         ixneighb=ixcell(i)+ixcentre
         call getperiodiccell(ixneighb,ncelldiv,xsub)
         cellnumberneighb=1+ixneighb+ncelldiv*(iyneighb+ncelldiv*izneighb)
c Check that this number has not already been included in the list
         itest=0
         found=.false.
         do while (itest.lt.ifound.and..not.found)
            itest=itest+1
            found=cellnumberneighb.eq.cellneighbours(itest)
         end do
         if (.not.found) then
            ifound=ifound+1
            cellneighbours(ifound)=cellnumberneighb
            xoffset(ifound)=xsub
            yoffset(ifound)=ysub
            zoffset(ifound)=zsub
            itypecell(ifound)=itype(i)
         endif
      end do

      return
      end subroutine getcellneighbours

c************************************************************************************************         

      subroutine getperiodiccell(neighb,ncelldiv,sub)

      implicit none
      integer neighb,ncelldiv,add,sub

c Calculate the bit to add to neighb to get the round off to work correctly for negative numbers

      add=min(0,sign(1,neighb))
      sub=(neighb/ncelldiv)+add
      neighb=neighb-ncelldiv*sub

      return
      end subroutine getperiodiccell

c************************************************************************************************         

      subroutine getdummycoordinates(cellnumber,cellnumberold
     1,x,y,z
     1,nfoundcell,ixcell,iycell,izcell
     1,nincell,xcell,ycell,zcell,incelllist
     1,ndummy,xdummy,ydummy,zdummy,indummylist
     1,cellneighbours,ncellatomnumbers,cellatomnumber
     2,ncelldiv,ncellatomnumberdiv,xoffset,yoffset,zoffset
     1,itype,itypecell,makegr)

      implicit none

c Input variables

      integer cellnumber,nincell,ndummy,cellneighbours(*)
     1,cellatomnumber(*),incelllist(*),indummylist(*)
     1,ncellatomnumbers(*)
     1,ncelldiv,ncellatomnumberdiv
     1,nfoundcell,ixcell(*),iycell(*),izcell(*),itype(*),itypecell(*)


      real x(*),y(*),z(*),xdummy(*),ydummy(*),zdummy(*),ran1
     1,xoffset(*),yoffset(*),zoffset(*),xcell(*),ycell(*),zcell(*)
c     1,xav,yav,zav

      logical makegr !Controls output if g(r) is being calculated

c Internal variables

      integer j,k,l,lmax,add,cellnumberk,cellnumberold
     1,nlcell

c Only get the coordinates if the cellnumber is different from the last time this routine
c was called

      if(cellnumber.ne.cellnumberold) then

         ndummy=0
         nincell=0

c Get the neighbours of this cell

         call getcellneighbours(cellnumber,nlcell,nfoundcell,ixcell,iycell,izcell
     1,cellneighbours,ncelldiv,xoffset,yoffset,zoffset,itype,itypecell)

c Step through the neighbours and add all the atoms to the neighbourlist

         k=0
         do while(k.lt.nlcell)
            k=k+1
            cellnumberk=cellneighbours(k)
            add=(cellnumberk-1)*ncellatomnumberdiv
            l=add
            lmax=ncellatomnumbers(cellnumberk)+add
            if(cellnumberk.ne.cellnumber) then

c Get the coordinates of atoms not in the same cell as the origin atom

               if(itypecell(k).eq.0.or.makegr) then

c For type 0 cells we gather the coordinates of the individual atoms

                  do while(l.lt.lmax)
                     l=l+1
                     j=cellatomnumber(l)
                     ndummy=ndummy+1
                     xdummy(ndummy)=x(j)+xoffset(k)
                     ydummy(ndummy)=y(j)+yoffset(k)
                     zdummy(ndummy)=z(j)+zoffset(k)
                     indummylist(ndummy)=1
                  end do
               else

c For type 1 cells we gather all the atoms into a single coordinate

                  if(ncellatomnumbers(cellnumberk).gt.0) then
c                     xav=0.0
c                     yav=0.0
c                     zav=0.0
c                     do while(l.lt.lmax)
c                        l=l+1
c                        j=cellatomnumber(l)
c                        xav=xav+x(j)
c                        yav=yav+y(j)
c                        zav=zav+z(j)
c                     end do
c Choose coordinate at random from within the list of atom numbers
                     l=add+ran1()*ncellatomnumbers(cellnumberk)+1
                     j=cellatomnumber(l)                  
                     ndummy=ndummy+1
c                     xdummy(ndummy)=xav/ncellatomnumbers(cellnumberk)+xoffset(k)
c                     ydummy(ndummy)=yav/ncellatomnumbers(cellnumberk)+yoffset(k)
c                     zdummy(ndummy)=zav/ncellatomnumbers(cellnumberk)+zoffset(k)
                     xdummy(ndummy)=x(j)+xoffset(k)
                     ydummy(ndummy)=y(j)+yoffset(k)
                     zdummy(ndummy)=z(j)+zoffset(k)
                     indummylist(ndummy)=ncellatomnumbers(cellnumberk)
                  endif
               endif
            else
               do while(l.lt.lmax)
                  l=l+1
                  j=cellatomnumber(l)
                  nincell=nincell+1
                  xcell(nincell)=x(j)
                  ycell(nincell)=y(j)
                  zcell(nincell)=z(j)
                  incelllist(nincell)=j
               end do
            endif
         end do
         cellnumberold=cellnumber
      endif

      return
      end subroutine getdummycoordinates

c************************************************************************************************         

      subroutine get_referencenumberlist(i,xorigin,yorigin,zorigin
     1,nincell,incell,xcell,ycell,zcell
     1,ndummy,xdummy,ydummy,zdummy,indummylist
     1,rsq,rcutoffsq,rsqstep,maprsq,nr
     1,atomnumberlist,referencenumberlist,weightlist
     1,natomnumberlist,natomnumberlistneighb)

      implicit none

c Input variables

      integer i,nincell,incell,ndummy,nr,maprsq(*),indummylist(*)
     1,atomnumberlist(*),referencenumberlist(*),weightlist(*)
     1,natomnumberlist,natomnumberlistneighb

      real xcell(*),ycell(*),zcell(*),xdummy(*),ydummy(*),zdummy(*)
     1,rsq(*),rcutoffsq,rsqstep
     1,xorigin,yorigin,zorigin

c Internal variables

      integer j,iref,ifind2

      real xdif,ydif,zdif,rdifsq

c i is the atom at the origin

      natomnumberlist=0

c Step through the neighbouring cells and add all the atoms to the neighbourlist if within range

      j=0
      do while(j.lt.ndummy)

         j=j+1
         xdif=xorigin-xdummy(j)
         ydif=yorigin-ydummy(j)
         zdif=zorigin-zdummy(j)
         rdifsq=(xdif*xdif+ydif*ydif+zdif*zdif)
         if(rdifsq.lt.rcutoffsq) then
            iref=ifind2(rdifsq,rsq,rsqstep,maprsq,nr)
            natomnumberlist=natomnumberlist+1
            atomnumberlist(natomnumberlist)=j
            referencenumberlist(natomnumberlist)=iref
            weightlist(natomnumberlist)=indummylist(j)
         endif
      end do
      natomnumberlistneighb=natomnumberlist

c Step through the current cell, picking off all the neighbours of atom i, but not including atom i

      j=0
      do while(j.lt.nincell)

         j=j+1
         if(j.ne.i) then
            xdif=xorigin-xcell(j)
            ydif=yorigin-ycell(j)
            zdif=zorigin-zcell(j)
            rdifsq=(xdif*xdif+ydif*ydif+zdif*zdif)
            if(rdifsq.lt.rcutoffsq) then
               iref=ifind2(rdifsq,rsq,rsqstep,maprsq,nr)
               natomnumberlist=natomnumberlist+1
               atomnumberlist(natomnumberlist)=j
               referencenumberlist(natomnumberlist)=iref
            endif
         else
c Save the index of the atom being moved
            incell=j
         endif
      end do

      return
      end subroutine get_referencenumberlist

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


      end module cellneighbourlist      
