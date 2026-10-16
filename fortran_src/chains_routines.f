      function match(string1,string2)

c Determines whether string1 is the same as string2 up any * that may be present in
c string2. If no * is present, the whole strings are matched.

      logical match
      character*3 string1
      character*4 string2

      nstar=index(string2,'*')

      if(nstar.gt.0) then

         nsearch=nstar-1
         if(nsearch.gt.0) then
            match=string1(1:nsearch).eq.string2(1:nsearch)
         else
            match=.true.
         endif
         return

      else

         match=string1.eq.string2
         return

      endif 

      end

      subroutine shortest_path(im1,nz,nmol,nneighb,neighb,rneighb1,rneighb2
     *,spvaluet,nspvaluet,mmol,mbond)

c Calculates the shortest (topological) path to im1 from all other molecules
c in the box, out to a 'diameter' of nz

      integer*4 nneighb(mmol),neighb(mmol,mbond)
c Variables for setting up the shortest path to the molecule at the origin
      integer rneighb1(mmol),rneighb2(mmol),spvaluet(mmol),nspvaluet(nmol,nz)
      integer, dimension(:), allocatable :: sourceatom1,sourceatom2
      integer intest,nneighb1,nneighb2
      logical found

      if(allocated(sourceatom1)) deallocate(sourceatom1,sourceatom2)
      allocate(sourceatom1(nmol),sourceatom2(nmol))

c Initialise the shortest topological path to other molecules larger than the number of
c molecules in the box.

      do im=1,nmol
         spvaluet(im)=nmol+1
         do icount=1,nz
            nspvaluet(im,icount)=0
         end do
      end do
c Set the shortest topological distances to other molecules in the vicinity of im1
      nneighb1=1
      rneighb1(nneighb1)=im1
      sourceatom1(nneighb1)=im1
c Distance to molecule at origin is zero
      spvaluet(im1)=0
      icount=0
      ineighb1=0
c limit the maximum value of icount needed
      do while(icount.lt.nz)
         icount=icount+1
c form a list of the neighbours of all the molecules listed in rneighb1
         nneighb2=0
         ineighb1save=ineighb1 !Need to save the current molecule in list to ensure we do not refer back to it
         do while(ineighb1.lt.nneighb1)
            ineighb1=ineighb1+1
            in=0
            do while(in.lt.nneighb(rneighb1(ineighb1)))
               in=in+1
               im2=neighb(rneighb1(ineighb1),in)
! Need to ensure im2 does not refer back to the molecule that precedes rneighb1(ineighb1) in the chain
               im2precede=sourceatom1(ineighb1)
               if(im2.ne.im2precede) then
!               write(6,*) icount,im1,im2,sourceatom1(ineighb1)
c set the value of shortest path if it is smaller than the current value
               if(icount.lt.spvaluet(im2)) then
                  spvaluet(im2)=icount
               end if
! Count number of times this chain length occurs.
               nspvaluet(im2,icount)=nspvaluet(im2,icount)+1
! Check whether molecule im2 is already in the list
               found=.false.
               intest=0
               do while(.not.found.and.intest.lt.nneighb1)
                  intest=intest+1
                  found=im2.eq.rneighb1(intest)
               end do
! Also check this molecule does not occur in the new list either
               intest=0
               do while(.not.found.and.intest.lt.nneighb2)
                  intest=intest+1
                  found=im2.eq.rneighb2(intest)
               end do
! If not found, set the new list
               if(.not.found.and.nneighb2.lt.mmol) then
                  nneighb2=nneighb2+1
                  rneighb2(nneighb2)=im2
                  sourceatom2(nneighb2)=rneighb1(ineighb1) !Save the source atom for the current set of neighbours
               end if
               end if
            end do
         end do
!Save the new values
         in=0
         do while(in.lt.nneighb2.and.nneighb1.lt.mmol)
            in=in+1
            nneighb1=nneighb1+1
            rneighb1(nneighb1)=rneighb2(in)
            sourceatom1(nneighb1)=sourceatom2(in)
         end do
!         write(6,*)
!         write(6,*) 'shortest_path> ',im1,icount,nneighb1,(rneighb1(in),spvaluet(rneighb1(in)),in=1,nneighb1)
      end do
!      if(im1.eq.1) then
!         write(6,*)
!         write(6,*) 'shortest_path> ',im1,icount,nneighb1
!         do in=1,nneighb1
!            write(6,*) rneighb1(in),spvaluet(rneighb1(in)
!         end do
!      end if
!      stop
      return

      end


