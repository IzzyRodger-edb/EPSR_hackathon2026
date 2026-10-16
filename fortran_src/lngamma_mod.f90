      module lngamma_mod

      integer                                   :: mgamma=0
      real(kind=8), dimension(:), allocatable   :: gamln

      contains

      subroutine fctrln(nterms)
      integer nterms,ncalc,i
      real(kind=8) g,a
      if(allocated(gamln)) deallocate(gamln)
      mgamma=nterms+1
      allocate(gamln(mgamma))
      ncalc=mgamma
      g=0.0
      a=0.0
      do i=1,ncalc
            gamln(i)=g
            a=a+1.0
            g=g+dlog(a)
      end do
      RETURN
      END subroutine fctrln

      function lngamma(ic)
      real(kind=8)                           :: lngamma
      integer                                :: ic
      if(ic.gt.mgamma.or..not.allocated(gamln)) then
         call fctrln(ic)
      endif
      lngamma=gamln(ic+1)
      return
      end function lngamma

      end module lngamma_mod


