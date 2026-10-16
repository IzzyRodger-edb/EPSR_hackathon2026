      module lngamma_mod
      integer, parameter :: mgamma=600
      real(kind=8), dimension(mgamma) :: gamln

      contains

      SUBROUTINE FCTRLN(nterms)
      real(kind=8) :: g,a
      nmax=mgamma
      ncalc=nterms+1
      if(ncalc.gt.nmax) then
            ncalc=nmax
            nterms=ncalc-1
      endif
      G=0.0
      GAMLN(1)=0.0
      A=0.0
      DO I=1,ncalc
            GAMLN(I)=G
            A=A+1.0
            g=g+dLOG(A)
      end do
      RETURN
      END subroutine fctrln

      end module lngamma_mod


