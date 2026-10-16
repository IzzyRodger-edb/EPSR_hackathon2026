   module pfroutines

      real, dimension(:), allocatable :: rratio,intenprobfunction

   contains

      subroutine getpfprobability(rfract)

      real, intent(in) :: rfract
      real :: probfunction
      real :: rstep,rsigma,gammacons,fullwidth
      integer ir,nr

      nr=50 !Must be greater than 1
      if(.not.allocated(intenprobfunction).or.size(intenprobfunction,1).lt.nr+1) then
         if(allocated(intenprobfunction)) deallocate(rratio,intenprobfunction)
         allocate(rratio(0:nr),intenprobfunction(0:nr))
      end if
      rsigma=max((1.0/(rfract*rfract)-1.0),0.0) !Keep the widths sensible: rfract cannot be > 1.
      gammacons=lngamma(rsigma) ! log of gamma(rsigma+1)
! Set up the probability function - assumed Poisson in form for r > 0.
      fullwidth=rsigma+1.0+10.0*sqrt(rsigma+1.0) !Full width of Poisson function in normalised units = position of max + 10 stds.
      rstep=fullwidth/real(nr)
      write(6,*) 'rsigma, gammacons, fullwidth, rstep: ', rsigma,gammacons,fullwidth,rstep
      rratio(0)=0.0
      intenprobfunction(0)=0.0
      do ir=1,nr
         r=real(ir)*rstep
         arg=rsigma*log(r)-r-gammacons
         if(abs(arg).lt.30.0) then
            probfunction=exp(arg)
         else
            probfunction=0.0
         end if
         rratio(ir)=r/(rsigma+1.0)
         intenprobfunction(ir)=intenprobfunction(ir-1)+probfunction*rstep!/(rratio(ir)**3)
      end do
! Final renormalisation to correct round-off errors
      do ir=1,nr
         intenprobfunction(ir)=intenprobfunction(ir)/intenprobfunction(nr)
      end do
      end subroutine getpfprobability

      function lngamma(argin) !log of gamma(argin+1.0)
      real argin,lngamma
      real*8 arg,gammaval,sumlngamma
      integer narg,i
      lngamma=0.0
      if(argin.lt.0.0) return !return 0 if illegal argument
      narg=int(argin) ! Integer part of arg
      arg=argin-real(narg)+1.0
      gammaval=gamma(arg)
      sumlngamma=log(gammaval)
!      write(6,*) 'lngamma> ',argin,narg,arg,gammaval,sumlngamma
      if(narg.gt.0) then
         do i=1,narg
            sumlngamma=sumlngamma+log(arg)
            arg=arg+1.0
         end do
      end if
!      write(6,*) 'lngamma> ',argin,narg,arg,gammaval,sumlngamma
      lngamma=real(sumlngamma)
      return
      end function lngamma

      function getrratio()

      real getrratio
      real ranrratio,grad,ran1
      logical found
      integer i

      ranrratio=ran1()
      i=0
      found=.false.
      do while(i.lt.size(rratio)-1.and..not.found)
         i=i+1
         found=ranrratio.lt.intenprobfunction(i)
      end do
      diff=intenprobfunction(i)-intenprobfunction(i-1)
      grad=0.0
      if(diff.gt.0.0) then
         grad=(rratio(i)-rratio(i-1))/diff
      end if
      getrratio=rratio(i-1)+(ranrratio-intenprobfunction(i-1))*grad

      end function getrratio

   end module pfroutines





            

