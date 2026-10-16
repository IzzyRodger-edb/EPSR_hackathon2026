      subroutine smoosharm(NQMIN,NQMAX,NSMOO,DIFF)
      DIMENSION DIFF(*)
C
C 3-POINT SMOOTHING keeping end points fixed
C
      n1=nqmin+1
      n2=nqmax-1
      diff(nqmax)=0.0 !Last value is set to zero exactly
      if((n2-n1).lt.2.or.nsmoo.eq.0) return
      DO IT=1,NSMOO
         sum1=0.25*(DIFF(NQMIN)+DIFF(N1))
         SM=DIFF(N1)
         DO I=N1,N2
            J=I+1
            SN=DIFF(J)
            sum2=0.25*(sm+sn)
            DIFF(I)=sum1+sum2
            SM=SN
            sum1=sum2
         end do
      end do
      RETURN
      END

      function quadinterp2d(grw,rrat,thrat,nro,nrp,nrm,ntho,nthp
     1,nthm,ntdim,nrdim)
      real, dimension(ntdim,nrdim), intent(in) :: grw
c
c quadratic interpolation along r
c
      delta=rrat
      dp=1.0+delta
      dm=1.0-delta
      dpm=dm*dp
      delta2=0.5*delta
      am=delta2*(dp*grw(nthm,nrp)-dm*grw(nthm,nrm))+dpm*grw(nthm,nro)
      ao=delta2*(dp*grw(ntho,nrp)-dm*grw(ntho,nrm))+dpm*grw(ntho,nro)
      ap=delta2*(dp*grw(nthp,nrp)-dm*grw(nthp,nrm))+dpm*grw(nthp,nro)
c
c now interpolate in theta
c
      delta=thrat
      dp=1.0+delta
      dm=1.0-delta
      dpm=dm*dp
      delta2=0.5*delta
      quadinterp2d=delta2*(dp*ap-dm*am)+dpm*ao
      return
      end

      function quadinterp(grw,rrat,thrat,phrat,nro,nrp,nrm,ntho,nthp
     1,nthm,npho,nphp,nphm,ntdim,npdim,nrdim)
      real, dimension(ntdim,npdim,nrdim), intent(in) :: grw
c
c quadratic interpolation along r
c
      delta=rrat
      dp=1.0+delta
      dm=1.0-delta
      dpm=dm*dp
      delta2=0.5*delta
      amm=delta2*(dp*grw(nthm,nphm,nrp)-dm*grw(nthm,nphm,nrm))
     1+dpm*grw(nthm,nphm,nro)
      amo=delta2*(dp*grw(nthm,npho,nrp)-dm*grw(nthm,npho,nrm))
     1+dpm*grw(nthm,npho,nro)
      amp=delta2*(dp*grw(nthm,nphp,nrp)-dm*grw(nthm,nphp,nrm))
     1+dpm*grw(nthm,nphp,nro)
      aom=delta2*(dp*grw(ntho,nphm,nrp)-dm*grw(ntho,nphm,nrm))
     1+dpm*grw(ntho,nphm,nro)
      aoo=delta2*(dp*grw(ntho,npho,nrp)-dm*grw(ntho,npho,nrm))
     1+dpm*grw(ntho,npho,nro)
      aop=delta2*(dp*grw(ntho,nphp,nrp)-dm*grw(ntho,nphp,nrm))
     1+dpm*grw(ntho,nphp,nro)
      apm=delta2*(dp*grw(nthp,nphm,nrp)-dm*grw(nthp,nphm,nrm))
     1+dpm*grw(nthp,nphm,nro)
      apo=delta2*(dp*grw(nthp,npho,nrp)-dm*grw(nthp,npho,nrm))
     1+dpm*grw(nthp,npho,nro)
      app=delta2*(dp*grw(nthp,nphp,nrp)-dm*grw(nthp,nphp,nrm))
     1+dpm*grw(nthp,nphp,nro)
!      if(nro.eq.51) then
!         write(6,*) nrm,nro,nrp,ntho,npho,delta2,dp,dm,dpm
!         write(6,*) aoo,rrat,grw(ntho,npho,nrm),grw(ntho,npho,nro),grw(ntho,npho,nrp)
!      end if
c
c now interpolate in theta
c
      delta=thrat
      dp=1.0+delta
      dm=1.0-delta
      dpm=dm*dp
      delta2=0.5*delta
      am=delta2*(dp*apm-dm*amm)+dpm*aom
      ao=delta2*(dp*apo-dm*amo)+dpm*aoo
      ap=delta2*(dp*app-dm*amp)+dpm*aop
c
c now interpolate in phi
c
      delta=phrat
      dp=1.0+delta
      dm=1.0-delta
      dpm=dm*dp
      delta2=0.5*delta
      quadinterp=delta2*(dp*ap-dm*am)+dpm*ao
!      if((abs(thrat).gt.0.5.and.ntho.gt.2.and.ntho.lt.18).or.abs(phrat).gt.0.5) then
!      if(nro.eq.51) then
!      if(abs(quadinterp-grw(ntho,npho,nro)).gt.2.0) then
!         write(6,*) nrm,nro,nrp,nthm,ntho,nthp,nphm,npho,nphp
!         write(6,*) rrat,thrat,phrat,quadinterp
!         write(6,*) aoo,grw(ntho,npho,nrm),grw(ntho,npho,nro),grw(ntho,npho,nrp)
!      end if
!      quadinterp=grw(ntho,npho,nro)
      return
      end

      subroutine sumoverm2(l1,l2,l,n1,m2f,m2l,n2,phil,ctl,phi,cz,chi,phs,sumc,sums)

      integer l1,l2,l,n1,m2f,m2l,n2,m2
      real cgcoef,dlm,phil,ctl,phi,cz,chi,phs,sumc,sums,phisum,cg,dl2m2n2,dlm0
! For negative n1, we need the complex conjugate (See Gray and Gubbins, eq 2.25). Sharm gives us only +ve n1, but all n2
      phssine=-real(sign(1,n1))
      m1=n1 !Imposed by requirement that orientation of molecule at origin is (000)
      do m2=m2f,m2l
         m=m1+m2
         cg=cgcoef(l1,l2,l,m1,m2)
         dl2m2n2=dlm(l2,m2,n2,cz)
         dlm0=dlm(l,m,0,ctl)
         phisum=m2*phi+n2*chi-m*phil !Reverse signs because of complex conjugate
         prod=phs*cg*dl2m2n2*dlm0
         sumc=sumc+prod*cosd(phisum)
         sums=sums+phssine*prod*sind(phisum)
!         write(6,'(7(i4),8(1x,e13.6))') l1,l2,l,n1,n2,m2,m,phs,phssine,cg,dl2m2n2,dlm0,phisum,prod 
      end do
      return
      end

      function suminterp(grw,rrat,thrat,nro,nrp,nrm,ntho,nthp
     1,nthm,np0,nphp,nphm,ntdim,npdim,nrdim)
      real, dimension(ntdim,npdim,nrdim), intent(in) :: grw
      real :: suminterp

      sum=0.0
      nph=nphp-nphm+1
      if(nph.gt.np0) then
         nphm=1
         nphp=np0
         nph=np0
      end if
      do npht=nphm,nphp
         npho=npht      
! Apply cyclic boundary conditions on phi
         do while (npho.gt.np0)
            npho=npho-np0
         end do
         do while (npho.lt.1)
            npho=npho+np0
         end do
c
c quadratic interpolation along r
c
         delta=rrat
         dp=1.0+delta
         dm=1.0-delta
         dpm=dm*dp
         delta2=0.5*delta
         amo=delta2*(dp*grw(nthm,npho,nrp)-dm*grw(nthm,npho,nrm))+dpm*grw(nthm,npho,nro)
         aoo=delta2*(dp*grw(ntho,npho,nrp)-dm*grw(ntho,npho,nrm))+dpm*grw(ntho,npho,nro)
         apo=delta2*(dp*grw(nthp,npho,nrp)-dm*grw(nthp,npho,nrm))+dpm*grw(nthp,npho,nro)
c
c now interpolate in theta
c
         delta=thrat
         dp=1.0+delta
         dm=1.0-delta
         dpm=dm*dp
         delta2=0.5*delta
         ao=delta2*(dp*apo-dm*amo)+dpm*aoo
         sum=sum+ao
      end do
      suminterp=sum/real(nph)
      return
      end


