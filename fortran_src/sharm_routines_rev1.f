      module sharm_mod

!Specification of l1, l2 and l values
      integer                                   :: nll1v,nll2v,nllv,n1step,n2step
      integer, dimension(:), allocatable        :: ll1v,ll2v,llv
!Arrays of l1,l2,l,n1,n2 values
      integer                                   :: nlmax,ndist,ldistinct
      integer, dimension(:), allocatable        :: ll1,ll2,ll,n1,n2
      real, dimension(:), allocatable           :: alsec
!Arrays of euler angles and coefficients
      integer                                   :: mdlmn,mz,mdist
      integer, dimension(:), allocatable        :: nm,il1rf,il2rf,ilrf,illrf !First or second molecule, reference integers
      real, dimension(:), allocatable           :: phim,cthetm,chim,xcent,ycent,zcent
      real, dimension(:), allocatable           :: philab
      real, dimension(:,:), allocatable         :: dlmnc
      real, dimension(:,:,:), allocatable       :: dlmnmc,chiphi
!Array of Clebsch-Gordan coefficients
      integer                                   :: mcg
      real, dimension(:), allocatable           :: termcg
!Arrays of angles and coefficients for testing
      real, dimension(:), allocatable           :: phimt1,cthetmt1,chimt1,phimt2,cthetmt2,chimt2
      real, dimension(:,:,:), allocatable       :: dlmnmct,chiphit
!Input and output arrays
      real, dimension(:), allocatable           :: z
      real, dimension(:,:), allocatable         :: gofz,gofzs,gzt,gzt2
      character(len=29), dimension(:), allocatable            :: titledat
!Logical variables to perform a test calculation and to output average or single molecules for plotting purposes
      logical                                   :: performtest,outputaveragemol

      end module sharm_mod

      subroutine get_l_and_n_values()

      use sharm_mod

c set up l and n values
! First count the number of terms
      nterm=0
c l1 values
      do ill1=1,nll1v
         il1=ll1v(ill1)
c l2 values
         do ill2=1,nll2v
            il2=ll2v(ill2)
            if((ldistinct.eq.0.and.il2.le.il1).or.ldistinct.eq.1) then
c l values - use closure relation on CG coefficients to determine range
               jlmax=iabs(il1+il2)
               jlmin=iabs(il1-il2)
               do ill=1,nllv
                  il=llv(ill)
                  if(il.ge.jlmin.and.il.le.jlmax) then
c n1 values
                     if(n1step.gt.0) then
                        n1max=il1+1
                        n1s=n1step
                     else
                        n1max=1
                        n1s=1
                     endif
                     do jn1=1,n1max,n1s
                        in1=jn1-1
c n2 values
                        if(n2step.gt.0) then
                           nts=il2/n2step
                           if(in1.gt.0) then !Avoid redundant coefficients
                              n2min=-nts*n2step
                           else
                              n2min=0
                           end if
                           n2max=nts*n2step
                           n2s=n2step
                        else
                           n2min=0
                           n2max=0
                           n2s=1
                        endif
                        do in2=n2min,n2max,n2s
                           nterm=nterm+1
                        end do
                     end do
                  endif
               end do
            endif
         end do
      end do
      ndist=nterm
      nlmax=nterm
      if(allocated(ll1)) deallocate(ll1,ll2,ll,n1,n2,il1rf,il2rf,ilrf,illrf,alsec)
      allocate(ll1(nlmax),ll2(nlmax),ll(nlmax),n1(nlmax),n2(nlmax),il1rf(nlmax),il2rf(nlmax),ilrf(nlmax),illrf(nlmax),alsec(nlmax))
      nterm=0
! save them!
c l1 values
      do ill1=1,nll1v
         il1=ll1v(ill1)
c l2 values
         do ill2=1,nll2v
            il2=ll2v(ill2)
            if((ldistinct.eq.0.and.il2.le.il1).or.ldistinct.eq.1) then
c l values - use closure relation on CG coefficients to determine range
               jlmax=iabs(il1+il2)
               jlmin=iabs(il1-il2)
               do ill=1,nllv
                  il=llv(ill)
                  if(il.ge.jlmin.and.il.le.jlmax) then
c n1 values
                     if(n1step.gt.0) then
                        n1max=il1+1
                        n1s=n1step
                     else
                        n1max=1
                        n1s=1
                     endif
                     do jn1=1,n1max,n1s
                        in1=jn1-1
! n2 values
                        if(n2step.gt.0) then
                           nts=il2/n2step
                           if(in1.gt.0) then !Avoid redundant coefficients
                              n2min=-nts*n2step
                           else
                              n2min=0
                           end if
                           n2max=nts*n2step
                           n2s=n2step
                        else
                           n2min=0
                           n2max=0
                           n2s=1
                        endif
                        do in2=n2min,n2max,n2s
                           if(nterm.lt.nlmax) then
                              nterm=nterm+1
                              ll1(nterm)=il1
                              ll2(nterm)=il2
                              ll(nterm)=il
                              n1(nterm)=in1
                              n2(nterm)=in2
                           endif
                        end do
                     end do
                  endif
               end do
            endif
         end do
      end do
      do i=1,ndist
         write(6,46) i,ll1(i),ll2(i),ll(i),n1(i),n2(i)
      end do
46    format(6(1x,i5))
! Form Clebsch-Gordan coefficients ln+ gamma values
      call form_cgcoeffs()
      return
      end

      subroutine form_cgcoeffs()

      use sharm_mod
      use lngamma_mod


      include 'epsr_dimension.inc'
      include 'sharm.inc'
      include 'sharm_sdf.inc'

      integer getrefint,getrefintnon
      integer il1,il2,il,in1,in2,m1,m2,m,m1min,m1max,m2min,m2max,mmin,mmax,mgammat

      do id=1,ndist
            alsec(id)=0.0
      end do
c
c set up the required CG coefficients and determine maximum number of ln gamma values required
c
! First count how many needed
      ic=0
      mgammat=0
      do id=1,ndist
         il1=ll1(id)
         in1=n1(id)
         il2=ll2(id)
         in2=n2(id)
         il=ll(id)
c
c initialize reference integer for molecule 1
c
         il1rf(id)=1+getrefint(il1,in1,ll1v,nll1v,n1step)+il1
c
c initialize reference integer for molecule 2
c
         il2rf(id)=1+getrefint(il2,in2,ll2v,nll2v,n2step)+il2
c
c initialize reference integer for laboratory function
c
         ilrf(id)=1+getrefintnon(il,llv,nllv)+il
         m1max=il1
         m1min=-il1
         do m1=m1min,m1max
            mgammat=max(mgammat,maxdlm(il1,m1,in1)) !Include the dlm's in the calculation
            mgammat=max(mgammat,maxdlm(il1,m1,-in1)) !Include the dlm's in the calculation
            m2max=min(il2,il-m1)
            m2min=max(-il2,-il-m1)
            do m2=m2min,m2max
               ic=ic+1
               m=m1+m2
               mgammat=max(mgammat,maxcgcoef(il1,il2,il,m1,m2))
               mgammat=max(mgammat,maxdlm(il2,m2,in2)) !Include the dlm's in the calculation
               mgammat=max(mgammat,maxdlm(il2,m2,-in2)) !Include the dlm's in the calculation
               mgammat=max(mgammat,maxdlm(il,m,0)) !Include the dlm's in the calculation
            end do
         end do
!         mmin=-il
!         mmax=il
!         do m=mmin,mmax
!            mgammat=max(mgammat,maxdlm(il,m,0))
!         end do
      end do
      write(6,*) 'form_cgcoeffs> Creating ',mgammat,' log gamma values'
      call fctrln(mgammat)
      mcg=ic
      if(allocated(termcg)) deallocate(termcg)
      allocate(termcg(mcg))
      ic=0
      do id=1,ndist
         il1=ll1(id)
         in1=n1(id)
         il2=ll2(id)
         in2=n2(id)
         il=ll(id)
         m1max=il1
         m1min=-il1
         do m1=m1min,m1max
            m2max=min(il2,il-m1)
            m2min=max(-il2,-il-m1)
            do m2=m2min,m2max
               ic=ic+1
               termcg(ic)=cgcoef(il1,il2,il,m1,m2)
               alsec(id)=alsec(id)+termcg(ic)*termcg(ic)
            end do
         end do
!         write(6,*) id,il1,in1,il1rf(id),il2,in2,il2rf(id),il,ilrf(id)
      end do
      write(6,*) 'form_cgcoeffs> ',mcg,' Clebsch-Gordan coefficients created'
      return
      end

      function getrefint(lreq,nreq,llv,nllv,nstep)
! Gets the reference integers for the Dlmn for specfied l and n values

      integer :: getrefint,lreq,nreq,nllv,nstep
      integer, dimension(*) :: llv
      integer :: id,sum,l,nts,ntmax

      sum=0
! Get the number of terms up to, but not including the current value of id
      id=1
      do while(id.lt.nllv.and.llv(id).lt.lreq)
         l=llv(id)
         if(nstep.gt.0) then
            nts=l/nstep+1
         else
            nts=1
         endif
! For each l and n value there will be 2l+1 m values. However since -n is treated at same time as +n, we treat -n same as +n
         sum=sum+nts*(2*l+1)
         id=id+1
      end do
! Add on the n values for the current id
      l=llv(id)
      ntmax=abs(nreq)
      if(nstep.gt.0) then
         nts=ntmax/nstep
         sum=sum+nts*(2*l+1)
      end if
      getrefint=sum
      return
      end

      function getrefintnon(lreq,llv,nllv)
! Gets the reference integers for the Dlmn for specfied l with n values set to zero
      integer :: getrefintnon,lreq
      integer, dimension(*) :: llv
      integer :: id,sum,l
      sum=0
! Get the number of terms up to, but not including the current value of id
      id=1
      do while(id.lt.nllv.and.llv(id).lt.lreq)
         l=llv(id)
! For each l and n value there will be 2l+1 m values. However since -n is treated at same time as +n, we treat -n same as +n
         sum=sum+2*l+1
         id=id+1
      end do
      getrefintnon=sum
      return
      end

      subroutine form_eulerangles(imol1,xt,yt,zt,naxis1,naxis2,ncz,ncy,lztyp,lytyp,alab,blab,clab,phim,cthetm,chim)
! Vectors alab, blab and clab must be defined before calling this routine

      integer                 :: imol1,ncz,ncy
      real, dimension(*)      :: xt,yt,zt
      real                    :: phim,cthetm,chim
      character(len=2)        :: naxis1,naxis2
      integer, dimension(*)   :: lztyp,lytyp

      real, dimension(3)      :: pt,at,bt,ct,a1,b1,c1,perp1,alab,blab,clab
      integer                 :: ierr,ifindz,ifindy,ia1,i
      real                    :: czx,czy,czz,cyx,cyy,cyz,oat,obt,oct,opt,oa1,ob1,oc1,operp1


! Get the vectors defining the two specified axes

      ierr=0
      ifindz=0
      ifindy=0
      czx=0.0
      czy=0.0
      czz=0.0
      cyx=0.0
      cyy=0.0
      cyz=0.0
      do iz=1,ncz
         ia1=lztyp(iz)
         ifindz=ifindz+1
         czx=czx+xt(ia1)
         czy=czy+yt(ia1)
         czz=czz+zt(ia1)
      end do
      do iy=1,ncy
         ia1=lytyp(iy)
         ifindy=ifindy+1
         cyx=cyx+xt(ia1)
         cyy=cyy+yt(ia1)
         cyz=cyz+zt(ia1)
      end do
      if(ifindz.gt.0) then
         ct(1)=czx/ifindz
         ct(2)=czy/ifindz
         ct(3)=czz/ifindz
      else
         ct(1)=0.0
         ct(2)=0.0
         ct(3)=1.0
      endif
      if(ifindy.gt.0) then
         pt(1)=cyx/ifindy
         pt(2)=cyy/ifindy
         pt(3)=cyz/ifindy
      else
         pt(1)=0.0
         pt(2)=0.0
         pt(3)=0.0
      endif
c     write(6,991) p1(1),p1(2),p1(3)
c     write(6,991) c1(1),c1(2),c1(3)
c renormalize these coordinates to unit vectors
c
      oct=sqrt(scapro(ct,ct))
c
c if oct is zero it means molecule is spherically symmetric, so no thetam can
c be calculated for this molecule
c
      if(oct.gt.0.0) then
         do i=1,3
            ct(i)=ct(i)/oct
         end do
         oct=1.0
      else
         oct=0.0
      endif
      opt=sqrt(scapro(pt,pt))
! if opt is zero it means molecule is linear so no chi dependence can be calculated.
      if(opt.gt.0.0) then
         do i=1,3
            pt(i)=pt(i)/opt
         end do
         opt=1.0
! vector product of pt with ct (z1) gives the at axis
         call vecpro(pt,ct,at)
         oat=sqrt(scapro(at,at))
! renormalize to give unit vector
         do i=1,3
            at(i)=at(i)/oat
         end do
         oat=1.0
! vector product of ct to at gives  bt
         call vecpro(ct,at,bt)
         obt=sqrt(scapro(bt,bt))
! renormalize to give unit vector
         do i=1,3
            bt(i)=bt(i)/obt
         end do
         obt=1.0
      else
         oat=0.0
         obt=0.0
      endif
      oa1=0.0
      ob1=0.0
      oc1=0.0
! put axis coordinates into correct arrays, depending on values of naxis1,naxis2
      if(naxis1.eq.'z'.and.naxis2.eq.'y') then
         do i=1,3
            a1(i)=at(i)
            b1(i)=bt(i)
            c1(i)=ct(i)
         end do
         oa1=oat
         ob1=obt
         oc1=oct
      else if(naxis1.eq.'-z'.and.naxis2.eq.'y') then
         do i=1,3
            a1(i)=at(i)
            b1(i)=bt(i)
            c1(i)=-ct(i)
         end do
         oa1=oat
         ob1=obt
         oc1=oct
      else if(naxis1.eq.'y'.and.naxis2.eq.'z') then
         do i=1,3
            a1(i)=-at(i)
            b1(i)=ct(i)
            c1(i)=bt(i)
         end do
         oa1=oat
         ob1=oct
         oc1=obt
      else if(naxis1.eq.'y'.and.naxis2.eq.'-z') then
         do i=1,3
            a1(i)=-at(i)
            b1(i)=ct(i)
            c1(i)=-bt(i)
         end do
         oa1=oat
         ob1=oct
         oc1=obt
      else if(naxis1.eq.'y'.and.naxis2.eq.'x') then
         do i=1,3
            a1(i)=bt(i)
            b1(i)=ct(i)
            c1(i)=at(i)
         end do
         oa1=obt
         ob1=oct
         oc1=oat
      else if(naxis1.eq.'x'.and.naxis2.eq.'y') then
         do i=1,3
            a1(i)=ct(i)
            b1(i)=bt(i)
            c1(i)=-at(i)
         end do
         oa1=oct
         ob1=obt
         oc1=oat
      else if(naxis1.eq.'z'.and.naxis2.eq.'x') then
         do i=1,3
            a1(i)=bt(i)
            b1(i)=-at(i)
            c1(i)=ct(i)
         end do
         oa1=obt
         ob1=oat
         oc1=oct
      else if(naxis1.eq.'-z'.and.naxis2.eq.'x') then
         do i=1,3
            a1(i)=bt(i)
            b1(i)=-at(i)
            c1(i)=-ct(i)
         end do
         oa1=obt
         ob1=oat
         oc1=oct
      else if(naxis1.eq.'x'.and.naxis2.eq.'z') then
         do i=1,3
            a1(i)=ct(i)
            b1(i)=at(i)
            c1(i)=bt(i)
         end do
         oa1=oct
         ob1=oat
         oc1=obt
      else if(naxis1.eq.'x'.and.naxis2.eq.'-z') then
         do i=1,3
            a1(i)=ct(i)
            b1(i)=at(i)
            c1(i)=-bt(i)
         end do
         oa1=oct
         ob1=oat
         oc1=obt
      endif
c
c now evaluate theta value for this molecule
c
      if(oc1.gt.0.0) then
         cthetm=scapro(clab,c1)
c
c vector product of clab to c1 gives the y' axis about which rotation thetam was
c made
c
         call vecpro(clab,c1,perp1)
c
c renormalize to give a unit vector
c
         operp1=sqrt(scapro(perp1,perp1))
         if(operp1.gt.0.0) then
            do i=1,3
               perp1(i)=perp1(i)/operp1
            end do
            operp1=1.0
            sphim=-scapro(perp1,alab)
            cphim=scapro(perp1,blab)
            phim=atan2(sphim,cphim)
         else
            phim=0.0
         endif
c now calculate chim
c
c scalar product of perp1 with a1 and b1 give chim rotation of molecule
c about molecules c-axis
c
         if(oa1.gt.0.0.and.ob1.gt.0.0) then
            schim=scapro(perp1,a1)
            cchim=scapro(perp1,b1)
            chim=atan2(schim,cchim)
         else
            chim=0.0
         endif
      else
         cthetm=0.0
         phim=0.0
         chim=0.0
      endif
      return
      end

      subroutine get_max_size_dlmn()

      use sharm_mod
      integer il,ilm,ntmax,nts,ntstep,mmin,mmax

      ilm=0
      do i=1,nll1v
         il=ll1v(i)
         ntstep=n1step
c n values
         if(ntstep.gt.0) then
            ntmax=il
            nts=ntstep
         else
            ntmax=0
            nts=1
         endif
         do in=0,ntmax,nts
            ilm=ilm+2*il+1
         end do
      end do
      mdlmn=ilm
      write(6,*) mdlmn,nll1v,ilm
      ilm=0
      do i=1,nll2v
         il=ll2v(i)
         ntstep=n2step
c n values
         if(ntstep.gt.0) then
            ntmax=il
            nts=ntstep
         else
            ntmax=0
            nts=1
         endif
         do in=0,ntmax,nts
            ilm=ilm+2*il+1
         end do
      end do
      mdlmn=max(mdlmn,ilm)
      write(6,*) mdlmn,nll2v,ilm
      ilm=0
      do ilt=1,nllv
         il=llv(ilt)
         ilm=ilm+2*il+1
      end do
      mdlmn=max(mdlmn,ilm)
      write(6,*) mdlmn,nllv,ilm
      return
      end

      subroutine form_dlmn(imol1,ierr)

      use reallocation_routines
      use sharm_mod

      include 'epsr_dimension.inc'
      include 'sharm.inc'
      include 'sharm_sdf.inc'
      integer                                :: nllt,ierr
      integer, dimension(:), allocatable     :: lltt
c
c calculate the D coefficients for this molecule
c
      ierr=0
      call reallocate1d_i(lltt,max(nll1v,nll2v))
      if(nm(imol1).eq.1) then
         nllt=nll1v
         do i=1,nll1v
            lltt(i)=ll1v(i)
         end do
         ntstep=n1step
      else
         nllt=nll2v
         do i=1,nll2v
            lltt(i)=ll2v(i)
         end do
         ntstep=n2step
      endif
c
c l values
c
      ilmnlm=0
      do ilt=1,nllt
         il=lltt(ilt)
c range of m values
         mmin=-il
         mmax=il
c n values
         if(ntstep.gt.0) then
            ntmax=il
            nts=ntstep
         else
            ntmax=0
            nts=1
         endif
         do in=0,ntmax,nts
            do m=mmin,mmax
               ilmnlm=ilmnlm+1
c     write(6,501) il,m,in,ilmnlm
c error exit if too many coeficients demanded
               if(ilmnlm.gt.mdlmn) then
                  ierr=ilmnlm-mdlmn
                  write(6,*) 'form_dlmn> Array reference ',ilmnlm,' larger than allowed ',mdlmn
                  return
               end if
               if(performtest) then
                  chiphi(1,imol1,ilmnlm)=m*phimt1(imol1)+in*chimt1(imol1)
                  dlmnmc(1,imol1,ilmnlm)=dlm(il,m,in,cthetmt1(imol1))
                  chiphit(1,imol1,ilmnlm)=m*phimt2(imol1)+in*chimt2(imol1)
                  dlmnmct(1,imol1,ilmnlm)=dlm(il,m,in,cthetmt2(imol1))
                  if(in.ne.0.and.(nm(imol1).eq.2.or.(nm(imol1).eq.1.and.ldistinct.eq.0))) then
                     chiphi(2,imol1,ilmnlm)=m*phimt1(imol1)-in*chimt1(imol1)
                     dlmnmc(2,imol1,ilmnlm)=dlm(il,m,-in,cthetmt1(imol1))
                     chiphit(2,imol1,ilmnlm)=m*phimt2(imol1)-in*chimt2(imol1)
                     dlmnmct(2,imol1,ilmnlm)=dlm(il,m,-in,cthetmt2(imol1))
                  end if
               else
                  chiphi(1,imol1,ilmnlm)=m*phim(imol1)+in*chim(imol1)
                  dlmnmc(1,imol1,ilmnlm)=dlm(il,m,in,cthetm(imol1))
                  if(in.ne.0.and.(nm(imol1).eq.2.or.(nm(imol1).eq.1.and.ldistinct.eq.0))) then
                     chiphi(2,imol1,ilmnlm)=m*phim(imol1)-in*chim(imol1)
                     dlmnmc(2,imol1,ilmnlm)=dlm(il,m,-in,cthetm(imol1))
                  end if
               end if
            end do
         end do
      end do
      if(performtest) then
         write(6,'(1x,i5,6(1x,e12.5))') imol1,phimt1(imol1),cthetmt1(imol1),chimt1(imol1),phimt2(imol1),cthetmt2(imol1),chimt2(imol1)
      else
         write(6,'(1x,i5,3(1x,e12.5))') imol1,phim(imol1),cthetm(imol1),chim(imol1)
      end if
      return
      end

      subroutine form_sum_prod_dlmn_cgcoeff(nmolsharm,ierr)

      use sharm_mod

      include 'epsr_dimension.inc'
      include 'sharm.inc'
      include 'sharm_sdf.inc'
      integer ierr,nsec,notype1,notype2,itimestotal,ntimestotal
      real phisum
      ierr=0
      nsec=0
      itimestotal=0
      gfac=4.0*pi/3.0
c
c initialize counters for number of molecules of each kind in the system
c
      if(ldistinct.eq.0) then
         if(nm(nmolsharm).eq.1.or.nm(nmolsharm).eq.2) then
            lsum1=1
            lsum2=1
         else
            lsum1=0
            lsum2=0
         endif
      else
         lsum1=0
         lsum2=0
         if(nm(nmolsharm).eq.1) lsum1=1
         if(nm(nmolsharm).eq.2) lsum2=1
      endif
! Count how many type 1 molecules there are
      notype1=0
      notype2=0
      do imol1=1,nmolsharm
         if(ldistinct.eq.0) then
            if(nm(imol1).eq.1.or.nm(imol1).eq.2) then
               notype1=notype1+1
               notype2=notype2+1
            endif
         else
            if(nm(imol1).eq.1) notype1=notype1+1
            if(nm(imol1).eq.2) notype2=notype2+1
         endif
      end do
      if(ldistinct.eq.0) then
         ntimestotal=notype1*(notype1-1)
      else
         ntimestotal=notype1*notype2
      end if
      write(6,*) nmolsharm,nm(nmolsharm),lsum1,lsum2
      do imol1=1,nmolsharm-1
         if(nm(imol1).ne.0) then
c     call system_clock(time3,rate,max)
c
c increment counters for number of molecules of types 1 and 2 in the system
c
            if(ldistinct.eq.0) then
               if(nm(imol1).eq.1.or.nm(imol1).eq.2) then
                  lsum1=lsum1+1
                  lsum2=lsum2+1
               endif
            else
               if(nm(imol1).eq.1) lsum1=lsum1+1
               if(nm(imol1).eq.2) lsum2=lsum2+1
            endif
            write(6,'(1x,i10,1x,i10,1x,i8,"% done",a)',advance='no') itimestotal,ntimestotal,itimestotal*100/ntimestotal,char(13)
            xc1=xcent(imol1)
            yc1=ycent(imol1)
            zc1=zcent(imol1)
c
c calc now look at second molecules
c
            nm2f=imol1+1
            do imol2=nm2f,nmolsharm
               if(imol1.ne.imol2.and.nm(imol2).ne.0) then
                  iflag1=0
                  iflag2=0
                  nm1=nm(imol1)
                  nm2=nm(imol2)
                  if(ldistinct.eq.0.and.nm1.eq.1.and.nm2.eq.1) then
                     iflag1=1
                     iflag2=1
                  endif
                  if(ldistinct.eq.1.and.nm1.eq.1.and.nm2.eq.2) then
                     iflag1=1
                     iflag2=0
                  endif
                  if(ldistinct.eq.1.and.nm1.eq.2.and.nm2.eq.1) then
                     iflag1=0
                     iflag2=1
                  endif
                  if((iflag1+iflag2).ne.0) then
                     itimestotal=itimestotal+iflag1+iflag2
                     xc2=xcent(imol2)
                     yc2=ycent(imol2)
                     zc2=zcent(imol2)
                     xdif=xc2-xc1
                     ydif=yc2-yc1
                     zdif=zc2-zc1
c apply periodic boundary conditions
                     call minimagec(xdif,ydif,zdif)
c                xc2=xc2+anint(xc1-xc2)
c                yc2=yc2+anint(yc1-yc2)
c                zc2=zc2+anint(zc1-zc2)
c calculate separation of molecular centres
c                pdif(1)=xc2-xc1
c                pdif(2)=yc2-yc1
c                pdif(3)=zc2-zc1
                     pdif(1)=xdif
                     pdif(2)=ydif
                     pdif(3)=zdif
c renormalize to unit vector
                     opdif=sqrt(scapro(pdif,pdif))
                     d12=opdif
                     do i=1,3
                        pdif(i)=pdif(i)/opdif
                     enddo
                     opdif=1.0
                     iref=int(d12/zstep)+1
c check that distance is within bounds prescribed for this distribution
                     if(iref.le.nz) then
c apply cosine rule to determine cos(theta) for intermolecular axis
                        cz=scapro(clab,pdif)
c vector perpendicular to c1 and pdif - gives the rotated y-axis
                        call vecpro(clab,pdif,perp1)
c normalize to give a unit vector
                        operp1=sqrt(scapro(perp1,perp1))
                        if(operp1.gt.0.0) then
                           do i=1,3
                              perp1(i)=perp1(i)/operp1
                           end do
                           operp1=1.0
c scalar product with a1 to give -sin(phi)
                           sphil=-scapro(alab,perp1)
c scalar product with b1 to give cos(phi)
                           cphil=scapro(blab,perp1)
c angle phi in lab coordinate system:
                           phil=atan2(sphil,cphil)
                        else
                           phil=0.0
                        endif
                        if(performtest) then
!                           call populate_disk(pi,phillab,thetallab,0.0,phinew,thetanew,chinew,thetaspread,phispread)
                           phil=phillab
                           cz=cos(thetallab)
                        end if
c
c calculate coefficients for laboratory coordinates
c
c
c l values
c
                        ilm=0
                        do ilt=1,nllv
                           il=llv(ilt)
! Don't bother with values disallowed by the CG coefficients
                           if(il.le.ll(ndist)) then
c range of m values
                              mmin=-il
                              mmax=il
                              do m=mmin,mmax
                                 ilm=ilm+1
c error exit if too many coeficients demanded
                                 philab(ilm)=m*phil
                                 dlmnc(1,ilm)=dlm(il,m,0,cz)
                                 dlmnc(2,ilm)=(-1)**m*dlm(il,m,0,-cz)
                              end do
                           endif
                        end do
c
c FINALLY we sum the products of coefficients
c
                        if(iflag1.eq.1) then !Molecule 1 is at the origin
                           ic=0
                           do id=1,ndist
                              il1=ll1(id)
                              in1=n1(id)
                              il2=ll2(id)
                              in2=n2(id)
                              il=ll(id)
c
c initialize reference integer for molecule 1
c
                              il1ref=il1rf(id)
c
c initialize reference integer for molecule 2
c
                              il2ref=il2rf(id)
c
c initialize reference integer for laboratory function
c
                              ilref=ilrf(id)
                              m1max=il1
                              m1min=-il1
                              do m1=m1min,m1max,1
                                 m2max=min(il2,il-m1)
                                 m2min=max(-il2,-il-m1)
!      if(il1.eq.3.and.il2.eq.3.and.il.eq.1.and.in1.eq.3.and.in2.eq.-3) then
!         ilmn1=il1ref+m1
!         ilmn2=il2ref+m2min
!         ilmn=ilref+m1+m2min
!         write (6,*) id,m1,m2min,m2max,il1ref+m1,il2ref,ilref
!         write(6,*) dlmnmc(1,imol1,ilmn1),dlmnmc(2,imol2,ilmn2),dlmnc(1,ilmn),termcg(ic+1)
!      end if
                                 do m2=m2min,m2max,1
                                    m=m1+m2
                                    ic=ic+1
                                    ilmn1=il1ref+m1
                                    ilmn2=il2ref+m2
                                    ilmn=ilref+m
                                    if(in2.lt.0) then
                                       if(performtest) then
                                          phisum=philab(ilmn)-chiphi(1,imol1,ilmn1)-chiphit(2,imol2,ilmn2)
                                          cosphi=cos(phisum)
                                          sinphi=sin(phisum)
                                          prod=dlmnmc(1,imol1,ilmn1)*dlmnmct(2,imol2,ilmn2)
                                       else
                                          phisum=philab(ilmn)-chiphi(1,imol1,ilmn1)-chiphi(2,imol2,ilmn2)
                                          cosphi=cos(phisum)
                                          sinphi=sin(phisum)
                                          prod=dlmnmc(1,imol1,ilmn1)*dlmnmc(2,imol2,ilmn2)
                                       end if
997   format(1x,6(1x,i4),3(1x,e12.5))
                                    else
                                       if(performtest) then
                                          phisum=philab(ilmn)-chiphi(1,imol1,ilmn1)-chiphit(1,imol2,ilmn2)
                                          cosphi=cos(phisum)
                                          sinphi=sin(phisum)
                                          prod=dlmnmc(1,imol1,ilmn1)*dlmnmct(1,imol2,ilmn2)
                                       else
                                          phisum=philab(ilmn)-chiphi(1,imol1,ilmn1)-chiphi(1,imol2,ilmn2)
                                          cosphi=cos(phisum)
                                          sinphi=sin(phisum)
                                          prod=dlmnmc(1,imol1,ilmn1)*dlmnmc(1,imol2,ilmn2)
                                       end if
                                    end if
                                    factor=prod*dlmnc(1,ilmn)*termcg(ic)
                                    add=factor*cosphi
                                    adds=factor*sinphi
                                    gofz(iref,id)=gofz(iref,id)+add
                                    gofzs(iref,id)=gofzs(iref,id)+adds
!      if(il1.eq.3.and.il2.eq.3.and.il.eq.1.and.in1.eq.3.and.in2.eq.-3) then
!      if(factor.gt.1.0e10.or.isnan(factor)) then
!         write(6,*) id,il1,il2,il,in1,in2,m1,m2,m
!         write(6,*) il1ref,il2ref,ilref,ilmn1,ilmn2,ilmn,ilm,size(dlmnc,1),size(dlmnc,2)
!         write(6,*) nllv,(llv(ilt),ilt=1,nllv)
!         write(6,*) dlmnmc(1,imol1,ilmn1),dlmnmc(2,imol2,ilmn2),dlmnc(1,ilmn),termcg(ic)
!         write(6,*) prod,cosphi,sinphi,add,adds,gofz(iref,id),gofzs(iref,id)
!      end if
                                 end do
                              end do
                           end do
                           nsec=nsec+1
                        endif
                        if(iflag2.eq.1.and..not.performtest) then !Molecule 2 is at the origin
                           ic=0
                           do id=1,ndist
                              il1=ll1(id)
                              in1=n1(id)
                              il2=ll2(id)
                              in2=n2(id)
                              il=ll(id)
c
c initialize reference integer for molecule 1
c
                              il1ref=il1rf(id)
c
c initialize reference integer for molecule 1
c
                              il2ref=il2rf(id)
c
c initialize reference integer for laboratory function
c
                              ilref=ilrf(id)
                              m1max=il1
                              m1min=-il1
                              do m1=m1min,m1max,1
                                 m2max=min(il2,il-m1)
                                 m2min=max(-il2,-il-m1)
                                 do m2=m2min,m2max,1
                                    m=m1+m2
                                    ic=ic+1
                                    ilmn1=il1ref+m1
                                    ilmn2=il2ref+m2
                                    ilmn=ilref+m
                                    if(in2.lt.0) then
                                       if(performtest) then
                                          phisum=philab(ilmn)-chiphi(1,imol2,ilmn1)-chiphit(2,imol1,ilmn2)
                                          cosphi=cos(phisum)
                                          sinphi=sin(phisum)
                                          prod=dlmnmc(1,imol2,ilmn1)*dlmnmct(2,imol1,ilmn2)
                                       else
                                          phisum=philab(ilmn)-chiphi(1,imol2,ilmn1)-chiphi(2,imol1,ilmn2)
                                          cosphi=cos(phisum)
                                          sinphi=sin(phisum)
                                          prod=dlmnmc(1,imol2,ilmn1)*dlmnmc(2,imol1,ilmn2)
                                       end if
                                    else
                                       if(performtest) then
                                          phisum=philab(ilmn)-chiphi(1,imol2,ilmn1)-chiphit(1,imol1,ilmn2)
                                          cosphi=cos(phisum)
                                          sinphi=sin(phisum)
                                          prod=dlmnmc(1,imol2,ilmn1)*dlmnmct(1,imol1,ilmn2)
                                       else
                                          phisum=philab(ilmn)-chiphi(1,imol2,ilmn1)-chiphi(1,imol1,ilmn2)
                                          cosphi=cos(phisum)
                                          sinphi=sin(phisum)
                                          prod=dlmnmc(1,imol2,ilmn1)*dlmnmc(1,imol1,ilmn2)
                                       end if
                                    end if
                                    factor=prod*dlmnc(2,ilmn)*termcg(ic)
                                    add=factor*cosphi
                                    adds=factor*sinphi
                                    gofz(iref,id)=gofz(iref,id)+add
                                    gofzs(iref,id)=gofzs(iref,id)+adds
!      if(il1.eq.3.and.il2.eq.3.and.il.eq.1.and.in1.eq.3.and.in2.eq.-3) then
!         write(6,*) id,m1,m2,ilmn1,ilmn2,ilmn
!         write(6,*) dlmnmc(2,imol1,ilmn1),dlmnmc(1,imol2,ilmn2),dlmnc(1,ilmn)
!         write(6,*) prod,cosphi,dlmnc(2,ilmn),termcg(ic),add,gofz(iref,id)
!      end if
                                 end do
                              end do
                           end do
                           nsec=nsec+1
                        endif
                     endif
                  endif
               endif
            end do
c     call system_clock(time4,rate,max)
c     write(6,*) rate,max,time4-time3,nsec
         endif
      end do
      write(6,'(1x,i10,1x,i10,1x,i8,"% done")') itimestotal,ntimestotal,itimestotal*100/ntimestotal
502   format(1x,i10,2(1x,i5))
500   format(1x,i6,5(1x,i3),3(1x,i5),1x,e12.5)
      do id=1,ndist
         write(6,500) id,ll1(id),ll2(id),ll(id),n1(id),n2(id),il1rf(id),il2rf(id),ilrf(id),alsec(id)
      end do
c renormalize g(r)'s to total number of atoms
c density of second atoms in sphere
      rho=float(lsum1)*float(lsum2)/vol
      write(6,193) rho
193   format(1x,'rho = ',e13.6)
      gfac1=gfac*rho
      do id=1,ndist
            z0=0.0
            z03=0.0
            do iz=1,nz
                  zm=z(iz)
                  zdiff=zm-z0
                  z1=zm+zdiff
                  z13=z1*z1*z1
c average number of particles in this volume element
                  rat=gfac1*(z13-z03)
                  gofz(iz,id)=(2*ll1(id)+1)*(2*ll2(id)+1)*gofz(iz,id)/rat
                  gofzs(iz,id)=(2*ll1(id)+1)*(2*ll2(id)+1)*gofzs(iz,id)/rat
                  if(id.eq.1) gofz(iz,id)=gofz(iz,id)-1.0
                  z0=z1
                  z03=z13
c           write(6,996) iz,zm,gofz(iz,id)
996   format(1x,i5,2(1x,e12.5))
            end do
      end do
      return
      end

      subroutine populate_disk(pi,phi0,theta0,chi0,phinew,thetanew,chinew,thetaspread,phispread)
! Randomly populates a tangential disk at (theta0,phi0) with a Gaussian-like distribution at constant radius,
! giving a new (theta,phi). Generates random chi about chi0.
      real pi,ctheta0,phi0,chi0,thetaspread,phispread,thetanew,phinew,chinew
      real x0,y0,z0,z02,deltax,deltay,stheta0,root1,root2,cthetanew,sthetanew,u1,u2,amp,cons

      ctheta0=cos(theta0)
      stheta0=sin(theta0)
      x0=stheta0*cos(phi0)
      y0=stheta0*sin(phi0)
      z0=ctheta0
      z02=z0*z0
      cons=-1.0
      do while(cons.lt.0.0)
         u1=ran1()
         u2=ran1()
         amp=sqrt(-2.0*log(u1))
         deltax=amp*cos(2.0*pi*u2)*phispread
         deltay=amp*sin(2.0*pi*u2)*thetaspread
         cons=z02-(deltax*(deltax+2.0*x0)+deltay*(deltay+2.0*y0))
      end do
! Choose the root closest to z0
      cthetanew=sqrt(cons)
      root1=acos(cthetanew)
      root2=acos(-cthetanew)
      dif1=abs(theta0-root1)/thetaspread
      dif2=abs(theta0-root2)/thetaspread
! If either or these differences is greater than 5 deviations, choose the other, otherwise toss a coin.
      if((dif2-dif1).gt.5.0) then
         thetanew=root1
      else if((dif1-dif2).gt.5.0) then
         thetanew=root2
      else !Toss a coin to choose between the two roots
         if(ran1().lt.0.5) then
            thetanew=root1
         else
            thetanew=root2
         end if
      end if
      sthetanew=sqrt(1.0-cthetanew*cthetanew)
      cphinew=(x0+deltax)/sthetanew
      sphinew=(y0+deltay)/sthetanew
      phinew=atan2(sphinew,cphinew)
      u1=u2
      u2=ran1()
      chinew=chi0+sqrt(-2*log(u1))*cos(2*pi*u2)*phispread
!      write(6,*) x0,y0,z0,deltax,deltay,root1,root2,dif1,dif2
!      write(6,*) phi0,theta0,chi0,phinew,thetanew,chinew
      return
      end


