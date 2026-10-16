      module containerforces_mod

           integer                                                :: numcontainerforces
           integer, dimension(:), allocatable                     :: numpairs,numcentraltype,numsecondarytype
           integer, dimension(:,:), allocatable                   :: pairrefstore
           integer, dimension(:,:), allocatable                   :: centraltype !List of atoms (q-type) at centres of "containers"
           integer, dimension(:,:), allocatable                   :: secondarytype !List of atom types to be contained (not q-type)

           real, dimension(:), allocatable                        :: ergold_container,ergnew_container

      contains

      subroutine get_containerforce_list(nmol,ic,rmaxex,it,jt,itype,qtype,natom,qtypea,ntype)

      use reallocation_routines

      real                                   :: rmaxex
      integer                                :: ia,ic,it,jt,im,natom,nt,ntype,numcentral,numsecondary,typecentral,typesecondary
      integer, dimension(:)                  :: itype,qtype,qtypea

      !initialise numcontainerforces if ic .eq. 1
      if(ic.eq.1) then
         numcontainerforces=0
         call reallocate1d_r(ergold_container,nmol)
         call reallocate1d_r(ergnew_container,nmol)
         do im=1,nmol
            ergold_container(im)=0.0
            ergnew_container(im)=0.0
         end do
      end if
      call reallocate1d_i(numpairs,ic)
      numpairs(ic)=0
      if (rmaxex.le.0.0) return ! Will only count types if rmaxex .gt. 0.0)
      numcontainerforces=numcontainerforces+1
! First determine the type of the central and secondary atoms
      do nt=1,ntype
         if(qtypea(nt).eq.0) then
            if(nt.eq.it) then
               typesecondary=it
            else if(nt.eq.jt) then
               typesecondary=jt
            end if
         else
            if(nt.eq.it) then
               typecentral=it
            else if(nt.eq.jt) then
               typecentral=jt
            end if
         end if
      end do
! Count the number of atoms of types it and jt
      numcentral=0
      numsecondary=0
      do ia=1,natom
         if(qtype(ia).eq.1.and.(itype(ia).eq.it.or.itype(ia).eq.jt)) then
            numcentral=numcentral+1
         else if(qtype(ia).eq.0.and.(itype(ia).eq.it.or.itype(ia).eq.jt)) then
            numsecondary=numsecondary+1
         end if
      end do
      numpairs(ic)=numcentral*numsecondary
      if(numpairs(ic).gt.0) then
         call reallocate2d_i(pairrefstore,2,numcontainerforces)
         call reallocate1d_i(numcentraltype,numcontainerforces)
         call reallocate1d_i(numsecondarytype,numcontainerforces)
         call reallocate2d_i(centraltype,numcentral,numcontainerforces)
! Store the secondary atom type number for this atom pair - this will be used to identify secondary atoms
         pairrefstore(1,numcontainerforces)=typecentral
         pairrefstore(2,numcontainerforces)=typesecondary
         numcentraltype(numcontainerforces)=numcentral
         numsecondarytype(numcontainerforces)=numsecondary
! Now assign the values
         numcentral=0
         do ia=1,natom
            if(qtype(ia).eq.1.and.(itype(ia).eq.it.or.itype(ia).eq.jt)) then
               numcentral=numcentral+1
               centraltype(numcentral,numcontainerforces)=ia
            end if
         end do
      else
! Just in case no pairs were found in fact - unlikely.
         numpairs(ic)=0
         if(numcontainerforces.gt.0) numcontainerforces=numcontainerforces-1
      end if

      end subroutine get_containerforce_list

      function en_containerforce(im,naf,nal,xadd,yadd,zadd,xa,ya,za &
         ,itype,mref,npairref,irefpcof,rmaxex,pcoefexmax,expecf,expecfpwr,dmax,harmonic)

! Calculates a quadratic energy if a specified atom pair has a minimum separation greater than specified rmaxex.

      real                                :: en_containerforce,expecf,dmax,mindist,gradlimit,loggradlimit,expecfuse
      real, dimension(:)                  :: xadd,yadd,zadd,xa,ya,za,rmaxex,pcoefexmax
      integer                             :: im,naf,nal,ia,icref,ic,jc,irefp
      integer, dimension(:)               :: itype,mref,irefpcof
      integer, dimension(:,:)             :: npairref
      logical harmonic

      en_containerforce=0.0
      if(numcontainerforces.eq.0) return
      sume=0.0
      gradlimit=expecfpwr !Limits the gradient on the container force to prevent the energy become too large
! If gradlimit is zero or less, then make the gradient limit very large - this is to make program compatible with previous versions
		if(gradlimit.gt.0.0) then
	      loggradlimit=log(gradlimit)
		else
			gradlimit=1.0
			loggradlimit=0.0
		end if
      do ic=1,numcontainerforces
         do ia=naf,nal
            if(itype(ia).eq.pairrefstore(2,ic)) then
               icref=npairref(itype(ia),pairrefstore(1,ic))
               irefp=irefpcof(icref)
               expecfuse=expecf*sqrt(abs(rmaxex(irefp)))
!               write(6,*) itype(ia),pairrefstore(1,ic),icref,irefp
               mindist=0.0
               x1=xadd(im)+xa(ia)
               y1=yadd(im)+ya(ia)
               z1=zadd(im)+za(ia)
! Need the minimum distance between ia and all the "central" atoms
               do jc=1,numcentraltype(ic)
                  ja=centraltype(jc,ic)
                  jm=mref(ja)
                  x2=xadd(jm)+xa(ja)
                  y2=yadd(jm)+ya(ja)
                  z2=zadd(jm)+za(ja)
                  xdif=x2-x1
                  ydif=y2-y1
                  zdif=z2-z1
! apply minimum image convention
                  call minimagec(xdif,ydif,zdif)
                  dsq=xdif*xdif+ydif*ydif+zdif*zdif
                  if(mindist.eq.0.0.or.dsq.lt.mindist) then
                     mindist=dsq
                  end if
               end do
               dtest=sqrt(mindist)*dmax
               if(.not.harmonic) then
                  rtest=rmaxex(irefp)
                  ddiff=(dtest-rtest)/expecfuse
                  if(gradlimit.gt.0.0.and.ddiff.ge.loggradlimit) then
                     sume=sume+abs(pcoefexmax(irefp))*gradlimit*(1.0+ddiff-loggradlimit) !Put linear term to prevent divergent exponential values
                  else
                     sume=sume+abs(pcoefexmax(irefp))*exp(ddiff)
                  end if
               else
                  rtest=rmaxex(irefp)-expecfuse
                  ddiff=(dtest-rtest)/expecfuse
                  if(ddiff.gt.0.0) then
							if(gradlimit.gt.0.0.and.ddiff.gt.gradlimit) then
								sume=sume+0.5*abs(pcoefexmax(irefp))*gradlimit*(2*ddiff-gradlimit)
							else
								sume=sume+0.5*abs(pcoefexmax(irefp))*ddiff*ddiff
							endif
						end if
               end if
!               if(isnan(sume)) then
!                  write(6,*) 'en_containerforce> ',im,itype(ia),irefp,rtest,ddiff,gradlimit,pcoefexmax(irefp)
!                  stop
!               end if
            end if
         end do
      end do
      en_containerforce=sume
         
      end function en_containerforce

      subroutine dist_containerforce(im,naf,nal,xadd,yadd,zadd,xa,ya,za &
         ,itype,mref,npairref,irefpcof,rmaxex,potexmax,dmax,addwav)

! For each requested container force, calculates whether atom pair exceeds maximum distance. If so, increments weighting on this term

      real                                :: expecf,dmax,mindist,addwav
      real, dimension(:)                  :: xadd,yadd,zadd,xa,ya,za,rmaxex,potexmax
      integer                             :: im,naf,nal,ia,icref,ic,jc,irefp
      integer, dimension(:)               :: itype,mref,irefpcof
      integer, dimension(:,:)             :: npairref

      if(numcontainerforces.eq.0) return
      do ic=1,numcontainerforces
!         addpcoef=addwav/real(numsecondarytype(ic))
         addpcoef=addwav
         do ia=naf,nal
            if(itype(ia).eq.pairrefstore(2,ic)) then
               icref=npairref(itype(ia),pairrefstore(1,ic))
               irefp=irefpcof(icref)
!               write(6,*) itype(ia),pairrefstore(1,ic),icref,irefp
               mindist=0.0
               x1=xadd(im)+xa(ia)
               y1=yadd(im)+ya(ia)
               z1=zadd(im)+za(ia)
! Need the minimum distance between ia and all the "central" atoms
               do jc=1,numcentraltype(ic)
                  ja=centraltype(jc,ic)
                  jm=mref(ja)
                  x2=xadd(jm)+xa(ja)
                  y2=yadd(jm)+ya(ja)
                  z2=zadd(jm)+za(ja)
                  xdif=x2-x1
                  ydif=y2-y1
                  zdif=z2-z1
! apply minimum image convention
                  call minimagec(xdif,ydif,zdif)
                  dsq=xdif*xdif+ydif*ydif+zdif*zdif
                  if(mindist.eq.0.0.or.dsq.lt.mindist) then
                     mindist=dsq
                  end if
               end do
               dtest=sqrt(mindist)*dmax
               if(dtest.gt.rmaxex(irefp)) potexmax(icref)=potexmax(icref)+addpcoef
            end if
         end do
      end do
         
      end subroutine dist_containerforce

      end module containerforces_mod

