      subroutine molfit(fadds,faddbs,nisum,intra,ireset)
      
      use timer
      use cellneighbourlist
      use ergpairreference
                    
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      real, dimension(:), allocatable         :: xr,yr,zr
      real stepmithis,steprithis
      
      temp=1.0/simt
      tol=tolout/temp
      do j=1,2
         do imoltype=1,nmoltype
            nmovei(j,imoltype)=0
            ninter(j,imoltype)=0
            nwhole(j,imoltype)=0
            nmovem(j,imoltype)=0
         end do
      end do

c Set the maximum translational step size

      trstepmax=1.0/real(natomr)
      trstepmax=trstepmax**(1.0/3.0)
      trstepmin=max(0.0001*trstepmax,1.0e-6) !Limit the smallest step size to that which will make a difference within the available precision.
      rotstepmax=1.0
      rotstepmin=1.0e-6
      fadd=0.0
      faddb=0.0
c
c if step size is non-zero, then do internal molecule structure refinement. This is to 
c prevent the build up of round-off errors from allowing the molecule to 
c become too distorted.
c
      if(nisum.ge.intra.and.stepmi.gt.0.0) then

         do im=1,nmol
            if(nato(im).gt.1) then
               call make_neighbour_list_from_mols(im)
!         write(6,*) im,nergpairm(im),molergpair(1,im),ergpair(ergpairref(1,im))
!         write(6,*) im,nnbref,ergold(1),ergnew(1)
               ! first save the coordinates of this molecule
               stepmithis=interntrans(mtype(im))
               if(stepmi.gt.0.0.and.stepmithis.eq.0.0) stepmithis=trstepmin
               interntrans(mtype(im))=stepmithis
               xs=xadd(im)
               ys=yadd(im)
               zs=zadd(im)
               if(.not.allocated(xr)) then
                   allocate(xr(nato(im)),yr(nato(im)),zr(nato(im)))
               else if (size(xr).lt.nato(im)) then
                   deallocate(xr,yr,zr)
                   allocate(xr(nato(im)),yr(nato(im)),zr(nato(im)))
               end if
               nat=naf(im)-1
               do ia=1,nato(im)
                  nat=nat+1
                  xr(ia)=xa(nat)
                  yr(ia)=ya(nat)
                  zr(ia)=za(nat)
               end do
               fdiff=0.0
               fdiffb=0.0
               
c move particles and calculate new chisq
! do not move the first atom of each molecule if this is the only
! atom in the molecule!
               nmovei(1,mtype(im))=nmovei(1,mtype(im))+1
               call movei(im,stepmithis)
!               call movea(nat,stepmithis)
               call en_change_intra(im,fdiffb)
               ediff=fdiffb*temp
               if(ediff.le.0.0) then
                  etest=10.0
               else if(ediff.lt.20.) then
                  etest=exp(-ediff)
               else
                  etest=0.
               end if
               rtest=ran1()
               if(rtest.lt.etest) then
                  faddb=faddb+fdiffb
                  call save_energies()
               else
                  nmovei(2,mtype(im))=nmovei(2,mtype(im))+1
                  nat=naf(im)-1
                  do ia=1,nato(im)
                     nat=nat+1
                     xa(nat)=xr(ia)
                     ya(nat)=yr(ia)
                     za(nat)=zr(ia)
                  end do
!                  xa(nat)=xr(ia)
!                  ya(nat)=yr(ia)
!                  za(nat)=zr(ia)
               end if
!               write(6,*) 'molfit> ',im,nat,ia,ediff
!               end do
c now try rotating internal groups if no of rotational groups is greater than 0
! This doesn't work because you need to include the interactions with neighbouring 
! molecules
c now adjust the coordinates so that the centre of mass has not changed
               call centreofmassi(im)
               !Save all the relevant pair energies
               call save_erg_pairs(im)
            else
               interntrans(mtype(im))=0.0
            endif
         end do
      else
!Step through cells and molecules 
         do ic=1,ncell
            if(nmolsincell(ic).gt.0) then
               do imref=1,nmolsincell(ic)
                  im=molsincell(imref,ic)
                  if(.not.allocated(xr)) then
                     allocate(xr(nato(im)),yr(nato(im)),zr(nato(im)))
                  else if (size(xr).lt.nato(im)) then
                     deallocate(xr,yr,zr)
                     allocate(xr(nato(im)),yr(nato(im)),zr(nato(im)))
                  end if
                  call system_clock(time1)
                  call make_neighbour_list_from_mols(im)
                  call system_clock(time2)
                  sumtime1=sumtime1+(time2-time1)/1000
                  if(stepri.gt.0.0.and.ngrot(im).gt.0) then
                     if(internrot(mtype(im)).eq.0.0) internrot(mtype(im))=rotstepmin
! Save the coordinates of this molecule
                     nat=naf(im)-1
                     do ia=1,nato(im)
                        nat=nat+1
                        xr(ia)=xa(nat)
                        yr(ia)=ya(nat)
                        zr(ia)=za(nat)
                     end do
                     fdiff=0.0
                     fdiffb=0.0
c
c perform internal rotations
c
                     ninter(1,mtype(im))=ninter(1,mtype(im))+1
                     do ig=nrf(im),nrl(im)
                        call rotati(ig,internrot(mtype(im)))
                     end do
 
c calculate change in intramolecular energy
c
                     call en_change_intra(im,fdiffb)
                     call en_change(im,fdiff)
c                        fdiff=0.0
                     ediff=(fdiff+fdiffb)*temp
                     if(ediff.le.0.0) then
                        etest=10.0
                     else if(ediff.lt.20.) then
                        etest=exp(-ediff)
                     else
                        etest=0.
                     endif
!                        write(6,*) 'Internal rotation: ',im,ediff
                     if(ran1().lt.etest) then
                        nstepreject(im)=0
                        fadd=fadd+fdiff
                        faddb=faddb+fdiffb
                        call save_energies()
                     else
! Move rejected - restore coordinates
                        nstepreject(im)=nstepreject(im)+1
                        do ia=1,nato(im)
                           nat=naf(im)+(ia-1)
                           xa(nat)=xr(ia)
                           ya(nat)=yr(ia)
                           za(nat)=zr(ia)
                        end do
                        ninter(2,mtype(im))=ninter(2,mtype(im))+1
                     end if
! Adjust the coordinates of centre of mass to correct for any change
                     call centreofmass(im)
! Save all the relevant pair energies
                     call save_erg_pairs(im)
                  else
                     internrot(mtype(im))=0.0
                  end if
! Rotate whole molecules and calculate new chisq. Only rotate molecules that have more than one atom
                  if(stepr.gt.0.0.and.nato(im).gt.1) then
! first save the coordinates of this molecule
                     nat=naf(im)-1
                     do ia=1,nato(im)
                        nat=nat+1
                        xr(ia)=xa(nat)
                        yr(ia)=ya(nat)
                        zr(ia)=za(nat)
                     end do
                     fdiff=0.0
                     fdiffb=0.0
                     if(wholerot(mtype(im)).le.0.0) wholerot(mtype(im))=rotstepmin
                     nwhole(1,mtype(im))=nwhole(1,mtype(im))+1
                     call rotatm(im,wholerot(mtype(im)))
! Calculate change in intramolecular and intermolecular energy
                     call en_change_intra(im,fdiffb)
                     call en_change(im,fdiff)
                     ediff=(fdiff+fdiffb)*temp
                     if(ediff.lt.0.0) then
                        etest=10.0
                     else if(ediff.lt.20.) then
                        etest=exp(-ediff)
                     else
                        etest=0.
                     endif
                     if(ran1().lt.etest) then
                        nstepreject(im)=0
                        fadd=fadd+fdiff
                        faddb=faddb+fdiffb
                        call save_energies()
                     else
! Move rejected - restore coordinates
                        nstepreject(im)=nstepreject(im)+1
                        do ia=1,nato(im)
                           nat=naf(im)+(ia-1)
                           xa(nat)=xr(ia)
                           ya(nat)=yr(ia)
                           za(nat)=zr(ia)
                        end do
                        nwhole(2,mtype(im))=nwhole(2,mtype(im))+1
                     end if
! Adjust the coordinates of centre of mass to correct for any change
                     call centreofmass(im)
! Save all the relevant pair energies
                     call save_erg_pairs(im)
                  else
                     wholerot(mtype(im))=0.0
                  end if
! Translate molecules
                  if(stepm.gt.0.0) then
! Save the coordinates of this molecule
                     xs=xadd(im)
                     ys=yadd(im)
                     zs=zadd(im)
                     fdiff=0.0
                     fdiffb=0.0
                     ! now translate each molecule whole
                     if(wholetrans(mtype(im)).le.0.0) wholetrans(mtype(im))=trstepmin
                     nmovem(1,mtype(im))=nmovem(1,mtype(im))+1
                     call movem(im,wholetrans(mtype(im)))
                     ! calculate change in intramolecular and intermolecular energy
                     call en_change_intra(im,fdiffb)
                     call en_change(im,fdiff)
                     ediff=(fdiff+fdiffb)*temp
                     if(ediff.lt.0.0) then
                        etest=10.0
                     else if(ediff.lt.20.) then
                        etest=exp(-ediff)
                     else
                        etest=0.
                     endif
                     if(ran1().lt.etest) then
                        nstepreject(im)=0
                        fadd=fadd+fdiff
                        faddb=faddb+fdiffb
                        call save_energies()
                     else
! Move rejected - restore coordinates
                        nstepreject(im)=nstepreject(im)+1
                        xadd(im)=xs
                        yadd(im)=ys
                        zadd(im)=zs
                        nmovem(2,mtype(im))=nmovem(2,mtype(im))+1
                     end if
! Adjust the coordinates of centre of mass to correct for any change
                     call centreofmass(im)
! Save all the relevant pair energies
                     call save_erg_pairs(im)
                  else
                     wholetrans(mtype(im))=0.0
                  end if
               end do
            end if
         end do
      end if

      do im=1,nmol
          call minimagec(xadd(im),yadd(im),zadd(im))
          if(tethered(im)) call minimagec(phix(im),phiy(im),phiz(im))
      end do

      if(nisum.ge.intra) then
         nisum=1
      else
         nisum=nisum+1
      end if

!Adjust step sizes to give required rejection rate

      if(ireset.ne.1) then

         do imoltype=1,nmoltype

            if(nmovei(2,imoltype).gt.0) then
               rat=rejrate*real(nmovei(1,imoltype))/real(nmovei(2,imoltype))
               if(rat.gt.1.1) rat=1.1
               if(rat.lt.0.9) rat=0.9
            else
               if(nmovei(1,imoltype).gt.0) then
                  rat=1.1
               else
                  rat=1.0
               endif
            end if
            interntrans(imoltype)=interntrans(imoltype)*rat
            if(interntrans(imoltype).gt.trstepmax) interntrans(imoltype)=trstepmax
            if(interntrans(imoltype).lt.trstepmin.and.interntrans(imoltype).ne.0.0) interntrans(imoltype)=trstepmin

            if(ninter(2,imoltype).gt.0) then
               rat=rejrate*real(ninter(1,imoltype))/real((ninter(2,imoltype)))
               if(rat.gt.1.1) rat=1.1
               if(rat.lt.0.9) rat=0.9
            else
               if(ninter(1,imoltype).gt.0) then
                  rat=1.1
               else
                  rat=1.0
               endif
            endif
            internrot(imoltype)=internrot(imoltype)*rat
            if(internrot(imoltype).gt.rotstepmax) internrot(imoltype)=rotstepmax
            if(internrot(imoltype).lt.rotstepmin.and.internrot(imoltype).ne.0.0) internrot(imoltype)=rotstepmin

            if(nwhole(2,imoltype).gt.0) then
               rat=rejrate*real(nwhole(1,imoltype))/real((nwhole(2,imoltype)))
               if(rat.gt.1.1) rat=1.1
               if(rat.lt.0.9) rat=0.9
            else
               if(nwhole(1,imoltype).gt.0) then
                  rat=1.1
               else
                  rat=1.0
               endif
            end if
            !Make sure a step size is set if stepr > 0
            
            wholerot(imoltype)=wholerot(imoltype)*rat
            if(wholerot(imoltype).gt.rotstepmax) wholerot(imoltype)=rotstepmax
            if(wholerot(imoltype).lt.rotstepmin.and.wholerot(imoltype).ne.0.0) wholerot(imoltype)=rotstepmin

            if(nmovem(2,imoltype).gt.0) then
               rat=rejrate*real(nmovem(1,imoltype))/real((nmovem(2,imoltype)))
               if(rat.gt.1.1) rat=1.1
               if(rat.lt.0.9) rat=0.9
            else
               if(nmovem(1,imoltype).gt.0) then
                  rat=1.1
               else
                  rat=1.0
               endif
            endif
            wholetrans(imoltype)=wholetrans(imoltype)*rat
            if(wholetrans(imoltype).gt.trstepmax) wholetrans(imoltype)=trstepmax
            if(wholetrans(imoltype).lt.trstepmin.and.wholetrans(imoltype).ne.0.0) wholetrans(imoltype)=trstepmin
!Adjust the whole rotation step size so that ratio with translation stays constant
            if(stepm.gt.0.0.and.wholetrans(imoltype).gt.0.0) then
               wholerot(imoltype)=wholetrans(imoltype)*stepr/(stepm*dmax)
            end if

         end do
      endif
 
      fadds=fadds+fadd/nmol
      faddbs=faddbs+faddb/nmol

      return
      end

      subroutine movei(im,stepmithis)

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      
      integer im,iref,ia
      real stepmithis
c 
c move the atoms
c
      do ia=naf(im),nal(im)
         call movea(ia,stepmithis)
      end do
      return
      end

      subroutine movea(ia,stepmithis)

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      real stepmithis

!      amass=sqrt(amassa(itype(ia)))
!      stepsz=stepmithis/amass
      stepsz=stepmithis
c move atoms by random amounts
      za(ia)=za(ia)+stepsz*(2.*ran1()-1)
      ya(ia)=ya(ia)+stepsz*(2.*ran1()-1)
      xa(ia)=xa(ia)+stepsz*(2.*ran1()-1)
      return
      end

      subroutine interm(im,fdiff,fdiffb)

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'

c
c perform internal rotations
c
      do ig=nrf(im),nrl(im)
         call rotati(ig,internrot(mtype(im)))
      end do
c
c calculate change in energy
c
      fdiff=0.0
      fdiffb=0.0
      call en_change_intra(im,fdiffb)
      call en_change(im,fdiff)
      return
      end

      subroutine rotati(ig,steprithis)

      use molecom_mod
		
		include 'epsr_dimension.inc'
      include 'molecom.inc'

c pvector is the vector about which the rotation occurs

      real pvector(3),steprithis

c xvector is the vector from centre of rotation to an individual atom
c xpvector is the component of this vector along the axis of rotation
c xrvector is the component of this vector perpendicular to the axis of rotation
c xrpvector is the vector perpendicular to both axis of rotation and xrvector

      real xvector(3),xpvector(3),xrvector(3),xrpvector(3)

      pi=3.141592653
c define random new orientation
      phi=pi*steprithis*(2*ran1()-1)
      cphi=cos(phi)
      sphi=sin(phi)
c define length of axis between the two points - the direction 1 -> 2 is
c treated as the positive going transformed z axis 
      pvector(1)=xa(labaxi(ig,2))-xa(labaxi(ig,1))
      pvector(2)=ya(labaxi(ig,2))-ya(labaxi(ig,1))
      pvector(3)=za(labaxi(ig,2))-za(labaxi(ig,1))
      rab2=scapro(pvector,pvector)
      if(rab2.gt.0) then
         rab=sqrt(rab2)
c Unit vector along axis of rotation
         do i=1,3
            pvector(i)=pvector(i)/rab
         end do
c step through atoms in group
         do ia=rotf(ig),rotl(ig)
            nat=labrot(ia)
c first calculate vector from centre of rotation to each atom in group
            nax=labaxi(ig,1)
            xvector(1)=xa(nat)-xa(nax)
            xvector(2)=ya(nat)-ya(nax)
            xvector(3)=za(nat)-za(nax)
c Vector parallel to axis of rotation
            frac=scapro(pvector,xvector)
            do i=1,3
               xpvector(i)=frac*pvector(i)
            enddo
c Vector perpendicular to axis of rotation
            do i=1,3
               xrvector(i)=xvector(i)-xpvector(i)
            enddo
c Vector perpendicular to plane defined by axis of rotation and xrvector
            call vecpro(pvector,xrvector,xrpvector)
c Generate rotated coordinates. The vector component parallel to the axis of rotation
c is not rotated
            xa(nat)=xa(nax)+xpvector(1)
            ya(nat)=ya(nax)+xpvector(2)
            za(nat)=za(nax)+xpvector(3)
            xa(nat)=xa(nat)+xrvector(1)*cphi+xrpvector(1)*sphi
            ya(nat)=ya(nat)+xrvector(2)*cphi+xrpvector(2)*sphi
            za(nat)=za(nat)+xrvector(3)*cphi+xrpvector(3)*sphi
         end do
      endif
      return
      end

      subroutine wholem(im)

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
c
c rotate entire molecule
c
      call rotatm(im,wholerot(mtype(im)))
c
c Calculate the change in energy
c
!      fdiff=0.0
!      fdiffb=0.0
!      call en_change(im,fdiff)
      return
      end

      subroutine rotatm(im,stepsz)
c routine to rotate entire molecule

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      pi=3.141592653
      iref1=naf(im)
      iref2=nal(im)
      if(iref2.gt.iref1) then
c choose an axis to rotate about
      phi=pi*stepsz*(2*ran1()-1)
      cphi=cos(phi)
      sphi=sin(phi)
c define rotations about x,y or z-axes
!      test=3.*ran1()
!      if(test.lt.1.0) then
c rotate about z axis
            do ia=iref1,iref2
                  xasave=xa(ia)*cphi-ya(ia)*sphi
                  ya(ia)=xa(ia)*sphi+ya(ia)*cphi
                  xa(ia)=xasave
            end do
!      else if(test.lt.2.) then
c rotate about y-axis
            do ia=iref1,iref2
                  zasave=za(ia)*cphi-xa(ia)*sphi
                  xa(ia)=za(ia)*sphi+xa(ia)*cphi
                  za(ia)=zasave
            end do
!      else
c rotate about x-axis
            do ia=iref1,iref2
                  yasave=ya(ia)*cphi-za(ia)*sphi
                  za(ia)=ya(ia)*sphi+za(ia)*cphi
                  ya(ia)=yasave
            end do
!      end if
      endif
      return
      end

      subroutine transm(im)

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
c
c move molecule
c
      call movem(im,wholetrans(mtype(im)))
c
c Calculate the change in energy
c
!      fdiff=0.0
!      fdiffb=0.0
!      call en_change(im,fdiff)

      return
      end

      subroutine movem(im,stepsz)

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      logical nochange

c store present coordinates
      xs=xadd(im)
      ys=yadd(im)
      zs=zadd(im)
c move atoms by random amounts
      zadd(im)=zadd(im)+stepsz*(2.*ran1()-1.0)
      yadd(im)=yadd(im)+stepsz*(2.*ran1()-1.0)
      xadd(im)=xadd(im)+stepsz*(2.*ran1()-1.0)
      nochange=xs.eq.xadd(im).and.ys.eq.yadd(im).and.zs.eq.zadd(im)
      do while (nochange)
          stepsz=stepsz*2.0
          zadd(im)=zadd(im)+stepsz*(2.*ran1()-1.0)
          yadd(im)=yadd(im)+stepsz*(2.*ran1()-1.0)
          xadd(im)=xadd(im)+stepsz*(2.*ran1()-1.0)
          nochange=xs.eq.xadd(im).and.ys.eq.yadd(im).and.zs.eq.zadd(im)
      end do
      return
      end

      subroutine en_change_intra(im,fdiffb)

      use ergpairreference
      
      real fdiffb,sumevib,sumeang,sumedih
      integer im

      sumevib=0.0
      sumeang=0.0
      sumedih=0.0
! Calculate the change in energy as a result of moving all the atoms in the molecule
      ergnew(1)=en_mole_intra(im,sumevib,sumeang,sumedih)
      fdiffb=fdiffb+ergnew(1)-ergold(1)
!      write(6,*) im,molergpair(1,im),ergpair(ergpairref(1,im)),ergold(1),ergnew(1),fdiffb

      return
      end                


      subroutine en_change(im,fdiff)
c
c calculates the change in energy as a result of moving all the atoms in a molecule
c
      use omp_lib
      use ergpairreference
      use timer

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      
      real,dimension(:), allocatable         :: x1temp,y1temp,z1temp
      real                                   :: x2temp,y2temp,z2temp,stepcl,sume
      integer                                :: icref,inn,it,jt,ia,ja,im,jm,iaadd
      
      if(.not.allocated(x1temp)) then
         allocate(x1temp(nato(im)),y1temp(nato(im)),z1temp(nato(im)))
      else if(nato(im).gt.size(x1temp)) then
         deallocate(x1temp,y1temp,z1temp)
         allocate(x1temp(nato(im)),y1temp(nato(im)),z1temp(nato(im)))
      end if     
!      write(6,*) 'en_change_whole> ',im
      stepcl=cellst/dmax
      call system_clock(time1)
c initialise the centre of mass coordinates and atom types which do not
c need to be updated for every atom - only when a molecule number changes
!Get new absolute coordinates of atoms on molecule im
      iaadd=naf(im)-1
      do iref=1,nato(im)
         ia=iref+iaadd
         x1temp(iref)=xa(ia)+xadd(im)
         y1temp(iref)=ya(ia)+yadd(im)
         z1temp(iref)=za(ia)+zadd(im)
      end do
!$omp parallel do private(inn)
      do inn=2,nnbref
         ergnew(inn)=0.0
      end do
!$omp end parallel do
!$omp parallel do private(inn,jm,sume,ja,x2temp,y2temp,z2temp,jt,iref,ia,it,icref) shared(ergnew)
!Step through list of neighbouring molecules
      do inn=2,nnbref
         jm=nbref(inn)
         !Don't bother if both molecules are tethered - this is dealt with separately
!         if(.not.(tethered(im).and.tethered(jm))) then
            sume=0.0
            do ja=naf(jm),nal(jm)
               x2temp=xa(ja)+xadd(jm)
               y2temp=ya(ja)+yadd(jm)
               z2temp=za(ja)+zadd(jm)
               jt=itype(ja)
!Step through moved atoms on molecule im
               do iref=1,nato(im)
                  ia=iref+iaadd
                  it=itype(ia)
                  icref=npairref(it,jt)
                  sume=sume+en_atom_inter_test(
     *x1temp(iref),y1temp(iref),z1temp(iref)
     *,x2temp,y2temp,z2temp,stepcl,xstar,dmax,icref,npt2)
               end do
            end do
            ergnew(inn)=ergnew(inn)+sume
!         endif
      end do
!$omp end parallel do
      do inn=2,nnbref
         fdiff=fdiff+ergnew(inn)-ergold(inn)
      end do
c Add on the effect of tethered molecules if present
      if(tethered(im)) call en_change_tethered(im,fdiff)
      call system_clock(time2)
      sumtime2=sumtime2+(time2-time1)/1000

      return
      end                

      function en_atom_intra(ia,ja,x1,y1,z1,x2,y2,z2,sumevib,sumeang)

c
c calculates the internal energy of atom pair ia,ja in a molecule
c
c it is assumed this pair have been preselected not to be the same atom
c
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'

      real en_atom_intra
      real x1,y1,z1,x2,y2,z2,xdif,ydif,zdif,dsq,dtest,dtdif,sume,stepcl
      real sumevib,sumeang
      integer ia,ja,it,jt,iz,icref,jref,image,ir,im,io1,io2
      logical notqtype,found

      en_atom_intra=0.0
!If ja is an image atom of ia, do nothing
      if(ja.eq.imageatom(ia)) return
      sume=0
c Calculate the internal energy for this atom pair
      stepcl=cellst/dmaxmol
c Set the atom type numbers
      it=itype(ia)
      jt=itype(ja)
      icref=npairref(it,jt)
      notqtype=qtype(ia).eq.0.and.qtype(ja).eq.0
c calculate separation of atoms within this molecule
      xdif=x2-x1
      ydif=y2-y1
      zdif=z2-z1
c first sum intramolecular harmonic potentials if needed
      jref=0
      iflag=0
      image=ia
      do while(jref.lt.nbond(image).and.iflag.eq.0)
         jref=jref+1
         if(ja.eq.iabond(image,jref)) iflag=jref
      end do
! If a bond was not found, then check the image atom for bonds
      if(iflag.eq.0) then
         image=0
         if(imageatom(ia).gt.0) then
            image=imageatom(ia)
            jref=0
            do while(jref.lt.nbond(image).and.iflag.eq.0)
               jref=jref+1
               if(ja.eq.iabond(image,jref)) iflag=jref
            end do
         end if
! Apply minimum image convention to these separations as they are likely to be across the cell boundary
         if(iflag.gt.0) call minimagec(xdif,ydif,zdif)
      end if
      if(iflag.gt.0.and.image.gt.0) then
         dsq=xdif*xdif+ydif*ydif+zdif*zdif
         dtest=sqrt(dsq)
c find the r value this corresponds to
         rat=dtest/stepcl
         iz=int(rat)+2
         wtms=amassr(icref)
         tb=abs(tbond(image,jref))
         wtm=wtms/tb
         if(tbond(image,jref).lt.0.0) wtm=wtm*tempang !Emphasize bond angles if requested.
         ddif=dtest-tb
         sume=wtm*ddif*ddif
         if(tbond(image,jref).gt.0.0) sumevib=sumevib+sume
         if(tbond(image,jref).lt.0.0) sumeang=sumeang+sume 
         if(refine_intra.and.iz.le.npt2.and.notqtype) then
            dtdif=dtest-xstar(iz)
            sume=sume+refine_intra_fac*(xavs(iz,icref)+xavs1(iz,icref)*dtdif)
         end if
      else
         !Check whether ja is listed as forming a dihedral angle with ia
         im=mref(ia)
         found=.false. 
! Step through the rotational groups for this molecule and pick out the dihedrals -
! Any atom pairs that are part of a dihedral are not used to calculate interatomic energies from the reference potential
         ir=nrf(im)-1
         do while(ir.lt.nrl(im).and..not.found)
            ir=ir+1
! Step through all the dihedral angles and calculate their contribution to the energy
            kd=0
            do while(kd.lt.ndihedralang(ir).and..not.found)
               kd=kd+1
               io1=labrotdih(1,kd,ir)
               io2=labrotdih(2,kd,ir)
               found=(io1.eq.ia.and.io2.eq.ja).or.(io1.eq.ja.and.io2.eq.ia)
            end do
         end do
! apply minimum image convention if atoms are not connected by bonds
         xdif=xdif*sizefactor
         ydif=ydif*sizefactor
         zdif=zdif*sizefactor
         call minimagec(xdif,ydif,zdif)
!         xdif=xdif/sizefactor
!         ydif=ydif/sizefactor
!         zdif=zdif/sizefactor
         dsq=xdif*xdif+ydif*ydif+zdif*zdif
         dtest=sqrt(dsq)
c find the r value this corresponds to
         rat=dtest/stepcl
         iz=int(rat)+2
         if(iz.le.npt2) then
            dtdif=dtest-xstar(iz)
c if no bonds were found for this atom pair apply molecular and empirical potential 
c in the usual way. However don't do the first of these things if there is a dihedral angle defined for this pair
c Q-type atoms are not included in this calculation since their positions are
c defined by reference to other atoms in the molecule.
            if(.not.found) then
               sume=xav(iz,icref)+xav1(iz,icref)*dtdif
            else
               sume=0.0
            end if
            if(refine_intra.and.notqtype) sume=sume
     1+refine_intra_fac*(xavs(iz,icref)+xavs1(iz,icref)*dtdif)

         end if
      end if
      en_atom_intra=sume
      return
      end

      function en_atom_inter(ia,ja,x1,y1,z1,x2,y2,z2)
c
c calculates the intermolecular energy between atoms ia and ja on distinct
c molecules im and jm. It is assumed the atoms have been preselected to occur 
c on distinct molecules.
c
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'epsrinput_par.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'

      real en_atom_inter,x1,y1,z1,x2,y2,z2,xdif,ydif,zdif,dsq,dtest,rat,dtdif,stepcl
      integer ia,ja,it,jt,iz,icref

      sume=0.0
      stepcl=cellst/dmax
!      write(6,*) ia,ja,x1,y1,z1,x2,y2,z2,stepcl
c Set the atom type numbers
      it=itype(ia)
      jt=itype(ja)
      icref=npairref(it,jt)
!      write(6,*) ia,ja,it,jt,icref
      xdif=x2-x1
      ydif=y2-y1
      zdif=z2-z1
c apply minimum image convention
      call minimagec(xdif,ydif,zdif)
      dsq=xdif*xdif+ydif*ydif+zdif*zdif
      dtest=sqrt(dsq)
!      write(6,*) ia,ja,dtest
c find the r value this corresponds to
      rat=dtest/stepcl
      iz=int(rat)+2
!      write(6,*) ia,ja,it,jt,iz,npt2
      if(iz.le.npt2) then
c sum the reference and empirical potentials
         dtdif=dtest-xstar(iz)
         sume=xavtot(iz,icref)+xavtot1(iz,icref)*dtdif
      else
         sume=0.0
      endif
      irefp=irefpcof(icref)
      if(rminex(irefp).lt.0.0.and.pcoefex(irefp).gt.0.0) then
         rtest=abs(rminex(irefp))*(1.0-expecf)
         dtesttest=dtest*dmax
         ddiff=dtesttest-rtest
         if(ddiff.gt.0.0) sume=sume+pcoefex(irefp)*ddiff*ddiff
      end if
      
      en_atom_inter=sume
      return
      end

      function en_atom_inter_test(x1,y1,z1,x2,y2,z2
     *,stepcl,xstar,dmax,icref,npt2)

      include 'epsr_dimension.inc'
      include 'epsrinput_par.inc'
      include 'potential_par.inc'
c
c calculates the intermolecular energy between atoms ia and ja on distinct
c molecules im and jm. It is assumed the atoms have been preselected to occur 
c on distinct molecules.
c
      real xstar(*),rtest
      real en_atom_inter_test,x1,y1,z1,x2,y2,z2
     *,xdif,ydif,zdif,dsq,dtest,rat,dtdif,stepcl
     *,dmax
      integer iz,icref,irefp

      sume=0.0
      xdif=x2-x1
      ydif=y2-y1
      zdif=z2-z1
c apply minimum image convention
      call minimagec(xdif,ydif,zdif)
      dsq=xdif*xdif+ydif*ydif+zdif*zdif
      dtest=sqrt(dsq)
c find the r value this corresponds to
      rat=dtest/stepcl
      iz=int(rat)+2
      if(iz.le.npt2) then
         dtdif=dtest-xstar(iz)
c sum the reference and empirical potentials
         sume=xavtot(iz,icref)+xavtot1(iz,icref)*dtdif
      else
         sume=0.0
      endif
      irefp=irefpcof(icref)
      if(rminex(irefp).lt.0.0.and.pcoefex(irefp).gt.0.0) then
         rtest=abs(rminex(irefp))*(1.0-expecf)
         dtesttest=dtest*dmax
         ddiff=dtesttest-rtest
         if(ddiff.gt.0.0) sume=sume+pcoefex(irefp)*ddiff*ddiff
      end if
      
      en_atom_inter_test=sume
      return
      end

      function en_mole_intra(im,sumevib,sumeang,sumedih)

!Calculates the intra-molecular energy of a molecule using
!the current coordinates.      

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      real sume
      real sumevib,sumeang,sumedih

      real en_mole_intra,en_atom_intra,en_mole_dihedral

      en_mole_intra=0.0
      ia=naf(im)-1
      do while(ia.lt.nal(im)-1)
         ia=ia+1
         ja=ia
         do while(ja.lt.nal(im))
            ja=ja+1
            sume=en_atom_intra(ia,ja,xa(ia),ya(ia),za(ia),xa(ja),ya(ja),za(ja),sumevib,sumeang)
            en_mole_intra=en_mole_intra+sume
         end do
      end do
      if(ngrot(im).gt.0.and.tempdih.gt.0.0) then
          sume=tempdih*simt*en_mole_dihedral(im)
          en_mole_intra=en_mole_intra+sume
          sumedih=sumedih+sume
      end if
      return
      end

      function en_mole_dihedral(im)

		use molecom_mod

!Calculates the intra- or inter- molecular energy of a pair of molecules using
!the current coordinates.      

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      real en_mole_dihedral
      real cd(3),co1(3),co2(3),cv1(3),cv2(3),cpar(3)
      integer im,ir,ia1,ia2,kd,io1,io2

      en_mole_dihedral=0.0
! Step through the rotational groups for this molecule and pick out the dihedrals
      do ir=nrf(im),nrl(im)
         if(rotgrp(ir).eq.'DIH') then
c
c define the separation vector cd between the two base atoms and normalize it to 
c unit value
c
            ia1=labaxi(ir,1)
            ia2=labaxi(ir,2)
            x1=xa(ia1)
            y1=ya(ia1)
            z1=za(ia1)
            x2=xa(ia2)
            y2=ya(ia2)
            z2=za(ia2)
            cd(1)=(x2-x1)
            cd(2)=(y2-y1)
            cd(3)=(z2-z1)
            dmag=sqrt(scapro(cd,cd))
            cd(1)=cd(1)/dmag
            cd(2)=cd(2)/dmag
            cd(3)=cd(3)/dmag
!            write(6,*) ir,ia1,ia2,ndihedralang(ir)
! Step through all the dihedral angles and calculate their contribution to the energy
            kd=0
            do while(kd.lt.ndihedralang(ir))
               kd=kd+1
c
c now determine the orientation vectors co for base atom 1 and base atom 2, using 
c the values of ltyp2 and ltyp3 
c
               io1=labrotdih(1,kd,ir)
               io2=labrotdih(2,kd,ir)
!               write(6,*) kd,io1,io2
               co1(1)=xa(io1)-x1
               co1(2)=ya(io1)-y1
               co1(3)=za(io1)-z1
               co2(1)=xa(io2)-x2
               co2(2)=ya(io2)-y2
               co2(3)=za(io2)-z2
c
c define new vectors normal to both cd and co
c
               call vecpro(co1,cd,cv1)
               call vecpro(co2,cd,cv2)
c
c renormalize these to unit vectors
c
               dmag=sqrt(scapro(cv1,cv1))
               cv1(1)=cv1(1)/dmag
               cv1(2)=cv1(2)/dmag
               cv1(3)=cv1(3)/dmag
               dmag=sqrt(scapro(cv2,cv2))
               cv2(1)=cv2(1)/dmag
               cv2(2)=cv2(2)/dmag
               cv2(3)=cv2(3)/dmag
c calculate the cosine of the angle between cv1 and cv2
               cosphi=scapro(cv1,cv2)
               if(cosphi.gt.1.0) cosphi=1.0
               if(cosphi.lt.-1.0) cosphi=-1.0
c
c calculate the vector product of cv1 and cv2 - this should be either parallel
c (phi positive) or antiparallel to cd (phi negative)
c
               call vecpro(cv1,cv2,cpar)
               sinphi=scapro(cpar,cd)
               cosdih=cos(dihedralang(kd,ir))
               sindih=sin(dihedralang(kd,ir))
               difc=cosdih-cosphi
               difs=sindih-sinphi
               diff=difc*difc+difs*difs
!      write(6,*) ir,kd,cosdih,cosphi,sindih,sinphi
               en_mole_dihedral=en_mole_dihedral+diff
            end do
         end if
      end do
      return

      end function en_mole_dihedral

      function en_mole_pair(im,jm)

      use omp_lib
      
!Calculates the inter- molecular energy of a pair of molecules using
!the current coordinates.      

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      real en_mole_pair,en_atom_inter,x1,y1,z1,x2,y2,z2,sume
      integer inn,ia,ja,im,jm
      

      en_mole_pair=0.0
!!$omp parallel do private(inn,ia,x1,y1,z1,ja,x2,y2,z2,sume) shared(en_mole_pair)
      do inn=1,1
      !Don't calculate pair energies if both molecules are tethered
!      if(.not.(tethered(im).and.tethered(jm))) then
         sume=0.0
         do ia=naf(im),nal(im)
            x1=xa(ia)+xadd(im)
            y1=ya(ia)+yadd(im)
            z1=za(ia)+zadd(im)
            do ja=naf(jm),nal(jm)
               x2=xa(ja)+xadd(jm)
               y2=ya(ja)+yadd(jm)
               z2=za(ja)+zadd(jm)
               sume=sume+en_atom_inter(ia,ja,x1,y1,z1,x2,y2,z2)
            end do
         end do   
         en_mole_pair=en_mole_pair+sume
!         write(6,*) omp_get_thread_num(),inn,im,jm,sume
!      end if
      end do
!!$omp end parallel do

      return
      end


