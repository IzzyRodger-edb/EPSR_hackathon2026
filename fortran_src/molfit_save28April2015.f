      subroutine molfit(fadds,faddbs,nisum,intra,rotfreq,ireset)
      
      use timer
      use cellneighbourlist
      use ergpairreference
      use molecom_mod
      use threebody_mod
                    
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      real, dimension(:), allocatable         :: xr,yr,zr
      real stepmithis,steprithis,en_change_intra
      integer iaadd,nisum,intra,rotfreq
      
      temp=1.0/simt
      tol=tolout/temp
!      tol=tolout/dmax
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
      rotstepmin=1.0e-5
      fadd=0.0
      faddb=0.0
! Ensure correct weights on angle and dihedral terms are set
      tempang=tempangwt
      tempdih=tempdihwt*tempb
c
c if step size is non-zero, then do internal molecule structure refinement. This is to 
c prevent the build up of round-off errors from allowing the molecule to 
c become too distorted.
c
      if(nisum.ge.intra) then

         do im=1,nmol
            call make_neighbour_list_from_mols(im)
            if(.not.allocated(xr)) then
                allocate(xr(nato(im)),yr(nato(im)),zr(nato(im)))
            else if (size(xr).lt.nato(im)) then
                deallocate(xr,yr,zr)
                allocate(xr(nato(im)),yr(nato(im)),zr(nato(im)))
            end if
            if(nato(im).gt.1.and.stepmi.gt.0.0) then
!         write(6,*) im,nergpairm(im),molergpair(1,im),ergpair(ergpairref(1,im))
!         write(6,*) im,nnbref,ergold(1),ergnew(1)
               ! first save the coordinates of this molecule
               stepmithis=interntrans(mtype(im))
               if(stepmi.gt.0.0.and.stepmithis.eq.0.0) stepmithis=trstepmin
               interntrans(mtype(im))=stepmithis
               xs=xadd(im)
               ys=yadd(im)
               zs=zadd(im)
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
!               write(6,*) 'molfit> 83 ',im,size(ergold),size(ergnew)
               fdiffb=en_change_intra(im,sumevib,sumeang,sumedih)
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
!               if(im.lt.11) write(6,*) 'molfit-trans> Molecular energy',im,ergpair(ergpairref(1,im))
!            else
!               interntrans(mtype(im))=0.0
            endif
!            if(im.lt.11) write(6,*) 'molfit-rot> ',im,stepri,ngrot(im),nrf(im),nrl(im)
! Now do the internal rotations
            if(stepri.gt.0.0.and.ngrot(im).gt.0) then
               if(internrot(mtype(im)).eq.0.0) internrot(mtype(im))=rotstepmin
c
c perform internal rotations
c
                  ninter(1,mtype(im))=ninter(1,mtype(im))+1
! Save the coordinates of this group
                  nat=naf(im)-1
                  do ia=1,nato(im)
                     nat=nat+1
                     xr(ia)=xa(nat)
                     yr(ia)=ya(nat)
                     zr(ia)=za(nat)
                  end do
                  fdiff=0.0
                  fdiffb=0.0
                  do ig=nrf(im),nrl(im)
                     call rotati(im,ig,internrot(mtype(im)),actualstep)
                  end do
!                  if(tethered(im)) call check_tethered(im,1) !Make sure move is allowed
c calculate change in intramolecular energy
!               write(6,*) 'molfit> 145 ',size(ergold),size(ergnew)
                  fdiffb=en_change_intra(im,sumevib,sumeang,sumedih)
                  fdiff=en_change(im)
!                  if(tethered(im)) fdiff=fdiff+en_change_tethered(im,tetheredatom(im))
!         if(im.lt.11) write(6,*) 'molfit> Vibrational, bond angle and dihedral energies: ',im,fdiffb,ergnew(1),ergold(1),sumevib,sumeang,sumedih
                  ediff=(fdiff+fdiffb)*temp
                  if(ediff.le.0.0) then
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
                     ergold_triplet(im)=ergnew_triplet(im)
                  else
! Move rejected - restore coordinates
                     nstepreject(im)=nstepreject(im)+1
! Restore the coordinates of this group
                     nat=naf(im)-1
                     do ia=1,nato(im)
                        nat=nat+1
                        xa(nat)=xr(ia)
                        ya(nat)=yr(ia)
                        za(nat)=zr(ia)
                     end do
                     ninter(2,mtype(im))=ninter(2,mtype(im))+1
                  end if
!                     write(6,*) 'molfit-rot> ',im,fadd,faddb,ninter(1,mtype(im)),ninter(2,mtype(im))
! Save all the relevant pair energies
               call save_erg_pairs(im)
!               if(im.lt.11) write(6,*) 'molfit-rot> Molecular energy',im,ergpair(ergpairref(1,im)),fdiffb
! Adjust the coordinates of centre of mass to correct for any change
               call centreofmass(im)
            else
               if(ngrot(im).eq.0) internrot(mtype(im))=0.0
            end if
!            if(im.lt.11) write(6,*) 'molfit> Vibrational, bond angle and dihedral energies: ',im,sumevib,sumeang,sumedih
         end do
      else
!Step through cells and molecules 
         do ic=1,ncell
!            write(6,*) 'molfit-rot> Cell ',ic,fadd,faddb
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
!                  open(25,file='ergold0.txt',status='unknown')
!                  do inn=1,nnbref
!                     write(25,*) im,nbref(inn),inn,ig,ergnew(inn),ergold(inn),ergnew(inn)-ergold(inn),fadd,faddb
!                  end do
!                  close(25)
                  if(stepri.gt.0.0.and.ngrot(im).gt.0.and.nisum.eq.rotfreq*int(nisum/rotfreq)) then !Only do this every rotfreq iterations
                     if(internrot(mtype(im)).eq.0.0) internrot(mtype(im))=rotstepmin
c
c perform internal rotations
c
!                        open(25,file='ergold1.txt',status='unknown')
!                        do inn=1,nnbref
!                           write(25,*) im,nbref(inn),inn,ig,ergnew(inn),ergold(inn),ergnew(inn)-ergold(inn),fadd,faddb
!                        end do
!                        close(25)
                        ninter(1,mtype(im))=ninter(1,mtype(im))+1
! Save the coordinates of this group
                        nat=naf(im)-1
                        do ia=1,nato(im)
                           nat=nat+1
                           xr(ia)=xa(nat)
                           yr(ia)=ya(nat)
                           zr(ia)=za(nat)
                        end do
                        fdiff=0.0
                        fdiffb=0.0
                        do ig=nrf(im),nrl(im)
                           call rotati(im,ig,internrot(mtype(im)),actualstep)
                        end do
!                        if(tethered(im)) call check_tethered(im,1) !Make sure move is allowed
 
c calculate change in intramolecular energy
c
!               write(6,*) 'molfit> 236 ',size(ergold),size(ergnew)
                        fdiffb=en_change_intra(im,sumevib,sumeang,sumedih)
                        fdiff=en_change(im)
!                        if(tethered(im)) fdiff=fdiff+en_change_tethered(im,tetheredatom(im))
!                        open(25,file='ergold2.txt',status='unknown')
!                        do inn=1,nnbref
!                           write(25,*) im,nbref(inn),inn,ig,ergnew(inn),ergold(inn),ergnew(inn)-ergold(inn),fadd,fdiff
!                        end do
!                        close(25)
                        ediff=(fdiff+fdiffb)*temp
!                        if(im.lt.11) then
!                           write(6,*) 'molfit-rot> ',im,ig,internrot(mtype(im)),fadd,faddb
!                           write(6,*) 'molfit-rot> ',fdiff,fdiffb,ediff,ninter(1,mtype(im)),ninter(2,mtype(im))
!                        end if
                        if(ediff.le.0.0) then
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
                           ergold_triplet(im)=ergnew_triplet(im)
                        else
! Move rejected - restore coordinates
                           nstepreject(im)=nstepreject(im)+1
! Restore the coordinates of this group
                           nat=naf(im)-1
                           do ia=1,nato(im)
                              nat=nat+1
                              xa(nat)=xr(ia)
                              ya(nat)=yr(ia)
                              za(nat)=zr(ia)
                           end do
                           ninter(2,mtype(im))=ninter(2,mtype(im))+1
                        end if
!                        open(25,file='ergold3.txt',status='unknown')
!                        do inn=1,nnbref
!                           write(25,*) im,nbref(inn),inn,ig,ergnew(inn),ergold(inn),ergnew(inn)-ergold(inn),fadd,fdiff
!                        end do
!                        close(25)
!                        write(6,'(a,i5,1x)',advance='no') 'molfit-rot3> Type 1 to continue, 0 to stop ', ig
!                        read(5,*) ngo
!                        if(ngo.eq.0) stop
!                     write(6,*) 'molfit-rot> ',im,fadd,faddb,ninter(1,mtype(im)),ninter(2,mtype(im))
! Save all the relevant pair energies
                     call save_erg_pairs(im)
! Adjust the coordinates of centre of mass to correct for any change
                     call centreofmass(im)
                  else
                     if(ngrot(im).eq.0) internrot(mtype(im))=0.0
                  end if
! Rotate and translate whole molecules and calculate energy change. Only rotate molecules that have more than one atom
! Rotate molecules
                  if(stepr.gt.0.0.and.nato(im).gt.1) then
! First save the coordinates of this molecule
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
!                     if(tethered(im)) call check_tethered(im,1) !Make sure move is allowed
! Calculate change in intramolecular and intermolecular energy
!               write(6,*) 'molfit> 308 ',size(ergold),size(ergnew)
                     fdiffb=en_change_intra(im,sumevib,sumeang,sumedih)
                     fdiff=en_change(im)
!                     if(tethered(im)) fdiff=fdiff+en_change_tethered(im,tetheredatom(im))
                        open(25,file='ergold4.txt',status='unknown')
!                        do inn=1,nnbref
!                           write(25,*) im,nbref(inn),inn,ergnew(inn),ergold(inn),ergnew(inn)-ergold(inn),fadd,fdiff
!                        end do
!                        close(25)
!                        write(6,'(a,i5,1x)',advance='no') 'molfit-wholerot> Type 1 to continue, 0 to stop ', im
!                        read(5,*) ngo
!                        if(ngo.eq.0) stop
                     ediff=(fdiff+fdiffb)*temp
                     if(ediff.le.0.0) then
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
                        ergold_triplet(im)=ergnew_triplet(im)
                     else
! Move rejected - restore coordinates
                        nstepreject(im)=nstepreject(im)+1
                        nat=naf(im)-1
                        do ia=1,nato(im)
                           nat=nat+1
                           xa(nat)=xr(ia)
                           ya(nat)=yr(ia)
                           za(nat)=zr(ia)
                        end do
                        nwhole(2,mtype(im))=nwhole(2,mtype(im))+1
                     end if
!                     write(6,*) 'molfit-wholerot> ',im,fadd,faddb,nwhole(1,mtype(im)),nwhole(2,mtype(im))
! Adjust the coordinates of centre of mass to correct for any change
                     call centreofmass(im)
! Save all the relevant pair energies
                     call save_erg_pairs(im)
                  else
                     if(nato(im).lt.2) wholerot(mtype(im))=0.0
                  end if
! Translate molecules
                  if(stepm.gt.0.0) then
! First save the coordinates of this molecule
                     xs=xadd(im)
                     ys=yadd(im)
                     zs=zadd(im)
                     fdiff=0.0
                     fdiffb=0.0
                     ! now translate each molecule whole
                     if(wholetrans(mtype(im)).le.0.0) wholetrans(mtype(im))=trstepmin
                     nmovem(1,mtype(im))=nmovem(1,mtype(im))+1
                     call movem(im,wholetrans(mtype(im)))
!                     if(tethered(im)) call check_tethered(im,1) !Make sure move is allowed
! Test tethered molecules if present
                     etest=10.0
                     if(tethered(im)) then
                        ediff=en_change_tethered(im,tetheredatom(im))*temp
                        if(ediff.le.0.0) then
                           etest=10.0
                        else if(ediff.lt.30.0) then
                           etest=exp(-ediff)
                        else
                           etest=0.0
                        end if
!      write(6,*) im,ediff,tetheredatom(im),temp,etest
                     end if
                     if(ran1().lt.etest) then
! OK now check intra- and inter-molecular energies
!               write(6,*) 'molfit> 366',366,size(ergold),size(ergnew)
                        if(nato(im).gt.1) fdiffb=en_change_intra(im,sumevib,sumeang,sumedih)
                        fdiff=en_change(im)
                        ediff=(fdiff+fdiffb)*temp
                        if(ediff.le.0.0) then
                           etest=10.0
                        else if(ediff.lt.30.0) then
                           etest=exp(-ediff)
                        else
                           etest=0.
                        endif
!      write(6,*) im,fdiffb,fdiff,ediff,etest
                        if(ran1().lt.etest) then
                           nstepreject(im)=0
                           fadd=fadd+fdiff
                           faddb=faddb+fdiffb
                           call save_energies()
                           ergold_triplet(im)=ergnew_triplet(im)
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
! Move rejected - restore coordinates
                        nstepreject(im)=nstepreject(im)+1
                        xadd(im)=xs
                        yadd(im)=ys
                        zadd(im)=zs
                        nmovem(2,mtype(im))=nmovem(2,mtype(im))+1
                     end if
                  end if
               end do
            end if
         end do
! Set the weighting on the dihedral energy
!         sumevib=0.0
!         sumeang=0.0
!         sumedih=0.0
!         tempang=1.0
!         tempdih=tempb
!         do im=1,nmol
!         !Get the intramolecular energy of this molecule
!            energy=en_mole_intra(im,sumevib,sumeang,sumedih)
!         write(6,*) im,nergpairm(im)
!         end do
!         write(6,*) 'molfit> Vibrational, bond angle, and dihedral energies: ',sumevib,sumeang,sumedih
!         if(sumevib.gt.0.0) then
!            ratioang=sumeang/sumevib
!            ratiodih=sumedih/sumevib
!         else
!            ratioang=1.0
!            ratiodih=1.0
!         end if
!         write(6,*) 'molfit> Relative bond angle and dihedral energies: ',ratioang,ratiodih
!         ratiomax=1.01
!         ratiomin=0.99
!         tempdihwt=tempdihreq
!         tempang=tempangwt
!         tempdih=tempdihwt*tempb
!         write(6,*) 'molfit> Bond angle and dihedral energy weights: ',tempang,tempdih
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
!      do ig=nrf(im),nrl(im)
!         call rotati(ig,internrot(mtype(im)))
!      end do
c
c calculate change in energy
c
!      fdiff=0.0
!      fdiffb=0.0
!      call en_change_intra(im,fdiffb)
!      call en_change(im,fdiff)
      return
      end

      subroutine rotati(im,ig,steprithis,actualstep)

      use molecom_mod
      
      include 'epsr_dimension.inc'
      include 'molecom.inc'

c pvector is the vector about which the rotation occurs

      real pvector(3),steprithis,actualstep

c xvector is the vector from centre of rotation to an individual atom
c xpvector is the component of this vector along the axis of rotation
c xrvector is the component of this vector perpendicular to the axis of rotation
c xrpvector is the vector perpendicular to both axis of rotation and xrvector

      real xvector(3),xpvector(3),xrvector(3),xrpvector(3)
      integer im,ig,i,ir,nat,ia,nax

      pi=3.141592653
c define random new orientation
      phi=pi*steprithis*(2*ran1()-1)
      actualstep=phi
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

      function en_change_intra(im,sumevib,sumeang,sumedih)

      use ergpairreference
      
      real en_change_intra
      integer im

      sumevib=0.0
      sumeang=0.0
      sumedih=0.0
! Calculate the change in energy as a result of moving all the atoms in the molecule
      ergnew(1)=en_mole_intra(im,sumevib,sumeang,sumedih)
      en_change_intra=ergnew(1)-ergold(1)
!      write(6,*) im,molergpair(1,im),ergpair(ergpairref(1,im)),ergold(1),ergnew(1),fdiffb

      return
      end                

      function en_change_intra_rot(im,ig,sumevib,sumeang,sumedih)

      use ergpairreference
      use molecom_mod

!Calculates the intra-molecular energy of a molecule using
!the current coordinates.      

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      real sumevib,sumeang,sumedih,en_change_intra_rot
      real en_mole_intra
      integer im,ig

      sumevib=0.0
      sumeang=0.0
      sumedih=0.0
! Calculate the change in energy as a result of moving the atoms in rotational group ig in the molecule
! The safest way - may not be the quickest - is simply calculate the full intramolecular energy of the revised molecule
      ergnew(1)=en_mole_intra(im,sumevib,sumeang,sumedih)
      en_change_intra_rot=ergnew(1)-ergold(1)
!      write(6,*) im,molergpair(1,im),ergpair(ergpairref(1,im)),ergold(1),ergnew(1),fdiffb
      return
      end


      function en_change(im)
c
c calculates the change in energy as a result of moving all the atoms in a molecule
c
      use omp_lib
      use ergpairreference
      use timer
      use threebody_mod
      use containerforces_mod

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'epsrinput_par.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      
      real                                   :: x1temp,y1temp,z1temp
      real                                   :: x2temp,y2temp,z2temp,stepcl,sume,fdiff
      integer                                :: icref,inn,it,jt,ia,ja,im,jm,iaadd
      logical bothtethered,bothqatoms,qtypeia,qtypeja
      
      stepcl=cellst/dmax
      call system_clock(time1)
!$omp parallel do private(inn)
      do inn=2,nnbref
         ergnew(inn)=0.0
      end do
!$omp end parallel do
!$omp parallel do private(inn,jm,sume,ia,ja,x1temp,y1temp,z1temp,x2temp,y2temp,z2temp,it,jt,icref,bothtethered,bothqatoms,qtupeia,qtypeja) shared(ergnew,nbref,xa,ya,za,xadd,yadd,zadd,itype,npairref)
!Step through list of neighbouring molecules
      do inn=2,nnbref
         jm=nbref(inn)
         !Don't bother if both molecules are tethered and both atoms are Q-type.
         bothtethered=tethered(im).and.tethered(jm)
         sume=0.0
!Step through moved atoms on molecule im
         do ia=naf(im),nal(im)
            it=itype(ia)
            x1temp=xa(ia)+xadd(im)
            y1temp=ya(ia)+yadd(im)
            z1temp=za(ia)+zadd(im)
            qtypeia=qtype(ia).eq.1.and.bothtethered
            do ja=naf(jm),nal(jm)
               x2temp=xa(ja)+xadd(jm)
               y2temp=ya(ja)+yadd(jm)
               z2temp=za(ja)+zadd(jm)
               jt=itype(ja)
               bothqatoms=qtype(ja).eq.1.and.qtypeia
               !Don't calculate energy if both molecules are tethered and both atoms are Q types
               if(.not.bothqatoms) then
                  icref=npairref(it,jt)
                  sume=sume+en_atom_inter_test(
     *x1temp,y1temp,z1temp,x2temp,y2temp,z2temp,stepcl,xstar,dmax,icref,npt2)
               end if
            end do
         end do
         ergnew(inn)=ergnew(inn)+sume
      end do
!$omp end parallel do
      fdiff=0.0
      do inn=2,nnbref
         fdiff=fdiff+ergnew(inn)-ergold(inn)
      end do
! Get the new triplet energy
      ergnew_triplet(im)=en_mole_triplets(im,xadd,yadd,zadd,xa,ya,za,mref)
      fdiff=fdiff+ergnew_triplet(im)-ergold_triplet(im)
! Get the new container force energy if requested
      if(numcontainerforces.gt.0) then
         ergnew_container(im)=en_containerforce(im,naf(im),nal(im),xadd,yadd,zadd,xa,ya,za
     1,itype,mref,npairref,irefpcof,rmaxex,pcoefexmax,expecf,expecfpwr,dmax,harmonic)
!         write(6,*) im,ergnew_container(im)
         fdiff=fdiff+ergnew_container(im)-ergold_container(im)
      end if
      en_change=fdiff
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

      real en_atom_intra,en_morse
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
         ddif=tb-dtest
         sume=en_morse(ddif,sqrt(wtm),tempang)*tempb*simt
!         sume=wtm*ddif*ddif*tempb*simt
         if(tbond(image,jref).gt.0.0) sumevib=sumevib+sume
         if(tbond(image,jref).lt.0.0) sumeang=sumeang+sume !Emphasize bond angles if requested.
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

      function en_morse(ddif,wtm,distort)
      real en_morse,ddif,wtm,distort,arg
! distort controls the amount of assymmetry in the bond. Larger means more symmetric, smaller less symmetric
      arg=wtm*ddif*distort
!      write(6,*) ddif,wtm,arg
      en_morse=(exp(arg)-1.0)/distort
      en_morse=en_morse*en_morse
      return
      end function en_morse

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
      do ia=naf(im),nal(im)-1
         do ja=ia+1,nal(im)
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
      real cd(3),co1(3),co2(3)
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
               en_mole_dihedral=en_mole_dihedral+en_dihedral(co1,co2,cd,dihedralang(kd,ir))
            end do
         end if
      end do
      return

      end function en_mole_dihedral

      function en_change_dihedral(im,ig,xr,yr,zr)

      use molecom_mod

!Calculates the intra- or inter- molecular energy of a pair of molecules using
!the current coordinates.      

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      real en_change_dihedral
      real xr(*),yr(*),zr(*)
      real cd(3),co1(3),co2(3),cv1(3),cv2(3),cpar(3)
      integer im,ig,ir,ia1,ia2,kd,io1,io2,iaadd

      en_change_dihedral=0.0
! Step through the rotational groups for this molecule and pick out the dihedrals
      ir=ig
      iaadd=1-naf(im)
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
!        write(6,*) ir,ia1,ia2,ndihedralang(ir)
! Step through all the dihedral angles and calculate their contribution to the energy
         kd=0
         do while(kd.lt.ndihedralang(ir))
            kd=kd+1
c
c now determine the orientation vectors co for base atom 1 and base atom 2, using 
c the values of ltyp2 and ltyp3 
c
! Energy of current dihedral
            io1=labrotdih(1,kd,ir)
            io2=labrotdih(2,kd,ir)
            co1(1)=xa(io1)-x1
            co1(2)=ya(io1)-y1
            co1(3)=za(io1)-z1
            co2(1)=xa(io2)-x2
            co2(2)=ya(io2)-y2
            co2(3)=za(io2)-z2
            en_change_dihedral=en_change_dihedral+en_dihedral(co1,co2,cd,dihedralang(kd,ir))
! Energy of old dihedral
            ir1=io1+iaadd
            ir2=io2+iaadd
            co1(1)=xr(ir1)-x1
            co1(2)=yr(ir1)-y1
            co1(3)=zr(ir1)-z1
            co2(1)=xr(ir2)-x2
            co2(2)=yr(ir2)-y2
            co2(3)=zr(ir2)-z2
            en_change_dihedral=en_change_dihedral-en_dihedral(co1,co2,cd,dihedralang(kd,ir))
         end do
      end if

      return
      end function en_change_dihedral

      function en_dihedral(co1,co2,cd,dihedralang)

! Calculates dihedral energy for two vectors, co1 and co2, with respect to the axis vector cd
! and specified angel dihedralang.

      real en_dihedral,dihedralang
      real co1(3),co2(3),cd(3),cv1(3),cv2(3),cpar(3)

c define new vectors normal to both cd and co
      call vecpro(co1,cd,cv1)
      call vecpro(co2,cd,cv2)
c renormalize these to unit vectors
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
c calculate the vector product of cv1 and cv2 - this should be either parallel
c (phi positive) or antiparallel to cd (phi negative)
      call vecpro(cv1,cv2,cpar)
      sinphi=scapro(cpar,cd)
      cosdih=cos(dihedralang)
      sindih=sin(dihedralang)
      difc=cosdih-cosphi
      difs=sindih-sinphi
      diff=difc*difc+difs*difs
!      write(6,*) ir,kd,cosdih,cosphi,sindih,sinphi
      en_dihedral=diff

      end function en_dihedral

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
      logical bothtethered,bothqatoms,qtypeia
      

      en_mole_pair=0.0
      !Don't calculate pair energies if both molecules are tethered and both atoms are Q-type
      bothtethered=tethered(im).and.tethered(jm)
!!$omp parallel do private(inn,ia,x1,y1,z1,ja,x2,y2,z2,sume,qtypeia,bothqatoms) shared(en_mole_pair,bothtethered)
      do inn=1,1
         sume=0.0
         do ia=naf(im),nal(im)
            x1=xa(ia)+xadd(im)
            y1=ya(ia)+yadd(im)
            z1=za(ia)+zadd(im)
            qtypeia=qtype(ia).eq.1.and.bothtethered
            do ja=naf(jm),nal(jm)
               bothqatoms=qtype(ja).eq.1.and.qtypeia
               if(.not.bothqatoms) then
                  x2=xa(ja)+xadd(jm)
                  y2=ya(ja)+yadd(jm)
                  z2=za(ja)+zadd(jm)
                  sume=sume+en_atom_inter(ia,ja,x1,y1,z1,x2,y2,z2)
               end if
            end do
         end do   
         en_mole_pair=en_mole_pair+sume
!         write(6,*) omp_get_thread_num(),inn,im,jm,sume
      end do
!!$omp end parallel do

      return
      end

      function en_change_tethered(im,ia)

      use molecom_mod

c With tethered molecules, for whole molecule moves or moves where part of a molecule is rotated
c we need
c to calculate the change in energy that derives from the change in position of centre of mass

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      real en_change_tethered,xdif,ydif,zdif,term
      integer ia,iaref

c Get the old energy of this molecule relative to the tether position

      xdif=xs-phix(im)
      ydif=ys-phiy(im)
      zdif=zs-phiz(im)
      if(ia.gt.0.and.ia.le.nato(im)) then
         iaref=naf(im)+ia-1
         xdif=xdif+xa(iaref)
         ydif=ydif+ya(iaref)
         zdif=zdif+za(iaref)
      end if
      call minimagec(xdif,ydif,zdif)
      term=xdif*xdif+ydif*ydif+zdif*zdif

c Get the new position of this molecule relative to the tether position

      xdif=xadd(im)-phix(im)
      ydif=yadd(im)-phiy(im)
      zdif=zadd(im)-phiz(im)
      if(ia.gt.0.and.ia.le.nato(im)) then
         iaref=naf(im)+ia-1
         xdif=xdif+xa(iaref)
         ydif=ydif+ya(iaref)
         zdif=zdif+za(iaref)
      end if

c Calculate the change in energy 

      call minimagec(xdif,ydif,zdif)
      en_change_tethered=tol*(xdif*xdif+ydif*ydif+zdif*zdif-term)
      return    
      end

      subroutine check_tethered(im,ia)

      use molecom_mod

! With tethered molecules, for whole molecule moves or moves where part of a molecule is rotated
! we need to check that the molecule has not moved more than the allowed
! distance from the tethering point.
! to calculate the change in energy that derives from the change in position of centre of mass

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      real xdif,ydif,zdif,dist,ratio
      integer ia,iaref

! Get the position of this molecule relative to the tether position

      iaref=naf(im)+ia-1
      xdif=xadd(im)-phix(im)
      ydif=yadd(im)-phiy(im)
      zdif=zadd(im)-phiz(im)
      if(ia.gt.0.and.ia.le.nato(im)) then
         iaref=naf(im)+ia-1
         xdif=xdif+xa(iaref)
         ydif=ydif+ya(iaref)
         zdif=zdif+za(iaref)
      end if
      call minimagec(xdif,ydif,zdif)
      dist=sqrt(xdif*xdif+ydif*ydif+zdif*zdif)
      ratio=tol/dist
!      write(6,*) im,dist,tol,ratio,dmax
      if(ratio.lt.1.0) then !Bring molecule back!
         xadd(im)=phix(im)+ratio*xdif
         yadd(im)=phiy(im)+ratio*ydif
         zadd(im)=phiz(im)+ratio*zdif
         if(ia.gt.0.and.ia.le.nato(im)) then
            iaref=naf(im)+ia-1
            xadd(im)=xadd(im)-xa(iaref)
            yadd(im)=yadd(im)-ya(iaref)
            zadd(im)=zadd(im)-za(iaref)
         end if
      end if
      return    
      end


