      subroutine molfit_epsr20(fadds,faddbs,nisum,intra,ireset)
      
      use timer
      use cellneighbourlist
      use ergpairreference
                    
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      temp=1.0/simt
      tol=tolout/temp
      do j=1,2
         nmovei(j)=0
         ninter(j)=0
         do imoltype=1,nmoltype
             nwhole(j,imoltype)=0
             nmovem(j,imoltype)=0
         end do
      end do

c Set the maximum translational step size

      trstepmax=1.0/real(natomr)
      trstepmax=trstepmax**(1.0/3.0)
      trstepmin=0.0001*trstepmax
      rotstepmax=1.0
      rotstepmin=0.00001*rotstepmax
      fadd=0.0
      faddb=0.0
!Initialise the number of internal coordinates. This is used to determine in en_change_epsr20
!whether the internal coordinates changed (intra-molecular atom move, internal rotational group
!or whole molecule rotation) or the molecular centre coordinates were changed
      nnatemp=0
      
      do ic=1,ncell
         if(nmolsincell(ic).gt.0) then
            do imref=1,nmolsincell(ic)
               im=molsincell(imref,ic)
               call system_clock(time1)
               call make_neighbour_list_from_mols(im)
               call system_clock(time2)
               sumtime1=sumtime1+(time2-time1)/1000
 
               if(nisum.ge.intra.and.stepmi.gt.0.0) then
c
c if step size is non-zero, then do internal molecule structure refinement. This is to 
c prevent the build up of round-off errors from allowing the molecule to 
c become too distorted.
c
c move particles and calculate new chisq
! do not move the first atom of each molecule if this is the only
! atom in the molecule!
                  if(nato(im).gt.1) then
                     iflag=0
                     nmovei(1)=nmovei(1)+1
                     call movei_epsr20(im,fdiffb)
                     ediff=(fdiff+fdiffb)*temp
                     if(ediff.le.0.0) then
                        etest=10.0
                     else if(ediff.lt.20.) then
                        etest=exp(-ediff)
                     else
                        etest=0.
                     end if
                     rtest=ran1()
                     if(rtest.le.etest) then
                        fadd=fadd+fdiff
                        faddb=faddb+fdiffb
                        call save_energies()
                     else
                        nmovei(2)=nmovei(2)+1
                        iaadd=naf(im)-1
                        do iref=1,nato(im)
                           ia=iref+iaadd
                           xa(ia)=xr(iref)
                           ya(ia)=yr(iref)
                           za(ia)=zr(iref)
                        end do
                     end if
                  endif
c now adjust the coordinates so that the centre of mass has not changed
                  call centreofmassi(im)
               endif

c now try rotating internal groups if no of rotational groups is greater than 0

!!This code needs updating. ig will refer to the wrong rotational groups in the current version

               if(ngrrot.gt.0.and.stepri.gt.0.) then
                  ninter(1)=0
                  ninter(2)=0
                  ig=0
                  if(ngrot(im).gt.0) then
!                     do igm=1,ngrot(im)
                        ninter(1)=ninter(1)+1
!                        ig=nrf(im)+igm
                        call interm_epsr20(im,fdiff,fdiffb)
c only allow rotation if it doesn't violate constraints
                        ediff=(fdiff+fdiffb)*temp
                        if(ediff.le.0.0) then
                           etest=10.0
                        else if(ediff.lt.20.) then
                           etest=exp(-ediff)
                        else
                           etest=0.
                        endif
                        if(ran1().le.etest) then
                           fadd=fadd+fdiff
                           faddb=faddb+fdiffb
c now adjust the coordinates of centre of mass to correct for any change
                           call centreofmass(im)
                           call save_energies()
                        else
                           ninter(2)=ninter(2)+1
                           do ia=1,nato(im)
                              nat=naf(im)+(ia-1)
                              xa(nat)=xr(ia)
                              ya(nat)=yr(ia)
                              za(nat)=zr(ia)
                           end do
                        endif
 !                    end do
                  endif
               endif

c now rotate each molecule whole
c Only rotate molecules that have more than one atom

               if(stepr.gt.0.0.and.nato(im).gt.1) then

!               write(6,*) 'wholem> ',im
                  nwhole(1,mtype(im))=nwhole(1,mtype(im))+1
                  call wholem_epsr20(im,fdiff,fdiffb)
                  if(abs(fdiff).gt.100.0) then
                     iflag=1
                  else
                     iflag=0
                  endif
                  ediff=(fdiff+fdiffb)*temp
                  if(ediff.le.0.0) then
                     etest=10.0
                  else if(ediff.lt.20.) then
                     etest=exp(-ediff)
                  else
                     etest=0.
                  endif
                  if(ran1().le.etest) then
                     fadd=fadd+fdiff
                     faddb=faddb+fdiffb
                     call save_energies()
                  else
                     nwhole(2,mtype(im))=nwhole(2,mtype(im))+1
                     do ia=1,nato(im)
                           nat=naf(im)+(ia-1)
                           xa(nat)=xr(ia)
                           ya(nat)=yr(ia)
                           za(nat)=zr(ia)
                     end do
                  endif
c
c now adjust the coordinates of centre of mass to correct for any change
c
                  call centreofmass(im)
   
               endif
c now translate each molecule 

               if(stepm.gt.0.0) then
!      write(6,*) 'transm> ',im,fadd
! Save the coordinates of these atoms
                  nmovem(1,mtype(im))=nmovem(1,mtype(im))+1
                  call transm_epsr20(im,fdiff,fdiffb)
                  ediff=(fdiff+fdiffb)*temp
                  if(ediff.le.0.0) then
                     etest=10.0
                  else if(ediff.lt.20.) then
                     etest=exp(-ediff)
                  else
                     etest=0.
                  endif
                  if(ran1().le.etest) then
                     fadd=fadd+fdiff
                     faddb=faddb+fdiffb
                     call save_energies()
                  else
                     nmovem(2,mtype(im))=nmovem(2,mtype(im))+1
                     xadd(im)=xs
                     yadd(im)=ys
                     zadd(im)=zs
                  endif
               endif
               !Save all the relevant pair energies
!               write(6,*) 'Save> ',im,fadd
               call save_erg_pairs(im)
!               write(6,*) 'Done> ',im,fadd
            end do
         end if
      end do

      do im=1,nmol
          call minimagec(xadd(im),yadd(im),zadd(im))
          if(tethered(im)) call minimagec(phix(im),phiy(im),phiz(im))
      end do

      if(nisum.ge.intra) then
         nisum=0
      else
         nisum=nisum+1
      end if

!Adjust step sizes to give required rejection rate

      if(ireset.ne.1) then

         if(nmovei(2).gt.0) then
            rat=rejrate*real(nmovei(1))/real(nmovei(2))
            if(rat.gt.1.1) rat=1.1
            if(rat.lt.0.9) rat=0.9
         else
            rat=1.1
         end if
         stepmi=stepmi*rat
         if(stepmi.gt.trstepmax) stepmi=trstepmax
         if(stepmi.lt.trstepmin) stepmi=trstepmin

         if(ninter(2).gt.0) then
            rat=rejrate*real(ninter(1))/real((ninter(2)))
            if(rat.gt.1.1) rat=1.1
            if(rat.lt.0.9) rat=0.9
         else
            rat=1.1
         endif
         stepri=stepri*rat
         if(stepri.gt.rotstepmax) stepri=rotstepmax
         if(stepri.lt.rotstepmin) stepri=rotstepmin

         do imoltype=1,nmoltype

            if(nwhole(2,imoltype).gt.0) then
               rat=rejrate*real(nwhole(1,imoltype))/real((nwhole(2,imoltype)))
               if(rat.gt.1.1) rat=1.1
               if(rat.lt.0.9) rat=0.9
            else
               rat=1.1
            end if
            wholerot(imoltype)=wholerot(imoltype)*rat
            if(wholerot(imoltype).gt.rotstepmax) wholerot(imoltype)=rotstepmax
            if(wholerot(imoltype).lt.rotstepmin) wholerot(imoltype)=rotstepmin

            if(nmovem(2,imoltype).gt.0) then
               rat=rejrate*real(nmovem(1,imoltype))/real((nmovem(2,imoltype)))
               if(rat.gt.1.1) rat=1.1
               if(rat.lt.0.9) rat=0.9
            else
               rat=1.1
            endif
            wholetrans(imoltype)=wholetrans(imoltype)*rat
            if(wholetrans(imoltype).gt.trstepmax) wholetrans(imoltype)=trstepmax
            if(wholetrans(imoltype).lt.trstepmin) wholetrans(imoltype)=trstepmin

         end do
      endif
 
      fadds=fadds+fadd/nmol
      faddbs=faddbs+faddb/nmol

      return
      end

      subroutine movei_epsr20(im,fdiffb)

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      
      integer im,iref,ia,iaadd
      real fdiffb
      
c store present coordinates
      iaadd=naf(im)-1
      do iref=1,nato(im)
         ia=iref+iaadd
         xr(iref)=xa(ia)
         yr(iref)=ya(ia)
         zr(iref)=za(ia)
      end do
c 
c move the atoms
c
      do ia=naf(im),nal(im)
         call movea_epsr20(ia)
      end do
      fdiffb=0
      call en_change_intra_epsr20(im,fdiffb)
      return
      end

      subroutine movea_epsr20(ia)

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      amass=amassa(itype(ia))
      stepsz=stepmi/amass
c move atoms by random amounts
      za(ia)=za(ia)+stepsz*(2.*ran1()-1)
      ya(ia)=ya(ia)+stepsz*(2.*ran1()-1)
      xa(ia)=xa(ia)+stepsz*(2.*ran1()-1)
      return
      end

      subroutine interm_epsr20(im,fdiff,fdiffb)

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'

c store present coordinates

      iaadd=naf(im)-1
      do iref=1,nato(im)
         ia=iref+iaadd
         xr(iref)=xa(ia)
         yr(iref)=ya(ia)
         zr(iref)=za(ia)
      end do
c
c perform internal rotations
c
      do ig=nrf(im),nrl(im)
         call rotati_epsr20(ig)
      end do
c
c calculate change in energy
c
      fdiff=0.0
      fdiffb=0.0
      call en_change_intra_epsr20(im,fdiffb)
      call en_change_epsr20(im,fdiff)
      return
      end

      subroutine rotati_epsr20(ig)

      include 'epsr_dimension.inc'
      include 'molecom.inc'

c pvector is the vector about which the rotation occurs

      real*4 pvector(3)

c xvector is the vector from centre of rotation to an individual atom
c xpvector is the component of this vector along the axis of rotation
c xrvector is the component of this vector perpendicular to the axis of rotation
c xrpvector is the vector perpendicular to both axis of rotation and xrvector

      real*4 xvector(3),xpvector(3),xrvector(3),xrpvector(3)

      pi=3.141592653
c define random new orientation
      phi=pi*stepri*(2*ran1()-1)
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
         do ia=1,nrot(ig)
            nat=labrot(ig,ia)
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

      subroutine wholem_epsr20(im,fdiff,fdiffb)

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
c
c rotate entire molecule
c
      call rotatm_epsr20(im,wholerot(mtype(im)))
c
c Calculate the change in energy
c
      fdiff=0.0
      fdiffb=0.0
      call en_change_epsr20(im,fdiff)
      return
      end

      subroutine rotatm_epsr20(im,stepsz)
c routine to rotate entire molecule

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      pi=3.141592653
c first save the coordinates of this molecule
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
      iref1=naf(im)
      iref2=nal(im)
      if(iref2.gt.iref1) then
c choose an axis to rotate about
      phi=pi*stepsz*(2*ran1()-1)
      cphi=cos(phi)
      sphi=sin(phi)
c define rotations about x,y or z-axes
      test=3.*ran1()
      if(test.lt.1.0) then
c rotate about z axis
            do ia=iref1,iref2
                  xasave=xa(ia)*cphi-ya(ia)*sphi
                  ya(ia)=xa(ia)*sphi+ya(ia)*cphi
                  xa(ia)=xasave
            end do
      else if(test.lt.2.) then
c rotate about y-axis
            do ia=iref1,iref2
                  zasave=za(ia)*cphi-xa(ia)*sphi
                  xa(ia)=za(ia)*sphi+xa(ia)*cphi
                  za(ia)=zasave
            end do
      else
c rotate about x-axis
            do ia=iref1,iref2
                  yasave=ya(ia)*cphi-za(ia)*sphi
                  za(ia)=ya(ia)*sphi+za(ia)*cphi
                  ya(ia)=yasave
            end do
      end if
      endif
      return
      end

      subroutine transm_epsr20(im,fdiff,fdiffb)

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
c
c move molecule
c
      call movem_epsr20(im,wholetrans(mtype(im)))
c
c Calculate the change in energy
c
      fdiff=0.0
      fdiffb=0.0
      call en_change_epsr20(im,fdiff)

      return
      end

      subroutine movem_epsr20(im,stepsz)

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

c store present coordinates
      xs=xadd(im)
      ys=yadd(im)
      zs=zadd(im)
c move atoms by random amounts
      zadd(im)=zadd(im)+stepsz*(2.*ran1()-1)
      yadd(im)=yadd(im)+stepsz*(2.*ran1()-1)
      xadd(im)=xadd(im)+stepsz*(2.*ran1()-1)
      return
      end

      subroutine en_change_intra_epsr20(im,fdiffb)

      use ergpairreference
      
      real fdiffb
      integer im

! Calculate the change in energy as a result of moving all the atoms in the molecule
!      write(6,*) im,molergpair(1,im),ergpair(1,im),ergold(1),ergnew(1)
      ergnew(1)=en_mole_intra(im)
      fdiffb=fdiffb+ergnew(1)-ergold(1)

      return
      end                


      subroutine en_change_epsr20(im,fdiff)
c
c calculates the change in energy as a result of moving all the atoms in a molecule
c
      use ergpairreference
      use timer

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      
      real :: xnew(near),ynew(near),znew(near)
      real :: x2temp,y2temp,z2temp
      real :: energynew
      integer :: icref,inn,it,jt,iref,iaadd
            
!      write(6,*) 'en_change_whole_epsr20> ',im
      stepcl=cellst/dmax
      call system_clock(time1)
c initialise the centre of mass coordinates and atom types which do not
c need to be updated for every atom - only when a molecule number changes
!Get new absolute coordinates of atoms on molecule im
      iaadd=naf(im)-1
      do iref=1,nato(im)
         ia=iref+iaadd
         xnew(iref)=xa(ia)+xadd(im)
         ynew(iref)=ya(ia)+yadd(im)
         znew(iref)=za(ia)+zadd(im)
      end do
!Step through list of neighbouring molecules
      do inn=1,nnbref
         jm=nbref(inn)
         !Don't bother if both molecules are tethered - this is dealt with separately
         if(jm.ne.im.and..not.(tethered(im).and.tethered(jm))) then
               energynew=0.0
               do ja=naf(jm),nal(jm)
!Get the absolute coordinates
                  x2temp=xa(ja)+xadd(jm)
                  y2temp=ya(ja)+yadd(jm)
                  z2temp=za(ja)+zadd(jm)
                  jt=itype(ja)
!Step through moved atoms on molecule im
                  do iref=1,nato(im)
                     ia=iref+iaadd
c Set the atom type numbers
                     it=itype(ia)
                     icref=npairref(it,jt)
                     energynew=energynew+en_atom_inter_test(
     *xnew(iref),ynew(iref),znew(iref),x2temp,y2temp,z2temp,stepcl,xstar,icref,npt2)
!                     energychange1=en_change_inter_test(
!     *xnew(nat),ynew(nat),znew(nat),xold(nat),yold(nat),zold(nat),
!     *x2temp,y2temp,z2temp,stepcl,xstar,icref,npt2)
                  end do
               end do
               ergnew(inn)=energynew
!               write(6,*) energynew,energyold,ergold(inn)
               fdiff=fdiff+energynew-ergold(inn)
         end if
      end do
c Add on the effect of tethered molecules if present
      if(tethered(im)) call en_change_tethered(im,fdiff)
      call system_clock(time2)
      sumtime2=sumtime2+(time2-time1)/1000

      return
      end                

      function en_atom_intra(ia,ja,x1,y1,z1,x2,y2,z2)
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
      integer ia,ja,it,jt,iz,icref,jref
      logical notqtype

      en_atom_intra=0.0
      sume=0
c Calculate the internal energy for this atom pair
      stepcl=cellst/dmax
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
      do while(jref.lt.nbond(ia).and.iflag.eq.0)
         jref=jref+1
         if(ja.eq.iabond(ia,jref)) iflag=jref
      end do
      if(iflag.gt.0) then

!Don't apply minimum image if atoms are connected by bonds.
         dsq=xdif*xdif+ydif*ydif+zdif*zdif
         dtest=sqrt(dsq)
c find the r value this corresponds to
         rat=dtest/stepcl
         iz=int(rat)+2
         dtdif=dtest-xstar(iz)
         wtms=amassr(icref)
         tb=tbond(ia,jref)
         wtm=wtms/tb
         ddif=dtest-tb
         sume=wtm*ddif*ddif
         if(refine_intra.and.iz.le.npt2.and.notqtype) sume=sume
     1+refine_intra_fac*(xavs(iz,icref)+xavs1(iz,icref)*dtdif)

      else 

! apply minimum image convention if atoms are not connected by bonds
         call minimagec(xdif,ydif,zdif)
         dsq=xdif*xdif+ydif*ydif+zdif*zdif
         dtest=sqrt(dsq)
c find the r value this corresponds to
         rat=dtest/stepcl
         iz=int(rat)+2
         if(iz.le.npt2) then
            dtdif=dtest-xstar(iz)
c if no bonds were found for this atom pair apply molecular repulsive and empirical potential 
c in the usual way.
c Q-type atoms are not included in this calculation since their positions are
c defined by reference to other atoms in the molecule.
            sume=xav(iz,icref)+xav1(iz,icref)*dtdif
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
      dtdif=dtest-xstar(iz)
      if(iz.le.npt2) then
c sum the reference and empirical potentials
         sume=xavtot(iz,icref)+xavtot1(iz,icref)*dtdif
      else
         sume=0.0
      endif
      en_atom_inter=sume
      return
      end

      function en_atom_inter_test(x1,y1,z1,x2,y2,z2,stepcl,xstar,icref,npt2)

      include 'epsr_dimension.inc'
      include 'potential_par.inc'
c
c calculates the intermolecular energy between atoms ia and ja on distinct
c molecules im and jm. It is assumed the atoms have been preselected to occur 
c on distinct molecules.
c
      real xstar(*)
      real en_atom_inter_test,x1,y1,z1,x2,y2,z2,xdif,ydif,zdif,dsq,dtest,rat,dtdif,stepcl
      integer iz,icref

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

      function en_mole_intra(im)
!Calculates the intra- or inter- molecular energy of a pair of molecules using
!the current coordinates.      
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      real en_mole_intra,en_atom_intra

      en_mole_intra=0.0
      ia=naf(im)-1
      do while(ia.lt.nal(im)-1)
         ia=ia+1
         ja=ia
         do while(ja.lt.nal(im))
            ja=ja+1
            en_mole_intra=en_mole_intra+en_atom_intra(ia,ja,xa(ia),ya(ia),za(ia),xa(ja),ya(ja),za(ja))
         end do
      end do

      return
      end

      function en_mole_pair(im,jm)
!Calculates the inter- molecular energy of a pair of molecules using
!the current coordinates.      
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'

      real en_mole_pair,en_atom_inter,x1,y1,z1,x2,y2,z2

      en_mole_pair=0.0
      !Don't calculate pair energies if both molecules are tethered
      if(.not.(tethered(im).and.tethered(jm))) then
         do ia=naf(im),nal(im)
            x1=xa(ia)+xadd(im)
            y1=ya(ia)+yadd(im)
            z1=za(ia)+zadd(im)
            do ja=naf(jm),nal(jm)
               x2=xa(ja)+xadd(jm)
               y2=ya(ja)+yadd(jm)
               z2=za(ja)+zadd(jm)
!               write(6,*) im,jm,ia,ja
               en_mole_pair=en_mole_pair+en_atom_inter(ia,ja,x1,y1,z1,x2,y2,z2)
            end do
         end do   
      end if

      return
      end


