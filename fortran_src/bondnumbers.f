      subroutine bondnumbers

      use molecom_mod
c
c to determine the correct atom numbers to assign to bonds in the input .ato
c file. This is to correct for the possibility that the atom numbers used to
c describe the intramolecular bonds may have got out sync. with the atom 
c numbering in this file
c
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
c
c The method is to step through all the atom numbers used to define bonds.
c The lowest of these is assumed to refer to the first atom in the molecule.
c The bond atom numbers are then revised to be in the correct relationship
c with the actual atom number of the first atom in the molecule
c
      ig=0
      do im=1,nmol
            namin=0
            do ia=naf(im),nal(im)
                  if(nbond(ia).gt.0) then
                        do ib=1,nbond(ia)
                            nt=iabond(ia,ib)
                            if(namin.eq.0) namin=nt
                            if(nt.lt.namin) namin=nt
                        end do
                  endif
            end do
c     write(6,*) im,itype(naf(im)),naf(im),namin
c
c reset the bond atom numbers if necessary
c
            if(namin.ne.0.and.namin.ne.naf(im)) then
                do ia=naf(im),nal(im)
                     if(imageatom(ia).gt.0) imageatom(ia)=imageatom(ia)-namin+naf(im)
                     if(nbond(ia).gt.0) then
                        do ib=1,nbond(ia)
                            nt=iabond(ia,ib)-namin
                            iabond(ia,ib)=naf(im)+nt
                        end do
                     endif
                end do
            endif
            if(ngrot(im).gt.0) then
                  do ih=1,ngrot(im)
                     ig=ig+1
c
c reset the internal rotation atom numbers if necessary
c
                     if(namin.ne.0.and.namin.ne.naf(im)) then
                        labaxi(ig,1)=labaxi(ig,1)-namin+naf(im)
                        labaxi(ig,2)=labaxi(ig,2)-namin+naf(im)
                        do j=rotf(ig),rotl(ig)
                           labrot(j)=labrot(j)-namin+naf(im)
                        end do
                        if(rotgrp(ig).eq.'DIH') then
                           do k=1,ndihedralang(ig)
                              do j=1,2
                                 labrotdih(j,k,ig)=labrotdih(j,k,ig)-namin+naf(im)
                              end do
                           end do
                        end if
                     endif
                  end do
            endif
      end do
      return
      end 
