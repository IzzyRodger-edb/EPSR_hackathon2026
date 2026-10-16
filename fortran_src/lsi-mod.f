      program standalone
      
c     Program to perform the EPSR setup command
      
c     2009/3/30 Attempt to make this a standalone program
      
c     Get the command from the command line arguments

      call get_command_from_args()

c     Run the command

      call run_lsi()
      
      stop
      end        


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine run_lsi()
c     
c     2019-01-04 Derived from chains_sp_criterion
c     This program calculates the distribution of local structure indices (LSI) for a set of bonded/non-bonded molecules.
c     The method is given by Shi, Russo, and Tanaka, J. Chem. Phys. 149, 224502 (2018); and Russo and Tanaka, 2014
c     NATURE COMMUNICATIONS DOI: 10.1038/ncomms4556
      
c     In this method a list of bonded and non-bonded near neighbours for each molecule is set up.
c     The LSI is defined as the difference in distance from the central molecule between
c     the nearest non-bonded and the furthest bonded molecule. 

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'lsi.inc'
      integer*4 ierrout
      character*20 text
      real*8 fc,fc1

      real, allocatable, dimension(:) :: avbondsave,navbondedsave,navnonbondedsave
      real, allocatable, dimension(:) :: avbondsave2,navbondedsave2,navnonbondedsave2
      pi=4.0*atan(1.0)
      piconv=pi/180.0
c     
c     setup the clusters symbols
c     
      call setup_symbols_lsi()
c     
c     now setup the input file
c     
      nin=11
      nout=6
      
      call setup_input_file(nin,nout,ierrout)
c     
c     if an error was encountered on setup, exit now
c     
      if(ierrout.ne.0) return
      
c     
c     read the input file
c     

c$$$  RH: Actually read the .dat parameters into symbols arrays
      write(6,*) "Ran setup_input_file ok."
      

      call read_input_file(nin,nout,ierrout)
      write(6,*) "Ran read_input_file ok."
c     
c     if an error was encountered on read, exit now
c     
      if(ierrout.ne.0) return
      

c$$$     RH:
c$$$     Using the variables from the input file (*CLUSTERS.dat), which
c$$$     a restored in the char array symbol_value. 
c$$$
c$$$     Values read: 
c$$$       nz -- the max cluster size
c$$$       nsumt -- number of accumulations
c$$$  and the values for each cluster distribuition

c     
c     set up the number of lsi bins
c     
      iref=2
      read(symbol_value(iref),*) nz
      if(nz.gt.mz) nz=mz
c     
c     set up min and max values of lsi
c     
      iref=iref+1
      read(symbol_value(iref),*) lsimin,lsimax
c     
c     set up current accumulator
c     
      iref=iref+1
      read(symbol_value(iref),*) nsumt
c     
c     check the number of distributions does not exceed the array size
c     


c$$$  RH: nfiles is declared in symbols.inc, and set in
c$$$  setup_symbols_clusters() in setup_symbols.f, and it is the number
c$$$  of distributions
      ndist=nfiles
      if(ndist.gt.mdist/3) ndist=mdist/3
      allocate(navbondedsave(ndist),navnonbondedsave(ndist),avbondsave(ndist))
      allocate(navbondedsave2(ndist),navnonbondedsave2(ndist),avbondsave2(ndist))
c
c Two sets of distributions are calculated: 
c 1) a histogram of bonded neighbours
c 2) a histogram of lsis
c
      ndist3=3*ndist
c
c setup the type numbers for the three atoms used to define a bond
c also the distance ranges 1 - 3 and 2 - 3 needed to define the bond
c
      do id=1,ndist
            iref=ireff(nsymboluptodata+1,id)
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
c Save the number of words read to test the other values
            nwordsave(id)=nwords
            do iw=1,nwords
                  atom1(id,iw)=sometext(ncf(iw):ncl(iw))
            end do
            iref=iref+1
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.ge.nwordsave(id)) then
               do iw=1,nwords
                  atom2(id,iw)=sometext(ncf(iw):ncl(iw))
               end do
            else
               write(6,200) id,symbol(iref)
               return
            endif
            iref=iref+1
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.ge.nwordsave(id)) then
               do iw=1,nwords
                  atom3(id,iw)=sometext(ncf(iw):ncl(iw))
               end do
            else
               write(6,200) id,symbol(iref)
               return
            endif
            iref=iref+1
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.ge.nwordsave(id)) then
               do iw=1,nwords
                  read(sometext(ncf(iw):ncl(iw)),*) dff1(id,iw)
               end do
            else
               write(6,200) id,symbol(iref)
               return
            endif
            iref=iref+1
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.ge.nwordsave(id)) then
               do iw=1,nwords
                  read(sometext(ncf(iw):ncl(iw)),*) dff2(id,iw)
               end do
            else
               write(6,200) id,symbol(iref)
               return
            endif
            iref=iref+1
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.ge.nwordsave(id)) then
               do iw=1,nwords
                  read(sometext(ncf(iw):ncl(iw)),*) dfs1(id,iw)
               end do
            else
               write(6,200) id,symbol(iref)
               return
            endif
            iref=iref+1
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.ge.nwordsave(id)) then
               do iw=1,nwords
                  read(sometext(ncf(iw):ncl(iw)),*) dfs2(id,iw)
               end do
            else
               write(6,200) id,symbol(iref)
               return
            endif
            iref=iref+1
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.ge.nwordsave(id)) then
               do iw=1,nwords
                  read(sometext(ncf(iw):ncl(iw)),*) dft1(id,iw)
               end do
            else
               write(6,200) id,symbol(iref)
               return
            endif
            iref=iref+1
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.ge.nwordsave(id)) then
               do iw=1,nwords
                  read(sometext(ncf(iw):ncl(iw)),*) dft2(id,iw)
               end do
            else
               write(6,200) id,symbol(iref)
               return
            endif
            !Read the minimum and maximum angles for defining a bond. Please note that since we
            !use the cosine rule, the minimum angle is actually the maximum in cosine, and the maximum angle
            !is the minimum in cosine. The ranges are therefore read in reverse order to that specified in the symbols table.
            iref=iref+1
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.ge.nwordsave(id)) then
               do iw=1,nwords
                  read(sometext(ncf(iw):ncl(iw)),*) ang213max(id,iw)
                  ang213max(id,iw)=cos(ang213max(id,iw)*piconv)
               end do
            else
               write(6,200) id,symbol(iref)
               return
            endif
            iref=iref+1
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.ge.nwordsave(id)) then
               do iw=1,nwords
                  read(sometext(ncf(iw):ncl(iw)),*) ang213min(id,iw)
                  ang213min(id,iw)=cos(ang213min(id,iw)*piconv)
               end do
            else
               write(6,200) id,symbol(iref)
               return
            endif
! Read current average number of bonded molecules per atom1
            iref=iref+1
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.ge.2) then
               read(symbol_value(iref),*) navbondedsave(id),navbondedsave2(id)
               navbondedsave2(id)=navbondedsave2(id)*navbondedsave2(id)
            else
               read(symbol_value(iref),*) navbondedsave(id)
               navbondedsave2(id)=0.0
            end if
! Read current average number of non-bonded molecules per atom1
            iref=iref+1
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.ge.2) then
               read(symbol_value(iref),*) navnonbondedsave(id),navnonbondedsave2(id)
               navnonbondedsave2(id)=navnonbondedsave2(id)*navnonbondedsave2(id)
            else
               read(symbol_value(iref),*) navnonbondedsave(id)
               navnonbondedsave2(id)=0.0
            end if
! Read current average bond length for this distribution
            iref=iref+1
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.ge.2) then
               read(symbol_value(iref),*) avbondsave(id),avbondsave2(id)
               avbondsave2(id)=avbondsave2(id)*avbondsave2(id)
            else
               read(symbol_value(iref),*) avbondsave(id)
               avbondsave2(id)=0.0
            end if
      end do
200   format('lsi> Not found enough items in distribution ',i2,', symbol: ',a11)
      write(6,*) nz,ndist
      zstep=(lsimax-lsimin)/float(nz-1)
      zhalf=0.5*zstep
      do j=1,nz
            z(j)=lsimin+real(j-1)*zstep
      end do

      call initgofz(nz,ndist3,gofz,mz,mdist)

c normalise atomic coordinates to box size

      call normalise_coordinates()
C
C calculate the average value and std dev for each quantity
c
      if(nsumt.gt.0) then
         fc=1.0/dble(nsumt)
      else
         fc=1.0
      end if
      fc1=1.0-fc

c set up the table of lsi indices

      call find_lsi(fc,fc1,navbondedsave,navnonbondedsave,avbondsave
     *,navbondedsave2,navnonbondedsave2,avbondsave2)

c axis title

      titler='lsi'
      idref3=0
      do id=1,ndist
         do iw=1,nwordsave(id)
            write(6,*) id,atom1(id,iw),' ',atom2(id,iw),' ',atom3(id,iw)
            write(6,*) dff1(id,iw),dff2(id,iw),dfs1(id,iw),dfs2(id,iw)
         end do
            idref1=idref3+1
            idref2=idref1+1
            idref3=idref2+1

c set the column title for this pair of atoms

            ilent1=index(atom1(id,1),' ')-1
            ilent2=index(atom2(id,1),' ')-1
            ilent3=index(atom3(id,1),' ')-1
            write(titledat(idref1),100) atom1(id,1)(1:ilent1),atom2(id,1)(1:ilent2)
     *,atom3(id,1)(1:ilent3)
100   format(4x,a,'-',a,'-',a)
            titledat(idref2)=titledat(idref1)
            titledat(idref3)=titledat(idref1)
      end do

      write(6,*) (titledat(id),' ',id=1,idref3)

c output the data

      ext='n'

      call outgofz(fnameinp,ext,nz,ndist3,nsumt,z,gofz,gzt,gzt2
     *,titler,titledat,mz,mdist)
      do id=1,ndist
! Calculate average bond length, average number of bonded molecules, an average number of non-bonded molecules
      end do
! Calculate uncertainties on these numbers
      do id=1,ndist
         navbondedsave2(id)=sqrt(navbondedsave2(id))
         navnonbondedsave2(id)=sqrt(navnonbondedsave2(id))
         avbondsave2(id)=sqrt(avbondsave2(id))
      end do           
c
c write input file
c
c
c update the current value of nsumt
c
      icount=0
      text='nsize'
      call search_for_symbol(text,icount)
      iref=ireff(icount,1)
      write(symbol_value(iref),*) nz
      iref=iref+2
      write(symbol_value(iref),*) nsumt
      iref=iref+1
      write(symbol_value(iref),*) ndist
      do id=1,ndist
            iref=ireff(nsymboluptodata+12,id)
            write(symbol_value(iref),'(2(1x,e12.5))') navbondedsave(id),navbondedsave2(id)
            iref=iref+1
            write(symbol_value(iref),'(2(1x,e12.5))') navnonbondedsave(id),navnonbondedsave2(id)
            iref=iref+1
            write(symbol_value(iref),'(2(1x,e12.5))') avbondsave(id),avbondsave2(id)
      end do

      call save_files_others(14)

      return
      end

      subroutine find_lsi(fc,fc1,navbondedsave,navnonbondedsave,avbondsave
     *,navbondedsave2,navnonbondedsave2,avbondsave2)

c     This program calculates the distribution of local structure indices (LSI) for a set of bonded/non-bonded molecules.
c     The method is given by Shi, Russo, and Tanaka, J. Chem. Phys. 149, 224502 (2018); and Russo and Tanaka, 2014
c     NATURE COMMUNICATIONS DOI: 10.1038/ncomms4556

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lsi.inc'

      logical match
      logical foundia1,foundib1,foundia2,foundd12,foundd13,foundd23,foundiw
      logical found,foundneighb

! If molecule is already used to define a value of zeta it is not used again for another value

      logical alreadyused(mmol)
   
      integer icount,icount1,icountmax,ia1save,ia2save,iwsave,nbonded,nbondedsum
      integer lenbondsum,lennonbondsum

      real mindist,anglefound,dnonbondedmin,dbondedav,dbondedsum
      real lsi(mmol),avbondmol(mmol),nneighbbond(mmol),nneighbnonbond(mmol)
      real navbonded,navnonbonded,avbond
      real navbondedsave(*),navnonbondedsave(*),avbondsave(*)
      real navbondedsave2(*),navnonbondedsave2(*),avbondsave2(*)
      real*8 fc,fc1

      alreadyused=.false. !Define the initial values of alreadyused. This stops the same molecule being used to define zeta for more than one central molecule.
      zstep=(lsimax-lsimin)/float(nz-1)
      pi=4.0*atan(1.0)
      piconv=pi/180.0
      dmaxsqr=dmax*dmax
      idref3=0
      do id=1,ndist
      idref1=idref3+1
      idref2=idref1+1
      idref3=idref2+1
      do iz=1,nz
         gofz(iz,idref1)=0
         gofz(iz,idref2)=0
         gofz(iz,idref3)=0
      end do
c
c First set up for each molecule a list of bonded neighbours.
c
      do i=1,nmol
            nneighb(i)=0
            dneighb(i,1)=0
            nneighbbond(i)=0
      end do
c step through all origin moleculesc

      nfsum=0
      do im1=1,nmol
         xc1=xadd(im1)
         yc1=yadd(im1)
         zc1=zadd(im1)
         do im2=1,nmol
! Need to find the minimum separation between atoms atom1 on molecule 1 and atom 3 on molecule 2
! The two atoms cannot be on the same molecule
            if(im1.ne.im2) then
               xc2=xadd(im2)
               yc2=yadd(im2)
               zc2=zadd(im2)
               mindist=dmax
               ia1save=0
               ia2save=0
               iwsave=0
               iw=0
!      write(6,*) im1,im2,id,nwordsave(id)
               do while(iw.lt.nwordsave(id))
                  iw=iw+1
                  ia1=naf(im1)-1
!      write(6,*) im1,im2,id,nwordsave(id),ia1,iw
                  do while(ia1.lt.nal(im1))
                     ia1=ia1+1
!      write(6,'(4(1x,i5),2(1x,a,1x,a))') im1,ia1,nwordsave(id),iw,type(itype(ia1)),atom1(id,iw)
                     foundia1=match(type(itype(ia1)),atom1(id,iw))
                     if(foundia1) then
                        xa1=xa(ia1)+xc1
                        ya1=ya(ia1)+yc1
                        za1=za(ia1)+zc1 
                        ia2=naf(im2)-1
                        do while(ia2.lt.nal(im2))
                           ia2=ia2+1
                           foundia2=match(type(itype(ia2)),atom3(id,iw))
                           if(foundia2) then
c Get the coordinates of the third atom
                              xa2=xa(ia2)+xc2
                              ya2=ya(ia2)+yc2
                              za2=za(ia2)+zc2 
c Calculate the distance between ia1 and ia2
                              xdif=xa2-xa1
                              ydif=ya2-ya1
                              zdif=za2-za1
c apply periodic boundary conditions
                              call minimagec(xdif,ydif,zdif)
                              dsq13=(xdif*xdif+ydif*ydif+zdif*zdif)*dmaxsqr
                              da13=sqrt(dsq13)
                              if(da13.lt.mindist) then
                                 mindist=da13
                                 ia1save=ia1
                                 ia2save=ia2
                                 iwsave=iw
                              end if
                           end if
                        end do
                     end if
                  end do
               end do
!      write(6,*) im1,im2,ia1save,ia2save,iwsave,mindist
                        
! Step through existing distances for im1 and if this distance is less than any existing distances, or the total number 
! of distances is less than mbond, insert the relevant value in order
               icount=0
               foundneighb=.false.
               do while (icount.lt.nneighb(im1).and..not.foundneighb)
                  icount=icount+1
                  if(mindist.lt.dneighb(im1,icount)) then
! Move all the existing values up by one, up to the limit mbond, but only if im2 is different from im1 - otherwise simply replace the current value
                     icountmax=min(mbond,nneighb(im1)+1)
                     icount2=icountmax
                     do while (icount2.gt.icount)
                        icount1=icount2-1
                        neighb(im1,icount2)=neighb(im1,icount1)
                        dneighb(im1,icount2)=dneighb(im1,icount1)
                        ia1min(im1,icount2)=ia1min(im1,icount1)
                        ia2min(im1,icount2)=ia2min(im1,icount1)
                        iwmin(im1,icount2)=iwmin(im1,icount1)
                        icount2=icount1
                     end do
! Finally insert the new value
                     nneighb(im1)=icountmax
                     neighb(im1,icount)=im2
                     dneighb(im1,icount)=mindist
                     ia1min(im1,icount)=ia1save
                     ia2min(im1,icount)=ia2save
                     iwmin(im1,icount)=iwsave
                     foundneighb=.true.
!                     icount=icountmax !Make sure it doesn't add
                  end if
               end do
! If no values were inserted, and number of bonds is less that maximum, then insert this value at the end.
               if(.not.foundneighb.and.nneighb(im1).lt.mbond) then
                  nneighb(im1)=nneighb(im1)+1
                  neighb(im1,nneighb(im1))=im2
                  dneighb(im1,nneighb(im1))=mindist
                  ia1min(im1,nneighb(im1))=ia1save
                  ia2min(im1,nneighb(im1))=ia2save
                  iwmin(im1,nneighb(im1))=iwsave
               end if
!      write(6,*) im1,icount,nneighb(im1),dneighb(im1,1)
            end if
         end do
!         write(6,'(1x,i5,1x,i5)') im1,nneighb(im1)
!         write(6,'(5(1x,i5,1x,f10.4))') (neighb(im1,j),dneighb(im1,j),j=1,5)
! Now find out which of the nearest neighbours of im1 are bonded to it
         do icount=1,nneighb(im1)
            iw=iwmin(im1,icount)
            da13=dneighb(im1,icount)
! First check that the distance of this pair is within the specified distance
            found=da13.ge.dfs1(id,iw).and.da13.le.dfs2(id,iw)
            foundia2=.false.
            if(found) then
               dsq13=da13*da13
               ia1=ia1min(im1,icount)
               im2=neighb(im1,icount)
               ia3=ia2min(im1,icount)
               xc3=xadd(im2)
               yc3=yadd(im2)
               zc3=zadd(im2)
               xa3=xc3+xa(ia3)
               ya3=yc3+ya(ia3)
               za3=zc3+za(ia3)
               ia2=naf(im1)-1
               do while (ia2.lt.nal(im1).and..not.foundia2)
                  ia2=ia2+1
! Make sure atom2 is the right type
                  foundia2=match(type(itype(ia2)),atom2(id,iw))
                  if(foundia2) then
!Check the separation of ia1 and ia2
                     xdif=xa(ia2)-xa(ia1)
                     ydif=ya(ia2)-ya(ia1)
                     zdif=za(ia2)-za(ia1)
                     dsq12=(xdif*xdif+ydif*ydif+zdif*zdif)*dmaxsqr
                     da12=sqrt(dsq12)
                     foundia2=da12.ge.dff1(id,iw).and.da12.le.dff2(id,iw)
                     if(foundia2) then
!Check the separation of ia2 and ia3
                        xa2=xc1+xa(ia2)
                        ya2=yc1+ya(ia2)
                        za2=zc1+za(ia2)
                        xdif=xa3-xa2
                        ydif=ya3-ya2
                        zdif=za3-za2
                        call minimagec(xdif,ydif,zdif)
                        dsq23=(xdif*xdif+ydif*ydif+zdif*zdif)*dmaxsqr
                        da23=sqrt(dsq23)
                        foundia2=da23.ge.dft1(id,iw).and.da23.le.dft2(id,iw)
!      write(6,*) im1,im2,ia1,ia2,ia3,da12,da13,da23,foundia2
                        if(foundia2) then
!Finally need to check the angle between da12 and da13
                           anglefound=(dsq12+dsq13-dsq23)/(2.0*da12*da13)
                           foundia2=anglefound.ge.ang213min(id,iw).and.anglefound.le.ang213max(id,iw)
                        end if
                     end if
                  end if
               end do
            end if
            bneighb(im1,icount)=found.and.foundia2
         end do
      end do
! For each non-bonded neighbour of a molecule, we need to check that that molecule is not bonded to the first molecule
      do im1=1,nmol
         do icount1=1,nneighb(im1)
            if(.not.bneighb(im1,icount1)) then
               im2=neighb(im1,icount1)
! Is im2 bonded to im1?
               icount2=0
               found=.false.
               do while (icount2.lt.nneighb(im2).and..not.found)
                  icount2=icount2+1
                  found=neighb(im2,icount2).eq.im1.and.bneighb(im2,icount2)
               end do
               bneighb(im1,icount1)=found
            end if
         end do                        
! Calculate the lsi for this molecule
         dnonbondedmin=dmax
         dbondedav=0.0
         nbonded=0
         im2save=0
         do icount1=1,nneighb(im1)
            if(bneighb(im1,icount1)) then
               nbonded=nbonded+1
               dbondedav=(real(nbonded-1)*dbondedav+dneighb(im1,icount1))/real(nbonded)
            else
               im2=neighb(im1,icount1)
               if(dneighb(im1,icount1).lt.dnonbondedmin.and..not.alreadyused(im2)) then
                  dnonbondedmin=dneighb(im1,icount1)
                  im2save=im2
               end if
            end if
         end do
         if(im2save.gt.0.and.dbondedav.gt.0.0) then
            lsi(im1)=dnonbondedmin-dbondedav
            alreadyused(im2save)=.true.
         else
            lsi(im1)=dmax !Ignore values which are undefined.
         endif
         write(6,'(1x,i5,1x,i5,1x,f10.4)') im1,nneighb(im1),lsi(im1)
         write(6,'(5(1x,i5,1x,f10.4,1x,l))') (neighb(im1,j),dneighb(im1,j),bneighb(im1,j),j=1,nneighb(im1))
      end do

!      do im1=1,nmol
!         write(6,'(1x,i5,1x,i5)') im1,nneighb(im1)
!         write(6,'(10(1x,i5,1x,f10.4))') (neighb(im1,j),dneighb(im1,j),j=1,nneighb(im1))
!      end do

c Calculate the number of bonds, non-bonds, and average bond length for each molecule

      do im1=1,nmol
         nneighbbond(im1)=0
         nneighbnonbond(im1)=0
         avbondmol(im1)=0.0
         if(nneighb(im1).ge.0) then
            do icount=1,nneighb(im1)
               if(bneighb(im1,icount)) then
                  nneighbbond(im1)=nneighbbond(im1)+1.0
                  avbondmol(im1)=avbondmol(im1)+dneighb(im1,icount)
               else
                  iw=iwmin(im1,icount)
                  da13=dneighb(im1,icount)
                  if(da13.ge.dfs1(id,iw).and.da13.le.dfs2(id,iw)) then
                     nneighbnonbond(im1)=nneighbnonbond(im1)+1.0
                  end if
               end if
            end do
         endif
         if(nneighbbond(im1).gt.0) avbondmol(im1)=avbondmol(im1)/nneighbbond(im1)
      end do
! Form averages
      navbonded=0
      navnonbonded=0
      avbond=0
      lensum=0
      do im1=1,nmol
         navbonded=navbonded+nneighbbond(im1)
         navnonbonded=navnonbonded+nneighbnonbond(im1)
         if(avbondmol(im1).gt.0.0) then
            avbond=avbond+avbondmol(im1)
            lensum=lensum+1
         end if
      end do
      navbonded=navbonded/real(nmol)
      navnonbonded=navnonbonded/real(nmol)
      if(lensum.gt.0) avbond=avbond/real(lensum)         
      write(6,994) navbonded,navnonbonded
      write(6,'("Average distance of bonded molecules: ",e12.5)') avbond
994   format(1x,' Average of ',e12.5,' bonds/mol and ',e12.5,' non bonds/mol (in first shell)')
      navbondedsave(id)=real(fc*dble(navbonded)+fc1*dble(navbondedsave(id)))
      navnonbondedsave(id)=real(fc*dble(navnonbonded)+fc1*dble(navnonbondedsave(id)))
      avbondsave(id)=real(fc*dble(avbond)+fc1*dble(avbondsave(id)))
! Form average deviations
      navbonded=0
      navnonbonded=0
      avbond=0
      lensum=0
      do im1=1,nmol
         diff=nneighbbond(im1)-navbondedsave(id)
         navbonded=navbonded+diff*diff
         diff=nneighbnonbond(im1)-navnonbondedsave(id)
         navnonbonded=navnonbonded+diff*diff
         if(avbondmol(im1).gt.0.0) then
            diff=avbondmol(im1)-avbondsave(id)
            avbond=avbond+diff*diff
            lensum=lensum+1
         end if
      end do
      navbonded=navbonded/real(nmol)
      navnonbonded=navnonbonded/real(nmol)
      if(lensum.gt.0) avbond=avbond/real(lensum)         
      navbondedsave2(id)=real(fc*dble(navbonded)+fc1*dble(navbondedsave2(id)))
      navnonbondedsave2(id)=real(fc*dble(navnonbonded)+fc1*dble(navnonbondedsave2(id)))
      avbondsave2(id)=real(fc*dble(avbond)+fc1*dble(avbondsave2(id)))

! Histogram number of bonds and lsi for each molecule
      do im1=1,nmol
         izref=1+nint(nneighbbond(im1))
         gofz(izref,idref1)=gofz(izref,idref1)+1.0
         izref=1+nint(nneighbnonbond(im1))
         gofz(izref,idref2)=gofz(izref,idref2)+1.0
         izref=nint((lsi(im1)-lsimin+zhalf)/zstep)+1
         if(izref.gt.0.and.izref.lt.nz+1) then
            gofz(izref,idref3)=gofz(izref,idref3)+1.0
         end if
      end do

      end do !end of id loop

      return
      end




