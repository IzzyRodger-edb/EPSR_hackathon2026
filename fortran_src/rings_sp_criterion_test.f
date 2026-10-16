      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_rings()

      stop
      end      

      subroutine run_rings()
c
c Program to estimate number of rings in a distribution
c of molecules according to the 'shortest path' method. 
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'chains.inc'
      integer*4 ierrout
      character*20 text
c
c setup the rings symbols
c
      call setup_symbols_rings()
c
c now setup the input file
c
      nin=11
      nout=6
      call setup_input_file(nin,nout,ierrout)
c
c if an error was encountered on setup, exit now
c
      if(ierrout.ne.0) return
c
c read the input file
c
      call read_input_file(nin,nout,ierrout)
c
c if an error was encountered on read, exit now
c
      if(ierrout.ne.0) return
c
c set up the number of ring length number bins
c
      iref=2
      read(symbol_value(iref),*) nz
      if(nz.gt.mring) nz=mring
c
c set up current accumulator
c
      iref=iref+1
      read(symbol_value(iref),*) nsumt
c
c check that the number of distributions is not too big for array sizes
c
      nfilesmax=mdist/2
      if(nfiles.gt.nfilesmax) nfiles=nfilesmax
      ndist=nfiles
c
c three sets of distributions are calculated: 
c 1) a histogram of bonded neighbours
c 2) a histogram of ring lengths
c
      ndist3=2*ndist
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
      end do
200   format('rings> Not found enough items in distribution ',i2,', symbol: ',a11)
      write(6,*) nz,ndist
      zstep=0.5/float(nz)
      zhalf=0.5*zstep
      do j=1,nz
            z(j)=real(j-1)
      end do

      call initgofz(nz,ndist3,gofz,mz,mdist)

c normalise atomic coordinates to box size

      call normalise_coordinates()

c set up the table of molecular bonds

      call find_rings()

c axis title

      titler='Ring_size'
      idref2=0
      do id=1,ndist
         do iw=1,nwordsave(id)
            write(6,*) id,atom1(id,iw),' ',atom2(id,iw),' ',atom3(id,iw)
            write(6,*) dff1(id,iw),dff2(id,iw),dfs1(id,iw),dfs2(id,iw)
         end do
            idref1=idref2+1
            idref2=idref1+1

c set the column title for this pair of atoms

            ilent1=index(atom1(id,1),' ')-1
            ilent2=index(atom2(id,1),' ')-1
            ilent3=index(atom3(id,1),' ')-1
            write(titledat(idref1),100) atom1(id,1)(1:ilent1),atom2(id,1)(1:ilent2)
     *,atom3(id,1)(1:ilent3)
100   format(4x,a,'-',a,'-',a)
            titledat(idref2)=titledat(idref1)
      end do

      write(6,*) (titledat(id),' ',id=1,idref2)

c output the data

      ext='n'

      call outgofz(fnameinp,ext,nz,ndist3,nsumt,z,gofz,gzt,gzt2
     *,titler,titledat,mz,mdist)

C
C calculate the average value and std dev for each quantity
c
      do id=1,ndist3

c Calculate normalising constant

         anorm=0.0
         do i=1,nz
            anorm=anorm+gzt(i,id)
         end do

         if(anorm.gt.0) then
            sum1=0.0
            sum2=0.0
            do i=1,nz
                  term=z(i)*gzt(i,id)/anorm     
                  sum1=sum1+term
                  sum2=sum2+z(i)*term
            end do
            chainlen(id,1)=sum1
            std=sum2-sum1*sum1
            if(std.gt.0.0) then
                  chainlen(id,2)=sqrt(std)
            else
                  chainlen(id,2)=0.0
            endif

         else

            chainlen(id,1)=0.0
            chainlen(id,2)=0.0

         endif
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
      iref=iref+1
      write(symbol_value(iref),*) nsumt
      iref=iref+1
      write(symbol_value(iref),*) ndist

      call save_files_others(14)

c Output mean values

      do id=1,ndist
            idfirst=2*id-1
            idlast=idfirst+1
            write(14,102) id,(chainlen(idr,1),chainlen(idr,2),idr=idfirst,idlast)
      end do
102   format(1x,i3,6(1x,e12.5))

      close(14)

      return
      end

      subroutine find_rings()

C Method used is based on the "Shortest Path" criterion of D S Franzblau, Phys Rev B, vol
c 44, p4925-4930 (1991). This seems a very reasonable criterion, applicable in a variety
c of situations.

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'chains.inc'

      logical match
      logical foundia1,foundib1,foundia2,foundd13,foundiw
      logical found

c nring is the total number of rings found
c lenring(*) gives the number of molecules in each ring

      integer*4 nring,lenring(mmol),id,iw

      zstep=0.5/float(nz)
      pi=3.141592653
      idref2=0

      do id=1,ndist

      idref1=idref2+1
      idref2=idref1+1
      nring=0
c
c First set up for each molecule a list of bonded neighbours.
c
      do i=1,nmol
            nneighb(i)=0
      end do

c step through all origin moleculesc

      nfsum=0
      nssum=0

      do im1=1,nmol
        ia1=naf(im1)
        ib1=naf(im1)
        do im2=1,nmol
          ia2=naf(im2)
c Need to check that these two molecules are not the same and that they 
c are not already included in each others lists

          found=im1.eq.im2
          imt=0
          do while(imt.lt.nneighb(im1).and..not.found)
            imt=imt+1
            found=im2.eq.neighb(im1,imt)
          end do

c Also need to check there is room for more neighbours in the arrays

          found=found.or.nneighb(im1).eq.mbond.or.nneighb(im2).eq.mbond

c Start the loop through sets of atoms

          foundiw=.false..or.found
          iw=0
          do while(.not.foundiw)
                     
            iw=iw+1
            foundia1=.false..or.found
            ia1=naf(im1)-1
            do while(.not.foundia1)

              ia1=ia1+1
              foundia1=match(type(itype(ia1)),atom1(id,iw))

              foundib1=.false..or.found
              ib1=naf(im1)-1
              do while(foundia1.and..not.foundib1)

                ib1=ib1+1
                foundib1=match(type(itype(ib1)),atom2(id,iw))
                if(foundib1) nfsum=nfsum+1

c Get the coordinates of the first two atoms

                xc1=xadd(im1)
                yc1=yadd(im1)
                zc1=zadd(im1)
                xa1=xa(ia1)+xc1
                ya1=ya(ia1)+yc1
                za1=za(ia1)+zc1
                xb1=xa(ib1)+xc1
                yb1=ya(ib1)+yc1
                zb1=za(ib1)+zc1

                foundia2=.false..or.found
                ia2=naf(im2)-1
                do while(foundib1.and..not.foundia2)

                  ia2=ia2+1
                  foundia2=match(type(itype(ia2)),atom3(id,iw))
                  if(foundia2) then
                 
c Get the coordinates of the third atom

                    xc2=xadd(im2)
                    yc2=yadd(im2)
                    zc2=zadd(im2)
                    xa2=xa(ia2)+xc2
                    ya2=ya(ia2)+yc2
                    za2=za(ia2)+zc2 

c Calculate the distance between ia1 and ia2

                    xdif=xa2-xa1
                    ydif=ya2-ya1
                    zdif=za2-za1
c apply periodic boundary conditions
                    call minimagec(xdif,ydif,zdif)
                    dsq12=xdif*xdif+ydif*ydif+zdif*zdif
                    da12=sqrt(dsq12)*dmax
c check that distance is within bounds prescribed for this distribution
                    foundd13=da12.ge.dff1(id,iw).and.da12.le.dff2(id,iw)

                    if(foundd13) then

c Calculate the distance between ib1 and ia2

                      xdif=xa2-xb1
                      ydif=ya2-yb1
                      zdif=za2-zb1
c apply periodic boundary conditions
                      call minimagec(xdif,ydif,zdif)
                      dsq12=xdif*xdif+ydif*ydif+zdif*zdif
                      da12=sqrt(dsq12)*dmax
                      found=da12.ge.dfs1(id,iw).and.da12.le.dfs2(id,iw)

                      if(found) then

c O.K. found a bond - save molecule numbers

                        nneighb(im1)=nneighb(im1)+1
                        nneighb(im2)=nneighb(im2)+1
                        neighb(im1,nneighb(im1))=im2
                        neighb(im2,nneighb(im2))=im1

                      endif
                    endif
                  endif
                  foundia2=ia2.ge.nal(im2).or.found
                enddo              
                foundib1=ib1.ge.nal(im1).or.found
              end do

c Ensure the loop will continue to the end in case no bonds have been found so far.

              foundia1=ia1.ge.nal(im1).or.found
            end do
            foundiw=iw.ge.nwordsave(id).or.found
          end do
          found=found.and.iw.gt.0
          if(found) write(6,*) id,iw,im1,im2,ia1,' ',type(itype(ia1)),' ',ib1,' ',type(itype(ib1))
     *,' ',ia2,' ',type(itype(ia2))
        end do
      end do

c     do im1=1,nmol
c       write(6,999) im1,nneighb(im1)
c     1,(neighb(im1,j),j=1,nneighb(im1))
c     end do
c999  format(10(1x,i5))

c Calculate the average number of bonds per molecule

      lensum=0
      nsum=0
      do im1=1,nmol
            if(nneighb(im1).ge.0) then
                  lensum=lensum+nneighb(im1)
                  nsum=nsum+1
            endif
      end do
      if(nsum.gt.0) then
            avbond=float(lensum)/float(nsum)
      else
            avbond=0.0
      endif
      write(6,994) nsum,avbond
994   format(1x,i5,' molecules, with average of ',e12.5,' bonds/mol.')
      write(6,*) nfsum

c maximum distance to calculate shortest path

      nspmax=nz/2

c Starting from each molecule in turn look for rings of length not exceeding nz.
c Cannot form rings containing less than three molecules.

      nring=0
      im1=0

      do while(im1.lt.nmol.and.nring.lt.mmol)

         im1=im1+1

c No chance of finding a ring if this molecule has less than 2 neighbours

         if(nneighb(im1).gt.1) then

c Initially allow all the molecules .ge. im1 to be included in the search

            do im=im1,nmol
               carryon(im)=.true.
            end do

c Set up the shortest paths to this molecule

            call shortest_path(im1,nspmax,nmol,nneighb,neighb
     *,rneighb1,rneighb2,spvaluet,mmol,mbond)

            do im=1,nmol
               spvalue(im,1)=spvaluet(im)
            end do

c Now begin the search for rings proper.

c Place molecule im1 in the stack

            nrlab(0)=im1
            nrlab(1)=im1
            spseq(0)=-1
            spseq(1)=0
            rneighb(1)=0

c Now search the neighbours of molecule nrlab(icount) if there are any left

            icount=1

c ibase will count the base to be used in the 'antipodal' shortest path check

            ibase=1

            do while (icount.gt.0.and.nring.lt.mmol)

c If we have searched all the neighbours of the current stack pointer need to backtrack 
c and try another molecule

c For the origin molecule, we do not need to use the last neighbour as all possible rings
c with this molecule will already have been found.

               isub=0
               if(icount.eq.1) isub=1
               if(rneighb(icount).ge.(nneighb(nrlab(icount))-isub)) then

c If this molecule is a neighbour of im1, ensure it will not be used again, so as to 
c avoid double counting rings.

                  if(icount.eq.2) carryon(nrlab(icount))=.false.
                  icount=icount-1
                  if(ibase.gt.1) ibase=ibase-1


               else

c Increment the neighbour counter for this member of the ring

                  rneighb(icount)=rneighb(icount)+1

c Get the molecule number of the new neighbour

                  im2=neighb(nrlab(icount),rneighb(icount))

c im2 must be allowed, otherwise the rings associated with it have already been counted.

                  if(carryon(im2)

c Also the number of bonds associated with im2 must be greater than 1

     *.and.nneighb(im2).gt.1

c Also im2 cannot be equal to the one immediately preceeding icount in the ring

     *.and.im2.ne.nrlab(icount-1)

c Also don't waste time looking at rings longer than nz 
               
     *.and.icount.lt.nz

c Also don't waste time if the ring counter has reached its maximum value

     *.and.nring.lt.mmol) then

c Check this molecule does not already exist in the ring (apart from the first one
c of course, in which case the ring is closed!!)

                     ifound=0
                     ic=1
                     do while(ifound.eq.0.and.ic.lt.icount)

                        ic=ic+1
                        if(im2.eq.nrlab(ic)) ifound=im2

                     end do

                     if(ifound.eq.0) then

c Appear to have found a valid neighbour - use it.

c     write(6,*) im1,im2,icount,rneighb(icount),nneighb(nrlab(icount))
c     *,(neighb(nrlab(icount),ic),ic=1,nneighb(nrlab(icount)))

c If im2 .eq. im1, the ring is closed.

                        if(im2.eq.im1.and.nring.lt.mmol) then

                           nring=nring+1
                           lenring(nring)=icount

      write(6,*)
      write(6,*) im1,nring,lenring(nring),(nrlab(ic),ic=1,icount)
      write(6,*) im1,(spseq(ic),ic=1,icount)
                  
                        else

c before we can include this molecule we need to check that it will create a unimodal
c distribution of SP distances

                           sptest=spvalue(im2,1)-spseq(icount)
                           spgrad=spseq(icount)-spseq(icount-1)

c If the old gradient is zero or negative, the new molecule can only be accepted if
c the new gradient is negative. If the old gradient is positive, the new molecule can be
c accepted under all circumstances. Also maximum increment + or - 1.

                           if(abs(sptest).lt.2
     *.and.spgrad.gt.0.or.(sptest.lt.0.and.spgrad.le.0)) then

c Before we can accept the ring need to check that antipodes are shortest paths

c If both previous and new gradients are positive, then not possible to check 
c the 'antipodes' to be shortest paths so we accept the ring at this stage.
c But we do need to calculate the shortest paths to other molecules

                              if(spgrad.gt.0.and.sptest.gt.0) then

                                 icount=icount+1
                                 nrlab(icount)=im2
                                 rneighb(icount)=0
                                 spseq(icount)=spvalue(im2,1)

c Set up the shortest paths to this new ring molecule

                                 call shortest_path(im2,nspmax,nmol
     *,nneighb,neighb,rneighb1,rneighb2,spvaluet,mmol,mbond)

                                 do im=1,nmol
                                    spvalue(im,icount)=spvaluet(im)
                                 end do

c Cannot check the shortest path between antinodes

                                 sptrue1=.false.
                                 sptrue2=sptrue1

                              else if(spgrad.gt.0.and.sptest.lt.0) then

c Turning point for even numbered rings. Check the antipode is also a shortest path

                                 nodd=.false.
                                 sptrue1=.true.

                              else if(spgrad.gt.0.and.sptest.eq.0) then

c Turning point for odd numbered rings. Check both antipodes are shortest path

                                 nodd=.true.
                                 sptrue1=.true.

                              endif

c Check antipode(s) is/are shortest path

                              if(sptrue1) then

c Distance of current ring vertex to current base

                                 rptest=icount-ibase
                                 sptrue1=rptest.eq.spvalue(nrlab(icount),ibase)
                                 if(nodd) then
                                    sptrue2=rptest.eq.spvalue(im2,ibase)
                                 else
                                    sptrue2=sptrue1
                                 endif
                              endif
                              if(sptrue1.and.sptrue2) then        

                                 icount=icount+1
                                 ibase=ibase+1
                                 nrlab(icount)=im2
                                 rneighb(icount)=0
                                 spseq(icount)=spvalue(im2,1)
c     write(6,*)
c     write(6,*) im1,im2,icount,ibase,spseq(icount)

                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            end do
         endif

c As all possible rings associated with this molecule have been counted, it cannot
c be used again

         carryon(im1)=.false.
      end do

      if(nring.gt.0) then
            lensum=0
            do in=1,nring
              lensum=lensum+lenring(in)
            end do
            avlen=float(lensum)/float(nring)
      else
            avlen=0.0
      endif
      write(6,996) nring,avlen
996   format(1x,i5,' rings of average size = ',e12.5,' molecules')

c form histogram of number of bonded neighbours

      do im1=1,nmol
        iref=nneighb(im1)+1
        if(iref.le.nz) gofz(iref,idref1)=gofz(iref,idref1)+1
      end do

c form histogram of ring lengths

      do ir=1,nring
        iref=lenring(ir)+1
        if(iref.le.nz) gofz(iref,idref2)=gofz(iref,idref2)+1
      end do

c normalise the distributions to probabilities
c
c     sum1=0.0
c     sum2=0.0
c     do i=1,nz
c        sum1=sum1+gofz(i,idref1)
c        sum2=sum2+gofz(i,idref2)
c     end do
c     do i=1,nz
c        if(nfsum.gt.0) gofz(i,idref1)=gofz(i,idref1)/nfsum
c        if(nfsum.gt.0) gofz(i,idref2)=gofz(i,idref2)/nfsum
c     end do

      end do
      return
      end


