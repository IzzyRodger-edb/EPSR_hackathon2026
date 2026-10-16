      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call readmole()

      stop
      end      

      subroutine readmole()
c
c reads a ghemical .gpr, or Jmol .jmol, or Cerius2 .ms or .msi, file  and creates a .mol template file
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'readmole.inc'
      include 'chemical_symbols.inc'

      logical test

      character*20 ext,formatlabel
      character*3 int_as_str1, int_as_str2, int_as_str3
      character*1 iaccept

      integer sizecounttype,sizechsymbol
      integer counttype(melements),labeltype,j,k
     +,ifoundneighb(near),natomtypelist(mcomp)
c Arrays to store vectors
      real a(3),a1(3),a2(3),b(3),c1perp(3),c2perp(3),csign(3)

      external cosd,sind
c      write(6,*) nwords,sometext(ncf(1):ncl(nwords))
c
c generate the extension name of the file to read
c
      ilenext=0
      ilenreadfile=0
      iext=0
      if(nwords.gt.1) then
            readfile=sometext(ncf(2):ncl(2))
            ilenreadfile=ncl(2)-ncf(2)+1

c Generate the extension for this file

            iext=indexlastdot(readfile(1:ilenreadfile))
            if(iext.gt.0) then
               ext=readfile(iext:ilenreadfile)
               ilenext=ilenreadfile-iext+1
               ilenreadfile=ilenreadfile-ilenext
            else
               ilenext=0
            endif
c Otherwise it is possible a third argument has been supplied to give the filename
            if(nwords.gt.2) then
               readfile=sometext(ncf(3):ncl(nwords))
               ilenreadfile=ncl(nwords)-ncf(3)+1
c Knockoff the bit from the last dot onwards if it exists
               iext=indexlastdot(readfile)
               ilenreadfile=ilenreadfile-iext
            endif
c If the start of the extension is after the first character, treat the first part of readfile
c as the file name
            if(ilenreadfile.gt.0) readfile=readfile(1:ilenreadfile)
      endif
      write(6,*) ilenreadfile,readfile(1:ilenreadfile),ilenext,ext(1:ilenext)   

c If no extension is specified so far then ask for it

      if(ilenext.eq.0) then
            write(6,510,advance='no')
510   format(/'readmole> specify type of file to read [gpr or jmol]: ')
            read(5,'(a)',iostat=ierr) ext
            if(ierr.ne.0) then
               write(6,512) 
512   format(/'readmole> Invalid extension specified.')
               return
            endif
            ilenext=index(ext,' ')
      endif
      

c Check this extension includes the dot

      if(ext(1:1).ne.'.') then
            ext='.'//ext
            ilenext=ilenext+1
      endif

c Check that valid extension has been specified

      if(ext(1:ilenext).ne.'.gpr'
     *.and.ext(1:ilenext).ne.'.jmol'
     *.and.ext(1:ilenext).ne.'.ms'
     *.and.ext(1:ilenext).ne.'.msi'
     *) then
            write(6,511) ext(1:ilenext)
511   format(/'readmole> Invalid extension specified: ',a)
            return
      endif

c Finally grab a file if it has not already been defined.

      if(ilenreadfile.eq.0) then
         nwords=0
         do while (nwords.le.0)
              call search_for_files(ext,ilenext,iaccept)
              call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
              if(nwords.le.0) then
                 write(6,499)
 499    format(/'readmole> No file specified! - routine will exit'/)
                 return
              endif
         end do
         readfile=sometext(ncf(1):ncl(nwords))
         ilenreadfile=ncl(nwords)-ncf(1)+1
      endif
c
c strip out any extension that was supplied - extension must be .ato, otherwise
c nothing is done.
c
      ilenuptodot=indexlastdot(readfile)-1
      if(ilenuptodot.gt.0) then
         readfile=readfile(1:ilenuptodot)
         ilenreadfile=ilenuptodot
      endif

c
c first stage is to read the molecular template file and generate
c the atoms and bonds reference integers
c
      fname_ato=readfile(1:ilenreadfile)//ext(1:ilenext)
      ilenfname_ato=ilenreadfile+ilenext
      fnameato=directory_name(1:nchardir)
     *//fname_ato(1:ilenfname_ato)
      ilenfnameato=nchardir+ilenfname_ato

      write(6,601) fnameato(1:ilenfnameato)
601   format(/'readmole> Reading from file: ',a)

!Assume the unit cell dimensions are initially undefined

      unitcelldefined=.false.
      if(ext(1:ilenext).eq.'.gpr') then

         call readgpr()
         iatomlabeladd=1

      else if(ext(1:ilenext).eq.'.jmol') then

         call readjmol()
         iatomlabeladd=0

      else if(ext(1:ilenext).eq.'.ms'.or.ext(1:ilenext).eq.'.msi') then

         call readcerius2()
         iatomlabeladd=0

      else

         return

      endif

      if(ierr.ne.0) return

c Generate new atom labels to reduce the number of labels

c Three options are available:-

c 1) Simply use chemical symbol to set up the label.

c 2) Use the bonding of each atom to set up the label. If the number of neighbours of a given atom and the atom types
c of those neighbours are the same as a previous atom of the same type, that atom is given the same label as the previous
c atom.

c 3) Use the same names as in Jmol; then use changelabel statements to change
c labels to the same as in option number 2.

      ilenatomtypelist=0
      test=.true.

c Set the counttypes to zero 

      do ia1=1,natom

         counttype(iatomicno(ia1))=0

      end do

      write(6,400,advance='no')
400   format(/'readmole> Specify whether to use only chemical symbols to define atom labels [1],'/
     +,' or use bonding patterns as well [2],'/
     +,' or use bonding patterns and a changelabel section [3]: ')

      read(5,*) labeltype
c      labeltype=2

      do ia1=1,natom

         write(6,*) ia1,' ',atomlabel(ia1)(1:ilenatomlabel(ia1)),' ',nbond(ia1),(iabond(ia1,j),j=1,nbond(ia1))

         ib1=0
         test=.false.

         do while (ib1.lt.ilenatomtypelist.and..not.test)

              ib1=ib1+1
c Check to see whether the atomic number matches one of those in the list
              test=atomicnolist(ib1).eq.iatomicno(ia1)
              if((test.and.labeltype.eq.2).or.(test.and.labeltype.eq.3)) then

c Now check, for labeltype=2 or labeltype=3, that this atom has the same number of bonds as that in the list

                  test=nbond(ia1).eq.natomlist(ib1)      
                  if(test) then

c If so, check the list of neighbours of the current atom and see if they each have an atomic number common to an 
c atom in the list.

                     do j=1,natomlist(ib1)
                        ifoundneighb(j)=0
                     end do

c Step through the neighbours. If the atomic numbers are the same, signal this but also mark that this atomic number
c has already been found, so cannot be used again.

                     nfound=0
                     do j=1,nbond(ia1)
                        k=1
                        test=iatomicno(iabond(ia1,j)).eq.iabondlist(ib1,k)
     +.and.ifoundneighb(k).eq.0
                        do while(k.lt.natomlist(ib1).and..not.test)
                           k=k+1
                           test=iatomicno(iabond(ia1,j)).eq.iabondlist(ib1,k)
     +.and.ifoundneighb(k).eq.0
                        end do
                        if(test) then
                            nfound=nfound+1
                            ifoundneighb(k)=1
                        endif
                     end do
c If the number of found neighbours is the same as in the list, then this atom is the same type as the one
c in the list. From above a 'found' neighbour is one which has the same atomic number as one in the list.
                     test=nfound.eq.natomlist(ib1)
                  endif
              endif     
         end do

c If an existing entry for this atom was not found, create a new one.

         if(.not.test) then

            if(ilenatomtypelist.lt.mcomp) then

               ib1=ilenatomtypelist+1
c Save the number of bonded atoms
               natomlist(ib1)=nbond(ia1)
               natomtypelist(ib1)=1
c Save the atomic number of this atom - this will be used to identify the atoms
               atomicnolist(ib1)=iatomicno(ia1)
c Save the list of atomic numbers in the neighbours of this atom
               do j=1,nbond(ia1)
                  iabondlist(ib1,j)=iatomicno(iabond(ia1,j))
               end do
c               chargelist(ib1)=acharge(ia1)
c Increment a counter for this atomic number
               counttype(iatomicno(ia1))=counttype(iatomicno(ia1))+1
c      write(6,'(1x,i5,1x,i5,1x,a)') iatomicno(ia1),counttype(iatomicno(ia1)),chsymbol(iatomicno(ia1))
               sizechsymbol=index(chsymbol(iatomicno(ia1)),' ')-1
               if(counttype(iatomicno(ia1)).lt.10) then
                  sizecounttype=1
               else if(counttype(iatomicno(ia1)).lt.100) then
                  sizecounttype=2
               else 
                  sizecounttype=3
               endif

c Ensure we don't generate an error when writing the atom label

               do while((sizechsymbol+sizecounttype).gt.3)

                  counttype(iatomicno(ia1))=counttype(iatomicno(ia1))-1
                  if(counttype(iatomicno(ia1)).lt.10) then
                     sizecounttype=1
                  else if(counttype(iatomicno(ia1)).lt.100) then
                     sizecounttype=2
                  else 
                     sizecounttype=3
                  endif

               end do

c            write(6,*) sizechsymbol,sizecounttype,counttype(iatomicno(ia1))

               write(formatlabel,201) sizechsymbol,sizecounttype
201   format('(a',i1,',i',i1,')')
c               write(6,*) formatlabel
               
               write(atomlabellist(ib1),formatlabel,iostat=ierr)
     +chsymbol(iatomicno(ia1)),counttype(iatomicno(ia1))

c If an error occurred in writing the label, then abort the program

               if(ierr.ne.0) then

                  write(6,506) counttype(iatomicno(ia1))
     +,chsymbol(iatomicno(ia1)) 
 506  format('readmole> Error creating atom label'
     +,i5,a2,/' Probably too many atom types -'
     +,' Will abort program')
                  return

               endif
               ilenatomtypelist=ib1
               atomlabel(ia1)=atomlabellist(ib1)
               atomtype(ia1)=ib1

            else

               write(6,500) mcomp
 500  format('readmole> When setting revised labels, more than ',i3
     +,' components required - program aborted')
                  return

            endif

         else

            atomlabel(ia1)=atomlabellist(ib1)
            atomtype(ia1)=ib1
            natomtypelist(ib1)=natomtypelist(ib1)+1

         endif

      end do

c Write out what has been found so far

      write(6,521)
521   format(/'readmole> Following atom types and number of each type found:-'/)
      do ib1=1,ilenatomtypelist

         write(6,522) ib1,atomlabellist(ib1),natomtypelist(ib1)
522   format(1x,i5,1x,a3,1x,i5)

      end do

c For Mopac to work correctly we need to ensure that the neighbours of each atom occur consecutively in the list.

      call reorder_atomlist()

c Run MOPAC on this molecule?

      write(6,520,advance='no')
520   format(/'readmole> Do you want to run MOPAC on this molecule? [y or n]? ')
      read(5,'(a)') sometext
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(sometext(ncf(1):ncl(1)).eq.'y'.or.sometext(ncf(1):ncl(1)).eq.'Y') then

         call run_mopac()

      endif

c Get the charge values in case these have changed

      do ib1=1,ilenatomtypelist
         nchargelist(ib1)=0
         chargelist(ib1)=0.0
      end do
      do ia1=1,natom
         ib1=0
         test=.false.
         do while(ib1.lt.ilenatomtypelist.and..not.test)
            ib1=ib1+1
            test=atomtype(ia1).eq.ib1
         enddo
         if(test) then
            write(6,*) 'chargelist> ',ib1,ia1,acharge(ia1)
            chargelist(ib1)=chargelist(ib1)+acharge(ia1)
            nchargelist(ib1)=nchargelist(ib1)+1
         endif
      end do
      do ib1=1,ilenatomtypelist
         write(6,*) ib1,nchargelist(ib1),chargelist(ib1)
!Average charge per atom type
         if(nchargelist(ib1).gt.0) then
            chargelist(ib1)=chargelist(ib1)/real(nchargelist(ib1))
         end if
      end do

c initial atom number density

      rho=0.001

c Set up the name of the output file

      fnameato=fnameato(1:indexlastdot(fnameato)-1)//'.mol'
      ilenfnameato=indexlastdot(fnameato)+3

      open(10,file=fnameato,status='unknown',iostat=ierr)
      if(ierr.ne.0) then

         write(6,501) fnameato(1:ilenfnameato)
 501  format('readmole> Output file ',a,' could not be opened')
         return

      endif

c Signal that this is a gmol file and that the atom numbers listed 
c are 1 less than EPSR atom numbers

      write(10,*) '.gmol ',iatomlabeladd


c Step through the atoms and set up bond distances.

        ia1=0
        do while (ia1.lt.natom)

           ia1=ia1+1
           ib1=0
           do while (ib1.lt.nbond(ia1))

              ib1=ib1+1
              xdif=xa(iabond(ia1,ib1))-xa(ia1)
              ydif=ya(iabond(ia1,ib1))-ya(ia1)
              zdif=za(iabond(ia1,ib1))-za(ia1)
              dsq=xdif*xdif+ydif*ydif+zdif*zdif
              tbond(ia1,ib1)=sqrt(dsq)

           end do
        end do

c Write out the atom coordinates and bond pairs,

       if(labeltype.ne.3) then
        do ia1=1,natom
           write(10,*) 'atom ',iatomlabel(ia1),' ',atomlabel(ia1),' '
     *,xa(ia1),ya(ia1),za(ia1),nbond(ia1)
     *,(iatomlabel(iabond(ia1,j)),j=1,nbond(ia1))
        end do
       else if(labeltype.eq.3) then
        do ia1=1,natom
           write(int_as_str1,'(I3.3)') iatomlabel(ia1)
           write(10,*) 'atom ',iatomlabel(ia1),' '//int_as_str1//' '
     *,xa(ia1),ya(ia1),za(ia1),nbond(ia1)
     *,(iatomlabel(iabond(ia1,j)),j=1,nbond(ia1))
        end do
       endif


c Set up the bond list. 
c Where multiple bonds for a given atom pair are defined,
c the mean bond distance is calculated

        nbondlist=0

        do ia1=1,natom

           ib1=0
           do while (ib1.lt.nbond(ia1))
              
              ib1=ib1+1
            ia2=iabond(ia1,ib1)
              if(ia2.gt.ia1) then

                 if(atomtype(ia1).le.atomtype(ia2)) then
                if(labeltype.ne.3) then
                bondlabel=atomlabel(ia1)//' '//atomlabel(ia2)
                else if(labeltype.eq.3) then
                write(int_as_str1,'(I3.3)') iatomlabel(ia1)
                write(int_as_str2,'(I3.3)') iatomlabel(ia2)
                bondlabel=int_as_str1//' '//int_as_str2
                end if

             else

                if(labeltype.ne.3) then
                bondlabel=atomlabel(ia2)//' '//atomlabel(ia1)
                else if(labeltype.eq.3) then
                write(int_as_str1,'(I3.3)') iatomlabel(ia1)
                write(int_as_str2,'(I3.3)') iatomlabel(ia2)
                bondlabel=int_as_str2//' '//int_as_str1
                end if

             endif

c Search the bond list to see if this label has already been defined!{{{

                 ib2=0
                 test=.true.
                 do while (ib2.lt.nbondlist.and.test)

                    ib2=ib2+1
                    test=bondlabellist(ib2).ne.bondlabel

                 enddo

                 if(test.and.nbondlist.lt.mstarlabel) then

                    ib2=nbondlist+1
                    bondlabellist(ib2)=bondlabel
                    bondlist(ib2)=0.0
                    nbondlabel(ib2)=0
                    nbondlist=ib2
                   
                 endif

                 bondlist(ib2)=bondlist(ib2)+abs(tbond(ia1,ib1))
                 nbondlabel(ib2)=nbondlabel(ib2)+1

              endif
           end do
        end do

c Write out the bonds

        ib1=0
        do while (ib1.lt.nbondlist)

           ib1=ib1+1
           write(10,*) 'bond ',bondlabellist(ib1),' '
     +,bondlist(ib1)/nbondlabel(ib1)

        end do

c Setup the bondangle list

c Step through the atoms and set up the angles between bonded pairs

        pi=4.0*atan(1.0)
        piconv=180/pi
        nbondanglelist=0

        ia1=0
        do while (ia1.lt.natom)

           ia1=ia1+1
           ib1=0
           do while (ib1.lt.nbond(ia1)-1)

              ib1=ib1+1
              iaref1=iabond(ia1,ib1)
              a(1)=xa(iaref1)-xa(ia1)
              a(2)=ya(iaref1)-ya(ia1)
              a(3)=za(iaref1)-za(ia1)
              ib2=ib1
              do while(ib2.lt.nbond(ia1))

                 ib2=ib2+1
                 iaref2=iabond(ia1,ib2)
                 b(1)=xa(iaref2)-xa(ia1)
                 b(2)=ya(iaref2)-ya(ia1)
                 b(3)=za(iaref2)-za(ia1)

                 angle=piconv*acos(scapro(a,b)
     */sqrt(scapro(a,a)*scapro(b,b)))!}}}

c Setup the bondangle label

             if(atomtype(iaref1).le.atomtype(iaref2)) then

                if(labeltype.ne.3) then
                bondanglelabel=atomlabel(iaref1)
     +//' '//atomlabel(ia1)
     +//' '//atomlabel(iaref2)
                else if(labeltype.eq.3) then
                write(int_as_str1,'(I3.3)') iatomlabel(iaref1)
                write(int_as_str2,'(I3.3)') iatomlabel(ia1)
                write(int_as_str3,'(I3.3)') iatomlabel(iaref2)
                bondanglelabel=int_as_str1//' '//int_as_str2//' '//int_as_str3
                end if

             else

                if(labeltype.ne.3) then
                bondanglelabel=atomlabel(iaref2)
     +//' '//atomlabel(ia1)
     +//' '//atomlabel(iaref1)
                else if(labeltype.eq.3) then
                write(int_as_str1,'(I3.3)') iatomlabel(iaref2)
                write(int_as_str2,'(I3.3)') iatomlabel(ia1)
                write(int_as_str3,'(I3.3)') iatomlabel(iaref1)
                bondanglelabel=int_as_str1//' '//int_as_str2//' '//int_as_str3
                end if

             endif

c Step through the existing list of labels to see if any exist already

             ib3=0
             test=.true.
             do while(ib3.lt.nbondanglelist.and.test)

                    ib3=ib3+1
                test=bondanglelabellist(ib3).ne.bondanglelabel

             end do

                 if(test.and.nbondanglelist.lt.mstarlabel) then

                    ib3=nbondanglelist+1
                    bondanglelabellist(ib3)=bondanglelabel
                bondanglelist(ib3)=0.0
                    nbondanglelabel(ib3)=0
                    nbondanglelist=ib3
                   
                 endif

                 bondanglelist(ib3)=bondanglelist(ib3)+angle
                 nbondanglelabel(ib3)=nbondanglelabel(ib3)+1

              end do
           end do
        end do

c Write out list of angles found

      ib1=0
      do while (ib1.lt.nbondanglelist)

         ib1=ib1+1
         write(10,*) 'angle ',bondanglelabellist(ib1)
     +,bondanglelist(ib1)/nbondanglelabel(ib1)
      
      end do

c Search for all possible dihedral angles

        ndihedrals=0

        ia1=0
        do while (ia1.lt.natom)

           ia1=ia1+1
           ib1=0
           do while (ib1.lt.nbond(ia1))

              ib1=ib1+1
              ia2=iabond(ia1,ib1)

c Define 2 axis atoms

              if(ia2.gt.ia1) then

                 b(1)=xa(ia2)-xa(ia1)
                 b(2)=ya(ia2)-ya(ia1)
                 b(3)=za(ia2)-za(ia1)

c Get neighbours of ia1

                 ib2=0
                 do while (ib2.lt.nbond(ia1))

                    ib2=ib2+1
                    ia1ref=iabond(ia1,ib2)
                    if(ia1ref.ne.ia2) then

                       a1(1)=xa(ia1ref)-xa(ia1)
                       a1(2)=ya(ia1ref)-ya(ia1)
                       a1(3)=za(ia1ref)-za(ia1)

c Calculate the vector perpendicular to the plane formed by these two vectors

                       call vecpro(a1,b,c1perp)

c Normalise this vector

                       amp=sqrt(scapro(c1perp,c1perp))
                       do i=1,3
                          c1perp(i)=c1perp(i)/amp
                       end do

c Now search the bonds of ia2 for the 4th atom

                       ib3=0
                       do while (ib3.lt.nbond(ia2))

                          ib3=ib3+1
                          ia2ref=iabond(ia2,ib3)
                          if(ia2ref.ne.ia1.and.ia2ref.ne.ia1ref) then

                             a2(1)=xa(ia2ref)-xa(ia2)
                             a2(2)=ya(ia2ref)-ya(ia2)
                             a2(3)=za(ia2ref)-za(ia2)

c Calculate the vector perpendicular to the plane formed by these two vectors

                             call vecpro(a2,b,c2perp)

c Normalise this vector

                             amp=sqrt(scapro(c2perp,c2perp))
                             do i=1,3
                                c2perp(i)=c2perp(i)/amp
                             end do

c Cosine of the angle between them

                             cosphi=scapro(c1perp,c2perp)
                             if(cosphi.gt.1.0) cosphi=1.0
                             if(cosphi.lt.-1.0) cosphi=-1.0

c Sine of the angle between them - see torangles for this

                             call vecpro(c1perp,c2perp,csign)
                             sinphi=scapro(csign,b)/sqrt(scapro(b,b))
                             angle=piconv*atan2(sinphi,cosphi)

                             ndihedrals=ndihedrals+1
                             if(ndihedrals.gt.mstarlabel) then
                                write(6,505) mstarlabel
 505                            format(
     +'readmole> No of dihedrals exceeds ',i5,' Program aborted')
                                return
                             endif
                             write(10,*) 'dihedral ',iatomlabel(ia1ref)
     *,iatomlabel(ia1),iatomlabel(ia2),iatomlabel(ia2ref),angle
                          endif
                       end do
                    endif
                 end do
              endif
           end do
        end do

c Search for all possible rotational groups

c        ia1=0
c        do while (ia1.lt.natom)

c           ia1=ia1+1

c Can only have a rotational group if the number of bonds
c to this atom is greater than 1

c           if(nbond(ia1).gt.1) then

c              ib1=0
c              do while (ib1.lt.nbond(ia1))

c                 ib1=ib1+1
c                 ia2=iabond(ia1,ib1)

c Define 2 axis atoms

c                 if(ia2.gt.ia1.and.nbond(ia2).gt.1) then

c                    write(10,*) 'rot',iatomlabel(ia1),iatomlabel(ia2)

c                 endif
c              end do
c           endif
c        end do

c Write the changelabel section for labeltype=3

       if(labeltype.eq.3) then
         do ia1=1,natom
           write(int_as_str1,'(I3.3)') iatomlabel(ia1)
           write(10,*) 'changelabel '//int_as_str1//' ',atomlabel(ia1)
         end do
       end if

c Write the initial potential parameters for each atom type

        do ia1=1,ilenatomtypelist

           if(chsymbol(atomicnolist(ia1)).eq.'H ') then

              epsilon=0
              sigma=0
              amass=2

           else if(chsymbol(atomicnolist(ia1)).eq.'O ') then

              epsilon=0.65
              sigma=3.1
              amass=chatwt(atomicnolist(ia1))

           else if(chsymbol(atomicnolist(ia1)).eq.'N ') then

              epsilon=0.7
              sigma=3.2
              amass=chatwt(atomicnolist(ia1))
           
           else if(chsymbol(atomicnolist(ia1)).eq.'C ') then

              epsilon=0.8
              sigma=3.7
              amass=chatwt(atomicnolist(ia1))
           
           else

              epsilon=0.8
              sigma=3.2
              amass=chatwt(atomicnolist(ia1))
           
           endif
           write(10,600) atomlabellist(ia1)
     *,epsilon,sigma,amass,chargelist(ia1)
     *,chsymbol(atomicnolist(ia1))
600   format(1x,'potential ',a3,1x,4(f10.5,1x),a3)

        end do

        write(10,*) 'temperature ',simtk
        write(10,*) 'vibtemp ', tempb
        write(10,*) 'angtemp ', tempang
        write(10,*) 'dihtemp ', tempdih
        write(10,*) 'ecoredcore ', ecore,dcore
        if(unitcelldefined) then
           write(10,*) 'unitcell ',alat,blat,clat,phib,thetac,phic
        else
           write(10,*) 'density ', rho
        endif
        close(10)

c Finally convert this to an ato file

        ilenfname_ato=indexlastdot(fname_ato)-1
        sometext='makemole '//fname_ato(1:ilenfname_ato)//'.mol'
        call parse(sometext,nwords,ncf,ncl,200,nstsize,nworddim)
        call makemole()

        return
        end

      subroutine readgpr()

c Subroutine to read a .gpr file from Ghemical

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'readmole.inc'
      include 'chemical_symbols.inc'

c
c read  file
c
      open(10,file=fnameato,status='old',iostat=ierr)

c Any kind of problems and the program return to calling routine
c and does nothing

      if(ierr.ne.0) then

         write(6,444) fnameato(1:ilenfnameato)
444   format(/'readmole> Could not open file : ',a)
         return

      endif

c initial temperature
      simtk=300
c initial vibrational temperature
      tempb=0.65e2
      tempang=1.0
      tempdih=10.0
c assumed values of ecore and dcore unless changed by the input. 
c These are used in the molecule generation program fmole.
      ecore=0.0
      dcore=0.0

      natom=0
      natomexpect=0
      nbondr=0
      nbondexpect=0

      do while (ierr.eq.0)

c Read a line of text

           read(10,100,iostat=ierr) sometext
100   format(a)

c parse the text into words. Any gaps of more than 200 spaces will 
c signal the end of the line

         call parse(sometext,nwords,ncf,ncl,200,nstsize,nworddim)
c
c check that the first word is not a control word

         first=sometext(ncf(1):ncl(1))
         ilenfirst=ncl(1)-ncf(1)+1
c        write(6,*) first(1:ilenfirst)
         ilenexclaim=index(first(1:ilenfirst),'!')
           if(ilenexclaim.gt.0) then

c command found

              command_text=first(ilenexclaim+1:ilenfirst)
              ncommand=ilenfirst-ilenexclaim

              if(command_text(1:ncommand).eq.'Atoms'.and.nwords.gt.1) then

                 read(sometext(ncf(2):ncl(2)),*,iostat=ierr) natomexpect

              else if(command_text(1:ncommand).eq.'Bonds'.and.nwords.gt.1) then

             read(sometext(ncf(2):ncl(2)),*,iostat=ierr) nbondexpect

            endif
           endif

           if(command_text(1:ncommand).eq.'Atoms'.and.nwords.eq.4) then

c Read the atom parameters

              if(natom.ge.mstarlabel) then
                 write(6,502) mstarlabel
 502             format('readgpr> No of atoms exceeds ',i5
     +,' readgpr aborted')
                 return
              endif
              natom=natom+1
            read(sometext(ncf(1):ncl(2)),*) iatomlabel(natom),iatomicno(natom)

c Generate the atom label

              atomlabel(natom)=sometext(ncf(1):ncl(1))
              ilenatomlabel(natom)=ncl(1)-ncf(1)+1
              atomlabel(natom)=atomlabel(natom)(1:ilenatomlabel(natom))
     *//chsymbol(iatomicno(natom))
              ilenchsymbol=index(chsymbol(iatomicno(natom)),' ')-1
              ilenatomlabel(natom)=ilenatomlabel(natom)+ilenchsymbol

c Set the number of bonds for this atom to zero

            nbond(natom)=0

         else if(command_text(1:ncommand).eq.'Bonds'.and.nwords.eq.3) then 

c Read the bonded pairs

              if(nbondr.ge.mstarlabel) then
                 write(6,503) mstarlabel
 503             format('readgpr> No of bonds exceeds ',i5
     +,' readgpr aborted')
                 return
              endif
              nbondr=nbondr+1
              read(sometext(ncf(1):ncl(2)),*) ib1,ib2

c It is assumed these numbers are 1 less than the reference integers within EPSR (which
c range from 1 to natom). 
c Set the bond numbers for each atom of the pair

            ib1=ib1+1
            ib2=ib2+1
            if(nbond(ib1).lt.near) then
             nbond(ib1)=nbond(ib1)+1
               iabond(ib1,nbond(ib1))=ib2
            endif
            if(nbond(ib2).lt.near) then
               nbond(ib2)=nbond(ib2)+1
               iabond(ib2,nbond(ib2))=ib1
            endif

         else if(command_text(1:ncommand).eq.'Coord'.and.nwords.eq.4) then 

c Read the bonded pairs

              read(sometext(ncf(1):ncl(nwords)),*) ia1,x,y,z

c Increment atom number by 1, and convert coordinates to A

              ia1=ia1+1
              xa(ia1)=x*10.0
              ya(ia1)=y*10.0
              za(ia1)=z*10.0

         else if(command_text(1:ncommand).eq.'PartialCharges'.and.nwords.eq.2) then 


c Read the bonded pairs

              read(sometext(ncf(1):ncl(nwords)),*) ia1,acharge(ia1+1)

         endif 
      end do

c End of reading .gpr file

      close(10)
      ierr=0
      return
      end

      subroutine readjmol()

c Subroutine to read a .jmol file from Jmol

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'readmole.inc'
      include 'chemical_symbols.inc'

      integer ilinelastatom,ilinelastbond,iref

      logical test
c
c read  file
c
      open(10,file=fnameato,status='old',iostat=ierr)

c Any kind of problems and the program return to calling routine
c and does nothing

      if(ierr.ne.0) then

         write(6,444) fnameato(1:ilenfnameato)
444   format(/'readmole> Could not open file : ',a)
         return

      endif

c initial temperature
      simtk=300
c initial vibrational temperature
      tempb=0.65e2
      tempang=1.0
      tempdih=10.0
c initial atom number density
      rho=0.001
c assumed values of ecore and dcore unless changed by the input. 
c These are used in the molecule generation program fmole.
      ecore=0.0
      dcore=0.0

      natom=0
      natomexpect=0
      nbondr=0
      nbondexpect=0
      ilinelastatom=0
      ilinelastbond=0
      iline=0

c No. of lines to ignore

      nignore=3

      do while (ierr.eq.0)

c Read a line of text

         read(10,100,iostat=ierr) sometext
100   format(a)

         if(ierr.eq.0) then

            iline=iline+1

c Ignore the first two lines

            if(iline.gt.nignore) then

c parse the text into words. Any gaps of more than 200 spaces will 
c signal the end of the line

               call parse(sometext,nwords,ncf,ncl,200,nstsize,nworddim)

c Third line has the number of atoms and bonds expected

               if(iline.eq.nignore+1) then

                  read(sometext,'(2i3)') natomexpect,nbondexpect
                  if(natomexpect.gt.999) then 
                     write(6,*) 'readjmol> Number of atoms ',natomexpect,' .gt. 999'
                     return
                  end if
                  ilinelastatom=iline+natomexpect
                  ilinelastbond=ilinelastatom+nbondexpect
c                  write(6,*) natomexpect,nbondexpect,ilinelastatom,ilinelastbond

               else if(iline.le.ilinelastatom) then

                  natom=iline-(nignore+1)
c Read the atom parameters

                  if(natom.ge.mstarlabel) then
                     write(6,502) mstarlabel
 502  format('readjmol> No of atoms exceeds ',i5,':  readjmol aborted.')
                     return
                  endif
                  read(sometext(ncf(1):ncl(3)),*) xa(natom),ya(natom),za(natom)
                  atomlabel(natom)=sometext(ncf(4):ncl(4)) 
                  ilenatomlabel(natom)=ncl(4)-ncf(4)+1

c Find the atomic number of this element

                  test=.false.
                  iref=0
                  do while(.not.test.and.iref.lt.melements)
                     iref=iref+1
                     test=atomlabel(natom)(1:ilenatomlabel(natom)).eq.chsymbol(iref)
                  end do
                  iatomicno(natom)=iref
c The next line is to fit with Ghemical atom numbering
                  iatomlabel(natom)=natom
                  write(sometext,*) iatomlabel(natom)
                  call parse(sometext,nwords,ncf,ncl,200,nstsize,nworddim)
                  atomlabel(natom)=sometext(ncf(1):ncl(1))
     *//chsymbol(iatomicno(natom))
                  ilenchsymbol=index(chsymbol(iatomicno(natom)),' ')-1
                  ilenatomlabel(natom)=ncl(1)-ncf(1)+1+ilenchsymbol
c      write(6,'("readjmol> ",i4,1x,a)') natom,atomlabel(natom)(1:ilenatomlabel(natom))
                  nbond(natom)=0
                  acharge(natom)=0

               else if(iline.le.ilinelastbond) then

c Read the bonded pairs

                  if(nbondr.ge.mstarlabel) then
                     write(6,503) mstarlabel
503   format('readjmol> No of bonds exceeds ',i5,': readjmol aborted')
                     return
                  endif
                  nbondr=nbondr+1
c                  read(sometext(ncf(1):ncl(2)),'(i3,i3)') ib1,ib2
                  read(sometext,'(i3,i3)') ib1,ib2
                  if(nbond(ib1).lt.near) then
                     nbond(ib1)=nbond(ib1)+1
                     iabond(ib1,nbond(ib1))=ib2
                  endif
                  if(nbond(ib2).lt.near) then
                     nbond(ib2)=nbond(ib2)+1
                     iabond(ib2,nbond(ib2))=ib1
                  endif
               endif
            endif
         endif
      end do

c End of reading .jmol file

      close(10)
      ierr=0
      return
      end

      subroutine readcerius2()

c Subroutine to read a .ms or .msi file from Cerius2 and output a corresponding .ato. In this version bonding information
c contained in the .ms file is ignored.

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'readmole.inc'
      include 'chemical_symbols.inc'

      integer iref,tempatlabel(mstarlabel),nleftbracket
     *,atom1(mstarlabel),atom2(mstarlabel)
      character*1 leftbracket

!Coordinates of a,b and c vectors
      real avec(3),bvec(3),cvec(3),aunit(3),bunit(3),cunit(3),perp(3),zaxis(3)

      logical test,foundatom1,foundatom2

      leftbracket='('
!Initialise the unit cell vectors in case they are not defined in the input file
      do i=1,3
         avec(i)=0.0
         bvec(i)=0.0
         cvec(i)=0.0
      end do
c
c read  file
c
      open(10,file=fnameato,status='old',iostat=ierr)

c Any kind of problems and the program return to calling routine
c and does nothing

      if(ierr.ne.0) then

         write(6,444) fnameato(1:ilenfnameato)
444   format(/'readmole> Could not open file : ',a)
         return

      endif

c initial temperature
      simtk=300
c initial vibrational temperature
      tempb=0.65e2
      tempang=1.0
      tempdih=10.0
      
c assumed values of ecore and dcore unless changed by the input. 
c These are used in the molecule generation program fmole.
      ecore=0.0
      dcore=0.0

      natom=0
      nbondr=0
      iline=0
      
c Read a line of text
      read(10,100,iostat=ierr) sometext
100   format(a)

      do while (ierr.eq.0)

         iline=iline+1
         call parse(sometext,nwords,ncf,ncl,200,nstsize,nworddim)
!         write(6,'(30(1x,a))') (sometext(ncf(i):ncl(i)),i=1,nwords)
! Look for lines containing the two words '(nnn Atom', where nnn is a number representing the atom
         if(nwords.eq.2.and.index(sometext(ncf(1):ncf(1)),leftbracket).eq.1.and.sometext(ncf(2):ncl(2)).eq.'Atom') then
            nleftbracket=1
            ncf(1)=ncf(1)+1
            natom=natom+1
            if(natom.ge.mstarlabel) then
               write(6,502) mstarlabel
502   format('readcerius2> No of atoms exceeds ',i5,':  readcerius2 aborted.')
               return
            endif
 !Get the listed number of this atom
            read(sometext(ncf(1):ncl(1)),*) tempatlabel(natom)
            write(6,'("Atom: ",2(1x,i5))') natom,tempatlabel(natom)
            !Need to read lines until the number of right brackets = number of left brackets
            read(10,100,iostat=ierr) sometext
            do while (ierr.eq.0.and.nleftbracket.gt.0)
               call parse_cerius2(10,nleftbracket,ierr)
!               write(6,'(30(1x,a))') (sometext(ncf(i):ncl(i)),i=1,nwords)
               if(sometext(ncf(3):ncl(3)).eq.'ACL'.and.nwords.ge.5) then !Read atom label
                  atomlabel(natom)=sometext(ncf(5):ncl(5)) 
                  ilenatomlabel(natom)=ncl(5)-ncf(5)+1

c Find the atomic number of this element

                  test=.false.
                  iref=0
                  do while(.not.test.and.iref.lt.melements)
                     iref=iref+1
                     test=atomlabel(natom)(1:ilenatomlabel(natom)).eq.chsymbol(iref)
                  end do
                  iatomicno(natom)=iref
c The next line is to fit with Ghemical atom numbering
                  iatomlabel(natom)=natom
                  write(sometext,*) iatomlabel(natom)
                  call parse(sometext,nwords,ncf,ncl,200,nstsize,nworddim)
                  atomlabel(natom)=sometext(ncf(1):ncl(1))
     *//chsymbol(iatomicno(natom))
                  ilenchsymbol=index(chsymbol(iatomicno(natom)),' ')-1
                  ilenatomlabel(natom)=ncl(1)-ncf(1)+1+ilenchsymbol
      write(6,'("readms> ",i4,1x,a)') natom,atomlabel(natom)(1:ilenatomlabel(natom))
                  nbond(natom)=0
                  acharge(natom)=0
               else if(sometext(ncf(3):ncl(3)).eq.'Charge'.and.nwords.ge.4) then !Read atom charge
                  read(sometext(ncf(4):ncl(4)),*) acharge(natom)
               else if(sometext(ncf(3):ncl(3)).eq.'XYZ'.and.nwords.ge.6) then
                  read(sometext(ncf(4):ncl(6)),*) xa(natom),ya(natom),za(natom) !Read atom xyz coordinates
               else if(sometext(ncf(3):ncl(3)).eq.'SymType') then !Ignore this atom if it is a symmetry atom
                  natom=natom-1
               end if

               read(10,100,iostat=ierr) sometext
            end do

! Look for lines containing the two words '(nnn Bond', where nnn is a number representing the bond

         else if(nwords.eq.2.and.index(sometext(ncf(1):ncf(1)),leftbracket).eq.1.and.sometext(ncf(2):ncl(2)).eq.'Bond') then
            nleftbracket=1
            ncf(1)=ncf(1)+1
            nbondr=nbondr+1
            if(nbondr.ge.mstarlabel) then
               write(6,503) mstarlabel
503   format('readms> No of bonds exceeds ',i5,': readms aborted')
               return
            endif
            write(6,'("Bond: ",2(1x,i5))') nbondr
            atom1(nbondr)=0
            atom2(nbondr)=0
            !Need to read lines until the number of right brackets = number of left brackets
            read(10,100,iostat=ierr) sometext
            do while (ierr.eq.0.and.nleftbracket.gt.0)
               call parse_cerius2(10,nleftbracket,ierr)
!               write(6,'(30(1x,a))') (sometext(ncf(i):ncl(i)),i=1,nwords)
               if(sometext(ncf(3):ncl(3)).eq.'Atom1'.and.nwords.ge.4) then !Read atom1 number label
                  read(sometext(ncf(4):ncl(4)),*) atom1(nbondr)
               else if(sometext(ncf(3):ncl(3)).eq.'Atom2'.and.nwords.ge.4) then
                  read(sometext(ncf(4):ncl(4)),*) atom2(nbondr) !Read atom2 number label
               else if(sometext(ncf(3):ncl(3)).eq.'SymType') then !Ignore this atom if it is a symmetry atom
                  nbondr=nbondr-1
               end if

               read(10,100,iostat=ierr) sometext
            end do
!Look for coordinates of a axis
         else if(nwords.ge.6.and.index(sometext(ncf(1):ncf(1)),leftbracket).eq.1.and.sometext(ncf(3):ncl(3)).eq.'A3') then
            nleftbracket=1
            call parse_cerius2(10,nleftbracket,ierr)
            read(sometext(ncf(4):ncl(6)),*) avec(1),avec(2),avec(3)
            read(10,100,iostat=ierr) sometext
!Look for coordinates of b axis
         else if(nwords.ge.6.and.index(sometext(ncf(1):ncf(1)),leftbracket).eq.1.and.sometext(ncf(3):ncl(3)).eq.'B3') then
            nleftbracket=1
            call parse_cerius2(10,nleftbracket,ierr)
            read(sometext(ncf(4):ncl(6)),*) bvec(1),bvec(2),bvec(3)
            read(10,100,iostat=ierr) sometext
!Look for coordinates of c axis
         else if(nwords.ge.6.and.index(sometext(ncf(1):ncf(1)),leftbracket).eq.1.and.sometext(ncf(3):ncl(3)).eq.'C3') then
            nleftbracket=1
            call parse_cerius2(10,nleftbracket,ierr)
            read(sometext(ncf(4):ncl(6)),*) cvec(1),cvec(2),cvec(3)
            read(10,100,iostat=ierr) sometext
         else
c Read a line of text
            read(10,100,iostat=ierr) sometext
         endif
      end do

c Set up the bond reference integers

      do ib=1,nbondr
         write(6,*) ib,atom1(ib),atom2(ib)
         !Find the atom numbers corresponding to atom1 and atom2
         ia=0
         foundatom1=.false.
         foundatom2=.false.
         do while(ia.lt.natom.and.(.not.foundatom1.or..not.foundatom2))
            ia=ia+1
            if(tempatlabel(ia).eq.atom1(ib)) then
               ib1=ia
               foundatom1=.true.
            end if
            if(tempatlabel(ia).eq.atom2(ib)) then
               ib2=ia
               foundatom2=.true.
            end if
         end do
         if(foundatom1.and.foundatom2) then
            if(nbond(ib1).lt.near) then
               nbond(ib1)=nbond(ib1)+1
               iabond(ib1,nbond(ib1))=ib2
            endif
            if(nbond(ib2).lt.near) then
               nbond(ib2)=nbond(ib2)+1
               iabond(ib2,nbond(ib2))=ib1
            endif
         end if
      end do

      write(6,'("readcerius2> No of lines, atoms and bonds found = ",i5,1x,i5,1x,i5)') iline,natom,nbondr
!      do ia=1,natom
!         write(6,*) ia,nbond(ia)
!      end do 

      !Convert lattice vectors to unit vectors if defined, calculate rotation on coordinates, and apply the rotation.

      alat=sqrt(scapro(avec,avec))
      blat=sqrt(scapro(bvec,bvec))
      clat=sqrt(scapro(cvec,cvec))
      if(alat.gt.0.0.and.blat.gt.0.0.and.clat.gt.0.0) then

         unitcelldefined=.true.
         do i=1,3
            aunit(i)=avec(i)
            bunit(i)=bvec(i)
            cunit(i)=cvec(i)
         end do
         call divvec(aunit,alat)
         call divvec(bunit,blat)
         call divvec(cunit,clat)
!Cosine and sine of a vector relative to z axis - theta
         ctheta=aunit(3)
         stheta=sqrt(aunit(1)*aunit(1)+aunit(2)*aunit(2))
         if(stheta.gt.0.0) then
            cphi=aunit(1)/stheta
            sphi=aunit(2)/stheta
         else
            cphi=1.0
            sphi=0.0
         end if
         write(6,*) ctheta,stheta,cphi,sphi
c Using cphi,sphi rotate all coordinates about the z axis. Rotation is in opposite direction to that of axis
         call rotateright(aunit(1),aunit(2),cphi,sphi)
         call rotateright(bunit(1),bunit(2),cphi,sphi)
         call rotateright(cunit(1),cunit(2),cphi,sphi)
         call rotateright(avec(1),avec(2),cphi,sphi)
         call rotateright(bvec(1),bvec(2),cphi,sphi)
         call rotateright(cvec(1),cvec(2),cphi,sphi)
         do ia=1,natom
            call rotateright(xa(ia),ya(ia),cphi,sphi)
         end do
! Using ctheta,stheta, rotate all coordinates about y axis. Note that in this case we want to make the a vector
! coincide with the x axis, but theta is measured wrt the z axis. Therefore we reverse the order of ctheta, stheta.
         call rotateright(aunit(1),aunit(3),stheta,ctheta)
         call rotateright(bunit(1),bunit(3),stheta,ctheta)
         call rotateright(cunit(1),cunit(3),stheta,ctheta)
         call rotateright(avec(1),avec(3),stheta,ctheta)
         call rotateright(bvec(1),bvec(3),stheta,ctheta)
         call rotateright(cvec(1),cvec(3),stheta,ctheta)
         do ia=1,natom
            call rotateright(xa(ia),za(ia),stheta,ctheta)
         end do
         write(6,*) (aunit(i),i=1,3)
         write(6,*) (bunit(i),i=1,3)
         write(6,*) (cunit(i),i=1,3)
! Finally we need to rotate about the x axis so that the b vector is in the xy plane
! First create a unit vector perpendicular to a and b
         call vecpro(aunit,bunit,perp)
         pmod=sqrt(scapro(perp,perp))
         if(pmod.gt.0.0) then
            call divvec(perp,pmod)
! z-axis vector
            zaxis(1)=0.0
            zaxis(2)=0.0
            zaxis(3)=1.0
            write(6,*) perp(1),perp(2),perp(3)
            ctheta=scapro(perp,zaxis)
            stheta=sqrt(1-ctheta*ctheta)
            write(6,*) ctheta,stheta
            call rotateright(aunit(3),aunit(2),ctheta,stheta)
            call rotateright(bunit(3),bunit(2),ctheta,stheta)
            call rotateright(cunit(3),cunit(2),ctheta,stheta)
            call rotateright(avec(3),avec(2),ctheta,stheta)
            call rotateright(bvec(3),bvec(2),ctheta,stheta)
            call rotateright(cvec(3),cvec(2),ctheta,stheta)
            do ia=1,natom
               call rotateright(za(ia),ya(ia),ctheta,stheta)
            end do
         end if
         write(6,*) (aunit(i),i=1,3)
         write(6,*) (bunit(i),i=1,3)
         write(6,*) (cunit(i),i=1,3)
!Convert unit cell vectors to form required by .ato file
!Moduli of vectors
         alat=sqrt(scapro(avec,avec))
         blat=sqrt(scapro(bvec,bvec))
         clat=sqrt(scapro(cvec,cvec))
!Angle between a and b
         angconv=180/(4.0*atan(1.0))
         phib=angconv*acos(scapro(aunit,bunit))
!Cosine and sine of c vector relative to z axis - theta
         ctheta=cunit(3)
         stheta=sqrt(cunit(1)*cunit(1)+cunit(2)*cunit(2))
         if(stheta.gt.0.0) then
            cphi=cunit(1)/stheta
            sphi=cunit(2)/stheta
         else
            cphi=1.0
            sphi=0.0
         end if
         thetac=angconv*acos(ctheta)
         phic=angconv*atan2(sphi,cphi)
         write(6,*) alat,blat,clat
         write(6,*) phib,thetac,phic
      end if



c End of reading .ms file

      close(10)
      ierr=0
      return
      end

      subroutine rotateright(x,y,c,s)

      !Assumes a right handed rotation of coordinate system. Left hand rotation
! is achieved by reversing order of x,y, or setting s to -s
      real x,y,c,s,xold,yold

      xold=x
      yold=y
      x=xold*c+yold*s
      y=yold*c-xold*s

      return
      end

      subroutine parse_cerius2(nin,nleftbracket,ierr)

      include 'commands.inc'

      !Removes brackets and apostrophies from Cerius input line. Also retrieves lines wrapped over several lines

      character*256 temptext
      character*1 leftbracket,rightbracket,doublequotes,backslash
      logical found,somefound
      integer nleftbracket,ierr,i,ndouble,nbackslash

      leftbracket='('
      rightbracket=')'
      doublequotes='"'
      backslash='\'
100   format(a)

      !First check if this line has a backslash (continuation character)
      nbackslash=index(sometext,backslash)
      do while (nbackslash.gt.0)
         !Get the text up to but not including the backslash
         read(nin,100,iostat=ierr) temptext
         sometext=sometext(1:nbackslash-1)//' '//temptext
         nbackslash=index(sometext,backslash)
      end do
      call parse(sometext,nwords,ncf,ncl,200,nstsize,nworddim)
      if(nwords.gt.0) then
         !Look for attribute lines - need to remove brackets and apostrophes
         ndouble=0
         do i = 1,nwords
            !brackets and apostrophes should occur only in the first and last character of each word
            !apostrophes will come before left brackets and after right brackets
            somefound=.true.
            do while(somefound)
               somefound=.false.
               found=index(sometext(ncf(i):ncf(i)),doublequotes).gt.0
               if(found) then
                  ndouble=ndouble+1
                  ncf(i)=ncf(i)+1
               end if
               somefound=somefound.or.found
               found=index(sometext(ncl(i):ncl(i)),doublequotes).gt.0
               if(found) then
                  ndouble=ndouble-1
                  ncl(i)=ncl(i)-1
               end if
               somefound=somefound.or.found
               found=index(sometext(ncf(i):ncf(i)),leftbracket).gt.0
               if(found) then
                  nleftbracket=nleftbracket+1
                  ncf(i)=ncf(i)+1
               end if
               somefound=somefound.or.found
               found=index(sometext(ncl(i):ncl(i)),rightbracket).gt.0
               if(found) then
                  nleftbracket=nleftbracket-1
                  ncl(i)=ncl(i)-1
               end if
               somefound=somefound.or.found
            end do
         end do
      end if

      return
      end

      subroutine run_mopac()

c Subroutine to run MOPAC7 on a supplied molecule. The program creates a MOPAC
c input file, then reads the resulting .OUT file to get the revised coordinates and atomic
c charges.

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'readmole.inc'
      include 'chemical_symbols.inc'

      character*256 command,mopacfile
      character*80 zmatrix(mstarlabel)
      character*5 symlabels(mstarlabel,3)
      character*5 tripletlabels(mstarlabel,3),doubletlabels(mstarlabel,2)
      integer ncharcommand,molcharge,ilenchargestring
     *,nextracharges,iextracharges,nextracoords,iextracoords,iref
      character*10 chargestring
      character*5 hamiltonian
      logical foundcharges,foundcoords,test,foundheat
c Integers that signify the bonding and angle relationships within the Mopac Z-matrix
      integer namop(mstarlabel),nbmop(mstarlabel),ncmop(mstarlabel)
c List of atom type triplets
      integer option,j,ntriplets,firsttripletatom(mstarlabel)
     +,secondtripletatom(mstarlabel),nooftriplets(mstarlabel)
     +,atomtriplet(mstarlabel),nlinez,idr
     +,ndoublets,noofdoublets(mstarlabel),atomdoublet(mstarlabel)
     +,firstdoubletatom(mstarlabel),seconddoubletatom(mstarlabel)
     +,norbitals,mopac_norbitals,mopac_natoms,nhyd,nnonhyd

c Initialise variables to stop gfortran complaining

      iextracharges=0
      iextracoords=0
      nextracharges=0
      nextracoords=0

c Check that the number of required atoms and orbitals is within the limit of the version of Mopac to be run

      sometext=system_mopac_norbitals
      read(sometext,*) mopac_norbitals
c Convert this to number of atoms
      mopac_natoms=2*(mopac_norbitals/5)
      if(natom.gt.mopac_natoms) then
         write(6,96) natom,mopac_natoms
96    format(/'run_mopac> Number of atoms, ',i4,' larger than allowed by Mopac: (',i4,').'
     +/' Either reduce the number of atoms or recompile Mopac.')
         return

      endif

c Count the number of hydrogen and non-hydrogen atoms

      nhyd=0
      nnonhyd=0
      do ia=1,natom
         if(iatomicno(ia).eq.1) then
            nhyd=nhyd+1
         else
            nnonhyd=nnonhyd+1
         end if
      end do
      write(6,55) nhyd,nnonhyd
55    format(/'run_mopac> Number of hydrogen atoms = ',i4,', number of other atoms = ',i4)
      norbitals=nhyd+4*nnonhyd
      if(norbitals.gt.mopac_norbitals) then

         write(6,97)
97    format(/'run_mopac> Too many orbitals requested (',i4,').'
     +' Either reduce the number of atoms or recompile Mopac.')
         return

      endif

      call build_zmatrix(nlinez,iatomicno,zmatrix)
c
c Get the required Hamiltonian
c
      write(6,91,advance='no')
91    format(/'run_mopac> Give the required Hamiltonian'
     +//'Options are: 1=AM1, 2=MNDO, 3=MNDO3, 4=PM3, 5=PM6 (if it exists): ')
      read(5,*) option

      if(option.eq.1) then

         hamiltonian='AM1'

      else if(option.eq.2) then

         hamiltonian='MNDO'

      else if(option.eq.3) then

         hamiltonian='MNDO3'

      else if(option.eq.5) then

         hamiltonian='PM6'

      else 

         hamiltonian='PM3'

      endif

c Specify the charge on this molecule (-2, -1, 0, 1, 2, etc.)

      write(6,104,advance='no')
104   format(/'run_mopac> Specify the charge on this molecule [0, 1, -1, etc.]: ')
      read(5,*) molcharge
c Only allow charges in the range -9 to +9!
      if(molcharge.gt.-1.and.molcharge.lt.10) then
         write(chargestring,"('CHARGE=',i1)") molcharge
         ilenchargestring=8
      else if(molcharge.gt.-10.and.molcharge.lt.0) then
         write(chargestring,"('CHARGE=',i2)") molcharge
         ilenchargestring=9
      else
         write(chargestring,"('CHARGE=0')")
         ilenchargestring=8
      endif

c Set all the bond reference integers to zero just in case there is a problem reading them

      do ia1=1,natom

         namop(ia1)=0
         nbmop(ia1)=0
         ncmop(ia1)=0

      end do

c Step through the lines of the z-matrix and get the bonding integers

      do iref=1,nlinez

         sometext=zmatrix(iref)
         call parse(sometext,nwords,ncf,ncl,200,nstsize,nworddim)

c Get the supplied bonding, bond angle and dihedral angle neighbours. This will be used
c to assign symmetries for the Mopac calculation. For
c present purpose, dihedral neighbours will be ignored since there is no easy way to define
c a symmetry for these atoms.

         read(sometext(ncf(8):ncl(8)),*) namop(iref)
         read(sometext(ncf(9):ncl(9)),*) nbmop(iref)
         read(sometext(ncf(10):ncl(10)),*) ncmop(iref)

      end do
      close(10)

c set up the bond symmetry labels to identify symmetric pairs

      do ia1=1,natom

         symlabels(ia1,1)=atomlabel(ia1)
         symlabels(ia1,2)='   '
         symlabels(ia1,3)='   '
         if(namop(ia1).gt.0) symlabels(ia1,2)=atomlabel(namop(ia1))
         if(nbmop(ia1).gt.0) symlabels(ia1,3)=atomlabel(nbmop(ia1))
c      write(6,'(3(1x,i4),3(1x,a3))') ia1,namop(ia1),nbmop(ia1),(symlabels(ia1,j),j=1,3)

      end do

c Now make a list of atom pairs based on their atom labels

      ndoublets=0
      do ia1=1,natom

         ib1=0
         test=.false.
         do while(ib1.lt.ndoublets.and..not.test)

            ib1=ib1+1
            test=
     +(symlabels(ia1,1).eq.doubletlabels(ib1,1)
     +.and.symlabels(ia1,2).eq.doubletlabels(ib1,2))
     +.or.
     +(symlabels(ia1,1).eq.doubletlabels(ib1,2)
     +.and.symlabels(ia1,2).eq.doubletlabels(ib1,1))

         end do

         if(test) then

c Increment the counter for this doublet and store the doublet number

            noofdoublets(ib1)=noofdoublets(ib1)+1
            atomdoublet(ia1)=ib1

c Otherwise generate a new doublet

         else

            ib1=ndoublets+1
            do j=1,2
               doubletlabels(ib1,j)=symlabels(ia1,j)
            end do
            noofdoublets(ib1)=1
            firstdoubletatom(ib1)=ia1
            seconddoubletatom(ib1)=namop(ia1)
            atomdoublet(ia1)=ib1
            ndoublets=ib1
         endif
      end do

c Write out the list of triplets

c      do ib1=1,ndoublets

c         write(6,'(4(1x,i4),3(1x,a3))') ib1,noofdoublets(ib1),firstdoubletatom(ib1)
c     +,seconddoubletatom(ib1),(doubletlabels(ib1,j),j=1,2)     

c      end do

c Now make list of the different triplets and assign each atom to a particular triplet. Also note how many previous
c atoms are in that triplet. If more than 1 then the atom should have the same bond length as the first atom. If more than 2 then
c the atom should have the same bond angle as the triplet.

      ntriplets=0
      do ia1=1,natom

         ib1=0
         test=.false.
         do while(ib1.lt.ntriplets.and..not.test)

            ib1=ib1+1
            test=symlabels(ia1,2).eq.tripletlabels(ib1,2)
            if(test) then

               test=
     +(symlabels(ia1,1).eq.tripletlabels(ib1,1)
     +.and.symlabels(ia1,3).eq.tripletlabels(ib1,3))
     +.or.
     +(symlabels(ia1,1).eq.tripletlabels(ib1,3)
     +.and.symlabels(ia1,3).eq.tripletlabels(ib1,1))

            endif

         end do
         if(test) then

c Increment the counter for this triplet and save the c Generate a new triplet

            nooftriplets(ib1)=nooftriplets(ib1)+1
            atomtriplet(ia1)=ib1

c Otherwise generate a new triplet

         else

            ib1=ntriplets+1
            do j=1,3
               tripletlabels(ib1,j)=symlabels(ia1,j)
            end do
            nooftriplets(ib1)=1
            firsttripletatom(ib1)=ia1
            secondtripletatom(ib1)=nbmop(ia1)
            atomtriplet(ia1)=ib1
            ntriplets=ib1
         endif
      end do

c Write out the list of triplets

c      do ib1=1,ntriplets

c         write(6,'(4(1x,i4),3(1x,a3))') ib1,nooftriplets(ib1),firsttripletatom(ib1)
c     +,secondtripletatom(ib1),(tripletlabels(ib1,j),j=1,3)     

c      end do
c
c Now prepare the MOPAC input file for the actual refinement
c
      ilenreadfile=indexlastdot(fnameato)-1
      readfile=fnameato(1:ilenreadfile)//'.dat'
      open(10,file=readfile,status='unknown')
      write(10,99) hamiltonian,chargestring(1:ilenchargestring)
99    format(1x,a5,' MMOK SYMMETRY PRECISE GEO-OK ',a)
c Write two information lines
      write(10,101) fnameato(1:indexlastdot(fnameato)-1)
      write(10,101) 'Run to minimize energy'
101   format(1x,a)
c Write the Z-matrix
      do ia1=1,nlinez
         write(10,'(a80)') zmatrix(ia1)
      end do
c Blank line
      write(10,*)
c Now write the required symmetries
c If for any atom types there is more than one atom of that type, we apply symmetry rules to these atoms
c This option is not effectivetripl when the atoms are specified via Cartesian coordinates however, so drop it for now.


      do ia1=1,natom

c Get the correspond triplet list for this atom

         idr=atomdoublet(ia1)

c Ensure the triplet has non-zero atom references

         if(seconddoubletatom(idr).gt.0
     +.and.ia1.gt.firstdoubletatom(idr)) then

c Will use this doublet to define a bond symmetry if the atom number is greater than firstdoubletatom

            write(10,106) firstdoubletatom(idr),ia1
106   format(1x,i5,1x,'1',1x,i5)

         endif

         itr=atomtriplet(ia1)
         if(secondtripletatom(itr).gt.0
     +.and.ia1.gt.firsttripletatom(itr)) then

c Will use this triplet to define a bond angle symmetry.

            write(10,107) firsttripletatom(itr),ia1
107   format(1x,i5,1x,'2',1x,i5)

         endif
      end do
c Blank line
      write(10,*)
      close(10)

c Next create the command to run Mopac with this file.
c The first time is run without refinement to create the z-matrix format

      command=system_mopac(1:ilensystem_mopac)//' '
     *//readfile(1:ilenreadfile)
      ncharcommand=ilensystem_mopac+1
     *+ilenreadfile

      write(6,103) command(1:ncharcommand)
103   format(/a)

c Submit the job

      call system(command,ierr)

c Now read the .out file to get the results

      mopacfile=readfile(1:ilenreadfile)//'.out'
      open(10,file=mopacfile,status='old',iostat=ierr)
      if(ierr.ne.0) then
         close(10)
         mopacfile=readfile(1:ilenreadfile)//'.OUT'
         open(10,file=mopacfile,status='old',iostat=ierr)
      end if
      write(6,*) ierr,mopacfile(1:70)

      foundheat=.false.
      foundcharges=.false.
      foundcoords=.false.

c Any error opening this file and the read is aborted

      do while(ierr.eq.0)

         read(10,'(a)',iostat=ierr) sometext
         if(ierr.eq.0) then

c Write out the heat of formation - this should be negative for a stable solution
            if(.not.foundheat) then
               foundheat=index(sometext,'FINAL HEAT OF FORMATION').gt.0
               if(foundheat) then

                  call parse(sometext,nwords,ncf,ncl,200,nstsize,nworddim)
                  write(6,90) sometext(ncf(1):ncl(nwords))
90    format(/'readmole>   -- MOPAC calculation --'//3x,a/)

               endif
            endif

c If the charges have not been found so far look for the specific text to signify this

            if(.not.foundcharges) then

               foundcharges=index(sometext,'NET ATOMIC CHARGES').gt.0
               if(foundcharges) then
c Set the number of extra lines to read before the list of charges is reached
                  iextracharges=0
                  nextracharges=2
               endif 

            else if(foundcharges.and..not.foundcoords
     *.and.iextracharges.lt.nextracharges) then

c ignore this line

               iextracharges=iextracharges+1

            else if(foundcharges.and..not.foundcoords
     *.and.iextracharges.ge.nextracharges
     *.and.iextracharges.lt.nextracharges+natom) then

c Split the string into its components. The third should be the atomic charge

               iextracharges=iextracharges+1
               call parse(sometext,nwords,ncf,ncl,200,nstsize,nworddim)
               iref=iextracharges-nextracharges
c      write(6,'(a)') sometext(ncf(1):ncl(nwords))
c      write(6,*) iextracharges,nextracharges,iref
               read(sometext(ncf(3):ncl(3)),*) acharge(iref)

            else if(foundcharges.and..not.foundcoords) then

c Search from the phrase signifying Cartesian Coordinates

               foundcoords=index(sometext,'CARTESIAN COORDINATES').gt.0
               if(foundcoords) then
                  iextracoords=0
                  nextracoords=3
               endif

            else if(foundcharges.and.foundcoords
     *.and.iextracoords.lt.nextracoords) then

c Ignore this line

               iextracoords=iextracoords+1

            else if(foundcharges.and.foundcoords
     *.and.iextracoords.ge.nextracoords
     *.and.iextracoords.lt.nextracoords+natom) then    

c Split the string into its components. The coordinates are given by the third, fourth and fifth values

               iextracoords=iextracoords+1
               call parse(sometext,nwords,ncf,ncl,200,nstsize,nworddim)
               iref=iextracoords-nextracoords
c      write(6,'(a)') sometext(ncf(1):ncl(nwords))
c      write(6,*) iextracoords,nextracoords,iref
               read(sometext(ncf(3):ncl(3)),*) xa(iref)
               read(sometext(ncf(4):ncl(4)),*) ya(iref)
               read(sometext(ncf(5):ncl(5)),*) za(iref)

            else if(foundcharges.and.foundcoords
     *.and.iextracharges.gt.nextracharges+natom
     *.and.iextracoords.gt.nextracoords+natom) then

c If all is read satisfactorily, no need to read any more.

               ierr=1

            endif
         endif
      end do
      close(10)

      if(.not.foundheat.or..not.foundcharges.or..not.foundcoords) then

         write(6,250)
250   format(/'run_mopac>   !! MOPAC calculation probably failed !!'
     +//'   Check the relevant .out file for more information.'
     +//'   You could try a different Hamiltonian, or'
     +/'   give the structure a shake with Jmol.')

      endif

      return
      end

      subroutine build_zmatrix(nlinez,iatomicno,zmatrix)

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'chemical_symbols.inc'

c Attempt to build the Mopac Z-matrix

      parameter(mstarlabel=10000)
      integer iatomicno(*)
      character*80 zmatrix(*)
      integer nlinez,atominzlist(mstarlabel)
      integer atomline(mstarlabel),bondline(mstarlabel),bangleline(mstarlabel),dihedline(mstarlabel)
      integer ia,ib,ic,id,ie,iaref,ibref,il
      logical test

      do ia=1,natom

         atominzlist(ia)=0
         atomline(ia)=0
         bondline(ia)=0
         bangleline(ia)=0
         dihedline(ia)=0
c      write(6,*) ia,nbond(ia),(iabond(ia,ib),ib=1,nbond(ia))

      end do

      ia=0
      ib=0
      nlinez=0
      do while(ia.lt.natom)

         ia=ia+1
         test=atominzlist(ia).eq.0
         if(test) then

c If not in the list, add this atom to the list

            nlinez=nlinez+1
            atomline(nlinez)=ia
            atominzlist(ia)=nlinez
c Step through the neighbours of ia. The first that is already in the list will form its bond
! This will only work if there at least 2 lines in the list
            if(nlinez.gt.1) then
            ib=0
            test=.false.
            do while(ib.lt.nbond(ia).and..not.test)

               ib=ib+1
               test=atominzlist(iabond(ia,ib)).gt.0
               if(test) bondline(nlinez)=iabond(ia,ib)

            end do
            end if

c Now search for an atom to make a bond angle with. Ideally this should be one of the neighbours of bondline
c that is already in the list, but obviously cannot be itself!
c If no suitable candidates are found, we step backwards down the list and use the first atom we find, providing
c it is not the current atom nor the bonded atom
! This will only work if there are at least 3 lines in the list

            if(nlinez.gt.2) then
            il=nlinez
            test=.false.
            do while(il.gt.0.and..not.test)

c Starting with the current line we see if any of the neighbours of bondline can be used to 
c form a bond angle. If not we go the next previous line and use the atom listed for that line

               if(il.eq.nlinez) then
                  iaref=bondline(il)
                  ib=0
                  test=.false.
                  do while(ib.lt.nbond(iaref).and..not.test)
                     ib=ib+1
                     test=ia.ne.iabond(iaref,ib).and.atominzlist(iabond(iaref,ib)).gt.0
                     if(test) bangleline(nlinez)=iabond(iaref,ib)
                  end do
               else
                  iaref=atomline(il)
                  test=iaref.ne.bondline(nlinez)
                  if(test) bangleline(nlinez)=iaref
               endif
               il=il-1

            end do
            end if

c For the dihedral angle we first look at the neighbours of bangleline. If none are suitable then
c use the nearest atom in the list that is not the present one, nor bonded to it, nor forms a dihedral angle with it.

            if(nlinez.gt.3) then

               il=nlinez
               test=.false.
               do while(il.gt.0.and..not.test)

                  if(il.eq.nlinez) then

                     iaref=bangleline(il)
                     test=iaref.eq.0
                     ib=0
                     do while(ib.lt.nbond(iaref).and..not.test)
                        ib=ib+1
                        ibref=iabond(iaref,ib)
                        test=ibref.ne.ia.and.atominzlist(ibref).gt.0.and.ibref.ne.bondline(nlinez)
                        if(test) dihedline(nlinez)=ibref
                     end do

                  else

                     iaref=atomline(il)
                     test=iaref.ne.bondline(nlinez).and.iaref.ne.bangleline(nlinez)
                     if(test) dihedline(nlinez)=atomline(il)
 
                  endif
                  il=il-1

               end do
            endif
         endif
      end do

      do ia=1,nlinez

         ib=atomline(ia)
         ic=bondline(ia)
         id=bangleline(ia)
         ie=dihedline(ia)
         if(ia.eq.1) then
            write(zmatrix(ia),97) chsymbol(iatomicno(ib)),bondlength(ib,ic),bangle(ib,ic,id),dihedangle(ib,ic,id,ie),ic,id,ie
97    format(1x,a3,3(1x,e13.6,' 0 '),3(1x,i4))
         else if(ia.eq.2) then
            write(zmatrix(ia),98) chsymbol(iatomicno(ib)),bondlength(ib,ic),bangle(ib,ic,id),dihedangle(ib,ic,id,ie),ic,id,ie
98    format(1x,a3,(1x,e13.6,' 1 '),2(1x,e13.6,' 0 '),3(1x,i4))
         else if(ia.eq.3) then
            write(zmatrix(ia),99) chsymbol(iatomicno(ib)),bondlength(ib,ic),bangle(ib,ic,id),dihedangle(ib,ic,id,ie),ic,id,ie
99    format(1x,a3,2(1x,e13.6,' 1 '),(1x,e13.6,' 0 '),3(1x,i4))
         else
            write(zmatrix(ia),100) chsymbol(iatomicno(ib)),bondlength(ib,ic),bangle(ib,ic,id),dihedangle(ib,ic,id,ie),ic,id,ie
100   format(1x,a3,3(1x,e13.6,' 1 '),3(1x,i4))
         endif
c         write(6,101) zmatrix(ia)
101   format(1x,a)

      end do
      return
      end

      function bondlength(ia1,ia2)

      include 'epsr_dimension.inc'
      include 'molecom.inc'

      integer ia1,ia2
      real xdif,ydif,zdif,bondlength

      if(ia1.eq.0.or.ia2.eq.0) then

         bondlength=0.0
         return

      endif

c Calculates the bond length between two atoms

      xdif=xa(ia2)-xa(ia1)
      ydif=ya(ia2)-ya(ia1)
      zdif=za(ia2)-za(ia1)
      bondlength=sqrt(xdif*xdif+ydif*ydif+zdif*zdif)

      return
      end

      function bangle(ia1,ia2,ia3)

      include 'epsr_dimension.inc'
      include 'molecom.inc'

      integer ia1,ia2,ia3
      real bangle,a(3),b(3)

      if(ia1.eq.0.or.ia2.eq.0.or.ia3.eq.0) then

         bangle=0.0
         return

      endif

c Calculates the angle between atoms ia1, ia2 and ia3, with ia2 in the centre

      pi=4.0*atan(1.0)
      piconv=180/pi

      a(1)=xa(ia1)-xa(ia2)
      a(2)=ya(ia1)-ya(ia2)
      a(3)=za(ia1)-za(ia2)
      b(1)=xa(ia3)-xa(ia2)
      b(2)=ya(ia3)-ya(ia2)
      b(3)=za(ia3)-za(ia2)

      bangle=piconv*acos(scapro(a,b)
     */sqrt(scapro(a,a)*scapro(b,b)))

      return
      end

      function dihedangle(ia1,ia2,ia3,ia4)

c Calculates the dihedral angle formed by ia1,ia2,ia3 and ia4. 
c The axis of rotation is defined by the vector b: ia2 > ia3.

      include 'epsr_dimension.inc'
      include 'molecom.inc'

      integer ia1,ia2,ia3,ia4
      real dihedangle,a1(3),a2(3),b(3),c1perp(3),c2perp(3),csign(3)

      if(ia1.eq.0.or.ia2.eq.0.or.ia3.eq.0.or.ia4.eq.0) then

         dihedangle=0.0
         return

      endif

      pi=4.0*atan(1.0)
      piconv=180/pi

c Define rotation axis atoms

      b(1)=xa(ia3)-xa(ia2)
      b(2)=ya(ia3)-ya(ia2)
      b(3)=za(ia3)-za(ia2)

c Vector from ia2 to ia1

      a1(1)=xa(ia1)-xa(ia2)
      a1(2)=ya(ia1)-ya(ia2)
      a1(3)=za(ia1)-za(ia2)

c Vector perpendicular to this plane

      call vecpro(a1,b,c1perp)

c Normalise this vector

      call normvec(c1perp)

c Vector from ia3 to ia4

      a2(1)=xa(ia4)-xa(ia3)
      a2(2)=ya(ia4)-ya(ia3)
      a2(3)=za(ia4)-za(ia3)

c Vector perpendicular to this plane

      call vecpro(a2,b,c2perp)

c Normalise this vector

      call normvec(c2perp)

c Cosine of the angle between them

      cosphi=scapro(c1perp,c2perp)
      if(cosphi.gt.1.0) cosphi=1.0
      if(cosphi.lt.-1.0) cosphi=-1.0

c Sine of the angle between them - see torangles for this

      call vecpro(c1perp,c2perp,csign)
      sinphi=scapro(csign,b)/sqrt(scapro(b,b))
      dihedangle=piconv*atan2(sinphi,cosphi)

      return
      end

      subroutine reorder_atomlist()

c Program to reorder the atom list so that neighbours occur consecutively.

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'readmole.inc'
      include 'swapatoms.inc'

      integer ia,ib,il,ilswap
      logical test

c First write out the current atom list and bond order

      do ia=1,natom

c         write(6,*) ia,atomlabel(ia),nbond(ia),(iabond(ia,ib),ib=1,nbond(ia))

c initialise the atomlist

         atominline(ia)=ia
         lineofatom(ia)=ia
         atominlinewanted(ia)=0
         lineofatomwanted(ia)=0

      end do

c Need to ensure each member of the list is bonded to at least one member already in the list.

      iline=0
      nlines=1
c atominline(iline) is the atom which occurs at line iline. We seed this with the first atom.
      atominlinewanted(nlines)=1
c lineofatom(ia) is the line at which atom ia occurs
      lineofatomwanted(1)=nlines
      do while(iline.lt.nlines)

         iline=iline+1
         ia=atominlinewanted(iline)
         test=lineofatomwanted(ia).eq.0
         if(test) then

c If not in the list, add this atom to the list

            nlines=nlines+1
            atominlinewanted(nlines)=ia
            lineofatomwanted(ia)=nlines

         endif

c Now step through the neighbours of ia. If not in the list, add them to it.

         ib=0
         do while(ib.lt.nbond(ia))

            ib=ib+1
            test=lineofatomwanted(iabond(ia,ib)).eq.0
            if(test) then

               nlines=nlines+1
               atominlinewanted(nlines)=iabond(ia,ib)
               lineofatomwanted(iabond(ia,ib))=nlines

            endif
         end do
      end do

c Write out the new list

c      write(6,*) natom,nlines
      do il=1,nlines

         if(atominline(il).ne.atominlinewanted(il)) then

c         write(6,*) il,atominline(il),atominlinewanted(il),lineofatom(il),lineofatomwanted(il)
c if the atom number at this line (atominline(il)) is not what is wanted, we swap the current line with the line at which
c that atom does occur
            ilswap=lineofatom(atominlinewanted(il))
            call swaplines(il,ilswap)
c         write(6,*) il,atominline(il),atominlinewanted(il),lineofatom(il),lineofatomwanted(il)

         endif
      end do

c Finally set the bond reference integers

      do il=1,nlines

         do ib=1,nbond(il)
            iabond(il,ib)=lineofatom(iabond(il,ib))
         end do
c         write(6,*) il,atomlabel(il),atominline(il),lineofatom(il),nbond(il),(iabond(il,ib),ib=1,nbond(il))
         write(6,*) il,atomlabel(il),nbond(il),(iabond(il,ib),ib=1,nbond(il))

      end do

      return
      end

      subroutine swaplines1(tia1,tia2)

c Swaps the atom numbers associated with two lines ia1 and ia2

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'readmole.inc'
      include 'swapatoms.inc'

c Input values

      integer tia1,tia2,atominlinet


c Setup the new line references

c      write(6,*) tia1,tia2,atominline(tia1),atominline(tia2),lineofatom(atominline(tia1)),lineofatom(atominline(tia2))
      lineofatom(atominline(tia1))=tia2
      lineofatom(atominline(tia2))=tia1
      atominlinet=atominline(tia1)
      atominline(tia1)=atominline(tia2)
      atominline(tia2)=atominlinet
c      lineofatom(atominline(tia1))=tia2
c      lineofatom(atominline(tia2))=tia1
c      write(6,*) tia1,tia2,atominline(tia1),atominline(tia2),lineofatom(atominline(tia1)),lineofatom(atominline(tia2))

      return
      end      

      subroutine swaplines(tia1,tia2)

c Swaps all labels and coordinates associated with two lines ia1 and ia2.
c BUT we don't swap iatomlabel, since these refer to the line order, not the atoms themselves

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'readmole.inc'
      include 'swapatoms.inc'

c Input values

      integer tia1,tia2,atominlinet

c Temporary values

      character*5 atomlabelt
      integer atomtypet,nbondt,iabondt(near),ilenatomlabelt,iatomicnot
      real xt,yt,zt
      real*8 acharget
      if(tia1.ne.tia2) then


c Setup the new line references

         lineofatom(atominline(tia1))=tia2
         lineofatom(atominline(tia2))=tia1
         atominlinet=atominline(tia1)
         atominline(tia1)=atominline(tia2)
         atominline(tia2)=atominlinet

c Save values in tia1

         atomlabelt=atomlabel(tia1)
         atomtypet=atomtype(tia1)
         nbondt=nbond(tia1)
         do ib=1,nbondt
            iabondt(ib)=iabond(tia1,ib)
         end do
         ilenatomlabelt=ilenatomlabel(tia1)
         iatomicnot=iatomicno(tia1)
         xt=xa(tia1)
         yt=ya(tia1)
         zt=za(tia1)
         acharget=acharge(tia1)
 
c Set new values for tia1

         atomlabel(tia1)=atomlabel(tia2)
         atomtype(tia1)=atomtype(tia2)
         nbond(tia1)=nbond(tia2)
         do ib=1,nbond(tia1)
            iabond(tia1,ib)=iabond(tia2,ib)
         end do
         ilenatomlabel(tia1)=ilenatomlabel(tia2)
         iatomicno(tia1)=iatomicno(tia2)
         xa(tia1)=xa(tia2)
         ya(tia1)=ya(tia2)
         za(tia1)=za(tia2)
         acharge(tia1)=acharge(tia2)

c Set new values for tia2


         atomlabel(tia2)=atomlabelt
         atomtype(tia2)=atomtypet
         nbond(tia2)=nbondt
         do ib=1,nbond(tia2)
            iabond(tia2,ib)=iabondt(ib)
         end do
         ilenatomlabel(tia2)=ilenatomlabelt
         iatomicno(tia2)=iatomicnot
         xa(tia2)=xt
         ya(tia2)=yt
         za(tia2)=zt
         acharge(tia2)=acharget

      end if

      return
      end      


