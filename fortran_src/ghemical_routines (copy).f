      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call readgpr()

      stop
      end      

      subroutine readgpr()
c
c reads a ghemical .gpr and creates a .mol template file
c
      parameter (mstarlabel=10000)
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
c     include 'clustercom.inc'

      logical test

      character*256 readfile,first
      character*4 ext

c Character string to represent the atom labels

      character*3 atomlabel(mstarlabel)
      character*3 atomlabellist(mcomp)

c Labels for bondlist

        character*7 bondlabel,bondlabellist(mstarlabel)
      character*11 bondanglelabel,bondanglelabellist(mstarlabel)

c List of atomic numbers used in this file

        integer*4 atomicnolist(mcomp),ilenatomtypelist

      integer*4 ilenreadfile,ilenfirst,ilenatomlabel(mstarlabel)

            integer*4 natomexpect,nbondexpect
     +,nbondlabel(mstarlabel),nbondlist
     +,nbondanglelabel(mstarlabel),nbondanglelist

c Reference integers for the atom label and atomic number

            integer*4 iatomlabel(mstarlabel),iatomicno(mstarlabel)

c Reference integers for atoms

      integer*4 ia1,ib1,ib2,atomtype(mstarlabel)

c Coordinates of individual atoms

        real*8 x,y,z,acharge(mstarlabel),chargelist(mcomp)

c List of bonds

      real*4 bondlist(mstarlabel),bondanglelist(mstarlabel)

c Arrays to store vectors

        real*4 a(3),a1(3),a2(3),b(3),c1perp(3),c2perp(3),csign(3)

      external cosd,sind

      include 'chemical_symbols.inc'
c
c initialise the type counters for each element
c
        do i=1,melements
           counttype(i)=0
        end do
c
c generate the name of the .gpr file
c
      ext='.gpr'
      ilenext=4
      if(nwords.lt.2) then
         nwords=0
         do while (nwords.le.0)
              call search_for_files(ext,ilenext,iaccept)
              call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
              if(nwords.le.0) then
                 write(6,499)
 499    format(/'makemole> No .mol file specified! - routine will exit'/)
                 return
              endif
         end do
         readfile=sometext(ncf(1):ncl(1))
         ilenreadfile=ncl(1)-ncf(1)+1
      else
         readfile=sometext(ncf(2):ncl(2))
         ilenreadfile=ncl(2)-ncf(2)+1
      endif
c
c strip out any extension that was supplied - extension must be .ato, otherwise
c nothing is done.
c
      ilenuptodot=index(readfile,'.')-1
      if(ilenuptodot.gt.0) then
         readfile=readfile(1:ilenuptodot)
         ilenreadfile=ilenuptodot
      endif

c
c first stage is to read the molecular template file and generate
c the atoms and bonds reference integers
c
      fnameato=directory_name(1:nchardir)
     *//readfile(1:ilenreadfile)//ext(1:ilenext)
      ilenfnameato=nchardir+ilenreadfile+ilenext
c
c read .gpr file
c
      open(10,file=fnameato,status='old',iostat=ierr)

c Any kind of problems and the program return to calling routine
c and does nothing

      if(ierr.ne.0) then

         write(6,444) fnameato(1:ilenfnameato)
444   format(/'readgpr> Could not open .gpr file : ',a)
         return

      endif

c initial temperature
      simtk=300
c initial vibrational temperature
      tempb=65.0
c initial atom number density
      rho=0.001
c assumed values of ecore and dcore unless changed by the input. 
c These are used in the molecule generation program fmole.
      ecore=1.0
      dcore=3.0

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
     +,' Program aborted')
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
     +,' Program aborted')
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
        close(10)

c End of reading .gpr file

c Set up the name of the output file

      fnameato=directory_name(1:nchardir)
     *//readfile(1:ilenreadfile)//'.mol'
      ilenfnameato=nchardir+ilenreadfile+4

        open(10,file=fnameato,status='unknown',iostat=ierr)
        if(ierr.ne.0) then

           write(6,501) fnameato(1:ilenfnameato)
 501       format('readgpr> Output file ',a,' could not be opened')
           return

        endif

c Signal that this is a gmol file and that the atom numbers listed 
c are 1 less than EPSR atom numbers

      iatomlabeladd=1
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

c Generate new atom labels to reduce the number of labels
c These are based on the chemical symbol and the charge on the atom
c Also set up the list of atomic numbers for potential parameter output

        ilenatomtypelist=0
        test=.true.

        do ia1=1,natom

           ib1=0
           test=.true.
           do while (ib1.lt.ilenatomtypelist.and.test)

              ib1=ib1+1
              test=atomicnolist(ib1).ne.iatomicno(ia1)
     +.or.chargelist(ib1).ne.acharge(ia1)

           end do

           if(test) then

              if(ilenatomtypelist.lt.mcomp) then

                 ib1=ilenatomtypelist+1
                 atomicnolist(ib1)=iatomicno(ia1)
                 chargelist(ib1)=acharge(ia1)
                 counttype(iatomicno(ia1))=counttype(iatomicno(ia1))+1
             if(index(chsymbol(iatomicno(ia1)),' ')-1.eq.2) then

                write(atomlabellist(ib1),201,iostat=ierr)
     + counttype(iatomicno(ia1)),chsymbol(iatomicno(ia1))
 201            format(i1,a2)

             else

                write(atomlabellist(ib1),202,iostat=ierr)
     + counttype(iatomicno(ia1)),chsymbol(iatomicno(ia1))
 202            format(i2,a1)

             endif

c If an error occurred in writing the label, then abort the program

             if(ierr.ne.0) then

                write(6,506) counttype(iatomicno(ia1))
     +,chsymbol(iatomicno(ia1)) 
 506            format('readgpr> Error creating atom label'
     +,i5,a2,/' Probably too many atom types -'
     +,' Will abort program')
                return

             endif
                 ilenatomtypelist=ib1

              else

                 write(6,500) mcomp
 500             format('readgpr> More than ',i3
     +,' components required - program aborted')
                 return

              endif
           endif
           atomlabel(ia1)=atomlabellist(ib1)
         atomtype(ia1)=ib1

        end do

c Write out the atom coordinates and bond pairs,

        do ia1=1,natom

           write(10,*) 'atom ',iatomlabel(ia1),' ',atomlabel(ia1),' '
     *,xa(ia1),ya(ia1),za(ia1),nbond(ia1)
     *,(iatomlabel(iabond(ia1,j)),j=1,nbond(ia1))

        end do

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

                bondlabel=atomlabel(ia1)//' '//atomlabel(ia2)

             else

                bondlabel=atomlabel(ia2)//' '//atomlabel(ia1)

             endif

c Search the bond list to see if this label has already been defined

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

                 bondlist(ib2)=bondlist(ib2)+tbond(ia1,ib1)
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
           do while (ib1.lt.nbond(ia1))

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
     */sqrt(scapro(a,a)*scapro(b,b)))

c Setup the bondangle label

             if(atomtype(iaref1).le.atomtype(iaref2)) then

                bondanglelabel=atomlabel(iaref1)
     +//' '//atomlabel(ia1)
     +//' '//atomlabel(iaref2)

             else

                bondanglelabel=atomlabel(iaref2)
     +//' '//atomlabel(ia1)
     +//' '//atomlabel(iaref1)

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
     +'readgpr> No of dihedrals exceeds ',i5,' Program aborted')
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

        ia1=0
        do while (ia1.lt.natom)

           ia1=ia1+1

c Can only have a rotational group if the number of bonds
c to this atom is greater than 1

           if(nbond(ia1).gt.1) then

              ib1=0
              do while (ib1.lt.nbond(ia1))

                 ib1=ib1+1
                 ia2=iabond(ia1,ib1)

c Define 2 axis atoms

                 if(ia2.gt.ia1.and.nbond(ia2).gt.1) then

                    write(10,*) 'rot',iatomlabel(ia1),iatomlabel(ia2)

                 endif
              end do
           endif
        end do

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
           write(10,*) 'potential ',atomlabellist(ia1),' '
     *,epsilon,sigma,amass,chargelist(ia1)
     *,' ',chsymbol(atomicnolist(ia1))

        end do

        write(10,*) 'temperature ',simtk
        write(10,*) 'vibtemp ', tempb
        write(10,*) 'density ', rho
        write(10,*) 'ecoredcore ', ecore,dcore

        close(10)


        return
        end


