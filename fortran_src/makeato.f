      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call makeato()

      stop
      end      

      subroutine makeato()

      use molecom_mod
c
c Makes a molecule or headgroup
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'

      character*256 readfile
      character*4 ext
      integer*4 ntypes(mcomp),ilenext
      external cosd,sind

c
c generate the name of the .ato file
c
      ext='.ato'
      ilenext=4
      if(nwords.lt.2) then
         nwords=0
         do while (nwords.le.0)
            write(6,100)
100   format(/'makeato> Type name of .ato file to create: ',$)
            read(5,101) sometext
101   format(a)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
499   format(/'makeato> No .ato file specified! - routine will exit'/)
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
      fnameato=directory_name(1:nchardir)//readfile(1:ilenreadfile)//'.ato'
      ilenfnameato=nchardir+ilenreadfile+4
      pi=3.141592653
      write(6,102)
102   format(/'makeato> MAKEATO will create a .ato file containing 1 molecule'
     *,/' which itself contains just 1 atom'/)
      natom=0
      nmol=1
      ntype=0
      natom1=1
      ntype1=1
!      write(6,103)
103   format(/'makeato> Give number of different atom types in molecule: ',$)
!      read(5,*) ntype
      ntype=1
      type(1)=readfile(1:3)
      do i=1,ntype
            if(i.eq.1) then
               write(6,1041) i,type(i)
1041   format(/'makeato> Atom type ',i2,' has label: ',a3)
            else
               write(6,104,advance='no') i
104   format(/'makeato> For atom type ',i2,' give it a label (up to 3 chars.): ')
               read(5,121) type(i)
121   format(a)
            end if
            
!            write(6,105)
105   format(/'makeato> How many atoms of this type are in the molecule? ',$)
!            read(5,*) ntypes(i)
            ntypes(i)=1
            write(6,106) type(i)
106   format(/'makeato> For atom type ',a3,':-') 
            write(6,107,advance='no')
107   format(/' Give its epsilon [kJ/mol]: ')
            read(5,*) epscor(i)
            write(6,108,advance='no')
108   format(/' Give its sigma [A]: ')
            read(5,*) sigcor(i)
            write(6,109,advance='no')
109   format(/' Give its atomic mass [amu]: ')
            read(5,*) amassa(i)
            write(6,110,advance='no')
110   format(/' Give its electronic charge [e] ')
            read(5,*) charge(i)
            write(6,111,advance='no')
111   format(/' Give its chemical symbol: ')
            read(5,*) atsymbol(i)
            qradius(i)=0.0
            iexchange(i)=0
            do j=1,ntypes(i)
                  natom=natom+1
                  alabel(natom)=type(i)
                  xa(natom)=0.0
                  ya(natom)=0.0
                  za(natom)=0.0
                  alabel(natom)=type(i)
                  itype(natom)=i
                  nbond(natom)=0
            end do
      end do
      write(6,112)
112   format(/'makeato> The following atoms have been defined:-'
     *//' At. no., type no., label, chemical symbol,'
     *,' (x,y,z) coordinates'/)
      do i=1,natom
         write(6,127) i,itype(i),alabel(i),atsymbol(i),xa(i),ya(i),za(i)
127   format(2(1x,i3),2(1x,a3),3(1x,f10.5))
      end do
c
c set up bond distances for pairs of atom types. If no bond distance is defined
c for a particular pair of atoms, then no constraint will be defined for
c this bond
c
!      write(6,113)
113   format(/'makeato> Set up bond distances for pairs of atoms.'/)
502   continue
!      write(6,114)
114   format(/'makeato> Give two atom numbers'
     *,' and the distance between them [A]'
     *//' (Type 0 0 0 when finished)'
     */' (N.B. if you get this wrong, or want to add or delete'
     *,' bonds,'
     */' you can make changes with changeato): ',$)
!      read(5,*) ia1,ia2,bond
      ia1=0
      ia2=0
      bond=0
      do while (ia1.gt.0.and.ia2.gt.0)
         nbond(ia1)=nbond(ia1)+1
         nbond(ia2)=nbond(ia2)+1
         iabond(ia1,nbond(ia1))=ia2
         iabond(ia2,nbond(ia2))=ia1
         tbond(ia1,nbond(ia1))=bond
         tbond(ia2,nbond(ia2))=bond
         write(6,115)
115   format(/'makeato> Give two atom numbers'
     *,' and the distance between them [A].'
     *//'(Type 0 0 0 when finished): ',$)
            read(5,*) ia1,ia2,bond
      end do
c
c write out the bonds defined for each atom
c
 !     do ia=1,natom
!            write(6,116) alabel(ia),nbond(ia)
!     1,(alabel(iabond(ia,i)),iabond(ia,i),tbond(ia,i),i=1,nbond(ia))
116   format(/1x,a3,1x,i2,/((1x,a3,1x,i2,1x,f10.5)))
!      end do
c
c set up bond angles for triplets of atom types. If no bond angle is defined
c for a particular triplet of atom types, then no constraint will be defined for
c this triplet.
c
!      write(6,117)
117   format(/'makeato> Type three atom numbers and angle between them'
     *//' (Middle atom is taken as centre).'
     */' (Type 0 0 0 0 when finished): ',$)
!      read(5,*) ia1,ia2,ia3,bangle
      ia1=0
      ia2=0
      ia3=0
      do while(ia1.gt.0.and.ia2.gt.0.and.ia3.gt.0)
c
c first check that atom ia2 has at least 2 bonds defined
c
            if(nbond(ia2).lt.2) then
                  write(6,118) ia2
118   format(/'makeato> Not enough bonds defined for atom ',i3)
                  goto 502
            endif
c
c Now check that both ia1 and ia3 have been defined as bonded to ia2
c
            nfound=0
            bond21=0.0
            bond23=0.0
            do ia=1,nbond(ia2)
                  if(iabond(ia2,ia).eq.ia1) then
                        nfound=nfound+1
                        bond21=tbond(ia2,ia)
                  else if(iabond(ia2,ia).eq.ia3) then
                        nfound=nfound+1
                        bond23=tbond(ia2,ia)
                  endif
            end do
c
c check that 2 bonds were found (if more than 2 were found it will take
c the last ones found
c
            if(nfound.lt.2) then
                  write(6,119) ia1,ia2,ia3
119   format(/'makeato> Suitable bonds between',3(1x,i3),' not found')
                  goto 502
            endif
c
c O.K. so set up distance between ia1 and ia3
c
            b13=bond21**2+bond23**2-2*bond21*bond23*cosd(bangle)
            b13=sqrt(b13)
            nbond(ia1)=nbond(ia1)+1
            nbond(ia3)=nbond(ia3)+1
            iabond(ia1,nbond(ia1))=ia3
            iabond(ia3,nbond(ia3))=ia1
            tbond(ia1,nbond(ia1))=b13
            tbond(ia3,nbond(ia3))=b13
            write(6,117)
            read(5,*) ia1,ia2,ia3,bangle
      end do
c
c write out the bonds defined for each atom
c
!      do ia=1,natom
!            write(6,116) alabel(ia),nbond(ia)
!     1,(alabel(iabond(ia,i)),iabond(ia,i),tbond(ia,i),i=1,nbond(ia))
!      end do
c
c determine if there are any rotational groups
c
      ig=0
      ngrrot=0
      ngrrotmax=0
      nrottot=0
      nrottotmax=0
! NOTE: If it is ever reinstated, the following code will need to be adjusted to take account
! of the revised method of indexing the rotational groups
!
!      write(6,120,advance='no')
!120   format(/'makeato> Type 2 atom numbers to form axis of rotational group'
!     *//' (Type 0 0 to exit when finished) ')
!      read(5,*) lab1,lab2
!      lab1=0
!      lab2=0

!      do while (lab1.gt.0.and.lab2.gt.0)
!            if(lab1.gt.natom.or.lab1.lt.1) then
 !                 write(6,126) lab1
!126   format(/'makeato> Specified atom ',a3,' does not exist!')
!                  goto 505
!            endif
!            if(lab2.gt.natom.or.lab2.lt.1) then
!                  write(6,121) lab2
!                  goto 505
!            endif
!            ig=ig+1
!            call testsize_rotgrps()
!            rotgrp(ig)='ROT'
!            labaxi(ig,1)=lab1
!            labaxi(ig,2)=lab2
!506         continue
!            write(6,122)
!122   format(/'makeato> Give no. atoms in rotational group, and their'
!     *' atom numbers: ',$)
 !           read(5,*) nrot(ig),(labrot(ig,j),j=1,nrot(ig))
c
c ensure these do not include the two atoms that form the axis!
c
!            do j=1,nrot(ig)
!                  if(labrot(ig,j).eq.labaxi(ig,1)
!     1.or.labrot(ig,j).eq.labaxi(ig,2)) then
!                        write(6,123)
!123   format(/'makeato> Atoms in rotational group cannot include'
!     1,' atoms used to define axis')
!                        go to 506
!                  endif
!            end do
!505         continue
!            write(6,120)
!            read(5,*) lab1,lab2
!      end do
      ngrot(1)=ig
c
c set some parameters ready for writing the .ato file
c
      nmol=1
      xadd(1)=0.0
      yadd(1)=0.0
      zadd(1)=0.0
      phix(1)=0.0
      phiy(1)=0.0
      phiz(1)=0.0
      nato(1)=natom
      naf(1)=1
      nal(1)=natom
      stepmi=0.1
      if(ngrot(1).gt.0) then
            stepri=1.0
      else
            stepri=0.0
      endif
      stepr=1.0
      stepm=2.0
      tempb=65
      write(6,124,advance='no')
124   format(/'makeato> Type the temperature [K]: ')
      read(5,*) simtk
      write(6,125,advance='no')
125   format(/'makeato> and atomic number density [per A**3]: ')
      read(5,*) rho
c     write(6,128,advance='no')
128   format(/'makeato> Specify ecore [kJ/mole] (used in fmole) [0.0]: ')
c     read(5,*) ecore
c     write(6,129,advance='no')
129   format(/'makeato> Specify dcore [A] (used in fmole) [0.0]: ')
c     read(5,*) dcore
!      write(6,130,advance='no')
130   format(/'makeato> Specify the vibrational weighting [65]: ')
!      read(5,*) tempb
      tempb=65
      ecore=0.0
      dcore=0.0
c
c convert density to box size
c
      volume=natom/rho
      dmax=volume**(1.0/3.0)
      cubic=.true.
c
c setup molecular size factor
c
      dmaxmol=dmax
      !Set up template file
      nmoltype=1
      moltemplate(nmoltype)=readfile(1:ilenreadfile)
      wholerot(nmoltype)=0.0
      wholetrans(nmoltype)=0.01
      call moloutc(fnameato)
      return
      end

