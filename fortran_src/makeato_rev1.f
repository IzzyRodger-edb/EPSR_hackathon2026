      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call makeato()

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
!      include 'lattice.inc'
      include 'readmole.inc'
      include 'chemical_symbols.inc'

      character*4 ext
      integer*4 ntypes(mcomp),ilenext,iref
      logical test
      external cosd,sind

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

c
c generate the name of the .mol and .ato file
c
      ext='.mol'
      ilenext=4
      if(nwords.lt.2) then
         nwords=0
         do while (nwords.le.0)
            write(6,100,advance="no")
100   format(/'makeato> Type name of .ato file to create: ')
            read(5,101) sometext
101   format(a)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
499   format(/'makeato> No .ato file specified! - routine will exit.'/)
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
c strip out any extension that was supplied - extension must be .mol.
c
      ilenuptodot=index(readfile,'.')-1
      if(ilenuptodot.gt.0) then
         readfile=readfile(1:ilenuptodot)
         ilenreadfile=ilenuptodot
      endif
      fnameato=directory_name(1:nchardir)//readfile(1:ilenreadfile)//'.mol'
      ilenfnameato=nchardir+ilenreadfile+4

      open(10,file=fnameato,status='unknown',iostat=ierr)
      if(ierr.ne.0) then

         write(6,501) fnameato(1:ilenfnameato)
 501  format('makeato> Output file ',a,' could not be opened')
         return

      endif

c Signal that this is a gmol file and that the atom numbers listed
c are 1 less than EPSR atom numbers

      write(10,*) '.gmol ',0

      write(6,102)
102   format(/'makeato> makeato will first create a .mol file containing 1 molecule'
     *,/' which itself contains just a single atom.'/
     *,/' It will then use this to create the corresponding .ato file using makemole.')
      natom=0
      ntype=1
      type(1)=readfile(1:min(3,ilenreadfile)) !The first (upto) three characters determine the atom label
      do i=1,ntype
            if(i.eq.1) then
               write(6,1041) i,type(i)
1041  format(/'makeato> Atom type ',i2,' has label: ',a3)
            else
               write(6,104,advance='no') i
104   format(/'makeato> For atom type ',i2,' give it a label (up to 3 chars.): ')
               read(5,121) type(i)
121   format(a)
            end if
            write(6,601) type(i),amassa(i)
601   format(/'makeato> Atomic weight found for atom type ',a3,' is: ',f10.5)
            ntypes(i)=1
            write(6,106) type(i)
106   format(/'makeato> For atom type ',a3,':-')
            write(6,107,advance='no')
107   format(/' Give its epsilon [kJ/mol]: ')
            read(5,*) epscor(i)
            write(6,108,advance='no')
108   format(/' Give its sigma [A]: ')
            read(5,*) sigcor(i)
!            write(6,109,advance='no')
109   format(/' Give its atomic mass [amu]: ')
!            read(5,*) amassa(i)
            write(6,110,advance='no')
110   format(/' Give its electronic charge [e] ')
           read(5,*) charge(i)
            test=.false.
            do while (.not.test)
               write(6,111,advance='no')
111   format(/' Give its chemical symbol: ')
               read(5,*) atsymbol(i)
!Get the atomic mass of this atom type
               iref=0
               do while(.not.test.and.iref.lt.melements)
                  iref=iref+1
                  test=atsymbol(i)(1:len_trim(atsymbol(i))).eq.chsymbol(iref)
               end do
               if(.not.test) then
                  write(6,602) atsymbol(i)
               end if
602   format(/'makeato> Chemical symbol ',a3,' not found. Try again.')
            end do
            amassa(i)=chatwt(iref)
            qradius(i)=0.0
            iexchange(i)=0
            do j=1,ntypes(i)
                  natom=natom+1
                  alabel(natom)=type(i)
                  xa(natom)=0.0
                  ya(natom)=0.0
                  za(natom)=0.0
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

      do ia1=1,natom
           write(10,*) 'atom ',ia1,' ',alabel(ia1),' '
     *,xa(ia1),ya(ia1),za(ia1),nbond(ia1)
      end do


      do ia1=1,natom
         write(10,600) alabel(ia1),epscor(ia1),sigcor(ia1),amassa(ia1)
     *,charge(ia1),atsymbol(ia1)
      end do
600   format(1x,'potential ',a3,1x,4(f10.5,1x),a3)

      write(10,*) 'temperature ',simtk
      write(10,*) 'vibtemp ', tempb
      write(10,*) 'angtemp ', tempang
      write(10,*) 'dihtemp ', tempdih
      write(10,*) 'ecoredcore ', ecore,dcore
      write(10,*) 'density ', rho
      close(10)

c Finally convert this to an ato file

      sometext='makemole '//readfile(1:ilenreadfile)//'.mol'
      call parse(sometext,nwords,ncf,ncl,200,nstsize,nworddim)
      call makemole()

      return
      end

