      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call xwts()

      end      

      subroutine xwts()
c
c Program xwts
c
c Sets up the X-ray form factors for the elements in the specified .ato file
c
c Form factor data are read from the following compilation:-
c
c  These data, computed for neutral atoms and ions are valid for the full  
c  range of sin(theta)/lambda from 0.0 to 6.0 A-1.This table have been     
c  developped from fits of o linear combination of five Gaussians to the   
c  values of the scattering factors tabulated in "International tables for 
c  Crystallography (1992) [Vol. C. Dordrecht: Kluwer Academic Publishers]. 
c  For a detailed discussion of the underlying validity of the             
c  information contained in these directories see, for example:            
c                                                                          
c          New Analytical Scattering Factor Functions for Free Atoms       
c                    and Ions for Free Atoms and Ions                      
c                       D. Waasmaier & A. Kirfel                           
c                    Acta Cryst. (1995). A51, 416-413                      
c                                                                          
c  fo the non-dispersive part of the atomic scattering factor is a          
c  function of the selected element and of sin(theta)/lambda, where         
c  lambda is the photon wavelengh and theta is incident angle.              
c  This function can be approximated by a function:                   
c                                                                          
c    f0[k] = c + [SUM a_i*EXP(-b_i*(k^2)) ]                           
c                i=1,5                                          
c                                                                          
c  where k = sin(theta) / lambda and c, a_i and b_i                   
c  are the coefficients tabulated in this file (in columns:                 
c  a1  a2  a3  a4  a5  c  b1  b2  b3  b4  b5                          
c
c Currently set up for 70 components (number of atom types)
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      character*256 wtsfile
      real*8 ax(5,mcomp),bx(5,mcomp),cx(mcomp)
     *,alphax(mcomp),deltax(mcomp)
      real*8 axt(5),bxt(5),cxt
      character*256 fname
      character*4 chsymtest,ffsym(mcomp)
      integer*4 nrho,normx,ncharwts
      write(6,100)

c Get parameters from input file
c
c setup the xwts symbols
c
      call setup_symbols_xwts()
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


      write(6,100)
100   format(/'xwts> Program to setup the X-ray form factors for'
     1,'for the specified .ato file.'/)
99    format(a80)
c
c Set up the base name for the weights file. This will be the same as the input file.

      iref=1
      wtsfile=fname_run(1:ilenfname_run)
      ncharwts=ilenfname_run
c
c set up whether to normalise to atoms or molecules.
c
      iref=iref+1
      read(symbol_value(iref),*) nrho

c Get type of normalisation employed

      iref=iref+1
      read(symbol_value(iref),*) normx

c Decide on whether to use IAFFs or MAFFs for normalisation

      iref=iref+1
      read(symbol_value(iref),*) iafformaff
c
c get the no. of atoms per molecule ratio in case it is needed
c
      write(6,202) (it,type(it),cfrac(it),atsymbol(it),it=1,ntype)
202   format(/'xwts> The following components were found in this file:-'
     1/1x,
     1/1x,'Component no., label,  atomic fraction, chemical symbol'
     1/1x,
     2/((5x,i5,8x,a3,5x,e12.5,5x,a3)))
203   continue
c
c norder = 1 means this is single dataset (X-ray form factors have no isotope 
c dependence).
c

c Check number of components in input file is same as atom types

      iref=iref+1
      read(symbol_value(iref),*) ncomponent
      if(ncomponent.ne.ntype) then

         write(6,235)
235   format(/'xwts> Number of components in input file different from number of atom types in .ato file')
         return

      endif

      write(6,223)
223   format(/'xwts> Get the form factor coefficients for all'
     1,' components in the sample'/)
      do it=1,ntype

c Skip over next entry in input file as it is be determined from the ato file

         iref=iref+1
c
c cannot see "Q" type atoms
c
         if(qtypea(it).eq.0) then

            iref=iref+1
            chsymtest=symbol_value(iref)
c
c test to see if this form factor symbol exists in the table
c
            ifind=0
            do while(ifind.eq.0)
               call form_fact(chsymtest,axt,bxt,cxt,ifind,5)
               if(ifind.eq.0) then
                  write(6,209,advance='no') chsymtest,type(it)
209   format(/'xwts> Did not recognise form factor symbol ',a4
     *,' for atom type ',a3//'xwts> Try typing it again: ')
                  read(5,*,iostat=ierr) chsymtest
               endif
            end do
            symbol_value(iref)=chsymtest
            do i=1,5
               ax(i,it)=axt(i)
               bx(i,it)=bxt(i)
            end do
            cx(it)=cxt
            ffsym(it)=chsymtest

c Read the alpha and delta values to define the modified atomic form factors (MAFF)
c 0 0 means use the independent atomic form factors (IAFF)

            iref=iref+1
            read(symbol_value(iref),*) alphax(it)
            iref=iref+1
            read(symbol_value(iref),*) deltax(it)

         else

            iref=iref+3
            do i=1,5
               ax(i,it)=0.0
               bx(i,it)=0.0
            end do
            alphax(it)=0.0
            deltax(it)=0.0
            ffsym(it)='nil'

         endif

      end do

c write out the total weights for the dataset(s)

      fname=directory_name(1:nchardir)//wtsfile(1:ncharwts)//'.wts'
      open(10,file=fname,status='unknown')

c signal whether these weights are per atom or per molecule (not used for
c X-ray weights)

      write(10,'(3(1x,i3))') nrho,normx,iafformaff
      do it=1,ntype
         write(10,114) type(it),ffsym(it)
114   format(1x,a3,1x,a4)

c 28/1/2006 we add alpha and delta for modified atomic form factors if needed.
c See equation 5, JCP, vol113, p 9149-9161 (2000), which however is incorrect!

         write(10,115) cx(it),(ax(i,it),bx(i,it),i=1,5)
     *,alphax(it),deltax(it)

      end do
115   format(11(1x,f11.6))
      close(10)
c
c Rewrite the input file in case anything was changed
c
      call save_files_others(14)

      return
      end

      subroutine form_fact(chsym,axt,bxt,cxt,ifind,ndim)
      include 'commands.inc'
      character*4 chsym
      character*15 fnameX
      real*8 axt(ndim),bxt(ndim),cxt

      fnameX='f0_WaasKirf.dat'

      ifind=0

c check the number of non-space characters in chsym

      ilenchsym=len_trim(chsym)
      if(ilenchsym.lt.1) return


c Open form factor file and skip some lines

      open(10,file=fnameX,status='old',iostat=ierr)

      if(ierr.ne.0) then
         write(6,301) fnameX
301   format(/'form_factor> X-ray file ',a15,' does not exist')
         return
      endif

      do while (ierr.eq.0.and.ifind.eq.0)

         read(10,'(a)',iostat=ierr) sometext

c Search sometext for the character string showing a valid data line

         if(index(sometext,'#S').gt.0) then

            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)

c The third word should contain the text being searched for

            if(sometext(ncf(3):ncl(3)).eq.chsym(1:ilenchsym)) then
               read(10,*,iostat=ierr)
               read(10,*,iostat=ierr)
               read(10,'(a)',iostat=ierr) sometext

c Check that the number of words is correct

               call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
               if(ierr.eq.0.and.nwords.eq.11) then
                  read(sometext,*) (axt(i),i=1,5),cxt,(bxt(i),i=1,5)
                  ifind=1
               endif
            endif
         endif
      end do
      close(10)
      return
      end
