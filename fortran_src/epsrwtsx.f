      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call epsrwtsx()

      stop
      end      

      subroutine epsrwtsx()
c
c Program EPSRWTSX
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
      character*256 readfile
      character*4 ext
      character*1 iaccept
      real*8 ax(5,mcomp),bx(5,mcomp),cx(mcomp)
     *,alphax(mcomp),deltax(mcomp)
      real*8 axt(5),bxt(5),cxt
      character*9 fnameb
      character*80 fname
      character*4 chsymtest,ffsym(mcomp)
      integer*4 nrho,normx
c
c generate the name of the .ato file
c
      ext='.ato'
      fnameb='X-ray.wts'
      ilenext=4
      if(nwords.lt.2) then
         nwords=0
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
499   format(/'epsrwtsx> No .ato file specified! - routine will exit'/)
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
      ilendot=index(readfile,'.')-1
      if(ilendot.gt.0) then
         readfile=readfile(1:ilendot)
         ilenreadfile=ilendot
      endif
c
c create the full name of the input file by adding on the directory name if not present
c
      call get_full_filename(readfile,ilenreadfile,fnameato,ilenfnameato)
      fnameato=fnameato(1:ilenfnameato)//'.ato'
      ilenfnameato=ilenfnameato+4
c
c initialise number of atom and molecule types to zero
c
      ntype=0
      nmoltype=0
c
c read .ato file
c
      sizefactor=1.0
      call minitc(fnameato,rmaxcell,cellmax,upol,rho,amubar,amubar2
     *,ndist,ierr)
      if(rho.eq.0.0.or.ierr.ne.0) then
         write(6,599) fnameato
599   format(/'epsrwtsx> Problems with .ato file ',a/'Routine will exit now'/)
         return
      endif
      write(6,*)
      write(6,100)
100   format(/'epsrwtsx> Program to setup the X-ray form factors for'
     1,'for the specified .ato file.'/)
99    format(a80)

      write(6,200)
200   format(/'epsrwtsx> Is the output to be per atom (1)'
     1,', or per molecule (2)? ',$)
      read(5,*) nrho

c Check that a sensible value of nrho is specified. If not it is assumed to be
c per molecule.

      if(nrho.lt.1.or.nrho.gt.2) nrho=2

c If nrho .eq. 2, then signal per molecule by making it -1

      if(nrho.eq.2) nrho=-1

c
c get the no. of atoms per molecule ratio in case it is needed
c
      atmolrat=float(natomr)/float(nmol)
      write(6,202) (it,type(it),cfrac(it),atsymbol(it),it=1,ntype)
202   format(/'epsrwtsx> The following components were found in this file:-'
     1/1x,
     1/1x,'Component no., label,  atomic fraction, chemical symbol'
     1/1x,
     2/((5x,i5,8x,a3,5x,e12.5,5x,a3)))
203   continue
      norder=1
c
c norder = 1 means this is single dataset (X-ray form factors have no isotope 
c dependence).
c
      write(6,223)
223   format(/'epsrwtsx> Get the form factor coefficients for all'
     1,' components in the sample'/)
      do it=1,ntype
c
c cannot see "Q" type atoms
c
         if(qtypea(it).eq.0) then
206         continue
            write(6,300) type(it)
300   format(/'epsrwtsx> For component ',a3,' Type form factor symbol: ',$)
            read(5,*,iostat=ierr) chsymtest
c
c test to see if this form factor symbol exists in the table
c
            ifind=0
            do while(ifind.eq.0)
               call form_fact(chsymtest,axt,bxt,cxt,ifind,5)
               if(ifind.eq.0) then
                  write(6,209) chsymtest,type(it)
209   format(/'epsrwtsx> Did not recognise form factor symbol ',a4
     *,' for atom type ',a3//'epsrwtsx> Try typing it again: ',$)
                  read(5,*,iostat=ierr) chsymtest
               endif
            end do
            do i=1,5
               ax(i,it)=axt(i)
               bx(i,it)=bxt(i)
            end do
            cx(it)=cxt
            ffsym(it)=chsymtest

c Read the alpha and delta values to define the modified atomic form factors (MAFF)
c 0 0 means use the independent atomic form factors (IAFF)

            write(6,210)
210   format(/'epsrwtsx> Specify q and delta for MAFF '
     *,'(0 0 for IAFF): ',$)
            read(5,*,iostat=ierr) alphax(it),deltax(it)

         else

            do i=1,5
               ax(i,it)=0.0
               bx(i,it)=0.0
            end do
            alphax(it)=0.0
            deltax(it)=0.0
            ffsym(it)='nil'

         endif

      end do

c Decide whether data are normalised to single atom scattering (1) or mean scattering
c length squared.

      write(6,224)
224   format(/'epsrwtsx> Are data normalised to single atom scattering (1)'
     1,' or mean form factor squared (2)? ',$)
      read(5,*) normx

c Decide whether to use the IAFFs or the MAFFs for the normalisation

      write(6,225)
225   format(/'epsrwtsx> Use the IAFFs (1) or MAFFs (2) for normalisation? ',$)
      read(5,*,iostat=ierr) iafformaff

c Based on the values of nrho, normx and iafformaff, calculate the data type.

      ndatatype=sign(1,nrho)*(normx+2*(iafformaff-1))

c write out the total weights for the dataset(s)

      nwords=0
      do while (nwords.eq.0)
         write(6,220)
220   format(/'epsrwtsx> Type basename of file to output weights to: ',$)
         read(5,111) sometext
111   format(a)
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         fname=directory_name(1:nchardir)//sometext(ncf(1):ncl(1))
      end do
c
c find out how many characters there are in this name
c
      nchar=nchardir+ncl(1)-ncf(1)+1
      nchart=nchar+9
c
c concatenate these two filenames
c
      fname=fname(1:nchar)//fnameb
      write(6,113) fname(1:nchart)
113   format(/'epsrwtsx> Writing X-ray weights to file ',a)
      open(10,file=fname,status='unknown')

c signal whether these weights are per atom or per molecule (not used for
c X-ray weights)

      write(10,*) nrho,normx,iafformaff
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
c write the .ato file in case anything was changed
c
c     call moloutc(fnameato)
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
