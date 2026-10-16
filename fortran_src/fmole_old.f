      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call fmole()

      stop
      end      

      subroutine fmole()

      use molecom_mod
      use cellneighbourlist
      
c perform a molecular "shake" on all the molecules in a cluster. This is repeated 
c the specified number of times
c
c 2006-11-12 Made a small change to initpt_fmole routine, so that truncation function 
c applied to potential, not gradient of potential.
c
      
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      character*256 readfile
      character*4 ext
      character*1 iaccept
      integer*4 nmoltypesave,ilenreadfile,num_threads
      character*256 moltemplatesave(mcomp)

c
c generate the name of the .ato file
c
      ext='.ato'
      ilenext=4
      if(nwords.lt.2) then
         nwords=0
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
499   format(/'fcluster> No .ato file specified! - routine will exit'/)
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
c determine the number of times to do the shake
c
      ntimes=0
      ierr=1
      if(nwords.gt.2) then
         read(sometext(ncf(3):ncl(3)),*,iostat=ierr) ntimes
      endif
      do while(ierr.ne.0)
         write(6,497)
497   format(/'fmole> Type the number of times to perform the shake: ',$)
         read(5,*,iostat=ierr) ntimes
      end do
c
c determine the frequency for neighbour list update
c
      nrenew=0
      ierr=1
      if(nwords.gt.3) then
         read(sometext(ncf(4):ncl(4)),*,iostat=ierr) nrenew
      endif
      do while(ierr.ne.0)
         write(6,496)
496   format(/'fmole> Type the frequency to update the neighbour list (0): ',$)
         read(5,*,iostat=ierr) nrenew
      end do
      if(nrenew.lt.1.or.nrenew.gt.ntimes) nrenew=ntimes
c
c strip out any extension that was supplied for the .ato file 
c - extension must be .ato, otherwise nothing is done.
c
      ilendot=index(readfile,'.')-1
      if(ilendot.gt.0) then
         readfile=readfile(1:ilendot)
         ilenreadfile=ilendot
      endif
      write(6,*) readfile(1:ilenreadfile)
c
c create the full name of the input file by adding on the directory name if not present
c
      call get_full_filename(readfile,ilenreadfile,fnameato,ilenfnameato)
      fnameato=fnameato(1:ilenfnameato)//'.ato'
      ilenfnameato=ilenfnameato+4
      write(6,*) readfile(1:ilenreadfile)
      write(6,*) fnameato(1:ilenfnameato)
      
c Save this for the second read

      readfile=fnameato
      ilenreadfile=ilenfnameato

c initialise number of atom and molecular types to zero

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
599   format(/'fmole> Problems with .ato file ',a/'Routine will exit now'/)
         return
      endif

c Save the names of the template files.

      nmoltypesave=nmoltype
      do it=1,nmoltype
         moltemplatesave(it)=moltemplate(it)
      end do

c Save the last values of tempb, tempdih, ecore, and dcore read in

      tempbsave=tempb
      tempdihsave=tempdih
      ecoresave=ecore
      dcoresave=dcore

c For each template, run makemole to ensure the corresponding .ato has the latest
c values.

      do it=1,nmoltypesave

         ilenmoltemplate=index(moltemplatesave(it),' ')-1
         if(ilenmoltemplate.gt.0) then

            sometext='makemole '//moltemplatesave(it)(1:ilenmoltemplate)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.gt.1) then
               write(6,595) sometext(ncf(1):ncl(2))
595   format(/'fmole> ',a/)
               call makemole()
            endif

         endif

      end do

c Makemole may have overwritten some values of the original ato file so now need 
c to read the original ato file again.

      fnameato=readfile
      ilenfnameato=ilenreadfile
      sizefactor=1.0
      ntype=0
      nmoltype=0
      call minitc(fnameato,rmaxcell,cellmax,upol,rho,amubar,amubar2
     *,ndist,ierr)

c Update this ato file from the template files

      call update_ato(ierr)

      call normalise_coordinates
c
c initialize a repulsive exponential potential based on ecore and dcore
c - this is to push non-bonded atoms as far apart as possible
c
      rminpt=4.0
      rmaxpt=6.0
      call initpt_fmole(rminpt,rmaxpt)
c
c if stepmi is set to zero, nothing will happen
c
      if(stepmi.le.0.0) then
         write(6,598)
598   format(/'fmole> Step size is set to zero!'/)
         return
      endif
c
c save the step sizes and set the others to zero
c
      stepms=stepm
      stepris=stepri
      steprs=stepr
      stepm=0.0
      stepr=0.0
      stepri=0.0
c
c calculate the initial intramolecular energy
c
      ftest=0.0
      ftestb=0.0

c return in an error occurred in fmolec1

c     if(ierr.ne.0) return
596   format(/'fmole> Starting intramolecular energy: ',(1x,e12.5)/)
c
c force the molecular "shake"
c
      intra=1
      refine_intra=.false.
      coul_intra=.false.
      refine_intra_fac=0.0
      coul_intra_fac=0.0
      !Set the requested rejection rate to 0.8
      rejrate=0.8
      nisum=intra
      ireset=0
      nmovei(1)=0
      nmovei(2)=0
      in=nrenew
      it=0

! set the maximum separation of atoms

      sepmax=0.0
      nmolcellin=10

! Initialise cell neighbour list

      call get_date_and_time('Start initialise_cellneighbourlist#',0)
      call initialise_cellneighbourlist(nmolcellin,xadd,yadd,zadd,nmol,sepmax)
      call get_date_and_time('Finish initialise_cellneighbourlist#',1)

      do while(it.lt.ntimes)
         it=it+1
         if(in.lt.nrenew) then
            in=in+1
         else
            in=1
c
c At beginning and every nrenew iterations we renew energy of distribution
c and reset the neighbour list
c
            num_threads=0
            call fmolec1(ftest,ftestb,ierr,num_threads)
            if(ierr.ne.0) return
         endif
         if(it.eq.1) write(6,596) ftestb
c         write(6,596) ftestb

c Ensure molecules stay within box.

         do im=1,nmol
            call minimagec(xadd(im),yadd(im),zadd(im))
            if(tethered(im)) call minimagec(phix(im),phiy(im),phiz(im))
         end do
         call molfit(ftest,ftestb,nisum,intra,ireset)
         if(100*(it/100).eq.it.or.it.eq.ntimes) then
            write(6,597) it,ftestb,nmovei(1),nmovei(2)
597   format(/'fmole> Iteration ',i5,' Intramolecular energy: ',(1x,e12.5)
     */1x,'No. of moves tried: ',i7,' No. of moves rejected: ',i7)
            nmovei(1)=0
            nmovei(2)=0
         end if
      end do
c
c restore step sizes
c
      stepm=stepms
      stepri=stepris
      stepr=steprs
      tempb=tempbsave
      tempdih=tempdihsave
      ecore=ecoresave
      dcore=dcoresave
c
c write out the new molecule coordinates
c
      call unnormalise_coordinates
      call moloutc(fnameato)

      write(6,*) 'Done fmole'

      return
      end

      subroutine initpt_fmole(rminpt,rmaxpt)
      
      use molecom_mod
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      pi=3.141592653
c
c maximum value of a sum - this it to prevent overflows
c
      summax=exp(25.0)
c
c generate x scae and renormalized x scale and potential damping function at low r
c
      npt=1000
      cellst=0.5*dmax/npt
      npt2=npt-1
      nptmin=1
      nptmax=npt2
      do j=1,npt
         x(j)=(j-1)*cellst
         if(x(j).le.rminpt) nptmin=j
         if(x(j).le.rmaxpt) nptmax=j
         xstar(j)=x(j)/dmax
      end do
      do j=1,npt2
         xx=x(j)
         wfac(j)=truncate(xx,rminpt,rmaxpt)
         wfacc(j)=truncatec(xx,rmaxpt)
      end do
      wfac(npt)=0.0
c
c set the empirical potentials to zero
c
      do ip=1,npairs
            do j=1,npt
                  xav(j,ip)=0.0
                  xavs(j,ip)=0.0
                  xav1(j,ip)=0.0
                  xavs1(j,ip)=0.0
                  epot(j,ip)=0.0
                  epot1(j,ip)=0.0
                  ecoul(j,ip)=0.0
                  ecoul1(j,ip)=0.0
            end do
            irefpcof(ip)=ip
            rminex(ip)=0.0
            pcoefex(ip)=0.0
      end do
      
! Generate molecular potential (based on ecore and dcore)     
      
      power=12.0
      call generate_mp()
      
      return
      end

      subroutine update_ato(ierrout)

      use molecom_mod

c Updates the .ato file, using the latest template files. Note that if an attempt is
c made to change any of the atom types or the number of different types, nothing will be 
c changed from the point at which the error was detected onwards for that molecule.

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      integer*4 ierrout
      character*20 ext
      integer*4 ilenext
      integer*4 icurrenttemplate,imoltype,ntypetemplate,natomtemplate,ilenmoltemplate

      character*3, dimension(:), allocatable        :: alabeltemplate,typetemplate,rotgrptemplate!(near)
      character*3                                            :: labeltest,symboltest

      integer, dimension(:), allocatable             :: nbondtemplate,nrottemplate!(near)
      integer, dimension(:,:), allocatable         :: iabondtemplate,labrottemplate!(near,near)
      integer, dimension(:,:), allocatable        :: labaxitemplate!(near,2)

      real, dimension(:,:), allocatable            :: tbondtemplate,dihedralangs!(near,near)

      integer, dimension(:), allocatable   :: ndihedralangs!(near,mmix)!(ma,mmix)
      integer, dimension(:,:,:), allocatable :: labrotdihs!(2,near,near,mmix)

      if(.not.allocated(typetemplate)) allocate(typetemplate(near))

c Label of the current template file

      icurrenttemplate=0

c Step through the molecules. Identify its type and read in the correct molecular template
c if it is different from the previous molecule

c Number the rotational groups

      ngrrot=0
      nrottot=0
      ig=0

c Extension of template file

      ext='.ato'
      ilenext=4

      do im=1,nmol

         imoltype=mtype(im)

c Only try to read a template file if one is specified

         if(imoltype.ne.icurrenttemplate) then

c Set the initial number of atoms in the template file to zero.
c This is to ensure nothing is done in the event that a template file does not exist.

            natomtemplate=0

c Assign the extension to the template file

            icurrenttemplate=imoltype
            ilenmoltemplate=index(moltemplate(imoltype),' ')-1

c Check that a template file is specified

            if(ilenmoltemplate.gt.0) then

c O.K so open it and get the molecular parameters

               fname_template=moltemplate(imoltype)(1:ilenmoltemplate)//
     *ext(1:ilenext)
               ilenfname_template=ilenmoltemplate+ilenext
      write(6,597) fname_template(1:ilenfname_template)
597   format(/'update_ato> ',a)


c Convert this to the full filename

               call get_full_filename(fname_template
     *,ilenfname_template,fnametemplate,ilenfnametemplate)

c Open it and read the relevant atom types, bond lengths and rotational groups.

               open(10,file=fnametemplate,status='old',iostat=ierr)

               if(ierr.eq.0) then

c Decide whether this is a cubic or a non-cubic simulation box

                  read(10,107,iostat=ierr) sometext
107   format(a)
                  call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
                  if(nwords.ne.3) then

c Not cubic, read 3 more lines

                     read(10,*,iostat=ierr)
                     read(10,*,iostat=ierr)
                     read(10,*,iostat=ierr)

                  else

! Otherwise read 1 more lines

                     read(10,*,iostat=ierr)

                  endif
                  ntypetemplate=0


c Number of atoms in the template molecule

                  read(10,*) natomtemplate
                        if(.not.allocated(alabeltemplate)) then
                            allocate(alabeltemplate(natomtemplate),nbondtemplate(natomtemplate)
     1,rotgrptemplate(natomtemplate),ndihedralangs(natomtemplate))
                            allocate(iabondtemplate(natomtemplate,near),tbondtemplate(natomtemplate,near))
                            allocate(nrottemplate(natomtemplate),labaxitemplate(natomtemplate,2),labrottemplate(natomtemplate,near))
                            allocate(dihedralangs(near,natomtemplate),labrotdihs(2,near,natomtemplate))
                        else if(natomtemplate.gt.size(alabeltemplate)) then
                            deallocate(alabeltemplate,nbondtemplate,iabondtemplate,tbondtemplate,rotgrptemplate
     1,nrottemplate,labaxitemplate,labrottemplate,ndihedralangs,dihedralangs,labrotdihs)
                            allocate(alabeltemplate(natomtemplate),nbondtemplate(natomtemplate)
     1,rotgrptemplate(natomtemplate),ndihedralangs(natomtemplate))
                            allocate(iabondtemplate(natomtemplate,near),tbondtemplate(natomtemplate,near))
                            allocate(nrottemplate(natomtemplate),labaxitemplate(natomtemplate,2),labrottemplate(natomtemplate,near))
                            allocate(dihedralangs(near,natomtemplate),labrotdihs(2,near,natomtemplate))
                        end if

c Run through the atoms in the template files and pick up the atom types and bond
c length parameters

                  do ia=1,natomtemplate

                     read(10,100,iostat=ierr) alabeltemplate(ia)
100   format(1x,a3,1x,a3)
                     it=0
                     ic=0

c Save the number of types in this molecule. Reminder it is not possible to change 
c the atom types when updating the molecule template.

                     do while (ic.lt.ntypetemplate.and.it.eq.0)
                        ic=ic+1
                        if(alabeltemplate(ia).eq.typetemplate(ic)) it=ic
                     end do

c If it.eq.0, then no atom type of this kind has been found so far in the template

                     if(it.eq.0) then
                        ntypetemplate=ntypetemplate+1
                        typetemplate(ntypetemplate)=alabeltemplate(ia)
                     endif
                     read(10,*)
                     read(10,*) nbondtemplate(ia),(iabondtemplate(ia,ja)
     *,tbondtemplate(ia,ja),ja=1,nbondtemplate(ia))

                  end do

c Read rotational groups

                  read(10,*) nrotgrptemplate
                  if(nrotgrptemplate.gt.0) then

                     do ir=1,nrotgrptemplate

                        read(10,100,iostat=ierr) rotgrptemplate(ir)
                        read(10,*) labaxitemplate(ir,1),labaxitemplate(ir,2)
                        if(rotgrptemplate(ir).eq.'DIH') then
                           read(10,*) ndihedralangs(ir)
                           do k=1,ndihedralangs(ir)
                              read(10,*) labrotdihs(1,k,ir),labrotdihs(2,k,ir),dihedralangs(k,ir)
                           end do
                        end if
                        read(10,*) nrottemplate(ir),(labrottemplate(ir,jr),jr=1,nrottemplate(ir)) 

                     end do

                  endif

c Read the potential parameters in case they have changed. It will be 
c assumed these occur in the same order as they are found in the molecule

                  do it=1,ntypetemplate

                     read(10,100,iostat=ierr) labeltest,symboltest

c Find the type number corresponding to this label

                     ic=0
                     is=0
                     do while(ic.lt.ntype.and.is.eq.0)

                        ic=ic+1
                        if(labeltest.eq.type(ic)) is=ic

                     end do

c Read the parameters if a matching type has been found. 

                     if(is.gt.0) then

                        atsymbol(is)=symboltest
                        read(10,*) epscor(is),sigcor(is)
     *,amassa(is),charge(is),qradius(is)

                     else

c If not ignore the next line

                        read(10,*)

                     endif
  
                  end do

               else

c Signal that the specified template file was not found

                  write(6,596) fnametemplate(1:ilenfnametemplate)
596   format(/'update_ato> Could not open template file :'
     *,a)

               endif
               close(10)

            endif

         endif

c O.K. so the template file has been read. Now try to apply it to the current molecule.

c First check that the number of atoms in the template file is equal to the number
c of atoms in the current molecule

         if(nato(im).eq.natomtemplate) then

c Now step through the atoms and if the labels correspond, set the neighbour numbers
c and bond lengths

            nadd=naf(im)-1

            do ia=naf(im),nal(im)

               itref=ia-nadd

c Only change the bond parameters if the labels match

               if(alabel(ia).eq.alabeltemplate(itref)) then

                  nbond(ia)=nbondtemplate(itref)
                  ja=0
                  do while(ja.lt.nbond(ia))
                     ja=ja+1
                     iabond(ia,ja)=iabondtemplate(itref,ja)+nadd
                     tbond(ia,ja)=tbondtemplate(itref,ja)
                  end do

               endif

            end do

c Rotational groups

            ngrot(im)=nrotgrptemplate
            ngrrot=ngrrot+ngrot(im)
            call testsize_rotgrps()
            ir=0
            do while(ir.lt.ngrot(im))

               ir=ir+1
               ig=ig+1
               rotgrp(ig)=rotgrptemplate(ir)
               labaxi(ig,1)=labaxitemplate(ir,1)+nadd
               labaxi(ig,2)=labaxitemplate(ir,2)+nadd
               nrot(ig)=nrottemplate(ir)
               rotf(ig)=nrottot+1
               nrottot=nrottot+nrot(ig)
               rotl(ig)=nrottot
               call testsize_labrot()
               do j=1,nrot(ig)
                  labref=rotf(ig)+j-1
                  labrot(labref)=labrottemplate(ir,j)+nadd
               end do
               if(rotgrp(ig).eq.'DIH') then
                  ndihedralang(ig)=ndihedralangs(ir)
                  do k=1,ndihedralang(ig)
                     dihedralang(k,ig)=dihedralangs(k,ir)
                     do j=1,2
                        labrotdih(j,k,ig)=labrotdihs(j,k,ir)+nadd
                     end do
                  end do
               else
                  ndihedralang(ig)=0
               end if

            end do

         else

c Remember to update the rotational groups reference integers, in case they are
c not accessed by the template file

            ngrrot=ngrrot+ngrot(im)
            ig=ngrrot

         endif

      end do
      
      return
      end


 


