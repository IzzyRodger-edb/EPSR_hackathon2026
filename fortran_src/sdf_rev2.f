      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_sdf()

      stop
      end      

      subroutine run_sdf()
c
c performs calculation of spherical harmonic coefficients from .ato file
c
c uses a specified set of atoms to define the central "molecule", then calculates the
c spatial density function (SDF) of a specified second "molecule" around this. 
c
      use sharm_mod
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'sharm.inc'
      character*50 formatout
      character*2 ifmt1,ifmt2,ifmt3,ifmt4,il1fmt,il2fmt,ilfmt,in1fmt,in2fmt
      real pi,piconv
c
c setup the SDF symbols
c
      call setup_symbols_sdf()
c
c now setup the input file
c
      nin=11
      nout=6
      call setup_input_file(nin,nout,ierrout)
c
c if an error occurred on setup, exit now
c
      if(ierrout.ne.0) return

      pi=4.0*atan(1.0)
      piconv=pi/180.0
c
c read the input file
c
      call read_input_file(nin,nout,ierrout)
c
c if an error occurred on read, exit now
c
      if(ierrout.ne.0) return
      iref=1
c
c set up number of radius values
c
      iref=iref+1
      read(symbol_value(iref),*) nz
c
c set up rmax of calculation
c
      iref=iref+1
      read(symbol_value(iref),*) zmax
c
c set up number of configurations already stored
c
      iref=iref+1
      read(symbol_value(iref),*) nsumt
c
c set up number of coefficients
c
      iref=iref+1
      read(symbol_value(iref),*) ndist
c
c set up l1 values
c
      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,1,256,nworddim)
      nll1v=nwords
      if(allocated(ll1v)) deallocate(ll1v)
      allocate(ll1v(nll1v))
      do i=1,nll1v
         read(symbol_value(iref)(ncf(i):ncl(i)),*,iostat=ierr) ll1v(i)
      end do
c
c set up l2 values
c
      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,1,256,nworddim)
      nll2v=nwords
      if(allocated(ll2v)) deallocate(ll2v)
      allocate(ll2v(nll2v))
      do i=1,nll2v
         read(symbol_value(iref)(ncf(i):ncl(i)),*,iostat=ierr) ll2v(i)
      end do
c
c set up l values
c
      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
      nllv=nwords
      if(allocated(llv)) deallocate(llv)
      allocate(llv(nllv))
      do i=1,nllv
         read(symbol_value(iref)(ncf(i):ncl(i)),*) llv(i)
      end do
c
c set up n1step
c
      iref=iref+1
      read(symbol_value(iref),*) n1step
c
c set up n2step
c
      iref=iref+1
      read(symbol_value(iref),*) n2step
c
c For the central "molecule" read the list of atom types to be used to define the 
c unit. The first in the list will be assigned to the centre of the molecule.
c
      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
      natofirst=0

c skip through the atom types to find out the type number of the specified atom

      ic=0
      do while (ic.lt.nwords.and.natofirst.lt.mcentre)
         ic=ic+1
         it=0
         ifound=0
         do while(it.lt.ntype.and.ifound.eq.0)
            it=it+1
            if(type(it).eq.symbol_value(iref)(ncf(ic):ncl(ic))) ifound=it
         end do
         if(ifound.le.0) then
            write(6,201) symbol_value(iref)(ncf(ic):ncl(ic))
     *,fname_ato(1:ilenfname_ato)
201   format(/'sdf> Atom type ',a,' does not exist in .ato file ',a)
            return
         else
            natofirst=natofirst+1
            lctyp1(natofirst)=ifound
            typefirst(natofirst)=ifound
         endif
      end do
      if(natofirst.le.0) then
         write(6,202)
202   format(/'sdf> None of the specified atom types were recognised. '
     */' - returning to main program')
         return
      else
         write(6,203) natofirst
203   format(/'sdf> No. of atoms used to define centre molecule = ',i2)
      endif
      ncent1=1

c For central molecule read the array of maximum distances between atoms in the order
c (1,2) (1,3) (1,4) ... (2,3) (2,4) ... (3,4) etc.
c There should be natofirst*(natofirst-1)/2 of them

      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
      npairfirstexpec=natofirst*(natofirst-1)/2
      if(nwords.lt.npairfirstexpec) then
         write(6,204) npairfirstexpec
204   format(/'sdf> Not enough pair distances defined '
     *,'for central molecule. Need: ',i5)
         return
      endif
      i=0
      do while (i.lt.npairfirstexpec)
         i=i+1
         read(symbol_value(iref)(ncf(i):ncl(i)),*) distfirst(i)
      end do
      npairfirst=i

c Variation of distance (half width) about the specified distance for an atom to be accepted.
c This is expressed as a fraction of the specified distance. For example fraction of 0.2
c for a bond distance of 2A means atoms within 2 +- 0.4A will be accepted

      iref=iref+1
      read(symbol_value(iref),*) fractfirst

c
c for molecule 1 read primary axis type (1=X, 2=Y, 3=Z), and
c the number and atom numbers of atoms used to define the molecule's nzt1
c axis
c
      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
      nzt1=symbol_value(iref)(ncf(1):ncl(1))
      ncz1=nwords-1
      if(ncz1.gt.mcentre) ncz1=mcentre
      do i=1,ncz1
         ic=i+1
         read(symbol_value(iref)(ncf(ic):ncl(ic)),*) lztyp1(i)
      end do
c
c for molecule 1 read secondary axis type (1=X, 2=Y, 3=Z), and
c the number and atom numbers of atoms used to define the
c molecule's nyt1 axis
c
      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
      nyt1=symbol_value(iref)(ncf(1):ncl(1))
      ncy1=nwords-1
      if(ncy1.gt.mcentre) ncy1=mcentre
      do i=1,ncy1
         ic=i+1
         read(symbol_value(iref)(ncf(ic):ncl(ic)),*) lytyp1(i)
      end do

c For second "molecule" read type of atom used to define the molecular
c centre. Only the first in the list is used.

      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)

c Skip through the specified atom types to find out the type number of the specified atom
c We also test to see if this atom is the same as the corresponding atom on the central
c molecule. If not the molecules are treated as distinct.

      natosecond=0
      ldistinct=0
      do ic=1,nwords
         it=0
         ifound=0
         do while(it.lt.ntype.and.ifound.eq.0
     *.and.natosecond.lt.mcentre)
            it=it+1
            if(type(it).eq.symbol_value(iref)(ncf(ic):ncl(ic))) ifound=it
         end do
         if(ifound.eq.0) then
            write(6,201) symbol_value(iref)(ncf(ic):ncl(ic))
     *,fname_ato(1:ilenfname_ato)
            return
         else
            natosecond=natosecond+1
            lctyp2(natosecond)=ifound

         endif
      end do
      if(natosecond.le.0) then
         write(6,202)
         return
      else
         write(6,206) natosecond
206   format(/'sdf> No. of atoms used to define second molecule = ',i2)
      endif
      ncent2=1

c test for identity of central and second molecules

      ldistinct=0
      if(natosecond.eq.natofirst) then
         ldistinct=0
         do it=1,natofirst
            if(lctyp2(it).ne.lctyp1(it)) ldistinct=1
         end do
      else
         ldistinct=1
      endif

c If the second and central molecules are identical, the rest of the input for the
c second molecule is identical to the first.

      if(ldistinct.ne.0) then

c For second molecule read the array of maximum distances between atoms in the order
c (1,2) (1,3) (1,4) ... (2,3) (2,4) ... (3,4) etc.
c There should be natosecond*(natosecond-1)/2 of them

         iref=iref+1
         call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
         npairsecondexpec=natosecond*(natosecond-1)/2
         if(nwords.lt.npairsecondexpec) then
            write(6,205) npairsecondexpec
205   format(/'sdf> Not enough pair distances defined '
     *,'for second molecule. Need: ',i5)
            return
         endif
         i=0
         do while (i.lt.npairsecondexpec)
            i=i+1
            read(symbol_value(iref)(ncf(i):ncl(i)),*) distsecond(i)
         end do
         npairsecond=i

c Variation of distance (half width) about the specified distance for an atom to be accepted.
c This is expressed as a fraction of the specified distance. For example fraction of 0.2
c for a bond distance of 2A means atoms within 2 +- 0.4A will be accepted

         iref=iref+1
         read(symbol_value(iref),*) fractsecond

c for molecule 2 read primary axis type (1=X, 2=Y, 3=Z), and
c the number and atom numbers of atoms used to define the molecule's nzt1
c axis

         iref=iref+1
         call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
         nzt2=symbol_value(iref)(ncf(1):ncl(1))
         ncz2=nwords-1
         do i=1,ncz2
            ic=i+1
            read(symbol_value(iref)(ncf(ic):ncl(ic)),*) lztyp2(i)
         end do

c for molecule 2 read secondary axis type (1=X, 2=Y, 3=Z), and
c the number and atom numbers of atoms used to define the
c molecule's nyt1 axis

         iref=iref+1
         call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
         nyt2=symbol_value(iref)(ncf(1):ncl(1))
         ncy2=nwords-1
         do i=1,ncy2
            ic=i+1
            read(symbol_value(iref)(ncf(ic):ncl(ic)),*) lytyp2(i)
         end do

      else

c For identical molecules set all the parameters for the second molecule
c equal to the first molecule.

         ncent2=ncent1
         natosecond=natofirst
         npairsecond=npairfirst
         do i=1,npairsecond
            distsecond(i)=distfirst(i)
         end do
         fractsecond=fractfirst
         nzt2=nzt1
         ncz2=ncz1
         do i=1,ncz2
            lztyp2(i)=lztyp1(i)
         end do
         nyt2=nyt1
         ncy2=ncy1
         do i=1,ncy2
            lytyp2(i)=lytyp1(i)
         end do
         iref=iref+4
      endif
      !Decide whether to use simulated phim, thetam, chim
      iref=iref+1
      call parse(symbol_value(iref),nwords,ncf,ncl,5,256,nworddim)
      performtest=symbol_value(iref)(ncf(1):ncl(1)).ne.'notest'.and.nwords.gt.7
      if(performtest) then
         read(symbol_value(iref)(ncf(1):ncl(nwords)),*,iostat=ierr) phim1,thetam1,chim1,phim2,thetam2,chim2,thetallab,phillab
         phim1=phim1*piconv
         thetam1=thetam1*piconv
         chim1=chim1*piconv
         phim2=phim2*piconv
         thetam2=thetam2*piconv
         chim2=chim2*piconv
         thetallab=thetallab*piconv
         phillab=phillab*piconv
      else
         phim1=0.0
         thetam1=0.0
         chim1=0.0
         phim2=0.0
         thetam2=0.0
         chim2=0.0
         thetallab=0.0
         phillab=0.0
      end if
      if(ierr.ne.0) performtest=.false.
      iref=iref+1
      outputaveragemol=symbol_value(iref).eq.'y'

c Only perform the calculation when nsum is zero or a multiple of nmult

      nmult=5
      nsum=nsumt/nmult
      nt=nsum*nmult
      if(nt.eq.nsumt) then
c
c normalise atomic coordinates to box size
c
         call normalise_coordinates()
         zmaxd=zmax/dmax
         if(zmaxd.gt.0.5) then
            zmaxd=0.5
            zmax=zmaxd*dmax
         endif
! Set up l and n values and allocate some needed arrays
         call get_l_and_n_values()
         mz=nz
         mdist=ndist
         if(allocated(z)) deallocate(z,gofz,gofzs,gzt,gzt2,titledat)
         allocate(z(mz),gofz(mz,mdist),gofzs(mz,mdist),gzt(mz,mdist),gzt2(mz,mdist),titledat(mdist))
c set up step size and z values
         zstep=zmaxd/real(nz)
         zhalf=0.5*zstep
         do j=1,nz
            z(j)=zstep*j-zhalf
         end do
         write(6,*) 'sdf> zmaxd: ',zmax,dmax,zmaxd,zstep
c zero gofz
         call initgofzs(nz,ndist,gofz,gofzs,mz,mdist)
         call sdf(nin,fname_inp)
c
c set the z values
c
         do i=1,nz
            z(i)=z(i)*dmax
         end do
c
c axis title
c
         titler='   r'
c
c possible formats for integers in title
c
         ifmt1='i1'
         ifmt2='i2'
         ifmt3='i3'
         ifmt4='i4'
c
c set the column titles
c
         do id=1,ndist
            il1=ll1(id)
            il2=ll2(id)
            il=ll(id)
            in1=n1(id)
            in2=n2(id)
c
c set up the format for the title
c
            if(il1.lt.10) then
               il1fmt=ifmt1
            else if(il1.lt.100) then
               il1fmt=ifmt2
            else if(il1.lt.1000) then
               il1fmt=ifmt3
            endif
            if(il2.lt.10) then
               il2fmt=ifmt1
            else if(il2.lt.100) then
               il2fmt=ifmt2
            else if(il2.lt.1000) then
               il2fmt=ifmt3
            endif
            if(il.lt.10) then
               ilfmt=ifmt1
            else if(il.lt.100) then
               ilfmt=ifmt2
            else if(il.lt.1000) then
               ilfmt=ifmt3
            endif
            if(in1.lt.10) then
               in1fmt=ifmt1
            else if(in1.lt.100) then
               in1fmt=ifmt2
            else if(in1.lt.1000) then
               in1fmt=ifmt3
            endif
            if(abs(in2).lt.10) then
               if(in2.lt.0) then
                  in2fmt=ifmt2
               else
                  in2fmt=ifmt1
               end if
            else if(abs(in2).lt.100) then
               if(in2.lt.0) then
                  in2fmt=ifmt3
               else
                  in2fmt=ifmt2
               end if
            else if(abs(in2).lt.1000) then
               if(in2.lt.0) then
                  in2fmt=ifmt4
               else
                  in2fmt=ifmt3
               end if
            endif
            formatout='(1h('//il1fmt//',1h,,'//il2fmt//',1h,,'//ilfmt//',1h;,'//in1fmt//',1h,,'//in2fmt//',1h))'
c
c write the title for this distribution
c
            write(titledat(id),fmt=formatout) il1,il2,il,in1,in2
         end do
         ext='h'
         call outgofzs_sharm(fnameinp,ext,nz,ndist,nsum,z,gofz,gofzs,gzt,gzt2,titler,titledat,mz,mdist)
      endif
      if(.not.performtest) nsumt=nsumt+1
c
c write input file
c
c
c make sure all current values are stored
c
      write(symbol_value(4),*) nsumt
      write(symbol_value(5),*) ndist
      write(6,*) ndist,' ',symbol_value(5)
      call save_files_others(14)
      close(14)
      if(performtest) write(6,*) 'sdf> Test performed.'
      return
      end

      subroutine sdf(nout,fname_inp)

      use reallocation_routines
      use sharm_mod

      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      include 'sharm.inc'
      include 'sharm_sdf.inc'
      real, dimension(:), allocatable           :: xat,yat,zat
      integer, dimension(:,:), allocatable      :: iasave !Store atom numbers of atoms used to define molecules, so as not to be used in calculating coefficients
      character(len=*)                       :: fname_inp
      character(len=256)                     :: fname_xyz,textline !Filename for xyz data, line of text
      integer                                :: nout !Unit number on which to write xyz file
      real                                   :: xt,yt,zt,xr,yr,zr
      real, dimension(mcentre)               :: xfirstc,yfirstc,zfirstc,xsecondc,ysecondc,zsecondc
      character(len=3), dimension(mcentre)   :: atsymbolfirstc,atsymbolsecondc,attypefirstc,attypesecondc
      integer :: ncentral,nsecond

      pi=4.0*atan(1.0)
      piconv=pi/180.0
      degconv=180.0/pi
      zstep=zmaxd/float(nz)
c volume of sphere of radius 0.5
      vol=(alatx*(blaty*clatz))
c set up counter for second atoms
      nsec=0
! Arrays to store Euler angles
      if(allocated(nm)) deallocate(nm,phim,cthetm,chim,xcent,ycent,zcent,dlmnmc,chiphi,dlmnc,philab)
      call get_max_size_dlmn()
      write(6,*) 'sdf> ',nmol,mdlmn
      allocate(nm(nmol),phim(nmol),cthetm(nmol),chim(nmol),xcent(nmol),ycent(nmol),zcent(nmol))
      allocate(dlmnmc(2,nmol,mdlmn),chiphi(2,nmol,mdlmn),dlmnc(2,mdlmn),philab(mdlmn))
      if(performtest) then
! Allocate suitable arrays if performing a test
         if (allocated(phimt1)) deallocate(phimt1,cthetmt1,chimt1,phimt2,cthetmt2,chimt2,dlmnmct,chiphit)
         allocate(phimt1(nmol),cthetmt1(nmol),chimt1(nmol),phimt2(nmol),cthetmt2(nmol),chimt2(nmol))
         allocate(dlmnmct(2,nmol,mdlmn),chiphit(2,nmol,mdlmn))
      endif
! Allocate arrays to store atom numbers used in molecule definitions
      if(allocated(iasave)) deallocate(iasave)
      allocate(iasave(max(natofirst,natosecond),nmol))

c set up coordinates of Laboratory X, Y and Z axes
c
      alab(1)=1.0
      alab(2)=0.0
      alab(3)=0.0
      blab(1)=0.0
      blab(2)=1.0
      blab(3)=0.0
      clab(1)=0.0
      clab(2)=0.0
      clab(3)=1.0
c
c set up coordinates of molecular centres and Euler angles for each molecule 
c in laboratory coordinate frame
c
      nsmmo1=0
      nsmmo2=0
! Save the label of the first centre molecule
      imolcentre=0

c Look for centre molecules. When found calculate their Euler angles.

      imol1=0
      ncentral=0
      nsecond=0

c step through the atoms and find those which correspond to ltyp1(1)
      do ia1=1,natom

c First set up the test conditions

c Signal that so far a true molecule has not been found

         iflag=0

         itmol=0

         do while(itmol.le.ldistinct)
            itmol=itmol+1
            if(itmol.eq.1) then
               natomol=natofirst
               npairmol=npairfirst
               do i=1,natomol
                  lctyp(i)=lctyp1(i)
               end do
               do i=1,npairmol
                  distcentre(i)=distfirst(i)/dmax
               end do
               fractcentre=fractfirst
            else
               natomol=natosecond
               npairmol=npairsecond
               do i=1,natomol
                  lctyp(i)=lctyp2(i)
               end do
               do i=1,npairmol
                  distcentre(i)=distsecond(i)/dmax
               end do
               fractcentre=fractsecond
            endif

c Test to see if this atom corresponds to the required centre atom.

            if(itype(ia1).eq.lctyp(1)) then

c Increment the molecule counter in case this is a valid molecule. It will later be
c decremented if this turns out to not to be valid.

               if(imol1.lt.nmol) imol1=imol1+1

c Set the label whether this is a first or second molecule

               nm(imol1)=itmol

c Now search for other atoms which satisfy the specified label and distance criteria

               nfound=1
               atomno(1)=ia1
               im1=mref(ia1)
               xo1=xadd(im1)+xa(ia1)
               yo1=yadd(im1)+ya(ia1)
               zo1=zadd(im1)+za(ia1)
c Temporary coordinates of this atom relative to the centre of the molecule
               xtemp(1)=0.0
               ytemp(1)=0.0
               ztemp(1)=0.0
               xcent(imol1)=xo1
               ycent(imol1)=yo1
               zcent(imol1)=zo1

c At this point there if there is only 1 atom in the molecule and it must be of the
c correct type, otherwise we wouldn't have got here, we signal immediately that a 
c molecule has been found, and there is no need to check intra-molecular distances.

               if(natomol.gt.nfound) then

c Set the other atom numbers for this molecule to zero. This will be used to signal whether
c or not an atom has been found at each index. If the correct set of atom types with the
c specified distances from neighbours cannot be found, this molecule is abandoned as a 
c possible centre molecule

                  do i=2,natomol
                     atomno(i)=0
                  end do

c Step through second atoms and see which if any satisfy the type and distance constraints

                  ia2=0
                  do while (ia2.lt.natom.and.nfound.lt.natomol)

                     ia2=ia2+1

c First check whether this atom already exists for this molecule

                     itest=1
                     do while (itest.le.natomol.and.atomno(itest).ne.ia2)
                        itest=itest+1
                     end do

c itest .gt. natomol means it was not found so far

                     if(itest.gt.natomol) then

c Now see whether it fits the type number for one of the atoms of this molecule.

                        i=1
                        ifound=0
                        do while(i.lt.natomol.and.ifound.eq.0)
                           i=i+1

c Check to see if ia2 matches the current label and that this slot has not
c already been assigned

                           if(lctyp(i).eq.itype(ia2).and.atomno(i).eq.0) then

c now measure its distance from the centre atom

                              im2=mref(ia2)
                              xo2=xadd(im2)+xa(ia2)
                              yo2=yadd(im2)+ya(ia2)
                              zo2=zadd(im2)+za(ia2)
                              xdiff=xo2-xo1
                              ydiff=yo2-yo1
                              zdiff=zo2-zo1
                              call minimagec(xdiff,ydiff,zdiff)
                              dmeas=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
                              dmeas=sqrt(dmeas)
                              dtestlow=distcentre(i-1)*(1.0-fractcentre)
                              dtesthigh=distcentre(i-1)*(1.0+fractcentre)
                              if(dmeas.ge.dtestlow.and.dmeas.le.dtesthigh) then

c atom found at correct distance from centre - save the atom number in atomno.

                                 ifound=1
                                 atomno(i)=ia2
                                 xtemp(i)=xdiff
                                 ytemp(i)=ydiff
                                 ztemp(i)=zdiff
                                 nfound=nfound+1
                              endif
                           endif
                        end do
                     endif
                  end do


c If atoms of the correct number and type and at the correct distance
c from the central atom have been found, then this group is a possible candidate for
c a molecule. Now check that the other distances are correct.

                  if(nfound.eq.natomol) then
                     iflag=1
                     iref=natomol-1
                     i1=2
                     i2=3
                     do while (iref.lt.npairmol.and.iflag.eq.1)
                        xdiff=xtemp(i2)-xtemp(i1)
                        ydiff=ytemp(i2)-ytemp(i1)
                        zdiff=ytemp(i2)-ztemp(i1)
                        call minimagec(xdiff,ydiff,zdiff)
                        dmeas=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
                        dmeas=sqrt(dmeas)
                        iref=iref+1
                        dtestlow=distcentre(iref)*(1.0-fractcentre)
                        dtesthigh=distcentre(iref)*(1.0+fractcentre)
                        if(dmeas.lt.dtestlow.or.dmeas.gt.dtesthigh) iflag=0
                        i2=i2+1
                        if(i2.gt.natomol) then
                           i1=i1+1
                           i2=i1+1
                        endif
                     end do
                  else
                     iflag=0
                  endif
               else
                  iflag=1
               endif
! If iflag is still 1 after this it means a valid "molecule" has been found
! so we determine its Euler angles
               if(iflag.eq.1) then
! Save the atom numbers that form this molecule
                  do i=1,natomol
                     iasave(i,imol1)=atomno(i)
                  end do
                  if(nm(imol1).eq.1) then
                     nsmmo1=nsmmo1+1
                  endif
                     if(nm(imol1).eq.2) then
                     nsmmo2=nsmmo2+1
                  endif
! Form atomic coordinates relative to this centre
                  call reallocate1d_r(xat,natomol)
                  call reallocate1d_r(yat,natomol)
                  call reallocate1d_r(zat,natomol)
                  do ia=1,natomol
                     xat(ia)=xtemp(ia)
                     yat(ia)=ytemp(ia)
                     zat(ia)=ztemp(ia)
                  end do
                  if(nm(imol1).eq.1) then
                     call form_eulerangles(imol1,xat,yat,zat,nzt1,nyt1,ncz1,ncy1,lztyp1,lytyp1,alab,blab,clab,phim(imol1),cthetm(imol1),chim(imol1))
                  else if(nm(imol1).eq.2) then
                     call form_eulerangles(imol1,xat,yat,zat,nzt2,nyt2,ncz2,ncy2,lztyp2,lytyp2,alab,blab,clab,phim(imol1),cthetm(imol1),chim(imol1))
                  end if
!If testing orientations, supply a fictitious but known distribution of angles
                  if(performtest) then
!                     write(6,499) imol1,phim1,thetam1,chim1,phim2,thetam2,chim2
!                     call populate_disk(pi,phim1,thetam1,chim1,phinew,thetanew,chinew,thetaspread,phispread)
! Rotate this molecule onto desired orientation
                     phimt1(imol1)=phim1
                     cthetmt1(imol1)=cos(thetam1)
                     chimt1(imol1)=chim1
!                     call populate_disk(pi,phim2,thetam2,chim2,phinew,thetanew,chinew,thetaspread,phispread)
                     phimt2(imol1)=phim2
                     cthetmt2(imol1)=cos(thetam2)
                     chimt2(imol1)=chim2
!                     write(6,499) imol1,phim(imol1),cthetm(imol1),chim(imol1),phimt(imol1),cthetmt(imol1),chimt(imol1)
499   format(1x,i5,6(1x,e12.5))
                  end if
                  call form_dlmn(imol1,ierr)
                  if(ierr.ne.0) then
                     write(6,*) 'sdf> Error when calculating dlm coefficients for molecule ',imol1
                     return
                  end if
! If this is a centre molecule, we form its average rotated coordinates for later output
                  if(nm(imol1).eq.1) then
                     if(ncentral.eq.0.or.outputaveragemol) then
                        ncentral=ncentral+1
                        do ia=1,natomol
! Rotate coordinates onto lab axes
! Temporary coordinates
                           xt=xtemp(ia)
                           yt=ytemp(ia)
                           zt=ztemp(ia)
! Need to rotate these to align with the laboratory axes - performed in reverse order about the laboratory axes
                           if(natofirst.gt.1) then
                              call rotateeuler(xt,yt,zt,-chim(imol1),-acos(cthetm(imol1)),-phim(imol1),xr,yr,zr)
                           else
                              xr=xt
                              yr=yt
                              zr=zt
                           end if
                           write(6,*) atomno(ia),xr*dmax,yr*dmax,zr*dmax
                           if(ncentral.eq.1) then
                              xfirstc(ia)=xr
                              yfirstc(ia)=yr
                              zfirstc(ia)=zr
                              atsymbolfirstc(ia)=atsymbol(itype(atomno(ia)))
                              attypefirstc(ia)=type(itype(atomno(ia)))
                           else
                              xfirstc(ia)=xfirstc(ia)+xr
                              yfirstc(ia)=yfirstc(ia)+yr
                              zfirstc(ia)=zfirstc(ia)+zr
                           end if
                        end do
                     end if
                  end if
! If this is a centre molecule, we form its average rotated coordinates for later output
                  if(nm(imol1).eq.2) then
                     if(nsecond.eq.0.or.outputaveragemol) then
                        nsecond=nsecond+1
                        do ia=1,natomol
! Rotate coordinates onto lab axes
! Temporary coordinates
                           xt=xtemp(ia)
                           yt=ytemp(ia)
                            zt=ztemp(ia)
! Need to rotate these to align with the laboratory axes - performed in reverse order about the laboratory axes
                           if(natosecond.gt.1) then
                              call rotateeuler(xt,yt,zt,-chim(imol1),-acos(cthetm(imol1)),-phim(imol1),xr,yr,zr)
                           else
                              xr=xt
                              yr=yt
                              zr=zt
                           end if
!                        write(6,*) atomno(ia),xr*dmax,yr*dmax,zr*dmax
                           if(nsecond.eq.1) then
                              xsecondc(ia)=xr
                              ysecondc(ia)=yr
                              zsecondc(ia)=zr
                              atsymbolsecondc(ia)=atsymbol(itype(atomno(ia)))
                              attypesecondc(ia)=type(itype(atomno(ia)))
                           else
                              xsecondc(ia)=xsecondc(ia)+xr
                              ysecondc(ia)=ysecondc(ia)+yr
                              zsecondc(ia)=zsecondc(ia)+zr
                           end if
                        end do
                     end if
                  end if
               else
c If this is not a valid molecule then we ignore this candidate
                  if(imol1.gt.0) imol1=imol1-1
               endif
            endif
         end do
      end do
      nmolsharm=imol1
! If ldistinct = 1 (molecules are distinct) we need to ensure that none of the atoms that occur in second molecules do not also occur in first molecules
      if(ldistinct.gt.0) then
         imol1=0
         do while(imol1.lt.nmolsharm)
            imol1=imol1+1
            if(nm(imol1).eq.1) then
! Test second molecules to see if there any atoms common with first molecule
               imol2=0
               do while(imol2.lt.nmolsharm)
                  imol2=imol2+1
                  iflag=0
                  if(imol1.ne.imol2.and.nm(imol2).eq.2) then
                     ia1=0
                     do while(ia1.lt.natofirst.and.iflag.eq.0) 
                        ia1=ia1+1
                        ia2=0
                        do while(ia2.lt.natosecond.and.iflag.eq.0)
                           ia2=ia2+1
                           if(iasave(ia1,imol1).eq.iasave(ia2,imol2)) iflag=1
                        end do
                     end do
                  end if
                  if(iflag.eq.1) then !This second molecule must be removed from list
                     write(6,'(a,i5)') 'sdf> Removing molecule: ',imol2
                     imol=imol2
                     do while (imol.lt.nmolsharm)
                        imolt=imol+1
                        nm(imol)=nm(imolt)
                        phim(imol)=phim(imolt)
                        cthetm(imol)=cthetm(imolt)
                        chim(imol)=chim(imolt)
                        if(nm(imol).eq.1) then
                           do ia=1,natofirst
                              iasave(ia,imol)=iasave(ia,imolt)
                           end do
                        else
                           do ia=1,natosecond
                              iasave(ia,imol)=iasave(ia,imolt)
                           end do
                        end if
                        imol=imolt
                     end do
                     nmolsharm=nmolsharm-1
                     nsmmo2=nsmmo2-1
                     if(imol1.ge.imol2) imol1=imol1-1 !Need to decrement imol1 if it is .gt. imol2, since imol2 has been removed from the list
                     imol2=imol2-1
                  end if
               end do
            end if
         end do
         !Write out lists of atoms to check that algorithm works
         do imol=1,nmolsharm
            if(nm(imol).eq.1) write(6,*) imol,nm(imol),(iasave(ia,imol),ia=1,natofirst)
         end do
         do imol=1,nmolsharm
            if(nm(imol).eq.2) write(6,*) imol,nm(imol),(iasave(ia,imol),ia=1,natosecond)
         end do
      end if
      write(6,333) nsmmo1,nsmmo2
333   format(1x,'sdf> Number of type 1 molecules = ',i10
     1/1x, 'sdf> Number of type 2 molecules = ',i10)

! Write out coordinates of central molecule

      fname_xyz=fname_inp(1:indexlastdot(fname_inp)-1)//'-first.xyz'
      open(nout,file=fname_xyz,status='unknown',iostat=ierr)
      if(ierr.eq.0.and.ncentral.gt.0) then
! Write number of atoms in molecule
         write(nout,'(i6)') natofirst
! Write information line - not used for plotting
         write(nout,'(a)') 'Central molecule coordinate file: '//fname_xyz(1:len_trim(fname_xyz))
! Shift coordinates to centre used by SHARM/SDF
         do ia=1,natofirst
            xr=xfirstc(ia)*dmax/real(ncentral)
            yr=yfirstc(ia)*dmax/real(ncentral)
            zr=zfirstc(ia)*dmax/real(ncentral)
            write(textline,'(a3,1x,f10.5,1x,f10.5,1x,f10.5,1x,a3)') atsymbolfirstc(ia),xr,yr,zr,attypefirstc(ia)
            textline=textline(1:len_trim(textline))//'  #CENTRALMOLECULE'
            write(nout,'(a)') textline(1:len_trim(textline))
         end do
         close(nout)
      else
         write(6,*) 'sdf> Error writing first molecule file'
         stop
      end if
! Write out coordinates of second molecule
      fname_xyz=fname_inp(1:indexlastdot(fname_inp)-1)//'-second.xyz'
      open(nout,file=fname_xyz,status='unknown',iostat=ierr)
      if(ierr.eq.0.and.ldistinct.eq.0.and.ncentral.gt.0) then
! Write number of atoms in molecule
         write(nout,'(i6)') natofirst
! Write information line - not used for plotting
         write(nout,'(a)') 'Central molecule coordinate file: '//fname_xyz(1:len_trim(fname_xyz))
! Shift coordinates to centre used by SHARM/SDF
         do ia=1,natofirst
            xr=xfirstc(ia)*dmax/real(ncentral)
            yr=yfirstc(ia)*dmax/real(ncentral)
            zr=zfirstc(ia)*dmax/real(ncentral)
            write(textline,'(a3,1x,f10.5,1x,f10.5,1x,f10.5,1x,a3)') atsymbolfirstc(ia),xr,yr,zr,attypefirstc(ia)
            textline=textline(1:len_trim(textline))//'  #CENTRALMOLECULE'
            write(nout,'(a)') textline(1:len_trim(textline))
         end do
         close(nout)
      else if(ierr.eq.0.and.ldistinct.eq.1.and.nsecond.gt.0) then !Second molecule is different from central molecule
! Write number of atoms in molecule
         write(nout,'(i6)') natosecond
! Write information line - not used for plotting
         write(nout,'(a)') 'Central molecule coordinate file: '//fname_xyz(1:len_trim(fname_xyz))
! Shift coordinates to centre used by SHARM/SDF
         do ia=1,natosecond
            xr=xsecondc(ia)*dmax/real(nsecond)
            yr=ysecondc(ia)*dmax/real(nsecond)
            zr=zsecondc(ia)*dmax/real(nsecond)
            write(textline,'(a3,1x,f10.5,1x,f10.5,1x,f10.5,1x,a3)') atsymbolsecondc(ia),xr,yr,zr,attypesecondc(ia)
            textline=textline(1:len_trim(textline))//'  #CENTRALMOLECULE'
            write(nout,'(a)') textline(1:len_trim(textline))
         end do
         close(nout)
      else
         write(6,*) 'sdf> Error writing second molecule file'
         stop
      end if
      call form_sum_prod_dlmn_cgcoeff(nmolsharm,ierr)
      return
      end

