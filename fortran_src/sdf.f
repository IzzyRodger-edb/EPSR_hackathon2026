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
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lmncom.inc'
      include 'sharm.inc'
      character*50 formatout
      character*2 ifmt1,ifmt2,ifmt3,il1fmt,il2fmt,ilfmt,in1fmt,in2fmt

      call fctrln(250)
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

c
c read the input file
c
      call read_input_file(nin,nout,ierrout)
c
c if an error occurred on read, exit now
c
      if(ierrout.ne.0) return
c
c set up number of radius values
c
      read(symbol_value(2),*) nz
c
c set up rmax of calculation
c
      read(symbol_value(3),*) zmax
c
c set up number of configurations already stored
c
      read(symbol_value(4),*) nsumt
c
c set up number of coefficients
c
      read(symbol_value(5),*) ndist
c
c set up l1 values
c
      call parse(symbol_value(6),nwords,ncf,ncl,1,256,nworddim)
      nll1v=nwords
      do i=1,nll1v
         read(symbol_value(6)(ncf(i):ncl(i)),*,iostat=ierr) ll1v(i)
      end do
c
c l2 values other than 0 are not used in SDF
c
      call parse(symbol_value(7),nwords,ncf,ncl,1,256,nworddim)
      nll2v=nwords
      do i=1,nll2v
         read(symbol_value(7)(ncf(i):ncl(i)),*,iostat=ierr) ll2v(i)
      end do
c
c set up l values
c For SDF these are set identical to l1 in the setup program since l2==0
      call parse(symbol_value(8),nwords,ncf,ncl,5,256,nworddim)
      nllv=nwords
      do i=1,nllv
         read(symbol_value(8)(ncf(i):ncl(i)),*) llv(i)
      end do
c
c set up n1step
c
      read(symbol_value(9),*) n1step
c
c set up n2step
c 
      read(symbol_value(10),*) n2step
c
c For the central "molecule" read the list of atom types to be used to define the 
c unit. The first in the list will be assigned to the centre of the molecule.
c
      call parse(symbol_value(11),nwords,ncf,ncl,5,256,nworddim)
      natofirst=0

c skip through the atom types to find out the type number of the specified atom

      ic=0
      do while (ic.lt.nwords.and.natofirst.lt.mcentre)
         ic=ic+1
         it=0
         ifound=0
         do while(it.lt.ntype.and.ifound.eq.0)
            it=it+1
            if(type(it).eq.symbol_value(11)(ncf(ic):ncl(ic))) ifound=it
         end do
         if(ifound.le.0) then
            write(6,201) symbol_value(11)(ncf(ic):ncl(ic))
     *,fname_ato(1:ilenfname_ato)
201   format(/'sdf> Atom type ',a,' does not exist in .ato file ',a)
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

      call parse(symbol_value(12),nwords,ncf,ncl,5,256,nworddim)
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
         read(symbol_value(12)(ncf(i):ncl(i)),*) distfirst(i)
      end do
      npairfirst=i

c Variation of distance (half width) about the specified distance for an atom to be accepted.
c This is expressed as a fraction of the specified distance. For example fraction of 0.2
c for a bond distance of 2A means atoms within 2 +- 0.4A will be accepted

      read(symbol_value(13),*) fractfirst

c
c for molecule 1 read primary axis type (1=X, 2=Y, 3=Z), and
c the number and atom numbers of atoms used to define the molecule's nzt1
c axis
c
      call parse(symbol_value(14),nwords,ncf,ncl,5,256,nworddim)
      nzt1=symbol_value(14)(ncf(1):ncl(1))
      ncz1=nwords-1
      if(ncz1.gt.mcentre) ncz1=mcentre
      do i=1,ncz1
         iref=i+1
         read(symbol_value(14)(ncf(iref):ncl(iref)),*) lztyp1(i)
      end do
c
c for molecule 1 read secondary axis type (1=X, 2=Y, 3=Z), and
c the number and atom numbers of atoms used to define the
c molecule's nyt1 axis
c
      call parse(symbol_value(15),nwords,ncf,ncl,5,256,nworddim)
      nyt1=symbol_value(15)(ncf(1):ncl(1))
      ncy1=nwords-1
      if(ncy1.gt.mcentre) ncy1=mcentre
      do i=1,ncy1
         iref=i+1
         read(symbol_value(15)(ncf(iref):ncl(iref)),*) lytyp1(i)
      end do

c For second "molecule" read type of atom used to define the molecular
c centre. Only the first in the list is used.

      call parse(symbol_value(16),nwords,ncf,ncl,5,256,nworddim)

c Skip through the specified atom types to find out the type number of the specified atom
c We also test to see if this atom is the same as the corresponding atom on the central
c molecule. If not the molecules are treated as distinct.

      natosecond=0
      ldsnct=0
      do ic=1,nwords
         it=0
         ifound=0
         do while(it.lt.ntype.and.ifound.eq.0
     *.and.natosecond.lt.mcentre)
            it=it+1
            if(type(it).eq.symbol_value(16)(ncf(ic):ncl(ic))) ifound=it
         end do
         if(ifound.eq.0) then
            write(6,201) symbol_value(16)(ncf(ic):ncl(ic))
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

      ldsnct=0
      if(natosecond.eq.natofirst) then
         ldsnct=0
         do it=1,natofirst
            if(lctyp2(it).ne.lctyp1(it)) ldsnct=1
         end do
      else
         ldsnct=1
      endif

c If the second and central molecules are identical, the rest of the input for the
c second molecule is identical to the first.

      if(ldsnct.ne.0) then

c For second molecule read the array of maximum distances between atoms in the order
c (1,2) (1,3) (1,4) ... (2,3) (2,4) ... (3,4) etc.
c There should be natosecond*(natosecond-1)/2 of them

         call parse(symbol_value(17),nwords,ncf,ncl,5,256,nworddim)
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
            read(symbol_value(17)(ncf(i):ncl(i)),*) distsecond(i)
         end do
         npairsecond=i

c Variation of distance (half width) about the specified distance for an atom to be accepted.
c This is expressed as a fraction of the specified distance. For example fraction of 0.2
c for a bond distance of 2A means atoms within 2 +- 0.4A will be accepted

         read(symbol_value(18),*) fractsecond

c for molecule 2 read primary axis type (1=X, 2=Y, 3=Z), and
c the number and atom numbers of atoms used to define the molecule's nzt1
c axis

         call parse(symbol_value(19),nwords,ncf,ncl,5,256,nworddim)
         nzt2=symbol_value(19)(ncf(1):ncl(1))
         ncz2=nwords-1
         do i=1,ncz2
            iref=i+1
            read(symbol_value(19)(ncf(iref):ncl(iref)),*) lztyp2(i)
         end do

c for molecule 2 read secondary axis type (1=X, 2=Y, 3=Z), and
c the number and atom numbers of atoms used to define the
c molecule's nyt1 axis

         call parse(symbol_value(20),nwords,ncf,ncl,5,256,nworddim)
         nyt2=symbol_value(20)(ncf(1):ncl(1))
         ncy2=nwords-1
         do i=1,ncy2
            iref=i+1
            read(symbol_value(20)(ncf(iref):ncl(iref)),*) lytyp2(i)
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

      endif

c Only perform the calculation when nsum is zero or a multiple of nmult

      nmult=5
      nsum=nsumt/nmult
      nt=nsum*nmult
      if(nt.eq.nsumt) then
c 
c set up l and n values
c
            nterm=0
c l1 values
            do ill1=1,nll1v
              il1=ll1v(ill1)
c l2 values
              do ill2=1,nll2v
                il2=ll2v(ill2)
                if((ldsnct.eq.0.and.il2.le.il1).or.ldsnct.eq.1) then
c l values - use closure relation on CG coefficients to determine range
                jlmax=iabs(il1+il2)
                jlmin=iabs(il1-il2)
                do ill=1,nllv
                  il=llv(ill)
                  if(il.ge.jlmin.and.il.le.jlmax) then
c n1 values
                  if(n1step.gt.0) then
                  n1max=il1+1
                    n1s=n1step
                  else
                  n1max=1
                  n1s=1
                  endif
                  do jn1=1,n1max,n1s
                    in1=jn1-1
c n2 values
                    if(n2step.gt.0) then
                    n2max=il2+1
                    n2s=n2step
                    else
                    n2max=1
                    n2s=1
                    endif
                  do jn2=1,n2max,n2s
                    in2=jn2-1
c test that this is allowed
                    ltest=il1+il2+il
                    ltest2=2*(ltest/2)
                    if(ltest.eq.ltest2.or.in1.ne.0.or.in2.ne.0) then
c save them!
                      if(nterm.lt.nlmax) then
                        nterm=nterm+1
                        ll1(nterm)=il1
                        ll2(nterm)=il2
                        ll(nterm)=il
                        n1(nterm)=in1
                        n2(nterm)=in2
                      endif
                    endif
                  end do
                  end do
                  endif
                end do
                endif
              end do
            end do
            ndist=nterm
            do i=1,ndist
                  write(6,46) i,ll1(i),ll2(i),ll(i),n1(i),n2(i)
            end do
46    format(6(1x,i5))
c
c normalise atomic coordinates to box size
c
            call normalise_coordinates()
            zmaxd=zmax/dmax
            if(zmaxd.gt.0.5) then
                  zmaxd=0.5
                  zmax=zmaxd*dmax
            endif
c
c set up step size and z values
c
c
c zero gofz
c
            zstep=zmaxd/float(nz)
            zhalf=0.5*zstep
            do j=1,nz
                  z(j)=zstep*j-zhalf
            end do
            call initgofz(nz,ndist,gofz,mz,mdist)

            call sdf(ndist,nin,fname_inp)
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
                  if(in2.lt.10) then
                     in2fmt=ifmt1
                  else if(in2.lt.100) then
                     in2fmt=ifmt2
                  else if(in2.lt.1000) then
                     in2fmt=ifmt3
                  endif
                  formatout='(1h('//il1fmt//',1h,,'//il2fmt//',1h,,'
     *//ilfmt//',1h;,'//in1fmt//',1h,,'//in2fmt//',1h))'
c
c write the title for this distribution
c
                  write(titledat(id),fmt=formatout) il1,il2,il,in1,in2
            end do

            ext='h'

            call outgofz_sharm(fnameinp,ext,nz,ndist,nsum,z,gofz,gzt,gzt2
     *,titler,titledat,mz,mdist,ll1,ll2,ll,n1,n2)

      endif

      nsumt=nsumt+1

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

      return
      end

      subroutine sdf(ndist,nout,fname_inp)
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lmncom.inc'
      include 'sharm.inc'
      include 'sharm_sdf.inc'

		character(len=*)                       :: fname_inp
		character(len=256)                     :: fname_xyz,textline !Filename for xyz data, line of text
		integer                                :: nout !Unit number on which to write xyz file
      real                                   :: xt,yt,zt,xr,yr,zr
      real, dimension(mcentre)               :: xfirstc,yfirstc,zfirstc,xsecondc,ysecondc,zsecondc
      character(len=3), dimension(mcentre)   :: atsymbolfirstc,atsymbolsecondc
      integer :: ncentral,nsecond

c     integer*4 rate,max,time1,time2,time3,time4
      pi=3.141592653
      zstep=zmaxd/float(nz)
c volume of sphere of radius 0.5
      gfac=4.0*pi/3.0
      vol=1.0
c set up counter for second atoms
      nsec=0
      do id=1,ndist
            alsec(id)=0.0
      end do
c
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
c set up the required CG coefficients
c
      ic=0
      do id=1,ndist
        il1=ll1(id)
        in1=n1(id)
        il2=ll2(id)
        in2=n2(id)
        il=ll(id)
c
c initialize reference integer for molecule 1
c
        ilmn1=1+(il1*(12+(il1-1)*(7+4*il1)))/6+in1*(2*il1+1)
        il1rf(id)=ilmn1
c
c initialize reference integer for molecule 1
c
        ilmn2=1+(il2*(12+(il2-1)*(7+4*il2)))/6+in2*(2*il2+1)
        il2rf(id)=ilmn2
c
c initialize reference integer for laboratory function
c
        ilmn=il*il+il+1
        ilrf(id)=ilmn
        nm1=2.0*il1+1
        m1min=-il1-1
        do mm1=1,nm1
          m1=mm1+m1min
          m2max=min0(il2,il-m1)
          m2min=max0(-il2,-il-m1)
          nm2=m2max-m2min+1
          m2min=m2min-1
          do mm2=1,nm2
            m2=mm2+m2min
            m=m1+m2
501   format(1x,8(1x,i5))
            ic=ic+1
c error exit if too many coefficients demanded
            if(ic.le.mcg) then
                  termcg(ic)=cgcoef(il1,il2,il,m1,m2)
                  alsec(id)=alsec(id)+termcg(ic)*termcg(ic)
            endif
          end do
        end do
      end do
      if(ic.le.mcg) then
        ncg=ic
        write(6,989) ncg
989   format(1x,'No. of CG coefficients = ',i10)
      else
        go to 999
      endif
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

        do while(itmol.le.ldsnct)
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

            if(imol1.lt.mmol) imol1=imol1+1

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
               if(nm(imol1).eq.1) then
                  nsmmo1=nsmmo1+1
                  naxis1=nzt1
                  naxis2=nyt1
               endif
               if(nm(imol1).eq.2) then
                  nsmmo2=nsmmo2+1
                  naxis1=nzt2
                  naxis2=nyt2
               endif

c find coordinates of arithmetic centre of lztyp1 or lztyp2 atoms and lytyp1 
c or lytyp2 atoms (relative to origin)

               ifindz=0
               ifindy=0
               czx=0.0
               czy=0.0
               czz=0.0
               cyx=0.0
               cyy=0.0
               cyz=0.0
               if(nm(imol1).eq.1) then
                  do iz=1,ncz1
                     ia=lztyp1(iz)
                     ifindz=ifindz+1
                     czx=czx+xtemp(ia)
                     czy=czy+ytemp(ia)
                     czz=czz+ztemp(ia)
                  end do
                  do iy=1,iabs(ncy1)
                     ia=lytyp1(iy)
                     ifindy=ifindy+1
                     cyx=cyx+xtemp(ia)
                     cyy=cyy+ytemp(ia)
                     cyz=cyz+ztemp(ia)
                  end do
               else
                  do iz=1,ncz2
                     ia=lztyp2(iz)
                     ifindz=ifindz+1
                     czx=czx+xtemp(ia)
                     czy=czy+ytemp(ia)
                     czz=czz+ztemp(ia)
                  end do
                  do iy=1,ncy2
                     ia=lytyp2(iy)
                     ifindy=ifindy+1
                     cyx=cyx+xa(ia)
                     cyy=cyy+ya(ia)
                     cyz=cyz+za(ia)
                  end do
               endif
               if(ifindz.gt.0) then
                     ct(1)=czx/ifindz
                     ct(2)=czy/ifindz
                     ct(3)=czz/ifindz
               else
                     ct(1)=0.0
                     ct(2)=0.0
                     ct(3)=1.0
               endif
               if(ifindy.gt.0) then
                     pt(1)=cyx/ifindy
                     pt(2)=cyy/ifindy
                     pt(3)=cyz/ifindy
               else
                     pt(1)=0.0
                     pt(2)=0.0
                     pt(3)=0.0
               endif
c     write(6,991) p1(1),p1(2),p1(3)
c     write(6,991) c1(1),c1(2),c1(3)

c renormalize these coordinates to unit vectors
c
               oct=sqrt(scapro(ct,ct))
c
c if oct is zero it means molecule is spherically symmetric, so no thetam can
c be calculated for this molecule
c
               if(oct.gt.0.0) then
                     do i=1,3
                        ct(i)=ct(i)/oct
                     end do
                     oct=1.0
               else
                     oct=0.0
               endif
               opt=sqrt(scapro(pt,pt))
c
c if opt is zero it means molecule is linear so no chi dependence can be 
c calculated.
c
               if(opt.gt.0.0) then
                     do i=1,3
                        pt(i)=pt(i)/opt
                     end do
                     opt=1.0
c
c vector product of pt with ct (z1) gives the at axis
c
                     call vecpro(pt,ct,at)
                     oat=sqrt(scapro(at,at))
c
c renormalize to give unit vector
c
                     do i=1,3
                        at(i)=at(i)/oat
                     end do
                     oat=1.0
c
c vector product of ct to at gives  bt
c
                     call vecpro(ct,at,bt)
                     obt=sqrt(scapro(bt,bt))
c
c renormalize to give unit vector
c
                     do i=1,3
                        bt(i)=bt(i)/obt
                     end do
                     obt=1.0
               else
                     oat=0.0
                     obt=0.0
               endif
               oa1=0.0
               ob1=0.0
               oc1=0.0
c
c put axis coordinates into correct arrays, depending on values of naxis1,naxis2
c
               if(naxis1.eq.'z'.and.naxis2.eq.'y') then
                     do i=1,3
                        a1(i)=at(i)
                        b1(i)=bt(i)
                        c1(i)=ct(i)  
                     end do
                     oa1=oat
                     ob1=obt
                     oc1=oct
               else if(naxis1.eq.'-z'.and.naxis2.eq.'y') then
                     do i=1,3
                        a1(i)=at(i)
                        b1(i)=bt(i)
                        c1(i)=-ct(i)       
                     end do
                     oa1=oat
                     ob1=obt
                     oc1=oct
               else if(naxis1.eq.'y'.and.naxis2.eq.'z') then
                     do i=1,3
                        a1(i)=-at(i)
                        b1(i)=ct(i)
                        c1(i)=bt(i)  
                     end do
                     oa1=oat
                     ob1=oct
                     oc1=obt
               else if(naxis1.eq.'y'.and.naxis2.eq.'-z') then
                     do i=1,3
                        a1(i)=-at(i)
                        b1(i)=ct(i)
                        c1(i)=-bt(i)       
                     end do
                     oa1=oat
                     ob1=oct
                     oc1=obt
               else if(naxis1.eq.'y'.and.naxis2.eq.'x') then
                     do i=1,3
                        a1(i)=bt(i)
                        b1(i)=ct(i)
                        c1(i)=at(i)  
                     end do
                     oa1=obt
                     ob1=oct
                     oc1=oat
               else if(naxis1.eq.'x'.and.naxis2.eq.'y') then
                     do i=1,3
                        a1(i)=ct(i)
                        b1(i)=bt(i)
                        c1(i)=-at(i)       
                     end do
                     oa1=oct
                     ob1=obt
                     oc1=oat
               else if(naxis1.eq.'z'.and.naxis2.eq.'x') then
                     do i=1,3
                        a1(i)=bt(i)
                        b1(i)=-at(i)
                        c1(i)=ct(i)  
                     end do
                     oa1=obt
                     ob1=oat
                     oc1=oct
               else if(naxis1.eq.'-z'.and.naxis2.eq.'x') then
                     do i=1,3
                        a1(i)=bt(i)
                        b1(i)=-at(i)
                        c1(i)=-ct(i)       
                     end do
                     oa1=obt
                     ob1=oat
                     oc1=oct
               else if(naxis1.eq.'x'.and.naxis2.eq.'z') then
                     do i=1,3
                        a1(i)=ct(i)
                        b1(i)=at(i)
                        c1(i)=bt(i)  
                     end do
                     oa1=oct
                     ob1=oat
                     oc1=obt
               else if(naxis1.eq.'x'.and.naxis2.eq.'-z') then
                     do i=1,3
                        a1(i)=ct(i)
                        b1(i)=at(i)
                        c1(i)=-bt(i)       
                     end do
                     oa1=oct
                     ob1=oat
                     oc1=obt
               endif
c
c now evaluate theta value for this molecule
c
               if(oc1.gt.0.0) then
                     cthetm(imol1)=scapro(clab,c1)
c
c vector product of clab to c1 gives the y' axis about which rotation thetam was
c made
c
                     call vecpro(clab,c1,perp1)
c
c renormalize to give a unit vector
c
                     operp1=sqrt(scapro(perp1,perp1))
                     if(operp1.gt.0.0) then
                        do i=1,3
                           perp1(i)=perp1(i)/operp1
                        end do
                        operp1=1.0
                        sphim=-scapro(perp1,alab)
                        cphim=scapro(perp1,blab)
                        phim(imol1)=atan2(sphim,cphim)
                     else
                        phim(imol1)=0.0
                     endif
c now calculate chim
c
c scalar product of perp1 with a1 and b1 give chim rotation of molecule 
c about molecules c-axis 
c
                     if(oa1.gt.0.0.and.ob1.gt.0.0) then
                        schim=scapro(perp1,a1)
                        cchim=scapro(perp1,b1)
                        chim(imol1)=atan2(schim,cchim)
                     else
                        chim(imol1)=0.0
                     endif
               else
                     cthetm(imol1)=0.0
                     phim(imol1)=0.0
                     chim(imol1)=0.0
               endif
! If this is a centre molecule, we form its average rotated coordinates for later output
               if(nm(imol1).eq.1) then
                  ncentral=ncentral+1
				      do ia=1,natomol
! Rotate coordinates onto lab axes
! Temporary coordinates
                     xt=xtemp(ia)
				         yt=ytemp(ia)
 				         zt=ztemp(ia)
! Need to rotate these to align with the laboratory axes - performed in reverse order about the laboratory axes
                     call rotateeuler(xt,yt,zt,-chim(imol1),-acos(cthetm(imol1)),-phim(imol1),xr,yr,zr)
      write(6,*) atomno(ia),xr*dmax,yr*dmax,zr*dmax
                     if(ncentral.eq.1) then
                        xfirstc(ia)=xr
                        yfirstc(ia)=yr
                        zfirstc(ia)=zr
                        atsymbolfirstc(ia)=atsymbol(itype(atomno(ia)))
                     else
                        xfirstc(ia)=xfirstc(ia)+xr
                        yfirstc(ia)=yfirstc(ia)+yr
                        zfirstc(ia)=zfirstc(ia)+zr
                     end if
                 end do
               end if
! If this is a centre molecule, we form its average rotated coordinates for later output
               if(nm(imol1).eq.2) then
                  nsecond=nsecond+1
				      do ia=1,natomol
! Rotate coordinates onto lab axes
! Temporary coordinates
                     xt=xtemp(ia)
				         yt=ytemp(ia)
 				         zt=ztemp(ia)
! Need to rotate these to align with the laboratory axes - performed in reverse order about the laboratory axes
                     call rotateeuler(xt,yt,zt,-chim(imol1),-acos(cthetm(imol1)),-phim(imol1),xr,yr,zr)
      write(6,*) atomno(ia),xr*dmax,yr*dmax,zr*dmax
                     if(nsecond.eq.1) then
                        xsecondc(ia)=xr
                        ysecondc(ia)=yr
                        zsecondc(ia)=zr
                        atsymbolsecondc(ia)=atsymbol(itype(atomno(ia)))
                     else
                        xsecondc(ia)=xsecondc(ia)+xr
                        ysecondc(ia)=ysecondc(ia)+yr
                        zsecondc(ia)=zsecondc(ia)+zr
                     end if
                 end do
               end if
c
c calculate the D coefficients for this molecule
c
               icount=0
               if(nm(imol1).eq.1) then
                     nllt=nll1v
                     do i=1,nll1v
                        lltt(i)=ll1v(i)
                     end do
                     ntstep=n1step
               else
                     nllt=nll2v
                     do i=1,nll2v
                        lltt(i)=ll2v(i)
                     end do
                     ntstep=n2step
               endif
c
c l values
c
               do ilt=1,nllt
                     il=lltt(ilt)
c
c initialize reference integer
c
                     ilmnl=(il*(6+(il-1)*(7+4*il)))/6
c
c range of m values
c
                     mmin=-il-1
                     mmax=2*il+1
c
c n values
c
                     if(ntstep.gt.0) then
                        ntmax=il+1
                        nts=ntstep
                     else
                        ntmax=1
                        nts=1
                     endif
                     do jn=1,ntmax,nts
                        in=jn-1
                        ilmnln=ilmnl+in*(2*il+1)
                        do im=1,mmax
                           m=im+mmin
                           ilmnlm=ilmnln+im
c     write(6,501) il,m,in,ilmnlm
c error exit if too many coeficients demanded
                           if(ilmnlm.gt.mdlmn) go to 999
                           chiphi(1,imol1,ilmnlm)=m*phim(imol1)+in*chim(imol1)
                           termd=dlm(il,m,in,cthetm(imol1))
                           dlmnmc(1,imol1,ilmnlm)=termd
                           if(in.ne.0.and.(nm(imol1).eq.2.or.(nm(imol1).eq.1.and.ldsnct.eq.0))) then
                              chiphi(2,imol1,ilmnlm)=m*phim(imol1)-in*chim(imol1)
                              dlmnmc(2,imol1,ilmnlm)=dlm(il,m,-in,cthetm(imol1))
                           end if
c                          alsec(ilmnlm)=alsec(ilmnlm)+termd*termd
                        end do
                     end do
               end do
!               write(6,499) imol1,nm(imol1),ia1,
!     *phim(imol1),cthetm(imol1),chim(imol1)
c     *,dlmnmc(imol1,1)
499   format(1x,3(i5,1x),4(1x,e12.5))
            else

c If this is not a valid molecule then we ignore this candidate

               if(imol1.gt.0) imol1=imol1-1

            endif
          endif
        end do
      end do
      nmolsharm=imol1
      write(6,333) nsmmo1,nsmmo2
333   format(1x,'Number of type 1 molecules = ',i10
     1/1x, 'Number of type 2 molecules = ',i10)

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
				write(textline,'(a3,1x,f10.5,1x,f10.5,1x,f10.5)') atsymbolfirstc(ia),xr,yr,zr
				textline=textline(1:len_trim(textline))//'  #CENTRALMOLECULE'
				write(nout,'(a)') textline(1:len_trim(textline))
			end do
			close(nout)
		end if
! Write out coordinates of second molecule
      fname_xyz=fname_inp(1:indexlastdot(fname_inp)-1)//'-second.xyz'
      open(nout,file=fname_xyz,status='unknown',iostat=ierr)
		if(ierr.eq.0.and.ldsnct.eq.0.and.ncentral.gt.0) then
! Write number of atoms in molecule
			write(nout,'(i6)') natofirst
! Write information line - not used for plotting
			write(nout,'(a)') 'Central molecule coordinate file: '//fname_xyz(1:len_trim(fname_xyz))
! Shift coordinates to centre used by SHARM/SDF
			do ia=1,natofirst
            xr=xfirstc(ia)*dmax/real(ncentral)
            yr=yfirstc(ia)*dmax/real(ncentral)
            zr=zfirstc(ia)*dmax/real(ncentral)
				write(textline,'(a3,1x,f10.5,1x,f10.5,1x,f10.5)') atsymbolfirstc(ia),xr,yr,zr
				textline=textline(1:len_trim(textline))//'  #CENTRALMOLECULE'
				write(nout,'(a)') textline(1:len_trim(textline))
			end do
			close(nout)
      else if(ierr.eq.0.and.ldscnt.eq.1.and.nsecond.gt.0) then !Second molecule is different from central molecule
! Write number of atoms in molecule
			write(nout,'(i6)') natosecond
! Write information line - not used for plotting
			write(nout,'(a)') 'Central molecule coordinate file: '//fname_xyz(1:len_trim(fname_xyz))
! Shift coordinates to centre used by SHARM/SDF
			do ia=1,natofirst
            xr=xsecondc(ia)*dmax/real(nsecond)
            yr=ysecondc(ia)*dmax/real(nsecond)
            zr=zsecondc(ia)*dmax/real(nsecond)
				write(textline,'(a3,1x,f10.5,1x,f10.5,1x,f10.5)') atsymbolsecondc(ia),xr,yr,zr
				textline=textline(1:len_trim(textline))//'  #CENTRALMOLECULE'
				write(nout,'(a)') textline(1:len_trim(textline))
			end do
			close(nout)
		end if
c
c initialize counters for number of molecules of each kind in the system
c
      if(ldsnct.eq.0) then
            if(nm(nmolsharm).eq.1.or.nm(nmolsharm).eq.2) then
                  lsum1=1
                  lsum2=1
            else
                  lsum1=0
                  lsum2=0
            endif
      else 
            lsum1=0
            lsum2=0
            if(nm(nmolsharm).eq.1) lsum1=1
            if(nm(nmolsharm).eq.2) lsum2=1
      endif
      write(6,501) nmolsharm,nm(nmolsharm),lsum1,lsum2
      do imol1=1,nmolsharm-1
        if(nm(imol1).ne.0) then
!            write(6,*) imol1
c     call system_clock(time3,rate,max)
c
c increment counters for number of molecules of types 1 and 2 in the system
c
            if(ldsnct.eq.0) then
                  if(nm(imol1).eq.1.or.nm(imol1).eq.2) then
                        lsum1=lsum1+1
                        lsum2=lsum2+1
                  endif
            else 
                  if(nm(imol1).eq.1) lsum1=lsum1+1
                  if(nm(imol1).eq.2) lsum2=lsum2+1
            endif
            xc1=xcent(imol1)
            yc1=ycent(imol1)
            zc1=zcent(imol1)
c
c calc now look at second molecules
c
            nm2f=imol1+1
            do imol2=nm2f,nmolsharm
              if(imol1.ne.imol2.and.nm(imol2).ne.0) then
              iflag1=0
              iflag2=0
              nm1=nm(imol1)
              nm2=nm(imol2)
              if(ldsnct.eq.0.and.nm1.eq.1.and.nm2.eq.1) then
                iflag1=1
                iflag2=1
              endif
              if(ldsnct.eq.1.and.nm1.eq.1.and.nm2.eq.2) then
                iflag1=1
                iflag2=0
              endif
              if(ldsnct.eq.1.and.nm1.eq.2.and.nm2.eq.1) then
                iflag1=0
                iflag2=1
              endif
              if((iflag1+iflag2).ne.0) then
                xc2=xcent(imol2)
                yc2=ycent(imol2)
                zc2=zcent(imol2)
c apply periodic boundary conditions
                xdif=xc2-xc1
                ydif=yc2-yc1
                zdif=zc2-zc1
c apply periodic boundary conditions
                call minimagec(xdif,ydif,zdif)
c                xc2=xc2+anint(xc1-xc2)
c                yc2=yc2+anint(yc1-yc2)
c                zc2=zc2+anint(zc1-zc2)
c calculate separation of molecular centres
c                pdif(1)=xc2-xc1
c                pdif(2)=yc2-yc1
c                pdif(3)=zc2-zc1
                pdif(1)=xdif
                pdif(2)=ydif
                pdif(3)=zdif
c renormalize to unit vector
                opdif=sqrt(scapro(pdif,pdif))
                d12=opdif
                do i=1,3
                  pdif(i)=pdif(i)/opdif
                enddo
                opdif=1.0
                iref=int(d12/zstep)+1
c check that distance is within bounds prescribed for this distribution
                if(iref.le.nz) then
c apply cosine rule to determine cos(theta) for intermolecular axis
                  cz=scapro(clab,pdif)
c vector perpendicular to c1 and pdif - gives the rotated y-axis
                  call vecpro(clab,pdif,perp1)
c normalize to give a unit vector
                  operp1=sqrt(scapro(perp1,perp1))
                  if(operp1.gt.0.0) then
                    do i=1,3
                      perp1(i)=perp1(i)/operp1
                    end do
                    operp1=1.0
c scalar product with a1 to give -sin(phi)
                    sphil=-scapro(alab,perp1)
c scalar product with b1 to give cos(phi)
                    cphil=scapro(blab,perp1)
c angle phi in molecule 1's coordinate system:
                    phil=atan2(sphil,cphil)
                  else
                    phil=0.0
                  endif
c
c calculate coefficients for laboratory coordinates
c
c
c l values
c
                  do ilt=1,nllv
                    il=llv(ilt)
c
c initialize reference integer
c
                    ilref=il*il
c
c range of m values
c
                    mmin=-il-1
                    mmax=2*il+1
                    do im=1,mmax
                      m=im+mmin
                      ilm=ilref+im
c error exit if too many coeficients demanded
                      if(ilm.gt.mdlmn) go to 999
                      philab(ilm)=m*phil
                      dlmnc(1,ilm)=dlm(il,m,0,cz)
                      dlmnc(2,ilm)=(-1)**m*dlm(il,m,0,-cz)
                    end do
                  end do
c
c FINALLY we sum the products of coefficients
c
                  if(iflag1.eq.1) then
                  ic=0
                  do id=1,ndist
                    il1=ll1(id)
                    in1=n1(id)
                    il2=ll2(id)
                    in2=n2(id)
                    il=ll(id)
c
c initialize reference integer for molecule 1
c
                    il1ref=il1rf(id)
c
c initialize reference integer for molecule 1
c
                    il2ref=il2rf(id)
c
c initialize reference integer for laboratory function
c
                    ilref=ilrf(id)
                    nm1=2.0*il1+1
                    m1min=-il1-1
                    do mm1=1,nm1
                      m1=mm1+m1min
                      m2max=min0(il2,il-m1)
                      m2min=max0(-il2,-il-m1)
                      nm2=m2max-m2min+1
                      m2min=m2min-1
                      do mm2=1,nm2
                        m2=mm2+m2min
                        m=m1+m2
                        ic=ic+1
                        ilmn1=il1ref+m1
                        ilmn2=il2ref+m2
                        ilmn=ilref+m
                        cosphi=cos(chiphi(1,imol1,ilmn1)
     1+chiphi(1,imol2,ilmn2)-philab(ilmn))
                        prod=dlmnmc(1,imol1,ilmn1)*dlmnmc(1,imol2,ilmn2)
                        add=prod*cosphi*dlmnc(1,ilmn)*termcg(ic)
997   format(1x,6(1x,i4),3(1x,e12.5))
                        gofz(iref,id)=gofz(iref,id)+add

                      end do
                    end do
                  end do
                  nsec=nsec+1
                  endif
c if molecules are identical, we add on the case when molecule2 is at the 
c origin
                  if(iflag2.eq.1) then
                  ic=0
                  do id=1,ndist
                    il1=ll1(id)
                    in1=n1(id)
                    il2=ll2(id)
                    in2=n2(id)
                    il=ll(id)
c
c initialize reference integer for molecule 1
c
                    il1ref=il1rf(id)
c
c initialize reference integer for molecule 1
c
                    il2ref=il2rf(id)
c
c initialize reference integer for laboratory function
c
                    ilref=ilrf(id)
                    nm1=2.0*il1+1
                    m1min=-il1-1
                    do mm1=1,nm1
                      m1=mm1+m1min
                      m2max=min0(il2,il-m1)
                      m2min=max0(-il2,-il-m1)
                      nm2=m2max-m2min+1
                      m2min=m2min-1
                      do mm2=1,nm2
                        m2=mm2+m2min
                        m=m1+m2
                        ic=ic+1
                        ilmn1=il1ref+m1
                        ilmn2=il2ref+m2
                        ilmn=ilref+m
                        cosphi=cos(chiphi(1,imol2,ilmn1)
     1+chiphi(1,imol1,ilmn2)-philab(ilmn))
                        prod=dlmnmc(1,imol2,ilmn1)*dlmnmc(1,imol1,ilmn2)
                        add=prod*cosphi*dlmnc(2,ilmn)*termcg(ic)
                        gofz(iref,id)=gofz(iref,id)+add
                      end do
                    end do
                  end do
                  nsec=nsec+1
                  endif
                endif
              endif
              endif
            end do
c     call system_clock(time4,rate,max)
c     write(6,*) rate,max,time4-time3,nsec
         endif
      end do
      write(6,502) nsec,lsum1,lsum2
502   format(1x,i10,2(1x,i5))
500   format(1x,5(i5,1x),5(e12.5,1x))
      do id=1,ndist
            write(6,500) ll1(id),ll2(id),ll(id),n1(id),n2(id)
     1,alsec(id)
      end do
c renormalize g(r)'s to total number of atoms
c density of second atoms in sphere
      rho=float(lsum1)*float(lsum2)/vol
      write(6,193) rho
193   format(1x,'rho = ',e13.6)
      gfac1=gfac*rho
      do id=1,ndist
            z0=0.0
            z03=0.0
            do iz=1,nz
                  zm=z(iz)
                  zdiff=zm-z0
                  z1=zm+zdiff
                  z13=z1*z1*z1
c average number of particles in this volume element
                  rat=gfac1*(z13-z03)
                  gofz(iz,id)=(2*ll1(id)+1)*(2*ll2(id)+1)
     1*gofz(iz,id)/rat
                  if(id.eq.1) gofz(iz,id)=gofz(iz,id)-1.0
                  z0=z1
                  z03=z13
c           write(6,996) iz,zm,gofz(iz,id)
996   format(1x,i5,2(1x,e12.5))
            end do
      end do
      return
999   write(6,998) ic,ilmnlm,ilm,mcg,mdlmn,mdlmn
998   format(1x,'Counters ic, ilmnlm or ilm in SDF',3(1x,i10)
     1/1x,' exceeded their limits', 3(1x,i10))
      stop
      end


