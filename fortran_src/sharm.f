      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_sharm()

      stop
      end      

      subroutine run_sharm()
c
c performs calculation of spherical harmonic coefficients from .ato file
c
c uses arithmetic centre of molecule rather than a particular atom to define
c molecular centre. 
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
      real pi,piconv
      call fctrln(250)
c
c setup the SHARM symbols
c
      call setup_symbols_sharm()
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
c set up l2 values
c
      call parse(symbol_value(7),nwords,ncf,ncl,1,256,nworddim)
      nll2v=nwords
      do i=1,nll2v
         read(symbol_value(7)(ncf(i):ncl(i)),*,iostat=ierr) ll2v(i)
      end do
c
c set up l values
c
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
c for molecule 1 read number and type of atoms used to define the molecular
c centre
c
      call parse(symbol_value(11),nwords,ncf,ncl,5,256,nworddim)
      ncent1=0
c
c skip through the atom types to find out the type number of the specified atom
c
      do ic=1,nwords
         it=0
         ifound=0
         do while(it.lt.ntype.and.ifound.eq.0)
            it=it+1
            if(type(it).eq.symbol_value(11)(ncf(ic):ncl(ic))) ifound=it
         end do
         if(ifound.le.0) then
            write(6,201) symbol_value(11)(ncf(ic):ncl(ic))
     *,fname_ato(1:ilenfname_ato)
201   format(/'sharm> Atom type ',a,' does not exist in .ato file ',a)
            return
         else
            if(ncent1.lt.mcentre) then
               ncent1=ncent1+1
               lctyp1(ncent1)=ifound
            endif
         endif
      end do
c
c for molecule 1 read primary axis type (1=X, 2=Y, 3=Z), and
c the number and atom numbers of atoms used to define the molecule's nzt1
c axis
c
      call parse(symbol_value(12),nwords,ncf,ncl,5,256,nworddim)
      nzt1=symbol_value(12)(ncf(1):ncl(1))
      ncz1=nwords-1
      if(ncz1.gt.mcentre) ncz1=mcentre
      do i=1,ncz1
         iref=i+1
         read(symbol_value(12)(ncf(iref):ncl(iref)),*) lztyp1(i)
      end do
c
c for molecule 1 read secondary axis type (1=X, 2=Y, 3=Z), and
c the number and atom numbers of atoms used to define the
c molecule's nyt1 axis
c
      call parse(symbol_value(13),nwords,ncf,ncl,5,256,nworddim)
      nyt1=symbol_value(13)(ncf(1):ncl(1))
      ncy1=nwords-1
      if(ncy1.gt.mcentre) ncy1=mcentre
      do i=1,ncy1
         iref=i+1
         read(symbol_value(13)(ncf(iref):ncl(iref)),*) lytyp1(i)
      end do

c for molecule 2 read number and type of atoms used to define the molecular
c centre. Also check whether this molecule is identical to the first

      call parse(symbol_value(14),nwords,ncf,ncl,5,256,nworddim)
      ncent2=0
c
c skip through the atom types to find out the type number of the specified atom
c
      do ic=1,nwords
         it=0
         ifound=0
         do while(it.lt.ntype.and.ifound.eq.0)
            it=it+1
            if(type(it).eq.symbol_value(14)(ncf(ic):ncl(ic))) ifound=it
         end do
         if(ifound.eq.0) then
            write(6,201) symbol_value(14)(ncf(ic):ncl(ic))
     *,fname_ato(1:ilenfname_ato)
            return
         else
            if(ncent2.lt.mcentre) then
               ncent2=ncent2+1
               lctyp2(ncent2)=ifound
            endif
         endif
      end do

c test for identity of central and second molecules

      ldsnct=0
      if(ncent2.eq.ncent1) then
         ldsnct=0
         do it=1,ncent1
            if(lctyp2(it).ne.lctyp1(it)) ldsnct=1
         end do
      else
         ldsnct=1
      endif

      write(6,*) ncent1,ncent2,ldsnct



      if(ldsnct.ne.0) then

c if molecules are different, then repeat this for the second molecule

c for molecule 2 read primary axis type (1=X, 2=Y, 3=Z), and
c the number and atom numbers of atoms used to define the molecule's nzt1
c axis

         call parse(symbol_value(15),nwords,ncf,ncl,5,256,nworddim)
         nzt2=symbol_value(15)(ncf(1):ncl(1))
         ncz2=nwords-1
         if(ncz2.gt.mcentre) ncz2=mcentre
         do i=1,ncz2
            iref=i+1
            read(symbol_value(15)(ncf(iref):ncl(iref)),*) lztyp2(i)
         end do

c for molecule 2 read secondary axis type (1=X, 2=Y, 3=Z), and
c the number and atom numbers of atoms used to define the
c molecule's nyt1 axis

         call parse(symbol_value(16),nwords,ncf,ncl,5,256,nworddim)
         nyt2=symbol_value(16)(ncf(1):ncl(1))
         ncy2=nwords-1
         if(ncy2.gt.mcentre) ncy2=mcentre
         do i=1,ncy2
            iref=i+1
            read(symbol_value(16)(ncf(iref):ncl(iref)),*) lytyp2(i)
         end do
      else

c otherwise set all the values to the same as the central molecule

         ncent2=ncent1
         nzt2=nzt1
         nyt2=nyt1
         if(ncent2.gt.0) then
            do ic=1,ncent1
               lctyp2(ic)=lctyp1(ic)
            end do
         end if
         ncz2=ncz1
         do i=1,ncz1
            lztyp2(i)=lztyp1(i)
         end do
         ncy2=ncy1
         do i=1,iabs(ncy1)
            lytyp2(i)=lytyp1(i)
         end do
      endif
c     write(6,*) ncent1,ncent2,ncz1,ncy1,ncz2,ncy2
      !Decide whether to use simulated phim, thetam, chim
      call parse(symbol_value(17),nwords,ncf,ncl,5,256,nworddim)
      performtest=symbol_value(17)(ncf(1):ncl(1)).ne.'notest'.and.nwords.gt.6
      if(performtest) then
         read(symbol_value(17)(ncf(1):ncl(7)),*,iostat=ierr) thetallab,phillab,phimmid,thetammid,chimmid,phispread,thetaspread
         thetallab=thetallab*piconv
         phillab=phillab*piconv
         phimmid=phimmid*piconv
         thetammid=thetammid*piconv
         chimmid=chimmid*piconv
         phispread=phispread*piconv
         thetaspread=thetaspread*piconv
      end if
      if(ierr.ne.0) performtest=.false.
c
c only perform the calculation when nsum is zero or a multiple of nmult
c
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

            call sharm(zmaxd,ndist,nin,fname_inp)
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
      if(performtest) write(6,*) 'sharm> Test performed.'

      return
      end

      subroutine sharm(zmaxd,ndist,nout,fname_inp)
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      include 'lmncom.inc'
      include 'sharm.inc'
      include 'sharm_sdf.inc'

c     integer*4 rate,max,time1,time2,time3,time4
      real amass,amasssum

		logical :: centralmolecule,found,foundatom
		character(len=*) :: fname_inp
		character(len=256) :: fname_xyz,textline !Filename for xyz data, line of text
		integer :: nout !Unit number on which to write xyz file
      real :: xt,yt,zt,xr,yr,zr
      real, dimension(:), allocatable           :: xacentral,yacentral,zacentral
      integer :: ncentral,natocentral
      character(len=3), dimension(:), allocatable  :: atsymbolcentral
      logical, dimension(:), allocatable :: centraldefinition
      real phimt(mmol),cthetmt(mmol),chimt(mmol),dlmnmct(mmol,mdlmn),chiphit(mmol,mdlmn)
      real pi,piconv,thetatest,u1,u2

      if(performtest) then
         open(10,file='ran1.dat',status='unknown')
         call initran1(10)
         close (10)
      end if

		pi=4.0*atan(1.0)
      piconv=pi/180.0
      degconv=180.0/pi

      zstep=zmaxd/float(nz)
c volume of sphere of radius 0.5
      gfac=4.0*pi/3.0
      vol=(alatx*(blaty*clatz))
c set up counter for second atoms
      nsec=0
      do id=1,ndist
            alsec(id)=0.0
            nmsec(id)=0
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
      do imol1=1,nmol
c define origin of this molecule
            xo1=0.0
            yo1=0.0
            zo1=0.0
            amasssum=0.0
            ncref=0
            nm(imol1)=0
c
c test to see if this is molecule 1 or molecule 2. Need to check that all the specified types occur - not just one
c
            found=.true.
            ic=0
            do while(ic.lt.ncent1.and.found)
               ic=ic+1
               foundatom=.false.
               ia1=naf(imol1)-1
               do while(ia1.lt.nal(imol1).and..not.foundatom)
                  ia1=ia1+1
                  foundatom=itype(ia1).eq.lctyp1(ic)
               end do
               if(ic.eq.1) then
                  found=foundatom
               else
                  found=found.and.foundatom
               end if
            end do
            if(found) then
               nm(imol1)=1
            else !molecule 1 not found, try molecule 2
               found=.true.
               ic=0
               do while(ic.lt.ncent2.and.found)
                  ic=ic+1
                  foundatom=.false.
                  ia1=naf(imol1)-1
                  do while(ia1.lt.nal(imol1).and..not.foundatom)
                     ia1=ia1+1
                     foundatom=itype(ia1).eq.lctyp2(ic)
                  end do
                  if(ic.eq.1) then
                     found=foundatom
                  else
                     found=found.and.foundatom
                  end if
               end do
               if(found) then
                  nm(imol1)=2
               endif
            endif
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
            if(nm(imol1).eq.1.and.ncent1.gt.0) then
              do ia1=naf(imol1),nal(imol1)
                do ic=1,ncent1
                  if(itype(ia1).eq.lctyp1(ic)) then
                  amass=amassa(itype(ia1))
                  xo1=xo1+xa(ia1)*amass
                  yo1=yo1+ya(ia1)*amass
                  zo1=zo1+za(ia1)*amass
                  amasssum=amasssum+amass
                  endif
                end do
              end do
            else if(nm(imol1).eq.2.and.ncent2.gt.0) then
              do ia1=naf(imol1),nal(imol1)
                do ic=1,ncent2
                  if(itype(ia1).eq.lctyp2(ic)) then
                  amass=amassa(itype(ia1))
                  xo1=xo1+xa(ia1)*amass
                  yo1=yo1+ya(ia1)*amass
                  zo1=zo1+za(ia1)*amass
                  amasssum=amasssum+amass
                  endif
                end do
              end do
            else 
              do ia1=naf(imol1),nal(imol1)
                  amass=amassa(itype(ia1))
                  xo1=xo1+xa(ia1)*amass
                  yo1=yo1+ya(ia1)*amass
                  zo1=zo1+za(ia1)*amass
                  amasssum=amasssum+amass
              end do
            endif
            if(amasssum.gt.0.0.and.nm(imol1).ne.0) then
              xo1=xo1/amasssum
              yo1=yo1/amasssum
              zo1=zo1/amasssum
c     write(6,991) xo1,yo1,zo1
              xcent(imol1)=xadd(imol1)+xo1
              ycent(imol1)=yadd(imol1)+yo1
              zcent(imol1)=zadd(imol1)+zo1
c
c step through all the atoms and ignore those which do not match the required
c ltype's
c
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
                naf1=naf(imol1)-1
                do iz=1,ncz1
                  ia1=naf1+lztyp1(iz)
                  ifindz=ifindz+1
                  czx=czx+xa(ia1)
                  czy=czy+ya(ia1)
                  czz=czz+za(ia1)
                end do
                ncy=ncy1
                do iy=1,iabs(ncy1)
                  ia1=naf1+lytyp1(iy)
                  ifindy=ifindy+1
                  cyx=cyx+xa(ia1)
                  cyy=cyy+ya(ia1)
                  cyz=cyz+za(ia1)
                end do
              else
                naf2=naf(imol1)-1
                do iz=1,ncz2
                  ia1=naf2+lztyp2(iz)
                  ifindz=ifindz+1
                  czx=czx+xa(ia1)
                  czy=czy+ya(ia1)
                  czz=czz+za(ia1)
                end do
                ncy=ncy2
                do iy=1,ncy2
                  ia1=naf2+lytyp2(iy)
                  ifindy=ifindy+1
                  cyx=cyx+xa(ia1)
                  cyy=cyy+ya(ia1)
                  cyz=cyz+za(ia1)
                end do
              endif
              if(ifindz.gt.0) then
                ct(1)=czx/ifindz-xo1
                ct(2)=czy/ifindz-yo1
                ct(3)=czz/ifindz-zo1
              else
                ct(1)=0.0
                ct(2)=0.0
                      ct(3)=1.0
              endif
              if(ifindy.gt.0) then
                pt(1)=cyx/ifindy-xo1
                pt(2)=cyy/ifindy-yo1
                pt(3)=cyz/ifindy-zo1
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
!If testing orientations, supply a fictitious but known distribution of angles
                  if(performtest) then
                     call populate_disk(pi,thetammid,phimmid,thetaspread,phispread,thetatemp,phimt(imol1))
                     cthetmt(imol1)=cos(thetatemp)
!                     u1=ran1()
!                     amp=-2.0*log(u1)
!                     u2=ran1()
                     chimt(imol1)=chimmid
!+sqrt(amp)*cos(u2*2.0*pi)*phispread
                     write(6,499) imol1,phimt(imol1),cthetmt(imol1),chimt(imol1)
                  end if
              else
                  cthetm(imol1)=0.0
                  phim(imol1)=0.0
                  chim(imol1)=0.0
              endif
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
                  chiphi(imol1,ilmnlm)=m*phim(imol1)+in*chim(imol1)
                  termd=dlm(il,m,in,cthetm(imol1))
                  dlmnmc(imol1,ilmnlm)=termd
                  if(performtest) then
                     chiphi(imol1,ilmnlm)=0.0
                     termd=dlm(il,m,in,1.0)
                     dlmnmc(imol1,ilmnlm)=termd
                     chiphit(imol1,ilmnlm)=m*phimt(imol1)+in*chimt(imol1)
                     termd=dlm(il,m,in,cthetmt(imol1))
                     dlmnmct(imol1,ilmnlm)=termd
                  end if
c                 alsec(ilmnlm)=alsec(ilmnlm)+termd*termd
                  end do
                end do
              end do
c           write(6,499) imol1,phim(imol1),cthetm(imol1),chim(imol1)
c     1,dlmnmc(imol1,1)
            else
              nm(imol1)=0
            endif
499   format(1x,i5,4(1x,e12.5))
      end do
      write(6,333) nsmmo1,nsmmo2
333   format(1x,'Number of type 1 molecules = ',i10
     1/1x, 'Number of type 2 molecules = ',i10)
!      stop
! Write out coordinates of centralmolecule. We will use the average coordinates of all these molecules.
		im=0
      ncentral=0
		centralmolecule=.false.
		do while(im.lt.nmol)
			im=im+1
			centralmolecule=nm(im).eq.1
			if(centralmolecule) then !Write out coordinates of this molecule
            ncentral=ncentral+1
            if(ncentral.eq.1) then !First time allocate the arrays to store summed coordinates
               if(allocated(xacentral)) deallocate(xacentral,yacentral,zacentral,atsymbolcentral,centraldefinition)
               natocentral=nato(im)
               allocate(xacentral(natocentral),yacentral(natocentral),zacentral(natocentral)
     *,atsymbolcentral(natocentral),centraldefinition(natocentral))
            end if
! Shift coordinates to centre used by SHARM/SDF
			   xshift=xadd(im)-xcent(im)
				yshift=yadd(im)-ycent(im)
				zshift=zadd(im)-zcent(im)
				do ia=naf(im),nal(im)
! Temporary coordinates
					xt=xa(ia)+xshift
					yt=ya(ia)+yshift
					zt=za(ia)+zshift
! Need to rotate these to align with the laboratory axes - performed in reverse order about the laboratory axes
               call rotateeuler(xt,yt,zt,-chim(im),-acos(cthetm(im)),-phim(im),xr,yr,zr)
               iaref=ia-naf(im)+1
               if(ncentral.eq.1) then
                  xacentral(iaref)=xr
                  yacentral(iaref)=yr
                  zacentral(iaref)=zr
                  atsymbolcentral(iaref)=atsymbol(itype(ia))
! Check whether this atom is part of the central molecule definition
						found=.false.
              		ic=0
						do while(.not.found.and.ic.lt.ncent1)
							ic=ic+1
							found=itype(ia).eq.lctyp1(ic)
	               end do
						iz=0
						do while (.not.found.and.iz.lt.ncz1)
							iz=iz+1
							found=iaref.eq.lztyp1(iz)
						end do
						iy=0
						do while (.not.found.and.iy.lt.ncy1)
							iy=iy+1
							found=iaref.eq.lytyp1(iy)
						end do
                  centraldefinition(iaref)=found
               else
                  xacentral(iaref)=xacentral(iaref)+xr
                  yacentral(iaref)=yacentral(iaref)+yr
                  zacentral(iaref)=zacentral(iaref)+zr
               end if
            end do
         end if
      end do
! Form average coordinates for central molecule
      if(ncentral.gt.0) then
         do ia=1,natocentral
            xacentral(ia)=xacentral(ia)*dmax/real(ncentral)
            yacentral(ia)=yacentral(ia)*dmax/real(ncentral)
            zacentral(ia)=zacentral(ia)*dmax/real(ncentral)
         end do
      	fname_xyz=fname_inp(1:indexlastdot(fname_inp)-1)//'-first.xyz'
   		open(nout,file=fname_xyz,status='unknown',iostat=ierr)
   		if(ierr.eq.0) then
! Write number of atoms in molecule
   			write(nout,'(i6)') natocentral
! Write information line - not used for plotting
   			write(nout,'(a)') 'Central molecule coordinate file: '//fname_xyz(1:len_trim(fname_xyz))
   			do ia=1,natocentral
               write(textline,'(a3,1x,f10.5,1x,f10.5,1x,f10.5)') atsymbolcentral(ia),xacentral(ia),yacentral(ia),zacentral(ia)
   				! Append a comment if this atom is part of the central molecule definition
   				if(centraldefinition(ia)) then
   					textline=textline(1:len_trim(textline))//'  #CENTRALMOLECULE'
   				end if
   				write(nout,'(a)') textline(1:len_trim(textline))
   			end do
   			close(nout)
   		end if
         if(ldscnt.eq.0) then !Also write this out as the second molecule
         	fname_xyz=fname_inp(1:indexlastdot(fname_inp)-1)//'-second.xyz'
      		open(nout,file=fname_xyz,status='unknown',iostat=ierr)
      		if(ierr.eq.0) then
! Write number of atoms in molecule
      			write(nout,'(i6)') natocentral
! Write information line - not used for plotting
      			write(nout,'(a)') 'Central molecule coordinate file: '//fname_xyz(1:len_trim(fname_xyz))
      			do ia=1,natocentral
                  write(textline,'(a3,1x,f10.5,1x,f10.5,1x,f10.5)') atsymbolcentral(ia),xacentral(ia),yacentral(ia),zacentral(ia)
      				! Append a comment if this atom is part of the central molecule definition
      				if(centraldefinition(ia)) then
      					textline=textline(1:len_trim(textline))//'  #CENTRALMOLECULE'
      				end if
      				write(nout,'(a)') textline(1:len_trim(textline))
      			end do
      			close(nout)
      		end if
         end if
      end if
! Now repeat this for second molecules if they are different from first molecules
! Write out coordinates of centralmolecule. We will use the average coordinates of all these molecules.
		im=0
      ncentral=0
		centralmolecule=.false.
		do while(im.lt.nmol)
			im=im+1
			centralmolecule=nm(im).eq.2
			if(centralmolecule) then !Write out coordinates of this molecule
            ncentral=ncentral+1
            if(ncentral.eq.1) then !First time allocate the arrays to store summed coordinates
               if(allocated(xacentral)) deallocate(xacentral,yacentral,zacentral,atsymbolcentral,centraldefinition)
               natocentral=nato(im)
               allocate(xacentral(natocentral),yacentral(natocentral),zacentral(natocentral)
     *,atsymbolcentral(natocentral),centraldefinition(natocentral))
            end if
! Shift coordinates to centre used by SHARM/SDF
			   xshift=xadd(im)-xcent(im)
				yshift=yadd(im)-ycent(im)
				zshift=zadd(im)-zcent(im)
				do ia=naf(im),nal(im)
! Temporary coordinates
					xt=xa(ia)+xshift
					yt=ya(ia)+yshift
					zt=za(ia)+zshift
! Need to rotate these to align with the laboratory axes - performed in reverse order about the laboratory axes
               call rotateeuler(xt,yt,zt,-chim(im),-acos(cthetm(im)),-phim(im),xr,yr,zr)
               iaref=ia-naf(im)+1
               if(ncentral.eq.1) then
                  xacentral(iaref)=xr
                  yacentral(iaref)=yr
                  zacentral(iaref)=zr
                  atsymbolcentral(iaref)=atsymbol(itype(ia))
! Check whether this atom is part of the central molecule definition
						found=.false.
              		ic=0
						do while(.not.found.and.ic.lt.ncent1)
							ic=ic+1
							found=itype(ia).eq.lctyp1(ic)
	               end do
						iz=0
						do while (.not.found.and.iz.lt.ncz1)
							iz=iz+1
							found=iaref.eq.lztyp1(iz)
						end do
						iy=0
						do while (.not.found.and.iy.lt.ncy1)
							iy=iy+1
							found=iaref.eq.lytyp1(iy)
						end do
                  centraldefinition(iaref)=found
               else
                  xacentral(iaref)=xacentral(iaref)+xr
                  yacentral(iaref)=yacentral(iaref)+yr
                  zacentral(iaref)=zacentral(iaref)+zr
               end if
            end do
         end if
      end do
! Form average coordinates for central molecule
      if(ncentral.gt.0) then
         do ia=1,natocentral
            xacentral(ia)=xacentral(ia)*dmax/real(ncentral)
            yacentral(ia)=yacentral(ia)*dmax/real(ncentral)
            zacentral(ia)=zacentral(ia)*dmax/real(ncentral)
         end do
      	fname_xyz=fname_inp(1:indexlastdot(fname_inp)-1)//'-second.xyz'
   		open(nout,file=fname_xyz,status='unknown',iostat=ierr)
   		if(ierr.eq.0) then
! Write number of atoms in molecule
   			write(nout,'(i6)') natocentral
! Write information line - not used for plotting
   			write(nout,'(a)') 'Central molecule coordinate file: '//fname_xyz(1:len_trim(fname_xyz))
   			do ia=1,natocentral
               write(textline,'(a3,1x,f10.5,1x,f10.5,1x,f10.5)') atsymbolcentral(ia),xacentral(ia),yacentral(ia),zacentral(ia)
   				! Append a comment if this atom is part of the central molecule definition
   				if(centraldefinition(ia)) then
   					textline=textline(1:len_trim(textline))//'  #CENTRALMOLECULE'
   				end if
   				write(nout,'(a)') textline(1:len_trim(textline))
   			end do
   			close(nout)
   		end if
      end if
c
c initialize counters for number of molecules of each kind in the system
c
      if(ldsnct.eq.0) then
            if(nm(nmol).eq.1.or.nm(nmol).eq.2) then
                  lsum1=1
                  lsum2=1
            else
                  lsum1=0
                  lsum2=0
            endif
      else 
            lsum1=0
            lsum2=0
            if(nm(nmol).eq.1) lsum1=1
            if(nm(nmol).eq.2) lsum2=1
      endif
      write(6,501) nmol,nm(nmol),lsum1,lsum2
      do imol1=1,nmol-1
        if(nm(imol1).ne.0) then
            write(6,*) imol1
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
            do imol2=nm2f,nmol
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
c angle phi in lab coordinate system:
                    phil=atan2(sphil,cphil)
                  else
                    phil=0.0
                  endif
                  if(performtest) then
                     call populate_disk(pi,thetallab,phillab,thetaspread,phispread,thetatemp,phil)
                     cz=cos(thetatemp)
                  end if
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
                  if(iflag1.eq.1) then !Molecule 1 is at the origin
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
c initialize reference integer for molecule 2
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
                        if(performtest) then
                           cosphi=cos(chiphi(imol1,ilmn1)
     1+chiphit(imol2,ilmn2)-philab(ilmn))
                           prod=dlmnmc(imol1,ilmn1)*dlmnmct(imol2,ilmn2)
                        else
                           cosphi=cos(chiphi(imol1,ilmn1)
     1+chiphi(imol2,ilmn2)-philab(ilmn))
                           prod=dlmnmc(imol1,ilmn1)*dlmnmc(imol2,ilmn2)
                        end if
                        add=prod*cosphi*dlmnc(1,ilmn)*termcg(ic)
997   format(1x,6(1x,i4),3(1x,e12.5))
                        gofz(iref,id)=gofz(iref,id)+add
                      end do
                    end do
                  end do
                  nsec=nsec+1
                  endif
                  if(iflag2.eq.1) then !Molecule 2 is at the origin
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
                        if(performtest) then
                           cosphi=cos(chiphi(imol2,ilmn1)
     1+chiphit(imol1,ilmn2)-philab(ilmn))
                           prod=dlmnmc(imol2,ilmn1)*dlmnmct(imol1,ilmn2)
                        else
                           cosphi=cos(chiphi(imol2,ilmn1)
     1+chiphi(imol1,ilmn2)-philab(ilmn))
                           prod=dlmnmc(imol2,ilmn1)*dlmnmc(imol1,ilmn2)
                        end if
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
     1,nmsec(id),alsec(id)
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
      if(performtest) then
         open(10,file='ran1.dat',status='unknown')
         call saveran1(10)
         close (10)
      end if
      return
999   write(6,998) ic,ilmnlm,ilm,mcg,mdlmn,mdlmn
998   format(1x,'Counters ic, ilmnlm or ilm in BANG',3(1x,i10)
     1/1x,' exceeded their limits', 3(1x,i10))
      stop
      end

      subroutine populate_disk(pi,theta0,phi0,thetaspread,phispread,thetanew,phinew)
! Randomly populates a tangential disk at (theta0,phi0) with a Gaussian-like distribution at constant radius,
! giving a new (theta,phi)
      real pi,ctheta0,phi0,thetaspread,phispread,thetanew,phinew
      real x0,y0,z0,z02,deltax,deltay,deltaz,stheta0,root1,root2,cthetanew,sthetanew,u1,u2,amp,cons

      ctheta0=cos(theta0)
      stheta0=sin(theta0)
      x0=stheta0*cos(phi0)
      y0=stheta0*sin(phi0)
      z0=ctheta0
      z02=z0*z0
      cons=-1.0
      do while(cons.lt.0.0)
         u1=ran1()
         u2=ran1()
         amp=sqrt(-2.0*log(u1))
         deltax=amp*cos(2.0*pi*u2)*phispread
         deltay=amp*sin(2.0*pi*u2)*thetaspread
         cons=z02-(deltax*(deltax+2.0*x0)+deltay*(deltay+2.0*y0))
      end do
! Choose the root closest to z0
      cthetanew=sqrt(cons)
      root1=acos(cthetanew)
      root2=acos(-cthetanew)
      dif1=abs(theta0-root1)/thetaspread
      dif2=abs(theta0-root2)/thetaspread
! If either or these differences is greater than 5 deviations, choose the other, otherwise toss a coin.
      if((dif2-dif1).gt.5.0) then
         thetanew=root1
      else if((dif1-dif2).gt.5.0) then
         thetanew=root2
      else !Toss a coin to choose between the two roots
         if(ran1(1).lt.0.5) then
            thetanew=root1
         else
            thetanew=root2
         end if
      end if
      sthetanew=sqrt(1.0-cthetanew*cthetanew)
      cphinew=(x0+deltax)/sthetanew
      sphinew=(y0+deltay)/sthetanew
      phinew=atan2(sphinew,cphinew)
!      write(6,*) x0,y0,z0,deltax,deltay,root1,root2,dif1,dif2
!      write(6,*) theta0,phi0,thetaspread,phispread,thetanew,phinew
      return
      end










