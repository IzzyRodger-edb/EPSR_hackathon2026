      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_coord()

      end      

      subroutine run_coord()
c
c programme to estimate coordination number histogram for specified atom type
c pairs
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'coord.inc'
      integer*4 ifound1,ifound2,ierrout
      character*20 text
c
c setup the coord symbols
c
      call setup_symbols_coord()
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
c
c set up the number of coordination number size bins
c
      iref=2
      read(symbol_value(iref),*) nz
      nz=min(nz,mz)
c
c set up current accumulator
c
      iref=iref+1
      read(symbol_value(iref),*) nsumt
!
! Set up the centre atom, atom-c, atom-z or atom-r (only 1 can be specified).
! This atom will be used to preselect the middle atom of the triangle to be within a 
! specified distance range of this atom.

      iref=iref+1
      ltypc=0
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      read(symbol_value(iref),'(a)',iostat=ierr) atomc
      !Check this symbol exists in the current file
      it=0
      do while (it.lt.ntype.and.ltypc.eq.0)
         it=it+1
         if(atomc.eq.type(it)) ltypc=it
      end do
c
c set up whether rho-z is to be calculated
c
      iref=iref+1
      dorhoz=.false.
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
! For atomz to work we need to specify an atom type to act as the origin. 
      read(symbol_value(iref),'(a)',iostat=ierr) atomz
      !Check this symbol exists in the current file
      it=0
      do while (it.lt.ntype.and.ltypc.eq.0)
         it=it+1
         if(atomz.eq.type(it)) then
            ltypc=it
            dorhoz=.true.
         end if
      end do
c
c set up whether rho-r is to be calculated. Note that if both atomz and atomr are set, 
! only the first will be used. 
c
      iref=iref+1
      dorhor=.false.
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
! For atomr to work we need to specify an atom type to act as the origin. 
      read(symbol_value(iref),'(a)',iostat=ierr) atomr
      !Check this symbol exists in the current file
      it=0
      do while (it.lt.ntype.and.ltypc.eq.0)
         it=it+1
         if(atomr.eq.type(it)) then
            ltypc=it
            dorhor=.true.
         end if
      end do
      write(6,'(l2,l2,1x,i4)') dorhoz,dorhor,ltypc
c
c check the number of distributions does not exceed the array size
c
      ndist=nfiles
      if(ndist.gt.mdist/2) ndist=mdist/2
      ndist2=2*ndist
c
c setup the type number of the first and second atoms used to define
c the distance range over which the coordination number is to be calculated
c
      do id=1,ndist
            it=0
            ifound1=0
            ifound2=0
            iref1=ireff(nsymboluptodata+1,id)
            iref2=iref1+1
            iref3=iref2+1
            iref4=iref3+1
            do while(it.lt.ntype.and.(ifound1.eq.0.or.ifound2.eq.0))
               it=it+1
               if(type(it).eq.symbol_value(iref1)) ifound1=it
               if(type(it).eq.symbol_value(iref2)) ifound2=it
            end do
            ltyp1(id)=ifound1
            ltyp2(id)=ifound2
            atom1(id)=type(ifound1)
            atom2(id)=type(ifound2)
            read(symbol_value(iref3),*) damin(id)
            read(symbol_value(iref4),*) damax(id)
            write(6,*) id,ltyp1(id),ltyp2(id),damin(id),damax(id)
            iref=iref4+1
            read(symbol_value(iref),*) rminc(id),rmaxc(id)
      end do
      write(6,*) nz,ndist
      zstep=0.5/float(nz)
      zhalf=0.5*zstep
c
c two sets of distributions are calculated: one of coordination number distribution
c the other of mean bond length for each coordination number
c
      do j=1,nz
            z(j)=float(j-1)
      end do

      call initgofz(nz,ndist2,gofz,mz,mdist)

c
c normalise atomic coordinates to box size
c
      call normalise_coordinates()
c
c set up the table of molecular bonds
c
      call find_coord()
c
c axis title
c
      titler='Coordination_number'
      idref2=0
      do id=1,nfiles
            idref1=idref2+1
            idref2=idref1+1
c
c set the column title for this pair of atoms
c
            ilent1=index(atom1(id),' ')-1
            ilent2=index(atom2(id),' ')-1
            write(6,*) ilent1,ilent2,ndist,ndist2
            write(6,100) atom1(id)(1:ilent1),atom2(id)(1:ilent2)
     *,damin(id),damax(id)
            write(titledat(idref1),100) atom1(id)(1:ilent1),atom2(id)(1:ilent2)
     *,damin(id),damax(id)
100   format(1x,a,'-',a,'(',e9.3,',',e9.3,')')
c
c title of corresponding column for bond length distribution
c
            titledat(idref2)=titledat(idref1)
      write(6,*) titledat(id)
      end do

c output the data

      ext='n'

      call outgofz(fnameinp,ext,nz,ndist2,nsumt,z,gofz,gzt,gzt2
     *,titler,titledat,mz,mdist)

c
c calculate average coordination number for each distribution
c
      idref2=0
      do id=1,ndist
            idref1=idref2+1
            idref2=idref1+1
            sum1=0.0
            sum2=0.0
            do i=1,nz
                  term=z(i)*gzt(i,idref1)
                  sum1=sum1+term
                  sum2=sum2+z(i)*term
            end do
            coordn(id,1)=sum1
            std=sum2-sum1*sum1
            if(std.gt.0.0) then
                  coordn(id,2)=sqrt(std)
            else
                  coordn(id,2)=0.0
            endif
      end do            

c
c write input file
c
c
c update the current value of nsumt
c
      icount=0
      text='nsize'
      call search_for_symbol(text,icount)
      iref=ireff(icount,1)
      write(symbol_value(iref),*) nz
      iref=iref+1
      write(symbol_value(iref),*) nsumt
      iref=iref+4
      write(symbol_value(iref),*) ndist

      call save_files_others(14)

      do id=1,ndist
c
c for each partial correlation function write the two atom types, 
c the average coordination number, and it rms deviation
c
            if (ltypc.eq.0) then
               write(14,101) atom1(id),atom2(id),damin(id),damax(id),coordn(id,1),coordn(id,2)
            else
               write(14,101) atom1(id),atom2(id),damin(id),damax(id),rminc(id),rmaxc(id),coordn(id,1),coordn(id,2)
            end if
      end do
101   format(1x,a3,1x,a3,4(1x,f10.5))

      close(14)

      return
      end

      subroutine find_coord()
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'coord.inc'
      real*4 bav(matom)
      integer*4 mbond(matom)
      logical carryon
      integer, dimension(:), allocatable  :: atom1list,atom2list
      integer                          :: natom1list,natom2list

      zstep=1
      pi=3.141592653
c volume of sphere of radius 0.5
      gfac=4.0*pi/3.0
      vol=gfac*0.125*dmax*dmax*dmax
      do id=1,ndist
c
c define the reference integers for gofz
c
         idref1=2*id-1
         idref2=2*id
!Make a list of all the atoms of type ltyp2(id) which are within the specified distance
! of atoms of type ltypc
         natom1list=0
         natom2list=0
         do ia2=1,natom
            if(itype(ia2).eq.ltyp1(id)) natom1list=natom1list+1
            if(itype(ia2).eq.ltyp2(id)) natom2list=natom2list+1
         end do
!         write(6,*) id,natom1list,natom2list,natom3list
         if(allocated(atom1list)) deallocate(atom1list)
         if(allocated(atom2list)) deallocate(atom2list)
         if(natom1list.gt.0) allocate(atom1list(natom1list))
         if(natom2list.gt.0) allocate(atom2list(natom2list))
!Repeat this saving the atom numbers                  
         iref1=0
         iref2=0
         do ia2=1,natom
            if(itype(ia2).eq.ltyp1(id)) then
!Need to select those middle atoms which are within specified distance range of atomc, if defined
               carryon=.true.
               if(ltypc.gt.0) then
                  carryon=.false.
                  im2=mref(ia2)
                  xc2=xadd(im2)
                  yc2=yadd(im2)
                  zc2=zadd(im2)
                  x2=xa(ia2)+xc2
                  y2=ya(ia2)+yc2
                  z2=za(ia2)+zc2    
                  ia4=0
                  do while(ia4.lt.natom.and..not.carryon)
                     ia4=ia4+1
                     if(itype(ia4).eq.ltypc) then
                        im4=mref(ia4)
                        xdif=x2-(xa(ia4)+xadd(im4))
                        ydif=y2-(ya(ia4)+yadd(im4))
                        zdif=z2-(za(ia4)+zadd(im4))
c apply periodic boundary conditions
                        call minimagec(xdif,ydif,zdif)
                        if(dorhoz) then
                           dsq12=zdif*zdif
                        else if (dorhor) then
                           dsq12=xdif*xdif+ydif*ydif
                        else
                           dsq12=xdif*xdif+ydif*ydif+zdif*zdif
                        end if
                        d12=sqrt(dsq12)*dmax
                        carryon=d12.ge.rminc(id).and.d12.lt.rmaxc(id)
!                    write(6,*) ia2,ia4,itype(ia4),carryon
                     end if
                  end do
               end if
               if(carryon) then
                  iref1=iref1+1
                  atom1list(iref1)=ia2
               end if
            end if
            if(itype(ia2).eq.ltyp2(id)) then
               iref2=iref2+1
               atom2list(iref2)=ia2
            endif
         end do
         natom1list=iref1
         natom2list=iref2
c
c set up for each molecule the number of bonded neighbours.
c
         do i=1,natom
            mbond(i)=0
            bav(i)=0.0
         end do
c
c number of origin atoms
c
         momol=0
c
c number of bonded atoms
c
         mbmol=0
c
c step through all origin molecules
c
         do iref1=1,natom1list
            ia1=atom1list(iref1)
            im1=mref(ia1)
            momol=momol+1
            xc1=xadd(im1)
            yc1=yadd(im1)
            zc1=zadd(im1)
            xa1=xa(ia1)+xc1
            ya1=ya(ia1)+yc1
            za1=za(ia1)+zc1
c now look at all the other molecules and determine which ones are 
c bonded
            do iref2=1,natom2list
                  ia2=atom2list(iref2)
                  im2=mref(ia2)
c
c calculate atom1 to atom2 distance, applying periodic boundary conditions
c
                  xc2=xadd(im2)
                  yc2=yadd(im2)
                  zc2=zadd(im2)
                  x2=xa(ia2)+xc2
                  y2=ya(ia2)+yc2
                  z2=za(ia2)+zc2    
                  xdif=x2-xa1
                  ydif=y2-ya1
                  zdif=z2-za1
c apply periodic boundary conditions
                  call minimagec(xdif,ydif,zdif)
                  dsq12=xdif*xdif+ydif*ydif+zdif*zdif
                  da12=sqrt(dsq12)*dmax
c check that distance is within bounds prescribed for this distribution
                  if(da12.ge.damin(id).and.da12.le.damax(id)) then
                      mbond(ia1)=mbond(ia1)+1
                      bav(ia1)=bav(ia1)+da12
                      mbmol=mbmol+1
                  endif 
            end do
         end do
c
c accumulate the histogram for this molecule
c
         do iref1=1,natom1list
            ia1=atom1list(iref1)
            if(mbond(ia1).gt.0) then
                  bav(ia1)=bav(ia1)/float(mbond(ia1))
            else
                  bav(ia1)=0.0
            endif
            iref=mbond(ia1)+1
            if(iref.gt.0.and.iref.lt.nz) then
                  gofz(iref,idref1)=gofz(iref,idref1)+1
                  gofz(iref,idref2)=gofz(iref,idref2)+bav(ia1)
            endif
         end do
         write(6,*) mbmol,momol
c
c calculate the average number of bonds per molecule and normalise histograms
c
         if(momol.gt.0) then
            avbond=float(mbmol)/float(momol)
            do iref=1,nz
                  if(gofz(iref,idref1).gt.0.0) then
                    gofz(iref,idref2)=gofz(iref,idref2)
     */gofz(iref,idref1)
                  else
                    gofz(iref,idref2)=0.0
                  endif
                  gofz(iref,idref1)=gofz(iref,idref1)/float(momol)
            end do
         else
            avbond=0.0
         endif
         write(6,994) ltyp2(id),ltyp1(id),avbond
994   format(1x,'Average coord. number for type ',i2,' around type '
     *,i2,' = ',f10.5)
      end do
      return
      end


