      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_partials()

      end

      subroutine run_partials()
c
c programme to estimate number of clusters in a distribution
c of molecules. Atoms are defined to be in a cluster when they are within a
c specified distance range of a neighbour
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      include 'partials.inc'
      integer*4 ierrout
      character*20 text
c
c setup the paritals symbols
c
      call setup_symbols_partials()
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
c set up the number of cluster size bins
c
      iref=2
      read(symbol_value(iref),*) nz
      nz=min(nz,mz)
! Set up range of calculation

      iref=iref+1
      rnormalise=1.0
      read(symbol_value(iref),*) rnormalise
      if(rnormalise.gt.1.0) rnormalise=1.0
c
c set up current accumulator
c
      iref=iref+1
      read(symbol_value(iref),*) nsumt
!
! Set up the centre atom, atom-c, atom-z or atom-r (only 1 can be specified).
! This atom will be used to preselect the centre atom about which densities will be calculated
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
      if(ndist.gt.mdist) ndist=mdist
      write(6,*) nz,ndist
c
c set up number of partial g(r)'s to calculate
c
      ndist=0
      do it=1,ntype
            do jt=1,ntype
               if(ndist.lt.mpairs) then
                  !Restricted to specified atom at origin
                  if(ltypc.gt.0) then
                     if(it.eq.ltypc) then
                        ndist=ndist+1
                        ltyp1(ndist)=it
                        ltyp2(ndist)=jt
                        atom1(ndist)=type(it)
                        atom2(ndist)=type(jt)
                     end if
                  else if(jt.ge.it) then
                     !Use all atom pairs
                     ndist=ndist+1
                     ltyp1(ndist)=it
                     ltyp2(ndist)=jt
                     atom1(ndist)=type(it)
                     atom2(ndist)=type(jt)
                  end if
               endif
            end do
      end do
      do id=1,ndist
         write(6,'(1x,i4,2(1x,a3))') id,atom1(id),atom2(id)
      end do
c
c setup step size and z values
c
c For dorhoz and dorhor set up factor to capture true distance range in each case
      if(dorhoz) rnormalise=rnormalise*clatz/dmax
      if(dorhor) rnormalise=rnormalise*min(alatx,blaty)/dmax
      zstep=0.5*rnormalise/float(nz)
      zhalf=0.5*zstep
      do j=1,nz
            z(j)=(zstep*j-zhalf)*dmax
      end do

      call initgofz(nz,ndist,gofz,mz,mdist)

c
c normalise atomic coordinates to box size
c
      call normalise_coordinates()
c
c set up the table of molecular bonds
c
      call find_partials()
c
c axis title
c
      titler='r'
      do id=1,ndist
c
c set the column title for this pair of atoms
c
            ilent1=index(atom1(id),' ')-1
            ilent2=index(atom2(id),' ')-1
            titledat(id)='   '//atom1(id)(1:ilent1)//'-'//atom2(id)(1:ilent2)
      end do

c output the data

      ext='g'

      call outgofz(fnameinp,ext,nz,ndist,nsumt,z,gofz,gzt,gzt2
     *,titler,titledat,mz,mdist)

c
c write input file
c

c update the current value of nsumt and ndist

      icount=0
      text='nz'
      call search_for_symbol(text,icount)
      iref=ireff(icount,1)
      write(symbol_value(iref),*) nz
      iref=iref+2
      write(symbol_value(iref),*) nsumt
      iref=iref+4
      write(symbol_value(iref),*) ndist

      call save_files_others(14)

      do id=1,ndist
c
c for each partial correlation function write the two atom type numbers needed.
c
            write(14,104) id,atom1(id),atom2(id)
      end do
104   format(1x,i5,5x,a3,1x,a3)

      close(14)

      return
      end

      subroutine find_partials()
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'partials.inc'
      include 'lattice.inc'
      integer*4 lsum1(mpairs),lsum2(mpairs)
!      integer, dimension(:), allocatable     :: xlist,ylist,zlist,xlistsave,ylistsave,zlistsave
      logical doloop
      pi=4.0*atan(1.0)
c volume of sphere of radius 0.5
      if(dorhoz) then
         gfac=1.0
         vol=2.0*alatx*blaty*dmax*dmax !The factor of 2 is needed because we symmetrise about z = 0
      else if(dorhor) then
         gfac=pi
         vol=clatz*dmax
      else
         gfac=4.0*pi/3.0
         vol=alatx*blaty*clatz*dmax*dmax*dmax
      endif
c
c step through all the atoms and ignore those which do not match the required
c ltyp1
c
c set up counters for number of origin and non-origin atoms for each partial
c
      do id=1,ndist
      lsum1(id)=0
      lsum2(id)=0
      do ia1=1,natom
          if(itype(ia1).eq.ltyp2(id)) lsum2(id)=lsum2(id)+1
          if(itype(ia1).eq.ltyp1(id)) then
            im1=mref(ia1)
            xc1=xadd(im1)
            yc1=yadd(im1)
            zc1=zadd(im1)
            x1=xa(ia1)+xc1
            y1=ya(ia1)+yc1
            z1=za(ia1)+zc1
            doloop=.true.
c            if(lsum1(id).gt.0) then
c               doloop=.not.((rhoz.eq.1.and.(abs(z1-z1s)).lt.0.05/dmax)
c     1.or.(rhoxy.eq.1.and.(abs(x1-x1s)+abs(y1-y1s)).lt.0.05/dmax))
c            end if
            if(doloop) then
               lsum1(id)=lsum1(id)+1
               do ia2=1,natom
                  im2=mref(ia2)
                  if(itype(ia2).eq.ltyp2(id).and.im1.ne.im2) then
c
c calculate atom1 to atom2 distance, applying periodic boundary conditions
c
                     xc2=xadd(im2)
                     yc2=yadd(im2)
                     zc2=zadd(im2)
                     x2=xa(ia2)+xc2
                     y2=ya(ia2)+yc2
                     z2=za(ia2)+zc2
                     xdif=x2-x1
                     ydif=y2-y1
                     zdif=z2-z1
c apply periodic boundary conditions
                     call minimagec(xdif,ydif,zdif)
                     if(dorhoz) then
                        d12=abs(zdif)
                     else if(dorhor) then
                        dsq12=xdif*xdif+ydif*ydif
                        d12=sqrt(dsq12)
                     else
                        dsq12=xdif*xdif+ydif*ydif+zdif*zdif
                        d12=sqrt(dsq12)
                     end if
                     iref=int(d12/zstep)+1
c check that distance is within bounds prescribed for this distribution
                     if(iref.gt.0.and.iref.le.nz) then
                        gofz(iref,id)=gofz(iref,id)+1.0
                     end if
                  end if
               end do
            end if
            !Save the coordinates of this atom. For rhoz=1 or rhoxy=1, if the next atom of this type
            !has barely changed values of z (rhoz=1) or (x,y) (rhoxy=1) then we don't need to include it
            !since the second atoms will already have been included in these cases
            x1s=x1
            y1s=y1
            z1s=z1
          endif
      end do
      end do
c renormalize g(r)'s to total number of atoms
      do id=1,ndist
c density of second atoms in sphere
            if(dorhoz.or.dorhor) then
               rho=real(lsum1(id))*vol
            else if(ltypc.gt.0) then
!               rho=1.0
               rho=real(lsum1(id))
            else
               rho=real(lsum1(id))*real(lsum2(id))/vol
            end if
            write(6,193) id,real(lsum1(id)),real(lsum2(id)),rho
c           do iz=100,nz
c                 write(6,193) iz,z(iz),gofz(iz,id)
c           end do
193   format(1x,i5,5(1x,e13.6))
            gfac1=gfac*rho
            z0=0.0
            z03=0.0
            do iz=1,nz
                  zm=z(iz)
                  zdiff=zm-z0
                  z1=zm+zdiff
                  if(dorhoz) then
                     z13=z1
                  else if(dorhor) then
                     z13=z1*z1
                  else
                     z13=z1*z1*z1
                  endif
c average number of particles in this volume element
                  rat=gfac1*(z13-z03)
                  gofz(iz,id)=gofz(iz,id)/rat
                  z0=z1
                  z03=z13
            end do
      end do
      return
      end



