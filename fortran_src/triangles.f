      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_triangles()

      stop
      end      

      subroutine run_triangles()
c
c programme to estimate the distribution of included angles
c between specified triplets of atoms at specified distance ranges
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'triangles.inc'
      integer*4 ifound1,ifound2,ierrout
      character*20 text
c
c setup the triangles symbols
c
      call setup_symbols_triangles()
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
c set up the number of cos(theta) bins
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
c set up the number of triangle distributions
c
      ndist=nfiles
      if(ndist.gt.mdist) ndist=mdist
c
c find the symbol numbers corresponding to atom1, atom2 and atom3
c
      icount1=0
      text='atom1'
      call search_for_symbol(text,icount1)
      icount2=icount1+1
      icount3=icount2+1 
c
c setup the type numbers of the three atoms used to define
c whether two molecules are bonded, and the maximum and minimum distances
c
      do id=1,nfiles
         iref1=ireff(icount1,id)
         iref2=ireff(icount2,id)
         iref3=ireff(icount3,id)
         ifound1=0
         ifound2=0
         ifound3=0
         it=0
         do while(it.lt.ntype)
            it=it+1
            if(type(it).eq.symbol_value(iref1)) ifound1=it
            if(type(it).eq.symbol_value(iref2)) ifound2=it
            if(type(it).eq.symbol_value(iref3)) ifound3=it
         end do
         if(ifound1.eq.0.or.ifound2.eq.0.or.ifound3.eq.0) then
            write(6,100) id
100   format(/'triangles> atoms not matched to '
     *,'.ato file for distribution ',i3)
            return
         endif    
         ltyp1(id)=ifound1
         ltyp2(id)=ifound2
         ltyp3(id)=ifound3
         atom1(id)=type(ifound1)
         atom2(id)=type(ifound2)
         atom3(id)=type(ifound3)
         iref=iref3+1
         read(symbol_value(iref),*) ltyp(id)
         iref=iref+1
         read(symbol_value(iref),*) damin(id,1)
         iref=iref+1
         read(symbol_value(iref),*) damax(id,1)
         iref=iref+1
         read(symbol_value(iref),*) damin(id,2)
         iref=iref+1
         read(symbol_value(iref),*) damax(id,2)
         iref=iref+1
         read(symbol_value(iref),*) rminc(id),rmaxc(id)
!         write(6,*) id,ltyp(id),ltyp1(id),ltyp2(id),ltyp3(id),ltypc
!         write(6,*) (damin(id,i),damax(id,i),i=1,2),rminc(id),rmaxc(id)
         do i=1,2
            damin(id,i)=damin(id,i)/dmax
            damax(id,i)=damax(id,i)/dmax
         end do
         rminc(id)=rminc(id)/dmax
         rmaxc(id)=rmaxc(id)/dmax
      end do
c
c define the z values
c
      zstep=180.0/float(nz)
      zhalf=0.5*zstep
      do j=1,nz
            z(j)=(zstep*j-zhalf)
      end do

c set gofz to zero

      call initgofz(nz,ndist,gofz,mz,mdist)

c normalise atomic coordinates to box size

      call normalise_coordinates()

c set up the distribution of tri-angles

      call find_triangles()

c set the axis title

      titler='theta[deg.]'

c set the column titles

      do id=1,ndist
            ilent1=index(atom1(id),' ')-1
            ilent2=index(atom2(id),' ')-1
            ilent3=index(atom3(id),' ')-1
            titledat(id)='   '//atom1(id)(1:ilent1)
     *//'-'//atom2(id)(1:ilent2)//'-'//atom3(id)(1:ilent3)
      end do

c set the file extension

      ext='c'

c write the accumulated results

      call outgofz(fnameinp,ext,nz,ndist,nsumt,z,gofz,gzt,gzt2
     *,titler,titledat,mz,mdist)

c update the current value of nsumt

      icount=0
      text='nsize'
      call search_for_symbol(text,icount)
      iref=ireff(icount,1)
      write(symbol_value(iref),*) nz
      iref=iref+1
      write(symbol_value(iref),*) nsumt
      iref=iref+4
      write(symbol_value(iref),*) ndist

c write input file

      call save_files_others(14)
      
      !Calculate q value for each distribution
      
      do id=1,ndist
      
         sum1=0.0
         sum2=0.0
         do iz=1,nz
            term=gzt(iz,id)*sind(z(iz))
            q=cosd(z(iz))
            if(ltyp(id).eq.1) then
               q=q+1.0/3.0
               q=q*q
            end if
            sum1=sum1+term
            sum2=sum2+q*term
         end do
         if(sum1.gt.0.0) then
            if(ltyp(id).eq.1) then
               write(14,*) id,0.5*(rminc(id)+rmaxc(id))*dmax,(1.0-(9.0*sum2)/(4.0*sum1))
            else
               write(14,*) id,0.5*(rminc(id)+rmaxc(id))*dmax,(sum2/sum1)
            endif
         else
            write(14,*) id,0.5*(rminc(id)+rmaxc(id))*dmax,0.0
         end if
      end do

      close(14)
      write(6,*) 'Triangles completed successfully'

      return
      end

      subroutine find_triangles()
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'triangles.inc'
      logical carryon
      integer, dimension(:), allocatable  :: atom1list,atom2list,atom3list
      integer                          :: natom1list,natom2list,natom3list
      zstep=180.0/float(nz)
      piconv=180.0/3.141592653
c
c step through all the atoms and ignore those which do not match the required
c ltyp1
c
      do id=1,ndist
!Make a list of all the atoms of type ltyp2(id) which are within the specified distance
! of atoms of type ltypc
         natom1list=0
         natom2list=0
         natom3list=0
         do ia2=1,natom
            if(itype(ia2).eq.ltyp1(id)) natom1list=natom1list+1
            if(itype(ia2).eq.ltyp2(id)) natom2list=natom2list+1
            if(itype(ia2).eq.ltyp3(id)) natom3list=natom3list+1
         end do
!         write(6,*) id,natom1list,natom2list,natom3list
         if(allocated(atom1list)) deallocate(atom1list)
         if(allocated(atom2list)) deallocate(atom2list)
         if(allocated(atom3list)) deallocate(atom3list)
         if(natom1list.gt.0) allocate(atom1list(natom1list))
         if(natom2list.gt.0) allocate(atom2list(natom2list))
         if(natom3list.gt.0) allocate(atom3list(natom3list))
!Repeat this saving the atom numbers                  
         iref1=0
         iref2=0
         iref3=0
         do ia2=1,natom
            if(itype(ia2).eq.ltyp1(id)) then
               iref1=iref1+1
               atom1list(iref1)=ia2
            endif
            if(itype(ia2).eq.ltyp3(id)) then
               iref3=iref3+1
               atom3list(iref3)=ia2
            end if
            carryon=.false.
            if(itype(ia2).eq.ltyp2(id)) then
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
                        d12=sqrt(dsq12)
                        carryon=d12.ge.rminc(id).and.d12.lt.rmaxc(id)
!                    write(6,*) ia2,ia4,itype(ia4),carryon
                     end if
                  end do
               end if
               if(carryon) then
                  iref2=iref2+1
                  atom2list(iref2)=ia2
               end if
            end if
         end do
         !If the type of the first atom is the same as that of the centre atom (ltypc) then we use the minimum and maximum distance for this atom as specified for this pair
         if(ltyp1(id).eq.ltypc) then
            damin(id,1)=rminc(id)
            damax(id,1)=rmaxc(id)
         end if
         natom1list=iref1
         natom2list=iref2
         natom3list=iref3
         write(6,*) id,natom1list,natom2list,natom3list
!Now do main calculation using these lists of ia1, ia2 and ia3 atoms
         do iref1=1,natom1list
            ia1=atom1list(iref1)
            xc1=xadd(mref(ia1))
            yc1=yadd(mref(ia1))
            zc1=zadd(mref(ia1))
            x1=xa(ia1)+xc1
            y1=ya(ia1)+yc1
            z1=za(ia1)+zc1
            do iref2=1,natom2list
               ia2=atom2list(iref2)
               im2=mref(ia2)
               xc2=xadd(im2)
               yc2=yadd(im2)
               zc2=zadd(im2)
               x2=xa(ia2)+xc2
               y2=ya(ia2)+yc2
               z2=za(ia2)+zc2    
c
c calculate atom1 to atom2 distance, applying periodic boundary conditions
c
               xdif=x2-x1
               ydif=y2-y1
               zdif=z2-z1
c apply periodic boundary conditions
               call minimagec(xdif,ydif,zdif)
               if(dorhoz.and.ltyp1(id).eq.ltypc) then
                  dsq12=zdif*zdif
               else if (dorhor.and.ltyp1(id).eq.ltypc) then
                  dsq12=xdif*xdif+ydif*ydif
               else
                  dsq12=xdif*xdif+ydif*ydif+zdif*zdif
               end if
               d12=sqrt(dsq12)
c check that distance is within bounds prescribed for this distribution
               if(d12.ge.damin(id,1).and.d12.le.damax(id,1)) then
c search through all the other atoms for the third atom to define
c the bond angle
                  if(ltyp(id).eq.1) then
                     do iref3=1,natom3list
                        ia3=atom3list(iref3)
                        if(ia3.ne.ia2.and.ia3.ne.ia1) then
                           im3=mref(ia3)
                           xc3=xadd(im3)
                           yc3=yadd(im3)
                           zc3=zadd(im3)
                           x3=xa(ia3)+xc3
                           y3=ya(ia3)+yc3
                           z3=za(ia3)+zc3    
                           xdif=x3-x2
                           ydif=y3-y2
                           zdif=z3-z2
                           call minimagec(xdif,ydif,zdif)
                           dsq23=xdif*xdif+ydif*ydif+zdif*zdif
                           d23=sqrt(dsq23)
                           if(d23.ge.damin(id,2).and.d23.le.damax(id,2)) then
                              !If centre atom is defined, ignore respective coordinates for dorhoz or dorhor
                              if(dorhoz.and.ltyp1(id).eq.ltypc) then
                                 xdif=x3-x2
                                 ydif=y3-y2
                                 zdif=z3-z1
                              else if (dorhor.and.ltyp1(id).eq.ltypc) then
                                 xdif=x3-x1
                                 ydif=y3-y1
                                 zdif=z3-z2
                              else
                                 xdif=x3-x1
                                 ydif=y3-y1
                                 zdif=z3-z1
                              end if
                              call minimagec(xdif,ydif,zdif)
                              dsq13=xdif*xdif+ydif*ydif+zdif*zdif
                              coss=(dsq12+dsq23-dsq13)/(2.0*d12*d23)
                              if(coss.gt.1.0) coss=1.0
                              if(coss.lt.-1.0) coss=-1.0
                              deg=piconv*acos(coss)
                              iref=int(deg/zstep)+1
                              gofz(iref,id)=gofz(iref,id)+1.0
                           endif
                        endif
                     end do
                  else
                     x3=0.0
                     y3=0.0
                     z3=0.0
                     nf=0
                     do ib=1,nbond(ia2)
                        ia3=iabond(ia2,ib)
                        if(itype(ia3).eq.ltyp3(id).and.ia3.ne.ia2) then
                           x3=x3+xa(ia3)
                           y3=y3+ya(ia3)
                           z3=z3+za(ia3)   
                           nf=nf+1
                        endif
                     end do
                     if(nf.gt.0) then
                        x3=x3/nf+xc2
                        y3=y3/nf+yc2
                        z3=z3/nf+zc2
                        xdif=x3-x2
                        ydif=y3-y2
                        zdif=z3-z2
                        call minimagec(xdif,ydif,zdif)
                        dsq23=xdif*xdif+ydif*ydif+zdif*zdif
                        d23=sqrt(dsq23)
                        !If centre atom is defined, ignore respective coordinates for dorhoz or dorhor
                        if(dorhoz.and.ltyp1(id).eq.ltypc) then
                           xdif=x3-x2
                           ydif=y3-y2
                           zdif=z3-z1
                        else if (dorhor.and.ltyp1(id).eq.ltypc) then
                           xdif=x3-x1
                           ydif=y3-y1
                           zdif=z3-z2
                        else
                           xdif=x3-x1
                           ydif=y3-y1
                           zdif=z3-z1
                        end if
                        call minimagec(xdif,ydif,zdif)
                        dsq13=xdif*xdif+ydif*ydif+zdif*zdif
                        coss=(dsq12+dsq23-dsq13)/(2.0*d12*d23)
                        if(coss.gt.1.0) coss=1.0
                        if(coss.lt.-1.0) coss=-1.0
                        deg=piconv*acos(coss)
                        iref=int(deg/zstep)+1
                        gofz(iref,id)=gofz(iref,id)+1.0
                     end if
                  end if
               end if 
            end do
         end do
      end do
      if(allocated(atom1list)) deallocate(atom1list)
      if(allocated(atom2list)) deallocate(atom2list)
      if(allocated(atom3list)) deallocate(atom3list)
c
c normalise each distribution function to have an area of unity
c
      do id=1,ndist
            sumg=0.0
            beta=-0.5*zstep
            do iz=1,nz
                  beta=beta+zstep
                  betar=beta/piconv
                  gofz(iz,id)=gofz(iz,id)/sin(betar)
                  sumg=sumg+gofz(iz,id)
            end do
            if(sumg.gt.0.0) then
                  do iz=1,nz
                        gofz(iz,id)=gofz(iz,id)/(sumg*zstep)
                  end do
            end if
      end do
      return
      end


      
