      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_torangles()

      stop
      end      

      subroutine run_torangles()
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
      include 'torangles.inc'
      integer*4 ifound,ierrout
      character*20 text
      character*50 formatout
      character*2 fmti1,fmti2,fmti3,ltypfmt,ltyp1fmt

c
c setup the triangles symbols
c
      call setup_symbols_torangles()
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
c set up the number of phi bins (between -180 and 180)
c
      iref=2
      read(symbol_value(iref),*) nz
      nz=min(nz,mz)
c
c set up current accumulator
c
      iref=iref+1
      read(symbol_value(iref),*) nsumt

c
c set up the number of triangle distributions
c
      ndist=nfiles
      if(ndist.gt.mdist) ndist=mdist
c
c find the symbol numbers corresponding to atom1
c
      icount=0
      text='atom1'
      call search_for_symbol(text,icount)
c
c setup the type numbers of the atom used to define
c which molecules to use
c
      do id=1,nfiles
         iref=ireff(icount,id)
         ifound=0
         it=0
         do while(it.lt.ntype.and.ifound.eq.0)
            it=it+1
            if(type(it).eq.symbol_value(iref)) ifound=it
         end do
         if(ifound.eq.0) then
            write(6,100) id
100   format(/'torangles> atom not matched to '
     *,'.ato file for distribution ',i3)
            return
         endif    
         nfirst(id)=ifound
         atom1(id)=type(ifound)
c
c for each torsion angle read the two base atom numbers needed to define the 
c bond 
c
         iref=iref+1
         read(symbol_value(iref),*) ltyp(id)
         iref=iref+1
         read(symbol_value(iref),*) ltyp1(id)
c
c search through the .ato file to find out the types of atoms ltyp and ltyp1
c
         im=0
         ifound=0
         do while(im.lt.nmol.and.ifound.eq.0)
            im=im+1
            if(atom1(id).eq.type(itype(naf(im)))) then
               ifound=im
               atom2(id)=type(itype(naf(im)+ltyp(id)-1))
               atom3(id)=type(itype(naf(im)+ltyp1(id)-1))
            endif
         end do
                           
c
c and for each base atom the two atom numbers needed to define the torsion angle.
c (These can be the same atom if desired so that only 1 atom at each end is needed
c to define the torsion angle.)
c
         iref=iref+1
         read(symbol_value(iref),*) (ltyp2(id,j),j=1,2)
         iref=iref+1
         read(symbol_value(iref),*) (ltyp3(id,j),j=1,2)
      end do
c
c define the z values
c
      zstep=360/float(nz)
      zhalf=0.5*zstep
      do j=1,nz
            z(j)=zstep*j-zhalf-180.0
      end do

c set gofz to zero

      call initgofz(nz,ndist,gofz,mz,mdist)

c normalise atomic coordinates to box size

      call normalise_coordinates()

c set up the distribution of torsional angles

      call find_torangles()

c set the axis title

      titler='Phi[deg]'

c set the column titles

      fmti1='i1'
      fmti2='i2'
      fmti3='i3'
      do id=1,ndist
            if(ltyp(id).lt.10) then
                  ltypfmt=fmti1
            else if(ltyp(id).lt.100) then
                  ltypfmt=fmti2
            else if(ltyp(id).lt.1000) then
                  ltypfmt=fmti3
            endif
            if(ltyp1(id).lt.10) then
                  ltyp1fmt=fmti1
            else if(ltyp1(id).lt.100) then
                  ltyp1fmt=fmti2
            else if(ltyp1(id).lt.1000) then
                  ltyp1fmt=fmti3
            endif
            write(formatout,332) ltypfmt,ltyp1fmt
332   format('(3x,a,2h:[,a,',a2,',1h-,a,',a2,',1h])')
      write(6,*) formatout
            ilent1=index(atom1(id),' ')-1
            ilent2=index(atom2(id),' ')-1
            ilent3=index(atom3(id),' ')-1
            write(titledat(id),fmt=formatout) atom1(id)(1:ilent1)
     *,atom2(id)(1:ilent2),ltyp(id),atom3(id)(1:ilent3),ltyp1(id)
      write(6,*) titledat(id)
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
      iref=iref+1
      write(symbol_value(iref),*) ndist

c write input file

      call save_files_others(14)

      close(14)
      write(6,*) 'Torangles completed successfully'

      return
      end

      subroutine find_torangles()
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'torangles.inc'
      real*4 cd(3),co1(3),co2(3),cv1(3),cv2(3)
     1,cpar(3)
      real*4 phisave(mdist)
      character*80 formatout
      pi=3.141592653
      zstep=360.0/float(nz)
c
c step through all the molecules
c
      ifindm=0
      open(29,file='toranglesvalues.txt',status='unknown',iostat=ierr)

c Setup the output format

      write(formatout,332) '(1x,i6,',ndist,'(1x,e14.7))'
332   format(a7,i4,a11)
      do im1=1,nmol
        do id=1,ndist
          phisave(id)=0.0
        end do
        do id=1,ndist
          if(itype(naf(im1)).eq.nfirst(id)) then
            ifindm=ifindm+1
            ia1=naf(im1)+ltyp(id)-1
            ia2=naf(im1)+ltyp1(id)-1
            x1=xa(ia1)
            y1=ya(ia1)
            z1=za(ia1)
            x2=xa(ia2)
            y2=ya(ia2)
            z2=za(ia2)
c
c define the separation vector cd between the two base atoms and normalize it to 
c unit value
c
            cd(1)=(x2-x1)
            cd(2)=(y2-y1)
            cd(3)=(z2-z1)
            dmag=sqrt(scapro(cd,cd))
            cd(1)=cd(1)/dmag
            cd(2)=cd(2)/dmag
            cd(3)=cd(3)/dmag
c
c now determine the orientation vectors co for base atom 1 and base atom 2, using 
c the values of ltyp2 and ltyp3 
c
            io1a=naf(im1)+ltyp2(id,1)-1
            io1b=naf(im1)+ltyp2(id,2)-1
            io2a=naf(im1)+ltyp3(id,1)-1
            io2b=naf(im1)+ltyp3(id,2)-1
            co1(1)=0.5*(xa(io1a)+xa(io1b))-x1
            co1(2)=0.5*(ya(io1a)+ya(io1b))-y1
            co1(3)=0.5*(za(io1a)+za(io1b))-z1
            co2(1)=0.5*(xa(io2a)+xa(io2b))-x2
            co2(2)=0.5*(ya(io2a)+ya(io2b))-y2
            co2(3)=0.5*(za(io2a)+za(io2b))-z2
c
c define new vectors normal to both cd and co
c
            call vecpro(co1,cd,cv1)
            call vecpro(co2,cd,cv2)
c
c renormalize these to unit vectors
c
            dmag=sqrt(scapro(cv1,cv1))
            cv1(1)=cv1(1)/dmag
            cv1(2)=cv1(2)/dmag
            cv1(3)=cv1(3)/dmag
            dmag=sqrt(scapro(cv2,cv2))
            cv2(1)=cv2(1)/dmag
            cv2(2)=cv2(2)/dmag
            cv2(3)=cv2(3)/dmag
c calculate the cosine of the angle between cv1 and cv2
            cosphi=scapro(cv1,cv2)
            if(cosphi.gt.1.0) cosphi=1.0
            if(cosphi.lt.-1.0) cosphi=-1.0
c
c calculate the vector product of cv1 and cv2 - this should be either parallel
c (phi positive) or antiparallel to cd (phi negative)
c
            call vecpro(cv1,cv2,cpar)
            sinphi=scapro(cpar,cd)
            phis=atan2(sinphi,cosphi)*180/pi
            phisave(id)=phis
c           write(6,991) ia1,ia2,io1a,io1b,io2a,io2b,phi,phis
            iref=int((phis+180.0)/zstep)+1
            gofz(iref,id)=gofz(iref,id)+1.0
991   format(1x,6(i4,1x),2(e13.6,1x))
          endif
        end do
        write(29,fmt=formatout) im1,(phisave(id),id=1,ndist)
      end do
c
c ensure final distributions are normalised to the sum of all terms and calculate average
c chain length
c
      do id=1,ndist
            sum0=0.0
            do i=1,nz
                  sum0=sum0+gofz(i,id)
            end do
            sum0=sum0/float(nz)
            if(sum0.ne.0.0) then
                  do i=1,nz
                        gofz(i,id)=gofz(i,id)/sum0
                  end do
            endif
      end do            
      write(6,945) ifindm
945   format(1x,'No of molecules found = ',i5)
      return
      end   

c     subroutine vecpro(a,b,c)
c     real*4 a(3),b(3),c(3)
c performs the vector product of c1 x c2 and puts result in c3
c     c(1)=a(2)*b(3)-a(3)*b(2)
c     c(2)=a(3)*b(1)-a(1)*b(3)
c     c(3)=a(1)*b(2)-a(2)*b(1)
c     return
c     end

c     function scapro(a,b)
c     real*4 a(3),b(3)
c     scapro=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
c     return
c     end



