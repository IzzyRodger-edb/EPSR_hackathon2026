      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call bondsclus()

      end      

      subroutine bondsclus()

      use reallocation_routines
c
c calculates average distance between each atom pair within each molecule
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
      character*3 chtype
      character*256 fnamebonds
      integer natommol
      character(len=3), dimension(:), allocatable   :: alabelt
      real, dimension(:,:), allocatable             :: bond
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
499   format(/'bonds> No .ato file specified! - routine will exit'/)
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
c generate name of output file  - .bonds
c
      fnamebonds=directory_name(1:nchardir)//readfile(1:ilenreadfile)//'.bonds'
      ilenfnameato=nchardir+ilenreadfile+6
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
599   format(/'Problems with .ato file ',a/'Routine will exit now'/)
         return
      endif
c
c open output file
c
      open(11,file=fnamebonds,status='unknown')
c
c zero bond array
c
      ibcnt=0
      chtype='   '
      iaref=0
      im=0
      do while(im.lt.nmol)
         im=im+1
         if(alabel(naf(im)).ne.chtype) then
c
c A new molecule has been found
c
c if ibcnt.ne.0 then write out previous bond lengths
c
            if(ibcnt.gt.0) then
               write(6,*) ibcnt,natommol,nsum
               write(11,104) (i,alabelt(i),i=1,natommol)
               write(11,*)
104   format(//7x,40(1x,i5,a3,2x))     
               do i=1,natommol-1
                  do j=i+1,natommol
                     bond(i,j)=bond(i,j)/real(nsum)
                     bond(j,i)=bond(j,i)/real(nsum)-bond(i,j)**2
                     if(bond(j,i).gt.0.0) then
                        bond(j,i)=sqrt(bond(j,i))
                     else
                        bond(j,i)=0.0
                     endif
                  end do
               end do
               do i=1,natommol
                  write(11,105) i,alabelt(i),(bond(i,j),j=1,natommol)
105   format(i3,a3,1x,(40(1x,e10.3)))
               end do
               write(11,*)
            endif
c
c save new values
c
            ibcnt=ibcnt+1
            chtype=alabel(naf(im))
            natommol=nato(im)
            if(allocated(alabelt)) deallocate(alabelt)
            if(allocated(bond)) deallocate(bond)
            allocate(alabelt(natommol))
            allocate(bond(natommol,natommol))
            do i=1,natommol
                  iref=i+naf(im)-1
                  alabelt(i)=alabel(iref)
            end do
c
c start new bond values
c
            nsum=1
            do i=1,natommol-1
               bond(i,i)=0.0
               iref=i+naf(im)-1
               x1=xa(iref)
               y1=ya(iref)
               z1=za(iref)
               do j=i+1,natommol
                  bond(j,j)=0.0
                  jref=j+naf(im)-1
                  x2=xa(jref)
                  y2=ya(jref)
                  z2=za(jref)
                  xd=x2-x1
                  yd=y2-y1
                  zd=z2-z1
                  dist=xd*xd+yd*yd+zd*zd
                  bond(i,j)=sqrt(dist)
                  bond(j,i)=dist
               end do
            end do
         else
            nsum=nsum+1
            do i=1,natommol-1
               iref=i+naf(im)-1
               x1=xa(iref)
               y1=ya(iref)
               z1=za(iref)
               do j=i+1,natommol
                  jref=j+naf(im)-1
                  x2=xa(jref)
                  y2=ya(jref)
                  z2=za(jref)
                  xd=x2-x1
                  yd=y2-y1
                  zd=z2-z1
                  dist=xd*xd+yd*yd+zd*zd
                  bond(i,j)=bond(i,j)+sqrt(dist)
                  bond(j,i)=bond(j,i)+dist
               end do
            end do
         endif
      end do
c
c write out the results from the last molecule
c
      write(6,*) ibcnt,natommol,nsum
      write(11,104) (i,alabelt(i),i=1,natommol)
      write(11,*)
      do i=1,natommol-1
         do j=i+1,natommol
            bond(i,j)=bond(i,j)/float(nsum)
            bond(j,i)=bond(j,i)/real(nsum)-bond(i,j)**2
            if(bond(j,i).gt.0.0) then
               bond(j,i)=sqrt(bond(j,i))
            else
               bond(j,i)=0.0
            endif
         end do
      end do
      do i=1,natommol
         write(11,105) i,alabelt(i),(bond(i,j),j=1,natommol)
      end do
      write(11,*)
      close (11)
      return
      end
