      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call getbragginten()

      stop
      end      

      subroutine getbragginten()
c
c Program bragginten
c
c Calculates the Bragg intensities from a specified .ato file.
c Intensities are separated into pairs of atom types.
c
c Currently set up for 70 components (number of atom types)
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      include 'weights.inc'
      include 'X-rayf0.inc'

      character(len=256) :: fname,outputfile,formatout
      integer :: multa,multb,multc !Number of unit cells along each direction a, b, c.
      integer :: nbragg,nin,nout,ierrout
      real, allocatable, dimension(:) :: qbragg       !Q scale for calculation
      real, allocatable, dimension(:,:) :: intentemp,intenbragg !Scattered intensity
      real, dimension(mcomp) :: scatlen
      real :: qstep,qmax,qwidth,sum,hklqmin,rmax
      integer :: nfile,ib,it,jt,ip,iref

c Get parameters from input file
c
c setup the xwts symbols
c
      call setup_symbols_bragginten()
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

      !Output file name
      outputfile=fname_run(1:ilenfname_run)
      ncharoutput=ilenfname_run
c
c normalise atomic coordinates to box size
c
      call normalise_coordinates()
      write(6,100)
100   format(/'bragginten> Program to calculate Bragg intensities'
     1,' for the specified .ato file.'/)
99    format(a80)
      iref=1
c
c Read underlying multiplicities along the a, b and c axes - this will have been defined when setting up the .ato file
c
      iref=iref+1
      read(symbol_value(iref),*) multa,multb,multc

c Minimum Q from (h,k,l) values

      iref=iref+1
      read(symbol_value(iref),*) hklqmin

c Maximum value of r - used for convolution function

      iref=iref+1
      read(symbol_value(iref),*) rmax

cc Get q step and q max

      iref=iref+1
      read(symbol_value(iref),*) qstep,qmax

c Get width of broadening function

      iref=iref+1
      read(symbol_value(iref),*) qwidth
      
      iref=iref+1

! Get the scattering lengths

      it=0
      do while(it.lt.ntype)
         it=it+1
         iref=iref+2
         read(symbol_value(iref),*) scatlen(it)
      end do
      
      !Set up the (neutron) weights
      ip=0
      do it=1,ntype
         do jt=it,ntype
            ip=ip+1
            cwt(1,ip)=scatlen(it)*scatlen(jt)
         end do
      end do
      !Set up Q scale
      nfile=1
      xfile(1)=0
      nbragg=nint(qmax/qstep)+1
      if(allocated(qbragg)) deallocate(qbragg)
      if(allocated(intentemp)) deallocate(intentemp)
      if(allocated(intenbragg)) deallocate(intenbragg)
      allocate(qbragg(nbragg),intentemp(npairs,nbragg),intenbragg(nfile,nbragg))
      do ib=1,nbragg
         qbragg(ib)=qstep*real(ib-1)
      end do
      call bragginten(qstep,qmax,qwidth,rmax,multa,multb,multc
     *,hklqmin,nbragg,qbragg,intentemp)
      !Combine the data according the supplied weights
      do im=1,nfile
         if(xfile(im).gt.0) then
            do ib=1,nbragg
               sum=0.0
               do ip=1,npairs
                  sum=sum+weightx(ib,ip,xfile(im))*intentemp(ip,ib)
               end do
               intenbragg(im,ib)=sum
            end do
         else
            do ib=1,nbragg
               sum=0.0
               do ip=1,npairs
                  sum=sum+cwt(im,ip)*intentemp(ip,ib)
               end do
               intenbragg(im,ib)=sum
            end do
         end if
      end do

c Create a filename for the output

      fname=directory_name(1:nchardir)//outputfile(1:ncharoutput)//'.inten'
      open(10,file=fname,status='unknown')
      write(formatout,103) nfile
103   format('(1x,f10.5,',i2.2,'(1x,e12.5))')
      write(6,*) formatout(1:index(formatout,' '))
      do ib=1,nbragg
         write(10,formatout) qbragg(ib),(intenbragg(im,ib),im=1,nfile)
      end do
      close(10)
      
c
c Rewrite the input file in case anything was changed
c
c      call save_files_others(14)

      if(allocated(qbragg)) deallocate(qbragg)
      if(allocated(intentemp)) deallocate(intentemp)
      if(allocated(intenbragg)) deallocate(intenbragg)
      return
      end


