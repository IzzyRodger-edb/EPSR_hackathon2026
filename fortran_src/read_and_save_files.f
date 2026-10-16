      subroutine read_input_file(nin,nout,ierrout)
c
c reads the EPSR input files
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      integer*4 ierrout,ifnameato,ifnamepcof,ilensymbol_value
      character*20 command
c
c don't do anything if there is an error in the input file - will start out with
c default values in this case
c
      if(ierrout.ne.0) then
         close(nin)
         return
      endif

c
c set up the symbol numbers for the fnameato and fnamepcof symbols
c
      ifnameato=0
      ifnamepcof=0
c
c set the number of potential coefficient pairs to zero - this is used to
c signal whether a valid .pcof file has been read in
c
      npairspcof=0
c
c signal whether any atoms have been read in
c
      ntype=0
c
c Read in commands from the already opened input file, and process each command
c in turn. At the end close this file.
c
      command='next'
      icount=1
      ifilecount=1
c
c we assume that the symbols are specified in order, but the
c process_commands routine will check each command in turn to ensure the correct
c values are updated. The last one should be a 'q' to signify the end of
c command input
c
      ierr=0
      do while(command.ne.'e'.and.ierr.eq.0)
c
c check the status of the file counter
c
         if(icount.gt.nsymboldata.and.ifilecount.lt.nfiles) then
            ifilecount=ifilecount+1
            icount=nsymboluptodata+1
         else if(icount.gt.nsymbol) then
            icount=1
         endif
c
c save the symbol number for fnameato and fnamepcof in case these files need to
c be read as well
c
         if(symbol(icount).eq.'fnameato') ifnameato=icount
         if(symbol(icount).eq.'fnamepcof') ifnamepcof=icount
c
c read the commands
c
         nwords=0
         read(nin,100,iostat=ierr) sometext
100   format(a)
         call parse(sometext,nwords,ncf,ncl,12,nstsize,nworddim)

c Need to ignore blank lines read from an input file

         if(nwords.gt.0) call process_commands(nin,nout,command,icount,ifilecount,ierrout,1)
!         write(6,*) 'read_input_file> ',icount,ifilecount
      end do

c
c finished input - close the input file
c
      close(nin)
      return
      end

      subroutine initgofz(nz,ndist,gofz,nxdim,nfdim)

c
c initialisation routine for histograms
c
      integer*4 nz,ndist,nxdim,nfdim
      real*4 gofz(nxdim,nfdim)
c
c zero gofz
c
      do j=1,nz
            do id=1,ndist
                  gofz(j,id)=0.
            end do
      end do
      return
      end

      subroutine initgofzs(nz,ndist,gofz,gofzs,nxdim,nfdim)

c
c initialisation routine for histograms
c
      integer*4 nz,ndist,nxdim,nfdim
      real*4 gofz(nxdim,nfdim),gofzs(nxdim,nfdim)
c
c zero gofz
c
      do j=1,nz
            do id=1,ndist
                  gofz(j,id)=0.
                  gofzs(j,id)=0.
            end do
      end do
      return
      end

      subroutine outgofz(fnameinp,ext,nz,ndist,nsum,z,gofz,gzt,gzt2
     *,titler,titledat,nxdim,nfdim)

c
c summing and output routine for histograms
c
      integer*4 nz,ndist,nsum,nxdim,nfdim,ierr,lptin,lpt,ilen
      real*4 z(nxdim),gofz(nxdim,nfdim),gzt(nxdim,nfdim),gzt2(nxdim,nfdim)
      character*14 titler
      character*29 titledat(nfdim)
      character*256 fnameinp,fname
      character*1 ext
      real*8 fc1,fc2,a,b,a1,a2,b1,b2,error

c
c set up the name of the output file with extension .n
c
c strip off the last . extension and add a new one
c
      ilen=indexlastdot(fnameinp)
      fname=fnameinp(1:ilen)//ext(1:1)
c
c read in data to be fitted if nsum is greater than zero
c
      ierr=0
      if(nsum.gt.0) then
         call r_int(lptin,fname,ndist,z,gzt,gzt2,nxdim,nfdim,ierr)
      else
         lptin=nz
      endif
c
c if an error occurs on opening this file, set nsum to -1 to start a new
c accumulation
c
      if(ierr.ne.0) nsum=-1
      if(nsum.gt.0) then
            do id=1,ndist
                  do j=1,nz
                        if(gzt2(j,id).lt.1.0e-20) gzt2(j,id)=0.0
                        if(abs(gzt(j,id)).lt.1.0e-20) gzt(j,id)=0.0
c check that number of points read in is not less than number expected
                        if(j.le.lptin) then
                          gzt2(j,id)=(gzt2(j,id)**2+gzt(j,id)**2)
                        else
                          gzt(j,id)=0.0
                          gzt2(j,id)=0.0
                        endif
                  end do
            end do
      else
            do id=1,ndist
                  do j=1,nz
                        gzt(j,id)=0.0
                        gzt2(j,id)=0.0
                  end do
            end do
      endif
      if(nsum.ge.0) then
            nsum=nsum+1
            fc1=1.0/float(nsum)
            fc2=1.0-fc1
      else
            fc1=1.0
            fc2=0.0
      endif
      lpt=nz                              !#points
      do id=1,ndist
            do i=1,nz
                  a2=fc2*gzt(i,id)
                  a1=fc1*gofz(i,id)
                  a=a1+a2
                  b2=fc2*gzt2(i,id)
                  b1=fc1*(gofz(i,id))**2
                  b=b1+b2
                  error=b-a**2
                  if(error.gt.0.) then
                        error=dsqrt(error)
                  else
                        error=0.
                  endif
                  gzt(i,id)=a
                  gzt2(i,id)=error
            end do
      end do
      call w_int(lpt,fname,ndist,z,gzt,gzt2,titler,titledat,nxdim,nfdim)      !read data
      return
      end

      subroutine outgofz_sharm(fnameinp,ext,nz,ndist,nsum,z,gofz,gzt,gzt2
     *,titler,titledat,nxdim,nfdim)

c
c summing and output routine for histograms
c
      integer*4 nz,ndist,nsum,nxdim,nfdim,ierr,lptin,lpt,ilen
      integer*4 ll1(nfdim),ll2(nfdim),ll(nfdim),n1(nfdim),n2(nfdim)
      real*4 z(nxdim),gofz(nxdim,nfdim),gzt(nxdim,nfdim),gzt2(nxdim,nfdim)
      character*14 titler
      character*29 titledat(nfdim)
      character*256 fnameinp,fname
      character*1 ext
      real*8 fc1,fc2,a,b,a1,a2,b1,b2,error

c
c set up the name of the output file with extension .n
c
c strip off the last . extension and add a new one
c
      ilen=indexlastdot(fnameinp)
      fname=fnameinp(1:ilen)//ext(1:1)
c
c read in data to be fitted if nsum is greater than zero
c
      ierr=0
      if(nsum.gt.0) then
         call r_int(lptin,fname,ndist,z,gzt,gzt2,nxdim,nfdim,ierr)
      else
         lptin=nz
      endif
c
c if an error occurs on opening this file, set nsum to -1 to start a new
c accumulation
c
      if(ierr.ne.0) nsum=-1
      if(nsum.gt.0) then
            do id=1,ndist
                  do j=1,nz
                        if(gzt2(j,id).lt.1.0e-20) gzt2(j,id)=0.0
                        if(abs(gzt(j,id)).lt.1.0e-20) gzt(j,id)=0.0
c check that number of points read in is not less than number expected
                        if(j.le.lptin) then
                          gzt2(j,id)=(gzt2(j,id)**2+gzt(j,id)**2)
                        else
                          gzt(j,id)=0.0
                          gzt2(j,id)=0.0
                        endif
                  end do
            end do
      else
            do id=1,ndist
                  do j=1,nz
                        gzt(j,id)=0.0
                        gzt2(j,id)=0.0
                  end do
            end do
      endif
      if(nsum.ge.0) then
            nsum=nsum+1
            fc1=1.0/float(nsum)
            fc2=1.0-fc1
      else
            fc1=1.0
            fc2=0.0
      endif
      lpt=nz                              !#points
      do id=1,ndist
            do i=1,nz
                  a2=fc2*gzt(i,id)
                  a1=fc1*gofz(i,id)
                  a=a1+a2
                  b2=fc2*gzt2(i,id)
                  b1=fc1*(gofz(i,id))**2
                  b=b1+b2
                  error=b-a**2
                  if(error.gt.0.) then
                        error=dsqrt(error)
                  else
                        error=0.
                  endif
                  gzt(i,id)=a
                  gzt2(i,id)=error
            end do
      end do
      call w_int(lpt,fname,ndist,z,gzt,gzt2,titler,titledat,nxdim,nfdim)      !read data
      return
      end

      subroutine outgofzs_sharm(fnameinp,ext,nz,ndist,nsum,z,gofz,gofzs,gzt,gzt2
     *,titler,titledat,nxdim,nfdim)

c
c summing and output routine for histograms
c
      integer*4 nz,ndist,nsum,nxdim,nfdim,ierr,lptin,lpt,ilen
      integer*4 ll1(nfdim),ll2(nfdim),ll(nfdim),n1(nfdim),n2(nfdim)
      real*4 z(nxdim),gofz(nxdim,nfdim),gofzs(nxdim,nfdim),gzt(nxdim,nfdim),gzt2(nxdim,nfdim)
      character*14 titler
      character*29 titledat(nfdim)
      character*256 fnameinp,fname
      character*1 ext
      real*8 fc1,fc2,a,b,a1,a2,b1,b2,error

c
c set up the name of the output file with extension .n
c
c strip off the last . extension and add a new one
c
      ilen=indexlastdot(fnameinp)
      fname=fnameinp(1:ilen)//ext(1:1)
c
c read in data to be fitted if nsum is greater than zero
c
      ierr=0
      if(nsum.gt.0) then
         call r_int(lptin,fname,ndist,z,gzt,gzt2,nxdim,nfdim,ierr)
      else
         lptin=nz
      endif
c
c if an error occurs on opening this file, set nsum to -1 to start a new
c accumulation
c
      if(ierr.ne.0) nsum=-1
      if(nsum.gt.0) then
            do id=1,ndist
                  do j=1,nz
c check that number of points read in is not less than number expected
                        if(j.gt.lptin) then
                          gzt(j,id)=0.0
                          gzt2(j,id)=0.0
                        endif
                  end do
            end do
      else
            do id=1,ndist
                  do j=1,nz
                        gzt(j,id)=0.0
                        gzt2(j,id)=0.0
                  end do
            end do
      endif
      if(nsum.ge.0) then
            nsum=nsum+1
            fc1=1.0/real(nsum)
            fc2=1.0-fc1
      else
            fc1=1.0
            fc2=0.0
      endif
      lpt=nz                              !#points
      do id=1,ndist
            do i=1,nz
                  a2=fc2*gzt(i,id)
                  a1=fc1*gofz(i,id)
                  b2=fc2*gzt2(i,id)
                  b1=fc1*gofzs(i,id)
                  gzt(i,id)=a1+a2
                  gzt2(i,id)=b1+b2
            end do
      end do
      call w_int(lpt,fname,ndist,z,gzt,gzt2,titler,titledat,nxdim,nfdim)      !read data
      return
      end

      subroutine save_files(nin,nout)
      use jmk_readline_module
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      integer*4 ilenuptodot
      character*256 promptstring
      integer*4 nin,nout,ierrfile,lenpromptstring
c
c first check whether the existing file name is O.K.
c
      ierrfile=1
      do while(ierrfile.ne.0)

         if(ilenfname_inp.gt.0) then

            write(nout,303) prompt(1:ilenprompt),fname_inp(1:ilenfname_inp)
303   format(/a,1x,'Current name of file is "',a,'"')

         else

            nwords=0
            do while(nwords.eq.0)
               write(nout,*)
               write(promptstring,302) prompt(1:ilenprompt),symbol_setup(1:ilensetup)
302   format(a,1x,' Specify the name of the ',a,' input file: ')
               lenpromptstring=len_trim(promptstring)
               call readline(sometext,promptstring(1:lenpromptstring))
               call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
               if(nwords.gt.0) then
                  call concatenate(1,nwords,fname_inp,ilenfname_inp
     *,nwords,ncf,ncl,sometext)
               endif
            end do

         endif
c
c check the name for the expected extension - if not present it is added on
c
c first strip out any extension that may have already been supplied
c
         ilenuptodot=index(fname_inp,'.')-1
         if(ilenuptodot.gt.0) then
            ilenfname_inp=ilenuptodot
         else
            ilenfname_inp=len_trim(fname_inp)
         endif
c
c if the expected extension has characters before the . ,then check if these are present
c - anything beyond these characters if found will be ignored
c
         if(ilenextinpdot.gt.1) then
            itest=index(fname_inp,extinp(1:(ilenextinpdot-1)))
            if(itest.gt.0) ilenfname_inp=itest-1
         endif
c
c reconstruct the input file name with the correct extension
c
         if(ilenfname_inp.gt.0) then
            fname_inp=fname_inp(1:ilenfname_inp)//extinp(1:ilenextinp)
         else
            fname_inp=extinp(1:ilenextinp)
         endif
         ilenfname_inp=ilenfname_inp+ilenextinp
         write(nout,304) prompt(1:ilenprompt),fname_inp(1:ilenfname_inp)
304   format(/a,1x,'Writing to input file "',a,'"')
c
c create the full name of the input file by adding on the directory name if not present

         call get_full_filename(fname_inp,ilenfname_inp,fnameinp,ilenfnameinp)
c
c open input file. For plot files we force it to accept the current filename
c
         if(symbol_setup.eq.'plot') then
            ierrfile=0
         else
            open(14,file=fnameinp,status='new',iostat=ierrfile)
         endif
         if(ierrfile.ne.0) then
            write(nout,301) prompt(1:ilenprompt),fname_inp(1:ilenfname_inp)
301   format(/a,1x,'File "',a,'" already exists.'//)
            nwords=0
            do while(nwords.eq.0)
               sometext=' '
               write(promptstring,305)
305   format(' Do you want to overwrite it (y or n)? ')
               lenpromptstring=len_trim(promptstring)
               call readline(sometext,promptstring(1:lenpromptstring))
               call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            end do
            if(sometext(ncf(1):ncl(1)).eq.'y') then
               ierrfile=0
            else
               ilenfname_inp=0
            endif
         endif
      end do
c
c so now go on and save the data!
c
      if(symbol_setup.eq.'changeato') then
         call save_files_changeato(14)
      else
         call save_files_others(14)
      endif

      close(14)
      return
      end

      subroutine force_save_files(nin,nout)
      use jmk_readline_module
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      integer*4 ilenuptodot
      character*256 promptstring
      integer*4 nin,nout,ierrfile,lenpromptstring
c
c first check whether the existing file name is O.K.
c
      ierrfile=1
      do while(ierrfile.ne.0)

         if(ilenfname_inp.gt.0) then

            write(nout,303) prompt(1:ilenprompt),fname_inp(1:ilenfname_inp)
303   format(/a,1x,'Current name of file is "',a,'"')

         else

            nwords=0
            do while(nwords.eq.0)
               write(nout,*)
               write(promptstring,302) prompt(1:ilenprompt),symbol_setup(1:ilensetup)
302   format(a,1x,' Specify the name of the ',a,' input file: ')
               lenpromptstring=len_trim(promptstring)
               call readline(sometext,promptstring(1:lenpromptstring))
               call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
               if(nwords.gt.0) then
                  call concatenate(1,nwords,fname_inp,ilenfname_inp
     *,nwords,ncf,ncl,sometext)
               endif
            end do

         endif
c
c check the name for the expected extension - if not present it is added on
c
c first strip out any extension that may have already been supplied
c
         ilenuptodot=index(fname_inp,'.')-1
         if(ilenuptodot.gt.0) then
            ilenfname_inp=ilenuptodot
         else
            ilenfname_inp=len_trim(fname_inp)
         endif
c
c if the expected extension has characters before the . ,then check if these are present
c - anything beyond these characters if found will be ignored
c
         if(ilenextinpdot.gt.1) then
            itest=index(fname_inp,extinp(1:(ilenextinpdot-1)))
            if(itest.gt.0) ilenfname_inp=itest-1
         endif
c
c reconstruct the input file name with the correct extension
c
         if(ilenfname_inp.gt.0) then
            fname_inp=fname_inp(1:ilenfname_inp)//extinp(1:ilenextinp)
         else
            fname_inp=extinp(1:ilenextinp)
         endif
         ilenfname_inp=ilenfname_inp+ilenextinp
         write(nout,304) prompt(1:ilenprompt),fname_inp(1:ilenfname_inp)
304   format(/a,1x,'Writing to input file "',a,'"')
c
c create the full name of the input file by adding on the directory name if not present

         call get_full_filename(fname_inp,ilenfname_inp,fnameinp,ilenfnameinp)
c
c open input file. For plot files we force it to accept the current filename
c
         ierrfile=0
      end do
c
c so now go on and save the data!
c
      if(symbol_setup.eq.'changeato') then
         call save_files_changeato(14)
      else
         call save_files_others(14)
      endif

      close(14)
      return
      end

      subroutine save_files_others(nout)
c
c writes the general input file
c
      use reallocation_routines
      use threebody_mod
      use coordconstraints_mod

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      include 'iresetnumbers.inc'
      integer*4 iqmin,npgausstest
      character*256 formatline,fullline,shortline

      piconv=4.0*atan(1.0)/180.0
c
c find the symbol reference for qmin
c
      if(symbol_setup.eq.'epsr') then
         iqmin=ireff(iqminref,nfiles)
      else
c
c ensure this value will not be used if it has not been defined
c
         iqmin=nsymbol
      endif
c
c open the new input file
c
      open(nout,file=fnameinp,status='old',iostat=ierr)
      if(ierr.ne.0) then !Try opening as an new file
         open(nout,file=fnameinp,status='new',iostat=ierr)
      end if
c Because gfortran gives an error when a file opened with status unknown does not exist!
!      ierr=0

      if(ierr.eq.0) then
c
c generate the run name - this should be the same as the input file but without
c the last .extension
c
         ilenfname_run=indexlastdot(fname_inp)-1
         fname_run=fname_inp(1:ilenfname_run)
c
c write new input file
c
         write(nout,229,iostat=ierr) fname_run(1:ilenfname_run)
229   format(a,15x,'Title of this file')
         icount=0
         ifilecount=1
         do while (icount.lt.nsymboldata)
            icount=icount+1
! For EPSR files skip over the ereq parameters which are now stored in the .pcof file
            if(symbol_setup.eq.'epsr') then
               if(icount.gt.iereqfirst.and.icount.le.iereqlast) then
                  icount=iereqlast+1
               else if(icount.gt.iqmin.and.icount.lt.nsymboluptodata) then
c for EPSR files skip over the bit related to the .pcof file - this is written
c to a separate file
                  icount=nsymboluptodata
               end if
            endif
            iref=ireff(icount,ifilecount)
            if(icount.eq.nsymboluptodata+1) then
               write(nout,204,iostat=ierr) dataprompt(1:ilendataprompt),ifilecount
204   format(/a,1x,i3/)
            endif
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            write(nout,101,iostat=ierr) symbol(icount)
     *,sometext(ncf(1):ncl(max(1,nwords))),symbol_text(icount)
101   format(a11,1x,a,15x,a)
            if(icount.ge.nsymboldata.and.ifilecount.lt.nfiles) then
               ifilecount=ifilecount+1
               icount=nsymboluptodata
            endif
         end do
c
c write a 'q' to signify the end of input (but save is not needed)
c
         write(nout,102,iostat=ierr)
102   format('q')

      endif

      if(ierr.ne.0) then

         write(6,103)
103   format(/'save_files> Could not open specified file for writing.')
         return

      endif

c
c if symbol_setup .eq. epsr then save the pcof and ato files
c
      if(symbol_setup.eq.'epsr') then

         close(nout)
         open(nout,file=fnamepcof,status='unknown',iostat=ierr)
         icount=iereqfirst+1 !Start from ereq
         do while (icount.lt.nsymboluptodata)
            iref=ireff(icount,ifilecount)
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            write(nout,101) symbol(icount),sometext(ncf(1):ncl(max(1,nwords)))
     *,symbol_text(icount)
            if(icount.eq.iereqlast) icount=iqmin !Skip to bit before data files
            icount=icount+1
         end do
         write(nout,102)
c
c now write the potential parameters
c
         icount=nsymboldata
         write(nout,*) npairspcof
         do ip=1,npairspcof
c
c write the atom symbols - these will be used to recognise this atom pair in other
c simulations
c
            write(nout,'(a3,1x,a3,1x,5x,i0,2x,f10.5)') tlab1(ip),tlab2(ip),ip,coulfac(ip)
c
c save the current values of rminex in case they have been changed
c
            icount=icount+1
            iref=ireff(icount,nfiles)
            sometext=' '
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)

c No. of Gaussians to write. We need three parameters for each Gaussian, its position
c its sigma and its height or depth. Any error in reading these and the data are ignored.

            iword=1
            read(sometext(ncf(iword):ncl(iword)),*,iostat=ierr) rminex(ip)
            if(iword.lt.nwords) then
                iword=iword+1
                read(sometext(ncf(iword):ncl(iword)),*,iostat=ierr) coordex(ip)
            else
                coordex(ip)=0.0
            end if
! Determine if any rmaxex has been defined.
            npgausstest=int((nwords-iword)/3)
            iresidual=nwords-iword-3*npgausstest
            if(iresidual.eq.1) then
               iword=iword+1
               read(sometext(ncf(iword):ncl(iword)),*,iostat=ierr) rmaxex(ip)
            else
               rmaxex(ip)=0.0
               pcoefexmax(ip)=0.0
            end if
            igauss=0
            do while(iword.lt.nwords.and.ierr.eq.0.and.igauss.lt.mpgauss)

               igauss=igauss+1
               iword=iword+1
               read(sometext(ncf(iword):ncl(iword)),*,iostat=ierr) rpgauss(igauss,ip)
               iword=iword+1
               read(sometext(ncf(iword):ncl(iword)),*,iostat=ierr) wpgauss(igauss,ip)
               iword=iword+1
               read(sometext(ncf(iword):ncl(iword)),*,iostat=ierr) hpgauss(igauss,ip)
               if(ierr.ne.0) then
                  igauss=igauss+1
               end if
            end do
            npgauss(ip)=igauss
!            write(6,*) ip,rminex(ip),pcoefex(ip)
c
c write the repulsive potential coefficients and Gaussian values if present
c
            if(npgauss(ip).eq.0) then
               if(rmaxex(ip).le.0.0) then
                  write(nout,*) rminex(ip),pcoefex(ip),coordex(ip)
               else
                  write(nout,*) rminex(ip),pcoefex(ip),coordex(ip),rmaxex(ip),pcoefexmax(ip)
               end if
            else
               if(rmaxex(ip).le.0.0) then
                  write(formatline,'("(3(1x,f10.5),",i1.1,"(1x,f9.4,1x,f9.4,1x,e11.4))")') npgauss(ip)
                  write(nout,formatline) rminex(ip),pcoefex(ip),coordex(ip)
     1,(rpgauss(i,ip),wpgauss(i,ip),hpgauss(i,ip),i=1,npgauss(ip))
               else
                  write(formatline,'("(5(1x,f10.5),",i1.1,"(1x,f9.4,1x,f9.4,1x,e11.4))")') npgauss(ip)
                  write(nout,formatline) rminex(ip),pcoefex(ip),coordex(ip),rmaxex(ip),pcoefexmax(ip)
     1,(rpgauss(i,ip),wpgauss(i,ip),hpgauss(i,ip),i=1,npgauss(ip))
               end if
            endif
c
c write the potential generation coefficients
c
            write(nout,*) (pcoefp(ic,ip),ic=1,ncoeffp)
         end do
! Write out any tb potential groups
         if(nthreeb.gt.0) then
            do ithreeb=1,nthreeb
c
c save the current values of the three-body potential in case they have been changed
c
               icount=icount+1
               iref=ireff(icount,nfiles)
               sometext=' '
               sometext=symbol_value(iref)
               call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
               if(nwords.gt.12) then
                  do i=1,3
                     read(sometext(ncf(i):ncl(i)),'(a3)') threebgrp(i,ithreeb)
                  end do
                  i=4
                  read(sometext(ncf(i):ncl(i)),*) threebamp(ithreeb)
                  i=i+1
                  read(sometext(ncf(i):ncl(i)),*) threebawidth(ithreeb)
                  threebawidth(ithreeb)=piconv*threebawidth(ithreeb)
                  i=i+1
                  read(sometext(ncf(i):ncl(i)),*) threebrwidth(ithreeb)
                  do j=1,2
                     i=i+1
                     read(sometext(ncf(i):ncl(i)),*) threebrmax(j,ithreeb)
                  end do
                  do j=1,3
                     i=i+1
                     read(sometext(ncf(i):ncl(i)),*) threebsamemol(j,ithreeb)
                  end do
                  nthreebangs(ithreeb)=(nwords-i)/2
                  if(nthreebangs(ithreeb).gt.9) nthreebangs(ithreeb)=9
                  if(nthreebangs(ithreeb).gt.nthreebangsmax) nthreebangsmax=nthreebangs(nthreeb)
                  call reallocate2d_r(threebangs,nthreebangsmax,nthreeb)
                  call reallocate2d_r(threebamps,nthreebangsmax,nthreeb)
                  do j=1,nthreebangs(ithreeb)
                     i=i+1
                     read(sometext(ncf(i):ncl(i)),*) threebangs(j,ithreeb)
                     threebangs(j,ithreeb)=piconv*threebangs(j,ithreeb)
                     i=i+1
                     read(sometext(ncf(i):ncl(i)),*) threebamps(j,ithreeb)
                  end do
               end if
               write(formatline,'("(a,3(1x,a3),1x,f10.2,2(1x,f7.2),2(1x,f10.5),3(1x,i1),",i1.1,"(1x,f7.2))")') 2*nthreebangs(ithreeb)
               write(6,'(a)') formatline(1:len_trim(formatline))
               write(nout,formatline)
     *'tb',(threebgrp(j,ithreeb),j=1,3),threebamp(ithreeb),threebawidth(ithreeb)/piconv,threebrwidth(ithreeb)
     *,(threebrmax(j,ithreeb),j=1,2),(threebsamemol(j,ithreeb),j=1,3),(threebangs(j,ithreeb)/piconv,threebamps(j,ithreeb),j=1,nthreebangs(ithreeb))
            end do
         end if
         if (ncoordconstraints.gt.0) then
            do icoord=1,ncoordconstraints
               write(formatline,'("(a,1x,a,1x,i3,2(1x,f12.3,1x,f8.3),",i2.2,"(1x,a,2(1x,f8.3)))")') coordconstraints(icoord)%ns
               write(6,'(a,1x,i5)') formatline(1:len_trim(formatline)),coordconstraints(icoord)%ns
               write(shortline,formatline) 'coordcon',coordconstraints(icoord)%atc,coordconstraints(icoord)%nc
     *,coordconstraints(icoord)%cctemp,coordconstraints(icoord)%ccdistort
     *,coordconstraints(icoord)%ccgaussweight,coordconstraints(icoord)%ccgausswidth
     *,(coordconstraints(icoord)%ats(i),coordconstraints(icoord)%ds(i),coordconstraints(icoord)%dsmax(i)
     *,i=1,coordconstraints(icoord)%ns)
               fullline=shortline
               if (coordconstraints(icoord)%hastb) then
                  write(formatline,'("(1x,a,2(1x,f12.3,1x,f8.3),",i2.2,"(1x,f8.3)))")')
     *,coordconstraints(icoord)%ntbangles
                  write(6,'(a,1x,i5)') formatline(1:len_trim(formatline)),coordconstraints(icoord)%ntbangles
                  write(shortline,formatline) 'angles',coordconstraints(icoord)%tbweight,coordconstraints(icoord)%tbdistort
     *,coordconstraints(icoord)%tbgaussweight,coordconstraints(icoord)%tbgausswidth
     *,(coordconstraints(icoord)%tbangles(i),i=1,coordconstraints(icoord)%ntbangles)
                  fullline=fullline(1:len_trim(fullline))//shortline
               end if
               if (coordconstraints(icoord)%hastbintra) then
                  write(formatline,'("(1x,a,2(1x,f12.3,1x,f8.3),",i2.2,"(1x,a3,1x,f8.3)))")')
     *,coordconstraints(icoord)%ntbintraangles
                  write(6,'(a,i12)') formatline(1:len_trim(formatline)),coordconstraints(icoord)%ntbintraangles
                  write(shortline,formatline) 'anglesintra',coordconstraints(icoord)%tbintraweight,coordconstraints(icoord)%tbintradistort
     *,coordconstraints(icoord)%tbintragaussweight,coordconstraints(icoord)%tbintragausswidth
     *,(coordconstraints(icoord)%tbintraats(i),coordconstraints(icoord)%tbintraangles(i),i=1,coordconstraints(icoord)%ntbintraangles)
                  fullline=fullline(1:len_trim(fullline))//shortline
               end if
               write(nout,'(a)') fullline
            end do
         end if
         close(nout)
c
c now write the specified .ato file
c but only write the .ato file if it contains some molecules!
c
         if(nmol.gt.0) call moloutc(fnameato)
      else
         close(nout)
      endif

      return
      end

      subroutine read_pcof(nin,nout,icount,ifilecount,ierrout)
c
c sets up the EPSR pcof values and symbols. This should always be done AFTER the
c .ato file has been read
c
      use reallocation_routines
      use threebody_mod
      use coordconstraints_mod

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      include 'iresetnumbers.inc'
      integer icount,ifilecount,ifound,ierrout,ierrread
      integer irefireset,irefiinit
      character*256 formatline,text
      character*20 ext,command
      character*1 iaccept
      integer ilenext,ilenlastdot,npgausstest,npcoeffread,ilendash
      real sumamp

      pi=4.0*atan(1.0)
      piconv=pi/180.0
c
c ensure the input stream number is different from the stream number where data are
c currently being read
c
      nin1=nin+2
! Initialise the three-body potential
      nthreeb=0
      nthreebangsmax=0
! Initialise the coordination constraints
      ncoordconstraints=0
c
c set up the integer corresponding to fnamepcof - it is assumed that symbol(icount)
c has already been checked before calling this routine
c
      ifound=ireff(icount,ifilecount)
!      write(6,*) 'read_pcof> ',icount,ifilecount,ifound
c
c get the symbol number for ireset and iinit, and read ireset. If not found
c then assume ireset = 1
c
      if(iresetref.gt.0.and.iresetref.le.nsymbol) then
         irefireset=ireff(iresetref,ifilecount)
         irefiinit=ireff(iinitref,ifilecount)
         read(symbol_value(irefireset),*) ireset
         read(symbol_value(irefiinit),*) iinit
      else
         ireset=1
         irefireset=0
         iinit=1
         irefiinit=0
      endif
c
c setup the expected extension for .pcof files
c
      ext='.pcof'
      ilenext=5
c
c open the .pcof file and read in the relevant parameters from there
c
      if(symbol_value(ifound).eq.'search'
     *.or.symbol_value(ifound).eq.'<undefined>') then
c
c search for a suitable pcof file
c
         nwords=0
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(nout,103) prompt(1:ilenprompt)
103   format(/a,1x,'No .pcof file specified! - routine will exit'/)
               if(symbol_save.ne.' ') then
                  symbol_value(ifound)=symbol_save
               else
                  symbol_value(ifound)=symbol_default(ifound)
               endif
               ierrout=0
               return
            endif
         end do
         call concatenate(1,nwords,fname_pcof,ilenfname_pcof
     *,nwords,ncf,ncl,sometext)
      else
         sometext=symbol_value(ifound)
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         call concatenate(1,nwords,fname_pcof,ilenfname_pcof
     *,nwords,ncf,ncl,sometext)
      endif

c ensure this file as the correct extension

      ilenlastdot=indexlastdot(fname_pcof)-1
      if(ilenlastdot.gt.0) then
         fname_pcof=fname_pcof(1:ilenlastdot)//ext
         ilenfname_pcof=ilenlastdot+ilenext
      else
         fname_pcof=fname_pcof(1:ilenfname_pcof)//ext
         ilenfname_pcof=ilenfname_pcof+ilenext
      endif
c
c convert this to a full path name if necessary
c
      call get_full_filename(fname_pcof,ilenfname_pcof,fnamepcof,ilenfnamepcof)
c
c open the file if possible and read the commands
c
      open(nin1,file=fnamepcof,status='old',iostat=ierr)
      if(ierr.eq.0) then
c
c reset the number of symbols to the number of data symbols
c
         nsymbol=nsymboldata
c
c get the symbol number for 'ereq' - this is now read in from pcof file, if present.
c
         icount=iereqfirst
c
c sum of all the read errors
c
         ierrread=0
c
c set up a store for the symbol number of ncoeffp
c
         incoeffpref=0
c
c initialise command
c
         command=" "
         ierr=0
         do while(command.ne.'q'.and.ierr.eq.0)
            icount=icount+1
c
c ignore lines with nothing in them
c
            nwords=0
            do while(nwords.eq.0.and.ierr.eq.0)
               read(nin1,100,iostat=ierr) sometext
100   format(a)
               ierrread=ierrread+abs(ierr)
               call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
            end do
            command=sometext(ncf(1):ncl(1))
c
c reflect the input text into the output stream
c
c           if(nwords.gt.0) write(nout,*) (sometext(ncf(i):ncl(i)),' ',i=1,nwords)

c
c don't waste time searching for a symbol which is not there for sure
c
            if(nwords.gt.1) then
c
c otherwise check the symbol number and save the new value
c
               call search_for_symbol(command,icount)
               if(icount.ne.0.and.icount.le.nsymbol) then
                  iref=ireff(icount,ifilecount)
!      write(6,*) 'read_pcof> ',icount,ifilecount,iref
                  call concatenate(2,nwords,symbol_value(iref),ilenvalue
     *,nwords,ncf,ncl,sometext)
c      write(6,*) command,' ',iref,icount,' ',symbol_value(iref)
               endif
c
c save the symbol number corresponding to ncoeffp
c
               if(command.eq.'ncoeffp') incoeffp=icount
               if(command.eq.'roverlap') iroverlap=icount
               if(command.eq.'rminfac') irminfac=icount

            endif
         end do

c
c Set up the array size for energy pair values
c
!         command='menerg'
!         call search_for_symbol(command,icount)
!         iref=ireff(icount,ifilecount)
!         read(symbol_value(iref),*) menerg

!         command='menergadd'
!         call search_for_symbol(command,icount)
!         iref=ireff(icount,ifilecount)
!         read(symbol_value(iref),*) menergadd
c
c get the number of empirical potential coefficients
c
         iref=ireff(incoeffp,ifilecount)
         read(symbol_value(iref),*,iostat=ierr) ncoeffp
c         write(6,*) ncoeffp
! Make sure this number does not exceed the allowed dimension
         if(ncoeffp.gt.mcoeff) then
            ncoeffp=mcoeff
         endif
c Write the new value to the symbol table
         write(sometext,*) ncoeffp
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         call concatenate(1,1,symbol_value(iref),ilensymbol_value,nwords,ncf,ncl,sometext)
         ierrread=ierrread+abs(ierr)
         iref=ireff(iroverlap,ifilecount)
         read(symbol_value(iref),*,iostat=ierr) roverlap
         ierrread=ierrread+abs(ierr)
         iref=ireff(irminfac,ifilecount)
         read(symbol_value(iref),*,iostat=ierr) rminfac
         ierrread=ierrread+abs(ierr)
c
c now get from this file, or create, the rminex values
c
         npairspcof=0
         read(nin1,*,iostat=ierr) npairspcof
c     write(6,*) npairspcof,ncoeffp
         ierrread=ierrread+abs(ierr)

         write(nout,104) fname_pcof(1:ilenfname_pcof)
104   format(/'.pcof file "',a,'"')

         do irefp=1,npairspcof
            efacp(irefp)=1.0 !Set the factor which controls whether this potential will be refined
            read(nin1,100,iostat=ierr) sometext
            ierrread=ierrread+abs(ierr)
            call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
            tlab1(irefp)=sometext(ncf(1):ncl(1))
            tlab2(irefp)=sometext(ncf(2):ncl(2))
            coulfac(irefp)=1.0
            if(nwords.gt.3.and.ireset.ne.1) then
               read(sometext(ncf(4):ncl(4)),*,iostat=ierr) coulfac(irefp)
            endif
c      write(6,*) irefp,npairspcof,nfiles,nwords,nsymbol
c
c generate the new symbol name
c
            nsymbol=nsymbol+1
            iref=ireff(nsymbol,nfiles)
            symbol(nsymbol)=sometext(ncf(1):ncl(1))//'-'
     *//sometext(ncf(2):ncl(2))
            iflag(nsymbol)=0
            symbol_default(nsymbol)=' '
            symbol_type(nsymbol)='r'
            read(nin1,100,iostat=ierr) sometext
            call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)

c No. of Gaussians to read. We need three parameters for each Gaussian, its position
c its sigma and its height or depth. Any error in reading these and the data are ignored.

            iword=1
! 2021-02-28 We  have added the expected coordination number, coordex(irefp) after pcoefex here
            read(sometext(ncf(iword):ncl(iword)),*,iostat=ierr) rminex(irefp)
            iword=iword+1
            read(sometext(ncf(iword):ncl(iword)),*,iostat=ierr) pcoefex(irefp)
! Make sure if the specified radius is zero or ireset = 2, the corresponding coefficient is zero
            if(rminex(irefp).eq.0.0.or.ireset.eq.2) pcoefex(irefp)=0.0
            if(iword.lt.nwords) then
                iword=iword+1
                read(sometext(ncf(iword):ncl(iword)),*,iostat=ierr) coordex(irefp)
            else
                coordex(irefp)=0.0
            end if
            rmaxex(irefp)=0.0
            pcoefexmax(irefp)=0.0
! Check whether rmaxex has been specified. If it has there will be 2+npgauss*3 words left.
! If not there will be only npgauss*3 words left.
            npgausstest=int((nwords-iword)/3)
            iwordresidual=nwords-iword-npgausstest*3
            if(iwordresidual.eq.2) then
               iword=iword+1
               read(sometext(ncf(iword):ncl(iword)),*,iostat=ierr) rmaxex(irefp)
               iword=iword+1
               read(sometext(ncf(iword):ncl(iword)),*,iostat=ierr) pcoefexmax(irefp)
               if(rmaxex(irefp).eq.0.0) pcoefexmax(irefp)=0.0
            end if
! Read Gaussians if present
            igauss=0
            do while(iword.lt.nwords.and.ierr.eq.0.and.igauss.lt.mpgauss)

               igauss=igauss+1
               iword1=iword+1
               iword2=iword+2
               iword3=iword+3
               read(sometext(ncf(iword1):ncl(iword1)),*,iostat=ierr)
     1rpgauss(igauss,irefp)
               read(sometext(ncf(iword2):ncl(iword2)),*,iostat=ierr)
     1wpgauss(igauss,irefp)
               read(sometext(ncf(iword3):ncl(iword3)),*,iostat=ierr)
     1hpgauss(igauss,irefp)

               if(ierr.ne.0) then
                  igauss=igauss-1
               endif

               iword=iword3

            end do
            npgauss(irefp)=igauss
c
c remove absurdly small values of pcoefex
c
            if(abs(pcoefex(irefp)).lt.1.0e-10) pcoefex(irefp)=0.0
            ierrread=ierrread+abs(ierr)
c
c check that the minimum allowed radius for this pair is correct
c However for Q or q atoms the minimum radius can be zero
c
            if(abs(rminex(irefp)).lt.roverlap
     *.and.sometext(ncf(1):ncf(1)).ne.'q'
     *.and.sometext(ncf(1):ncf(1)).ne.'Q'
     *.and.sometext(ncf(2):ncf(2)).ne.'q'
     *.and.sometext(ncf(2):ncf(2)).ne.'Q'

     *) rminex(irefp)=sign(roverlap,rminex(irefp))

c
c generate the appropriate symbol value
c
            sometext=' '
            if(npgauss(irefp).eq.0) then
               if(rmaxex(irefp).le.0.0) then
                  write(sometext,'(1x,f10.5,a,f10.5)',iostat=ierr) rminex(irefp),' ',coordex(irefp)
               else
                  write(sometext,'((1x,f10.5,a,f10.5,1x,f10.5))',iostat=ierr) rminex(irefp),' ',coordex(irefp),rmaxex(irefp)
               end if
            else
               if(rmaxex(irefp).le.0.0) then
                  write(formatline,'("((1x,f10.5,a,f10.5),",i1.1,"(1x,f9.4,1x,f9.4,1x,e11.4))")') npgauss(irefp)
                  write(sometext,formatline,iostat=ierr) rminex(irefp),' ',coordex(irefp)
     1,(rpgauss(igauss,irefp),wpgauss(igauss,irefp),hpgauss(igauss,irefp)
     1,igauss=1,npgauss(irefp))
               else
                  write(formatline,'("((1x,f10.5,a,f10.5,1x,f10.5),",i1.1,"(1x,f9.4,1x,f9.4,1x,e11.4))")') npgauss(irefp)
                  write(sometext,formatline,iostat=ierr) rminex(irefp),' ',coordex(irefp),rmaxex(irefp)
     1,(rpgauss(igauss,irefp),wpgauss(igauss,irefp),hpgauss(igauss,irefp)
     1,igauss=1,npgauss(irefp))
               end if
            endif
            ierrread=ierrread+abs(ierr)
            call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
            call concatenate(1,nwords,symbol_value(iref),ilensymbol
     *,nwords,ncf,ncl,sometext)
            symbol_text(nsymbol)='Minimum (and maximum) separation and P,W,H of Gaussians for this atom pair'
            write(nout,*) symbol(nsymbol)(1:10),symbol_value(iref)(1:ilensymbol)
c
c read the potential coefficients or reset them, depending on the value of ireset
c
            read(nin1,*,iostat=ierr) (pcoefp(ic,irefp),ic=1,ncoeffp) !The first one is not used by the program and will be used to signal a factor for the remainder
            ierrread=ierrread+abs(ierr)
            if(ireset.ne.0) then
! Reset them to zero if requested. The first one is not used by the program and will be used to signal a factor for the remainder
               if(ireset.lt.2) pcoefp(1,irefp)=1.0
               do ic=2,ncoeffp
                  pcoefp(ic,irefp)=0.0
               end do
            endif
            efacp(irefp)=pcoefp(1,irefp)
! Check whether any of the additional terms requires the EP to be held fixed: to achieve this the position or width needs to be set < 0
            do ig=1,npgauss(irefp)
               if(rpgauss(ig,irefp).lt.0.0.or.wpgauss(ig,irefp).lt.0.0) efacp(irefp)=0.0
            end do
            pcoefp(1,irefp)=efacp(irefp)
         end do
c
c try to match each pair of atoms in the current .ato file with one of the
c sets of coefficients read in. If not then create a new set of coefficients.
c
         do ip=1,npairs
            it=n1extr(ip)
            jt=n2extr(ip)
            irefpfound=0
            irefp=0
            do while(irefp.lt.npairspcof.and.irefpfound.eq.0)
               irefp=irefp+1
               if(
     *(tlab1(irefp).eq.type(it).and.tlab2(irefp).eq.type(jt))
     *.or.
     *(tlab1(irefp).eq.type(jt).and.tlab2(irefp).eq.type(it))
     *) irefpfound=irefp
            end do
c
c if the correct pair was found, then save it
c
            if(irefpfound.gt.0) then
               irefpcof(ip)=irefpfound
               write(6,'(a,2(1x,i5),2(1x,a))') 'read_pcof> Pair found: ',ip,irefpcof(ip),tlab1(irefpfound),tlab2(irefpfound)
            else
c
c otherwise create a new set of coefficients
c
               npairspcof=npairspcof+1
               irefpcof(ip)=npairspcof
               tlab1(npairspcof)=type(it)
               tlab2(npairspcof)=type(jt)
               pcoefp(1,npairspcof)=1.0
               do ic=2,ncoeffp
                  pcoefp(ic,npairspcof)=0.0
               end do

c Do not set the minimum distance for qtype atoms

               if(type(it)(1:1).ne.'q'.and.type(it)(1:1).ne.'Q'
     *.and.type(jt)(1:1).ne.'q'.and.type(jt)(1:1).ne.'Q') then

                  rminex(npairspcof)=rminfac*dcor(ip)
                  if(rminex(npairspcof).lt.roverlap)
     *rminex(npairspcof)=roverlap

               else

                  rminex(npairspcof)=0.0

               endif
               pcoefex(npairspcof)=0.0
               rmaxex(npairspcof)=0.0
               pcoefexmax(npairspcof)=0.0
               npgauss(npairspcof)=0
c
c generate the new symbol name
c
               write(sometext,*) type(it),' ',type(jt)
               call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
               nsymbol=nsymbol+1
               iref=ireff(nsymbol,nfiles)
               symbol(nsymbol)=sometext(ncf(1):ncl(1))//'-'//sometext(ncf(2):ncl(2))
               iflag(nsymbol)=0
               symbol_default(nsymbol)=' '
               symbol_type(nsymbol)='r'
c
c generate the appropriate symbol value
c
               sometext=' '
               write(sometext,*,iostat=ierr) rminex(npairspcof)
               ierrread=ierrread+abs(ierr)
               call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
               call concatenate(1,1,symbol_value(iref),ilensymbol,nwords,ncf,ncl,sometext)
               symbol_text(nsymbol)='Minimum separation and P,W,H of Gaussians for this atom pair'
               write(nout,*) symbol(nsymbol)(1:10),symbol_value(iref)(1:ilensymbol)
            endif

         end do
! Get any three body terms or coordination constraints that have been defined
         do while(ierr.eq.0)
            read(nin1,'(a)',iostat=ierr) sometext
c
c parse this text into words
c
            if(ierr.eq.0) then
               call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
               if(nwords.gt.0) then
!                  write(6,'(a,i3,1x,a)') 'minitc> ',nwords,sometext(ncf(1):ncl(nwords))
! Trap the program here in case a three body potential group has been read in
                  if(sometext(ncf(1):ncl(1)).eq.'tb'.and.nwords.gt.13) then
! Three body group to be read in
                     nthreeb=nthreeb+1
                     call reallocate2d_c(threebgrp,3,3,nthreeb)
                     call reallocate1d_r(threebamp,nthreeb)
                     call reallocate2d_r(threebrmax,2,nthreeb)
                     call reallocate1d_r(threebawidth,nthreeb)
                     call reallocate1d_r(threebrwidth,nthreeb)
                     call reallocate1d_i(nthreebangs,nthreeb)
                     call reallocate2d_i(threebsamemol,3,nthreeb)
                     iref=1
                     do i=1,3
                         iref=iref+1
                         read(sometext(ncf(iref):ncl(iref)),'(a3)') threebgrp(i,nthreeb)
                     end do
! Ensure third atom is same as first (for time being), otherwise this three body group is ignored
!                     if(threebgrp(1,nthreeb).eq.threebgrp(3,nthreeb)) then
                         iref=iref+1
                         read(sometext(ncf(iref):ncl(iref)),*) threebamp(nthreeb)
                         iref=iref+1
                         read(sometext(ncf(iref):ncl(iref)),*) threebawidth(nthreeb)
                         threebawidth(nthreeb)=1.0
                         iref=iref+1
                         read(sometext(ncf(iref):ncl(iref)),*) threebrwidth(nthreeb)
                         do i=1,2
                             iref=iref+1
                             read(sometext(ncf(iref):ncl(iref)),*) threebrmax(i,nthreeb)
                         end do
                         do i=1,3
                             iref=iref+1
                             read(sometext(ncf(iref):ncl(iref)),*) threebsamemol(i,nthreeb)
                         end do
                         nthreebangs(nthreeb)=(nwords-iref)/2
                         write(6,*) 'read_pcof> No. of tb angles: ',nthreeb,nthreebangs(nthreeb)
                         if(nthreebangs(nthreeb).gt.9) nthreebangs(nthreeb)=9
                         if(nthreebangs(nthreeb).gt.nthreebangsmax) nthreebangsmax=nthreebangs(nthreeb)
                         call reallocate2d_r(threebangs,nthreebangsmax,nthreeb)
                         call reallocate2d_r(threebamps,nthreebangsmax,nthreeb)
                         sumamp=0.0
                         do i=1,nthreebangs(nthreeb)
                             iref=iref+1
                             read(sometext(ncf(iref):ncl(iref)),*) threebangs(i,nthreeb)
                             threebangs(i,nthreeb)=cos(piconv*threebangs(i,nthreeb))
                             iref=iref+1
                             read(sometext(ncf(iref):ncl(iref)),*) threebamps(i,nthreeb)
                             threebamps(i,nthreeb)=abs(threebamps(i,nthreeb))
                             sumamp=sumamp+threebamps(i,nthreeb)
                         end do
                         !Normalise individual terms to total sum of amplitudes
                         if(sumamp.gt.0.0) then
!                            do i=1,nthreebangs(nthreeb)
! 2016-02-09 Cut this out for compatibility with earlier simulations
!                               threebamps(i,nthreeb)=threebamps(i,nthreeb)/sumamp
!                            end do
                         else
                            threebamp(nthreeb)=0.0
                         end if
!                         threebawidth(nthreeb)=threebawidth(nthreeb)*threebangs(1,nthreeb)/pi
!                         threebangs(1,nthreeb)=cos(threebangs(1,nthreeb))+threebawidth(nthreeb)*log(2.0)
!                     else
!                         nthreeb=nthreeb-1
!                     end if
c
c generate the new symbol name
c
                     nsymbol=nsymbol+1
                     iref=ireff(nsymbol,nfiles)
                     symbol(nsymbol)=sometext(ncf(1):ncl(1))
                     iflag(nsymbol)=0
                     symbol_default(nsymbol)=' '
                     symbol_type(nsymbol)='c'
                     call concatenate(2,nwords,symbol_value(iref),ilensymbol,nwords,ncf,ncl,sometext)
                     symbol_text(nsymbol)='Three atom labels, amplitude, width, angle width, two dist., and angles'
                     write(nout,*) symbol(nsymbol)(1:10),symbol_value(iref)(1:ilensymbol)
                  else if (sometext(ncf(1):ncl(1)).eq.'coordcon'.and.nwords.gt.8) then
                     ncoordconstraints=ncoordconstraints+1
                     call reallocate1d_cc(coordconstraints,ncoordconstraints)
                     text=sometext(ncf(2):ncl(nwords))
      write(6,*) 'read_pcof> No. coord. constraints: ',ncoordconstraints,' ',text(1:len_trim(text))
                     call init_constraint(coordconstraints(ncoordconstraints),text)
!                  call find_constraints(coordconstraints(ncoordconstraints),nmol,natom,xadd,yadd,zadd,xa,ya,za,mref
!     1,itype,ntype,type)
                  end if
               end if
            end if
         end do
c
c if any error occurred during input, then ireset is set to 1
c
         if(ierrread.ne.0) then
            ireset=1
            iinit=1
            ierrout=ierrread
         else
            ierrout=0
         endif
         close(nin1)
      else
         write(6,101) prompt(1:ilenprompt),fname_pcof(1:ilenfname_pcof)
101   format(/a,1x,'Problems with specified .pcof file: ',a
     */' - will use defaults or existing values')
         write(6,*) npairspcof
c
c get the default parameters needed to set up the potential coefficients
c
         iref=ireff(incoeffp,ifilecount)
         read(symbol_value(iref),*,iostat=ierr) ncoeffp

c Make sure this number does not exceed the allowed dimension

         if(ncoeffp.gt.mcoeff) ncoeffp=mcoeff
         iref=ireff(iroverlap,ifilecount)
         read(symbol_value(iref),*,iostat=ierr) roverlap
         iref=ireff(irminfac,ifilecount)
         read(symbol_value(iref),*,iostat=ierr) rminfac
c
c the value of npairspcof will be used to signify whether a set of coefficients
c has already been defined. If so we leave well alone, but keep the new name
c for the .pcof file so this will be written when the save is done.
c
         if(npairspcof.le.0) then
c
c set ireset to 1 to ensure the .pcof file and other parameters are initialised correctly
c on saving the input files and running the simulation
c
            ireset=1
            iinit=1
            ierrread=1
         else
            ierrout=0
            ierrread=0
            iref=ireff(nsymboldata,nfiles)
            do ip=1,npairspcof
               iref=iref+1
               write(6,*) rminex(ip),' ',symbol_value(iref)(1:20)
            end do
         endif

      endif
c
c if there was an error on reading the .pcof file we generate default values
c
      if(ierrread.ne.0.or.ireset.eq.1) then
c
c set the number of symbols to the current maximum
c
         nsymbol=nsymboldata
c
c set up minimum overlap distances for each atom pair. These are based on the values of
c the L-J sigma for each atom, but can be altered as required.
c
         npairspcof=npairs
         write(nout,104) fname_pcof(1:ilenfname_pcof)
         do ip=1,npairs
            it=n1extr(ip)
            jt=n2extr(ip)

c Do not set the minimum distance for qtype atoms

            if(type(it)(1:1).ne.'q'.and.type(it)(1:1).ne.'Q'
     *.and.type(jt)(1:1).ne.'q'.and.type(jt)(1:1).ne.'Q') then

               rminex(ip)=rminfac*dcor(ip)
               if(rminex(ip).lt.roverlap) rminex(ip)=roverlap

            else

               rminex(ip)=0.0

            endif
            pcoefex(ip)=0.0
            coordex(ip)=0.0
            rmaxex(ip)=0.0
            pcoefexmax(ip)=0.0
            irefpcof(ip)=ip
            coulfac(irefpcof(ip))=1.0
            tlab1(ip)=type(it)
            tlab2(ip)=type(jt)
            npgauss(ip)=0
c
c generate the new symbol name
c
            write(sometext,*) tlab1(ip),' ',tlab2(ip)
            call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
            nsymbol=nsymbol+1
            iref=ireff(nsymbol,nfiles)
            symbol(nsymbol)=sometext(ncf(1):ncl(1))//'-'//sometext(ncf(2):ncl(2))
            iflag(nsymbol)=0
            symbol_default(nsymbol)=' '
c
c generate the appropriate symbol value
c
            sometext=' '
            write(sometext,*) rminex(ip),' ',coordex(ip)
            call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
            call concatenate(1,1,symbol_value(iref),ilensymbol,nwords,ncf,ncl,sometext)
            symbol_text(nsymbol)='Minimum (coordination number and maximum) separation and P,W,H of Gaussians for this atom pair'
            write(nout,*) symbol(nsymbol)(1:10),symbol_value(iref)(1:ilensymbol)
c
c Initialize the potential coefficients
c
            pcoefp(1,ip)=1.0
            do ic=2,ncoeffp
               pcoefp(ic,ip)=0.0
            end do
         end do

         ierrout=0

      endif

      symbol_value(ifound)=fname_pcof(1:ilenfname_pcof)
c
c write the new values of ireset and iinit in case they have been changed
c
      if(iresetref.ne.0) then
         write(sometext,*) ireset
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         call concatenate(1,1,symbol_value(irefireset),ilensymbol_value,nwords,ncf,ncl,sometext)
      endif
      if(iinitref.ne.0) then
         write(sometext,*) iinit
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         call concatenate(1,1,symbol_value(irefiinit),ilensymbol_value,nwords,ncf,ncl,sometext)
      endif

      return
      end

      subroutine read_ato(nin,nout,icount,ifilecount,ierrout)
c
c reads the EPSR .ato file
c
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      integer*4 icount,ifilecount,ifound,ierrout
      integer*4 isizefactor,irho
      character*20 ext,text
      character*1 iaccept
      integer*4 ilenext
c
c ensure the input stream number is different from the stream number where data are
c currently being read
c
      nin1=nin+2
c
c set up the integer corresponding to fnameato - it is assumed that symbol(icount)
c has already been checked before calling this routine
c
      ifound=ireff(icount,ifilecount)
c
c get the symbol number for sizefactor and read sizefactor. If not found then assume
c sizefactor is 1.0
c
      text='sizefactor'
      isizefactor=0
      call search_for_symbol(text,isizefactor)
      if(isizefactor.gt.0) then
         iref=ireff(isizefactor,ifilecount)
         read(symbol_value(iref),*) sizefactor
      else
         sizefactor=1.0
      endif

c set up the expected extension of the .ato file

      ext='.ato'
      ilenext=4

c open the .ato file and read in the relevant parameters from there

      if(symbol_value(ifound).eq.'search'
     *.or.symbol_value(ifound).eq.'<undefined>') then

c search for a suitable ato file

         nwords=0
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(nout,103) prompt(1:ilenprompt)
103   format(/a,1x,'No .ato file specified! - routine will exit'/)
               if(symbol_save.ne.' ') then
                  symbol_value(ifound)=symbol_save
               else
                  symbol_value(ifound)=symbol_default(ifound)
               endif
               ierrout=0
               return
            endif
         end do
         call concatenate(1,nwords,fname_ato,ilenfname_ato
     *,nwords,ncf,ncl,sometext)

      else

         sometext=symbol_value(ifound)
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         call concatenate(1,nwords,fname_ato,ilenfname_ato
     *,nwords,ncf,ncl,sometext)

      endif

c ensure this file as the correct extension

      ilenlastdot=indexlastdot(fname_ato)-1
      if(ilenlastdot.gt.0) then
         fname_ato=fname_ato(1:ilenlastdot)//ext
         ilenfname_ato=ilenlastdot+ilenext
      else
         fname_ato=fname_ato(1:ilenfname_ato)//ext
         ilenfname_ato=ilenfname_ato+ilenext
      endif

c
c convert this to a full path name if necessary
c
      call get_full_filename(fname_ato,ilenfname_ato,fnameato,ilenfnameato)

      write(nout,202) fnameato(1:ilenfnameato)
202   format(/'Attempting to read file: ',a/)
c
c now attempt to read the specified .ato file
c
c
c initialise number of atom and molecule types to zero, after saving them
c
      ntypesave=ntype
      nmoltypesave=nmoltype
      ntype=0
      nmoltype=0
      call minitc(fnameato,rmaxcell,cellmax,upol,rho,amubar,amubar2
     *,ndist,ierrout)

      if(ierrout.ne.0) then
c
c the file did not exist or had errors, but we will leave any existing
c atom parameters as they are and write them on saving to the new file
c
         write(nout,31)
31    format('No atoms were read from this file')
c
c restore the number of atom and molecule types
c
         ntype=ntypesave
         nmoltype=nmoltypesave
         ierrout=0
      else
c
c write out the list of atom types in this file
c
!         write(nout,32) ntype
!32    format(/'There are ',i2,' types of atom in this file'/)
!         do it=1,ntype
!            write(nout,33) it,type(it)
!33    format(1x,'Atom type ',i2,' has label ',a3)
!         end do
      endif

c save the current fnameato symbol_value

      symbol_value(ifound)=fname_ato(1:ilenfname_ato)
c
c update the value of rho and save it
c
      irho=0
      text='rho'
      call search_for_symbol(text,irho)
      if(irho.ne.0) then
         iref=ireff(irho,ifilecount)
         write(symbol_value(iref),*) rho
      endif

      return
      end


