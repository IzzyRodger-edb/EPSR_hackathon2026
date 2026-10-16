      subroutine get_command_from_args()

c Get the command into sometext from the command line arguments

      use String_Utility
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'version.inc'

      character*256 text
      integer*4 ilensometext

c Get the current version

      version=date
      ilenversion=ilendate
c
c get the system commands for the particular operating system in use
c
      call get_system_commands()

c Add on the path separator to the system_binaries.

      system_binaries=system_binaries(1:ilensystem_binaries)
     *//system_termination(1:ilensystem_termination)
      ilensystem_binaries=ilensystem_binaries+ilensystem_termination

      write(6,600) system_binaries(1:ilensystem_binaries)
600   format(/'Binaries folder: ',a)
c
c set initial directory name - this is called the home directory
c
      call get_directory_name()
      ncharhomedir=nchardir
      home_directory=directory_name

c Get input from any supplied arguments

      narg=iargc()
      iarg=0
      sometext=' '
      ilensometext=0
      do while (iarg.lt.narg)

         iarg=iarg+1
         call getarg(iarg,text)

c Parse the text

         call parse(text,nwords,ncf,ncl,5,nstsize,nworddim)

c Set the current directory - should be the first argument

         if(iarg.eq.1) then

            directory_name = text(ncf(1):ncl(nwords))
            nchardir = ncl(nwords)-ncf(1)+1
            write(6,100) directory_name(1:nchardir)
100   format(/'Current folder for input files is: ',a)

         else if(iarg.eq.2) then

! The second argument will be normally be symbol_setup: these are always lower case.

            sometext=text(ncf(1):ncl(nwords))
            ilensometext=ncl(nwords)-ncf(1)+1
            symbol_setup = StrLowCase(sometext(1:ilensometext))
            ilensetup=ilensometext

         else

c Otherwise accumulate the arguments into a single text line separated by spaces

            sometext=sometext(1:ilensometext)//' '//text(ncf(1):ncl(nwords))
            ilensometext=ilensometext + ncl(nwords)-ncf(1)+2

         endif

      end do

      write(6,102) sometext(1:ilensometext)
102   format(a)

c Read the command

      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      write(6,*) nwords

      return
      end

      subroutine get_directory_name()
c
c sets the directory where input will be read from.
c
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      ierr=getcwd(sometext)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)

c concatenate the words into the directory name

      nlast=nwords
      call concatenate(1,nlast,directory_name,nchardir
     *,nwords,ncf,ncl,sometext)

c To cope with spaces in folder names we need to insert the directory name between quotes.
c These cannot include the first (if present) and last system termination character.

c Check whether the first character is a system termination character. 

c      if(directory_name(1:ilensystem_termination).eq.system_termination(1:ilensystem_termination)) then

c         directory_name=system_termination(1:ilensystem_termination)
c     *//'"'
c     *//directory_name(ilensystem_termination+1:nchardir)

c      else

c         directory_name='"'//directory_name(ilensystem_termination+1:nchardir)

c      endif
c      nchardir=nchardir+1

c check that the last characters are the correct termination sequence. If not add it.

      ntest=nchardir-ilensystem_termination+1
      if(directory_name(ntest:nchardir)
     *.ne.system_termination(1:ilensystem_termination)) then
         directory_name=directory_name(1:nchardir)
c    *   //'"'
     *   //system_termination(1:ilensystem_termination)
         nchardir=nchardir+ilensystem_termination
      else
         directory_name=directory_name(1:ntest-1)
c     *   //'"'
     *   //system_termination(1:ilensystem_termination)
      endif
c      nchardir=nchardir+1
   
      write(6,100) directory_name(1:nchardir)
100   format(/'Home folder is: ',a)

      return
      end

      subroutine set_directory_name()
c
c sets the current directory
c
c if no directory is given the current directory is shown
c
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
c Temporary character store
      character*256 temp_name,stringcheck,dir_temp,dir_test
      character*20 testseq
      integer*4 nchartemp,ilentestseq,ipos,ilenstringcheck,nchardirtemp,nchardirtest
      logical fileexists
c
c Store the current directory in a temporary variable
    
      dir_temp=directory_name(1:nchardir)
      nchardirtemp=nchardir
c
c check that a name was specified - if not do nothing
c
      if(nwords.gt.1) then
         if(sometext(ncf(2):ncl(2)).eq.'home'
     *.or.sometext(ncf(2):ncl(2)).eq.'.') then

c If the second word is 'home'or '.' 
c the directory name is reset to the home directory
c (which is the directory where EPSRshell was first started)

            directory_name=home_directory
            nchardir=ncharhomedir
         
         else if(sometext(ncf(2):ncl(2)).eq.'..') then

c If the second word is '..' the directory name goes one level back up the tree

c Working from the end of the directory name search for the first occurrence of 
c a system termination character (but ignoring the very last one). If none are found,
c do nothing.

            ntest=nchardir-ilensystem_termination
            ipos=ntest+ilensystem_termination-1
            do while(directory_name(ntest:ipos)
     *.ne.system_termination(1:ilensystem_termination).and.ntest.gt.0)
               ntest=ntest-1
               ipos=ipos-1
            end do
            if(ntest.gt.0) then
               nchardir=ipos
               directory_name=directory_name(1:nchardir)
            endif

         else

            nlast=nwords

c Place the remaining words into a store

            call concatenate(2,nlast,temp_name,nchartemp
     *,nwords,ncf,ncl,sometext)

c Search for system_termination characters

            icount=0
            testseq=system_termination(1:ilensystem_termination)
            ilentestseq=ilensystem_termination
            ifirst=1
            iadd=index(temp_name(ifirst:nchartemp),testseq(1:ilentestseq))

            do while (iadd.gt.0)

               ilast=ifirst+iadd-2

c Check what is the string before this. If ilast.lt.ifirst there is nothing before this
c so start a new directory name.

c Otherwise check whether the characters before the system_termination are '.' or '..'

               if(ilast.ge.ifirst) then

                  stringcheck=temp_name(ifirst:ilast)
                  ilenstringcheck=ilast-ifirst+1

                  if(temp_name(ifirst:ilast).eq.'.') then

c Single '.' means backup to home directory

                     directory_name=home_directory
                     nchardir=ncharhomedir
         
                  else if(temp_name(ifirst:ilast).eq.'..') then

c Double '..' means backup one directory

c Working from the end of the directory name search for the first occurrence of 
c a system termination character (but ignoring the very last one). If none are found,
c do nothing.

                     ntest=nchardir-ilensystem_termination
                     ipos=ntest+ilensystem_termination-1
                     do while(directory_name(ntest:ipos)
     *.ne.system_termination(1:ilensystem_termination).and.ntest.gt.0)
                        ntest=ntest-1
                        ipos=ipos-1
                     end do
                     if(ntest.gt.0) then
                        nchardir=ipos
                        directory_name=directory_name(1:nchardir)
                     endif

                  else

c Otherwise add the specified string to the existing directory name

                     directory_name=directory_name(1:nchardir)
     *//stringcheck(1:ilenstringcheck)//system_termination(1:ilensystem_termination)
                     nchardir=nchardir+ilenstringcheck+ilensystem_termination

                  endif

               else if(ifirst.eq.1) then

c In the event that we encounter a system_termination string at the beginnning, 
c we start a new directory name

                  directory_name=system_termination(1:ilensystem_termination)
                  nchardir=ilensystem_termination

               else

c Do nothing

               endif

               icount=icount+1
               ifirst=ilast+ilensystem_termination+1 
c     write(6,*) icount,ifirst,ilast

               if(ifirst.le.nchartemp) then

                  iadd=index(temp_name(ifirst:nchartemp),testseq(1:ilentestseq))

c If there is still text in the word, we need one more pass to process it

                  if(iadd.eq.0) iadd=nchartemp-ifirst+2

               else

                  iadd=0

               endif

            end do

c In case there is still text in the word, we add it to the current directory name

            if(ifirst.le.nchartemp) then

               directory_name=directory_name(1:nchardir)
     *//temp_name(ifirst:nchartemp)//system_termination(1:ilensystem_termination)
               nchardir=nchardir+nchartemp-ifirst+1+ilensystem_termination

            endif

         endif

c check that the last characters are the correct termination sequence. If not add it.

         ntest=nchardir-ilensystem_termination+1
         if(directory_name(ntest:nchardir)
     *.ne.system_termination(1:ilensystem_termination)) then
            directory_name=directory_name(1:nchardir)
     *//system_termination(1:ilensystem_termination)
            nchardir=nchardir+ilensystem_termination
         endif
      endif

c Now check that this folder exists

      dir_test=directory_name(1:nchardir)//'.'
      nchardirtest=nchardir+1
c Check whether this file exists
      inquire(file=dir_test(1:nchardirtest),exist=fileexists)
c      write(6,*) dir_test(1:nchardirtest),' ',fileexists
      if(.not.fileexists) then
         write(6,101) directory_name(1:nchardir)
101   format(/'No such folder! : ',a)
         directory_name=dir_temp(1:nchardirtemp)
         nchardir=nchardirtemp
      endif

      write(6,100) directory_name(1:nchardir)
100   format(/'Current folder for input files: ',a)

      return
      end

      subroutine list_directory()
c
c shows the files in the current input directory
c
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      character*256 command,text
      integer*4 nchartext
      if(nwords.gt.1) then
c
c concatenate the remaining words into the filename list
c
         nlast=nwords
         call concatenate(2,nlast,text,nchartext
     *,nwords,ncf,ncl,sometext)
         command=system_ls(1:ilensystem_ls)
     *//' '//system_quotes(1:ilensystem_quotes)
     *//directory_name(1:nchardir)//text(1:nchartext)
     *//system_quotes(1:ilensystem_quotes)

      else

         command=system_ls(1:ilensystem_ls)
     *//' '//system_quotes(1:ilensystem_quotes)
     *//directory_name(1:nchardir)
     *//system_quotes(1:ilensystem_quotes)

      endif
      call system(command,ierr)
      return
      end

      subroutine make_directory()
c
c shows the files in the current input directory
c
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      character*256 command,text,fnamefull
      integer*4 nchartext,ilenfnamefull
      if(nwords.gt.1) then
c
c concatenate the remaining words into the filename list
c
         nlast=nwords
         call concatenate(2,nlast,text,nchartext
     *,nwords,ncf,ncl,sometext)
      else
         write(6,100)
100   format(/'make_directory> Type the directory name to create: ',$)
         read(5,101) sometext
101   format(a)
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         call concatenate(1,nwords,text,nchartext
     *,nwords,ncf,ncl,sometext)
      endif

      call get_full_filename(text,nchartext,fnamefull,ilenfnamefull)
      command=system_md(1:ilensystem_md)
     *//' '
     *//system_quotes(1:ilensystem_quotes)
     *//fnamefull(1:ilenfnamefull)
     *//system_quotes(1:ilensystem_quotes)

      call system(command,ierr)

      return
      end

      subroutine edit_file()
c
c shows the files in the current input directory
c
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      character*256 command,text,fnamefull
      integer*4 nchartext,ilenfnamefull
      if(nwords.gt.1) then
c
c concatenate the remaining words into the filename list
c
         nlast=nwords
         call concatenate(2,nlast,text,nchartext
     *,nwords,ncf,ncl,sometext)
         call get_full_filename(text,nchartext,fnamefull,ilenfnamefull)
      else
         fnamefull=' '
         ilenfnamefull=0
      endif

      write(6,100) fnamefull(1:ilenfnamefull)
100   format(/'edit_file> ',a)

      if(ilenfnamefull.gt.0) then
         command=system_edit(1:ilensystem_edit)
     *//' '
     *//system_quotes(1:ilensystem_quotes)
     *//fnamefull(1:ilenfnamefull)
     *//system_quotes(1:ilensystem_quotes)

         ilencommand=ilensystem_edit+ilenfnamefull+1
     *+2*ilensystem_quotes

      else

         command=system_edit(1:ilensystem_edit)
         ilencommand=ilensystem_edit

      endif

      write(6,100) command(1:ilencommand)

      call system(command,ierr)

      return
      end

      subroutine delete_file()
c
c shows the files in the current input directory
c
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      character*256 command,text,fnamefull
      integer*4 nchartext,ilenfnamefull
      if(nwords.gt.1) then
c
c concatenate the remaining words into the filename list
c
         nlast=nwords
         call concatenate(2,nlast,text,nchartext
     *,nwords,ncf,ncl,sometext)
         call get_full_filename(text,nchartext,fnamefull,ilenfnamefull)
         command=system_del(1:ilensystem_del)
     *//' '
     *//system_quotes(1:ilensystem_quotes)
     *//fnamefull(1:ilenfnamefull)
     *//system_quotes(1:ilensystem_quotes)

         call system(command,ierr)

      else

         write(6,100)
100   format(/'delete_file> No files specified to delete. ')

      endif

      return
      end

      subroutine get_full_filename(fnameshort,ilenshort,fnamelong,ilenlong)
c
c checks to see whether a directory path is present in fnameshort, and if not
c adds it. Output in fnamelong.
c
      include 'directory_name.inc'
      include 'system_commands.inc'
      character*256 fnameshort,fnamelong
      integer*4 ilenshort,ilenlong,ilentermination
c
c search for directory termination characters in fnameshort
c
      ilentermination=index(fnameshort
     *,system_termination(1:ilensystem_termination))

      if(nchardir.gt.0.and.ilentermination.eq.0) then
         fnamelong=directory_name(1:nchardir)
     *//fnameshort(1:ilenshort)
         ilenlong=nchardir+ilenshort
      else
         fnamelong=fnameshort(1:ilenshort)
         ilenlong=ilenshort
      endif

      return
      end

      subroutine get_nterms(nin,fname,extin,nterms,ierr)
c
c returns the number of columns in a file. Lines containing a # are ignored.
c It is assumed the maximum number of columns in a file is 201 (EPSR or Gudrun
c files). If this number is reached the nextfile in the sequence is opened if
c possible
c
      include 'commands.inc'
      character*256 fname
      character*2 ext
      character*10 extin
      integer*4 ilenextin,nterms

      ilenextin=index(extin,' ')-1
      write(6,*) fname,' ',extin(1:ilenextin)
      ntest=100
      nterms=0
      ifilecount=0
      ierr=0
      do while(ntest.ge.100.and.ierr.eq.0)
         ifilecount=ifilecount+1
c
c generate the version number
c
         write(ext,100) ifilecount
100   format(i2.2)
c
c strip any existing extension off the input filename
c
         ilenfname=indexlastdot(fname)-1
         if(ilenfname.le.0) ilenfname=index(fname,' ')-1
c
c add the expected extension back on
c
         fname=fname(1:ilenfname)//extin(1:ilenextin)
         ilenfname=ilenfname+ilenextin
c
c strip the last two characters off the filename *assume it is of the form
c xxxxx.y01, etc.) and create the sequence of versions in turn
c
         ilenfname=ilenfname-2
         fname=fname(1:ilenfname)//ext(1:2)
c
      write(6,*) fname(1:ilenfname+2)
c
c try to open this file
c
         open(nin,file=fname,status='old',iostat=ierr)
         if(ierr.eq.0) then
c
c read lines of text until one occurs without a # in it
c
            nhash=1
            ierr=0
            do while(nhash.ne.0.and.ierr.eq.0)
               read(nin,101,iostat=ierr) sometext
101   format(a)
               nhash=index(sometext,'#')
            end do
c
c find out how many words this line contains
c
            if(nhash.eq.0) then
               call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
c
c number of terms = (ncolumns-1)/2
c
               ntest=(nwords-1)/2
               nterms=nterms+ntest
            endif
            close(nin)
         endif
      end do
      return
      end

      subroutine setup_input_file(nin,nout,ierrout)
c
c sets up the EPSR input files
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
c     include 'molecom.inc'
c     include 'clustercom.inc'
c input and output stream number
      integer*4 nin,nout
      integer*4 ierrout,ilenext,ilenuptodot
      character*1 iaccept
      character*20 ext
      character*256 text
c
c generate the extension from the last dot onwards
c
      text=extinp
      ilenextinpdot=indexlastdot(text)
      extdot=extinp(ilenextinpdot:ilenextinp)
      ilenextdot=ilenextinp-ilenextinpdot+1
c     write(6,*) extdot(1:ilenextdot),ilenextinp,ilenextinpdot,ilenextdot
c
c write out the value of the input class and symbol_setup to ensure we are
c indeed accessing the class of file requested
c
      write(nout,101) symbol_setup(1:ilensetup),extinp(1:ilenextinp)
101   format(/'setup_input_file> File class: "',a
     *,'"; file extension: "',a,'"')
c
c parse the supplied text
c
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
c
c generate the name of the input file
c
      iw=1
      do while(iw.le.nwords.and.sometext(ncf(iw):ncl(iw))
     *.ne.symbol_setup(1:ilensetup))
         iw=iw+1
      end do
      if(iw.ge.nwords) then
c
c no filename found so prompt for it from the user
c
         ext=extinp
         ilenext=ilenextinp
c     write(6,*) ext(1:ilenext)
         nwords=0
         sometext=' '
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(nout,100)
100   format(/'setup input file> No input file specified - routine will exit'/)
               ierrout=1
               return
            endif
         end do
         iw=0
      endif
      iw=iw+1
      call concatenate(iw,nwords,fname_inp,ilenfname_inp
     *,nwords,ncf,ncl,sometext)
c
c if this filename has a .ext we do nothing and will attempt to read it as is,
c warts and all!
c
      text=fname_inp
      ilenuptodot=index(fname_inp,'.')-1
c     write(6,*) fname_inp(1:ilenfname_inp),ilenfname_inp,ilenuptodot
      if(ilenuptodot.lt.0) then
c
c .ext not supplied - therefore we need to add on the correct one
c
         ilenuptodot=ilenfname_inp
         fname_inp=fname_inp(1:ilenuptodot)//extinp(1:ilenextinp)
         ilenfname_inp=ilenuptodot+ilenextinp
      endif
c
c create the full name of the input file by adding on the directory name if not present
c
      call get_full_filename(fname_inp,ilenfname_inp,fnameinp,ilenfnameinp)
      write(6,105) fnameinp(1:ilenfnameinp)
105   format(/'setup_input_file> Full filename = ',a)
c
c open this file if possible
c
      open(nin,file=fnameinp,status='old',iostat=ierrout)
      if(ierrout.eq.0) then
c
c  If O.K. read in the data from the specified input file
c
            write(nout,102) fname_inp(1:ilenfname_inp)
102   format(/'setup_input_file> Reading input file: "',a,'"')
c
c find the number of characters in this filename up to the last dot
c
            text=fname_inp
            ilenuptodot=indexlastdot(text)-1
c
c read the name of the run - this should correspond to fname_inp but without the 
c last .extension. It will be used to construct the names of input
c and output files (except .ato and .pcof files which may be named independently
c from the name of the run)
c
            read(nin,*,iostat=ierrout) fname_run
            ilenfname_run=len_trim(fname_run)
c
c we don't do the following two checks for plot input files
c
c
c remove the .extension if supplied
c
            if(symbol_setup(1:4).ne.'plot') then
c
c check that this name corresponds to filename (without the .extension) specified. 
c If not signal the error but carry on, using the correct name when writing the file
c
               if(fname_run(1:ilenfname_run)
     *.ne.fname_inp(1:ilenuptodot)) then
                  write(nout,104) 
104   format(/'setup_input_file> Run name in input file is different from '
     *,'filename specified')
                  fname_run=fname_inp(1:ilenuptodot)
                  ilenfname_run=ilenuptodot
               endif
            endif
c
c generate the long form of fname_run
c
            call get_full_filename(fname_run,ilenfname_run
     *,fnamerun,ilenfnamerun)
c
c signal that if an error occurred the attempt to open the input file
c was unsuccessful, but we can still read commands from the terminal
c
            if(ierrout.ne.0) ierrout=-1
            return
      else
            fname_run=fname_inp(1:ilenuptodot)
            ilenfname_run=ilenuptodot
            write(nout,103) fname_inp(1:ilenfname_inp)
103   format(/'setup_input_file> Problems with specified input file: ',a
     *,//'   - will use default values')
c
c generate the long form of fname_run
c
            call get_full_filename(fname_run,ilenfname_run
     *,fnamerun,ilenfnamerun)
c
c signify that the input file could not be opened or read, but a run
c name has been successfully defined with the default values
c
            ierrout=-1
            return
      endif
      end

      subroutine search_for_files(ext,ilenext,iaccept)
c
c list any files with the required extension in the current directory
c
      include 'commands.inc'
      include 'system_commands.inc'
      include 'directory_name.inc'
      character*256 fname
      character(len=*) ext
      character(len=*) iaccept
      integer*4 ilenext,ilenstar,ilenextin
      integer*4 ifirstchar,idirpath
c
c first set up the command to do a directory listing into a file
c
      sometext=' '
c      write(6,*) sometext(1:len_trim(sometext))
      sometext=system_ls(1:ilensystem_ls)//' '
     *//system_quotes(1:ilensystem_quotes)
     *//directory_name(1:nchardir)//'*'//ext(1:ilenext)
     *//system_quotes(1:ilensystem_quotes)
     *//' > dirlist.txt'
c      write(6,*) sometext(1:len_trim(sometext))
c
c get the directory listing into a file, then asks the user to select a file
c if found then iaccept is set to 'y'. If none found or none specified then iaccept is set to 'n'
c
      call system(sometext,ierr)
c
c open this file and list it line by line
c
      if(ierr.eq.0) then
            open(10,file='dirlist.txt',status='old',iostat=ierr)
      endif
c
c strip out the * in ext if present
c
      ilenstar=index(ext,'*')
      if(ilenstar.gt.0) then
         ilenextin=ilenstar-1
      else
         ilenextin=ilenext
      endif
      iaccept=' '
      ifound=0
c
c number of lines read
c
      nline=0
      do while (ierr.eq.0.and.iaccept.ne.'y'.and.iaccept.ne.'e')
         sometext=' '
         nwords=0
         read(10,101,iostat=ierr) sometext
101   format(a)
         nline=nline+1
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         if(nwords.gt.0) then
            ifound=ifound+1
            call concatenate(1,nwords,fname,ilenfname
     *,nwords,ncf,ncl,sometext)
c           write(6,*) ifound,nwords

c Strip out everything up to the last directory path character sequence

            idirpath=1
            ifirstchar=1
            do while(idirpath.gt.0.and.ifirstchar.le.ilenfname)
               idirpath=index(fname(ifirstchar:ilenfname)
     *,system_termination(1:ilensystem_termination))
               if(idirpath.gt.0) ifirstchar=ifirstchar+idirpath
     *+ilensystem_termination-1
            end do
            fname=fname(ifirstchar:ilenfname)
            ilenfname=ilenfname-ifirstchar+1
            write(6,123) fname(1:ilenfname)
c
c save the filename in case it is accepted
c
            sometext=fname(1:ilenfname)
123   format(/'Filename: ',a,'? (Type y to accept, u to go back, e to exit) ',$)
            read(5,101) iaccept
!      write(6,*) 'search_for_files> ',iaccept(1:len_trim(iaccept))
            if(iaccept.eq.'u') then
c     write(6,*) nline
c
c request to backspace through file
c
               backspace(10)
               backspace(10)
               nline=nline-2
c
c If nline is .lt. 0 it means we have attempted to skip past the beginning of file.
c Therefore read to end of file, counting the number of lines and reposition before last 
c line
c
               if(nline.lt.0) then
c
c ensure file is positioned at beginning
c
                  rewind(10)
                  nline=0
                  ierr=0
                  do while (ierr.eq.0)
                     sometext=' '
                     nwords=0
                     read(10,101,iostat=ierr) sometext
                     nline=nline+1
                     call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
                  end do
c     write(6,*) ierr,nline,nwords,sometext(ncf(1):ncl(1))
                  backspace(10)
                  nline=nline-1
c
c If the line that caused ierr to be non-zero contains no words, we need to backspace
c again.
c
                  if(nwords.le.0) then
                     backspace(10)
                     nline=nline-1
                  endif
                  ierr=0
               endif
            endif
         endif
         if(ifound.gt.0.and.ierr.ne.0) then
            rewind(10)
            ifound=0
            ierr=0
            nline=0
         endif
      end do
      close(10)
      if(ifound.eq.0) then
         write(6,125) ext(1:ilenext),directory_name(1:nchardir)
125   format(/1x,'No files of extension "',a,'" found in directory "',a,'"'/)
         iaccept='n'
      endif
      if(iaccept.ne.'y') then
         write(6,124)
124   format(/1x,'No files selected...'/
     */1x,'Type the required filename (without extension): ',$)
         sometext=' '
         read(5,101,iostat=ierr) sometext
c
c strip out any leading spaces
c
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         if(nwords.gt.0) then
            iaccept='y'
            call concatenate(1,nwords,fname,ilenfname
     *,nwords,ncf,ncl,sometext)
            sometext=fname(1:ilenfname)
         else
            iaccept='n'
         endif          
      endif
      return
      end

      subroutine get_a_line(nin,nspcword
     *,nchartext,nfirst,nlast,ifilecount)
c
c gets a line of text from the input file and parses it into words
c ncf and ncl contain the first and last character numbers of each word
c nchartext is the maximum number of characters in a sometext,
c nworddim is dimension of the
c word array, nspcword is the number of spaces between words needed to decide that the 
c end of text has been reached.
c
c The words are read into the symbols if they exist, otherwise the defaults are adopted
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      integer*4 nspcword,nchartext
      integer*4 nfirst,nlast,ifilecount
c
c make sure allowed limits are not exceeded
c
      if(nfirst.lt.1) nfirst=1
      if(nlast.gt.nsymbol) nlast=nsymbol
c
c read a line of text
c
      read(nin,100,iostat=ierr) sometext
100   format(a)
c
c parse it into words
c
      call parse(sometext,nwords,ncf,ncl,nspcword,nchartext
     *,nworddim)
      i=0
      icount=nfirst-1
      do while(icount.lt.nlast)
         i=i+1
         icount=icount+1
c
c set up the reference value for this symbol
c
         iref=ireff(icount,ifilecount)
         if(i.le.nwords) then
            symbol_value(iref)=sometext(ncf(i):ncl(i))
         else
            symbol_value(iref)=symbol_default(icount)
         endif
      end do
      return
      end

      function ireff(icount,ifilecount)
      include 'epsr_dimension.inc'
      include 'symbols.inc'
c
c calculates the reference integer for the symbol table
c
c
c set up the reference value for this symbol
c
      if(icount.gt.nsymboldata) then
         ireff=icount+(nfiles-1)*nsymbolplus
!         ifilecount=nfiles
      else if(icount.gt.nsymboluptodata) then
         ireff=icount+(ifilecount-1)*nsymbolplus
      else
         ireff=icount
!         ifilecount=0
      endif
      return
      end

      subroutine search_for_symbol(text,ifound)
c
c returns the symbol number corresponding to a specified text string
c
      include 'epsr_dimension.inc'
      include 'symbols.inc'
      integer*4 icount,ifound
      character*20 text

c
c first check that the current value of ifound does not already correspond
c to the specified symbol
c
c     write(6,101) text
101   format(/'search_for_symbol> searching for symbol ',a)

      if(ifound.gt.0.and.ifound.le.nsymbol) then
         if(symbol(ifound).eq.text) return
      endif
      ifound=0
      icount=0
      do while(icount.lt.nsymbol.and.ifound.eq.0)
         icount=icount+1
         if(symbol(icount).eq.text) ifound=icount
      end do
      if(ifound.eq.0) then
c        write(6,100) text
100   format(/'search_for_symbol> text ',a,' not found in symbol table')
      else
c        write(6,102) symbol(ifound),ifound
102   format(/'search_for_symbol> found symbol ',a,' number ',i3)
      endif

      return
      end

