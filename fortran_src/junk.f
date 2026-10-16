      subroutine get_system_commands()
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'system_commands.inc'
c     write(6,*) 'Got to get_system_commands'
c
c initialise the system commands
c
      system_ls=' '
      ilensystem_ls=0
      system_cd=' '
      ilensystem_cd=0
      system_md=' '
      ilensystem_md=0
      system_edit=' '
      ilensystem_edit=0
      system_del=' '
      ilensystem_del=0
c
c termination character(s) for the end of a directory name
c This will be added to the end of directory names if it is not present
c
      system_termination=' '
      ilensystem_termination=0
c
c characters to join commands in one line
c
      system_join=' '
      ilensystem_join=0
c
c command to invoke gnuplot
c
      system_gnu=' '
      ilensystem_gnu=0

c Command to invoke jmol

      system_jmol=' '
      ilensystem_jmol=0

c System binaries location

      system_binaries=' '
      ilensystem_binaries=0
c
c commands to do 2d and 3d plotting
c
      system_plot3d=' '
      ilensystem_plot3d=0
      system_plot2d=' '
      ilensystem_plot2d=0
      system_plotato=' '
      ilensystem_plotato=0
      system_pgout=' '
      ilensystem_pgout=0
      system_quotes=' '
      ilensystem_quotes=0
      gnu_pointtype=0
      gnu_pointsize=0
c
c open the system_command file in the current directory
c If a fault occurs, the Windows Command Prompt commands are assumed
c
      open(10,file='system_commands.txt',status='old',iostat=ierr)
      do while (ierr.eq.0)
         sometext=' '
         read(10,100,iostat=ierr) sometext
100   format(a)
c
c parse this text to see if contains any words
c
           call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
c     write(6,*) nwords,(' ',sometext(ncf(i):ncl(i)),i=1,nwords)
         if(nwords.gt.1) then
            if(sometext(ncf(1):ncl(1)).eq.'system_ls') then
               nlast=nwords
               call concatenate(2,nlast,system_ls,ilensystem_ls
     *,nwords,ncf,ncl,sometext)
c              write(6,102) ilensystem_ls,system_ls(1:ilensystem_ls)
102   format(1x,i3,1x,a)
            else if(sometext(ncf(1):ncl(1)).eq.'system_cd') then
               nlast=nwords
               call concatenate(2,nlast,system_cd,ilensystem_cd
     *,nwords,ncf,ncl,sometext)
            else if(sometext(ncf(1):ncl(1)).eq.'system_md') then
               nlast=nwords
               call concatenate(2,nlast,system_md,ilensystem_md
     *,nwords,ncf,ncl,sometext)
            else if(sometext(ncf(1):ncl(1)).eq.'system_edit') then
               nlast=nwords
               call concatenate(2,nlast,system_edit,ilensystem_edit
     *,nwords,ncf,ncl,sometext)
c              write(6,102) ilensystem_edit,system_edit(1:ilensystem_edit)
            else if(sometext(ncf(1):ncl(1)).eq.'system_del') then
               nlast=nwords
               call concatenate(2,nlast,system_del,ilensystem_del
     *,nwords,ncf,ncl,sometext)
            else if(sometext(ncf(1):ncl(1)).eq.'system_termination') then
               nlast=nwords
               call concatenate(2,nlast,system_termination
     *,ilensystem_termination,nwords,ncf,ncl,sometext)
c              write(6,102) ilensystem_termination
c     *,system_termination(1:ilensystem_termination)
            else if(sometext(ncf(1):ncl(1)).eq.'system_join') then
               nlast=nwords
               call concatenate(2,nlast,system_join
     *,ilensystem_join,nwords,ncf,ncl,sometext)
            else if(sometext(ncf(1):ncl(1)).eq.'system_gnu') then
               nlast=nwords
               call concatenate(2,nlast,system_gnu
     *,ilensystem_gnu,nwords,ncf,ncl,sometext)
            else if(sometext(ncf(1):ncl(1)).eq.'system_jmol') then
               nlast=nwords
               call concatenate(2,nlast,system_jmol
     *,ilensystem_jmol,nwords,ncf,ncl,sometext)
            else if(sometext(ncf(1):ncl(1)).eq.'system_binaries') then
               nlast=nwords
               call concatenate(2,nlast,system_binaries
     *,ilensystem_binaries,nwords,ncf,ncl,sometext)
            else if(sometext(ncf(1):ncl(1)).eq.'system_plot3d') then
               nlast=nwords
               call concatenate(2,nlast,system_plot3d
     *,ilensystem_plot3d,nwords,ncf,ncl,sometext)
            else if(sometext(ncf(1):ncl(1)).eq.'system_plot2d') then
               nlast=nwords
               call concatenate(2,nlast,system_plot2d
     *,ilensystem_plot2d,nwords,ncf,ncl,sometext)
            else if(sometext(ncf(1):ncl(1)).eq.'system_plotato') then
               nlast=nwords
               call concatenate(2,nlast,system_plotato
     *,ilensystem_plotato,nwords,ncf,ncl,sometext)
            else if(sometext(ncf(1):ncl(1)).eq.'system_pgout') then
               nlast=nwords
               call concatenate(2,nlast,system_pgout
     *,ilensystem_pgout,nwords,ncf,ncl,sometext)
            else if(sometext(ncf(1):ncl(1)).eq.'system_quotes') then
               nlast=nwords
               call concatenate(2,nlast,system_quotes
     *,ilensystem_quotes,nwords,ncf,ncl,sometext)
            else if(sometext(ncf(1):ncl(1)).eq.'gnu_pointtype') then
               read(sometext(ncf(2):ncl(nwords)),*,iostat=ierr) gnu_pointtype
            else if(sometext(ncf(1):ncl(1)).eq.'gnu_pointsize') then
               read(sometext(ncf(2):ncl(nwords)),*,iostat=ierr) gnu_pointsize
            endif
         endif
      end do
      close(10)
      if(gnu_pointtype.eq.0) then
         gnu_pointtype=8
      endif
      if(gnu_pointsize.eq.0.0) then
         gnu_pointsize=0.3
      endif
      return
      end


