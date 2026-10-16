      subroutine plot_new()
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      integer*4 ierrout
c
c set up plot symbols
c
      call setup_symbols_plot()
c
c now attempt to open the file containing the default plot values
c
      sometext='plot '//extinp(1:ilenextinp)
      nin=11
      nout=6
      call setup_input_file(nin,nout,ierrout)
c
c now read it
c
      if(ierrout.eq.0) call read_plot_input_file(nin,nout,ierrout)
c
c ierrout .lt. 0 signifies the file could not be opened or on reading encountered an 
c error, but an input file has been successfully defined using the default values,
c so these values are now available for viewing and changing.
c

      call process_plot_symbols(5,6,ierrout,0)
      return
      end

      subroutine setup_symbols_plot()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
      include 'plot_pars.inc'
c
c sets the symbols and their default values for input to torangles.
c
c iflag decides whether this parameter can be changed
c a parameter can only be changed if ireset >= iflag
c

      iref=0

      iref=iref+1
      symbol(iref)='l'
      symbol_type(iref)='c'
      symbol_text(iref)='Lists available plot types'
      iflag(iref)=0
      symbol_default(iref)='0'

      iref=iref+1
      symbol(iref)='p'
      symbol_type(iref)='i'
      symbol_text(iref)='Plot using the current or specified plot type'
      iflag(iref)=0
      symbol_default(iref)=' '
      icountp=iref

      iref=iref+1
      symbol(iref)='npt'
      symbol_type(iref)='i'
      symbol_text(iref)='Number of types of plot'
      iflag(iref)=0
      symbol_default(iref)='1'

      nsymboluptodata=iref

      iref=iref+1
      symbol(iref)='type'
      symbol_type(iref)='c'
      symbol_text(iref)='Type of plot'
      iflag(iref)=1
      symbol_default(iref)= '<undefined>'

      iref=iref+1
      symbol(iref)='ext'
      symbol_type(iref)='c'
      symbol_text(iref)='File extension(s) used to search for plot files'
      iflag(iref)=0
      symbol_default(iref)='<undefined>'

      iref=iref+1
      symbol(iref)='f'
      symbol_type(iref)='c'
      symbol_text(iref)='File name to plot'
      iflag(iref)=0
      symbol_default(iref)='<undefined>'
      icountf=iref

      iref=iref+1
      symbol(iref)='b'
      symbol_type(iref)='c'
      symbol_text(iref)='Block numbers to plot (e.g. 1 2 - 5 9 - 6)'
      iflag(iref)=0
      symbol_default(iref)='<undefined>'
      icountb=iref

      iref=iref+1
      symbol(iref)='bspace'
      symbol_type(iref)='i'
      symbol_text(iref)='Spacing between blocks in plot file'
      iflag(iref)=0
      symbol_default(iref)='2'

      iref=iref+1
      symbol(iref)='boffset'
      symbol_type(iref)='i'
      symbol_text(iref)='Column number(s) for first block(s) in plot file(s)'
      iflag(iref)=0
      symbol_default(iref)='2 2'

      iref=iref+1
      symbol(iref)='xmin'
      symbol_type(iref)='r'
      symbol_text(iref)='Minimum value on x-axis'
      iflag(iref)=0
      symbol_default(iref)='0'

      iref=iref+1
      symbol(iref)='xmax'
      symbol_type(iref)='r'
      symbol_text(iref)='Maximum value on x-axis'
      iflag(iref)=0
      symbol_default(iref)='1'

      iref=iref+1
      symbol(iref)='ymin'
      symbol_type(iref)='r'
      symbol_text(iref)='Minimum value on y-axis'
      iflag(iref)=0
      symbol_default(iref)='0'

      iref=iref+1
      symbol(iref)='ymax'
      symbol_type(iref)='r'
      symbol_text(iref)='Maximum value on y-axis'
      iflag(iref)=0
      symbol_default(iref)='1'

      iref=iref+1
      symbol(iref)='ydel'
      symbol_type(iref)='r'
      symbol_text(iref)='Spacing between plots'
      iflag(iref)=0
      symbol_default(iref)='1'

      iref=iref+1
      symbol(iref)='yfact'
      symbol_type(iref)='r'
      symbol_text(iref)='Factor(s) outside y values'
      iflag(iref)=0
      symbol_default(iref)='1.0 1.0'

      iref=iref+1
      symbol(iref)='xf'
      symbol_type(iref)='r'
      symbol_text(iref)='Fractional x coordinate of labels'
      iflag(iref)=0
      symbol_default(iref)='0.75'

      iref=iref+1
      symbol(iref)='yf'
      symbol_type(iref)='r'
      symbol_text(iref)='Fractional y coordinate of labels'
      iflag(iref)=0
      symbol_default(iref)='0.25'

      iref=iref+1
      symbol(iref)='title'
      symbol_type(iref)='c'
      symbol_text(iref)='Title of plot'
      iflag(iref)=0
      symbol_default(iref)='<undefined>'

      iref=iref+1
      symbol(iref)='xlabel'
      symbol_type(iref)='c'
      symbol_text(iref)='x-axis label'
      iflag(iref)=0
      symbol_default(iref)='<undefined>'

      iref=iref+1
      symbol(iref)='ylabel'
      symbol_type(iref)='c'
      symbol_text(iref)='y-axis label'
      iflag(iref)=0
      symbol_default(iref)='<undefined>'

      iref=iref+1
      symbol(iref)='logx'
      symbol_type(iref)='c'
      symbol_text(iref)='log scale x (yes or no)'
      iflag(iref)=0
      symbol_default(iref)='no'

      iref=iref+1
      symbol(iref)='logy'
      symbol_type(iref)='c'
      symbol_text(iref)='log scale y (yes or no)'
      iflag(iref)=0
      symbol_default(iref)='no'

      iref=iref+1
      symbol(iref)='xcolumn'
      symbol_type(iref)='i'
      symbol_text(iref)='column number for x values [normally 1 or 0]'
      iflag(iref)=0
      symbol_default(iref)='1'

      iref=iref+1
      symbol(iref)='histogram'
      symbol_type(iref)='c'
      symbol_text(iref)='histogram data (yes or no)'
      iflag(iref)=0
      symbol_default(iref)='no'

      nsymboldata=iref
      nsymbolplus=nsymboldata-nsymboluptodata

      nsymbol=iref
c
c initialise the symbol values
c
      do icount=1,nsymbol
         symbol_value(icount)=symbol_default(icount)
      end do
c
c initialise number of plot types
c
      nfiles=1
c
c setup the output prompt
c
      prompt='plot>'
      ilenprompt=index(prompt,'>')
c
c set up the data prompt
c
      dataprompt='pt'
      dataprompttext='Sets the specified plot type '
      ilendataprompt=index(dataprompt,' ')-1
c
c setup the extension for the input file
c
      extinp='plot_defaults.txt'
      ilenextinp=index(extinp,' ')-1

      return
      end

      subroutine process_plots(nin,nout,icount,ifilecount,ierr)
c
c plots a specified set of results from EPSR or its associated routines
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'plot_pars.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      character*20 command
      character*256 input_text
      character*1 yorn
      integer*4 ierr,icount,ifilecount,icountwr,ifile

c
c first determine the requested command
c
      command=symbol(icount)
c
c save the associated command text
c
      iref=ireff(icount,ifilecount)
      input_text=symbol_value(iref)
c     write(6,*) iref,' ',input_text(1:50)
      if(command.eq.'l') then
c
c save the reference integer
c
         irefl=iref
         write(nout,100) prompt(1:ilenprompt)
100   format(/a,1x,'Available plot types are:-'//)
c
c set up the reference integer to get the symbol_values
c
         icountwr=nsymboluptodata+1
         ifile=0
         do while(ifile.lt.nfiles)
            ifile=ifile+1
            iref=ireff(icountwr,ifile)
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
c
c attempt to read the first word as a type number
c
            if(nwords.gt.0) then
               read(sometext(ncf(1):ncl(1)),*,iostat=ierr) iplottype
               if(ierr.ne.0) then
c
c if it isn't set a plot type number
c
                  iplottype=ifile
                  ibegin=1
               else if(nwords.gt.1) then
                  ibegin=2
               else
                  ibegin=0
               endif
            else
               ibegin=0
            end if
            if(ibegin.gt.0) then
               write(nout,101) iplottype
     *,sometext(ncf(ibegin):ncl(max(ibegin,nwords)))
101   format(1x,i3,'  ',a)
            else
               write(nout,102) iplottype
            endif
102   format(1x,i3)
         end do
         symbol_value(irefl)=' '
      else if(command.eq.'f') then
         call setup_plot_filenames(nout,ierr)
         if(ierr.ne.0) then
            write(nout,104) prompt(1:ilenprompt)
104   format(/a,1x,' Try again (y or n)? ',$)
            read(nin,*) yorn
            if(yorn.eq.'y') then
               symbol_value(icountf)='search'
               call setup_plot_filenames(nout,ierr)
            endif
         endif
      else if(command.eq.'b') then
         call setup_blocknumbers(nin,nout,ierr)
      else if(command.eq.'p') then
c
c perform the plot. The plot type is determined from the value of
c symbol_value(icountp)
c
         read(symbol_value(icountp),*) plottype
         ifilecount=plottype
         icount=nsymboluptodata+1
         iref=ireff(icount,ifilecount)
         symbol_save=symbol_value(iref)
         call gnuplot(nin,nout,ierr)
      else
c
c if none of the above we set the plottype to ifilecount
c
         if(ifilecount.gt.0.and.ifilecount.le.nfiles) then
            plottype=ifilecount
         else
            plottype=1
         endif
         write(symbol_value(icountp),*) plottype
      endif
      return
      end

      subroutine setup_plot_filenames(nout,ierr)
c
c sets up the filenames for plotting and determines the number
c of available columns
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'plot_pars.inc'
      integer*4 ierr
      character*1 iaccept
c
c set the number of plots to zero - this will be used to signal
c whether valid plotting files have been defined
c
      nplot=0
c
c get the active plot type
c
      read(symbol_value(icountp),*,iostat=ierr) plottype
      if(ierr.ne.0.or.plottype.lt.1.or.plottype.gt.nfiles) then
         write(nout,102) 'setup_plot_filenames>',1,nfiles
102   format(/a,1x,'Specify a valid plot type between'
     *,i2,' and ',i2)
         ierr=1
         return
      endif
c
c get the extensions associated with this plottype
c
      itest=nsymboluptodata+2
      iref=ireff(itest,plottype)
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.gt.0) then
         ext1=sometext(ncf(1):ncl(1))
         ilenext1=ncl(1)-ncf(1)+1
         nplot=1
         if(nwords.gt.1) then
            ext2=sometext(ncf(2):ncl(2))
            ilenext2=ncl(2)-ncf(2)+1
            nplot=2
         endif
      else
         write(nout,103) 'setup_plot_filenames>',plottype
103   format(/a,1x,'No extensions supplied for plottype ',i3)
         ierr=1
         nplot=0
         return
      endif
c
c get the values of bspace and boffset - spacing between blocks in file
c and offset for first block
c
      iref=iref+1
      read(symbol_value(iref),*,iostat=ierr1) bspace
      iref=iref+1

c With this version we introduce the possibility of having different offsets
c for each dataset when two are plotted together. We have to ensure this does
c not crash older versions of plot_defaults.txt which will only have one entry
c for boffset

      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.gt.1) then
         read(symbol_value(iref),*,iostat=ierr2) boffset,boffset2
      else
         read(symbol_value(iref),*,iostat=ierr2) boffset
         boffset2=boffset
      endif
      if(ierr1.ne.0.or.ierr2.ne.0) then
         write(nout,104) 'setup_plot_filenames>'
104   format(/a,1x,'Block offset and spacing not defined correctly')
         nplot=0
         ierr=1
         return
      endif
c
c identify a file to be plotted, ignoring any characters after the first
c dot.
c
      sometext=symbol_value(icountf)
      ilenuptodot=index(sometext,'.')-1
      if(ilenuptodot.gt.0) sometext=sometext(1:ilenuptodot)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.lt.1.or.symbol_value(icountf).eq.'search'
     *.or.symbol_value(icountf).eq.'<undefined>') then
         call search_for_files(ext1,ilenext1,iaccept)
c
c remove any leading or trailing spaces in the filename
c
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      else
c
c Add the correct extension back in.
c
         sometext=sometext(ncf(1):ncl(nwords))//ext1
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      endif
c
c regenerate the filename into a single word
c
      call concatenate(1,nwords,fnameplot,ilenfnameplot
     *,nwords,ncf,ncl,sometext)
c
c build the full path plot filename(s)
c
      call get_full_filename(fnameplot,ilenfnameplot
     *,fname1,ilenfname1)

c
c strip out anything beyond the extension that has been supplied
c
      ilenuptoext=index(fnameplot,ext1)-1
c
c open the requested files and determine the number of columns
c
      call find_ncolumn(fname1)
      if(nblock.eq.0) then
         write(nout,125)
125   format(/'setup_plot_filenames> This file is empty!')
         ierr=1
         nplot=0
         return
      else
         write(nout,126) nblock,fname1(1:ilenfname1)
126   format(/'setup_plot_filenames> There are ',i3,' blocks in the file ',a)
      endif
      ierr=0
c
c save the good value of fnameplot
c
      symbol_value(icountf)=fnameplot(1:ilenuptoext)
c
c generate extension for corresponding data file
c
      if(nplot.gt.1) then
c
c add on the new extension
c
         fnameplot=fnameplot(1:ilenuptoext)//ext2(1:ilenext2)
         ilenfnameplot=ilenuptoext+ilenext2
c
c create the full file name
c
         call get_full_filename(fnameplot,ilenfnameplot
     *,fname2,ilenfname2)
         write(nout,127) 'setup_plot_filenames> ',fname2(1:ilenfname2)
127   format(/'setup_plot_filenames> ',a)
      endif
      return
      end

      subroutine find_ncolumn(fname)
c
c this routine applies strictly to EPSR generated files and is used by the plotting
c command. A more general purpose routine is get_nterms.
c
      include 'epsr_dimension.inc'
      include 'symbols.inc'
      include 'commands.inc'
      include 'plot_pars.inc'
      character*256 fname
      character*20 text
      integer*4 bdiv,boff
c
c get the value of xcolumn
c
      text='xcolumn'
      call search_for_symbol(text,ifound)
      iref=ireff(ifound,plottype)
      read(symbol_value(iref),*) xcolumn

c     write(6,*) fname
      ilenfname=index(fname,' ')-1
      ifilefind=0
      open(10,file=fname,status='old',iostat=ierr)
      if(ierr.ne.0) then
            write(6,*)
            write(6,*) 'Problems with file ',fname(1:ilenfname)
            nblock=0
            return
      endif
c read a line of text
      read(10,100,iostat=ierr) sometext
100   format(a)
      if(ierr.eq.0) then
c find out how many words it contains
         call parse(sometext,nwords,ncf,ncl,30,nstsize,nworddim)
         ntitle=0
c
c if the first word has a # in it then the following words are 
c the x axis title, and n column titles
c
         ilent=0
         ilent=index(sometext(ncf(1):ncl(1)),'#')
         if(ilent.gt.0) then
c
c A label, if it exists at all, will appear at the head of each column
c or pair of columns
c
c
c ignore the first word which contains the # and second word which contains the
c x-axis label (if xcolumn.gt.0)
c
            if(xcolumn.gt.0) then
               boff=3
            else
               boff=3
            endif
            iw=boff-1
            do while (iw.lt.nwords)
               iw=iw+1
               ntitle=ntitle+1
               plabel(ntitle)=sometext(ncf(iw):ncl(iw))
            end do
c
c now read another line of text and use this to determine the number of blocks
c
            read(10,100,iostat=ierr) sometext
            call parse(sometext,nwords,ncf,ncl,30,nstsize,nworddim)
         endif
c
c set up the block numbers
c
         iw=boffset
         nblock=0
         do while(iw.le.nwords)
            nblock=nblock+1
c
c if ntitle is still zero the block numbers are used to invent suitable titles
c
            if(ntitle.lt.nblock) then
               write(plabel(nblock),*) 'Block',nblock
               ntitle=nblock
            endif
            iw=iw+bspace
         end do
c
c the ratio ntitle/nblock gives the number of titles per block. If this ratio
c is greater than 1, it means there are more titles than we need, so we need
c to pick out the appropriate ones for the present value of boffset
c
         write(6,125) ntitle,nblock
125   format(/'find_ncolumn> ',2(1x,i5))
         if(nblock.gt.0) then
            bdiv=ntitle/nblock
            if(bdiv.gt.1) then
               boff=boffset-xcolumn
               it=boff
               ib=1
               do while(it.le.ntitle)
                  plabel(ib)=plabel(it)
                  it=it+bdiv
                  ib=ib+1
               end do
            endif
         endif
      else
         nblock=0
      endif
      close(10)
126   format(/'find_ncolumn> ',6(1x,i5))
      write(6,126) plottype,nwords,boffset,boffset2,bspace,nblock
      return
      end

      subroutine setup_blocknumbers(nin,nout,ierr)

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'plot_pars.inc'
c
c first check that a valid filename has been defined - any error and the 
c routine exits
c
      call setup_plot_filenames(nout,ierr)
      if(ierr.ne.0) then
         write(nout,103)
103   format(/'setup_blocknumbers> A valid filename must be defined before '
     *,'block numbers can be set.')
         return
      endif
c
c set up the plot columns
c
c     write(6,*) icountb,' ',symbol_value(icountb)(1:50)
      sometext=symbol_value(icountb)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.le.0) then
         write(nout,100) nblock
100   format(/'setup_blocknumbers> There are ',i3,' blocks. Specify which ones to plot'
     *//'(Use spaces between numbers, start with a number,'
     */' and use n1 - n2 to define a range): ',$)
         read(nin,101) sometext
101   format(a)
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         do while(nwords.eq.0)
            write(nout,102) prompt(1:ilenprompt)
102   format(/a,1x,'No column numbers specified - try again: ',$)
            read(nin,101) sometext
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         end do
      endif
c
c From the input list determine the block numbers to be plotted.
c The first word should always be a number - therefore
c read the words until a number is found
c
      iw=1
      npcolumn=0
      read(sometext(ncf(iw):ncl(iw)),*,iostat=ierr) nfirst
      iw=iw+1
c
c now step through the remaining words and find out the other column numbers
c
      do while(ierr.eq.0.and.iw.le.nwords)
         if(sometext(ncf(iw):ncl(iw)).eq.'-') then
c
c a range has been specified
c
            iw=iw+1
            read(sometext(ncf(iw):ncl(iw)),*,iostat=ierr) nlast
            if(nlast.gt.nfirst) then
               idel=1
            else if(nlast.lt.nfirst) then
               idel=-1
            else
               idel=0
            endif
            do while(nfirst.ne.nlast)
               if(nfirst.gt.0.and.nfirst.le.nblock) then
                  npcolumn=npcolumn+1
                  icolumn(npcolumn)=nfirst
               endif
               nfirst=nfirst+idel
            end do
         else
            if(nfirst.gt.0.and.nfirst.le.nblock) then
               npcolumn=npcolumn+1
               icolumn(npcolumn)=nfirst
            endif
            read(sometext(ncf(iw):ncl(iw)),*,iostat=ierr) nfirst
         endif
         iw=iw+1
      end do
c
c save the last value
c
      if(nfirst.gt.0.and.nfirst.le.nblock) then
         npcolumn=npcolumn+1
         icolumn(npcolumn)=nfirst
      endif

      if(npcolumn.gt.0) then

         write(6,235) npcolumn
235   format(/'setup_blocknumbers> Number of plotting columns: ',i5)

      else

         ierr=1
         write(6,236)
236   format(/'setup_blocknumbers> Need to define valid block numbers for this file')

      endif

c
c if the command is successful, save it
c
      if(ierr.eq.0) then
         call concatenate(1,nwords,symbol_value(icountb),ilen
     *,nwords,ncf,ncl,sometext)
      endif

      return
      end

      subroutine gnuplot(nin,nout,ierr)
c
c sets up a GNUplot input file and run GNU plot
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'plot_pars.inc'
      include 'system_commands.inc'
      character*256 plot_output,plot_command
      character*2 typelabel
      real*4 labelx,labely
c
c now check that block numbers have been defined - any error and the 
c routine exits
c
      call setup_blocknumbers(nin,nout,ierr)
      if(ierr.ne.0) then
         write(nout,122)
122   format(/'GNUplot> Block numbers need be defined before '
     *,'data can be plot.')
         return
      endif
c
c get the other information needed to form the plot
c
c
c get the extensions associated with this plottype
c
      itest=nsymboluptodata+5
      iref=ireff(itest,plottype)
      read(symbol_value(iref),*) xmin
      iref=iref+1
      read(symbol_value(iref),*) xmax
      iref=iref+1
      read(symbol_value(iref),*) ymin
      iref=iref+1
      read(symbol_value(iref),*) ymax
      iref=iref+1
      read(symbol_value(iref),*) ydel

c factor for y values

      iref=iref+1
      yfact1=1.0
      yfact2=1.0
      read(symbol_value(iref),*,iostat=ierr) yfact1,yfact2
      write(6,*) symbol_value(iref)(1:20), yfact1,yfact2
c
c get fractional positions of labels
c
      iref=iref+1
      read(symbol_value(iref),*) xfrac
      iref=iref+1
      read(symbol_value(iref),*) yfrac
c
c plot title
c
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.gt.0) then
         call concatenate(1,nwords,title,ilentitle
     *,nwords,ncf,ncl,sometext)
      else
         title=' '
         ilentitle=1
      endif
c
c set up a default title in case it is undefined
c
      if(title(1:ilentitle).eq.'<undefined>') then
         ilenuptodot=index(fnameplot,'.')-1
         title=directory_name(1:nchardir)
     *//fnameplot(1:ilenuptodot)
         ilentitle=nchardir+ilenuptodot
      endif
c
c x-axis label
c
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.gt.0) then
         call concatenate(1,nwords,xlabel,ilenxlabel
     *,nwords,ncf,ncl,sometext)
      else
         xlabel=' '
         ilenxlabel=1
      endif
c
c y-axis label
c
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.gt.0) then
         call concatenate(1,nwords,ylabel,ilenylabel
     *,nwords,ncf,ncl,sometext)
      else
         ylabel=' '
         ilenylabel=1
      endif
c
c log x-axis
c
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.gt.0) then
         call concatenate(1,1,logx,ilenlogx
     *,nwords,ncf,ncl,sometext)
      else
         logx='no'
         ilenlogx=2
      endif
c
c if this axis is log x, check that the minimum x value is > 0
c
      if(logx(1:1).eq.'y'.and.xmin.le.0.0) xmin=1.0e-12

c
c log y-axis
c
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.gt.0) then
         call concatenate(1,1,logy,ilenlogy
     *,nwords,ncf,ncl,sometext)
      else
         logy='no'
         ilenlogy=2
      endif
c
c if this axis is log y, check that the minimum y value is > 0
c
      if(logy(1:1).eq.'y'.and.ymin.le.0.0) ymin=1.0e-12
      write(6,*) logx(1:1),logy(1:1)
c
c column number for x values. 0 means use the item number as x value
c
      iref=iref+1
      read(symbol_value(iref),*) xcolumn
c
c decide if this is a histogram or not
c
      iref=iref+1
      sometext=symbol_value(iref)
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      if(nwords.gt.0) then
         call concatenate(1,1,histogram,ilenhistogram
     *,nwords,ncf,ncl,sometext)
      else
         histogram='no'
         ilenhistogram=2
      endif
c
c set the dimensions of plot
c
      yup=ymax+(npcolumn-1)*ydel
c
c generate the type label
c
      write(typelabel,150) plottype
150   format(i2.2)
c
c generate the name of the output file
c
c
c first strip out the extension in fnameplot
c
c
c strip out anything beyond the . in fnameplot
c
      ilenuptodot=index(fnameplot,'.')-1
      plot_output=directory_name(1:nchardir)//'plot'
     *//typelabel(1:2)//fnameplot(1:ilenuptodot)//'.gnu'
      ilenplot_output=nchardir+ilenuptodot+10
      open(11,file=plot_output,status='unknown',iostat=ierr)
c
c write the title and x and y labels
c
         write(11,101) title(1:ilentitle)
101   format('reset',/,'set title ',1h',a,1h')
         write(11,102) xlabel(1:ilenxlabel)
102   format('set xlabel ',1h',a,1h')
         write(11,103) ylabel(1:ilenylabel)
103   format('set ylabel ',1h',a,1h')
         write(11,104) xmin,xmax
104   format('set xrange [',e10.3,':',e10.3,']')
         write(11,105) ymin,yup
105   format('set yrange [',e10.3,':',e10.3,']')
c
c log or linear axes
c
      if(logx(1:1).eq.'y') then
         write(11,113) 'x'
113   format('set logscale ',a1)
      else
         write(11,114) 'x'
114   format('set nologscale ',a1)
      endif
      if(logy(1:1).eq.'y') then
         write(11,113) 'y'
      else
         write(11,114) 'y'
      endif

c
c set the labels and linestyles
c
c x coordinate of data label
c
      labelx=xmin+xfrac*(xmax-xmin)
      do ip=1,npcolumn
c
c y coordinate of data label
c
         labely=yfrac+(ip-1)*ydel
         iplot=icolumn(ip)
         ilent=index(plabel(iplot),' ')-1
         write(11,106) ip,plabel(iplot)(1:29),labelx,labely
106   format('set label ',i3,1x,1h',a,1h',' at ',e10.3,',',e10.3)
         write(11,107) ip,ip,gnu_pointtype,gnu_pointsize
107   format('set style line',i3,' lt ',i3,' pt ',i3,' ps ',f7.2)
      end do
      
c decide whether a histogram is needed

      if(histogram(1:1).eq.'y') then
         style='histeps'
         ilenstyle=7
      else
         style='lines'
         ilenstyle=5
      endif
c
c now write the plot command
c
      write(11,108)
      write(6,*) boffset,boffset2
108   format('plot \\')
      do ip=1,npcolumn
         yadd=(ip-1)*ydel
      write(6,*)ip,boffset,bspace,icolumn(ip),iplot
         iplot=boffset+bspace*(icolumn(ip)-1)
         write(11,109) fname1(1:ilenfname1),xcolumn,yfact1,iplot
     *,yadd,style(1:ilenstyle),ip
109   format(1h',a,1h',' u ',i3,':(',f10.5,'*column(',i3,')+',f8.3
     *,') notitle w ',a,' ls ',i3,$)
c
c add the continuation line if necessary
c
         if(ip.ne.npcolumn.or.nplot.ne.1) write(11,110)
110   format(', \\')
         if(nplot.eq.2) then
            iplot=boffset2+bspace*(icolumn(ip)-1)
            write(11,111) fname2(1:ilenfname2),xcolumn,iplot,yfact2,iplot
     *,yadd,ip
111   format(1h',a,1h',' u ',i3,':(abs(column(',i3,'))>=1.0e-5?',f10.5
     *,'*column(',i3,')+',f8.3
     *,':1/0) notitle w p ls ',i3,$)
            if(ip.ne.npcolumn) then
               write(11,110)
            endif
         endif
      end do
c
c finish off the last line
c
      write(11,*)
c     write(11,*) 'pause -1'
      close(11)
      plot_command=system_gnu(1:ilensystem_gnu)
     *//' '
     *//system_quotes(1:ilensystem_quotes)
     *//plot_output(1:ilenplot_output)
     *//system_quotes(1:ilensystem_quotes)
     *//' -'
      call system(plot_command,ierr)
      return
      end

      subroutine read_plot_input_file(nin,nout,ierrout)
c
c reads the PLOT input file
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      integer*4 ierrout
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
c read the commands
c
         nwords=0
         read(nin,100,iostat=ierr) sometext
100   format(a)
         call parse(sometext,nwords,ncf,ncl,12,nstsize,nworddim)
c     write(6,*) sometext(ncf(1):ncl(nwords))

c Need to ignore blank lines read from an input file

         if(nwords.gt.0) 
     1call process_plot_commands(nin,nout,command,icount,ifilecount,ierrout,1)
      end do

c
c finished input - close the input file
c
      close(nin)
      return
      end

      subroutine process_plot_symbols(nin,nout,ierrout,ignore)    
c
c processes commands coming in on the input stream, nin. Diagnostics are written
c to the output stream nout. An 'e' on the input stream signifies the end of input.
c 'data' <number> signifies the following entries relate to dataset <number>. If no
c number is specified, then dataset 1 is assumed.
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      character*20 command
      integer*4 nin,nout,ierrout,ignore
c
c set an error flag so that calling routine can tell if an error has occurred
c
      ierrout=0
      command='next'
      icount=0
c
c set up a counter to say which dataset is being looked at
c
      ifilecount=1
c
c signal whether current values have been saved or not
c
      isave=0
c
c setup an integer for the number of not recognised commands
c
      inotrecog=0
      ierr=0
      do while (command.ne.'e'.and.ierr.eq.0)
         if(icount.gt.0.and.icount.le.nsymbol) then
c
c get the reference value for the symbol
c
            iref=ireff(icount,ifilecount)
c     write(nout,*) icount,iref,ifilecount,nfiles
c
c sort out the words of symbol_value
c
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            ilensymbol=index(symbol(icount),' ')-1
            write(nout,203) prompt(1:ilenprompt),symbol(icount)(1:ilensymbol)
     *,symbol_text(icount)(1:max(1,index(symbol_text(icount),'  ')))
203   format(/a,' ',a,' - ',a)
            write(nout,202) symbol(icount)(1:ilensymbol)
     *,sometext(ncf(1):ncl(max(1,nwords)))
202   format(/a,': ',a,' ? ',$)
         else
c
c reset the file counter if you have got here
c
c           ifilecount=1
            write(nout,201) prompt(1:ilenprompt)
201   format(/a,' ',$)
         endif

         sometext=' '
         read(nin,100,iostat=ierr) sometext
100   format(a)
         call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
         call process_plot_commands(nin,nout,command,icount,ifilecount
     *,ierrout,ignore)
c        write(6,*) command,isave
      end do

      if(isave.ne.1) then
c
c plot files are always saved without question
c
         if(symbol_setup.ne.'plot') then
            nwords=0
            sometext=' '
            write(nout,401) prompt(1:ilenprompt)
401   format(/a,' Current data have not been saved.'
     *//'Press <Enter> to save, or q to exit without saving: ',$)
            read(nin,100,iostat=ierr) sometext
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.gt.0) then
               command=sometext(ncf(1):ncl(1))
            else
               command=' '
            endif
            if(command.ne.'q') then
               call save_plot_files(nin,nout)
            endif
         else
            call save_plot_files(nin,nout)
         endif
      endif   
      return
      end

      subroutine process_plot_commands(nin,nout,command,icount,ifilecount
     *,ierrout,ignore)  
c
c processes command in the input stream, nin. Diagnostics are written
c to the output stream nout. An 'e' on the input stream signifies the end of input.
c 'data' <number> signifies the following entries relate to dataset <number>. If no
c number is specified, then dataset 1 is assumed.
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      character*20 command,text
      integer*4 nout,icount,ifilecount,ifilecountold,ierrout,ignore
c
c set an error flag so that calling routine can tell if an error has occurred
c
      ierrout=0
c
c if nwords .eq. 0 it means the previous change (if any) has been accepted.
c if however symbol_save is already blank, it means this was a blank line
c in which case the symbol counter should not be incremented
c
      if(nwords.gt.0) then
         command=sometext(ncf(1):ncl(1))
      else
         command='d'
      endif
      if(command.eq.'d'.or.command.eq.'u') then
         if(symbol_save.ne.' ') then
            symbol_save=' '
         endif
         ifilecountold=ifilecount

         if(command.eq.'d') then

            icount=icount+1
c
c check the status of the file counter
c
            if(icount.gt.nsymboldata.and.ifilecount.lt.nfiles) then
               ifilecount=ifilecount+1
               icount=nsymboluptodata+1
            else if(icount.gt.nsymbol) then
               icount=0
               ifilecount=1
            endif
c
c signal that you want to write out the file counter
c
            if(icount.eq.nsymboluptodata+1) then
               ifilecountold=0
            else
               ifilecountold=ifilecount
            endif
c     write(6,*) icount,ifilecount,ifilecountold,nsymboluptodata,nsymboldata,nsymbol
         else if(command.eq.'u') then

            icount=icount-1
c
c check the status of the file counter
c
            if(icount.lt.0) then
               icount=nsymbol
               ifilecount=nfiles
            endif
            if(icount.le.nsymboluptodata.and.ifilecount.gt.1) then
               icount=nsymboldata
               ifilecount=ifilecount-1
            endif
c
c signal that you want to write out the file counter
c
            if(icount.eq.nsymboldata) then
               ifilecountold=0
            else
               ifilecountold=ifilecount
            endif
         endif
c
c set up the data text if ifilecount has changed
c
         if(ifilecount.ne.ifilecountold) then
            write(sometext,250) ifilecount
250   format(i3)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.gt.0) then
               datatext=sometext(ncf(1):ncl(1))
               ilendatatext=ncl(1)-ncf(1)+1
            else
               datatext=' '
               ilendatatext=1
            endif
c
c and write it!
c
            if(nin.eq.5) write(nout,204) 
     *dataprompt(1:ilendataprompt),datatext(1:ilendatatext)
     *,dataprompttext(1:max(1,index(dataprompttext,'  ')))
204   format(//a,' ',a,' - ',a/)
         endif
         command='next'
         return
      endif
c
c perform any necessary translations
c
      if(command.eq.'bye') command='e'
c
c reflect the input text into the output stream
c
c     write(nout,*) sometext(ncf(1):ncl(nwords))
c
c first check whether the value of command is one of the special values 
c
      if(command.eq.'default') then
            inotrecog=0
            is=0
            ifilecount=1
            do while(is.lt.nsymboluptodata)
               is=is+1
               symbol_value(is)=symbol_default(is)
            end do
            do while(ifilecount.lt.nfiles)
               is=is+1
               iref=is+nsymbolplus*(ifilecount-1)
               symbol_value(iref)=symbol_default(is)
               if(is.ge.nsymboldata.and.ifilecount.lt.nfiles) then
                  is=nsymboluptodata
                  ifilecount=ifilecount+1
               endif
            end do
            command='next'
      else if(command.eq.'ls'.or.command.eq.'dir') then
            inotrecog=0
            call list_directory()
            command='next'
      else if(command.eq.'cd'.or.command.eq.'pwd') then
            inotrecog=0
            call set_directory_name()
            command='next'
      else if(command.eq.dataprompt) then
            inotrecog=0
c
c try to decide which data set you wish to access
c
            ifound=0
            if(nwords.gt.1) then
               read(sometext(ncf(2):ncl(2)),*,iostat=ierr) ndata
               if(ierr.eq.0) ifound=1
            endif
            if(ifound.eq.0) ndata=1
            if(ndata.lt.1) ndata=1
            if(ndata.gt.nfiles) ndata=nfiles
            ifilecount=ndata
            icount=nsymboluptodata+1
            write(sometext,250) ifilecount
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.gt.0) then
               datatext=sometext(ncf(1):ncl(1))
               ilendatatext=ncl(1)-ncf(1)+1
            else
               datatext=' '
               ilendatatext=1
            endif
c
c and write it!
c
            if(nin.eq.5) write(nout,204) 
     *dataprompt(1:ilendataprompt),datatext(1:ilendatatext)
     *,dataprompttext(1:max(1,index(dataprompttext,'  ')))
c
c for plotting we need to process this command 
c
            if(symbol_setup.eq.'plot'.and.ignore.eq.0) then
               iref=ireff(icount,ifilecount)
c allow this command through the nochanges filter in check_symbol
               symbol_save=symbol_value(iref)
               command='store'
            else
               command='next'
            endif
      else if(command.eq.'help'.or.command.eq.'?') then
            inotrecog=0
c
c write out the list of commands
c
            write(nout,207)
207   format(/'List of parameters:-'//)
            do is=1,nsymbol
               ilensymbol=index(symbol(is),' ')-1
               write(nout,209) symbol(is),symbol_text(is)(1:65)
209   format(a10,' - ',a65)
            end do
            write(nout,210)
210   format(/'Other commands:-'/)
            if(ilendataprompt.gt.0) then
               text=dataprompt(1:ilendataprompt)//' <number> '
               write(nout,211) text(1:20),dataprompttext
211   format(a20,'-',1x,a45)
            endif
            write(nout,208)
208   format(
     *'ls or dir <filename> - Lists files in the current directory'
     */'cd <directory name>  - Changes the current directory'
     */'pwd                  - Show the current directory'
     */'save                 - Saves the current values to file'
     */'default              - Sets numerical values to their defaults'
     */'help or ?            - Show this list'
     */'<value>              - New value for the listed parameter'
     */'<parameter> <value>  - Specified parameter and its new value'
     */'n                    - Revert to the previous value'
     */'<Enter> or d            - Accept current parameter and go forwards'
     */'u                    - Accept current parameter and '
     *,'go backwards'
     */'q                    - Return to EPSR shell without saving'
     */'bye or e             - Return to EPSR shell'
     */)
            command='next'          
      else if(command.eq.'save') then
            inotrecog=0
            call save_plot_files(nin,nout)
            isave=1
            command='next'
      else if(command.eq.'q') then
            inotrecog=0
            isave=1
            command='e'
      else if(command.ne.'e'.and.command.ne.'next'
     *.and.command.ne.'n') then
c
c have not recognised this command so far
c
c
c first check if a valid symbol has been called
c
c (we add 1 to the value of icount since that is the most likely
c value if commands are being read from a file - this is to save search
c time)
c
            ifound=icount+1
            call search_for_symbol(command,ifound)
            if(ifound.gt.0) then
c
c set up the correct reference integer depending on the value of ifilecount
c
               iref=ireff(ifound,ifilecount)
c
c save the current value of this symbol
c
               symbol_save=symbol_value(iref)
c
c can only store the new value if it is defined! (nwords.gt.1)
c
c concatenate the text from the words, starting from the second, and place it
c in the appropriate symbol_value
c
c     write(6,*) ifound,ifilecount,iref,command,' ',symbol_value(iref)(1:30)
               if(nwords.gt.1) then
                  call concatenate(2,nwords,symbol_value(iref),ilenfound
     *,nwords,ncf,ncl,sometext)
                  icount=ifound
                  command='store'
               else

c otherwise go to this parameter and ask for the new value
c if no new value was specified we go to this new value. Normally we do nothing
c at that point, unless these are plot symbols, when action with the current value
c may need to be processed

                  icount=ifound
                  if(symbol_setup.eq.'plot') then
                     command='store'
                  else
                     command='next'
                  endif
               endif
               if(ifound.gt.nsymboldata) ifilecount=nfiles
            else
c
c if icount is greater than zero
c save the text in the appropriate symbol value
c however if the input channel is not 5, it probably means the value is
c spurious since values from a file would be preceded by a recognisable symbol
c
               if(icount.gt.0.and.nin.eq.5) then
c
c set up the correct reference integer depending on the value of ifilecount
c
                  iref=ireff(icount,ifilecount)
c
c save the current value of the symbol
c
                  symbol_save=symbol_value(iref)
                  call concatenate(1,nwords,symbol_value(iref),ilenfound
     *,nwords,ncf,ncl,sometext)
                  command='store'
               else
100   format(/'process_commands> Command not recognised: ',a,': ',a)
                  write(nout,100) command,' ',sometext(1:40)
                  inotrecog=inotrecog+1
                  ierrout=1
                  command='next'
               endif
            endif
      endif
c     write(6,*) command,isave
c
c if command is 'store' it means we need to check the new symbol_value
c and set it back to the saved value if there is an error, or make any other changes
c to the symbol values that are required. See check_symbol for details.
c
      if(command.eq.'store') then

c
c signal that the new value will not have been written to the input file
c
            inotrecog=0
            isave=0

            call check_plot_symbol(nin,nout,icount,ifilecount,ierr,ignore)

c
c set the value of ifilecount if icount does not refer to a dataset
c
c           if(icount.le.nsymboluptodata) ifilecount=1
c           if(icount.gt.nsymboldata) ifilecount=nfiles
            if(ierr.eq.0) then
               command='next'
            else
               command='n'
            endif
c     write(6,*) iref,ierr,symbol_value(iref)(1:30),command
      endif

      if(command.eq.'n') then

c
c change of value was not accepted
c

            inotrecog=0
c
c set up the correct reference integer depending on the value of ifilecount
c
            if(symbol_save.ne.' '.and.icount.gt.0) then
c
c the new value was not accepted so revert to the previous value
c
               iref=ireff(icount,ifilecount)
               symbol_value(iref)=symbol_save

            else
c
c previous saved value was accepted - nothing to change
c
               icount=icount+1
            endif
            command='next'
      endif

      return
      end

      subroutine check_plot_symbol(nin,nout,icount,ifilecount,ierr,ignore)

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'

      integer*4 icount,ifilecount,ierr,ignore,iref
      integer*4 icheck
      real*4 rcheck
      character*20 ext
      character*256 text
      character*1 iaccept
      integer*4 ilenext
c     write(6,*) symbol_value(iinitref)(1:10),iinitref

c
c make sure icount does not exceed the maximum allowed value
c
      if(icount.lt.1) icount=1
      if(icount.gt.nsymbol) icount=nsymbol
c     if(icount.gt.nsymboldata) ifilecount=nfiles
      iref=ireff(icount,ifilecount)
c
c set error flag to zero
c
      ierr=0
c
c check that it is allowed to change this value - if a change has been requested
c
      if(iflag(icount).gt.0.and.ignore.eq.0
     *.and.symbol_value(iref).ne.symbol_save) then
            write(nout,200) prompt(1:ilenprompt)
200   format(/a,1x,'Cannot change this value!')
            symbol_value(iref)=symbol_save
            return
      endif
c
c now check that the specified value(s) is(are) of the correct type if it is not 
c a character variable
c
      if(symbol_type(icount).ne.'c') then
         sometext=symbol_value(iref)
c     write(6,*) symbol_value(iref),sometext(1:50)
         call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
c     write(6,*) symbol_value(iref),sometext(1:50),nwords,ncf(1),ncl(1)
         iw=0
         if(symbol_type(icount).eq.'i') then
            do while(iw.lt.nwords.and.ierr.eq.0)
               iw=iw+1
               read(sometext(ncf(iw):ncl(iw)),*,iostat=ierr) icheck
               if(ierr.ne.0) then
                  write(6,201) prompt(1:ilenprompt),sometext(ncf(iw):ncl(iw))
201   format(/a,1x,'Value ',a,' is not an integer - try again!')
                  symbol_value(iref)=symbol_save
                  ierr=0
                  return
               endif
            end do
         else
c     write(6,*) nwords,iw,sometext(ncf(1):ncl(1))
            do while(iw.lt.nwords.and.ierr.eq.0)
               iw=iw+1
               read(sometext(ncf(iw):ncl(iw)),*,iostat=ierr) rcheck
c     write(6,*) ierr,rcheck
               if(ierr.ne.0) then
                  write(nout,202) prompt(1:ilenprompt),sometext(ncf(iw):ncl(iw))
202   format(/a,1x,'Value ',a,' is not a real number - try again!')
                  symbol_value(iref)=symbol_save
                  ierr=0
                  return
               endif
            end do
         endif
      endif
c
c if icount .eq. nsymboluptodata the number of files may have been
c changed. This means the values of symbol_value above nsymboluptodata will 
c need to be resequenced
c
      if(icount.eq.nsymboluptodata) then
c
c read the new value of nfiles into a temporary location
c
         read(symbol_value(iref),*,iostat=ierr) nfilet

         if(ierr.eq.0) then
c
c First reassign the symbol_values above nsymboldata
c
            nfadd=(nfilet-nfiles)*nsymbolplus
c
c If new number of files is greater than old this must be done in reverse order
c to avoid overwriting existing values
c
            if(nfadd.gt.0) then
               isymbol=nsymbol+1
               ivalue=ireff(nsymbol,nfiles)+1
               do while(isymbol.gt.nsymboldata+1)
                  isymbol=isymbol-1
                  ivalue=ivalue-1
                  ivaluenew=ivalue+nfadd
c     write(6,*) nsymboluptodata,nsymboldata,nsymbol,isymbol,ivalue,ivaluenew
                  if(ivaluenew.le.msymbol) then
                     symbol_value(ivaluenew)=symbol_value(ivalue)
                  else
                     write(nout,203) prompt(1:ilenprompt)
203   format(/a,1x,'ERROR - tried to generate too many symbol values')
                     ierr=1
                     return
                  endif
               end do
            else if(nfadd.lt.0) then
               isymbol=nsymboldata
               ivalue=ireff(nsymboldata,nfiles)
               do while(isymbol.lt.nsymbol)
                  isymbol=isymbol+1
                  ivalue=ivalue+1
                  ivaluenew=ivalue+nfadd
c     write(6,*) nsymboluptodata,nsymboldata,nsymbol,isymbol,ivalue,ivaluenew
                  if(ivaluenew.le.msymbol) then
                     symbol_value(ivaluenew)=symbol_value(ivalue)
                  else
                     write(nout,203) prompt(1:ilenprompt)
                     ierr=1
                     return
                  endif
               end do
c
c no need to do anything if nfadd = 0
c
            endif
            ifiletemp=nfiles
            nfiles=nfilet
c
c if new number of files .gt. old number of files
c set the new values of symbol_value to the default values
c
            if(nfadd.gt.0) then
               do while(ifiletemp.lt.nfiles)
                  ifiletemp=ifiletemp+1
                  isymbol=nsymboluptodata
                  do while(isymbol.lt.nsymboldata)
                     isymbol=isymbol+1
                     iref=ireff(isymbol,ifiletemp)
                     symbol_value(iref)=symbol_default(isymbol)
                  end do
               end do
            endif
         endif
      endif

c
c plotting commands are processed by their own routine, but not during reading them in
c
      if(symbol_setup.eq.'plot'.and.ignore.eq.0) then
         call process_plots(nin,nout,icount,ifilecount,ierr)
         return
      endif
c
c changes specific to 2D and 3D plotting
c

      if(symbol_setup.eq.'plot2d'.or.symbol_setup.eq.'plot3d') then

         if(symbol(icount).eq.'shcoeffs') then
c
c changing the name of the coefficients file - read the corresponding .dat file
c
            text=symbol_value(ireff(nsymbol,nfiles))
            nwords=0
            call parse(text,nwords,ncf,ncl,5,256,nworddim)
            ext=text(ncf(1):ncl(1))
            ilenext=ncl(1)-ncf(1)+1
            nwords=0
c
c if symbol_value(1) is 'search' we need to find a valid SHARM coefficients file
c
            if(symbol_value(iref).eq.'search'
     *.or.symbol_value(iref).eq.'<undefined>') then
               call search_for_files(ext,ilenext,iaccept)
               call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            else
               sometext=symbol_value(iref)
               call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
c
c add the extension if not supplied
c
               ilendot=index(sometext(ncf(1):ncl(1)),'.')
               if(ilendot.le.0) then
                  sometext=sometext(ncf(1):ncl(1))//ext
                  ncl(1)=ncl(1)+ilenext
               endif
            endif
            if(nwords.le.0) then
               write(nout,206) prompt(1:ilenprompt)
206   format(/a,1x,'No coefficients file specified! - routine will exit'/)
               symbol_value(iref)=symbol_save
               ierr=0
               return
            else
               symbol_value(iref)=sometext(ncf(1):ncl(1))
            endif
            call find_plot_par(12,nout,ierr)
         endif
c
c if setone has been changed, make sure new value is either 0 or 1
c
         if(symbol(icount).eq.'setone') then
            read(symbol_value(iref),*,iostat=ierr) icheck
            if(ierr.ne.0.or.icheck.ne.0.or.icheck.ne.1) then
               write(nout,207) prompt(1:ilenprompt)
207   format(/a,1x,'setone must be 0 or 1!')
               symbol_value(iref)=symbol_save
               ierr=0
               return
            endif
         endif          
c
c if nvary has changed (plot2d, plot3d) need to change symbol_text and symbol_value
c
         if(symbol(icount).eq.'nvary') then
c
c read the new value of nvary
c
            iref=ireff(icount,ifilecount)
            read(symbol_save,*,iostat=ierr1) nvary
            read(symbol_value(iref),*,iostat=ierr2) nvaryt
            if(ierr1.eq.0.and.ierr2.eq.0.and.nvaryt.ne.nvary) then
               iref=iref+1
               if(symbol_setup.eq.'plot3d') then
                  if(nvaryt.eq.1) then
                     symbol_text(iref)='phim, thetam, chim '
     *//'(degrees - separated by spaces)'
                  else if(nvaryt.eq.2) then
                     symbol_text(iref)='thetal, phil, chim '
     *//'(degrees - separated by spaces)'
                  else if(nvaryt.eq.3) then
                     symbol_text(iref)='thetal, phil, phim '
     *//'(degrees - separated by spaces)'
                  endif
               else if(symbol_setup.eq.'plot2d') then
                  if(nvaryt.eq.1) then
                     symbol_text(iref)='phil, phim, thetam, chim '
     *//'(degrees - separated by spaces)'
                  else if(nvaryt.eq.2) then
                     symbol_text(iref)='thetal, phim, thetam, chim '
     *//'(degrees - separated by spaces)'
                  else if(nvaryt.eq.3) then
                     symbol_text(iref)='thetal, phil, thetam, chim '
     *//'(degrees - separated by spaces)'
                  else if(nvaryt.eq.4) then
                     symbol_text(iref)='thetal, phil, phim, chim '
     *//'(degrees - separated by spaces)'
                  else if(nvaryt.eq.5) then
                     symbol_text(iref)='thetal, phil, phim, thetam '
     *//'(degrees - separated by spaces)'
                  endif
               endif
               symbol_value(iref)=symbol_default(icount+1)
            else
               ierr=1
            endif
         endif
      endif

      return
      end

      subroutine read_plot2d(nin,nout,ierrout)
c
c sets up the plot2d input files
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      integer*4 nin,nout,ierr,ierrout
c
c return if an error occurred on trying to open the input file
c
      if(ierrout.ne.0) then
            close(nin)
            return
      endif
      write(6,*) fname_run(1:ilenfname_run)
c
c  If O.K. read in the data from the specified input file
c
      symbol_value(1)=fname_run(1:ilenfname_run)
c
c commands to read the rest of the input file
c
      icount=1
      ifilecount=1
      ierr=0
      do while(icount.lt.nsymbol.and.ierr.eq.0)
         icount=icount+1
         iref=ireff(icount,ifilecount)
         symbol_save=symbol_value(iref)
         read(nin,100,iostat=ierr) sometext
         if(ierr.eq.0) then
100   format(a)
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         if(nwords.gt.0) then
            if(icount.eq.nsymboluptodata+1) then
               call concatenate(1,3,symbol_value(iref),ilen
     *,nwords,ncf,ncl,sometext)
               iref1=iref+1
               call concatenate(4,4,symbol_value(iref1),ilen1
     *,nwords,ncf,ncl,sometext)
               icount=icount+1
            else
               call concatenate(1,nwords,symbol_value(iref),ilen
     *,nwords,ncf,ncl,sometext)
               if(icount.eq.nsymboluptodata) then
                  read(symbol_value(iref),*) nfiles
               endif

c if symbol_value .eq. nvary has changed need to change symbol_text for the 
c next line.

               if(symbol(icount).eq.'nvary') then
                  read(symbol_value(iref),*,iostat=ierr2) nvaryt
                  if(nvaryt.eq.1) then
                     symbol_text(iref+1)='phim, thetam, chim '
     *//'(degrees - separated by spaces)'
                  else if(nvaryt.eq.2) then
                     symbol_text(iref+1)='thetal, phil, chim '
     *//'(degrees - separated by spaces)'
                  else if(nvaryt.eq.3) then
                     symbol_text(iref+1)='thetal, phil, phim '
     *//'(degrees - separated by spaces)'
                  endif
               endif

            endif
         else
            symbol_value(iref)=' '
         endif
         if(icount.ge.nsymboldata.and.ifilecount.lt.nfiles) then
            ifilecount=ifilecount+1
            icount=nsymboluptodata
         endif
         endif
      end do
      close(nin)
c
c determine the number of coefficients in the coefficients file and whether the
c molecules are identical or distinct
c
      call find_plot_par(nin,nout,ierrout)
      return
      end

      subroutine read_plot3d(nin,nout,ierrout)
c
c sets up the plot3d input files
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      integer*4 nin,nout,ierr,ierrout
c
c return if an error occurred on trying to open the input file
c
      if(ierrout.ne.0) then
            close(nin)
            return
      endif
      write(6,*) fname_run(1:ilenfname_run)
c
c  If O.K. read in the data from the specified input file
c
      symbol_value(1)=fname_run(1:ilenfname_run)
c
c commands to read the rest of the input file
c
      icount=1
      ifilecount=1
      ierr=0
      do while(icount.lt.nsymbol.and.ierr.eq.0)
         icount=icount+1
         iref=ireff(icount,ifilecount)
         symbol_save=symbol_value(iref)
         read(nin,100,iostat=ierr) sometext
         if(ierr.eq.0) then
100   format(a)
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         if(nwords.gt.0) then
            if(icount.eq.nsymboluptodata+1) then
               call concatenate(1,1,symbol_value(iref),ilen
     *,nwords,ncf,ncl,sometext)
               iref1=iref+1
               call concatenate(2,4,symbol_value(iref1),ilen1
     *,nwords,ncf,ncl,sometext)
               iref2=iref+2
               call concatenate(5,7,symbol_value(iref2),ilen2
     *,nwords,ncf,ncl,sometext)
               icount=icount+2
            else
               call concatenate(1,nwords,symbol_value(iref),ilen
     *,nwords,ncf,ncl,sometext)
               if(icount.eq.nsymboluptodata) then
                  read(symbol_value(iref),*) nfiles
               endif

c if symbol_value .eq. nvary has changed need to change symbol_text for the 
c next line.

               if(symbol(icount).eq.'nvary') then
                  read(symbol_value(iref),*,iostat=ierr2) nvaryt
                  if(nvaryt.eq.1) then
                     symbol_text(iref+1)='phim, thetam, chim '
     *//'(degrees - separated by spaces)'
                  else if(nvaryt.eq.2) then
                     symbol_text(iref+1)='thetal, phil, chim '
     *//'(degrees - separated by spaces)'
                  else if(nvaryt.eq.3) then
                     symbol_text(iref+1)='thetal, phil, phim '
     *//'(degrees - separated by spaces)'
                  endif
               endif

            endif
         else
            symbol_value(iref)=' '
         endif
         if(icount.ge.nsymboldata.and.ifilecount.lt.nfiles) then
            ifilecount=ifilecount+1
            icount=nsymboluptodata
         endif
         endif
      end do
      close(nin)
c
c determine the number of coefficients in the coefficients file and whether the
c molecules are identical or distinct
c
      call find_plot_par(nin,nout,ierrout)
      return
      end

      subroutine find_plot_par(nin,nout,ierr)
c
c finds certain plotting parameters from the SHARM input file
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'symbols.inc'
      integer*4 nin,nout,ierr,ilensymbol_value
      character*256 fnametemp
      character*256 centrec,centres
      character*10 text
      character*20 moretext

c search for the symbol "shcoeffs"

      moretext='shcoeffs'
      ifound=0
      call search_for_symbol(moretext,ifound)
      iref=ireff(ifound,nfiles)
c
c generate the SHARM input filename
c
      ilensymbol_value=indexlastdot(symbol_value(iref))-1
      if(ilensymbol_value.le.0) then
         ilensymbol_value=index(symbol_value(iref),' ')-1
      endif
      fnametemp=directory_name(1:nchardir)
     *//symbol_value(iref)(1:ilensymbol_value)//'.dat'
      write(nout,*) fnametemp(1:(nchardir+ilensymbol_value+4))
      open(nin,file=fnametemp,status='old',iostat=ierr)
      if(ierr.ne.0) then
         write(nout,100) prompt(1:ilenprompt)
     *,symbol_value(iref)(1:ilensymbol_value)//'.dat'
100   format(/a,' SHARM input file ',a,' does not exist'
     *//'Will try to get number of coefficients directly'
     *,' from coefficients file')
         fnametemp=directory_name(1:nchardir)
     *//symbol_value(iref)(1:ilensymbol_value)//'.h01'
         text='.h01'
         nin1=nin+1
         call get_nterms(nin1,fnametemp,text,ncoeffs)
c
c write this to symbol_value(iref+1)
c
         write(symbol_value(iref+1),*) ncoeffs
         return
      endif
c
c read lines of text from this file, looking for specific keywords
c
101   format(a)
      ierr=0
      centrec=' '
      centres=' '
      do while(ierr.eq.0.and.(centrec.eq.' '.or.centres.eq.' '))
         read(nin,101,iostat=ierr) sometext
         call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
c
c check for keyword
c
102   format(a,1x,a)
         if(sometext(ncf(1):ncl(1)).eq.'ncoeffs') then
            symbol_value(iref+1)=sometext(ncf(2):ncl(2))
         else if(sometext(ncf(1):ncl(1)).eq.'atom-c') then
            centrec=sometext(ncf(2):ncl(nwords))
         else if(sometext(ncf(1):ncl(1)).eq.'atom-s') then
            centres=sometext(ncf(2):ncl(nwords))
         endif
      end do
      close(nin)
c
c if centrec or centres have not been set it probably means there is an error
c with the .dat file
c
      if(centrec.eq.' '.or.centres.eq.' ') then
         write(nout,103) prompt(1:ilenprompt)
     *,symbol_value(iref)(1:ilensymbol_value)//'.dat'
103   format(/a,' Problems with file ',a
     *//'Will try to get number of coefficients directly'
     *,' from coefficients file')
         fnametemp=directory_name(1:nchardir)
     *//symbol_value(iref)(1:ilensymbol_value)//'.h01'
         text='.h01'
         nin1=nin+1
         call get_nterms(nin1,fnametemp,text,ncoeffs)
c
c write this to symbol_value(iref+1)
c
         write(symbol_value(iref+1),*) ncoeffs
         ierr=0
      else
         ierr=0
      endif
c
c decide whether molecules were the same or different
c
      if(centrec.ne.centres) then
         symbol_value(iref+2)='1'
      else
         symbol_value(iref+2)='0'
      endif

      return
      end

      subroutine save_plot_files(nin,nout)
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      integer*4 ilenuptodot
c
c first check whether the existing file name is O.K.
c
      ierrfile=1
      do while(ierrfile.ne.0)

         if(ilenfname_inp.gt.0) then

            write(nout,303) prompt(1:ilenprompt),fname_inp(1:ilenfname_inp)
303   format(/a,1x,'Current name of file is "',a,'"')
c     *//'Specify the new filename or press <Enter> to accept this name: '
c     *,$)
c           read(nin,100,iostat=ierr) sometext
100   format(a)
c           call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
c           if(nwords.gt.0) then
c              call concatenate(1,nwords,fname_inp,ilenfname_inp
c     *,nwords,ncf,ncl,sometext)
c           endif

         else

            nwords=0
            do while(nwords.eq.0) 
               write(nout,302) prompt(1:ilenprompt),symbol_setup(1:ilensetup)
302   format(/a,1x,' Specify the name of the ',a,' input file: ',$)
               read(nin,100,iostat=ierr) sometext
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
            ilenfname_inp=index(fname_inp,'  ')-1
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
            nwords=0
            do while(nwords.eq.0)
               sometext=' '
               write(6,301) prompt(1:ilenprompt),fname_inp(1:ilenfname_inp)
301   format(/a,1x,'File "',a,'" already exists.'
     *//' Do you want to overwrite it (y or n)? ',$)
               read(5,100,iostat=ierr) sometext
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
      if(symbol_setup.eq.'plot2d') then
         call save_files_plot2d(14)
      else if(symbol_setup.eq.'plot3d') then
         call save_files_plot3d(14)
      else
         call save_files_plots(14)
      endif

      close(14)
      return
      end

      subroutine save_files_plots(nout)
c
c writes the general input file
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
c
c open the new input file
c
      open(nout,file=fnameinp,status='unknown',iostat=ierr)
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
103   format(/'save_files_plots> Could not open specified file for writing.')
         return

      endif
      return
      end

      subroutine save_files_plot2d(nout)
c
c writes the plot2d input file
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
c
c open the new input file
c
      open(nout,file=fnameinp,status='unknown',iostat=ierr)
c
c write .txt file
c
      icount=0
      ifilecount=1
      do while (icount.lt.nsymbol)
         icount=icount+1
         iref=ireff(icount,ifilecount)
c
c lines containing sphere parameters has to be printed as a single line
c
         if(icount.eq.nsymboluptodata+1) then
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            call concatenate(1,nwords,symbol_value(iref),ilen
     *,nwords,ncf,ncl,sometext)
            iref1=iref+1
            sometext=symbol_value(iref1)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            call concatenate(1,nwords,symbol_value(iref1),ilen1
     *,nwords,ncf,ncl,sometext)
            sometext=symbol_value(iref)(1:ilen)
     *//' '//symbol_value(iref1)(1:ilen1)
            write(nout,102) sometext(1:(ilen+ilen1+1))
c102  format(a,11x,'(x,y,z) of circle, and circle radius')
102   format(a,11x)
            icount=icount+1
         else
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
c
c character values cannot be printed with the corresponding text
c
            if(symbol_type(icount).eq.'c') then
               write(nout,103) sometext(ncf(1):ncl(max(1,nwords)))
103   format(a)
            else
               write(nout,101) sometext(ncf(1):ncl(max(1,nwords)))
     *,symbol_text(icount)
            endif
101   format(a,15x,a)
         endif
         if(icount.ge.nsymboldata.and.ifilecount.lt.nfiles) then
            ifilecount=ifilecount+1
            icount=nsymboluptodata
         endif
      end do

      return
      end

      subroutine save_files_plot3d(nout)
c
c writes the EPSR .inp and .pcof files
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
c
c open the new input file
c
      open(nout,file=fnameinp,status='unknown',iostat=ierr)
c
c write .txt file
c
      icount=0
      ifilecount=1
      do while (icount.lt.nsymbol)
         icount=icount+1
         iref=ireff(icount,ifilecount)
c
c lines containing sphere parameters has to be printed as a single line
c
         if(icount.eq.nsymboluptodata+1) then
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            call concatenate(1,nwords,symbol_value(iref),ilen
     *,nwords,ncf,ncl,sometext)
            iref1=iref+1
            sometext=symbol_value(iref1)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            call concatenate(1,nwords,symbol_value(iref1),ilen1
     *,nwords,ncf,ncl,sometext)
            iref2=iref+2
            sometext=symbol_value(iref2)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            call concatenate(1,nwords,symbol_value(iref2),ilen2
     *,nwords,ncf,ncl,sometext)
            sometext=symbol_value(iref)(1:ilen)
     *//' '//symbol_value(iref1)(1:ilen1)
     *//' '//symbol_value(iref2)(1:ilen2)
            write(nout,102) sometext(1:(ilen+ilen1+ilen2+2))
102   format(a,11x,'sphere radius, (r,theta,phi), (r,g,b colour indices)')
            icount=icount+2
         else
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
c
c character values cannot be printed with the corresponding text
c
            if(symbol_type(icount).eq.'c') then
               write(nout,103) sometext(ncf(1):ncl(max(1,nwords)))
103   format(a)
            else
               write(nout,101) sometext(ncf(1):ncl(max(1,nwords)))
     *,symbol_text(icount)
            endif
101   format(a,15x,a)
         endif
         if(icount.ge.nsymboldata.and.ifilecount.lt.nfiles) then
            ifilecount=ifilecount+1
            icount=nsymboluptodata
         endif
      end do

      return
      end




