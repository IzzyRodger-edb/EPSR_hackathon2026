      subroutine process_symbols(nin,nout,ierrout,ignore)   
c
c processes commands coming in on the input stream, nin. Diagnostics are written
c to the output stream nout. An 'e' on the input stream signifies the end of input.
c 'data' <number> signifies the following entries relate to dataset <number>. If no
c number is specified, then dataset 1 is assumed.
c
      use jmk_readline_module
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      character*20 command
      character*256 promptstring
      integer*4 nin,nout,ierrout,ignore,lenpromptstring
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
c
c sort out the words of symbol_value
c
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            ilensymbol=index(symbol(icount),' ')-1
            write(nout,203) prompt(1:ilenprompt),symbol(icount)(1:ilensymbol)
     *,symbol_text(icount)(1:max(1,index(symbol_text(icount),'  ')))
203   format(/a,' ',a,' - ',a/)
            write(promptstring,'(a,": ",a," ? ")') symbol(icount)(1:ilensymbol)
     *,sometext(ncf(1):ncl(max(1,nwords)))
202   format(/a,': ',a,' ? ')
         else
c
c reset the file counter if you have got here
c
c           ifilecount=1
            write(*,*) !Blank line
            write(promptstring,'(a,1x)') prompt(1:ilenprompt)
201   format(/a,' ')
         endif

         lenpromptstring=len_trim(promptstring)
         sometext=' '
         call readline(sometext,promptstring(1:lenpromptstring))
!         read(nin,100,iostat=ierr) sometext
100   format(a)
         call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
c         write(nout,*) nwords,' ',sometext(ncf(1):ncl(nwords)
         call process_commands(nin,nout,command,icount,ifilecount
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
            write(nout,401,advance='no') prompt(1:ilenprompt)
401   format(/a,' Current data have not been saved.'
     *//'Press <Enter> to save, or q to exit without saving: ')
            read(nin,100,iostat=ierr) sometext
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.gt.0) then
               command=sometext(ncf(1):ncl(1))
            else
               command=' '
            endif
            if(command.ne.'q') then
               call save_files(nin,nout)
            endif
         else
            call save_files(nin,nout)
         endif
      endif
      call save_history_list()   
      return
      end

      subroutine process_commands(nin,nout,command,icount,ifilecount
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
     */'forcesave            - Saves the current values to file directly without checking if OK to do so'
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
            call save_files(nin,nout)
            isave=1
            command='next'
      else if(command.eq.'forcesave') then
            inotrecog=0
            call force_save_files(nin,nout,.true.)
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
               if(ifound.lt.nsymboluptodata+1) ifilecount=1
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

            call check_symbol(nin,nout,icount,ifilecount,ierr,ignore)

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
c
c for some modules this may mean we need to reset more than just the symbol_value
c
               if(symbol_setup.eq.'changeato') then
                  call check_symbol(nin,nout,icount,ifilecount,ierr,ignore)
               endif

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


      subroutine check_symbol(nin,nout,icount,ifilecount,ierr,ignore)

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'iresetnumbers.inc'
      include 'potential_par.inc'

      integer*4 icount,ifilecount,ierr,ignore,iref
      integer*4 icheck
      real*4 rcheck
      character*20 ext
      character*256 text
      character*1 iaccept
      integer*4 ilenext,ilentext
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
c if the name of the ato file has changed, we need to read in the
c new values
c
      if(symbol(icount).eq.'fnameato') then
c     write(6,*) 'Got to fnameato check_symbol'
         call read_ato(nin,nout,icount,ifilecount,ierrout)
      endif

c
c if the name of the pcof file has changed, we need to read in the
c new values
c
      if(symbol(icount).eq.'fnamepcof') then
c     write(6,*) 'Got to fnamepcof check_symbol'
         icountsave=icount
         call read_pcof(nin,nout,icount,ifilecount,ierrout)
         icount=icountsave
      endif

c
c ATOM CHECK
c
c If the first four characters of symbol are 'atom', it means the associated variable
c is a list of one or more atom types. We need to check that these atom types actually
c exist in the current .ato file.
c
      if(symbol(icount)(1:4).eq.'atom') then
c
c First, set up the list of atoms
c
         sometext=symbol_value(iref)
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)

c now check that the atom types actually exist in the current .ato file

         iw=0
         do while(iw.lt.nwords)
            iw=iw+1

c Search for a * in the atom type: atom types will only be matched up to the *

            ncstar=index(sometext(ncf(iw):ncl(iw)),'*')
            if(ncstar.gt.0) then
               ncsearch=ncstar-1
               ncluse=ncf(iw)+ncsearch-1
            else
               ncsearch=ncl(iw)-ncf(iw)+1
               ncluse=ncl(iw)
            endif
c
c step through the atom types and check that one exists with the same type
c as that specified
c
            it=0
            ifound=0
            do while(it.lt.ntype.and.ifound.eq.0) 
               it=it+1

c If number of characters being searched is .gt. 0, then check the atom type

               if(ncsearch.gt.0) then

                  if(type(it)(1:ncsearch).eq.sometext(ncf(iw):ncluse)) ifound=it

               else

c Otherwise always accept this atom type

                  ifound=it

               endif
            end do
            if(ifound.eq.0) then
               write(nout,205) prompt(1:ilenprompt)
     *,sometext(ncf(iw):ncl(iw)),fname_ato(1:ilenfname_ato)
205   format(/a,1x,'WARNING! Atom type ',a
     *,' not found in .ato file ',a)
               symbol_value(iref)=symbol_save
            endif
         end do
      endif

c
c changes specific to SHARM
c

      if(symbol_setup.eq.'sharm') then

         

c reformat the symbol_values (11) and (14) so there is only 1 space between words

         if(iref.eq.11.or.iref.eq.14) then
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            call concatenate(1,nwords,symbol_value(iref),ilenfound
     *,nwords,ncf,ncl,sometext)

         endif
c
c if centres is the same as centrec set the axis definition for the second molecule 
c to be the same as the central molecule
c
         if(iref.ge.14.and.symbol_value(11).eq.symbol_value(14)) then
c set axiss1 same as axisc1
            symbol_value(15)=symbol_value(12)
c set axiss2 same as axisc2
            symbol_value(16)=symbol_value(13)
         endif
      endif

c
c changes specific to SDF
c

      if(symbol_setup.eq.'sdf') then

c reformat the symbol_valueS (11) and (16) so there is only 1 space between words

         if(iref.eq.11.or.iref.eq.16) then
            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            call concatenate(1,nwords,symbol_value(iref),ilenfound
     *,nwords,ncf,ncl,sometext)

         else if(iref.eq.12) then

c Check that the correct number of bonds have been read in

            call parse(symbol_value(11),nwords,ncf,ncl,5,256,nworddim)
            nbondexpect=nwords*(nwords-1)/2
            call parse(symbol_value(12),nwords,ncf,ncl,5,256,nworddim)
            if(nwords.lt.nbondexpect) write(6,501) prompt(1:ilenprompt)
     *,nbondexpect
501   format(/a,1x,'WARNING! Number of bonds specified is smaller '
     *,'than expected. Should be: ',i3)

         else if(iref.ge.16) then

c For second molecule if it only has one atom, ensure that no atoms are specified for
c the coordinate axes

            call parse(symbol_value(16),nwords,ncf,ncl,5,256,nworddim)
            naxisexpect=nwords
            call parse (symbol_value(19),nwords,ncf,ncl,5,256,nworddim)
            if(naxisexpect.eq.1.and.nwords.gt.1) then
               write(6,502) prompt(1:ilenprompt)
502   format(/a,1x,'WARNING! No need for axis atoms in this case')
               symbol_value(19)=symbol_value(19)(ncf(1):ncl(1))
            endif 
            call parse (symbol_value(20),nwords,ncf,ncl,5,256,nworddim)
            if(naxisexpect.eq.1.and.nwords.gt.1) then
               write(6,502) prompt(1:ilenprompt)
               symbol_value(20)=symbol_value(20)(ncf(1):ncl(1))
            endif

         else if(iref.eq.17) then

c Check that the correct number of bonds have been read in

            call parse(symbol_value(16),nwords,ncf,ncl,5,256,nworddim)
            nbondexpect=nwords*(nwords-1)/2
            call parse(symbol_value(17),nwords,ncf,ncl,5,256,nworddim)
            if(nwords.lt.nbondexpect) write(6,501) prompt(1:ilenprompt)
     *,nbondexpect

         endif

         if(symbol_value(11).eq.symbol_value(16)) then

c if atom-s is the same as atom-c set the axis definition for the second molecule 
c to be the same as the central molecule

c set bond distances to be the same for both molecules
            symbol_value(17)=symbol_value(12)
c set acceptance fractions to be the same
            symbol_value(18)=symbol_value(13)
c set axiss1 same as axisc1
            symbol_value(19)=symbol_value(14)
c set axiss2 same as axisc2
            symbol_value(20)=symbol_value(15)
         endif
        
      endif

c if requested search for a suitable data or wts file

      if((symbol(icount).eq.'datafile'
     *.or.symbol(icount).eq.'wtsfile')
     *.and.symbol_value(iref).eq.'search') then
         if(symbol(icount).eq.'datafile') then
            nwords=0
            do while(nwords.le.0)
c
c specify the extension to search on
c
               write(nout,120,advance='no')
120   format(/'Specify the extension to search on: ')
               read(nin,100) sometext
100   format(a)
               call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            end do
            call concatenate(1,1,text,ilentext
     *,nwords,ncf,ncl,sometext)
            if(ilentext.le.20) then
               ext=text(1:ilentext)
               ilenext=ilentext
            else
               ext=' '
               ilenext=1
            endif
         else
            ext='.wts'
            ilenext=4
         endif
         call search_for_files(ext,ilenext,iaccept)
         if(iaccept.ne.'y') then
            write(nout,121) prompt(1:ilenprompt)
121   format(/a,1x,'No file specified - value will not be changed'/)
            symbol_value(iref)=symbol_save
            ierr=0
         else
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            call concatenate(1,nwords,symbol_value(iref),ilen
     *,nwords,ncf,ncl,sometext)
            ierr=0
         endif
      endif
! Changes specific to mapgrentropy
      if(symbol_setup.eq.'mapgrentropy'.or.symbol_setup.eq.'mapgr') then

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
      write(6,'(a)') ext(1:ilenext)
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
         end if
      end if
! Other changes specific to mapgr
      if(symbol_setup.eq.'mapgr') then
c
c if nvary has changed mapgr needs to change symbol_text for all the files.
c
         if(symbol(icount).eq.'nvary') then
c
c read the new value of nvary
c
            iref=ireff(icount,ifilecount)
            read(symbol_save,*,iostat=ierr1) nvary
            read(symbol_value(iref),*,iostat=ierr2) nvaryt
            if(ierr1.eq.0.and.ierr2.eq.0.and.nvaryt.ne.nvary) then
               do i=1,nfiles
                  iref=ireff(nsymboluptodata,i)
                  if(nvaryt.eq.1) then
                     iref=iref+1
                     symbol_text(iref)='thetal and phil of view (degrees - separated by spaces)'
                     iref=iref+1
                     symbol_text(iref)='phim, thetam, chim '
     *//'(degrees - separated by spaces)'
                  else if(nvaryt.eq.2) then
                     iref=iref+1
                     symbol_text(iref)='thetam and phim of view (degrees - separated by spaces)'
                     iref=iref+1
                     symbol_text(iref)='thetal, phil, chim '
     *//'(degrees - separated by spaces)'
                  else if(nvaryt.eq.3) then
                     iref=iref+1
                     symbol_text(iref)='thetam and chim of view (degrees - separated by spaces)'
                     iref=iref+1
                     symbol_text(iref)='thetal, phil, phim '
     *//'(degrees - separated by spaces)'
                  endif
               end do
            end if
         end if
      end if
c
c if nvary has changed (plot2d, splot2d, plot3d) need to change symbol_text and symbol_value
c
c
c changes specific to EPSR
c
      if(symbol_setup.eq.'epsr') then

c
c if a parameter with a negative iflag is changed, we need to set iinit or ireset to 1
c (but not when reading from .pcof file - otherwise ireset will always be set to 1!)
c
         irefiinit=ireff(iinitref,ifilecount)
         irefireset=ireff(iresetref,ifilecount)
         if(ignore.eq.0) then
            if(iflag(iref).lt.0) then
               symbol_value(irefiinit)='1'
            endif
            if(iflag(iref).eq.-2) then
               symbol_value(irefireset)='1'
            else if(iflag(iref).lt.-2) then
               symbol_value(irefireset)='2'
            endif
         endif
c
c if ireset has been set, need to set iinit to 1 as well
c
         if(symbol_value(irefireset).ne.'0') then
            symbol_value(irefiinit)='1'
         endif

c Ensure the number of potential coefficients does not exceed allowed dimension

         if(symbol(icount).eq.'ncoeffp') then

            irefncoeffp=ireff(incoeffp,ifilecount)
            read(symbol_value(irefncoeffp),*,iostat=ierr) ncoeffp

c Make sure this number does not exceed the allowed dimension

            if(ncoeffp.gt.mcoeff) then

               write(6,300) mcoeff
300   format(/'check_symbol> Too many potential coefficients defined.'
     1,' Will reduce to maximum allowed (',i3,')')
               ncoeffp=mcoeff
               write(sometext,*,iostat=ierr) ncoeffp
               call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
               call concatenate(1,nwords,symbol_value(irefncoeffp),ilenfound
     *,nwords,ncf,ncl,sometext)

            endif

         endif

      endif

c Changes specific to nwts

      if(symbol_setup.eq.'nwts') then

c         write(6,*) iref,icount,ifilecount

c If fnameato has been set then get the number of components

         if(symbol(icount).eq.'fnameato'.and.symbol_value(iref).ne.symbol_default(iref)) then

            nfilet=ntype
            if(nfilet.lt.10) then
               write(symbol_value(nsymboluptodata),'(i1)') nfilet
            else if(nfilet.lt.100) then
               write(symbol_value(nsymboluptodata),'(i2)') nfilet
            else
               write(symbol_value(nsymboluptodata),'(i3)') nfilet
            endif

            ifiletemp=nfiles
            nfiles=nfilet
      write(6,*) nfilet,ifiletemp,nfiles
c
c if new number of files .gt. old number of files
c set the new values of symbol_value to the default values
c
            do while(ifiletemp.lt.nfiles)
               ifiletemp=ifiletemp+1
               isymbol=nsymboluptodata
      write(6,*) ifiletemp,isymbol
               do while(isymbol.lt.nsymboldata)
                  isymbol=isymbol+1
                  iref=ireff(isymbol,ifiletemp)
                  symbol_value(iref)=symbol_default(isymbol)
               end do
            end do

c Now set the specific values given by the .ato file

            il=nsymboluptodata
            if1=0
c      write(6,*) il,if1
            do while(if1.lt.nfiles)
               if1=if1+1
c      write(6,*) if1,il,type(if1),' ',atsymbol(if1)
               symbol_value(il+1)=type(if1)
               symbol_value(il+2)=atsymbol(if1)
               if(symbol_value(il+2).ne.'H  ') then
                  write(symbol_value(il+3),'(i1)') 0
               endif 
               il=il+nsymbolplus
            end do

c Check that output is 1 or 2

         else if(symbol(icount).eq.'output') then

            read(symbol_value(iref),*) itest
            if(itest.lt.1) itest=1
            if(itest.gt.2) itest=2
            write(symbol_value(iref),'(i1)') itest

c Get the number of samples

         else if(symbol(icount).eq.'nsamples') then

            read(symbol_value(iref),*) nsamples
            if(nsamples.gt.3) nsamples=3
            if(nsamples.lt.1) nsamples=1
            write(symbol_value(iref),'(i1)') nsamples

c xfrac can only be set if number of samples is 3

            if(nsamples.eq.3) then

               iflag(icount+1)=0

            else

               iflag(icount+1)=1
               symbol_value(iref+1)=symbol_default(icount+1)

            endif

c nsamples is 2 or 3, then isotope substitution is allowed, otherwise not               

            if(nsamples.lt.2) then

c Set the flags on the individual components to prevent isotope substitution flags being set

               iflag(nsymboluptodata+5)=1
               iflag(nsymboluptodata+6)=1
               il=nsymboluptodata
               do if1=1,nfiles

                  symbol_value(il+5)=symbol_default(nsymboluptodata+5)
                  symbol_value(il+6)=symbol_default(nsymboluptodata+6)
                  il=il+nsymbolplus

               end do

            else

               iflag(nsymboluptodata+5)=0
               iflag(nsymboluptodata+6)=0

            endif

c Make sure normtot has an allowed value

         else if(symbol(icount).eq.'normtot') then

            read(symbol_value(iref),*) itest
            if(itest.lt.0) itest=0
            if(itest.gt.2) itest=2
            write(symbol_value(iref),*) itest

c Make sure normdif has an allowed value

         else if(symbol(icount).eq.'normdif') then

            read(symbol_value(iref),*) itest
            if(itest.lt.0) itest=0
            if(itest.gt.1) itest=1
            write(symbol_value(iref),*) itest

c Only 'H' atoms are allowed to exchange

         else if(symbol(icount).eq.'iexchange') then

            if(symbol_value(iref-1).ne.'H  ') then

               write(symbol_value(iref),'(i1)') 0

            endif
            read(symbol_value(iref),*) itest
            if(itest.lt.0) itest=0
            if(itest.gt.1) itest=1
            write(symbol_value(iref),*) itest

c For abundances need to ensure there are an even number entries

         else if(symbol(icount).eq.'abundances'.or.symbol(icount).eq.'abund2') then

            sometext=symbol_value(iref)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            do while(nwords.ne.2*(nwords/2))
               write(6,252,advance='no')
252   format(/'check_symbol> Need to specify mass numbers and abundances in pairs: ')
               read(5,'(a)') sometext
               call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            end do
            symbol_value(iref)=sometext(ncf(1):ncl(nwords))

         endif
      endif

c Changes specific to xwts

      if(symbol_setup.eq.'xwts') then

c         write(6,*) iref,icount,ifilecount

c If fnameato has been set then get the number of components

         if(symbol(icount).eq.'fnameato'.and.symbol_value(iref).ne.symbol_default(iref)) then

            nfilet=ntype
            if(nfilet.lt.10) then
               write(symbol_value(nsymboluptodata),'(i1)') nfilet
            else if(nfilet.lt.100) then
               write(symbol_value(nsymboluptodata),'(i2)') nfilet
            else
               write(symbol_value(nsymboluptodata),'(i3)') nfilet
            endif

            ifiletemp=nfiles
            nfiles=nfilet
c      write(6,*) nfilet,ifiletemp,nfiles
c
c if new number of files .gt. old number of files
c set the new values of symbol_value to the default values
c
            do while(ifiletemp.lt.nfiles)
               ifiletemp=ifiletemp+1
               isymbol=nsymboluptodata
c     write(6,*) ifiletemp,isymbol
               do while(isymbol.lt.nsymboldata)
                  isymbol=isymbol+1
                  iref=ireff(isymbol,ifiletemp)
                  symbol_value(iref)=symbol_default(isymbol)
               end do
            end do

c Now set the specific values given by the .ato file

            il=nsymboluptodata
            if1=0
c      write(6,*) il,if1
            do while(if1.lt.nfiles)
               if1=if1+1
c      write(6,*) if1,il,type(if1),' ',atsymbol(if1)
               symbol_value(il+1)=type(if1)
               if(symbol_value(il+2).eq.symbol_default(nsymboluptodata+2)) symbol_value(il+2)=atsymbol(if1)
              il=il+nsymbolplus
            end do

c Check that output is 1 or 2

         else if(symbol(icount).eq.'output') then

            read(symbol_value(iref),*) itest
            if(itest.lt.1) itest=1
            if(itest.gt.2) itest=2
            write(symbol_value(iref),'(i1)') itest

c Normtot cannot be less than 0 or greater than 2

         else if(symbol(icount).eq.'normtot') then

            read(symbol_value(iref),*) itest            
      write(6,*) itest
            if(itest.lt.0) itest=0
            if(itest.gt.2) itest=2
            write(symbol_value(iref),'(i1)') itest

         endif

      endif
c
c changes specific to changeato
c

      if(symbol_setup.eq.'changeato') then

            call change_symbol_changeato(icount,ifilecount,ierr)

      endif
      
! Changes specific to bragginten

c Changes specific to xwts

      if(symbol_setup.eq.'bragginten') then

c         write(6,*) iref,icount,ifilecount

c If fnameato has been set then get the number of components

         if(symbol(icount).eq.'fnameato'.and.symbol_value(iref).ne.symbol_default(iref)) then

            nfilet=ntype
            if(nfilet.lt.10) then
               write(symbol_value(nsymboluptodata),'(i1)') nfilet
            else if(nfilet.lt.100) then
               write(symbol_value(nsymboluptodata),'(i2)') nfilet
            else
               write(symbol_value(nsymboluptodata),'(i3)') nfilet
            endif

            ifiletemp=nfiles
            nfiles=nfilet
c      write(6,*) nfilet,ifiletemp,nfiles
c
c if new number of files .gt. old number of files
c set the new values of symbol_value to the default values
c
            do while(ifiletemp.lt.nfiles)
               ifiletemp=ifiletemp+1
               isymbol=nsymboluptodata
c     write(6,*) ifiletemp,isymbol
               do while(isymbol.lt.nsymboldata)
                  isymbol=isymbol+1
                  iref=ireff(isymbol,ifiletemp)
                  symbol_value(iref)=symbol_default(isymbol)
               end do
            end do

c Now set the specific values given by the .ato file

            il=nsymboluptodata
            if1=0
c      write(6,*) il,if1
            do while(if1.lt.nfiles)
               if1=if1+1
c      write(6,*) if1,il,type(if1),' ',atsymbol(if1)
               symbol_value(il+1)=type(if1)
               !Q type atoms cannot scatter, therefore set scattering length for these atoms to zero
               if(qtypea(if1).gt.0) then
                  symbol_value(il+2)='0.0'
               end if
               il=il+nsymbolplus
            end do
         end if
      end if



      return
      end

      subroutine save_files_changeato(nout)
      include 'epsr_dimension.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      integer*4 nout
c
c close the opened file as it is not needed
c
      close(nout)
c
c write out the revised .ato file
c
      call moloutc(fnameinp)
      return
      end

      subroutine write_symbols_changeato(nfirst,nlast)
c
c writes the values read in from the .ato file to the appropriate symbol_value
c
      use molecom_mod
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'

      integer*4 icount,ifilecount
c
c ensure nfirst and nlast are within the allowed range
c
      if(nfirst.lt.1) nfirst=1
      if(nfirst.le.nsymboldata) then
         ifilecount=1
      else
         ifilecount=nfiles
      endif
      if(nlast.gt.nsymbol) nlast=nsymbol
c
c do all the symbols that can be defined from the .ato file - the others
c are set to their default values
c
      icount=nfirst-1
      do while (icount.le.nlast)
            if(icount.ge.nsymboldata.and.ifilecount.lt.nfiles) then
               ifilecount=ifilecount+1
               icount=nsymboluptodata
            endif
            icount=icount+1
            iref=ireff(icount,ifilecount)
            if(symbol(icount).eq.'bond') then
               symbol_value(iref)=symbol_default(iref)
            else if(symbol(icount).eq.'label') then
               symbol_value(iref)=symbol_default(iref)
            else if(symbol(icount).eq.'rho') then
               write(symbol_value(iref),*) rho
            else if(symbol(icount).eq.'temp') then
               write(symbol_value(iref),*) simtk
            else if(symbol(icount).eq.'tol') then
               write(symbol_value(iref),*) tolout
            else if(symbol(icount).eq.'stepmi') then
               write(symbol_value(iref),*) stepmi
            else if(symbol(icount).eq.'stepri') then
               write(symbol_value(iref),*) stepri
            else if(symbol(icount).eq.'stepr') then
               write(symbol_value(iref),*) stepr
            else if(symbol(icount).eq.'stepm') then
               write(symbol_value(iref),*) stepm
            else if(symbol(icount).eq.'vibtemp') then
               write(symbol_value(iref),*) tempb
            else if(symbol(icount).eq.'angtemp') then
               write(symbol_value(iref),*) tempangreq
            else if(symbol(icount).eq.'dihtemp') then
               write(symbol_value(iref),*) tempdihreq
            else if(symbol(icount).eq.'ecore') then
               write(symbol_value(iref),*) ecore
            else if(symbol(icount).eq.'dcore') then
               write(symbol_value(iref),*) dcore
            else if(symbol(icount).eq.'box') then
               if(cubic) then
                  write(symbol_value(iref),*) 'cubic'
               else
                  write(symbol_value(iref),*) 'noncubic'
               endif
            else if(symbol(icount).eq.'addbond') then
               symbol_value(iref)=symbol_default(iref)
            else if(symbol(icount).eq.'removebond') then
               symbol_value(iref)=symbol_default(iref)
            else if(symbol(icount).eq.'natomtypes') then
               write(symbol_value(iref),*) ntype
               nfiles=ntype
            else if(symbol(icount).eq.'atlabel') then
               write(symbol_value(iref),*) type(ifilecount)
            else if(symbol(icount).eq.'epsilon') then
               write(symbol_value(iref),*) epscor(ifilecount)
            else if(symbol(icount).eq.'sigma') then
               write(symbol_value(iref),*) sigcor(ifilecount)
            else if(symbol(icount).eq.'atmass') then
               write(symbol_value(iref),*) amassa(ifilecount)
            else if(symbol(icount).eq.'charge') then
               write(symbol_value(iref),*) charge(ifilecount)
            else if(symbol(icount).eq.'qradius') then
               write(symbol_value(iref),*) qradius(ifilecount)
            else if(symbol(icount).eq.'atsymbol') then
               write(symbol_value(iref),*) atsymbol(ifilecount)
            else if(symbol(icount).eq.'iexchange') then
               write(symbol_value(iref),*) iexchange(ifilecount)
            else if(symbol(icount).eq.'natomsph') then
               write(symbol_value(iref),*) natomsphere(ifilecount)
            else if(symbol(icount).eq.'radminsph') then
               write(symbol_value(iref),*) radminsphere(ifilecount)
            else if(symbol(icount).eq.'radmaxsph') then
               write(symbol_value(iref),*) radmaxsphere(ifilecount)
            else if(symbol(icount).eq.'cgmolfac') then
               write(symbol_value(iref),*) pseudofacmol(ifilecount)
            endif

      end do

      return
      end

      subroutine change_symbol_changeato(icount,ifilecount,ierr)
c
c resets the value of a specified parameter in the .ato file if requested 
c by calling program
c
      use molecom_mod
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      integer*4 ifilecount,ierr,icount,itemp
      real :: alpha,beta,gam,newalpha,newbeta,newgam
      real*4 rtemp,afac,bfac,cfac,newphib,newthetac,newphic
      integer*4 massno
      character*3 ctemp
      character(len=256) :: tokens
      integer :: ntf(nworddim),ntl(nworddim),ntokens

c
c initial error signal
c
      ierr=0
            
      iref=ireff(icount,ifilecount)

      if(symbol(icount).eq.'bond') then

         if(symbol_value(iref).eq.'y') call change_bond()
         symbol_value(iref)=symbol_save

      else if(symbol(icount).eq.'label') then

         if(symbol_value(iref).eq.'y') call change_label()
         symbol_value(iref)=symbol_save

      else if(symbol(icount).eq.'tether') then

! Get the command and break into tokens
         tokens=symbol_value(iref)
         call parse(tokens,ntokens,ntf,ntl,11,nstsize,nworddim)
         if(tokens(ntf(1):ntl(1)).eq.'y') call change_tether(tokens,ntokens,ntf,ntl)
         symbol_value(iref)=symbol_save

      else if(symbol(icount).eq.'remove') then

         if(symbol_value(iref).eq.'y') call remove_molecule_types()
         symbol_value(iref)=symbol_save

      else if(symbol(icount).eq.'rho') then

         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0.and.rtemp.gt.0.0) then
            rhoold=natomrcg/(alatx*blaty*clatz)
            rho=rtemp
            vol=rhoold/rho
            dmaxfac=vol**(1.0/3.0)
            afac=dmaxfac
            bfac=dmaxfac
            cfac=dmaxfac
            newphib=phib
            newthetac=thetac
            newphic=phic
            call change_boxsize(afac,bfac,cfac,newphib,newthetac,newphic)
c      write(6,*) rhoold,rho,dmaxfac,alat,blat,clat

         endif

      else if(symbol(icount).eq.'temp') then

         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            simtk=rtemp
         endif

      else if(symbol(icount).eq.'tol') then

         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            tolout=rtemp
         endif

      else if(symbol(icount).eq.'stepmi') then

         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            stepmi=rtemp
         endif

      else if(symbol(icount).eq.'stepri') then

         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            stepri=rtemp
         endif

      else if(symbol(icount).eq.'stepr') then

         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            stepr=rtemp
         endif

      else if(symbol(icount).eq.'stepm') then

         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            stepm=rtemp
         endif

      else if(symbol(icount).eq.'vibtemp') then

         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            tempb=rtemp
         endif

      else if(symbol(icount).eq.'angtemp') then

         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            tempangreq=rtemp
         endif

      else if(symbol(icount).eq.'dihtemp') then

         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            tempdihreq=rtemp
         endif

      else if(symbol(icount).eq.'ecore') then

         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            ecore=rtemp
         endif

      else if(symbol(icount).eq.'dcore') then

         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            dcore=rtemp
         endif

      else if(symbol(icount).eq.'addbond') then

         if(symbol_value(iref).eq.'y') call add_bond()
         symbol_value(iref)=symbol_save

      else if(symbol(icount).eq.'removebond') then

         if(symbol_value(iref).eq.'y') call remove_bond()
         symbol_value(iref)=symbol_save

      else if(symbol(icount).eq.'box') then

! Get the command and break into tokens
         tokens=symbol_value(iref)
         call parse(tokens,ntokens,ntf,ntl,11,nstsize,nworddim)
c Only need to do anything if symbol_value is 'noncubic' or if cubic = .false. and symbol_value is 'cubic'
c Note that to change shape of a cubic box you need to type 'noncubic'.
c Otherwise there is no change in box shape, so nothing to do.
         if((cubic.and.tokens(ntf(1):ntl(1)).eq.'noncubic').or.tokens(ntf(1):ntl(1)).eq.'change') then
c First deal with case of going from cubic box to non-cubic box, or wanting to change a non-cubic box.
            write(6,102) alat,blat,clat
102   format(/'changeato> Current a, b, and c axes: ',3(1x,f10.5))
c Get the new axes
            if(ntokens.gt.3) then
               sometext=tokens(ntf(2):ntl(4))
            else
               write(6,100,advance='no')
100   format(/'changeato> Type a factor for each axis, or one factor for all axes,'
     */' (Use -ve for actual values, or press <Enter> to leave unchanged: ')
               sometext=' '
               read(5,'(a)',iostat=ierr) sometext
            end if
            call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
            if(nwords.eq.1) then
               read(sometext(ncf(1):ncl(1)),*) afac
               bfac=afac
               cfac=afac
            else if(nwords.eq.3) then
               read(sometext(ncf(1):ncl(3)),*) afac,bfac,cfac
            else
               afac=1.0
               bfac=1.0
               cfac=1.0
            endif
            if(afac.lt.0.0) afac=abs(afac)/alat
            if(bfac.lt.0.0) bfac=abs(bfac)/blat
            if(cfac.lt.0.0) cfac=abs(cfac)/clat
            
!            write(6,'(a,3(1x,f10.5))') 'changeato> Current polar angles phib, thetac and phic: ',phib,thetac,phic
            call convertphibthetacphic(phib,thetac,phic,alpha,beta,gam)
            write(6,'(a,3(1x,f10.5))') 'changeato> Current crystallographic angles alpha, beta, gamma: ',alpha,beta,gam
            if(ntokens.gt.6) then
               sometext=tokens(ntf(5):ntl(7))
            else
               write(6,101,advance='no')
101   format(/'changeato> Type new alpha, beta, gamma values or press <Enter> to leave unchanged: ')
               sometext=' '
               read(5,'(a)',iostat=ierr) sometext
            end if
            call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
            if(nwords.eq.3) then
               read(sometext(ncf(1):ncl(3)),*) newalpha,newbeta,newgam
               call convertalphabetagamma(newalpha,newbeta,newgam,newphib,newthetac,newphic)
               write(6,'(a,3(1x,f10.5))') 'changeato> New polar angles phib, thetac and phic: ',newphib,newthetac,newphic
            else
               newphib=phib
               newthetac=thetac
               newphic=phic
            end if

            call change_boxsize(afac,bfac,cfac,newphib,newthetac,newphic)
            
         else if((.not.cubic.and.tokens(ntf(1):ntl(1)).eq.'cubic')) then

c Now for the case when going from non-cubic to cubic

c Get the volume of the current simulation box

            volold=alatx*blaty*clatz

c Get the new axes

            alat=1.0
            blat=1.0
            clat=1.0
            phib=90.0
            thetac=0.0
            phic=0.0

c set up lattice parameters for applying the minimum image convention

            alatx=alat
            blatx=blat*cosd(phib)
            blaty=blat*sind(phib)
            clatx=clat*sind(thetac)*cosd(phic)
            claty=clat*sind(thetac)*sind(phic)
            clatz=clat*cosd(thetac)

c Get the volume of this new box

            vol=volold/(alatx*blaty*clatz)

c Scale the box dimensions to give the save volume as previously

            dmaxfac=vol**(1.0/3.0)
            alat=alat*dmaxfac
            alatx=alatx*dmaxfac
            blat=blat*dmaxfac
            blatx=blatx*dmaxfac
            blaty=blaty*dmaxfac
            clat=clat*dmaxfac
            clatx=clatx*dmaxfac
            claty=claty*dmaxfac
            clatz=clatz*dmaxfac
c
c define the Cartesian coordinates of the lattice vectors
c
            alatc(1)=alatx
            alatc(2)=0.0
            alatc(3)=0.0
            blatc(1)=blatx
            blatc(2)=blaty
            blatc(3)=0.0
            clatc(1)=clatx
            clatc(2)=claty
            clatc(3)=clatz
            bxax=blatx/alatx
            cxax=clatx/alatx
            cyby=claty/blaty

c
c calculate the perpendicular distance between the faces of the lattice
c
            perpab=perpd(alatc,blatc,clatc)
            perpbc=perpd(blatc,clatc,alatc)
            perpca=perpd(clatc,alatc,blatc)
            write(6,*) perpab,perpbc,perpca
c
c set up maximum radius to be used for Fourier transforms and to truncate
c the potential. This is defined as the minimum perpendicular distance between
c planes. 
c
c G(r) and potentials will be calculated to rmax, or maximum allowed distance
c based on number of points and step size, whichever is smaller.
c
            dmax=amin1(perpab,perpbc,perpca)
            write(6,*) nmol,dmax
            cell=dmax
            cubic=.true.

         endif

!Calculate the new density and assign the value to the appropriate symbol

         rho=natomr/(alatx*blaty*clatz)
         write(symbol_value(irefrho),*) rho

c Make sure a sensible value is left, in case the typed value was nonsense

         if(cubic) then
            write(symbol_value(iref),*) 'cubic'
         else
            write(symbol_value(iref),*) 'noncubic'
         endif
         
      else if(symbol(icount).eq.'findmols') then
      
         if(symbol_value(iref).eq.'y') call get_molecules()
         symbol_value(iref)=symbol_save

      else if(symbol(icount).eq.'getmolfs') then
      
         if(symbol_value(iref).eq.'y') call make_mol_files()
         symbol_value(iref)=symbol_save

      else if(symbol(icount).eq.'add_w_rot') then
      
! Get the command and break into tokens
         tokens=symbol_value(iref)
         call parse(tokens,ntokens,ntf,ntl,11,nstsize,nworddim)
         if(tokens(ntf(1):ntl(1)).eq.'y') call add_whole_rot(tokens,ntokens,ntf,ntl)
         symbol_value(iref)=symbol_save

      else if(symbol(icount).eq.'axisshift') then
      
! Get the command and break into tokens
         tokens=symbol_value(iref)
         call parse(tokens,ntokens,ntf,ntl,11,nstsize,nworddim)
         if(tokens(ntf(1):ntl(1)).eq.'y') call shift_molecules(tokens,ntokens,ntf,ntl)
         symbol_value(iref)=symbol_save

      else if(symbol(icount).eq.'axisrotate') then
      
! Get the command and break into tokens
         tokens=symbol_value(iref)
         call parse(tokens,ntokens,ntf,ntl,11,nstsize,nworddim)
         if(tokens(ntf(1):ntl(1)).eq.'y') call rotate_molecules(tokens,ntokens,ntf,ntl)
         symbol_value(iref)=symbol_save

      else if(symbol(icount).eq.'alignmole') then
      
! Get the command and break into tokens
         tokens=symbol_value(iref)
         call parse(tokens,ntokens,ntf,ntl,11,nstsize,nworddim)
         if(tokens(ntf(1):ntl(1)).eq.'y') call align_molecules(tokens,ntokens,ntf,ntl)
         symbol_value(iref)=symbol_save

      else if(symbol(icount).eq.'setimage') then
      
! Get the command and break into tokens
         tokens=symbol_value(iref)
         call parse(tokens,ntokens,ntf,ntl,11,nstsize,nworddim)
         if(tokens(ntf(1):ntl(1)).eq.'y') call add_image_atom(.true.,tokens,ntokens,ntf,ntl)
         symbol_value(iref)=symbol_save

      else if(symbol(icount).eq.'unsetimage') then
      
! Get the command and break into tokens
         tokens=symbol_value(iref)
         call parse(tokens,ntokens,ntf,ntl,11,nstsize,nworddim)
         if(tokens(ntf(1):ntl(1)).eq.'y') call add_image_atom(.false.,tokens,ntokens,ntf,ntl)
         symbol_value(iref)=symbol_save

      else if(symbol(icount).eq.'molscalefc') then

! Get the command and break into tokens
         tokens=symbol_value(iref)
         call parse(tokens,ntokens,ntf,ntl,11,nstsize,nworddim)
         if(tokens(ntf(1):ntl(1)).eq.'y') call set_molecule_scale_factors(tokens,ntokens,ntf,ntl)
         symbol_value(iref)=symbol_save

      else if(symbol(icount).eq.'atlabel') then

         read(symbol_value(iref),*,iostat=ierr) ctemp
         if(ierr.eq.0) then
c
c step through all the atoms and make sure any with the old
c type are changed to the new
c
            do ia=1,natom
               if(itype(ia).eq.ifilecount) alabel(ia)=ctemp
            end do
c
c finally change the type itself
c
            type(ifilecount)=ctemp

         endif

      else if(symbol(icount).eq.'epsilon') then
         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            epscor(ifilecount)=rtemp
         endif

      else if(symbol(icount).eq.'sigma') then
         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            sigcor(ifilecount)=rtemp
         endif

      else if(symbol(icount).eq.'atmass') then
         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            amassa(ifilecount)=rtemp
         endif

      else if(symbol(icount).eq.'charge') then
         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            charge(ifilecount)=rtemp
         endif

      else if(symbol(icount).eq.'qradius') then
         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            qradius(ifilecount)=rtemp
         endif

      else if(symbol(icount).eq.'atsymbol') then
         read(symbol_value(iref),*,iostat=ierr) ctemp
         if(ierr.eq.0) then
c
c need to check that this chemical symbol exists in the Periodic Table
c
c
c test to see if this chemical symbol exists in the table
c
            massno=0
            btest=scatt_len(ctemp,massno,ifind)
            if(ifind.eq.0) then
c
c reference integer for atom label
c
               irefl=ireff(nsymboluptodata+1,ifilecount)
               write(6,209) ctemp,symbol_value(irefl)
209   format(/'changeato> Did not recognise chemical symbol ',a3
     *,' for atom type ',a3)
               ierr=1
            else
               atsymbol(ifilecount)=ctemp
            endif
         endif

      else if(symbol(icount).eq.'iexchange') then
         read(symbol_value(iref),*,iostat=ierr) itemp
         if(ierr.eq.0) then
            iexchange(ifilecount)=itemp
         endif

      else if(symbol(icount).eq.'natomsph') then
         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            natomsphere(ifilecount)=rtemp
         endif

      else if(symbol(icount).eq.'radminsph') then
         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            radminsphere(ifilecount)=rtemp
         endif

      else if(symbol(icount).eq.'radmaxsph') then
         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0) then
            radmaxsphere(ifilecount)=rtemp
         endif

      else if(symbol(icount).eq.'cgmolfac') then
         read(symbol_value(iref),*,iostat=ierr) rtemp
         if(ierr.eq.0.and.rtemp.ge.1.0) then
            pseudofacmol(ifilecount)=rtemp
         else
            pseudofacmol(ifilecount)=1.0
         endif

      endif

c
c reset the value of the symbol if an error occurred on trying to set it
c
      if(ierr.ne.0) then
         symbol_value(iref)=symbol_save
         symbol_save=' '
         ierr=0
      endif

      return
      end

      subroutine change_bond()
      include 'epsr_dimension.inc'

      include 'molecom.inc'
      include 'clustercom.inc'
      write(6,31,advance='no')
31    format(/1x,'Give 2 atom type numbers, 2 atom numbers, and new bond: ')
      read(5,*) ncha1,ncha2,nat1,nat2,bonnew
      do while (ncha1.ne.0)
            do ia=1,natom
                  mmm=mref(ia)
                  na1=ia+1-naf(mmm)
                  if(nbond(ia).gt.0.and.(na1.eq.nat1.or.na1.eq.nat2)
     1) then
                        do ib=1,nbond(ia)
                              na2=iabond(ia,ib)+1-naf(mmm)
                              if(itype(ia).eq.ncha1
     1.and.itype(iabond(ia,ib)).eq.ncha2
     1.and.na1.eq.nat1.and.na2.eq.nat2) tbond(ia,ib)=bonnew
                              if(itype(ia).eq.ncha2
     1.and.itype(iabond(ia,ib)).eq.ncha1
     1.and.na1.eq.nat2.and.na2.eq.nat1) tbond(ia,ib)=bonnew
                        end do
                  end if
            end do
            write(6,31)
            read(5,*) ncha1,ncha2,nat1,nat2,bonnew
      end do
      return
      end

      subroutine change_label()

      use molecom_mod
      include 'epsr_dimension.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      character*3 a1,ats
      integer findtype(mcomp),iexch,it,ifound,nthistype,nchange,ichange,itest
      integer, dimension(:), allocatable      ::changeselect,changethis
      real natomsph,rhosph,radminsph,radmaxsph,pseudofacsph
      write(6,35)
35    format(/'change_label> The following atom types are currently defined:-')
      do it=1,ntype
         write(6,*) it,type(it)
      end do
      write(6,'(/a)') 'change_label> The following molecule types are currently defined:-'
      do it=1,nmoltype
         write(6,'(i5,1x,a3)') it,moltype(it)
      end do
      write(6,36,advance='no')
36    format(/'change_label> Give atom type number'
     1,' and molecule type number to change (0 0 to exit): ')
      read(5,*) it,iamol
      if(it.lt.1.or.it.gt.ntype.or.iamol.lt.1.or.iamol.gt.nmoltype) return
      !First find out how many atoms of this type exist in the file
      nthistype=0
      do im=1,nmol
         if(mtype(im).eq.iamol) then
            do ia=naf(im),nal(im)
               if(itype(ia).eq.it) nthistype=nthistype+1
            end do
         end if
      end do
      if(nthistype.le.0) then
         write(6,1784) type(it),iamol
1784     format(/'change_label> Warning! No atoms of type ',a3,' at atom number ',i4)
         return
      else
         write(6,1788) nthistype,type(it),moltype(iamol)
1788     format(/'change_label> There are ',i7,' atoms of type ',a3,' in molecules of type ',a3)
      end if
!Grab the list of atoms to (possibly) change
      if(allocated(changeselect)) deallocate(changeselect)
      allocate(changeselect(nthistype))
      nthistype=0
      do im=1,nmol
         if(mtype(im).eq.iamol) then
            do ia=naf(im),nal(im)
               if(itype(ia).eq.it) then
                  nthistype=nthistype+1
                  changeselect(nthistype)=ia
               end if
            end do
         end if
      end do
      write(6,37,advance='no') type(it)
37    format(/'change_label> Current label for this type is ',a3
     1//,' Give new label (a3): ')
      read(5,38) a1
38    format(a)
      write(6,39,advance='no')
39    format(/'change_label> Give fraction of these atoms to change (1.0 for all): ')
      read(5,*) fracchange
      notochange=min(nthistype,int(fracchange*real(nthistype)))
      if(allocated(changethis)) deallocate(changethis)
      allocate(changethis(notochange))
!Finally pick out the atoms to change at random from this list, remembering to reduce the list once it is changed
      !Get the list of atoms to change
      ichange=0
      icount=0
      do while(ichange.lt.notochange)
         icount=icount+1
         if(icount.gt.nthistype) icount=1
! Pick atom at random
         if((ran1().lt.fracchange.or.fracchange.ge.1.0).and.changeselect(icount).gt.0) then
            ichange=ichange+1
            changethis(ichange)=changeselect(icount)
            changeselect(icount)=0 !Set this to zero to indicate this atom has already been selected
         end if
      end do
      nchange=ichange
      if(nchange.le.0) then
         write(6,1811) 
1811     format(/'change_label> Warning: No valid atoms to change labels found!!'/)
         return
      else
         write(6,1814) nchange 
1814     format(/'change_label> Changing labels on ',i7,' atoms'/)
      end if
      ichange=0
      do while(ichange.lt.nchange)
         ichange=ichange+1
         ia=changethis(ichange)
         im=mref(ia)
         iaref=ia-naf(im)+1
         epsit=epscor(it)
         sigit=sigcor(it)
         amass=amassa(it)
         qt=charge(it)
         rqt=qradius(it)
         ats=atsymbol(it)
         iexch=iexchange(it)
         alabel(ia)=a1
         natomsph=natomsphere(it)
         rhosph=rhosphere(it)
         radminsph=radminsphere(it)
         radmaxsph=radmaxsphere(it)
         pseudofacsph=pseudofacmol(it)
         call ftype(ia)
         call ctype(a1,ats,iexch,epsit,sigit,amass,qt,rqt,natomsph,rhosph,radminsph,radmaxsph,pseudofacsph)
      end do
      !Clean up the type lists
      call cleanup_types()
c
c now check that all symbol_values have been redefined correctly for any type changes
c that have been made
c
      call write_symbols_changeato(nsymboluptodata,nsymbol) 

      return
      end

      subroutine change_tether(tokens,ntokens,ntf,ntl)

!Tethers specified molecule types

         include 'epsr_dimension.inc'
         include 'commands.inc'
         include 'symbols.inc'
         include 'molecom.inc'
         include 'clustercom.inc'
         integer it,moltyperead,moltetheredatomread
         logical moltetheredread
         character(len=*), intent(in) :: tokens
         integer, intent(in) :: ntokens
         integer, dimension(*), intent(in) :: ntf,ntl

!Make a list of all current molecule types and their tethered status

         write(6,100)
100   format(/1x,'molecule type number',2x,'label',2x,'tethered status (T or F)',2x,'tethered atom'/)
         do it=1,nmoltype
            write(6,*) it,moltype(it),moltethered(it),moltetheredatom(it)
101   format(10x,i4,10x,a3,5x,l1,5x,i5)
         end do
         if(ntokens.gt.4) then
            sometext=tokens(ntf(2):ntl(4))
         else
            write(6,102)
102   format(/'change_tether> Give type number, new status (T or F), and tethering atom (0 for C.M.) ( to quit).')
            read(5,'(a)') sometext
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         end if
         read(sometext(ncf(1):ncl(1)),*,iostat=ierr) moltyperead
         if(ierr.eq.0.and.nwords.gt.1) then
            moltetheredread=sometext(ncf(2):ncl(2)).eq.'T'.or.sometext(ncf(2):ncl(2)).eq.'t'
            if(moltetheredread.and.nwords.gt.2) then
               read(sometext(ncf(3):ncl(3)),*) moltetheredatomread
            else
               moltetheredatomread=0
            end if
            !Apply the new values
            do im=1,nmol
               if(mtype(im).eq.moltyperead) then
! If this molecule is not tethered and we are making it tethered we need to set the tethering coordinates
                  if(.not.tethered(im).and.moltetheredread) then
                     phix(im)=xadd(im)
                     phiy(im)=yadd(im)
                     phiz(im)=zadd(im)
                  end if
                  tethered(im)=moltetheredread
                  tetheredatom(im)=moltetheredatomread
                  !If not tethered then the thethering coordinates should be set to zero
                  if(.not.tethered(im)) then
                     phix(im)=0.0
                     phiy(im)=0.0
                     phiz(im)=0.0
                     tetheredatom(im)=0
                  end if
               endif
            end do
            moltethered(moltyperead)=moltetheredread
            moltetheredatom(moltyperead)=moltetheredatomread
         end if
c
c now check that all symbol_values have been redefined correctly for any type changes
c that have been made
c
         call write_symbols_changeato(nsymboluptodata,nsymbol) 
         return
      end

      subroutine shift_molecules(tokens,ntokens,ntf,ntl)

!Applies a constant shift along a specified axis to all molecules of a specified type
!

         include 'epsr_dimension.inc'
         include 'commands.inc'
         include 'symbols.inc'
         include 'molecom.inc'
         include 'clustercom.inc'
         integer it,moltyperead
         character*1 axisread
         real axisshift
         character(len=*), intent(in) :: tokens
         integer, intent(in) :: ntokens
         integer, dimension(*), intent(in) :: ntf,ntl

!Make a list of all current molecule types and their tethered status

         write(6,100)
100   format(/1x,'molecule type number',2x,'label'/)
         do it=1,nmoltype
            write(6,101) it,moltype(it)
101   format(10x,i4,10x,a3)
         end do
         if(ntokens.gt.3) then
            sometext=tokens(ntf(2):ntl(4))
         else
            write(6,102,advance='no')
102   format(/'shift_molecules> Give type number, axis to shift along (x, y, or z),'/
     *,' and amount to shift. (q to quit): ')
            read(5,'(a)') sometext
         end if
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         read(sometext(ncf(1):ncl(1)),*,iostat=ierr) moltyperead
         if(ierr.eq.0.and.nwords.gt.2) then
            axisread=sometext(ncf(2):ncl(2))
            read(sometext(ncf(3):ncl(3)),*) axisshift
            !Apply the new values
            do im=1,nmol
               if(mtype(im).eq.moltyperead) then
                  if(axisread.eq.'x') xadd(im)=xadd(im)+axisshift
                  if(axisread.eq.'y') yadd(im)=yadd(im)+axisshift
                  if(axisread.eq.'z') zadd(im)=zadd(im)+axisshift
                  if(tethered(im)) then
                     if(axisread.eq.'x') phix(im)=phix(im)+axisshift
                     if(axisread.eq.'y') phiy(im)=phiy(im)+axisshift
                     if(axisread.eq.'z') phiz(im)=phiz(im)+axisshift
                  end if
               endif
            end do
         end if
c
c now check that all symbol_values have been redefined correctly for any type changes
c that have been made
c
         call write_symbols_changeato(nsymboluptodata,nsymbol) 
         return
      end

      subroutine rotate_molecules(tokens,ntokens,ntf,ntl)

!Applies a rotation of a molecule about a specified axis to all molecules of a specified type
!

         include 'epsr_dimension.inc'
         include 'commands.inc'
         include 'symbols.inc'
         include 'molecom.inc'
         include 'clustercom.inc'
         integer it,molnumberread
         character*3 moltyperead
         character*1 axisread
         real axisshift,pi,degconv
         character(len=*), intent(in) :: tokens
         integer, intent(in) :: ntokens
         integer, dimension(*), intent(in) :: ntf,ntl

         pi=4.0*atan(1.0)
         degconv=pi/180.0

!Make a list of all current molecule types and their tethered status

         write(6,100) nmol
100   format(/'rotate_molecules> Number of molecules: ',i0/,'molecule types:- '/)
         do it=1,nmoltype
            write(6,101) it,moltype(it)
101   format(10x,i4,10x,a3)
         end do
         if(ntokens.gt.3) then
            sometext=tokens(ntf(2):ntl(4))
         else
            write(6,102,advance='no')
102   format(/'rotate_molecules> Give molecule number or type, axis to rotate about (x, y, or z),'/
     *,' and amount to rotate [deg.]. (q to quit): ')
            read(5,'(a)') sometext
         end if
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         read(sometext(ncf(1):ncl(1)),*,iostat=ierr) molnumberread
         if(ierr.ne.0) then !Assume this is a molecule type, not a number
            moltyperead=sometext(ncf(1):ncl(1))
            molnumberread=0
         else
            moltyperead=' '
         end if
         if(ierr.eq.0.and.nwords.gt.2) then
            axisread=sometext(ncf(2):ncl(2))
            read(sometext(ncf(3):ncl(3)),*) axisshift
            phi=degconv*axisshift
            cphi=cos(phi)
            sphi=sin(phi)
            !Apply the new values
            do im=1,nmol
               if(moltype(mtype(im)).eq.moltyperead.or.im.eq.molnumberread) then
                  if(axisread.eq.'z') then
c rotate about z axis
                     do ia=naf(im),nal(im)
                        xasave=xa(ia)*cphi-ya(ia)*sphi
                        ya(ia)=xa(ia)*sphi+ya(ia)*cphi
                        xa(ia)=xasave
                     end do
                  else if(axisread.eq.'y') then
c rotate about y-axis
                     do ia=naf(im),nal(im)
                        zasave=za(ia)*cphi-xa(ia)*sphi
                        xa(ia)=za(ia)*sphi+xa(ia)*cphi
                        za(ia)=zasave
                     end do
                  else if(axisread.eq.'x') then
c rotate about x-axis
                     do ia=naf(im),nal(im)
                        yasave=ya(ia)*cphi-za(ia)*sphi
                        za(ia)=ya(ia)*sphi+za(ia)*cphi
                        ya(ia)=yasave
                     end do
                  end if
               end if
            end do
         end if
c
c now check that all symbol_values have been redefined correctly for any type changes
c that have been made
c
         call write_symbols_changeato(nsymboluptodata,nsymbol) 
         return
      end

      subroutine align_molecules(tokens,ntokens,ntf,ntl)

!Rotates all molecules of a specified type to align with a specified axis.
!

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      integer it,moltyperead,axisatom1,axisatom2
      character*1 axisread
      real labaxis(3),molaxis(3),labmod,molmod
      character(len=*), intent(in) :: tokens
      integer, intent(in) :: ntokens
      integer, dimension(*), intent(in) :: ntf,ntl

!Make a list of all current molecule types
      write(6,100)
100   format(/1x,'molecule type number',2x,'label'/)
      do it=1,nmoltype
         write(6,101) it,moltype(it)
101   format(10x,i4,10x,a3)
      end do
      if(ntokens.gt.4) then
         sometext=tokens(ntf(2):ntl(5))
      else
         write(6,102,advance='no')
102   format(/'align_molecules> Give type number, axis to align to (x,y,z,a,b,c),'/
     *,' and 2 atom numbers to define molecular axis (q to quit): ')
         read(5,'(a)') sometext
      end if
      call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
      read(sometext(ncf(1):ncl(1)),*,iostat=ierr) moltyperead
      if(ierr.eq.0.and.nwords.gt.3) then
         axisread=sometext(ncf(2):ncl(2))
         read(sometext(ncf(3):ncl(4)),*) axisatom1,axisatom2
! Set up lab axis vector
         if(axisread.eq.'x'.or.axisread.eq.'a') then
            labaxis(1)=1.0
            labaxis(2)=0.0
            labaxis(3)=0.0
         else if(axisread.eq.'y') then
            labaxis(1)=0.0
            labaxis(2)=1.0
            labaxis(3)=0.0
         else if(axisread.eq.'z') then 
            labaxis(1)=0.0
            labaxis(2)=0.0
            labaxis(3)=1.0
         else if(axisread.eq.'b') then
            labaxis(1)=blatx
            labaxis(2)=blaty
            labaxis(3)=0.0
         else if(axisread.eq.'c') then 
            labaxis(1)=clatx
            labaxis(2)=claty
            labaxis(3)=clatz
         end if
! Ensure this vector is a unit vector
         call normvec(labaxis)
! Rotate specified molecules
         do im=1,nmol
            if(mtype(im).eq.moltyperead) then
               call rotate_mol(im,labaxis,axisatom1,axisatom2)
            endif
         end do
      end if
c
c now check that all symbol_values have been redefined correctly for any type changes
c that have been made
c
      call write_symbols_changeato(nsymboluptodata,nsymbol) 
      return
      end

      
      subroutine rotate_mol(im,labaxis,axisatom1,axisatom2)
      
      ! Rotates a molecule so the axis going from axisatom1 to axisatom2 coincides with 
      ! labaxis. Labaxis is converted to a unit vector if it is not so already.      
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      
      integer im,axisatom1,axisatom2,nax1,nax2
      real labaxis(3),molaxis(3)

c pvector is the vector about which the rotation occurs

      real pvector(3)

c xvector is the vector from centre of rotation to an individual atom
c xpvector is the component of this vector along the axis of rotation
c xrvector is the component of this vector perpendicular to the axis of rotation
c xrpvector is the vector perpendicular to both axis of rotation and xrvector

      real xvector(3),xpvector(3),xrvector(3),xrpvector(3),cphi,sphi

      pi=3.141592653
      call normvec(labaxis)      
! First check specified axis atoms exist in this molecule.
! If not nothing is done.

!      write(6,*) im,axisatom1,axisatom2
      if(axisatom1.lt.1.or.axisatom1.gt.nato(im)) return
      if(axisatom2.lt.1.or.axisatom2.gt.nato(im)) return
      if(axisatom1.eq.axisatom2) return
      nax1=axisatom1+naf(im)-1
      nax2=axisatom2+naf(im)-1
      !Axis defined by axisatom1 and axisatom2
      molaxis(1)=xa(nax2)-xa(nax1)
      molaxis(2)=ya(nax2)-ya(nax1)
      molaxis(3)=za(nax2)-za(nax1)
      call normvec(molaxis)
!      write(6,*) im,(molaxis(j),j=1,3),(labaxis(j),j=1,3)
      !Cosine of required rotation
      cphi=scapro(labaxis,molaxis)
      sphi=sqrt(1.0-cphi*cphi)
      write(6,*) im,cphi,sphi
      !Axis of rotation is at right angles to the plane defined by molaxis and lab axis
      call vecpro(molaxis,labaxis,pvector)
      rab2=scapro(pvector,pvector)
      if(rab2.gt.0) then
         rab=sqrt(rab2)
c Unit vector along axis of rotation
         do i=1,3
            pvector(i)=pvector(i)/rab
         end do
c step through atoms in group
         do ia=1,nato(im)
c            if(ia.ne.axisatom1.and.ia.ne.axisatom2) then
               nat=ia+naf(im)-1
c first calculate vector from centre of rotation to each atom in group
               xvector(1)=xa(nat)-xa(nax1)
               xvector(2)=ya(nat)-ya(nax1)
               xvector(3)=za(nat)-za(nax1)
c Vector parallel to axis of rotation
               frac=scapro(pvector,xvector)
               do i=1,3
                  xpvector(i)=frac*pvector(i)
               enddo
c Vector perpendicular to axis of rotation
               do i=1,3
                  xrvector(i)=xvector(i)-xpvector(i)
               enddo
c Vector perpendicular to plane defined by axis of rotation and xrvector
               call vecpro(pvector,xrvector,xrpvector)
c Generate rotated coordinates. The vector component parallel to the axis of rotation
c is not rotated
               xa(nat)=xa(nax1)+xpvector(1)
               ya(nat)=ya(nax1)+xpvector(2)
               za(nat)=za(nax1)+xpvector(3)
               xa(nat)=xa(nat)+xrvector(1)*cphi+xrpvector(1)*sphi
               ya(nat)=ya(nat)+xrvector(2)*cphi+xrpvector(2)*sphi
               za(nat)=za(nat)+xrvector(3)*cphi+xrpvector(3)*sphi
c            end if
         end do
      endif
      return
      end

      subroutine add_bond()
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      write(6,41,advance='no')
41    format(/1x,'Give 2 atom type numbers, 2 atom numbers,'
     1,' and bond length: ')
      read(5,*) ncha1,ncha2,nat1,nat2,bonnew
      do while (ncha1.gt.0)
            do ia=1,natom
                  mmm=mref(ia)
                  na1=ia+1-naf(mmm)
                  if(itype(ia).eq.ncha1.and.na1.eq.nat1) then
                        nbond(ia)=nbond(ia)+1
                        ib=nbond(ia)
                        iref=naf(mmm)+nat2-1
                        iabond(ia,ib)=iref
                        tbond(ia,ib)=bonnew
                  endif
                  if (itype(ia).eq.ncha2.and.na1.eq.nat2) then
                        nbond(ia)=nbond(ia)+1
                        ib=nbond(ia)
                        iref=naf(mmm)+nat1-1
                        iabond(ia,ib)=iref
                        tbond(ia,ib)=bonnew
                  end if
            end do
            write(6,41)
            read(5,*) ncha1,ncha2,nat1,nat2,bonnew
      end do
      return
      end

      subroutine remove_bond()
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      write(6,42)
42    format(1x,'Give 2 atom type numbers and 2 atom numbers'
     *,' of bond to delete')
      read(5,*) ncha1,ncha2,nat1,nat2
      do while (ncha1.gt.0)
            do ia=1,natom
                  mmm=mref(ia)
                  na1=ia-naf(mmm)+1
                  if(itype(ia).eq.ncha1.and.na1.eq.nat1
     *.and.nbond(ia).gt.0) then
                      do ib=1,nbond(ia)
                        ibsav=ib
                        ibref=iabond(ia,ib)
                        na2=ibref-naf(mmm)+1
                        if(itype(ibref).eq.ncha2
     *.and.na2.eq.nat2) go to 43
                      end do
                      go to 44
43                    continue
                      nbond(ia)=nbond(ia)-1
                      if(nbond(ia).ge.ibsav) then
                        do ib=ibsav,nbond(ia)
                            ibref=ib+1
                            iabond(ia,ib)=iabond(ia,ibref)
                            tbond(ia,ib)=tbond(ia,ibref)
                        end do
                      endif
44                    continue
                  else if(itype(ia).eq.ncha2.and.na1.eq.nat2
     *.and.nbond(ia).gt.0) then
                      do ib=1,nbond(ia)
                        ibsav=ib
                        ibref=iabond(ia,ib)
                        na2=ibref-naf(mmm)+1
                        if(itype(ibref).eq.ncha1
     *.and.na2.eq.nat1) go to 45
                      end do
                      go to 46
45                    continue
                      nbond(ia)=nbond(ia)-1
                      if(nbond(ia).ge.ibsav) then
                        do ib=ibsav,nbond(ia)
                            ibref=ib+1
                            iabond(ia,ib)=iabond(ia,ibref)
                            tbond(ia,ib)=tbond(ia,ibref)
                        end do
                      endif
46                    continue
                  end if
            end do
            write(6,42)
            read(5,*) ncha1,ncha2,nat1,nat2
      end do
      return
      end
      
      subroutine remove_molecule_types()
      !Removes a specified type of molecule from the ato file.
      
         include 'epsr_dimension.inc'
         include 'commands.inc'
         include 'symbols.inc'
         include 'molecom.inc'
         include 'clustercom.inc'
         integer it,mtyperead,ierr,nremove,ntoremove,findtype(mcomp)
         logical readagain
         character(len=3) :: moltyperead

!Determine the number of molecules of each type

      do it=1,nmoltype
         findtype(it)=0
      end do
      do im=1,nmol
         imoltype=0
         ifound=0
         do while(imoltype.lt.nmoltype.and.ifound.eq.0)
            imoltype=imoltype+1
            if(alabel(naf(im)).eq.moltype(imoltype)) then
               ifound=imoltype
            end if
         end do
         if(ifound.gt.0) then
            findtype(imoltype)=findtype(imoltype)+1
         endif
      end do

!Make a list of all current molecule types and their tethered status

         readagain=.true.
         do while(readagain)
            write(6,100)
100   format(/1x,'molecule type number',2x,'label',2x,'number',2x,'tethered status (T or F)'/)
            do it=1,nmoltype
               write(6,101) it,moltype(it),findtype(it),moltethered(it)
101   format(10x,i4,10x,a3,5x,i7,1x,l1)
            end do
            write(6,102)
102   format(/'remove_molecule_types> Give type number of molecule(s) to remove (q to quit).'
     */' and number to remove'//' If no number is given, all will be removed')
            read(5,'(a)',iostat=ierr) sometext
            if(ierr.eq.0.and.len_trim(sometext).gt.0) then
               call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
               readagain=sometext(ncf(1):ncl(1)).ne.'q'
               if(readagain) then
                  read(sometext(ncf(1):ncl(1)),*) mtyperead
                  ntoremove=-1
                  if(nwords.gt.1) then
                     read(sometext(ncf(2):ncl(2)),*,iostat=ierr) ntoremove
                  end if
                  if(ntoremove.lt.0) ntoremove=findtype(mtyperead)
                  moltyperead=moltype(mtyperead)
                  im=1
                  nremove=0
                  do while (im.le.nmol)
                     if(alabel(naf(im)).eq.moltyperead.and.nremove.lt.ntoremove) then
                        write(6,*) im,alabel(naf(im)),moltyperead
                        call remove_molecule(im)
                        nremove=nremove+1
                        findtype(mtyperead)=findtype(mtyperead)-1
                     else
                        im=im+1
                     end if
                  end do
               end if
!            else
!               readagain=.false.
            endif
         end do
c
c now check that all symbol_values have been redefined correctly for any type changes
c that have been made
c
         call write_symbols_changeato(nsymboluptodata,nsymbol) 
         return
      end
      
      subroutine remove_molecule(iremove)

		use molecom_mod
      !Removes the specified molecule
      
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      
      integer iremove,im,ia,ig,ib,ic,ifound,sumrotgrp
      
      !Obviously can't remove a non-existent molecule. Also we will always leave at least one molecule
      !in file, otherwise the .ato file becomes unreadable.
      if(iremove.gt.0.and.iremove.le.nmol.and.nmol.gt.1) then
         !Get the rotational group number up to this point
         ig=0
         im=1
         do while (im.lt.iremove)
            ig=ig+ngrot(im)
            im=im+1
         end do
         natom=natom-nato(im)
         ia=naf(im)-1
         do while(im.lt.nmol)
            naf(im)=ia+1
            nato(im)=nato(im+1)
            xadd(im)=xadd(im+1)
            yadd(im)=yadd(im+1)
            zadd(im)=zadd(im+1)
            phix(im)=phix(im+1)
            phiy(im)=phiy(im+1)
            phiz(im)=phiz(im+1)
            tethered(im)=tethered(im+1)
            ic=0
            iaref=naf(im+1)-1
            do while (ic.lt.nato(im))
               ic=ic+1
               ia=ia+1
               iaref=iaref+1
               alabel(ia)=alabel(iaref)
               xa(ia)=xa(iaref)
               ya(ia)=ya(iaref)
               za(ia)=za(iaref)
               nbond(ia)=nbond(iaref)
               do ib=1,nbond(ia)
                  iabond(ia,ib)=iabond(iaref,ib)
                  tbond(ia,ib)=tbond(iaref,ib)
               end do
            end do
            nal(im)=ia
            !Now transfer the rotational groups
            igref=ig+ngrot(im)
            ngrot(im)=ngrot(im+1)
            ib=0
            do while(ib.lt.ngrot(im))
               ib=ib+1
               ig=ig+1
               igref=igref+1
               rotgrp(ig)=rotgrp(igref)
               labaxi(ig,1)=labaxi(igref,1)
               labaxi(ig,2)=labaxi(igref,2)
               nrot(ig)=nrot(igref)
               rotf(ig)=rotf(igref)
               rotl(ig)=rotl(igref)
            end do
            im=im+1
         end do
         nmol=nmol-1
      end if
      !Finally clean up the type lists
      call cleanup_types()
      return
      end

      subroutine cleanup_types()
      !Renews atom and molecule type lists, removing anything that is no longer needed

      use molecom_mod
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      
      integer findtype(mcomp),it,im,ia,ifound,imoltype
c
c Check if any of the previous types no longer exist - if so remove them from the list
c
      do it=1,ntype
         findtype(it)=0
      end do
      do im=1,nmol
         do ia=naf(im),nal(im)
            ifound=0
            it=0
            do while(it.lt.ntype.and.ifound.eq.0)
               it=it+1
               if(alabel(ia).eq.type(it)) ifound=it
            end do
            if(ifound.gt.0) findtype(ifound)=findtype(ifound)+1
         end do
      end do
c
c remove any types that are no longer needed. Start from the end to make sure we catch them all
c
      it=ntype
      do while (it.gt.0)
         if(findtype(it).le.0) then
c
c remove this type
c
            jt=it
            do while (jt.lt.ntype)
               type(jt)=type(jt+1)
               epscor(jt)=epscor(jt+1)
               sigcor(jt)=sigcor(jt+1)
               amassa(jt)=amassa(jt+1)
               charge(jt)=charge(jt+1)
               qradius(jt)=qradius(jt+1)
               atsymbol(jt)=atsymbol(jt+1)
               iexchange(jt)=iexchange(jt+1)
               rhosphere(jt)=rhosphere(jt+1)
               natomsphere(jt)=natomsphere(jt+1)
               radminsphere(jt)=radminsphere(jt+1)
               radmaxsphere(jt)=radmaxsphere(jt+1)
               jt=jt+1
            end do
            ntype=ntype-1
         endif
         it=it-1
      end do
!Finally reset the type numbers for each atom
      do ia=1,natom
         it=0
         ifound=0
         do while(it.lt.ntype.and.ifound.eq.0)
            it=it+1
            if(alabel(ia).eq.type(it)) ifound=it
         end do
         itype(ia)=ifound
      end do
      
!Repeat this for the list of molecule types

      do it=1,nmoltype
         findtype(it)=0
      end do
      do im=1,nmol
         imoltype=0
         ifound=0
         do while(imoltype.lt.nmoltype.and.ifound.eq.0)
            imoltype=imoltype+1
            if(alabel(naf(im)).eq.moltype(imoltype)) then
               ifound=imoltype
            end if
         end do
         if(ifound.gt.0) then
            findtype(imoltype)=findtype(imoltype)+1
         endif
      end do
      do imoltype=1,nmoltype
         write(6,'(i4,1x,a3,1x,i7)') imoltype,moltype(imoltype),findtype(imoltype)
      end do
      !remove any types no longer needed. Start from end to make sure we catch all of them.
      imoltype=nmoltype
      do while(imoltype.gt.0)
         if(findtype(imoltype).le.0) then
            jt=imoltype
            do while (jt.lt.nmoltype)
               moltype(jt)=moltype(jt+1)
               moltethered(jt)=moltethered(jt+1)
               moltemplate(jt)=moltemplate(jt+1)
               wholerot(jt)=wholerot(jt+1)
               wholetrans(jt)=wholetrans(jt+1)
               jt=jt+1
            end do
            nmoltype=nmoltype-1
         end if
         imoltype=imoltype-1
      end do
      !Finally reset the type numbers for each molecule
      do im=1,nmol
         imoltype=0
         ifound=0
         do while(imoltype.lt.nmoltype.and.ifound.eq.0)
            imoltype=imoltype+1
            if(alabel(naf(im)).eq.moltype(imoltype)) then
               ifound=imoltype
            end if
         end do
         mtype(im)=ifound
      end do
      
      return
      end

      subroutine add_whole_rot(tokens,ntokens,ntf,ntl)

      use reallocation_routines
      use molecom_mod

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      character(len=3) :: ctemp
      integer rotatom1,rotatom2,rotatom1this,rotatom2this
      character(len=*), intent(in) :: tokens
      integer, intent(in) :: ntokens
      integer, dimension(*), intent(in) :: ntf,ntl

      !This routine will add a rotational group which will rotate the whole molecule about a specified pair of atoms
      if(ntokens.gt.1) then
         sometext=tokens(ntf(2):ntl(2))
      else
         write(6,101,advance='no')
101   format(/'changeato> Type atom label of first atom in molecule (to identify molecule),'
     */'<Enter> to leave unchanged: ')
         sometext=' '
         read(5,'(a)',iostat=ierr) sometext
      end if
      call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
      if(nwords.gt.0) then
         read(sometext(ncf(1):ncl(1)),*) ctemp
         !Find the molecule type this corresponds to
         imoltype=0
         ifound=0
         do while(imoltype.lt.nmoltype.and.ifound.eq.0)
            imoltype=imoltype+1
            if(ctemp.eq.moltype(imoltype)) then
               ifound=imoltype
            end if
         end do
         if(ifound.gt.0) then
            rotatom1=0
            rotatom2=0
            if(ntokens.gt.3) then
               read(tokens(ntf(3):ntl(4)),*,iostat=ierr) rotatom1,rotatom2
            else
               write(6,102,advance='no')
102   format(/'changeato> Type two atom numbers (within molecule) to form rotational group: ')
               read(5,*,iostat=ierr) rotatom1,rotatom2
            end if
            !Only proceed if both atom numbers are .gt. 0 and not the same value
            if(rotatom1.gt.0.and.rotatom2.gt.0.and.rotatom1.ne.rotatom2) then
               !Step through each molecule and stop if the molecule is the same type as that specified
               ig=0
               do im=1,nmol
                  ig=ig+ngrot(im)
                  if(mtype(im).eq.ifound.and.rotatom1.le.nato(im).and.rotatom2.le.nato(im)) then
                     !Set the rotational axis atom to the actual atom number for this molecule
                     rotatom1this=rotatom1+naf(im)-1
                     rotatom2this=rotatom2+naf(im)-1
                     !Create a new rotational group. First need to shift all the existing groups beyond ig up by 1
                     ig=ig+1
                     ih=ngrrot
                     ih1=ngrrot+1
                     ngrrot=ih1
                     call testsize_rotgrps()
                     do while(ih.ge.ig)
                        rotgrp(ih1)=rotgrp(ih)
                        labaxi(ih1,1)=labaxi(ih,1)
                        labaxi(ih1,2)=labaxi(ih,2)
                        nrot(ih1)=nrot(ih)
                        rotf(ih1)=rotf(ih)
                        rotl(ih1)=rotl(ih)
! Do any dihedral groups.
                        if(rotgrp(ih).eq.'DIH') then
                           ndihedralang(ih1)=ndihedralang(ih)
                           do k=1,ndihedralang(ih1)
                              dihedralang(k,ih1)=dihedralang(k,ih)
                              do j=1,2
                                 labrotdih(j,k,ih1)=labrotdih(j,k,ih)
                              end do
                           end do
                        end if
                        ih1=ih
                        ih=ih-1
                     end do
                     !The rotational group will be all the atoms in the molecule except rotatom1 and rotatom2
                     ngrot(im)=ngrot(im)+1
!                     call testsize_rotgrps()
                     rotgrp(ig)='ROT'
                     labaxi(ig,1)=rotatom1this
                     labaxi(ig,2)=rotatom2this
                     rotf(ig)=nrottot+1
! Number of atoms in rotational group is equal to number of atoms in molecule except the two forming the axis
                     nrot(ig)=nato(im)-2
                     nrottot=nrottot+nrot(ig)
                     rotl(ig)=nrottot
                     call testsize_labrot()
                     labref=rotf(ig)-1
                     write(6,*) ig,rotf(ig),rotl(ig),rotatom1this,rotatom2this,naf(im),nal(im)
                     do ia=naf(im),nal(im)
                        if(ia.ne.rotatom1this.and.ia.ne.rotatom2this) then
                           labref=labref+1
                           labrot(labref)=ia
                           write(6,*) labref,ia,labrot(labref)
                        end if
                     end do
                  end if
               end do
            end if
         end if
      end if
      return
      end

      subroutine add_image_atom(ido,tokens,ntokens,ntf,ntl)

      use reallocation_routines
      use molecom_mod

      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      character(len=3) :: ctemp
      character(len=1) :: yorno
      integer rotatom1,rotatom2,rotatom1s,rotatom2s,ifound
      logical ido,found
      character(len=*), intent(in) :: tokens
      integer, intent(in) :: ntokens
      integer, dimension(*), intent(in) :: ntf,ntl

      !This routine will add an image atom, joining one end of a long molecule to its tail on the opposite side of the simulation box.
      if(ntokens.gt.1) then
         sometext=tokens(ntf(2):ntl(2))
      else
         write(6,101,advance='no')
101   format(/'changeato> Type atom label of first atom in molecule (to identify molecule),'
     */'<Enter> to leave unchanged: ')
         sometext=' '
         read(5,'(a)',iostat=ierr) sometext
      end if
      call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
      if(nwords.gt.0) then
         read(sometext(ncf(1):ncl(1)),*) ctemp
         !Find the molecule type this corresponds to
         imoltype=0
         ifound=0
         do while(imoltype.lt.nmoltype.and.ifound.eq.0)
            imoltype=imoltype+1
            if(ctemp.eq.moltype(imoltype)) then
               ifound=imoltype
            end if
         end do
         if(ifound.gt.0.and.ido) then
            rotatom1s=-1
            rotatom2s=-1
            do while(rotatom1s.ne.0.and.rotatom2s.ne.0)
               if(ntokens.gt.3) then
                  sometext=tokens(ntf(3):ntl(4))
               else
                  write(6,102,advance='no')
102   format(/'changeato> Type atom number and image atom number to set (within same molecule) (<Enter> to quit): ')
                  read(5,'(a)',iostat=ierr) sometext
               end if
               call parse(sometext,nwords,ncf,ncl,11,nstsize,nworddim)
               if(nwords.gt.1) then
                  read(sometext,*) rotatom1s,rotatom2s
               else
                  rotatom1s=0
                  rotatom2s=0
               end if
               !Only proceed if both atom numbers are .gt. 0 and not the same value
               if(rotatom1s.gt.0.and.rotatom2s.gt.0.and.rotatom1s.ne.rotatom2s) then
                  !Step through each molecule and stop if the molecule is the same type as that specified
                  ig=0
                  do im=1,nmol
                     if(mtype(im).eq.ifound.and.rotatom1s.le.nato(im).and.rotatom2s.le.nato(im)) then
                        !Set the rotational axis atom to the actual atom number for this molecule
                        rotatom1=rotatom1s+naf(im)-1
                        rotatom2=rotatom2s+naf(im)-1
      write(6,*) im,mtype(im),ifound,rotatom1s,rotatom2s,rotatom1,rotatom2,naf(im)
                        imageatom(rotatom1)=rotatom2
                        imageatom(rotatom2)=rotatom1
                     end if
                  end do
               end if
            end do
         else if (ifound.gt.0.and..not.ido) then
!Make a list of the image atoms for this molecular type
            im=0
            found=.false.
            do while(im.lt.nmol.and..not.found)
               im=im+1
               found=mtype(im).eq.ifound
            end do
            if(found) then
               nadd=1-naf(im)
               do ia=naf(im),nal(im)
                  if(imageatom(ia).gt.0) then
                     iaremove=ia+nadd
                     iaremoveimage=imageatom(ia)+nadd
!                     write(6,103,advance='no') iaremove,iaremoveimage
103   format('changeato> Remove connection between atom ',i5,' and its image ',i5,' (y or n)? ')
!                     read(5,'(a)',iostat=ierr) yorno
                     yorno='y'
                     if(yorno.eq.'y'.or.yorno.eq.'Y') then
                        !Do this for all the other molecules of this type. im should be the first of its type
                        !found in the .ato file
                        im1=im-1
                        do while (im1.lt.nmol)
                           im1=im1+1
                           if(mtype(im).eq.ifound) then
                              nadd1=naf(im1)-1
                              imageatom(iaremove+nadd1)=0
                              imageatom(iaremoveimage+nadd1)=0
                           end if
                        end do
                     end if
                  end if
               end do
            end if
         end if
      end if
      return
      end

      subroutine make_mol_files()
      use reallocation_routines
c
c resets the value of a specified parameter in the .ato file if requested 
c by calling program
c
      include 'epsr_dimension.inc'
      include 'epsr16.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      include 'chemical_symbols.inc'
      character(len=3)                                :: int_as_str1,int_as_str2,int_as_str3
! Variables for checking
      character(len=256)                              :: moltemplatecheck,moltemplateshort
      integer                                         :: lenmoltemplatecheck,lenmoltemplateshort,nbondsfound,ndihedralsfound
      real, dimension(3)                              :: a,a1,a2,b,c1perp,c2perp,csign
      real                                            :: pi,piconv,angle
      logical                                         :: found

      pi=4.0*atan(1.0)
      piconv=180/pi
      nbondsfound=0
      ndihedralsfound=0
      do imoltype=1,nmoltype
         lenmoltemplateshort=indexlastdot(moltemplate(imoltype))-1
         if(lenmoltemplateshort.le.0) lenmoltemplateshort=len_trim(moltemplate(imoltype))
         moltemplateshort=moltemplate(imoltype)(1:lenmoltemplateshort)//'.mol'
         lenmoltemplateshort=lenmoltemplateshort+4
         call get_full_filename(moltemplateshort,lenmoltemplateshort
     *,moltemplatecheck,lenmoltemplatecheck)
         write(6,'(/1x,a)') moltemplatecheck(1:lenmoltemplatecheck)
! Try to open this file. If it exists, we won't do anything.
         open (10,file=moltemplatecheck,status='old',iostat=ierr)
         if(ierr.ne.0) then
! Create a new template file based on .ato filename plus label of first atom in molecule
            moltemplate(imoltype)=fname_inp(1:indexlastdot(fname_inp)-1)
     *//moltype(imoltype)(1:len_trim(moltype(imoltype)))
            lenmoltemplateshort=indexlastdot(fname_inp)-1+len_trim(moltype(imoltype))
            write(6,'(/1x,a,a)') 'make_mol_files> molecule template name: ',moltemplate(imoltype)(1:lenmoltemplateshort)
            moltemplatecheck=fnameinp(1:indexlastdot(fnameinp)-1)
     *//moltype(imoltype)(1:len_trim(moltype(imoltype)))//'.mol'
            lenmoltemplatecheck=indexlastdot(fnameinp)-1+len_trim(moltype(imoltype))+4
            write(6,'(/1x,"make_mol_files> creating .mol file: ",a)') moltemplatecheck(1:lenmoltemplatecheck)
            im=0
            found=.false.
            do while (im.lt.nmol.and..not.found)
               im=im+1
               found=alabel(naf(im)).eq.moltype(imoltype)
            end do
            if(found) then !Use this molecule to form the template
! The next bit is a replica of ghemical_routines_busch.f (readmole) with appropriate changes related to atom labelling.
               open (10,file=moltemplatecheck,status='unknown',iostat=ierr)
! Signal that this is a .gmol file and that the atom numbers listed are the EPSR atom numbers
               if(ierr.eq.0) then
                  write(10,*) '.gmol ',0
                  iasub=1-naf(im)
! Write atoms and their bonding arrangements.
                  do ia1=naf(im),nal(im)
                     ia=ia1+iasub
                     write(int_as_str1,'(i3.3)') ia
                     write(10,*) ' atom ',ia,' '//int_as_str1//' '
     *,xa(ia1),ya(ia1),za(ia1),nbond(ia1),((iabond(ia1,j)+iasub),j=1,nbond(ia1))
                  end do
! Write the bond list. 
                  do ia1=naf(im),nal(im)
                     write(int_as_str1,'(i3.3)') ia1+iasub
                     ib1=0
                     do while (ib1.lt.nbond(ia1))
                        ib1=ib1+1
                        ia2=iabond(ia1,ib1)
                        if(ia2.gt.ia1) then
                           nbondsfound=nbondsfound+1
                           write(int_as_str2,'(i3.3)') ia2+iasub
                           write(10,'(" bond ",a3,1x,a3,1x,f10.5)') int_as_str1,int_as_str2,tbond(ia1,ib1)
                        end if
                     end do
                  end do
! Write the bond angle list
                  do ia1=naf(im),nal(im)
                     write(int_as_str1,'(i3.3)') ia1+iasub
                     do ib1=1,nbond(ia1)
                        ia2=iabond(ia1,ib1)
                        write(int_as_str2,'(i3.3)') ia2+iasub
                        a(1)=xa(ia2)-xa(ia1)
                        a(2)=ya(ia2)-ya(ia1)
                        a(3)=za(ia2)-za(ia1)
                        ia2ref=ia2+iasub
                        do ib2=1,nbond(ia1)
                           ia3=iabond(ia1,ib2)
                           if(ia3.gt.ia2) then
                              write(int_as_str3,'(i3.3)') ia3+iasub
                              b(1)=xa(ia3)-xa(ia1)
                              b(2)=ya(ia3)-ya(ia1)
                              b(3)=za(ia3)-za(ia1)
                              angle=piconv*acos(scapro(a,b)/sqrt(scapro(a,a)*scapro(b,b)))
                              write(10,'(a,a,1x,f10.5)') ' angle '
     *,int_as_str2//' '//int_as_str1//' '//int_as_str3,angle
                           end if
                        end do
                     end do
                  end do
! Set up all possible dihedral angles
                  do ia1=naf(im),nal(im)
                     ib1=0
                     do while (ib1.lt.nbond(ia1))
                        ib1=ib1+1
                        ia2=iabond(ia1,ib1)
! Define 2 axis atoms
                        if(ia2.gt.ia1) then
                           b(1)=xa(ia2)-xa(ia1)
                           b(2)=ya(ia2)-ya(ia1)
                           b(3)=za(ia2)-za(ia1)
! Get neighbours of ia1
                           ib2=0
                           do while (ib2.lt.nbond(ia1))
                              ib2=ib2+1
                              ia1ref=iabond(ia1,ib2)
                              if(ia1ref.ne.ia2) then
                                 a1(1)=xa(ia1ref)-xa(ia1)
                                 a1(2)=ya(ia1ref)-ya(ia1)
                                 a1(3)=za(ia1ref)-za(ia1)
! Calculate the vector perpendicular to the plane formed by these two vectors
                                 call vecpro(a1,b,c1perp)
! Normalise this vector
                                 amp=sqrt(scapro(c1perp,c1perp))
                                 do i=1,3
                                    c1perp(i)=c1perp(i)/amp
                                 end do
! Now search the bonds of ia2 for the 4th atom
                                 ib3=0
                                 do while (ib3.lt.nbond(ia2))
                                    ib3=ib3+1
                                    ia2ref=iabond(ia2,ib3)
                                    if(ia2ref.ne.ia1.and.ia2ref.ne.ia1ref) then
                                       a2(1)=xa(ia2ref)-xa(ia2)
                                       a2(2)=ya(ia2ref)-ya(ia2)
                                       a2(3)=za(ia2ref)-za(ia2)
! Calculate the vector perpendicular to the plane formed by these two vectors
                                       call vecpro(a2,b,c2perp)
! Normalise this vector
                                       amp=sqrt(scapro(c2perp,c2perp))
                                       do i=1,3
                                          c2perp(i)=c2perp(i)/amp
                                       end do
! Cosine of the angle between them
                                       cosphi=scapro(c1perp,c2perp)
                                       if(cosphi.gt.1.0)cosphi=1.0
                                       if(cosphi.lt.-1.0) cosphi=-1.0
! Sine of the angle between them - see torangles for this
                                       call vecpro(c1perp,c2perp,csign)
                                       sinphi=scapro(csign,b)/sqrt(scapro(b,b))
                                       angle=piconv*atan2(sinphi,cosphi)
                                       write(10,*) 'dihedral ',ia1ref,ia1,ia2,ia2ref,angle
                                       ndihedralsfound=ndihedralsfound+1
                                    end if
                                 end do
                              end if
                           end do
                        endif
                     end do
                  end do
! Write the changelabel commands
                  do ia1=naf(im),nal(im)
                     write(int_as_str1,'(i3.3)') ia1+iasub
                     write(10,*) 'changelabel '//int_as_str1//' ',alabel(ia1)
                  end do
! Write potential parameters
                  do it=1,ntype
! Check that this type exists in the molecule
                     found=.false.
                     ia=naf(im)-1
                     do while (ia.lt.nal(im).and..not.found)
                        ia=ia+1
                        found=alabel(ia).eq.type(it)
                     end do
                     if(found) write(10,600) type(it)
     *,epscor(it),sigcor(it),amassa(it),charge(it),atsymbol(it)
600   format(1x,'potential ',a3,1x,4(f10.5,1x),a3)
                  end do
                  write(10,*) 'temperature ',300
                  write(10,*) 'vibtemp ', 65.0
                  write(10,*) 'angtemp ', 1.0
                  write(10,*) 'dihtemp ', 10.0
                  write(10,*) 'ecoredcore ', 0.0, 0.0
                  write(10,*) 'density ', rho
! Done. 
                  close(10)
                  write(6,'(/1x,a)') 'make_mol_files> *** Remember to run fmole on the new .ato file. ***'
               end if
            end if
         end if
      end do
      if(nbondsfound.gt.0) then
         stepmi=0.1
      else
         stepmi=0.0
      end if
      if(ndihedralsfound.gt.0) then
         stepri=0.1
      else
         stepri=0.0
      end if
! Ensure the bond angle and dihedral weights have been set
      tempangreq=1.0
      tempdihreq=10.0
      tempangwt=1.0
      tempdihwt=10.0
      return
      end

      subroutine set_molecule_scale_factors(tokens,ntokens,ntf,ntl)
      !Sets the molecular scale factor for a particular molecule type
      
      use pfroutines
      use molecom_mod

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      integer it,mtyperead,ierr,nset,ntoset,findtype(mcomp)
      character(len=3) :: moltyperead
      real rfract,thrd !Relative width of Poisson function. 0.0 means all molecules have the same size, scale factor = 1.0
      character(len=*), intent(in) :: tokens
      integer, intent(in) :: ntokens
      integer, dimension(*), intent(in) :: ntf,ntl

      thrd=1.0/3.0
!Determine the number of molecules of each type
      do it=1,nmoltype
         findtype(it)=0
      end do
      do im=1,nmol
         imoltype=0
         ifound=0
         do while(imoltype.lt.nmoltype.and.ifound.eq.0)
            imoltype=imoltype+1
            if(alabel(naf(im)).eq.moltype(imoltype)) then
               ifound=imoltype
            end if
         end do
         if(ifound.gt.0) then
            findtype(imoltype)=findtype(imoltype)+1
         endif
      end do

!Make a list of all current molecule types and their current scale factor
      write(6,100)
100   format(/1x,'molecule type number',2x,'label',2x,'number',2x,'relative width'/)
      do it=1,nmoltype
         write(6,101) it,moltype(it),findtype(it),rfractmol(it)
101   format(10x,i4,10x,a3,5x,i0,1x,f10.5)
      end do
      if(ntokens.gt.2) then
         sometext=tokens(ntf(2):ntl(3))
         ierr=0
      else
         write(6,102)
102   format(/'set_molecule_scale_factors> Give type number of molecule(s) to set (q to quit).'
     */' and relative width of Poisson function (0 - 1)')
         read(5,'(a)',iostat=ierr) sometext
      end if
      if(ierr.eq.0.and.len_trim(sometext).gt.0) then
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         read(sometext(ncf(1):ncl(1)),*,iostat=ierr) mtyperead
         if(nwords.gt.1.and.ierr.eq.0) then
            read(sometext(ncf(2):ncl(2)),*,iostat=ierr) rfract
            if(rfract.gt.1.0) rfract=1.0
            if(rfract.gt.0.0) call getpfprobability(rfract)
            rfractmol(mtyperead)=rfract
            ntoset=findtype(mtyperead)
            moltyperead=moltype(mtyperead)
            im=0
            nset=0
            do while (im.lt.nmol)
               im=im+1
               if(alabel(naf(im)).eq.moltyperead.and.nset.lt.ntoset) then
                  if(rfract.gt.0.0) then
                     scalefacmol(im)=(getrratio())**thrd !Assume the MASS of the particle is distributed normally
                  else
                     scalefacmol(im)=1.0
                  end if
                  write(6,*) im,alabel(naf(im)),' ',moltyperead,scalefacmol(im)
                  nset=nset+1
               end if
            end do
         end if
      endif
      return
      end
      

