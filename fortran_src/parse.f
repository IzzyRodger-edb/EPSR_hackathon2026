      subroutine parse(sometext,nwords,ncf,ncl,nspcword
     *,nchartext,nworddim)
c
c parses the text in sometext to separate words with the assumed delimiter a space
c ncf and ncl contain the first and last character numbers of each word
c ncharword is the maximum number of characters in a word, nworddim is dimension of the
c word array, nspcword is the number of spaces between words needed to decide that the 
c end of text has been reached.
c
      character(len=*) sometext
      integer nwords,ncf(*),ncl(*),nspcword,nchartext,nworddim

!         write(6,*) sometext(1:len_trim(sometext))
c      write(6,'(a60)') sometext(1:60)
c
c step through the characters and identify words (separated by spaces). Extra spaces
c are ignored unless the number of spaces exceeds nspcword, in which case the search
c for words terminates
c
c First bypass any blank space before the first non-space character
c
      nspace=0
      nwords=0
      i=1
      do while(sometext(i:i).eq.' '.and.i.lt.nchartext)
         i=i+1
      end do
      nfirst=i
      nlast=i
      do while(i.lt.nchartext.and.nspace.le.nspcword)
         i=i+1
c Check whether the current letter is a space
         if(sometext(i:i).eq.' ') then
c A space is found
c Check whether this is the end of a word
            if(nlast.ge.nfirst) then
c New word found
               if(nwords.lt.nworddim) then
                  nwords=nwords+1
                  ncf(nwords)=nfirst
                  ncl(nwords)=nlast
               end if
               nspace=1
            else
               nspace=nspace+1
            endif
            nfirst=i+1
         else
            nlast=i
         endif
      end do
c Get the final word in case the text does not end with a space
c First check whether the last letter is a space
      if(sometext(i:i).ne.' ') then
         if(nlast.ge.nfirst) then
c New word found
            if(nwords.lt.nworddim) then
               nwords=nwords+1
               ncf(nwords)=nfirst
               ncl(nwords)=nlast
            endif
         endif
      endif
c      write(6,'(1x,i5,/((2(i5,1x),a)))') nwords,(ncf(j),ncl(j),sometext(ncf(j):ncl(j)),j=1,nwords)

      return
      end

      subroutine concatenate(nfirst,nlast,return_string,ncharstring
     *,nwords,ncf,ncl,input_string)
c
c concatenates a series of words into a single text string separated by single spaces. Maximum length of string
c is 256 characters. The string begins with nfirst) and ends with nlast.
c
      character(len=*) input_string,return_string
      character(len=256) output_string
      integer nwords,ncf(*),ncl(*),ncharword,ncharstring
c
c maximum number of characters in a word
c
      ncharword=256
      ncharstring=0
      if(nlast.gt.nwords) nlast=nwords
      if(nfirst.lt.1) nfirst=1
c
c set a flag to say whether the character limit has been reached
c
      iflag=0
      iw=nfirst
      do while(iw.le.nlast.and.ncharstring.lt.ncharword)
         ncharadd=ncl(iw)-ncf(iw)+1
         nchartest=ncharstring+ncharadd
c Need to allow for the space between words if there is more than one word in the string
         if(iw.gt.nfirst) nchartest=nchartest+1
         if(nchartest.gt.ncharword) then
            ncharadd=ncharword-ncharstring
c Again need to allow for space between words if there is more than one word already
            if(iw.gt.nfirst) ncharadd=ncharadd-1
            iflag=1
         endif
         if(iw.eq.nfirst) then
            output_string=input_string(ncf(iw):min(ncl(iw),ncf(iw)+ncharadd-1))
            ncharstring=ncharadd
         else
c
c test the final number of characters in the text
c
            if(ncharadd.gt.0) output_string=output_string(1:ncharstring)
     *//' '//input_string(ncf(iw):min(ncl(iw),ncf(iw)+ncharadd-1))
            ncharstring=ncharstring+ncharadd+1
         endif
         iw=iw+1
      end do
! Make we don't exceed the length of return_string
      if(ncharstring.gt.len(return_string)) ncharstring=len(return_string)
      return_string=output_string(1:ncharstring)
      return
      end


