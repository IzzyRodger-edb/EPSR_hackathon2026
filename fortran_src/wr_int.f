      subroutine r_int(lpt,fname,nfile,xd,yd,ed,nxdim,nfdim,ierr) !read data
      include 'wrcom.inc'
      real*4      xd(nxdim),yd(nxdim,nfdim),ed(nxdim,nfdim)
      character*256 fname,sometext
      integer*4 ierr
      real, parameter :: limit=1.0e-20
      nt=(nfile+(nwrmax-1))/nwrmax
      ilen=indexlastdot(fname)+1
      icount=0
      do it=1,nt
c
c number of elements to read for each x value
c
            iwr=nfile-(it-1)*nwrmax
            if(iwr.gt.nwrmax) iwr=nwrmax
            write(ext,331) it
331   format(i2.2)
            fname=fname(1:ilen)//ext(1:2)
            open(8,file=fname,status='old',iostat=ierr)
            if(ierr.ne.0) return
c
c check to see if the first line is a comment line
c
            read(8,100) sometext
100   format(a)
            ilenhash=index(sometext,'#')
            if(ilenhash.le.0) then
c
c hash not found, therefore assume the line is a line of numbers -
c start the file from the beginning
c
               rewind(8)
            endif
c
c READ it
c
            i=1
200         continue
                  read(8,*,end=201) xo(i),(yo(i,ir),eo(i,ir),ir=1,iwr)
! Remove absurdly small values
                  do ir=1,iwr
                     if(abs(yo(i,ir)).lt.limit) yo(i,ir)=0.0
                     if(abs(eo(i,ir)).lt.limit) eo(i,ir)=0.0
                  i=i+1
            go to 200
201         lpt=i-1
            write(6,101) lpt,fname(1:ilen+2)
101         format(1x,'Read ',i4,' points from ',a)
c
c save the data for this file
c
            do ir=1,iwr
                  icount=icount+1
                  do i=1,lpt
                        yd(i,icount)=yo(i,ir)
                        ed(i,icount)=eo(i,ir)
                  end do
            end do
            close(8)
      end do
c
c copy the x values
c
      do i=1,lpt
            xd(i)=xo(i)
      end do
      return
      end
      subroutine w_int(lpt,fname,nfile,xd,yd,ed
     *,titlex,titled,nxdim,nfdim)
      include 'wrcom.inc'
      real*4 xd(nxdim),yd(nxdim,nfdim),ed(nxdim,nfdim)
      character fname*256,formatout*30
      character formattitle*30
      character*14 titlex
      character*29 titled(nfdim)
      nt=(nfile+(nwrmax-1))/nwrmax
      ilen=indexlastdot(fname)+1
      icount=0
c
c copy the x values
c
      do i=1,lpt
            xo(i)=xd(i)
      end do
      do it=1,nt
c
c number of elements to write for each x value
c
            iwr=nfile-(it-1)*nwrmax
            if(iwr.gt.nwrmax) iwr=nwrmax
c
c setup output format specification
c
            nblock=2*iwr+1
            write(formattitle,333) '(1h#,1x,a14,',iwr,'(1x,a29))'
333   format(a12,i3.3,a9)
            write(formatout,332) '(',nblock,'(1x,e14.7))'
            lrec=nblock*15+1
332   format(a1,i3.3,a11)
            write(6,*) formatout
            write(ext,331) it
331   format(i2.2)
            fname=fname(1:ilen)//ext(1:2)
            write(6,101) lpt,fname(1:ilen+2)
101         format(1x,'Writing ',i4,' points to ',a)
            open(8,file=fname,status='unknown',form='formatted'
     1,recl=lrec)
c
c get the data for this file
c
            do ir=1,iwr
                  icount=icount+1
                  titleo(ir)=titled(icount)
                  do i=1,lpt
                        yo(i,ir)=yd(i,icount)
                        eo(i,ir)=ed(i,icount)
                  end do
            end do
c
c write it
c
            write(8,fmt=formattitle,iostat=ierr) titlex,(titleo(ir),ir=1,iwr)
            do i=1,lpt
                  write(8,fmt=formatout,iostat=ierr) xo(i),(yo(i,ir),eo(i,ir),ir=1,iwr)
            end do
            close(8)
      end do
100   format(13(1x,e14.7))
      return
      end
      subroutine w_intr(lpt1,lpt2,fname,nfile,xd,yd,ed
     *,titlex,titled,nxdim,nfdim)
      include 'wrcom.inc'
      real*4 xd(nxdim),yd(nxdim,nfdim),ed(nxdim,nfdim)
      character fname*256,formatout*30
      character formattitle*30
      character*14 titlex
      character*29 titled(nfdim)
      nt=(nfile+(nwrmax-1))/nwrmax
      ilen=indexlastdot(fname)+1
      icount=0
      lpt=lpt2-lpt1+1
c
c copy the x values
c
      do i=lpt1,lpt2
            xo(i)=xd(i)
      end do
      do it=1,nt
c
c number of elements to write for each x value
c
            iwr=nfile-(it-1)*nwrmax
            if(iwr.gt.nwrmax) iwr=nwrmax
c
c setup output format specification
c
            nblock=2*iwr+1
            write(formattitle,333) '(1h#,1x,a14,',nblock,'(1x,a29))'
333   format(a12,i3.3,a9)
            write(formatout,332) '(',nblock,'(1x,e14.7))'
            lrec=nblock*15+1
332   format(a1,i3.3,a11)
            write(6,*) formatout
            write(ext,331) it
331   format(i2.2)
            fname=fname(1:ilen)//ext(1:2)
            write(6,101) lpt,fname(1:ilen+2)
101         format(1x,'Writing ',i4,' points to ',a)
            open(8,file=fname,status='unknown',form='formatted'
     1,recl=lrec)
c
c get the data for this file
c
            do ir=1,iwr
                  icount=icount+1
                  titleo(ir)=titled(icount)
                  do i=lpt1,lpt2
                        yo(i,ir)=yd(i,icount)
                        eo(i,ir)=ed(i,icount)
                  end do
            end do
c
c write it
c
            write(8,fmt=formattitle) titlex,(titleo(ir),ir=1,iwr)
            do i=lpt1,lpt2
                  write(8,fmt=formatout) xo(i),(yo(i,ir),eo(i,ir),ir=1,iwr)
            end do
            close(8)
      end do
100   format(13(1x,e14.7))
      return
      end


