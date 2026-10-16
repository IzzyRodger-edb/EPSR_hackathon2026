      subroutine initran1(nchan)
      include 'ran1com.inc'
      integer nchan,ierr
c
c open file containing random number parameters, or else initialise
c
            read(nchan,*,iostat=ierr) mcall,idum,iy,(iv(i),i=1,ntab)
c
c if the number of calls to ran1 exceeds 1.0E8, the sequence is reset, as this
c is approximately 1/20th of the period of the ran1 generator
c
c           if(mcall.gt.100000000) then
c                 idum=0
c                 iy=0
c                 mcall=0
c           endif

c if error or end of file occurs, initialise values

            if(ierr.ne.0) then
 
               idum=0
               iy=0
               mcall=0

            endif

            return

      end

      subroutine saveran1(nchan)
      include 'ran1com.inc'
      integer nchan
c
c open file containing random number parameters, or else initialise
c
      write(nchan,*) mcall,idum,iy,(iv(i),i=1,ntab)
      return
      end

      FUNCTION ran1()
      include 'ran1com.inc'
c
c \Minimal" random number generator of Park and Miller with Bays-Durham shue and
c added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
c the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
c alter idum between successive deviates in a sequence. RNMX should approximate the largest
c floating value that is less than 1.
c
c
c Initialize if idum.le.0
c
      if (idum.le.0.or.iy.eq.0) then 
c
c Be sure to prevent idum = 0.
c
            idum=max(-idum,1) 
c
c Load the shue table (after 8 warm-ups).
c
            do j=NTAB+8,1,-1 
                  k=idum/IQ
                  idum=IA*(idum-k*IQ)-IR*k
                  if (idum.lt.0) idum=idum+IM
                  if (j.le.NTAB) iv(j)=idum
            enddo
            iy=iv(1)
      endif
c
c Start here when not initializing.
c
      k=idum/IQ 
c
c Compute idum=mod(IA*idum,IM) without overflows by Schrage's method.
c
      idum=IA*(idum-k*IQ)-IR*k 
      if (idum.lt.0) idum=idum+IM
c
c Will be in the range 1:NTAB.
c
      j=1+iy/NDIV 
      iy=iv(j) 
c
c Output previously stored value and rell the shue table.
c
      iv(j)=idum
c
c Because users don't expect endpoint values.
c
      ran1=min(AM*iy,RNMX) 
      mcall=mcall+1
      return
      END
