      INTEGER IA,IM,IQ,IR,NTAB,NDIV
      REAL AM,EPS,RNMX  
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     * NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum,iv(NTAB),iy,mcall
      common /ran1com/ mcall,idum,iy,iv
