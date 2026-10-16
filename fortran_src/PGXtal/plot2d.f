CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C***<plot>**************************************************************
C
      SUBROUTINE PLOT(X,NX,Y,NY,Z,N1,N2,W,SIZE,IWIDTH,XLBL,YLBL,TITL)
C     ---------------------------------------------------------------
C
      REAL          X(*),Y(*),Z(N1,N2),W(*)
      CHARACTER*(*) XLBL,YLBL,TITL
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C Purpose
C      This subroutine plots "data" defined on a regularly-spaced 
C   rectangular grid of points Z(I,J). With the default choice for the 
C   PGCELL routine that is linked, the output is a linearly-interpolated 
C   map (rather than coarse rectangular boxes).
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    Z        R*4    I    N1 x N2   The recatangular "data"-array.
C    Z        R*4    O    N1 x N2   A scaled, and clipped, version of 
C                                   the input array(!).
C    N1       I*4    I       -      The first dimension of array Z.
C    N2       I*4    I       -      The second dimension of array Z.
C    X        R*4    I      NX      Array of X-coordinates.
C    NX       I*4    I       -      Number of X-pixels to be plotted 
C                                  (usually = N1, but must be <= N1).
C    Y        R*4    I      NY      Array of Y-coordinates.
C    NY       I*4    I       -      Number of Y-pixels to be plotted 
C                                  (usually = N2, but must be <= N2).
C    W        R*4    I      NX      Work array, at least NX long.
C    SIZE     R*4    I       -      Character-size for plot (try 1.5).
C    IWIDTH   I*4    I       -      Line-width for plot (try 2).
C    XLBL     A*1    I     *(*)     Label for X-axis.
C    YLBL     A*1    I     *(*)     Label for Y-axis.
C    TITL     A*1    I     *(*)     Title for plot.
C
C Globals
C    COLTABS.INC
C
C History
C   Initial release.                                    DSS:  3 Jul 1992
C   Minor changes to conform with new PGCELL.           DSS:  6 Feb 1995
C   Put in option to over-lay contours.                 DSS: 21 Feb 1995
C   No longer uses non-standard Q-format for defaults.  DSS: 27 Aug 1997
C   Fortran made LINUX-friendly!                        DSS: 15 Sep 1997
C-----------------------------------------------------------------------
C
      LOGICAL      OVRLAY
      CHARACTER*32 STRING
C
      WRITE(*,*)
      WRITE(*,*)'                (0) Contour'
      WRITE(*,*)'                (1) Surface'
      WRITE(*,*)'                (2) Colour: Grey-Scale'
      WRITE(*,*)'                (3) Colour: Heat'
      WRITE(*,*)'                (4) Colour: Rainbow Spectrum'
      WRITE(*,*)'                (5) Colour: BGYRW'
      WRITE(*,*)'                (6) Colour: Serpent'
      WRITE(*,*)'                (7) Colour: Read in from file'
      WRITE(*,*)
   1  IPLOT=0
      WRITE(*,100)
 100  FORMAT(' PLOT>  Type ?  : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=1) IPLOT
      IF (IPLOT.EQ.0) THEN
        CALL CONTOR(Z,NX,NY,X,Y,N1,N2,SIZE,IWIDTH,.FALSE.)
      ELSEIF (IPLOT.EQ.1) THEN
        CALL FMXMN(Z,N1*N2,YMAX1,YMIN)
        YINC=(YMAX1-YMIN)/FLOAT(NY-1)
        YMAX=2.0*YMAX1-YMIN
        XINC=(X(NX)-X(1))/FLOAT(NY-1)
        XMAX=2.0*X(NX)-X(1)
        CALL PGBEGIN(0,'?',1,1)
        CALL PGPAPER(0.0,1.0)
        CALL PGSCH(SIZE)
        CALL PGSLW(IWIDTH)
        CALL PGENV(X(1),XMAX,YMIN,YMAX,0,-2)
        CALL PGHI2D(Z,N1,N2,1,NX,1,NY,X,NX/NY,YINC,.TRUE.,W)
      ELSEIF (IPLOT.LE.7) THEN
        CALL GREY(Z,NX,NY,X,Y,N1,N2,IPLOT,SIZE,IWIDTH)
      ELSE
        GOTO 1
      ENDIF
      CALL PGLABEL(XLBL,YLBL,TITL)
      IF (IPLOT.GT.1) THEN
        WRITE(*,*)
        CALL LOGQYN(' PLOT> Over-lay contours ?','N',OVRLAY)
        IF (OVRLAY) CALL CONTOR(Z,NX,NY,X,Y,N1,N2,SIZE,IWIDTH,OVRLAY)
      ENDIF
      CALL PGEND
      END
C
      SUBROUTINE FMXMN(Y,N,YMAX,YMIN)
C     -------------------------------
C
      REAL         Y(*)
      CHARACTER*32 STRING
C
      SUM=0.0
      YMIN1=+1.0E25
      YMAX1=-1.0E25
      DO 10 I=1,N
        YI=Y(I)
        YMIN1=MIN(YMIN1,YI)
        YMAX1=MAX(YMAX1,YI)
        SUM=SUM+YI
  10  CONTINUE
      WRITE(*,*)
      WRITE(*,*)' Min, Max & Sum = ',YMIN1,YMAX1,SUM
      WRITE(*,*)
   1  YMIN=YMIN1
      YMAX=YMAX1
      WRITE(*,100) YMIN,YMAX
 100  FORMAT(' Surface>  Zmin & Zmax for plot ? (def=',1pe10.3,
     *       ',',E10.3,')  : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=1) YMIN,YMAX
      END
C
C***<contour>***********************************************************
C
      SUBROUTINE CONTOR(MAP,NX,NY,X,Y,N1,N2,SIZE,IWIDTH,OVRLAY)
C     ---------------------------------------------------------
C
      REAL    MAP(*),X(*),Y(*)
      REAL    CONT(25,2),TR(6)
      INTEGER NCONT(2),LCOLOR(2),LSTYLE(2),LWIDTH(2)
      LOGICAL OVRLAY
      DATA    LCOLOR,LSTYLE,LWIDTH,NWIND,IDP /1,3,1,1,1,2,1,0/
C
      CALL INIT(X,Y,NX,NY,TR,MAP,DMIN,DMAX,N1*N2)
      IF (OVRLAY) THEN
        CALL AUTCNT(NCONT,CONT,DMIN,DMAX,OVRLAY)
        CALL PGSLW(1)
        CALL PGCONT(MAP,N1,N2,1,NX,1,NY,CONT,NCONT,TR)
        CALL PGSLW(IWIDTH)
        RETURN
      ENDIF
      CALL ASKCNT(NCONT,CONT,DMIN,DMAX)
      CALL PGBEGIN(0,'?',1,1)
      CALL PGPAPER(0.0,1.0)
      CALL PGSCH(SIZE)
      CALL PGSLW(IWIDTH)
      CALL PGENV(X(1),X(NX),Y(1),Y(NY),0,0)
      DO 20 I=1,2
        CALL PGSCI(LCOLOR(I))
        CALL PGSLS(LSTYLE(I))
        LWDTH=2*LWIDTH(I)
        IF (LWDTH.GT.7) LWDTH=7
        CALL PGSLW(LWDTH)
        CALL PGCONT(MAP,N1,N2,1,NX,1,NY,CONT(1,I),NCONT(I),TR)
  20  CONTINUE
      CALL PGSCI(1)
      CALL PGSLW(IWIDTH)
      END
C
      SUBROUTINE INIT(X,Y,NX,NY,TR,MAP,DMIN,DMAX,NMAP)
C     ------------------------------------------------
C
      REAL X(*),Y(*),TR(*),MAP(*)
C
      TR(1)=X(1)-(X(NX)-X(1))/FLOAT(NX-1)
      TR(2)=(X(NX)-X(1))/FLOAT(NX-1)
      TR(3)=0.0
      TR(4)=Y(1)-(Y(NY)-Y(1))/FLOAT(NY-1)
      TR(5)=0.0
      TR(6)=(Y(NY)-Y(1))/FLOAT(NY-1)
      DMIN=+1.0E+20
      DMAX=-1.0E+20
      DTOT=0.0
      DO 10 I=1,NMAP
        F=MAP(I)
        DTOT=DTOT+F
        IF (F.GT.DMAX) DMAX=F
        IF (F.LT.DMIN) DMIN=F
  10  CONTINUE
      WRITE(*,*)
      WRITE(*,*)' Minimum value = ',DMIN
      WRITE(*,*)' Maximum value = ',DMAX
      WRITE(*,*)' Total flux    = ',DTOT
      WRITE(*,*)
      END
C
      SUBROUTINE ASKCNT(NCONT,CONT,DMIN,DMAX)
C     ---------------------------------------
C
      REAL           CONT(25,*),C(25)
      INTEGER        NCONT (*)
      LOGICAL        AUTO
      CHARACTER*132 STRING
C
      WRITE(*,*)
      CALL LOGQYN(' Contours>  Autoscale (linear) ?','Y',AUTO)
      IF (AUTO) THEN
        CALL AUTCNT(NCONT,CONT,DMIN,DMAX,.FALSE.)
        RETURN
      ENDIF
      WRITE(*,*)
      WRITE(*,*)'          ***** Contour Values  *****'
      WRITE(*,*)
   1  WRITE(*,100)
 100  FORMAT(' Thin>  ',$)
      READ(*,200,ERR=1) STRING
 200  FORMAT(A)
      CALL FINDNC(STRING,132,NCONT(1))
      READ(STRING,*,ERR=1) (CONT(I,1),I=1,NCONT(1))
      WRITE(*,*)
   2  WRITE(*,110)
 110  FORMAT(' Thick>  ',$)
      READ(*,200,ERR=2) STRING
      CALL FINDNC(STRING,132,NCONT(2))
      READ(STRING,*,ERR=2) (CONT(I,2),I=1,NCONT(2))
      END
C
      SUBROUTINE FINDNC(STRING,NCHARS,NC)
C     -----------------------------------
C
      CHARACTER STRING(*)
C
      DO 10 I=1,NCHARS
  10    IF (STRING(I).NE.' ') GOTO 1
   1  IMIN=I
      DO 20 I=NCHARS,1,-1
  20    IF (STRING(I).NE.' ') GOTO 2
   2  IMAX=I
      IF (IMIN.LE.IMAX) THEN
        NC=1
        J=IMIN
        DO 30 I=IMIN,IMAX
          IF (J.GT.IMAX) RETURN
          IF (STRING(J).EQ.' ' .OR. STRING(J).EQ.',') THEN
            NC=NC+1
   3        IF (STRING(J+1).EQ.' ' .OR. STRING(J+1).EQ.',') THEN
              J=J+1
              GOTO 3
            ENDIF
          ENDIF
          J=J+1
  30    CONTINUE
      ELSE
        NC=0
      ENDIF
      END
C
      SUBROUTINE AUTCNT(NCONT,CONT,DMIN,DMAX,OVRLAY)
C     ----------------------------------------------
C
      REAL         CONT(25,*)
      INTEGER      NCONT(*)
      LOGICAL       OVRLAY
      CHARACTER*32 STRING
C
      WRITE(*,*)
   1  WRITE(*,100) DMIN,DMAX
 100  FORMAT(' Autoscale> Range ? (def=',1pe10.3,' to ',e10.3,')  : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=1) DMIN1,DMAX1
      IF (NNN.NE.0) THEN
        DMIN=DMIN1
        DMAX=DMAX1
      ENDIF
   2  N=10
      WRITE(*,110) N
 110  FORMAT(' Autoscale>  No. of contours ?  (def=',I2,') : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=2) N
      N=MAX(MIN(N,50),2)
      NCONT(1)=(N+1)/2
      NCONT(2)=N-NCONT(1)
      IF (OVRLAY) THEN
        IF (N.GT.25) N=25
        NCONT(1)=N
        NCONT(2)=0
      ENDIF
      XINC=(DMAX-DMIN)/FLOAT(N)
      X=DMIN+0.5*XINC
      DO 10 I=1,NCONT(1)
        CONT(I,1)=X
        X=X+XINC
  10  CONTINUE
      DO 20 I=1,NCONT(2)
        CONT(I,2)=X
        X=X+XINC
  20  CONTINUE
      END
C
C***<grey scale>********************************************************
C
      SUBROUTINE GREY(MAP,NX,NY,X,Y,N1,N2,IPLOT,SIZE,IWIDTH)
C     ------------------------------------------------------
C
      REAL MAP(*),X(*),Y(*)
      REAL TR(6),TRCOL(6),COLBAR(2,256),LUTUSR(3,256)
      REAL RED(256),GR(256),BL(256)
      DATA TRCOL /0.0,1.0,0.0,0.0,0.0,1.0/
C
      NCOLS=256
      CALL INIT(X,Y,NX,NY,TR,MAP,DMIN,DMAX,N1*N2)
      CALL STGREY(MAP,N1*N2,DMIN,DMAX,COLBAR,NCOLS)
      IF (IPLOT.EQ.7) CALL LUTIN(LUTUSR,NCOLS,NCLUSR,IPLOT)
      CALL PGBEGIN(0,'?',1,1)
      CALL PGPAPER(0.0,1.0)
      CALL PGSCH(SIZE)
      CALL PGSLW(IWIDTH)
      CALL PGVPORT(0.89,0.93,0.2,0.90)
      CALL PGWINDOW(1.0,2.0,1.0,256.0)
      CALL SETCOL(IPLOT,NCOLS,LUTUSR,RED,GR,BL)
      CALL PGBOX('BC',0.0,0,'BC',0.0,0)
      CALL PGCELL(COLBAR,2,256,1,2,1,256,1.0,0.0,TRCOL,NCOLS,RED,GR,BL)
      CALL PGBOX('BC',0.0,0,'BC',0.0,0)
      CALL PGVPORT(0.12,0.82,0.2,0.90)
      XMIN=TR(1)+TR(2)
      XMAX=TR(1)+FLOAT(NX)*TR(2)
      YMIN=TR(4)+TR(6)
      YMAX=TR(4)+FLOAT(NY)*TR(6)
      CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
      CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
      CALL PGCELL(MAP,N1,N2,1,NX,1,NY,1.0,0.0,TR,NCOLS,RED,GR,BL)
      CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
      END
C
      SUBROUTINE STGREY(MAP,NMAP,FMIN,FMAX,COLBAR,NCOLS)
C     --------------------------------------------------
C
      REAL         MAP(*),COLBAR(2,*)
      CHARACTER*32 STRING
C
   1  WRITE(*,100) FMIN,FMAX
 100  FORMAT(' >>  Zmin & Zmax for plot ? (def=',1pe10.3,
     *       ',',E10.3,')  : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=1) FMIN2,FMAX2
      IF (NNN.NE.0) THEN
        FMIN=FMIN2
        FMAX=FMAX2
      ENDIF
      IF (ABS(FMAX-FMIN).LT.1.0E-20) GOTO 1
      FNORM=1.0/(FMAX-FMIN)
   2  C=1.0
      WRITE(*,110) C
 110  FORMAT(' >>  Contrast factor ?  (def=',F3.1,')  : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=2) C
      C=MAX(MIN(C,10.0),0.01)
      DO 10 I=1,NMAP
        F=(MAP(I)-FMIN)*FNORM
	IF (F.LE.0.0) THEN
          MAP(I)=0.0
        ELSEIF (F.GE.1.0) THEN
          MAP(I)=1.0
	ELSE
          MAP(I)=F**C
	ENDIF
  10  CONTINUE
      DCOL=0.999/FLOAT(NCOLS-1)
      COL=0.0
      IF (FMAX.LT.FMIN) THEN
        COL=1.0
        DCOL=-DCOL
      ENDIF
      DO 20 I=1,NCOLS
        COLBAR(1,I)=COL**C
        COLBAR(2,I)=COL**C
        COL=COL+DCOL
  20  CONTINUE
      END
C
      SUBROUTINE LUTIN(LUTUSR,NCOLS,NCLUSR,IPLOT)
C     -------------------------------------------
C
      CHARACTER*72 FILNAM
      REAL LUTUSR(3,NCOLS)
C
   1  WRITE(*,100)
 100  FORMAT(' INPUT> Filename for user colour-table ?  : ',$)
      READ(*,200,ERR=1) FILNAM
 200  FORMAT(A)
      OPEN(UNIT=17,FILE=FILNAM,STATUS='OLD',FORM='FORMATTED',ERR=1)
      DO 10 J=1,NCOLS
  10    READ(17,*,ERR=2,END=2) (LUTUSR(I,J),I=1,3)
   2  CLOSE(UNIT=17)
      NCLUSR=J-1
      WRITE(*,*)' No. of colour indicies read in = ',NCLUSR
      IF (NCLUSR.LE.1) THEN
        IPLOT=2
      ELSE
        NCOLS=NCLUSR
      ENDIF
      END
C
      SUBROUTINE SETCOL(IPLOT,NCOLS,LUTUSR,R,G,B)
C     -------------------------------------------
C
      INCLUDE 'COLTABS.INC'
      REAL     LUTUSR(3,*),R(*),G(*),B(*)
C
      IF (IPLOT.EQ.2) THEN
        Z=0.0
        DZ=0.999/FLOAT(NCOLS-1)
        DO 10 I=1,NCOLS
          R(I)=Z
          G(I)=Z
          B(I)=Z
          Z=Z+DZ
  10    CONTINUE
      ELSEIF (IPLOT.EQ.3) THEN
        CALL STCOL1(HEAT,R,G,B,NCOLS)
      ELSEIF (IPLOT.EQ.4) THEN
        CALL STCOL1(SPECTRUM,R,G,B,NCOLS)
      ELSEIF (IPLOT.EQ.5) THEN
        CALL STCOL1(BGYRW,R,G,B,NCOLS)
      ELSEIF (IPLOT.EQ.6) THEN
        CALL STCOL1(SERP,R,G,B,NCOLS)
      ELSE
        CALL STCOL1(LUTUSR,R,G,B,NCOLS)
      ENDIF
      END
C
      SUBROUTINE STCOL1(LUT,R,G,B,N)
C     ------------------------------
C
      REAL LUT(3,*),R(*),G(*),B(*)
C
      DO 10 I=1,N
        R(I)=LUT(1,I)
        G(I)=LUT(2,I)
        B(I)=LUT(3,I)
  10  CONTINUE
      END
C
C***<some untility routines>********************************************
C
      SUBROUTINE VCOPY(X,Y,N)
C     -----------------------
C
      REAL X(*),Y(*)
C
      DO 10 I=1,N
  10    Y(I)=X(I)
      END
C
      SUBROUTINE VRFILL(X,A,N)
C     ------------------------
C
      REAL X(*)
C
      DO 10 I=1,N
  10    X(I)=A
      END
C
      SUBROUTINE LOGQYN(S,D,L)
C     ------------------------
C
      LOGICAL      L,L2
      CHARACTER*1  D,D2,A
      CHARACTER*45 STRING
      CHARACTER    S(*)
C
      IF (D.EQ.'Y') THEN
        L=.TRUE.
        D2='N'
        L2=.FALSE.
      ELSEIF (D.EQ.'N') THEN
        L=.FALSE.
        D2='Y'
        L2=.TRUE.
      ELSE
        WRITE(*,*)' Default should be Y or N !'
        RETURN
      ENDIF
      CALL STPACK(STRING,S,45)
   1  WRITE(*,100) STRING,D
 100  FORMAT(A45,' Y/N (def=',A1,')  : ',$)
      CALL FORMTQ(A,1,NNN)
      IF (NNN.EQ.0) RETURN
      IF (A.EQ.'y' .OR. A.EQ.'T' .OR. A.EQ.'t') A='Y'
      IF (A.EQ.'n' .OR. A.EQ.'F' .OR .A.EQ.'f') A='N'
      IF (A.EQ.D) THEN
        RETURN
      ELSEIF (A.EQ.D2) THEN
        L=L2
        RETURN
      ENDIF
      GOTO 1
      END
C
      SUBROUTINE STPACK(STRING,S,N)
C     -----------------------------
C
      CHARACTER STRING(*),S(*)
C
      DO 10 I=1,N
        STRING(I)=S(I)
        IF (S(I).EQ.'?') GOTO 20
  10  CONTINUE
  20  DO 30 J=I+1,N
  30    STRING(J)=' '
      END
C
      SUBROUTINE FORMTQ(STRING,NCHAR,NNN)
C     -----------------------------------
C
      CHARACTER*(*) STRING
C
      READ(*,200) STRING
 200  FORMAT(A)
      NNN=0
      DO 10 I=1,NCHAR
  10    IF (STRING(I:I).NE.' ') NNN=NNN+1
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
