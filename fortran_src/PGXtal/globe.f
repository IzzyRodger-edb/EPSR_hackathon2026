CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C***<sphere-specific stuff>*********************************************
C
      SUBROUTINE GLOBE(Z,THTMIN,PHIMIN,NTHT,NPHI,DTHT,DPHI,SIZE,IWIDTH,
     *                 XLBL,YLBL,TITLE)
C     -----------------------------------------------------------------
C
      COMMON /TRCOM/ THT1,PHI1,N2,N1,DCLAT,DLONG,RXY2,DEYE,DSTR,
     *               DSCRN,X0(3),Y0(3),Z0(3)
      REAL           Z(NPHI,NTHT)
      CHARACTER*(*)  XLBL,YLBL,TITLE
      DATA           RADANG /57.2957795/
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      This is a "front-end" for the globe-plotting routine PLTSPH. All
C    the angle parameters should be passed down in degrees.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    Z        R*4    I  NPHI x NTHT The recatangular "data"-array.
C    Z        R*4    O  NPHI x NTHT A scaled, and clipped, version of 
C                                   the input array(!).
C    THTMIN   R*4    I       -      The min. value of co-latitude theta.
C    PHIMIN   R*4    I       -      The min. value of longitude phi.
C    NPHI     I*4    I       -      The first dimension of array Z.
C    NTHT     I*4    I       -      The second dimension of array Z.
C    DTHT     R*4    I       -      The theta increment for array Z.
C    DPHI     R*4    I       -      The phi increment for array Z.
C    SIZE     R*4    I       -      Character-size for plot.
C    IWIDTH   I*4    I       -      Line-width for plot.
C    XLBL     C*(*)  I       -      The x-axis label for the plot.
C    YLBL     C*(*)  I       -      The y-axis label for the plot.
C    TITLE    C*(*)  I       -      The title for the plot.
C
C History
C   D. S. Sivia      21 Feb 1995  Initial release.
C   D. S. Sivia      27 Feb 1995  Corrected minor problem with contours.
C   D. S. Sivia      27 Aug 1997  No longer use non-standard Q-format.
C   D. S. Sivia      15 Sep 1997  Fortran made LINUX-friendly!
C   D. S. Sivia      24 Sep 1997  Put in extra square-root protection!
C-----------------------------------------------------------------------
C
      THT1=THTMIN/RADANG
      PHI1=PHIMIN/RADANG
      N1=NPHI
      N2=NTHT
      DCLAT=DTHT/RADANG
      DLONG=DPHI/RADANG
      CALL TRINIT(RXY,RXY2,DEYE,DSTR,DSCRN,X0,Y0,Z0,RADANG)
      CALL PLTSPH(Z,NPHI,NTHT,1,NPHI,1,NTHT,RXY,SIZE,IWIDTH,
     *            XLBL,YLBL,TITLE)
      END
C
      SUBROUTINE TRINIT(RXY,RXY2,DEYE,DSTR,DSCRN,X0,Y0,Z0,RADANG)
C     -----------------------------------------------------------
C
      REAL         X0(*),Y0(*),Z0(*)
      CHARACTER*64 STRING
C
      WRITE(*,*)
   1  DEYE=10.0
      WRITE(*,100)
 100  FORMAT(' Eye >  Distance from unit sphere ?  (def=10.0) : ',$)
      CALL FORMTQ(STRING,64,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=1) DEYE
      IF (DEYE.LT.1.0E-5) GOTO 1
      DEYE=DEYE+1.01
      DSCRN=1.0/DEYE
      DSTR=DEYE-DSCRN
      RXY2=1.0-1.0/DEYE**2
      RXY=SQRT(ABS(RXY2)+1.0E-20)
   2  IPHI1=90
      ITHT=60
      IPHI2=135
      WRITE(*,110)
 110  FORMAT(' Euler >  Phi1, theta, Phi2 ?   (def=90,60,135) : ',$)
      CALL FORMTQ(STRING,64,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=2) IPHI1,ITHT,IPHI2
      IF (ITHT.LT.0 .OR. ITHT.GT.180) GOTO 2
      IF (IPHI1.LT.0 .OR. IPHI1.GT.360) GOTO 2
      IF (IPHI2.LT.0 .OR. IPHI2.GT.360) GOTO 2
      THT=FLOAT(ITHT)/RADANG
      PHI1=FLOAT(IPHI1)/RADANG
      PHI2=FLOAT(IPHI2)/RADANG
      SINTHT=-SIN(THT)
      COSTHT=COS(THT)
      SINPH1=-SIN(PHI1)
      COSPH1=COS(PHI1)
      SINPH2=SIN(PHI2)
      COSPH2=COS(PHI2)
      CALL ROTY(COSPH2,SINPH2,0.0,UX,VX,WX,SINTHT,COSTHT)
      CALL ROTY(-SINPH2,COSPH2,0.0,UY,VY,WY,SINTHT,COSTHT)
      CALL ROTY(0.0,0.0,1.0,UZ,VZ,WZ,SINTHT,COSTHT)
      CALL ROTZ(UX,VX,WX,X0(1),X0(2),X0(3),SINPH1,COSPH1)
      CALL ROTZ(UY,VY,WY,Y0(1),Y0(2),Y0(3),SINPH1,COSPH1)
      CALL ROTZ(UZ,VZ,WZ,Z0(1),Z0(2),Z0(3),SINPH1,COSPH1)
      END
C
      SUBROUTINE ROTY(X,Y,Z,U,V,W,S,C)
C     --------------------------------
C
      U=X*C-Z*S
      V=Y
      W=X*S+Z*C
      END
C
      SUBROUTINE ROTZ(X,Y,Z,U,V,W,S,C)
C     --------------------------------
C
      U=X*C+Y*S
      V=-X*S+Y*C
      W=Z
      END
C
      SUBROUTINE TRNL(XW,YW,XI,YJ)
C     ----------------------------
C
      COMMON /TRCOM/ THTMIN,PHIMIN,NTHT,NPHI,DTHT,DPHI,RXY2,DEYE,DSTR,
     *               DSCRN,X0(3),Y0(3),Z0(3)
      DATA           SMALL,SAFE,TWOPI /1.0E-10,0.9999,6.283185/
      DATA           SMALL2 /1.0E-20/
C
      R2=XW*XW+YW*YW
      IF (R2.GT.RXY2) THEN
        XI=-1.0
        YJ=-1.0
        RETURN
      ELSEIF (R2.LT.SMALL) THEN
        X=0.0
        Y=0.0
        Z=1.0
      ELSE
        DSTR2=DSTR*DSTR
        COSBT2=DSTR2/(R2+DSTR2)
        SINBT2=1.0-COSBT2
        B=DEYE*SINBT2
        Z=B+SQRT(ABS(B*(B-DEYE)+COSBT2)+SMALL2)
        ALPHA=SQRT(ABS((1.0-Z*Z)/R2)+SMALL2)
        X=ALPHA*XW
        Y=ALPHA*YW
      ENDIF
      U=X*X0(1)+Y*X0(2)+Z*X0(3)
      V=X*Y0(1)+Y*Y0(2)+Z*Y0(3)
      W=X*Z0(1)+Y*Z0(2)+Z*Z0(3)
      IF (W.GE.SAFE) THEN
        XI=1.001
        YJ=1.001
      ELSEIF (W.LE.-SAFE) THEN
        XI=1.001
        YJ=FLOAT(NTHT)+0.999
      ELSE     
        PHI=ATAN2(V,U)
        IF (PHI.LT.0.0) PHI=PHI+TWOPI
        THT=ACOS(W)
        XI=1.0+(PHI-PHIMIN)/DPHI
        YJ=1.0+(THT-THTMIN)/DTHT
      ENDIF
      END
C
      SUBROUTINE TRNLB(XW,YW,XI,YJ)
C     -----------------------------
C
      COMMON /TRCOM/ THTMIN,PHIMIN,NTHT,NPHI,DTHT,DPHI,RXY2,DEYE,DSTR,
     *               DSCRN,X0(3),Y0(3),Z0(3)
C
      PHI=PHIMIN+DPHI*(XI-1.0)
      THT=THTMIN+DTHT*(YJ-1.0)
      SINPHI=SIN(PHI)
      COSPHI=COS(PHI)
      SINTHT=SIN(THT)
      COSTHT=COS(THT)
      CALL TRNLC(XW,YW,SINPHI,COSPHI,SINTHT,COSTHT)
      END
C
      SUBROUTINE TRNLC(XW,YW,SINPHI,COSPHI,SINTHT,COSTHT)
C     ---------------------------------------------------
C
      COMMON /TRCOM/ THTMIN,PHIMIN,NTHT,NPHI,DTHT,DPHI,RXY2,DEYE,DSTR,
     *               DSCRN,X0(3),Y0(3),Z0(3)
C
      UX=SINTHT*COSPHI
      UY=SINTHT*SINPHI
      UZ=COSTHT
      Z=UX*X0(3)+UY*Y0(3)+UZ*Z0(3)
      IF (Z.LT.DSCRN) THEN
        XW=100.0
        YW=100.0
      ELSE
        X=UX*X0(1)+UY*Y0(1)+UZ*Z0(1)
        Y=UX*X0(2)+UY*Y0(2)+UZ*Z0(2)
        BETA=DSTR/(DEYE-Z)
        XW=X*BETA
        YW=Y*BETA
      ENDIF
      END
C
C***<plot2dx>***********************************************************
C
      SUBROUTINE PLTSPH(Z,N1,N2,I1,I2,J1,J2,RADIUS,SIZE,IWIDTH,
     *                  XLBL,YLBL,TITLE)
C     ---------------------------------------------------------
C
      REAL          Z(N1,N2)
      REAL          W(2050),X(361),Y(361),SINA(361),COSA(361)
      CHARACTER*(*) XLBL,YLBL,TITLE
      LOGICAL       OVRLAY
      CHARACTER*32  STRING
      EXTERNAL      TRNL,TRNLB,TRNLC
      DATA          NANGLE /361/
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      This is a globe-specific version of PLT2DX. It takes "data" 
C    defined on a rectangular grid of points Z(I,J), corresponding to
C    latitude and longitude, and plots out a linearly-interpolated map 
C    on the surface of a sphere.
C      The appropriate Jacobian term should be factored into the array 
C    Z(I,J) before it is passed down (by pre-multilication). 
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    Z        R*4    I    N1 x N2   The recatangular "data"-array.
C    Z        R*4    O    N1 x N2   A scaled, and clipped, version of 
C                                   the input array(!).
C    N1       I*4    I       -      The first dimension of array Z.
C    N2       I*4    I       -      The second dimension of array Z.
C    I1,I2    I*4    I       -      The inclusive range of the first
C                                   index (I) to be plotted.
C    J1,J2    I*4    I       -      The inclusive range of the second
C                                   index (J) to be plotted.
C    RADIUS   R*4    I       -      Radius for plot, in output units.
C    SIZE     R*4    I       -      Character-size for plot.
C    IWIDTH   I*4    I       -      Line-width for plot.
C    XLBL     C*(*)  I       -      The x-axis label for the plot.
C    YLBL     C*(*)  I       -      The y-axis label for the plot.
C    TITLE    C*(*)  I       -      The title for the plot.
C
C Globals
C    COLTABS.INC
C
C History
C   D. S. Sivia      21 Feb 1995  Modified PLT2DX to PLTSPH.
C-----------------------------------------------------------------------
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
        CALL CONTOR(Z,N1,N2,I1,I2,J1,J2,RADIUS,SIZE,IWIDTH,.FALSE.)
      ELSEIF (IPLOT.EQ.1) THEN
        WRITE(*,*)
        WRITE(*,*)' Sorry, but this option is not (yet) operational!'
        WRITE(*,*)
        GOTO 1
      ELSEIF (IPLOT.LE.7) THEN
        CALL GREY(Z,N1,N2,I1,I2,J1,J2,RADIUS,IPLOT,SIZE,IWIDTH)
      ELSE
        GOTO 1
      ENDIF
      CALL PGBBUF
      CALL CIRCL(RADIUS,X,Y,SINA,COSA,NANGLE)
      CALL PGSLW(1)
      CALL LATTS(RADIUS,X,Y,SINA,COSA,NANGLE,TRNLC)
      CALL LONGS(RADIUS,X,Y,SINA,COSA,NANGLE,TRNLC)
      CALL PGSLW(IWIDTH)
      CALL PGEBUF
      CALL PGLABEL(XLBL,YLBL,TITLE)
      IF (IPLOT.GT.1) THEN
        WRITE(*,*)
        CALL LOGQYN(' PLOT> Over-lay contours ?','N',OVRLAY)
        IF (OVRLAY) 
     *    CALL CONTOR(Z,N1,N2,I1,I2,J1,J2,RADIUS,SIZE,IWIDTH,OVRLAY)
      ENDIF
      CALL PGEND
      END
C
C***<contour>***********************************************************
C
      SUBROUTINE CONTOR(MAP,N1,N2,I1,I2,J1,J2,RADIUS,SIZE,IWIDTH,OVRLAY)
C     ------------------------------------------------------------------
C
      REAL     MAP(*)
      REAL     CONT(25,2)
      INTEGER  NCONT(2),LCOLOR(2),LSTYLE(2),LWIDTH(2)
      LOGICAL  OVRLAY
      EXTERNAL PLOT,TRNLB
      DATA     LCOLOR,LSTYLE,LWIDTH,NWIND,IDP /1,3,1,1,1,2,1,0/
C
      CALL INIT(MAP,DMIN,DMAX,N1*N2)
      IF (OVRLAY) THEN
        CALL AUTCNT(NCONT,CONT,DMIN,DMAX,OVRLAY)
        CALL PGSLW(1)
        CALL PGCONX(MAP,N1,N2,I1,I2,J1,J2,CONT,NCONT,PLOT)
        CALL PGSLW(IWIDTH)
        RETURN
      ENDIF
      CALL ASKCNT(NCONT,CONT,DMIN,DMAX)
      CALL PGBEGIN(0,'?',1,1)
      CALL PGPAPER(0.0,1.0)
      CALL PGSCH(SIZE)
      CALL PGSLW(IWIDTH)
      CALL PGENV(-RADIUS,RADIUS,-RADIUS,RADIUS,0,-2)
      DO 20 I=1,2
        CALL PGSCI(LCOLOR(I))
        CALL PGSLS(LSTYLE(I))
        LWDTH=2*LWIDTH(I)
        IF (LWDTH.GT.7) LWDTH=7
        CALL PGSLW(LWDTH)
        CALL PGCONX(MAP,N1,N2,I1,I2,J1,J2,CONT(1,I),NCONT(I),PLOT)
  20  CONTINUE
      CALL PGSCI(1)
      CALL PGSLW(IWIDTH)
      END
C
      SUBROUTINE PLOT(IVISIB,XI,YJ,ZIJ)
C     ---------------------------------
C
      EXTERNAL TRNLB
      LOGICAL  LOUT1,LOUT2
C
      LOUT2=.FALSE.
      CALL TRNLB(XW,YW,XI,YJ)
      IF (XW.GT.1.0 .OR. YW.GT.1.0) LOUT2=.TRUE.
      IF (IVISIB.EQ.0 .OR. LOUT1 .OR. LOUT2) THEN
        CALL PGMOVE(XW,YW)
      ELSE
        CALL PGDRAW(XW,YW)
      ENDIF
      LOUT1=LOUT2
      END
C
      SUBROUTINE INIT(MAP,DMIN,DMAX,NMAP)
C     -----------------------------------
C
      REAL MAP(*)
C
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
      SUBROUTINE GREY(MAP,N1,N2,I1,I2,J1,J2,RADIUS,IPLOT,SIZE,IWIDTH)
C     ---------------------------------------------------------------
C
      REAL     MAP(*)
      REAL     TRCOL(6),COLBAR(2,256),LUTUSR(3,256)
      REAL     RED(256),GR(256),BL(256)
      EXTERNAL TRNL
      DATA     TRCOL /0.0,1.0,0.0,0.0,0.0,1.0/
C
      NCOLS=256
      CALL INIT(MAP,DMIN,DMAX,N1*N2)
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
      CALL PGWINDOW(-RADIUS,RADIUS,-RADIUS,RADIUS)
      CALL PGCELX(MAP,N1,N2,I1,I2,J1,J2,1.0,0.0,TRNL,NCOLS,RED,GR,BL)
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
C***<circlular border, latitudes and longitudes>************************
C
      SUBROUTINE CIRCL(R,X,Y,SINA,COSA,N)
C     -----------------------------------
C
      REAL X(*),Y(*),SINA(*),COSA(*)
C
      DANG=2.0*3.141592654/FLOAT(N-1)
      ANG=0.0
      DO 10 I=1,N
        SINA(I)=SIN(ANG)
        COSA(I)=COS(ANG)
        X(I)=R*COSA(I)
        Y(I)=R*SINA(I)
        ANG=ANG+DANG
  10  CONTINUE
      CALL PGLINE(N,X,Y)
      END
C
      SUBROUTINE LATTS(R,X,Y,SINA,COSA,N,TRNLC)
C     -----------------------------------------
C
      REAL     X(*),Y(*),SINA(*),COSA(*)
      EXTERNAL TRNLC
C
      DTHT=3.141592654/6.0
      THT=DTHT
      DO 20 J=1,5
        SINTHT=SIN(THT)
        COSTHT=COS(THT)
        K=0
        DO 10 I=1,N
          CALL TRNLC(XW,YW,SINA(I),COSA(I),SINTHT,COSTHT)
          IF (XW.GT.R) THEN
            IF (K.GT.1) CALL PGLINE(K,X,Y)
            K=0
          ELSE
            K=K+1
            X(K)=XW
            Y(K)=YW
          ENDIF
  10    CONTINUE
        IF (K.GT.1) CALL PGLINE(K,X,Y)
        THT=THT+DTHT
  20  CONTINUE
      END
C
      SUBROUTINE LONGS(R,X,Y,SINA,COSA,N,TRNLC)
C     -----------------------------------------
C
      REAL     X(*),Y(*),SINA(*),COSA(*)
      EXTERNAL TRNLC
C
      DPHI=3.141592654/6.0
      PHI=0.0
      DO 20 J=1,6
        SINPHI=SIN(PHI)
        COSPHI=COS(PHI)
        K=0
        DO 10 I=1,N
          CALL TRNLC(XW,YW,SINPHI,COSPHI,SINA(I),COSA(I))
          IF (XW.GT.R) THEN
            IF (K.GT.1) CALL PGLINE(K,X,Y)
            K=0
          ELSE
            K=K+1
            X(K)=XW
            Y(K)=YW
          ENDIF
  10    CONTINUE
        IF (K.GT.1) CALL PGLINE(K,X,Y)
        PHI=PHI+DPHI
  20  CONTINUE
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
