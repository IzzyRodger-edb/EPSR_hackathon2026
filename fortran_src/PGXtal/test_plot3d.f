CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      PROGRAM TEST_PLOT3D
C     -------------------
C
C  A little program to illustrate the visualisation of a 3D density 
C  density distribution: f=f(x,y,z).
C
      INCLUDE   'COLTABS.INC'
      PARAMETER (NX=20,NY=20,NZ=20)
      REAL       F(-NX:NX,-NY:NY,-NZ:NZ)
      REAL       EYE(3),LATICE(3,4,2)
      REAL       SLNORM(3),APOINT(3)
      REAL       LIGHT(3),RGBBKG(3),RGBSRF(3)
      LOGICAL    ANOTHR
      DATA       A,B,C /10.0,10.0,10.0/
      DATA       ALFA,BETA,GAMA /90.0,90.0,90.0/
      DATA       EYE /0.0,0.0,1000.0/
      DATA       LIGHT /-1.0,-1.0,-1.0/
      DATA       RGBBKG /0.25,0.25,0.25/
      DATA       RGBSRF /1.00,0.00,0.00/
      DATA       IPHI1,ITHT,IPHI2 /90,70,75/
      DATA       PI /3.141592654/
C
      RADANG=PI/180.0
      CALL FILARY(F,NX,NY,NZ,FMIN,FMAX)
      CALL AXINIT(A,B,C,ALFA,BETA,GAMA,LATICE,RADANG,ZSHIFT)
      N1=2*NX
      N2=2*NY
      N3=2*NY
C
      FSURF=FMIN+0.21*(FMAX-FMIN)
      F1=FMIN
      F2=MIN(2.0*FSURF,FMAX)
   1  CALL VIEW3A(IPHI1,ITHT,IPHI2,LATICE,RADANG,ZSHIFT)
      CALL VIEW3B(LATICE,EYE,WMAX)
      CALL VIEW3C(LATICE,SLNORM,APOINT)
C
      CALL PGBEGIN(0,'?',1,1)
      CALL PGPAPER(0.0,1.0)
      CALL PGVPORT(0.0,1.0,0.0,1.0)
      CALL PGWINDOW(-WMAX,WMAX,-WMAX,WMAX)
      CALL PGBOX('BC',0.0,0,'BC',0.0,0)
      CALL SBFINT(RGBBKG,4,1,1,MAXBUF)
      CALL EDGES(EYE,LATICE,1)
      CALL COLINT(RGBSRF,16,39,0.5,0.0,1.0)
      CALL COLTAB(SPECTRUM,256,0.7,40,63)
      CALL SBSURF(EYE,LATICE,F,N1,N2,N3,FSURF,16,39,LIGHT,.FALSE.)
      CALL SBSLIC(EYE,LATICE,F,N1,N2,N3,F1,F2,40,63,SLNORM,APOINT,7)
      CALL SBFCLS(1)
      CALL PGBOX('BC',0.0,0,'BC',0.0,0)
      CALL PGEND
C
      CALL LOGQYN(' Plot> Another view ?','Y',ANOTHR)
      IF (ANOTHR) GOTO 1
      END
C
      SUBROUTINE EDGES(EYE,LATICE,ICOL)
C     ---------------------------------
C
      REAL EYE(*),LATICE(3,*),END1(3),END2(3)
C
      II=0
      DO 10 L=1,12,4
        I=II+2
        J=MOD(II+1,3)+2
        K=MOD(II+2,3)+2
        CALL SBLINE(EYE,LATICE,LATICE(1,K),ICOL,.FALSE.)
        END1(1)=LATICE(1,I)+LATICE(1,J)-LATICE(1,1)
        END1(2)=LATICE(2,I)+LATICE(2,J)-LATICE(2,1)
        END1(3)=LATICE(3,I)+LATICE(3,J)-LATICE(3,1)
        CALL SBLINE(EYE,LATICE(1,I),END1,ICOL,.FALSE.)
        CALL SBLINE(EYE,LATICE(1,J),END1,ICOL,.FALSE.)
        END2(1)=END1(1)+LATICE(1,K)-LATICE(1,1)
        END2(2)=END1(2)+LATICE(2,K)-LATICE(2,1)
        END2(3)=END1(3)+LATICE(3,K)-LATICE(3,1)
        CALL SBLINE(EYE,END1,END2,ICOL,.FALSE.)
        II=II+1
  10  CONTINUE
      END
C
      SUBROUTINE FILARY(F,NX,NY,NZ,FMIN,FMAX)
C     ---------------------------------------
C
      REAL F(-NX:NX,-NY:NY,-NZ:NZ)
      DATA SMALL,RMAX,THREE /1.0E-10,10.0,3.0/
C
      FMAX=0.0
      FMIN=1.0E20
      DX=RMAX/FLOAT(NX)
      DY=RMAX/FLOAT(NY)
      DZ=RMAX/FLOAT(NZ)
      Z=-RMAX
      DO 30 K=-NZ,NZ
        ZZ=Z**2
        Y=-RMAX
        DO 20 J=-NY,NY
          YY=Y**2
          X=-RMAX
          DO 10 I=-NX,NX
            RR=X**2+YY+ZZ+SMALL
            R=SQRT(RR)
            PSI2=EXP(-R)*(THREE*ZZ-RR)**2
            F(I,J,K)=PSI2
            FMAX=MAX(PSI2,FMAX)
            FMIN=MIN(PSI2,FMIN)
            X=X+DX
  10      CONTINUE
          Y=Y+DY
  20    CONTINUE
        Z=Z+DZ
  30  CONTINUE
      END
C
      SUBROUTINE AXINIT(A,B,C,ALFA,BETA,GAMA,LATICE,RADANG,ZSHIFT)
C     ------------------------------------------------------------
C
      REAL LATICE(3,4,*)
C
      ALF=RADANG*ALFA
      BET=RADANG*BETA
      GAM=RADANG*GAMA
      AX=A
      AY=0.0
      AZ=0.0
      BX=B*COS(GAM)
      BY=B*SIN(GAM)
      BZ=0.0
      CX=C*COS(BET)
      CY=(C*COS(ALF)-CX*COS(GAM))/SIN(GAM)
      CZ=SQRT(C**2-CX**2-CY**2)
      XSHIFT=-0.5*(AX+BX+CX)
      YSHIFT=-0.5*(AY+BY+CY)
      ZSHIFT=-0.5*(AZ+BZ+CZ)
      LATICE(1,1,2)=XSHIFT
      LATICE(2,1,2)=YSHIFT
      LATICE(3,1,2)=ZSHIFT
      LATICE(1,2,2)=AX+XSHIFT
      LATICE(2,2,2)=AY+YSHIFT
      LATICE(3,2,2)=AZ+ZSHIFT
      LATICE(1,3,2)=BX+XSHIFT
      LATICE(2,3,2)=BY+YSHIFT
      LATICE(3,3,2)=BZ+ZSHIFT
      LATICE(1,4,2)=CX+XSHIFT
      LATICE(2,4,2)=CY+YSHIFT
      LATICE(3,4,2)=CZ+ZSHIFT
      ZSHIFT=-MAX(A,B,C)
      END
C
      SUBROUTINE VIEW3A(JPHI1,JTHT,JPHI2,LATICE,RADANG,ZSHIFT)
C     --------------------------------------------------------
C
      REAL      LATICE(3,4,*)
      CHARACTER STRING*60
C
      WRITE(*,*)
   1  WRITE(*,100) JPHI1,JTHT,JPHI2
 100  FORMAT(' View> Euler angles ?  (def=',I2,',',I2,',',I2,
     *       ' deg) : ',$)
      CALL FORMTQ(STRING,60,NNN)
      IF (NNN.NE.0) THEN
        READ(STRING,*,ERR=1) IPHI1,ITHT,IPHI2
        IF (ITHT.LT.0 .OR. ITHT.GT.180) GOTO 1
        IF (IPHI1.LT.0 .OR. IPHI1.GT.360) GOTO 1
        IF (IPHI2.LT.0 .OR. IPHI2.GT.360) GOTO 1
        JTHT=ITHT
        JPHI1=IPHI1
        JPHI2=IPHI2
      ENDIF
      THT=FLOAT(JTHT)*RADANG
      PHI1=FLOAT(JPHI1)*RADANG
      PHI2=FLOAT(JPHI2)*RADANG
      SINTHT=-SIN(THT)
      COSTHT=COS(THT)
      SINPH1=-SIN(PHI1)
      COSPH1=COS(PHI1)
      SINPH2=SIN(PHI2)
      COSPH2=COS(PHI2)
      DO 10 J=1,4
        CALL ROTZ(LATICE(1,J,2),LATICE(2,J,2),LATICE(3,J,2),U,V,W,
     *       SINPH2,COSPH2)
        CALL ROTY(U,V,W,X,Y,Z,SINTHT,COSTHT)
        CALL ROTZ(X,Y,Z,LATICE(1,J,1),LATICE(2,J,1),LATICE(3,J,1),
     *       SINPH1,COSPH1)
        LATICE(3,J,1)=LATICE(3,J,1)+ZSHIFT
  10  CONTINUE
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
      SUBROUTINE VIEW3B(LATICE,EYE,WMAX)
C     ----------------------------------
C
      REAL LATICE(3,*),EYE(*)
C
      WMAX=1.0E-10
      DO 10 J=1,4
        CALL SBLIN1(EYE,LATICE(1,J),LATICE(2,J),LATICE(3,J),XW,YW)
        WMAX=MAX(ABS(XW),ABS(YW),WMAX)
  10  CONTINUE
      WMAX=1.05*WMAX
      END
C
      SUBROUTINE VIEW3C(LATICE,SLNORM,APOINT)
C     ---------------------------------------
C
      REAL LATICE(3,*),SLNORM(*),APOINT(*)
C
      DO 10 I=1,3
        APOINT(I)=0.5*(LATICE(I,4)+LATICE(I,3)+LATICE(I,2)-LATICE(I,1))
        SLNORM(I)=LATICE(I,2)-LATICE(I,1)
  10  CONTINUE
      END
C
C***<logical question utitlity>*****************************************
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
