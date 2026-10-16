 c from James Stewart:
 c MOPAC has a subroutine XYZINT that can form the core of a program
 c to perform this conversion.  Feel free to use it.
 c from Frank Jensen (frj \\at// dou.dk):
 c Dr. Jensen has put together a short program to perform this conversion
 c by assembling subroutines from MOPAC.  The text follows:
       PROGRAM carint
       IMPLICIT REAL*8(A-H,O-Z)
 C     PROGRAM TO CONVERT cartesian coordinates to z-matrix
       character*4 test
       DIMENSION xyz(3,100),na(100),nb(100),nc(100),geo(3,100),iz(100)
       DEGREE = 57.29578d+00
 c     open(unit=5,form='formatted',status='old',file='xx')
 c     open(unit=6,form='formatted',status='unknown')
       i=0
 30    read(5,*,end=40,err=40)x,y,z
       i=i+1
       xyz(1,i)=x
       xyz(2,i)=y
       xyz(3,i)=z
       goto 30
 40    continue
       numat=i
       call XYZINT(XYZ,NUMAT,NA,NB,NC,DEGREE,GEO)
       do 50 i=1,3
 50    write(6,*)' '
       do 60 i=1,numat
 	 if (geo(3,i).gt.180.d0) geo(3,i)=geo(3,i)-360.d0
          write(6,103)na(i),geo(1,i),nb(i),geo(2,i),nc(i),geo(3,i)
 60    continue
       write(6,*)' '
 110   format(a4,3f15.10)
 103   format(i5,f12.5,i5,f12.5,i5,f12.5)
 200   stop
       END
       SUBROUTINE DIHED(XYZ,I,J,K,L,ANGLE)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION XYZ(3,*)
 *********************************************************************
 *
 *      DIHED CALCULATES THE DIHEDRAL ANGLE BETWEEN ATOMS I, J, K,
 *            AND L.  THE CARTESIAN COORDINATES OF THESE ATOMS
 *            ARE IN ARRAY XYZ.
 *
 *     DIHED IS A MODIFIED VERSION OF A SUBROUTINE OF THE SAME NAME
 *           WHICH WAS WRITTEN BY DR. W. THEIL IN 1973.
 *
 *********************************************************************
       XI1=XYZ(1,I)-XYZ(1,K)
       XJ1=XYZ(1,J)-XYZ(1,K)
       XL1=XYZ(1,L)-XYZ(1,K)
       YI1=XYZ(2,I)-XYZ(2,K)
       YJ1=XYZ(2,J)-XYZ(2,K)
       YL1=XYZ(2,L)-XYZ(2,K)
       ZI1=XYZ(3,I)-XYZ(3,K)
       ZJ1=XYZ(3,J)-XYZ(3,K)
       ZL1=XYZ(3,L)-XYZ(3,K)
 C      ROTATE AROUND Z AXIS TO PUT KJ ALONG Y AXIS
       DIST= SQRT(XJ1**2+YJ1**2+ZJ1**2)
       COSA=ZJ1/DIST
       IF(COSA.GT.1.0D0) COSA=1.0D0
       IF(COSA.LT.-1.0D0) COSA=-1.0D0
       DDD=1.0D0-COSA**2
       IF(DDD.LE.0.0) GO TO 10
       YXDIST=DIST* SQRT(DDD)
       IF(YXDIST.GT.1.0D-9) GO TO 20
    10 CONTINUE
       XI2=XI1
       XL2=XL1
       YI2=YI1
       YL2=YL1
       COSTH=COSA
       SINTH=0.D0
       GO TO 30
    20 COSPH=YJ1/YXDIST
       SINPH=XJ1/YXDIST
       XI2=XI1*COSPH-YI1*SINPH
       XJ2=XJ1*COSPH-YJ1*SINPH
       XL2=XL1*COSPH-YL1*SINPH
       YI2=XI1*SINPH+YI1*COSPH
       YJ2=XJ1*SINPH+YJ1*COSPH
       YL2=XL1*SINPH+YL1*COSPH
 C      ROTATE KJ AROUND THE X AXIS SO KJ LIES ALONG THE Z AXIS
       COSTH=COSA
       SINTH=YJ2/DIST
    30 CONTINUE
       YI3=YI2*COSTH-ZI1*SINTH
       YL3=YL2*COSTH-ZL1*SINTH
       CALL DANG(XL2,YL3,XI2,YI3,ANGLE)
       IF (ANGLE .LT. 0.) ANGLE=6.2831853D0+ANGLE
       IF (ANGLE .GE. 6.2831853D0 ) ANGLE=0.D0
       RETURN
       END
       SUBROUTINE BANGLE(XYZ,I,J,K,ANGLE)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION XYZ(3,*)
 *********************************************************************
 *
 * BANGLE CALCULATES THE ANGLE BETWEEN ATOMS I,J, AND K. THE
 *        CARTESIAN COORDINATES ARE IN XYZ.
 *
 *********************************************************************
       D2IJ = (XYZ(1,I)-XYZ(1,J))**2+
      1       (XYZ(2,I)-XYZ(2,J))**2+
      2       (XYZ(3,I)-XYZ(3,J))**2
       D2JK = (XYZ(1,J)-XYZ(1,K))**2+
      1       (XYZ(2,J)-XYZ(2,K))**2+
      2       (XYZ(3,J)-XYZ(3,K))**2
       D2IK = (XYZ(1,I)-XYZ(1,K))**2+
      1       (XYZ(2,I)-XYZ(2,K))**2+
      2       (XYZ(3,I)-XYZ(3,K))**2
       XY = SQRT(D2IJ*D2JK)
       TEMP = 0.5D0 * (D2IJ+D2JK-D2IK) / XY
       IF (TEMP .GT. 1.0D0) TEMP=1.0D0
       IF (TEMP .LT. -1.0D0) TEMP=-1.0D0
       ANGLE = ACOS( TEMP )
       RETURN
       END
       SUBROUTINE DANG(A1,A2,B1,B2,RCOS)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 **********************************************************************
 *
 *    DANG  DETERMINES THE ANGLE BETWEEN THE POINTS (A1,A2), (0,0),
 *          AND (B1,B2).  THE RESULT IS PUT IN RCOS.
 *
 **********************************************************************
       PI=2.0D0* ASIN(1.0D00)
       ZERO=1.0D-6
       IF( ABS(A1).LT.ZERO.AND. ABS(A2).LT.ZERO) GO TO 10
       IF( ABS(B1).LT.ZERO.AND. ABS(B2).LT.ZERO) GO TO 10
       ANORM=1.0D0/ SQRT(A1**2+A2**2)
       BNORM=1.0D0/ SQRT(B1**2+B2**2)
       A1=A1*ANORM
       A2=A2*ANORM
       B1=B1*BNORM
       B2=B2*BNORM
       SINTH=(A1*B2)-(A2*B1)
       COSTH=A1*B1+A2*B2
       IF(COSTH.GT.1.0D0) COSTH=1.0D0
       IF(COSTH.LT.-1.0D0) COSTH=-1.0D0
       RCOS= ACOS(COSTH)
       IF( ABS(RCOS).LT.4.0D-4) GO TO 10
       IF(SINTH.GT.0.D0) RCOS=6.2831853D0-RCOS
       RCOS=-RCOS
       RETURN
    10 RCOS=0.0D0
       RETURN
       END
       SUBROUTINE XYZGEO(XYZ,NUMAT,NA,NB,NC,DEGREE,GEO)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION XYZ(3,*), NA(*), NB(*), NC(*), GEO(3,*)
 ***********************************************************************
 *
 *   XYZGEO CONVERTS COORDINATES FROM CARTESIAN TO INTERNAL.
 *
 *     ON INPUT XYZ  = ARRAY OF CARTESIAN COORDINATES
 *              NUMAT= NUMBER OF ATOMS
 *              NA   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
 *                     BY DISTANCE
 *              NB   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
 *                     BY ANGLE
 *              NC   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
 *                     BY DIHEDRAL
 *
 *    ON OUTPUT GEO  = INTERNAL COORDINATES IN ANGSTROMS, RADIANS,
 *                     AND RADIANS
 *
 ***********************************************************************
       DO 10 I=2,NUMAT
          J=NA(I)
          K=NB(I)
          L=NC(I)
          IF(I.LT.3) GOTO 10
          II=I
          CALL BANGLE(XYZ,II,J,K,GEO(2,I))
          GEO(2,I)=GEO(2,I)*DEGREE
          IF(I.LT.4) GOTO 10
          CALL DIHED(XYZ,II,J,K,L,GEO(3,I))
          GEO(3,I)=GEO(3,I)*DEGREE
    10 GEO(1,I)= SQRT((XYZ(1,I)-XYZ(1,J))**2+
      1                   (XYZ(2,I)-XYZ(2,J))**2+
      2                   (XYZ(3,I)-XYZ(3,J))**2)
       GEO(1,1)=0.D0
       GEO(2,1)=0.D0
       GEO(3,1)=0.D0
       GEO(2,2)=0.D0
       GEO(3,2)=0.D0
       GEO(3,3)=0.D0
       RETURN
       END
       SUBROUTINE XYZINT(XYZ,NUMAT,NA,NB,NC,DEGREE,GEO)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION XYZ(3,*), NA(*), NB(*), NC(*), GEO(3,*)
 ***********************************************************************
 *
 * XYZINT WORKS OUT THE INTERNAL COORDINATES OF A MOLECULE.
 *        THE "RULES" FOR THE CONNECTIVITY ARE AS FOLLOWS:
 *        ATOM I IS DEFINED AS BEING AT A DISTANCE FROM THE NEAREST
 *        ATOM J, ATOM J ALREADY HAVING BEEN DEFINED.
 *        ATOM I MAKES AN ANGLE WITH ATOM J AND THE ATOM K, WHICH HAS
 *        ALREADY BEEN DEFINED, AND IS THE NEAREST ATOM TO J
 *        ATOM I MAKES A DIHEDRAL ANGLE WITH ATOMS J, K, AND L. L HAVING
 *        BEEN DEFINED AND IS THE NEAREST ATOM TO K
 *
 *        NOTE THAT GEO AND XYZ MUST NOT BE THE SAME IN THE CALL.
 *
 *   ON INPUT XYZ    = CARTESIAN ARRAY OF NUMAT ATOMS
 *            DEGREE = 1 IF ANGLES ARE TO BE IN RADIANS
 *            DEGREE = 57.29578 IF ANGLES ARE TO BE IN RADIANS
 *
 ***********************************************************************
       NAI1=0
       NAI2=0
       DO 20 I=1,NUMAT
          NA(I)=2
          NB(I)=3
          NC(I)=4
          IM1=I-1
          IF(IM1.EQ.0)GOTO 20
          SUM=100.D0
          DO 10 J=1,IM1
             R=(XYZ(1,I)-XYZ(1,J))**2+
      1          (XYZ(2,I)-XYZ(2,J))**2+
      2          (XYZ(3,I)-XYZ(3,J))**2
             IF(R.LT.SUM.AND.NA(J).NE.J.AND.NB(J).NE.J) THEN
                SUM=R
                K=J
             ENDIF
    10    CONTINUE
 C
 C   ATOM I IS NEAREST TO ATOM K
 C
          NA(I)=K
          IF(I.GT.2)NB(I)=NA(K)
          IF(I.GT.3)NC(I)=NB(K)
 C
 C   FIND ANY ATOM TO RELATE TO NA(I)
 C
    20 CONTINUE
       NA(1)=0
       NB(1)=0
       NC(1)=0
       NB(2)=0
       NC(2)=0
       NC(3)=0
       CALL XYZGEO(XYZ,NUMAT,NA,NB,NC,DEGREE,GEO)
       RETURN
       END
      SUBROUTINE GMETRY(GEO,COORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1NA(NUMATM),NB(NUMATM),NC(NUMATM)
     2       /EULER / TVEC(3,3), ID
      COMMON /REACTN/ STEP, GEOA(3,NUMATM), GEOVEC(3,NUMATM),COLCST
      COMMON /GEOOK/ IGEOOK
      COMMON /NUMCAL/ NUMCAL
      DIMENSION GEO(3,*),COORD(3,*)
      CHARACTER *15 NDIMEN(4)
      LOGICAL GEOOK
      SAVE ICALCN, NDIMEN, GEOOK
      DATA ICALCN/0/
      DATA NDIMEN/' MOLECULE     ',' POLYMER       ',
     1'LAYER STRUCTURE',' SOLID         '/
C***********************************************************************
C
C    GMETRY  COMPUTES COORDINATES FROM BOND-ANGLES AND LENGTHS.
C *** IT IS ADAPTED FROM THE PROGRAM WRITTEN BY M.J.S. DEWAR.
C
C     THREE SEPARATE OPTIONS EXIST WITHIN GMETRY. THESE ARE:
C    (A) IF NA(1) IS EQUAL TO 99 (IMPOSSIBLE UNDER NORMAL CIRCUMSTANCES)
C        THEN GEO IS ASSUMED TO BE IN CARTESIAN RATHER THAN INTERNAL
C        COORDINATES, AND COORD IS THEN SET EQUAL TO GEO.
C    (B) IF STEP IS NON-ZERO (THIS IS THE CASE WHEN "SADDLE" IS USED)
C        THEN GEO IS FIRST MODIFIED BY SHIFTING THE INTERNAL COORDINATES
C        ALONG A RADIUS FROM GEOA TO PLACE GEO AT ADISTANCESTEPFROMGEOA.
C    (C) NORMAL CONVERSION FROM INTERNAL TO CARTESIAN COORDINATESISDONE.
C
C  ON INPUT:
C         GEO    = ARRAY OF INTERNAL COORDINATES.
C         NATOMS = NUMBER OF ATOMS, INCLUDING DUMMIES.
C         NA     = ARRAY OF ATOM LABELS FOR BOND LENGTHS.
C
C  ON OUTPUT:
C         COORD  = ARRAY OF CARTESIAN COORDINATES
C
C***********************************************************************
C                                     OPTION (B)
      GEOOK=(IGEOOK.EQ.99)
      IF(ABS(STEP) .GT. 1.D-4) THEN
         SUM=0.D0
         DO 10 J=1,3
C$DOIT VBEST
            DO 10 I=1,NATOMS
               GEOVEC(J,I)=GEO(J,I)-GEOA(J,I)
   10    SUM=SUM+GEOVEC(J,I)**2
         SUM=SQRT(SUM)
         ERROR=(SUM-STEP)/SUM
      ELSE
         ERROR=0.D0
      ENDIF
      DO 20 J=1,3
C$DOIT VBEST
         DO 20 I=1,NATOMS
   20 GEO(J,I)=GEO(J,I)-ERROR*GEOVEC(J,I)
C                                     OPTION (A)
      IF(NA(1).EQ.99) THEN
         DO 30 I=1,3
C$DOIT VBEST
            DO 30 J=1,NATOMS
   30    COORD(I,J)=GEO(I,J)
         GOTO 100
      ENDIF
C                                     OPTION (C)
      COORD(1,1)=0.0D00
      COORD(2,1)=0.0D00
      COORD(3,1)=0.0D00
      COORD(1,2)=GEO(1,2)
      COORD(2,2)=0.0D00
      COORD(3,2)=0.0D00
      IF(NATOMS.EQ.2) GOTO 100
      CCOS=COS(GEO(2,3))
      IF(NA(3).EQ.1)THEN
         COORD(1,3)=COORD(1,1)+GEO(1,3)*CCOS
      ELSE
         COORD(1,3)=COORD(1,2)-GEO(1,3)*CCOS
      ENDIF
      COORD(2,3)=GEO(1,3)*SIN(GEO(2,3))
      COORD(3,3)=0.0D00
      DO 90 I=4,NATOMS
         COSA=COS(GEO(2,I))
         MB=NB(I)
         MC=NA(I)
         XB=COORD(1,MB)-COORD(1,MC)
         YB=COORD(2,MB)-COORD(2,MC)
         ZB=COORD(3,MB)-COORD(3,MC)
         RBC=XB*XB+YB*YB+ZB*ZB
         IF(RBC.LT.1.D-16)THEN
C
C     TWO ATOMS ARE COINCIDENT.  A FATAL ERROR.
C
            WRITE(6,'(A,I4,A,I4,A)')' ATOMS',MB,' AND',MC,' ARE COINCIDE
     1NT'
            WRITE(6,'(A)')' THIS IS A FATAL ERROR, RUN STOPPED IN GMETRY
     1'
            STOP
         ELSE
            RBC=1.0D00/SQRT(RBC)
         ENDIF
         MA=NC(I)
         XA=COORD(1,MA)-COORD(1,MC)
         YA=COORD(2,MA)-COORD(2,MC)
         ZA=COORD(3,MA)-COORD(3,MC)
C
C     ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.  IF XYB IS
C     TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
C
         XYB=SQRT(XB*XB+YB*YB)
         K=-1
         IF (XYB.GT.0.1D00) GO TO 40
         XPA=ZA
         ZA=-XA
         XA=XPA
         XPB=ZB
         ZB=-XB
         XB=XPB
         XYB=SQRT(XB*XB+YB*YB)
         K=+1
C
C     ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
C
   40    COSTH=XB/XYB
         SINTH=YB/XYB
         XPA=XA*COSTH+YA*SINTH
         YPA=YA*COSTH-XA*SINTH
         SINPH=ZB*RBC
         COSPH=SQRT(ABS(1.D00-SINPH*SINPH))
         ZQA=ZA*COSPH-XPA*SINPH
C
C     ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE.
C
         YZA=SQRT(YPA**2+ZQA**2)
         IF(YZA.LT.1.D-4)GOTO 60
         IF(YZA.LT.2.D-2 .AND. .NOT.GEOOK)THEN
            WRITE(6,'(//20X,'' CALCULATION ABANDONED AT THIS POINT'')')
            WRITE(6,'(//10X,'' THREE ATOMS BEING USED TO DEFINE THE'',/
     110X,'' COORDINATES OF A FOURTH ATOM, WHOSE BOND-ANGLE IS'')')
            WRITE(6,'(10X,'' NOT ZERO OR 180 DEGREEES, ARE '',
     1''IN AN ALMOST STRAIGHT'')')
            WRITE(6,'(10X,'' LINE.  THERE IS A HIGH PROBABILITY THAT THE
     1'',/10X,'' COORDINATES OF THE ATOM WILL BE INCORRECT.'')')
            WRITE(6,'(//20X,''THE FAULTY ATOM IS ATOM NUMBER'',I4)')I
            CALL GEOUT(1)
            WRITE(6,'(//20X,''CARTESIAN COORDINATES UP TO FAULTY ATOM'')
     1')
            WRITE(6,'(//5X,''I'',12X,''X'',12X,''Y'',12X,''Z'')')
            DO 50 J=1,I
   50       WRITE(6,'(I6,F16.5,2F13.5)')J,(COORD(K,J),K=1,3)
            WRITE(6,'(//6X,'' ATOMS'',I3,'','',I3,'', AND'',I3,
     1'' ARE WITHIN'',F7.4,'' ANGSTROMS OF A STRAIGHT LINE'')')
     2MC,MB,MA,YZA
            STOP
         ENDIF
         COSKH=YPA/YZA
         SINKH=ZQA/YZA
         GOTO 70
   60    CONTINUE
C
C   ANGLE TOO SMALL TO BE IMPORTANT
C
         COSKH=1.D0
         SINKH=0.D0
   70    CONTINUE
C
C     COORDINATES :-   A=(???,YZA,0),   B=(RBC,0,0),  C=(0,0,0)
C     NONE ARE NEGATIVE.
C     THE COORDINATES OF I ARE EVALUATED IN THE NEW FRAME.
C
         SINA=SIN(GEO(2,I))
         SIND=-SIN(GEO(3,I))
         COSD=COS(GEO(3,I))
         XD=GEO(1,I)*COSA
         YD=GEO(1,I)*SINA*COSD
         ZD=GEO(1,I)*SINA*SIND
C
C     TRANSFORM THE COORDINATES BACK TO THE ORIGINAL SYSTEM.
C
         YPD=YD*COSKH-ZD*SINKH
         ZPD=ZD*COSKH+YD*SINKH
         XPD=XD*COSPH-ZPD*SINPH
         ZQD=ZPD*COSPH+XD*SINPH
         XQD=XPD*COSTH-YPD*SINTH
         YQD=YPD*COSTH+XPD*SINTH
         IF (K.LT.1) GO TO 80
         XRD=-ZQD
         ZQD=XQD
         XQD=XRD
   80    COORD(1,I)=XQD+COORD(1,MC)
         COORD(2,I)=YQD+COORD(2,MC)
         COORD(3,I)=ZQD+COORD(3,MC)
   90 CONTINUE
C
C *** NOW REMOVE THE TRANSLATION VECTORS, IF ANY, FROM THE ARRAY COOR
C
  100 CONTINUE
      K=NATOMS
  110 IF(LABELS(K).NE.107) GOTO 120
      K=K-1
      GOTO 110
  120 K=K+1
      IF(K.GT.NATOMS) GOTO 170
C
C   SYSTEM IS A SOLID, OF DIMENSION NATOMS+1-K
C
      L=0
      DO 130 I=K,NATOMS
         L=L+1
         MC=NA(I)
         TVEC(1,L)=COORD(1,I)-COORD(1,MC)
         TVEC(2,L)=COORD(2,I)-COORD(2,MC)
         TVEC(3,L)=COORD(3,I)-COORD(3,MC)
  130 CONTINUE
      ID=L
      IF (ICALCN.NE.NUMCAL) THEN
         ICALCN=NUMCAL
         WRITE(6,140)NDIMEN(ID+1)
  140    FORMAT(/10X,'    THE SYSTEM IS A ',A15,/)
         IF(ID.EQ.0) GOTO 170
         WRITE(6,150)
         WRITE(6,160)(I,(TVEC(J,I),J=1,3),I=1,ID)
  150    FORMAT(/,'                UNIT CELL TRANSLATION VECTORS',/
     1/,'              X              Y              Z')
  160    FORMAT('    T',I1,' = ',F11.7,'    ',F11.7,'    ',F11.7)
      ENDIF
  170 CONTINUE
      J=0
      DO 190 I=1,NATOMS
         IF (LABELS(I).EQ.99.OR.LABELS(I).EQ.107) GO TO 190
         J=J+1
C$DOIT ASIS
         DO 180 K=1,3
  180    COORD(K,J)=COORD(K,I)
  190 CONTINUE
      RETURN
      END
*COMDECK SIZES                                                          
************************************************************************
*   THIS FILE CONTAINS ALL THE ARRAY SIZES FOR USE IN MOPAC.            
*                                                                       
*     THERE ARE ONLY 5 PARAMETERS THAT THE PROGRAMMER NEED SET:          
*     MAXHEV = MAXIMUM NUMBER OF HEAVY ATOMS (HEAVY: NON-HYDROGEN ATOMS)
*     MAXLIT = MAXIMUM NUMBER OF HYDROGEN ATOMS.                        
*     MAXTIM = DEFAULT TIME FOR A JOB. (SECONDS)                        
*     MAXDMP = DEFAULT TIME FOR AUTOMATIC RESTART FILE GENERATION (SECS)
*     ISYBYL = 1 IF MOPAC IS TO BE USED IN THE SYBYL PACKAGE, =0 OTHERWISE
*     SEE ALSO NMECI, NPULAY AND MESP AT THE END OF THIS FILE
*                                                                       
      PARAMETER (MAXHEV=60, MAXLIT=60)
      PARAMETER (MAXTIM=3600, MAXDMP=3600)
      PARAMETER (ISYBYL=0)
*                                                                       
************************************************************************
*                                                                       
*   THE FOLLOWING CODE DOES NOT NEED TO BE ALTERED BY THE PROGRAMMER    
*                                                                       
************************************************************************
*                                                                       
*    ALL OTHER PARAMETERS ARE DERIVED FUNCTIONS OF THESE TWO PARAMETERS 
*                                                                       
*      NAME                   DEFINITION                                
*     NUMATM         MAXIMUM NUMBER OF ATOMS ALLOWED.                   
*     MAXORB         MAXIMUM NUMBER OF ORBITALS ALLOWED.                
*     MAXPAR         MAXIMUM NUMBER OF PARAMETERS FOR OPTIMISATION.     
*     N2ELEC         MAXIMUM NUMBER OF TWO ELECTRON INTEGRALS ALLOWED.  
*     MPACK          AREA OF LOWER HALF TRIANGLE OF DENSITY MATRIX.     
*     MORB2          SQUARE OF THE MAXIMUM NUMBER OF ORBITALS ALLOWED.  
*     MAXHES         AREA OF HESSIAN MATRIX                             
*     MAXALL         LARGER THAN MAXORB OR MAXPAR.
************************************************************************
      PARAMETER (VERSON=7.01D0)
      PARAMETER (NUMATM=MAXHEV+MAXLIT)
      PARAMETER (MAXORB=4*MAXHEV+MAXLIT)
      PARAMETER (MAXPAR=3*NUMATM)
      PARAMETER (MAXBIG=MAXORB*MAXORB*2)
      PARAMETER (N2ELEC=(50*MAXHEV*(MAXHEV-1)+10*MAXHEV*MAXLIT
     +                     +(MAXLIT*(MAXLIT-1))/2))
      PARAMETER (MAXHES=(MAXPAR*(MAXPAR+1))/2,MORB2=MAXORB**2)
      PARAMETER (MPACK=(MAXORB*(MAXORB+1))/2)
      PARAMETER (MAXPR=6*MAXORB)
      PARAMETER (MAXALL=4*MAXHEV+3*MAXLIT)
      PARAMETER (NMECI=8, NPULAY=MPACK, MMCI=60)
      PARAMETER (MESP=50000)
      PARAMETER (LENABC=400)
      PARAMETER (LENAB2=LENABC*(LENABC+5))
      PARAMETER (NPPA = 1082, MAXNSS = 500)
      PARAMETER (MAXDEN=10*MAXHEV+MAXLIT)
*
      PARAMETER (IGREEN=600000)
************************************************************************
*DECK MOPAC 
