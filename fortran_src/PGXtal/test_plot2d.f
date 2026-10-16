CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      PROGRAM TEST_PLOT2D
C     -------------------
C
C  A little program to test the "PLOT" subroutine for plotting a 
C  2-dimensional array with a simple rectangular mapping.
C
      REAL X(65),Y(65),Z(65,65),W(65)
C
      NG=40
      NDAT=2*(NG/2)+1
      N2=NDAT/2+1
      SIG2=15.0/FLOAT(N2**2)
      NX=NDAT
      NY=NDAT
      CALL FILARY(NX,NY,N2,X,Y,Z,SIG2,NDAT)
      CALL PLOT(X,NX,Y,NY,Z,NDAT,NDAT,W,1.5,2,'X-axis','Y-axis',
     *         'An example of 2-D plotting with PLOT')
      END
C
      SUBROUTINE FILARY(NX,NY,N2,X,Y,Z,SIG2,NDAT)
C     -------------------------------------------
C
      REAL X(*),Y(*),Z(NDAT,NDAT)
C
      RAD0=2.0/SQRT(SIG2)
      X0=-RAD0/2.0
      Y0=-X0/2.0
      SIG2B=3.0*SIG2
      DO 20 J=1,NY
        Y(J)=FLOAT(J-N2)
        YY=Y(J)**2
        DO 10 I=1,NX
          X(I)=FLOAT(I-N2)
          RAD2=(X(I)-X0)**2+(Y(J)-Y0)**2
          Z(I,J)=EXP(-SIG2*(X(I)**2+YY))
     *          +0.2*EXP(-SIG2B*(SQRT(RAD2)-RAD0)**2)
  10    CONTINUE
  20  CONTINUE
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
