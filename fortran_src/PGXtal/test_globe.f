CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      PROGRAM TEST_GLOBE
C     ------------------
C
C  A little program to test the "GLOBE" subroutine for plotting a 
C  2-dimensional array onto the surface of a sphere.
C
      REAL    Z(91,91)
      LOGICAL ANOTHR
      DATA    NTHT,NPHI /91,91/
C
      CALL FILARY(Z,THTMIN,PHIMIN,NTHT,NPHI,DTHT,DPHI)
      CALL GLOBE(Z,THTMIN,PHIMIN,NTHT,NPHI,DTHT,DPHI,1.5,2,' ',' ',
     *           'An example of 2-D plotting with GLOBE')
      END
C
      SUBROUTINE FILARY(Z,THTMIN,PHIMIN,NTHT,NPHI,DTHT,DPHI)
C     ------------------------------------------------------
C
      REAL Z(NPHI,NTHT)
C
      THTMIN=0.0
      THTMAX=90.0
      DTHT=(THTMAX-THTMIN)/FLOAT(NTHT-1)
      PHIMIN=0.0
      PHIMAX=90.0
      DPHI=(PHIMAX-PHIMIN)/FLOAT(NPHI-1)
      NY2=NTHT/2+1
      NX2=NPHI/2+1
      N2=(NTHT+NPHI)/4+1
      SIG2=15.0/FLOAT(N2**2)
      RAD0=2.0/SQRT(SIG2)
      X0=-RAD0/2.0
      Y0=-X0/2.0
      SIG2B=3.0*SIG2
      DO 20 J=1,NTHT
        YJ=FLOAT(J-NY2)
        YY=YJ**2
        DO 10 I=1,NPHI
          XI=FLOAT(I-NX2)
          RAD2=(XI-X0)**2+(YJ-Y0)**2
          Z(I,J)=EXP(-SIG2*(XI**2+YY))
     *          +0.2*EXP(-SIG2B*(SQRT(RAD2)-RAD0)**2)+0.02
  10    CONTINUE
  20  CONTINUE
      Z(1,1)=0.0
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
