      SUBROUTINE SBSURF(EYE,LATICE,DENS,N1,N2,N3,DSURF,IC1,IC2,LIGHT,
     *                  LSHINE)
C     ---------------------------------------------------------------
C
      REAL             EYE(*),LATICE(3,*),DENS(0:N1,0:N2,0:N3),LIGHT(*)
      LOGICAL          LSHINE
C
      REAL*8           XLNORM,PNEYE(3),DSMAL2
      REAL             BAS(3,3),MTRX
      REAL             XYZ(3),DXYZ(3,3),FRCXYZ(12),DDXYZ(3,12,2)
      REAL             DLOCAL(8),VERT(3,12),GRDSCL(3)
      LOGICAL          LPS,LCOLOR,LEMPTY
      INTEGER          IVERT(8)
      COMMON  /SRFCOM/ GRDCUB(3,8),MTRX(3,3),ORIG(3),XL2,COL0,COLSCL
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      This subroutine plots an iso-surface through a unit-cell of
C    density. All (x,y,z) values are taken to be given in world 
C    coordinates. The z-component of the eye-poisition should be 
C    positive and that of all the lattice-vertices should be negative; 
C    the viewing-screen is fixed at z=0.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    EYE      R*4    I       3      (x,y,z) coordinate of eye-position.
C    LATICE   R*4    I     3 x 4    (x,y,z) coordinates of the origin
C                                   and the a, b & C lattice-vertices.
C    DENS     R*4    I     (N1+1)   The density at regular points within
C                        x (N2+1)   the unit cell, wrapped around so
C                        x (N3+1)   that DENS(0,J,K)=DENS(N1,J,K) etc..
C    N1,N2,N3 I*4    I       -      The dimensions of the unit-cell grid.
C    DSURF    R*4    I       -      Density for the iso-surface.
C    IC1,IC2  I*4    I       -      Lowest & highest colour-index to be
C                                   used for the shading.
C    LIGHT    R*4    I       3      (x,y,z) direction of flood-light.
C    LSHINE   L*1    I       -      Shiny surface if TRUE, else diffuse.
C
C Globals 
C    SFTBUF
C    SRFCOM
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     SBSRF0     A quick check in case there are no iso-surafces.
C     SBSRF1     Anlayses a 2-d box in the surface of the unit-cell.
C     SBSRF2     Paints a 2-d box in the surface of the unit-cell.
C     SBSRF3     Analyses a 3-d box within the unit-cell.
C     SBSRF4     Initialises the gradients for a 3-d box.
C     SBSRF5     Breaks up the iso-surface in a 3-d box into triangles.
C     SBSRF6     Paints a triangular patch of an iso-surface.
C
C History
C   D. S. Sivia       3 May 1995  Initial release.
C   D. S. Sivia       7 Jul 1995  Fixed bug in determinant calculation.
C   D. S. Sivia      20 Oct 1995  Speeded up computations slightly.
C   D. S. Sivia      13 Dec 1995  A bit more tinkering for speed.
C   D. S. Sivia      14 Jun 1996  Completely new algorithm!
C   D. S. Sivia      24 Oct 1997 "Safe-guarded" some rounding errors.
C-----------------------------------------------------------------------
C
C Carry out some initial checks.
C
      SMALL=1.0E-10
      SMALL2=SMALL**2
      DSMAL2=DBLE(SMALL2)
      IF (EYE(3).LE.SMALL) RETURN
      IF (N1.LT.1 .OR. N2.LT.1 .OR. N3.LT.1) RETURN
      XL2=LIGHT(1)**2+LIGHT(2)**2+LIGHT(3)**2
      IF (XL2.LT.SMALL) RETURN
      IF (LATICE(3,1).GE.0.0) RETURN
      ZFAR=LATICE(3,2)+LATICE(3,3)+LATICE(3,4)-2.0*LATICE(3,1)
      IF (ZFAR.GE.0.0) RETURN
      DO 10 J=1,3
        IF (LATICE(3,J+1).GE.0.0) RETURN
        BAS(1,J)=LATICE(1,J+1)-LATICE(1,1)
        BAS(2,J)=LATICE(2,J+1)-LATICE(2,1)
        BAS(3,J)=LATICE(3,J+1)-LATICE(3,1)
        IF ((ZFAR-BAS(3,J)).GE.0.0) RETURN
        BAS2J=BAS(1,J)**2+BAS(2,J)**2+BAS(3,J)**2
        IF (BAS2J.LT.SMALL2) RETURN
  10  CONTINUE
      NTOT=(N1+1)*(N2+1)*(N3+1)
      CALL SBSRF0(DENS,NTOT,DSURF,LEMPTY)
      IF (LEMPTY) RETURN
C
C Set up matrix for real-space to lattice-index transformation.
C
      XN1=0.99999*FLOAT(N1)
      XN2=0.99999*FLOAT(N2)
      XN3=0.99999*FLOAT(N3)
      DET=BAS(1,1)*BAS(2,2)*BAS(3,3)+BAS(1,2)*BAS(2,3)*BAS(3,1)
     *   +BAS(1,3)*BAS(2,1)*BAS(3,2)-BAS(3,1)*BAS(2,2)*BAS(1,3)
     *   -BAS(3,2)*BAS(2,3)*BAS(1,1)-BAS(3,3)*BAS(2,1)*BAS(1,2)
      IF (ABS(DET).LT.SMALL2) RETURN
      DETNRM=1.0/DET
      MTRX(1,1)=XN1*DETNRM*(BAS(2,2)*BAS(3,3)-BAS(2,3)*BAS(3,2))
      MTRX(1,2)=XN2*DETNRM*(BAS(2,3)*BAS(3,1)-BAS(2,1)*BAS(3,3))
      MTRX(1,3)=XN3*DETNRM*(BAS(2,1)*BAS(3,2)-BAS(2,2)*BAS(3,1))
      MTRX(2,1)=XN1*DETNRM*(BAS(3,2)*BAS(1,3)-BAS(3,3)*BAS(1,2))
      MTRX(2,2)=XN2*DETNRM*(BAS(3,3)*BAS(1,1)-BAS(3,1)*BAS(1,3))
      MTRX(2,3)=XN3*DETNRM*(BAS(3,1)*BAS(1,2)-BAS(3,2)*BAS(1,1))
      MTRX(3,1)=XN1*DETNRM*(BAS(1,2)*BAS(2,3)-BAS(1,3)*BAS(2,2))
      MTRX(3,2)=XN2*DETNRM*(BAS(1,3)*BAS(2,1)-BAS(1,1)*BAS(2,3))
      MTRX(3,3)=XN3*DETNRM*(BAS(1,1)*BAS(2,2)-BAS(1,2)*BAS(2,1))
      CALL SBRCOP(LATICE,ORIG,3)
C
C Some general initialisations.
C
      DDSURF=MAX(ABS(DSURF),SMALL)
      IF (DSURF.LT.0.0) DDSURF=-DDSURF
      GRDSCL(1)=-0.5/(DDSURF*FLOAT(N1))
      GRDSCL(2)=-0.5/(DDSURF*FLOAT(N2))
      GRDSCL(3)=-0.5/(DDSURF*FLOAT(N3))
      COL0=FLOAT(IC1)
      COLSCL=FLOAT(IC2-IC1)
      DO 30 I=1,3
        DXYZ(I,1)=BAS(I,1)/FLOAT(N1)
        DXYZ(I,2)=BAS(I,2)/FLOAT(N2)
        DXYZ(I,3)=BAS(I,3)/FLOAT(N3)
        DDXYZ(I,1,1)=0.0
        DDXYZ(I,1,2)=DXYZ(I,1)
        DDXYZ(I,2,1)=DDXYZ(I,1,1)+DDXYZ(I,1,2)
        DDXYZ(I,2,2)=DXYZ(I,2)
        DDXYZ(I,3,1)=DDXYZ(I,2,1)+DDXYZ(I,2,2)
        DDXYZ(I,3,2)=-DXYZ(I,1)
        DDXYZ(I,4,1)=DDXYZ(I,3,1)+DDXYZ(I,3,2)
        DDXYZ(I,4,2)=-DXYZ(I,2)
        DO 20 J=1,4
          DDXYZ(I,J+4,1)=DDXYZ(I,J,1)
          DDXYZ(I,J+4,2)=DXYZ(I,3)
          DDXYZ(I,J+8,1)=DDXYZ(I,J,1)+DXYZ(I,3)
          DDXYZ(I,J+8,2)=DDXYZ(I,J,2)
  20    CONTINUE
  30  CONTINUE
C
C First paint the edges of the lattice.
C
      DO 60 IFACE=1,3
        I=IFACE
        J=MOD(IFACE,3)+1
        K=MOD(J,3)+1
        IF (IFACE.EQ.1) THEN
          IN=N1
          JN=N2
        ELSEIF (IFACE.EQ.2) THEN
          IN=N2
          JN=N3
        ELSE
          IN=N3
          JN=N1
        ENDIF
        KK=0
        XN=BAS(2,J)*BAS(3,I)-BAS(2,I)*BAS(3,J)
        YN=BAS(3,J)*BAS(1,I)-BAS(3,I)*BAS(1,J)
        ZN=BAS(1,J)*BAS(2,I)-BAS(1,I)*BAS(2,J)
        DNRM=SQRT(XN**2+YN**2+ZN**2+SMALL2)
        PNEYE(1)=DBLE(EYE(1)-0.5*(LATICE(1,I+1)+LATICE(1,J+1)))
        PNEYE(2)=DBLE(EYE(2)-0.5*(LATICE(2,I+1)+LATICE(2,J+1)))
        PNEYE(3)=DBLE(EYE(3)-0.5*(LATICE(3,I+1)+LATICE(3,J+1)))
        DEYE=SNGL(DSQRT(PNEYE(1)**2+PNEYE(2)**2+PNEYE(3)**2+DSMAL2))
        XLNORM=DBLE(XN)*PNEYE(1)+DBLE(YN)*PNEYE(2)+DBLE(ZN)*PNEYE(3)
        COSSEE=SNGL(XLNORM)/(DEYE*DNRM)
        IF (COSSEE.LT.0.001) THEN
          KK=N3
          IF (IFACE.EQ.2) KK=N1
          IF (IFACE.EQ.3) KK=N2
          PNEYE(1)=PNEYE(1)+DBLE(BAS(1,K))
          PNEYE(2)=PNEYE(2)+DBLE(BAS(2,K))
          PNEYE(3)=PNEYE(3)+DBLE(BAS(3,K))
          DEYE=SNGL(DSQRT(PNEYE(1)**2+PNEYE(2)**2+PNEYE(3)**2+DSMAL2))
          XLNORM=DBLE(XN)*PNEYE(1)+DBLE(YN)*PNEYE(2)+DBLE(ZN)*PNEYE(3)
          COSSEE=-SNGL(XLNORM)/(DEYE*DNRM)
        ENDIF
        IF (COSSEE.GT.0.001) THEN
          XYZ1=FLOAT(KK)*DXYZ(1,K)+LATICE(1,1)
          XYZ2=FLOAT(KK)*DXYZ(2,K)+LATICE(2,1)
          XYZ3=FLOAT(KK)*DXYZ(3,K)+LATICE(3,1)
          DO 50 J1=1,JN
            J0=J1-1
            DO 40 I1=1,IN
              I0=I1-1
              IF (IFACE.EQ.1) THEN
                DLOCAL(1)=DENS(I0,J0,KK)-DSURF
                DLOCAL(2)=DENS(I1,J0,KK)-DSURF
                DLOCAL(3)=DENS(I1,J1,KK)-DSURF
                DLOCAL(4)=DENS(I0,J1,KK)-DSURF
              ELSEIF (IFACE.EQ.2) THEN
                DLOCAL(1)=DENS(KK,I0,J0)-DSURF
                DLOCAL(2)=DENS(KK,I1,J0)-DSURF
                DLOCAL(3)=DENS(KK,I1,J1)-DSURF
                DLOCAL(4)=DENS(KK,I0,J1)-DSURF
              ELSE
                DLOCAL(1)=DENS(J0,KK,I0)-DSURF
                DLOCAL(2)=DENS(J0,KK,I1)-DSURF
                DLOCAL(3)=DENS(J1,KK,I1)-DSURF
                DLOCAL(4)=DENS(J1,KK,I0)-DSURF
              ENDIF
              CALL SBSRF1(DLOCAL,IBSIDE,FRCXYZ)
              IF (IBSIDE.NE.0) THEN
                XYZ(1)=XYZ1+DXYZ(1,I)*FLOAT(I0)
                XYZ(2)=XYZ2+DXYZ(2,I)*FLOAT(I0)
                XYZ(3)=XYZ3+DXYZ(3,I)*FLOAT(I0)
                CALL SBSRF2(XYZ,DXYZ(1,I),DXYZ(1,J),IBSIDE,FRCXYZ,VERT,
     *                      EYE,LIGHT,LSHINE)
              ENDIF
  40        CONTINUE
            XYZ1=XYZ1+DXYZ(1,J)
            XYZ2=XYZ2+DXYZ(2,J)
            XYZ3=XYZ3+DXYZ(3,J)
  50      CONTINUE
        ENDIF
  60  CONTINUE
C
C Step through each "cube" in the lattice, and paint any isosurfaces
C found therein.
C
      X00K=LATICE(1,1)
      Y00K=LATICE(2,1)
      Z00K=LATICE(3,1)
      DO 90 K1=1,N3
        K0=K1-1
        DO 80 J1=1,N2
          J0=J1-1
          X0JK=X00K+DXYZ(1,2)*FLOAT(J0)
          Y0JK=Y00K+DXYZ(2,2)*FLOAT(J0)
          Z0JK=Z00K+DXYZ(3,2)*FLOAT(J0)
          DO 70 I1=1,N1
            I0=I1-1
            DLOCAL(1)=DENS(I0,J0,K0)-DSURF
            DLOCAL(2)=DENS(I1,J0,K0)-DSURF
            DLOCAL(3)=DENS(I1,J1,K0)-DSURF
            DLOCAL(4)=DENS(I0,J1,K0)-DSURF
            DLOCAL(5)=DENS(I0,J0,K1)-DSURF
            DLOCAL(6)=DENS(I1,J0,K1)-DSURF
            DLOCAL(7)=DENS(I1,J1,K1)-DSURF
            DLOCAL(8)=DENS(I0,J1,K1)-DSURF
            CALL SBSRF3(DLOCAL,IVERT,FRCXYZ,ISUMV,ISUMF)
            IF (ISUMV.NE.0) THEN
              XYZ(1)=X0JK+DXYZ(1,1)*FLOAT(I0)
              XYZ(2)=Y0JK+DXYZ(2,1)*FLOAT(I0)
              XYZ(3)=Z0JK+DXYZ(3,1)*FLOAT(I0)
              CALL SBSRF4(DENS,N1,N2,N3,I0,J0,K0,GRDSCL,BAS,GRDCUB)
              CALL SBSRF5(XYZ,DDXYZ,ISUMV,ISUMF,IVERT,FRCXYZ,VERT,EYE,
     *                    LIGHT,LSHINE)
            ENDIF
  70      CONTINUE
  80    CONTINUE
        X00K=X00K+DXYZ(1,3)
        Y00K=Y00K+DXYZ(2,3)
        Z00K=Z00K+DXYZ(3,3)
  90  CONTINUE
      END
C
      SUBROUTINE SBSRF0(DENS,NTOT,DSURF,LEMPTY)
C     -----------------------------------------
C
      REAL    DENS(*)
      LOGICAL LEMPTY
C
      LEMPTY=.TRUE.
      DO 10 I=1,NTOT
        IF (DENS(I).GT.DSURF) THEN
          LEMPTY=.FALSE.
          RETURN
        ENDIF
  10  CONTINUE
      END
C
      SUBROUTINE SBSRF1(D,IB,DF)
C     --------------------------
C
      REAL D(*),DF(*)
      DATA SMALL /1.0E-20/
C
      IB=0
      IF (D(1).GE.0.0) IB=1
      IF (D(2).GE.0.0) IB=IB+2
      IF (D(3).GE.0.0) IB=IB+4
      IF (D(4).GE.0.0) IB=IB+8
      IF (IB.EQ.0 .OR. IB.EQ.15) RETURN
      DO 10 I=1,4
        J=1+MOD(I,4)
        IF (D(I)*D(J).LT.-SMALL) THEN
          DI=ABS(D(I))
          DF(I)=DI/(DI+ABS(D(J)))
        ENDIF
  10  CONTINUE
      END
C
      SUBROUTINE SBSRF2(XYZ,D1,D2,IB,FRC,VERT,EYE,LIGHT,LSHINE)
C     ---------------------------------------------------------
C
      REAL    XYZ(*),D1(*),D2(*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      LOGICAL LSHINE
C
      IF (IB.EQ.15) THEN
        DO 10 I=1,3
          VERT(I,1)=XYZ(I)
          VERT(I,2)=XYZ(I)+D1(I)
          VERT(I,3)=VERT(I,2)+D2(I)
          VERT(I,4)=XYZ(I)+D2(I)
  10    CONTINUE
        CALL SBSRF6(EYE,4,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.1) THEN
        DO 20 I=1,3
          VERT(I,1)=XYZ(I)
          VERT(I,2)=VERT(I,1)+FRC(1)*D1(I)
          VERT(I,3)=XYZ(I)+(1.0-FRC(4))*D2(I)
  20    CONTINUE
        CALL SBSRF6(EYE,3,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.2) THEN
        DO 30 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)
          VERT(I,2)=VERT(I,1)+FRC(2)*D2(I)
          VERT(I,3)=XYZ(I)+FRC(1)*D1(I)
  30    CONTINUE
        CALL SBSRF6(EYE,3,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.4) THEN
        DO 40 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)+D2(I)
          VERT(I,2)=VERT(I,1)-FRC(3)*D1(I)
          VERT(I,3)=VERT(I,1)-(1.0-FRC(2))*D2(I)
  40    CONTINUE
        CALL SBSRF6(EYE,3,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.8) THEN
        DO 50 I=1,3
          VERT(I,1)=XYZ(I)+D2(I)
          VERT(I,2)=VERT(I,1)-FRC(4)*D2(I)
          VERT(I,3)=VERT(I,1)+(1.0-FRC(3))*D1(I)
  50    CONTINUE
        CALL SBSRF6(EYE,3,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.7) THEN
        DO 60 I=1,3
          VERT(I,1)=XYZ(I)
          VERT(I,2)=VERT(I,1)+D1(I)
          VERT(I,3)=VERT(I,2)+D2(I)
          VERT(I,4)=VERT(I,3)-FRC(3)*D1(I)
          VERT(I,5)=XYZ(I)+(1.0-FRC(4))*D2(I)
  60    CONTINUE
        CALL SBSRF6(EYE,5,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.14) THEN
        DO 70 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)
          VERT(I,2)=VERT(I,1)+D2(I)
          VERT(I,3)=XYZ(I)+D2(I)
          VERT(I,4)=XYZ(I)+(1.0-FRC(4))*D2(I)
          VERT(I,5)=XYZ(I)+FRC(1)*D1(I)
  70    CONTINUE
        CALL SBSRF6(EYE,5,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.13) THEN
        DO 80 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)+D2(I)
          VERT(I,2)=XYZ(I)+D2(I)
          VERT(I,3)=XYZ(I)
          VERT(I,4)=XYZ(I)+FRC(1)*D1(I)
          VERT(I,5)=XYZ(I)+D1(I)+FRC(2)*D2(I)
  80    CONTINUE
        CALL SBSRF6(EYE,5,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.11) THEN
        DO 90 I=1,3
          VERT(I,1)=XYZ(I)+D2(I)
          VERT(I,2)=XYZ(I)
          VERT(I,3)=XYZ(I)+D1(I)
          VERT(I,4)=VERT(I,3)+FRC(2)*D2(I)
          VERT(I,5)=VERT(I,1)+(1.0-FRC(3))*D1(I)
  90    CONTINUE
        CALL SBSRF6(EYE,5,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.3) THEN
        DO 100 I=1,3
          VERT(I,1)=XYZ(I)
          VERT(I,2)=XYZ(I)+D1(I)
          VERT(I,3)=VERT(I,2)+FRC(2)*D2(I)
          VERT(I,4)=XYZ(I)+(1.0-FRC(4))*D2(I)
 100    CONTINUE
        CALL SBSRF6(EYE,4,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.6) THEN
        DO 110 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)
          VERT(I,2)=VERT(I,1)+D2(I)
          VERT(I,3)=VERT(I,2)-FRC(3)*D1(I)
          VERT(I,4)=XYZ(I)+FRC(1)*D1(I)
 110    CONTINUE
        CALL SBSRF6(EYE,4,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.12) THEN
        DO 120 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)+D2(I)
          VERT(I,2)=XYZ(I)+D2(I)
          VERT(I,3)=VERT(I,2)-FRC(4)*D2(I)
          VERT(I,4)=XYZ(I)+D1(I)+FRC(2)*D2(I)
 120    CONTINUE
        CALL SBSRF6(EYE,4,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.9) THEN
        DO 130 I=1,3
          VERT(I,1)=XYZ(I)+D2(I)
          VERT(I,2)=XYZ(I)
          VERT(I,3)=XYZ(I)+FRC(1)*D1(I)
          VERT(I,4)=VERT(I,1)+(1.0-FRC(3))*D1(I)
 130    CONTINUE
        CALL SBSRF6(EYE,4,VERT,LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.5) THEN
        DO 140 I=1,3
          VERT(I,1)=XYZ(I)
          VERT(I,2)=VERT(I,1)+FRC(1)*D1(I)
          VERT(I,3)=XYZ(I)+(1.0-FRC(4))*D2(I)
          VERT(I,4)=XYZ(I)+D1(I)+D2(I)
          VERT(I,5)=VERT(I,4)-FRC(3)*D1(I)
          VERT(I,6)=VERT(I,4)-(1.0-FRC(2))*D2(I)
 140    CONTINUE
        CALL SBSRF6(EYE,3,VERT,LSHINE,LIGHT,0) 
        CALL SBSRF6(EYE,3,VERT(1,4),LSHINE,LIGHT,0) 
      ELSEIF (IB.EQ.10) THEN
        DO 150 I=1,3
          VERT(I,1)=XYZ(I)+D1(I)
          VERT(I,2)=VERT(I,1)+FRC(2)*D2(I)
          VERT(I,3)=XYZ(I)+FRC(1)*D1(I)
          VERT(I,4)=XYZ(I)+D2(I)
          VERT(I,5)=VERT(I,4)-FRC(4)*D2(I)
          VERT(I,6)=VERT(I,4)+(1.0-FRC(3))*D1(I)
 150    CONTINUE
        CALL SBSRF6(EYE,3,VERT,LSHINE,LIGHT,0) 
        CALL SBSRF6(EYE,3,VERT(1,4),LSHINE,LIGHT,0) 
      ENDIF
      END
C
      SUBROUTINE SBSRF3(D,IVERT,DF,ISUMV,ISUMF)
C     -----------------------------------------
C
      REAL    D(*),DF(*)
      INTEGER IVERT(*),IC(8)
      DATA    SMALL /1.0E-20/
C
      ISUMV=0
      DO 10 I=1,8
        IF (D(I).LT.0.0) THEN
          IC(I)=0
        ELSE
          IC(I)=1
          ISUMV=ISUMV+1
        ENDIF
  10  CONTINUE
      IF (ISUMV.EQ.0 .OR. ISUMV.EQ.8) THEN
        ISUMV=0
        RETURN
      ENDIF
      IF (ISUMV.GT.4) THEN
        ISUMV=8-ISUMV
        DO 20 I=1,8
  20      IC(I)=MOD(IC(I)+1,2)
      ENDIF
      J=0
      DO 30 I=1,8
        IF (IC(I).EQ.1) THEN
          J=J+1
          IVERT(J)=I
        ENDIF
  30  CONTINUE
      ISUMF=0
      DO 40 I=1,4
        J=1+MOD(I,4)
        IF (D(I)*D(J).LT.-SMALL) THEN
          DI=ABS(D(I))
          DF(I)=DI/(DI+ABS(D(J)))
          ISUMF=ISUMF+1
        ENDIF
        K=I+4
        IF (D(I)*D(K).LT.-SMALL) THEN
          DI=ABS(D(I))
          DF(K)=DI/(DI+ABS(D(K)))
          ISUMF=ISUMF+1
        ENDIF
        L=J+4
        IF (D(K)*D(L).LT.-SMALL) THEN
          DK=ABS(D(K))
          DF(I+8)=DK/(DK+ABS(D(L)))
          ISUMF=ISUMF+1
        ENDIF
  40  CONTINUE
      END
C
      SUBROUTINE SBSRF4(DENS,N1,N2,N3,I0,J0,K0,GRDSCL,BAS,GRD)
C     --------------------------------------------------------
C
      REAL DENS(0:N1,0:N2,0:N3),GRDSCL(*),BAS(3,*),GRD(3,*)
      REAL G(3,8)
C
      IM=I0-1
      IF (IM.LT.0) IM=N1
      JM=J0-1
      IF (JM.LT.0) JM=N2
      KM=K0-1
      IF (KM.LT.0) KM=N3
      I1=I0+1
      J1=J0+1
      K1=K0+1
      IP=I1+1
      IF (IP.GT.N1) IP=0
      JP=J1+1
      IF (JP.GT.N2) JP=0
      KP=K1+1
      IF (KP.GT.N3) KP=0
      G(1,1)=GRDSCL(1)*(DENS(I1,J0,K0)-DENS(IM,J0,K0))
      G(2,1)=GRDSCL(2)*(DENS(I0,J1,K0)-DENS(I0,JM,K0))
      G(3,1)=GRDSCL(3)*(DENS(I0,J0,K1)-DENS(I0,J0,KM))
      G(1,2)=GRDSCL(1)*(DENS(IP,J0,K0)-DENS(I0,J0,K0))
      G(2,2)=GRDSCL(2)*(DENS(I1,J1,K0)-DENS(I1,JM,K0))
      G(3,2)=GRDSCL(3)*(DENS(I1,J0,K1)-DENS(I1,J0,KM))
      G(1,3)=GRDSCL(1)*(DENS(IP,J1,K0)-DENS(I0,J1,K0))
      G(2,3)=GRDSCL(2)*(DENS(I1,JP,K0)-DENS(I1,J0,K0))
      G(3,3)=GRDSCL(3)*(DENS(I1,J1,K1)-DENS(I1,J1,KM))
      G(1,4)=GRDSCL(1)*(DENS(I1,J1,K0)-DENS(IM,J1,K0))
      G(2,4)=GRDSCL(2)*(DENS(I0,JP,K0)-DENS(I0,J0,K0))
      G(3,4)=GRDSCL(3)*(DENS(I0,J1,K1)-DENS(I0,J1,KM))
      G(1,5)=GRDSCL(1)*(DENS(I1,J0,K1)-DENS(IM,J0,K1))
      G(2,5)=GRDSCL(2)*(DENS(I0,J1,K1)-DENS(I0,JM,K1))
      G(3,5)=GRDSCL(3)*(DENS(I0,J0,KP)-DENS(I0,J0,K0))
      G(1,6)=GRDSCL(1)*(DENS(IP,J0,K1)-DENS(I0,J0,K1))
      G(2,6)=GRDSCL(2)*(DENS(I1,J1,K1)-DENS(I1,JM,K1))
      G(3,6)=GRDSCL(3)*(DENS(I1,J0,KP)-DENS(I1,J0,K0))
      G(1,7)=GRDSCL(1)*(DENS(IP,J1,K1)-DENS(I0,J1,K1))
      G(2,7)=GRDSCL(2)*(DENS(I1,JP,K1)-DENS(I1,J0,K1))
      G(3,7)=GRDSCL(3)*(DENS(I1,J1,KP)-DENS(I1,J1,K0))
      G(1,8)=GRDSCL(1)*(DENS(I1,J1,K1)-DENS(IM,J1,K1))
      G(2,8)=GRDSCL(2)*(DENS(I0,JP,K1)-DENS(I0,J0,K1))
      G(3,8)=GRDSCL(3)*(DENS(I0,J1,KP)-DENS(I0,J1,K0))
      DO 20 J=1,8
        DO 10 I=1,3
  10      GRD(I,J)=G(1,J)*BAS(I,1)+G(2,J)*BAS(I,2)+G(3,J)*BAS(I,3)
  20  CONTINUE
      END
C
      SUBROUTINE SBSRF5(XYZ,DXYZ,ISV,ISF,IV,FRC,VERT,EYE,LIGHT,LSHINE)
C     ----------------------------------------------------------------
C
      REAL    XYZ(*),DXYZ(3,12,*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      INTEGER IV(*),IV4MAP(12)
      LOGICAL LSHINE
      DATA    IV4MAP /12,8,4,3,11,7,6,2,10,9,5,1/
C
      IF (ISV.EQ.1) THEN
        CALL SBSF5A(XYZ,DXYZ,FRC,VERT,IV(1),EYE,LSHINE,LIGHT)
      ELSEIF (ISV.EQ.2) THEN
        IF (ISF.EQ.6) THEN
          CALL SBSF5A(XYZ,DXYZ,FRC,VERT,IV(1),EYE,LSHINE,LIGHT)
          CALL SBSF5A(XYZ,DXYZ,FRC,VERT,IV(2),EYE,LSHINE,LIGHT)
        ELSE
          IJDIF=IV(2)-IV(1)
          IF (IV(1).LE.4) THEN
            IF (IV(2).LE.4) THEN
              K2=IV(1)
              IF (IJDIF.EQ.3) K2=IV(2)
            ELSE
              K2=IV(2)
            ENDIF
          ELSE
            K2=IV(1)+4
            IF (IJDIF.EQ.3) K2=IV(2)+4
          ENDIF
          CALL SBSF5B(XYZ,DXYZ,FRC,VERT,K2,EYE,LSHINE,LIGHT,1)
        ENDIF
      ELSEIF (ISV.EQ.3) THEN
        IF (ISF.EQ.9) THEN
          DO 10 I=1,3
  10        CALL SBSF5A(XYZ,DXYZ,FRC,VERT,IV(I),EYE,LSHINE,LIGHT)
        ELSEIF (ISF.EQ.6) THEN
          DO 20 I1=1,3
            I2=1+MOD(I1,3)
            I=MIN(I1,I2)
            J=MAX(I1,I2)
            K2=0
            IJDIF=IV(J)-IV(I)
            IF (IV(I).LE.4) THEN
              IF (IV(J).LE.4) THEN
                IF (IJDIF.EQ.1) THEN
                  K2=IV(I)
                ELSEIF (IJDIF.EQ.3) THEN
                  K2=IV(J)
                ENDIF
              ELSE
                IF (IJDIF.EQ.4) K2=IV(J)
              ENDIF
            ELSE
              IF (IJDIF.EQ.1) THEN
                K2=IV(I)+4
              ELSEIF (IJDIF.EQ.3) THEN
                K2=IV(J)+4
              ENDIF
            ENDIF
            IF (K2.GT.0) GOTO 1
  20      CONTINUE
   1      CALL SBSF5B(XYZ,DXYZ,FRC,VERT,K2,EYE,LSHINE,LIGHT,1)
        ELSE
          K3=IV(1)+IV(2)+IV(3)-5+2*(IV(1)/5+2*(IV(2)/5)+IV(3)/5)
          CALL SBSF5C(XYZ,DXYZ,FRC,VERT,K3,EYE,LSHINE,LIGHT)
        ENDIF
      ELSE
        IF (ISF.EQ.12) THEN
          DO 30 I=1,4
  30        CALL SBSF5A(XYZ,DXYZ,FRC,VERT,IV(I),EYE,LSHINE,LIGHT)
        ELSEIF (ISF.EQ.4) THEN
          K4=(IV(1)+IV(2)+IV(3)+IV(4)-6)/4
          IF ((IV(2)-IV(1)).EQ.3) K4=6
          CALL SBSF5B(XYZ,DXYZ,FRC,VERT,K4,EYE,LSHINE,LIGHT,2)
        ELSEIF (ISF.EQ.6) THEN
          IF (IV(3).LE.4) THEN
            K3=IV(1)+IV(2)+IV(3)-6
            K4=MOD((IV(4)+K3),4)+3*K3
          ELSE
            IF (IV(2).GE.5) THEN
              K3=IV(2)+IV(3)+IV(4)-18
              K4=IV4MAP(MOD((IV(1)+K3),4)+3*K3)
            ELSE
              K4=12+IV(3)-IV(2)
              IF ((IV(1)+IV(2)+IV(3)+IV(4)).EQ.22) K4=29-K4
            ENDIF
          ENDIF
          CALL SBSF5D(XYZ,DXYZ,FRC,VERT,K4,EYE,LSHINE,LIGHT)
        ELSE
          K4=IV(1)+IV(2)+IV(3)+IV(4)
          IF (K4.EQ.16 .OR. K4.EQ.20) THEN          
            CALL SBSF5B(XYZ,DXYZ,FRC,VERT,IV(3),EYE,LSHINE,LIGHT,1)
            CALL SBSF5B(XYZ,DXYZ,FRC,VERT,IV(4),EYE,LSHINE,LIGHT,1)
          ELSEIF (K4.EQ.18) THEN
            K4A=IV(1)
            IF ((IV(2)-K4A).EQ.3) K4A=4
            K4B=9+MOD(K4A+1,4)
            CALL SBSF5B(XYZ,DXYZ,FRC,VERT,K4A,EYE,LSHINE,LIGHT,1)
            CALL SBSF5B(XYZ,DXYZ,FRC,VERT,K4B,EYE,LSHINE,LIGHT,1)
          ELSE
            IF (K4.EQ.14) THEN
              K4A=IV(4)
              K3=IV(1)+IV(2)+IV(3)-5+2*(IV(1)/5+2*(IV(2)/5)+IV(3)/5)
            ELSEIF (K4.EQ.22) THEN
              K4A=IV(1)
              K3=IV(2)+IV(3)+IV(4)-5+2*(IV(2)/5+2*(IV(3)/5)+IV(4)/5)
            ELSE
              IF (MOD((IV(1)+IV(2)),2).EQ.0) THEN
                IF (IV(4).EQ.6 .OR. (IV(3)-IV(2)).EQ.2) THEN
                  K4A=IV(2)
                  K3=IV(1)+IV(3)+IV(4)-5+2*(IV(1)/5+2*(IV(3)/5)+IV(4)/5)
                ELSE
                  K4A=IV(1)
                  K3=IV(2)+IV(3)+IV(4)-5+2*(IV(2)/5+2*(IV(3)/5)+IV(4)/5)
                ENDIF
              ELSE
                IF (IV(1).EQ.3 .OR. (IV(3)-IV(2)).EQ.2) THEN
                  K4A=IV(3)
                  K3=IV(1)+IV(2)+IV(4)-5+2*(IV(1)/5+2*(IV(2)/5)+IV(4)/5)
                ELSE
                  K4A=IV(4)
                  K3=IV(1)+IV(2)+IV(3)-5+2*(IV(1)/5+2*(IV(2)/5)+IV(3)/5)
                ENDIF
              ENDIF
            ENDIF
            CALL SBSF5A(XYZ,DXYZ,FRC,VERT,K4A,EYE,LSHINE,LIGHT)
            CALL SBSF5C(XYZ,DXYZ,FRC,VERT,K3,EYE,LSHINE,LIGHT)
          ENDIF
        ENDIF
      ENDIF
      END
C
      SUBROUTINE SBSF5A(XYZ,DXYZ,FRC,VERT,IV,EYE,LSHINE,LIGHT)
C     ---------------------------------------------------------
C
      REAL    XYZ(*),DXYZ(3,12,*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      LOGICAL LSHINE
C
      IF (IV.LE.4) THEN
        J=IV
        K=1+MOD(IV+2,4)
        L=IV+4
      ELSE
        J=IV+4
        K=9+MOD(IV-2,4)
        L=IV
      ENDIF
      DO 10 I=1,3
        VERT(I,1)=XYZ(I)+DXYZ(I,J,1)+FRC(J)*DXYZ(I,J,2)
        VERT(I,2)=XYZ(I)+DXYZ(I,K,1)+FRC(K)*DXYZ(I,K,2)
        VERT(I,3)=XYZ(I)+DXYZ(I,L,1)+FRC(L)*DXYZ(I,L,2)
  10  CONTINUE
      CALL SBSRF6(EYE,3,VERT,LSHINE,LIGHT,1) 
      END  
C
      SUBROUTINE SBSF5B(XYZ,DXYZ,FRC,VERT,KK,EYE,LSHINE,LIGHT,LL)
C     ------------------------------------------------------------
C
      REAL    XYZ(*),DXYZ(3,12,*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      INTEGER IVL(4,12,2)
      LOGICAL LSHINE
      DATA    IVL /5,6,2,4,6,7,3,1,7,8,4,2,8,5,1,3,9,1,4,12,10,2,1,9,
     *     11,3,2,10,12,4,3,11,12,5,6,10,9,6,7,11,10,7,8,12,11,8,5,9,
     *     5,6,7,8,4,2,10,12,1,3,11,9,4,2,10,12,5,6,7,8,1,3,11,9,24*0/
C
      J=IVL(1,KK,LL)
      K=IVL(2,KK,LL)
      L=IVL(3,KK,LL)
      M=IVL(4,KK,LL)
      DO 10 I=1,3
        VERT(I,1)=XYZ(I)+DXYZ(I,J,1)+FRC(J)*DXYZ(I,J,2)
        VERT(I,2)=XYZ(I)+DXYZ(I,K,1)+FRC(K)*DXYZ(I,K,2)
        VERT(I,3)=XYZ(I)+DXYZ(I,L,1)+FRC(L)*DXYZ(I,L,2)
        VERT(I,4)=XYZ(I)+DXYZ(I,M,1)+FRC(M)*DXYZ(I,M,2)
        VERT(I,5)=VERT(I,1)
        VERT(I,6)=0.25*(VERT(I,1)+VERT(I,2)+VERT(I,3)+VERT(I,4))
  10  CONTINUE
      DO 20 I=1,4
        CALL SBRCOP(VERT(1,I),VERT(1,7),6)
        CALL SBSRF6(EYE,3,VERT(1,6),LSHINE,LIGHT,1)
  20  CONTINUE
      END  
C
      SUBROUTINE SBSF5C(XYZ,DXYZ,FRC,VERT,K3,EYE,LSHINE,LIGHT)
C     ---------------------------------------------------------
C
      REAL    XYZ(*),DXYZ(3,12,*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      INTEGER IV3(5,24)
      LOGICAL LSHINE
      DATA    IV3 /5,6,7,3,4,8,5,6,2,3,7,8,5,1,2,6,7,8,4,1,
     *             12,4,2,6,9,4,2,10,9,5,9,1,3,8,12,9,1,3,7,10,
     *             1,3,11,10,6,1,3,11,12,5,4,2,10,11,8,2,4,12,11,7,
     *             10,12,4,1,6,12,10,2,1,5,11,9,1,4,8,1,9,11,7,2,
     *             3,11,9,6,2,3,11,9,5,4,2,10,12,8,3,4,12,10,7,3,
     *             5,6,7,11,12,8,5,6,10,11,7,8,5,9,10,6,7,8,12,9/
C
      J=IV3(1,K3)
      K=IV3(2,K3)
      L=IV3(3,K3)
      M=IV3(4,K3)
      N=IV3(5,K3)
      DO 10 I=1,3
        VERT(I,1)=XYZ(I)+DXYZ(I,J,1)+FRC(J)*DXYZ(I,J,2)
        VERT(I,2)=XYZ(I)+DXYZ(I,K,1)+FRC(K)*DXYZ(I,K,2)
        VERT(I,3)=XYZ(I)+DXYZ(I,L,1)+FRC(L)*DXYZ(I,L,2)
        VERT(I,4)=XYZ(I)+DXYZ(I,M,1)+FRC(M)*DXYZ(I,M,2)
        VERT(I,5)=XYZ(I)+DXYZ(I,N,1)+FRC(N)*DXYZ(I,N,2)
        VERT(I,6)=VERT(I,1)
        VERT(I,7)=0.2*(VERT(I,1)+VERT(I,2)+VERT(I,3)+VERT(I,4)+
     *                 VERT(I,5))
  10  CONTINUE
      DO 20 I=1,5
        CALL SBRCOP(VERT(1,I),VERT(1,8),6)
        CALL SBSRF6(EYE,3,VERT(1,7),LSHINE,LIGHT,1)
  20  CONTINUE
      END
C
      SUBROUTINE SBSF5D(XYZ,DXYZ,FRC,VERT,K4,EYE,LSHINE,LIGHT)
C     ---------------------------------------------------------
C
      REAL    XYZ(*),DXYZ(3,12,*),FRC(*),VERT(3,*),EYE(*),LIGHT(*)
      INTEGER IV4(6,16)
      LOGICAL LSHINE
      DATA    IV4 /12,9,6,7,3,4,10,9,5,4,3,7,11,10,6,5,4,3,
     *             11,12,5,6,2,3,9,12,8,3,2,6,10,9,5,8,3,2,
     *             10,11,8,5,1,2,12,11,7,2,1,5,9,12,8,7,2,1,
     *             9,10,7,8,4,1,11,10,6,1,4,8,12,11,7,6,1,4,
     *             12,10,6,1,3,8,12,10,7,3,1,5,11,9,6,2,4,8,
     *             11,9,5,4,2,7/
      DATA    VNORM /0.1666666667/
C
      CALL SBRFIL(VERT(1,8),0.0,3)
      DO 20 J=1,6
        K=IV4(J,K4)
        DO 10 I=1,3
          VERT(I,J)=XYZ(I)+DXYZ(I,K,1)+FRC(K)*DXYZ(I,K,2)
          VERT(I,8)=VERT(I,8)+VERT(I,J)
  10    CONTINUE
  20  CONTINUE
      CALL SBRCOP(VERT,VERT(1,7),3)
      DO 30 I=1,3
  30    VERT(I,8)=VERT(I,8)*VNORM
      DO 40 I=1,6
        CALL SBRCOP(VERT(1,I),VERT(1,9),6)
        CALL SBSRF6(EYE,3,VERT(1,8),LSHINE,LIGHT,1)
  40  CONTINUE
      END
C
      SUBROUTINE SBSRF6(EYE,NV,VERT,LSHINE,LIGHT,INSIDE)
C     --------------------------------------------------
C
      REAL             EYE(*),VERT(3,*),LIGHT(*)
      LOGICAL          LSHINE
C
      REAL*8           XLNORM,ZZ,DZZ
      REAL             XW(20),YW(20),MTRX
      LOGICAL          LPS,LCOLOR
      COMMON  /SRFCOM/ GRDCUB(3,8),MTRX(3,3),ORIG(3),XL2,COL0,COLSCL
      COMMON  /SFTBUF/ SBBUFF(2000000),NXP,NYP,IBFMOD,KSTART,
     *                 IR(0:255),IG(0:255),IB(0:255),LPS,LCOLOR,
     *                 XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                 YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC
C
C Carry out some initial checks and calculate the coordinates of the 
C projected triangle.
C
      IF (NV.LT.3 .OR. NV.GT.10) RETURN
      SMALL=1.0E-10
      XMIN=+1.0E20
      XMAX=-1.0E20
      YMIN=+1.0E20
      YMAX=-1.0E20
      DO 10 I=1,NV
        CALL SBLIN1(EYE,VERT(1,I),VERT(2,I),VERT(3,I),XW(I),YW(I))
        IF (XW(I).LT.XMIN) THEN
          XMIN=XW(I)
          ILEFT=I
        ENDIF
        IF (YW(I).LT.YMIN) THEN
          YMIN=YW(I)
          JBOTOM=I
        ENDIF
        XMAX=MAX(XW(I),XMAX)
        YMAX=MAX(YW(I),YMAX)
  10  CONTINUE
      IF (XMIN.GE.XTRC .OR. XMAX.LE.XBLC) RETURN
      IF (YMIN.GE.YTRC .OR. YMAX.LE.YBLC) RETURN
C
C Find the outward normal seen by the eye.
C
      AX=VERT(1,2)-VERT(1,1)
      AY=VERT(2,2)-VERT(2,1)
      AZ=VERT(3,2)-VERT(3,1)
      BX=VERT(1,1)-VERT(1,NV)
      BY=VERT(2,1)-VERT(2,NV)
      BZ=VERT(3,1)-VERT(3,NV)
      XN=BY*AZ-AY*BZ
      YN=BZ*AX-AZ*BX
      ZN=BX*AY-AX*BY
      TEN=XN*(EYE(1)-VERT(1,1))+YN*(EYE(2)-VERT(2,1))
     *   +ZN*(EYE(3)-VERT(3,1))
      IF (TEN.LT.0.0) THEN
        XN=-XN
        YN=-YN
        ZN=-ZN
        TEN=-TEN
      ENDIF
C
C Plot the projected triangle.
C
      XLNORM=DBLE(TEN)
      EYENRM=XN*EYE(1)+YN*EYE(2)+ZN*EYE(3)
      DX=FLOAT(NXP-1)/(XTRC-XBLC)
      DY=FLOAT(NYP-1)/(YTRC-YBLC)
      DYJ=1.0/DY
      DXI=1.0/DX
      SAFER=0.0001
      IF ((XMAX-XMIN).GT.(YMAX-YMIN)) THEN
        JYMIN=INT((YMIN-YBLC)*DY)+2
        JYMAX=MIN(INT((YMAX-YBLC)*DY)+1,NYP)
        IF (JYMIN.GT.JYMAX) RETURN
        YJ=YBLC+(FLOAT(JYMIN-1)+SAFER)*DYJ
        NVL2=JBOTOM
        NVR2=JBOTOM
        J1=JYMIN
        DO 40 IVERT=1,NV
          IF (YJ.GT.YW(NVL2)) THEN
   1        NVL1=NVL2
            NVL2=NVL1-1
            IF (NVL2.LT.1) NVL2=NV
            IF (NVL2.EQ.JBOTOM) RETURN
            IF (YJ.GT.YW(NVL2)) GOTO 1
            YDIFL=YW(NVL2)-YW(NVL1)
            IF (ABS(YDIFL).LT.SMALL) YDIFL=SMALL
            GRADL=(XW(NVL2)-XW(NVL1))/YDIFL
          ENDIF
          IF (YJ.GT.YW(NVR2)) THEN
   2        NVR1=NVR2
            NVR2=NVR1+1
            IF (NVR2.GT.NV) NVR2=1
            IF (NVR2.EQ.JBOTOM) RETURN
            IF (YJ.GT.YW(NVR2)) GOTO 2
            YDIFR=YW(NVR2)-YW(NVR1)
            IF (ABS(YDIFR).LT.SMALL) YDIFR=SMALL
            GRADR=(XW(NVR2)-XW(NVR1))/YDIFR
          ENDIF
          IF (YW(NVL2).LT.YW(NVR2)) THEN
            J2=MIN(INT((YW(NVL2)-YBLC)*DY)+1,JYMAX)
          ELSE
            J2=MIN(INT((YW(NVR2)-YBLC)*DY)+1,JYMAX)
          ENDIF
          DO 30 J=J1,J2
            IF (J.GE.1) THEN
              XL=XW(NVL1)+GRADL*(YJ-YW(NVL1))
              XR=XW(NVR1)+GRADR*(YJ-YW(NVR1))
              ISTEP=1
              IX1=MAX(INT((XL-XBLC)*DX)+2,1)
              IX2=MIN(INT((XR-XBLC)*DX)+1,NXP)
              IF (IX1.GT.IX2) THEN
                ISTEP=-1
                IX1=MIN(IX1-1,NXP)
                IX2=MAX(IX2+1,1)
              ENDIF
              XI=XBLC+FLOAT(IX1-1)*DXI
              SDXI=FLOAT(ISTEP)*DXI
              DZZ=DBLE(SDXI*XN)
              ZZ=DBLE(EYENRM-XI*XN-YJ*YN)
              K=(J-1)*NXP+IX1
              DO 20 I=IX1,IX2,ISTEP
                XLAMDA=SNGL(XLNORM/ZZ)
                Z=EYE(3)*(1.0-XLAMDA)
                IF (Z.GT.SBBUFF(K)) THEN
                  SBBUFF(K)=Z
                  X=EYE(1)+XLAMDA*(XI-EYE(1))
                  Y=EYE(2)+XLAMDA*(YJ-EYE(2))
                  IF (INSIDE.EQ.0) THEN
                    GX=XN
                    GY=YN
                    GZ=ZN
                  ELSE
                    CALL SBSF6A(X,Y,Z,ORIG,MTRX,GRDCUB,GX,GY,GZ)
                  ENDIF
                  CALL SBSF6B(EYE,X,Y,Z,GX,GY,GZ,LIGHT,XL2,LSHINE,CLR)
                  SBBUFF(KSTART+K)=COL0+COLSCL*CLR
                ENDIF
                XI=XI+SDXI
                ZZ=ZZ-DZZ
                K=K+ISTEP
  20          CONTINUE
            ENDIF
            YJ=YJ+DYJ
  30      CONTINUE
          J1=J2+1
          IF (J1.GT.JYMAX) RETURN
  40    CONTINUE
      ELSE
        IXMIN=INT((XMIN-XBLC)*DX)+2
        IXMAX=MIN(INT((XMAX-XBLC)*DX)+1,NXP)
        IF (IXMIN.GT.IXMAX) RETURN
        XI=XBLC+(FLOAT(IXMIN-1)+SAFER)*DXI
        NVL2=ILEFT
        NVR2=ILEFT
        I1=IXMIN
        DO 70 IVERT=1,NV
          IF (XI.GT.XW(NVL2)) THEN
   3        NVL1=NVL2
            NVL2=NVL1-1
            IF (NVL2.LT.1) NVL2=NV
            IF (NVL2.EQ.ILEFT) RETURN
            IF (XI.GT.XW(NVL2)) GOTO 3
            XDIFL=XW(NVL2)-XW(NVL1)
            IF (ABS(XDIFL).LT.SMALL) XDIFL=SMALL
            GRADL=(YW(NVL2)-YW(NVL1))/XDIFL
          ENDIF
          IF (XI.GT.XW(NVR2)) THEN
   4        NVR1=NVR2
            NVR2=NVR1+1
            IF (NVR2.GT.NV) NVR2=1
            IF (NVR2.EQ.ILEFT) RETURN
            IF (XI.GT.XW(NVR2)) GOTO 4
            XDIFR=XW(NVR2)-XW(NVR1)
            IF (ABS(XDIFR).LT.SMALL) XDIFR=SMALL
            GRADR=(YW(NVR2)-YW(NVR1))/XDIFR
          ENDIF
          IF (XW(NVL2).LT.XW(NVR2)) THEN
            I2=MIN(INT((XW(NVL2)-XBLC)*DX)+1,IXMAX)
          ELSE
            I2=MIN(INT((XW(NVR2)-XBLC)*DX)+1,IXMAX)
          ENDIF
          DO 60 I=I1,I2
            IF (I.GE.1) THEN
              YL=YW(NVL1)+GRADL*(XI-XW(NVL1))
              YR=YW(NVR1)+GRADR*(XI-XW(NVR1))
              ISTEP=1
              JY1=MAX(INT((YL-YBLC)*DY)+2,1)
              JY2=MIN(INT((YR-YBLC)*DY)+1,NYP)
              IF (JY1.GT.JY2) THEN
                ISTEP=-1
                JY1=MIN(JY1-1,NYP)
                JY2=MAX(JY2+1,1)
              ENDIF
              YJ=YBLC+FLOAT(JY1-1)*DYJ
              SDYJ=FLOAT(ISTEP)*DYJ
              DZZ=DBLE(SDYJ*YN)
              ZZ=DBLE(EYENRM-YJ*YN-XI*XN)
              K=(JY1-1)*NXP+I
              KSTEP=ISTEP*NXP
              DO 50 J=JY1,JY2,ISTEP
                XLAMDA=SNGL(XLNORM/ZZ)
                Z=EYE(3)*(1.0-XLAMDA)
                IF (Z.GT.SBBUFF(K)) THEN
                  SBBUFF(K)=Z
                  X=EYE(1)+XLAMDA*(XI-EYE(1))
                  Y=EYE(2)+XLAMDA*(YJ-EYE(2))
                  IF (INSIDE.EQ.0) THEN
                    GX=XN
                    GY=YN
                    GZ=ZN
                  ELSE
                    CALL SBSF6A(X,Y,Z,ORIG,MTRX,GRDCUB,GX,GY,GZ)
                  ENDIF
                  CALL SBSF6B(EYE,X,Y,Z,GX,GY,GZ,LIGHT,XL2,LSHINE,CLR)
                  SBBUFF(KSTART+K)=COL0+COLSCL*CLR
                ENDIF
                YJ=YJ+SDYJ
                ZZ=ZZ-DZZ
                K=K+KSTEP
  50          CONTINUE
            ENDIF
            XI=XI+DXI
  60      CONTINUE
          I1=I2+1
          IF (I1.GT.IXMAX) RETURN
  70    CONTINUE
      ENDIF
      END
C
      SUBROUTINE SBSF6A(X,Y,Z,ORIG,MTRX,GRD,XN,YN,ZN)
C     -----------------------------------------------
C
      REAL ORIG(*),MTRX(3,*),GRD(3,*)
      DATA ZERO,ONE /0.00001,0.99999/
C
      X0=X-ORIG(1)
      Y0=Y-ORIG(2)
      Z0=Z-ORIG(3)
      XI=X0*MTRX(1,1)+Y0*MTRX(2,1)+Z0*MTRX(3,1)
      YJ=X0*MTRX(1,2)+Y0*MTRX(2,2)+Z0*MTRX(3,2)
      ZK=X0*MTRX(1,3)+Y0*MTRX(2,3)+Z0*MTRX(3,3)
      DX=MIN(ONE,MAX(ZERO,XI-FLOAT(INT(XI))))
      DY=MIN(ONE,MAX(ZERO,YJ-FLOAT(INT(YJ))))
      DZ=MIN(ONE,MAX(ZERO,ZK-FLOAT(INT(ZK))))
      XN1=(1.0-DX)*(GRD(1,1)+DY*(GRD(1,4)-GRD(1,1)))
     *        +DX*(GRD(1,2)+DY*(GRD(1,3)-GRD(1,2)))
      XN2=(1.0-DX)*(GRD(1,5)+DY*(GRD(1,8)-GRD(1,5)))
     *        +DX*(GRD(1,6)+DY*(GRD(1,7)-GRD(1,6)))
      XN=XN1+DZ*(XN2-XN1)
      YN1=(1.0-DX)*(GRD(2,1)+DY*(GRD(2,4)-GRD(2,1)))
     *        +DX*(GRD(2,2)+DY*(GRD(2,3)-GRD(2,2)))
      YN2=(1.0-DX)*(GRD(2,5)+DY*(GRD(2,8)-GRD(2,5)))
     *        +DX*(GRD(2,6)+DY*(GRD(2,7)-GRD(2,6)))
      YN=YN1+DZ*(YN2-YN1)
      ZN1=(1.0-DX)*(GRD(3,1)+DY*(GRD(3,4)-GRD(3,1)))
     *        +DX*(GRD(3,2)+DY*(GRD(3,3)-GRD(3,2)))
      ZN2=(1.0-DX)*(GRD(3,5)+DY*(GRD(3,8)-GRD(3,5)))
     *        +DX*(GRD(3,6)+DY*(GRD(3,7)-GRD(3,6)))
      ZN=ZN1+DZ*(ZN2-ZN1)
      END
C
      SUBROUTINE SBSF6B(EYE,X,Y,Z,XN,YN,ZN,LIGHT,XL2,LSHINE,COLOUR)
C     -------------------------------------------------------------
C
      REAL    EYE(*),LIGHT(*)
      LOGICAL LSHINE
      DATA    SMALL2 /1.0E-20/
C
      COLOUR=0.0
      XNL=XN*LIGHT(1)+YN*LIGHT(2)+ZN*LIGHT(3)
      IF (XNL.GE.0.0) RETURN
      XN2=XN**2+YN**2+ZN**2+SMALL2
      IF (LSHINE) THEN
        RFNORM=2.0*XNL/XN2
        RX=LIGHT(1)-XN*RFNORM
        RY=LIGHT(2)-YN*RFNORM
        RZ=LIGHT(3)-ZN*RFNORM
        VX=EYE(1)-X
        VY=EYE(2)-Y
        VZ=EYE(3)-Z
        XRV=RX*VX+RY*VY+RZ*VZ
        IF (XRV.LT.0.0) RETURN
        V2=VX**2+VY**2+VZ**2
        COLOUR=MIN(XRV**2/(ABS(XL2*V2)+SMALL2),1.0)
      ELSE
        COLOUR=MIN(-XNL/SQRT(ABS(XL2*XN2)+SMALL2),1.0)
      ENDIF
      END
C
C
C***<vector utility routines>*******************************************
C
      SUBROUTINE SBRFIL(X,A,N)
C     ------------------------
C
      REAL X(*)
C
      DO 10 I=1,N
  10    X(I)=A
      END
C
      SUBROUTINE SBRCOP(X,Y,N)
C     ------------------------
C
      REAL X(*),Y(*)
C
      DO 10 I=1,N
  10    Y(I)=X(I)
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
