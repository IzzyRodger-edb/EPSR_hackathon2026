CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE PGCELL(A,IDIM,JDIM,I1,I2,J1,J2,FG,BG,TR,NCOLS,R,G,B)
C     ---------------------------------------------------------------
C
      REAL         A(IDIM,JDIM),TR(6)
      REAL         R(0:NCOLS-1),G(0:NCOLS-1),B(0:NCOLS-1)
      INTEGER      IDIM,JDIM,I1,I2,J1,J2,NCOLS
C
      INCLUDE     'grpckg1.inc'
      INTEGER      IR(0:255),IG(0:255),IB(0:255)
      LOGICAL      LPS,LCOLOR
      CHARACTER*16 TYPE,CHR
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      This subroutine is designed to do the job of CELL_ARRAY in GKS;
C   that is, it shades elements of a rectangular array with the 
C   appropriate colours passed down in the RGB colour-table. Essentially,
C   it is a colour version of PGGRAY. 
C      The colour-index used for particular array pixel is given by:
C         Colour Index = NINT{[A(i,j)-BG/(FG-BG)]*FLOAT(NCOLS-1)} ,
C   with truncation at 0 and NCOLS-1, as necessary.
C      The transform matrix TR is used to calculate the (bottom left) 
C   world coordinates of the cell which represents each array element:
C         X = TR(1) + TR(2)*I + TR(3)*J
C         Y = TR(4) + TR(5)*I + TR(6)*J  .
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    A        R*4    I   IDIMxJDIM  The array to be plotted.
C    IDIM     I*4    I       -      The first dimension of array A.
C    JDIM     I*4    I       -      The second dimension of array A.
C    I1,I2    I*4    I       -      The inclusive range of the first
C                                   index (I) to be plotted.
C    J1,J2    I*4    I       -      The inclusive range of the second
C                                   index (J) to be plotted.
C    FG       R*4    I       -      The array value which is to appear
C                                   with shade 1 ("foreground").
C    BG       R*4    I       -      The array value which is to appear
C                                   with shade 0 ("background").
C    TR       R*4    I       6      Transformation matrix between array
C                                   grid and world coordinates.
C    NCOLS    I*4    I       -      Number of colours in colour-table.
C    R        R*4    I      NCOLS   Red intensity for colour-table.
C    G        R*4    I      NCOLS   Green intensity for colour-table.
C    B        R*4    I      NCOLS   Blue intensity for colour-table.
C
C Globals
C    GRPCKG1.INC
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     PGQCOL     Inquires about the colour capability.
C     PGQCI      Inquires about the current colour index.
C     PGQINF     Inquires about general PGPLOT information.
C     PGBBUF     Recommended initial call (to start a PGPLOT buffer).
C     PGEBUF     Recommended final call (to end a PGPLOT buffer).
C     PGCLPX     Pixel-device support subroutine for PGCELL.
C     PGCLPS     PostScript support subroutine for PGCELL.
C     GREXEC     Dispatches command to appropriate device driver.
C     DSQINF     Inquires about viewport and window dimensions.
C
C History
C   D. S. Sivia       3 Jul 1992  Initial release.
C   D. S. Sivia       7 Feb 1995  Uses GRGRAY approach instead of PGPOLY.
C   D. S. Sivia       6 Mar 1995  Slight changes for Postscript output.
C   D. S. Sivia       1 Aug 1996  Replaced pgplot.inc with DSQINF!
C   D. S. Sivia      21 Oct 1997  Made slightly friendlier for NT.
C-----------------------------------------------------------------------
C
C Find out device-type. If not Postscript, then (i) return if less than
C 16 shades available; (ii) save initial colour index (and hope it's 
C less than ICLOW).
C
      NC=256
      LPS=.TRUE.
      CALL PGQINF('TYPE',TYPE,LCHR)
      IF (TYPE.EQ.'PS' .OR. TYPE.EQ.'VPS') THEN
        LCOLOR=.FALSE.
      ELSEIF (TYPE.EQ.'CPS' .OR. TYPE.EQ.'VCPS') THEN
        LCOLOR=.TRUE.
      ELSE
        LPS=.FALSE.
        CALL PGQCOL(IC1,IC2)
        NC=IC2-IC1+1
        IF (NC.LT.16) THEN
          WRITE(*,*)' *** Not enough colours available on this device!'
          RETURN
        ELSE
          ICLOW=4
          IF (NC.GE.96) ICLOW=16
          IC1=IC1+ICLOW
          NC=NC-ICLOW
        ENDIF
        CALL PGQCI(ICSAVE)
        IF (ICSAVE.GE.IC1) ICSAVE=IC1+1
      ENDIF
      CALL PGBBUF
C
C Activate the colour table. If NCOLS is less than the number of colours
C available, simply assign NCOLS; otherwise, use a linear interpolation.
C
      IF (NCOLS.LE.NC) THEN
        NC=NCOLS
        DO 10 I=0,NCOLS-1
          CALL PGSCR(IC1+I,R(I),G(I),B(I))
          IR(I)=NINT(R(I)*255.0)
          IG(I)=NINT(G(I)*255.0)
          IB(I)=NINT(B(I)*255.0)
  10    CONTINUE
      ELSE
        COL=0.0
        DCOL=0.999*FLOAT(NCOLS-1)/FLOAT(NC-1)
        DO 20 I=0,NC-1
          ICOL=INT(COL)
          DICOL=COL-FLOAT(ICOL)
          RL=R(ICOL)+DICOL*(R(ICOL+1)-R(ICOL))
          GL=G(ICOL)+DICOL*(G(ICOL+1)-G(ICOL))
          BL=B(ICOL)+DICOL*(B(ICOL+1)-B(ICOL))
          CALL PGSCR(IC1+I,RL,GL,BL)
          IR(I)=NINT(RL*255.0)
          IG(I)=NINT(GL*255.0)
          IB(I)=NINT(BL*255.0)
          COL=COL+DCOL
  20    CONTINUE
      ENDIF
      ASCALE=FLOAT(NC-1)/(FG-BG)
C
C Check to see whether a pixel device or Postscript is being used, and 
C call the appropriate PGCELL support subrotuine.
C
      IF (LPS) THEN
        CALL PGCLPS(A,IDIM,JDIM,I1,I2,J1,J2,BG,TR,ASCALE,IR,IG,IB,NC,
     *              LCOLOR)
      ELSE
        NBUF=0
        LCHR=LEN(CHR)
        CALL GREXEC(GRGTYP,4,RBUF,NBUF,CHR,LCHR)
        IF (CHR(7:7).EQ.'P') THEN
          CALL PGCLPX(A,IDIM,JDIM,I1,I2,J1,J2,BG,TR,ASCALE,IC1,NC)
        ELSE
          WRITE(*,*)' Sorry, PGCELL does not support this device!'
        ENDIF
      ENDIF
C
C Reset the initial colour index.
C
      IF (.NOT. LPS) CALL PGSCI(ICSAVE)
      CALL PGEBUF
      END
C
      SUBROUTINE PGCLPX(A,IDIM,JDIM,I1,I2,J1,J2,BG,TR,ASCALE,IC1,NC)
C     --------------------------------------------------------------
C
C Light-up the device pixels, with colours determined by array A.
C
      INCLUDE  'grpckg1.inc'
      REAL      A(IDIM,*),TR(*),BUFFER(2050)
      CHARACTER CHR*16
C
      CALL DSQINF(XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *            YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC)
      IX1=NINT(XOFF)
      IX2=NINT(XOFF+XLEN)
      JY1=NINT(YOFF)
      JY2=NINT(YOFF+YLEN)
      DET=TR(2)*TR(6)-TR(3)*TR(5)
      TR11=+TR(6)/DET
      TR12=-TR(3)/DET
      TR21=-TR(5)/DET
      TR22=+TR(2)/DET
      DX=1.0/XSCALE
      DY=1.0/YSCALE
      X1=(XOFF-XORG)*DX-TR(1)
      Y1=(YOFF-YORG)*DY-TR(4)
      DXX=TR11*DX
      DXY=TR12*DY
      DYX=TR21*DX
      DYY=TR22*DY
      XI1=TR11*X1+TR12*Y1
      YJ1=TR21*X1+TR22*Y1
      DO 40 JY=JY1,JY2
        XI=XI1
        YJ=YJ1
        BUFFER(2)=FLOAT(JY)
        NPIX=2
        DO 30 IX=IX1,IX2
          I=INT(XI)
          J=INT(YJ)
          IF (I.GE.I1 .AND. I.LE.I2 .AND. J.GE.J1 .AND. J.LE.J2) THEN
            IF (NPIX.EQ.2) BUFFER(1)=FLOAT(IX)
            K=NINT((A(I,J)-BG)*ASCALE)
            IF (K.LT.0) K=0
            IF (K.GE.NC) K=NC-1
            NPIX=NPIX+1
            BUFFER(NPIX)=FLOAT(K+IC1)
          ENDIF
          XI=XI+DXX
          YJ=YJ+DYX
  30    CONTINUE
        CALL GREXEC(GRGTYP,26,BUFFER,NPIX,CHR,LCHR)
        XI1=XI1+DXY
        YJ1=YJ1+DYY
  40  CONTINUE
      END
C
      SUBROUTINE PGCLPS(A,IDIM,JDIM,I1,I2,J1,J2,BG,TR,ASCALE,IR,IG,IB,
     *                  NC,LCOLOR)
C     ----------------------------------------------------------------
C
C Postscript support subroutine for PGCELL.
C
      INCLUDE  'grpckg1.inc'
      REAL      A(IDIM,*),TR(*)
      INTEGER   IR(0:*),IG(0:*),IB(0:*)
      INTEGER   VALUE(33)
      LOGICAL   LCOLOR
      CHARACTER INLINE*80
C
C Set clipping rectangle in device.
C
      CALL DSQINF(XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *            YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC)
      WRITE (INLINE,100) NINT(XOFF),NINT(YOFF),NINT(XLEN),NINT(YLEN),
     *          -NINT(XLEN)
 100  FORMAT(I6,I6,' moveto ',I6, ' 0 rlineto  0 ',I6,' rlineto ',
     *I6,' 0 rlineto')
      CALL GRTERM
      CALL GRESC(' newpath ')
      CALL GRESC(INLINE)
      CALL GRESC(' closepath clip ')
C
C Work out the nunmber of X and Y pixels for PS image, with NDOTS per 
C inch, and build an image transformation matrix.
C
      NDOTS=100
C      IF (LCOLOR) NDOTS=50
      NXP=NINT(FLOAT(NDOTS)*XLEN/XPERIN)
      NYP=NINT(FLOAT(NDOTS)*YLEN/YPERIN)
      DXI=(FLOAT(I2-I1)+0.999)/FLOAT(NXP)
      DYJ=(FLOAT(J2-J1)+0.999)/FLOAT(NYP)
      DET=TR(2)*TR(6)-TR(3)*TR(5)
      TR11=+TR(6)/DET
      TR12=-TR(3)/DET
      TR21=-TR(5)/DET
      TR22=+TR(2)/DET
      AT=TR11/(XSCALE*DXI)
      BT=TR21/(XSCALE*DYJ)
      CT=TR12/(YSCALE*DXI)
      DT=TR22/(YSCALE*DYJ)
      TX=-(TR11*(XORG/XSCALE+TR(1))+TR12*(YORG/YSCALE+TR(4))+I1)/DXI
      TY=-(TR21*(XORG/XSCALE+TR(1))+TR22*(YORG/YSCALE+TR(4))+J1)/DYJ
C
C Use a PostScript "image" operator.
C
      WRITE (INLINE, '(A,I5,A)') '/picstr ',NXP,' string def'
      CALL GRESC(INLINE)
      WRITE (INLINE,110) NXP,NYP,AT,BT,CT,DT,TX,TY
 110  FORMAT(2I4,' 8 [',6(1PE10.3,' '),']')
      CALL GRESC(INLINE)
      CALL GRESC('{ currentfile picstr readhexstring pop}')
      IF (LCOLOR) THEN
        CALL GRESC('  false 3 colorimage')
      ELSE
        CALL GRESC('  image')
      ENDIF
C
C Write out the image array in hexadecimal.
C
      ASCALE=ASCALE*255.0/FLOAT(NC-1)
      YJ=FLOAT(J1)
      DO 20 JP=1,NYP
        J=INT(YJ)
        XI=FLOAT(I1)
        L=0
        DO 10 IP=1,NXP
          L=L+1
          I=INT(XI)
          IC=NINT((A(I,J)-BG)*ASCALE)
          IF (IC.LT.0) IC=0
          IF (IC.GE.NC) IC=NC-1
          IF (LCOLOR) THEN
            VALUE(L)=IR(IC)
            VALUE(L+1)=IG(IC)
            VALUE(L+2)=IB(IC)
            L=L+2
          ELSE
            VALUE(L)=(IR(IC)+IG(IC)+IB(IC))/3
          ENDIF
          IF (L.EQ.33) THEN
            WRITE(INLINE,120) (VALUE(K),K=1,33)
 120        FORMAT(33Z2.2)
            CALL GRESC(INLINE(1:66))
            L=0
          ENDIF
          XI=XI+DXI
   10   CONTINUE
        IF (L.NE.0) THEN
          WRITE(INLINE,120) (VALUE(K),K=1,L)
          CALL GRESC(INLINE(1:2*L))
        ENDIF
        YJ=YJ+DYJ
   20 CONTINUE
      CALL GRESC(' newpath initclip ')
      CALL GRTERM
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
