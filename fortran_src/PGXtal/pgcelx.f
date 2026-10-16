CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE PGCELX(A,IDIM,JDIM,I1,I2,J1,J2,FG,BG,TRNL,NCOLS,R,G,B)
C     -----------------------------------------------------------------
C
      REAL         A(IDIM,JDIM)
      REAL         R(0:NCOLS-1),G(0:NCOLS-1),B(0:NCOLS-1)
      INTEGER      IDIM,JDIM,I1,I2,J1,J2,NCOLS
      EXTERNAL     TRNL
C
      INCLUDE     'grpckg1.inc'
      INTEGER      IR(0:255),IG(0:255),IB(0:255)
      LOGICAL      LPS,LCOLOR
      CHARACTER*16 TYPE,CHR
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      This subroutine is a non-linear version of PGCELL. That is to 
C   say, it paints a linearly-interpolated version of the rectangular
C   grid A(i,j), where the colour-index used for particular array pixel
C   is given by:
C         Colour Index = NINT{[A(i,j)-BG/(FG-BG)]*FLOAT(NCOLS-1)} ,
C   with truncation at 0 and NCOLS-1, as necessary, according to the 
C   user-provided non-linear transform TRNL. The subrouitne TRNL will
C   be called with four arguments:
C         SUBROUTINE TRNL(XW,YW,XI,YJ)
C   and should return the real (fractional) array indicies XI and YJ
C   which correspond to the (real) world coordinates XW and YW. The 
C   appropriate Jacobian term should be factored into the array A(i,j)
C   before it is passed down (by pre-multilication). 
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
C    TRNL     R*4    I       -      Subroutine which performs the inverse
C                                   of the desired non-linear transform.
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
C     PGCXPX     Pixel-device support subroutine for PGCELX.
C     PGCXPS     PostScript support subroutine for PGCELX.
C     GREXEC     Dispatches command to appropriate device driver.
C     DSQINF     Inquires about viewport and window dimensions.
C
C History
C   D. S. Sivia      14 Feb 1995  Initial release.
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
C call the appropriate PGCELX support subrotuine.
C
      IF (LPS) THEN
        CALL PGCXPS(A,IDIM,JDIM,I1,I2,J1,J2,BG,TRNL,ASCALE,IR,IG,IB,
     *              NC,LCOLOR)
      ELSE
        NBUF=0
        LCHR=LEN(CHR)
        CALL GREXEC(GRGTYP,4,RBUF,NBUF,CHR,LCHR)
        IF (CHR(7:7).EQ.'P') THEN
          CALL PGCXPX(A,IDIM,JDIM,I1,I2,J1,J2,BG,TRNL,ASCALE,IC1,NC)
        ELSE
          WRITE(*,*)' Sorry, PGCELX does not support this device!'
        ENDIF
      ENDIF
C
C Reset the initial colour index.
C
      IF (.NOT. LPS) CALL PGSCI(ICSAVE)
      CALL PGEBUF
      END
C
      SUBROUTINE PGCXPX(A,IDIM,JDIM,I1,I2,J1,J2,BG,TRNL,ASCALE,IC1,NC)
C     ----------------------------------------------------------------
C
C Light-up the device pixels, with colours determined by a linearly
C interpolating array A.
C
      INCLUDE  'grpckg1.inc'
      REAL      A(IDIM,*),BUFFER(2050)
      CHARACTER CHR*16
      EXTERNAL  TRNL
C
      CALL DSQINF(XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *            YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC)
      IX1=NINT(XOFF)
      IX2=NINT(XOFF+XLEN)
      JY1=NINT(YOFF)
      JY2=NINT(YOFF+YLEN)
      DX=1.0/XSCALE
      DY=1.0/YSCALE
      XW0=(XOFF-XORG)*DX
      YW=(YOFF-YORG)*DY
      DO 20 JY=JY1,JY2
        XW=XW0
        BUFFER(2)=FLOAT(JY)
        NPIX=2
        DO 10 IX=IX1,IX2
          CALL TRNL(XW,YW,XI,YJ)
          I=INT(XI)
          J=INT(YJ)
          IF (I.GE.I1 .AND. I.LT.I2 .AND. J.GE.J1 .AND. J.LT.J2) THEN
            IF (NPIX.EQ.2) BUFFER(1)=FLOAT(IX)
            II=I+1
            JJ=J+1
            X=XI-FLOAT(I)
            Y=YJ-FLOAT(J)
            AXY=(1.0-X)*(A(I,J)+Y*(A(I,JJ)-A(I,J)))+
     *               X*(A(II,J)+Y*(A(II,JJ)-A(II,J)))
            K=NINT((AXY-BG)*ASCALE)
            IF (K.LT.0) K=0
            IF (K.GE.NC) K=NC-1
            NPIX=NPIX+1
            BUFFER(NPIX)=FLOAT(K+IC1)
          ELSE
            IF (NPIX.GT.2) THEN
              CALL GREXEC(GRGTYP,26,BUFFER,NPIX,CHR,LCHR)
              NPIX=2
            ENDIF
          ENDIF
          XW=XW+DX
  10    CONTINUE
        IF (NPIX.GT.2) CALL GREXEC(GRGTYP,26,BUFFER,NPIX,CHR,LCHR)
        YW=YW+DY
  20  CONTINUE
      END
C
      SUBROUTINE PGCXPS(A,IDIM,JDIM,I1,I2,J1,J2,BG,TRNL,ASCALE,IR,IG,IB,
     *                  NC,LCOLOR)
C     ------------------------------------------------------------------
C
C Postscript support subroutine for PGCELX.
C
      INCLUDE  'grpckg1.inc'
      REAL      A(IDIM,*)
      INTEGER   IR(0:*),IG(0:*),IB(0:*)
      INTEGER   IBUFER(2050)
      LOGICAL   LCOLOR
      EXTERNAL  TRNL
C
C Work out the nunmber of X and Y pixels for PS image, with NDOTS per 
C inch, and build an image transformation matrix.
C
      CALL DSQINF(XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *            YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC)
      NDOTS=100
C      IF (LCOLOR) NDOTS=50
      NXP=NINT(FLOAT(NDOTS)*XLEN/XPERIN)
      NYP=NINT(FLOAT(NDOTS)*YLEN/YPERIN)
      DX=XLEN/FLOAT(NXP)
      DY=YLEN/FLOAT(NYP)
      DXW=DX/XSCALE
      DYW=DY/YSCALE
      XW0=(XOFF-XORG)/XSCALE+DXW/2.0
      YW=(YOFF-YORG)/YSCALE+DYW/2.0
      AT=1.0/DX
      BT=0.0
      CT=0.0
      DT=1.0/DY
      TX0=-XOFF/DX
      TY=-YOFF/DY
C
      ASCALE=ASCALE*255.0/FLOAT(NC-1)
      YY=YOFF
      DO 20 JY=1,NYP
        XW=XW0
        XX=XOFF
        TX=TX0
        DXX=0.0
        NPIX=0
        DO 10 IX=1,NXP
          CALL TRNL(XW,YW,XI,YJ)
          I=INT(XI)
          J=INT(YJ)
          IF (I.GE.I1 .AND. I.LT.I2 .AND. J.GE.J1 .AND. J.LT.J2) THEN
            DXX=DXX+DX
            NPIX=NPIX+1
            II=I+1
            JJ=J+1
            X=XI-FLOAT(I)
            Y=YJ-FLOAT(J)
            AXY=(1.0-X)*(A(I,J)+Y*(A(I,JJ)-A(I,J)))+
     *               X*(A(II,J)+Y*(A(II,JJ)-A(II,J)))
            IC=NINT((AXY-BG)*ASCALE)
            IF (IC.LT.0) IC=0
            IF (IC.GE.NC) IC=NC-1
            IBUFER(NPIX)=FLOAT(IC)
          ELSE
            IF (NPIX.GT.0) THEN
              CALL PGCXP1(XX,DXX,YY,DY,NPIX,IBUFER,LCOLOR,IR,IG,IB,
     *                    AT,BT,CT,DT,TX,TY)
              TX=TX-FLOAT(NPIX)
              NPIX=0
            ENDIF
            TX=TX-1.0
            XX=XX+DXX+DX
            DXX=0.0
          ENDIF
          XW=XW+DXW
  10    CONTINUE
        IF (NPIX.GT.0) CALL PGCXP1(XX,DXX,YY,DY,NPIX,IBUFER,LCOLOR,
     *                             IR,IG,IB,AT,BT,CT,DT,TX,TY)
        YW=YW+DYW
        YY=YY+DY
        TY=TY-1.0
  20  CONTINUE
      CALL GRESC(' newpath ')
      CALL GRTERM
      END
C
      SUBROUTINE PGCXP1(XX,DXX,YY,DY,NPIX,IBUFER,LCOLOR,IR,IG,IB,
     *                  AT,BT,CT,DT,TX,TY)
C     -----------------------------------------------------------
C
C Postscript support subroutine for PGCXPS.
C
      INTEGER   IBUFER(*),IR(0:*),IG(0:*),IB(0:*)
      INTEGER   VALUE(33)
      LOGICAL   LCOLOR
      CHARACTER INLINE*80
C
C Set the path.
C
      WRITE (INLINE,100) NINT(XX),NINT(YY),NINT(DXX),NINT(DY),
     *          -NINT(DXX)
 100  FORMAT(I6,I6,' moveto ',I6, ' 0 rlineto  0 ',I6,' rlineto ',
     *I6,' 0 rlineto')
      CALL GRTERM
      CALL GRESC(' newpath ')
      CALL GRESC(INLINE)
      CALL GRESC(' closepath ')
C
C Use a PostScript "image" operator.
C
      WRITE (INLINE, '(A,I5,A)') '/picstr ',NPIX,' string def'
      CALL GRESC(INLINE)
      WRITE (INLINE,110) NPIX,AT,BT,CT,DT,TX,TY
 110  FORMAT(I4,' 1 8 [',6(1PE10.3,' '),']')
      CALL GRESC(INLINE)
      CALL GRESC('{ currentfile picstr readhexstring pop}')
      IF (LCOLOR) THEN
        CALL GRESC('  false 3 colorimage')
      ELSE
        CALL GRESC('  image')
      ENDIF
C
C Paint an array of pixels.
C
      L=0
      DO 10 I=1,NPIX
        L=L+1
        IC=IBUFER(I)
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
 120      FORMAT(33Z2.2)
          CALL GRESC(INLINE(1:66))
          L=0
        ENDIF
  10  CONTINUE
      IF (L.NE.0) THEN
        WRITE(INLINE,120) (VALUE(K),K=1,L)
        CALL GRESC(INLINE(1:2*L))
      ENDIF
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
