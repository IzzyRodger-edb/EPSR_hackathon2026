CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C***<plot>**************************************************************
C
      SUBROUTINE PLT2DZ(ZMIN,ZMAX,XMIN,XMAX,YMIN,YMAX,SIZE,IWIDTH,XLBL,
     *                  YLBL,TITL)
C     -----------------------------------------------------------------
C
      CHARACTER*(*) XLBL,YLBL,TITL
      CHARACTER*32  STRING
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      This subroutine plots a 2D colour map given the users own
C    (peculiar) mapping. To this end, the following routine must be 
C    supplied by the user:
C        SUBROUTINE TRNL(XW,YW,ZCOLOR)
C    should return the colour-scale, between 0.0 (for ZMIN) and 1.0
C    (for ZMAX), which corresponds to the (real) world coordinates XW 
C    and YW.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    ZMIN     R*4    I       -      Lowest "intensity" for the plot.
C    ZMAX     R*4    I       -      Highest "intensity" for the plot.
C    XMIN     R*4    I       -      Xmin for plot, in output units.
C    XMAX     R*4    I       -      Xmax for plot, in output units.
C    YMIN     R*4    I       -      Ymin for plot, in output units.
C    YMIN     R*4    I       -      Ymax for plot, in output units.
C    SIZE     R*4    I       -      Character-size for plot.
C    IWIDTH   I*4    I       -      Line-width for plot.
C    XLBL     A*1    I     *(*)     Label for X-axis.
C    YLBL     A*1    I     *(*)     Label for Y-axis.
C    TITL     A*1    I     *(*)     Title for plot.
C
C Globals
C    COLTABS.INC
C
C History
C   D. S. Sivia      15 Feb 1995  Initial release.
C   D. S. Sivia      21 Feb 1995  Put in option to over-lay contours.
C   D. S. Sivia       8 Aug 1995  Modified PLOT2DBX to PLOT2DBZ.
C   D. S. Sivia      28 Sep 1995  Fixed "bug" with the contrast-factor.
C   D. S. Sivia      27 Aug 1997  No longer use non-standard Q-format.
C   D. S. Sivia      15 Sep 1997  Fortran made LINUX-friendly!
C-----------------------------------------------------------------------
C
      WRITE(*,*)
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
      IF (IPLOT.GE.2 .AND. IPLOT.LE.7) THEN
        CALL GREY(ZMIN,ZMAX,XMIN,XMAX,YMIN,YMAX,IPLOT,SIZE,IWIDTH)
      ELSE
        GOTO 1
      ENDIF
      CALL PGLABEL(XLBL,YLBL,TITL)
      CALL PGEND
      END
C
C***<grey scale>********************************************************
C
      SUBROUTINE GREY(DMIN,DMAX,XMIN,XMAX,YMIN,YMAX,IPLOT,SIZE,IWIDTH)
C     ----------------------------------------------------------------
C
      REAL     TRCOL(6),COLBAR(2,256),LUTUSR(3,256)
      REAL     RED(256),GR(256),BL(256)
      EXTERNAL TRNL
      DATA     TRCOL /0.0,1.0,0.0,0.0,0.0,1.0/
C
      NCOLS=256
      CALL STGREY(DMIN,DMAX,COLBAR,NCOLS,CNTRST)
      TRCOL(6)=(DMAX-DMIN)/255.0
      TRCOL(4)=DMIN-TRCOL(6)
      IF (IPLOT.EQ.7) CALL LUTIN(LUTUSR,NCOLS,NCLUSR,IPLOT)
      CALL PGBEGIN(0,'?',1,1)
      CALL PGPAPER(0.0,1.0)
      CALL PGSCH(SIZE)
      CALL PGSLW(IWIDTH)
      CALL PGVPORT(0.86,0.90,0.2,0.90)
      CALL PGWINDOW(1.0,2.0,DMIN,DMAX)
      CALL SETCOL(IPLOT,NCOLS,LUTUSR,RED,GR,BL,CNTRST)
      CALL PGSCH(0.50*SIZE)
      CALL PGBOX('B',0.0,0,'B',0.0,0)
      CALL PGBOX('C',0.0,0,'CMTSVI',0.0,0)
      CALL PGCELL(COLBAR,2,256,1,2,1,256,1.0,0.0,TRCOL,NCOLS,RED,GR,BL)
      CALL PGBOX('B',0.0,0,'B',0.0,0)
      CALL PGBOX('C',0.0,0,'CMTSVI',0.0,0)
      CALL PGSCH(SIZE)
      CALL PGVPORT(0.12,0.82,0.2,0.90)
      CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
      CALL PGBOX('C',0.0,0,'C',0.0,0)
      CALL PGBOX('BNSTI',0.0,0,'BNSTI',0.0,0)
      CALL PGCELX(TRNL,NCOLS,RED,GR,BL)
      CALL PGBOX('C',0.0,0,'C',0.0,0)
      CALL PGBOX('BNSTI',0.0,0,'BNSTI',0.0,0)
      END
C
      SUBROUTINE STGREY(FMIN,FMAX,COLBAR,NCOLS,CNTRST)
C     ------------------------------------------------
C
      REAL         COLBAR(2,*)
      CHARACTER*32 STRING
C
      WRITE(*,*)
   2  CNTRST=1.0
      WRITE(*,110)
 110  FORMAT(' >>  Contrast factor ?  (def=1.0)  : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=2) CNTRST
      CNTRST=MAX(MIN(CNTRST,10.0),0.01)
      DCOL=0.999/FLOAT(NCOLS-1)
      COL=0.0
      IF (FMAX.LT.FMIN) THEN
        COL=1.0
        DCOL=-DCOL
      ENDIF
      DO 20 I=1,NCOLS
        COLBAR(1,I)=COL
        COLBAR(2,I)=COL
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
      SUBROUTINE SETCOL(IPLOT,NCOLS,LUTUSR,R,G,B,C)
C     ---------------------------------------------
C
      INCLUDE 'COLTABS.INC'
      REAL     LUTUSR(3,*),R(*),G(*),B(*)
C
      IF (IPLOT.EQ.2) THEN
        Z=0.0
        DZ=0.999/FLOAT(NCOLS-1)
        DO 10 I=1,NCOLS
          ZC=Z**C
          R(I)=ZC
          G(I)=ZC
          B(I)=ZC
          Z=Z+DZ
  10    CONTINUE
      ELSEIF (IPLOT.EQ.3) THEN
        CALL STCOL1(HEAT,R,G,B,NCOLS,C)
      ELSEIF (IPLOT.EQ.4) THEN
        CALL STCOL1(SPECTRUM,R,G,B,NCOLS,C)
      ELSEIF (IPLOT.EQ.5) THEN
        CALL STCOL1(BGYRW,R,G,B,NCOLS,C)
      ELSEIF (IPLOT.EQ.6) THEN
        CALL STCOL1(SERP,R,G,B,NCOLS,C)
      ELSE
        CALL STCOL1(LUTUSR,R,G,B,NCOLS,C)
      ENDIF
      END
C
      SUBROUTINE STCOL1(LUT,R,G,B,N,C)
C     --------------------------------
C
      REAL LUT(3,*),R(*),G(*),B(*)
C
      DO 10 I=1,N
        R(I)=LUT(1,I)**C
        G(I)=LUT(2,I)**C
        B(I)=LUT(3,I)**C
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
C
      SUBROUTINE PGCELX(TRNL,NCOLS,R,G,B)
C     -----------------------------------
C
      REAL         R(0:NCOLS-1),G(0:NCOLS-1),B(0:NCOLS-1)
      INTEGER      NCOLS
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
C      This subroutine is a non-linear version of PGCELL modified
C   specifically for PLOT2DBX2! 
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
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
C   D. S. Sivia       8 Aug 1995  Modified for PLOT2DBX2!
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
      ASCALE=FLOAT(NC-1)
C
C Check to see whether a pixel device or Postscript is being used, and 
C call the appropriate PGCELX support subrotuine.
C
      IF (LPS) THEN
        CALL PGCXPS(TRNL,ASCALE,IR,IG,IB,NC,LCOLOR)
      ELSE
        NBUF=0
        LCHR=LEN(CHR)
        CALL GREXEC(GRGTYP,4,RBUF,NBUF,CHR,LCHR)
        IF (CHR(7:7).EQ.'P') THEN
          CALL PGCXPX(TRNL,ASCALE,IC1,NC)
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
      SUBROUTINE PGCXPX(TRNL,ASCALE,IC1,NC)
C     -------------------------------------
C
C Light-up the device pixels, with colours determined by a linearly
C interpolating array A.
C
      INCLUDE  'grpckg1.inc'
      REAL      BUFFER(2050)
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
          CALL TRNL(XW,YW,XCOLOR)
          IF (XCOLOR.GE.0.0) THEN
            IF (NPIX.EQ.2) BUFFER(1)=FLOAT(IX)
            K=NINT(XCOLOR*ASCALE)
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
      SUBROUTINE PGCXPS(TRNL,ASCALE,IR,IG,IB,NC,LCOLOR)
C     -------------------------------------------------
C
C Postscript support subroutine for PGCELX.
C
      INCLUDE  'grpckg1.inc'
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
          CALL TRNL(XW,YW,XCOLOR)
          IF (XCOLOR.GE.0.0) THEN
            DXX=DXX+DX
            NPIX=NPIX+1
            IC=NINT(XCOLOR*ASCALE)
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
