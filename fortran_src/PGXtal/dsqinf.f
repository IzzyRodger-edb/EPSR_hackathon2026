CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE DSQINF(XOFF,XLEN,XORG,XSCALE,XPERIN,XBLC,XTRC,
     *                  YOFF,YLEN,YORG,YSCALE,YPERIN,YBLC,YTRC)
C     ---------------------------------------------------------
C
      DATA BIG,SMALL /1.0E+20,1.0E-20/
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C Purpose
C      Obtain some viewport and window information about the current 
C    PGPLOT device, without directly accessing the common blocks in
C    pgplot.inc.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    XPERIN   R*4    O       -      Plot X scale in dots/inch.
C    YPERIN   R*4    O       -      Plot Y scale in dots/inch.
C    XOFF     R*4    O       -      Absolute coord of blc of viewport.
C    YOFF     R*4    O       -      Absolute coord of blc of viewport.
C    XLEN     R*4    O       -      Width of viewport in absolute coord.
C    YLEN     R*4    O       -      Height of viewport in absolute coord.
C    XORG     R*4    O       -      Absolute coord of world X=0.
C    YORG     R*4    O       -      Absolute coord of world Y=0.
C    XSCALE   R*4    O       -      Absolute units per world coord in X.
C    YSCALE   R*4    O       -      Absolute units per world coord in Y.
C    XBLC     R*4    O       -      World X coord at blc of window.
C    XTRC     R*4    O       -      World X coord at trc of window.
C    YBLC     R*4    O       -      World Y coord at blc of window.
C    YTRC     R*4    O       -      World Y coord at trc of window.
C
C Globals
C     None.
C
C External Calls
C   SUBROUTINE   DESCRIPTION
C     PGQVP      Inquires about viewport dimensions.
C     PGQWIN     Inquires about world coords of window.
C
C History
C   D. S. Sivia       1 Aug 1996  Initial release.
C-----------------------------------------------------------------------
C
      CALL PGQWIN(XBLC,XTRC,YBLC,YTRC)
      CALL PGQVP(1,XI1,XI2,YI1,YI2)
      CALL PGQVP(3,XOFF,XP2,YOFF,YP2)
      XLEN=ABS(XP2-XOFF)
      YLEN=ABS(YP2-YOFF)
      XPERIN=XLEN/(ABS(XI2-XI1)+SMALL)
      YPERIN=YLEN/(ABS(YI2-YI1)+SMALL)
      XWDIF=XTRC-XBLC
      YWDIF=YTRC-YBLC
      AXWDIF=BIG
      AYWDIF=BIG
      IF (ABS(XWDIF).GT.SMALL) AXWDIF=1.0/XWDIF
      IF (ABS(YWDIF).GT.SMALL) AYWDIF=1.0/YWDIF
      XSCALE=XLEN*AXWDIF
      YSCALE=YLEN*AYWDIF
      XORG=(XOFF*XTRC-XP2*XBLC)*AXWDIF
      YORG=(YOFF*YTRC-YP2*YBLC)*AYWDIF
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
