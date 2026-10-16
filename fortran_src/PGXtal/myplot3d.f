CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C***<plot>**************************************************************
C
      SUBROUTINE PLOT3d(nin,dens,n1,n2,n3,size,iwidth,icoll,side,nspher
     *,rdcoor,thcoor,phcoor,radius,rgbbal,xtr,ytr,csize,titl,dsurf)
C     ---------------------------------------------------------------
C
      CHARACTER*(*) TITL
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
!      INCLUDE  'COLTABS.INC'
      REAL dens(0:n1,0:n2,0:n3),RGB(3,3),EYE(3),LIGHT(3),LATICE(3,4)
     *,centre(3),rgbobj(3),rgbbak(3),rdcoor(*),thcoor(*),phcoor(*)
     *,rgbbal(3,*),radius(*)
 	integer ifcol(50),ilcol(50)
      CHARACTER STRING*32,TYPE*16,CHR*16
      LOGICAL   OVRLAY,LSHIN
      DATA      EYE,LIGHT /0.0,0.0,1000.0,-1.0,-1.0,-1.0/
      DATA      RGB /0.0,0.0,1.0,0.35,0.35,0.35,1.0,1.0,1.0/
c        CALL FMXMN(dens,n1,n2,n3,DHIGH,DLOW)
	write(6,487)
487	format(1x,' Type 3 RGB Colour indices for background')
	read(nin,*) rgbbak(1),rgbbak(2),rgbbak(3)
	write(6,485)
485	format(1x,' Type ishade (1-8): 0 means no shading')
	read(nin,*) ishade
	write(6,489)
489	format(1x,' Type 3 RGB Colour indices for object')
	read(nin,*) rgbobj(1),rgbobj(2),rgbobj(3)
	write(6,488)
488	format(1x,' Type 3 (x,y,z) coordinates for light source')
	read(nin,*) xl,yl,zl
	write(6,486)
486	format(1x,' Type fade factor (0 = no fading, 1=full fading)')
	read(nin,*) fadefc
	write(6,484)
484	format(1x,' Type transparency of object (0=0%,1=25%,2=50%,3=75%)')
	read(nin,*) itrans
c these are converted to coordinate system of plotting frame, and reversed
c to show the direction of the light
c
	light(1)=-yl
	light(2)=-zl
	light(3)=-xl
	write(6,490)
490	format(1x,' Type difuse, shine, polish and contrast')
	read(nin,*) difus,shin,polish,alfa
	lshin=.false.
        CALL EULER(side,LATICE,sina,cosa,sinb,cosb)
        write(6,*)
        do j=1,4
           write(6,492) (latice(i,j),i=1,3)
492	format(3(1x,e13.6))
        end do
        CALL PGQCOL(ICMIN,ICMAX)
        CALL PGSCH(SIZE)
        CALL PGSLW(IWIDTH)
	xwmin=-side
	xwmax=side
	ywmin=-side
	ywmax=side
        CALL PGWINDOW(xwmin,xwmax,ywmin,ywmax)
        CALL PGSCI(icoll)
c        CALL PGBOX('BC',0.0,0,'BC',0.0,0)
        CALL PGQINF('TYPE',TYPE,LCHR)
c        IF (TYPE.EQ.'PS' .OR. TYPE.EQ.'VPS' .OR. TYPE.EQ.'CPS' .OR.
c     *      TYPE.EQ.'VCPS') THEN
c	        CALL PGSCI(0)
c       ELSE
c	        CALL PGSCI(1)
c        ENDIF
        CALL SBFINT(RGBbak,16,1,1,MAXBUF)
	lastc=16
	nclbal=16
	ictest=(icmax-16)/(nspher+2)
	if(ictest.lt.nclbal) nclbal=ictest
	if(iabs(ishade).gt.0) then
	  ifbak=lastc+1
	  ilbak=lastc+nclbal
	  lastc=ilbak
	  CALL COLINT(RGBbak,ifbak,ilbak,DIFUS,SHIN,POLISH)
	  if(ishade.gt.0) call sbfbkg(ifbak,ilbak,ishade)
	  if(ishade.lt.0) call sbfbkg(ilbak,ifbak,iabs(ishade))
	endif
	if(nspher.gt.0) then
	  do is=1,nspher
		zo=rdcoor(is)*sin(thcoor(is))*cos(phcoor(is))
		xo=rdcoor(is)*sin(thcoor(is))*sin(phcoor(is))
		yo=rdcoor(is)*cos(thcoor(is))
		CALL ROTY(xo,yo,zo,U,V,W,SINA,COSA)
		CALL ROTX(U,V,W,centre(1),centre(2),centre(3)
     *,SINB,COSB)
		centre(3)=centre(3)-2.0*side
		ifcol(is)=lastc+1
		ilcol(is)=lastc+nclbal
	        CALL COLINT(RGBbal(1,is),ifcol(is),ilcol(is),DIFUS,SHIN,POLISH)
		call sbball(eye,centre,radius(is)
     *,ifcol(is),ilcol(is),light,lshin,xo,yo,zo)
		lastc=ilcol(is)
	  end do
	  icobj=lastc+1
	else
	  icobj=17
	end if
	call colint(rgbobj,icobj,icmax,difus,shin,polish)
	if(itrans.gt.0.and.itrans.lt.4) then
	      CALL SBTSUR(EYE,LATICE,dens,n1,n2,n3,dsurf,icobj,icmax
     *,light,lshin,itrans,fadefc)
	else
	      CALL SBSURF(EYE,LATICE,dens,n1,n2,n3,dsurf,icobj,icmax
     *,light,lshin,fadefc)
	endif
!	write(6,*) 'myplot3d line 114 '
        CALL SBFCLS(1)
c	xst1=x(2)-x(1)
c	nxst1=int(x(1)/xst1)
c	xmin=nxst1*xst1
c	yst1=y(2)-y(1)
c	nyst1=int(y(1)/yst1)
c	ymin=nyst1*yst1
c	write(6,599) xmin,ymin
599	format(2(1x,e13.6))
        CALL AXES3D(EYE,LATICE,side,size)
c	xti=xwmin+xtr*(xwmax-xwmin)
c	yti=ywmin+ytr*(ywmax-ywmin)
        CALL PGSCH(cSIZE)
        CALL PGmTEXT('T',ytr,xtr,0.0,TITL)
	return
      END
C
C***<3d-surface>********************************************************
C
      SUBROUTINE FMXMN(Y,n1,n2,n3,YMAX,YMIN)
C     --------------------------------------
C
      REAL         Y(0:n1,0:n2,0:n3)
      CHARACTER*32 STRING
C
      YMIN1=y(0,0,0)
      YMAX1=y(0,0,0)
	do k=0,n3
		do j=0,n2
			do i=0,n1
        			if(y(i,j,k).lt.ymin1) ymin1=y(i,j,k)
				if(y(i,j,k).gt.ymax1) ymax1=y(i,j,k)
 			end do
		end do
	end do
1     YMIN=YMIN1
      YMAX=YMAX1
      WRITE(*,100) YMIN,YMAX
 100  FORMAT(' Surface>  Zmin & Zmax for plot ? (def= ',e10.3,
     *       ',',E10.3,')  : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=1) YMIN,YMAX
      IF (YMAX.LE.YMIN) GOTO 1
	if(ymin.lt.0.0) ymin=0.0
      END
C
      SUBROUTINE EULER(side,LATICE,sina,cosa,sinb,cosb)
C     ------------------------
C
      CHARACTER*80 STRING
      REAL         LATICE(3,4),IA,IB
      DATA         PIRAD /0.01745329252/
C
   1  WRITE(*,100)
 100  FORMAT(' Surface> Rotation and tilt ?  (def=0,0 deg) : ',$)
      CALL FORMTQ(STRING,80,NNN)
      IF (NNN.EQ.0) THEN
        SINA=0
        COSA=1
        SINB=0
        COSB=1
      ELSE
        READ(STRING,*,ERR=1) IA,IB
        SINA=-SIN(IA*PIRAD)
        COSA=COS(IA*PIRAD)
        SINB=-SIN(IB*PIRAD)
        COSB=COS(IB*PIRAD)
      ENDIF
c initial coordinates of centre of object
	xo=-side
	yo=-side
	zo=-side
c initial coordinates of a-axis
	xa=-side
	ya=-side
	za=side
c initial coordinates of b-axis
	xb=side
	yb=-side
	zb=-side
c initial coordinates of c-axis
	xc=-side
	yc=side
	zc=-side
c perform rotation about X and Y
      CALL ROTY(xo,yo,zo,U,V,W,SINA,COSA)
      CALL ROTX(U,V,W,LATICE(1,1),LATICE(2,1),latice(3,1),SINB,COSB)
      CALL ROTY(xa,ya,za,U,V,W,SINA,COSA)
      CALL ROTX(U,V,W,LATICE(1,2),LATICE(2,2),latice(3,2),SINB,COSB)
      CALL ROTY(xb,yb,zb,U,V,W,SINA,COSA)
      CALL ROTX(U,V,W,LATICE(1,3),LATICE(2,3),latice(3,3),SINB,COSB)
      CALL ROTY(xc,yc,zc,U,V,W,SINA,COSA)
      CALL ROTX(U,V,W,LATICE(1,4),LATICE(2,4),latice(3,4),SINB,COSB)
	zshift=2*side
	do i=1,4
c
c shift z coordinates to behind the screen
c
		latice(3,i)=latice(3,i)-zshift
	end do
      END
C
      SUBROUTINE ROTX(X,Y,Z,U,V,W,S,C)
C     --------------------------------
C
      U=X
      V=Y*C+Z*S
      W=-Y*S+Z*C
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
      SUBROUTINE SRFCOL(RGB,NCBAND,ICTAB,DIF,SHIN,POLISH,LSHIN,IC1,IC2)
C     -----------------------------------------------------------------
C
      REAL         RGB(*)
      LOGICAL      LSHIN
      CHARACTER*32 STRING
C
      POLISH=1.0
   1  ICTAB=0
      WRITE(*,100) IC1,IC2,ICTAB
 100  FORMAT(' Surface> Colour table ?  (',I1,'-',I1,',def=',I1,
     *       ') : ',$)
      CALL FORMTQ(STRING,32,NNN)
      IF (NNN.NE.0) READ(STRING,*,ERR=1) ICTAB
      IF (ICTAB.GT.IC2) THEN
        GOTO 1
      ELSEIF (ICTAB.LT.IC1) THEN
        NCBAND=1
   2    RGB(1)=0.0
        RGB(2)=1.0
        RGB(3)=0.0
        WRITE(*,110) (RGB(I),I=1,3)
 110    FORMAT(' Surface> RGB colour ?  (def=',F3.1,',',F3.1,',',
     *    F3.1,') : ',$)
        CALL FORMTQ(STRING,32,NNN)
        IF (NNN.NE.0) READ(STRING,*,ERR=2) (RGB(I),I=1,3)
        IF (RGB(1)+RGB(2)+RGB(3).LE.0.05) GOTO 2 
      ELSE
   3    NCBAND=8
        WRITE(*,120) NCBAND
 120    FORMAT(' Surface> No. of colour-bands ?  (def=',I1,') : ',$)
        CALL FORMTQ(STRING,32,NNN)
        IF (NNN.NE.0) READ(STRING,*,ERR=3) NCBAND
        NCBAND=MAX(MIN(NCBAND,64),1)
      ENDIF
      CALL LOGQYN(' Surface> A shiny gloss ?','N',LSHIN)
      IF (LSHIN) THEN
        SHIN=1.0
        DIF=0.0
      ELSE
        SHIN=0.0
   4    DIF=0.7
        WRITE(*,130) DIF
 130    FORMAT(' Surface> Diffusiveness ?  (def=',F3.1,') : ',$)
        CALL FORMTQ(STRING,32,NNN)
        IF (NNN.NE.0) READ(STRING,*,ERR=4) DIF
        DIFUSE=MAX(MIN(DIF,1.0),0.1)
      ENDIF
      END
C
      SUBROUTINE AXES3D(EYE,LATICE,side,csize)
C     ----------------------------------------------------------------
C
      REAL          EYE(*),LATICE(3,*)
      REAL          pivx(3),pivy(3),pivz(3),LATCAB(3,4),orx(3,2)
     *,ory(3,2),orz(3,2)
      CHARACTER*10 XLBL,YLBL,zlbl
      AX=LATICE(1,2)-LATICE(1,1)
      AY=LATICE(2,2)-LATICE(2,1)
      AZ=LATICE(3,2)-LATICE(3,1)
      BX=LATICE(1,3)-LATICE(1,1)
      BY=LATICE(2,3)-LATICE(2,1)
      BZ=LATICE(3,3)-LATICE(3,1)
      cX=LATICE(1,4)-LATICE(1,1)
      cY=LATICE(2,4)-LATICE(2,1)
      cZ=LATICE(3,4)-LATICE(3,1)
	latcab(1,1)=0.0
	latcab(2,1)=0.0
	latcab(3,1)=-2.0*side
      LATCAB(1,2)=latcab(1,1)+0.2*aX
      LATCAB(2,2)=latcab(2,1)+0.2*aY
      LATCAB(3,2)=latcab(3,1)+0.2*aZ
      LATCAB(1,3)=latcab(1,1)+0.2*BX
      LATCAB(2,3)=latcab(2,1)+0.2*BY
      LATCAB(3,3)=latcab(3,1)+0.2*BZ
      LATCAB(1,4)=latcab(1,1)+0.2*cX
      LATCAB(2,4)=latcab(2,1)+0.2*cY
      LATCAB(3,4)=latcab(3,1)+0.2*cZ
c
c draw axes
c
	call pgsls(2)
      CALL SBLINE(EYE,LATcab(1,1),LATcab(1,2),1,.FALSE.)
      CALL SBLINE(EYE,LATcab(1,1),LATcab(1,3),1,.FALSE.)
      CALL SBLINE(EYE,LATcab(1,1),LATcab(1,4),1,.FALSE.)
	call pgsls(1)
c
c draw labels at end of each axis
c
	xlbl='\fix'
	ylbl='\fiy'
	zlbl='\fiz'
c
c orientation of xlbl
c
c	orient(1,1)=ax
c	orient(2,1)=ay
c	orient(3,1)=az
c	orient(1,1)=bx
c	orient(2,1)=by
c	orient(3,1)=bz
c	call sbtext(eye,xlbl,1,latcab(1,2),0,orient,csize)
c
c orientation of xlbl
c
c	orient(1,1)=ax
c	orient(2,1)=ay
c	orient(3,1)=az
c	orient(1,1)=bx
c	orient(2,1)=by
c	orient(3,1)=bz
c	call sbtext(eye,ylbl,1,latcab(1,3),0,orient,csize)
c
c orientation of ylbl
c
c	orient(1,1)=cx
c	orient(2,1)=cy
c	orient(3,1)=cz
c	orient(1,1)=bx
c	orient(2,1)=by
c	orient(3,1)=bz
c	call sbtext(eye,zlbl,1,latcab(1,4),0,orient,csize)
	scla=csize*0.2
      XSIGN=+1.0
      XSCL=-SCLA
      IF (ax.lT.0.0) THEN
        XSIGN=-1.0
      ENDIF
      ORX(1,1)=1.0
      ORX(2,1)=0.0
      ORX(3,1)=0.0
      ORX(1,2)=cX
      ORX(2,2)=cY
      ORX(3,2)=cZ
      PIVX(1)=LATcab(1,2)+0.02*ax
      PIVX(2)=LATcab(2,2)+0.02*ay
      PIVX(3)=LATcab(3,2)+0.02*az
      CALL SBTEXT(EYE,XLBL,1,PIVX,0.5,ORX,SCLA)
      YSIGN=+1.0
      IF (bx.lt.0.0) THEN
        YSIGN=-1.0
      ENDIF
      ORY(1,1)=1.0
      ORY(2,1)=0.0
      ORY(3,1)=0.0
      ORY(1,2)=cX
      ORY(2,2)=cY
      ORY(3,2)=cZ
      PIVY(1)=LATcab(1,3)+0.02*bx
      PIVY(2)=LATcab(2,3)+0.02*by
      PIVY(3)=LATcab(3,3)+0.02*bz
      CALL SBTEXT(EYE,YLBL,1,PIVY,0.5,ORY,SCLA)
      ORz(1,1)=1.0
      ORz(2,1)=0.0
      ORz(3,1)=0.0
      ORz(1,2)=cX
      ORz(2,2)=cY
      ORz(3,2)=cZ
      PIVz(1)=LATcab(1,4)+0.02*cx
      PIVz(2)=LATcab(2,4)+0.02*cy
      PIVz(3)=LATcab(3,4)+0.02*cz
      CALL SBTEXT(EYE,zLBL,1,PIVz,0.5,ORz,SCLA)
      END
C
      SUBROUTINE AXNUMS(EYE,XMIN,XMAX,PIVX,ORX,SCLA,XSIGN)
C     ----------------------------------------------------
C
      REAL      EYE(*),PIVX(*),ORX(3,*)
      REAL      END1(3),END2(3),PIVOT(3)
      CHARACTER NLBL*20
      DATA      FRTICK,FRNUM /0.02,0.10/
C
      XR=PGRND(XMAX-XMIN,NSUB)
      DX=XR/FLOAT(NSUB)
      IF (DX.LE.1.0E-20) RETURN
   1  XJ=DX*FLOAT(1+INT(XMIN/DX))
      IF ((XJ+DX).GE.XMAX) THEN
        DX=DX/2.0
        NSUB=NSUB*2
        GOTO 1
      ENDIF
      IF (XMIN.LT.0.0) XJ=XJ-DX
      XN=XSIGN/(XMAX-XMIN)
      XH=0.5*(XMIN+XMAX)
      DO 20 J=1,NSUB
        IF (XJ.GT.XMAX) RETURN
        XF=XN*(XJ-XH)
        DO 10 I=1,3
          END1(I)=PIVX(I)+XF*ORX(I,1)+SCLA*ORX(I,2)
          END2(I)=END1(I)-FRTICK*ORX(I,2)
          PIVOT(I)=END1(I)-FRNUM*ORX(I,2)
  10    CONTINUE
        IPOWER=INT(LOG10(ABS(XJ)+1.0E-10))-5
        IF (XJ.LT.1.0) IPOWER=IPOWER-1
        X=XJ/(10.0**IPOWER)
        IMANTS=NINT(X)
        CALL PGNUMB(IMANTS,IPOWER,0,NLBL,NC)
        CALL SBLINE(EYE,END1,END2,1,.FALSE.)
        CALL SBTEXT(EYE,NLBL,1,PIVOT,0.5,ORX,SCLA*0.15)
        XJ=XJ+DX
  20  CONTINUE
      END
C
      SUBROUTINE AZNUMS(EYE,ZMIN,ZMAX,LATZ,SCLA)
C     ------------------------------------------
C
      REAL      EYE(*),LATZ(3,*)
      REAL      END1(3),END2(3),PIVOT(3),ORIENT(3,3)
      CHARACTER NLBL*20
      DATA      FRTICK,FRNUM /-0.07,-0.10/
C
      ORIENT(1,1)=SCLA
      ORIENT(2,1)=0.0
      ORIENT(3,1)=0.0
      ZSIGN=+1.0
      IF ((LATZ(2,2)-LATZ(2,1)).LT.0.0) ZSIGN=-1.0
      DO 10 I=1,3
        ORIENT(I,3)=LATZ(I,2)-LATZ(I,1)
        ORIENT(I,2)=ZSIGN*SCLA*ORIENT(I,3)
  10  CONTINUE
      ZR=PGRND(ZMAX-ZMIN,NSUB)
      DZ=ZR/FLOAT(NSUB)
      IF (DZ.LE.1.0E-20) RETURN
   1  ZJ=DZ*FLOAT(1+INT(ZMIN/DZ))
      IF ((ZJ+DZ).GE.ZMAX) THEN
        DZ=DZ/2.0
        NSUB=NSUB*2
        GOTO 1
      ENDIF
      IF (ZMIN.LE.0.0) ZJ=ZJ-DZ
      ZN=1.0/(ZMAX-ZMIN)
      DO 30 J=1,NSUB
        IF (ZJ.GT.ZMAX) RETURN
        ZF=ZN*(ZJ-ZMIN)
        DO 20 I=1,3
          END1(I)=LATZ(I,1)+ZF*ORIENT(I,3)
          END2(I)=END1(I)+FRTICK*ORIENT(I,1)
          PIVOT(I)=END1(I)+FRNUM*ORIENT(I,1)-FRTICK*ORIENT(I,2)
  20    CONTINUE
        IPOWER=INT(LOG10(ABS(ZJ)+1.0E-10))-5
        IF (ZJ.LT.1.0) IPOWER=IPOWER-1
        Z=ZJ/(10.0**IPOWER)
        IMANTS=NINT(Z)
        CALL PGNUMB(IMANTS,IPOWER,0,NLBL,NC)
        CALL SBLINE(EYE,END1,END2,1,.FALSE.)
        CALL SBTEXT(EYE,NLBL,1,PIVOT,1.0,ORIENT,SCLA*0.12)
        ZJ=ZJ+DZ
  30  CONTINUE
      END
C
C***<contour>***********************************************************
C
      SUBROUTINE CONTOR(MAP,NX,NY,X,Y,N1,N2,SIZE,IWIDTH,OVRLAY)
C     ---------------------------------------------------------
C
      REAL    MAP(*),X(*),Y(*)
      REAL    CONT(25,2),TR(6)
      INTEGER NCONT(2),LCOLOR(2),LSTYLE(2),LWIDTH(2)
      LOGICAL OVRLAY
      DATA    LCOLOR,LSTYLE,LWIDTH,NWIND,IDP /1,3,1,1,1,2,1,0/
C
      CALL INIT(X,Y,NX,NY,TR,MAP,DMIN,DMAX,N1*N2)
      IF (OVRLAY) THEN
        CALL AUTCNT(NCONT,CONT,DMIN,DMAX,OVRLAY)
        CALL PGSLW(1)
        CALL PGCONT(MAP,N1,N2,1,NX,1,NY,CONT,NCONT,TR)
        CALL PGSLW(IWIDTH)
        RETURN
      ENDIF
      CALL ASKCNT(NCONT,CONT,DMIN,DMAX)
      CALL PGBEGIN(0,'?',1,1)
      CALL PGPAPER(0.0,1.0)
      CALL PGSCH(SIZE)
      CALL PGSLW(IWIDTH)
      CALL PGENV(X(1),X(NX),Y(1),Y(NY),0,0)
      DO 20 I=1,2
        CALL PGSCI(LCOLOR(I))
        CALL PGSLS(LSTYLE(I))
        LWDTH=2*LWIDTH(I)
        IF (LWDTH.GT.7) LWDTH=7
        CALL PGSLW(LWDTH)
        CALL PGCONT(MAP,N1,N2,1,NX,1,NY,CONT(1,I),NCONT(I),TR)
  20  CONTINUE
      CALL PGSCI(1)
      CALL PGSLW(IWIDTH)
      END
C
      SUBROUTINE INIT(X,Y,NX,NY,TR,MAP,DMIN,DMAX,NMAP)
C     ------------------------------------------------
C
      REAL X(*),Y(*),TR(*),MAP(*)
C
      TR(1)=X(1)-(X(NX)-X(1))/FLOAT(NX-1)
      TR(2)=(X(NX)-X(1))/FLOAT(NX-1)
      TR(3)=0.0
      TR(4)=Y(1)-(Y(NY)-Y(1))/FLOAT(NY-1)
      TR(5)=0.0
      TR(6)=(Y(NY)-Y(1))/FLOAT(NY-1)
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
      READ(10,200,ERR=1) STRING
 200  FORMAT(A)
      CALL FINDNC(STRING,132,NCONT(1))
      READ(STRING,*,ERR=1) (CONT(I,1),I=1,NCONT(1))
      WRITE(*,*)
   2  WRITE(*,110)
 110  FORMAT(' Thick>  ',$)
      READ(10,200,ERR=2) STRING
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
      SUBROUTINE GREY(MAP,NX,NY,X,Y,N1,N2,IPLOT,SIZE,IWIDTH)
C     ------------------------------------------------------
C
      REAL MAP(*),X(*),Y(*)
      REAL TR(6),TRCOL(6),COLBAR(2,256),LUTUSR(3,256)
      REAL RED(256),GR(256),BL(256)
      DATA TRCOL /0.0,1.0,0.0,0.0,0.0,1.0/
C
      NCOLS=256
      CALL INIT(X,Y,NX,NY,TR,MAP,DMIN,DMAX,N1*N2)
      CALL STGREY(MAP,N1*N2,DMIN,DMAX,COLBAR,NCOLS)
      TRCOL(6)=(DMAX-DMIN)/255.0
      TRCOL(4)=DMIN-TRCOL(6)
      IF (IPLOT.EQ.7) CALL LUTIN(LUTUSR,NCOLS,NCLUSR,IPLOT)
      CALL PGBEGIN(0,'?',1,1)
      CALL PGPAPER(0.0,1.0)
      CALL PGSCH(SIZE)
      CALL PGSLW(IWIDTH)
      CALL PGVPORT(0.86,0.90,0.2,0.90)
      CALL PGWINDOW(1.0,2.0,DMIN,DMAX)
      CALL SETCOL(IPLOT,NCOLS,LUTUSR,RED,GR,BL)
      CALL PGSCH(0.50*SIZE)
      CALL PGBOX('B',0.0,0,'B',0.0,0)
      CALL PGBOX('C',0.0,0,'CMTSVI',0.0,0)
      CALL PGCELL(COLBAR,2,256,1,2,1,256,1.0,0.0,TRCOL,NCOLS,RED,GR,BL)
      CALL PGBOX('B',0.0,0,'B',0.0,0)
      CALL PGBOX('C',0.0,0,'CMTSVI',0.0,0)
      CALL PGSCH(SIZE)
      CALL PGVPORT(0.12,0.82,0.2,0.90)
      XMIN=TR(1)+TR(2)
      XMAX=TR(1)+FLOAT(NX)*TR(2)
      YMIN=TR(4)+TR(6)
      YMAX=TR(4)+FLOAT(NY)*TR(6)
      CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
      CALL PGBOX('C',0.0,0,'C',0.0,0)
      CALL PGBOX('BNSTI',0.0,0,'BNSTI',0.0,0)
      CALL PGCELL(MAP,N1,N2,1,NX,1,NY,1.0,0.0,TR,NCOLS,RED,GR,BL)
      CALL PGBOX('C',0.0,0,'C',0.0,0)
      CALL PGBOX('BNSTI',0.0,0,'BNSTI',0.0,0)
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
      read(10,200,ERR=1) FILNAM
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
      READ(10,200) STRING
 200  FORMAT(A)
      NNN=0
      DO 10 I=1,NCHAR
  10    IF (STRING(I:I).NE.' ') NNN=NNN+1
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
