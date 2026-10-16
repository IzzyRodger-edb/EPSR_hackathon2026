c
c filenames and variables common to SHARM and SDF
c
      character*2 naxis1,naxis2

      real*4 pi,piconv
     *,a1(3),b1(3),c1(3)
     *,pdif(3),perp1(3),alab(3),blab(3),clab(3)
     *,at(3),bt(3),ct(3),pt(3)
     *,phim1,thetam1,chim1,phim2,thetam2,chim2,thetallab,phillab,phispread,thetaspread
     *,phinew,thetanew,chinew

      integer*4 atomno(mcentre)

      common/SHARM_SDF_c/naxis1,naxis2

      common/SHARM_SDF_r/pi,piconv
     *,a1,b1,c1
     *,pdif,perp1,alab,blab,clab
     *,at,bt,ct,pt
     *,phim1,thetam1,chim1,phim2,thetam2,chim2,thetallab,phillab,phispread,thetaspread
     *,phinew,thetanew,chinew

      common/SHARM_SDF_i/atomno

      
