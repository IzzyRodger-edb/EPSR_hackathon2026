c
c filenames and variables used within coord
c

      parameter(mdist=100,mz=1000)

      character*14 titler
      character*29 titledat(mdist)
      character*1 ext
      character*4 atom1(mdist),atom2(mdist)
      character*3 atomc,atomz,atomr
      logical dorhoz,dorhor

      real*4 damin(mdist),damax(mdist),coordn(mdist,2),rminc(mdist),rmaxc(mdist)
      real*4 zstep,zhalf,z(mz),gofz(mz,mdist),gzt(mz,mdist),gzt2(mz,mdist)

      integer*4 nz,nsumt,ltyp1(mdist),ltyp2(mdist),ltypc,ndist,ndist2

      common/coord_c/titler,titledat,ext,atom1,atom2,atomc,atomz,atomr

      common/coord_r/damin,damax,rminc,rmaxc,coordn,zstep,zhalf,z,gofz,gzt,gzt2
     *,fract2

      common/coord_i/nz,nsumt,ltyp1,ltyp2,ltypc,ndist,ndist2

      common/coord_l/dorhoz,dorhor

