c
c filenames and variables used within torangles
c

      parameter(mdist=mfile,mz=200)

      character*14 titler
      character*29 titledat(mdist)
      character*1 ext
      character*4 atom1(mdist),atom2(mdist),atom3(mdist)


      real*4 damin(mdist,2),damax(mdist,2)
      real*4 zstep,zhalf,z(mz),gofz(mz,mdist),gzt(mz,mdist),gzt2(mz,mdist)


      integer*4 nz,nsumt,ltyp(mdist),ltyp1(mdist),ltyp2(mdist,2),ltyp3(mdist,2)
     *,ndist,nfirst(mdist)

      common/torangles_c/titler,titledat,ext,atom1,atom2,atom3

      common/torangles_r/damin,damax,zstep,zhalf,z,gofz,gzt,gzt2

      common/torangles_i/nz,nsumt,ltyp,ltyp1,ltyp2,ltyp3,ndist,nfirst
