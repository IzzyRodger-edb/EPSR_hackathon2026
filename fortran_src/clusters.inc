c
c filenames and variables used within clusters
c

      parameter(mdist=mfile,mz=5000,mpix=201)

      character*14 titler
      character*29 titledat(mdist)
      character*1 ext
      character*4 atom1(mdist),atom2(mdist)

      real*4 damin(mdist),damax(mdist)
      real*4 zstep,zhalf,z(mz),gofz(mz,mdist),gzt(mz,mdist),gzt2(mz,mdist)
c
c fraction of molecules in clusters of 2 or more molecules
c
      real*4 fract2(mdist),voidmax


      integer*4 nz,nsumt,ltyp(mdist),ltyp1(mdist),ndist,npix

      common/clusters_c/titler,titledat,ext,atom1,atom2

      common/clusters_r/damin,damax,zstep,zhalf,z,gofz,gzt,gzt2
     *,fract2,voidmax

      common/clusters_i/nz,nsumt,ltyp,ltyp1,ndist,npix

