c
c filenames and variables used within partials
c

      parameter(mdist=mpairs,mz=1000)

      character*14 titler
      character*29 titledat(mdist)
      character*1 ext
      character*4 atom1(mdist),atom2(mdist)
      character*3 atomc,atomz,atomr
      logical dorhoz,dorhor

      real*4 zstep,zhalf,z(mz),gofz(mz,mdist),gzt(mz,mdist),gzt2(mz,mdist),rnormalise

      integer*4 nz,nsumt,ltyp1(mdist),ltyp2(mdist),ltypc,ndist

      common/partials_c/titler,titledat,ext,atom1,atom2,atomc,atomz,atomr

      common/partials_r/zstep,zhalf,z,gofz,gzt,gzt2,rnormalise
 
      common/partials_i/nz,nsumt,ltyp1,ltyp2,ltypc,ndist
      
      common/partials_l/dorhoz,dorhor

