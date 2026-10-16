c
c filenames and variables used within triangles
c

      parameter(mdist=mfile,mz=200)

      character*14 titler
      character*29 titledat(mdist)
      character*1 ext
      character*4 atom1(mdist),atom2(mdist),atom3(mdist)
      character*3 atomc,atomz,atomr
      logical dorhoz,dorhor


      real damin(mdist,2),damax(mdist,2),rminc(mdist),rmaxc(mdist)
      real zstep,zhalf,z(mz),gofz(mz,mdist),gzt(mz,mdist),gzt2(mz,mdist)


      integer nz,nsumt,ltyp(mdist),ltyp1(mdist),ltyp2(mdist),ltyp3(mdist),ltypc,ndist

      common/triangles_c/titler,titledat,ext,atom1,atom2,atom3,atomc,atomz,atomr

      common/triangles_r/damin,damax,rminc,rmaxc,zstep,zhalf,z,gofz,gzt,gzt2

      common/triangles_i/nz,nsumt,ltyp,ltyp1,ltyp2,ltyp3,ltypc,ndist

      common/triangles_l/dorhoz,dorhor

