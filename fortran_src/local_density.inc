c
c filenames and variables used within clusters
c

        parameter(mdist=mfile,mz=1000,mpix=201,maxncbonds=10)

c$$$    RH: let maximum number of bond definitions per cluster
c$$$    distribution = 10
        
        character*14 titler
        character*29 titledat(mdist)
        character*1 ext
        character*4 atom1(mdist,maxncbonds),atom2(mdist,maxncbonds)
        character*4 atom1s(mdist,maxncbonds),atom2s(mdist,maxncbonds)
        real*4 damin(mdist),damax(mdist)
        real*4 zstep,zhalf,z(mz),gofz(mz,mdist),gzt(mz,mdist),
     $       gzt2(mz,mdist)
        real rmins(mdist,maxncbonds)
        real rmaxs(mdist,maxncbonds)

c
c fraction of molecules in clusters of 2 or more molecules
c
        real*4 fract2(mdist),voidmax

        integer*4 nz,nsumt,ltyp(mdist,maxncbonds),
     $       ltyp1(mdist,maxncbonds),ndist,npix, numbdefs(mdist)
c$$$        RH: numbdefs() contains the number of bond definitions for
c$$$        each cluster distribution

        common/clusters_c/titler,titledat,ext,atom1,atom2,
     $       atom1s,atom2s

        common/clusters_r/zstep,zhalf,z,gofz,gzt,gzt2
     *,fract2,voidmax,rmins,rmaxs,damin,damax

        common/clusters_i/nz,nsumt,ltyp,ltyp1,ndist,npix,numbdefs

