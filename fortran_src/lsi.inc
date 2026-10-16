c
c filenames and variables used within chains
c

      parameter(mdist=30,mz=200,mchain=200,mbond=20,mring=25)

      character*14 titler
      character*29 titledat(mdist)
      character*1 ext
      character*4 atom1(mdist,mdist),atom2(mdist,mdist),atom3(mdist,mdist)

      real dff1(mdist,mdist),dff2(mdist,mdist),dfs1(mdist,mdist),dfs2(mdist,mdist),chainlen(mdist,2)
      real dft1(mdist,mdist),dft2(mdist,mdist)
      real zstep,zhalf,z(mz),gofz(mz,mdist),gzt(mz,mdist),gzt2(mz,mdist)
      real ang213min(mdist,mdist),ang213max(mdist,mdist)
      real lsimin,lsimax,rtrimax

      integer nz,nsumt,lftype(mdist),lstype(mdist),lttype(mdist),ndist,ndist3

c nneighb(im) stores the number of neighbours of molecule im
c neighb(im,*) lists the molecule numbers of the neighbours of molecule im

      integer neighb(mmol,mbond),nneighb(mmol),ia1min(mmol,mbond),ia2min(mmol,mbond),iwmin(mmol,mbond)

c List of distances to all neighbours. This will be in order of increasing distance

      real dneighb(mmol,mbond)

c Whether molecule in list is bonded or not, in order of increasing distance.

      logical bneighb(mmol,mbond)

c Stores the number criteria for defining bonds for each distribution

      integer nwordsave(mdist)


      common/lsi_c/titler,titledat,ext,atom1,atom2,atom3

      common/lsi_r/dff1,dff2,dfs1,dfs2,dft1,dft2
     *,chainlen,zstep,zhalf,z,gofz,gzt,gzt2,ang213min,ang213max
     *,lsimin,lsimax,rtrimax,dneighb
 
      common/lsi_i/nz,nsumt,lftype,lstype,lttype,ndist,ndist3
     *,neighb,nneighb,ia1min,ia2min,iwmin,nwordsave


