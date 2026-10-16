c
c filenames and variables used within chains
c

      parameter(mdist=20,mz=200,mchain=200,mbond=30,mring=25)

      character*14 titler
      character*29 titledat(mdist)
      character*1 ext
      character*4 atom1(mdist,mdist),atom2(mdist,mdist),atom3(mdist,mdist)

      real dff1(mdist,mdist),dff2(mdist,mdist),dfs1(mdist,mdist),dfs2(mdist,mdist),chainlen(mdist,2)
      real dft1(mdist,mdist),dft2(mdist,mdist)
      real zstep,zhalf,z(mz),gofz(mz,mdist),gzt(mz,mdist),gzt2(mz,mdist)
      real ang213min(mdist,mdist),ang213max(mdist,mdist)

      integer nz,nsumt,lftype(mdist),lstype(mdist),lttype(mdist),ndist,ndist3

c nneighb(im) stores the number of neighbours of molecule im
c neighb(im,*) lists the molecule numbers of the neighbours of molecule im

      integer neighb(mmol,mbond),nneighb(mmol)

c icount is a stack counter to be used for the current chain or ring
c nrlab(mdist) is a list of the molecules in sequence for the current chain or ring 
c rneighb(mdist) stores the neighbour number of the current trial neighbour molecule
c - when this exceeds the number of neighbours of molecule icount, this particular molecule
c is abandoned and a new one is sought.

      integer icount,nrlab(0:mz),rneighb(mz)

c spvalue(mmol,mring) is the shortest topological distances (shortest paths)
c of each molecule to the one at the begining at successive places along the chain or ring.
c spvaluet(mmol) is a temporary array of shortest paths
c spseq(*) is the list of shortest topological distances in the chain or ring. For a shortest path
c chain or ring this sequence should be 'unimodal'
c nspvaluet will carry the number of times the shortest path occurs between this molecule and that at the origin

      integer spvalue(mmol,mring),spvaluet(mmol)
     *,spseq(0:mchain),sptest,spgrad,spmin,spmax,rptest

c Variables for setting up the shortest path to the molecule at the origin

      integer rneighb1(mmol),rneighb2(mmol)

c Stores the number criteria for defining bonds for each distribution

      integer nwordsave(mdist)

c carryon will be used to decide whether to proceed down the chain or ring or not

      logical nodd,sptrue1,sptrue2,carryon(mmol)



      common/chains_c/titler,titledat,ext,atom1,atom2,atom3

      common/chains_r/dff1,dff2,dfs1,dfs2,dft1,dft2
     *,chainlen,zstep,zhalf,z,gofz,gzt,gzt2,ang213min,ang213max
 
      common/chains_i/nz,nsumt,lftype,lstype,lttype,ndist,ndist3
     *,neighb,nneighb,icount,nrlab,rneighb,spvalue,spvaluet,spseq
     *,sptest,spgrad,spmin,spmax,rptest,rneighb1,rneighb2,nwordsave

      common/chains_l/nodd,sptrue1,sptrue2,carryon

