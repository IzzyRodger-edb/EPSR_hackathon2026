c
c labels of molecules in current file
c
      character*3 moltype(mcomp)
      character*256 moltemplate(mcomp)

c The following integers will store the number of attempted moves and the
c number of rejected moves. These numbers are used to assess the step size
c for subsequent moves. The whole molecule translations and rotations
c numbers are kept separate for different molecule types to ensure that roughly the same
c fraction of attempted moves for different molecules are accepted. Used in molfit.

      integer nmol,nato(mmol),ngrot(mmol),naf(mmol),nal(mmol),nrf(mmol),nrl(mmol)
     *,mtype(mmol),nmoltype
     *,nmovei(2,mcomp),ninter(2,mcomp),nwhole(2,mcomp),nmovem(2,mcomp)
     *,nmoveit(2),nintert(2),nwholet(2),nmovemt(2)
     *,nstepreject(mmol),tetheredatom(mmol),moltetheredatom(mcomp)
     
      logical tethered(mmol),moltethered(mcomp)

c The following store the different rotation and translation steps for the different
c molecule types

      real xadd(mmol),yadd(mmol),zadd(mmol),phix(mmol)
     *,phiy(mmol),phiz(mmol),dmax,volume,amasst(mmol),moldiam(mmol)
     *,amasstsqrt(mmol)
     *,interntrans(mcomp),internrot(mcomp),wholerot(mcomp),wholetrans(mcomp)

      common/clusterc/moltype,moltemplate

      common/clusteri/nmol,nato,ngrot,naf,nal,nrf,nrl
     *,mtype,nmoltype
     *,nmovei,ninter,nwhole,nmovem
     *,nmoveit,nintert,nwholet,nmovemt
     *,nstepreject,tetheredatom,moltetheredatom
     
      common/clusterl/tethered,moltethered

      common/clusterr/xadd,yadd,zadd,phix
     *,phiy,phiz,dmax,volume,amasst,moldiam
     *,amasstsqrt
     *,interntrans,internrot,wholerot,wholetrans
