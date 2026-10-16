      integer*4 nq,naccess(mpairs)

      real*4 delfofq(mpts,mfile),datadif(mpts,mfile),chisqfac,ratfofq
     *,cprodsave(mpairs),nofr(mpz,mpairs)

      common/delta_chisqr/delfofq,datadif,chisqfac,ratfofq,cprodsave,nofr

      common/delta_chisqi/nq,naccess