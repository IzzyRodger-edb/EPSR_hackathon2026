      parameter (mqinv=15)
      real*8 cwt(mfile,mpairs),ciwt(mfile,mpairs)
     *,cwtpot(mpairs,mpairs),cwtt(mpairs,mpairs+mfile)
     *,carray(mpairs+mfile),carrayold(mpairs+mfile,mpairs+mfile)
     *,cwttq(mqinv,mpairs,mpairs+mfile)
      real cwtav
      integer nxfile,xfile(mfile),nqinv,qinv(mqinv)
      common/weights_d/cwt,ciwt,cwtpot,cwtt,carray,carrayold,cwttq
      common/weights_r/cwtav
      common/weights_i/nxfile,xfile,nqinv,qinv
