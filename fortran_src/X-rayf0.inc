c
c Store the X-ray form factors for each atom type
c

      integer mxfile
      parameter (mxfile=2)
      real*8 f0X(mpts,mcomp,mxfile),fm0x(mpts,mcomp,mxfile)
     *,weightx(mpts,mpairs,mxfile),backx(mpts,mxfile)
      integer datapoint(mxfile),datanorm(mxfile),normxx(mxfile),iafformaffx(mxfile)

c     common/Xrayf0_c/

      common/Xrayf0_r/f0x,fm0x,weightx,backx
 
      common/Xrayf0_i/datapoint,datanorm,normxx,iafformaffx

