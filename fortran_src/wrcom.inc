      parameter (nwrmax=100,nform=(nwrmax*2+1))
      real*4 xo(10000),yo(10000,nwrmax),eo(10000,nwrmax)    
      character*2 ext
      character*29 titleo(nwrmax)
      common /wrcom/xo,yo,eo
      common /wrcomc/ext,titleo
