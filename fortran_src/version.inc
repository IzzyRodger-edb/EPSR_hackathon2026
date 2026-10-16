c
c Stores the current version of EPSR
c
      character*20 date,version
      integer*4 ilendate,ilenversion
      parameter (date='26: 2017-03-10',ilendate=14)

      common/version_c/version
      common/version_i/ilenversion
