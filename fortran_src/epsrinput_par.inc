c
c This include file includes parameters needed in EPSR16/17
c which have been read in by subroutine INITA and will
c be written out in subroutine NORMA
c
      integer NISUM,NESUM,NESUMT,nsum,nsums
      real roverlap,rminfac,RMINPT,RMAXPT,widthptc,RBROAD,EXPECF,expecfpwr,RCHARGE
     *,ENSUM,ENSUM2,PRSUM,PRSUM2,lorchqwidth,lorchqwidthq
     *,rminbragg,rmaxgr

       common/epsrinput_par/nisum,nesum,nesumt,nsum,nsums
     *,roverlap,rminfac,rminpt,rmaxpt,widthptc,rbroad,expecf,expecfpwr,rcharge
     *,ENSUM,ENSUM2,PRSUM,PRSUM2,lorchqwidth,lorchqwidthq
     *,rminbragg,rmaxgr
