      integer mstarlabel
      parameter (mstarlabel=10000)

c Character string to represent the atom labels

      character*5 atomlabel(mstarlabel)
      character*3 atomlabellist(mcomp)

c Labels for bondlist

      character*11 bondlabel,bondlabellist(mstarlabel)
      character*17 bondanglelabel,bondanglelabellist(mstarlabel)
      character*256 first,readfile

c List of atomic numbers used in this file

      integer atomicnolist(mcomp),ilenatomtypelist

      integer ilenreadfile,ilenfirst,ilenatomlabel(mstarlabel)

      integer natomexpect,nbondexpect
     +,nbondlabel(mstarlabel),nbondlist
     +,nbondanglelabel(mstarlabel),nbondanglelist

c Reference integers for the atom label and atomic number

      integer iatomlabel(mstarlabel),iatomicno(mstarlabel)

c Reference integers for atoms

      integer ia1,ib1,ib2,atomtype(mstarlabel),ierr,iext
     +,iabondlist(mcomp,near),natomlist(mcomp),nchargelist(mcomp)


c Coordinates of individual atoms

      real*8 x,y,z,acharge(mstarlabel),chargelist(mcomp)

c List of bonds

      real bondlist(mstarlabel),bondanglelist(mstarlabel),alat,blat,clat,phib,thetac,phic

      logical unitcelldefined

      common/readmole_c/atomlabel,atomlabellist
     +,bondlabel,bondlabellist
     +,bondanglelabel,bondanglelabellist,first,readfile

      common/readmole_i/atomicnolist,ilenatomtypelist
     +,ilenreadfile,ilenfirst,ilenatomlabel
     +,natomexpect,nbondexpect
     +,nbondlabel,nbondlist
     +,nbondanglelabel,nbondanglelist
     +,iatomlabel,iatomicno
     +,ia1,ib1,ib2,atomtype,ierr,iext
     +,iabondlist,natomlist

      common/readmole_r/x,y,z,acharge,chargelist
     +,bondlist,bondanglelist,alat,blat,clat,phib,thetac,phic

      common/readmole_l/unitcelldefined

