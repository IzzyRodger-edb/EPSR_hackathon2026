c
c filenames and variables used within SHARM
c
      parameter(mcentre=20,mcentre2=mcentre*(mcentre-1)/2)

      character*14 titler
      character*1 ext
      character*2 nzt1,nzt2,nyt1,nyt2

      real vol,zstep,zhalf,zmax,zmaxd

      integer*4 nz,nsumt,ncoeffs
     *,ncent1,ncent2,ncz1,ncy1,ncz2,ncy2
     *,lztyp1(mcentre),lztyp2(mcentre),lytyp1(mcentre)
     *,lytyp2(mcentre),lctyp1(mcentre)
     *,lctyp2(mcentre),lctyp(mcentre)

c variables specific to # SDF

c Number of atoms at centre and type number of each atom

      integer*4 natofirst,npairfirst,typefirst(mcentre)
     *,natosecond,npairsecond,typesecond(mcentre)
     *,natomol,npairmol

c Array of distances between centre atoms in the order 
c   (1,2) (1,3) (1,4)...(2,3) (2,4)...(3,4) etc.

      real*4 distfirst(mcentre2),fractfirst
     *,distsecond(mcentre2),fractsecond,distcentre(mcentre2)
     *,fractcentre

c Temporary coordinates of atoms in a molecule

     *,xtemp(mcentre),ytemp(mcentre),ztemp(mcentre)


      common/SHARM_c/titler,ext,nzt1,nzt2,nyt1,nyt2


      common/SHARM_r/vol,zstep,zhalf,zmax,zmaxd
     *,distfirst,fractfirst,distsecond,fractsecond
     *,distcentre,fractcentre
     *,xtemp,ytemp,ztemp
    


      common/SHARM_i/nz,nsumt,ncoeffs
     *,ncent1,ncent2,ncz1,ncy1,ncz2,ncy2
     *,lztyp1,lztyp2,lytyp1,lytyp2,lctyp1,lctyp2,lctyp
     *,natofirst,npairfirst,typefirst
     *,natosecond,npairsecond,typesecond
     *,natomol,npairmol



