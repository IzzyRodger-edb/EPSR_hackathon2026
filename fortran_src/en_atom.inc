c
c common variables - temporary coordinates
c     used en_atom_intra and en_atom_inter
c
!Centre of molecule positions
      real xaddold,yaddold,zaddold,ergold(matom),ergnew(matom)
      integer it,jt,imold,jmold

      integer iergref,nnbref,nbref(matom)
!Number and atom indices of intra-molecular atoms used to change energy
      integer nnatemp,natemp(near)

      common/en_atom_r/xaddold,yaddold,zaddold,ergold,ergnew

      common/en_atom_i/it,jt,imold,jmold
     *,iergref,nnbref,nbref,nnatemp,natemp

