c
c parameters for accessing the neighbour list.
c This is stored as a bit array to save space. Each word of the array contains 32 bits, 
c and each bit refers to a particular atom pair. A bit set to 1 means this pair is 
c included in the neighbour list
c
c needs to be preceded by EPSR_dimension.inc
c
      parameter(mlength=32,mterms=(matom*(matom-1)/(2*mlength))
     *,mtermsmol=(mmol*(mmol-1)/(2*mlength)))


c neighb is top triangular matrix which stores a bit map. Each bit indicates whether or not
c a pair of atoms is in each other's neighbour list. Likewise for the neighbour list of molecules.
c neighbref will represent the number of set bits in each of the rows 
c of the top triangular matrix, neighb. 

      integer*4 neighb(0:mterms)
c     *,neighbmol(0:mtermsmol)
      integer*4 neighbref(0:mterms+1)

c ibref and iwref store internal references to the bit number and word number of neighb
c corresponding to a specified pair i,j

      integer*4 ibref,iwref,nwordsref,ibadd

      common/neighbour_list_i1/neighb
c     *,neighbmol

      common/neighbour_list_i2/neighbref

      common/neighbour_list_i4/ibref,iwref,nwordsref,ibadd

