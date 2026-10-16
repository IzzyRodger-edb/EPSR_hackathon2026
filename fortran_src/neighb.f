      subroutine init_neighb_list(natom,nmol)
c
c initialise the neighbour list
c
      include 'epsr_dimension.inc'
      include 'neighbour_list.inc'

c
c number of terms in neighbour list
c
      nterms=natom*(natom-1)/2
c
c number of words required
c
      nwords=nterms/mlength
      nwordsref=nwords+1
c
c zero the neighbour list
c
      do i=0,nwords
         neighb(i)=0
      end do

c zero the sums of set bits for each word of neighb

      do i=0,nwordsref
         neighbref(i)=0
      end do

c initialise the index for the reference array

      nref=0

c
c now do the molecule neighbour list
c
c     nterms=nmol*(nmol-1)/2
c     nwords=nterms/mlength
c     do i=0,nwords
c        neighbmol(i)=0
c     end do
c
      return
      end

      subroutine set_neighb_list(i,j,natom)
c
c subroutine to set a bit in the neighbour list
c
      include 'epsr_dimension.inc'
      include 'neighbour_list.inc'
      integer*4 i,j,natom
c
c generate the bit reference and set the bit
c
      call get_bit_ref(i,j,natom)
      if(ibadd.ge.0) then
         neighb(iwref)=IBSet(neighb(iwref),ibref)
         neighbref(iwref+1)=neighbref(iwref+1)+1
      endif
      return
      end

c     subroutine set_neighbmol_list(i,j,nmol)
c
c subroutine to set a bit in the molecule neighbour list
c
c     include 'epsr_dimension.inc'
c     include 'neighbour_list.inc'
c     integer*4 i,j,nmol
c
c generate the bit reference and set the bit
c
c     call get_bit_ref(i,j,nmol)
c     if(ibadd.ge.0) neighbmol(iwref)=IBSet(neighbmol(iwref),ibref)
c     
c     return
c     end

      subroutine form_cummulative_sums()
      include 'epsr_dimension.inc'
      include 'neighbour_list.inc'
      i=0
      do while(i.lt.nwordsref)
         i=i+1
         neighbref(i)=neighbref(i)+neighbref(i-1)
      end do

      return
      end   

      function test_neighb_list(neighb_ref)
c
c function to test whether an atom pair occurs in the neighbour list
c
c if found, neighb_list contains the sequence number of the set bit
c otherwise it is set to zero
c
      include 'epsr_dimension.inc'
      include 'neighbour_list.inc'
      logical test_neighb_list
      integer*4 neighb_ref,ibitsum

      neighb_ref=0

c Test the bit

      test_neighb_list=BTest(neighb(iwref),ibref)

c If a valid pair is found get the reference number (i.e. the number of set bits
c up to this point in the neighbour list

      if(test_neighb_list) then
         neighb_ref=ibitsum()
      endif
      return
      end

c     function test_neighbmol_list(i,j,nmol)
c
c logical function to test whether an atom pair occurs in the neighbour list
c
c     include 'epsr_dimension.inc'
c     include 'neighbour_list.inc'
c     logical test_neighbmol_list
c     integer*4 i,j,nmol
c
c
c test the bit
c
c     test_neighbmol_list=BTest(neighbmol(iwref),ibref)
c
c     return
c     end

      subroutine get_bit_ref(i,j,natom)

      include 'epsr_dimension.inc'
      include 'neighbour_list.inc'
      integer*4 i,j,natom

c setup the bit reference integer
c subtract 1 from the bit reference since bit numbering starts at 0

      if(j.gt.i.and.j.le.natom) then
         ibadd=(2*natom-i)*(i-1)/2+j-i-1
      else if(i.gt.j.and.i.le.natom) then
         ibadd=(2*natom-j)*(j-1)/2+i-j-1
      else
         ibadd=-1
         return
      endif

c
c set up the word reference
c
      iwref=ibadd/mlength
c
c get the bit reference for this word
c
      ibref=ibadd-mlength*iwref

      return
      end

      subroutine make_neighbour_list(ia,jafirst,jalast,natom
     *,nnbref,nbref,nergref)

      include 'epsr_dimension.inc'
      include 'neighbour_list.inc'
      integer*4 ia,ja,jalast,natom,nnbref,nbref(*),nergref(*)

c Sets up the list of neighbours of atom ia and their corresponding
c energy table reference numbers

      nnbref=0

      ja=jafirst

      call get_bit_ref(ia,jafirst,natom)

      do while(ja.le.jalast)

c Test the current bit to see if it is set

         if(BTest(neighb(iwref),ibref)) then

c bit is set. Therefore return the number of set bits up to this
c point in the neighbour table: this is used to reference the energy table

            nnbref=nnbref+1
            nbref(nnbref)=ja
            nergref(nnbref)=ibitsum()

         endif

c Go to the next atom.

         ja=ja+1
         if(ja.lt.ia) then
               ibadd=ibadd+natom-ja
         else if(ja.gt.ia) then
               ibadd=ibadd+1
         else
c atoms are the same - definitely no neighbour reference for this pair!
               ibadd=ibadd+natom-ja+1
               ja=ja+1
         endif

c set up the word reference

         iwref=ibadd/mlength

c get the bit reference for this word

         ibref=ibadd-mlength*iwref

      end do

      return
      end

      function ibitsum()
c
c Gets the sequence number of the set bits in neighb for a specified i,j
c
      include 'epsr_dimension.inc'
      include 'neighbour_list.inc'
      integer*4 ibitsum

c search forwards or backwards depending on value of ib to save time

c     write(6,*) ibref,iwref,neighbref(iwref),neighbref(iwref+1)
      if(ibref.lt.mlength/2) then

         ib=-1
         ibitsum=neighbref(iwref)
         do while(ib.lt.ibref)
            ib=ib+1     
            ibitsum=ibitsum+IBits(neighb(iwref),ib,1)
         end do

      else

         ib=mlength
         ibitsum=neighbref(iwref+1)

c When going backwards don't include the bit corresponding to ibref since
c it must be included in the sum!

         do while(ib.gt.ibref+1)
            ib=ib-1
            ibitsum=ibitsum-IBits(neighb(iwref),ib,1)
         end do

      endif

      return
      end   
