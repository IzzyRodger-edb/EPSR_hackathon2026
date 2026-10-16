      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call run_clusters()

      end      

      subroutine run_clusters()
c
c programme to estimate number of clusters in a distribution
c of molecules. Atoms are defined to be in a cluster when they are within a
c specified distance range of a neighbour 

c 2007-01-25 Modified and simplified to correct for a small counting error when estimating 
c which molecules belong to which cluster.

c 2007-01-26 Modified to estimate the percentage of clusters which are percolating, 
c i.e. spanning the box in one or more dimensions.

      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'directory_name.inc'
      include 'system_commands.inc'
      include 'symbols.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'clusters.inc'
      integer*4 ifound1,ifound2,ierrout
      character*20 text
c
c setup the clusters symbols
c
      call setup_symbols_clusters()
c
c now setup the input file
c
      nin=11
      nout=6
      call setup_input_file(nin,nout,ierrout)
c
c if an error was encountered on setup, exit now
c
      if(ierrout.ne.0) return

c
c read the input file
c
      call read_input_file(nin,nout,ierrout)
c
c if an error was encountered on read, exit now
c
      if(ierrout.ne.0) return
c
c set up the number of cluster size bins
c
      iref=2
      read(symbol_value(iref),*) nz
      if(nz.gt.mz) nz=mz
c
c set up current accumulator
c
      iref=iref+1
      read(symbol_value(iref),*) nsumt
c
c check the number of distributions does not exceed the array size
c
      ndist=nfiles
      if(ndist.gt.mdist) ndist=mdist
c
c setup the type number of the first and second atoms used to define
c whether two molecules are bonded, and the maximum and minimum distances
c
      do id=1,ndist
            it=0
            ifound1=0
            ifound2=0
            iref1=ireff(nsymboluptodata+1,id)
            iref2=iref1+1
            iref3=iref2+1
            iref4=iref3+1
            do while(it.lt.ntype.and.(ifound1.eq.0.or.ifound2.eq.0))
               it=it+1
               if(type(it).eq.symbol_value(iref1)) ifound1=it
               if(type(it).eq.symbol_value(iref2)) ifound2=it
            end do
            ltyp(id)=ifound1
            ltyp1(id)=ifound2
            atom1(id)=type(ifound1)
            atom2(id)=type(ifound2)
            read(symbol_value(iref3),*) damin(id)
            read(symbol_value(iref4),*) damax(id)
            write(6,*) id,ltyp(id),ltyp1(id),damin(id),damax(id)
      end do
      write(6,*) nz,ndist
c
c setup step size and z values
c
      zstep=0.5/float(nz)
      zhalf=0.5*zstep
      do j=1,nz
            z(j)=float(j-1)
      end do

      call initgofz(nz,ndist,gofz,mz,mdist)

c
c normalise atomic coordinates to box size
c
      call normalise_coordinates()
c
c set up the table of molecular bonds
c
      call find_clusters()
c
c axis title
c
      titler='Cluster_size'
      do id=1,ndist
c
c set the column title for this pair of atoms
c
            ilent1=index(atom1(id),' ')-1
            ilent2=index(atom2(id),' ')-1
            titledat(id)='   '//atom1(id)(1:ilent1)//'-'//atom2(id)(1:ilent2)
      end do

c output the data

      ext='n'

      call outgofz(fnameinp,ext,nz,ndist,nsumt,z,gofz,gzt,gzt2
     *,titler,titledat,mz,mdist)

c
c For each distribution calculate the fraction of molecules in clusters of 2 or more
c
      do id=1,ndist
            sum1=0.0
            sum2=0.0
            do i=1,nz
                  sum1=sum1+z(i)*gzt(i,id)
                  if(z(i).gt.1.0) sum2=sum2+z(i)*gzt(i,id)
            end do
            fract2(id)=sum2/sum1
      end do

c
c write input file
c

c update the current value of nsumt and ndist

      icount=0
      text='nsize'
      call search_for_symbol(text,icount)
      iref=ireff(icount,1)
      write(symbol_value(iref),*) nz
      iref=iref+1
      write(symbol_value(iref),*) nsumt
      iref=iref+1
      write(symbol_value(iref),*) ndist

      call save_files_others(14)
      
      do id=1,ndist
c
c for each cluster distribution indicate the mean number of molecules per cluster
c
            write(14,101) id,fract2(id)
      end do
101   format(/'Cluster ',i3,': fraction containing 2 or more molecules'
     *,f10.5)

      close(14)

      return
      end

      subroutine find_clusters()
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'clusters.inc'
      integer*4 neighb(mmol,mdist),mbond(mmol),labcluster(mmol)
     1,lencluster(mmol),ncllab(mmol)
      nfdim=mdist
      zstep=0.5/float(nz)
      pi=3.141592653
c
c look at each requested pair of atoms
c
      do id=1,ndist
            write(6,*) ltyp(id),ltyp1(id),damin(id),damax(id)
c
c set up for each molecule a list of bonded neighbours.
c
c mbond keeps track of how many molecules each molecule is bonded to. -1 
c indicates that no specified atom pairs were found for this molecule, so that 
c it is not included in the cluster analysis.
c
      do i=1,nmol
            mbond(i)=-1
      end do
c
c step through all atoms 
c
      do im1=1,nmol-1
         do ia1=naf(im1),nal(im1)
c
c now look at all atoms in other (distinct) molecules and determine which 
c ones are bonded
c
            do im2=im1+1,nmol
               do ia2=naf(im2),nal(im2)
c
c find out if this pair of atoms match the specified types
c
                  if(
     *(itype(ia1).eq.ltyp(id).and.itype(ia2).eq.ltyp1(id))
     *.or.
     *(itype(ia2).eq.ltyp(id).and.itype(ia1).eq.ltyp1(id))
     *) then
c
c O.K. so we have found a pair of molecules with atoms of the correct labels
c
c zero the bond number accumulator if this has not already been set
c
                     if(mbond(im1).lt.0) mbond(im1)=0
                     if(mbond(im2).lt.0) mbond(im2)=0
c
c calculate atom1 to atom2 distance, applying periodic boundary conditions
c
                     xc1=xadd(im1)
                     yc1=yadd(im1)
                     zc1=zadd(im1)
                     xa1=xa(ia1)+xc1
                     ya1=ya(ia1)+yc1
                     za1=za(ia1)+zc1
                     xc2=xadd(im2)
                     yc2=yadd(im2)
                     zc2=zadd(im2)
                     xa2=xa(ia2)+xc2
                     ya2=ya(ia2)+yc2
                     za2=za(ia2)+zc2      
                     xdif=xa2-xa1
                     ydif=ya2-ya1
                     zdif=za2-za1
c
c apply periodic boundary conditions
c
                     call minimagec(xdif,ydif,zdif)
                     dsq12=xdif*xdif+ydif*ydif+zdif*zdif
                     da12=sqrt(dsq12)*dmax
c
c check that distance is within bounds prescribed for this distribution
c
                     if(da12.ge.damin(id).and.da12.le.damax(id))
     * then
c
c O.K. found a bond - save molecule numbers
c
                        mbond(im1)=mbond(im1)+1
                        mbond(im2)=mbond(im2)+1
                        neighb(im1,mbond(im1))=im2
                        neighb(im2,mbond(im2))=im1
c     write(6,*) im1,im2,mbond(im1),mbond(im2)
                     endif    
                  endif
               end do
            end do
         end do
c        write(6,*) im1,mbond(im1),(neighb(im1,im2),im2=1,mbond(im1))
      end do
c
c Generate initial cluster labels - ncllab(im1). Each valid molecule is put in its own 
c unique cluster. If a molecule is not valid (i.e. does not satisfy the type requirements
c so that mbond = -1) it is given a fictitious cluster label nmol+1. Such a cluster label
c can never exist because there cannot be more clusters than the number of molecules in 
c the box.

      ncluster=0
      do im1=1,nmol
         if(mbond(im1).ge.0) then
            ncluster=ncluster+1
            ncllab(im1)=ncluster
            labcluster(ncluster)=ncluster
         else
            ncllab(im1)=nmol+1
         endif
      end do
      nclusterold=0
c
c begin loop to look for clusters
c
      do while (nclusterold.ne.ncluster)

c save the previous number of clusters

         nclusterold=ncluster

c step through molecules looking at the bonds to other molecules

         do im1=1,nmol

c Ignore molecules of the wrong type

            if(mbond(im1).ge.0) then
               nseek=ncllab(im1)
               in=0
               do while(in.lt.mbond(im1))

c This molecule is bonded to at least one neighbour

c Find out the lowest cluster label (if any) of bonded neighbours

                  in=in+1
                  im2=neighb(im1,in)
                  if(ncllab(im2).lt.nseek) nseek=ncllab(im2)
               end do

c Set all these bonded atoms into the same cluster

               ncllab(im1)=nseek
               in=0
               do while(in.lt.mbond(im1))
                  in=in+1
                  im2=neighb(im1,in)
                  ncllab(im2)=nseek
               end do
            endif
         end do

c Now count how many distinct clusters there are

         ncluster=0

         do im1=1,nmol
            if(mbond(im1).ge.0) then

c Run through the existing clusters and check whether the cluster this molecule belongs to
c already exists.

               iflag=0
               ic=0
               do while(ic.lt.ncluster.and.iflag.eq.0) 

                  ic=ic+1
                  if(ncllab(im1).eq.labcluster(ic)) iflag=ic

               end do

c If this cluster does not exist, we start a new one.

               if(iflag.eq.0) then

                  ncluster=ncluster+1
                  labcluster(ncluster)=ncllab(im1)

               endif
            endif
         end do

         write(6,*) ncluster,nclusterold
                  
      end do

c now work out the length of each cluster

      write(6,*) ncluster
      do ic=1,ncluster
         sum=0
         do im=1,nmol
            if(ncllab(im).eq.labcluster(ic)) sum=sum+1
         end do
         lencluster(ic)=sum
c         write(6,*) ic,labcluster(ic),lencluster(ic)
      end do

c calculate the average size of each cluster

      lensum=0
      do ic=1,ncluster
            lensum=lensum+lencluster(ic)
      end do
      if(ncluster.gt.0) then
            avlen=real(lensum)/real(ncluster)
      else  
            avlen=0.0
      endif
      write(6,995) id,ncluster,avlen
995   format(/'For cluster ',i2,' there are ',i5
     *,' clusters of average size = ',e12.5,' molecules')

c Look to see which clusters are percolating, i.e. span the simulation box in at least
c one dimension.

      npercolate=0
      do ic=1,ncluster

c First find a molecule that is in this cluster

         im1=0
         iflag=0
         do while(im1.lt.nmol.and.iflag.eq.0)

            im1=im1+1
            if(ncllab(im1).eq.labcluster(ic)) iflag=im1

         end do

c Record the maximum separation of all the other molecules in this cluster from this 
c starting molecule, applying periodic boundary conditions

         xmin=0.0
         ymin=0.0
         zmin=0.0
         xmax=0.0
         ymax=0.0
         zmax=0.0
         do im2=1,nmol
            if(ncllab(im2).eq.labcluster(ic)) then

               xdif=xadd(im2)-xadd(im1)
               ydif=yadd(im2)-yadd(im1)
               zdif=zadd(im2)-zadd(im1)
               call minimagec(xdif,ydif,zdif)

c Store the largest and smallest separations

               if(xdif.gt.xmax) xmax=xdif
               if(xdif.lt.xmin) xmin=xdif
               if(ydif.gt.ymax) ymax=ydif
               if(ydif.lt.ymin) ymin=ydif
               if(zdif.gt.zmax) zmax=zdif
               if(zdif.lt.zmin) zmin=zdif

            endif
         end do
c        write(6,*) ic,xmin,xmax,ymin,ymax,zmin,zmax

c Apply periodic boundary conditions to the difference between maximum and minimum 
c coordinates

         xdif=xmax-xmin
         ydif=ymax-ymin
         zdif=zmax-zmin
         diftest=1.0-damax(id)/dmax

c Test for overlap. If any one of these differences is less than the maximum separation
c for this cluster, it is regarded as percolating.

c        write(6,*) ic,xdif,ydif,zdif,diftest

         if(xdif.ge.diftest.or.ydif.ge.diftest.or.zdif.ge.diftest) then

c Cluster is percolating.

            npercolate=npercolate+lencluster(ic)

         endif

      end do

c Fraction of molecules in percolating clusters:-

      fracpercolate=real(npercolate)/real(lensum)

      write(6,994) id,fracpercolate
994   format(/'Clusters> For cluster ',i2
     1,' fraction of molecules which percolate is ',f10.5)

c form histogram of number of bonded neighbours

c form histogram of cluster sizes

      do ic=1,ncluster
            iref=lencluster(ic)+1
            if(iref.le.nz) then
                  gofz(iref,id)=gofz(iref,id)+1
            else
                  gofz(nz,id)=gofz(nz,id)+1
            endif
      end do

c normalise the distributions to the number of clusters

      do ic=2,nz
            gofz(ic,id)=gofz(ic,id)/real(ncluster)
      end do

c save the percolation fraction

      gofz(1,id)=fracpercolate
c
c repeat for other distributions
c
      end do
      return
      end



      
