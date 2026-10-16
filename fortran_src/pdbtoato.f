      program pdbtoato
c
c this is a modified version of minitc to convert a cubic pdb file to EPSR ato format
c
      include 'epsr_dimension.inc'
      include 'commands.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'potential_par.inc'
      include 'lattice.inc'
      character*256 fname1,fname2
      character*20 oldmolnum,newmolnum
      character*2 version
      character*3 oldmolname,newmolname,search1,search2
      external cosd,sind

c Initialise random number generator

      call initran1(10)

      write(6,600)
600      format(/'pdbtoato> Give name of .pdb file to convert: ',$)
      read(5,598) fname1

c Set up the distance to be used to characterise bonded atoms

      write(6,597)
597      format(/'pdbtoato> Give maximum distance to classify bonds')
      read(5,*) bmax

c Set up a search sequence for molecules which are distinct in the pdb file
c but are paired in practice

      write(6,596)
596      format(/'pdbtoato> Give two molecule names to join (a3, a3)')
      read(5,*) search1,search2

      open(12,file=fname1,status='old',iostat=ierr)
c
c if there is a problem with this file, print out an error message and exit
c
      if(ierr.ne.0) then
         ilenfname=index(fname1,'.')
         ilenfname=ilenfname+index(fname1(ilenfname:256),' ')-1
         write(6,599) fname1(1:ilenfname)
599      format(/'pdbtoato> Problems with file ',a/)
         stop
      endif

c Read the pdb file line by line and search for particular sequences

c Set some values not defined by .pdb file

      simtk=300
      tol=0.0
      stepmi=0.1
      stepri=1.0
      stepr=1.0
      stepm=0.3
      tempb=65.0
      phib=90.0
      thetac=0.0
      phic=0.0
      tolout=tol
      icount=0

c Set up look to generate multiple ato files from a single pdb file

      do while(ierr.eq.0)

      ia=0
      ig=0
      im=0
      nmol=0
      natom=0
      natomr=0
      ngrrot=0
      ntype=0
      nmoltype=0
      oldmolname=' '
      oldmolnum=' '
      iter=0

c Read lines from the .pdb file and look for particular sequences

      do while (ierr.eq.0.and.iter.eq.0.and.im.lt.mmol.and.ia.lt.matom)

         read(12,598,iostat=ierr) sometext
598      format(a)

c Split this text into its component words

         call parse(sometext,nwords,ncf,ncl,10,nstsize,nworddim)

c      write(6,*) nwords, sometext(ncf(1):ncl(1))



c Look for the word "CRYST1".

         if(sometext(ncf(1):ncl(1)).eq.'CRYST1'.and.nwords.ge.4.and.ierr.eq.0) then

c Next word is the lattice constant

            read(sometext(ncf(2):ncl(2)),*) alat
            read(sometext(ncf(3):ncl(3)),*) blat
            read(sometext(ncf(4):ncl(4)),*) clat
            dmax=alat

         else if(sometext(ncf(1):ncl(1)).eq.'ATOM'.and.nwords.ge.9.and.ierr.eq.0) then

c An atom is found. There must be at least 9 words on this line for us to be able
c to use it.

c Find the molecule label and increment the molecule counter if this is different from the 
c previous molecule.

c Set the type of the new molecule

            newmolname=sometext(ncf(4):ncl(4))

c Likewise new molecule number

            newmolnum=sometext(ncf(nwords-5):ncl(nwords-5))

            if(newmolnum.ne.oldmolnum
     *.and.oldmolname.ne.search1.and.newmolname.ne.search2) then

c If the molecule number has changed, close the current molecule and calculate its centre of mass
c and adjusted internal coordinates.

               if(im.gt.0) then

                nato(im)=ia-naf(im)+1
                natom=natom+nato(im)
                nal(im)=ia
                xsum=0.0
                ysum=0.0
                zsum=0.0
                do ib=naf(im),nal(im)
                   xsum=xsum+xa(ib)
                   ysum=ysum+ya(ib)
                   zsum=zsum+za(ib)
                end do
                xadd(im)=xsum/nato(im)
                yadd(im)=ysum/nato(im)
                zadd(im)=zsum/nato(im)
                do ib=naf(im),nal(im)
                   xa(ib)=xa(ib)-xadd(im)
                   ya(ib)=ya(ib)-yadd(im)
                   za(ib)=za(ib)-zadd(im)
                end do
                phix(im)=0.0
                phiy(im)=0.0
                phiz(im)=0.0

c Check to see whether this molecule type already exists

                  imoltype=0
                  ifound=0
                  do while(imoltype.lt.nmoltype.and.ifound.eq.0)
                   imoltype=imoltype+1
                   if(oldmolname.eq.moltype(imoltype)) ifound=imoltype
                   mtype(im)=ifound
                  end do
                  if(ifound.eq.0) then
                   nmoltype=nmoltype+1
                   moltype(nmoltype)=oldmolname
c
c initialise the template filename
c
                   moltemplate(nmoltype)=oldmolname
                   mtype(im)=nmoltype

c Initialise the rotational and translational step sizes

                   wholerot(nmoltype)=stepr
                   wholetrans(nmoltype)=stepm
                  endif
c Set the number of internal rotation groups
                ngrot(im)=0

             endif

c Start a new molecule

               im=im+1
             nmol=im
             nato(im)=0
             naf(im)=ia+1

            endif

            oldmolname=newmolname
            oldmolnum=newmolnum

            ia=ia+1
            nato(im)=nato(im)+1
            mref(ia)=im
            read(sometext(ncf(3):ncl(3)),*) alabel(ia)

c set up 'type' of atom - this is used in assigning neutron weights

            call ftype(ia)
c
c do not count atoms which begin with 'Q' or 'q' - these are assumed 
c invisible to the scattering radiation and experiment (although they
c can be used to calculate other structural quantities if necessary
c
            if(qtype(ia).ne.1) then
             natomr=natomr+1
            else
             iflag=1
            endif
c get its starting coordinates and sigmas
            read(sometext(ncf(nwords-4):ncl(nwords-2)),*) xa(ia),ya(ia),za(ia)
c read the number of atoms to which it is bound, the atom numbers of
c those neighbours and the anticipated bond length 
            nbond(ia)=0

         else if(sometext(ncf(1):ncl(1)).eq.'ENDMDL'.and.ierr.eq.0) then

c Process the last molecule and force an exit

            if(im.gt.0) then

             natom=natom+nato(im)
             nal(im)=ia
             xsum=0.0
             ysum=0.0
             zsum=0.0
             do ib=naf(im),nal(im)
                xsum=xsum+xa(ib)
                ysum=ysum+ya(ib)
                zsum=zsum+za(ib)
             end do
             xadd(im)=xsum/nato(im)
             yadd(im)=ysum/nato(im)
             zadd(im)=zsum/nato(im)
             do ib=naf(im),nal(im)
                xa(ib)=xa(ib)-xadd(im)
                ya(ib)=ya(ib)-yadd(im)
                za(ib)=za(ib)-zadd(im)
             end do
             phix(im)=0.0
             phiy(im)=0.0
             phiz(im)=0.0

c Check to see whether this molecule type already exists

               imoltype=0
               ifound=0
               do while(imoltype.lt.nmoltype.and.ifound.eq.0)
                imoltype=imoltype+1
                if(oldmolname.eq.moltype(imoltype)) ifound=imoltype
                mtype(im)=ifound
               end do
               if(ifound.eq.0) then
                nmoltype=nmoltype+1
                moltype(nmoltype)=oldmolname
c
c initialise the template filename
c
                moltemplate(nmoltype)=oldmolname
                mtype(im)=nmoltype

c Initialise the rotational and translational step sizes

                wholerot(nmoltype)=stepr
                wholetrans(nmoltype)=stepm
               endif
c Set the number of internal rotation groups
             ngrot(im)=0

            endif

c Force the exit

            iter=1
            icount=icount+1

         endif
      end do

c For each atom in each molecule, make a list of the bonded atoms

      do im=1,nmol
         if(nato(im).gt.1) then
            do ia1=naf(im),nal(im)-1
               do ia2=ia1+1,nal(im)
                xdiff=(xa(ia1)-xa(ia2))/dmax
                ydiff=(ya(ia1)-ya(ia2))/dmax
                zdiff=(za(ia1)-za(ia2))/dmax

c Minimum image convention

                xdiff=xdiff-nint(xdiff)
                ydiff=ydiff-nint(ydiff)
                zdiff=zdiff-nint(zdiff)

                dtest=dmax*sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff)
                if(dtest.le.bmax) then
                  nbond(ia1)=nbond(ia1)+1
                  nbond(ia2)=nbond(ia2)+1
                  iabond(ia1,nbond(ia1))=ia2
                  iabond(ia2,nbond(ia2))=ia1
                  tbond(ia1,nbond(ia1))=dtest
                  tbond(ia2,nbond(ia2))=dtest
                endif
             end do
            end do
         endif
      end do
c make a list of the atom types read in from the current file

      if(ierr.eq.0) then
      ntypefile=0
      ia=0
      do while(ia.lt.natom)
         ia=ia+1
         itypefile=0
         ifound=0
         do while(itypefile.lt.ntypefile.and.ifound.eq.0)
            itypefile=itypefile+1
            if(alabel(ia).eq.typefile(itypefile)) ifound=itypefile
         end do
         if(ifound.eq.0) then
            ntypefile=ntypefile+1
            typefile(ntypefile)=alabel(ia)
         endif
      end do

c Write the atom types out as a check

      write(6,204)
204      format(/'pdbtoato> Following atom types found:-'/)

      do it=1,ntype

         write(6,202) it,type(it)

      end do

c Write the molecule types out as a check

      write(6,201)
201      format(/'pdbtoato> Following molecule types found:-'/)

      do imoltype=1,nmoltype

         write(6,202) imoltype,moltype(imoltype),wholerot(imoltype)
     *,wholetrans(imoltype)
202      format(1x,i3,1x,a3,2(1x,e12.5))

      end do
c
c read core energy and diameter - used in fmole for non-bonded atoms
c
      ecore=1.0
      dcore=3.0

c set up lattice parameters for applying the minimum image convention

      alatx=alat
      blatx=blat*cosd(phib)
      blaty=blat*sind(phib)
      clatx=clat*sind(thetac)*cosd(phic)
      claty=clat*sind(thetac)*sind(phic)
      clatz=clat*cosd(thetac)
c
c define the Cartesian coordinates of the lattice vectors
c
      alatc(1)=alatx
      alatc(2)=0.0
      alatc(3)=0.0
      blatc(1)=blatx
      blatc(2)=blaty
      blatc(3)=0.0
      clatc(1)=clatx
      clatc(2)=claty
      clatc(3)=clatz
      bxax=blatx/alatx
      cxax=clatx/alatx
      cyby=claty/blaty
c
c calculate the perpendicular distance between the faces of the lattice
c
      perpab=perpd(alatc,blatc,clatc)
      perpbc=perpd(blatc,clatc,alatc)
      perpca=perpd(clatc,alatc,blatc)
c      write(6,*) perpab,perpbc,perpca
c
c set up maximum radius to be used for Fourier transforms and to truncate
c the potential. This is defined as the minimum perpendicular distance between
c planes. 
c
c G(r) and potentials will be calculated to rmax, or maximum allowed distance
c based on number of points and step size, whichever is smaller.
c
      dmax=amin1(perpab,perpbc,perpca)
      cell=dmax
      write(6,*) cell
c
c convert temperature in K to kJ/mole
c
      simt=0.008314*simtk
c
c calculate volume of box and atomic number density
c
      volume=alatx*blaty*clatz
      rho=natomr/volume
      write(6,*) volume,rho
c
c calculate atomic fractions
c
      do it=1,ntype
            cfrac(it)=0
      end do
      do ia=1,natom
            iref=itype(ia)
            cfrac(iref)=cfrac(iref)+1
      end do
      do it=1,ntype
            cfrac(it)=cfrac(it)/natomr
      end do
      do it=1,ntype
            write(6,916) it,type(it),cfrac(it)
      end do
916      format(1x,i3,1x,'Atomic fraction of ',a3,' = ',e12.5)
c
c setup molecular size factor
c
      dmaxmol=dmax*sizefactor
c
c calculate mass of each molecule
c
      do im=1,nmol
            amasst(im)=0.0
            do ia=naf(im),nal(im)
                  if(qtype(ia).eq.0) then
                     amass=amassa(itype(ia))
                     amasst(im)=amasst(im)+amass
                  endif
            end do
            amasstsqrt(im)=sqrt(amasst(im))
c
c set maximum distance of each atom from centre of mass
c
            moldiam(im)=0.0
            do ia=naf(im),nal(im)
                  sum=xa(ia)*xa(ia)+ya(ia)*ya(ia)+za(ia)*za(ia)
                  if(sum.gt.moldiam(im)) moldiam(im)=sum
            end do
            moldiam(im)=sqrt(moldiam(im))
      end do
      ip=0
      do it=1,ntype
            do jt=it,ntype
                  ip=ip+1
                  n1extr(ip)=it
                  n2extr(ip)=jt
                  npairref(it,jt)=ip
                  npairref(jt,it)=ip
                  ect=sqrt(epscor(it)*epscor(jt))
                  ecor(ip)=ect
                  qcor(ip)=charge(it)*charge(jt)*1390
                  dct=0.5*(sigcor(it)+sigcor(jt))
                  dcor(ip)=dct
                  dscor(ip)=dct*dct
                  am1=amassa(it)
                  am2=amassa(jt)
                  amassr(ip)=tempb*simt*sqrt(am1*am2/(am1+am2))
            end do
      end do
      npairs=ip
      ndist=npairs

      write(6,917) nmol,natomr,ntype,npairs
917      format(4(1x,i5))

c Generate the name of the output file

      write(version,918) icount
918      format(i2.2)
      fname2=fname1(1:index(fname1,'.')-1)//version//'.ato'
      ilenfname2=index(fname2,'.')+3

      write(6,915) fname2(1:ilenfname2)
915      format(/'pdbtoato> Output filename is ',a)

c Write the new .ato file

      call moloutc(fname2)

      endif

c Read the next file if it exists

      end do

      close (12)

      stop
      end


      subroutine ftype(ia)
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      character*1 qtypet
      character*3 ats,typetest
      integer*4 iexch
c
c is this atom of "Q"-type?
c
      qtypet=alabel(ia)(1:1)
      if(qtypet.eq.'Q'.or.qtypet.eq.'q') then
            qtype(ia)=1
      else
            qtype(ia)=0
      endif
      ifound=0
      it=0
      do while(it.lt.ntype.and.ifound.eq.0)
         it=it+1
         if(alabel(ia).eq.type(it)) ifound=it
      end do
      if(ifound.gt.0) then
         itype(ia)=ifound
      else if(ntype.lt.mcomp) then
         ntype=ntype+1
         type(ntype)=alabel(ia)
         qtypea(ntype)=qtype(ia)
         itype(ia)=ntype

c Open labels.txt and see if this atom type has already been defined

         open(20,file='labels.txt',status='old',iostat=ierr)
         ifound=0
         do while(ierr.eq.0.and.ifound.eq.0)
            read(20,*,iostat=ierr) typetest,epsit,sigt,amass,qt,ats
            if(typetest.eq.type(ntype)) ifound=1
         end do
         if(ifound.eq.0) then
            write(6,108) type(ntype)
108      format(/'pdbtoato> Give epsilon, sigma, mass, charge and atomic symbol for atom type: ',a3)
            read(5,*) epsit,sigt,amass,qt,ats
         endif

         rqt=0.0
         iexch=-1
         atsymbol(ntype)=ats
         iexchange(ntype)=iexch
         epscor(ntype)=epsit
         sigcor(ntype)=sigt
         amassa(ntype)=amass
         charge(ntype)=qt
         qradius(ntype)=rqt

c Write the new values at the end if no values for this atom exist

         if(ifound.eq.0) then
            write(20,101) type(ntype),epscor(ntype),sigcor(ntype),amassa(ntype)
     *,charge(ntype),atsymbol(ntype)
101      format(1x,a3,4(1x,e13.6),1x,a3)
         endif
         close(20)
      else
         write(6,100) mcomp
100      format(1x,'Tried to generate more than ',i3,' types of atom')
         stop
      endif
      return
      end

      function perpd(a,b,c)
      real*4 a(3),b(3),c(3),t(3)
c
c calculates the perpendicular projection of vector c onto the direction normal
c to the plane defined by vectors a and b
c
c form vector product a X b
c
      call vecpro(a,b,t)

c
c normalise it to a unit vector
c
      amp=sqrt(scapro(t,t))
      do i=1,3
            t(i)=t(i)/amp
      end do
c
c project c onto this normal vector
c
      perpd=scapro(c,t)
      return
      end      

      subroutine vecpro(a,b,c)
      real*4 a(3),b(3),c(3)
c performs the vector product of c1 x c2 and puts result in c3
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end

      function scapro(a,b)
      real*4 a(3),b(3)
      scapro=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      return
      end

      subroutine parse(sometext,nwords,ncf,ncl,nspcword
     *,nchartext,nworddim)
c
c parses the text in sometext to separate words with the assumed delimiter a space
c ncf and ncl contain the first and last character numbers of each word
c ncharword is the maximum number of characters in a word, nworddim is dimension of the
c word array, nspcword is the number of spaces between words needed to decide that the 
c end of text has been reached.
c
      character sometext
      integer*4 nwords,ncf(*),ncl(*)
     *,nspcword,nchartext,nworddim
c
c step through the characters and identify words (separated by spaces). Extra spaces
c are ignored unless the number of spaces exceeds nspcword, in which case the search
c for words terminates
c
      nwords=0
      nfirst=1
      nlast=0
      nspace=0
      i=0
      do while(i.lt.nchartext.and.nspace.le.nspcword)

         i=i+1

c       
c       check the number of characters in the current word - if equal to ncharword we save 
c       the current word
c       
         if(sometext(i:i).eq.' ') then
c       
c       A space is found
c       
c       check that this is not the end of a word
c       
            if(nlast.ge.nfirst) then
c       
c       new word found
c       
             if(nwords.lt.nworddim) then
                nwords=nwords+1
                ncf(nwords)=nfirst
                ncl(nwords)=nlast
                nspace=0
             else
                return
             end if
            endif
            nspace=nspace+1
            nfirst=i+1
         else
            nlast=i
         endif
      end do
      return
      end

      subroutine moloutc(fname1)
      include 'epsr_dimension.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      include 'lattice.inc'
      character*256 fname1
      open(10,file=fname1,status='unknown',iostat=ierr)
c
c if there is a problem with this file, print out an error message and exit
c
      if(ierr.ne.0) then
         ilenfname=index(fname1,'.')
         ilenfname=ilenfname+index(fname1(ilenfname:256),' ')-1
         write(6,599) fname1(1:ilenfname)
599      format(/'moloutc> Problems with file ',a/)
         nmol=0
         return
      endif

c      Write a warning message if ROT groups are defined, but internal rotations are
c       switched off.

      if(ngrrot.gt.0.and.stepri.le.0.0) then
         write(6,598)
598      format(/'moloutc> WARNING! ROT groups defined BUT step size is ZERO!')
      endif

      ig=0
c
c write no. of molecules and tolerance on bond distances
c
c      write(10,207) nmol,simtk
c      write(10,208) alat,blat,clat
c      write(10,208) phib,thetac,phic
      write(10,207) nmol,dmax,simtk
      write(10,206) tol,stepmi,stepri,stepr
     1,stepm,tempb
      do im=1,nmol
      write(10,200) nato(im),xadd(im),yadd(im)
     1,zadd(im),phix(im),phiy(im),phiz(im),im
c
c integer to convert from absolute atom numbers to relative atom numbers
c
      nadd=1-naf(im)
      do ia=naf(im),nal(im)
c for each atom write its label
            write(10,201) alabel(ia), (ia-naf(im)+1)
c write its starting coordinates and sigmas
            write(10,202) xa(ia),ya(ia),za(ia)
c
c write the number of atoms to which it must be bound, the atom numbers of
c those neighbours and the anticipated bond length 
c
            write(10,203) nbond(ia),((iabond(ia,j)+nadd),tbond(ia,j)
     1,j=1,nbond(ia))
      end do
c write number of internal rotation groups
      write(10,204) ngrot(im)
      if(ngrot(im).gt.0) then
            do ih=1,ngrot(im)
                  ig=ig+1
                  write(10,201) rotgrp(ig)
c for rotational groups read atom numbers of two atoms to form axis of rotation
                  write(10,204) (labaxi(ig,1)+nadd),(labaxi(ig,2)+nadd)
c write number of atoms in group and their numbers
                  write(10,204) nrot(ig),((labrot(ig,j)+nadd),j=1,nrot(ig))
            end do
      endif
      end do
      do it=1,ntype
            write(10,211) type(it),atsymbol(it),iexchange(it)
            write(10,202) epscor(it),sigcor(it),amassa(it),charge(it)
     1,qradius(it)
      end do
c write core energy and diameter which represents the closest overlap of ANY 
c non-bonded atom
      write(10,202) ecore,dcore
c
c save random number generator parameters
c
      call saveran1(10)
c
c write out the molecular template files if defined
c
      do imoltype=1,nmoltype

         ilen=index(moltemplate(imoltype),' ')-1
         write(10,212) imoltype,moltemplate(imoltype)(1:ilen)
     *,wholerot(imoltype),wholetrans(imoltype)

      end do
      close(10)
      return
200      format(1x,i4,6(1x,e13.6),1x,i5)
207      format(1x,i5,7(1x,e13.6))
208      format(1x,6(1x,e13.6))
206      format(1x,6(1x,e10.3))
201      format(1x,a3,5x,i5)
202      format(1x,6(e12.5,1x))
203      format(1x,i2,(5(1x,i4,1x,e9.3)))
204      format(1x,14(i4,1x))
205      format(1x,i32)
211      format(1x,a3,1x,a3,1x,i2)
212      format(1x,i3,1x,a,2(1x,e13.6))
      end




            



