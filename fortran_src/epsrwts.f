      program standalone

c Program to perform the EPSR setup command

c 2009/3/30 Attempt to make this a standalone program

c Get the command from the command line arguments

      call get_command_from_args()

c Run the command

      call epsrwts()

      stop
      end      

      subroutine epsrwts()
c
c Program EPSRWTS
c
c Calculates weights for a given input data (S(Q)) file or set of files, based on a specified
c .ATO file
c
c Currently set up for 70 components (number of atom types)
c
      include 'epsr_dimension.inc'
      parameter(misotope=5)
      include 'commands.inc'
      include 'symbols.inc'
      include 'directory_name.inc'
      include 'molecom.inc'
      include 'clustercom.inc'
      character*256 readfile
      character*4 ext
      character*1 iaccept
      real*4 batom(mcomp,misotope,2),abund(mcomp,misotope,2)
     1,cprod(mpairs)
      real*4 binter(mcomp,3)
      real*4 wintra(mpairs,3),winter(mpairs,3),wtsum(3)
      integer*4 iflagsubs(mcomp),nisotope(mcomp,2)
      integer*4 massno
      character*8 partlabel(mpairs),fnameb
      character*80 fname
      character*3 chsymtest
c
c generate the name of the .ato file
c
      ext='.ato'
      ilenext=4
      if(nwords.lt.2) then
         nwords=0
         do while (nwords.le.0)
            call search_for_files(ext,ilenext,iaccept)
            call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
            if(nwords.le.0) then
               write(6,499)
499   format(/'epsrwts> No .ato file specified! - routine will exit'/)
               return
            endif
         end do
         readfile=sometext(ncf(1):ncl(1))
         ilenreadfile=ncl(1)-ncf(1)+1
      else
         readfile=sometext(ncf(2):ncl(2))
         ilenreadfile=ncl(2)-ncf(2)+1
      endif
c
c strip out any extension that was supplied - extension must be .ato, otherwise
c nothing is done.
c
      ilendot=index(readfile,'.')-1
      if(ilendot.gt.0) then
         readfile=readfile(1:ilendot)
         ilenreadfile=ilendot
      endif
c
c create the full name of the input file by adding on the directory name if not present
c
      call get_full_filename(readfile,ilenreadfile,fnameato,ilenfnameato)
      fnameato=fnameato(1:ilenfnameato)//'.ato'
      ilenfnameato=ilenfnameato+4
c
c initialise number of atom and molecule types to zero
c
      ntype=0
      nmoltype=0
c
c read .ato file
c
      sizefactor=1.0
      call minitc(fnameato,rmaxcell,cellmax,upol,rho,amubar,amubar2
     *,ndist,ierr)
      if(rho.eq.0.0.or.ierr.ne.0) then
         write(6,599) fnameato
599   format(/'epsrwts> Problems with .ato file ',a/'Routine will exit now'/)
         return
      endif
      write(6,*)
      write(6,100)
100   format(/'epsrwts> Program to calculate inter- and intra-molecular'
     1,' weightings '
     1,/1x,'for DCS, 1st- or 2nd-order difference data.'/)
99    format(a80)
      write(6,200)
200   format(/'epsrwts> Is the output to be per atom (1), or per molecule (2)? ',$)
      read(5,*) nrho
c
c get the no. of atoms per molecule ratio in case it is needed
c
      atmolrat=float(natomr)/float(nmol)
      write(6,202) (it,type(it),cfrac(it),atsymbol(it),it=1,ntype)
202   format(/'epsrwts> The following components were found in this file'
     1/1x,
     1/1x,'Component no., label,  atomic fraction, chemical symbol'
     1/1x,
     2/((5x,i5,8x,a3,5x,e12.5,5x,a3)))
203   continue
      write(6,219)
219   format(/'epsrwts> How many samples (1,2, or 3)? (0 to quit) ',$)
      read(5,*)  norder
c
c the program assumes
c norder = 1 means this is single dataset
c norder = 2 means this is a first-order difference experiment
c norder = 3 means this is a second-order difference experiment
c
      if(norder.lt.1) return
      if(norder.gt.3) go to 203
c
c set flags to say which (if any) atoms have been substituted
c
      do it=1,ntype
            iflagsubs(it)=0
      end do
      nsubs=0
c
c Get the scattering lengths for the first sample
c
      if(norder.gt.1) then
            write(6,207)
207   format(/'epsrwts> For sample 1, get the scattering lengths for'
     1,' all components in the sample'/)
      else
            write(6,223)
223   format(/'epsrwts> Get the scattering lengths for all'
     1,' components in the sample'/)
      end if
      do it=1,ntype
c
c cannot see "Q" type atoms
c
         if(qtypea(it).eq.0) then
206         continue
            chsymtest=atsymbol(it)
            massno=0
c
c test to see if this chemical symbol exists in the table
c
            btest=scatt_len(chsymtest,massno,ifind)
            if(ifind.eq.0) then
               write(6,209) chsymtest,type(it)
209   format(/'epsrwts> Did not recognise chemical symbol ',a3
     *,' for atom type ',a3//'epsrwts> Try typing it again: ',$)
               read(5,*,iostat=ierr) atsymbol(it)
               go to 206
            endif
c
c for hydrogen atoms we need to check whether they exchange with other hydrogens
c
            if(atsymbol(it).eq.'H') then
               iexchange(it)=-1
               do while(iexchange(it).ne.0.and.iexchange(it).ne.1)
                  write(6,103) type(it)
103   format(/'epsrwts> For hydrogen component ',a3,':-'
     *//' If it exchanges with atoms on other molecules type 1,'
     *,' if not type 0: ',$)
                  read(5,*,iostat=ierr) iexchange(it)
               end do
               if(iexchange(it).ne.1) iexchange(it)=0
            else
               iexchange(it)=0
            endif
            asum=0
            nsum=1
205         continue
            write(6,211) type(it)
211   format(/'epsrwts> For component ',a3
     *//' Type 0 for a natural isotope or mass number'
     *,' for a specific isotope: ',$)
            read(5,*,err=205) massno
            write(6,224)
224   format(/' and its abundance (0.0-1.0): ',$)
            read(5,*,err=205)  abund(it,nsum,1)
            if(massno.eq.0) then
               batom(it,nsum,1)=btest
            else
               batom(it,nsum,1)=scatt_len(chsymtest,massno,ifind)
               if(ifind.eq.0) then
                  write(6,212) atsymbol(it)
212   format(/'epsrwts> Did not find this massno for chemical'
     1,' symbol ',a3,' - try typing again'/)
                  go to 205
               endif
            endif
c     write(6,*) it,nsum,batom(it,nsum,1)
            asum=asum+abund(it,nsum,1)
            if(asum.lt.1.0) then
               nsum=nsum+1
               go to 205
            else
               do isum=1,nsum
                  abund(it,nsum,1)=abund(it,nsum,1)/asum
               end do
            endif 
            nisotope(it,1)=nsum
         endif
      end do
c
c we now need to check that components which
c are stated to be exchangeable through the sample do indeed all have the 
c the same scattering length.
c
c Step through the components and if iexchange is 1 form the average scattering 
c for all the atoms with the same chemical symbol.
c
      bscatsum=0.0
      cfracsum=0.0
      do it=1,ntype
            if(iexchange(it).eq.1) then
                  cfracsum=cfracsum+cfrac(it)
                  do ii=1,nisotope(it,1)
                        term=cfrac(it)*abund(it,ii,1)*batom(it,ii,1)
                        bscatsum=bscatsum+term
                  end do
            endif
      end do
      if(cfracsum.gt.0.0) then
            bscatsum=bscatsum/cfracsum
      else
            bscatsum=0.0
      endif
c
c go through each atomic type and if iexchange is 1, then replace the individual
c scattering lengths with a single average scattering length
c
      do it=1,ntype
            if(iexchange(it).eq.1) then
                  nisotope(it,1)=1
                  abund(it,1,1)=1.0
                  batom(it,1,1)=bscatsum
            endif
      end do
c
c now for the second sample, if appropriate identify those atoms
c which are substituted and given the relevant scattering lengths
c for the substituted atoms.
c
      if(norder.gt.1) then
204         continue
            write(6,213)
213   format(/'epsrwts> For sample 2:-'
     *//' Give the component number of '
     *,'a substituted atom, (0 to exit) ',$)
c
c obviously it is important to ensure you do not specify an atom
c which previously was treated as if only 1 isotope exists. The
c programme does not do any checking on this.
c
            read(5,*,err=204) isubs
            if(isubs.ge.1.and.isubs.le.ntype) then
                  if(qtypea(isubs).ne.0) then
                     write(6,155)
155   format(/'epsrwts> Cannot substitute on a Q atom!')
                     goto 204
                  endif
                  nsubs=nsubs+1
                  iflagsubs(isubs)=1
c
c read in mass number and abundance
c
                  asum=0
                  nsum=1
208         continue
                  write(6,211) type(isubs)
                  read(5,*,err=208)  massno
                  write(6,224)
                  read(5,*,err=208) abund(isubs,nsum,2)
                  chsymtest=atsymbol(isubs)
                  batom(isubs,nsum,2)=scatt_len(chsymtest,massno,ifind)
                  if(ifind.eq.0) then
                        write(6,212) atsymbol(isubs)
                        go to 208
                  endif
                  asum=asum+abund(isubs,nsum,2)
                  if(asum.lt.1.0) then
                     nsum=nsum+1
                     go to 208
                  else
                     do isum=1,nsum
                        abund(isubs,nsum,2)=abund(isubs,nsum,2)/asum
                     end do
                  endif 
                  nisotope(isubs,2)=nsum
                  go to 204
            endif
c
c for all the NON-substituted atoms, set the relevant scattering
c lengths and parameters to those of 1st sample
c
            do it=1,ntype
                  if(iflagsubs(it).eq.0) then
                        nisotope(it,2)=nisotope(it,1)
                        do ib=1,nisotope(it,2)
                              batom(it,ib,2)=batom(it,ib,1)
                              abund(it,ib,2)=abund(it,ib,1)
                        end do
                  endif
            end do
c
c we now need to check that components which have the same chemical symbol and which
c are stated to be exchangeable through the sample do indeed all have the 
c the same scattering length.
c
c Step through the components and if iexchange is 1 form the average scattering 
c for all the atoms with the same chemical symbol.
c
            bscatsum=0.0
            cfracsum=0.0
            do it=1,ntype
                  if(iexchange(it).eq.1) then
                        cfracsum=cfracsum+cfrac(it)
                        do ii=1,nisotope(it,2)
                              term=cfrac(it)*abund(it,ii,2)*batom(it,ii,2)
                              bscatsum=bscatsum+term
                        end do
                  endif
            end do
            if(cfracsum.gt.0.0) then
                  bscatsum=bscatsum/cfracsum
            else
                  bscatsum=0.0
            endif
c
c go through each atomic type and if iexchange is 1, then replace the individual
c scattering lengths with a single average scattering length
c
            do it=1,ntype
                  if(iexchange(it).eq.1) then
                        nisotope(it,2)=1
                        abund(it,1,2)=1.0
                        batom(it,1,2)=bscatsum
                  endif
            end do
      endif
c
c if this is a 2nd-order difference experiment, then third sample is
c treated as a mixture of the first two. If this is not the case then
c you can either use this programme to set up the weights for each DCS file
c separately, and so fit each total DCS, or else write another version!
c
      if(norder.eq.3) then
            write(6,214) 
214   format(/'epsrwts> For sample 3, specify the mixing fraction'
     *//' (This should be the same as that used'
     *,' in the .hdb_dat file): ',$)
            read(5,*) fracmix
            frac1=fracmix
            frac2=(1.0-fracmix)
      else
            frac1=1.0
            frac2=0.0
      endif
c           
c for each sample calculate both intra-molecular and inter-molecular weights
c for the total differential cross-sections
c
c total number of unique pairs expected
c
      npart=ntype*(ntype+1)/2
c
c calculate the weights for the first one or two samples
c
c weights for the third sample if present will be calculated
c as an average of the first two samples
c
      if(norder.gt.2) then
            ndata=2
      else
            ndata=norder
      endif
c
c go through all the components and calculate the isotope averaged scattering
c length for datasets 1 and 2
c
      do idata=1,ndata
            do it1=1,ntype
               if(qtypea(it1).eq.0) then
c
c average scattering length of this component
c
                  binter(it1,idata)=0.0
                  do ib=1,nisotope(it1,idata)
                        term1=abund(it1,ib,idata)*batom(it1,ib,idata)
                        binter(it1,idata)=binter(it1,idata)+term1
                  end do
               endif
            end do
c
c step through all component pairs and write down the intra- and inter-
c molecular weights
c
c check the number of partials generated
c
            ipart=0
c
c sum of weights
c
            wtsum(idata)=0.0
c
c first components
c
            do it1=1,ntype
                  binter1=binter(it1,idata)
c
c step through second components
c
                  do it2=it1,ntype
                     ipart=ipart+1
                     if(qtypea(it1).eq.0.and.qtypea(it2).eq.0) then
c
c if idata=1 write a label for this pair of components
c
                        if(idata.eq.1) then
                              write(partlabel(ipart),108) 
     1 type(it1),type(it2)
108   format(1x,a3,1x,a3)
c
c atomic fraction product for this pair
c
                              if(it1.eq.it2) then
                                    cprod(ipart)=cfrac(it1)**2
                              else
                                    cprod(ipart)=2.0*cfrac(it1)*cfrac(it2)
                              endif
                        endif
c
c average scattering length of second component
c
                        binter2=binter(it2,idata)
c
c scattering length product for intermolecular terms
c
                        winter(ipart,idata)=binter1*binter2
c
c for intra-molecular terms involving H if both components are non-exchanging
c and no. of isotopes is the same for both components, need to sum
c the scattering length product for each isotope.
c
c Note that no checking has been done to see if this component pair
c occurs on the same molecule. Strictly if it does not the intramolecular
c weight should be set to zero. However EPSR will never find
c such pairs.
c
                        if(iexchange(it1).eq.0.and.iexchange(it2).eq.0
     1.and.nisotope(it1,idata).eq.nisotope(it2,idata)
     1.and.atsymbol(it1).eq.'H '.and.atsymbol(it2).eq.'H ') then
                              bintra=0.0
                              do ib=1,nisotope(it1,idata)
c also the abundances have to be identical, but we won't check
c since if they are not, then this is probably not a valid
c correlation in any case
                                    term=abund(it1,ib,idata)
     1*batom(it1,ib,idata)*batom(it2,ib,idata)
                                    bintra=bintra+term
                              end do
                              wintra(ipart,idata)=bintra
                        else
                              wintra(ipart,idata)=winter(ipart,idata)
                        endif
                        term=cprod(ipart)*winter(ipart,idata)
                        wtsum(idata)=wtsum(idata)+term
                     endif
                  end do
            end do
      end do
c
c now do third sample if present
c
      if(norder.eq.3) then
            ipart=0
            wtsum(3)=0
            do it1=1,ntype
               if(qtypea(it1).eq.0) then
c
c average scattering length for this component
c
                  binter(it1,3)=frac1*binter(it1,1)+frac2*binter(it1,2)
               endif
            end do
            do it1=1,ntype
               if(qtypea(it1).eq.0) then
                  do it2=it1,ntype
                     if(qtypea(it2).eq.0) then
                        ipart=ipart+1
c
c for the intermolecular terms, we just have the random mixing result
c
                        winter(ipart,3)=binter(it1,3)*binter(it2,3)
c
c for intramolecular terms, we have to mix the intra-molecular terms in proportion
c unless the atoms are exchangeable, in which you get just the random mixing result
c
                        if(iexchange(it1).eq.0.and.iexchange(it2).eq.0) then
                              wintra(ipart,3)=frac1*wintra(ipart,1)
     1+frac2*wintra(ipart,2)
                        else
c
c otherwise we just have the random mixing result
c
                              wintra(ipart,3)=winter(ipart,3)
                        endif
                        term=cprod(ipart)*winter(ipart,3)
                        wtsum(3)=wtsum(3)+term
                     endif
                  end do
               endif
            end do
      endif
c
c write out the total weights for the dataset(s)
c
      nwords=0
      do while (nwords.eq.0)
         write(6,220)
220   format(/'epsrwts> Type basename of file to output weights to: ',$)
         read(5,111) sometext
111   format(a)
         call parse(sometext,nwords,ncf,ncl,5,nstsize,nworddim)
         fname=directory_name(1:nchardir)//sometext(ncf(1):ncl(1))
      end do
c
c find out how many characters there are in this name
c
      nchar=nchardir+ncl(1)-ncf(1)+1
      do idata=1,norder
            write(6,201) idata
201   format(/'epsrwts> For total data ',i2,' has data been normalised'
     1,' (1) or not (0)? ',$)
c
c ("Normalised" means divided by the sum of the weight factors)
c
            read(5,*) nstruct
c
c generate file name to write to
c
            if(norder.gt.1) then
                  write(fnameb,112) idata
112   format('tot',i1.1,'.wts')
                  nchart=nchar+8
            else
                  write(fnameb,109)
109   format('tot.wts')
                  nchart=nchar+7
            endif
c
c concatenate these two filenames
c
            fname=fname(1:nchar)//fnameb
            write(6,113) fname(1:nchart)
113   format(/'epsrwts> Writing TOTAL weights to file ',a)
            open(10,file=fname,status='unknown')
c
c signal whether these weights are per atom or per molecule
c
            if(nrho.eq.2) then
                  if(nstruct.eq.1) then
                        rhofac=1.0/atmolrat
                  else
                        rhofac=atmolrat
                  endif
            else
                  rhofac=1.0
            endif
            write(10,*) 1.0
            do ipart=1,npart
               it1=n1extr(ipart)
               it2=n2extr(ipart)
               if(qtypea(it1).eq.0.and.qtypea(it2).eq.0) then
                  winterwr=winter(ipart,idata)
                  wintrawr=wintra(ipart,idata)
c
c relative contribution to total
c
                  relativewt=winter(ipart,idata)*cprod(ipart)
c
c If required normalise the weights to the total sum of weights
c
                  if(nstruct.eq.1) then
                        winterwr=winterwr/wtsum(idata)
                        wintrawr=wintrawr/wtsum(idata)
                        relativewt=relativewt/wtsum(idata)
                  endif
                  write(10,114) partlabel(ipart)
114   format(a8)
                  write(10,*) wintrawr*rhofac,winterwr*rhofac,relativewt
               endif
            end do
            close(10)
      end do
c
c now if norder .eq. 2, take the first order difference between samples 1 and 2
c and write out the weights (but not if both intra and inter weights are zero
c
      if(norder.eq.2) then
            write(6,301) 
301   format(/'epsrwts> Has first order difference data data been normalised'
     1,' (1) or not (0)? ',$)
c
c ("Normalised" means divided by the sum of the weight factors)
c
            read(5,*) nstruct
c
c generate file name to write to
c
            write(fnameb,115) 
115   format('dif.wts')
c
c concatenate these two filenames
c
            fname=fname(1:nchar)//fnameb
            write(6,116) fname(1:nchar+7)
116   format(/'epsrwts> Writing 1st order difference weights to file ',a)
            open(10,file=fname,status='unknown')
c
c signal whether these weights are per atom or per molecule
c
            if(nrho.eq.2) then
                  if(nstruct.eq.1) then
                        rhofac=1.0/atmolrat
                  else
                        rhofac=atmolrat
                  endif
            else
                  rhofac=1.0
            endif
            write(10,*) 1.0
            wtsumdif=wtsum(1)-wtsum(2)
            do ipart=1,npart
               it1=n1extr(ipart)
               it2=n2extr(ipart)
               if(qtypea(it1).eq.0.and.qtypea(it2).eq.0) then
c
c only write out those partials which involve a substituted atom
c
               if(iflagsubs(it1).eq.1.or.iflagsubs(it2).eq.1) then
                  winterwr=winter(ipart,1)-winter(ipart,2)
                  wintrawr=wintra(ipart,1)-wintra(ipart,2)
c
c relative contribution to total
c
                  relativewt=cprod(ipart)*winterwr
c
c If required normalise the weights to the total sum of weights
c
                  if(nstruct.eq.1.and.wtsumdif.ne.0.0) then
                        winterwr=winterwr/wtsumdif
                        wintrawr=wintrawr/wtsumdif
                        relativewt=relativewt/wtsumdif
                  endif
                  write(10,114) partlabel(ipart)
                  write(10,*) wintrawr*rhofac,winterwr*rhofac,relativewt
               endif
               endif
            end do
            close(10)
      endif
c
c Do second order difference if appropriate
c
      if(norder.eq.3) then
c
c it will be assumed here that partial structure factors have
c been normalised
c
            nstruct=1
c
c generate 3 filenames to write to
c
            write(fnameb,117) 
117   format('xx.wts')
c
c concatenate these two filenames
c
            fname=fname(1:nchar)//fnameb
            write(6,*)
            write(6,118) fname(1:nchar+6)
118   format(1x,'Writing XX weights to file ',a)
            open(10,file=fname,status='unknown')
c
c signal whether these weights are per atom or per molecule
c
            if(nrho.eq.2) then
                  rhofac=1.0/atmolrat
            else
                  rhofac=1.0
            endif
            write(10,*) 1.0
            write(fnameb,119) 
119   format('xy.wts')
c
c concatenate these two filenames
c
            fname=fname(1:nchar)//fnameb
            write(6,*)
            write(6,120) fname(1:nchar+6)
120   format(1x,'Writing XY (=XH) weights to file ',a)
            open(11,file=fname,status='unknown')
c
c signal whether these weights are per atom or per molecule
c
            write(11,*) 1.0
            write(fnameb,121) 
121   format('yy.wts')
c
c concatenate these two filenames
c
            fname=fname(1:nchar)//fnameb
            write(6,*)
            write(6,122) fname(1:nchar+6)
122   format(1x,'Writing YY (=HH) weights to file ',a)
            open(12,file=fname,status='unknown')
c
c signal whether these weights are per atom or per molecule
c
            write(12,*) 1.0
c
c set up the total atomic fractions for X (non-substituted)
c and Y (substituted) atoms and their respective average
c scattering lengths
c
            cx=0
            cy=0
            cxbx=0
            cyby1=0
            cyby2=0
            do it=1,ntype
               if(qtypea(it).eq.0) then
                  if(iflagsubs(it).eq.1) then
                        cy=cy+cfrac(it)
                        cyby1=cyby1+cfrac(it)*binter(it,1)
                        cyby2=cyby2+cfrac(it)*binter(it,2)
                  else
                        cx=cx+cfrac(it)
                        cxbx=cxbx+cfrac(it)*binter(it,1)
                  endif
               endif
            end do
c
c write down the sum of weights for each of the YY, XY
c and XX terms
c
            wtsumyy=frac1*wtsum(1)+frac2*wtsum(2)-wtsum(3)
            wtsumxy=2.0*cxbx*(cyby1-cyby2)
            wtsumxx=cxbx*cxbx
c
c also write down the overall weighting of the YY term in the
c first order difference sample 1 - sample 2
c
            wtsumyydif=cyby1*cyby1-cyby2*cyby2
c 
c write down the weights of the YY and XY terms in the 1st sample
c
            wtsumyy1=cyby1*cyby1
            wtsumxy1=2.0*cxbx*cyby1
c
c now step through all the partials and write down their
c weights. It is assumed the data have been normalised.
c
            do ipart=1,npart
c
c get atom type numbers for this pair
c
               it1=n1extr(ipart)
               it2=n2extr(ipart)
               if(qtypea(it1).eq.0.and.qtypea(it2).eq.0) then
c
c YY terms
c
               if(iflagsubs(it1).eq.1.AND.iflagsubs(it2).eq.1) then
                  winteryy=(frac1*winter(ipart,1)+frac2*winter(ipart,2)
     1-winter(ipart,3))/wtsumyy
                  if(iexchange(it1).eq.1.or.iexchange(it2).eq.1) then
                     wintrayy=(frac1*wintra(ipart,1)+frac2*wintra(ipart,2)
     1-wintra(ipart,3))/wtsumyy
                  else
                     wintrayy=0.0
                  endif
               else
                  winteryy=0.0
                  wintrayy=0.0
               endif
               relativeyy=cprod(ipart)*winteryy
c
c XY terms
c
               if(iflagsubs(it1).NE.iflagsubs(it2)) then
                  winterxy=(winter(ipart,1)-winter(ipart,2)
     1-wtsumyydif*winteryy)/wtsumxy
                  wintraxy=(wintra(ipart,1)-wintra(ipart,2)
     1-wtsumyydif*wintrayy)/wtsumxy
               else if(iflagsubs(it1).eq.1.AND.iflagsubs(it2).eq.1
     1.and.iexchange(it1).eq.0.and.iexchange(it2).eq.0) then
                  winterxy=0.0
                  wintraxy=(wintra(ipart,1)-wintra(ipart,2)
     1-wtsumyydif*wintrayy)/wtsumxy
               else
                  winterxy=0.0
                  wintraxy=0.0
               endif
               relativexy=cprod(ipart)*winterxy
c
c XX terms
c
               if(iflagsubs(it1).eq.0.AND.iflagsubs(it2).eq.0) then
                  winterxx=(winter(ipart,1)-wtsumxy1*winterxy
     1-wtsumyy1*winteryy)/wtsumxx
                  wintraxx=(wintra(ipart,1)-wtsumxy1*wintraxy
     1-wtsumyy1*wintrayy)/wtsumxx
               else if(iflagsubs(it1).eq.1.AND.iflagsubs(it2).eq.1
     1.and.iexchange(it1).eq.0.and.iexchange(it2).eq.0) then
                  winterxx=0.0
                  wintraxx=(wintra(ipart,1)-wtsumxy1*wintraxy
     1-wtsumyy1*wintrayy)/wtsumxx
               else
                  winterxx=0.0
                  wintraxx=0.0
               endif
               relativexx=cprod(ipart)*winterxx
c
c ignore small contributions
c
               if(winteryy.ne.0.0.or.wintrayy.ne.0.0) then
                  write(12,114) partlabel(ipart)
                  write(12,*) wintrayy*rhofac,winteryy*rhofac,relativeyy
               endif
               if(winterxy.ne.0.0.or.wintraxy.ne.0.0) then
                  write(11,114) partlabel(ipart)
                  write(11,*) wintraxy*rhofac,winterxy*rhofac,relativexy
               endif
               if(winterxx.ne.0.0.or.wintraxx.ne.0.0) then
                  write(10,114) partlabel(ipart)
                  write(10,*) wintraxx*rhofac,winterxx*rhofac,relativexx
               endif
               endif
            end do
            close(10)
            close(11)
            close(12)
      endif
c
c write the .ato file in case anything was changed
c
      call moloutc(fnameato)
      return
      end

