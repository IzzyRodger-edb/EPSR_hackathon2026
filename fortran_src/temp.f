c Set up the name of the output file

      fnameato=fnameato(1:indexlastdot(fnameato)-1)//'.mol'
      ilenfnameato=indexlastdot(fnameato)+3

      open(10,file=fnameato,status='unknown',iostat=ierr)
      if(ierr.ne.0) then

         write(6,501) fnameato(1:ilenfnameato)
 501  format('readmole> Output file ',a,' could not be opened')
         return

      endif

c Signal that this is a gmol file and that the atom numbers listed 
c are the EPSR atom numbers

      write(10,*) '.gmol ',0


c Step through the atoms and set up bond distances.

        ia1=0
        do while (ia1.lt.natom)

           ia1=ia1+1
           ib1=0
           do while (ib1.lt.nbond(ia1))

              ib1=ib1+1
              xdif=xa(iabond(ia1,ib1))-xa(ia1)
              ydif=ya(iabond(ia1,ib1))-ya(ia1)
              zdif=za(iabond(ia1,ib1))-za(ia1)
              dsq=xdif*xdif+ydif*ydif+zdif*zdif
              tbond(ia1,ib1)=sqrt(dsq)

           end do
        end do

c Write out the atom coordinates and bond pairs,

       if(labeltype.ne.3) then
        do ia1=1,natom
           write(10,*) 'atom ',iatomlabel(ia1),' ',atomlabel(ia1),' '
     *,xa(ia1),ya(ia1),za(ia1),nbond(ia1)
     *,(iatomlabel(iabond(ia1,j)),j=1,nbond(ia1))
        end do
       else if(labeltype.eq.3) then
        do ia1=1,natom
           write(int_as_str1,'(I3.3)') iatomlabel(ia1)
           write(10,*) 'atom ',iatomlabel(ia1),' '//int_as_str1//' '
     *,xa(ia1),ya(ia1),za(ia1),nbond(ia1)
     *,(iatomlabel(iabond(ia1,j)),j=1,nbond(ia1))
        end do
       endif


c Set up the bond list. 
c Where multiple bonds for a given atom pair are defined,
c the mean bond distance is calculated

        nbondlist=0

        do ia1=1,natom

           ib1=0
           do while (ib1.lt.nbond(ia1))
              
              ib1=ib1+1
            ia2=iabond(ia1,ib1)
              if(ia2.gt.ia1) then

                 if(atomtype(ia1).le.atomtype(ia2)) then
                if(labeltype.ne.3) then
                bondlabel=atomlabel(ia1)//' '//atomlabel(ia2)
                else if(labeltype.eq.3) then
                write(int_as_str1,'(I3.3)') iatomlabel(ia1)
                write(int_as_str2,'(I3.3)') iatomlabel(ia2)
                bondlabel=int_as_str1//' '//int_as_str2
                end if

             else

                if(labeltype.ne.3) then
                bondlabel=atomlabel(ia2)//' '//atomlabel(ia1)
                else if(labeltype.eq.3) then
                write(int_as_str1,'(I3.3)') iatomlabel(ia1)
                write(int_as_str2,'(I3.3)') iatomlabel(ia2)
                bondlabel=int_as_str2//' '//int_as_str1
                end if

             endif

c Search the bond list to see if this label has already been defined!{{{

                 ib2=0
                 test=.true.
                 do while (ib2.lt.nbondlist.and.test)

                    ib2=ib2+1
                    test=bondlabellist(ib2).ne.bondlabel

                 enddo

                 if(test.and.nbondlist.lt.mstarlabel) then

                    ib2=nbondlist+1
                    bondlabellist(ib2)=bondlabel
                    bondlist(ib2)=0.0
                    nbondlabel(ib2)=0
                    nbondlist=ib2
                   
                 endif

                 bondlist(ib2)=bondlist(ib2)+abs(tbond(ia1,ib1))
                 nbondlabel(ib2)=nbondlabel(ib2)+1

              endif
           end do
        end do

c Write out the bonds

        ib1=0
        do while (ib1.lt.nbondlist)

           ib1=ib1+1
           write(10,*) 'bond ',bondlabellist(ib1),' '
     +,bondlist(ib1)/nbondlabel(ib1)

        end do

c Setup the bondangle list

c Step through the atoms and set up the angles between bonded pairs

        pi=4.0*atan(1.0)
        piconv=180/pi
        nbondanglelist=0

        ia1=0
        do while (ia1.lt.natom)

           ia1=ia1+1
           ib1=0
           do while (ib1.lt.nbond(ia1)-1)

              ib1=ib1+1
              iaref1=iabond(ia1,ib1)
              a(1)=xa(iaref1)-xa(ia1)
              a(2)=ya(iaref1)-ya(ia1)
              a(3)=za(iaref1)-za(ia1)
              ib2=ib1
              do while(ib2.lt.nbond(ia1))

                 ib2=ib2+1
                 iaref2=iabond(ia1,ib2)
                 b(1)=xa(iaref2)-xa(ia1)
                 b(2)=ya(iaref2)-ya(ia1)
                 b(3)=za(iaref2)-za(ia1)

                 angle=piconv*acos(scapro(a,b)
     */sqrt(scapro(a,a)*scapro(b,b)))!}}}

c Setup the bondangle label

             if(atomtype(iaref1).le.atomtype(iaref2)) then

                if(labeltype.ne.3) then
                bondanglelabel=atomlabel(iaref1)
     +//' '//atomlabel(ia1)
     +//' '//atomlabel(iaref2)
                else if(labeltype.eq.3) then
                write(int_as_str1,'(I3.3)') iatomlabel(iaref1)
                write(int_as_str2,'(I3.3)') iatomlabel(ia1)
                write(int_as_str3,'(I3.3)') iatomlabel(iaref2)
                bondanglelabel=int_as_str1//' '//int_as_str2//' '//int_as_str3
                end if

             else

                if(labeltype.ne.3) then
                bondanglelabel=atomlabel(iaref2)
     +//' '//atomlabel(ia1)
     +//' '//atomlabel(iaref1)
                else if(labeltype.eq.3) then
                write(int_as_str1,'(I3.3)') iatomlabel(iaref2)
                write(int_as_str2,'(I3.3)') iatomlabel(ia1)
                write(int_as_str3,'(I3.3)') iatomlabel(iaref1)
                bondanglelabel=int_as_str1//' '//int_as_str2//' '//int_as_str3
                end if

             endif

c Step through the existing list of labels to see if any exist already

             ib3=0
             test=.true.
             do while(ib3.lt.nbondanglelist.and.test)

                    ib3=ib3+1
                test=bondanglelabellist(ib3).ne.bondanglelabel

             end do

                 if(test.and.nbondanglelist.lt.mstarlabel) then

                    ib3=nbondanglelist+1
                    bondanglelabellist(ib3)=bondanglelabel
                     bondanglelist(ib3)=0.0
                    nbondanglelabel(ib3)=0
                    nbondanglelist=ib3
                   
                 endif

                 bondanglelist(ib3)=bondanglelist(ib3)+angle
                 nbondanglelabel(ib3)=nbondanglelabel(ib3)+1

              end do
           end do
        end do

c Write out list of angles found

      ib1=0
      do while (ib1.lt.nbondanglelist)

         ib1=ib1+1
         write(10,*) 'angle ',bondanglelabellist(ib1)
     +,bondanglelist(ib1)/nbondanglelabel(ib1)
      
      end do

c Search for all possible dihedral angles

        ndihedrals=0

        ia1=0
        do while (ia1.lt.natom)

           ia1=ia1+1
           ib1=0
           do while (ib1.lt.nbond(ia1))

              ib1=ib1+1
              ia2=iabond(ia1,ib1)

c Define 2 axis atoms

              if(ia2.gt.ia1) then

                 b(1)=xa(ia2)-xa(ia1)
                 b(2)=ya(ia2)-ya(ia1)
                 b(3)=za(ia2)-za(ia1)

c Get neighbours of ia1

                 ib2=0
                 do while (ib2.lt.nbond(ia1))

                    ib2=ib2+1
                    ia1ref=iabond(ia1,ib2)
                    if(ia1ref.ne.ia2) then

                       a1(1)=xa(ia1ref)-xa(ia1)
                       a1(2)=ya(ia1ref)-ya(ia1)
                       a1(3)=za(ia1ref)-za(ia1)

c Calculate the vector perpendicular to the plane formed by these two vectors

                       call vecpro(a1,b,c1perp)

c Normalise this vector

                       amp=sqrt(scapro(c1perp,c1perp))
                       do i=1,3
                          c1perp(i)=c1perp(i)/amp
                       end do

c Now search the bonds of ia2 for the 4th atom

                       ib3=0
                       do while (ib3.lt.nbond(ia2))

                          ib3=ib3+1
                          ia2ref=iabond(ia2,ib3)
                          if(ia2ref.ne.ia1.and.ia2ref.ne.ia1ref) then

                             a2(1)=xa(ia2ref)-xa(ia2)
                             a2(2)=ya(ia2ref)-ya(ia2)
                             a2(3)=za(ia2ref)-za(ia2)

c Calculate the vector perpendicular to the plane formed by these two vectors

                             call vecpro(a2,b,c2perp)

c Normalise this vector

                             amp=sqrt(scapro(c2perp,c2perp))
                             do i=1,3
                                c2perp(i)=c2perp(i)/amp
                             end do

c Cosine of the angle between them

                             cosphi=scapro(c1perp,c2perp)
                             if(cosphi.gt.1.0) cosphi=1.0
                             if(cosphi.lt.-1.0) cosphi=-1.0

c Sine of the angle between them - see torangles for this

                             call vecpro(c1perp,c2perp,csign)
                             sinphi=scapro(csign,b)/sqrt(scapro(b,b))
                             angle=piconv*atan2(sinphi,cosphi)

                             ndihedrals=ndihedrals+1
                             if(ndihedrals.gt.mstarlabel) then
                                write(6,505) mstarlabel
 505                            format(
     +'readmole> No of dihedrals exceeds ',i5,' Program aborted')
                                return
                             endif
                             write(10,*) 'dihedral ',iatomlabel(ia1ref)
     *,iatomlabel(ia1),iatomlabel(ia2),iatomlabel(ia2ref),angle
                          endif
                       end do
                    endif
                 end do
              endif
           end do
        end do

c Write the initial potential parameters for each atom type

        do ia1=1,ilenatomtypelist

            

           if(chsymbol(atomicnolist(ia1)).eq.'H ') then

              epsilon=0
              sigma=0
              amass=2

           else if(chsymbol(atomicnolist(ia1)).eq.'O ') then

              epsilon=0.65
              sigma=3.1
              amass=chatwt(atomicnolist(ia1))

           else if(chsymbol(atomicnolist(ia1)).eq.'N ') then

              epsilon=0.7
              sigma=3.2
              amass=chatwt(atomicnolist(ia1))
           
           else if(chsymbol(atomicnolist(ia1)).eq.'C ') then

              epsilon=0.8
              sigma=3.7
              amass=chatwt(atomicnolist(ia1))
           
           else

              epsilon=0.8
              sigma=3.2
              amass=chatwt(atomicnolist(ia1))
           
           endif
           write(10,600) atomlabellist(ia1)
     *,epsilon,sigma,amass,chargelist(ia1)
     *,chsymbol(atomicnolist(ia1))
600   format(1x,'potential ',a3,1x,4(f10.5,1x),a3)

        end do

        write(10,*) 'temperature ',simtk
        write(10,*) 'vibtemp ', tempb
        write(10,*) 'angtemp ', tempang
        write(10,*) 'dihtemp ', tempdih
        write(10,*) 'ecoredcore ', ecore,dcore
        if(unitcelldefined) then
           write(10,*) 'unitcell ',alat,blat,clat,phib,thetac,phic
        else
           write(10,*) 'density ', rho
        endif
        close(10)

c Finally convert this to an ato file

        ilenfname_ato=indexlastdot(fname_ato)-1
        sometext='makemole '//fname_ato(1:ilenfname_ato)//'.mol'
        call parse(sometext,nwords,ncf,ncl,200,nstsize,nworddim)
        call makemole()


