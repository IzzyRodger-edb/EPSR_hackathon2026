c
c The following bit was used to control the amplitude of the potential coefficients.
c Now not used.
c

c Now adjust the potential amplitude if so requested

c We attempt to modify potential coefficients only if requrested and we are a reasonable no. of iterations
c from the last intramolecular disordering - otherwise the energy/pressure are not at their equilibrium values
c and we may distort the EP incorrectly.

      if(controlep.and.nisum.ge.nisumm) then

c Set a limit to the epfactor

         epfactorlimit=0.1*potfac
         epfactorupperlimit=1.0
         epfactorlowerlimit=1.0-epfactorlimit
         if(epfactorupperlimit.gt.1.1) epfactorupperlimit=1.1
         if(epfactorlowerlimit.lt.0.9) epfactorlowerlimit=0.9

c Energy and pressure from reference plus repulsive potentials alone

         energrp=energd-energep
         pressrp=pressd-pressep
         epchangerp=energrp*requestew+pressrp*requestpw

c Total weighted sum of requested energy and pressure

         epchangerq=requeste*requestew+requestp*requestpw

c Difference between these quantities
      
         epchangedif=epchangerq-epchangerp
 
c Actual difference that IS delivered by current EP

         epchangeep=energep*requestew+pressep*requestpw
         write(6,*) nisum,nisumm,epchangerp,epchangerq,epchangedif,epchangeep
         epfactor=1.0

c Can only change coefficients by a factor if the sign of the requested change and that of the current EP 
c are the same

         if(abs(epchangeep).gt.0.0.and.sign(1.0,epchangeep).eq.sign(1.0,epchangedif)) then

            epfactor=(epchangedif)/(epchangeep)
            if(epfactor.gt.epfactorupperlimit) epfactor=epfactorupperlimit
            if(epfactor.lt.epfactorlowerlimit) epfactor=epfactorlowerlimit
            write(6,*) ' Potential factor = ',epfactor

c In case there is no chance of energy/pressure from EP being of the correct sign, or else if the coefficients
c are not changed, allow perturbation of EP

            if(epfactor.eq.1.0) then
               changedpcof=.false.
            else
               changedpcof=.true.
            endif     

c Process individual site-site terms depending on the sign of the current EP energy and pressure

            do ip=1,npairs

               ipref=irefpcof(ip)
c               epfact(ipref)=1.0
c               epchangeept=energf(ipref,2)*requestew+pressf(ipref,2)*requestpw
c Enhance coefficients that will improve the energy or pressure, and reduce those which won't.
c               if(sign(1.0,epchangeept).eq.sign(1.0,epchangerq)) then
               epfact(ipref)=epfactor
c              else
c                  epfact(ipref)=1.0/epfactor
c               endif
            end do
            call multiply_pot_coeffs(epfact)
         endif
      endif

