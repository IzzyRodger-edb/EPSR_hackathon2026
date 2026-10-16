      parameter (mcoeff=1000,mpgauss=5)

      character*3 tlab1(mpairs),tlab2(mpairs)
      character*20 addpotential

      integer ncoeffp,npcoeff(mfile),npits,npitss,npfreq
     *,nzero,nnonzero,npairspcof,irefpcof(mpairs)
c Number of entries in energy table
     *,nergref
c number of Gaussians to be added to EP
     *,npgauss(mpairs)



      real pdecay(mfile),pstepd(mfile),pstepp(mfile)
     *,pwt(mcoeff),pwtold(mcoeff),expgr(mpz*mcoeff),expsq(mpz*mcoeff)
     *,pdmin,pdstep,pdmax,paccept,psigma2,psigma2read,psigmafac,pcoefd(mcoeff,mfile)
     *,pcoefp(mcoeff,mpairs),rminex(mpairs),coordex(mpairs)
     *,rmaxex(mpairs),pcoefex(mpairs),pcoefexmax(mpairs)
     *,potex(mpairs),potexmax(mpairs)
     *,weight(2,mfile,mpairs),uwt(mpairs),uwait(2,mfile+mpairs,mpairs)
     *,xav(mpz,mpairs),xav1(mpz,mpairs),wsum(mfile)
     *,xavs(mpz,mpairs),xavs1(mpz,mpairs),wfac(mpz),wfacc(mpz)
     *,xavtot(mpz,mpairs),xavtot1(mpz,mpairs)
     *,epot(mpz,mpairs),epot1(mpz,mpairs)
     *,ecoul(mpz,mpairs),ecoul1(mpz,mpairs)
     *,feedback,feedbackold,efacm(mfile),efacold(mfile),efacp(mpairs)
     *,rpgauss(mpgauss,mpairs),wpgauss(mpgauss,mpairs),hpgauss(mpgauss,mpairs)
     *,absolute_energy,potfacamp,potfac,refpotfac,ereqstep,erequnit
     *,coulfac(mpairs)

      logical gaussian,harmonic

      common/potential_c/tlab1,tlab2,addpotential

      common/potential_i/ncoeffp,npcoeff,npits,npitss,npfreq
     *,nzero,nnonzero,npairspcof,irefpcof,nergref,nneighb,npgauss

      common/potential_r/pdecay,pstepd,pstepp,pwt,pwtold
     *,expgr,expsq,pdmin,pdstep,pdmax,paccept
     *,psigma2,psigma2read,psigmafac,pcoefd,pcoefp,rminex,coordex
     *,rmaxex,pcoefex,pcoefexmax
     *,potex,potexmax,weight,uwt,uwait
     *,xav,xav1,wsum,xavs,xavs1,wfac,wfacc,epot,epot1,ecoul,ecoul1
     *,xavtot,xavtot1
     *,feedback,feedbackold,efacm,efacold,efacp
     *,rpgauss,wpgauss,hpgauss,absolute_energy,potfacamp,potfac,refpotfac,ereqstep,erequnit
     *,coulfac

      common/potential_l/gaussian,harmonic
