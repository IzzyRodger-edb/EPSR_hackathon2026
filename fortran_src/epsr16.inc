      real x(mpz),xstar(mpz),y(mpz),ycoh(mpz),entwt(mpz)
     *,foqols(mpz,mpts),xsav(mpz),ysav(mpz),xavbar,rlog
     *,soqinter(mpts,mpairs),sq2inter(mpts,mpairs)
     *,soqintra(mpts,mpairs),sq2intra(mpts,mpairs)
     *,add,cellst,sigma1(mfile),sigma2(mfile),addcon(mfile)
     *,rhofac(mfile),xwt(mfile),cell,datfc2,sm1sav(mpz)
     *,rmaxcell,cellmax,upol,rho,amubar,amubar2,rejrate
     *,qminm,qminp
     
      common/dista/x,xstar,y,ycoh,entwt
     *,foqols,xsav,ysav,xavbar,rlog
     *,soqinter,sq2inter
     *,soqintra,sq2intra
     *,add,cellst,sigma1,sigma2,addcon
     *,rhofac,xwt,cell,datfc2,sm1sav
     *,rmaxcell,cellmax,upol,rho,amubar,amubar2,rejrate
     *,qminm,qminp

      real gorinter(mpz,mpairs),gr2inter(mpz,mpairs)
     *,gorintra(mpz,mpairs),gr2intra(mpz,mpairs)
     *,energyref(mpairs),energyemp(mpairs),energytot(mpairs)

      common/xaver/gorinter,gr2inter
     *,gorintra,gr2intra
     *,energyref,energyemp,energytot
     
      integer npt,npt1,npt2,nmove,nfile,ndist,nqmin(mfile)
     *,nqmax(mfile),nqpts(mpts),nptmin,nptmax
     
      common/idista/npt,npt1,npt2,nmove,nfile,ndist,nqmin
     *,nqmax,nqpts,nptmin,nptmax

      real :: entsig(mpz,mpairs)

      common/xaver8/entsig


