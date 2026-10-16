      character*3 alabel(matom),type(mcomp),atsymbol(mcomp)
c type labels of atoms in current file
     *,typefile(mcomp)

      real tol,tolout,step,tempf,tempb,tempc,stepf,stepd,steps
     +,tempdih,tempang,tempdihreq,tempangreq,tempangwt,tempdihwt
     +,stepr,stepri,stepm,stepmi,xa(matom),ya(matom),za(matom)
     +,sxa(matom),sya(matom),sza(matom),amassa(mcomp),amassr(mpairs)
     +,tbond(matom,near),bound(matom),eatomt(matom),eatomv(matom)
     +,xs,ys,zs,fold,trial,trialt(matom),trialv(matom)
     +,simt,simtk,charge(mcomp),qradius(mcomp),qcor(mpairs)
     +,dcor(mpairs),sigcor(mcomp),epscor(mcomp)
     +,res,sig1m,sig2m,scaleb,distob,dscor(mpairs),ecor(mpairs)
     +,xbar(matom),ybar(matom),zbar(matom),dscorm(mpairs)
     +,x2bar(matom),y2bar(matom),z2bar(matom),sig1b,sig2b,sig3b
     +,xid(matom),yid(matom),zid(matom),scalef,distor
     +,sig1b2,sig2b2,sig3b2,scal2b,dist2b,cfrac(mcomp),cfraccg(mcomp) !The latter takes account of any coarse graining 
     +,dcore,ecore,dcorm,dcorsq,dcorsm,power
     +,sizefactor,decline_rate,dmaxmol,refine_intra_fac,coul_intra_fac
     +,sumcharge
     +,dihedralang(near,matom) !Records the required dihedral angle (in deg.) for rotational groups with the label 'DIH'
     +,maxmolsize !Stores the maximum dimension of all molecules (reduced coordinates)

      integer labrotdih(2,near,matom),ndihedralang(matom)
     +,nearn(matom),nbondp(matom),nbonds(matom),nbond(matom),iabond(matom,near)
     +,ntype,itype(matom)
     +,npair,iadd(matom),ntry,ntryo,nbar,iok(matom),natom,natomr
     +,mref(matom)
     +,qtype(matom),qtypea(mcomp),npairs,npairref(mcomp,mcomp)
     +,n1extr(mpairs),n2extr(mpairs),iexchange(mcomp),ntypefile,iatomicnotype(mcomp)
     +,natomsmoved,atomsmoved(near) !The number of atoms within a molecule moved when changing internal energy
     +,imageatom(matom) !Image atom for long molecules that need to connect across box boundary

      logical refine_intra,coul_intra

      common/molecom_r/tol,tolout,step,tempf,tempb,tempc,stepf,stepd,steps
     +,tempdih,tempang,tempdihreq,tempangreq,tempdihwt,tempangwt
     +,stepr,stepri,stepm,stepmi,xa,ya,za
     +,sxa,sya,sza,amassa,amassr
     +,tbond,bound,eatomt,eatomv
     +,xs,ys,zs,fold,trial,trialt,trialv
     +,simt,simtk,charge,qradius,qcor
     +,dcor,sigcor,epscor
     +,res,sig1m,sig2m,scaleb,distob,dscor,ecor
     +,xbar,ybar,zbar,dscorm
     +,x2bar,y2bar,z2bar,sig1b,sig2b,sig3b
     +,xid,yid,zid,scalef,distor
     +,sig1b2,sig2b2,sig3b2,scal2b,dist2b,cfrac,cfraccg
     +,dcore,ecore,dcorm,dcorsq,dcorsm,power
     +,sizefactor,decline_rate,dmaxmol,refine_intra_fac,coul_intra_fac
     +,sumcharge
     +,dihedralang
     +,maxmolsize

      common/molecom_c/alabel,type,atsymbol,typefile

      common/molecom_i/labrotdih,ndihedralang
     +,nearn,nbondp,nbonds,nbond,iabond
     +,ntype,itype
     +,npair,iadd,ntry,ntryo,nbar,iok,natom,natomr
     +,mref
     +,qtype,qtypea,npairs,npairref
     +,n1extr,n2extr,iexchange,ntypefile,iatomicnotype
     +,natomsmoved,atomsmoved,imageatom

      common/molecom_l/refine_intra,coul_intra

