      parameter (ibar2=1602,jbar2=290,meshx=20,meshy=20,
     &           nobd=200,nfrsrfd=20,nprd=200)
      parameter (mshx=meshx+1, mshy=meshy+1, nxy=ibar2*jbar2)

      character*80 prbname
      character*8 dat,jnm,ochn,tim
      character*10 fdir
      character*2 idt
      logical sym

      real*8 nxo,mxo,nyo,myo,nxf,mxf,nyf,myf,jcb,mod1,
     &		 nxp,mxp,nyp,myp

      common /begini/ aa,h0,h0r,xcenter,alpha1,c1,
     1		xstart,tstart,tfinish,ybc(2000,40),ubc(2000,40),
     1      tstart_a,tfinish_a,prtdt_a,prtdt_t,twprt_a,
     2      vbc(2000,40),xxxf,ttend,xxl,xxt,xxk,
     3		xxlr,xxkr,c1r,crest,energyp0,
     3		segma,fxstart,t_pad(2000),u_pad(2000),h_pad(2000),
     4		x_pad(2000),ulinear,usecond,areturn,stn

	  common /addition2/ det_k(100),ytrough,mod1

	  common /addition3/ j_det,istrflag,nextra,
     1	nrs,nloc,nanimation,nout(100)

      common /addition4/ xout(100)

      common /addition14/ xpol(10),ypol(10),conc(10)

	  common /addition0/ ttmax,yjet,ujet,vjet,wjet,pjet,ylow,yup,xxc,
     3          time_jet,yout,uout,vout,wout,ylowout,yupout,
     4		ustar,xdis,eta0old,eta1old,eta2old,eta2aold,eta0x,eta2x,
     5		etaimax,c01,c02,cim1,crbc,xforce,yforce,
     5		eta0,eta1,eta1a,eta2,uleft1,vleft1,xmass0,xmass_d,
     5		xsaka,xmass_a,xmass_a1,xmass_a2,cnf,petit,
     5		umean,ratio,uturb,xovertop,overtopmass,
     5		ssource,tsource,aawave(100),twave(100),cwave(100),
     5		swave(100),nwave,isources,isourcee,jsources,jsourcee,
     5		nsource,islip,ninflow,ibg,ieg,jbg,jeg,
     6		it,nnn,istable,interx,intery,nexit,
     7		istablemx,kemodel,npollutant,nweakref,nopen,novertop
c
      common /fldvarr/ un(nxy),vn(nxy),pn(nxy),fn(nxy),ar(nxy),at(nxy),
     4				dissipturb(nxy),dissipmole(nxy),productturb(nxy),
     5				xmnf(nxy),xmnu(nxy),xmnv(nxy),xmnk(nxy),uxmb(10),
     6				vymb(10),xnewloc,ynewloc,acmin,axmb,aymb,tpmb,nmb

	  common /addition1/ ac(nxy),cvol(nxy),tss(nxy),
     5                 tauxx(nxy),tauyy(nxy),tauxy(nxy),ttots(nxy),
     5                 tauxxturb(nxy),tauyyturb(nxy),tauxyturb(nxy),
     6                 ftilde(nxy),u(nxy),v(nxy),p(nxy),f(nxy),xpn(nxy),
     7                 xnut(nxy),xnutn(nxy),xk(nxy),xkn(nxy),xp(nxy),
     8				   xnuty(nxy),xnutyn(nxy),xep(nxy),xepn(nxy),
     9		           xprod(nxy),amp1,amp2,amp3,amp4,amp5,realize,
     2				   eddycoef,c_mu_bc,XE,XEpor,c1e,c2e,c_mu,sege,segk,
     9				   xmaxxnut,ymaxxnut,c1ea,c1eb,ticf,roughness,xmnt
c
      common /fldvari/ nf(nxy),nfold(nxy),nfs(nxy)
c
      common /geomr/ jcb(nxy),alp(nxy),gam(nxy)
c
      common /meshr/ x(ibar2),xi(ibar2),delx(ibar2),etah(ibar2),
     1             etatmp(ibar2),rdx(ibar2),rx(ibar2),y(jbar2),
     2			   yj(jbar2),dely(jbar2),rdy(jbar2),xl(mshx),
     3             xc(meshx),dxmn(meshx),yl(mshy),yc(meshy),
     4             dymn(meshy),r(ibar2),ri(ibar2),
     5             xb(5),yb(5)
c
      common /meshi/ nxl(meshx),nxr(meshx),nyl(meshy),nyr(meshy)
c
      common /intvarr/ delt,t,autot,prtdt,twprt,pltdt,twplt,
     1             sfdt,sfprt,
     1             twfin,xnu,dtend,dmpdt, rhof, vchgt, ttl,
     2             sigma, cyl, gx, gy, ui, vi, utop,
     3             alpha, beta, flgc, xmin, xmax, ymin,
     4             omcyl, ymax,tmv,tke,power,xsponge,adamp,
     5             dtvis, dudr,dudl,dudt,dudb,dvdr,dvdl,dvdt,dvdb,
     6             tquit,tbeg,dxmin,dymin,psat,tseti,tsetf,
     7             dtmax,uinf(4),vinf(4),con,fcvlim,
     8             erriccg,xmu,opalp,omalp,frctn,twdmp

      common /intvari/ ibar,jbar,imax,jmax,im1,jm1,nkx,nky,ncyc,
     1             icyl,kl,kr,kb,kt,iter,nocon,nflgc,
     3             i,j,liter,npack,idiv,ibcflg0,
     4             ibcflg,ibcfinal,itc,jtc,iotty,
     5             ivis,jvis,isft,jsft,itmxiccg,sym

      common /constr/ emf,emf1,em6,em10,ep10,pi,tpi,rpd,
     1                em6p1,em61,fact
c
      common /obstr/ oa2(nobd),oa1(nobd),ob2(nobd),ob1(nobd),
     1               oc2(nobd),oc1(nobd),xobs(nxy),
     2               yobs(nxy),od1(nobd),nxo(nobd),od2(nobd),
     3               mxo(nobd),oe1(nobd),nyo(nobd),oe2(nobd),myo(nobd)

      common /porousmd/ pa2(nprd),pa1(nprd),pb2(nprd),pb1(nprd),
     1           pc2(nprd),pc1(nprd),
     2           pd1(nprd),nxp(nprd),pd2(nprd),ipr(nprd),
     3           mxp(nprd),pe1(nprd),nyp(nprd),pe2(nprd),
     4		     myp(nprd),xporosity(100),porousa(nxy),porousb(nxy),
     8			 porousc(nxy),porousd(nxy),porousp(nxy),
     9			 porstart,porheight,porlength(50),porslope(50),
     5		     xalpha(100),xbeta(100),xgamma(100),d50(100),
     6		     xa(100),xxb(100),gc(100),npc(nxy),nporous(100),
     7		     nmovbd(nxy),nportype,npor,nporsection

      common /obsti/ nmean,istart,iend,iinterval,
     1		nobstype,nobs(100),noc(nxy),ioh(nobd),ijobs(nxy)
c
      common /mean/ tmstart,tmend,tinterval,
     1	eflux(ibar2),efluxt(ibar2),disp(ibar2),ep(ibar2),ep0(ibar2),
     1	ek(ibar2),et(ibar2),etaa(ibar2)
c
      common /frsrfr/ fa2(nfrsrfd),fa1(nfrsrfd),
     1                fb2(nfrsrfd),fb1(nfrsrfd),fc2(nfrsrfd),
     2                fc1(nfrsrfd),ifh(nfrsrfd),fd1(nfrsrfd),
     3                nxf(nfrsrfd),fd2(nfrsrfd),mxf(nfrsrfd),
     4                fe1(nfrsrfd),nyf(nfrsrfd),fe2(nfrsrfd),
     5                myf(nfrsrfd),flht,frsurf
c
      common /frsrfi/ nfrsrf
c
      common /ghostr/ pbc(4),pbcl,pbcr,pbct,pbcb
c
      common /labeli/ prbname,dat,jnm,ochn,tim,idt,fdir
c
      common /iccgr/ a(nxy),b(nxy),c(nxy),d(nxy),e(nxy),
     1               af(nxy),bf(nxy),cf(nxy),df(nxy),ef(nxy),
     2               ag(nxy),bg(nxy),cg(nxy),dg(nxy),eg(nxy),
     3               ah(nxy),bh(nxy),ch(nxy),dh(nxy),eh(nxy),
     4               soln(nxy),srce(nxy),residu(nxy),sc1(nxy),
     6               sc2(nxy)

c.....For Youngs model
      logical periodic, four_sweep, tracer, slic
      integer tracers

      common /intrfcer/ tol, cutvof, eps, epsmach,
     &                  avgiter, vol_error,
     &                  advect_vol(nxy),utt(nxy),vrr(nxy),ftmp(nxy),
     &                  grad_vof(2,nxy),
     &                  flux_pt(2,nxy), dvol(nxy),
     &                  xi0(nxy), xi1(nxy), yi0(nxy), yi1(nxy),
     &                  cell_vol(nxy), dvol_tot(2,nxy)

      common /intrfcei/ itmax0, nmix, 
     &                  idir, periodic, iface0(nxy),
     &                  iface1(nxy), ilook(4), ijmix(nxy), four_sweep,
     &                  isweeps, frac, io_interval, tracers, tracer,
     &                  slic,nfree





