      subroutine rinput
c
c ======================================================================
c
c   Purpose -
c     read the input data; initialize selected and default variables
c
c   RINPUT is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c           SETUP   
c
c
c   RINPUT calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c         SETARRY    ripple
c
c ======================================================================
c
c##############################################################
      implicit real*8 (a-h,o-z)
      
c##############################################################
c
c############
      include  "comdk2.h" 
c############
c
      external ck,ce,cn
	dimension tmp(50)
c
c.... Namelist input
c
      namelist /numparam/ delt,twfin,prtdt,pltdt,sfdt,
     1                    alpha,kl,kr,kt,kb,autot,
     2                    npack,con,dmpdt,
     3                    dtmax,idiv,
     4                    erriccg,fcvlim,
     5                    frctn,
     6                    itmxiccg,sym
      namelist /fldparam/ xnu,icyl,gx,gy,ui,vi,utop,psat,
     1                    rhof,uinf,vinf,pbc
      namelist /mesh/ nkx,xl,xc,nxl,nxr,dxmn,nky,yl,yc,nyl,nyr,dymn
      namelist /obstcl/ nobs,oa2,oa1,ob2,ob1,oc2,oc1,ioh,od1,od2,
     &                  oe1,oe2,nxo,mxo,nyo,myo
C.....add to deal with multiple porous medias
      namelist /porous/ nporous,pa2,pa1,pb2,pb1,pc2,pc1,ipr,pd1,pd2,
     &                  pe1,pe2,nxp,mxp,nyp,myp,nportype,xporosity,
     &                  xalpha,xbeta,d50,xgamma
      namelist /freesurf/ nfrsrf,fa2,fa1,fb2,fb1,fc2,fc1,ifh,fd1,fd2,
     &                    fe1,fe2,nxf,mxf,nyf,myf,
     &                    flht
c
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... Branch around defaults if restarting
c
	character*2 AAA,achac(30)
        fdir='./Results/'
        achac(1)='1'
        achac(2)='2'
        achac(3)='3'
        achac(4)='4'
        achac(5)='5'
        achac(6)='6'
        achac(7)='7'
        achac(8)='8'
	  achac(9)='9'
	  achac(10)='10'
        achac(11)='11'
        achac(12)='12'
        achac(13)='13'
        achac(14)='14'
        achac(15)='15'
        achac(16)='16'
        achac(17)='17'
        achac(18)='18'
        achac(19)='19'
        achac(20)='20'
        achac(21)='21'
        achac(22)='22'
        achac(23)='23'
        achac(24)='24'
        achac(25)='25'
        achac(26)='26'
        achac(27)='27'
        achac(28)='28'
        achac(29)='29'
        achac(30)='30'

	rsflg=0.0
	nnn=0
	fact=1.0
	t0=0.0

      if (rsflg.ne.0.0) go to 60
      rsflg = 1.0d0
c
c.... Set defaults for input variables
c
c.... Namelist /numparam/
c
      con = 0.30d0
      dtmax = 0.5d0
      alpha = 0.3d0
      autot = 1.0d0
      erriccg = 1.0d-08
      fcvlim = 0.39d0
      frctn = 5.0d-03
      npack = 0
	dmpdt=0.0
      idiv = 1
      kl = 1
      kr = 1
      kt = 1
      kb = 1
      itmxiccg = 2000
      sym = .true.
	prtdt=0.05
	pltdt=0.05
c
c.... Namelist /fldparam/
c
       rhof = 1.0d0
       psat = 0.0d0
       xnu = 0.0d0
       gx = 0.0d0
       gy = 0.0d0
       ui = 0.0d0
       vi = 0.0d0
       utop=0.0d0
       icyl = 0
c
c.... Namelist /obstcl/
c
      nobstype = 0
c
c.... Namelist /freesurf/
c
      flht = 0.0d0
      nfrsrf = 0
c
c.... Initialize constants
c
      emf = 1.0d-06
      em6 = 1.0d-06
      em10 = 1.0d-10
      ep10 = 1.0d+10
      tquit = 60.d0
      em6p1 = 1.000001d0
      pi = 3.14159265359d0
      rpd = 0.0174532925d0
      tbeg = 0.d0
      tpi = 6.283185307d0
      em61 = 0.999999d0
c
c.... Initialize selected variables
c
      t = 0.d0
      flgc = 0.d0
      twdmp = 0.d0
      twplt = 0.d0
      twprt = 0.d0
      sfprt = 0.d0
      vchgt = 0.d0
      pbcl = 0.d0
      pbcr = 0.d0
      pbct = 0.d0
      pbcb = 0.d0
      iter = 0
      ncyc = 0
      nflgc = 0
      nocon = 0
      nvar = 18
      ibcflg = 0
	petit=0.0
	nextra=0
	umax=1.0e6
	vmax=1.0e6
	istablemx=10
      j_det=0
      xdis=0.0
      xmass_a1=0.0
      xmass_a2=0.0
      xmass_d=0.0
      eta0x=0.0
	xsaka=0.0
C.....ratio of specifying initial turbulence level
	ticf=1.0
	crest=1.0d16
C.....nonlinear eddy viscosity coef for turbulence realizability (default: disable)
	c1ea=0.0
	c1eb=0.0
      realize=2.0
C.....transient turbulence adjustment coef. (if large, no adjustment)
C	eddycoef=3.4
	eddycoef=50.0

	XE=9.0d0
	XEpor=9.0d0
	kemodel=4
	nweakref=0
	islip=0
	xovertop=1.0e16

C.....Default value for free surface tracking method 
C.....nfree=0: no free surface updating; nfree=1: donar-acceptor; nfree=4: Youngs method
	nfree=1

C.....moving boundary flag
	nmb=1
	xnewloc=0.0
	ynewloc=0.0
	acmin=0.1
	do n=1,10
		uxmb(n)=0.0
		vymb(n)=0.0
	end do
c
c.... set default array values
c
c     ---------------------------------------
      call setarry (uinf(1),0.0d0,4)
      call setarry (vinf(1),0.0d0,4)
      call setarry (pbc(1),0.0d0,4)
      call setarry (oa1(1),0.0d0,nobd)
      call setarry (oa2(1),0.0d0,nobd)
      call setarry (ob1(1),0.0d0,nobd)
      call setarry (ob2(1),0.0d0,nobd)
      call setarry (oc1(1),0.0d0,nobd)
      call setarry (oc2(1),0.0d0,nobd)
      call setarry (od1(1),0.0d0,nobd)
      call setarry (od2(1),0.0d0,nobd)
      call setarry (oe1(1),0.0d0,nobd)
      call setarry (oe2(1),0.0d0,nobd)
      call setarry (nxo(1),0.0d0,nobd)
      call setarry (mxo(1),0.0d0,nobd)
      call setarry (nyo(1),0.0d0,nobd)
      call setarry (myo(1),0.0d0,nobd)
      call setarry (fa1(1),0.0d0,nfrsrfd)
      call setarry (fa2(1),0.0d0,nfrsrfd)
      call setarry (fb1(1),0.0d0,nfrsrfd)
      call setarry (fb2(1),0.0d0,nfrsrfd)
      call setarry (fc1(1),0.0d0,nfrsrfd)
      call setarry (fc2(1),0.0d0,nfrsrfd)
      call setarry (fd1(1),0.0d0,nfrsrfd)
      call setarry (fd2(1),0.0d0,nfrsrfd)
      call setarry (fe1(1),0.0d0,nfrsrfd)
      call setarry (fe2(1),0.0d0,nfrsrfd)
      call setarry (nxf(1),0.0d0,nfrsrfd)
      call setarry (mxf(1),0.0d0,nfrsrfd)
      call setarry (nyf(1),0.0d0,nfrsrfd)
      call setarry (myf(1),0.0d0,nfrsrfd)
c     ---------------------------------------
c
c.... write out run time info
c
   60 read(5,400) prbname
      write(9,25)
      write(13,25)
      write(9,50) dat,tim,ochn
      write(13,50) dat,tim,ochn
      write(9,100) prbname
      write(13,100) prbname
      write(9,150)
c
c...  read initial input data; write it out to paper and film
c
        read(5,*)
        read(5,*)delt,twfin,prtdt,autot
        read(5,*)kl,kr,kt,kb
        read(5,*)
        write(13,numparam)
	  dtmax=min(dtmax,delt*2.5)

        read(5,*)
        read(5,*)xnu,gx,gy,ui,vi,rhof
        read(5,*)
        write(13,fldparam)

        read(5,*)
        read(5,*)nkx,(xl(i),i=1,nkx+1),(xc(i),i=1,nkx),(nxl(i),
     &  i=1,nkx),(nxr(i),i=1,nkx)
        read(5,*)(dxmn(i),i=1,nkx)
        read(5,*)nky,(yl(i),i=1,nky+1),(yc(i),i=1,nky),(nyl(i),
     &  i=1,nky),(nyr(i),i=1,nky)
        read(5,*)(dymn(i),i=1,nky)
        read(5,*)
        write(13,mesh)
      
	  read(5,*)
	  read(5,*)nobstype
	  do n=1,max(1,nobstype)
          read(5,*)nobs(n)
c	    read(5,*)(oa2((n-1)*20+i),i=1,nobs(n)),
c     &		(oa1((n-1)*20+i),i=1,nobs(n)),(ob2((n-1)*20+i),
c     &		i=1,nobs(n)),(ob1((n-1)*20+i),i=1,nobs(n))
c          read(5,*)(oc2((n-1)*20+i),i=1,nobs(n)),
c     &		(oc1((n-1)*20+i),i=1,nobs(n)),(ioh((n-1)*20+i),
c     &		i=1,nobs(n))
c          read(5,*)(od2((n-1)*20+i),i=1,nobs(n)),
c     &		(od1((n-1)*20+i),i=1,nobs(n)),(oe2((n-1)*20+i),
c     &		i=1,nobs(n)),(oe1((n-1)*20+i),i=1,nobs(n))
c	    read(5,*)(nxo((n-1)*20+i),i=1,nobs(n)),
c     &		(mxo((n-1)*20+i),i=1,nobs(n)),(nyo((n-1)*20+i),
c     &		i=1,nobs(n)),(myo((n-1)*20+i),i=1,nobs(n))

        do nn=1,nobs(n)
         read(5,*)oa2((n-1)*20+nn),
     &  oa1((n-1)*20+nn),ob2((n-1)*20+nn),
     &  ob1((n-1)*20+nn)
         read(5,*)oc2((n-1)*20+nn),
     &  oc1((n-1)*20+nn),ioh((n-1)*20+nn)
c remove higher power terms, fyshi
c            read(5,*)od2((n-1)*20+nn),
c     &  od1((n-1)*20+nn),oe2((n-1)*20+nn),
c     &  oe1((n-1)*20+nn)
c           read(5,*)nxo((n-1)*20+nn),
c     &  mxo((n-1)*20+nn),nyo((n-1)*20+nn),
c     &  myo((n-1)*20+nn)
c assign zero for those terms, fyshi
         od2((n-1)*20+nn)=0.
         od1((n-1)*20+nn)=0.
         oe2((n-1)*20+nn)=0.
         oe1((n-1)*20+nn)=0.
         nxo((n-1)*20+nn)=0.
         mxo((n-1)*20+nn)=0.
         nyo((n-1)*20+nn)=0.
         myo((n-1)*20+nn)=0.
        enddo
	
        end do
        read(5,*)
        write(13,obstcl)

C.......read porous material information
	  read(5,*)
	  read(5,*)npor
	  if (npor.ne.0.and.npor.ne.1.and.npor.ne.10) then
		write(*,*)'you must specify npor=0 (no), 1 (porous) or 10'
		stop
	  endif
	  write(9,*)'porous switch is', npor
        xporosity(1)=1.0d0
        gc(1)=1.0d0
        xa(1)=0.0d0
        xxb(1)=0.0d0
        d50(1)=1.0d0
C.......define porous struture (npor=1)
        if (npor.eq.1) then
          read(5,*)nportype
          do 700 n=2,nportype+1
            read(5,*)d50(n),xporosity(n),xalpha(n),xbeta(n),xgamma(n)
c make alpha beta and gamma readable, fyshi    
c                xalpha(n)=300.0
c                xbeta(n)=1.1
c                xgamma(n)=0.34
                xxnu=max(1.0d-6,xnu)
                xa(n)=xalpha(n)*(1.0-xporosity(n))**2
     &                  /xporosity(n)**3*xxnu/abs(gy)/d50(n)**2
                xxb(n)=xbeta(n)*(1.0-xporosity(n))/xporosity(n)
     &                  **3/abs(gy)/d50(n)
                gc(n)=(1.0+xgamma(n)*(1.0-xporosity(n))
     &                  /xporosity(n))/xporosity(n)
        		write(9,*)'porous media',n,'a=',xa(n)*abs(gy)
     &		/gc(n),'b=',xxb(n)*(1.0+7.5/(0.01*1.0
     &		/xporosity(n)/d50(n)))*abs(gy)/gc(n)
            read(5,*)nporous(n-1)
           do i=1,nporous(n-1)
            read(5,*)pa2((n-2)*20+i),
     &	     pa1((n-2)*20+i),
     &	     pb2((n-2)*20+i),
     &	     pb1((n-2)*20+i)
            read(5,*)pc2((n-2)*20+i),
     &	    pc1((n-2)*20+i),
     &	    ipr((n-2)*20+i)
           enddo
c	      read(5,*)(pd2((n-2)*20+i),i=1,nporous(n-1)),
c     &	    (pd1((n-2)*20+i),i=1,nporous(n-1)),
c     &	    (pe2((n-2)*20+i),i=1,nporous(n-1)),
c     &	    (pe1((n-1)*20+i),i=1,nporous(n-1)) 
c            read(5,*)(nxp((n-2)*20+i),i=1,nporous(n-1)),
c     &	    (mxp((n-2)*20+i),i=1,nporous(n-1)),
c     &	    (nyp((n-2)*20+i),i=1,nporous(n-1)),
c     &	    (myp((n-2)*20+i),i=1,nporous(n-1))
            write(13,porous)
700       continue
        endif
C.......define simple porous layer (npor=10; allow partial cell)
        if (npor.eq.10) then
	    read(5,*)d50(2),xporosity(2)
          xalpha(2)=300.0
          xbeta(2)=1.1
          xgamma(2)=0.34
          read(5,*)porstart,porheight,nporsection,
     &		(porlength(n),porslope(n),n=1,nporsection)
	  endif
	  read(5,*)
	  write(9,*)'porous media finished'
	  write(9,*)

	  read(5,*)
        read(5,*)nfrsrf
        read(5,*)(fa2(i),i=1,nfrsrf),(fa1(i),i=1,nfrsrf)
	print*,(fa2(i),i=1,nfrsrf),(fa1(i),i=1,nfrsrf)
        read(5,*)(fb2(i),i=1,nfrsrf),(fb1(i),i=1,nfrsrf)
        read(5,*)(fc2(i),i=1,nfrsrf),(fc1(i),i=1,nfrsrf),(ifh(i),
     &	i=1,nfrsrf)
        read(5,*)(fd2(i),i=1,nfrsrf),(fd1(i),i=1,nfrsrf)
        read(5,*)(fe2(i),i=1,nfrsrf),(fe1(i),i=1,nfrsrf)
        read(5,*)(nxf(i),i=1,nfrsrf),(mxf(i),i=1,nfrsrf)
        read(5,*)(nyf(i),i=1,nfrsrf),(myf(i),i=1,nfrsrf)
        read(5,*)flht
        read(5,*)
        write(13,freesurf)

C.......read wave parameter
	  read(5,*)
	  read(5,*)aa,h0,h0r,ninflow
      if (kl.ne.6.and.ninflow.ne.7.and.ninflow.ne.8.and.ninflow.ne.70.
     &	and.ninflow.ne.80.and.ninflow.ne.71.and.ninflow.ne.81.and.
     &	ninflow.ne.100.and.ninflow.ne.0.and.ninflow.ne.10) then
	write(*,*)'inconsistence between kl and ninflow; respecify'
	stop
	  endif
        write(9,*)'wave height is',aa,'inflow water depth is',h0,
     &	'outflow water depth is',h0r,'type of inflow is',ninflow
        if (ninflow.eq.100) then
	read(5,*)isources,isourcee,jsources,jsourcee,nsource,tsource
		xxt=tsource
	  endif

C.......if pressure-driven wavemaker is used (currently only linear wave
C.......with ninflow=54), the matrix will become asymmetric; enforce
C.......sym=.false. to ensure correct pressure computation
        if (ninflow.eq.54) sym=.false.

        if (ninflow.eq.100) then
C.........irregular wave
          if (nsource.eq.44) then
	      read(5,*)nwave
		  do nw=1,nwave
		    read(5,*)aawave(nw),twave(nw)
		  end do
	      do nw=1,nwave
C...........calculate the wave number based on dispersion relation using
C...........Newton Ralphson method
              segmasqr=(2.0d0*pi/twave(nw))**2
              ctmp=sqrt(-gy*h0)
              xltmp=ctmp*twave(nw)
              xxk=2.0d0*pi/xltmp
              n=0
355           fk=-gy*xxk*tanh(xxk*h0)-segmasqr
              if (abs(fk).le.1.0e-6.or.n.gt.1000) goto 365
              fkdif=-gy*xxk*h0*(1.0d0-tanh(xxk*h0)**2)-gy*tanh(xxk*h0)
              xxk=xxk-fk/fkdif
              n=n+1
              goto 355
365           write(9,*)'nsource=',nsource,'xxk=',xxk
              xxl=2.0d0*pi/xxk
              cwave(nw)=xxl/twave(nw)
	      end do
	    endif

C.........fifth order Stokes wave 
	    if (nsource.eq.14) then
	      call stokes(aa,h0,tsource,xxk,amp1,amp2,amp3,
     &		amp4,amp5,cnf,gy,pi)
	      write(9,*)'nsource=',nsource,'xxk=',xxk,'cnf=',cnf
	      write(9,*)'a1=',amp1,'a2=',amp2,'a3=',amp3
	      write(9,*)'a4=',amp4,'a5=',amp5		
            xxl=2.0d0*pi/xxk
            segma=2.*pi/tsource
            c1=xxl/tsource
	    endif

          if (nsource.eq.34.or.nsource.eq.4) then
C...........calculate the wave number based on dispersion relation using
C...........Newton Ralphson method
            segmasqr=(2.0d0*pi/tsource)**2
            ctmp=sqrt(-gy*h0)
            xltmp=ctmp*tsource
            xxk=2.0d0*pi/xltmp
            n=0
155         fk=-gy*xxk*tanh(xxk*h0)-segmasqr
            if (abs(fk).le.1.0e-6.or.n.gt.1000) goto 165
            fkdif=-gy*xxk*h0*(1.0d0-tanh(xxk*h0)**2)-gy*tanh(xxk*h0)
            xxk=xxk-fk/fkdif
            n=n+1
            goto 155
165         write(9,*)'nsource=',nsource,'n=',n,'xxk=',xxk,'f(k)=',fk
            xxl=2.0d0*pi/xxk
            segma=2.*pi/tsource 
            c1=xxl/tsource
C...........if it's stokes wave, find zero for second-order waves;
	      if (nsource.eq.4) then
              atmp=aa/2
              btmp=aa**2*xxk/16.*cosh(xxk*h0)/sinh(xxk*h0)**3
     &          *(2.+cosh(2.*xxk*h0))
              xtmp=(-atmp+sqrt(atmp**2+8.0d0*btmp**2))/4.0d0/btmp
              cnf=asin(xtmp)
	      write(9,*)'a=',atmp,'b=',btmp,'cnf=',cnf
	      if (0.3*atmp.lt.btmp) then
		  write(9,*)'Long wave; Stokes wave is not applicable'
		  stop
	      endif
	      endif
          endif
	  
	    if (nsource.eq.5) then
	      read(5,*)fxstart
C...........reset fxstart=6.0
C	      fxstart=6.0
            c1=sqrt(-gy*h0*(1.+aa/h0))
            xstart=fxstart*h0/sqrt(aa/h0)
            write(9,*)'nsource=',nsource,'c1=',c1
	    endif
          
		if (nsource.eq.24) then 
            call cnoidal(aa,h0,tsource,xxl,xxc,ytrough,mod1)
            if (xxl.gt.1.0e6) then
              write(9,*)'no cnoidal wave exists for given parameters.'
              write(9,*)'xxl=',xxl,'yt=',ytrough,'xxc=',xxc,'mod1=',mod1
              write(iotty,*)
     &          'no cnoidal wave exists for given parameters.'
              stop
            endif
            write(9,*)'xxl=',xxl,'yt=',ytrough,'xxc=',xxc,'mod1=',mod1
            c1=xxc
            xxk=2*pi/xxl
C...........iterative method to find zero of cnoidal wave
                zup=1.0
                zlow=0.0
                zmid=(zup+zlow)/2.0d0
                nzero=0
3000            nzero=nzero+1
                zero1=ytrough+aa*cn(zmid*0.5*ck(mod1),mod1)**2
                write(9,*)nzero,zero1,zmid
                if (abs(zero1).le.1.0e-6) goto 3100
                if (nzero.gt.1000) then
                        write(9,*)'too many iterations; stop'
                        stop
                endif
                if (zero1.lt.0.0) then
                        zup=zmid
                        zmid=(zup+zlow)/2.0d0
                        goto 3000
                else
                        zlow=zmid
                        zmid=(zup+zlow)/2.0d0
                        goto 3000
                endif
3100            continue
                cnf=zmid
                write(9,*)'cnf=',cnf
          endif
	  endif

C.......for periodic wave (second-order Stokes or linear) 
        if (ninflow.eq.4.
     &    or.ninflow.eq.34.or.ninflow.eq.54) then
          read(5,*)xxt,areturn
C.........calculate the wave number based on dispersion relation using
C.........Newton Ralphson method
          segmasqr=(2.0d0*pi/xxt)**2
          ctmp=sqrt(-gy*h0)
          xltmp=ctmp*xxt
          xxk=2.0d0*pi/xltmp
          n=0
55        fk=-gy*xxk*tanh(xxk*h0)-segmasqr
          if (abs(fk).le.1.0e-6.or.n.gt.1000) goto 65
          fkdif=-gy*xxk*h0*(1.0d0-tanh(xxk*h0)**2)-gy*tanh(xxk*h0)
          xxk=xxk-fk/fkdif
          n=n+1
          goto 55
65        write(9,*)'ninflow=',ninflow,'n=',n,'xxk=',xxk,'f(k)=',fk
          xxl=2.0d0*pi/xxk
          segma=2.*pi/xxt
	    c1=xxl/xxt

C.........if it's stokes wave, find zero for second-order waves;
          if (ninflow.eq.4) then
              atmp=aa/2
              btmp=aa**2*xxk/16.*cosh(xxk*h0)/sinh(xxk*h0)**3
     &          *(2.+cosh(2.*xxk*h0))
              xtmp=(-atmp+sqrt(atmp**2+8.0d0*btmp**2))/4.0d0/btmp
              cnf=asin(xtmp)
              write(9,*)'a=',atmp,'b=',btmp,'cnf=',cnf
              if (atmp/3.0d0.lt.btmp) then
                write(9,*)'Large wave; Stokes not valid; try cnoidal'
                stop
	        endif
C.............calculate the mass transport due to periodic wave
			ulinear=abs(gy)*aa**2/c1/h0/8.0
              usecond=abs(gy)*
     &          (2*aa**2*xxk/16.*cosh(xxk*h0)/sinh(xxk*h0)**3
     &          *(2.+cosh(2.*xxk*h0)))**2/c1/h0/8.0
              write(9,*)'ulinear=',ulinear,'usecond=',usecond
          endif

C.........if it's linear wave, only find the mass transport
          if (ninflow.eq.34) then
              ulinear=abs(gy)*aa**2/c1/h0/8.0
              write(9,*)'ulinear=',ulinear,'usecond=',usecond
          endif
        endif

C.......parameters for solitary waves
        if (ninflow.eq.5) then
C	    read(5,*)fxstart
C.........reset fxstart=6.0
          fxstart=4.0
          alpha1=aa/h0
          c1=sqrt(-gy*h0*(1.+aa/h0))
          xstart=fxstart*h0/sqrt(aa/h0)
        endif

C.......for cnoidal wave only
        if (ninflow.eq.24) then
          read(5,*)xxt,areturn
          call cnoidal(aa,h0,xxt,xxl,xxc,ytrough,mod1)
	    if (xxl.gt.1.0e6) then
            write(9,*)'no cnoidal wave exists for given parameters.'
            write(9,*)'xxl=',xxl,'yt=',ytrough,'xxc=',xxc,'mod1=',mod1
            write(iotty,*)
     &          'no cnoidal wave exists for given parameters.'
            write(iotty,*)
     &          'xxl=',xxl,'yt=',ytrough,'xxc=',xxc,'mod1=',mod1
	      stop
	    endif	
          write(9,*)'xxl=',xxl,'yt=',ytrough,'xxc=',xxc,'mod1=',mod1
          c1=xxc
          xxk=2*pi/xxl
C.........calculate the mass transport due to periodic wave
          if (ninflow.eq.24) then
                ulinear=abs(gy)*aa**2/c1/h0/8.0
                write(9,*)'ulinearc=',ulinear
C...............iterative method to find zero of cnoidal wave
                zup=1.0
                zlow=0.0
                zmid=(zup+zlow)/2.0d0
                nzero=0
2000            nzero=nzero+1
                zero1=ytrough+aa*cn(zmid*0.5*ck(mod1),mod1)**2
                write(9,*)nzero,zero1,zmid
                if (abs(zero1).le.1.0e-6) goto 2100
                if (nzero.gt.1000) then
                        write(9,*)'too many iterations; stop'
                        stop
                endif
                if (zero1.lt.0.0) then
                        zup=zmid
                        zmid=(zup+zlow)/2.0d0
                        goto 2000
                else
                        zlow=zmid
                        zmid=(zup+zlow)/2.0d0
                        goto 2000
                endif
2100            continue
                cnf=zmid
                write(9,*)'cnf=',cnf
          endif
        endif

C.......for plunging jet or plane jet only
	  if ((ninflow.eq.7.or.ninflow.eq.70.or.ninflow.eq.71).or.
     &	(ninflow.eq.8.or.ninflow.eq.80.or.ninflow.eq.81)) then
C.........ninflow=7 is the horizontal jet and ninflow=8 is the vertical jet
	    read(5,*)ujet,vjet,wjet,yjet,time_jet
	    yup=yjet+wjet/2
	    ylow=yjet-wjet/2
       if (ninflow.eq.71.or.ninflow.eq.81) ustar=sqrt(ujet**2+vjet**2)
     &	  /(5.75*log10(wjet/2/0.0001)+6.0)
	  endif

C ... for plunging jet and take out fyshi
          if (ninflow.eq.10)then
            read(5,*)ujet,vjet,wjet,yjet,time_jet
            read(5,*)uout,vout,wout,yout
            yup=yjet+wjet/2.
            ylow=yjet-wjet/2.
            yupout=yout+wout/2.
            ylowout=yout-wout/2.
          endif

C.......To specif boundary condition from another numerical results
        if (ninflow.eq.9) then
          open(51,file='kuang_bim',status="unknown")
	    open(52,file='kuang_uvw',status="unknown")
	    open(53,file='kuang_det',status="unknown")
	    do 21 j=1,50
		  read(53,*)det_k(j)
21	    continue
	    j_det=1
	    prtdt=det_k(j_det)
          it=1
	    itotal=536
	    iinter=399
	    tpd=1.0
          do 69 j=1,itotal
            read(51,*)t_pad(j),x_pad(j),h_pad(j)
            h_pad(j)=h_pad(j)/1000.0
            x_pad(j)=x_pad(j)/1000.0
            a_pad=0.09
            omega=2*pi*0.6
            u_pad(j)=(x_pad(j)-x_pad(j-1))/
     &          max(1.0d-6,t_pad(j)-t_pad(j-1))
69        continue
          do 79 j=itotal-iinter,itotal
            t_pad(j+iinter+1)=t_pad(j)+tpd
            x_pad(j+iinter+1)=x_pad(j)
            h_pad(j+iinter+1)=h_pad(j)
79        continue
          do 791 j=itotal+1,itotal+1+iinter
            t_pad(j+iinter+1)=t_pad(j)+tpd
            x_pad(j+iinter+1)=x_pad(j)
            h_pad(j+iinter+1)=h_pad(j)
791       continue
          do 792 j=itotal+iinter+2,itotal+2+iinter+iinter
            t_pad(j+iinter+1)=t_pad(j)+tpd
            x_pad(j+iinter+1)=x_pad(j)
            h_pad(j+iinter+1)=h_pad(j)
792       continue
	    j=1
	    k=0
689       read(52,'(a2)')AAA
          if (AAA.ne.' T') then
                backspace(52)
		  k=k+1
                read(52,*)t1,ybc(j,k),ubc(j,k),vbc(j,k)
          else
		  j=j+1
		  k=0
          endif
	    if (j.le.itotal) goto 689
	    do 679 k=1,35
	      ybc(1,k)=ybc(2,k)
	      ubc(1,k)=0.0
	      vbc(1,k)=0.0
679	    continue
          do 687 k=1,35
            do 681 j=itotal-iinter,itotal
                ybc(j+iinter+1,k)=ybc(j,k)
                ubc(j+iinter+1,k)=ubc(j,k)
                vbc(j+iinter+1,k)=vbc(j,k)
681         continue
            do 682 j=itotal+1,itotal+1+iinter
                ybc(j+iinter+1,k)=ybc(j,k)
                ubc(j+iinter+1,k)=ubc(j,k)
                vbc(j+iinter+1,k)=vbc(j,k)
682         continue
            do 683 j=itotal+2+iinter,itotal+2+iinter+iinter
                ybc(j+iinter+1,k)=ybc(j,k)
                ubc(j+iinter+1,k)=ubc(j,k)
                vbc(j+iinter+1,k)=vbc(j,k)
683         continue
687       continue
	  endif

        xxkr=0.0
	  if (h0r.eq.h0) then
		xxkr=xxk
		c1r=c1
		goto 275
	  endif
	
	  if (nsource.eq.44) goto 275
        if ((nopen.eq.1.or.nopen.eq.11).and.kr.eq.3) then
          if (ninflow.eq.100) then
		if (nsource.eq.34.or.nsource.eq.4.or.nsource.eq.14) 
     &		then
	                ttmp=tsource
		else
			if (nsource.eq.5) then
		    	  c1r=sqrt(-gy*(h0r+aa))
			  goto 275
			endif
			if (nsource.eq.24) then
			  c1r=sqrt(-gy*(h0r+aa+ytrough))
			  goto 275
			endif
		endif
          else
              if (ninflow.eq.4.
     &          or.ninflow.eq.34) then
                ttmp=xxt
		    else
			  if (ninflow.eq.24) then
				c1r=sqrt(-gy*(h0r+aa+ytrough))
				goto 275
			  else
	                if (ninflow.eq.5) then
	                          c1r=sqrt(-gy*(h0r+aa))
	                          goto 275
					endif
					c1r=sqrt(-gy*h0r)
					goto 275
			  endif	
              endif
          endif
          segmasqr=(2.0d0*pi/ttmp)**2
          ctmp=sqrt(-gy*h0r)
          xltmp=ctmp*ttmp
          xxkr=2.0d0*pi/xltmp
          n=0
255       fk=-gy*xxkr*tanh(xxkr*h0r)-segmasqr
          if (abs(fk).le.1.0e-6.or.n.gt.1000) goto 265
          fkdif=-gy*xxkr*h0r*(1.0d0-tanh(xxkr*h0r)**2)-
     &		gy*tanh(xxkr*h0r)
          xxkr=xxkr-fk/fkdif
          n=n+1
          goto 255
265       xxlr=2.0d0*pi/xxkr
          segma=2.*pi/ttmp
          c1r=xxlr/ttmp
	  endif

275	  write(9,*)'c1=',c1,'c1r=',c1r
	  write(9,*)'inflow and outflow boundary condition finished'
	  
	  read(5,*)
	  write(9,*)'wave parameter finished'
	  write(9,*)

C.......read output format
	  read(5,*)
	  read(5,*)lout,nanimation,nmean
        write(9,*)'gauge switch is',lout,'  animation switch is',
     &	nanimation,'  mean flow switch is',nmean
        twprt=tstart
        twplt=0.0
        write(9,*)'output format finished'
	  write(9,*)
C.......output location
        if (lout.eq.1) then
          read(5,*)nloc,(xout(i),i=1,nloc),xxxf,ttend,prtdt_t
          do 88 n=1,nloc
            open(700+5*(n-1)+1,file=fdir//
     &      'gage'//achac(n),status="unknown")
88        continue
        endif
        if (prtdt.lt.0.0) then
                read(5,*)ndt,(det_k(i),i=1,ndt)
                prtdt=det_k(1)
                j_det=1
        endif
C.......animation data set
        if (nanimation.ne.0) read(5,*)tstart_a,tfinish_a,prtdt_a
        twprt_a=tstart_a
C.......controling parameter for mean quantities output
        if (nmean.eq.1) then
          read(5,*)istart,iend,iinterval,tmstart,tmend,tinterval
          open(70,file='eflux',status="unknown")
          open(72,file='energy',status="unknown")
          open(73,file='eta',status="unknown")
          open(74,file='efluxt',status="unknown")
        endif

	  if (nmean.eq.2) then
		xmnt=0.0
		call setarry (xmnf(1),0.0d0,imax*jmax)
		call setarry (xmnu(1),0.0d0,imax*jmax)
		call setarry (xmnv(1),0.0d0,imax*jmax)
		call setarry (xmnk(1),0.0d0,imax*jmax)
	  endif

	  read(5,*)

C.......read other parameters
	  read(5,*)
	  read(5,*)kemodel,roughness,nopen,npollutant,nrs,novertop,nfree,
     &	islip

	  if (nopen.eq.11) then
		xsponge=1.5*xxl
	    adamp=200.0
	    power=10.0
		write(*,*)'sponge layer parameter',xsponge,adamp,power
	  endif

C.......Check validaty range of kemodel
	  if (kemodel.ne.0.and.kemodel.ne.1.and.kemodel.ne.4) then
		write(*,*) 'you must specify kemodel=0 (laminar), 1 (linear), 
     &	or 4 (nonlinear)'
		stop
	  endif

C.....For turbulent flow, based on roughness to justify whether the surface is smooth or rough
	if (kemodel.gt.0) then
	  if (roughness.ge.0.0000100000001) then
C.........rough thus XE<0
		XE=-1.0
	  else
C.........smooth and define XE=9.0
		XE=9.0
	  endif
	endif

	if (kemodel.gt.0) then
C          open(27,file='k',status="unknown")
C          open(28,file='eps',status="unknown")
C          open(30,file='eddy',status="unknown")
C		 open(71,file='disp',status="unknown")
	endif
      if (kemodel.gt.0.and.xnu.eq.0.0) then
            xnu=1.0e-6
            write(13,*)'warning: xnu=0.0 when turbulence model on.'
            write(iotty,*)'warning: xnu=0.0 when turbulence model on.'
      endif

	if (npollutant.ne.0) then
           open(33,file='pol',status="unknown")
           read(5,*)(xpol(n),ypol(n),conc(n),n=1,npollutant)
	endif

	if (novertop.eq.1) read(5,*)xovertop

	write(9,*)'turbulence switch is',kemodel,'  roughness is',
     &	roughness,'coefficient for smooth surface is',XE,
     &	'open boundary type is',nopen,'  pollutant switch is',
     &	npollutant,'  restart switch is',nrs,'  overtopping
     &	switch is',novertop,'  with xovertop=',xovertop
	write(9,*)'other parameters finished'
	write(9,*)
	read(5,*)

	read(5,*)
	read(5,*)nmb,acmin
	if (nmb.eq.2.or.nmb.eq.3) read(5,*)uxmb(1),vymb(1)
	if (nmb.eq.4.or.nmb.eq.5) read(5,*)axmb,aymb,tpmb

C.....parameters for solitary waves
      if (nmb.eq.25) then
C	    read(5,*)fxstart
		fxstart=4.0d0
          c1=sqrt(-gy*h0*(1.+aa/h0))
          xstart=fxstart*h0/sqrt(aa/h0)
      endif

C.....for cnoidal wave only
	if (nmb.eq.24) then
          read(5,*)xxt
          call cnoidal(aa,h0,xxt,xxl,xxc,ytrough,mod1)
	    if (xxl.gt.1.0e6) then
            write(9,*)'no cnoidal wave exists for given parameters.'
            write(9,*)'xxl=',xxl,'yt=',ytrough,'xxc=',xxc,'mod1=',mod1
            write(iotty,*)
     &          'no cnoidal wave exists for given parameters.'
            write(iotty,*)
     &          'xxl=',xxl,'yt=',ytrough,'xxc=',xxc,'mod1=',mod1
	      stop
	    endif	
          write(9,*)'xxl=',xxl,'yt=',ytrough,'xxc=',xxc,'mod1=',mod1
          c1=xxc
          xxk=2*pi/xxl
C.........iterative method to find zero of cnoidal wave
          zup=1.0
          zlow=0.0
          zmid=(zup+zlow)/2.0d0
          nzero=0
2001      nzero=nzero+1
          zero1=ytrough+aa*cn(zmid*0.5*ck(mod1),mod1)**2
          write(9,*)nzero,zero1,zmid
          if (abs(zero1).le.1.0e-6) goto 2101
          if (nzero.gt.1000) then
                        write(9,*)'too many iterations; stop'
                        stop
          endif
          if (zero1.lt.0.0) then
                        zup=zmid
                        zmid=(zup+zlow)/2.0d0
                        goto 2001
          else
                        zlow=zmid
                        zmid=(zup+zlow)/2.0d0
                        goto 2001
         endif
2101     continue
         cnf=zmid
         write(9,*)'cnf=',cnf
      endif

	write(9,*)'nmb=',nmb,'acmin=',acmin,'uxmb=',
     &		uxmb(1),'vymb=',vymb(1)
	write(9,*)'moving body parameters finished'
	read(5,*)

c
C.......if ninflow=100 (source function used), force idiv=0
	if (ninflow.eq.100) idiv=0

C.......set flht=h0 if there is no obstacle
	if (nobstype.eq.0) flht=h0

c
c.... Get problem name and print it out
      iotty=6
      jnm=" RIPPLE "
      write (iotty,500) prbname
c
c.... Compute constant terms and initialize necessary variables
c
      pi=acos(-1.0d0)
      if (icyl.eq.0) tpi=1.0d0
      cyl=float(icyl)
      omcyl=1.-cyl
      emf1=1.d00-emf
c.... flux limiter
      if (con*1.3.le.fcvlim) go to 10
      write (13,200) con,fcvlim
      write (iotty,200) con,fcvlim
      fcvlim=1.3*con
   10 continue
c.... minimum time step
      dtend=0.0001*delt
      if (dmpdt.eq.0.0) dmpdt=twfin
      twdmp=dmpdt
      frsurf=frctn*rhof
      omalp=1.-alpha
      opalp=1.+alpha
      if (kl.eq.5) pbcl=1.0d0
      if (kr.eq.5) pbcr=1.0d0
      if (kt.eq.5) pbct=1.0d0
      if (kb.eq.5) pbcb=1.0d0
      xmu=xnu*rhof
      if (alpha.gt.1.0d0) alpha=1.0d0
c
9999  return
   25 format(32x," RIPPLE ")
   50 format(15x,"Date: ",a8,2x,"Time: ",a8,2x,"Machine: ",a8)
  100 format (/,5x,a80,/)
  150 format(/,15x,"* * * * * * * * * * * * * * * * * * * *",
     &       /,15x,"           ERROR DIAGNOSTICS           ",
     &       /,15x,"* * * * * * * * * * * * * * * * * * * *",/)
  175 format(/,15x,"* * * * * * * * * * * * * * * * * * * *",
     &       /,15x,"               RESTART                 ",
     &       /,15x,"* * * * * * * * * * * * * * * * * * * *",/)
  200 format (1x,4hcon=,1pe10.3,1x,11hand fcvlim=,e10.3,1x,
     &         41hare incompatible. setting fcvlim=1.3*con.)
  400 format(a80)
  500 format (" RIPPLE:  ",a80,/," Processing input data . . .")
      end
  
 
