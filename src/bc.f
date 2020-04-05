      subroutine bc
c
c ======================================================================
c
c   Purpose -
c     set boundary conditions
c
c   BC is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE     ACCEL   CONVECT  CONVECTC     SETUP    SRFFRCE
c	   VTILDE     VOFADV  MAC
c
c
c   BC calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c            none
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
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
	external ck,ce,cn

C.....automatically set kl=3 when solitary wave ends
	if (ninflow.eq.5.and.kl.eq.6) then
	  if (t.ge.2.0*xstart/c1) kl=3
	endif

      if (iter.eq.0) then
          eta0=0.0
          eta2=0.0
          eta1=0.0
          eta1a=0.0  
          eta0old=0.0
          eta2old=0.0  
          eta1old=0.0 
          eta2aold=0.0
      endif

	if (kl.ne.6) goto 600

	if (ninflow.eq.5) then 
C.......Specify Special Boundary Condition: solitary at left
	  if (iter.eq.0) then
		tex=sqrt(0.75*aa/h0**3)
                eta0=aa/cosh(tex*(xstart-delx(1)/2.-c1*t))**2
		eta2=aa/cosh(tex*(xstart+delx(2)/2.-c1*t))**2
		eta1=aa/cosh(tex*(xstart-c1*t))**2
                eta1a=aa/cosh(tex*(xstart-c1*(t-delt/2.0d0)))**2
		eta0x=eta0
		eta2x=eta2
C...............testing the weakly reflected boundary condition 
		if (nweakref.eq.1.and.ncyc.ge.1) then
                  eta0old=aa/
     &		    cosh(tex*(xstart-delx(1)/2.-c1*(t-delt)))**2
                  eta2old=aa/
     &		    cosh(tex*(xstart+delx(2)/2.-c1*(t-delt)))**2
                  eta1old=aa/cosh(tex*(xstart-c1*(t-delt)))**2
                  eta2aold=aa/
     &              cosh(tex*(xstart+delx(2)-c1*(t-delt)))**2
		endif  
	  endif
	  goto 600
	endif

	if (ninflow.eq.4) then
C.......Specify Special Boundary Condition: periodic wave from left
C.......Stokes second order
        if (iter.eq.0) then
                eta0=aa/2.*cos(pi/2.-xxk*(delx(1)/2.)-segma*t-cnf)
     &		+aa**2*xxk/16.*cosh(xxk*h0)/sinh(xxk*h0)**3
     &		*(2.+cosh(2.*xxk*h0))
     &		*cos(2.*(pi/2.-xxk*delx(1)/2-segma*t-cnf))
                eta2=aa/2.*cos(pi/2.+xxk*(delx(2)/2.)-segma*t-cnf)
     &          +aa**2*xxk/16*cosh(xxk*h0)/sinh(xxk*h0)**3
     &          *(2+cosh(2*xxk*h0))
     &		*cos(2*(pi/2.+xxk*delx(2)/2-segma*t-cnf))
                eta1=aa/2.*cos(pi/2.-segma*t-cnf)
     &          +aa**2*xxk/16*cosh(xxk*h0)/sinh(xxk*h0)**3
     &          *(2+cosh(2*xxk*h0))*cos(2*(pi/2.-segma*t-cnf))
                eta1a=aa/2.*cos(pi/2.-segma*(t-delt/2.0d0)-cnf)
     &          +aa**2*xxk/16*cosh(xxk*h0)/sinh(xxk*h0)**3
     &          *(2+cosh(2*xxk*h0))*
     &		cos(2*(pi/2.-segma*(t-delt/2.0)-cnf))
                eta0x=eta0
			eta2x=eta2
C...............testing the weakly reflected boundary condition
                if (nweakref.eq.1.and.ncyc.ge.1) then
                  eta0old=aa/2.0
     &		  *cos(pi/2.0-xxk*delx(1)/2.0-segma*(t-delt))
     &            +aa**2*xxk/16.*cosh(xxk*h0)/sinh(xxk*h0)**3
     &            *(2.+cosh(2.*xxk*h0))
     &		  *cos(2.*(pi/2.-segma*(t-delt)-xxk*delx(1)/2.))
                  eta2old=aa/2.
     &		  *cos(pi/2.0+xxk*delx(2)/2.0-segma*(t-delt))
     &            +aa**2*xxk/16*cosh(xxk*h0)/sinh(xxk*h0)**3
     &            *(2.+cosh(2*xxk*h0))
     &		  *cos(2.*(-segma*(t-delt)+pi/2.+xxk*delx(2)/2.))
                  eta1old=aa/2.
     &		  *cos(pi/2.-segma*(t-delt))
     &            +aa**2*xxk/16*cosh(xxk*h0)/sinh(xxk*h0)**3
     &            *(2.+cosh(2.*xxk*h0))
     &		  *cos(2.*(-segma*(t-delt)+pi/2.))
                  eta2aold=aa/2.0
     &            *cos(pi/2.0+xxk*delx(2)-segma*(t-delt))
     &            +aa**2*xxk/16.0*cosh(xxk*h0)/sinh(xxk*h0)**3
     &            *(2.0+cosh(2.0*xxk*h0))
     &            *cos(2.0*(pi/2.0-segma*(t-delt)+xxk*delx(2)))
		endif
        else
        endif
	  goto 600
	endif

	if (ninflow.eq.9) then
C.....special b.c. from specified condition
237     if (t.gt.t_pad(it+1)) then
          it=it+1
          goto 237
        else
          if (t.lt.t_pad(it)) then
            eta1=h_pad(it)
          else
            eta1=(h_pad(it)*(t_pad(it+1)-t)+h_pad(it+1)*(t-t_pad(it)))
     &          /(t_pad(it+1)-t_pad(it))
          endif
          eta0=eta1
          eta2=eta1
          eta0x=eta0
	    eta2x=eta2
          uleft1=eta1*1.30/(eta1+h0)
          vleft1=0.0
        endif
        goto 700
	endif

	if (ninflow.eq.34) then
C.......Specify Special Boundary Condition: linear periodic wave from left
        if (iter.eq.0) then
                eta0=aa/2.*cos(pi/2.-xxk*(delx(1)/2.)-segma*t)
                eta2=aa/2.*cos(pi/2.+xxk*(delx(2)/2.)-segma*t)
                eta1=aa/2.*cos(pi/2.-segma*t)
                eta1a=aa/2.*cos(pi/2.-segma*(t-delt/2.0d0))
C...............consider to compensate reflected wave (for vertical wall)
                if (nweakref.eq.10) then
C.................constant water depth
C                 distotal=(x(im1)-x(1))*2.0d0
C                 refcoef=1.0d0
C                 tref=distotal/c1
C.................wave on slope; specify the tref and refcorf manually
C.................from experiment. In Sakakiyama's case, assuming the
C.................reflected wave has a=1.3cm from exp
                  tref=10.09
                  refcoef=1.3/10.5
                  if (t.ge.tref) then
				  eta0=eta0+refcoef*
     &		      aa/2.*cos(pi/2.-xxk*(delx(1)/2.)-segma*(t-tref))
                    eta2=eta2+refcoef*
     &              aa/2.*cos(pi/2.+xxk*(delx(1)/2.)-segma*(t-tref))
                    eta1=eta1+refcoef*
     &              aa/2.*cos(pi/2.-segma*(t-tref))
                  endif
                endif

                eta0x=eta0
                eta2x=eta2
C...............testing the weakly reflected boundary condition
                if (nweakref.eq.1.and.ncyc.ge.1) then
                  eta0old=aa/2.0
     &            *cos(pi/2.0-xxk*delx(1)/2.0-segma*(t-delt))
                  eta2old=aa/2.
     &            *cos(pi/2.0+xxk*delx(2)/2.0-segma*(t-delt))
                  eta1old=aa/2.
     &            *cos(pi/2.-segma*(t-delt))
                  eta2aold=aa/2.
     &		  *cos(pi/2.+xxk*delx(2)-segma*(t-delt))
                endif
        else
        endif
	  goto 600
	endif

        if (ninflow.eq.24) then
C.........Specify Special Boundary Condition: cnoidal wave from left
          if (iter.eq.0) then
          	eta0=ytrough+aa*cn(cnf*0.5*ck(mod1)+2.0*ck(mod1)*
     &          (-delx(1)/2/xxl-t/xxt),mod1)**2
          	eta2=ytrough+aa*cn(cnf*0.5*ck(mod1)+2.0*ck(mod1)*
     &          (+delx(2)/2/xxl-t/xxt),mod1)**2
          	eta1=ytrough+aa*cn(cnf*0.5*ck(mod1)+2.0*ck(mod1)*
     &          (-t/xxt),mod1)**2
                eta1a=ytrough+aa*cn(cnf*0.5*ck(mod1)+2.0*ck(mod1)*
     &          (-(t-delt/2.0d0)/xxt),mod1)**2
                eta0x=eta0
                eta2x=eta2
C...............testing the weakly reflected boundary condition
          	if (nweakref.eq.1.and.ncyc.ge.1) then
            	  eta0old=ytrough+aa*cn(cnf*0.5*ck(mod1)+2.0*ck(mod1)*
     &            (-delx(1)/2/xxl-(t-delt)/xxt),mod1)**2
            	  eta2old=ytrough+aa*cn(cnf*0.5*ck(mod1)+2.0*ck(mod1)*
     &            (+delx(2)/2/xxl-(t-delt)/xxt),mod1)**2
            	  eta1old=ytrough+aa*cn(cnf*0.5*ck(mod1)+2.0*ck(mod1)*
     &            (-(t-delt)/xxt),mod1)**2
            	  eta2aold=ytrough+aa*cn(cnf*0.5*ck(mod1)+2.0*ck(mod1)*
     &            (+delx(2)/xxl-(t-delt)/xxt),mod1)**2
	  		endif
          endif
	   endif

600	  continue

        if (((nopen.eq.1.or.nopen.eq.11.or.nopen.eq.2).and.kl.eq.3.and
     &	.ncyc.ge.1.and.iter.eq.0.and.ibcflg0.ne.0).or.
C.........testing the weakly reflected boundary condition
     &    (nweakref.eq.1.and.ncyc.ge.1.and.iter.eq.0)) then
	    if (nopen.eq.2) then
	      eta0x=flht-h0
		  goto 720
	    endif

		eta00old=0.0
	    do 555 j=2,jm1
	      ij00=(j-1)*imax+1
		  eta00old=eta00old+ac(ij00)*dely(j)*fn(ij00)
	      if (fn(ij00).gt.0.0.and.fn(ij00+imax).eq.0.0) then
C     			eta00old=y(j-1)+fn(ij00)*dely(j)-flht
			eta00old=eta00old-h0
			goto 565
	      endif
555	    continue
565	    continue

          eta02old=0.0
		do 557 j=2,jm1
            ij02=(j-1)*imax+2
		  eta02old=eta02old+ac(ij02)*dely(j)*fn(ij02)
            if (fn(ij02).gt.0.0.and.fn(ij02+imax).eq.0.0) then
C                eta02old=y(j-1)+fn(ij02)*dely(j)-flht
			eta02old=eta02old-h0
		goto 558
	    endif
557       continue
558	    continue

	    eta01old=(eta00old+eta02old)/2.0

C.........use dispersion relationship to calculate phase velocity
	    if (xxk.gt.0.0) then
C              c01=sqrt(-gy/xxk*tanh(xxk*(h0+eta01old)))
	      c01=sqrt(-gy/xxk*tanh(xxk*h0))
C.........It's found the smaller C tends to give better (larger) mass conservation
C.........for nonlinear waves; tentative treatment; further study needed
C		  c01=c01/(1.0+aa/h0)
          else
C.........use averaged wave to calculate phase speed
              c01=sqrt(-gy*h0)
C	      c01=sqrt(-gy*(h0+eta01old))
          endif

C.........use dispersion relationship to calculate phase velocity
          if (xxk.gt.0.0) then
C     	      c02=sqrt(-gy/xxk*tanh(xxk*(h0+eta02old)))
	      c02=sqrt(-gy/xxk*tanh(xxk*h0))
C.........It's found the smaller C tends to give better (larger) mass conservation
C.........for nonlinear waves; tentative treatment; further study needed
C		  c02=c02/(1.0+aa/h0)
          else
C.........use averaged wave to calculate phase speed
              c02=sqrt(-gy*h0)
C	      c02=sqrt(-gy*(h0+eta02old))
          endif

C.........find the actual free surface displacement at the current time step
	    eta00=c01*delt/((delx(1)+delx(2))/2.0)*((eta02old-eta2old)
     &		-(eta00old-eta0old))+(eta00old-eta0old)+eta0
C.........added a RHS from Petit et al (1994)
          if (petit.ne.0.0.and.max(xxt,tsource).gt.0.0)
     &          eta00=eta00-delt*2.0*pi/max(xxt,tsource)
     &          *petit*eta02old
C.........assign the value of eta00 to calculated free surface displacement
          eta0x=eta00
	  endif
720	  continue

        if ((nopen.eq.1.or.nopen.eq.11.or.nopen.eq.2).and.kr.eq.3.and
     &	.ncyc.ge.1.and.iter.eq.0.and.ibcflg0.ne.0) then
C        if (nopen.eq.1.and.kr.eq.3.and.ncyc.ge.1) then
C.......the first line works well for non-optimizing compiling; bugs
	    if (nopen.eq.2) then
		  etaimax=flht-h0r
		  goto 700
		endif

		etaim1old=0.0
          do 556 j=2,jm1
            ijim1=(j-1)*imax+im1
		  etaim1old=etaim1old+fn(ijim1)*dely(j)*ac(ijim1)
            if (fn(ijim1).gt.0.0.and.fn(ijim1+imax).eq.0.0) then
C                etaim1old=y(j-1)+fn(ijim1)*dely(j)-flht
			etaim1old=etaim1old-h0r
	        goto 656
	    endif
556	    continue
656	    continue

		etaimaxold=0.0
		do 756 j=2,jm1
            ijimax=(j-1)*imax+imax
		  etaimaxold=etaimaxold+fn(ijimax)*dely(j)*ac(ijimax)
            if (fn(ijimax).gt.0.0.and.fn(ijimax+imax).eq.0.0) then
C                etaimaxold=y(j-1)+fn(ijimax)*dely(j)-flht
			etaimaxold=etaimaxold-h0r
	        goto 856
	    endif
756       continue
856	  continue

          etarbcold=(etaim1old+etaimaxold)/2.0

C.........use dispersion relationship to calculate phase velocity
          if (xxkr.gt.0.0) then
C              cim1=sqrt(-gy/xxkr*tanh(xxkr*(h0r+etaim1old)))
	      cim1=sqrt(-gy/xxkr*tanh(xxkr*h0r))
C.........It's found the smaller C tends to give better (larger) mass conservation
C.........for nonlinear waves; tentative treatment; further study needed
C		  cim1=cim1/(1.0+aa/h0)
          else
C.........use shallow water relation to calculate phase speed
              cim1=sqrt(-gy*h0r)
C	      cim1=sqrt(-gy*(h0r+etaim1old))
          endif

C.........use dispersion relationship to calculate phase velocity
	  if (xxkr.gt.0.0.and.nsource.ne.44) then
C              crbc=sqrt(-gy/xxkr*tanh(xxkr*(h0r+etarbcold)))
	      crbc=sqrt(-gy/xxkr*tanh(xxkr*h0r))
C.........It's found the smaller C tends to give better (larger) mass conservation
C.........for nonlinear waves; tentative treatment; further study needed
C		  crbc=crbc/(1.0+aa/h0)
          else
C.........use shallow water relation to calculate phase speed
C              crbc=sqrt(-gy*(h0r+etarbcold))
	       crbc=sqrt(-gy*h0r)
          endif

          etaimax=-crbc*delt/((delx(im1)+delx(imax))/2.0)*
     &		(etaimaxold-etaim1old)+etaimaxold
	endif

700	continue

      ijl=1
      ijr=imax
      do 100 j=1,jmax
C.......modified to account for the impact
        if (nf(ijl+1).ne.2) then
		f(ijl)=f(ijl+1)
	    p(ijl)=p(ijl+1)
	    nf(ijl)=nf(ijl+1)	
	  else
		f(ijl)=0.0
          p(ijl)=0.0
          nf(ijl)=6
	  endif
	  if (kl.eq.1) nf(ijl)=0
        if (nf(ijr-1).ne.1) then
		f(ijr)=f(ijr-1)
	    p(ijr)=p(ijr-1)
		nf(ijr)=nf(ijr-1)
	  else
		f(ijr)=0.0
          p(ijr)=0.0
          nf(ijr)=6
	  endif
	  if (kr.eq.1) nf(ijr)=0

        go to (10,20,30,40,35,5,6), kl

   35	  continue
	  p(ijl)=pbc(3)+gy*yj(j)
	  go to 50

    6	  continue
C.......special boundary for grid turbulence
	  u(ijl)=1.0
	  xk(ijl)=0.05
	  xep(ijl)=0.045
	  v(ij)=0.0
	  goto 50
C.......end of modification

    5	  continue
C.......This is for any specified flows from left boundary
C.......We first start with free surface flow whether both u,v, and etc are given
C.......for periodic waves of second order Stokes
	  if (ninflow.eq.4) then
          u(ijl)=aa/2.*segma*cosh(xxk*yj(j))/sinh(xxk*h0)*
     &        cos(pi/2.-segma*t-cnf)
     &		+3./32.*2.*xxk*aa**2*segma*cosh(2.*xxk*yj(j))/
     &		sinh(xxk*h0)**4*cos(2.*(pi/2.-segma*t-cnf))
C...............substract the mean mass transport velocity to simulate the
C...............wave tank situation where no pure mass transport exists
C...............U=E/rho Ch (see Dean & Dalrymple, 1991) and E=1/8rho g H^2
     &          -(ulinear+usecond)*areturn
C.........minus the mass accumalation each time step
          if (nweakref.gt.1) u(ijl)=u(ijl)-xmass_d/delt/(h0+eta1)
          v(ijl)=aa/2.*segma*sinh(xxk*y(j))/sinh(xxk*h0)*
     &          sin(pi/2.-segma*t-cnf-xxk*(delx(1)/2.))
     &		+3./32.*2.*xxk*aa**2*segma*sinh(2.*xxk*y(j))
     &		/sinh(xxk*h0)**4
     &		*sin(pi/2-segma*t-cnf-xxk*(delx(1)/2.))
C.........testing the weakly reflected boundary condition
          if (nweakref.eq.1.and.ncyc.ge.1) then
            uleft1=u(ijl)
            uleft1old=aa/2.*segma*cosh(xxk*yj(j))/sinh(xxk*h0)*
     &          cos(pi/2.-segma*(t-delt)-cnf)
     &          +3./32.*2.*xxk*aa**2*segma*cosh(2.*xxk*yj(j))/
     &          sinh(xxk*h0)**4*cos(2.*(pi/2.-segma*(t-delt)-cnf))
     &          -(ulinear+usecond)*areturn
            uleft2old=aa/2.*segma*cosh(xxk*yj(j))/sinh(xxk*h0)*
     &          cos(pi/2.-segma*(t-delt)-cnf+xxk*delx(2))
     &          +3./32.*2.*xxk*aa**2*segma*cosh(2.*xxk*yj(j))/
     &          sinh(xxk*h0)**4
     &		*cos(2.*(-segma*(t-delt)-cnf+pi/2.+xxk*delx(2)))
     &          -(ulinear+usecond)*areturn
            u(ijl)=c02*delt/delx(2)*((un(ijl+1)-uleft2old)-
     &          (un(ijl)-uleft1old))+
     &          (un(ijl)-uleft1old)+uleft1
C...............added a RHS from Petit et al (1994)
     &          -delt*2.0*pi/xxt*petit*un(ijl+1)

            vleft1=v(ijl)
            vleft1old=aa/2.*segma*sinh(xxk*y(j))/sinh(xxk*h0)*
     &          sin(pi/2.-segma*(t-delt)-cnf-xxk*(delx(1)/2.))
     &          +3./32.*2.*xxk*aa**2*segma*sinh(2.*xxk*y(j))
     &          /sinh(xxk*h0)**4
     &          *sin(pi/2-segma*(t-delt)-cnf-xxk*(delx(1)/2.))
            vleft2old=aa/2.*segma*sinh(xxk*y(j))/sinh(xxk*h0)*
     &          sin(pi/2.-segma*(t-delt)-cnf+xxk*(delx(2)/2.))
     &          +3./32.*2.*xxk*aa**2*segma*sinh(2.*xxk*y(j))
     &          /sinh(xxk*h0)**4
     &          *sin(pi/2-segma*(t-delt)-cnf+xxk*(delx(2)/2.))
            v(ijl)=c01*delt/((delx(1)+delx(2))/2)*
     &          ((vn(ijl+1)-vleft2old)-(vn(ijl)-vleft1old))
     &          +(vn(ijl)-vleft1old)+vleft1
	    endif
	  endif

C.......for linear periodic waves
        if (ninflow.eq.34) then
          u(ijl)=aa/2.*segma*cosh(xxk*yj(j))/sinh(xxk*h0)*
     &          cos(pi/2.-segma*t)
     &          -ulinear*areturn
C.........minus the mass accumalation each time step
          if (nweakref.gt.1) u(ijl)=u(ijl)-xmass_d/delt/(h0+eta1)

          v(ijl)=aa/2.*segma*sinh(xxk*y(j))/sinh(xxk*h0)*
     &          sin(pi/2.-segma*t-xxk*(delx(1)/2.))
C.........testing the weakly reflected boundary condition
          if (nweakref.eq.1.and.ncyc.ge.1) then
            uleft1=u(ijl)
            uleft1old=aa/2.*segma*cosh(xxk*yj(j))/sinh(xxk*h0)*
     &          cos(pi/2.-segma*(t-delt))
     &          -ulinear*areturn
            uleft2old=aa/2.*segma*cosh(xxk*yj(j))/sinh(xxk*h0)*
     &          cos(pi/2.-segma*(t-delt)+xxk*delx(2))
     &          -ulinear*areturn
            u(ijl)=c02*delt/delx(2)*((un(ijl+1)-uleft2old)-
     &          (un(ijl)-uleft1old))+
     &          (un(ijl)-uleft1old)+uleft1
C...............added a RHS from Petit et al (1994)
     &          -delt*2.0*pi/xxt*petit*un(ijl+1)

            vleft1=v(ijl)
            vleft1old=aa/2.*segma*sinh(xxk*y(j))/sinh(xxk*h0)*
     &          sin(pi/2.-segma*(t-delt)-xxk*(delx(1)/2.))
            vleft2old=aa/2.*segma*sinh(xxk*y(j))/sinh(xxk*h0)*
     &          sin(pi/2.-segma*(t-delt)+xxk*(delx(2)/2.))
            v(ijl)=c01*delt/((delx(1)+delx(2))/2)*
     &          ((vn(ijl+1)-vleft2old)-(vn(ijl)-vleft1old))
     &          +(vn(ijl)-vleft1old)+vleft1
          endif
        endif

C.......for cnoidal wave; after Wiegel (1960) JFM
        if (ninflow.eq.24) then
	  ccn1=cn(cnf*0.5*ck(mod1)+2.0*ck(mod1)*(-t/xxt),mod1)
	  ccn2=cn(cnf*0.5*ck(mod1)+2.0*ck(mod1)*(-delx(1)/2/xxl-t/xxt),
     &		mod1)
	  csn1=sqrt(1.0-ccn1**2)
	  cdn1=sqrt(1.0-mod1*(1-ccn1**2))
	  csn2=sqrt(1.0-ccn2**2)
	  cdn2=sqrt(1.0-mod1*(1-ccn2**2))
C.........it is found Wiegel's vertical velocity solution only valid for 2K(k)
C.........modification made to correct it
          ftest=mod(cnf*0.5*ck(mod1)+2.0*ck(mod1)*(-delx(1)/2/xxl-
     &		t/xxt),4*ck(mod1))
          if (ftest.ge.2.0*ck(mod1).or.(ftest.ge.-2.0*ck(mod1).and.
     &		ftest.le.0.0)) cdn2=-cdn2

          u(ijl)=sqrt(-gy*h0)*(-5./4.+3.*(h0+ytrough)/2./h0-
     &	    (h0+ytrough)**2/4./
     &	    h0**2+(3.*aa/2./h0-(h0+ytrough)*aa/2./h0**2)*ccn1**2
     &	    -aa**2/4/h0**2*
     &	    ccn1**4-8*aa*ck(mod1)**2/xxl**2*(h0/3-yj(j)**2/2/h0)*
     &	    (-mod1*csn1**2*ccn1**2+ccn1**2*cdn1**2-csn1**2
     &	    *cdn1**2))
     &      -ulinear*areturn
C......minus the mass accumalation each time step
          if (nweakref.gt.1) u(ijl)=u(ijl)-xmass_d/delt/(h0+eta1)

          v(ijl)=sqrt(-gy*h0)*y(j)*2*aa*ck(mod1)/xxl/h0*(1+
     &	     (h0+ytrough)/h0+aa/h0*ccn2**2+32*ck(mod1)**2/3/xxl**2*
     &       (h0**2-y(j)**2/2)*(mod1*csn2**2-mod1*ccn2**2
     &	     -cdn2**2))*csn2*ccn2*cdn2

C.........testing the weakly reflected boundary condition
          if (nweakref.eq.1.and.ncyc.ge.1) then
            ccn1=cn(cnf*0.5*ck(mod1)+2.0*ck(mod1)*(-(t-delt)/xxt),mod1)
            ccn2=cn(cnf*0.5*ck(mod1)+2.0*ck(mod1)*
     &		(-delx(1)/2/xxl-(t-delt)/xxt),mod1)
            csn1=sqrt(1.0-ccn1**2)
            cdn1=sqrt(1.0-mod1*(1-ccn1**2))
            csn2=sqrt(1.0-ccn2**2)
            cdn2=sqrt(1.0-mod1*(1-ccn2**2))
            ccn3=cn(cnf*0.5*ck(mod1)+2.0*ck(mod1)*
     &		(+delx(2)/xxl-(t-delt)/xxt),mod1)
            ccn4=cn(cnf*0.5*ck(mod1)+2.0*ck(mod1)*
     &          (+delx(2)/2/xxl-(t-delt)/xxt),mod1)
            csn3=sqrt(1.0-ccn3**2)
            cdn3=sqrt(1.0-mod1*(1-ccn3**2))
            csn4=sqrt(1.0-ccn4**2)
            cdn4=sqrt(1.0-mod1*(1-ccn4**2))

C...........it is found Wiegle's velocity solution only valid for 2K(k)
C...........modification made to correct it
            ftest2=mod(cnf*0.5*ck(mod1)+2.0*ck(mod1)*(-delx(1)/2/xxl-
     &          (t-delt)/xxt),4*ck(mod1))
            if (ftest2.ge.2.0*ck(mod1).or.(ftest2.ge.-2.0*ck(mod1).and.
     &          ftest2.le.0.0)) cdn2=-cdn2
            ftest4=mod(cnf*0.5*ck(mod1)+2.0*ck(mod1)*(+delx(1)/2/xxl-
     &          (t-delt)/xxt),4*ck(mod1))
            if (ftest4.ge.2.0*ck(mod1).or.(ftest4.ge.-2.0*ck(mod1).and.
     &          ftest4.le.0.0)) cdn4=-cdn4

            uleft1=u(ijl)
            uleft1old=sqrt(-gy*h0)*(-5./4.+3.*(h0+ytrough)/2./h0-
     &      (h0+ytrough)**2/4./
     &      h0**2+(3.*aa/2./h0-(h0+ytrough)*aa/2./h0**2)*ccn1**2
     &      -aa**2/4/h0**2*
     &      ccn1**4-8*aa*ck(mod1)**2/xxl**2*(h0/3-yj(j)**2/2/h0)*
     &      (-mod1*csn1**2*ccn1**2+ccn1**2*cdn1**2-csn1**2
     &      *cdn1**2))
     &      -ulinear*areturn
            uleft2old=sqrt(-gy*h0)*(-5./4.+3.*(h0+ytrough)/2./h0-
     &      (h0+ytrough)**2/4./
     &      h0**2+(3.*aa/2./h0-(h0+ytrough)*aa/2./h0**2)*ccn3**2
     &      -aa**2/4/h0**2*
     &      ccn3**4-8*aa*ck(mod1)**2/xxl**2*(h0/3-yj(j)**2/2/h0)*
     &      (-mod1*csn3**2*ccn3**2+ccn3**2*cdn3**2-csn3**2
     &      *cdn3**2))
     &      -ulinear*areturn
            u(ijl)=c02*delt/delx(2)*((un(ijl+1)-uleft2old)-
     &          (un(ijl)-uleft1old))+
     &          (un(ijl)-uleft1old)+uleft1
C...............added a RHS from Petit et al (1994)
     &          -delt*2.0*pi/xxt*petit*un(ijl+1)

            vleft1=v(ijl)
            vleft1old=sqrt(-gy*h0)*y(j)*2*aa*ck(mod1)/xxl/h0*(1+
     &       (h0+ytrough)/h0+aa/h0*ccn2**2+32*ck(mod1)**2/3/xxl**2*
     &       (h0**2-y(j)**2/2)*(mod1*csn2**2-mod1*ccn2**2
     &       -cdn2**2))*csn2*ccn2*cdn2
            vleft2old=sqrt(-gy*h0)*y(j)*2*aa*ck(mod1)/xxl/h0*(1+
     &       (h0+ytrough)/h0+aa/h0*ccn4**2+32*ck(mod1)**2/3/xxl**2*
     &       (h0**2-y(j)**2/2)*(mod1*csn4**2-mod1*ccn4**2
     &       -cdn4**2))*csn4*ccn4*cdn4
            v(ijl)=c01*delt/((delx(1)+delx(2))/2)*
     &          ((vn(ijl+1)-vleft2old)-(vn(ijl)-vleft1old))
     &          +(vn(ijl)-vleft1old)+vleft1
	  endif
C.........after Dean and Dalrymple (1991) (proven to be wrong)
C	  u(ijl)=xxc*(eta1/h0-(eta1/h0)**2+3./2./mod1*(1./6.-
C     &		(yj(j)/h0)**2)*((1.+sqrt(mod1))*aa**2/h0**2+2.*
C     &		(2.*sqrt(mod1)-1.)
C     &		*(eta1/h0)**2-3.*mod1*(eta1/h0)**4*h0**2))
C	  dedx=-aa*sqrt(3.*aa/h0)*xxl/2/pi/h0*ccn2*sqrt(1.-ccn2**2)
C     &		*sqrt(1.+sqrt(mod1)*(ccn2**2-1))
C	  v(ijl)=-xxc*(dedx*y(j)/h0*(-2.*eta1/h0+1)+3./h0/mod1*
C     &	 	(1.-(y(j)/h0)**2)*((2*sqrt(mod1)-1.)*aa/h0-3.*ck(mod1)
C     &		*sqrt(mod1)*eta1))

C.........shallow water assumption
C	  u(ijl)=xxc*eta1/(eta1+h0)
C	  v(ijl)=0.0
	endif

C.......for specified boundary condition
        if (ninflow.eq.9) then
	    kbc=1
	    mbc=1
689       continue
          if (ybc(it,kbc+1).eq.0.0) then
            if (y(j-1).le.(eta1+h0)) then
                u(ijl)=ubc(it,kbc)
            else
                u(ijl)=0.0
            endif
            goto 688
          endif
          if (yj(j).lt.ybc(it,kbc)+h0) then
                uit=ubc(it,kbc)
                uitp=ubc(it+1,kbc)
                if (t.lt.t_pad(it)) then
                  u(ijl)=uit
                else
                  u(ijl)=(uit*(t_pad(it+1)-t)+uitp*(t-t_pad(it)))
     &            /(t_pad(it+1)-t_pad(it))
                endif
                goto 688
          endif
          if (yj(j).gt.ybc(it,kbc)+h0.and.yj(j).le.ybc(it,kbc+1)+h0)
     &          then
                uit=(ubc(it,kbc)*(ybc(it,kbc+1)+h0-yj(j))+ubc(it,kbc+1)
     &          *(yj(j)-ybc(it,kbc)-h0))/(ybc(it,kbc+1)-ybc(it,kbc))
                uitp=(ubc(it+1,kbc)*(ybc(it+1,kbc+1)+h0-yj(j))+
     &          ubc(it+1,kbc+1)*(yj(j)-ybc(it+1,kbc)-h0))
     &          /(ybc(it+1,kbc+1)-ybc(it+1,kbc))
                if (yj(j).gt.ybc(it+1,kbc)+h0.and.yj(j).le.ybc(it+1,
     &          kbc+1)+h0.and.t.ge.t_pad(it)) then
                   u(ijl)=(uit*(t_pad(it+1)-t)+uitp*(t-t_pad(it)))
     &          /(t_pad(it+1)-t_pad(it))
                else
                  u(ijl)=uit
                endif
                goto 688
          else
                kbc=kbc+1
                goto 689
          endif
688       continue

789       continue
          if (ybc(it,mbc+1).eq.0.0) then
            if (y(j-1).le.(eta1+h0)) then
                v(ijl)=vbc(it,mbc)
            else
                v(ijl)=0.0
            endif
            goto 788
          endif
          if (y(j).lt.ybc(it,mbc)+h0) then
                vit=vbc(it,mbc)
                vitp=vbc(it+1,mbc)
                if (t.lt.t_pad(it)) then
                  v(ijl)=vit
                else
                  v(ijl)=(vit*(t_pad(it+1)-t)+vitp*(t-t_pad(it)))
     &            /(t_pad(it+1)-t_pad(it))
                endif
                goto 788
          endif
          if (y(j).gt.ybc(it,mbc)+h0.and.y(j).le.ybc(it,mbc+1)+h0)
     &          then
                vit=(vbc(it,mbc)*(ybc(it,mbc+1)+h0-y(j))+vbc(it,mbc+1)
     &          *(y(j)-ybc(it,mbc)-h0))/(ybc(it,mbc+1)-ybc(it,mbc))
                vitp=(vbc(it+1,mbc)*(ybc(it+1,mbc+1)+h0-y(j))+
     &          vbc(it+1,mbc+1)*(y(j)-ybc(it+1,mbc)-h0))
     &          /(ybc(it+1,mbc+1)-ybc(it+1,mbc))
                if (y(j).gt.ybc(it+1,mbc)+h0.and.y(j).le.ybc(it+1,
     &          mbc+1)+h0.and.t.ge.t_pad(it)) then
                  v(ijl)=(vit*(t_pad(it+1)-t)+vitp*(t-t_pad(it)))
     &          /(t_pad(it+1)-t_pad(it))
                else
                  v(ijl)=vit
                endif
                goto 788
          else
                mbc=mbc+1
                goto 789
          endif
788       continue
        endif
	  
C.......Solitary wave from Boussinesq equation (Lee et al, 1982)
	  if (ninflow.eq.5) then
  	    atmp=sqrt(3.0/4.0*aa/h0**3)
          xxx1=xstart-c1*t
	    c2=sqrt(-gy*h0)
	    d2edx2=2.0*aa*atmp**2*(2.0*cosh(atmp*xxx1)**2-3)
     &		/cosh(atmp*xxx1)**4
	    u(ijl)=c2*eta1/h0*(1.0-1.0/4.0*eta1/h0+h0/3.0*(h0/eta1)*
     &	  (1.0-3.0/2.0*yj(j)**2/h0**2)*d2edx2)

          xxx2=xstart-delx(1)/2.0-c1*t
	    dedx=-2.0*aa*sinh(atmp*xxx2)*atmp/cosh(atmp*xxx2)**3
	    d3edx3=-8.0*aa*sinh(atmp*xxx2)*atmp**3*
     &	  (cosh(atmp*xxx2)**2-3.0)/cosh(atmp*xxx2)**5
	    v(ijl)=-c2*y(j)/h0*((1.0-1.0/2.0*eta0/h0)*dedx+1.0/3.0*h0**2*
     &	  (1.0-1.0/2.0*y(j)**2/h0**2)*d3edx3)

C.........testing the weakly reflected boundary condition
          if (nweakref.eq.1.and.ncyc.ge.1) then
            xxx1=xstart-c1*(t-delt)
            d2edx2=2.0*aa*atmp**2*(2.0*cosh(atmp*xxx1)**2-3)
     &          /cosh(atmp*xxx1)**4
            xxx1a=xstart+delx(2)-c1*(t-delt)
            d2edx2a=2.0*aa*atmp**2*(2.0*cosh(atmp*xxx1a)**2-3)
     &          /cosh(atmp*xxx1a)**4
	      uleft1=u(ijl)
            uleft1old=c2*eta1old/h0*(1.0-1.0/4.0*eta1old/h0+
     &	    h0/3.0*(h0/eta1old)*
     &      (1.0-3.0/2.0*yj(j)**2/h0**2)*d2edx2)
            uleft2old=c2*eta2aold/h0*(1.0-1.0/4.0*eta2aold/h0+
     &      h0/3.0*(h0/eta2aold)*
     &      (1.0-3.0/2.0*yj(j)**2/h0**2)*d2edx2a)
            u(ijl)=c02*delt/delx(2)*((un(ijl+1)-uleft2old)-
     &		(un(ijl)-uleft1old))+
     &          (un(ijl)-uleft1old)+uleft1

            xxx2=xstart-delx(1)/2.0-c1*(t-delt)
            xxx2a=xstart+delx(2)/2.0-c1*(t-delt)
            dedx=-2.0*aa*sinh(atmp*xxx2)*atmp/cosh(atmp*xxx2)**3
            dedxa=-2.0*aa*sinh(atmp*xxx2)*atmp/cosh(atmp*xxx2a)**3
            d3edx3=-8.0*aa*sinh(atmp*xxx2)*atmp**3*
     &      (cosh(atmp*xxx2)**2-3.0)/cosh(atmp*xxx2)**5
            d3edx3a=-8.0*aa*sinh(atmp*xxx2a)*atmp**3*
     &      (cosh(atmp*xxx2a)**2-3.0)/cosh(atmp*xxx2a)**5
	      vleft1=v(ijl)
            vleft1old=-c2*y(j)/h0*((1.0-1.0/2.0*eta0old/h0)*dedx+1.0
     &	    /3.0*h0**2*(1.0-1.0/2.0*y(j)**2/h0**2)*d3edx3)
            vleft2old=-c2*y(j)/h0*((1.0-1.0/2.0*eta2old/h0)*dedxa+1.0
     &      /3.0*h0**2*(1.0-1.0/2.0*y(j)**2/h0**2)*d3edx3a)
            v(ijl)=c01*delt/((delx(1)+delx(2))/2)*
     &          ((vn(ijl+1)-vleft2old)-(vn(ijl)-vleft1old))
     &		+(vn(ijl)-vleft1old)+vleft1
          endif
C.........it is found defining v(ijl+1) always cause trouble to satisfy
C.........continuity problem, though it is from analytical solution
C          xxx0=xstart+delx(2)/2.0-c1*t
C          dedxa=-2.0*aa*sinh(atmp*xxx0)*atmp/cosh(atmp*xxx0)**3
C          d3edx3a=-8.0*aa*sinh(atmp*xxx0)*atmp**3*
C     &	   (cosh(atmp*xxx0)**2-3.0)/cosh(atmp*xxx0)**5
C          v(ijl+1)=-c2*y(j)/h0*((1.0-1.0/2.0*eta2/h0)*dedxa+1.0/3.0*
C     &	   h0**2*(1.0-1.0/2.0*y(j)**2/h0**2)*d3edx3a)
	  endif

C.......specify the VOF based on the given or calculated free surface
	  if (y(j).le.eta0x+h0) then
		f(ijl)=1.0
	  else 
	     if (y(j-1).ge.eta0x+h0) then
	     	f(ijl)=0.00
			if (f(ijl+1).eq.0.0) u(ijl)=0.00
			v(ijl)=0.00
			p(ijl)=0.0
	     else
			f(ijl)=(eta0x+h0-y(j-1))/dely(j)
	     endif
	  endif

C.......specify the inflow turbulence information 
        if (kemodel.gt.0) then
                xk(ijl)=0.5*uturb**2
                xnut(ijl)=5.00*xnu*ratio**2
                xep(ijl)=C_mu*xk(ijl)**2/xnut(ijl)
        endif

        go to 50

   10   continue
	  u(ijl)=0.0
        v(ijl)=v(ijl+1)
	  go to 50

   20   u(ijl)=0.0
        v(ijl)=-v(ijl+1)*delx(1)/delx(2)
        go to 50

   30   continue
C	  if (iter.gt.0.and.nnn.eq.0) go to 50
C.......should be able to specify better open flow boundary condition using
C.......extrapolation. The currently one is more like the zero-gradient 
C.......(solid wall) that would cause excessive reflection
	  if (nopen.eq.0) then
            u(ijl)=u(ijl+1)*(x(2)*rx(1)*cyl+1.0-cyl)
            v(ijl)=v(ijl+1)
		  f(ijl)=f(ijl+1)
		  if (f(ijl).le.1.0e-6) u(ijl)=0.0
	  endif
C.......open boundary condition
C.......Noted that for plane jet problem, the left boundary must be set to
C.......nopen=1; the conventional nopen=0 may accumalate k to very high
C.......level and cause instability (6/21/99) 
        if ((nopen.eq.1.or.nopen.eq.11.or.nopen.eq.2).and.ncyc.ge.1)then
	    if (nopen.eq.1.or.nopen.eq.11) then
            u(ijl)=c02*delt/delx(2)*
     &          (un(ijl+1)-un(ijl))+un(ijl)
            v(ijl)=c01*delt/((delx(1)+delx(2))/2.0)*
     &          (vn(ijl+1)-vn(ijl))+vn(ijl)
	      if (kemodel.gt.0.and.f(ijl).ne.0.0) then
              xk(ijl)=xk(ijl+1)
              xep(ijl)=xep(ijl+1)
              xnut(ijl)=xnut(ijl+1)
	      endif
		else
	      u(ijl)=u(ijl+1)
		  v(ijl)=v(ijl+1)
		  if (kemodel.gt.0.and.f(ijl).ne.0.0) then
              xk(ijl)=xk(ijl+1)
              xep(ijl)=xep(ijl+1)
              xnut(ijl)=xnut(ijl+1)
	      endif
		endif

          if (y(j).le.eta0x+flht) then
                f(ijl)=1.0
          else
              if (y(j-1).ge.eta0x+flht) then
                f(ijl)=0.00
              else
                f(ijl)=(eta0x+flht-y(j-1))/dely(j)
              endif
          endif
		if (nopen.eq.2) f(ijl+1)=f(ijl)
        endif
        go to 50

   40   continue
        f(ijl)=f(ijr-2)
  	  u(ijl)=u(ijr-2)
        v(ijl)=v(ijr-2)
	  p(ijl)=p(ijr-2)

   50   continue
C.......start the specified flow without free surface (e.g., pipe flows and jet)
C.......Since it can be combined with any kl=1,2,3,and even 6 (not done yet)
C.......it is put at the end of boundary condition

C.......specify plane jet (uniform velocity) from the left wall
        if (ninflow.eq.7) then
         if(t.le.time_jet)then
          if (yj(j).gt.ylow.and.yj(j).le.yup) then
            f((j-1)*imax+1)=1.0
            u((j-1)*imax+1)=ujet
            v((j-1)*imax+1)=vjet
            v((j-1)*imax+1-imax)=vjet
          endif
          endif
	  endif
C ... jet and take out
        if (ninflow.eq.10) then
         if(t.le.time_jet)then
          if (yj(j).gt.ylowout.and.yj(j).le.yupout) then
c left
            f((j-1)*imax+1)=1.0
            u((j-1)*imax+1)=uout
            v((j-1)*imax+1)=vout
            v((j-1)*imax+1-imax)=vout
c right should be right boundary condition below fyshi
c            f((j-1)*imax+imax)=1.0
c            u((j-1)*imax+imax)=-uout
c            v((j-1)*imax+imax)=vout
c            v((j-1)*imax)=vout
          endif
          endif
          endif


C.......specify plane jet (parabolic shape) from the left wall
        if (ninflow.eq.70) then
          if (yj(j).gt.ylow.and.yj(j).le.yup) then
		  ijl=(j-1)*imax+1
            f(ijl)=1.0
            u(ijl)=-6.0*ujet/wjet**2*(yj(j)-yjet)**2+1.5*ujet
            v(ijl)=-6.0*vjet/wjet**2*(y(j)-yjet)**2+1.5*vjet
            v(ijl-imax)=-6.0*vjet/wjet**2*(y(j-1)-yjet)**2+1.5*vjet
          endif
	  endif
C.......specify plane jet (log-law profile) from the left wall
        if (ninflow.eq.71) then
          if (yj(j).gt.ylow.and.yj(j).le.yup) then
	      ijl=(j-1)*imax+1
            f(ijl)=1.0
	      if ((wjet/2-abs(yj(j)-yjet))/0.0001.le.1.0) then
		    u(ijl)=0.0
		  else
			 u(ijl)=ustar*(5.75*log10((wjet/2-abs(yj(j)-yjet))/0.0001)
     &			+8.5)*ujet/sqrt(ujet**2+vjet**2)
            endif
		  if ((wjet/2-abs(y(j)-yjet))/0.0001.le.1.0) then
			  v(ijl)=0.0
		  else
			  v(ijl)=ustar*(5.75*log10((wjet/2-abs(y(j)-yjet))/0.0001)
     &			+8.5)*vjet/sqrt(ujet**2+vjet**2)
		  endif
		  if ((wjet/2-abs(y(j-1)-yjet))/0.0001.le.1.0) then
			  u(ijl-imax)=0.0
			else
			  v(ijl-imax)=ustar*(5.75*log10((wjet/2-abs(y(j-1)-yjet))
     &		  /0.0001)+8.5)*vjet/sqrt(ujet**2+vjet**2)
		  endif
		endif
	  endif

	  go to (60,70,80,90,85,55), kr
   85	  continue
	  p(ijr)=pbc(4)+gy*yj(j)
	  go to 95
   55   u(ijr-1)=uinf(4)
        v(ijr)=v(ijr-1)
        go to 95
   60   u(ijr-1)=0.0
        v(ijr)=v(ijr-1)
        go to 95
   70   u(ijr-1)=0.0
        v(ijr)=-v(ijr-1)*delx(imax)/delx(im1)
        go to 95
   80   continue
	  if (nopen.eq.0) then
          u(ijr-1)=u(ijr-2)*(x(ibar)*rx(im1)*cyl+1.0-cyl)
          v(ijr)=v(ijr-1)
		f(ijr)=f(ijr-1)
		if (f(ijr).le.1.0e-6) u(ijr-1)=0.0
	  endif
C.......testing the open boundary condition
        if ((nopen.eq.1.or.nopen.eq.11.or.nopen.eq.2).and.ncyc.ge.1)then
	    if (nopen.eq.1.or.nopen.eq.11) then
            u(ijr-1)=-cim1*delt/delx(im1)*
     &          (un(ijr-1)-un(ijr-2))+un(ijr-1)
            v(ijr)=-crbc*delt/((delx(im1)+delx(imax))/2.0)*
     &          (vn(ijr)-vn(ijr-1))+vn(ijr)
C.......it is found (12/4/98) that when very large domain (2000 horizontal 
C.......grid) is used, the right open boundary condition becomes unstable
C.......after certain time even before the wave train arrives. It is
C.......suspected that it is caused by the errors in pressure solver.
C.......Currently using small delt to restrict the instability (either
C.......constant delt or small dtmax).
		  if (kemodel.gt.0.and.f(ijr).ne.0.0) then
              xk(ijr)=xk(ijr-1)  
              xep(ijr)=xep(ijr-1) 
              xnut(ijr)=xnut(ijr-1) 
	      endif
		else
	      u(ijr-1)=u(ijr-2)
		  v(ijr)=v(ijr-1)
	      if (kemodel.gt.0.and.f(ijr).ne.0.0) then
              xk(ijr)=xk(ijr-1)  
              xep(ijr)=xep(ijr-1) 
              xnut(ijr)=xnut(ijr-1)
	      endif
		endif
          if (y(j).le.etaimax+flht) then
              f(ijr)=1.0
          else
              if (y(j-1).ge.etaimax+flht) then
                f(ijr)=0.00
              else
                f(ijr)=(etaimax+flht-y(j-1))/dely(j)
              endif
          endif
	  endif
	  if (nopen.eq.2) f(ijr-1)=f(ijr) 
        go to 95
   90   continue
        f(ijr-1)=f(ijl+1)
        f(ijr)=f(ijl+2)
	  u(ijr-1)=u(ijl+1)
        v(ijr-1)=v(ijl+1)
        p(ijr-1)=p(ijl+1)
        u(ijr)=u(ijl+2)
        v(ijr)=v(ijl+2)
	  p(ijr)=p(ijl+2)
   95   continue
	  ijl=ijl+imax
        ijr=ijr+imax

C ... jet and take out
        if (ninflow.eq.10) then
         if(t.le.time_jet)then
          if (yj(j).gt.ylowout.and.yj(j).le.yupout) then
c left should be left condition above fyshi
c            f((j-1)*imax+1)=1.0
c            u((j-1)*imax+1)=uout
c            v((j-1)*imax+1)=vout
c            v((j-1)*imax+1-imax)=vout
c right
            f((j-1)*imax+imax-1)=1.0
            u((j-1)*imax+imax-1)=-uout
            v((j-1)*imax+imax)=vout
            v((j-1)*imax)=vout
          endif
          endif
          endif



  100 continue
c
      ijb=1
      ijt=imax*jm1+1
      do 200 i=1,imax
C.......modified to account for the impact
        if (nf(ijb+imax).ne.4) then
		f(ijb)=f(ijb+imax)
	    p(ijb)=p(ijb+imax)
	    nf(ijb)=nf(ijb+imax)
	  else
		f(ijb)=0.0
          p(ijb)=0.0
          nf(ijb)=6
	  endif
	  if (kb.eq.1) nf(ijb)=0
        if (nf(ijt-imax).ne.3) then
		f(ijt)=f(ijt-imax)
	    p(ijt)=p(ijt-imax)
		nf(ijt)=nf(ijt-imax)
	  else
		f(ijt)=0.0
          p(ijt)=0.0
          nf(ijt)=6
	  endif
	  if (kt.eq.1) nf(ijt)=0

        go to (110,120,130,140,135,105,106), kt
  135	  continue
	  p(ijt)=pbc(2)
	  go to 150
  105   v(ijt-imax)=vinf(2)
        u(ijt)=u(ijt-imax)
        go to 150
  110   v(ijt-imax)=0.0
        u(ijt)=u(ijt-imax)
        go to 150
  120   v(ijt-imax)=0.0
        u(ijt)=-u(ijt-imax)*dely(jmax)/dely(jm1)
        go to 150
  130   continue
        v(ijt-imax)=v(ijt-2*imax)
        u(ijt)=u(ijt-imax)
        go to 150
  140   continue
        f(ijt-imax)=f(ijb+imax)
        f(ijt)=f(ijb+2*imax)
	  if (iter.gt.0.and.nnn.eq.0) go to 150
	  v(ijt-imax)=v(ijb+imax)
        u(ijt-imax)=u(ijb+imax)
        p(ijt-imax)=p(ijb+imax)
        u(ijt)=u(ijb+2*imax)
	  go to 150        
C.......add the slip boundary condition for cavity case
  106   if (t.gt.h0) then
                xmu=0.0
                xnu=0.0
        else
        endif
        if (t.gt.aa) then   
                kt=1
                goto 110
        else
                v(ijt-imax)=0.0
		u(ijt)=utop
                goto 150
        endif
  150   continue
C.......start the specified flow without free surface (e.g., pipe flows and jet)
C.......Since it can be combined with any kt=1,2,3,it is put at the end of boundary condition

C.......specify plane jet (uniform velocity) from the top wall
        if (ninflow.eq.8) then
          if (xi(i).gt.ylow.and.xi(i).le.yup) then
		    f((jmax-1)*imax+i)=1.0
              u((jmax-1)*imax+i)=ujet
	  	    u((jmax-1)*imax+i-1)=ujet
              v((jmax-2)*imax+i)=vjet
	    endif
	  endif
C ... jet and take out fyshi
        if (ninflow.eq.10) then
          if (xi(i).gt.ylow.and.xi(i).le.yup) then
                    f((jmax-1)*imax+i)=1.0
              u((jmax-1)*imax+i)=ujet
                    u((jmax-1)*imax+i-1)=ujet
              v((jmax-2)*imax+i)=vjet
            endif
          endif
C.......specify plane jet (parabolic shape) from the top wall
        if (ninflow.eq.80) then
          if (xi(i).gt.ylow.and.xi(i).le.yup) then
	  	    f((jmax-1)*imax+i)=1.0
		    u((jmax-1)*imax+i)=-6.0*ujet/wjet**2*(x(i)-yjet)**2+1.5*ujet
		    u((jmax-1)*imax+i-1)=-6.0*ujet/wjet**2*(x(i-1)-yjet)**2
     &		  +1.5*ujet
              v((jmax-2)*imax+i)=-6.0*vjet/wjet**2*(xi(i)-yjet)**2
     &	 	  +1.5*vjet
          endif
	  endif
C.......specify plane jet (log-law profile) from the top wall
        if (ninflow.eq.81) then
          if (xi(i).gt.ylow.and.xi(i).le.yup) then
		    f((jmax-1)*imax+i)=1.0
		    if ((wjet/2-abs(x(i)-yjet))/0.0001.le.1.0) then
			  u((jmax-1)*imax+i)=0.0
		    else
			  u((jmax-1)*imax+i)=ustar*(5.75*log10((wjet/2-abs(x(i)
     &			-yjet))/0.0001)+8.5)*ujet/sqrt(ujet**2+vjet**2)
              endif
		    if ((wjet/2-abs(x(i-1)-yjet))/0.0001.le.1.0) then
			  u((jmax-1)*imax+i-1)=0.0
			else
			  u((jmax-1)*imax+i-1)=ustar*(5.75*log10((wjet/2-
     &		 abs(x(i-1)-yjet))/0.0001)+8.5)*ujet/sqrt(ujet**2+vjet**2)
			endif
		    if ((wjet/2-abs(xi(i)-yjet))/0.0001.le.1.0) then
			  v((jmax-2)*imax+i)=0.0
			else
			  v((jmax-2)*imax+i)=ustar*(5.75*log10((wjet/2
     &		 -abs(xi(i)-yjet))/0.0001)+8.5)*vjet/sqrt(ujet**2+vjet**2)
			endif
		endif
	  endif

	  go to (160,170,180,190,185,155,156), kb
  185	  continue
	  p(ijb)=pbc(1)
	  go to 195
  156   u(ijb+imax)=0.0
	  v(ijb)=0.0
	  go to 195 
  155   v(ijb)=vinf(1)
        u(ijb)=u(ijb+imax)
        go to 195
  160   v(ijb)=0.0
        u(ijb)=u(ijb+imax)
        go to 195
  170   v(ijb)=0.0
        u(ijb)=-u(ijb+imax)*dely(1)/dely(2)
        go to 195
  180   continue
        v(ijb)=v(ijb+imax)
        u(ijb)=u(ijb+imax)
        go to 195
  190   continue
        f(ijb)=f(ijt-2*imax)
	  if (iter.gt.0.and.nnn.eq.0) go to 195
	  v(ijb)=v(ijt-2*imax)
        u(ijb)=u(ijt-2*imax)
  195   ijt=ijt+1
        ijb=ijb+1
  200   continue

  850 continue
c
      if (ibcflg.eq.0) go to 9999
c.... free surface and sloped boundary conditions
c
      do 421 i=2,im1
        do 420 j=2,jm1
c
	  ij=(j-1)*imax+i
          ijm=ij-imax
          imj=ij-1
          ipj=ij+1
          ijp=ij+imax
          imjm=ij-1-imax
          imjp=ij-1+imax
          ipjm=ij-imax+1
          ipjp=ij+imax+1

C.........no free surface BC imposed on obstacle cell, which could be in motion
C	    if (ac(ij)*ar(ij)*ar(imj)*at(ij)*at(ijm).lt.1.0d0-em6)goto 420
		if (ac(ij).le.em6) goto 420
c
C.........added for impact purpose (it's found from the test of dam break
C.........that even with a thin layer of water on the right corner, this
C.........condition won't take effect; it plays the role only when the 
C.........corner is originally dry).
	    if (i.eq.2.and.nf(ij).eq.2.and.ibcfinal.eq.1.and.
     &		nfold(ij).ne.6) goto 420
          if (i.eq.im1.and.nf(ij).eq.1.and.ibcfinal.eq.1.and.
     &		nfold(ij).ne.6) goto 420
          if (j.eq.2.and.nf(ij).eq.4.and.ibcfinal.eq.1.and.
     &		nfold(ij).ne.6) goto 420
          if (j.eq.jm1.and.nf(ij).eq.3.and.ibcfinal.eq.1.and.
     &		nfold(ij).ne.6) goto 420
c
C.........special attention should be made to take care the cell jump
C.........from empty to full due to breaking
          if ((nf(ij).eq.0.and.nfold(ij).ne.6).
     &		or.(nf(ij).gt.5.and.nfold(ij).gt.5)) go to 420
          if (i.eq.2.and.kl.eq.5) go to 420
          if (i.eq.im1.and.kr.eq.5) go to 420
          if (j.eq.jm1.and.kt.eq.5) go to 420
          if (j.eq.2.and.kb.eq.5) go to 420
c	
          if (nnn.eq.0)  goto 235

C.........consider the new cell filling treatment
          if (nfold(ij).gt.5) then
C.........for isolate new filled-cell
		if (nf(ij).eq.5.or.nfs(ij).eq.5) then
			if (v(ij).eq.0.0.and.v(ijm).eq.0.0) then
			  if (u(ij).lt.0.0) then
				if (at(ij).gt.em6) v(ij)=vn(ipj)*
     &					at(ipj)/at(ij)
				xk(ij)=xk(ipj)
				xep(ij)=xep(ipj)
				xnut(ij)=xnut(ipj)
				xp(ij)=xp(ipj)
			  else
			     if (at(ij).gt.em6) v(ij)=vn(imj)*
     &				at(imj)/at(ij)
                   xk(ij)=xk(imj)
                   xep(ij)=xep(imj)
                   xnut(ij)=xnut(imj)
                   xp(ij)=xp(imj)
			  endif
			else
			  if (v(ij).lt.0.0) then
				if (ar(ij).gt.em6) u(ij)=un(ijp)*
     &				ar(ijp)/ar(ij)
                  xk(ij)=xk(ijp)
                  xep(ij)=xep(ijp)
                  xnut(ij)=xnut(ijp)
                  xp(ij)=xp(ijp)
			  else
				if (ar(ij).gt.em6) u(ij)=un(ijm)*
     &				ar(ijm)/ar(ij)
                  xk(ij)=xk(ijm)
                  xep(ij)=xep(ijm)
                  xnut(ij)=xnut(ijm)
                  xp(ij)=xp(ijm)
			  endif
			endif
              if (xk(ij).le.0.0.and.kemodel.gt.0) then
                                goto 75
              else
                                goto 234
              endif
		endif

C.........for surface cell
		if (nf(ij).gt.0.and.nf(ij).lt.5) then
C.................for vertical interface (consider volume effect)
                  if (nfold(ijp).gt.5.and.at(ij).gt.em6)
     &		  	v(ij)=cvmgt(v(ipj)*at(ipj)/at(ij),
     &		  	v(imj)*at(imj)/at(ij),
     &		  	fn(ipj)+fn(ipjp).gt.fn(imj)+fn(imjp))

                  if (nfold(ijm).gt.5.and.at(ijm).gt.em6) v(ijm)=
     &		  	cvmgt(v(ipjm)*at(ipjm)/at(ijm),
     &		  	v(imjm)*at(imjm)/at(ijm),
     &			fn(ipjm)+fn(ipj).gt.fn(imjm)+fn(imj)) 

C.................for horizontal interface (consider volume effect)
                  if (nfold(ipj).gt.5.and.ar(ij).gt.em6) u(ij)=
     &			cvmgt(u(ijp)*ar(ijp)/ar(ij),
     &			u(ijm)*ar(ijm)/ar(ij),
     &			fn(ijp)+fn(ipjp).gt.fn(ijm)+fn(ipjm))

                  if (nfold(imj).gt.5.and.ar(imj).gt.em6) u(imj)=
     &			cvmgt(u(imjp)*ar(imjp)/ar(imj),
     &			u(imjm)*ar(imjm)/ar(imj),
     &			fn(imjp)+fn(ijp).gt.fn(imjm)+fn(ijm))
	
                  if (kemodel.eq.0) goto 234
                  go to (71,72,73,74), nf(ij)
71                if (xk(imj).eq.0.0) goto 72
                  xk(ij)=xk(imj)
                  xep(ij)=xep(imj)
                  xnut(ij)=xnut(imj)
                  xp(ij)=xp(imj)
                  goto 234
72                if (xk(ipj).eq.0.0) goto 73
                  xk(ij)=xk(ipj)
                  xep(ij)=xep(ipj)
                  xnut(ij)=xnut(ipj)
                  xp(ij)=xp(ipj)
                  goto 234
73                if (xk(ijm).eq.0.0) goto 74
                  xk(ij)=xk(ijm)
                  xep(ij)=xep(ijm)
                  xnut(ij)=xnut(ijm)
                  xp(ij)=xp(ijm)
                  goto 234
74                if (xk(ijp).eq.0.0) goto 75
                  xk(ij)=xk(ijp)
                  xep(ij)=xep(ijp)
                  xnut(ij)=xnut(ijp)
                  xp(ij)=xp(ijp)
                  goto 234
		endif

C.........for cell jump from empty to full (consider volume effect)
		if (nf(ij).eq.0) then
C.................for vertical interface (consider volume effect)
                  if (nfold(ijp).gt.5.and.at(ij).gt.em6)
     &		  	v(ij)=cvmgt(v(ipj)*at(ipj)/at(ij),
     &		  	v(imj)*at(imj)/at(ij),
     &		  	fn(ipj)+fn(ipjp).gt.fn(imj)+fn(imjp))

                  if (nfold(ijm).gt.5.and.at(ijm).gt.em6) v(ijm)=
     &		  	cvmgt(v(ipjm)*at(ipjm)/at(ijm),
     &		  	v(imjm)*at(imjm)/at(ijm),
     &			fn(ipjm)+fn(ipj).gt.fn(imjm)+fn(imj)) 

C.................for horizontal interface (consider volume effect)
                  if (nfold(ipj).gt.5.and.ar(ij).gt.em6) u(ij)=
     &			cvmgt(u(ijp)*ar(ijp)/ar(ij),
     &			u(ijm)*ar(ijm)/ar(ij),
     &			fn(ijp)+fn(ipjp).gt.fn(ijm)+fn(ipjm))

                  if (nfold(imj).gt.5.and.ar(imj).gt.em6) u(imj)=
     &			cvmgt(u(imjp)*ar(imjp)/ar(imj),
     &			u(imjm)*ar(imjm)/ar(imj),
     &			fn(imjp)+fn(ijp).gt.fn(imjm)+fn(ijm))


				if (kemodel.gt.0.and.npc(ij).eq.1) then
				xk(ij)=max(max(xk(ipj),xk(imj)),
     &			max(xk(ijp),xk(ijm)))
				xep(ij)=max(max(xep(ipj),xep(imj)),
     &			max(xep(ijp),xep(ijm)))
				xnut(ij)=C_mu*xk(ij)**2/xep(ij)
                    xp(ij)=max(max(xp(ipj),xp(imj)),
     &			max(xp(ijp),xp(ijm)))
				goto 234
				endif
		endif

75        continue
		if (kemodel.gt.0.and.npc(ij).eq.1) then
  		  xk(ij)=max(max(xk(ipj),xk(imj)),
     &			max(xk(ijp),xk(ijm)),0.5*uturb**2)
                  xnut(ij)=max(max(xnut(ipj),xnut(imj)),
     &			max(xnut(ijp),xnut(ijm)),5.00*xnu*ratio**2)
                  xep(ij)=C_mu*xk(ij)**2/xnut(ij)
                  xp(ij)=max(max(xp(ipj),xp(imj)),
     &			max(xp(ijp),xp(ijm)))
                  if (xk(ij).eq.0.0) then
                    xk(ij)=max(max(max(xk(ipj),xk(imj)),
     &				max(xk(ijp),xk(ijm))),
     &                max(max(xk(ijp+1),xk(ijp-1)),
     &				max(xk(ijm+1),xk(ijm-1))))
                    xep(ij)=max(max(max(xep(ipj),xep(imj)),
     &				max(xep(ijp),xep(ijm))),
     &                max(max(xep(ijp+1),xep(ijp-1)),
     &				max(xep(ijm+1),xep(ijm-1))))
                    xnut(ij)=C_mu*xk(ij)**2/xep(ij)
                    xp(ij)=max(max(max(xp(ipj),xp(imj)),
     &				max(xp(ijp),xp(ijm))),
     &                max(max(xp(ijp+1),xp(ijp-1)),
     &				max(xp(ijm+1),xp(ijm-1))))
                  endif
		endif

234		continue
C.........when porous media is present, it is possible to create the 
C.........small pressure distrubance near the shoreline which
C.........generate the flow nearby; this treatment accounts for
C.........this special case
		if (npor.ne.0.and.kemodel.gt.0) then
		  if (min(ac(ijp),ac(ijm),ac(ipj),ac(imj)).le.em6) then
                    xk(ij)=max(xk(ipj),xk(imj),xk(ijp),xk(ijm),
     &                xk(ijp+1),xk(ijp-1),xk(ijm+1),xk(ijm-1),
     &				0.5*uturb**2)
                    xnut(ij)=max(xnut(ipj),xnut(imj),xnut(ijp),
     &              xnut(ijm),xnut(ijp+1),xnut(ijp-1),
     &              xnut(ijm+1),xnut(ijm-1),5.00*xnu*ratio**2)
                    xep(ij)=C_mu*xk(ij)**2/xnut(ij)
		  endif
		endif

          	if ((xk(ij).le.1.0e-16.or.xep(ij).le.1.0e-16).
     &			and.kemodel.gt.0.and.npc(ij).eq.1) then
                  	write(9,*)t,i,j,'new filled cell failed',
     &			xk(ij),xep(ij),f(ij),f(ijp),f(ijm),f(ipj),
     &			f(imj),nf(ij),nf(ijp),nf(ijm),nf(ipj),nf(imj),
     &			nfold(ij),nfold(ijp),nfold(ijm),nfold(ipj),
     &			nfold(imj),nfs(ij),nfs(ijp),nfs(ijm),nfs(ipj),
     &           	nfs(imj),xk(ijp),xk(ijm),xk(ipj),xk(imj),
     &                  f(ijp+1),f(ijp-1),f(ijm+1),f(ijm-1),
     &                  xk(ijp+1),xk(ijp-1),xk(ijm+1),xk(ijm-1)

                  	xk(ij)=0.5*uturb**2
                  	xnut(ij)=5.00*xnu*ratio**2
                  	xep(ij)=C_mu*xk(ij)**2/xnut(ij)
          	endif
		goto 235
	  else
	  endif

C.......for newly emptied cell
	  if (nfold(ij).le.5.and.nf(ij).gt.5) then
	    if (nf(ipj).gt.5) then
			u(ij)=0.0
              if (j.eq.2) u(ijm)=0.0
              if (j.eq.jm1) u(ijp)=0.0
	    endif
	    if (nf(imj).gt.5) then
			u(imj)=0.0
			if (j.eq.2) u(imjm)=0.0
			if (j.eq.jm1) u(imjp)=0.0
	    endif
	    if (nf(ijp).gt.5) then
			v(ij)=0.0
              if (i.eq.2) v(imj)=0.0
              if (i.eq.im1) v(ipj)=0.0
	    endif
	    if (nf(ijm).gt.5) then
			v(ijm)=0.0
              if (i.eq.2) v(imjm)=0.0
              if (i.eq.im1) v(ipjm)=0.0
	    endif
          if (i.eq.2) then
                xk(imj)=xk(ij)
                xep(imj)=xep(ij)
                xnut(imj)=xnut(ij)
                xp(imj)=xp(ij)
          endif
          if (i.eq.im1) then
                xk(ipj)=xk(ij)
                xep(ipj)=xep(ij)
                xnut(ipj)=xnut(ij)
                xp(ipj)=xp(ij)
          endif

          if (j.eq.jm1) then
                xk(ijp)=xk(ij)
                xep(ijp)=xep(ij)
                xnut(ijp)=xnut(ij)
                xp(ijp)=xp(ij)
          endif
          if (j.eq.2) then
                xk(ijm)=xk(ij)
                xep(ijm)=xep(ij)
                xnut(ijm)=xnut(ij)
                xp(ijm)=xp(ij)
          endif

	    goto 420
	  endif
235     continue

C.........skip if the cell is previously fluid cell
C          if (nnn.eq.1.and.nfold(ij).eq.0) goto 250

          if (nf(ipj).lt.6.or.(ar(ij).lt.em6)) goto 260
          if (nf(imj)+nf(imjp).eq.0.and.nf(ijp).ne.0.and.nf(ipjp).eq.6.
     &		and.(at(ij).gt.1.0-em6).and.(at(imj).gt.1.0-em6)) then
          	if (nf(ijp)*nf(ij).eq.0) goto 2600
C...............zero tangential gradient; close to zero stress
              v(ij)=v(imj)*cvmgt(1.0d0,at(imj),at(imj).le.em6)
     &		/cvmgt(1.0d0,at(ij),at(ij).le.em6)
		endif	
2600      if (nf(imj)+nf(imjm).eq.0.and.nf(ijm).ne.0.and.nf(ipjm).eq.6.
     &		and.(at(ijm).gt.1.0-em6).and.(at(imjm).gt.1.0-em6)) then
              if (nf(ijm)*nf(ij).eq.0) goto 260
C...............zero tangential gradient; close to zero stress
              v(ijm)=v(imjm)*cvmgt(1.0d0,at(imjm),at(imjm).le.em6)
     &		/cvmgt(1.0d0,at(ijm),at(ijm).le.em6)
          endif
C          if (nf(imj).ne.6) 
C     		u(ij)=u(imj)*cvmgt(1.0d0,ar(imj),ar(imj).le.em6)
C     &		*r(i-1)/(cvmgt(1.0d0,ar(ij),ar(ij).le.em6)*r(i))
C		u(ij)=uxmb(nmovbd(ij))+(u(imj)-uxmb(nmovbd(ij)))*ar(imj)
C     &		/ar(ij)
C.........By such doing, it will slow down momentum transfer on slope; physically happens
		if (nf(imj).eq.6) then
			u(ij)=cvmgt(u(ij),u(imj),abs(u(ij)).lt.abs(u(imj)))
		else
			u(ij)=u(imj)
		endif

  260     continue
          if(nf(ijp).lt.6.or.(at(ij).lt.em6)) goto 270
          if (nf(ijm)+nf(ipjm).eq.0.and.nf(ipj).ne.0.and.nf(ipjp).eq.6.
     &		and.(ar(ij).gt.1.0-em6).and.(ar(ijm).gt.1.0-em6)) then
              if (nf(ipj)*nf(ij).eq.0) goto 2700
C...............zero tangential gradient; close to zero stress
              u(ij)=u(ijm)*cvmgt(1.0d0,ar(ijm),ar(ijm).le.em6)
     &		/cvmgt(1.0d0,ar(ij),ar(ij).le.em6)
          endif
2700      if (nf(ijm)+nf(imjm).eq.0.and.nf(imj).ne.0.and.nf(imjp).eq.6.
     &		and.(ar(imj).gt.1.0-em6).and.(ar(imjm).gt.1.0-em6)) then
              if (nf(imj)*nf(ij).eq.0) goto 270
C...............zero tangential gradient; close to zero stress
              u(imj)=u(imjm)*cvmgt(1.0d0,ar(imjm),ar(imjm).le.em6)
     &		/cvmgt(1.0d0,ar(imj),ar(imj).le.em6)
          endif
C          if (nf(ijm).ne.6) 
C     		v(ij)=v(ijm)*cvmgt(1.0d0,at(ijm),at(ijm).le.em6)
C     &		/cvmgt(1.0d0,at(ij),at(ij).le.em6)
C     		v(ij)=vymb(nmovbd(ij))+(v(ijm)-vymb(nmovbd(ij)))*at(ijm)
C     &		/at(ij)
		if (nf(ijm).eq.6) then
     			v(ij)=cvmgt(v(ij),v(ijm),abs(v(ij)).lt.abs(v(ijm)))
		else
			v(ij)=v(ijm)
		endif

  270     continue
          if (nf(imj).lt.6.or.(ar(imj).lt.em6)) goto 240
		if (nf(ipj)+nf(ipjp).eq.0.and.nf(ijp).ne.0.and.nf(imjp).eq.6.
     &		and.(at(ij).gt.1.0-em6).and.(at(ipj).gt.1.0-em6)) then
              if (nf(ijp)*nf(ij).eq.0) goto 2400
C...............zero tangential gradient; close to zero stress
     			v(ij)=v(ipj)*cvmgt(1.0d0,at(ipj),at(ipj).le.em6)
     &		/cvmgt(1.0d0,at(ij),at(ij).le.em6)
		endif
2400		if (nf(ipj)+nf(ipjm).eq.0.and.nf(ijm).ne.0.and.nf(imjm).eq.6.
     &		and.(at(ijm).gt.1.0-em6).and.(at(ipjm).gt.1.0-em6)) then 
          	if (nf(ijm)*nf(ij).eq.0) goto 240
C...............zero tangential gradient; close to zero stress
     			v(ijm)=v(ipjm)*cvmgt(1.0d0,at(ipjm),at(ipjm).le.em6)
     &		/cvmgt(1.0d0,at(ijm),at(ijm).le.em6)
		endif
C		if (nf(ipj).ne.6) 
C     		u(imj)=u(ij)*cvmgt(1.0d0,ar(ij),ar(ij).le.em6)
C     &		*r(i)/(cvmgt(1.0d0,ar(imj),ar(imj).le.em6)*r(i-1))
C		u(imj)=uxmb(nmovbd(ij))+(u(ij)-uxmb(nmovbd(ij)))*ar(ij)
C     &		/ar(imj)
		if (nf(ipj).eq.6) then
			u(imj)=cvmgt(u(ij),u(imj),abs(u(ij)).lt.abs(u(imj)))
		else
			u(imj)=u(ij)
		endif

  240     continue
          if(nf(ijm).lt.6.or.(at(ijm).lt.em6)) goto 250
		if (nf(ijp)+nf(ipjp).eq.0.and.nf(ipj).ne.0.and.nf(ipjm).eq.6.
     &		and.(ar(ij).gt.1.0-em6).and.(ar(ijp).gt.1.0-em6)) then
          	if (nf(ipj)*nf(ij).eq.0) goto 2500
C...............zero tangential gradient; close to zero stress
     			u(ij)=u(ijp)*cvmgt(1.0d0,ar(ijp),ar(ijp).le.em6)
     &		/cvmgt(1.0d0,ar(ij),ar(ij).le.em6)
          endif
2500		if (nf(ijp)+nf(imjp).eq.0.and.nf(imj).ne.0.and.nf(imjm).eq.6.
     &		and.(ar(imj).gt.1.0-em6).and.(ar(imjp).gt.1.0-em6)) then
          	if (nf(imj)*nf(ij).eq.0) goto 250
C...............zero tangential gradient; close to zero stress
     			u(imj)=u(imjp)*cvmgt(1.0d0,ar(imjp),ar(imjp).le.em6)
     &		/cvmgt(1.0d0,ar(imj),ar(imj).le.em6)
		endif
C		if (nf(ijp).ne.6) 
C     		v(ijm)=v(ij)*cvmgt(1.0d0,at(ij),at(ij).le.em6)
C     &		/cvmgt(1.0d0,at(ijm),at(ijm).le.em6)
C     		v(ijm)=vymb(nmovbd(ij))+(v(ij)-vymb(nmovbd(ij)))*at(ij)
C     &		/at(ijm)
     		if (nf(ijp).eq.6) then
			v(ijm)=cvmgt(v(ij),v(ijm),abs(v(ij)).lt.abs(v(ijm)))
		else
			v(ijm)=v(ij)
		endif

  250     continue
  420   continue
  421 continue

      do 5421 i=2,im1
        do 5420 j=2,jm1
c
	  ij=(j-1)*imax+i
          ijm=ij-imax
          imj=ij-1
          ipj=ij+1
          ijp=ij+imax
          imjm=ij-1-imax
          imjp=ij-1+imax
          ipjm=ij-imax+1
          ipjp=ij+imax+1

C.........no free surface BC imposed on obstacle cell, which could be in motion
C	    if(ac(ij)*ar(ij)*ar(imj)*at(ij)*at(ijm).lt.1.0d0-em6)goto 5420 
c
		if (ac(ij).le.em6) goto 5420
C.........added for impact purpose (it's found from the test of dam break
C.........that even with a thin layer of water on the right corner, this
C.........condition won't take effect; it plays the role only when the 
C.........corner is originally dry).
	    if (i.eq.2.and.nf(ij).eq.2.and.ibcfinal.eq.1.and.
     &		nfold(ij).ne.6) goto 380
          if (i.eq.im1.and.nf(ij).eq.1.and.ibcfinal.eq.1.and.
     &		nfold(ij).ne.6) goto 380
          if (j.eq.2.and.nf(ij).eq.4.and.ibcfinal.eq.1.and.
     &		nfold(ij).ne.6) goto 380
          if (j.eq.jm1.and.nf(ij).eq.3.and.ibcfinal.eq.1.and.
     &		nfold(ij).ne.6) goto 380
C.........skip the newly empty cell
     		if (nfold(ij).le.5.and.nf(ij).gt.5) goto 5420	
C.........for cell jumping from empty to fluid cell, div(u)<>0 in general and
C.........no boundary condition can help resolve this problem

c
C.........special attention should be made to take care the cell jump
C.........from empty to full due to breaking
          if ((nf(ij).eq.0.and.nfold(ij).ne.6).
     &		or.(nf(ij).gt.5.and.nfold(ij).gt.5)) go to 5420
          if (i.eq.2.and.kl.eq.5) go to 5420
          if (i.eq.im1.and.kr.eq.5) go to 5420
          if (j.eq.jm1.and.kt.eq.5) go to 5420
          if (j.eq.2.and.kb.eq.5) go to 5420

C.....block the original condition
C          dij=rdx(i)*(cvmgt(1.0d0,ar(ij),ar(ij).le.em6)
C     &	    *r(i)*u(ij)-cvmgt(1.0d0,ar(imj),ar(imj).le.em6)
C     &	    *r(i-1)*u(imj))+
C     &        rdy(j)*(cvmgt(1.0d0,at(ij),at(ij).le.em6)
C     &	    *ri(i)*v(ij)-cvmgt(1.0d0,at(ijm),at(ijm).le.em6)
C     &	    *ri(i)*v(ijm))
C.....Apply moving body
          dij=rdx(i)*(ar(ij)
     &	    *r(i)*(u(ij)-uxmb(nmovbd(ij)))-ar(imj)
     &	    *r(i-1)*(u(imj)-uxmb(nmovbd(ij))))+
     &        rdy(j)*(at(ij)
     &	    *ri(i)*(v(ij)-vymb(nmovbd(ij)))-at(ijm)
     &	    *ri(i)*(v(ijm)-vymb(nmovbd(ij))))

C.........If two neighbors are all open to air, redistribute dij based on ar and at
  310     if(nf(ipj).lt.6.or.(ar(ij).lt.em6)) go to 330
		if (nf(ijp).eq.6.and.at(ij).ge.em6) then
			u(ij)=u(ij)-delx(i)*dij*(ar(ij)/(ar(ij)+at(ij)))
     &			/(cvmgt(1.0d0,ar(ij),ar(ij).le.em6)*r(i))
			v(ij)=v(ij)-dely(j)*dij*(at(ij)/(ar(ij)+at(ij)))
     &			/(cvmgt(1.0d0,at(ij),at(ij).le.em6)*ri(i))
		else
			u(ij)=u(ij)-delx(i)*dij
     &			/(cvmgt(1.0d0,ar(ij),ar(ij).le.em6)*r(i))
		endif
          go to 380

  330     if(nf(ijp).lt.6.or.(at(ij).lt.em6)) go to 320
		if (nf(imj).eq.6.and.ar(imj).ge.em6) then
			v(ij)=v(ij)-dely(j)*dij*(at(ij)/(ar(imj)+at(ij)))
     &			/(cvmgt(1.0d0,at(ij),at(ij).le.em6)*ri(i))
			u(imj)=u(imj)+delx(i)*dij*(ar(imj)/(ar(imj)+at(ij)))
     &			/(cvmgt(1.0d0,ar(imj),ar(imj).le.em6)*r(i-1))
		else
			v(ij)=v(ij)-dely(j)*dij
     &			/(cvmgt(1.0d0,at(ij),at(ij).le.em6)*ri(i))
		endif
		go to 380

  320     if(nf(imj).lt.6.or.(ar(imj).lt.em6)) go to 340
		if (nf(ijm).eq.6.and.at(ijm).ge.em6) then
			u(imj)=u(imj)+delx(i)*dij*(ar(imj)/(ar(imj)+at(ijm)))
     &			/(cvmgt(1.0d0,ar(imj),ar(imj).le.em6)*r(i-1))
			v(ijm)=v(ijm)+dely(j)*dij*(at(ijm)/(ar(imj)+at(ijm)))
     &			/(cvmgt(1.0d0,at(ijm),at(ijm).le.em6)*ri(i))
		else
			u(imj)=u(imj)+delx(i)*dij
     &			/(cvmgt(1.0d0,ar(imj),ar(imj).le.em6)*r(i-1))
		endif
          go to 380

  340     if(nf(ijm).lt.6.or.(at(ijm).lt.em6)) go to 380
          v(ijm)=v(ijm)+dely(j)*dij/(cvmgt(1.0d0,at(ijm),
     &		at(ijm).le.em6)*ri(i))
          go to 380

c
c....     set velocities in empty cells
c         adjacent to partial fluid cells
c
  380     continue
C.........supply boundary condition
	  if (nfold(ij).eq.6) then
                if (i.eq.2) v(imj)=v(ij)
                if (i.eq.2) v(imjm)=v(ijm)
                if (i.eq.im1) v(ipj)=v(ij)
                if (i.eq.im1) v(ipjm)=v(ijm)
                if (j.eq.2) u(ijm)=u(ij)
                if (j.eq.2) u(imjm)=u(imj)
                if (j.eq.jm1) u(ijp)=u(ij)
                if (j.eq.jm1) u(ij+imax-1)=u(imj)
	  endif

 5420   continue
 5421 continue
c
 9999 continue
c
c.... Special velocity boundary conditions
c
      return
      end
