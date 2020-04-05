      subroutine k_epsilon

c##############################################################
      implicit real*8 (a-h,o-z)      
c##############################################################
c
c############
      include  "comdk2.h"
c############
c
	  xkappa=0.41d0
	  gamma=1.0d0
C.......it's found that for plane jet problem, instability occurs when
C.......gamma<1.0 (even when gamma=0.9)
C	  gamma=con
	  cmass=0.01
C       cmass=0.00
	  xmaxxnut=xnu*1.0e-6
	  prodmax=0.0
	  crest=1.0d16

C.......use the calculate stress to compute the production term
c------------------------------
	  call stress
c------------------------------

	  do 10 j=2,jmax-1
	  do 10 i=2,imax-1
	    ij=(j-1)*imax+i
          ijm=ij-imax
          imj=ij-1
          ipj=ij+1
          ijp=ij+imax
          imjm=ij-1-imax
          imjp=ij-1+imax
          ipjm=ij-imax+1
          ipjp=ij+imax+1

C.........for solid, porous media, or empty cell, set zero and jump to next cell
	    if (ac(ij).le.em6.or.npc(ij).ne.1.or.nf(ij).gt.5) then
              xk(ij)=0.0
              xep(ij)=0.0
              xnut(ij)=0.0
			goto 10
	    endif

C.........enforce the low turbulence level around the source region
          if (ninflow.eq.100) then
          	if (i.ge.isources-10.and.i.le.isourcee+10) then
                  xk(ij)=0.5*uturb**2
                  xnut(ij)=5.00*xnu*ratio**2
                  xep(ij)=c_mu*xk(ij)**2/xnut(ij)
				goto 10
          	endif 
          endif

C.........special way to compute the cell at wall boundary
C.........left and right wall boundary or nearly vertical obstacle
	    if ((i.eq.2.and.(kl.eq.1.or.kl.eq.2)).or.
     &		(i.eq.imax-1.and.(kr.eq.1.or.kr.eq.2)).or.
     &		(ac(ij).lt.1.0.and.max(at(ij),at(ijm)).lt.
     &		max(ar(ij),ar(imj))).or.
     &		(ac(ij).eq.1.0.and.(ac(ipj)+ac(imj)).lt.
     &		(ac(ijp)+ac(ijm)).
     &		and.(ac(ipj)+ac(imj)).lt.2.0)) then
		
		  Vt=(v(ij)*at(ij)+v(ijm)*at(ijm))/2.0
C...........rough wall
		  if (XE.le.0.0) then
			if (abs(Vt).lt.1.0e-6) then
				xnut(ij)=xnu
				goto 10
			endif
			fricv=max(abs(Vt)/(5.75*log10(min(0.005,delx(i))*0.5
     &			*ac(ij)**1/roughness)+8.5),
     &			sqrt(xnu/delx(i)*abs(Vt))*1.0e-6)	
		  else
			if (abs(Vt).lt.1.0e-6) then
				xnut(ij)=xnu
				goto 10
			endif
C.............iteration to find out the friction velocity using
C.............Newton-Raphson for smooth wall
	    	Xpp=0.05d0
			niter=0
			XAA=log(XE*abs(Vt)*min(0.005,delx(i))*0.5*ac(ij)**1/xnu)
5			Xn=Xpp+(xkappa-Xpp*(XAA+log(Xpp)))/(1.0+xkappa/Xpp)
			if (abs((Xn-Xpp)/Xpp).gt.1.0e-3.and.niter.le.10.
     &			and.Xn.gt.0.0) then
				Xpp=Xn
				niter=niter+1
				goto 5
			else
 			  fricv=max(Xn*abs(Vt),sqrt(xnu/delx(i)*abs(Vt))*1.0e-6)
			endif
		  endif
		
		  xk(ij)=fricv**2/sqrt(C_mu_bc)
C...........damping function works for low turbulence
C...........(Launder & Spalding, 1972; Hinze, 1972; Lemos, 1991)
		  xep(ij)=fricv**3/(xkappa*delx(i)*0.5*(1.0-
     &			exp(-fricv*delx(i)*0.5/26.0/xnu)))
		  xnut(ij)=C_mu_bc*xk(ij)**2/xep(ij)
            goto 10
          endif

C.........bottom and top wall boundary or nearly horizontal obstacle
          if ((j.eq.2.and.(kb.eq.1.or.kb.eq.2)).or.
     &        (j.eq.jmax-1.and.(kt.eq.1.or.kt.eq.2)).or.
     &		(ac(ij).lt.1.0.and.max(at(ij),at(ijm)).ge.
     &		max(ar(ij),ar(imj))).or.
     &		(ac(ij).eq.1.0.and.(ac(ipj)+ac(imj)).gt.
     &		(ac(ijp)+ac(ijm)).and.(ac(ijp)+
     &		ac(ijm)).lt.2.0)) then   

		  Ut=(u(ij)*ar(ij)+u(imj)*ar(imj))/2.0
C...........rough wall
		  if (XE.le.0.0) then
			if (abs(Ut).le.1.0e-6) then
                      xnut(ij)=xnu
					goto 10
			endif
			fricv=max(abs(Ut)/(5.75*log10(min(0.005,dely(j))*0.5
     &			*ac(ij)**1/roughness)+8.5),
     &			sqrt(xnu/dely(j)*abs(Ut))*1.0e-6)
		  else
			if (abs(Ut).le.1.0e-6) then
                      xnut(ij)=xnu
					goto 10
			endif
C.............iteration to find out the friction velocity using
C.............Newton-Raphson
              Xpp=0.05
              niter=0
              XAA=log(XE*abs(Ut)*min(0.005,dely(j))*0.5*ac(ij)**1/xnu)
6             Xn=Xpp+(xkappa-Xpp*(XAA+log(Xpp)))/(1.0+xkappa/Xpp)
              if (abs((Xn-Xpp)/Xpp).gt.1.0e-3.and.niter.le.10.and.
     &			Xn.gt.0.0) then
                  Xpp=Xn
                  niter=niter+1
				goto 6
              else
                fricv=max(Xn*abs(Ut),sqrt(xnu/dely(j)*abs(Ut))*1.0e-6)
              endif
		  endif
            xk(ij)=fricv**2/sqrt(C_mu_bc)
            xep(ij)=fricv**3/(xkappa*dely(j)*0.5*(1.0-
     &			exp(-fricv*dely(j)*0.5/26.0/xnu)))
            xnut(ij)=C_mu_bc*xk(ij)**2/xep(ij)
		  goto 10
          endif

C.........special way to compute the cell adjacent to porous media cell
C.........nearly vertical wall
          if (npc(ipj).ne.1) then
		  ijtmp=ipj
C...........note that the direction of normal velocity is important;
C...........positive when it goes out and vise verse.
		  Vnn=-u(ij)*ar(ij)
	    else
		  if (npc(imj).ne.1) then
			ijtmp=imj
			Vnn=u(imj)*ar(imj)
		  else
			goto 555
		  endif
	    endif
C.........iteration to find out the new friction velocity using
C.........Newton-Raphson
          niter=0
          Vt=(v(ij)*at(ij)+v(ijm)*at(ijm))/2.0
          Vte=Vt-(v(ijtmp)*at(ijtmp)+v(ijtmp-imax)*at(ijtmp-imax))/2.0
	    Vte=abs(Vte)
	    Vt=abs(Vt)
C.........tentative treatment to make Vte=Vt and Vnn=0.0 so that it
C.........treats porous bed as impermeable plane
          Xpp=0.05d0*Vte
C.........it's found that the normal flux has strong effect to the total
C.........stress; when negative, it will make tau_m >> tau_w if not
C.........restricted; tentatively skip the process if -tau_m > tau_w to
C.........avoid negative value in sqrt. Should seek more appropriate
C.........equation than in Ilegbusi (1989) J. Heat Transfer Vol. 32.
	    if (Xpp**2/4.0.lt.-cmass*Vnn*Vt/(1.0+Xpp/xkappa/Vt)) then
              fricv=sqrt(xnu/delx(i)/0.5/ac(ij)*abs(Vt))
			goto 550
	    endif
          if (abs(Vte).lt.1.0e-6) then
             	xnut(ij)=xnu
              goto 10
          endif
15        XAA=-Vte/(Xpp**2)-(2.0*Xpp-cmass*
     &		Vnn/xkappa/(1.0+Xpp/xkappa/Vt)
     &		**2)/xkappa/2.0/(Xpp**2+Vt*cmass*
     &		Vnn/(1.0+Xpp/xkappa/Vt))
	    ffx=Vte/Xpp-log(XEpor*delx(i)/2.0/xnu*
     &		sqrt(Xpp**2+cmass*Vnn*Vt/(1.0+Xpp/xkappa/Vt)))/xkappa
          Xn=Xpp-ffx/XAA
          if (abs((Xn-Xpp)/Xpp).gt.1.0e-3.and.niter.le.10.
     &          and.Xn.gt.0.0) then 
              Xpp=Xn
              niter=niter+1
              goto 15
          else
 		fricv=max(
     &		sqrt(max(0.0d0,Xn**2+cmass*Vnn*Vt/(1.0+Xn/xkappa*Vt))),
     &		sqrt(xnu/delx(i)/0.5*abs(Vt)))
          endif
550       xk(ij)=fricv**2/sqrt(C_mu_bc)
          xep(ij)=fricv**3/(xkappa*delx(i)*0.5*(1.0-
     &                  exp(-fricv*delx(i)*0.5/26.0/xnu)))
          xnut(ij)=C_mu_bc*xk(ij)**2/xep(ij)
          goto 10

555	    continue

C.........special way to compute the cell adjacent to porous media cell
C.........nearly horizontal wall
          if (npc(ijp).ne.1) then
                ijtmp=ijp
                Vnn=-v(ij)*at(ij)
          else   
                if (npc(ijm).ne.1) then
                        ijtmp=ijm
                        Vnn=v(ijm)*at(ijm)
                else
                        goto 565
                endif
          endif
C.........iteration to find out the new friction velocity using
C.........Newton-Raphson
          niter=0
          Vt=(u(ij)*ar(ij)+u(imj)*ar(imj))/2.0
          Vte=Vt-(u(ijtmp)*ar(ijtmp)+u(ijtmp-1)*ar(ijtmp-1))/2.0
	    Vte=abs(Vte)
          Vt=abs(Vt)
          Xpp=0.05d0*Vte
          if (Xpp**2/2.0.lt.-cmass*Vnn*Vt/(1.0+Xpp/xkappa/Vt)) then
                fricv=sqrt(xnu/dely(j)/0.5*abs(Vt))
                goto 560
          endif
          if (abs(Vte).lt.1.0e-6) then
                xnut(ij)=xnu
                goto 10
          endif
16        XAA=-Vte/(Xpp**2)-(2.0*Xpp-cmass*
     &		Vnn/xkappa/(1.0+Xpp/xkappa/Vt)
     &		**2)/xkappa/2.0/(Xpp**2+Vt*cmass*
     &		Vnn/(1.0+Xpp/xkappa/Vt))
	    ffx=Vte/Xpp-log(XEpor*dely(j)/2.0/xnu*
     &		sqrt(Xpp**2+cmass*Vnn*Vt/(1.0+Xpp/xkappa/Vt)))/xkappa
          Xn=Xpp-ffx/XAA
          if (abs((Xn-Xpp)/Xpp).gt.1.0e-3.and.niter.le.10.
     &          and.Xn.gt.0.0) then 
			Xpp=Xn
			niter=niter+1
              goto 16
          else 
              fricv=max(
     &		sqrt(max(0.0d0,Xn**2+cmass*Vnn*Vt/(1.0+Xn/xkappa*Vt))),
     &		sqrt(xnu/dely(j)/0.5*abs(Vt)))
          endif
560       xk(ij)=fricv**2/sqrt(C_mu_bc)
          xep(ij)=fricv**3/(xkappa*dely(j)*0.5*(1.0-
     &                  exp(-fricv*dely(j)*0.5/26.0/xnu)))
          xnut(ij)=C_mu_bc*xk(ij)**2/xep(ij)
          goto 10

565	    continue

C.........start the main computation for cell not adjacent to any type of 
C.........boundaries (solid, porous, and computational boundaries)
		frt=(f(ij)+f(ipj)+f(ijp)+f(ipjp))/4.0
		frb=(f(ij)+f(ipj)+f(ijm)+f(ipjm))/4.0
		flt=(f(ij)+f(imj)+f(ijp)+f(imjp))/4.0
		flb=(f(ij)+f(imj)+f(ijm)+f(imjm))/4.0

	    dxmin=(delx(i)+delx(i-1))/2.0
	    dxplus=(delx(i)+delx(i+1))/2.0
	    dymin=(dely(j)+dely(j-1))/2.0
	    dyplus=(dely(j)+dely(j+1))/2.0
	    uct=(u(ij)*ar(ij)+u(imj)*ar(imj))/2.0
	    vct=(v(ij)*at(ij)+v(ijm)*at(ijm))/2.0
          delxa=dxmin+dxplus+gamma*sign(1.0d0,uct)*(dxplus-dxmin)
          delya=dymin+dyplus+gamma*sign(1.0d0,vct)*(dyplus-dymin)

C.........calculate effective eddy viscosity
	    xnre=(xnutn(ij)*delx(i+1)+xnutn(ipj)*delx(i))/
     &		 (delx(i)+delx(i+1))/sege+xnu
          if (((kr.eq.7.or.kr.eq.3).and.i.eq.im1).or.nf(ipj).ge.6)
     &		xnre=xnutn(ij)/sege+xnu
	    xnle=(xnutn(ij)*delx(i-1)+xnutn(imj)*delx(i))/
     &		 (delx(i)+delx(i-1))/sege+xnu
	    if (((kl.eq.7.or.kl.eq.3).and.i.eq.2).or.nf(imj).ge.6) 
     &		xnle=xnutn(ij)/sege+xnu
	    xnte=(xnutn(ij)*dely(j+1)+xnutn(ijp)*dely(j))/
     &		 (dely(j)+dely(j+1))/sege+xnu
          if (((kt.eq.7.or.kt.eq.3).and.j.eq.jm1).or.nf(ijp).ge.6) 
     &		xnte=xnutn(ij)/sege+xnu
	    xnbe=(xnutn(ij)*dely(j-1)+xnutn(ijm)*dely(j))/
     &		 (dely(j)+dely(j-1))/sege+xnu
          if (((kb.eq.7.or.kb.eq.3).and.j.eq.2).or.nf(ijm).ge.6) 
     &		xnbe=xnutn(ij)/sege+xnu
          xnrk=(xnutn(ij)*delx(i+1)+xnutn(ipj)*delx(i))/
     &           (delx(i)+delx(i+1))/segk+xnu
          if (((kr.eq.7.or.kr.eq.3).and.i.eq.im1).or.nf(ipj).ge.6) 
     &		xnrk=xnutn(ij)/segk+xnu
          xnlk=(xnutn(ij)*delx(i-1)+xnutn(imj)*delx(i))/
     &           (delx(i)+delx(i-1))/segk+xnu
          if (((kl.eq.7.or.kl.eq.3).and.i.eq.2).or.nf(imj).ge.6) 
     &		xnlk=xnutn(ij)/segk+xnu
          xntk=(xnutn(ij)*dely(j+1)+xnutn(ijp)*dely(j))/
     &           (dely(j)+dely(j+1))/segk+xnu
          if (((kt.eq.7.or.kt.eq.3).and.j.eq.jm1).or.nf(ijp).ge.6) 
     &		xntk=xnutn(ij)/segk+xnu
          xnbk=(xnutn(ij)*dely(j-1)+xnutn(ijm)*dely(j))/
     &           (dely(j)+dely(j-1))/segk+xnu
          if (((kb.eq.7.or.kb.eq.3).and.j.eq.2).or.nf(ijm).ge.6) 
     &		xnbk=xnutn(ij)/segk+xnu

C.........calculate gradient of k and epsilon
          dkdl=(xkn(ij)-xkn(imj))/dxmin
          dkdr=(xkn(ipj)-xkn(ij))/dxplus
          dkdb=(xkn(ij)-xkn(ijm))/dymin
          dkdt=(xkn(ijp)-xkn(ij))/dyplus
	    dedl=(xepn(ij)-xepn(imj))/dxmin
	    dedr=(xepn(ipj)-xepn(ij))/dxplus
	    dedb=(xepn(ij)-xepn(ijm))/dymin
	    dedt=(xepn(ijp)-xepn(ij))/dyplus

C.........impose zero gradient on free surface
C.........By so doing, for diffusion, it represent mirror condition on free surface
C.........For advection, since the upwind scheme is used, it means no air effect will
C.........by taken if flow is going towards the water (back of wave crest) 
C.........and water effects will be fully taken if flow is leaving surface (front crest)
C.........It is found that turbulence level can be very high near the free surface 
C.........after running many waves; We need further check whether this is because
C.........of the turbulence model or improper boundary condition applied;
C.........Numerical experiment has been taken by setting dk/dn(surface)=dk/dn(immediate
C.........interior cell but instability developed by having very small or negative k & epsilon
	    if (fn(ipj).eq.0.0) then
			dkdr=0.0
			dedr=0.0
	    endif
          if (fn(imj).eq.0.0) then
              dkdl=0.0
              dedl=0.0
          endif
          if (fn(ijp).eq.0.0) then
              dkdt=0.0
              dedt=0.0
          endif
          if (fn(ijm).eq.0.0) then
              dkdb=0.0
              dedb=0.0
          endif

C.........Update epsilon equation
C.........Compute the advective term
	    fex=uct/delxa*(dxplus*dedl+dxmin*dedr+gamma*sign(1.0d0,uct)*
     &			(dxplus*dedl-dxmin*dedr))
	    fey=vct/delya*(dyplus*dedb+dymin*dedt+gamma*sign(1.0d0,vct)*
     &			(dyplus*dedb-dymin*dedt))

C.........Compute viscous term
	    vise=(dedr*xnre-dedl*xnle)/delx(i)+(dedt*xnte-dedb*xnbe)
     &		 /dely(j)
C.........term to correct the negative diffusion by FTCS scheme
     &		 +0.5*(1.0-gamma)*delt*uct**2/delx(i)*((xepn(ipj)
     &		 -xepn(ij))/dxplus-(xepn(ij)-xepn(imj))/dxmin)
     &           +0.5*(1.0-gamma)*delt*vct**2/dely(j)*((xepn(ijp)
     &           -xepn(ij))/dyplus-(xepn(ij)-xepn(ijm))/dymin)
C.........see Lemos p96-97

C.........Compute the production term
	    dudx=(u(ij)*ar(ij)-u(imj)*ar(imj))/delx(i)
	    dvdy=(v(ij)*at(ij)-v(ijm)*at(ijm))/dely(j)
	    uyvxrt=(u(ijp)*ar(ijp)-u(ij)*ar(ij))
     &		/((dely(j)+dely(j+1))/2.0)+
     &		(v(ipj)*at(ipj)-v(ij)*at(ij))
     &		/((delx(i)+delx(i+1))/2.0)
          uyvxrb=(u(ij)*ar(ij)-u(ijm)*ar(ijm))
     &          /((dely(j)+dely(j-1))/2.0)+
     &          (v(ipjm)*at(ipjm)-v(ijm)*at(ijm))
     &          /((delx(i)+delx(i+1))/2.0)
          uyvxlt=(u(ijp-1)*ar(ijp-1)-u(imj)*ar(imj))
     &          /((dely(j)+dely(j+1))/2.0)+
     &          (v(ij)*at(ij)-v(imj)*at(imj))
     &          /((delx(i)+delx(i-1))/2.0)
          uyvxlb=(u(imj)*ar(imj)-u(imjm)*ar(imjm))
     &          /((dely(j)+dely(j-1))/2.0)+
     &          (v(ijm)*at(ijm)-v(imjm)*at(imjm))
     &          /((delx(i)+delx(i-1))/2.0)

C.........special effort to make shear contribution=0 on free surface
	    if (nf(ij).ne.0) then
C...........stability problem occurs when violent free surface overturning
C...........happends (wave plunging) and Youngs model is used. It's found
C...........instability always starts from free surface where the
C...........production terms are unrealistically calculated. Tentatively
C...........restrict prodnew=0.0 on free surface for this case.
              if (ninflow.eq.7.or.ninflow.eq.8.or.ninflow.eq.10) then
                        dudx=0.0
                        dvdy=0.0
              endif
			uyvxrt=0.0
			uyvxrb=0.0
			uyvxlt=0.0
			uyvxlb=0.0
	    endif

C.........It's likely to have unrealistically large shear on free
C.........surface when impact or plunging occurs; force it to be zero
C.........when adjacent to free surface
          if (nf(ijp).ne.0.or.nf(ipj).ne.0.or.nf(ipjp).
     &          ne.0) uyvxrt=0.0
          if (nf(ijm).ne.0.or.nf(ipj).ne.0.or.nf(ipjm).
     &          ne.0) uyvxrb=0.0
          if (nf(ijp).ne.0.or.nf(imj).ne.0.or.nf(imjp).
     &          ne.0) uyvxlt=0.0
          if (nf(ijm).ne.0.or.nf(imj).ne.0.or.nf(imjm).
     &          ne.0) uyvxlb=0.0

C.........calculating production terms
	    prodnew=(tauxxturb(ij)*dudx
     &	  +tauyyturb(ij)*dvdy
     &	  +(tauxyturb(ij)*uyvxrt+tauxyturb(imj)*uyvxlt
     &	  +tauxyturb(ijm)*uyvxrb+tauxyturb(imjm)*uyvxlb)
     &	  /4.0d0)/rhof*f(ij)**2

C.........do not account for negative production near free surface
	    if (prodnew.lt.0.0.and.(nf(ipj).ne.0.or.nf(imj).ne.0.
     &		or.nf(ijp).ne.0.or.nf(ijm).ne.0)) prodnew=0.0

C.........do not account for negative production for cell with very small
C.........turbulence where molecular viscous effect is important and
C.........turbulence calculation is inaccurate
          if (prodnew.lt.0.0.and.xk(ij).lt.0.5*(uturb*5.0)**2)
     &  	prodnew=0.0

	    if (prodnew.lt.0.0) then
C.........for k equation, it is linear decay; for epsilon equation,
C.........it is exponential decay! calculated differently
              crest=min(crest,xk(ij)/(-prodnew+xep(ij))/2.0,
     &		xk(ij)/(-c1e*prodnew+c2e*xep(ij)))
	    endif

C.........calculate time scale of turbulence over time scale of mean strain
C.........only used when RNG is used
	    ttots(ij)=xk(ij)/xep(ij)*
     &		dsqrt(2.0*dudx**2+2.0*dvdy**2+
     &		((uyvxrt+uyvxrb+uyvxlt+uyvxlb)/4.0)**2)
C.........calculate total mean strain (only used for check)
          tss(ij)=dsqrt(2.0*dudx**2+2.0*dvdy**2+
     &          ((uyvxrt+uyvxrb+uyvxlt+uyvxlb)/4.0)**2)

C.........check the abnormal situation for xkn(ij)
	    txep=xep(ij)
	    if (xkn(ij).lt.1.0e-16) then
			write(9,*)'strange',t,i,j,f(ij),ac(ij),nf(ij),
     &		nf(ipj),nf(imj),nf(ijp),nf(ijm),
     &		nfold(ij),nfold(ipj),nfold(imj),nfold(ijp),
     &		nfold(ijm),f(ipj),f(imj),f(ijp),f(ijm),
     &		xkn(ij),xkn(ipj),xkn(imj),xkn(ijp),xkn(ijm)
			xep(ij)=txep
			goto 233
	    endif 

C.........compute xep(ij)
          c1ee=c1e+c1ea-c1ea/(1.0+c1eb*ttots(ij))
		c2ee=c2e
C.........when RNG theory is used, c2e has chance to go negative when
C.........S*k/epsilon is large; Adapt explicit scheme to avoid
C.........numerical instability when this (rarely) happens.
          if (c2ee.ge.0.0) then
              xep(ij)=(xepn(ij)/delt-fex-fey+vise+c1ee*xepn(ij)/xkn(ij)
     &          *prodnew)/(1.0/delt+c2ee*xepn(ij)/xkn(ij))
          else
              xep(ij)=xepn(ij)+delt*(-fex-fey+vise+c1ee*xepn(ij)/xkn(ij)
     &          *prodnew-c2ee*xepn(ij)**2/xkn(ij))
          endif

C.........record abnormal situation
	    if (xep(ij).lt.1.0e-16) then
			write(9,*)'eps too small',t,i,j,xep(ij),xepn(ij),vise,
     &		xkn(ij),dedt,xnte,dedb,xnbe,dedr,xnre,dedl,xnle,
     &		prodnew,u(ij),u(imj),v(ij),
     &		v(ijm),f(ij),f(ipj),f(imj),f(ijp),f(ijm)
			xep(ij)=txep
	    endif
	    if (xep(ij).gt.1.0e16) then 
			write(9,*)'unstable; eps too large',t,i,j,xep(ij)
			xep(ij)=txep
	    endif

233	    continue

C.........update K equation
	    fkx=uct/delxa*(dxplus*dkdl+dxmin*dkdr+gamma*sign(1.0d0,uct)*
     &			(dxplus*dkdl-dxmin*dkdr))
	    fky=vct/delya*(dyplus*dkdb+dymin*dkdt+gamma*sign(1.0d0,vct)*
     &			(dyplus*dkdb-dymin*dkdt))

          visx=(dkdr*xnrk-dkdl*xnlk)/delx(i)
     &           +0.5*(1.0-gamma)*delt*uct**2/delx(i)
     &           *((xkn(ipj)-xkn(ij))/dxplus-(xkn(ij)-xkn(imj))/dxmin)
          visy=(dkdt*xntk-dkdb*xnbk)/dely(j)
     &           +0.5*(1.0-gamma)*delt*vct**2/dely(j)
     &           *((xkn(ijp)-xkn(ij))/dyplus
     &           -(xkn(ij)-xkn(ijm))/dymin)
          visk=visx+visy

C.........record earlier number
	    txk=xk(ij)
C.........for output only
		xprod(ij)=prodnew
	    
		xk(ij)=xkn(ij)+delt*(-fkx-fky+visk+prodnew-xep(ij))

C.........if RNG is used, eps can grow very rapicdly due to negative c2ee
C.........use semi-implicit method to calculate k to ensure stability
	    if (xk(ij).le.0.0.and.c2ee.le.0.0) then
	      xk(ij)=(xkn(ij)/delt
     &	      -fkx-fky+visk+prodnew)/(1.0/delt+c_mu*xkn(ij)/xnutn(ij))
	      write(9,*)'rng, negative k using explicit method'
	    endif

C.........record abnormal situation
	    if (xk(ij).lt.1.0e-16) then
              write(9,*)'k too small',t,i,j,xk(ij),xkn(ij),delt,
     &        visk,fkx,fky,xep(ij),xepn(ij),prodnew,
     &		xkn(ipj),xkn(imj),xkn(ijp),xkn(ijm),
     &		dkdr,xnrk,dkdl,xnlk,dkdt,xntk,dkdb,xnbk,
     &		u(ij),u(imj),v(ij),v(ijm),f(ij),
     &		f(ipj),f(imj),f(ijp),f(ijm)
	 		xk(ij)=txk
	    endif
	    if (xk(ij).le.1.0e-16.or.xep(ij).le.1.0e-16) then
			write(9,*)'***check for scheme'
			xk(ij)=0.5*uturb**2
              xnut(ij)=5.00*xnu*ratio**2
              xep(ij)=C_mu*xk(ij)**2/xnut(ij)
          endif

C.........update eddy viscosity
	    xnut(ij)=C_mu*xk(ij)**2/xep(ij)
C.........enforcing damping function so that the diffusion of k, eps
C.........will be calculated in consistant with stress tau calculation
	    xnut(ij)=xnut(ij)*
     &                  xnut(ij)/(xnu*eddycoef)*
     &                  (1.0-exp(-eddycoef*xnu/max(1.0d-20,xnut(ij))))

C.........record maximum value of eddy viscosity for stability condition
	    if (nf(ij).eq.0) then
                if (xnut(ij).ge.xmaxxnut) then
                       xmaxxnut=xnut(ij)
                       xmaxxk=xk(ij)
                       xmaxxep=xep(ij)
					 prodmax=prodnew
                       ni=i
                       nj=j
					 xx=xi(i)
					 yy=yj(j)
                endif
	    endif
C.........record total dissipation for output
          xdis=xdis+delt*xep(ij)
     &      *(x(i)-x(i-1))*(y(j)-y(j-1))*f(ij)*cvmgt(porousp(ij),
     &	    ac(ij),npc(ij).ne.1.and.npor.ne.0)

10	  continue


CCC	    special boundary on wall for inflow and outflow
        il=2
        ir=im1
        do 65 j=1,jmax
          	ijl=(j-1)*imax+il
          	ijr=(j-1)*imax+ir
          	if (kl.ne.6) then
	    	  xk(ijl-1)=xk(ijl)
	    	  xep(ijl-1)=xep(ijl)
	    	  xnut(ijl-1)=xnut(ijl)
	  		endif
	  		xk(ijr+1)=xk(ijr)
	  		xep(ijr+1)=xep(ijr)
          	xnut(ijr+1)=xnut(ijr)
65      continue
	
	  jb=2
	  jt=jm1
	  do 55 i=1,imax
	  		ijb=(jb-1)*imax+i
	  		ijt=(jt-1)*imax+i
	  		xk(ijb-imax)=xk(ijb)
	  		xep(ijb-imax)=xep(ijb)
	  		xnut(ijb-imax)=xnut(ijb)
	  		xk(ijt+imax)=xk(ijt)
	  		xep(ijt+imax)=xep(ijt)
	  		xnut(ijt+imax)=xnut(ijt)
55	  continue

        write(9,'(3i5,3f8.3,4e12.4)')ncyc,ni,nj,t,xx,yy,xmaxxk,
     &	    xmaxxep,prodmax,xmaxxnut

	  return
	  end

	
