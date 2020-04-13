      program MONGOOSE_MULTIPHASE
c
c ======================================================================
c This version combines
c 1) the original version of mongoose (ripple + multi-scale les 
c    + porous media, Shi et al., 2004)
c 2) k-e turbulence closure (Lin and Liu, 1998)
c 3) bubble evolution using prescribed bubble entrainment
c    (Shi et al., 2008)
c 4) sand transport (not completed)
c ======================================================================
c       
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c#######################################################################
      implicit real*8 (a-h,o-z)
      include  "comdk2.h" 
      include  "bubble1.h"
c############
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... establish diagnostics
c
c     ------------------
      call begin
c     ------------------
c
c.... read input and do the problem setup
c
c     ------------------
      call setup
c     ------------------
c
c
c.... print initial input data
c
c     ------------------
      call prtplt2 (1)
c     ------------------
c
c.... set initial boundary conditions
c
c     ------------------
      call bc
c     ------------------
c
c.... end of problem generation,
c     write out misc info
c
      tset=tsetf-tseti
      write (13,100) tset,t,ncyc
      write (iotty,100) tset,t,ncyc
c
      if (nrs.eq.1) then
          read(10)t,delt,twprt_tmp,twplt,flgc,ncyc,
     &          iter,nflgc,nanimation,twprt_a,xmass0,xnewloc,ynewloc,
     &		  oa1,ob1,oc1,vchgt
          do 1111 j=1,jmax
		  read(10)(u((j-1)*imax+i),i=1,imax)
		  read(10)(v((j-1)*imax+i),i=1,imax)
		  read(10)(p((j-1)*imax+i),i=1,imax)
		  read(10)(f((j-1)*imax+i),i=1,imax)
		  read(10)(nf((j-1)*imax+i),i=1,imax)
		  read(10)(xk((j-1)*imax+i),i=1,imax)
		  read(10)(xep((j-1)*imax+i),i=1,imax)
		  read(10)(xnut((j-1)*imax+i),i=1,imax)
            read(10)(ac((j-1)*imax+i),i=1,imax)
            read(10)(at((j-1)*imax+i),i=1,imax)
            read(10)(ar((j-1)*imax+i),i=1,imax)
		  read(10)(nmovbd((j-1)*imax+i),i=1,imax)
1111      continue
	    read(10)(ep0(i),i=1,imax),uxmb,vymb
      endif

      twprt=tstart
	twplt=twprt
      twprt_a=tstart_a

c initial uvz
        if(ninflow.eq.101) then
         call initial_uvz
        endif
c
c.... Begin cycle
c
c     ------------
   10 call newcyc
c     ------------
c
      if (nexit.eq.9999.or.nexit.eq.9998) then
C.....have a dump file ready for check when the program is terminated normally or abnormally
          write(8)t,delt,twprt,twplt,flgc,ncyc,
     &          iter,nflgc,nanimation,twprt_a,xmass0,xnewloc,ynewloc,
     &			oa1,ob1,oc1,vchgt
          do j=1,jmax
                write(8)(u((j-1)*imax+i),i=1,imax)
                write(8)(v((j-1)*imax+i),i=1,imax)
                write(8)(p((j-1)*imax+i),i=1,imax)
                write(8)(f((j-1)*imax+i),i=1,imax)
                write(8)(nf((j-1)*imax+i),i=1,imax)
                write(8)(xk((j-1)*imax+i),i=1,imax)
                write(8)(xep((j-1)*imax+i),i=1,imax)
                write(8)(xnut((j-1)*imax+i),i=1,imax)
                write(8)(ac((j-1)*imax+i),i=1,imax)
                write(8)(at((j-1)*imax+i),i=1,imax)
                write(8)(ar((j-1)*imax+i),i=1,imax)
			  write(8)(nmovbd((j-1)*imax+i),i=1,imax)
          end do
          write(8) (ep0(i),i=1,imax),uxmb,vymb
          goto 777
      endif
	iter=0
c
c.... Compute the tilde velocities
c
c     -------------
      call vtilde
c     -------------
c
c
c.... Directly solve the pressure Poisson equation
c
c     -------------
      call implctp
c     -------------
c
c

c.... turbulent model
c
c     --------------
      if (kemodel.gt.0) call k_epsilon
c	--------------
c
c
c.... re-compute without having to touch vofadv if courant condition is violated
c
c     --------------
	if (flgc.eq.2.0) goto 10
c	--------------
c
c.....skip VOFADV if KR=5 (applied pressure for no free surface case)
c
	if (kr.eq.5) goto 10
c
c.... Advect the VOF function
c
c     -------------
      if (nfree.eq.1) then
        call vofadv
      else
        if (nfree.eq.4) then
		call Youngs
        endif
      endif
c     -------------
c
c
      go to 10
c
  550 format (110f8.3)
  100 format (" End of problem generation at ",1pe12.5," s ",/,
     &        " Entering main loop at time ",1pe12.5," , ",
     &        " cycle ",i5)

777   continue

	close(9)
	close(11)
	close(12)
	close(13)

      end

      subroutine initial_uvz
c##############################################################
      implicit real*8 (a-h,o-z)
      
c##############################################################
c
c############
      include  "comdk2.h"
      include  "bubble1.h" 
      real eta_1d(imax),u_1d(imax)
c############

      print*,'read u v eta, specified in main.f subroutine'
      print*, 'initial_uvz'
      print*, 'imax,jmax = ', imax, jmax
      
      open(1,file='eta_ripple.txt')
       do i=1,imax
         read(1,*) eta_1d(i)
       enddo
       close(1)
      open(1,file='u_ripple.txt')
       do i=1,imax
         read(1,*) u_1d(i)
       enddo
       close(1)

      do i=1,imax
       do 20 j=1,jmax

         eta10=eta_1d(i)
         ij1=(j-1)*imax+i

         f(ij1)=(eta10-y(j))/dely(j)
         if(f(ij1).le.em6) f(ij1)=0.0
         if(f(ij1).ge.(1.0-em6)) f(ij1)=1.0

         p(ij1)=(eta10-yj(j))*9.8*rhof
         
         u(ij1)=u_1d(i)
         v(ij1)=0.0

20    continue

      enddo

      return
      end

      subroutine newcyc
c
c ======================================================================
c
c   Purpose -
c     begin cycle:  provide monitor print, optional graphics
c                   and/or field variable print, and test for run
c                   termination.  if continuing, increment time
c                   and cycle, move advance time arrays into
c                   time-n arrays, and adjust time step.
c
c   NEWCYC is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE
c
c
c   NEWCYC calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c          PRTPLT    ripple    SECOND    system      KILL    ripple
c          TAPOUT    ripple   DELTADJ    ripple     EXITA    system
c         PLOTOUT    ripple
c
c ======================================================================
c
c##############################################################
      implicit real*8 (a-h,o-z)
      
c##############################################################
c
c############
      include  "comdk2.h"
      include  "bubble1.h" 
c############
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      dimension uc(nxy),vc(nxy)
c
c.... save the Lagrangian velocities
c
      do 5 ij=1,nxy
        uc(ij)=un(ij)
        vc(ij)=vn(ij)
    5 continue

c.... if vof function convection limit exceeded,
c     decrease the timestep and try again
c
C      if (flgc.gt.0.5.and.autot.ne.0.0) go to 100

      if (nloc.ne.0) then
        if (t.gt.xxxf.and.t.lt.ttend) then
          do 88 n=1,nloc
                write(700+5*(n-1)+1,'(100e14.6)')t,etah(nout(n)),
     &          (f((j-1)*imax+nout(n)),j=2,jmax-1),
     &          ((u((j-1)*imax+nout(n))*ar((j-1)*imax+nout(n))
     &		   +u((j-1)*imax+nout(n)-1)*ar((j-1)*imax+nout(n)-1))/2,
     &			j=2,jmax-1),
     &          ((v((j-1)*imax+nout(n))*at((j-1)*imax+nout(n))
     &			+v((j-2)*imax+nout(n))*at((j-2)*imax+nout(n)))/2,
     &			j=2,jmax-1),
     &          (p((j-1)*imax+nout(n)),j=2,jmax-1),
     &          (sqrt(2.0*xk((j-1)*imax+nout(n))),j=2,jmax-1),
     &          (xnut((j-1)*imax+nout(n))*1000.0,j=2,jmax-1),
     &          ((tauxy((j-1)*imax+nout(n))+tauxy((j-1)*imax+nout(n)-1)
     &          +tauxy((j-2)*imax+nout(n))+tauxy((j-2)*imax+nout(n)-1))
     &			/4,j=2,jmax-1)
88        continue
          xxxf=xxxf+prtdt_t
          endif
        endif
c
c.... print time and cycle data on paper
c
c     ----------------
      call prtplt2 (2)
c     ----------------
c
        xmass=0.0
      do 181 i=2,imax-1
        do 181 j=2,jmax-1
	  if (kr.eq.4.and.i.eq.imax-1) goto 181
          ij=(j-1)*imax+i
          imj=ij-1
          ijm=ij-imax
          xmass=xmass+f(ij)*cvmgt(porousp(ij),ac(ij),npc(ij).
     &	ne.1.and.npor.ne.0)*(x(i)-x(i-1))*(y(j)-y(j-1))
  181 continue
C.......for sakakiyama's case only
C.......the whole wave tank is 20m and the computational domain is 7m long;
C.......thus the overtopping water M will reduce the mean d by M/(20/7)
	if (nweakref.eq.10) xmass=xmass+xsaka/3.0
        if (t.eq.0.0) xmass0=xmass
        xmass=xmass-xmass0

C.......analytical solution for total mass
        if (nweakref.gt.1.and.ninflow.eq.34) then
                xmass_a=aa*xxl/4.0d0/pi*(1.0-cos(2.0d0*pi*t/xxt))
C.......consider to compensate the reflected wave (for vertical wall)
C.......noted when this alternative is used, consistent change should be
C.......made for free surface
		if (nweakref.eq.10) then
C.......constant water depth
C		  distotal=(x(im1)-x(1))*2.0d0
C		  refcoef=1.0d0
C		  tref=distotal/c1
C.......wave on slope; specify the tref and refcorf manually from experiment.
C.......sakakiyama's case; assume the reflected wave has a=1.3cm from exp
		  tref=10.09
		  refcoef=1.3/10.5
		  if (t.ge.tref) then
		    xmass_a=xmass_a-refcoef*
     &		    aa*xxl/4.0d0/pi*(1.0-cos(2.0d0*pi*(t-tref)/xxt))
		  endif
		endif

		if (ncyc.eq.0) then
	                xmass_a1=xmass_a1+eta1*delt*c1
		else
	                xmass_a1=xmass_a1+eta1a*delt*c1
		endif
                xmass_d=xmass-xmass_a
        end if

        if (nweakref.gt.1.and.ninflow.eq.24) then
                if (ncyc.eq.0) then
                        xmass_a1=xmass_a1+eta1*delt*c1
                else
                        xmass_a1=xmass_a1+eta1a*delt*c1
                endif
                xmass_d=xmass-xmass_a1
        end if

        dissipt=0.0
        dissipm=0.0
	  productt=0.0
        do 271 j=2,jmax-1
        do 271 i=2,imax-1
          if (kr.eq.4.and.i.eq.imax-1) goto 271
          ij=(j-1)*imax+i
          imj=ij-1
          ipj=ij+1
          ijm=ij-imax
          ijp=ij+imax
          imjp=ijp-1
          ipjm=ijm+1
          imjm=ijm-1
          disturb=xep(ij)*f(ij)
          dissipturb(ij)=dissipturb(ij)+delt*disturb
          dissipt=dissipt+disturb
     &          *delx(i)*dely(j)*cvmgt(porousp(ij),
     &          ac(ij),npc(ij).ne.1.and.npor.ne.0)
		prodturb=xprod(ij)*f(ij)
		productturb(ij)=productturb(ij)+delt*prodturb
		productt=productt+prodturb
     &		*delx(i)*dely(j)*cvmgt(porousp(ij),
     &		ac(ij),npc(ij).ne.1.and.npor.ne.0)
          uupr=(u(ij)*dely(j+1)+u(ijp)*dely(j))/(dely(j)+dely(j+1))
          uupl=(u(imj)*dely(j+1)+u(imjp)*dely(j))/(dely(j)+dely(j+1))
          uup=(uupr+uupl)/2
          ubtr=(u(ij)*dely(j-1)+u(ijm)*dely(j))/(dely(j)+dely(j-1))
          ubtl=(u(imj)*dely(j-1)+u(imjm)*dely(j))/(dely(j)+dely(j-1))
          ubt=(ubtr+ubtl)/2
          vrtu=(v(ij)*delx(i+1)+v(ipj)*delx(i))/(delx(i)+delx(i+1))
          vrtb=(v(ijm)*delx(i+1)+v(ipjm)*delx(i))/(delx(i)+delx(i+1))
          vrt=(vrtu+vrtb)/2
          vlfu=(v(ij)*delx(i-1)+v(imj)*delx(i))/(delx(i)+delx(i-1))
          vlfb=(v(ijm)*delx(i-1)+v(imjm)*delx(i))/(delx(i)+delx(i-1))
          vlf=(vlfu+vlfb)/2
          gradsquare=((u(ij)-u(imj))/delx(i))**2
     &          +((v(ij)-v(ijm))/dely(j))**2
     &          +((uup-ubt)/dely(j))**2+((vrt-vlf)/delx(i))**2
          dismole=xnu*gradsquare*f(ij)
          dissipmole(ij)=dissipmole(ij)+delt*dismole
          dissipm=dissipm+dismole
     &          *delx(i)*dely(j)*cvmgt(porousp(ij),
     &          ac(ij),npc(ij).ne.1.and.npor.ne.0)
271	continue

        if (t+em6.lt.twplt.and.t.lt.twfin) go to 60
        twplt=twplt+pltdt
c
C========================================================================
           
        energyk=0.0
        energyp=0.0
	  energyt=0.0
	  xpollutant=0.0
	  xforce=0.0
	  yforce=0.0
 
        do 171 j=2,jmax-1
        do 171 i=2,imax-1
	  if (kr.eq.4.and.i.eq.imax-1) goto 171
          ij=(j-1)*imax+i
          imj=ij-1
          ijm=ij-imax
          energyk=energyk+(((u(ij)+u(imj))/2.)**2/2.0+
     &  	((v(ij)+v(ijm))/2.)**2/2.0)
     &  	*delx(i)*dely(j)*f(ij)*cvmgt(porousp(ij),
     &		ac(ij),npc(ij).ne.1.and.npor.ne.0)
          energyt=energyt+xk(ij)
     &  	*delx(i)*dely(j)*f(ij)*cvmgt(porousp(ij),
     &		ac(ij),npc(ij).ne.1.and.npor.ne.0)
          energyp=energyp+(y(j-1)+f(ij)*dely(j)/2.0-h0)*(-gy)*
     &  	delx(i)*dely(j)*f(ij)*cvmgt(porousp(ij),
     &		ac(ij),npc(ij).ne.1.and.npor.ne.0)
	  if (ncyc.eq.0) energyp0=energyp
          xpollutant=xpollutant+xp(ij)*cvmgt(porousp(ij),ac(ij),
     &		npc(ij).ne.1.and.npor.ne.0)*delx(i)*dely(j)
C.......calculate force in x-direction
	  if (ac(ij+1).ne.1.0.or.ac(ij-1).ne.1.0) then
C.......  estimate distance between central point to obstacle surface
C.......  for smooth boundary (partial cell treatment)
          disco=((at(ij)+at(ij-imax))/2.0-0.5)*delx(i)
C.......  for stair type boundary
          if ((ar(ij)*ar(ij-1).eq.0.0.and.ar(ij)+ar(ij-1).eq.1.0).and.
     &      (at(ij)*at(ij-imax).eq.0.0.and.at(ij)+at(ij-imax).eq.1.0))
     &      disco=0.5*delx(i)
C.......  obstacle on the right 	
	  if (ac(ij+1).lt.ac(ij-1)) then
C.......    calculate slope of pressure
	    slopep=(p(ij)-p(ij-1))/((delx(i)+delx(i-1))/2.0)
C.......  obstacle on the left
	  else
C.......    calculate slope of pressure
            slopep=(p(ij)-p(ij+1))/((delx(i)+delx(i+1))/2.0)
	  endif
C.......  estimate p on obstacle based on linear inter(extra)polation
          pobstacle=p(ij)+slopep*disco
  	  xforce=xforce+(ar(ij-1)-ar(ij))*pobstacle*dely(j)
	  endif

C.......calculate force in y-direction
          if (ac(ij+imax).ne.1.0.or.ac(ij-imax).ne.1.0) then
C.......  estimate distance between central point to obstacle surface
C.......  for smooth boundary (partial cell treatment)
          disco=((ar(ij)+ar(ij-1))/2.0-0.5)*dely(j)
C.......  for stair type boundary
          if ((ar(ij)*ar(ij-1).eq.0.0.and.ar(ij)+ar(ij-1).eq.1.0).and.
     &      (at(ij)*at(ij-imax).eq.0.0.and.at(ij)+at(ij-imax).eq.1.0))
     &      disco=0.5*dely(j)
C.......  obstacle on the top
          if (ac(ij+imax).lt.ac(ij-imax)) then
C.......    calculate slope of pressure
            slopep=(p(ij)-p(ij-imax))/((dely(j)+dely(j-1))/2.0)
C.......  obstacle on the bottom
          else
C.......    calculate slope of pressure
            slopep=(p(ij)-p(ij+imax))/((dely(j)+dely(j+1))/2.0)
          endif
C.......  estimate p on obstacle based on linear inter(extra)polation
          pobstacle=p(ij)+slopep*disco
	  yforce=yforce+(at(ij-imax)-at(ij))*pobstacle*delx(i)
	  endif
  171   continue

	  write(12,550)t,xmass,energyk,energyp-energyp0,energyt,dissipt,
     &	dissipm,productt,xforce,yforce,xsaka,xnewloc,ynewloc
550     format(f8.4,12e15.7)
C==========================================================================
c
   60 continue

        if (nmean.eq.1) then
	  if (ncyc.eq.0.or.(t.ge.tmstart.and.t.le.tmend+1.0e-6)) then
            do 191 i=istart,iend
              eflux(i)=0.0
              efluxt(i)=0.0
              disp(i)=0.0
              ek(i)=0.0
              et(i)=0.0
              ep(i)=0.0
              etaa(i)=0.0
              do 192 j=jm1,2,-1
                ij=(j-1)*imax+i
                imj=ij-1
                ijm=ij-imax
C...............energy flux
                eflux(i)=eflux(i)+(p(ij)+rhof*
     &                  gy*(h0-(y(j-1)+f(ij)*dely(j)/2.0)))
     &                  *(u(ij)+u(imj))/2.
     &                  *dely(j)*f(ij)*cvmgt(porousp(ij),
     &                  ac(ij),npc(ij).ne.1.and.npor.ne.0)
                efluxt(i)=efluxt(i)+p(ij)
     &                  *(u(ij)+u(imj))/2.
     &                  *dely(j)*f(ij)*cvmgt(porousp(ij),
     &                  ac(ij),npc(ij).ne.1.and.npor.ne.0)
C...............dissipation
                if (kemodel.gt.0) then
                  disp(i)=disp(i)+rhof*xep(ij)*delx(i)
     &                  *dely(j)*f(ij)*cvmgt(porousp(ij),
     &                  ac(ij),npc(ij).ne.1.and.npor.ne.0)
                endif
C...............internal energy
                ek(i)=ek(i)+rhof*(((u(ij)+u(imj))/2.)**2/2.0+
     &                  ((v(ij)+v(ijm))/2.)**2/2.0)
     &                  *delx(i)*dely(j)
     &                  *f(ij)*cvmgt(porousp(ij),
     &                  ac(ij),npc(ij).ne.1.and.npor.ne.0)
                if (kemodel.gt.0) then
                  et(i)=et(i)+rhof*xk(ij)
     &                  *delx(i)*dely(j)
     &                  *f(ij)*cvmgt(porousp(ij),
     &                  ac(ij),npc(ij).ne.1.and.npor.ne.0)
                endif
                ep(i)=ep(i)+rhof*(y(j-1)+f(ij)*(y(j)-y(j-1))/2.0-h0)
     &                  *(-gy)*delx(i)*dely(j)
     &                  *f(ij)*cvmgt(porousp(ij),
     &                  ac(ij),npc(ij).ne.1.and.npor.ne.0)
C...............mean free surface
                etaa(i)=etaa(i)+f(ij)*dely(j)
192           continue
	      if (ncyc.eq.0) ep0(i)=ep(i)
              etaa(i)=etaa(i)-h0
              ep(i)=ep(i)-ep0(i)
191         continue
	      if (t.lt.tmstart) go to 193
            write(70,'(300e14.6)')t,(eflux(i),i=istart,iend,iinterval)
C            if (kemodel.gt.0) write(71,'(3000e14.6)')
C     &          t,(disp(i),i=istart,iend)
            write(72,'(300e14.6)')t,((ek(i)+et(i)+ep(i)),i=istart,iend)
            write(73,'(300e14.6)')t,(etaa(i),i=istart,iend)
            write(74,'(300e14.6)')t,(efluxt(i),i=istart,iend,iinterval)
            tmstart=tmstart+tinterval
	  endif
        endif

193	continue

	if (nmean.eq.2) then
	if (t.ge.twfin-xxt*2.0) then
	 if (t.le.twfin) then
	  xmnt=xmnt+delt
	  do i=2,im1
	  do j=2,jm1
		ij=(j-1)*imax+i
		xmnf(ij)=xmnf(ij)+f(ij)*delt
		xmnu(ij)=xmnu(ij)+(u(ij)*ar(ij)+u(ij-1)*ar(ij-1))*delt*f(ij)
     &		/2.0
		xmnv(ij)=xmnv(ij)+(v(ij)*at(ij)+v(ij-imax)*at(ij-imax))
     &		*delt*f(ij)/2.0
		xmnk(ij)=xmnk(ij)+xk(ij)*delt*f(ij)
	  end do
	  end do
	 else
	  write(*,*)'xmnt=',xmnt
	  write(9,*)'xmnt=',xmnt
	  open(3000,file='mnf',status="unknown")
	  open(3001,file='mnu',status="unknown")
	  open(3002,file='mnv',status="unknown")
	  open(3003,file='mnk',status="unknown")
	  do j=2,jm1
	    write(3000,'(3000f8.3)') (xmnf((j-1)*imax+i)/xmnt,i=2,im1)
		write(3001,'(3000f8.3)') (xmnu((j-1)*imax+i)/xmnt,i=2,im1)
		write(3002,'(3000f8.3)') (xmnv((j-1)*imax+i)/xmnt,i=2,im1)
		write(3003,'(3000f8.3)') (xmnk((j-1)*imax+i)/xmnt,i=2,im1)
	  end do
	 endif
	endif
	endif

c output surface
      if (t+em6.lt.sfprt) go to 172
       sfprt=sfprt+sfdt
       call surface
172    continue      


       goto 1909
c - put bubble fyshi
      if(t+em6.le.delt)then
       print*,'initial bubble ...',t
      call init_bubble
      else
       print*,'call bubble'
      call bubble_source
      call bubble_all
      endif
1909  continue      

      if (t+em6.lt.twprt) go to 72
      twprt=twprt+prtdt
      if (ninflow.eq.9.or.j_det.ne.0) then
        j_det=j_det+1
        prtdt=det_k(j_det)
      endif
c
c.... print field variable data on paper
      print*,'prtdt=',prtdt,t
c
      goto 1911
      do k_bub=1,mbub
      call prtplt_bubble(k_bub)
      enddo
      call prtplt_any 
1911  continue      
c     ----------------
   70 call prtplt2 (3)
c     ----------------
c
72    if (nanimation.eq.0) go to 75
      if (t+em6.lt.twprt_a) go to 75
      twprt_a=twprt_a+prtdt_a
c
c.... print field variable data on paper
c
c     ----------------
c      call prtplt (4)
c     ----------------
c
   75 if (t+em6.lt.twdmp) go to 80
      twdmp=twdmp+dmpdt
c
c.... check to see if problem finish time surpassed
c
   80 if (t.gt.twfin) go to 130

	if (kt.eq.7.and.t.gt.aa) then
		xnu=0.0
		xmu=0.0
	endif
	if (kt.eq.7.and.t.gt.h0) kt=1
c
c.... set the advance time arrays into the time-n arrays
c
C==========================================================================
C.....store the old time step value
      do 90 j=1,jmax
        do 85 i=1,imax
          ij=(j-1)*imax+i
          un(ij)=u(ij)
          vn(ij)=v(ij)
	    xnutn(ij)=xnut(ij)
	    xnutyn(ij)=xnuty(ij)
	    xpn(ij)=xp(ij)
	    xkn(ij)=xk(ij)
	    xepn(ij)=xep(ij)
          pn(ij)=p(ij)
          fn(ij)=f(ij)
          nfold(ij)=nf(ij)
   85 continue
   90  continue

C==========================================================================
c
c.... adjust the time step (delt)
c
c     ----------------
  100 call deltadj
c     ----------------
c
c.... advance time
c
      t=t+delt
c
c.... check for vof function convection limit,
c     terminate if it has been exceeded
c
      if (nflgc.lt.10000) go to 110
c
      write (9,170) ncyc,t
      write (iotty,170) ncyc,t
c
c.... exit from newcyc and terminate
c
c     --------------------------
c	routine="NEWCYC"
c      call kill (iotty,ncyc,routine,nexit)
c     --------------------------
c
c.... check for pressure convergence failure,
c     terminate if number of failed cycles
c     (noncon) is nonzero
c
  110 if (nocon.lt.5) go to 120
c
      write (9,160) ncyc,t
      write (iotty,160) ncyc,t
c
c.... exit from newcyc and terminate
c
c     --------------------------------
c	routine="NEWCYC"
c      call kill (iotty,ncyc,routine,nexit)
c     --------------------------------
c
c.... everything ok, advance cycle
c
  120 ncyc=ncyc+1
c
c
	time=0.0
      tleft=ttl-time-tquit
c
      if (mod(ncyc,25).ne.0) go to 9999
c
	tend=0.0
	tendd=0.0
      grind=1000.*(tend+tend-tendd-tbeg)/float(ibar*jbar)
      write (13,140) t,ncyc,grind,iter
      write (iotty,140) t,ncyc,grind,iter
      go to 9999
c
c.... normal termination (t > twfin)
c
130	continue
      write (iotty,180)
      write (13,180)
c.... exit and destroy the dropfile
c
c     --------------
      nexit=9999
c     --------------
c
c.... check timestep: if too small (delt < dtend),
c     terminate run on next cycle
c
 9999 if(delt.gt.dtend) go to 20

	write(9,*)delt,dtend,twfin,t

      write (iotty,150) ncyc,t
      write (9,150) ncyc,t
	
      twfin=t*0.999
c
   20 continue
c	iter=0
c
      return
c                    * * * * error section * * * *
  140 format (2x,"t = ",1pe12.5,"  cycle = ",i5,2x,"grind = ",1pe12.5,
     &        " ms","  iter = ",i3)
  150 format (1x,29hdelt less than dtend on cycle,i6,1x,3ht =,1pe15.7)
  160 format (//1x,25htoo many pressit failures,3x,5hncyc=,i7,1x,2ht=,1p
     1         e14.6//)
  170 format (//1x,24htoo many vofadv failures,3x,5hncyc=,i7,1x,2ht=,1pe
     1         14.6//)
  180 format (/,18x,"* * * * Normal Termination * * * *")
      end

