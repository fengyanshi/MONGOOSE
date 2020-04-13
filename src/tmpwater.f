
        subroutine accel
c
c ======================================================================
c
c   Purpose -
c     compute Lagrangian velocities caused by the implicit
c     pressure and explicit viscous forces
c
c   ACCEL is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c         IMPLCTP
c
c
c   ACCEL calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c              BC    ripple
c
c ======================================================================
c
c##############################################################
      implicit real*8 (a-h,o-z)
c##############################################################
c
c############
      include  "comdk2.h"      
      include   "bubble1.h"
c############
c
      data tiny /1.0d-25/, zero /0.0d0/
      dimension dpdx_t(10,250),dpdy_t(10,250)
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      do 20 j=1,jm1
        m=1
        n=1
        do 20 i=1,im1
c
          ij=(j-1)*imax+i
          ipj=ij+1
          imj=ij-1
          ijp=ij+imax
          ijm=ij-imax
          ipjp=ipj+imax
          ipjm=ipj-imax
          imjp=imj+imax
          imjm=imj-imax
c
          rhoij=f(ij)*rhof
          rhoipj=f(ipj)*rhof
          rhoijp=f(ijp)*rhof
          rhorc=(delx(i+1)*rhoij+delx(i)*rhoipj)/(delx(i)+delx(i+1))
          rhotc=(dely(j+1)*rhoij+dely(j)*rhoijp)/(dely(j)+dely(j+1))
          rhox=rhorc+tiny
          rhoy=rhotc+tiny
c
          if (j.eq.1.and.kb.eq.3) goto 15
c
c....     check for surrounding fluid
          rhobar=rhox
          if (ar(ij).lt.em6) go to 10
c
          rhox=rhobar*(delx(i+1)+delx(i))

          if (npor.eq.0) then
            gctmpra=1.0
            gctmpr=1.0
          endif
          if (npor.ne.0) then
            if (ar(ij).gt.em6.and.(npc(ij).ne.1.or.npc(ipj).ne.1)) 
     &		then
                gctmpra=(porousc(ij)+porousc(ipj))/2.0
                if (npc(ij).ne.1) then
                  xpstmp=porousp(ij)
                  xatmp=porousa(ij)
                  xxbtmp=porousb(ij)
                  d50tmp=porousd(ij)
                  xkc=max(abs(un(ij)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict0=(xatmp+xbb
     &			*sqrt(un(ij)**2+vn(ij)**2))*abs(gy)
                else
                  xfrict0=0.0
                endif
                if (npc(ipj).ne.1) then
                  xpstmp=porousp(ipj)
                  xatmp=porousa(ipj)
                  xxbtmp=porousb(ipj)
                  d50tmp=porousd(ipj)
                  xkc=max(abs(un(ij)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict1=(xatmp+xbb
     &			*sqrt(un(ij)**2+vn(ipj)**2))*abs(gy)
                else
                  xfrict1=0.0
                endif
			  xfriction=(xfrict0+xfrict1)/2.0
                gctmpr=1.0/(gctmpra+xfriction*delt)
	    else
                gctmpra=1.0
			  gctmpr=1.0
	    endif
          endif

C.......adding sponge layer (Wei and Kirby, 1995)
	  if (nopen.eq.11) then
		if (kl.eq.3) then
		  if (x(i).le.xsponge) then
			udamp=(exp(((xsponge-x(i))/xsponge)**power)-1.0)
     &			/(exp(1.0)-1.0)
C.....re-define gctmpr [1/(c+dt*f)]
			fricold=(1./gctmpr-gctmpra)/delt
			fricnew=fricold+adamp/xxt*udamp
			gctmpr=1.0/(gctmpra+delt*fricnew)
		  endif
		endif
		if (kr.eq.3) then
		  if (x(i).ge.x(im1)-xsponge) then
			udamp=(exp(((-x(im1)+x(i)+xsponge)/xsponge)**power)-1.0)
     &			/(exp(1.0)-1.0)
C.....re-define gctmpr [1/(c+dt*f)]
			fricold=(1./gctmpr-gctmpra)/delt
			fricnew=fricold+adamp/xxt*udamp
			gctmpr=1.0/(gctmpra+delt*fricnew)
		  endif
		endif
	  endif

          u(ij)=gctmpra*gctmpr*u(ij)+delt*gctmpr*(p(ij)-p(ipj))*2.0
     &		/rhox*fact+delt*gctmpr*gx

   10     continue
	    u(ij)=cvmgt(zero,u(ij),((rhorc.lt.frsurf).AND.
     &          ((nf(ij).gt.5).and.(nf(ipj).gt.5))).or.
     &		(ar(ij).lt.em6))
C.........since frsurf is normally smaller than em6, which is used to
C.........clear out small fluid cell, there's a small chance that a cell
C.........may have finite f(ij) but velocity equals to zero based on the
C.........following criteria. However, it is found when pressure-driven
C.........wavemaker is used, a large velocity can occur at surface cell
C.........if the above condition is used; inconsistency calls further work
	    if (ninflow.eq.54) then
            u(ij)=cvmgt(zero,u(ij),((rhorc.lt.frsurf).OR.
     &          ((nf(ij).gt.5).and.(nf(ipj).gt.5))).or.
     &          (ar(ij).lt.em6))
	  endif
c
15	  continue
        if (i.eq.1.and.(kl.eq.3.or.kl.eq.6.or.kl.eq.4)) goto 20
c
c....   reset y-velocity and check for surrounding fluid
        rhobar=rhoy
        if (at(ij).lt.em6) go to 25
c
        rhoy=rhobar*(dely(j+1)+dely(j))

	  if (npor.eq.0) then
            gctmpta=1.0
            gctmpt=1.0
	  endif
        if (npor.ne.0) then
            if (at(ij).gt.em6.and.(npc(ij).ne.1.or.npc(ijp).ne.1)) 
     &		then
                gctmpta=(porousc(ij)+porousc(ijp))/2.0
                if (npc(ij).ne.1) then
                  xpstmp=porousp(ij)
                  xatmp=porousa(ij)
                  xxbtmp=porousb(ij)
                  d50tmp=porousd(ij)
                  xkc=max(abs(vn(ij)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict0=(xatmp+xbb
     &			*sqrt(un(ij)**2+vn(ij)**2))*abs(gy)
                else
                  xfrict0=0.0
                endif
                if (npc(ijp).ne.1) then
                  xpstmp=porousp(ijp)
                  xatmp=porousa(ijp)
                  xxbtmp=porousb(ijp)
                  d50tmp=porousd(ijp)
                  xkc=max(abs(vn(ij)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict1=(xatmp+xbb
     &			*sqrt(un(ijp)**2+vn(ij)**2))*abs(gy)
                else
                  xfrict1=0.0
                endif
			  xfriction=(xfrict0+xfrict1)/2.0
                gctmpt=1.0/(gctmpta+xfriction*delt)
	    else
                gctmpta=1.0
			  gctmpt=1.0
          endif
        endif

c
C.......adding sponge layer
	  if (nopen.eq.11) then
		if (kl.eq.3) then
		  if (xi(i).le.xsponge) then
			vdamp=(exp(((xsponge-xi(i))/xsponge)**power)-1.0)
     &			/(exp(1.0)-1.0)
C.....re-define gctmpt [1/(c+dt*f)]
			fricold=(1./gctmpt-gctmpta)/delt
			fricnew=fricold+adamp/xxt*vdamp
			gctmpt=1.0/(gctmpta+delt*fricnew)
		  endif
		endif
		if (kr.eq.3) then
		  if (xi(i).ge.x(im1)-xsponge) then
			vdamp=(exp(((-x(im1)+xi(i)+xsponge)/xsponge)**power)-1.0)
     &			/(exp(1.0)-1.0)
C.....re-define gctmpt [1/(c+dt*f)]
			fricold=(1./gctmpt-gctmpta)/delt
			fricnew=fricold+adamp/xxt*vdamp
			gctmpt=1.0/(gctmpta+delt*fricnew)
		  endif
		endif
	  endif

	  v(ij)=gctmpta*gctmpt*v(ij)+delt*gctmpt*(p(ij)-p(ijp))*2.0
c     &		/rhoy*fact+delt*gctmpt*gy
c bubble buoyoncy added fyshi
     &          /rhoy*fact+delt*gctmpt*gy*(1.-alpha_g(ij))
   25   continue
	  v(ij)=cvmgt(zero,v(ij),((rhotc.lt.frsurf).AND.
     &          ((nf(ij).gt.5).and.(nf(ijp).gt.5))).or.
     &		(at(ij).lt.em6))
	  if (ninflow.eq.54) then
            v(ij)=cvmgt(zero,v(ij),((rhotc.lt.frsurf).OR.
     &          ((nf(ij).gt.5).and.(nf(ijp).gt.5))).or. 
     &          (at(ij).lt.em6))
	  endif

   20 continue
c
C	do j=3,jm1-1
C	do i=3,im1-1
C	ij=(j-1)*imax+i
C	ijm=ij-imax
C	imj=ij-1
C	if (f(ij).gt.0.0d0.and.ac(ij).gt.0.0.and.i.gt.50) then
C	if ((u(ij)*ar(ij)-u(imj)*ar(imj))/delx(i)+(v(ij)*at(ij)
C     &	-v(ijm)*at(ijm))/dely(j).gt.1.0e-4) then
C	write(9,*)'!!!1',i,j,f(ij),f(ij+1),f(imj),f(ij+imax),f(ijm),
C     &	f(ij+1+imax),f(imj+imax),f(ijm+1),f(ijm-1),
C     &	nf(ij),nf(ij+1),nf(imj),nf(ij+imax),nf(ijm),
C     &	nf(ij+1+imax),nf(imj+imax),nf(ijm+1),nf(ijm-1),
C     &	ac(ij),ar(ij),ar(imj),at(ij),at(ijm),
C     &	u(ij),u(imj),v(ij),v(ijm),
C     &	(u(ij)*ar(ij)-u(imj)*ar(imj))/delx(i)+(v(ij)*at(ij)
C     &	-v(ijm)*at(ijm))/dely(j)
C	endif
C	endif
C	end do
C	end do
c
c.... update the boundary conditions
c
C.....It is tricky here whether we should adjust free surface velocity 
C.....at this moment (set ibcflg=1). The original RIPPLE does not. It is
C.....found that by seting ibcflg=1, mass increases slightly but free
C.....surface computation is better; without doing this, mass drops a
C.....little and the free surface is less smooth. Considering that
C.....the free surface boundary condition is approximate which is basically
C.....to ensure div(u)=0, we decided not to enforce it at this moment.
      ibcflg=1
C.....added to account for the impact
C     ibcfinal=1
c     ---------
      call bc
c     ---------
      ibcflg=0
C     ibcfinal=0
c
C	stop
C	do j=3,jm1-1
C	do i=3,im1-1
C	ij=(j-1)*imax+i
C	ijm=ij-imax
C	imj=ij-1
C	if (f(ij).gt.0.0d0.and.ac(ij).gt.0.0.and.i.gt.50) then
C	if ((u(ij)*ar(ij)-u(imj)*ar(imj))/delx(i)+(v(ij)*at(ij)
C     &	-v(ijm)*at(ijm))/dely(j).gt.1.0e-5) then
C	write(9,*)'!!!2',i,j,f(ij),f(ij+1),f(imj),f(ij+imax),f(ijm),
C     &	f(ij+1+imax),f(imj+imax),f(ijm+1),f(ijm-1),
C     &	nf(ij),nf(ij+1),nf(imj),nf(ij+imax),nf(ijm),
C     &	nf(ij+1+imax),nf(imj+imax),nf(ijm+1),nf(ijm-1),
C     &	ac(ij),ar(ij),ar(imj),at(ij),at(ijm),
C     &	u(ij),u(imj),v(ij),v(ijm),
C     &	(u(ij)*ar(ij)-u(imj)*ar(imj))/delx(i)+(v(ij)*at(ij)
C     &	-v(ijm)*at(ijm))/dely(j)
C	endif
C	endif
C	end do
C	end do

      return
      end

       real*8 function dasum(n,dx,incx)                       
c
c***purpose  sum of magnitudes of d.p. vector components              
c***description                                                      
c                                                                   
c                b l a s  subprogram                               
c    description of parameters                                    
c                                                                
c     --input--                                                 
c        n  number of elements in input vector(s)              
c       dx  double precision vector with n elements           
c     incx  storage spacing between elements of dx           
c                                                           
c     --output--                                           
c    dasum  double precision result (zero if n .le. 0)    
c                                                        
c     returns sum of magnitudes of double precision dx. 
c     dasum = sum from 0 to n-1 of dabs(dx(1+i*incx))  
c***references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,     
c                 *basic linear algebra subprograms for fortran usage*, 
c                 algorithm no. 539, transactions on mathematical       
c                 software, volume 5, number 3, september 1979, 308-323 
c***routines called  (none)                                             
c***end prologue  dasum                                                 
c                                                                       
      real*8 dx(1)                                         
c***first executable statement  dasum                               
      dasum = 0.d0                                                 
      if(n.le.0)return                                            
      if(incx.eq.1)goto 20                                       
c                                                               
c        code for increments not equal to 1.                   
c                                                             
      ns = n*incx                                            
          do 10 i=1,ns,incx                                 
          dasum = dasum + dabs(dx(i))                      
   10     continue                                        
      return                                             
c                                                       
c        code for increments equal to 1.               
c                                                     
c                                                    
c        clean-up loop so remaining vector length is a multiple of 6.   
c                                                                      
   20 m = mod(n,6)                                                    
      if( m .eq. 0 ) go to 40                                          
      do 30 i = 1,m                                                   
         dasum = dasum + dabs(dx(i))                                 
   30 continue                                                      
      if( n .lt. 6 ) return                                        
   40 mp1 = m + 1                                                 
      do 50 i = mp1,n,6                                          
         dasum = dasum + dabs(dx(i)) + dabs(dx(i+1)) + dabs(dx(i+2))  
     1   + dabs(dx(i+3)) + dabs(dx(i+4)) + dabs(dx(i+5))             
   50 continue                                                      
      return                                                       
      end                                                         
  
 
      subroutine aset(ntype)
c
c ======================================================================
c
c   Purpose -
c     Initialize any interior obstacles with the function f(x,y);
c     for f(x,y) < 0, the point (x,y) is either open or closed
c     to the flow for ioh=0 or 1, respectively
c
c   ASET is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c           SETUP
c
c
c   ASET calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c             OBS    ripple   SETARRY    ripple       FXY    ripple
c
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
      dimension iflg(5),dis(4),xm(5),ym(5)
c
      data zero /0.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      xmode=2.*pi/(x(im1)-x(1))
      ymode=2.*pi/(y(jm1)-y(1))
      rad2o2=sqrt(2.0d0)/2.0d0
c
c     -------------------------------------
C      call setarry (a(1),0.0d0,nxy)
c     -------------------------------------
c
      if (nobs(ntype).le.0) go to 240
c
      do 230 k=(ntype-1)*20+1,(ntype-1)*20+nobs(ntype)
        do 221 j=1,jmax
          do 220 i=1,imax
            ij=(j-1)*imax+i
            if (noc(ij).eq.1) goto 220
            rdxdy=1.0/(delx(i)*dely(j))
            do 60 m=1,4
              go to (10,20,30,40), m
   10         x1=x(i)
              y1=cvmgt(y(j)-dely(j),y(j-1),j.eq.1)
              dis(1)=dely(j)
              go to 50
   20         y1=y(j)
              x1=x(i)
              dis(2)=delx(i)
              go to 50
   30         x1=cvmgt(x(i)-delx(i),x(i-1),i.eq.1)
              y1=y(j)
              dis(3)=dely(j)
              go to 50
   40         y1=cvmgt(y(j)-dely(j),y(j-1),j.eq.1)
              x1=cvmgt(x(i)-delx(i),x(i-1),i.eq.1)
              dis(4)=delx(i)
   50         iflg(m)=0
              fconic=fxy(x1,y1,oa1(k),oa2(k),ob1(k),ob2(k),oc1(k),
     &                   oc2(k),od1(k),nxo(k)*xmode,od2(k),mxo(k)*xmode,
     &                   oe1(k),nyo(k)*ymode,oe2(k),myo(k)*ymode)
              if (abs(fconic).le.1.0d-12) fconic=0.0
              if (fconic.le.0.0) iflg(m)=1
              xm(m)=x1
              ym(m)=y1
   60       continue
c
            iflg(5)=iflg(1)
            xm(5)=xm(1)
            ym(5)=ym(1)
            iflgs=0
c
            do 70 m=1,4
              iflgs=iflgs+iflg(m)
   70       continue
c
            brij=0.0
            btij=0.0
            if (iflgs.eq.0) go to 220
            if (iflgs.lt.4) go to 80
            bij=1.0d0
            brij=1.0d0
            btij=1.0d0
            go to 200
   80       if (iflg(1).eq.1.and.iflg(2).eq.1) brij=1.0d0
            if (iflg(2).eq.1.and.iflg(3).eq.1) btij=1.0d0
c
            do 160 m=1,4
              if (iflg(m).eq.iflg(m+1)) go to 160
              x1=xm(m)
              y1=ym(m)
              x2=xm(m+1)
              y2=ym(m+1)
              if (iflg(m).eq.0) go to 90
              x2=xm(m)
              y2=ym(m)
              x1=xm(m+1)
              y1=ym(m+1)
   90         epsif=0.001d0*(abs(x2-x1)+abs(y2-y1))
              smn=0.0
              fmn=fxy(x2,y2,oa1(k),oa2(k),ob1(k),ob2(k),oc1(k),
     &                oc2(k),od1(k),nxo(k)*xmode,od2(k),mxo(k)*xmode,
     &                oe1(k),nyo(k)*ymode,oe2(k),myo(k)*ymode)
              smx=1.0d0
              fmx=fxy(x1,y1,oa1(k),oa2(k),ob1(k),ob2(k),oc1(k),
     &                oc2(k),od1(k),nxo(k)*xmode,od2(k),mxo(k)*xmode,
     &                oe1(k),nyo(k)*ymode,oe2(k),myo(k)*ymode)
              s=0.5
  100         xt=s*x1+(1.0-s)*x2
              yt=s*y1+(1.0-s)*y2
              fs=fxy(xt,yt,oa1(k),oa2(k),ob1(k),ob2(k),oc1(k),
     &                oc2(k),od1(k),nxo(k)*xmode,od2(k),mxo(k)*xmode,
     &                oe1(k),nyo(k)*ymode,oe2(k),myo(k)*ymode)
              if (abs(fs).lt.epsif) go to 130
              if (fs.ge.0.0) go to 110
              fden=abs(fs-fmn)+1.0d-10
              se=s-fs*(s-smn)/fden
              if (se.gt.smx) se=smx
              fmn=fs
              smn=s
              go to 120
  110         fden=abs(fmx-fs)+1.0d-10
              se=s-fs*(smx-s)/fden
              if (se.lt.smn) se=smn
              fmx=fs
              smx=s
  120         si=s-fs*(smx-smn)/(fmx-fmn)
              s=0.5*(se+si)
              go to 100
  130         dis(m)=sqrt((xt-x2)**2+(yt-y2)**2)
              go to (140,150,160,160), m
  140         brij=dis(1)/dely(j)
              go to 160
  150         btij=dis(2)/delx(i)
  160       continue
c
            m=0
            bij=0.0
  170       continue
            m=m+1
            if (m.eq.5) go to 190
            if (iflg(m).eq.0) go to 170
            mp1=m+1
            if (mp1.eq.5) mp1=1
            mm1=m-1
            if (mm1.eq.0) mm1=4
            bij=bij+dis(m)*dis(mm1)
            if (iflg(mp1).eq.1) go to 180
            dis2=dis(m)
  180       continue
            if (iflg(mm1).eq.1) go to 170
            dis1=dis(mm1)
            go to 170
  190       continue
            if (iflgs.eq.3) bij=bij-dis1*dis2
            bij=0.5*bij*rdxdy
            if (bij.gt.1.0) bij=1.0d0
  200       continue
            if (ioh(k).eq.0) go to 210
            bij=-bij
            brij=-brij
            btij=-btij
  210       ac(ij)=ac(ij)+bij
            if(ac(ij).gt.0.9999d0) ac(ij)=1.0d0
            if(ac(ij).lt.0.01d0) ac(ij)=0.0
            ar(ij)=ar(ij)+brij
            if(ar(ij).gt.0.9999d0) ar(ij)=1.0d0
            if(ar(ij).lt.0.01d0) ar(ij)=0.0
            at(ij)=at(ij)+btij
            if(at(ij).gt.0.9999d0) at(ij)=1.0d0
            if(at(ij).lt.0.01d0) at(ij)=0.0
C...........modified by lin to reduce the chance of very large velocity
55	      continue
C	      goto 56
C...........another way to make nearly blocked cell completely blocked
C...........gaurentee to avoid instability by partial cell; spurious reflection
C            if (ac(ij).le.0.20) then
C...........reduce instability possibilities; minimize spurious reflection
            if (ac(ij).le.acmin) then
	 		ac(ij)=0.0
			ar(ij)=0.0
			at(ij)=0.0
			ar(ij-1)=0.0
			at(ij-imax)=0.0
	      endif
56	    continue
  220     continue
  221   continue
  230 continue
c
c
c.... get the physical location of the obstacle boundary
c
c     ---------
C      call obs
c     ---------
c
  240 continue
c
      do 311 j=2,jm1
        do 310 i=2,im1
          ij=(j-1)*imax+i
          imj=ij-1
          ijm=ij-imax
          if (ac(ij).gt.em6) go to 305
          ar(ij)=0.0
          ar(imj)=0.0
          at(ij)=0.0
          at(ijm)=0.0
  305		continue
		if (ac(ij)*ar(ij)*ar(imj)*at(ij)*at(ijm).lt.0.999999.
     &	and.noc(ij)+noc(imj)+noc(ij+1)+noc(ijm)+noc(ij+imax).eq.0) 
     &	nmovbd(ij)=ntype
  310   continue
  311 continue
c
      do 450 n=1,nxy
C.......works always but can cause stair-type edge 
C.......if multi-structures are not properly set in sequence 
        if (ac(n).eq.0.0) noc(n)=1
C.......create smooth edge but introduce gap between two structures 
C.......if structure does not align with cell boundary
C	  if (ac(n).ne.1.0) noc(n)=1
450   continue
c
c.... set ar and at values for inflow and outflow boundary segments
c     here as update modification depending on application
c
      return
      end


      subroutine bdycell(nx,ny,wl,wr,wt,wb,bcl,bcr,bct,bcb,
     &                    bound,p,xi,yj,gx,gy)
c
c ======================================================================
c
c   Purpose -
c     fills the ghost cells with the boundary quantities
c
c   BDYCELL is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c         IMPLCTP  
c
c
c   BDYCELL calls the following subroutines and functions -
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
c -- I enlarge the dimension for p(nx+ny)
      dimension p((nx+1)*(ny+1)),bound(nx+ny),xi(nx+1),yj(ny+1)
      integer wl,wr,wt,wb
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      nxp=nx+1
      nyp=ny+1
      ijt=ny*nxp+2-nxp
      ijb=2+nxp
      ijl=2
      ijr=nx
 
c.... bottom boundary
c
      ij=2
      do 300 i=2,nx
        p(ij)=bcb*bound(1)+(1.-bcb)*p(ij+nx+1)
        ij=ij+1
        ijt=ijt+1
  300 continue
c
c.... top boundary
c
      ij=ny*(nx+1)+2
      do 400 i=2,nx
        p(ij)=bct*bound(2)+(1.-bct)*p(ij-nx-1)
        ij=ij+1
        ijb=ijb+1
  400 continue
c
c.... left boundary
c
      ij=1
      do 100 j=1,ny+1
        p(ij)=bcl*(bound(3)+gy*yj(j))+(1.-bcl)*p(ij+1)
        ij=ij+nx+1
        ijr=ijr+nx+1
  100 continue
c
c.... right boundary
c
      ij=nx+1
      do 200 j=1,ny+1
        p(ij)=bcr*(bound(4)+gy*yj(j))+(1.-bcr)*p(ij-1)
        ij=ij+nx+1
        ijl=ijl+nx+1
  200 continue
c
      return
      end


      subroutine begin
c
c ======================================================================
c
c   Purpose -
c     get job identification, run-time limit, date, time of day,
c     machine, and initialize graphics
c     tape5:  Input file (INPUT)
c     tape8:  Binary restart dump file
c     tape9:  Error message output file
c     tape10:  Binary restart read file
c     tape13:  General edit output file
c     6 or 59:  Standard output
c
c
c   BEGIN is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE
c
c
c   BEGIN calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c            MACH    system    SECOND    system    GETJTL    system
c         SETFLSH    system     GPLOT    system   LIBDISP    system
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
c.... open the input and output files
c
c     ------------------------------------------------
      open (unit=5,file="input",form="formatted",status="unknown")
      open (unit=8,file="tape8",form="unformatted",status="unknown")
      open (unit=9,file="errors",form="formatted",status="unknown")
      open (unit=10,file="tape10",form="unformatted",status="unknown")
	open (12,file='conserve',form="formatted",status="unknown")
      open (unit=13,file="edit",form="formatted",status="unknown")
C	open (21,file='vof',status="unknown")
C	open (22,file='ux',status="unknown")
C	open (23,file='vy',status="unknown")
C	open (24,file='ps',status="unknown")
C      open (29,file='vortex',status="unknown")
	open (31,file='xc',status="unknown")
	open (32,file='yc',status="unknown")
c
      return
      end



      subroutine cgitj(a,b,c,d,e,af,bf,cf,df,ef,
     1     n,m,x,y,r,up,aux,istop,ercg,
     2     xratio,ratio,xnorm,ynorm)
c
c***purpose
c
c      given a 9-point matrix multiply routine called matmul9
c      and an approximate factorization a = l*d*lt, this routine
c      provides a conjugate gradient acceleration to iterate
c      to an approximate solution to a*x=y.
c***description
c
c       cgitj requires that the coefficients of the ith
c       equation be stored in the ith element of the diagonal
c       arrays. five diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c       the ldlt factors are stored similarly in five
c      lower diagonal arrays.
c
c     on entry
c
c        a,b,c,d,e   are 1-d arrays of length
c            at least nxm. the arrays contain, from left to right,
c            the upper diagonals of a as one sweeps a left to right.
c
c        el,dl,cl,bl,af   are 1-d arrays of length
c            at least nxm. the arrays contain from left to right
c            the lower diagonals of l,d as they are swept left to right.
c
c        n     is the size  or order of each tridiagonal block if
c              a is viewed as a block tridiagonal system of equations.
c
c        m     is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c
c        x     is an estimate to the solution vector of size nxm.
c
c        y     is an r.h.s. vector of dimenison at least nxm.
c
c        r,up,aux  are work variables at least nxm in length.
c               the final residual is left array r.
c
c        istop  is the maximum number of iterations allowed.
c               upon exit it contains the actual number of
c               iterations required to converge to accuracy
c               specified.
c
c        ercg   is the error tolerance for convegence
c
c     on return
c
c        x     is the solution vector of dimension at least
c              nxm containing the solution vector.
c
c        xratio   is the ratio of dxnorm to xnorm.
c
c        ratio    is the ratio of rnorm to ynorm
c
c        xnorm    is the l2 norm of the solution vector x.
c
c        ynorm    is the l2 norm of the r.h.s.
c
c***routines called  ldlt9f matmul9


      implicit real*8 (a-h,o-z)
      real*8 e(m*n),d(m*n),c(m*n),b(m*n),a(m*n),ef(m*n),df(m*n),cf(m*n),
     &       bf(m*n),af(m*n),x(m*n),y(m*n),r(m*n),up(m*n),aux(m*n),
     &	   el(m*n),dl(m*n),cl(m*n),bl(m*n)
c  el,dl,cl,bl is added by Fengyan

      mn = m*n
      e3=ercg

c --- the linux f77 compiler can not make b(0) with the defination b(m*n)
c      so use nfun instead of 0 

	nfun=0


      do  3 i=1,n+1
        ef(i) = 0.
        df(i) = 0.
        cf(i) = 0.
    3 continue
      do 5 i=1,mn-n-1
        ef(i+n+1) = e(i)
        cf(i+n) = c(i+1)
    5 continue
      do 10 i= 1,mn-n
        df(i+n) = d(i)
   10 continue
      do 15 i=1,mn-1
        bf(i+1) = b(i)
   15 continue
      do 20 i=1,mn
        af(i) = a(i)
   20 continue

c added by Fengyan
	do i=1,mn-1
		bl(i+1) = b(i)
	enddo
	do i=1,mn-n+1
		cl(i+n-1) = c(i)
	enddo
	do i=1,mn-n
		dl(i+n) = d(i)
	enddo
	do i=1,mn-n-1
		el(i+n+1) = e(i)
	enddo

c
ccc periodic block copy
c
      do 22 i=2,n
        ef(i) = e(mn-n+i-1)
   22 continue
      do 25 i=1,n
        df(i) = d(mn-n+i)
   25 continue
      do 28 i=1,n-1
        cf(i) = c(mn-n+i+1)
   28 continue
      call ldlt9f(n,m,ef,df,cf,bf,af)
c
c initialize the residual r
c
C      call matmul9(e(-n),d(1-n),c(2-n),b(0),a,b,c,d,e,n,m,x,r)
C      call matmul9(e(-n),d(1-n),c(2-n),0.0d0,a,b,c,d,e,n,m,x,r)
CC      call matmul9(e,d,c,b(0),a,b,c,d,e,n,m,x,r)
c Fengyan changed the first call to following call
      call matmul9(el,dl,cl,bl,a,b,c,d,e,n,m,x,r)


      do 30 i=1,mn
        r(i) = y(i)-r(i)
   30 continue
      if (dasum(mn,r,1) .eq. 0.)     then
        i = 0
        go to 200
      end if
      ynorm = dasum(mn,y,1)
c
ccc solve l * d * l(trans)*aux = r0
c
      call ldlt(n,m,ef,df,cf,bf,af,bf(2),cf(n),df(n+1),ef(n+2),
     1          up,r)
c
ccc dot r with aux
c
      rdot = ddot(mn,r,1,up,1)
c
c put up into aux (p0)
c
      do 40 i4=1,mn
 40     aux(i4) = up(i4)
c p is stored in aux and mp or (ldlt)-1(r) are in up.
c
ccc begin main loop
c
      do 100 i = 1,istop
c
c---- find m * p and store it in up
c
c      call matmul9(e(-n),d(1-n),c(2-n),b(0),a,b,c,d,e,n,m,aux,up)
C      call matmul9(e(-n),d(1-n),c(2-n),0.0d0,a,b,c,d,e,n,m,aux,up)
CC      call matmul9(e,d,c,b(0),a,b,c,d,e,n,m,aux,up)
C Fengyan changed
      call matmul9(el,dl,cl,bl,a,b,c,d,e,n,m,aux,up)


c
ccc compute alpha(i)
c
      alpha = ddot(mn,up,1,aux,1)
      alpha = rdot/alpha
c
c
ccc compute x(i+1)
c
      do 80 i8 = 1,mn
 80     x(i8) = x(i8) + alpha*aux(i8)
c
ccc compute r(i+1)
c
      do 90 i9=1,mn
        r(i9) = r(i9) - alpha*up(i9)
 90   continue
c
ccc termination test
ccc
ccc stop whenever norm(aux) / norm (x) .lt. e3
c
      rnorm = dasum(mn,r,1)
      dxnorm=dasum(mn,aux,1)*dabs(alpha)
      xnorm=dasum(mn,x,1)
c***
ccc output ratio
c
      ratio = rnorm/ynorm
      xratio=dxnorm/xnorm
c
c***
      if ( xratio .le. e3.and.ratio.le.e3 ) go to 200
c
ccc compute beta(i)
c
      call ldlt(n,m,ef,df,cf,bf,af,bf(2),cf(n),df(n+1),ef(n+2),
     1          up,r)
      alpha = ddot(mn,r,1,up,1)
      beta = alpha/rdot
      rdot = alpha
c get new p (stored in aux)
      do 95  i11=1,mn
 95     aux(i11)=beta*aux(i11)+up(i11)
  100 continue
C      write (59,102 ) xratio,ratio
      write (6,102 ) xratio,ratio
 102  format (1h0,'cg terminates because of too many iterations. ', /,
     x         'xratio = ',e12.4, 'ratio = ',e12.4)
  200 istop=i

      return
      end

 
      subroutine convect
c
c ======================================================================
c
c   Purpose -
c     compute the convective transport terms needed for
c     changing the velocities from a Lagrangian to an
c     Eulerian frame
c
c   CONVECT is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          VTILDE
c
c
c   CONVECT calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c              BC    ripple     DVCAL    ripple
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
      dimension uc(nxy),vc(nxy)
c
      data tiny /1.0d-25/, zero /0.0d0/, one /1.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... save the Lagrangian velocities
c
      do 5 ij=1,nxy
        uc(ij)=un(ij)
        vc(ij)=vn(ij)
    5 continue
c
c.... compute the convective transport terms and
c     transform the Lagrangian velocities to Eulerian
c
      do 20 j=2,jm1
        m=1
        n=1
        do 20 i=2,im1
c
          ij=(j-1)*imax+i
          ipj=ij+1
          imj=ij-1
          ijp=ij+imax
          ijm=ij-imax
          ipjp=ipj+imax
          ipjm=ipj-imax
          imjp=imj+imax
          imjm=imj-imax
c
c....     get the velocity fluxes
c         -----------------
          call dvcal(uc,vc,i,j)
c         -----------------
c
          rhoij=f(ij)*rhof
          rhoipj=f(ipj)*rhof
          rhoijp=f(ijp)*rhof
          rhorc=(delx(i+1)*rhoij+delx(i)*rhoipj)/(delx(i)+delx(i+1))
          rhotc=(dely(j+1)*rhoij+dely(j)*rhoijp)/(dely(j)+dely(j+1))
	  
C          aut=(delx(i+1)*at(ij)+delx(i)*at(ipj))/(delx(i)+delx(i+1))
C          aub=(delx(i+1)*at(ijm)+delx(i)*at(ipjm))/(delx(i)+delx(i+1))
C          avr=(dely(j+1)*ar(ij)+dely(j)*ar(ijp))/(dely(j)+dely(j+1))
C          avl=(dely(j+1)*ar(imj)+dely(j)*ar(imjp))/(dely(j)+dely(j+1))
c
C          dudr=cvmgt(ac(ipj),1.0d0,ac(ipj).ge.em6)*dudr
C          dudl=cvmgt(ac(ij),1.0d0,ac(ij).ge.em6)*dudl

C          dudt=aut*dudt
C          dudb=aub*dudb

C          dvdt=cvmgt(ac(ijp),1.0d0,ac(ijp).ge.em6)*dvdt
C          dvdb=cvmgt(ac(ij),1.0d0,ac(ij).ge.em6)*dvdb

C          dvdr=avr*dvdr
C          dvdl=avl*dvdl

54	  continue
c
c....     reset x-velocity and check for surrounding fluid
          rdelx=1.0/(delx(i)+delx(i+1))
          rdely=1.0/(dely(j)+dely(j+1))
          if (ar(ij).lt.em6) go to 10
c
          if (i.eq.im1.and.kr.eq.5) goto 15

c....     compute the x-direction convective flux
          sgu=sign(one,uc(ij))
          rdxa=delx(i)+delx(i+1)+alpha*sgu*(delx(i+1)-delx(i))
          rdxa=1.0/rdxa
          fux=rdxa*uc(ij)*(delx(i)*dudr+delx(i+1)*dudl+
     &         alpha*sgu*(delx(i+1)*dudl-delx(i)*dudr))
          vbt=(delx(i)*vc(ipj)+delx(i+1)*vc(ij))*rdelx
          vbb=(delx(i)*vc(ipjm)+delx(i+1)*vc(ijm))*rdelx
          vav=0.5*(vbt+vbb)
          dyt=0.5*(dely(j)+dely(j+1))
          dyb=0.5*(dely(j-1)+dely(j))
          sgv=sign(one,vav)
          dya=dyt+dyb+alpha*sgv*(dyt-dyb)
          fuy=(vav/dya)*(dyb*dudt+dyt*dudb+alpha*sgv*
     &         (dyt*dudb-dyb*dudt))
c
	  if (npor.ne.0.and.(npc(ij).ne.1.or.npc(ipj).ne.1)) then
	    if (npc(ij).ne.1) then
		xpstmp=porousp(ij)
		gctmp=porousc(ij)
	    else
                xpstmp=porousp(ipj)
                gctmp=porousc(ipj)
	    endif
            u(ij)=u(ij)-fact*delt*(fux+fuy)/xpstmp**2/gctmp
	  else
            u(ij)=u(ij)-fact*delt*(fux+fuy)
	  endif

   10     u(ij)=cvmgt(zero,u(ij),((rhorc.lt.frsurf).and.
     &          ((nf(ij).gt.5).and.(nf(ipj).gt.5))).or.
     &		(ar(ij).lt.em6))

c
c....     reset y-velocity and check for surrounding fluid
          if (at(ij).lt.em6) go to 25
c
15	  continue
          if (j.eq.jm1.and.kt.eq.5) goto 20

c....     compute the y-direction convective flux
          ubr=(dely(j+1)*uc(ij)+dely(j)*uc(ijp))*rdely
          ubl=(dely(j+1)*uc(imj)+dely(j)*uc(imjp))*rdely
          uav=0.5*(ubr+ubl)
          dxr=0.5*(delx(i)+delx(i+1))
          dxl=0.5*(delx(i)+delx(i-1))
          sgu=sign(one,uav)
          dxa=dxr+dxl+alpha*sgu*(dxr-dxl)
          fvx=(uav/dxa)*(dxl*dvdr+dxr*dvdl+alpha*sgu*
     &         (dxr*dvdl-dxl*dvdr))
          sgv=sign(one,vc(ij))
          dya=dely(j+1)+dely(j)+alpha*sgv*(dely(j+1)-dely(j))
          fvy=(vc(ij)/dya)*(dely(j)*dvdt+dely(j+1)*dvdb+
     &          alpha*sgv*(dely(j+1)*dvdb-dely(j)*dvdt))
c
	  if (npor.ne.0.and.(npc(ij).ne.1.or.npc(ijp).ne.1)) then
            if (npc(ij).ne.1) then
                xpstmp=porousp(ij)
                gctmp=porousc(ij)
            else
                xpstmp=porousp(ijp)
                gctmp=porousc(ijp)
            endif
            v(ij)=v(ij)-fact*delt*(fvx+fvy)/xpstmp**2/gctmp
	  else
            v(ij)=v(ij)-fact*delt*(fvx+fvy)
	  endif

   25     v(ij)=cvmgt(zero,v(ij),((rhotc.lt.frsurf).and.
     &          ((nf(ij).gt.5).and.(nf(ijp).gt.5))).or.
     &		(at(ij).lt.em6))
   20 continue
c
c.... update the boundary conditions
c
C      ibcflg=1
	ibcflg0=1
c     ---------
   45 call bc
c     ---------
C      ibcflg=0
	ibcflg0=0
c
      return
      end


      real*8 function cvmgt(x1,x2,x3)
c
c ======================================================================
c
c   purpose -
c     Cray vectorization function:
c     Return x1 if x3 is true, otherwise return x2
c
c   cvmgt is called by -
c
c   cvmgt calls the following subroutines and functions -
c ======================================================================
c
      implicit real*8 (a-h,o-z)
      logical x3
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      if (x3) then
		cvmgt=x1
      else
		cvmgt=x2
      endif
c
      return
      end
  

C**********************************************************************
C**                  FUNCTION DDOT                                   **
C**********************************************************************
      real*8 function ddot(n,dx,incx,dy,incy)
      real*8 dx(1),dy(1)
      ddot = 0.d0
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
C
C.....CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
         ddot = ddot + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
c----------------------------------------------------------------DEBUG
c     ddot = 9.87654321
c----------------------------------------------------------------DEBUG
      return
C
C.....CODE FOR BOTH INCREMENTS EQUAL TO 1
C     CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
         ddot = ddot + dx(i)*dy(i)
   30 continue
c----------------------------------------------------------------DEBUG
c     ddot = 9.87654321
c----------------------------------------------------------------DEBUG
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
         ddot = ddot + dx(i)*dy(i) + dx(i+1)*dy(i+1) +
     1   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
c----------------------------------------------------------------DEBUG
c     ddot = 9.87654321
c----------------------------------------------------------------DEBUG
      return
C
C.....CODE FOR POSITIVE EQUAL INCREMENTS .NE.1
   60 continue
      ns = n*incx
          do 70 i=1,ns,incx
          ddot = ddot + dx(i)*dy(i)
   70     continue
c----------------------------------------------------------------DEBUG
c     ddot = 9.87654321
c----------------------------------------------------------------DEBUG
      return
C2345678901234567890123456789012345678901234567890123456789012345678901
C**********************************************************************
C**  Modification history: Created by C.L. Lawson (jpl), R.J. Hanson **
C**                        (snla), D.R. Kinkaid (u. of texas), F.T.  **
C**                        Krogh (jpl) OCT79. Revised DEC86          **
C**                        Adapted from NASA-VOF2D by MCW, JAN90     **
C**                                                                  **
C**   Purpose:  d.p. inner product of d.p. vectors. See NASA-VOF2D   **
C**             for further details.                                 **
C**                                                                  **
C**   Called from:         CGITJ,ILUCGJ                              **
C**   Calls to:                                                      **
C**   External fctns:                                                **
C**                                                                  **
C**********************************************************************
      end
  
 
      subroutine deltadj
c
c ======================================================================
c
c   Purpose -
c     time step adjustment
c
c   DELTADJ is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          NEWCYC 
c
c
c   DELTADJ calls the following subroutines and functions -
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
      data itmin /20/
      data itcjr / 100 /
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      itmost=itmxiccg
      deltn=delt
      itc=0
      jtc=0
      if (flgc.lt.0.5) go to 20
      if (autot.eq.0.0) go to 20
C.....modified to skip the rewind step, which is incompatible with moving body
	delt=0.8*delt
	goto 20
C.....end of modification


      t=t-delt
      ncyc=ncyc-1
      it=min(0,it-10)
      delt=0.8*delt
c
      do 11 i=1,imax
        do 10 j=1,jmax
          ij=(j-1)*imax+i
          xdis=xdis-delt*xep(ij)
     &    *(x(i)-x(i-1))*(y(j)-y(j-1))*f(ij)*cvmgt(porousp(ij),
     &	  ac(ij),npc(ij).ne.1.and.npor.ne.0)
          p(ij)=pn(ij)
          f(ij)=fn(ij)
	    nf(ij)=nfold(ij)
          u(ij)=0.0
          v(ij)=0.0
	    xp(ij)=xpn(ij)
	    xk(ij)=xkn(ij)
		xep(ij)=xepn(ij)
		xnut(ij)=xnutn(ij)
		xnuty(ij)=xnutyn(ij)
   10   continue
   11 continue
c
      flgc=0.0
      nflgc=nflgc+1
   20 continue
      dumx=em10
      dvmx=em10
C.....added to account xnut(ij)
      dtvis=1.0d+10
      if (autot.eq.0.0) go to 40
c
      ijx=0
      ijy=0
      do 31 i=1,im1
        do 30 j=2,jm1
          ij=(j-1)*imax+i
          udm=abs(un(ij))*ar(ij)/(xi(i+1)-xi(i))
          vdm=abs(vn(ij))*at(ij)/(yj(j+1)-yj(j))
C.........For small AC cell, it is equivalent to say that it has small DX & DY
C.........Therefore, dt needs be reduced accordingly; Here, an approximate method
C.........is used to estimate AC effects so that it will maintain a balance between
C.........the computational efficiency and stability of VOF & velocity computation
C.........e.g., if obstacle is horizontal & vertical, ac effects not considered
C.........on slope, we need to carefully reduce dt to attain stable solution
		udm=min(abs(un(ij))/(xi(i+1)-xi(i)),udm*cvmgt(ac(ij)/
     &		(ac(ij+1)+em6),ac(ij+1)/(ac(ij)+em6),ac(ij).ge.ac(ij+1)))
		vdm=min(abs(vn(ij))/(yj(j+1)-yj(j)),vdm*cvmgt(ac(ij)/
     &		(ac(ij+imax)+em6),ac(ij+imax)/(ac(ij)+em6),
     &		ac(ij).ge.ac(ij+imax)))
          dumx=max(dumx,udm)
          dvmx=max(dvmx,vdm)
	    if (dumx.eq.udm) ijx = ij
          if (dvmx.eq.vdm) ijy = ij
C.........added to account xnut(ij)
          dxsq=delx(i)**2
          dysq=dely(j)**2
          rdsq=dxsq*dysq/(dxsq+dysq)
          rdsq=rdsq/(3.0*(xnu+xnut(ij))+1.0d-60)
          dtvis=min(dtvis,rdsq)
	    if (dtvis.eq.rdsq) ijxnut=ij
   30   continue
   31 continue
c
      jx=ijx/imax + 1
      ix=ijx-(jx-1)*imax
      jy=ijy/imax + 1
      iy=ijy-(jy-1)*imax
C.....added to account xnut(ij)
      jvis=ijxnut/imax + 1
      ivis=ijxnut-(jvis-1)*imax
c
      dtmp=1.025d0
      if(iter.lt.itmin) dtmp=1.05d0
      if(iter.gt.itmost.and.liter.gt.itmost) dtmp=0.99d0
      delto=delt*dtmp
      delto=min(delto,dtmax)
      delt=min(delto,con/dumx,con/dvmx,dtvis,crest)

      if (delt.eq.con/dumx) then
        itc=ix
        jtc=jx
      endif
      if (delt.eq.con/dvmx) then
        itc=iy
        jtc=jy

      endif
      if (delt.eq.dtvis) then
        itc=ivis
        jtc=jvis
      endif
   40 if (delt.eq.deltn) go to 60
c
   60 continue
      liter=iter
c.... time step growing
      if ((delt.eq.delto).and.(dtmp.gt.1.0d0)) idt="g"
c.... time step decaying
      if ((delt.eq.delto).and.(dtmp.lt.1.0d0)) idt="d"
c.... time step unchanged
      if ((delt.eq.deltn).or.(dtmp.eq.1.0d0)) idt="f"
c.... x Courant limit
      if (delt.eq.con/dumx) idt="cx"
c.... y Courant limit
      if (delt.eq.con/dvmx) idt="cy"
c.... viscous limit
      if (delt.eq.dtvis) idt="v"
c.... maximum allowed time step
      if (delt.eq.dtmax) idt="m"
c
C.... limited by negative production
C      if (delt.eq.1.0d0/xxmax) idt="ps"
C.... limited by negative production
      if (delt.eq.crest) idt="np"

C.....add to stop the program if delt is too small
      if (delt.lt.1.0e-4*dtmax) then
        write(13,*)'delt',delt, 'is less than 1.0e-4*dtmax',
     &          dtmax*1.0e-4
        write(iotty,*)'delt',delt, 'is less than 1.0e-4*dtmax', 
     &		dtmax*1.0e-4
        nexit=9999 
      endif

      return
      end
 

C**********************************************************************
C**                 FUNCTION DNRM2                                   **
C**********************************************************************
      real*8 function dnrm2(n,dx,incx)
      save cutlo, cuthi, zero, one
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      integer next
      real*8 dx(1), cutlo, cuthi, hitest, sum, xmax, zero, one
      data zero, one /0.0d0, 1.0d0/
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
      if(n .gt. 0) go to 10
      dnrm2 = zero
      go to 300
   10 assign 30 to next
      sum = zero
      nn = n * incx
C
C.....BEGIN MAIN LOOP
      i = 1
   20 go to next,(30, 50, 70, 110)
   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
C
C.....PHASE 1:  SUM IS ZERO
   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85
C
C.....PREPARE FOR PHASE 2
      assign 70 to next
      go to 105
C
C.....PREPARE FOR PHASE 4
  100 i = j
      assign 110 to next
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115
C
C.....PHASE 2:  SUM IS SMALL
C     SCALE TO AVOID DESTRUCTIVE UNDERFLOW
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
C
C.....COMMON CODE FOR PHASES 2 AND 4
C     IN PHASE 4 SUM IS LARGE SO SCALE TO AVOID OVERFLOW
  110 if( dabs(dx(i)) .le. xmax ) go to 115
      sum = one + sum * (xmax / dx(i))**2
      xmax = dabs(dx(i))
      go to 200
C
  115 sum = sum + (dx(i)/xmax)**2
      go to 200
C
C.....PREPARE FOR PHASE 3
   75 sum = (sum * xmax) * xmax
C
C.....FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
   85 hitest = cuthi/float( n )
C
C.....PHASE 3:  SUM IS MID-RANGE THEREFORE NO SCALING
      do 95 j =i,nn,incx
      if(dabs(dx(j)) .ge. hitest) go to 100
   95 sum = sum + dx(j)**2
      dnrm2 = dsqrt( sum )
      go to 300
C
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
C
C.....END OF MAIN LOOP
C     COMPUTE SQUARE ROOT AND ADJUST FOR SCALING
      dnrm2 = xmax * dsqrt(sum)
  300 continue
      return
C2345678901234567890123456789012345678901234567890123456789012345678901
C**********************************************************************
C**  Modification history: Created by C.L. Lawson (jpl), R.J. Hanson **
C**                        (snla), D.R. Kincaid (u. of texas), F.T.  **
C**                        Krogh (jpl) OCT79.  Revised DEC86         **
C**                        Adapted from NASA-VOF2D by MCW, JAN90     **
C**                                                                  **
C**    Purpose: Euclidean length (l2 norm) of d.p. vector.  See      **
C**             NASA-VOF2D for further details.                      **
C**                                                                  **
C**    Called from:        ILUCGJ                                    **
C**    Calls to:                                                     **
C**    External fctns:                                               **
C**                                                                  **
C**********************************************************************
      end
  
 
      subroutine dvcal (uc,vc,i5,j5)
c
c ======================================================================
c
c   Purpose -
c     compute the x- and y-velocity fluxes for cell ij
c
c   DVCAL is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c         CONVECT
c
c
c   DVCAL calls the following subroutines and functions -
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
c.... provide free-slip-like (islip=0) or standard free-slip (islip=1)
c.....or no-slip (islip=2)
c.....boundary conditions for all obstacle surface flow cells
c
C      data islip / 0 /
      dimension uc(1),vc(1)
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c                                dvdt
c                                 x
c
c
c
c                               vc(ij)
c                   dvdl x--------*--------x dudt,dvdr
c                        |                 |
c                        |                 |
c                        |                 |
c                        |    dudl,dvdb    |      dudr
c                 uc(imj)*        x        *uc(ij)  x
c                        |                 |
c                        |                 |
c                        |                 |
c                        |                 |
c                        x--------*--------x dudb
c                               vc(ijm)
c
      ij=(j-1)*imax+i
      ipj=ij+1
      ijp=ij+imax
      ipjp=ipj+imax
      ipjm=ipj-imax
      ijm=ij-imax
      imj=ij-1
      imjp=imj+imax
      imjm=imj-imax
c
C.....add the f(ij) as the weighting factor to reduce surface stress (2004/9/9)
C      dudr=((uc(ipj)-uxmb(nmovbd(ipj)))*ar(ipj)
C     &		-(uc(ij)-uxmb(nmovbd(ipj)))*ar(ij))*rdx(i+1)
C     &		*f(ipj)*ac(ipj)
C      dudl=((uc(ij)-uxmb(nmovbd(ij)))*ar(ij)
C     &		-(uc(imj)-uxmb(nmovbd(ij)))*ar(imj))*rdx(i)
C     &		*f(ij)*ac(ij)
C      dudt=((uc(ijp)-uxmb(nmovbd(ij)))*ar(ijp)
C     &		-(uc(ij)-uxmb(nmovbd(ij)))*ar(ij))*2.0
C     &		/(dely(j)+dely(j+1))*(f(ij)*ac(ij)+f(ipj)*ac(ipj)
C     &		+f(ijp)*ac(ijp)+f(ipjp)*ac(ipjp))/4.0d0
C      dudb=((uc(ij)-uxmb(nmovbd(ij)))*ar(ij)
C     &		-(uc(ijm)-uxmb(nmovbd(ij)))*ar(ijm))*2.0
C     &		/(dely(j)+dely(j-1))*(f(ij)*ac(ij)+f(ipj)*ac(ipj)
C     &		+f(ijm)*ac(ijm)+f(ipjm)*ac(ipjm))/4.0d0
C      dvdr=((vc(ipj)-vymb(nmovbd(ij)))*at(ipj)
C     &		-(vc(ij)-vymb(nmovbd(ij)))*at(ij))*2.0
C     &		/(delx(i)+delx(i+1))*(f(ij)*ac(ij)+f(ipj)*ac(ipj)
C     &		+f(ijp)*ac(ijp)+f(ipjp)*ac(ipjp))/4.0d0
C      dvdl=((vc(ij)-vymb(nmovbd(ij)))*at(ij)
C     &		-(vc(imj)-vymb(nmovbd(ij)))*at(imj))*2.0
C     &		/(delx(i)+delx(i-1))*(f(ij)*ac(ij)+f(imj)*ac(imj)
C     &		+f(ijp)*ac(ijp)+f(imjp)*ac(imjp))/4.0d0
C      dvdt=((vc(ijp)-vymb(nmovbd(ijp)))*at(ijp)
C     &		-(vc(ij)-vymb(nmovbd(ijp)))*at(ij))*rdy(j+1)
C     &		*f(ijp)*ac(ijp)
C      dvdb=((vc(ij)-vymb(nmovbd(ij)))*at(ij)
C     &		-(vc(ijm)-vymb(nmovbd(ij)))*at(ijm))*rdy(j)
C     &		*f(ij)*ac(ij)
c
C.....It is found that when f(ij),ar(ij),at(ij), ac(ij) effects are included in
C.....the calculation of dvcal, it always gives earlier collapse of breaking wave
C.....front & overestimation of runup (e.g., front moves faster); For this reason,
C.....simple formula based on velocity only is used to evaluate dvcal 
      dudr=(uc(ipj)-uc(ij))*rdx(i+1)
      dudl=(uc(ij)-uc(imj))*rdx(i)
      dudt=(uc(ijp)-uc(ij))*2.0/(dely(j)+dely(j+1))
      dudb=(uc(ij)-uc(ijm))*2.0/(dely(j)+dely(j-1))
      dvdr=(vc(ipj)-vc(ij))*2.0/(delx(i)+delx(i+1))
      dvdl=(vc(ij)-vc(imj))*2.0/(delx(i)+delx(i-1))
      dvdt=(vc(ijp)-vc(ij))*rdy(j+1)
      dvdb=(vc(ij)-vc(ijm))*rdy(j)

C=======================================================================
C.....eliminate the velocity gradient at the free surface
      if (f(ipj).le.emf) dudr=0.0
      if (f(ij).le.emf) dudl=0.0
C.....modified to give the zero tangential gradient on free surface
      if ((f(ijp).le.emf.and.f(ipjp).le.emf)) dudt=0.0
      if ((f(ijm).le.emf.and.f(ipjm).le.emf)) dudb=0.0
      if ((f(ipj).le.emf.and.f(ipjp).le.emf)) dvdr=0.0
      if ((f(imj).le.emf.and.f(imjp).le.emf)) dvdl=0.0
C	if (nf(ij).eq.6.or.nf(ipj).eq.6.or.nf(ijp).eq.6.or.nf(ipjp).eq.6)
C     &	then
C		dudt=0.0
C		dvdr=0.0
C	endif
C	if (nf(ij).eq.6.or.nf(ipj).eq.6.or.nf(ijm).eq.6.or.nf(ipjm).eq.6)
C     &	dudb=0.0
C	if (nf(ij).eq.6.or.nf(imj).eq.6.or.nf(ijp).eq.6.or.nf(imjp).eq.6)
C     &	dvdl=0.0
      if (f(ijp).le.emf) dvdt=0.0
      if (f(ij).le.emf) dvdb=0.0
c
C=======================================================================

      if (islip.ne.0) go to 20
C.....upwind concept is applied to decide whether the obstacle information should be included
      if (ar(ipj).lt.em6.and.u(ij).gt.0.0) dudr=0.0
      if (ar(imj).lt.em6.and.u(ij).lt.0.0) dudl=0.0
      if (ar(ijp).lt.em6) dudt=0.0
      if (ar(ijm).lt.em6) dudb=0.0
      if (at(ipj).lt.em6) dvdr=0.0
      if (at(imj).lt.em6) dvdl=0.0
      if (at(ijp).lt.em6.and.v(ij).gt.0.0) dvdt=0.0
      if (at(ijm).lt.em6.and.v(ij).lt.0.0) dvdb=0.0

      go to 9999

   20 if (islip.ne.1) go to 9999
      if (ar(ijp).lt.em6) dudt=0.0
      if (ar(ijm).lt.em6) dudb=0.0
      if (at(ipj).lt.em6) dvdr=0.0
      if (at(imj).lt.em6) dvdl=0.0
c
 9999 return
      end

 
      function fxy(x,y,a1,a2,b1,b2,c1,c2,d1,
     &                              kxc,d2,kxs,e1,kyc,e2,kys)
c
c ======================================================================
c
c   Purpose -
c     solve the function f(x,y) for a point (x,y), where:
c     f(x,y) = a1*x + a2*x**2 + b1*y + b2*y**2 + c1 + c2*x*y
c            +d1*cos(kxc*x)+d2*sin(kxs*x)+e1*cos(kyc*y)+e2*sin(kys*y)
c
c   FXY is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c            ASET  INITVOFF
c
c
c   FXY calls the following subroutines and functions -
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
      real*8 kxc,kxs,kyc,kys
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      fxy=a2*x*x+a1*x+b2*y*y+b1*y+c2*x*y+c1
     &      +d1*cos(kxc*x)+d2*sin(kxs*x)
     &      +e1*cos(kyc*y)+e2*sin(kys*y)
c
      return
      end

 
      subroutine geom
c
c ======================================================================
c
c   Purpose -
c     calculate the metric for the mesh
c
c   GEOM is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c         MESHSET
c
c
c   GEOM calls the following subroutines and functions -
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
      real*8 jacob
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... zero-out the volume arrays
c
c     ---------------------------------
      call setarry (cvol(1),0.0d0,nxy)
c     ---------------------------------
c
c.... Compute the metric coefficients used
c     in the implicit pressure solution
c
      do 100 j=2,jm1
        js=j-1
        do 100 i=2,im1
          is=i-1
          ij=(j-1)*imax+i
          x1=0.5*(x(is+1)+x(is+1)-x(is)-x(is))
          x2=0.5*(x(is+1)+x(is)-x(is)-x(is+1))
          y1=0.5*(y(js)+y(js+1)-y(js)-y(js+1))
          y2=0.5*(y(js+1)+y(js+1)-y(js)-y(js))
          jacob=x1*y2-x2*y1
c
          jcb(ij)=jacob
          alp(ij)=(x2**2 + y2**2)/jacob
          gam(ij)=(x1**2 + y1**2)/jacob
          cvol(ij)=ri(i)*delx(i)*dely(j)
c
  100 continue
c
      jb=1
      jt=jmax
      ijb=2
      ijt=jm1*imax+2
      do 101 i=2,im1
        alp(ijt)=alp(ijt-imax)
        gam(ijt)=gam(ijt-imax)
        jcb(ijt)=jcb(ijt-imax)
        cvol(ijt)=cvol(ijt-imax)
        alp(ijb)=alp(ijb+imax)
        gam(ijb)=gam(ijb+imax)
        jcb(ijb)=jcb(ijb+imax)
        cvol(ijb)=cvol(ijb+imax)
        ijb=ijb+1
        ijt=ijt+1
  101 continue
c
      il=1
      ir=imax
      ijl=1
      ijr=imax
      do 102 j=1,jmax
        alp(ijr)=alp(ijr-1)
        gam(ijr)=gam(ijr-1)
        jcb(ijr)=jcb(ijr-1)
        cvol(ijr)=cvol(ijr-1)
        alp(ijl)=alp(ijl+1)
        gam(ijl)=gam(ijl+1)
        jcb(ijl)=jcb(ijl+1)
        cvol(ijl)=cvol(ijl+1)
        ijl=ijl+imax
        ijr=ijr+imax
  102 continue
c
      return
      end


      subroutine ilucgj(e,d,c,b,a,bu,cu,du,eu,
     1      el,dl,cl,bl,af,bf,cf,df,ef,
     2      n,m,x,y,r,up,aux,istop,ercg,
     3      xratio,ratio,xnorm,ynorm)
c
c***purpose
c
c      given a 9-point matrix multiply routine called matmul9
c      and an approximate factorization a = l*d*u, this routine
c      provides a conjugate gradient acceleration to iterate
c      to an approximate solution to a*x=y.
c
c***description
c
c       ilucgj requires that the coefficients of the ith
c       equation be stored in the ith element of the diagonal
c       arrays. nine diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c       the ldu factors are stored similarly in nine
c       diagonal arrays.
c
c     on entry
c
c        e,d,c,b,a,bu,cu,du,eu,   are 1-d arrays of length
c            at least nxm. the arrays contain, from left to right,
c            the diagonals of a as one sweeps a left to right.
c
c        el,dl,cl,bl,af,bf,cf,df,ef   are 1-d arrays of length
c            at least nxm. the arrays contain from left to right
c            the diagonals of l,d,u as they are swept left to right.
c
c        n     is the size  or order of each tridiagonal block if
c              a is viewed as a block tridiagonal system of equations.
c
c        m     is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c
c        x     is an estimate to the solution vector of size nxm.
c
c        y     is an r.h.s. vector of dimenison at least nxm.
c
c        r,up,aux  are work variables at least nxm in length.
c               the final residual is left array r.
c
c        istop  is the maximum number of iterations allowed.
c               upon exit it contains the actual number of
c               iterations required to converge to accuracy
c               specified.
c
c        ercg  is the error tolerance for convegence
c
c     on return
c
c        x     is the output vector of dimension at least
c              nxm containing the solution vector.
c
c        xratio   is the ratio of dxnorm to xnorm.
c
c        ratio    is the ratio of rnorm to ynorm
c
c        xnorm    is the l2 norm of the solution vector x.
c
c        ynorm    is the l2 norm of the r.h.s.
c
c***routines called   lu9p matmul9 ldlt utdu rypax ymax9p
c
      implicit real*8 (a-h,o-z)
      real*8 e(m*n),d(m*n),c(m*n),b(m*n),a(m*n),bu(m*n),
     1	     cu(m*n),du(m*n),eu(m*n),el(m*n),dl(m*n),cl(m*n),bl(m*n),
     2       af(m*n),bf(m*n),cf(m*n),df(m*n),ef(m*n),
     3       x(m*n),y(m*n),r(m*n),up(m*n),aux(m*n)
      mn = m*n
      e3=ercg
c
ccc copy a into factor storage b
c
      do 10 i=1,mn
        el(i) = e(i)
        dl(i) = d(i)
        cl(i) = c(i)
        bl(i) = b(i)
        af(i) = a(i)
        bf(i) = bu(i)
        cf(i) = cu(i)
        df(i) = du(i)
        ef(i) = eu(i)
   10 continue
c
ccc factor b into incomplete ldu
c
      call lu9p(n,m,el,dl,cl,bl,af,bf,cf,df,ef)
c initialize the residual r
      ynorm = dnrm2(mn,y,1)
c
ccc compute r0
c
      call matmul9(e(1),d(1),c(1),b(1),a(1),
     1             bu(1),cu(1),du(1),eu(1),n,m,x,r)
      do 30 i=1,mn
        r(i) = y(i)-r(i)
   30 continue
      if (dnrm2(mn,r,1) .eq. 0.)     return
c
ccc solve l * l(trans)*aux = r0
c
      call slve9(n,m,el(1),dl(1),cl(1),bl(1),af(1),
     1         bl(2),cl(n),dl(n+1),el(n+2),aux,r)
c
ccc dot r with ldlt[-1]r
c
      rdot = ddot(mn,r,1,aux,1)
c
ccc compute q(0)
c
      call rypax(eu(-n),du(1-n),cu(2-n),bu(0),a(1),
     1     b(2),c(n),d(n+1),e(n+2),n,m,aux,up,0.d0)
C      call rypax(eu(-n),du(1-n),cu(2-n),0.0d0,a(1),
C     1     b(2),c(n),d(n+1),e(n+2),n,m,aux,up,0.d0)
CC      call rypax(eu,du,cu,bu(0),a(1),
CC     1     b(2),c(n),d(n+1),e(n+2),n,m,aux,up,0.)

c
ccc compute p0
c
      call slve9(n,m,ef(-n),df(1-n),cf(2-n),bf(0),af(1),
     1             bf(1),cf(1),df(1),ef(1),aux,up)
C      call slve9(n,m,ef(-n),df(1-n),cf(2-n),0.0d0,af(1),
C     1             bf(1),cf(1),df(1),ef(1),aux,up)
CC      call slve9(n,m,ef,df,cf,bf(0),af(1),
CC     1             bf(1),cf(1),df(1),ef(1),aux,up)

c
ccc begin main loop
c
      do 100 i100 = 1,istop
c
ccc compute alpha
c
      alpha = ddot(mn,up,1,aux,1)
      alpha = rdot/alpha
c
ccc multiply aux by alpha
c
      do 70 i7 = 1,mn
 70   aux(i7) = alpha*aux(i7)
c
ccc compute new x
c
      do 80 i8 = 1,mn
 80   x(i8) = x(i8) + aux(i8)
c
ccc r = r-ax
c
      call ymax9p(e(1),d(1),c(1),b(1),a(1),
     1             bu(1),cu(1),du(1),eu(1),n,m,aux,r)
c
ccc termination test
ccc
ccc stop whenever norm(aux) / norm (x) .lt. e3
c
      rnorm = dnrm2(mn,r,1)
      dxnorm = dnrm2(mn,aux,1)
      xnorm = dnrm2(mn,x,1)
c***
ccc output ratio
c
      ratio = rnorm/ynorm
      xratio=dxnorm/xnorm
c
c***
      if ( xratio .le. e3.and.ratio.le.e3 ) go to 200
c
ccc compute beta
c
      call slve9(n,m,el(1),dl(1),cl(1),bl(1),af(1),
     1         bl(2),cl(n),dl(n+1),el(n+2),aux,r)
      alpha = ddot(mn,r,1,aux,1)
      beta = alpha/rdot
      rdot = alpha
c
ccc compute q
c
      call rypax(eu(-n),du(1-n),cu(2-n),bu(0),a(1),
     1     b(2),c(n),d(n+1),e(n+2),n,m,aux,up,beta)
C      call rypax(eu(-n),du(1-n),cu(2-n),0.0d0,a(1),
C     1     b(2),c(n),d(n+1),e(n+2),n,m,aux,up,beta)
CC      call rypax(eu,du,cu,bu(0),a(1),
CC     1     b(2),c(n),d(n+1),e(n+2),n,m,aux,up,beta)

c
ccc compute p
c
      call slve9(n,m,ef(-n),df(1-n),cf(2-n),bf(0),af(1),
     1             bf(1),cf(1),df(1),ef(1),aux,up)
C      call slve9(n,m,ef(-n),df(1-n),cf(2-n),0.0d0,af(1),
C     1             bf(1),cf(1),df(1),ef(1),aux,up)
CC      call slve9(n,m,ef,df,cf,bf(0),af(1),
CC     1             bf(1),cf(1),df(1),ef(1),aux,up)

  100 continue
      write (6,102 ) xratio,ratio
 102  format (1h0,'cg terminates because of too many iterations. ', /,
     x         'xratio = ',e12.4, 'ratio = ',e12.4)
  200 istop=i100
      return
      end


      subroutine implctp
c
c ======================================================================
c
c   Purpose -
c     calculate the implicit pressure
c
c   IMPLCTP is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE 
c
c
c   IMPLCTP calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c           ACCEL    ripple    ILUCGJ      iccg   BDYCELL    ripple
c           CGITJ      iccg    STRAIN    ripple
c
c ======================================================================
c
c##############################################################
      implicit real*8 (a-h,o-z)
      
c##############################################################
c
c###########
      include  "comdk2.h" 
       
c###########
c
      dimension ap(100),bp(100),bfp(100),dp(100),dfp(100),sp(100)
      data tiny /1.0d-25/, zerod /0.0d0/, oned /1.0d0/
      external ck,ce,cn
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... get the strain and rotation tensors
c
c     -----------------------------------------------------------------
      call strain(ar,at,ac,x,y,ri,r,u,v,im1,jm1,nxy,
     &            tauxx,tauyy,tauxy,nf,uxmb,vymb,nmovbd)
c     -----------------------------------------------------------------
c
c.... reflect pressures to ghost cells
c
c     --------------------------------------------------------
      call bdycell(im1,jm1,1,1,1,1,pbcl,pbcr,pbct,pbcb,pbc,p,
     &	xi,yj,gx,gy)
c     --------------------------------------------------------
c
c.... load the iccg arrays, which are the diagonals of
c     the matrix to be inverted.  The arrays represent
c     the coefficients for cell ij and each of its 8
c     neighbors, as shown below
c
c                       ----------------------------
c                       |        |        |        |
c                       |        |        |        |
c                       |  imjp  |  ijp   |  ipjp  |
c                       |        |        |        |
c                       |  (C)   |  (D)   |  (E)   |
c                       ----------------------------
c                       |        |        |        |
c                       |        |        |        |
c                       |  imj   |   ij   |  ipj   |
c                       |        |        |        |
c                       |  (BF)  |  (A)   |  (B)   |
c                       ----------------------------
c                       |        |        |        |
c                       |        |        |        |
c                       | imjm   |  ijm   |  ipjm  |
c                       |        |        |        |
c                       |  (EF)  |  (DF)  |  (CF)  |
c                       ----------------------------
c

      nm=1
      coef=4.*xmu/3.
	nirrgwvct=0
      do 4001 j=2,jm1
        do 4000 i=2,im1
          ij=(j-1)*imax+i
          ipj=ij+1
          ijp=ij+imax
          ipjp=ipj+imax
          imjp=ijp-1
          imj=ij-1
          imjm=imj-imax
          ijm=ij-imax
          ipjm=ipj-imax
          if ((i.eq.2).and.(cyl.eq.1.0d0)) imjm=ijm
          if ((i.eq.2).and.(cyl.eq.1.0d0)) imj=ij
          if ((i.eq.2).and.(cyl.eq.1.0d0)) imjp=ijp
          p(ij)=pn(ij)-psat
c
          rhoij=f(ij)*rhof
          rhoipj=f(ipj)*rhof
          rhoijp=f(ijp)*rhof
          rhoimj=f(imj)*rhof
          rhoijm=f(ijm)*rhof
          zero=cvmgt(zerod,oned,
     &                (rhoij.lt.em6).or.(ac(ij).lt.emf))
c
          rl=omcyl+cyl*0.5*(x(i-1)+x(i-1))
          rhol=(delx(i-1)*rhoij+delx(i)*rhoimj)/(delx(i)+delx(i-1))
          rrhol=1./(rhol+tiny)
          rr=omcyl+cyl*0.5*(x(i)+x(i))
          rhor=(delx(i+1)*rhoij+delx(i)*rhoipj)/(delx(i)+delx(i+1))
          rrhor=1./(rhor+tiny)
          rt=omcyl+cyl*0.5*(x(i-1)+x(i))
          rhot=(dely(j+1)*rhoij+dely(j)*rhoijp)/(dely(j)+dely(j+1))
          rrhot=1./(rhot+tiny)
          rb=omcyl+cyl*0.5*(x(i-1)+x(i))
          rhob=(dely(j-1)*rhoij+dely(j)*rhoijm)/(dely(j)+dely(j-1))
          rrhob=1./(rhob+tiny)
c
          if (npor.eq.0) then
            gctmpra=1.0
            gctmpr=1.0
            gctmpla=1.0
            gctmpl=1.0
          endif
          if (npor.ne.0) then
            if (ar(ij).gt.em6.and.(npc(ij).ne.1.or.npc(ipj).ne.1)) 
     &		then
                gctmpra=(porousc(ij)+porousc(ipj))/2.0
                if (npc(ij).ne.1) then
                  xpstmp=porousp(ij)
                  xatmp=porousa(ij)
                  xxbtmp=porousb(ij)
                  d50tmp=porousd(ij)
                  xkc=max(abs(un(ij)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict0=(xatmp+xbb
     &			*sqrt(un(ij)**2+vn(ij)**2))*abs(gy)
                else
                  xfrict0=0.0
                endif
                if (npc(ipj).ne.1) then
                  xpstmp=porousp(ipj)
                  xatmp=porousa(ipj)
                  xxbtmp=porousb(ipj)
                  d50tmp=porousd(ipj)
                  xkc=max(abs(un(ij)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict1=(xatmp+xbb
     &			*sqrt(un(ij)**2+vn(ipj)**2))*abs(gy)
                else
                  xfrict1=0.0
                endif
		xfriction=(xfrict0+xfrict1)/2.0
                gctmpr=1.0/(gctmpra+xfriction*delt)
	    else
                gctmpra=1.0
		gctmpr=1.0
	    endif
            if (ar(imj).gt.em6.and.(npc(imj).ne.1.or.npc(ij).ne.1)) 
     &		then
                gctmpla=(porousc(imj)+porousc(ij))/2.0
                if (npc(ij).ne.1) then
                  xpstmp=porousp(ij)
                  xatmp=porousa(ij)
                  xxbtmp=porousb(ij)
                  d50tmp=porousd(ij)
                  xkc=max(abs(un(imj)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict0=(xatmp+xbb
     &			*sqrt(un(imj)**2+vn(ij)**2))*abs(gy)
                else
                  xfrict0=0.0
                endif
                if (npc(imj).ne.1) then
                  xpstmp=porousp(imj)
                  xatmp=porousa(imj)
                  xxbtmp=porousb(imj)
                  d50tmp=porousd(imj)
                  xkc=max(abs(un(imj)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict1=(xatmp+xbb
     &			*sqrt(un(imj)**2+vn(imj)**2))*abs(gy)
                else
                  xfrict1=0.0
                endif
                xfriction=(xfrict0+xfrict1)/2.0
		gctmpl=1.0/(gctmpla+xfriction*delt)
	    else
                gctmpla=1.0
		gctmpl=1.0
	    endif
          endif
	    if (i.gt.2.and.i.lt.im1) then
            gctmprb=gctmpr
            gctmplb=gctmpl
          else
            gctmprb=0.0
            gctmplb=0.0
          endif

          alfr=gctmpr*ar(ij)*rr*rrhor*(2.0*dely(j)/(delx(i)+delx(i+1)))
          alfl=gctmpl*ar(imj)*rl*rrhol*(2.0*dely(j)/(delx(i)+delx(i-1)))

	  if (npor.eq.0) then
            gctmpta=1.0
            gctmpt=1.0
            gctmpba=1.0
            gctmpb=1.0
	  endif
          if (npor.ne.0) then
            if (at(ij).gt.em6.and.(npc(ij).ne.1.or.npc(ijp).ne.1)) 
     &		then
                gctmpta=(porousc(ij)+porousc(ijp))/2.0
		if (npc(ij).ne.1) then
                  xpstmp=porousp(ij)
                  xatmp=porousa(ij)
                  xxbtmp=porousb(ij)
                  d50tmp=porousd(ij)
                  xkc=max(abs(vn(ij)),0.01d0)*max(1.0d0,xxt)
     &          	/xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
		  xfrict0=(xatmp+xbb
     &			*sqrt(un(ij)**2+vn(ij)**2))*abs(gy)
		else
		  xfrict0=0.0
		endif
                if (npc(ijp).ne.1) then
                  xpstmp=porousp(ijp)
                  xatmp=porousa(ijp)
                  xxbtmp=porousb(ijp)
                  d50tmp=porousd(ijp)
                  xkc=max(abs(vn(ij)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict1=(xatmp+xbb
     &			*sqrt(un(ijp)**2+vn(ij)**2))*abs(gy)
                else
                  xfrict1=0.0
                endif
                xfriction=(xfrict0+xfrict1)/2.0
		gctmpt=1.0/(gctmpta+xfriction*delt)
	    else
                gctmpta=1.0
		gctmpt=1.0
            endif
            if (at(ijm).gt.em6.and.(npc(ijm).ne.1.or.npc(ij).ne.1)) 
     &		then
                gctmpba=(porousc(ijm)+porousc(ij))/2.0
                if (npc(ij).ne.1) then
                  xpstmp=porousp(ij)
                  xatmp=porousa(ij)
                  xxbtmp=porousb(ij)
                  d50tmp=porousd(ij)
                  xkc=max(abs(vn(ijm)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict0=(xatmp+xbb
     &			*sqrt(un(ij)**2+vn(ijm)**2))*abs(gy)
                else
                  xfrict0=0.0
                endif
                if (npc(ijm).ne.1) then
                  xpstmp=porousp(ijm)
                  xatmp=porousa(ijm)
                  xxbtmp=porousb(ijm)
                  d50tmp=porousd(ijm)
                  xkc=max(abs(vn(ijm)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict1=(xatmp+xbb
     &			*sqrt(un(ijm)**2+vn(ijm)**2))*abs(gy)
                else
                  xfrict1=0.0
                endif
                xfriction=(xfrict0+xfrict1)/2.0
		gctmpb=1.0/(gctmpba+xfriction*delt)
	    else
                gctmpba=1.0
	  	gctmpb=1.0
	    endif
          endif

	  if (j.gt.2.and.j.lt.jm1) then
            gctmptb=gctmpt
            gctmpbb=gctmpb
          else
            gctmptb=0.0
            gctmpbb=0.0
          endif

          gamt=gctmpt*at(ij)*rt*rrhot*(2.0*delx(i)/(dely(j)+dely(j+1)))
          gamb=gctmpb*at(ijm)*rb*rrhob*(2.0*delx(i)/(dely(j)+dely(j-1)))
c

C===========================================================================
C.........add srceold and redefined srce

C.........in the original RIPPLE, srce(nm) is obtained by multiplying 
C.........ac(ij) to each term; This is because the fact that when div(ij)
C.........is computed, is has been divided by ac(ij). Here the multiplying 
C.........of ac(ij) is just to compensate it. Therefor, it should not
C.........be carried into the gravity computation.
          divergence=((u(ij)-uxmb(nmovbd(ij)))*gctmpra*gctmpr*ar(ij)
     &		-(u(imj)-uxmb(nmovbd(ij)))*gctmpla*gctmpl*ar(imj))/delx(i)
     &          +((v(ij)-vymb(nmovbd(ij)))*gctmpta*gctmpt*at(ij)
     &		-(v(ijm)-vymb(nmovbd(ij)))*gctmpba*gctmpb*at(ijm))/dely(j)


          srce(nm)=-cvol(ij)*divergence/delt/fact
     &      -gx*(gctmprb*ar(ij)-gctmplb*ar(imj))*dely(j)
     &	    -gy*(gctmptb*at(ij)-gctmpbb*at(ijm))*delx(i)

C.........added source term
	  if (ninflow.eq.100) then
	    if (i.ge.isources.and.i.le.isourcee.and.
     &	    j.ge.jsources.and.j.le.jsourcee) then
	      if (nsource.eq.34) then
          	srce(nm)=srce(nm)+ac(ij)*cvol(ij)/delt/fact*
     &		ssource*sin(2.0*pi/tsource*t)
	      endif
            if (nsource.eq.44) then
			if (nirrgwvct.eq.0) then
			  sswave=0.0
C.........here is the place where random phase angle can be included if needed
	          do nw=1,nwave
                  sswave=sswave+swave(nw)*sin(2.0*pi/twave(nw)*t)
	          end do
			  nirrgwvct=1
			end if
              srce(nm)=srce(nm)+ac(ij)*cvol(ij)/delt/fact*sswave			 
            endif
            if (nsource.eq.4) then
              srce(nm)=srce(nm)+ac(ij)*cvol(ij)/delt/fact*
     &          ssource*(cos(pi/2-2.0*pi/tsource*t-cnf)+
     &			aa*xxk/8.0d0*cosh(xxk*h0)/sinh(xxk*h0)**3
     &			*(2.+cosh(2.*xxk*h0))
     &			*cos(2.*(pi/2.-2.0*pi/tsource*t-cnf)))
            endif
	      if (nsource.eq.14) then
              srce(nm)=srce(nm)+ac(ij)*cvol(ij)/delt/fact*
     &          ssource*(cos(pi/2-2.0*pi/tsource*t-cnf)
     &          +amp2/amp1*cos(2.*(pi/2.-2.0*pi/tsource*t-cnf))
     &          +amp3/amp1*cos(3.*(pi/2.-2.0*pi/tsource*t-cnf))
     &          +amp4/amp1*cos(4.*(pi/2.-2.0*pi/tsource*t-cnf))
     &          +amp5/amp1*cos(5.*(pi/2.-2.0*pi/tsource*t-cnf)))
	      endif 
            if (nsource.eq.5) then
              srce(nm)=srce(nm)+ac(ij)*cvol(ij)/delt/fact*
     &          ssource/cosh(sqrt(0.75*aa/h0**3)*(xstart-c1*t))**2
            endif
            if (nsource.eq.24) then
              srce(nm)=srce(nm)+ac(ij)*cvol(ij)/delt/fact*
     &          ssource*(ytrough/aa+cn(cnf*0.5*ck(mod1)+2.0*ck(mod1)*
     &          (-t/tsource),mod1)**2)
            endif
C...........Add overtopping mass back into the flume to ensure mass conservation
C	      if (novertop.eq.1) then
C          	area=(x(isourcee)-x(isources-1))*
C     &  			(y(jsourcee)-y(jsources-1))
C			srce(nm)=srce(nm)+ac(ij)*cvol(ij)/delt/fact*
C     &			2.0d0*overtopmass/delt/area
C	      endif
	    endif
	  endif

C.........end of modification
C===========================================================================

          betr=0.0
          betl=0.0
          bett=0.0
          betb=0.0
          a(nm)=(alfr+alfl+gamb+gamt)
          b(nm)=-(alfr+0.25*(betb-bett))
          bf(nm)=-(alfl+0.25*(bett-betb))
          c(nm)=-0.25*(betl+bett)
          cf(nm)=-0.25*(betr+betb)
          d(nm)=-(gamt+0.25*(betl-betr))
          df(nm)=-(gamb+0.25*(betr-betl))
          e(nm)=0.25*(betr+bett)
          ef(nm)=0.25*(betl+betb)
          soln(nm)=p(ij)
c
C.........modified to treat obstacle as porous media; zero=0 for obstacle or air but
C.........zero=1 for valid fluid cell
	  a(nm)=cvmgt(oned/tiny,a(nm),(rhoij.lt.em6).or.(ac(ij).lt.emf))
          b(nm)=zero*b(nm)
          bf(nm)=zero*bf(nm)
          c(nm)=zero*c(nm)
          cf(nm)=zero*cf(nm)
          d(nm)=zero*d(nm)
          df(nm)=zero*df(nm)
          e(nm)=zero*e(nm)
          ef(nm)=zero*ef(nm)
          srce(nm)=zero*srce(nm)
          soln(nm)=zero*soln(nm)
c
C.........apply pressure-driven wavemaker (from Armenio, 1998; Ocean Eng)
          if (ninflow.eq.54.and.xi(i).le.xxl/2.0d0.
     &          and.rhoij.lt.em6) then
            a(nm)=1.0d0
            srce(nm)=aa*(-gy)/2.0*(1.0+cos(2.0*xi(i)*pi/xxl))
     &          *sin(2.0*pi/xxt*t)
          endif 
c
          nm=nm+1
c
 4000   continue
 4001 continue

      do 400 j=3,jm1-1
        do 400 i=3,im1-1
c
          ij=(j-1)*imax+i
          if (ac(ij).lt.em6) go to 400
          nm=(j-2)*ibar+i-1
          ipj=ij+1
          ijp=ij+imax
          ijm=ij-imax
          imj=ij-1
c
          if (ac(ijp).lt.em6.and.j.lt.jm1) then
            a(nm)=a(nm)+d(nm)
            d(nm)=0.0
C            srce(nm)=srce(nm)-gy*delx(i)*at(ij)
          endif
          if (ac(ijm).lt.em6.and.j.gt.2) then
            a(nm)=a(nm)+df(nm)
            df(nm)=0.0
C	    srce(nm)=srce(nm)-gy*delx(i)*at(ijm)
          endif
          if (ac(imj).lt.em6.and.i.gt.2) then
            a(nm)=a(nm)+bf(nm)
            bf(nm)=0.0
C            srce(nm)=srce(nm)-gx*dely(j)*ar(imj)
          endif
          if (ac(ipj).lt.em6.and.i.lt.im1) then
            a(nm)=a(nm)+b(nm)
            b(nm)=0.0
C            srce(nm)=srce(nm)-gx*dely(j)*ar(ipj)
          endif
c
  400 continue
c
c.... Left boundary
c
      if (kl.eq.4) go to 115
      wsl=0.0
      if (kl.eq.5) wsl=1.0d0
      imj=imax+1
      imjm=imj-imax
      imjp=imj+imax
      nm=1
      do 110 j=2,jm1
        a(nm)=a(nm)+(1.-wsl)*bf(nm)
        d(nm)=d(nm)+(1-wsl)*c(nm)
        df(nm)=df(nm)+(1-wsl)*ef(nm)
        srce(nm)=srce(nm)
     &       -wsl*(c(nm)*p(imjp)+bf(nm)*p(imj)
     &            +ef(nm)*p(imjm))
        bf(nm)=0.0
        ef(nm)=0.0
        c(nm)=0.0
	  if (kl.eq.5) goto 2004
	  if (nf(imj+1).eq.6) goto 2004

          if (npor.ne.0) then
            if (npc(imj+1).ne.1) then
                gctmpla=(porousc(imj+1)+porousc(imj+2))/2.0
                if (npc(imj+1).ne.1) then
                  xpstmp=porousp(imj+1)
                  xatmp=porousa(imj+1)
                  xxbtmp=porousb(imj+1)
                  d50tmp=porousd(imj+1)
                  xkc=max(abs(un(imj+1)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
				xfrict0=(xatmp+xbb
     &			*sqrt(un(imj+1)**2+vn(imj+1)**2))*abs(gy)
                else
                  xfrict0=0.0
                endif
                if (npc(imj+2).ne.1) then
                  xpstmp=porousp(imj+2)
                  xatmp=porousa(imj+2)
                  xxbtmp=porousb(imj+2)
                  d50tmp=porousd(imj+2)
                  xkc=max(abs(un(imj+1)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict1=(xatmp+xbb
     &			*sqrt(un(imj+1)**2+vn(imj+2)**2))*abs(gy)
                else
                  xfrict1=0.0
                endif
                xfriction=(xfrict0+xfrict1)/2.0
                gctmpl=1.0/(gctmpla+xfriction*delt)
                srce(nm)=srce(nm)-gx*dely(j)*gctmpl*ar(imj+1)
            else
                srce(nm)=srce(nm)-gx*dely(j)*ar(imj+1)
            endif
          else
            srce(nm)=srce(nm)-gx*dely(j)*ar(imj+1)
          endif
        
2004	  continue

112     imj=imj+imax
        imjm=imjm+imax
        imjp=imjp+imax
        nm=nm+im1-1
  110 continue
c
c.... Top boundary
c
  115 continue
CCC	if (kt.eq.4) go to 125
      wst=0.0
      if (kt.eq.5) wst=1.0d0
      ij=jm1*imax+2
      imj=ij-1
      ipj=ij+1
      nm=(jm1-2)*(im1-1)+1
      do 120 i=2,im1
        a(nm)=a(nm)+(1.-wst)*d(nm)
        bf(nm)=bf(nm)+(1-wst)*c(nm)
        b(nm)=b(nm)+(1-wst)*e(nm)
        srce(nm)=srce(nm)-
     &           wst*(c(nm)*p(imj)+d(nm)*p(ij)
     &           +e(nm)*p(ipj))
        c(nm)=0.0
        d(nm)=0.0
        e(nm)=0.0
	  if (kt.eq.5) goto 2001
	  if (nf(ij-imax).eq.6) goto 2001

        if (npor.ne.0) then
            if (npc(ij-imax).ne.1) then
                gctmpta=(porousc(ij-2*imax)+porousc(ij-imax))/2.0
                if (npc(ij-imax).ne.1) then
                  xpstmp=porousp(ij-imax)
                  xatmp=porousa(ij-imax)
                  xxbtmp=porousb(ij-imax)
                  d50tmp=porousd(ij-imax)
                  xkc=max(abs(vn(ij-2*imax)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict0=(xatmp+xbb
     &		    *sqrt(un(ij-imax)**2+vn(ij-2*imax)**2))*abs(gy)
                else
                  xfrict0=0.0
                endif
                if (npc(ij-2*imax).ne.1) then
                  xpstmp=porousp(ij-2*imax)
                  xatmp=porousa(ij-2*imax)
                  xxbtmp=porousb(ij-2*imax)
                  d50tmp=porousd(ij-2*imax)
                  xkc=max(abs(vn(ij-2*imax)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict1=(xatmp+xbb
     &		    *sqrt(un(ij-2*imax)**2+vn(ij-2*imax)**2))*abs(gy)
                else
                  xfrict1=0.0
                endif
			  xfriction=(xfrict0+xfrict1)/2.0
                gctmpt=1.0/(gctmpta+xfriction*delt)
                srce(nm)=srce(nm)+gy*delx(i)*gctmpt*at(ij-2*imax)
            else
                srce(nm)=srce(nm)+gy*delx(i)*at(ij-2*imax)
            endif
        else
            srce(nm)=srce(nm)+gy*delx(i)*at(ij-2*imax)
        endif
2001	  continue

116     nm=nm+1
        imj=imj+1
        ij=ij+1
        ipj=ipj+1
  120 continue
c
c.... Right boundary
c
  125 continue
      if (kr.eq.4) go to 135
      wsr=0.0
      if (kr.eq.5) wsr=1.0d0
      nm=im1-1
      ij=2*imax
      ijp=ij+imax
      ijm=ij-imax
      do 130 j=2,jm1
        a(nm)=a(nm)+(1.-wsr)*b(nm)
        d(nm)=d(nm)+(1-wsr)*e(nm)
        df(nm)=df(nm)+(1-wsr)*cf(nm)
        srce(nm)=srce(nm)-wsr*(b(nm)*p(ij)+e(nm)*p(ijp)
     &           +cf(nm)*p(ijm))
        b(nm)=0.0
        e(nm)=0.0
        cf(nm)=0.0
	  if (kr.eq.5) goto 2003
	  if (nf(ij-1).eq.6) goto 2003

        if (npor.ne.0) then
            if (npc(ij-1).ne.1) then
                gctmpra=(porousc(ij-1)+porousc(ij-2))/2.0
                if (npc(ij-1).ne.1) then
                  xpstmp=porousp(ij-1)
                  xatmp=porousa(ij-1)
                  xxbtmp=porousb(ij-1)
                  d50tmp=porousd(ij-1)
                  xkc=max(abs(un(ij-2)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict0=(xatmp+xbb
     &			*sqrt(un(ij-2)**2+vn(ij-1)**2))*abs(gy)
                else
                  xfrict0=0.0
                endif
                if (npc(ij-2).ne.1) then
                  xpstmp=porousp(ij-2)
                  xatmp=porousa(ij-2)
                  xxbtmp=porousb(ij-2)
                  d50tmp=porousd(ij-2)
                  xkc=max(abs(un(ij-2)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict1=(xatmp+xbb
     &			*sqrt(un(ij-2)**2+vn(ij-2)**2))*abs(gy)
                else
                  xfrict1=0.0
                endif
		xfriction=(xfrict0+xfrict1)/2.0
                gctmpr=1.0/(gctmpra+xfriction*delt)
                srce(nm)=srce(nm)+gx*dely(j)*gctmpr*ar(ij-2)
            else
                srce(nm)=srce(nm)+gx*dely(j)*ar(ij-2)
            endif
        else
            srce(nm)=srce(nm)+gx*dely(j)*ar(ij-2)
        endif
2003	  continue

113     ij=ij+imax
        ijm=ijm+imax
        ijp=ijp+imax
        nm=nm+im1-1
  130 continue
c
c.... Bottom boundary
c
  135 continue
CCC	if (kb.eq.4) go to 145
      wsb=0.0
      if (kb.eq.5) wsb=1.0d0
      ijm=2
      imjm=1
      ipjm=3
      nm=1
      do 140 i=2,im1
        a(nm)=a(nm)+(1.-wsb)*df(nm)
        bf(nm)=bf(nm)+(1-wsb)*ef(nm)
        b(nm)=b(nm)+(1-wsb)*cf(nm)
        srce(nm)=srce(nm)-wsb*
     &           (ef(nm)*p(imjm)+df(nm)*p(ijm)
     &          +cf(nm)*p(ipjm))
        ef(nm)=0.0
        df(nm)=0.0
        cf(nm)=0.0
	  if (kb.eq.5) goto 2002
	  if (nf(ijm+imax).eq.6) goto 2002 

        if (npor.ne.0) then
            if (npc(ijm+imax).ne.1) then
                gctmpba=(porousc(ijm+2*imax)+porousc(ijm+imax))/2.0
                if (npc(ijm+imax).ne.1) then
                  xpstmp=porousp(ijm+imax)
                  xatmp=porousa(ijm+imax)
                  xxbtmp=porousb(ijm+imax)
                  d50tmp=porousd(ijm+imax)
                  xkc=max(abs(vn(ijm+imax)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict0=(xatmp+xbb
     &		    *sqrt(un(ijm+imax)**2+vn(ijm+imax)**2))*abs(gy)
                else
                  xfrict0=0.0
                endif
                if (npc(ijm+2*imax).ne.1) then
                  xpstmp=porousp(ijm+2*imax)
                  xatmp=porousa(ijm+2*imax)
                  xxbtmp=porousb(ijm+2*imax)
                  d50tmp=porousd(ijm+2*imax)
                  xkc=max(abs(vn(ijm+imax)),0.01d0)*max(1.0d0,xxt)
     &                  /xpstmp/d50tmp
                  xbb=xxbtmp*(1.0+7.5/xkc)
                  xfrict1=(xatmp+xbb
     &		    *sqrt(un(ijm+2*imax)**2+vn(ijm+imax)**2))*abs(gy)
                else
                  xfrict1=0.0
                endif
		xfriction=(xfrict0+xfrict1)/2.0
                gctmpb=1.0/(gctmpba+xfriction*delt)
		srce(nm)=srce(nm)-gy*delx(i)*gctmpb*at(ijm+imax)
            else
                srce(nm)=srce(nm)-gy*delx(i)*at(ijm+imax)
            endif
        else
            srce(nm)=srce(nm)-gy*delx(i)*at(ijm+imax)
        endif

2002	  continue

126     nm=nm+1
        ijm=ijm+1
        imjm=imjm+1
        ipjm=ipjm+1
  140 continue
c
  145 continue

      itmax=itmxiccg
      its=itmxiccg
      error=erriccg
      nxm=ibar
      nym=jbar

	if (t.le.1000000.0) goto 147
        n=1
        do 143 j=2,jmax-1
        do 143 i=2,imax-1
	ap(n)=a(n)
	bp(n)=b(n)
	bfp(n)=bf(n)
	dp(n)=d(n)
	dfp(n)=df(n)
	sp(n)=srce(n)
	write(9,'(2i3,6e10.3,1e20.13)')i,j,a(n),b(n),bf(n),
     &	d(n),df(n),
     &		soln(n),srce(n)
        n=n+1
143     continue
147	continue

      if (sym) then
c       ------------------------------------------------------------
        call cgitj(a,b,c,d,e,af,bf,cf,df,ef,nxm,nym,soln,srce,residu,
     &              sc1,sc2,its,error,xratio,yratio,xnorm,ynorm)
c       ------------------------------------------------------------
      else
c       ------------------------------------------------------------
        call ilucgj(ef,df,cf,bf,a,b,c,d,e,eh,dh,ch,bh,ag,bg,cg,dg,
     &              eg,nxm,nym,soln,srce,residu,sc1,sc2,its,error,
     &              xratio,yratio,xnorm,ynorm)
c       ------------------------------------------------------------
      endif
c
      iter=its
      if (iter .gt. itmax) then
        write(iotty,150) ncyc,t,its
        write(13,150) ncyc,t,its
        nocon=nocon+1
c	call prtplt(3)
      endif
c
      nm=1
      do 200 j=2,jm1
        do 200 i=2,im1
          ij=(j-1)*imax+i
          p(ij)=soln(nm)+psat
C.........attempt to apply variable pressure; fail
C	   p(ij)=soln(nm)+0.05*x(i)
C.........end
          nm=nm+1
  200 continue
c

      if (iter.gt.itmax) nexit=9999

c.... Set pressures in the ghost cells
c
c     --------------------------------------------------------
      call bdycell(im1,jm1,1,1,1,1,pbcl,pbcr,pbct,pbcb,pbc,p,
     &	xi,yj,gx,gy)
c     --------------------------------------------------------
c
c     -----------
 9999 continue
	call accel
c     -----------
c
      return
  150 format("IMPLCTP pressure solution for cycle ",i5," , time ",
     &        1pe12.5,/,"did not converge after ",i5," iterations")
      end


      subroutine initreg
c
c ======================================================================
c
c   Purpose -
c     initialize region quantities
c
c   INITREG is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c           SETUP
c
c
c   INITREG calls the following subroutines and functions -
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
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... default initial velocities are zero
c
c     -------------------------------
      call setarry (u(1),0.0d0,nxy)
      call setarry (v(1),0.0d0,nxy)
c     -------------------------------
c
c.... set initial velocity field into u and v arrays
c
      do 100 i=1,imax
	etah(i)=0.0
        do 100 j=1,jmax
          ij=(j-1)*imax+i
          ijp=ij+imax
          ipj=ij+1
          if ((at(ij).gt.em6).and.((f(ij).gt.emf).or.
     &          (f(ijp).gt.emf))) v(ij)=vi
          if ((ar(ij).gt.em6).and.((f(ij).gt.emf).or.
     &          (f(ipj).gt.emf))) u(ij)=ui
          un(ij)=u(ij)
          vn(ij)=v(ij)
  100 continue
c
      return
      end

 
      subroutine initvoff
c
c ======================================================================
c
c   Purpose -
c     initialize any free surfaces (via the volume of fluid function)
c     with the function f(x,y); for f(x,y) < 0, the point (x,y) is
c     either within the fluid or beyond the free surface for ifh=0
c     or 1, respectively
c
c   INITVOFF is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c           SETUP
c
c
c   INITVOFF calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c             FXY    ripple
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
      dimension iflg(5), dis(4), xm(5), ym(5)
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      xmode=2.*pi/(x(im1)-x(1))
      ymode=2.*pi/(y(jm1)-y(1))
c
      do 230 k=1,nfrsrf
        do 221 j=2,jm1
          do 220 i=2,im1
            ij=(j-1)*imax+i
            rdxdy=1.0/(delx(i)*dely(j))
            do 60 m=1,4
              go to (10,20,30,40), m
   10         x1=x(i)
              y1=y(j-1)
              dis(1)=dely(j)
              go to 50
   20         y1=y(j)
              x1=x(i)
              dis(2)=delx(i)
              go to 50
   30         x1=x(i-1)
              y1=y(j)
              dis(3)=dely(j)
              go to 50
   40         y1=y(j-1)
              x1=x(i-1)
              dis(4)=delx(i)
   50         iflg(m)=0
c             ---------------------------------------------------------
              fconic=fxy(x1,y1,fa1(k),fa2(k),fb1(k),fb2(k),fc1(k),
     &                   fc2(k),fd1(k),nxf(k)*xmode,fd2(k),mxf(k)*xmode,
     &                   fe1(k),nyf(k)*ymode,fe2(k),myf(k)*ymode)
c             ---------------------------------------------------------
              if (abs(fconic).le.1.0d-12) fconic=0.0
              if (fconic.le.0.0) iflg(m)=1
              xm(m)=x1
              ym(m)=y1
   60       continue
c
            iflg(5)=iflg(1)
            xm(5)=xm(1)
            ym(5)=ym(1)
            iflgs=0
c
            do 70 m=1,4
              iflgs=iflgs+iflg(m)
   70       continue
c
            brij=0.0
            btij=0.0
            if (iflgs.eq.0) go to 220
            if (iflgs.lt.4) go to 80
            bij=1.0d0
            brij=1.0d0
            btij=1.0d0
            go to 200
   80       if (iflg(1).eq.1.and.iflg(2).eq.1) brij=1.0d0
            if (iflg(2).eq.1.and.iflg(3).eq.1) btij=1.0d0
c
            do 160 m=1,4
              if (iflg(m).eq.iflg(m+1)) go to 160
              x1=xm(m)
              y1=ym(m)
              x2=xm(m+1)
              y2=ym(m+1)
              if (iflg(m).eq.0) go to 90
              x2=xm(m)
              y2=ym(m)
              x1=xm(m+1)
              y1=ym(m+1)
   90         epsif=0.001*(abs(x2-x1)+abs(y2-y1))
              smn=0.0
c             -------------------------------------------------------
              fmn=fxy(x2,y2,fa1(k),fa2(k),fb1(k),fb2(k),fc1(k),
     &                fc2(k),fd1(k),nxf(k)*xmode,fd2(k),mxf(k)*xmode,
     &                fe1(k),nyf(k)*ymode,fe2(k),myf(k)*ymode)
c             -------------------------------------------------------
              smx=1.0d0
c             -------------------------------------------------------
              fmx=fxy(x1,y1,fa1(k),fa2(k),fb1(k),fb2(k),fc1(k),
     &                fc2(k),fd1(k),nxf(k)*xmode,fd2(k),mxf(k)*xmode,
     &                fe1(k),nyf(k)*ymode,fe2(k),myf(k)*ymode)
c             -------------------------------------------------------
              s=0.5
  100         xt=s*x1+(1.0-s)*x2
              yt=s*y1+(1.0-s)*y2
c             -------------------------------------------------------
              fs=fxy(xt,yt,fa1(k),fa2(k),fb1(k),fb2(k),fc1(k),
     &                fc2(k),fd1(k),nxf(k)*xmode,fd2(k),mxf(k)*xmode,
     &                fe1(k),nyf(k)*ymode,fe2(k),myf(k)*ymode)
c             -------------------------------------------------------
              if (abs(fs).lt.epsif) go to 130
              if (fs.ge.0.0) go to 110
              fden=abs(fs-fmn)+1.0d-10
              se=s-fs*(s-smn)/fden
              if (se.gt.smx) se=smx
              fmn=fs
              smn=s
              go to 120
  110         fden=abs(fmx-fs)+1.0d-10
              se=s-fs*(smx-s)/fden
              if (se.lt.smn) se=smn
              fmx=fs
              smx=s
  120         si=s-fs*(smx-smn)/(fmx-fmn)
              s=0.5*(se+si)
              go to 100
  130         dis(m)=sqrt((xt-x2)**2+(yt-y2)**2)
              go to (140,150,160,160), m
  140         brij=dis(1)/dely(j)
              go to 160
  150         btij=dis(2)/delx(i)
  160       continue
c
            m=0
            bij=0.0
  170       continue
            m=m+1
            if (m.eq.5) go to 190
            if (iflg(m).eq.0) go to 170
            mp1=m+1
            if (mp1.eq.5) mp1=1
            mm1=m-1
            if (mm1.eq.0) mm1=4
            bij=bij+dis(m)*dis(mm1)
            if (iflg(mp1).eq.1) go to 180
            dis2=dis(m)
  180       continue
            if (iflg(mm1).eq.1) go to 170
            dis1=dis(mm1)
            go to 170
  190       continue
            if (iflgs.eq.3) bij=bij-dis1*dis2
            bij=0.5*bij*rdxdy
            if (bij.gt.1.0d0) bij=1.0d0
  200       continue
            if (ifh(k).eq.0) go to 210
            bij=-bij
            brij=-brij
            btij=-btij
  210       f(ij)=f(ij)+bij
            if(f(ij).gt.0.9999) f(ij)=1.0d0
            if(f(ij).lt.0.0001) f(ij)=0.0
            fn(ij)=f(ij)
  220     continue
  221   continue
  230 continue

	if (novertop.eq.1) then
	  do j=2,jm1
          do i=2,im1
            ij=(j-1)*imax+i
		  if (x(i).gt.xovertop) f(ij)=0.0
		end do
	  end do
	endif
c
      return
      end

      subroutine kill (iotty,ncyc,routine,nexit)
c
c ======================================================================
c
c   Purpose -
c     graceful termination
c
c   KILL is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c           EQUIB   MESHSET    NEWCYC     TAPIN    TAPOUT
c
c
c
c   KILL calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c          FEMPTY    system     EXITA    system
c
c ======================================================================
c
c##############################################################
      implicit real*8 (a-h,o-z)
      
c##############################################################
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      write (9,10) ncyc,routine
      write (13,10) ncyc,routine
      write (iotty,10) ncyc,routine
c
c     --------------
	nexit=9998
c     --------------
c
      return
   10 format(5x," *  *  *  *  *  *  *  *  *  *  *  *  *  *",/,
     1       " RIPPLE terminated on cycle ",i5,
     2       " in subroutine ",a8,/,
     3       5x," *  *  *  *  *  *  *  *  *  *  *  *  *  *")
      end


      subroutine ldlt(n,m,e,d,c,b,a,bu,cu,du,eu,x,y)
c
c***purpose
c
c      solve (incomplete) ax=(l*d*u)x=y, where the matrix a is
c      a 9-diagonal matrix associated with a 9-point difference
c      operator. factors are generated lu9p.
c
c***description
c
c       ldlt requires that the coefficients of the ith
c       equation be stored in the ith element of the diagonal
c       arrays. nine diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c       the incomplete factors of a are stored over the input matrix.
c
c     on entry
c
c        e,d,c,b,a,bu,cu,du,eu,   are 1-d arrays of length
c            at least nxm. the arrays contain, from left to right,
c            the diagonals of a as one sweeps a left to right.
c
c        n     is the size  or order of each tridiagonal block if
c              a is viewed as a block tridiagonal system of equations.
c
c        m     is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c
c        y     is an input vector of dimenisonat least nxm.
c
c     on return
c
c        x     is the output vector of dimension at least
c              nxm containing the solution vector.
c
c***routines called   none
c
      implicit real*8 (a-h,o-z)
      dimension e(n,m),d(n,m),c(n,m),b(n,m),a(n,m),
     1          bu(n,m),cu(n,m),du(n,m),eu(n,m),x(n,m),y(n,m)
c*** first executable statement   ldlt
c
ccc   copy y into x
c
      do 1 j=1,m
        do 1 i=1,n
          x(i,j) = y(i,j)
    1 continue
c
ccc   forward solve
c
      do 100 j=1,m
        if (j .ne. 1)    then
          do 15 i=1,n-1
            x(i,j) = x(i,j)-c(i,j)*x(i+1,j-1)
   15     continue
          do 20 i=1,n
            x(i,j) = x(i,j)-d(i,j)*x(i,j-1)
   20     continue
          do 25 i=2,n
            x(i,j) = x(i,j)-e(i,j)*x(i-1,j-1)
   25     continue
c
ccc vertical periodic equations
c
        if (j .eq. m)    then
          do 16 i=1,n-1
            x(i,m) = x(i,m)-c(i,1)*x(i+1,1)
   16     continue
          do 21 i=1,n
            x(i,m) = x(i,m)-d(i,1)*x(i,1)
   21     continue
          do 26 i=2,n
            x(i,m) = x(i,m)-e(i,1)*x(i-1,1)
   26     continue
        end if
        end if
c
ccc   van der worst approximation
c
cdir$ ivdep
        do 5 i=n,2,-1
          x(i,j) = x(i,j)-b(i,j)*x(i-1,j)
    5   continue
cdir$ ivdep
        do 10 i=n,3,-1
          x(i,j) = x(i,j)+b(i-1,j)*b(i,j)*x(i-2,j)
   10   continue
  100 continue
c
ccc   diagonal solve
c
      do 30 j=1,m
        do 30 i=1,n
          x(i,j) = a(i,j)*x(i,j)
   30 continue
c
ccc  backward solve
c
      do 200 j=m,1,-1
        if (j .ne. m)    then
          do 45 i=2,n
            x(i,j) = x(i,j)-cu(i,j)*x(i-1,j+1)
   45     continue
          do 50 i=1,n
            x(i,j) = x(i,j)-du(i,j)*x(i,j+1)
   50     continue
          do 55 i=1,n-1
            x(i,j) = x(i,j)-eu(i,j)*x(i+1,j+1)
   55     continue
c
ccc vertical periodic equations
c
        if (j .eq. 1)   then
          do 46 i=2,n
            x(i,1) = x(i,1)-c(i-1,1)*x(i-1,m)
   46     continue
          do 51 i=1,n
            x(i,1) = x(i,1)-d(i,1)*x(i,m)
   51     continue
          do 56 i=1,n-1
            x(i,1) = x(i,1)-e(i+1,1)*x(i+1,m)
   56     continue
        end if
        end if
c
ccc van der worst approximation
c
cdir$ ivdep
        do 35 i=1,n-1
          x(i,j) = x(i,j)-bu(i,j)*x(i+1,j)
   35   continue
cdir$ ivdep
        do 40 i=1,n-2
          x(i,j) = x(i,j)+bu(i+1,j)*bu(i,j)*x(i+2,j)
   40   continue
  200 continue
      return
      end


      subroutine ldlt9f(n,m,e,d,c,b,a)
c
c***purpose
c
c      factor (incomplete) a=l*d*lt, where the matrix a is
c      a 5-diagonal matrix associated with a 9-point difference
c      s.p.d, operator.
c
c***description
c
c       ldlt9f requires that the coefficients of the ith
c       equation be stored in the ith element of the diagonal
c       arrays. five diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c       the incomplete factors of a are stored over the input matrix.
c
c     on entry
c
c        e,d,c,b,a   are 1-d arrays of length
c            at least nxm. the arrays contain, from left to right,
c            the lower diagonals of a as one sweeps a left to right.
c
c        n    is the size  or order of each tridiagonal block if
c              a is viewed as a block tridiagonal system of equations.
c
c        m    is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c
c     on return
c
c        e,d,c,b,a    contain the factors
c            ldlt. the reciprocal of d is stored in array a.
c
c***routines called   none
c
c*** first executable statement   ldlt9f
      implicit real*8 (a-h,o-z)
      dimension e(n,m),d(n,m),c(n,m),b(n,m),a(n,m)
c
ccc   form a = l*v
c
      do 100 j=1,m
c
ccc factor central block
c
        a(1,j) = 1./a(1,j)
        do 5 i=2,n
          a(i,j) = 1./(a(i,j)-b(i,j)*(a(i-1,j)*b(i,j)))
    5   continue
        if (j .ne. m)    then
c
ccc   form lower diagonals
c
        do 10 i=2,n
          d(i,j+1) = d(i,j+1)-(a(i-1,j)*e(i,j+1))*b(i,j)
   10   continue
        do 15 i=1,n-1
          c(i,j+1) = c(i,j+1)-(a(i,j)*d(i,j+1))*b(i+1,j)
   15   continue
c
ccc factor lower periodic diagonals
c
        if (j .eq. 1)    then
        do 11 i=2,n
          d(i,1) = d(i,1)-(a(i-1,1)*e(i,1))*b(i,1)
   11   continue
        do 16 i=1,n-1
          c(i,1) = c(i,1)-(a(i,1)*d(i,1))*b(i+1,1)
   16   continue
        end if
c
ccc   form central diagonals
c
        do 20 i=2,n
          b(i,j+1) = b(i,j+1)-(a(i-1,j)*e(i,j+1))*d(i-1,j+1)
     1            -(a(i,j)*d(i,j+1))*c(i-1,j+1)
   20   continue
          a(1,j+1) = a(1,j+1)-(a(1,j)*d(1,j+1))*d(1,j+1)
     1             -(a(2,j)*c(1,j+1))*c(1,j+1)
        do 30 i=2,n-1
          a(i,j+1) = a(i,j+1)-(a(i-1,j)*e(i,j+1))*e(i,j+1)
     1            -(a(i,j)*d(i,j+1))*d(i,j+1)
     2            -(a(i+1,j)*c(i,j+1))*c(i,j+1)
   30   continue
          a(n,j+1) = a(n,j+1)-(a(n-1,j)*e(n,j+1))*e(n,j+1)
     1            -(a(n,j)*d(n,j+1))*d(n,j+1)
c
ccc factor central periodic block
c
        if (j .eq. 1)    then
        do 21 i=2,n
          b(i,m) = b(i,m)-(a(i-1,1)*e(i,1))*d(i-1,1)
     1            -(a(i,1)*d(i,1))*c(i-1,1)
   21   continue
          a(1,m) = a(1,m)-(a(1,1)*d(1,1))*d(1,1)
     1             -(a(2,1)*c(1,1))*c(1,1)
        do 31 i=2,n-1
          a(i,m) = a(i,m)-(a(i-1,1)*e(i,1))*e(i,1)
     1            -(a(i,1)*d(i,1))*d(i,1)
     2            -(a(i+1,1)*c(i,1))*c(i,1)
   31   continue
          a(n,m) = a(n,m)-(a(n-1,1)*e(n,1))*e(n,1)
     1            -(a(n,1)*d(n,1))*d(n,1)
        end if
c
ccc scale diagonals
c
          do 40 i=2,n
          e(i,j+1) = (a(i-1,j)*e(i,j+1))
   40   continue
        do 50 i=1,n
          d(i,j+1) = (a(i,j)*d(i,j+1))
   50   continue
        do 60 i=1,n-1
          c(i,j+1) = (a(i+1,j)*c(i,j+1))
   60   continue
c
ccc scale periodic block
c
        if (j .eq. 1)    then
          do 41 i=2,n
          e(i,1) = (a(i-1,1)*e(i,1))
   41   continue
        do 51 i=1,n
          d(i,1) = (a(i,1)*d(i,1))
   51   continue
        do 61 i=1,n-1
          c(i,1) = (a(i+1,1)*c(i,1))
   61   continue
        end if
        end if
        do 70 i=2,n
          b(i,j) = a(i-1,j)*b(i,j)
   70 continue
  100 continue
      return
      end
       

      subroutine lu9p(n,m,e,d,c,b,a,bu,cu,du,eu)
c
c***purpose
c
c      factor (incomplete) a=l*d*u, where the matrix a is
c      a 9-diagonal matrix associated with a 9-point difference
c      operator.
c
c***description
c
c       lu9p requires that the coefficients of the ith
c       equation be stored in the ith element of the diagonal
c       arrays. nine diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c       the incomplete factors of a are stored over the input matrix.
c
c     on entry
c
c        e,d,c,b,a,bu,cu,du,eu,   are 1-d arrays of length
c            at least nxm. the arrays contain, from left to right,
c            the diagonals of a as one sweeps a left to right.
c
c        n    is the size  or order of each tridiagonal block if
c              a is viewed as a block tridiagonal system of equations.
c
c        m    is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c
c     on return
c
c        e,d,c,b,a,bu,cu,du,eu    contain the factors
c            ldu. the reciprocal of d is stored in array a.
c
c***routines called   none
c
c*** first executable statement   lu9p
      implicit real*8 (a-h,o-z)
      dimension e(n,m),d(n,m),c(n,m),b(n,m),a(n,m),
     1        bu(n,m),cu(n,m),du(n,m),eu(n,m)
c
ccc   form a = l*v
c
      do 100 j=1,m
c
ccc factor central block
c
        a(1,j) = 1./a(1,j)
        do 5 i=2,n
          b(i,j) = a(i-1,j)*b(i,j)
          a(i,j) = 1./(a(i,j)-b(i,j)*bu(i-1,j))
    5   continue
        if (j .ne. m)    then
c
ccc   form upper diagonals
c
          do 10 i=2,n
            du(i,j) = du(i,j)-b(i,j)*eu(i-1,j)
   10     continue
          do 15 i=2,n
            cu(i,j) = cu(i,j)-b(i,j)*du(i-1,j)
   15     continue
c
ccc   form lower diagonals
c
          do 20 i=2,n
           e(i,j+1) = a(i-1,j)*e(i,j+1)
   20     continue
          d(1,j+1) = a(1,j)*d(1,j+1)
          do 25 i=2,n
            d(i,j+1) = a(i,j)*(d(i,j+1)-e(i,j+1)*bu(i-1,j))
   25     continue
          do 30 i=1,n-1
            c(i,j+1) = a(i+1,j)*(c(i,j+1)-d(i,j+1)*bu(i,j))
   30     continue
c
ccc   form central diagonals
c
          do 35 i=2,n
            b(i,j+1) = b(i,j+1)-e(i,j+1)*du(i-1,j)-d(i,j+1)*cu(i,j)
   35     continue
          do 40 i=1,n-1
            bu(i,j+1) = bu(i,j+1)-d(i,j+1)*eu(i,j)-c(i,j+1)*du(i+1,j)
   40     continue
          a(1,j+1) = a(1,j+1)-d(1,j+1)*du(1,j)-c(1,j+1)*cu(2,j)
          do 45 i=2,n-1
            a(i,j+1) = a(i,j+1)-e(i,j+1)*eu(i-1,j)
     1                -d(i,j+1)*du(i,j)-c(i,j+1)*cu(i+1,j)
   45     continue
          a(n,j+1) = a(n,j+1)-e(n,j+1)*eu(n-1,j)-d(n,j+1)*du(n,j)
        end if
  100 continue
      do 110 j=1,m
        do 110 i=1,n
          bu(i,j) = a(i,j)*bu(i,j)
  110 continue
c
ccc   form u = d[-1]*v
c
      do 120 j=1,m-1
        do 120 i=1,n
          cu(i,j) = a(i,j)*cu(i,j)
          du(i,j) = a(i,j)*du(i,j)
          eu(i,j) = a(i,j)*eu(i,j)
  120 continue
      return
      end


      subroutine matmul9(e,d,c,b,a,bu,cu,du,eu,nh,nv,x,y)
c
c***purpose
c
c      form the matrix product a*x=y,where the matrix a is
c      a 9-diagonal matrix associated with a 9-point difference
c      operator.
c
c***description
c
c       matmul9 requires that the coefficients of the ith
c       equation be stored in the ith element of the diagonal
c       arrays. nine diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c
c     on entry
c
c        e,d,c,b,a,bu,cu,du,eu,   are 1-d arrays of length
c            at least nhxnv. the arrays contain, from left to right,
c            the diagonals of a as one sweeps a left to right.
c
c        nh    is the size  or order of each tridiagonal block if
c              a is viewed as a block tridiagonal system of equations.
c
c        nv    is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c
c        x     is the input vector of dimenison at least nhxnv.
c
c     on return
c
c        y     is the output vector of dimension at least
c              nhxnv containing a*x.
c
c
c***routines called   none
c
      implicit real*8 (a-h,o-z)
      real*8 e(*),d(*),c(*),b(*),a(*),bu(*),cu(*),du(*),eu(*),
     &     x(*),y(*)
c
c*** first executable statement   matmul9
      n = nh*nv
      nhm1 = nh-1
      nhp1 = nh+1
      i = 1
      y(i)=a(i)*x(i)+bu(i)*x(i+1)
     1      +du(i)*x(i+nh)+eu(i)*x(i+nhp1)
      do 10 i=2,nh
      y(i)=b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1)
     1        +(cu(i)*x(i+nhm1)+(du(i)*x(i+nh)+(eu(i)*x(i+nhp1)
     2        )))))
   10 continue
      i = nhp1
      y(i)=d(i)*x(i-nh)+c(i)*x(i-nhm1)
     1        +b(i)*x(i-1)+a(i)*x(i)+bu(i)*x(i+1)
     2        +cu(i)*x(i+nhm1)+du(i)*x(i+nh)+eu(i)*x(i+nhp1)
      do 20 i=nh+2,n-nhp1
      y(i)=e(i)*x(i-nhp1)+(d(i)*x(i-nh)+(c(i)*x(i-nhm1)
     1        +(b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1)
     2        +(cu(i)*x(i+nhm1)+(du(i)*x(i+nh)+(eu(i)*x(i+nhp1)
     3        ))))))))
   20 continue
      i=n-nh
      y(i)=e(i)*x(i-nhp1)+d(i)*x(i-nh)+c(i)*x(i-nhm1)
     1        +b(i)*x(i-1)+a(i)*x(i)+bu(i)*x(i+1)
     2        +cu(i)*x(i+nh-1)+du(i)*x(i+nh)
      do 30 i=n-nhm1,n-1
      y(i)=e(i)*x(i-nhp1)+(d(i)*x(i-nh)+(c(i)*x(i-nhm1)
     1        +(b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1)
     2        )))))
   30 continue
      i = n
      y(i)=e(i)*x(i-nhp1)+d(i)*x(i-nh)
     1      +b(i)*x(i-1)+a(i)*x(i)
c
ccc periodic equations
c
      l = n-nh
      k = l-1
      m = l+1
      i = 1
      y(i) = y(i)+du(i+l)*x(i+l)+eu(i+l)*x(i+m)
      y(i+l) = y(i+l)+du(i+l)*x(i)+cu(i+m)*x(i+1)
      do 40 i=2,nhm1
        y(i) = y(i)+cu(i+l)*x(i+k)+du(i+l)*x(i+l)
     1        +eu(i+l)*x(i+m)
        y(i+l) = y(i+l)+eu(i+k)*x(i-1)+du(i+l)*x(i)
     1        +cu(i+m)*x(i+1)
   40 continue
      i = nh
      y(i) = y(i)+cu(i+l)*x(i+k)+du(i+l)*x(i+l)
      y(i+l) = y(i+l)+eu(i+k)*x(i-1)+du(i+l)*x(i)
      return
      end
  
 
      subroutine meshset
c
c ======================================================================
c
c   Purpose -
c     mesh generator
c
c   MESHSET is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c           SETUP
c
c
c   MESHSET calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c            KILL    ripple      GEOM    ripple
c
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
      data one /1.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      i=1
      j=1
      x(1)=xl(1)
      y(1)=yl(1)
c
      do 40 k=1,nkx
        if (nxl(k).eq.0) go to 20
        dxml=(xc(k)-xl(k))/nxl(k)
        nt=nxl(k)
        tn=nt
        tn=max(tn,one+em6)
        dxmn1=min(dxmn(k),dxml)
        cmc=(xc(k)-xl(k)-tn*dxmn1)*tn/(tn-1.0)
        if (nt.eq.1) cmc=0.0
        bmc=xc(k)-xl(k)-cmc
        do 10 l=1,nt
          i=i+1
          rln=(float(l)-tn)/tn
   10   x(i)=xc(k)+bmc*rln-cmc*rln*rln
   20   if (nxr(k).eq.0) go to 40
        dxmr=(xl(k+1)-xc(k))/nxr(k)
        nt=nxr(k)
        tn=nt
        tn=max(tn,one+em6)
        dxmn1=min(dxmn(k),dxmr)
        cmc=(xl(k+1)-xc(k)-tn*dxmn1)*tn/(tn-1.0)
        if (nt.eq.1) cmc=0.0
        bmc=xl(k+1)-xc(k)-cmc
        do 30 l=1,nt
          i=i+1
          rln=float(l)/tn
   30   x(i)=xc(k)+bmc*rln+cmc*rln*rln
   40 continue
c
      if (kr.ne.4) go to 50
      i=i+1
      x(i)=x(i-1)+x(2)-x(1)
   50 continue
c
      do 90 k=1,nky
        if (nyl(k).eq.0) go to 70
        dyml=(yc(k)-yl(k))/nyl(k)
        nt=nyl(k)
        tn=nt
        tn=max(tn,one+em6)
        dymn1=min(dymn(k),dyml)
        cmc=(yc(k)-yl(k)-tn*dymn1)*tn/(tn-1.0)
        if (nt.eq.1) cmc=0.0
        bmc=yc(k)-yl(k)-cmc
        do 60 l=1,nt
          j=j+1
          rln=(float(l)-tn)/tn
   60   y(j)=yc(k)+bmc*rln-cmc*rln*rln
   70   if (nyr(k).eq.0) go to 90
        dymr=(yl(k+1)-yc(k))/nyr(k)
        nt=nyr(k)
        tn=nt
        tn=max(tn,one+em6)
        dymn1=min(dymn(k),dymr)
        cmc=(yl(k+1)-yc(k)-tn*dymn1)*tn/(tn-1.0)
        if (nt.eq.1) cmc=0.0
        bmc=yl(k+1)-yc(k)-cmc
        do 80 l=1,nt
          j=j+1
          rln=float(l)/tn
   80   y(j)=yc(k)+bmc*rln+cmc*rln*rln
   90 continue
c
      if (kt.ne.4) go to 100
      j=j+1
      y(j)=y(j-1)+y(2)-y(1)
  100 continue
      numx=i
      numy=j
      numxm1=numx-1
      numym1=numy-1
      numxp1=numx+1
      numyp1=numy+1
      ibar=numx-1
      jbar=numy-1
      imax=ibar+2
      jmax=jbar+2
      im1=imax-1
      jm1=jmax-1
c
c.... calculate values needed for variable mesh
c
      do 120 i=1,numx
        if (x(i).eq.0.0) go to 110
        rx(i)=1.0/x(i)
        go to 120
  110   rx(i)=0.0
  120 continue
c
      do 130 i=2,numx
        xi(i)=0.5*(x(i-1)+x(i))
        delx(i)=x(i)-x(i-1)
  130 rdx(i)=1.0/delx(i)
c
      delx(1)=delx(2)
      xi(1)=xi(2)-delx(2)
      rdx(1)=1.0/delx(1)
      delxa=delx(numx)
      if (kr.eq.4) delxa=delx(3)
      delx(numxp1)=delxa
      xi(numxp1)=xi(numx)+delxa
      x(numxp1)=xi(numxp1)+0.5*delx(numxp1)
      rdx(numxp1)=1.0/delx(numxp1)
c
      do 140 i=2,numy
        yj(i)=0.5*(y(i-1)+y(i))
        dely(i)=y(i)-y(i-1)
        rdy(i)=1.0/dely(i)
  140 continue
c
      dely(1)=dely(2)
      rdy(1)=1.0/dely(1)
      yj(1)=yj(2)-dely(2)
      delya=dely(numy)
      if (kt.eq.4) delya=dely(3)
      dely(numyp1)=delya
      yj(numyp1)=yj(numy)+delya
      y(numyp1)=yj(numyp1)+0.5*dely(numyp1)
      rdy(numyp1)=1.0/dely(numyp1)
c
c.... set r and ri array for plane or cylindrical geometry
c
      do 145 i=1,imax
        r(i)=x(i)
        ri(i)=xi(i)
        if(icyl.eq.1) go to 145
        r(i)=1.0d0
        ri(i)=1.0d0
  145 continue
c
c.... get the metric for the mesh
c
c     ----------
      call geom
c     ----------
c
c.... write out mesh information
c
      write (13,210)
      do 190 i=1,numxp1
        write (13,220) i,x(i),i,rx(i),i,delx(i),i,rdx(i),i,xi(i)

	write(31,'(2f12.5)')x(i)
  190 continue
      close(31)      
c
      write (13,210)
      do 200 i=1,numyp1
        write (13,230) i,y(i),i,dely(i),i,rdy(i),i,yj(i)
	write(32,'(2f12.5)')y(i)
  200 continue
      close(32)
c
c.... test array size
c
      if (imax.le.ibar2.and.jmax.le.jbar2) go to 9999
c
c                    * * * * error section * * * *
c
      write (9,240)
c
c     ---------------------------------
c	routine="MESHSET"
c      call kill (iotty,ncyc,routine,nexit)
c     ---------------------------------
c
 9999 return
c
  210 format (1h1)
  220 format (1x,2hx(,i2,2h)=,1pe12.5,2x,3hrx(,i2,2h)=,1pe12.5,2x,
     &        5hdelx(,i2,2h)=,1pe12.5,1x,4hrdx(,i2,2h)=,1pe12.5,2x,
     &        3hxi(,i2,2h)=,1pe12.5)
  230 format (1x,2hy(,i2,2h)=,1pe12.5,3x,5hdely(,i2,2h)=,1pe12.5,3x,
     &        4hrdy(,i2,2h)=,1pe12.5,3x,3hyj(,i2,2h)=,1pe12.5)
  240 format (45h mesh size inconsistent with array dimensions)
      end

 
      subroutine obs
c
c ======================================================================
c
c   Purpose -
c     find the physical location of any obstacles
c
c   OBS is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c            ASET
c
c
c   OBS calls the following subroutines and functions -
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
      npt=1
      do 171 i=2,im1
        atr=1.0-em6
        atl=1.0-em6
        atcc=1.0-em6
        do 170 j=2,jm1
          ij=(j-1)*imax+i
          imj=ij-1
          ijm=ij-imax
          if (ac(ij).lt.em6) go to 170
          afr=1.0d0
          aft=1.0d0
          afl=1.0d0
          afb=1.0d0
          if (ar(ij).lt.atr) afr=ar(ij)/atr
          if (at(ij).lt.atcc) aft=at(ij)/atcc
          if (ar(imj).lt.atl) afl=ar(imj)/atl
          if (at(ijm).lt.atcc) afb=at(ijm)/atcc
          go to 7
    7     continue
          if (ac(ij).ge.atcc) go to 120
          if ((aft+afb).lt.em6.or.(afl+afr).lt.em6) go to 170
          m=1
          amn=afb+afr
          if ((afr+aft).gt.amn) go to 10
          m=2
          amn=afr+aft
   10     if ((aft+afl).gt.amn) go to 20
          m=3
          amn=aft+afl
   20     if ((afl+afb).gt.amn) go to 30
          m=4
   30     go to (40,60,80,100), m
   40     x1=x(i-1)+aft*delx(i)
          y1=y(j)
          if (aft.lt.1.0d0) go to 50
          y1=y1-afr*dely(j)
   50     x2=x(i-1)
          y2=y(j)-afl*dely(j)
          if (afl.lt.1.0d0) go to 160
          x2=x2+afb*delx(i)
          go to 160
   60     x1=x(i-1)
          y1=y(j-1)+afl*dely(j)
          if (afl.lt.1.0d0) go to 70
          x1=x1+aft*delx(i)
   70     x2=x(i-1)+afb*delx(i)
          y2=y(j-1)
          if (afb.lt.1.0d0) go to 160
          y2=y2+afr*dely(j)
          go to 160
   80     x1=x(i)-afb*delx(i)
          y1=y(j-1)
          if (afb.lt.1.0d0) go to 90
          y1=y1+afl*dely(j)
   90     x2=x(i)
          y2=y(j-1)+afr*dely(j)
          if (afr.lt.1.0d0) go to 160
          x2=x2-aft*delx(i)
          go to 160
  100     x1=x(i)
          y1=y(j)-afr*dely(j)
          if (afr.lt.1.0d0) go to 110
          x1=x1-afb*delx(i)
  110     x2=x(i)-aft*delx(i)
          y2=y(j)
          if (aft.lt.1.0d0) go to 160
          y2=y2-afl*dely(j)
          go to 160
  120     if (afr.gt.em6) go to 130
          x1=x(i)
          y1=y(j-1)
          x2=x1
          y2=y(j)
          xobs(npt)=x1
          xobs(npt+1)=x2
          yobs(npt)=y1
          yobs(npt+1)=y2
          ijobs(ij)=npt
          npt=npt+2
  130     if (aft.gt.em6) go to 140
          x1=x(i-1)
          y1=y(j)
          x2=x(i)
          y2=y1
          xobs(npt)=x1
          xobs(npt+1)=x2
          yobs(npt)=y1
          yobs(npt+1)=y2
          ijobs(ij)=npt
          npt=npt+2
  140     if (afl.gt.em6) go to 150
          x1=x(i-1)
          y1=y(j)
          x2=x1
          y2=y(j-1)
          xobs(npt)=x1
          xobs(npt+1)=x2
          yobs(npt)=y1
          yobs(npt+1)=y2
          ijobs(ij)=npt
          npt=npt+2
  150     if (afb.gt.em6) go to 170
          x1=x(i-1)
          y1=y(j-1)
          x2=x(i)
          y2=y1
          xobs(npt)=x1
          xobs(npt+1)=x2
          yobs(npt)=y1
          yobs(npt+1)=y2
          ijobs(ij)=npt
          npt=npt+2
  160     xobs(npt)=x1
          xobs(npt+1)=x2
          yobs(npt)=y1
          yobs(npt+1)=y2
          ijobs(ij)=npt
          npt=npt+2
  170   continue
  171 continue
      nobspt=npt-1
c
      return
      end
  

      subroutine pollutant

c##############################################################
      implicit real*8 (a-h,o-z)  
c##############################################################
c
c############
      include  "comdk2.h"
c############
c
	do 10 j=2,jmax-1
	  do 10 i=2,imax-1
	    ij=(j-1)*imax+i
          ijm=ij-imax
          imj=ij-1
          ipj=ij+1
          ijp=ij+imax

	    if (ac(ij).le.em6) then
			goto 10
	    endif
	    if (nf(ij).gt.5) then
			xp(ij)=0.0
			goto 10
	    endif
          if (xpn(ij).eq.0.0.and.xpn(ipj).eq.0.0.and.xpn(imj).eq.
     &		0.0.and.xpn(ijp).eq.0.0.and.xpn(ijm).eq.0.0) then
              xp(ij)=0.0
              goto 10
          endif

C	it is expected that the no-flux boundary condition should be
C	satisfied automatically because the flux is represented by
C	pu\theta; on the boundary, \theta=0 and thus no flux is allowed

	    dxmin=(delx(i)+delx(i-1))/2.0
	    dxplus=(delx(i)+delx(i+1))/2.0
	    dymin=(dely(j)+dely(j-1))/2.0
	    dyplus=(dely(j)+dely(j+1))/2.0

C.......prepare the pollutant at the boundary
	    polr=(xpn(ij)*ac(ij)*delx(i+1)+xpn(ipj)*ac(ipj)*
     &		delx(i))/(delx(i)+delx(i+1))
            poll=(xpn(ij)*ac(ij)*delx(i-1)+xpn(imj)*ac(imj)*
     &		delx(i))/(delx(i)+delx(i-1))
            polt=(xpn(ij)*ac(ij)*dely(j+1)+xpn(ijp)*ac(ijp)*
     &		dely(j))/(dely(j)+dely(j+1))
            polb=(xpn(ij)*ac(ij)*dely(j-1)+xpn(ijm)*ac(ijm)*
     &		dely(j))/(dely(j)+dely(j-1))

            xnrp=(xnutn(ij)*delx(i+1)+xnutn(ipj)*delx(i))/
     &           (delx(i)+delx(i+1))+xnu
            if (((kr.eq.7.or.kr.eq.3).and.i.eq.im1).or.nf(ipj).ge.6)
     &          xnrp=xnutn(ij)+xnu
            xnlp=(xnutn(ij)*delx(i-1)+xnutn(imj)*delx(i))/
     &           (delx(i)+delx(i-1))+xnu
            if (((kl.eq.7.or.kl.eq.3).and.i.eq.2).or.nf(imj).ge.6)
     &          xnlp=xnutn(ij)+xnu
            xntp=(xnutn(ij)*dely(j+1)+xnutn(ijp)*dely(j))/
     &           (dely(j)+dely(j+1))+xnu
            if (((kt.eq.7.or.kt.eq.3).and.j.eq.jm1).or.nf(ijp).ge.6)
     &          xntp=xnutn(ij)+xnu
            xnbp=(xnutn(ij)*dely(j-1)+xnutn(ijm)*dely(j))/
     &           (dely(j)+dely(j-1))+xnu
            if (((kb.eq.7.or.kb.eq.3).and.j.eq.2).or.nf(ijm).ge.6)
     &          xnbp=xnutn(ij)+xnu

            dpdl=(xpn(ij)*ac(ij)-xpn(imj)*ac(imj))/dxmin
            dpdr=(xpn(ipj)*ac(ipj)-xpn(ij)*ac(ij))/dxplus
            dpdb=(xpn(ij)*ac(ij)-xpn(ijm)*ac(ijm))/dymin
            dpdt=(xpn(ijp)*ac(ijp)-xpn(ij)*ac(ij))/dyplus

C...........make sure the flux equal zero on the free surface and bottom
	    if (fn(ipj).eq.0.0.or.ar(ij).eq.0.0) then
		dpdr=0.0
		polr=0.0
	    endif
            if (fn(imj).eq.0.0.or.ar(imj).eq.0.0) then
                dpdl=0.0
		poll=0.0
            endif
            if (fn(ijp).eq.0.0.or.at(ij).eq.0.0) then
                dpdt=0.0
		polt=0.0
            endif
            if (fn(ijm).eq.0.0.or.at(ijm).eq.0.0) then
                dpdb=0.0
		polb=0.0
            endif

C........update pollutant equation

C........using conservative Lax-Wendroff method
            fpx=(u(ij)*ar(ij)*polr-u(imj)*ar(imj)*poll)/delx(i)
C........Lax-Wendroff compensation terms
     &		-(dpdr*(u(ij)*ar(ij))**2*delt/2.0
     &		-dpdl*(u(imj)*ar(imj))**2*delt/2.0)/delx(i)

            fpy=(v(ij)*at(ij)*polt-v(ijm)*at(ijm)*polb)/dely(j)
C........Lax-Wendroff compensation terms
     &		-(dpdt*(v(ij)*at(ij))**2*delt/2.0
     &		-dpdb*(v(ijm)*at(ijm))**2*delt/2.0)/dely(j)

            vispx=(dpdr*xnrp-dpdl*xnlp)/delx(i)

            vispy=(dpdt*xntp-dpdb*xnbp)/dely(j)

            visp=vispx+vispy

	    txp=xp(ij)
	    xp(ij)=xpn(ij)+delt*(-fpx-fpy+visp)/ac(ij)
10	continue

C.......special way to redistribute the pollutant when occasional the
C.......negative diffusion occurs to cause the negative value
	do 15 j=2,jmax-1
	  do 15 i=2,imax-1
		ij=(j-1)*imax+i
          ijm=ij-imax
          imj=ij-1
          ipj=ij+1
          ijp=ij+imax
		if (xp(ij).lt.0.0) then
		  xpmax=max(abs(xp(ipj)),abs(xp(imj)),
     &			abs(xp(ijp)),abs(xp(ijm)))
		  if (xpmax.eq.abs(xp(ipj))) then
			xp(ipj)=xp(ipj)+xp(ij)*ac(ij)*delx(i)
     &			/(ac(ipj)*delx(i+1))
			xp(ij)=0.0
			goto 15
		  endif
                  if (xpmax.eq.abs(xp(imj))) then
                        xp(imj)=xp(imj)+xp(ij)*ac(ij)*delx(i)
     &			/(ac(imj)*delx(i-1))
                        xp(ij)=0.0
                        goto 15
                  endif
                  if (xpmax.eq.abs(xp(ijp))) then
                        xp(ijp)=xp(ijp)+xp(ij)*ac(ij)*dely(j)
     &			/(ac(ijp)*dely(j+1))
                        xp(ij)=0.0
                        goto 15
                  endif
                  if (xpmax.eq.abs(xp(ijm))) then
                        xp(ijm)=xp(ijm)+xp(ij)*ac(ij)*dely(j)
     &			/(ac(ijm)*dely(j-1))
                        xp(ij)=0.0
                        goto 15
                  endif
		endif
15	continue

CCC	special boundary on wall for inflow and outflow
        il=2
        ir=im1
        do 65 j=1,jmax
          ijl=(j-1)*imax+il
          ijr=(j-1)*imax+ir
          if (kl.ne.6) then
	    xp(ijl-1)=xp(ijl)
	  endif
	  xp(ijr+1)=xp(ijr)
65      continue
	
	jb=2
	jt=jm1
	do 55 i=1,imax
	  ijb=(jb-1)*imax+i
	  ijt=(jt-1)*imax+i
	  xp(ijb-imax)=xp(ijb)
	  xp(ijt+imax)=xp(ijt)
55	continue

700   continue

	return
	end
	

      subroutine prtplt (n)
c
c ======================================================================
c
c   Purpose -
c     plot and provide formatted writes to paper and film,
c     where    = 1: initial mesh and problem data
c            n = 2: time step and cycle info
c              = 3: field variables
c
c   PRTPLT is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE    NEWCYC
c
c
c   PRTPLT calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c          ripple
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
	dimension vortex(nxy),xy(nxy),xx(nxy),yy(nxy)
c
      go to (20,120,140,145), n
c
c.... prtplt (1) write out mesh data
c
   20 write (13,*)
      go to 9999
c
c.... prtplt (2)  write time step, cycle information
c

  120 write (13,330) iter,t,idt,itc,jtc,delt,ncyc,vchgt
c
      go to 9999
c
c.... prtplt (3)  write field variables to paper
c
  140	write(13,*)
        if (t+em6.lt.tstart) goto 9999
        if (t-em6.gt.tfinish) goto 9999

        do 170 i=1,imax
        do 170 j=1,jmax
          ij=(j-1)*imax+i
          imj=ij-1
          ijm=ij-imax
	    dij=0.0
	    vortex(ij)=0.0
          if (j.eq.1.or.i.eq.1.or.j.eq.jmax.or.i.eq.imax) go to 160
	    vortex(ij)=(v(ij+1)*at(ij+1)-v(ij)*at(ij))/
     &			(delx(i)+delx(i+1))*2.0
     &			-(u(ij+imax)*ar(ij+imax)-u(ij)*ar(ij))
     &			/(dely(j)+dely(j+1))*2.0
          if (nf(ij).gt.5.or.(nf(ij+1).gt.5.and.nf(ij+1+imax).gt.5).
     &    or.(nf(ij+imax).gt.5.and.nf(ij+1+imax).gt.5)) vortex(ij)=0.0
          if (ac(ij).eq.0..or.(ac(ij+1).eq.0..and.ac(ij+1+imax).eq.0.).
     &    or.(ac(ij+imax).eq.0..and.ac(ij+1+imax).eq.0.)) vortex(ij)=0.
          dij=rdx(i)*(ar(ij)*r(i)*u(ij)-ar(imj)*r(i-1)*u(imj))+
     1         rdy(j)*(at(ij)*ri(i)*v(ij)-at(ijm)*ri(i)*v(ijm))
          if (ac(ij).gt.0.0) dij=dij/(ac(ij)*ri(i))
  160     continue
  170    continue
c
	 do 175 j=jbg,jeg,intery
C	   write(21,550)(f((j-1)*imax+i),i=ibg,ieg,interx)
C	   write(22,550)(u((j-1)*imax+i)*ar((j-1)*imax+i),i=ibg,ieg,
C     &		interx)
C	   write(23,550)(v((j-1)*imax+i)*at((j-1)*imax+i),i=ibg,ieg,
C     &		interx)
C	   write(24,550)(p((j-1)*imax+i),i=ibg,ieg,interx)
C         write(29,560)(vortex((j-1)*imax+i),i=ibg,ieg,interx)
C         write(29,560)(tss((j-1)*imax+i),i=ibg,ieg,interx)
	   if (kemodel.gt.0) then
C 	   write(27,550)(sqrt(2.0*xk((j-1)*imax+i)),
C     &		i=ibg,ieg,interx)
C         write(28,550)(xep((j-1)*imax+i)*1000.0,i=ibg,ieg,interx)
C	   write(30,550)(xnut((j-1)*imax+i)*1000.0,i=ibg,ieg,interx)
C         write(30,550)(ttots((j-1)*imax+i)*1.0,i=ibg,ieg,interx)
	   endif
C         if (npollutant.eq.1)
C     &          write(33,550)(xp((j-1)*imax+i),i=1,imax)
175	 continue
	 goto 9999

145	continue
        if (t+em6.lt.tstart_a) goto 9999
        if (t-em6.gt.tfinish_a) goto 9999

        do 180 i=2,imax-1
        do 180 j=2,jmax-1
          ij=(j-1)*imax+i
          imj=ij-1
          ijm=ij-imax
          vortex(ij)=(v(ij+1)*at(ij+1)-v(ij)*at(ij))/
     &                  (delx(i)+delx(i+1))*2.0
     &                  -(u(ij+imax)*ar(ij+imax)-u(ij)*ar(ij))
     &                  /(dely(j)+dely(j+1))*2.0
          if (nf(ij).gt.5.or.(nf(ij+1).gt.5.and.nf(ij+1+imax).gt.5).
     &    or.(nf(ij+imax).gt.5.and.nf(ij+1+imax).gt.5)) vortex(ij)=0.0
          if (ac(ij).eq.0..or.(ac(ij+1).eq.0..and.ac(ij+1+imax).eq.0.).
     &    or.(ac(ij+imax).eq.0..and.ac(ij+1+imax).eq.0.)) vortex(ij)=0.
          dij=rdx(i)*(ar(ij)*r(i)*u(ij)-ar(imj)*r(i-1)*u(imj))+
     1         rdy(j)*(at(ij)*ri(i)*v(ij)-at(ijm)*ri(i)*v(ijm))
          if (ac(ij).gt.0.0) dij=dij/(ac(ij)*ri(i))
180	continue

          nanima=(nanimation-1)/1000
          nanimb=(nanimation-nanima*1000-1)/100
          nanimc=(nanimation-nanima*1000-nanimb*100-1)/10
          nanimd=nanimation-nanima*1000-nanimb*100-nanimc*10-1
          open(1000,file='f'//char(nanima+48)//char(nanimb+48)//
     &		char(nanimc+48)//char(nanimd+48),status="unknown")
          open(1001,file='k'//char(nanima+48)//char(nanimb+48)//
     &          char(nanimc+48)//char(nanimd+48),status="unknown")
          open(1002,file='p'//char(nanima+48)//char(nanimb+48)//
     &          char(nanimc+48)//char(nanimd+48),status="unknown")
	    if (npollutant.eq.1) 
     &    open(1003,file='c'//char(nanima+48)//char(nanimb+48)//
     &          char(nanimc+48)//char(nanimd+48),status="unknown")
          open(1004,file='u'//char(nanima+48)//char(nanimb+48)//
     &          char(nanimc+48)//char(nanimd+48),status="unknown")
          open(1005,file='w'//char(nanima+48)//char(nanimb+48)//
     &          char(nanimc+48)//char(nanimd+48),status="unknown")
C          open(1006,file='d'//char(nanima+48)//char(nanimb+48)//
C     &          char(nanimc+48)//char(nanimd+48),status="unknown")
          open(1007,file='v'//char(nanima+48)//char(nanimb+48)//
     &          char(nanimc+48)//char(nanimd+48),status="unknown")
C		open(1008,file='z'//char(nanima+48)//char(nanimb+48)//
C     &		char(nanimc+48)//char(nanimd+48),status="unknown")
     
          do 185 j=jbg+1,jeg-1,intery
	    if (kr.eq.5) then
              write(1000,550)(xnut((j-1)*imax+i),
     &		i=ibg+1,ieg-1,interx)
		else
			write(1000,550)(f((j-1)*imax+i),i=ibg+1,ieg-1,interx)
	    endif
	    write(1001,590)(sqrt(2.0*xk((j-1)*imax+i)),
     &          i=ibg+1,ieg-1,interx)
          write(1002,560)(p((j-1)*imax+i),i=ibg+1,ieg-1,interx)
C            write(1002,560)(vortex((j-1)*imax+i),i=ibg+1,ieg-1,interx)
          if (npollutant.eq.1)
     &      write(1003,550)(xp((j-1)*imax+i),i=ibg+1,ieg-1,interx)
          write(1004,590)(cvmgt((u((j-1)*imax+i)*ar((j-1)*imax+i)
     &		+u((j-1)*imax+i-1)*ar((j-1)*imax+i-1))/2,
     &          0.0d0,f((j-1)*imax+i).gt.0.0),i=ibg+1,ieg-1,interx)
          write(1005,590)(cvmgt((v((j-1)*imax+i)*at((j-1)*imax+i)
     &		+v((j-2)*imax+i)*at((j-2)*imax+i))/2,
     &          0.0d0,f((j-1)*imax+i).gt.0.0),i=ibg+1,ieg-1,interx)
C            write(1006,590)((dissipturb((j-1)*imax+i)
C     &		+dissipmole((j-1)*imax+i)),i=ibg+1,ieg-1,interx)
CC            write(1007,590)((
CC     &          vortex((j-1)*imax+i)),i=ibg+1,ieg-1,interx)
            write(1007,590)((
     &          ac((j-1)*imax+i)),i=ibg+1,ieg-1,interx)
C		  write(1008,590)((
C     &		  productturb((j-1)*imax+i)),i=ibg+1,ieg-1,interx)
185	  continue
	  nanimation=nanimation+1

	  close(1000)
	  close(1001)
	  close(1002)
	  if (npollutant.eq.1) close(1003)
	  close(1004)
	  close(1005)
C	  close(1006)
C	  close(1007)
C	  close(1008)	
c
 9999 return
c
  270 format (1h1,1x,43hnf field (incl. fictitious cells) for cycle,i6,
     1         2x,2ht=,1pe14.7,2x,5hdelt=,e12.5//5x,32i4)
  275 format (1h1,1x,43hnw field (incl. fictitious cells) for cycle,i6,
     1         2x,2ht=,1pe14.7,2x,5hdelt=,e12.5//5x,32i4)
  280 format (1x,i3,1x,63i2/(5x,63i2))
  290 format (1h1)
  300 format (a80)
  310 format (4x,1hi,5x,1hj,9x,1hu,14x,1hv,15x,1hp,15x,1hd,12x,
     1         1hf,11x,2hnf)
  320 format (2x,i3,3x,i3,6(3x,1pe12.5),3x,i3)
  330 format (1x," iter = ",i5," time = ",1pe12.5," dt",a2,"(",
     &        i4,",",i4,") = ",1pe12.5," cycle = ",i6," vchgt = ",
     &        1pe12.5)
  350 format (1h0)
  360 format (1h ,18x,a80,1x,a10,2(1x,a8))
  370 format (2x,5hnkx= ,i4)
  380 format (2x,8hmesh-x= ,i4,3x,4hxl= ,1pe12.5,3x,4hxc= ,e12.5,3x,
     1       4hxr= ,e12.5,3x,5hnxl= ,i4,3x,5hnxr= ,i4,3x,6hdxmn= ,e12.5)
  390 format (2x,8hmesh-y= ,i4,3x,4hyl= ,1pe12.5,3x,4hyc= ,e12.5,3x,
     1       4hyr= ,e12.5,3x,5hnyl= ,i4,3x,5hnyr= ,i4,3x,6hdymn= ,e12.5)
  400 format (2x,5hnky= ,i4)
  430 format (13x,i3,2x,1pe12.5,3x,e12.5)

  550 format (3000f8.3)
  560 format (3000f9.3)
  570 format (3000f9.4)
  590 format (300e12.4)
	end

 
      subroutine pset(ntype)
c
c ======================================================================
c
c   Purpose -
c     Initialize any interior porous media with the function f(x,y);
c     for f(x,y) < 0, the point (x,y) is either open or closed
c     to the flow for ipr=0 or 1, respectively
c
c   ASET is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c           SETUP
c
c
c   PSET calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c                    ripple   SETARRY    ripple       FXY    ripple
c
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
      dimension iflg(5),dis(4),xm(5),ym(5)
c
      data zero /0.0d0/
	dimension npctmp(nxy),actmp(nxy)
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
C.....note sometimes the boundary cell could be miss-flaged; increase the
C.....border a little (0.01%) bit larger to solve the problem; bug in future
      do 420 j=1,jmax
      do 420 i=1,imax
		ij=(j-1)*imax+i
		npctmp(ij)=npc(ij)	
		actmp(ij)=1.0
420	continue

      xmode=2.*pi/(x(im1)-x(1))
      ymode=2.*pi/(y(jm1)-y(1))
      rad2o2=sqrt(2.0d0)/2.0d0
c
c     -------------------------------------
      call setarry (a(1),0.0d0,nxy)
c     -------------------------------------
c
      do 230 k=(ntype-1)*20+1,(ntype-1)*20+nporous(ntype)
        do 221 j=1,jmax
          do 220 i=1,imax
            ij=(j-1)*imax+i
	    if (ac(ij).le.0.5.or.npc(ij).ne.1) goto 220
            rdxdy=1.0/(delx(i)*dely(j))
            do 60 m=1,4
              go to (10,20,30,40), m
   10         x1=x(i)
              y1=cvmgt(y(j)-dely(j),y(j-1),j.eq.1)
              dis(1)=dely(j)
              go to 50
   20         y1=y(j)
              x1=x(i)
              dis(2)=delx(i)
              go to 50
   30         x1=cvmgt(x(i)-delx(i),x(i-1),i.eq.1)
              y1=y(j)
              dis(3)=dely(j)
              go to 50
   40         y1=cvmgt(y(j)-dely(j),y(j-1),j.eq.1)
              x1=cvmgt(x(i)-delx(i),x(i-1),i.eq.1)
              dis(4)=delx(i)
   50         iflg(m)=0
              fconic=fxy(x1,y1,pa1(k),pa2(k),pb1(k),pb2(k),pc1(k),
     &                   pc2(k),pd1(k),nxp(k)*xmode,pd2(k),mxp(k)*xmode,
     &                   pe1(k),nyp(k)*ymode,pe2(k),myp(k)*ymode)
              if (abs(fconic).le.1.0d-12) fconic=0.0
              if (fconic.le.0.0) iflg(m)=1
              xm(m)=x1
              ym(m)=y1
   60       continue
c
            iflg(5)=iflg(1)
            xm(5)=xm(1)
            ym(5)=ym(1)
            iflgs=0
c
            do 70 m=1,4
              iflgs=iflgs+iflg(m)
   70       continue
c
            brij=0.0
            btij=0.0
            if (iflgs.eq.0) go to 220
            if (iflgs.lt.4) go to 80
            bij=1.0d0
            brij=1.0d0
            btij=1.0d0
            go to 200
   80       if (iflg(1).eq.1.and.iflg(2).eq.1) brij=1.0d0
            if (iflg(2).eq.1.and.iflg(3).eq.1) btij=1.0d0
c
            do 160 m=1,4
              if (iflg(m).eq.iflg(m+1)) go to 160
              x1=xm(m)
              y1=ym(m)
              x2=xm(m+1)
              y2=ym(m+1)
              if (iflg(m).eq.0) go to 90
              x2=xm(m)
              y2=ym(m)
              x1=xm(m+1)
              y1=ym(m+1)
   90         epsif=0.001d0*(abs(x2-x1)+abs(y2-y1))
              smn=0.0
              fmn=fxy(x2,y2,pa1(k),pa2(k),pb1(k),pb2(k),pc1(k),
     &                pc2(k),pd1(k),nxp(k)*xmode,pd2(k),mxp(k)*xmode,
     &                pe1(k),nyp(k)*ymode,pe2(k),myp(k)*ymode)
              smx=1.0d0
              fmx=fxy(x1,y1,pa1(k),pa2(k),pb1(k),pb2(k),pc1(k),
     &                pc2(k),pd1(k),nxp(k)*xmode,pd2(k),mxp(k)*xmode,
     &                pe1(k),nyp(k)*ymode,pe2(k),myp(k)*ymode)
              s=0.5
  100         xt=s*x1+(1.0-s)*x2
              yt=s*y1+(1.0-s)*y2
              fs=fxy(xt,yt,pa1(k),pa2(k),pb1(k),pb2(k),pc1(k),
     &                pc2(k),pd1(k),nxp(k)*xmode,pd2(k),mxp(k)*xmode,
     &                pe1(k),nyp(k)*ymode,pe2(k),myp(k)*ymode)
              if (abs(fs).lt.epsif) go to 130
              if (fs.ge.0.0) go to 110
              fden=abs(fs-fmn)+1.0d-10
              se=s-fs*(s-smn)/fden
              if (se.gt.smx) se=smx
              fmn=fs
              smn=s
              go to 120
  110         fden=abs(fmx-fs)+1.0d-10
              se=s-fs*(smx-s)/fden
              if (se.lt.smn) se=smn
              fmx=fs
              smx=s
  120         si=s-fs*(smx-smn)/(fmx-fmn)
              s=0.5*(se+si)
              go to 100
  130         dis(m)=sqrt((xt-x2)**2+(yt-y2)**2)
              go to (140,150,160,160), m
  140         brij=dis(1)/dely(j)
              go to 160
  150         btij=dis(2)/delx(i)
  160       continue
c
            m=0
            bij=0.0
  170       continue
            m=m+1
            if (m.eq.5) go to 190
            if (iflg(m).eq.0) go to 170
            mp1=m+1
            if (mp1.eq.5) mp1=1
            mm1=m-1
            if (mm1.eq.0) mm1=4
            bij=bij+dis(m)*dis(mm1)
            if (iflg(mp1).eq.1) go to 180
            dis2=dis(m)
  180       continue
            if (iflg(mm1).eq.1) go to 170
            dis1=dis(mm1)
            go to 170
  190       continue
            if (iflgs.eq.3) bij=bij-dis1*dis2
            bij=0.5*bij*rdxdy
            if (bij.gt.1.0) bij=1.0d0
  200       continue
            if (ipr(k).eq.0) go to 210
            bij=-bij
            brij=-brij
            btij=-btij
  210       actmp(ij)=actmp(ij)+bij
            if(actmp(ij).gt.0.9999d0) actmp(ij)=1.0d0
            if(actmp(ij).lt.0.0001d0) actmp(ij)=0.0d0
            if(actmp(ij).lt.0.5d0) then
		npctmp(ij)=ntype+1
	    else
		npctmp(ij)=npc(ij)
	    endif
  220     continue
  221   continue
  230 continue

      do 450 n=1,nxy
		npc(n)=npctmp(n)
450   continue
c
C.....reset ac, ar, at = 0.0 in the vicinity of porous cells
C.....but do not correct the cell who is not inside porous media (29/09/02)
      do i=2,im1
        do j=2,jm1
          ij=(j-1)*imax+i
          ijm=ij-imax
          imj=ij-1
          ipj=ij+1
          ijp=ij+imax
          if (npc(ij).gt.1) then
            if (ac(ipj).gt.0.5.and.ac(ipj).lt.1.0) then
			if (npc(ij+1+imax).gt.1.or.npc(ij+1-imax).gt.1) then
                ac(ipj)=1.0d0
                at(ipj)=1.0
                ar(ipj)=1.0
			endif
            endif
            if (ac(ipj).gt.0.0.and.ac(ipj).le.0.5) then
			if (npc(ij+1+imax).gt.1.or.npc(ij+1-imax).gt.1) then
                ac(ipj)=0.0   
                at(ipj)=0.0   
                ar(ipj)=0.0   
                ar(ij)=0.0
                at(ij+1-imax)=0.0
			endif
            endif

            if (ac(imj).gt.0.5.and.ac(imj).lt.1.0) then
              if (npc(ij-1+imax).gt.1.or.npc(ij-1-imax).gt.1) then 
			  ac(imj)=1.0d0
                at(imj)=1.0
                ar(imj)=1.0
			endif
            endif
            if (ac(imj).gt.0.0.and.ac(imj).le.0.5) then
              if (npc(ij-1+imax).gt.1.or.npc(ij-1-imax).gt.1) then
			  ac(imj)=0.0
                at(imj)=0.0   
                ar(imj)=0.0   
                ar(ij-2)=0.0   
                at(ij-1-imax)=0.0
			endif
            endif
            
            if (ac(ijp).gt.0.5.and.ac(ijp).lt.1.0) then
              if (npc(ij+imax+1).gt.1.or.npc(ij+imax-1).gt.1) then  
			  ac(ijp)=1.0d0
                at(ijp)=1.0
                ar(ijp)=1.0
			endif
            endif
            if (ac(ijp).gt.0.0.and.ac(ijp).le.0.5) then
              if (npc(ij+imax+1).gt.1.or.npc(ij+imax-1).gt.1) then
			  ac(ijp)=0.0
                at(ijp)=0.0
                ar(ijp)=0.0
                ar(ij+imax-1)=0.0
                at(ij)=0.0
			endif
            endif

            if (ac(ijm).gt.0.5.and.ac(ijm).lt.1.0) then
              if (npc(ij-imax+1).gt.1.or.npc(ij-imax-1).gt.1) then
			  ac(ijm)=1.0d0
                at(ijm)=1.0
                ar(ijm)=1.0
			endif      
		  endif
            if (ac(ijm).gt.0.0.and.ac(ijm).le.0.5) then
              if (npc(ij-imax+1).gt.1.or.npc(ij-imax-1).gt.1) then
			  ac(ijm)=0.0
                at(ijm)=0.0
                ar(ijm)=0.0
                ar(ij-imax-1)=0.0
                at(ij-2*imax)=0.0
			endif
            endif
          endif
        end do
      end do
c
      do i=2,im1
          do j=2,jm1
            ij=(j-1)*imax+i
            if (ac(ij).eq.0.0) then
                ar(ij)=0.0
                at(ij)=0.0
                ar(ij-1)=0.0
                at(ij-imax)=0.0
            endif 
          end do
      end do
c
      return
      end
 

	subroutine cnoidal(height,depth,period,l,c,x2,mod)

	implicit none

	real*8 mod,ck,c,l,ce,cn,x,xa,xb,height,depth,period,grav
	real*8 moda,modb,x2,xtemp
	integer n

	external ck,ce,cn

	grav=9.8

	mod=0.75
	moda=0.5
	modb=1.0

	n=0

40	n=n+1	

	x=period-dsqrt(depth/grav)*4.0*ck(mod)*dsqrt(mod*depth/3.0/
     &	height)/dsqrt(1.0+height/depth/mod*(-mod+2.0-3.0*ce(mod)/
     &	ck(mod)))
c	write(1,*)'~',mod,x

	xa=period-dsqrt(depth/grav)*4.0*ck(moda)*dsqrt(moda*depth/3.0/
     &	height)/dsqrt(1.0+height/depth/moda*(-moda+2.0-3.0*ce(moda)/
     &	ck(moda)))

	xb=period-dsqrt(depth/grav)*4.0*ck(modb)*dsqrt(modb*depth/3.0/
     &	height)/dsqrt(1.0+height/depth/modb*(-modb+2.0-3.0*ce(modb)/
     &	ck(modb)))

	if (abs(x).lt.1.0e-16.or.n.gt.100) goto 50

	if (x*xa.lt.0.0) then 

	   modb=mod
	   mod=(moda+mod)/2.0
	else
	   moda=mod
	   mod=(mod+modb)/2.0
	endif

	goto 40

50	continue
C.......after sobet el at (1987, J. Waterway)
	l=4.0*ck(mod)*depth*dsqrt(mod*depth/height/3.0)
C.......after mei (1983) or simply c=L/T
	xtemp=-mod+2.0-3.0*ce(mod)/ck(mod)
        c=sqrt(grav*depth*(1.0+height/depth/mod*xtemp))
	x2=height/mod*(1.0-mod-ce(mod)/ck(mod))

	return
	end

c	..................................................................
	
	function ck(mod)
	real*8 ck,mod,mod1,a0,a1,b0,b1,c0,c1

	mod1=1.0-mod
	ck=0.0
	a0=1.0
	b0=dsqrt(mod1)
	c0=dsqrt(mod)
	n=1

15	if (abs(c0).lt.1.0e-15.or.n.gt.1000) then
	   goto 20
	else
	   n=n+1
	   a1=(a0+b0)/2.0
	   b1=dsqrt(a0*b0)
	   c1=(a0-b0)/2.0

	   a0=a1
	   b0=b1
	   c0=c1

	   goto 15

	endif

20	ck=3.1415926535897932384626/2.0/a0

	return
	end

c	................................................................

	function ce(mod)
	
	real*8 sum,mod,mod1,a0,a1,b0,b1,c0,c1,ck,c(1000),ce

	ce=0.0
	mod1=1.0-mod
	a0=1.0
	b0=dsqrt(mod1)
	c0=dsqrt(mod)
	n=1
	c(n)=c0

15	if (abs(c0).lt.1.0e-15.or.n.gt.1000) then
	   goto 20
	else
	   n=n+1
	   a1=(a0+b0)/2.0
	   b1=dsqrt(a0*b0)
	   c1=(a0-b0)/2.0

	   a0=a1
	   b0=b1
	   c0=c1
	   c(n)=c0

	   goto 15

	endif

20	ck=3.1415926535897932384626/2.0/a0

	sum=0.0
	do 30 k1=1,n
	   sum=sum+2.0**(k1-2)*c(k1)**2
30	continue
	
	ce=ck*(1.0-sum)

	return
	end

c	..............................................................

	function cn(u,mod)

	real*8 mod,mod1,u,a0,a1,b0,b1,c0,c1,c(1000),a(1000),y(1000),cn

	mod1=1.0-mod	
	a0=1.0
	b0=dsqrt(mod1)
	c0=dsqrt(mod)
	n=1
	a(n)=a0
	c(n)=c0

15	if (abs(c0).lt.1.0e-15.or.n.gt.1000) then
	   goto 20
	else
	   n=n+1
	   a1=(a0+b0)/2.0
	   b1=dsqrt(a0*b0)
	   c1=(a0-b0)/2.0

	   a0=a1
	   b0=b1
	   c0=c1

	   a(n)=a0
	   c(n)=c0

	   goto 15

	endif

20	y(n)=2.0**(n-1)*a(n)*u

	do 10 i=n-1,1,-1
	   y(i)=(y(i+1)+dasin(c(i+1)/a(i+1)*dsin(y(i+1))))/2.0
10	continue

	cn=dcos(y(1))

C.....it is found that in MS Developer Studio, the debug and release
C.....modes give different results; by adding a dummy close, such problem
C.....is gone; should be the bug in the compiler
	close(9991)

	return 
	end


       subroutine rypax(e,d,c,b,a,bu,cu,du,eu,nh,nv,x,y,r)
c
c***purpose
c
c      form the matrix product y=r*y+a[t]*x,where the matrix a is
c      a 9-diagonal matrix associated with a 9-point difference
c      operator.
c
c***description
c
c       rypax requires that the coefficients of the ith
c       equation be stored in the ith element of the diagonal
c       arrays. nine diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c
c     on entry
c
c        e,d,c,b,a,bu,cu,du,eu,   are 1-d arrays of length
c            at least nhxnv. the arrays contain, from left to right,
c            the diagonals of a as one sweeps a left to right.
c
c        nh    is the size  or order of each tridiagonal block if
c              a is viewed as a block tridiagonal system of equations.
c
c        nv    is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c
c        x     is an input vector of dimenison at least nhxnv.
c
c        y     is an input vector of dimension at least nhxnv.
c
c        r     is a scalar multiplier of the input vector y.
c
c     on return
c
c        y     is the output vector of dimension at least
c              nhxnv containing r*y+a[t]*x.
c
c***routines called   none
c
      implicit real*8 (a-h,o-z)
      real*8 e(*),d(*),c(*),b(*),a(*),bu(*),cu(*),du(*),eu(*),
     &       x(*),y(*)
c
c*** first executable statement   rypax
      n = nh*nv
      nhm1 = nh-1
      nhp1 = nh+1
      i = 1
      y(i) = r*y(i)+(a(i)*x(i)+bu(i)*x(i+1)
     1      +du(i)*x(i+nh)+eu(i)*x(i+nhp1))
      do 10 i=2,nh
      y(i) = r*y(i)+(b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1)
     1        +(cu(i)*x(i+nhm1)+(du(i)*x(i+nh)+(eu(i)*x(i+nhp1))
     2        )))))
   10 continue
      i = nhp1
      y(i) = r*y(i)+(d(i)*x(i-nh)+c(i)*x(i-nhm1)
     1        +b(i)*x(i-1)+a(i)*x(i)+bu(i)*x(i+1)
     2        +cu(i)*x(i+nhm1)+du(i)*x(i+nh)+eu(i)*x(i+nhp1))
      do 20 i=nh+2,n-nhp1
      y(i) = r*y(i)+(e(i)*x(i-nhp1)+(d(i)*x(i-nh)+(c(i)*x(i-nhm1)
     1        +(b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1)
     2        +(cu(i)*x(i+nhm1)+(du(i)*x(i+nh)+(eu(i)*x(i+nhp1))
     3        ))))))))
   20 continue
      i=n-nh
      y(i) = r*y(i)+(e(i)*x(i-nhp1)+d(i)*x(i-nh)+c(i)*x(i-nhm1)
     1        +b(i)*x(i-1)+a(i)*x(i)+bu(i)*x(i+1)
     2        +cu(i)*x(i+nh-1)+du(i)*x(i+nh))
      do 30 i=n-nhm1,n-1
      y(i) = r*y(i)+(e(i)*x(i-nhp1)+(d(i)*x(i-nh)+(c(i)*x(i-nhm1)
     1        +(b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1))
     2        )))))
   30 continue
      i = n
      y(i) = r*y(i)+(e(i)*x(i-nhp1)+d(i)*x(i-nh)
     1      +b(i)*x(i-1)+a(i)*x(i))
      return
      end


      subroutine setarry (a,c,n)
c
c ======================================================================
c
c   Purpose -
c     sets the contents of array a of length n equal to c
c
c   SETARRY is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c            ASET      GEOM   INITREG   PLOTOUT    RINPUT
c           SETUP    STRAIN   TENSION    VOFADV    
c
c
c   SETARRY calls the following subroutines and functions -
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
      dimension a(1)
c
c  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      do 10 i=1,n
        a(i)=c
   10 continue
c
      return
      end

      subroutine setarryi (iia,iic,n)
c
c ======================================================================
c
c   Purpose -
c     sets the contents of array a of length n equal to c
c
c   SETARRY is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c            ASET      GEOM   INITREG   PLOTOUT    RINPUT
c           SETUP    STRAIN   TENSION    VOFADV    
c
c
c   SETARRY calls the following subroutines and functions -
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
      dimension iia(1)
c
c  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      do 10 i=1,n
        iia(i)=iic
   10 continue
c
      return
      end
 
      subroutine setnf
c
c ======================================================================
c
c   Purpose -
c     determine the surface cell type nf(i,j)
c
c   SETNF is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c           SETUP    VOFADV
c
c
c   SETNF calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c      
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
C==========================================================================
C.....add the secondary orientation of free
C.....surface, measure by nfs
C==========================================================================

c.... Set the default value of nf
c
      do 100 ij=1,nxy
        nf(ij)=0
	  nfs(ij)=0
  100 continue
c
      do 751 i=2,im1
        do 750 j=2,jm1
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

c
c....     if cell is obstacle cell skip to end of loops
          if (ac(ij).lt.em6) go to 750
c
c....     declare empty cell to be a void cell
          if (f(ij).lt.1.0e-6) then
		  nf(ij)=6
		  nfs(ij)=6
		  f(ij)=0.0
		  goto 750
	    endif
c
c....     if cell is empty (or full but psat = 0.0) skip to end of loops;
c....     cells will have default values
          if ((f(ij).gt.emf1.and.psat.eq.0.0))
     &       go to 750
c
c....     four tests to see whether one of the four neighbor cells is both
c         empty and open to flow from (i,j) cell; if so, enter main do
c         loops thru 190
          if (f(ipj).lt.emf.and.ar(ij).gt.em6) go to 190
          if (f(imj).lt.emf.and.ar(imj).gt.em6) go to 190
          if (f(ijp).lt.emf.and.at(ij).gt.em6) go to 190
          if (f(ijm).lt.emf.and.at(ijm).gt.em6) go to 190
c
c....     cell is not a surface cell, obstacle cell,empty cell or partic-
c         ular type of full cell. if it satisfies pressure test, set nf=5
c         for isolated cell and skip  to end of loops; otherwise cell is
c         fluid cell and we skip to end of loops with default values
          if (p(ij).le.psat*em6p1.and.psat.gt.0.0) nf(ij)=5
          go to 750
c
c....     we now enter calculational parts of main do loops
  190     continue
c
c....     calculate the partial derivatives of f
c
c....     distances from midpoint of cell to midpoint of neighbor cells
c         distance to right and left neighbors
          dxr=0.5*(delx(i)+delx(i+1))
          dxl=0.5*(delx(i)+delx(i-1))
c
c....     distance to top and bottom neighbors
          dyt=0.5*(dely(j)+dely(j+1))
          dyb=0.5*(dely(j)+dely(j-1))
c
c....     denominators for finite difference formulas for partial
c         derivatives in x and y directions
          rxden=1.0/(dxr*dxl*(dxr+dxl))
          ryden=1.0/(dyt*dyb*(dyt+dyb))
c
c....     fofm (and fofp) indicate whether cells with lesser (greater) in-
c         dices contribute to average fluid heights in three cell array
c         fofm=1.0 when cell contributes;=0.0 otherwise
c         index is j for vertical heights; i for horizontal heights
c         obstacle cell does not contribute if no fluid in neighbor cell of
c         three cell array
          fofm=1.0d0
          if (ac(ipjm).le.em6.and.f(ijm).lt.emf) fofm=0.0
          fofp=1.0d0
          if (ac(ipjp).le.em6.and.f(ijp).lt.emf) fofp=0.0
c
          avfr=cvmgt(porousp(ipjm),ac(ipjm),npc(ipjm).ne.1)
     &	        *f(ipjm)*fofm*dely(j-1)+
     1          cvmgt(porousp(ipj),ac(ipj),npc(ipj).ne.1)
     &		*f(ipj)*dely(j)+
     2          cvmgt(porousp(ipjp),ac(ipjp),npc(ipjp).ne.1)
     &		*f(ipjp)*fofp*dely(j+1)
C.........end

          fofm=1.0d0
          if (ac(imjm).le.em6.and.f(ijm).lt.emf) fofm=0.0
          fofp=1.0d0
          if (ac(imjp).le.em6.and.f(ijp).lt.emf) fofp=0.0
c
          avfl=cvmgt(porousp(imjm),ac(imjm),npc(imjm).ne.1)
     &		*f(imjm)*fofm*dely(j-1)+
     1         	cvmgt(porousp(imj),ac(imj),npc(imj).ne.1)
     &		*f(imj)*dely(j)+
     2          cvmgt(porousp(imjp),ac(imjp),npc(imjp).ne.1)
     &		*f(imjp)*fofp*dely(j+1)
C.........end

          fofm=1.0d0
          if (ac(ijm).le.em6.and.f(imjm)+f(ipjm).lt.emf)
     &       fofm=0.0
          fofp=1.0d0
          if (ac(ijp).le.em6.and.f(imjp)+f(ipjp).lt.emf)
     &       fofp=0.0
c
c....     y fluid height in central cells = avfcx
          avfcx=cvmgt(porousp(ijm),ac(ijm),npc(ijm).ne.1)
     &		*f(ijm)*fofm*dely(j-1)+
     1          cvmgt(porousp(ij),ac(ij),npc(ij).ne.1)
     &		*f(ij)*dely(j)+
     2         	cvmgt(porousp(ijp),ac(ijp),npc(ijp).ne.1)
     &		*f(ijp)*fofp*dely(j+1)
C.........end of modification

          fofm=1.0d0
          if (ac(imjp).le.em6.and.f(imj).lt.emf) fofm=0.0
          fofp=1.0d0
          if (ac(ipjp).le.em6.and.f(ipj).lt.emf) fofp=0.0
c
c....     x fluid width in cells above = avft
          avft=cvmgt(porousp(imjp),ac(imjp),npc(imjp).ne.1)
     &		*f(imjp)*fofm*delx(i-1)+
     1          cvmgt(porousp(ijp),ac(ijp),npc(ijp).ne.1)
     &		*f(ijp)*delx(i)+
     2          cvmgt(porousp(ipjp),ac(ipjp),npc(ipjp).ne.1)
     &		*f(ipjp)*fofp*delx(i+1)
C.........end

          fofm=1.0d0
          if (ac(imjm).le.em6.and.f(imj).lt.emf) fofm=0.0
          fofp=1.0d0
          if (ac(ipjm).le.em6.and.f(ipj).lt.emf) fofp=0.0
c
c....     x fluid width in cells below = avfb
          avfb=cvmgt(porousp(imjm),ac(imjm),npc(imjm).ne.1)
     &		*f(imjm)*fofm*delx(i-1)+
     1          cvmgt(porousp(ijm),ac(ijm),npc(ijm).ne.1)
     &		*f(ijm)*delx(i)+
     2          cvmgt(porousp(ipjm),ac(ipjm),npc(ipjm).ne.1)
     &		*f(ipjm)*fofp*delx(i+1)
C.........end

          fofm=1.0d0
          if (ac(imj).le.em6.and.f(imjm)+f(imjp).lt.emf)
     &       fofm=0.0
          fofp=1.0d0
          if (ac(ipj).le.em6.and.f(ipjp)+f(ipjm).lt.emf)
     &    fofp=0.0
c
c....     x fluid width in central cells = avfcy
          avfcy=cvmgt(porousp(imj),ac(imj),npc(imj).ne.1)
     &		*f(imj)*fofm*delx(i-1)+
     1         	cvmgt(porousp(ij),ac(ij),npc(ij).ne.1)
     &		*f(ij)*delx(i)+
     2          cvmgt(porousp(ipj),ac(ipj),npc(ipj).ne.1)
     &		*f(ipj)*fofp*delx(i+1)
C.........end
c
c....     avfl set by convention in first column of cells in x and y directions
c....     obstacles are placed on floor from which distances are measured;
c....     fluid is then above obstacles
c....     if nf = 2 or 4 floor is at top of cell; at i+1 or j+1 respectively
c
c....     three point finite difference formulas for surface slopes with
c         variable mesh sizes; formulas exact for quadratics
c         slope dh/dx for almost horizontal fluid = pfx
          pfx=rxden*((avfr-avfcx)*dxl**2+(avfcx-avfl)*dxr**2)
c
c....     slope dw/dy for almost vertical fluid = pfy
          pfy=ryden*((avft-avfcy)*dyb**2+(avfcy-avfb)*dyt**2)
c
c....     pf = sum of squares of tangents; used as flag to differentiate
c         surface cells from isolated cells (nf=5)
          pf=pfx**2+pfy**2
c
c....     if pf very small, cell is declared isolated (instead of surface)
c         and we continue; otherwise go to 660
CC		what about the surrounding is all fluids ??
CC          if (pf.gt.em10) go to 660
	go to 660
c
c....     set nf(i,j) and p(i,j) for isolated cell; bypass determination
c         of pressure interpolation cell,calculation of surface pressure
  655     nf(ij)=5
	    nfs(ij)=5
c
c....     having set nf for the isolated cell,
c         we now skip to end of main loops
          go to 750
  660     continue
c
c....     for surface cells we pick up calculations of main loops; having
c         determined slopes and some auxiliary quantities,
c         we now get nf
c....     in order to set flags to be used later, we sum the f's in cols.
c         to right and left and in rows above and below the (i,j) cell
c
c         sfim = sum of f's in (i-1) col.
c         sfip = sum of f's in (i+1) col.
c         sfjp = sum of f's in (j+1) row
c         sfjm = sum of f's in (j-1) row
c         sfic = sum of f's in i col.
c         sfjc = sum of f's in j row
c
          sfim=f(imjp)+f(imj)+f(imjm)
          sfic=f(ijp)+f(ij)+f(ijm)
          sfip=f(ipjp)+f(ipj)+f(ipjm)
          sfjp=f(ipjp)+f(ijp)+f(imjp)
          sfjc=f(ipj)+f(ij)+f(imj)
          sfjm=f(ipjm)+f(ijm)+f(imjm)
c
c....     if there is little fluid in three by three array of cells: set
c         cell pressure as for an isolated cell and go to end of main loop
c         flags are initially set = 0

CC          if(sfim+sfic+sfip+sfjm+sfjc+sfjp.lt.0.10) go to 655
C

	if (f(ipj)+f(ijp)+f(imj)+f(ijm).lt.em10) goto 655
C.....	end of modification
c
c....     if any cell face is completely closed to flow: skip row and
c         column tests
          if(ar(ij).lt.em6.or.ar(imj).lt.em6.or.at(ij).lt.em6.
     &         or.at(ijm).lt.em6) go to 670
          iflgx=0
          jflgy=0
c
c....     if either col. to left or col. to right is empty or is full:
c         iflgx = 1
          if (sfim.lt.emf.or.sfim.gt.3.0-emf) iflgx=1
          if (sfip.lt.emf.or.sfip.gt.3.0-emf) iflgx=1
c
c....     if either row above or row below is empty or is full: jflgy = 1
          if (sfjp.lt.emf.or.sfjp.gt.3.0-emf) jflgy=1
          if (sfjm.lt.emf.or.sfjm.gt.3.0-emf) jflgy=1
c
c....     if both flags = 1: continue execution at 670 without intervention
          if (iflgx.eq.1.and.jflgy.eq.1) go to 670
c
c....     if exactly one flag = 1: change the corresponding slope
          if (iflgx.eq.1) pfx=1.0d10*pfx
          if (jflgy.eq.1) pfy=1.0d10*pfy
  670     continue
c
c....     we have concluded slope increases
c
c....     determine the pressure interpolation cell nf
c
c         algorithm guarantees that a neighboring fluid cell (one or more
c         always exist) always lies at floor of the surface cell. the fluid
c         cell at the floor is used as the interpolation neighbor cell in
c         presit or prescr
c
c....     get absolute value of slopes; minimum will determine whether
c         surface has near horizontal or near vertical orientation
          abpfx=abs(pfx)
          abpfy=abs(pfy)
c
c....     set default values of the indices of the interpolation neighbor
c         cell (l,m): (l,m) = (i,j)
          l=i
          m=j
c
c....     if surface is more nearly horizontal: go to 680. otherwise
c         surface is more nearly vertical and we set nf = 2 or 1
          if (abpfy.ge.abpfx) go to 680
c
c....     set minimum slope, l index, and nf consistently with nf = 2
          pfmn=pfy
          nf(ij)=2
          l=i+1
c
c....     compute length intervals needed for evaluating the pressure
c         interpolation
          dmx=delx(i)
          dmin=0.5*(dmx+delx(i+1))
          dnbr=delx(i+1)
c
c....     Sign of larger absolute slope (pfx) determines nf value.
c         if pfx > 0.0 then nf = 2 and we skip to end of nf routine at 690
c         otherwise nf = 1 and we go on, repeating immediately previous
c         steps appropriately for nf = 1
          if (pfx.gt.0.0) go to 672
          nf(ij)=1
          pfmn=-pfy
          l=i-1
          dmx=delx(i)
          dmin=0.5*(dmx+delx(i-1))
          dnbr=delx(i-1)
c

672	  continue

c....     having completed calculations for nf = 1, we skip to end of nf
c         routine holding nf = 1 data
c
          pfmn=-pfx
          nfs(ij)=4
          m=j+1
          dmx=dely(j)
          dmin=0.5*(dmx+dely(j+1))
          dnbr=dely(j+1)

          if (pfy.gt.0.0) go to 690
          nfs(ij)=3
          pfmn=pfx
          m=j-1
          dmx=dely(j)
          dmin=0.5*(dmx+dely(j-1))
          dnbr=dely(j-1)

          go to 690

  680     continue

c
c....     we are now in surface more nearly horizontal case.
c         we start with nf = 4
c         set minimum slope, l index and nf consistently with nf = 4
          pfmn=-pfx
          nf(ij)=4
          m=j+1
          dmx=dely(j)
          dmin=0.5*(dmx+dely(j+1))
          dnbr=dely(j+1)
c
c....     Sign of larger slope (pfy) determines nf value.
c         if pfy > 0.0 then nf = 4 and we skip to end of nf routine
c         otherwise nf = 3 and we go on, repeating immediately previous
c         steps appropriately for nf = 3
          if (pfy.gt.0.0) go to 675
          nf(ij)=3
          pfmn=pfx
          m=j-1
          dmx=dely(j)
          dmin=0.5*(dmx+dely(j-1))
          dnbr=dely(j-1)

675	  continue
c
c....     set minimum slope, l index, and nf consistently with nf = 2
          pfmn=pfy
          nfs(ij)=2
          l=i+1
c
c....     compute length intervals needed for evaluating the pressure
c         interpolation
          dmx=delx(i)
          dmin=0.5*(dmx+delx(i+1))
          dnbr=delx(i+1)
          if (pfx.gt.0.0) go to 690
          nfs(ij)=1
          pfmn=-pfy
          l=i-1
          dmx=delx(i)
          dmin=0.5*(dmx+delx(i-1))
          dnbr=delx(i-1)
c
  690     continue
c
  750   continue
  751 continue
c

C.....sweep the second time to identify the fluid cell with F<1
      do 800 j=2,jm1
	do 800 i=2,im1
	  ij=(j-1)*imax+i
	  if (nf(ij).ge.5) goto 800
	  if (nf(ij+1).ne.6.and.nf(ij-1).ne.6.and.nf(ij+imax).ne.6.
     &		and.nf(ij-imax).ne.6) nf(ij)=0
	  if (nf(ij).eq.0.and.nf(ij+1).eq.6) nf(ij)=1
        if (nf(ij).eq.0.and.nf(ij-1).eq.6) nf(ij)=2
        if (nf(ij).eq.0.and.nf(ij+imax).eq.6) nf(ij)=3
        if (nf(ij).eq.0.and.nf(ij-imax).eq.6) nf(ij)=4
800   continue

      return
      end

 
      subroutine setup
c
c ======================================================================
c
c   Purpose -
c     do the problem setup
c
c   SETUP is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE
c
c
c   SETUP calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c         INITREG    ripple  INITVOFF    ripple      ASET    ripple
c              BC    ripple    RINPUT    ripple   SETARRY    ripple
c           SETNF    ripple     TAPIN    ripple
c           EQUIB    ripple   MESHSET    ripple
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
      

C
	dimension uu(250,100),vv(250,100),hh(250)
	dimension ipol(10),jpol(10)
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... read the input data
c
c     ------------------
      call rinput
c     ------------------
c
c.... If the restart flag is greater than zero,
c     read a restart dump and skip the rest of the setup routine.
c
c.... generate the computational mesh
c
c     -------------
   40 call meshset
c     -------------
c

C.....define output data range
	tstart=0.0
	tfinish=twfin
	ibg=1
	ieg=imax
	jbg=1
	jeg=jmax
	interx=1
	intery=1

C.....added to check the boundary resolution
      write(9,*)'aa/2*sin(k*dx)=',aa/2.0d0*sin(xxk*delx(2)),
     &	'dely/2=',0.5d0*dely(jm1)
      if (aa/2.0d0*sin(xxk*delx(2)).gt.0.5d0*dely(jm1)) then
          write(9,*)'please use eiher smaller delx or larger dely'
          write(9,*)'at inflow boundary to provide stable solution'
          write(9,*)'program will continue; accuracy unknown'
          write(9,*)'expect to have false breaking in constant water'
      endif
c
C.....calculate the strength of the source and sink
      if (ninflow.eq.100) then
          area=(x(isourcee)-x(isources-1))*
     &	(y(jsourcee)-y(jsources-1))
	  if (nsource.eq.34.or.nsource.eq.4) then
	    ssource=aa*c1/area
            write(9,*)'ssource=',ssource
	  endif
          if (nsource.eq.14) then
            ssource=amp1*2.0*c1/area
            write(9,*)'ssource=',ssource
          endif
          if (nsource.eq.44) then
	    do nw=1,nwave
              swave(nw)=aawave(nw)*cwave(nw)/area
              write(9,*)'swave(',nw,'=',swave(nw)
	    end do
          endif
          if (nsource.eq.5) then
            ssource=2.0d0*aa*c1/area
            write(9,*)'ssource=',ssource
          endif
          if (nsource.eq.24) then
            ssource=2.0d0*aa*c1/area
            write(9,*)'ssource=',ssource
          endif
      endif

c.... set default cell array values
c     ----------------------------------------------
C.....set default porous type as 1 (non-porous) and obstacle 0 (open)
      do n=1,imax*jmax
	  npc(n)=1
	  noc(n)=0
      end do
      call setarry (ac(1),1.0d0,imax*jmax)
      call setarry (ar(1),1.0d0,imax*jmax)
      call setarry (at(1),1.0d0,imax*jmax)
	call setarryi (nmovbd(1),10,imax*jmax)
      call setarry (f(imax+1),1.0d0,imax*(jmax-1))
      call setarry (dissipturb(1),0.0d0,imax*jmax)
      call setarry (dissipmole(1),0.0d0,imax*jmax)
	call setarry (productturb(1),0.0d0,imax*jmax)
c     ----------------------------------------------
c
c.... generate any interior obstacles or special boundary conditions
c.... set the corresponding values of the flow arrays ar,ac,at
c
      if (nobstype.ne.0) then
	  do n=1,nobstype
c     ----------
          call aset(n)
c     ----------
        end do
	  open (25,file='obs',status="unknown")
	  do j=1,jmax
      	write(25,550)(ac((j-1)*imax+i),i=1,imax)
	  end do
	  close(25)
	endif
c
c.... generate porous media
c
      if (npor.eq.1) then
	  do 700 n=1,nportype
c     ----------
		call pset(n)
c     ----------
700	  continue
	  do n=1,nxy
		porousa(n)=xa(npc(n))
		porousb(n)=xxb(npc(n))
		porousc(n)=gc(npc(n))
		porousd(n)=d50(npc(n))
		porousp(n)=xporosity(npc(n))
	  end do
	  open (26,file='porous',status="unknown")
	  do j=1,jmax
     	    write(26,550)(porousp((j-1)*imax+i),i=1,imax)
	  end do
	  close(26)
      endif

      if (npor.eq.10) then
	  dtotal=0.0
	  nsec=1
	  xdepth=porheight
	  xstart=porstart
	  xend=porstart+porlength(nsec)
        do i=1,imax
		if (xi(i).gt.xend) then
		  nsec=nsec+1
	      xstart=xend
		  if (nsec.gt.nporsection) then
			xend=1.0e16
			xdepth=-1.0e16
			porslope(nsec)=0.0
		  else
		    xend=xend+porlength(nsec)
	        xdepth=dtotal
		  endif
	    endif
	    if (xi(i).gt.xstart) then
		  dtotal=xdepth+(xi(i)-xstart)*porslope(nsec)
	    else
		  dtotal=-1.0e16
		endif
		do j=1,jmax
		  ij=(j-1)*imax+i
		  if (y(j).le.dtotal) then
			npc(ij)=2
			porousd(ij)=d50(2)
			porousp(ij)=xporosity(2)
              porousa(ij)=xalpha(2)*(1.0-xporosity(2))**2
     &           /xporosity(2)**3*max(1.0d-6,xnu)/abs(gy)/d50(2)**2
              porousb(ij)=xbeta(2)*(1.0-xporosity(2))/xporosity(2)
     &           **3/abs(gy)/d50(2)
              porousc(ij)=(1.0+xgamma(2)*(1.0-xporosity(2))
     &           /xporosity(2))/xporosity(2)
		  else
			if (dtotal.gt.y(j-1)) then
			  fraction=(dtotal-y(j-1))/(y(j)-y(j-1))
			  npc(ij)=2
			  porousd(ij)=d50(2)
			  porousp(ij)=xporosity(2)*fraction+(1-fraction)
                porousa(ij)=xalpha(2)*(1.0-xporosity(2))**2
     &            /xporosity(2)**3*max(1.0d-6,xnu)/abs(gy)/d50(2)**2
                porousb(ij)=xbeta(2)*(1.0-xporosity(2))/xporosity(2)
     &            **3/abs(gy)/d50(2)
                porousc(ij)=(1.0+xgamma(2)*(1.0-xporosity(2))
     &            /xporosity(2))/xporosity(2)
			else
			  npc(ij)=1
			  porousd(ij)=d50(1)
			  porousp(ij)=xporosity(1)
                porousa(ij)=xa(1)
                porousb(ij)=xxb(1)
                porousc(ij)=gc(1)
			endif
		  endif
          end do
        end do
	  open (26,file='porous',status="unknown")
	  do j=1,jmax
     	    write(26,550)(porousp((j-1)*imax+i),i=1,imax)
	  end do
	  close(26)
      endif

c
      if (nfrsrf.gt.0) go to 90
c
c.... Compute the vof functions associated according to fluid height
c
      do 120 i=1,imax
        do 110 j=2,jmax
          jj=j
          ij=(jj-1)*imax+i
          f(ij)=1.0d0
          if (flht.gt.y(jj-1).and.flht.lt.y(jj)) then
            f(ij)=rdy(jj)*(flht-y(jj-1))
          endif
          if (y(jj-1).ge.flht) f(ij)=0.0
  110   continue
        ijb=imax+i
        f(ijb-imax)=f(ijb)
  120 continue
      go to 95        
c
c.... compute initial vof functions characterized
c     by the user-input free-surface functions
c
c     --------------
   90 call initvoff
c     --------------
c
   95 continue
c
c.... initialize region quantities
c
c     -------------
      call initreg
c     -------------
c
C.....setup the k-eps model's constants
      if (kemodel.gt.0) then
        	c1e=1.44
        	c2e=1.92
        	c_mu=0.09
        	sege=1.3
        	segk=1.0
	    C_mu_bc=0.09

C.......	setup the initial condition for k and eps
C.......	low initial turbulent specification
C.......	in order to maintain the stability, we require when
C.......	initial xk is small, both xnut and xep should be
C.......	proportionally small; The original relation is obtained
C.......	for umean=O(1 m/s); if we have other situation (umean
C.......	large for real scale) or we deliberately require
C.......	initial k is smaller, we use ratio=a*umean to adjust xnut
C.......	as well.
		umean=max(1.0d0,sqrt(-gy*aa**2/(1.0e-6+h0+aa)))
		ratio=ticf
		uturb=0.0025*umean*ratio
		do 35 j=1,jmax
		  do 35 i=1,imax
		   	ij=(j-1)*imax+i
	  	   	if (ac(ij).le.em6) then
	    		  xk(ij)=0.0
	    		  xep(ij)=0.0
	    		  xnut(ij)=0.0
	    		  goto 35
	  	   	endif
	  	    xk(ij)=0.5*uturb**2
	  	    xnut(ij)=5.00*xnu*ratio**2
	  	    xep(ij)=c_mu*xk(ij)**2/xnut(ij)
35		continue
      endif	  

      ibcflg=1
c     ----------------
      call bc
c     ----------------
      ibcflg=0

c.... Set the nf flag
c
c     ----------------
      call setnf
c     ----------------
c
C.....output location
      if (nloc.ne.0) then
          kout=1
          if (xout(kout).lt.0.0) then
                nout(kout)=1
                kout=kout+1
                if (kout.gt.nloc) goto 89
          endif
          do 88 i=1,imax-1
            if (xout(kout).ge.x(i).and.xout(kout).lt.x(i+1)) then
                nout(kout)=i+1
                kout=kout+1
                if (kout.gt.nloc) goto 89
            endif
88        continue
89        continue
      endif

      if (npollutant.ne.0) then
           do 84 n=1,npollutant
             do 86 i=2,imax
86           if (xpol(n).ge.xi(i).and.xpol(n).lt.xi(i+1)) ipol(n)=i
             do 87 j=2,jmax
87           if (ypol(n).ge.yj(j).and.ypol(n).lt.yj(j+1)) jpol(n)=j
             xp(jpol(n)*imax+ipol(n))=conc(n)
84         continue
           write(9,*)'&&',xpol(n),ipol(n),ypol(n),jpol(n),
     &          xp(jpol(n)*imax+ipol(n))
      endif

C.....call Youngs' initial if nfree=4
      if (nfree.eq.4) then
          call input
      endif

 9998 continue
c
 9999 return
  300 format (1x,2hk=,1pe12.4,2x,3hxi=,e12.4,2x,4hper=,e12.4)
  330 format (a80)
  340 format (10x,i10)
  350 format (10x,e20.6)
  370 format (10x,7i10)
  380 format (1x,4hcon=,1pe10.3,1x,11hand fcvlim=,e10.3,1x,
     &         41hare incompatible. setting fcvlim=1.3*con.)
  420 format (1h1)
  440 format (2x,i3,3x,i3,6(3x,1pe12.5))
  550 format (3000f8.3)
      end


       subroutine slve9(n,m,e,d,c,b,a,bu,cu,du,eu,x,y)
c
c***purpose
c
c      solve (incomplete) ax=(l*lt)x=y, where the matrix a is
c      a 9-diagonal matrix associated with a 9-point difference
c      operator. factors are generated lu9p.
c
c***description
c
c       ldlt requires that the coefficients of the ith
c       equation be stored in the ith element of the diagonal
c       arrays. nine diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c       the incomplete factors of a are stored over the input matrix.
c
c     on entry
c
c        e,d,c,b,a,bu,cu,du,eu,   are 1-d arrays of length
c            at least nxm. the arrays contain, from left to right,
c            the diagonals of a as one sweeps a left to right.
c
c        n     is the size  or order of each tridiagonal block if
c              a is viewed as a block tridiagonal system of equations.
c
c        m     is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c
c        y     is an input vector of dimenison at least nxm.
c
c     on return
c
c        x     is the output vector of dimension at least
c              nxm containing the solution vector.
c
c***routines called   none
c
      implicit real*8 (a-h,o-z)
      dimension e(n,m),d(n,m),c(n,m),b(n,m),a(n,m),
     1          bu(n,m),cu(n,m),du(n,m),eu(n,m),x(n,m),y(n,m)
c
ccc   copy y into x
c
      do 1 j=1,m
        do 1 i=1,n
          x(i,j) = y(i,j)
    1 continue
c
ccc   forward solve
c
      do 100 j=1,m
        if (j .ne. 1)    then
        x(1,j) = x(1,j)-c(1,j)*x(2,j-1)-d(1,j)*x(1,j-1)
cdir$ ivdep
        do 20 i=2,n-1
          x(i,j) = x(i,j)-c(i,j)*x(i+1,j-1)
     1           -d(i,j)*x(i,j-1)-e(i,j)*x(i-1,j-1)
   20   continue
        x(n,j) = x(n,j)-d(n,j)*x(n,j-1)-e(n,j)*x(n-1,j-1)
        end if
c
ccc   van der worst approximation
c
cdir$ ivdep
        do 5 i=n,2,-1
          x(i,j) = x(i,j)-b(i,j)*x(i-1,j)
    5   continue
cdir$ ivdep
        do 10 i=n,3,-1
          x(i,j) = x(i,j)+b(i-1,j)*b(i,j)*x(i-2,j)
   10   continue
  100 continue
c
ccc diagonal solve
c
      do 150 j=1,m
        do 150 i=1,n
          x(i,j) = a(i,j)*x(i,j)
  150 continue
c
ccc  backward solve
c
      do 200 j=m,1,-1
        if (j .ne. m)    then
        x(n,j) = x(n,j)-cu(n,j)*x(n-1,j+1)-du(n,j)*x(n,j+1)
cdir$ ivdep
        do 50 i=2,n-1
          x(i,j) = x(i,j)-cu(i,j)*x(i-1,j+1)
     1           -du(i,j)*x(i,j+1)-eu(i,j)*x(i+1,j+1)
   50   continue
        x(1,j) = x(1,j)-du(1,j)*x(1,j+1)-eu(1,j)*x(2,j+1)
        end if
c
ccc van der worst approximation
c
cdir$ ivdep
        do 35 i=1,n-1
          x(i,j) = x(i,j)-bu(i,j)*x(i+1,j)
   35   continue
cdir$ ivdep
        do 40 i=1,n-2
          x(i,j) = x(i,j)+bu(i+1,j)*bu(i,j)*x(i+2,j)
   40   continue
  200 continue
      return
      end
  
 
      subroutine strain(ar,at,ac,x,y,ri,r,u,v,im1,jm1,
     &                  nxy,sxx,syy,sxy,nf,uxmb,vymb,nmovbd)
c
c ======================================================================
c
c   Purpose -
c     compute the strain rate tensor
c
c   STRAIN is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c         IMPLCTP    STRESS
c
c
c   STRAIN calls the following subroutines and functions -
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
      dimension ar(1),at(1),ac(1),x(1),y(1),r(1),
     &          ri(1),u(1),v(1),sxx(1),sxy(1),syy(1),
     &		  uxmb(1),vymb(1),nmovbd(1),nf(1)
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... reset the pertinent arrays
c
c     --------------------------------
      call setarry (sxx(1),0.0d0,nxy)
      call setarry (syy(1),0.0d0,nxy)
      call setarry (sxy(1),0.0d0,nxy)
c     --------------------------------
c
      nxp=im1+1
      do 100 i=1,im1
c
        delxp=x(i+1)-x(i)
	  if (i.eq.1) then
		delx=x(i+1)-x(i)	
	  else
		delx=x(i)-x(i-1)
	  endif
c
        do 100 j=1,jm1
c
          ij=(j-1)*nxp+i
          ipj=ij+1
          ijp=ij+nxp
c      I removed ac(ipjp) that causes array bounds excceeded, Fengyan
          if (ac(ij).lt.1.0d-06.and.ac(ijp).lt.1.0d-06.and.
c     &	ac(ipj).lt.1.0d-06.and.ac(ipjp).lt.1.0d-06) go to 100
     &	ac(ipj).lt.1.0d-06) go to 100
		imj=ij-1
		ijm=ij-nxp
          ipjp=ijp+1
          ipjm=ijm+1
          imjp=ijp-1
          imjm=ijm-1

c - I added the following max -Fengyan
          imj=max(1,imj)
	    ijm=max(1,ijm)
	    imjm=max(1,imjm)
c
          delyp=y(j+1)-y(j)
	    if (j.eq.1) then
		  dely=y(j+1)-y(j)
	    else
		  dely=y(j)-y(j-1)
	    endif
c     I added max(1,imj) for ar and max(1,ijm) for at - Fengyan
          dudxc=((u(ij)-uxmb(nmovbd(ij)))*ar(ij)
     &		-(u(max(1,imj))-uxmb(nmovbd(ij)))*ar(max(1,imj)))/delx
          dvdyc=((v(ij)-vymb(nmovbd(ij)))*at(ij)
     &		-(v(max(1,ijm))-vymb(nmovbd(ij)))*at(max(1,ijm)))/dely
          dudytr=((u(ijp)-uxmb(nmovbd(ijp)))*ar(ijp)
     &		-(u(ij)-uxmb(nmovbd(ij)))*ar(ij))
     &		/(0.5*(dely+delyp))
          dvdxtr=((v(ipj)-vymb(nmovbd(ipj)))*at(ipj)
     &		-(v(ij)-vymb(nmovbd(ij)))*at(ij))
     &		/(0.5*(delx+delxp))

C.........Note that on free surface, the strain of shear and thus shear stress are always zero
C.........This is different from solid, where the strain is not specially treated here; rather
C.........shear stress is calculated specially at stress.f
C		if ((nf(ij).eq.6.and.nf(ipj).eq.6).or.(nf(ijp).eq.6.and.
C     &		nf(ipjp).eq.6)) dudytr=0.0
C		if ((nf(ij).eq.6.and.nf(ijp).eq.6).or.(nf(ipj).eq.6.and.
C     &		nf(ipjp).eq.6)) dvdxtr=0.0      
C.....modified to give the zero tangential gradient on free surface
		if (nf(ijp).eq.6.or.nf(ipjp).eq.6.or.nf(ij).eq.6.
     &		or.nf(ipj).eq.6) then
     			dudytr=0.0
     			dvdxtr=0.0
		endif
		if (nf(ij).eq.6) then
			dudxc=0.0
			dvdyc=0.0
		endif
c
          sxx(ij)=dudxc
          syy(ij)=dvdyc
          sxy(ij)=0.5*(dudytr+dvdxtr)
c
  100 continue
c
      return
      end


	subroutine stokes(aa,h0,xxt,xxk,amp1,amp2,amp3,amp4,
     &		amp5,cnf,gy,pi)

        implicit real*8 (a-h,o-z)

	  dimension b(2,2)

C.......this program calculate the wave amplitude and wave number using 
C.......given wave height (aa), water depth (h0), and wave period, based
C.......fifth-order Stokes wave theory as proposed by Skjelbreia and
C.......Hendrickson (1967) (5-th conf. coastal engng.)

C.......The wave amplitude (amp) and wave number (xxk) are obtained by
C.......solving two nonlinear algebriac equations (nonlinear wave
C.......amplitude equation (f1) and dispersion relationship (f2)).
C.......The Newton's method is used to solve the nonlinear equations
C.......system.
C.......f1=pi*aa/h0-2*pi/(xxk*h0)*(amp+amp**3*b33+amp**5*(b35+b55))
C.......f2=(2*pi*h0)/(-gy*xxt**2)
C     &	-xxk*h0/(2*pi)*tanh(xxk*h0)*(1+amp**2*c1+amp**4*c2)

C.......first step of Newton's method: guess the value of amp and xxk.

	  xxk=2*pi/(xxt*sqrt(-gy*h0))
	  amp=aa/2.0*xxk
	  n=0

C.......second step: form matrix B=partial (f1,f2)/partial (xxk,amp)
C.......(Jacobian matrix)

C.......prepare the most frequently used coefficients


800	   continue
        c=cosh(xxk*h0)
        s=sinh(xxk*h0)
        ck=h0*s
        sk=h0*c
        b33=(3.0*(8.0*c**6+1.0))/(64.0*s**6)
        b35=(88128.0*c**14-208224.0*c**12+70848.0*c**10+54000.0*c**8
     &          -21816.0*c**6+6264.0*c**4-54.0*c**2-81.0)
     &          /(12288.0*s**12*(6.0*c**2-1.0))
        b55=(192000.0*c**16-262720.0*c**14+83680.0*c**12+20160.0*c**10
     &  -7280.0*c**8+7160.0*c**6-1800.0*c**4-1050.0*c**2+225.0)
     &  /(12288.0*s**10*(6.0*c**2-1.0)*(8.0*c**4-11.0*c**2+3.0))

        b33k=9.0*c**5*ck/(4.0*s**6)-
     &		 (9.0*(8.0*c**6+1.0))/(32.0*s**7)*sk
	  b35k=(14.0*88128.0*c**13*ck-12.0*208224.0*c**11*ck
     &		+10.0*70848.0*c**9*ck+8.0*54000.0*c**7*ck
     &          -6.0*21816.0*c**5*ck+4.0*6264.0*c**3*ck
     &		  -2.0*54.0*c**1*ck)
     &          /(12288.0*s**12*(6.0*c**2-1.0))
     &	  -(88128.0*c**14-208224.0*c**12+70848.0*c**10+54000.0*c**8
     &          -21816.0*c**6+6264.0*c**4-54.0*c**2-81.0)*12.0
     &          /(12288.0*s**13*(6.0*c**2-1.0))*sk
     &	  -(88128.0*c**14-208224.0*c**12+70848.0*c**10+54000.0*c**8
     &          -21816.0*c**6+6264.0*c**4-54.0*c**2-81.0)*12.0*c*ck
     &          /(12288.0*s**12*(6.0*c**2-1.0)**2)
	  b55k=(16.0*192000.0*c**15*ck-14.0*262720.0*c**13*ck
     &		+12.0*83680.0*c**11*ck+10.0*20160.0*c**9*ck
     &  	-8.0*7280.0*c**7*ck+6.0*7160.0*c**5*ck
     &		-4.0*1800.0*c**3*ck-2.0*1050.0*c**1*ck)
     &  	/(12288.0*s**10*(6.0*c**2-1.0)*(8.0*c**4-11.0*c**2+3.0))
     &	  -(192000.0*c**16-262720.0*c**14+83680.0*c**12+20160.0*c**10
     &      -7280.0*c**8+7160.0*c**6-1800.0*c**4-1050.0*c**2+225.0)*10.0
     &      /(12288.0*s**11*(6.0*c**2-1.0)*(8.0*c**4-11.0*c**2+3.0))*sk
     &    -(192000.0*c**16-262720.0*c**14+83680.0*c**12+20160.0*c**10
     &      -7280.0*c**8+7160.0*c**6-1800.0*c**4-1050.0*c**2+225.0)
     &      *12.0*c*ck
     &      /(12288.0*s**10*(6.0*c**2-1.0)**2*(8.0*c**4-11.0*c**2+3.0))
     &    -(192000.0*c**16-262720.0*c**14+83680.0*c**12+20160.0*c**10
     &      -7280.0*c**8+7160.0*c**6-1800.0*c**4-1050.0*c**2+225.0)
     &      *(32.0*c**3-22.0*c)*ck
     &      /(12288.0*s**10*(6.0*c**2-1.0)*(8.0*c**4-11.0*c**2+3.0)**2)

	  c1=(8.0*c**4-8.0*c**2+9.0)/(8.0*s**4)
	  c2=(3840.0*c**12-4096.0*c**10+2592.0*c**8-1008.0*c**6+
     &	  5944.0*c**4-1830.0*c**2+147.0)/(512.0*s**10*(6.0*c**2-1.0))

	  c1k=(4.0*8.0*c**3*ck-2.0*8.0*c**1*ck)/(8.0*s**4)
     &	  -(8.0*c**4-8.0*c**2+9.0)*4.0*sk/(8.0*s**5)
	  c2k=(12.0*3840.0*c**11*ck-10.0*4096.0*c**9*ck
     &		+8.0*2592.0*c**7*ck-6.0*1008.0*c**5*ck+
     &		4.0*5944.0*c**3*ck-2.0*1830.0*c**1*ck)
     &		/(512.0*s**10*(6.0*c**2-1.0))
     &	  -(3840.0*c**12-4096.0*c**10+2592.0*c**8-1008.0*c**6+
     &		5944.0*c**4-1830.0*c**2+147.0)*10.0*sk/
     &		(512.0*s**11*(6.0*c**2-1.0))
     &	  -(3840.0*c**12-4096.0*c**10+2592.0*c**8-1008.0*c**6+
     &		5944.0*c**4-1830.0*c**2+147.0)*12.0*c*ck
     &		/(512.0*s**10*(6.0*c**2-1.0)**2)

C.......calculate B coefficient
C.......partial f1/partial xxk
	  b(1,1)=2.0*pi/(xxk**2*h0)*(amp+amp**3*b33+amp**5*(b35+b55))
     &		-2.0*pi/(xxk*h0)*(amp**3*b33k+amp**5*(b35k+b55k))
C.......partial f1/partial amp
	  b(1,2)=-2.0*pi/(xxk*h0)*
     &		(1.0+3.0*amp**2*b33+5.0*amp**4*(b35+b55))
C.......partial f2/partial xxk
	  b(2,1)=-h0/(2.0*pi)*tanh(xxk*h0)*(1.0+amp**2*c1+amp**4*c2)
     &		-xxk*h0/(2.0*pi)*(1.0-(tanh(xxk*h0))**2)*h0*
     &			(1.0+amp**2*c1+amp**4*c2)
     &		-xxk*h0/(2.0*pi)*tanh(xxk*h0)
     &			*(amp**2*c1k+amp**4*c2k)
C.......partial f2/partial amp
	  b(2,2)=-xxk*h0/(2.0*pi)*tanh(xxk*h0)
     &		*(2.0*amp*c1+4.0*amp**3*c2)

C.......calculate f1 and f2
	  f1=pi*aa/h0-2.0*pi/(xxk*h0)*(amp+amp**3*b33+amp**5*(b35+b55))
	  f2=(2.0*pi*h0)/(-gy*xxt**2)-xxk*h0/(2.0*pi)*tanh(xxk*h0)
     &		*(1.0+amp**2*c1+amp**4*c2)

C.......judge if criteria are satisfied:
	  if ((abs(f1).lt.1.0e-12.and.abs(f2).lt.1.0e-12).or.n.gt.100)
     &	goto 1000

C.......solve the two-equation system Bx=-F
C....... b(1,1)*xk + b(1,2)*xa = -f1
C....... b(2,1)*xk + b(2,2)*xa = -f2 	   

	  xa=(f1*b(2,1)-f2*b(1,1))/(b(1,1)*b(2,2)-b(2,1)*b(1,2))
	  xk=(f2*b(1,2)-f1*b(2,2))/(b(1,1)*b(2,2)-b(2,1)*b(1,2))

C.......update xxk and amp

	  xxk=xxk+xk
	  amp=amp+xa

	n=n+1
	goto 800

1000	continue
C	write(*,*)'n=',n,'f1=',f1,'f2=',f2,'xxk=',xxk,'amp=',amp
	   
C.......find the other coefficients
	b22=(2.0*c**2+1.0)*c/(4.0*s**3)
	b24=(c*(272.0*c**8-504.0*c**6-192.0*c**4+322.0*c**2+21.0))
     &		/(384.0*s**9)
	b44=(c*(768.0*c**10-448.0*c**8-48.0*c**6+48.0*c**4+106.0*c**2
     &		-21.0))/(384.0*s**9*(6.0*c**2-1.0))
	amp1=amp/xxk
	amp2=(amp**2*b22+amp**4*b24)/xxk
	amp3=(amp**3*b33+amp**5*b35)/xxk
	amp4=amp**4*b44/xxk
	amp5=amp**5*b55/xxk
C	write(2,200)b22,b24,b44
C	write(2,200)amp1,amp2,amp3,amp4,amp5

C.......find the zero using Newton's process
C.......guess the zero point
	cnf=pi/4.0
	m=0
1800	continue
	coef=amp1*sin(pi/2.0-cnf)+amp2*2.0*sin(2.0*(pi/2.0-cnf))
     &	 +amp3*3.0*sin(3.0*(pi/2.0-cnf))+amp4*4.0*sin(4.0*(pi/2.0-cnf))
     &	 +amp5*5.0*sin(5.0*(pi/2.0-cnf))
	f=amp1*cos(pi/2.0-cnf)+amp2*cos(2.0*(pi/2.0-cnf))
     &	 +amp3*cos(3.0*(pi/2.0-cnf))+amp4*cos(4.0*(pi/2.0-cnf))
     &	 +amp5*cos(5.0*(pi/2.0-cnf))
	if (abs(f).lt.1.0e-12.or.m.gt.100) goto 2000
	x=-f/coef
	cnf=cnf+x
	goto 1800

2000	continue
C	write(2,*)'cnf=',cnf
200	format(6f13.4)
	return
	end


      subroutine stress
c
c ======================================================================
c
c   Purpose -
c     compute the viscous stress tensor
c
c   STRESS is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          VTILDE	 K_Epsilon
c
c
c   STRESS calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c          STRAIN    ripple
c
c ======================================================================
c
c##############################################################
      implicit real*8 (a-h,o-z)
c##############################################################
c
c###########
      include  "comdk2.h" 
c###########
c
	dimension sxy(nxy)
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... get the strain and rotation tensors
c
c     -----------------------------------------------------------------
      call strain(ar,at,ac,x,y,ri,r,u,v,im1,jm1,nxy,
     &            tauxx,tauyy,tauxy,nf,uxmb,vymb,nmovbd)
c     -----------------------------------------------------------------
c
c.....In the following computation for turbulent flow outside porous media,
c.....both Reynolds stress (tau**turb) and total stress are calculated
c.....The former is used in k_epsilon to calculate turbulence production
c.....and the latter is for velocity calculation

	do n=1,nxy
	  sxy(n)=tauxy(n)
	end do

      do 100 j=1,jmax
        do 100 i=1,imax

        ij=(j-1)*imax+i
	  ipj=ij+1
	  imj=ij-1
	  ijp=ij+imax
	  ijm=ij-imax

	imj=max(1,imj)
	ijm=max(1,ijm)

	  ipjp=ijp+1
	  imjp=ijp-1
	  ipjm=ijm+1
	  imjm=ijm-1

c I added the above and below max -Fengyan

	imjp=max(1,imjp)
	imjm=max(1,imjm)

C.......For laminar flows
	  if (kemodel.eq.0) then
            tauxx(ij)=2.*(xmu)*tauxx(ij)
            tauyy(ij)=2.*(xmu)*tauyy(ij)
            tauxy(ij)=2.*(xmu)*sxy(ij)
		  goto 100
	  endif 

C.......For turbulent flows
	  if (kemodel.gt.0) then
C.........when kemodel=1 & 4, turbulence production is calculated by the
C.........product of stress and strain. It is noted the stress here
C.........contains Reynolds stress only based on the derivation. The 
C.........inclusion of the molecular stress may contaminate results when
C.........turbulence is low, though it is insignificant for large 
C.........turbulence flow. Thus, it is more accurate to introduce Reynolds
C.........stress (turbulence) and total stress (mean flow) separatedly.

C.........in certain case, very small xep (e.g., 1.0e-160) may occur during
C.........computation; the third power will make it become zero.
C          if (abs(xep(ij)).lt.1.0e-32.and.ac(ij).gt.1.0e-6) then
C              tauxxturb(ij)=0.0
C              tauxyturb(ij)=0.0
C              tauyyturb(ij)=0.0
C              tauxx(ij)=0.0
C              tauxy(ij)=0.0
C              tauyy(ij)=0.0
C              goto 100
C          endif

	    goto (130,130,130,160), kemodel
C.........standard k-e model
130	    continue
C.........only calculate non-solid cell; for solid cell adjacent to fluid
C.........or porous, use boundary condition to calculate bed shear stress

          if (ac(ij)+ac(ipj)+ac(ijp)+ac(ipjp).lt.1.0d-06) go to 100

C.........porous media calculation
	    if (npc(ij).ne.1.and.npor.ne.0) then
            tauxx(ij)=2.*(xmu)*tauxx(ij)
            tauyy(ij)=2.*(xmu)*tauyy(ij)
            tauxy(ij)=2.*(xmu)*sxy(ij)
C...........specify bed shear stress
C...........boundary on left
		  if ((ac(imj).le.1.0e-6.and.ac(imjp).le.1.0e-6).
     &		or.(i.eq.2.and.kl.eq.1)) then
		    tauxy(imj)=2*xmu*sxy(imj)
		  endif
C...........boundary on bottom
		  if ((ac(ijm).le.1.0e-6.and.ac(ipjm).le.1.0e-6).
     &		or.(j.eq.2.and.kb.eq.1)) then
		    tauxy(ijm)=2*xmu*sxy(ijm)
		  endif
	      goto 100
		endif
		 
		c_d=c_mu*xnut(ij)/(xnu*eddycoef)*
     &           (1.0-exp(-eddycoef*xnu/max(1.0d-20,xnut(ij))))
	    c_dr=c_mu*xnut(ipj)/(xnu*eddycoef)*
     &           (1.0-exp(-eddycoef*xnu/max(1.0d-20,xnut(ipj))))
	    c_dt=c_mu*xnut(ijp)/(xnu*eddycoef)*
     &		(1.0-exp(-eddycoef*xnu/max(1.0d-20,xnut(ijp))))
	    c_drt=c_mu*xnut(ipjp)/(xnu*eddycoef)*
     &           (1.0-exp(-eddycoef*xnu/max(1.0d-20,xnut(ipjp))))
          tauxxturb(ij)=2.*(c_d*xk(ij)**2/(xep(ij)+1.0e-16)
     &		*rhof)*tauxx(ij)-2./3.*rhof*xk(ij)
          tauxx(ij)=2.*(xmu+c_d*xk(ij)**2/(xep(ij)+1.0e-16)
     &		*rhof)*tauxx(ij)-2./3.*rhof*xk(ij)
          tauyyturb(ij)=2.*(c_d*xk(ij)**2/(xep(ij)+1.0e-16)
     &		*rhof)*tauyy(ij)-2./3.*rhof*xk(ij)
          tauyy(ij)=2.*(xmu+c_d*xk(ij)**2/(xep(ij)+1.0e-16)
     &		*rhof)*tauyy(ij)-2./3.*rhof*xk(ij)

C.........Calculate tauxy on wall (wall or bed shear stress) by using 
C.........Three cases are treated: 1. computational boundary 
C.........2. solid adjacent to fluid 3. porous adjacent to fluid

C.........boundary on right
		if (((ac(ipj).le.1.0e-6.and.ac(ipjp).le.1.0e-6).
     &		or.(npc(ipj).ne.1.and.npc(ipjp).ne.1.and.npor.ne.0).
     &		or.(i.eq.im1.and.kr.eq.1)).and.at(ij).gt.em6) then
		  tauxy(ij)=rhof*sqrt(c_mu_bc)*(xk(ij)+xk(ijp))/2
     &		*(-sign(1.0d0,v(ij)))+2*xmu*sxy(ij)
		  tauxyturb(ij)=rhof*sqrt(c_mu_bc)*(xk(ij)+xk(ijp))/2
     &		*(-sign(1.0d0,v(ij)))
		  goto 100
		endif
C.........boundary on top
		if (((ac(ijp).le.1.0e-6.and.ac(ipjp).le.1.0e-6).
     &		or.(npc(ijp).ne.1.and.npc(ipjp).ne.1.and.npor.ne.0).
     &		or.(j.eq.jm1.and.kt.eq.1)).and.ar(ij).gt.em6) then
		  tauxy(ij)=rhof*sqrt(c_mu_bc)*(xk(ij)+xk(ipj))/2
     &		*(-sign(1.0d0,u(ij)))+2*xmu*sxy(ij)
		  tauxyturb(ij)=rhof*sqrt(c_mu_bc)*(xk(ij)+xk(ipj))/2
     &		*(-sign(1.0d0,u(ij)))
		  goto 100
		endif
C.........boundary on left
		if (((ac(imj).le.1.0e-6.and.ac(imjp).le.1.0e-6).
     &		or.(npc(imj).ne.1.and.npc(imjp).ne.1.and.npor.ne.0).
     &		or.(i.eq.2.and.kl.eq.1)).and.at(ij).gt.em6) then
		  tauxy(imj)=rhof*sqrt(c_mu_bc)*(xk(ij)+xk(ijp))/2
     &		*(sign(1.0d0,v(ij)))+2*xmu*sxy(imj)
		  tauxyturb(imj)=rhof*sqrt(c_mu_bc)*(xk(ij)+xk(ijp))/2
     &		*(sign(1.0d0,v(ij)))
		endif
C.........boundary on bottom
		if (((ac(ijm).le.1.0e-6.and.ac(ipjm).le.1.0e-6).
     &		or.(npc(ijm).ne.1.and.npc(ipjm).ne.1.and.npor.ne.0).
     &		or.(j.eq.2.and.kb.eq.1)).and.ar(ij).gt.em6) then
		  tauxy(ijm)=rhof*sqrt(c_mu_bc)*(xk(ij)+xk(ipj))/2
     &		*(sign(1.0d0,u(ij)))+2*xmu*sxy(ijm)
		  tauxyturb(ijm)=rhof*sqrt(c_mu_bc)*(xk(ij)+xk(ipj))/2
     &		*(sign(1.0d0,u(ij)))
		endif

C	    if (xep(ij).eq.0.0.or.xep(ipj).eq.0.0.or.xep(ijp).
C     &		eq.0.0.or.xep(ipjp).eq.0.0) then
C            	tauxyturb(ij)=2.0*(c_d*xk(ij)**2/xep(ij)*rhof)*sxy(ij)
C              tauxy(ij)=2.0*(xmu+c_d*xk(ij)**2/xep(ij)*rhof)*sxy(ij)
C	if ((i.eq.28.or.i.eq.30).and.j.eq.16) write(9,*)'^^1',i,j,
C     &tauxy(ij),sxy(ij),xnutc,rhof
C	    else
			xnutr=(c_d*xk(ij)**2/(xep(ij)+1.0e-16)*delx(i+1)+
     &		c_dr*xk(ipj)**2/(xep(ipj)+1.0e-16)*delx(i))
     &      	/(delx(i)+delx(i+1))
            	xnutru=(c_dt*xk(ijp)**2/(xep(ijp)+1.0e-16)*delx(i+1)
     &		+c_drt*xk(ipjp)**2/(xep(ipjp)+1.0e-16)*delx(i))
     &      	/(delx(i)+delx(i+1))
            	xnutc=(xnutr*dely(j+1)+xnutru*dely(j))/(dely(j)+dely(j+1))
              tauxyturb(ij)=2.*(xnutc*rhof)*sxy(ij)
              tauxy(ij)=2.*(xmu+xnutc*rhof)*sxy(ij)
C	    endif
	    goto 100

C.........Nonlinear eddy viscosity model with realizability
C.........Noted that tauxx and tauyy are defined at the center and
C.........tauxy is defined at the top-right vertice of cell
160       continue
C.........only calculate non-solid cell; for solid cell adjacent to fluid
C.........or porous, use boundary condition to calculate bed shear stress
          if (ac(ij)+ac(ipj)+ac(ijp)+ac(ipjp).lt.1.0d-06) go to 100

C.........porous media calculation
	    if (npc(ij).ne.1.and.npor.ne.0) then
            tauxx(ij)=2.*(xmu)*tauxx(ij)
            tauyy(ij)=2.*(xmu)*tauyy(ij)
            tauxy(ij)=2.*(xmu)*sxy(ij)
C...........specify bed shear stress
C...........boundary on left
		  if ((ac(imj).le.1.0e-6.and.ac(imjp).le.1.0e-6).
     &		or.(i.eq.2.and.kl.eq.1)) then
		    tauxy(imj)=2*xmu*sxy(imj)
		  endif
C...........boundary on bottom
		  if ((ac(ijm).le.1.0e-6.and.ac(ipjm).le.1.0e-6).
     &		or.(j.eq.2.and.kb.eq.1)) then
		    tauxy(ijm)=2*xmu*sxy(ijm)
		  endif
	      goto 100
		endif 

C.........non-porous part calculation
C.........eveluate du/dy at the cell center
		if (nf(ij).eq.6) then
			dudyc=0.0
			dvdxc=0.0
			dudxc=0.0
			dvdyc=0.0
		else
c      I added max(1,imj)  - Fengyan
			ucc=.5*(u(ij)*ar(ij)+u(max(1,imj))*ar(max(1,imj)))
			utc=.5*(u(ijp)*ar(ijp)+u(imjp)*ar(imjp))
			if (nf(ijp).eq.6) utc=ucc
			uuc=.5*(u(max(1,ijm))*ar(max(1,ijm))+u(max(1,imjm))
     &          *ar(max(1,imjm)))
			if (nf(max(1,ijm)).eq.6) uuc=ucc
			upc=(ucc*dely(j+1)+utc*dely(j))/(dely(j+1)+dely(j))
			umc=(ucc*dely(max(1,j-1))+uuc*dely(j))/(dely(max(1,j-1))
     &           +dely(j))
			dudyc=(upc-umc)/dely(j)

C.............evaluate dv/dx at the cell center
			vcc=0.5*(v(ij)*at(ij)+v(max(1,ijm))*at(max(1,ijm)))
			vrc=0.5*(v(max(1,ipj))*at(max(1,ipj))+
     &             v(max(1,ipjm))*at(max(1,ipjm)))
			if (nf(ipj).eq.6) vrc=vcc
			vlc=0.5*(v(max(1,imj))*at(max(1,imj))+v(max(1,imjm))
     &         *at(max(1,imjm)))
			if (nf(max(1,imj)).eq.6) vlc=vcc
			vpc=(vcc*delx(i+1)+vrc*delx(i))/(delx(i+1)+delx(i))
	vmc=(vcc*delx(max(1,i-1))+vlc*delx(i))/(delx(max(1,i-1))+delx(i))
			dvdxc=(vpc-vmc)/delx(i)

C.............define dudxc and dvdyc
			dudxc=tauxx(ij)
			dvdyc=-tauxx(ij)
		endif
C.........define coefficients
	    smax=abs(dudxc*xk(ij)/(xep(ij)+1.0e-16))
C	    smax=0.0
	    c_d=1./3./(3.7+smax)
C.........from Lemos (1992), pp. 27, Harlow & Nakayama, 1967 (Physics
C.........of Fluids), C_d is modified locally to reduce effective eddy
C.........viscosity for low turbulence.
	    c_d=c_d*xnut(ij)/(xnu*eddycoef)*
     &          (1.0-exp(-eddycoef*xnu/max(1.0d-20,xnut(ij))))

	    smax1=max(smax,abs(dudyc*xk(ij)/(xep(ij)+1.0e-16)),
     &		abs(dvdxc*xk(ij)/(xep(ij)+1.0e-16)))

C.........evaluate dudx at the top-right vertice
		if (nf(ij).eq.6.or.nf(ipj).eq.6.or.nf(ijp).eq.6.
     &		or.nf(ipjp).eq.6) then
			dudxrt=0.0
			dvdyrt=0.0
			dudyrt=0.0
			dvdxrt=0.0
		else
		  uct=(u(ij)*ar(ij)*dely(j+1)+u(ijp)*ar(ijp)*dely(j))
     &		/(dely(j)+dely(j+1))
		  ult=(u(max(1,imj))*ar(max(1,imj))*dely(j+1)+u(imjp)*
     &		ar(imjp)*dely(j))
     &		/(dely(j)+dely(j+1))
		  urt=(u(ipj)*ar(ipj)*dely(j+1)+u(ipjp)*
     &		ar(ipjp)*dely(j))
     &		/(dely(j)+dely(j+1))
	      if (i.eq.imax) then
			dudxrt=(uct-ult)/delx(i)
	      else
			if (i.eq.1) then
			  dudxrt=(urt-uct)/delx(i)
			else
	    	  dudxrt=((urt-uct)*delx(i)/delx(i+1)+(uct-ult)*delx(i+1)
     &			/delx(i))/(delx(i)+delx(i+1))
	        endif
		  endif

C...........evaluate dvdy at the top-right vertice
            dvdyrt=-dudxrt

C...........define dudyrt and dvdxrt
	      dudyrt=(u(ijp)*ar(ijp)-u(ij)*ar(ij))*2
     &		/(dely(j)+dely(j+1))
	      dvdxrt=(v(ipj)*at(ipj)-v(ij)*at(ij))*2
     &		/(delx(i)+delx(i+1))
		endif

C.........Shear term is the one which most likely causes instability
C.........problem. For two dimensional problem, the nonlinear shear stress
C.........terms can be lumped into one (c2-c3)*dudx*(dvdx-dudy); In other
C.........words, it is the vorticity rather than shear strain  contributing
C.........to the shear stress; thus, it is sensible to also include the
C.........vorticity into consideration when realizability is enforced.
C.........Negative production can be induced if sxy~0 but omega and dudx
C.........is large so that nonlinear effect dominates.
          smax1=max(smax1,abs((dudyc-dvdxc)*xk(ij)/(xep(ij)+1.0e-16)))

C.........It's difficult to decide the value of realize. When it is
C.........large, it reduces the contribution of nonlinear terms of
C.........strong vortices, but it will also affect coefficients when
C.........shear stress is low; should choose it between 2/3 to 20/3;
C.........the latter (larger) value should be used when instability
C.........is caused by too strong nonlinear contribution.
          c_1=2./3./(123.5+realize*smax1**2)*
     &          xnut(ij)/(xnu*eddycoef)*
     &          (1.0-exp(-eddycoef*xnu/max(1.0d-20,xnut(ij))))
          c_2=-2./3./(39.2+realize*smax1**2)*
     &          xnut(ij)/(xnu*eddycoef)*
     &          (1.0-exp(-eddycoef*xnu/max(1.0d-20,xnut(ij))))
          c_3=2./3./(246.9+realize*smax1**2)*
     &          xnut(ij)/(xnu*eddycoef)*
     &          (1.0-exp(-eddycoef*xnu/max(1.0d-20,xnut(ij))))

C.........calculate tauxx
	    xxl1turb=2.*(c_d*xk(ij)**2/(xep(ij)+1.0e-16)*rhof)*dudxc-2./3.
     &		  *rhof*xk(ij)
          xxl1=2.*(xmu+c_d*xk(ij)**2/(xep(ij)+1.0e-16)*rhof)*dudxc-2./3.
     &          *rhof*xk(ij)
	    xxnl1=rhof*xk(ij)**3/(xep(ij)+1.0e-16)**2*
     &          c_1*(2*(dudxc**2+dudyc*dvdxc)
     &          -2./3.*(dudxc**2+2*dudyc*dvdxc+dvdyc**2))
 	    xxnl2=rhof*xk(ij)**3/(xep(ij)+1.0e-16)**2*
     &          c_2*(dudxc**2+dudyc**2
     &          -1./3.*(dudxc**2+dudyc**2+dvdxc**2+dvdyc**2))
	    xxnl3=rhof*xk(ij)**3/(xep(ij)+1.0e-16)**2*
     &          c_3*(dudxc**2+dvdxc**2
     &          -1./3.*(dudxc**2+dudyc**2+dvdxc**2+dvdyc**2))
          tauxxturb(ij)=xxl1turb+xxnl1+xxnl2+xxnl3
	    tauxx(ij)=xxl1+xxnl1+xxnl2+xxnl3

C.........calculate tauyy for normal cell
          yyl1turb=2.*(c_d*xk(ij)**2/(xep(ij)+1.0e-16)*rhof)*dvdyc
     &          -2./3.*rhof*xk(ij)
	    yyl1=2.*(xmu+c_d*xk(ij)**2/(xep(ij)+1.0e-16)*rhof)*dvdyc
     &          -2./3.*rhof*xk(ij)
	    yynl1=rhof*xk(ij)**3/(xep(ij)+1.0e-16)**2*
     &          c_1*(2*(dvdyc**2+dudyc*dvdxc)
     &          -2./3.*(dudxc**2+2*dudyc*dvdxc+dvdyc**2))
	    yynl2=rhof*xk(ij)**3/(xep(ij)+1.0e-16)**2*
     &          c_2*(dvdxc**2+dvdyc**2
     &          -1./3.*(dudxc**2+dudyc**2+dvdxc**2+dvdyc**2))
	    yynl3=rhof*xk(ij)**3/(xep(ij)+1.0e-16)**2*
     &          c_3*(dudyc**2+dvdyc**2
     &          -1./3.*(dudxc**2+dudyc**2+dvdxc**2+dvdyc**2))
          tauyyturb(ij)=yyl1turb+yynl1+yynl2+yynl3
          tauyy(ij)=yyl1+yynl1+yynl2+yynl3

C.........Calculate tauxy on wall (wall or bed shear stress) by using 
C.........Three cases are treated: 1. computational boundary 
C.........2. solid adjacent to fluid 3. porous adjacent to fluid

C.........boundary on left
		if (((ac(max(1,imj)).le.1.0e-6.and.ac(imjp).le.1.0e-6).
     &		or.(npc(imj).ne.1.and.npc(imjp).ne.1.and.npor.ne.0).
     &		or.(i.eq.2.and.kl.eq.1)).and.at(ij).gt.em6) then
		  tauxy(imj)=rhof*sqrt(c_mu_bc)*(xk(ij)+xk(ijp))/2
     &		*(sign(1.0d0,v(ij)))+2*xmu*sxy(imj)
		  tauxyturb(imj)=rhof*sqrt(c_mu_bc)*(xk(ij)+xk(ijp))/2
     &		*(sign(1.0d0,v(ij)))
		endif
C.........boundary on bottom
		if (((ac(ijm).le.1.0e-6.and.ac(ipjm).le.1.0e-6).
     &		or.(npc(ijm).ne.1.and.npc(max(1,ipjm)).ne.1.and.npor.ne.0).
     &		or.(j.eq.2.and.kb.eq.1)).and.ar(ij).gt.em6) then
		  tauxy(ijm)=rhof*sqrt(c_mu_bc)*(xk(ij)+xk(ipj))/2
     &		*(sign(1.0d0,u(ij)))+2*xmu*sxy(ijm)
		  tauxyturb(ijm)=rhof*sqrt(c_mu_bc)*(xk(ij)+xk(ipj))/2
     &		*(sign(1.0d0,u(ij)))
		endif
C.........boundary on right
		if (((ac(ipj).le.1.0e-6.and.ac(ipjp).le.1.0e-6).
     &		or.(npc(ipj).ne.1.and.npc(ipjp).ne.1.and.npor.ne.0).
     &		or.(i.eq.im1.and.kr.eq.1)).and.at(ij).gt.em6) then
		  tauxy(ij)=rhof*sqrt(c_mu_bc)*(xk(ij)+xk(ijp))/2
     &		*(-sign(1.0d0,v(ij)))+2*xmu*sxy(ij)
		  tauxyturb(ij)=rhof*sqrt(c_mu_bc)*(xk(ij)+xk(ijp))/2
     &		*(-sign(1.0d0,v(ij)))
		  goto 100
		endif
C.........boundary on top
		if (((ac(ijp).le.1.0e-6.and.ac(ipjp).le.1.0e-6).
     &		or.(npc(ijp).ne.1.and.npc(ipjp).ne.1.and.npor.ne.0).
     &		or.(j.eq.jm1.and.kt.eq.1)).and.ar(ij).gt.em6) then
		  tauxy(ij)=rhof*sqrt(c_mu_bc)*(xk(ij)+xk(ipj))/2
     &		*(-sign(1.0d0,u(ij)))+2*xmu*sxy(ij)
		  tauxyturb(ij)=rhof*sqrt(c_mu_bc)*(xk(ij)+xk(ipj))/2
     &		*(-sign(1.0d0,u(ij)))
		  goto 100
		endif



C		if (nf(ij).eq.6.or.nf(ipj).eq.6.or.nf(ijp).eq.6.
C     &		or.nf(ipjp).eq.6) then
C			tauxy(ij)=0.0
C			tauxyturb(ij)=0.0
C          else
		    xkr=(xk(ij)*delx(i+1)+xk(ipj)*delx(i))
     &			/(delx(i)+delx(i+1))
			xkru=(xk(ijp)*delx(i+1)+xk(ipjp)*delx(i))
     &			/(delx(i)+delx(i+1))
			xkrt=(xkr*dely(j+1)+xkru*dely(j))/(dely(j)+dely(j+1))
			xepr=((xep(ij)+1.0e-16)*delx(i+1)+
     &			(xep(ipj)+1.0e-16)*delx(i))/(delx(i)+delx(i+1))
			xepru=((xep(ijp)+1.0e-16)*delx(i+1)
     &			+(xep(ipjp)+1.0e-16)*delx(i))/(delx(i)+delx(i+1))
			xeprt=(xepr*dely(j+1)+xepru*dely(j))/(dely(j)+dely(j+1))
			xnutr=(xnut(ij)*delx(i+1)+xnut(ipj)*delx(i))
     &			/(delx(i)+delx(i+1))
			xnutru=(xnut(ijp)*delx(i+1)+xnut(ipjp)*delx(i))
     &			/(delx(i)+delx(i+1))
			xnutrt=(xnutr*dely(j+1)+xnutru*dely(j))
     &			/(dely(j)+dely(j+1))

			smaxrt=abs(dudxrt*xkrt/xeprt)
			c_drt=1./3./(3.7+smaxrt)
			c_drt=c_drt*xnutrt/(xnu*eddycoef)*
     &			(1.0-exp(-eddycoef*xnu/max(1.0d-20,xnutrt)))
			xnutr=c_drt*xkrt**2/xeprt
		      
			xk3xep2rt=xkrt**3/xeprt**2

			smax1rt=max(smaxrt,abs(dudyrt*xkrt/xeprt),
     &		  abs(dvdxrt*xkrt/xeprt),abs((dudyrt-dvdxrt)*xkrt/xeprt))
			c_1rt=2./3./(123.5+realize*smax1rt**2)*
     &			xnutrt/(xnu*eddycoef)*
     &			(1.0-exp(-eddycoef*xnu/max(1.0d-20,xnutrt)))
			c_2rt=-2./3./(39.2+realize*smax1rt**2)*
     &			xnutrt/(xnu*eddycoef)*
     &			(1.0-exp(-eddycoef*xnu/max(1.0d-20,xnutrt)))
			c_3rt=2./3./(246.9+realize*smax1rt**2)*
     &			xnutrt/(xnu*eddycoef)*
     &			(1.0-exp(-eddycoef*xnu/max(1.0d-20,xnutrt)))
                
			xyl1turb=2.0*(xnutrt*rhof)*sxy(ij)
              xyl1=2.0*(xmu+xnutrt*rhof)*sxy(ij)
              xynl2=rhof*xk3xep2rt*c_2rt*(dudxrt*dvdxrt+dudyrt*dvdyrt)
              xynl3=rhof*xk3xep2rt*c_3rt*(dudxrt*dudyrt+dvdxrt*dvdyrt)
              tauxyturb(ij)=xyl1turb+xynl2+xynl3
              tauxy(ij)=xyl1+xynl2+xynl3
C          endif    
	  endif
c
  100 continue
c

C	if (abs((tauxy(ij)+tauxy(ijs-1))/(tauxy(ij)+1.0e-10)).gt.0.01.and.
C     &	abs(tauxy(ij)).gt.0.00001.and.i.ne.imax/2) then
C		write(9,*)'**tauxy',i,j,ncyc,'f',f(ij),f(ij+1),f(ij-1),
C     &	f(ij+imax),f(ij-imax),'ac',ac(ij),ac(ij+1),ac(ij-1),
C     &	ac(ij+imax),ac(ij-imax),ac(ij+1+imax),
C     &	'tauxy',tauxy(ij),tauxy(ijs)
C		stop
C	endif
C			tauxx(ij)=tauxx(ijs)
C	tauyy(ij)=tauyy(ijs)
C	if (i.eq.imax/2) then
C	tauxy(ij)=0
C	else
C	tauxy(ij)=-tauxy(ijs-1)
C	endif

C	write(9,'(1000f8.3)')(tauxy(

C	do i=1,imax
C		do j=1,jmax
C		ij=(j-1)*imax+i
C
C			tauxx(ij)=0.0
C	tauyy(ij)=0.0
C	tauxy(ij)=0.0
C	end do
C	end do
      return
      end
  

      subroutine vofadv
c
c ======================================================================
c
c   Purpose -
c     compute volume of fluid advection fluxes from the newly
c     determined velocity field, updating the volume of fluid
c     function f(i,j)
c
c   VOFADV is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE
c
c
c   VOFADV calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c              BC    ripple   SETARRY    ripple     SETNF    ripple
c          VOFCOR    ripple    VOFERR    ripple   VOFPACK    ripple
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
      dimension umom(nxy),vmom(nxy)

      data tiny /1.0d-25/, zero /0.0d0/, tenth /0.10d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

c
c.... reset the pertinent arrays
c     ---------------------------------
      call setarry (umom(1),0.0d0,nxy)
      call setarry (vmom(1),0.0d0,nxy)
c     ---------------------------------
c
c.... Compute cell momenta
c
      do 200 j=1,jm1
        do 200 i=1,im1
          ij=(j-1)*imax+i
          ip=i+1
          jp=j+1
          ipj=ij+1
          ijp=ij+imax
          rhoij=f(ij)*rhof
          rhoipj=f(ipj)*rhof
          rhoijp=f(ijp)*rhof
          rhorf=(delx(ip)*rhoij+delx(i)*rhoipj)/(delx(i)+delx(ip))
          rhotf=(dely(jp)*rhoij+dely(j)*rhoijp)/(dely(j)+dely(jp))
          umom(ij)=rhorf*u(ij)
          vmom(ij)=rhotf*v(ij)
  200 continue
c
      if (ncyc.lt.1) go to 100
c
      flgc=0.0
C.....generally speaking, how to sweep does not affect the results because most of
C.....computation uses information of fn(ij). The only exception is epsd which is based
C.....on current f(ij) to prevent the negative f(ij) calculation. Because of this,
C.....the order of sweep can causes difference when a cell is to be empty. The following
C.....is to maintain horizontal symmetry be ensure left and right have the same sweep order
      do 50 j=1,jm1
        do 50 icyc=1,im1
		if (mod(icyc,2).eq.1) then
			i=icyc/2+1
		else
			i=imax-icyc/2
		endif
          ij=(j-1)*imax+i
          vx=u(ij)*delt
          vy=v(ij)*delt
          abvx=abs(vx)*ar(ij)
          abvy=abs(vy)*at(ij)
          if ((abvx.le.fcvlim*delx(i).or.j.eq.1).and.(abvy.
     &         le.fcvlim*dely(j).or.i.eq.1)) go to 5

          if (ac(ij).gt.em6.and.autot.ne.0.0) then
            flgc=1.
          endif

c
C    5     if (mod(ncyc,2).eq.0) go to 30
    5     if (mod(icyc,2).eq.1) go to 30
   10     continue
C          if ((ar(ij).lt.em6).and.(mod(ncyc,2).eq.0)) 
C     &		go to 55
          if ((ar(ij).lt.em6).and.(mod(icyc,2).eq.1)) 
     &		go to 55
          if (ar(ij).lt.em6) go to 30
          ia=i+1
          id=i
          idm=max0(i-1,1)
          ijdm=(j-1)*imax+idm
          ardm=ar(ijdm)
          rb=ar(ij)*r(i)+tiny
C          ra=cvmgt(ac(ij+1),porousp(ij+1),npc(ij+1).eq.1)
C     &		*ri(i+1)+tiny
C          rd=cvmgt(ac(ij),porousp(ij),npc(ij).eq.1)
C     &		*ri(i)+tiny
C.....for small ac, it requires extremely small dt to avoid interior air,
C.....even when div(u*theta)=0, due to the restriction imposed above
C.....In the modification below, such restriction will not be applied
C.....in the interior region, where cell is supposed to be full all time
          ra=cvmgt(cvmgt(ac(ij+1),1.0d0,nf(ij+1).ne.0.
     &		or.fn(ij+1).le.0.9),
     &		porousp(ij+1),npc(ij+1).eq.1)*ri(i+1)+tiny
          rd=cvmgt(cvmgt(ac(ij),1.0d0,nf(ij).ne.0.or.fn(ij).le.0.9),
     &		porousp(ij),npc(ij).eq.1)*ri(i)+tiny
          incf=1
          incu=0
          if (vx.ge.0.0) go to 20
          ia=i
          id=i+1
          idm=min0(i+2,imax)
          ijdm=(j-1)*imax+idm
          ardm=ar(ijdm-1)
C          ra=cvmgt(ac(ij),porousp(ij),
C     &		npc(ij).eq.1)*ri(i)+tiny
C          rd=cvmgt(ac(ij+1),
C     &		porousp(ij+1),npc(ij+1).eq.1)*ri(i+1)+tiny
          ra=cvmgt(cvmgt(ac(ij),1.0d0,nf(ij).ne.0.or.fn(ij).le.0.9),
     &		porousp(ij),npc(ij).eq.1)*ri(i)+tiny
          rd=cvmgt(cvmgt(ac(ij+1),1.0d0,nf(ij+1).ne.0.
     &		or.fn(ij+1).le.0.9),
     &		porousp(ij+1),npc(ij+1).eq.1)*ri(i+1)+tiny
          incf=-1
          incu=-1
c
   20     continue
          ija=(j-1)*imax+ia
          ijd=(j-1)*imax+id
          ijdu=(j-1)*imax+id+incu
          iad=ia
          if (nf(ijd).eq.3.or.nf(ijd).eq.4) iad=id
C.........since the single void with f>0 in the interior region might
C.........be identified as the full cell, we need special scheme to
C.........recognize this case and set the scheme to donar mathod without
C.........corrector so that the void can move instead of be static
          if (fn(ijd).le.0.999.and.nf(ijd+1).eq.0.and.nf(ijd-1).eq.0.
     &          and.nf(ijd+imax).eq.0.and.nf(ijd-imax).eq.0) iad=id
          if (fn(ija).lt.emf.or.fn(ijdm).lt.emf) iad=ia
          ijad=(j-1)*imax+iad
          fdm=max(fn(ijdm),fn(ijd),1.0d-4)
          if(ardm.lt.em6) fdm=1.0d0
C.........It it found that by setting fdm=1.0d0 for empty donor cell, it can prevent
C.........free surface front from spreading out by slowing down the actual motion;
C.........useful to maintain stability during breaking
		if (nf(ija).eq.6.and.nf(ijdm).ne.6) fdm=1.0d0
C		if (nf(ija).eq.6) fdm=1.0d0
          if (iad.eq.ia) then
            fx1=fn(ijad)*abs(vx)+max((fdm-fn(ijad))*abs(vx)-
     &         (fdm-fn(ijd))*delx(id),zero)
          else
            fx1=fn(ijad)*abs(vx)
          endif
          fx=min(fx1,fn(ijd)*delx(id)*rd/rb)
          delfd=fx*rdx(id)*(rb/rd)
          epsd=f(ijd)-delfd
          delfd=cvmgt(delfd+epsd,delfd,epsd.lt.zero)
          delfa=rdx(ia)*rd*delfd/(rdx(id)*ra)
          f(ijd)=f(ijd)-delfd
          f(ija)=f(ija)+delfa
          delma=delfa*rhof*cvol(ijd)
          delmd=delfd*rhof*cvol(ijd)
          umom(ija)=umom(ija)+delma*u(ijdu)
          umom(ijd)=umom(ijd)-delmd*u(ijdu)
C          if (mod(ncyc,2).eq.0) go to 55
          if (mod(icyc,2).eq.1) go to 55
c
   30     if (at(ij).lt.em6) go to 45
          ja=j+1
          jd=j
          jdm=max0(j-1,1)
          ijdm=(jdm-1)*imax+i
          atdm=at(ijdm)
          rb=at(ij)+tiny
C          ra=cvmgt(ac(ij+imax),
C     &		porousp(ij+imax),npc(ij+imax).eq.1)+tiny
C          rd=cvmgt(ac(ij),
C     &		porousp(ij),npc(ij).eq.1)+tiny
          ra=cvmgt(cvmgt(ac(ij+imax),1.0d0,nf(ij+imax).ne.0.
     &		or.fn(ij+imax).le.0.9),
     &		porousp(ij+imax),npc(ij+imax).eq.1)+tiny
          rd=cvmgt(cvmgt(ac(ij),1.0d0,nf(ij).ne.0.or.fn(ij).le.0.9),
     &		porousp(ij),npc(ij).eq.1)+tiny
          incf=1
          incv=0
          if (vy.ge.0.0) go to 40
          ja=j
          jd=j+1
          jdm=min0(j+2,jmax)
          ijdm=(jdm-1)*imax+i
          atdm=at(ijdm-imax)
C          ra=cvmgt(ac(ij),
C     &		porousp(ij),npc(ij).eq.1)+tiny
C          rd=cvmgt(ac(ij+imax),
C     &		porousp(ij+imax),npc(ij+imax).eq.1)+tiny
          ra=cvmgt(cvmgt(ac(ij),1.0d0,nf(ij).ne.0.or.fn(ij).le.0.9),
     &		porousp(ij),npc(ij).eq.1)+tiny
          rd=cvmgt(cvmgt(ac(ij+imax),1.0d0,nf(ij+imax).ne.0.
     &		or.fn(ij+imax).le.0.9),
     &		porousp(ij+imax),npc(ij+imax).eq.1)+tiny
          incf=-1
          incv=-1
c
   40     continue
          ija=(ja-1)*imax+i
          ijd=(jd-1)*imax+i
          ijdv=(jd-1+incv)*imax+i
          jad=ja
          if (nf(ijd).eq.1.or.nf(ijd).eq.2) jad=jd
C.........since the single void with f>0 in the interior region might
C.........be identified as the full cell, we need special scheme to
C.........recognize this case and set the scheme to donar mathod without
C.........corrector so that the void can move instead of be static
          if (fn(ijd).le.0.999.and.nf(ijd+1).eq.0.and.nf(ijd-1).eq.0.
     &          and.nf(ijd+imax).eq.0.and.nf(ijd-imax).eq.0) jad=jd
          if (fn(ija).lt.emf.or.fn(ijdm).lt.emf) jad=ja
          ijad=(jad-1)*imax+i
          fdm=max(fn(ijdm),fn(ijd),1.0d-4)
          if(atdm.lt.em6) fdm=1.0d0
		if (nf(ija).eq.6.and.nf(ijdm).ne.6) fdm=1.0d0
C		if (nf(ija).eq.6) fdm=1.0d0
          if (jad.eq.ja) then
            fy1=fn(ijad)*abs(vy)+max((fdm-fn(ijad))*abs(vy)-
     &         (fdm-fn(ijd))*dely(jd),zero)
          else
            fy1=fn(ijad)*abs(vy)
          endif
          fy=min(fy1,fn(ijd)*dely(jd)*rd/rb)
          delfd=fy*rdy(jd)*(rb/rd)
          epsd=f(ijd)-delfd
          delfd=cvmgt(delfd+epsd,delfd,epsd.lt.zero)
          delfa=rdy(ja)*rd*delfd/(rdy(jd)*ra)
          f(ijd)=f(ijd)-delfd
          f(ija)=f(ija)+delfa
          delma=delfa*rhof*cvol(ijd)
          delmd=delfd*rhof*cvol(ijd)
          vmom(ija)=vmom(ija)+delma*v(ijdv)
          vmom(ijd)=vmom(ijd)-delmd*v(ijdv)
C   45     if (mod(ncyc,2).eq.0) go to 10
   45     if (mod(icyc,2).eq.1) go to 10

   55 continue

	if (ninflow.eq.100) then
          if (i.ge.isources.and.i.le.isourcee.and.
     &    j.ge.jsources.and.j.le.jsourcee) f((j-1)*imax+i)=1.0
	endif

   50 continue
c
c.... defoam it by packing vertically
c
c       -------------------------------------------
C        call vofpack(2,im1,2,jbar,imax,at,nf,f,ac,dely,nxy,npack)
c       -------------------------------------------
c
c
c.... divergence correction
c
      if (idiv.ne.0)
c       -------------------------------------------
     &  call vofcor(2,im1,2,jm1,imax,nf,fn,f,
     &              ac,at,ar,rdx,rdy,r,ri,u,v,
     &              delt,kt,kb,kl,kr,uxmb,vymb,nmovbd)
c       -------------------------------------------
c
c
  100	continue
c
c.... make sure the new vof function is
c     above the floor and below the ceiling
c
c     --------------------------------------------------
      call voferr
c     --------------------------------------------------
c
c
c.... Update boundary conditions
c
      nnn=1
      ibcflg=1
C.....added to account for the impact
      ibcfinal=1
c     -------------
      call setnf
      call bc
c     -------------
      ibcfinal=0
      ibcflg=0
      nnn=0

C	do j=2,jm1
C	do i=2,im1
C	ij=(j-1)*imax+i
C	ijm=ij-imax
C	imj=ij-1
C	if (f(ij).gt.0.0d0.and.ac(ij).gt.0.0.and.i.gt.50) then
C	if (nfold(ij).eq.6.and.nf(ij).eq.0) goto 111
C	if ((u(ij)*ar(ij)-u(imj)*ar(imj))/delx(i)+(v(ij)*at(ij)
C     &	-v(ijm)*at(ijm))/dely(j).gt.1.0e-5) then
C	write(9,*)'!!!$2',i,j,f(ij),f(ij+1),f(imj),f(ij+imax),f(ijm),
C     &	f(ij+1+imax),f(imj+imax),f(ijm+1),f(ijm-1),
C     &	nf(ij),nf(ij+1),nf(imj),nf(ij+imax),nf(ijm),
C     &	nf(ij+1+imax),nf(imj+imax),nf(ijm+1),nf(ijm-1),
C     &	nfold(ij),nfold(ij+1),nfold(imj),nfold(ij+imax),nfold(ijm),
C     &	nfold(ij+1+imax),nfold(imj+imax),nfold(ijm+1),nfold(ijm-1),
C     &	ac(ij),ar(ij),ar(imj),at(ij),at(ijm),
C     &	u(ij),u(imj),v(ij),v(ijm),
C     &	(u(ij)*ar(ij)-u(imj)*ar(imj))/delx(i)+(v(ij)*at(ij)
C     &	-v(ijm)*at(ijm))/dely(j)
C	endif
C	endif
C111	continue
CC	end do
c
 9999 return
c
  150 format (1x,12hVOFADV ERROR/1x,5hncyc=,i7,1x,2ht=,1pe14.6,1x,
     1         5hdelt=,e12.4,1x,4hi,j=,2i4,1x,7hfcvlim=,e11.3/3x,
     2         5habvx=,e12.4,1x,5hdelx=,e12.4,1x,5habvy=,e12.4,1x,
     3         5hdely=,e12.4,2e12.4)
      end

 
      subroutine vofcor(i1,i2,j1,j2,nqj,nf,fn,f,ac,at,ar,
     &           rdx,rdy,r,ri,u,v,dt,kt,kb,kl,kr,uxmb,vymb,nmovbd)
c
c ======================================================================
c
c   Purpose -
c     correct the vof function for velocity divergence errors
c
c   VOFCOR is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          VOFADV
c
c
c   VOFCOR calls the following subroutines and functions -
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
      dimension at(1),nf(1),f(1),fn(1),ac(1),ar(1),
     &   rdx(1),rdy(1),r(1),ri(1),u(1),v(1),uxmb(1),vymb(1),nmovbd(1)
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
	imax=nqj
      do 80 j=j1,j2
        do 80 i=i1,i2
          ij=(j-1)*nqj+i
          imj=ij-1
          ijm=ij-nqj
          if (ac(ij).le.1.0e-6) go to 80
C.........modified to eliminate the empty cell adjacent to surface
C.........where divergence.ne.0.0
	    if (fn(ij).le.0.0) then
C			f(ij)=0.0d0
			go to 80
		endif
          divergence=(rdx(i)*(ar(ij)*r(i)*(u(ij)-uxmb(nmovbd(ij)))-
     1             ar(imj)*r(i-1)*(u(imj)-uxmb(nmovbd(ij))))+rdy(j)*
     2             (at(ij)*ri(i)*(v(ij)-vymb(nmovbd(ij)))-at(ijm)*ri(i)*
     3             (v(ijm)-vymb(nmovbd(ij)))))/(ri(i))
     4		   /cvmgt(ac(ij),1.0d0,nf(ij).ne.0.or.fn(ij).le.0.9)
          f(ij)=f(ij)+dt*fn(ij)*divergence
   80 continue
c
      return
      end

 
      subroutine voferr
c
c ======================================================================
c
c   Purpose -
c     Correct the VOF function if it lies below the floor
c     (emf) or above the ceiling (emf1)
c
c   VOFERR is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          VOFADV
c
c
c   VOFERR calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c            none
c
c ======================================================================
c
c##############################################################
      implicit real*8 (a-h,o-z)
      
      include  "comdk2.h"
c##############################################################
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
C	Add correction term for moving body and set the upper and lower limit
	do 80 i=2,im1
	do 80 j=2,jm1
		ij=(j-1)*imax+i
		ipj=ij+1
		imj=ij-1
		ijp=ij+imax
		ijm=ij-imax
		vchg=0.0
		if (ac(ij).le.em6) goto 80
C.........modified to eliminate the empty cell adjacent to surface
C.........where divergence.ne.0.0
	    if (fn(ij).le.0.0) then
C			f(ij)=0.0d0
C			go to 80
			goto 85
		endif

	    f(ij)=f(ij)
     &	+fn(ij)*delt*(uxmb(nmovbd(ij))*(ar(ij)-ar(ij-1))/delx(i)
     &	+vymb(nmovbd(ij))*(at(ij)-at(ij-imax))/dely(j))
     &	/cvmgt(ac(ij),1.0d0,nf(ij).ne.0.or.fn(ij).le.0.9)

   85		continue
	    if (f(ij).gt.0.999999) then
		  vchg=1.0d0-f(ij)
C...........The acceptor method works perfectly for 1D case when a surface cell is
C...........to be filled; however, for 2D case, when a surface cell is to be filled,
C...........there is a large chance that it will be over-filled from another dimension
C...........where donor method is used; in this case f>1; The example include f(ij)=0.7
C...........f(ipj)=0.5, f(imj)=0.9, f(ijp)=0, f(ijm)=1.0 and Crx=Cry=0.3; In the vertical
C...........direction, the acceptor method will be used and it will fill up f(ij)=1.0;
C...........In the horizontal direction, the donor method will be used and it will further
C...........add in fluid in this cell; In the following treatment,
C...........this overfilling is re-distributed into the smallest f cell (near empty),
C...........as it is supposed to be if the VOFADV scheme is perfect
C...........A better VOF scheme could possibly avoid this problem
		  if (vchg.lt.-1.0d-6) then
			if (f(ipj).lt.f(imj)) then
				ija=ipj
			else
				ija=imj
			endif
			if (f(ijp).lt.f(ijm)) then
				ijb=ijp
			else
				ijb=ijm
			endif
			if (f(ija).lt.f(ijb)) then
				if (ac(ija).gt.em6) 
     &			f(ija)=min(1.0d0,f(ija)-vchg*ac(ij)/ac(ija))
			else
				if (ac(ijb).gt.em6)
     &			f(ijb)=min(1.0d0,f(ijb)-vchg*ac(ij)/ac(ijb))
			endif
		  endif
	      f(ij)=1.0d0
	    endif

	    if (f(ij).lt.0.0001) then
		  vchg=-f(ij)
		  if (vchg.gt.1.0d-4) then
			if (f(ipj).gt.f(imj)) then
				ija=ipj
			else
				ija=imj
			endif
			if (f(ijp).gt.f(ijm)) then
				ijb=ijp
			else
				ijb=ijm
			endif
			if (f(ija).gt.f(ijb)) then
				if (ac(ija).gt.em6) 
     &			f(ija)=max(0.0d0,f(ija)-vchg*ac(ij)/ac(ija))
			else
				if (ac(ijb).gt.em6)
     &			f(ijb)=max(0.0d0,f(ijb)-vchg*ac(ij)/ac(ijb))
			endif
		  endif
	      f(ij)=0.0d0
	    endif

C.........eliminate small isolated & splashed fluid and fluid above abstacle
	    if (nf(ij).ne.0.and.f(ij).gt.0.0.and.
     &	f(ij)*ac(ij)+f(ipj)*ac(ipj)+f(imj)*ac(imj)
     &	+f(ijp)*ac(ijp)+f(ijm)*ac(ijm).lt.0.1) then
		  vchg=-f(ij)
	      f(ij)=0.0d0
	    endif

		vchgt=vchgt+vchg*delx(i)*dely(j)*
     &	  cvmgt(porousp(ij),ac(ij),npc(ij).ne.1)
   80	continue

	do i=2,im1
	do j=2,jm1
	  ij=(j-1)*imax+i
	  ipj=ij+1
	  imj=ij-1
	  ijp=ij+imax
	  ijm=ij-imax
	  if (ac(ij).lt.em6) then
C.........define the proper f(ij) in obstacle adjacent to fluid cell
	    if (ac(ipj)+ac(imj)+ac(ijp)+ac(ijm).lt.em6) then
			f(ij)=1.0d0
		else
			if (ac(ipj).gt.ac(imj)) then
				ija=ipj
			else
				ija=imj
			endif
			if (ac(ijp).gt.ac(ijm)) then
				ijb=ijp
			else
				ijb=ijm
			endif
			if (ac(ija).gt.ac(ijb)) then
				f(ij)=f(ija)
			else
				f(ij)=f(ijb)
			endif
		endif
	  endif
	end do
	end do

      return
      end


      subroutine vofpack(i1,i2,j1,j2,nqj,at,nf,f,ac,dy,nxy,npack)
c
c ======================================================================
c
c   Purpose -
c     pack the vof function vertically downward
c
c   VOFPACK is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          VOFADV
c
c
c   VOFPACK calls the following subroutines and functions -
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
      dimension at(1),nf(1),f(1),ac(nxy),dy(nxy)
      data emf /1.0d-06/, emf1 /0.999999d0/, one /1.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c

C.....Do a real packing down job; no porous media only
      do 70 i=i1,i2
        tvof=0.0d0
C.......integrate vol vertically
        do 50 j=j1,j2
          ij=(j-1)*nqj+i
          tvof=tvof+ac(ij)*f(ij)*dy(j) 
   50   continue
C.......redistribute the vof nicely if npack=1 (forced packing not suitable for breaking wave or thin layer front)
	  if (npack.eq.1) then
          do 60 j=j1,j2
          ij=(j-1)*nqj+i
	    if (ac(ij).lt.1.0e-6) then
			f(ij)=1.0
		else
			f(ij)=max(min(tvof/dy(j)/ac(ij),1.0d0),0.0d0)
		endif
          tvof=tvof-f(ij)*ac(ij)*dy(j)
   60     continue
	  endif
   70 continue  

      return
      end


      subroutine vtilde
c
c ======================================================================
c
c   Purpose -
c     compute the incremental velocity change due to
c     gravitational acceleration
c
c   VTILDE is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE
c
c   VTILDE calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c              BC    ripple   CONVECT    ripple  CONVECTC    ripple
c          STRESS    ripple
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
      data zero /0.0d0/
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      do 1 j=1,jmax
        do 1 i=1,imax
          ij=(j-1)*imax+i
          u(ij)=un(ij)
          v(ij)=vn(ij)
    1 continue
c
c.... get the stress tensor
c
c     ------------
      call stress
c     ------------
c
c     --------------
      call convect
c     --------------
c
      visx=0.0
      visy=0.0
      tiny=1.0d-25
      do 20 j=2,jm1
        m=1
        n=1
        do 20 i=2,im1
c
          ij=(j-1)*imax+i
          ipj=ij+1
          imj=ij-1
          ijp=ij+imax
          ijm=ij-imax
          ipjp=ipj+imax
          ipjm=ipj-imax
          imjp=imj+imax
          imjm=imj-imax
c
          rhoij=f(ij)*rhof
          rhoipj=f(ipj)*rhof
          rhoijp=f(ijp)*rhof
          rhorc=(delx(i+1)*rhoij+delx(i)*rhoipj)/(delx(i)+delx(i+1))
          rhotc=(dely(j+1)*rhoij+dely(j)*rhoijp)/(dely(j)+dely(j+1))
          rhox=rhorc+tiny
          rhoy=rhotc+tiny
c
c....     check for surrounding fluid
          rhobar=rhox
          if (ar(ij).lt.em6) go to 10
c
          if (i.eq.im1.and.kr.eq.5) goto 35

c....     compute the x-direction viscous stress
c
          tauxyy=(tauxy(ij)-tauxy(ijm))/dely(j)
          tauxxx=(ri(i+1)*tauxx(ipj)-ri(i)*tauxx(ij))
     &                /(r(i)*0.5*(delx(i+1)+delx(i)))
C		if (nf(ij).eq.6.or.nf(ipj).eq.6) tauxxx=0.0
          visx=(tauxxx+tauxyy)/rhof

C.........enforce vanishing stress gradient on free surface 
	    if (fn(ij).ge.6.or.fn(ipj).ge.6) visx=0.0
C.........end
c
    5     continue
c....     explicitly update the x-velocity with the
c         x-direction convective flux, viscous stress,
c         body force, and pressure force
          rhox=rhobar*(delx(i+1)+delx(i))
          u(ij)=u(ij)+delt*visx
   10     u(ij)=cvmgt(zero,u(ij),((rhorc.lt.frsurf).and.
     &          ((nf(ij).gt.5).and.(nf(ipj).gt.5))).or.
     &		(ar(ij).lt.em6))
c
c....     reset y-velocity and check for surrounding fluid
          rhobar=rhoy
          if (at(ij).lt.em6) go to 25
c
35        if (j.eq.jm1.and.kt.eq.5) goto 20

c....     compute the y-direction viscous stress
c
CC          tauxyx=(r(i)*tauxy(ipjp)-r(i-1)*tauxy(ijp))/(ri(i)*delx(i))
          tauxyx=(r(i)*tauxy(ij)-r(i-1)*tauxy(imj))/(ri(i)*delx(i))
C.........end of modification
          tauyyy=2.*(tauyy(ijp)-tauyy(ij))/(dely(j+1)+dely(j))
C		if (nf(ij).eq.6.or.nf(ijp).eq.6) tauyyy=0.0
          visy=(tauyyy+tauxyx)/rhof
C.........enforce vanishing stress gradient on free surface
          if (fn(ij).ge.6.or.fn(ijp).ge.6) visy=0.0
C.........end
c
   15     continue
c....     explicitly update the y-velocity with the
c         y-direction convective flux, viscous stress,
c         body force, and pressure force
          rhoy=rhobar*(dely(j+1)+dely(j))
          v(ij)=v(ij)+delt*visy

17        continue

   25     v(ij)=cvmgt(zero,v(ij),((rhotc.lt.frsurf).and.
     &          ((nf(ij).gt.5).and.(nf(ijp).gt.5))).or.
     &		(at(ij).lt.em6))
c
   20 continue
c
	ibcflg=1
C	ibcflg0=1
c     ---------
      call bc
c     ---------
C	ibcflg0=0
	ibcflg=0
c
      return
      end


       subroutine ymax9p(e,d,c,b,a,bu,cu,du,eu,nh,nv,x,y)
c
c***purpose
c
c      form the matrix product y=y-a*x,where the matrix a is
c      a 9-diagonal matrix associated with a 9-point difference
c      operator.
c
c***description
c      
c       ymax9p requires that the coefficients of the ith 
c       equation be stored in the ith element of the diagonal
c       arrays. nine diagonal arrays contain coefficients
c       as one sweeps accross the matrix a from left to right.
c   
c     on entry
c
c        e,d,c,b,a,bu,cu,du,eu,   are 1-d arrays of length
c            at least nhxnv. the arrays contain, from left to right,
c            the diagonals of a as one sweeps a left to right.
c       
c              a is viewed as a block tridiagonal system of equations.
c      
c        nv    is the block order if the matrix is viewed
c              as a block tridiagonal matrix.
c    
c        x     is an input vector of dimenison at least nhxnv.
c            
c        y     is an input vector of dimension at least nhxnv.
c
c     on return
c              
c        y     is the output vector of dimension at least
c              nhxnv containing y-a*x.
c
c
c***routines called   none
c
      implicit real*8 (a-h,o-z)
      real*8 e(*),d(*),c(*),b(*),a(*),bu(*),cu(*),du(*),eu(*),
     &       x(*),y(*)
c
c*** first executable statement   ymax9p
      n = nh*nv
      nhm1 = nh-1
      nhp1 = nh+1
      i = 1
      y(i) = y(i)-(a(i)*x(i)+bu(i)*x(i+1)
     1      +du(i)*x(i+nh)+eu(i)*x(i+nhp1))
      do 10 i=2,nh
      y(i) = y(i)-(b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1)
     1        +(cu(i)*x(i+nhm1)+(du(i)*x(i+nh)+(eu(i)*x(i+nhp1))
     2        )))))
   10 continue
      i = nhp1
      y(i) = y(i)-(d(i)*x(i-nh)+c(i)*x(i-nhm1)
     1        +b(i)*x(i-1)+a(i)*x(i)+bu(i)*x(i+1)
     2        +cu(i)*x(i+nhm1)+du(i)*x(i+nh)+eu(i)*x(i+nhp1))
      do 20 i=nh+2,n-nhp1
      y(i) = y(i)-(e(i)*x(i-nhp1)+(d(i)*x(i-nh)+(c(i)*x(i-nhm1)
     1        +(b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1)
     2        +(cu(i)*x(i+nhm1)+(du(i)*x(i+nh)+(eu(i)*x(i+nhp1))
     3        ))))))))
   20 continue
      i=n-nh
      y(i) = y(i)-(e(i)*x(i-nhp1)+d(i)*x(i-nh)+c(i)*x(i-nhm1)
     1        +b(i)*x(i-1)+a(i)*x(i)+bu(i)*x(i+1)
     2        +cu(i)*x(i+nh-1)+du(i)*x(i+nh))
      do 30 i=n-nhm1,n-1
      y(i) = y(i)-(e(i)*x(i-nhp1)+(d(i)*x(i-nh)+(c(i)*x(i-nhm1)
     1        +(b(i)*x(i-1)+(a(i)*x(i)+(bu(i)*x(i+1))
     2        )))))
   30 continue
      i = n
      y(i) = y(i)-(e(i)*x(i-nhp1)+d(i)*x(i-nh)
     1      +b(i)*x(i-1)+a(i)*x(i))
      return
      end


C=================Youngs' model for tracking free surface more accurately=================
C=================The current version works well for non-obstacle case====================
C=================For irregular obstacle (e.g., slope), further testing is needed=========

      subroutine Youngs

      implicit real*8 (a-h,o-z)
      integer dir_flag

      include  "comdk2.h" 

c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C.......not ready for wave runup on dry beach

	do ntmp=1,imax*jmax
	  ftmp(ntmp)=f(ntmp)
	end do

c     ...Cycle through the advection loop

      	flgc=0.0
	do 1000 j=1,jm1
	 do 1000 i=1,im1

          ij=(j-1)*imax+i
          vx=u(ij)*delt
          vy=v(ij)*delt
          abvx=abs(vx)*ar(ij)
          abvy=abs(vy)*at(ij)
          if (abvx.le.fcvlim*delx(i).and.abvy.
     &         le.fcvlim*dely(j)) go to 1000
          if (ac(ij).gt.em6) then
            flgc=1.
C            write (9,150) ncyc,t,delt,i,j,fcvlim,abvx,delx(i),
C     &                        abvy,dely(j)
          endif

  150 format (1x,12hVOFADV ERROR/1x,5hncyc=,i7,1x,2ht=,1pe14.6,1x,
     1         5hdelt=,e12.4,1x,4hi,j=,2i4,1x,7hfcvlim=,e11.3/3x,
     2         5habvx=,e12.4,1x,5hdelx=,e12.4,1x,5habvy=,e12.4,1x,
     3         5hdely=,e12.4,2e12.4)

1000	continue

   10 id_order = mod(ncyc,2)

      do idir = 1,isweeps

        if (four_sweep) then
           dir_flag = id_order + mod(idir,2) +
     &                mod(idir,3) + mod(idir,4)
        else
           dir_flag = id_order + mod(idir,2)
        end if

        do 2000 j=1,jm1
         do 2000 i=1,im1
          ij=(j-1)*imax+i
	  if (i.eq.1.and.j.eq.1) then
		utt(ij)=u(ij)
		vrr(ij)=v(ij)
		goto 1200
	  endif
	  if (i.gt.1) then
          	utt(ij)=(u(ij)+u(ij-1)+u(ij+imax)+u(ij-1+imax))/4.0
          	if (nf(ij+imax).ge.6) utt(ij)=(u(ij)+u(ij-1))/2.0
	  else
		utt(ij)=(u(ij)+u(ij+imax))/2.0
                if (nf(ij+imax).ge.6) utt(ij)=u(ij)
	  endif
          if (j.gt.1) then
          	vrr(ij)=(v(ij)+v(ij+1)+v(ij-imax)+v(ij+1-imax))/4.0
          	if (nf(ij+1).ge.6) vrr(ij)=(v(ij)+v(ij-imax))/2.0
	  else
		vrr(ij)=(v(ij)+v(ij+1))/2.0
                if (nf(ij+1).ge.6) vrr(ij)=v(ij)
	  endif
1200      flux_pt(1,ij) = x(i) - frac*u(ij)*ar(ij)*delt
          flux_pt(2,ij) = y(j) - frac*v(ij)*at(ij)*delt

	  if (i.gt.1) then
		dx=x(i)-x(i-1)
          	ravg = cyl*0.5*(x(i)+x(i-1))+1.-cyl
	  else
		dx=x(2)-x(1)
		ravg = cyl*0.5*(x(1)+x(1))+1.-cyl	
	  endif
	  if (j.gt.1) then
		dy=y(j)-y(j-1)
	  else
		dy=y(2)-y(1)
	  endif
          cell_vol(ij) = ravg*dx*dy
          advect_vol(ij) = cell_vol(ij)

          dvol_tot(1,ij) = u(ij)*ar(ij)*frac*delt*dy
     &                      *(cyl*x(i) + (1.0 - cyl))
          if (i.gt.1) then
            if (f(ij).eq.0.0.and.f(ij-1).gt.0.0)
     &          dvol_tot(1,ij) = u(ij-1)*ar(ij-1)*frac*delt*
     &                  (y(j)-y(j-1))*(cyl*x(i) + (1.0 - cyl))
          endif
          dvol_tot(2,ij) = v(ij)*at(ij)*frac*delt*dx
     &                      *(cyl*xi(i) + (1.0 - cyl))
          if (j.gt.1) then
            if (f(ij).eq.0.0.and.f(ij-imax).gt.0.0)   
     &          dvol_tot(2,ij) = v(ij-imax)*at(ij-imax)*frac*delt*
     &                  (x(i)-x(i-1))*(cyl*xi(i) + (1.0 - cyl))
          endif
2000	continue
c     ...Boundary conditions

      j = 1          ! Bottom
      do i = 1,im1
          ij = (j-1)*imax + i
	  if (i.gt.1) then
		dx=x(i)-x(i-1)
	  else
		dx=x(2)-x(1)
	  endif
          dvol_tot(1,ij) = 0.0
          dvol_tot(2,ij) = v(ij)*at(ij)*frac*delt*dx
     &                     *(cyl*xi(i) + (1.0 - cyl))
      end do
      i = 1          ! Left
      do j = 1,jm1
          ij = (j-1)*imax + i
	  if (j.gt.1) then
	  	dy=y(j)-y(j-1)
	  else
		dy=y(2)-y(1)
	  endif
          dvol_tot(1,ij) = u(ij)*ar(ij)*frac*delt*dy
     &                      *(cyl*x(i) + (1.0 - cyl))
          dvol_tot(2,ij) = 0.0
      end do

c       ------------------------------------------------------------
        call normal(1,im1,1,jm1,imax,x,y,f,epsmach,cutvof,grad_vof,
     &          ibar2,jbar2,nxy)
c       ------------------------------------------------------------

c       -----------------
        call reconstruct
c       -----------------

        if (mod(dir_flag,2) .eq. 0) then
          call advect_xvof
        else
          call advect_yvof
        endif

	  if (idir.eq.isweeps) then
	    do ntmp=1,imax*jmax
			f(ntmp)=ftmp(ntmp)
	    end do

          if (ninflow.eq.100) then
            do i=isources,isourcee
            do j=jsources,jsourcee
                f((j-1)*imax+i)=1.0
            end do
            end do
		endif

              nnn=1
              ibcflg=1
c            -------------
              call setnf
              call bc
c            -------------
              ibcflg=0
              nnn=0
	  endif

c       ...Check for any volume advection errors

	  total_volume = 0.0
        do 900 j = 1,jm1
          do 900 i = 1,im1
 
            ij = (j-1)*imax + i

	    if (f(ij).gt.1.0) vof_error=f(ij)-1.0
	    if (f(ij).lt.0.0) vof_error=-f(ij)

C            vof_error = dim(f(ij),1.0) + dim(0.0,f(ij))

C.........to empty the small isolate drop or surface cells not adjacent
C.........to the full cell; mass may lose
          if (ac(ij).lt.em6) goto 960
          if (f(ij).gt.1.0e-3) goto 950
	  if (idir.lt.isweeps) goto 950
          if (nf(ij).eq.0.or.nf(ij).eq.6) goto 950
          if ((nf(ij+1).eq.0.and.ac(ij+1).gt.em6).or.(nf(ij-1).eq.0.and.
     &          ac(ij-1).gt.em6).or.(nf(ij+imax).eq.0.and.ac(ij+imax).
     &          gt.em6).or
     &          .(nf(ij-imax).eq.0.and.ac(ij-imax).gt.em6)) goto 950
          f(ij)=0.0
          nf(ij)=6
          goto 900

950	  continue

C.........claim that small fluid are empty cell
          if (f(ij).gt.1.0e-6) goto 955
          f(ij)=0.0
          if (idir.eq.isweeps) nf(ij)=6
C	  nf(ij)=6
          goto 900

C.........claim that nearly full cells are full
955	  if (f(ij).lt.1.0-1.0e-6) goto 960
          f(ij)=1.0
          goto 900

960	  continue

            if (vof_error.gt.cutvof) then
               if (f(ij) .gt. 1.0) then
CC                  write (iotty,50) vof_error,i,j
                  write (19,50) vof_error,i,j
 50               format (' WARNING: VOF - 1 = ',1pe13.5,
     &                    ' in cell (',i4,',',i4,')')
               else if (f(ij) .lt. 0.0) then
CC                  write (iotty,55) vof_error,i,j
                  write (19,55) vof_error,i,j
 55               format (' WARNING: VOF = -',1pe13.5,
     &                    ' in cell (',i4,',',i4,')')
               endif
            endif
            vol_error = vol_error + advect_vol(ij)*vof_error
 
            f(ij) = min(1.0d0,max(0.0d0,f(ij)))
            total_volume = total_volume + f(ij)*cell_vol(ij)

C.......assign velocity to new filled cell at intermediate step

C.......it is found that if calling bc.F and setnf.F for each sweep, 
C.......velocity at surface cell might be altered in the middle way
C.......of VOF advection, which is not reasonable and actually causes
C.......instability for jet problem (near jet-surf intersection).
C.......Efforts should be made to modify bc.F to correct the problem.
C.......some simple trials are made and the problem remains.

C.......An alternative way for this problem is to call bc.F only after the
C.......final sweep is made. However, for the new-filled cell where
C.......velocity is not defined on some faces, this method could be 
C.......inaccurate (though works). To compensate this shortcoming, the
C.......simple treatment is made to assign the velocity for these cells.

	if (idir.lt.isweeps.and.nf(ij).eq.6.and.f(ij).gt.1.0e-6) then
	    if (u(ij).eq.0.0) then
		if (u(ij-1).ne.0.0) then
			u(ij)=u(ij-1)
		else
			if (u(ij+imax).ne.0.0) then
				u(ij)=u(ij+imax)
			else
				u(ij)=u(ij-imax)
			endif
		endif
	    endif
            if (u(ij-1).eq.0.0) then
                if (u(ij).ne.0.0) then
                        u(ij-1)=u(ij)
                else
                        if (u(ij-1+imax).ne.0.0) then
                                u(ij-1)=u(ij-1+imax)
                        else
                                u(ij-1)=u(ij-1-imax)
                        endif
                endif
            endif
            if (v(ij).eq.0.0) then
                if (v(ij-imax).ne.0.0) then
                        v(ij)=v(ij-imax)
                else
                        if (v(ij+1).ne.0.0) then
                                v(ij)=v(ij+1)
                        else
                                v(ij)=v(ij-1)
                        endif
                endif
            endif
            if (v(ij-imax).eq.0.0) then
                if (v(ij).ne.0.0) then
                        v(ij-imax)=v(ij)
                else
                        if (v(ij+1-imax).ne.0.0) then
                                v(ij-imax)=v(ij+1-imax)   
                        else
                                v(ij-imax)=v(ij-1-imax)   
                        endif
                endif
            endif	
	endif
900	continue

      enddo          ! Go back for next sweep

C	idir=0

      write (19,*) ' Cycle ',ncyc,': total fluid volume = ',
     &              total_volume

 100     format (/,'#  Interfaces at cycle ',i4,/)

c     ...Reinitialize the pre-advection volumes

      do ij = 1,nxy
         advect_vol(ij) = cell_vol(ij)
      enddo

      write (19,*) ' Fluid volume advection error = ', vol_error

      return
      end


      SUBROUTINE ADVECT_XVOF

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

      include  "comdk2.h"

      dimension xint(4), yint(4), xs(10), ys(10), pre_vol(nxy)
      integer iface(4)

c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

      do ij = 1,nxy
	 dvol(ij) = 0.0
      end do

      do j = 1,jm1
         do i = 1,im1
            ij = (j-1)*imax + i
            pre_vol(ij) = advect_vol(ij)
         end do
      end do

      do 10 j = 1,jm1
        do 10 i = 1,im1

          ij = (j-1)*imax + i

          isgn = int(sign(1.0d0,u(ij)))
	  ijd = ij + (1 - isgn)/2
          id = i + (1 - isgn)/2
          idm = max(1,id - 1)
          is = max(1,i - isgn)

          if ((f(ijd).lt.cutvof).or.(u(ij) .eq. 0.0)) go to 10

          dvol_xy = (x(i) - flux_pt(1,ij))*(y(j) - y(j-1))
          r_factor = 0.5*(flux_pt(1,ij)+x(i))
          dvol(ij) = dvol_xy*(1.0 + cyl*(r_factor - 1.0))

          if (f(ijd).gt. 1.0 - cutvof) go to 10

          if ((flux_pt(1,ij)-xi0(ijd))*
     &            (xi1(ijd)-flux_pt(1,ij)).ge.0.0) then

	    do n=1,4
	      xint(n) = 0.0
	      yint(n) = 0.0
	      iface(n) = 0.0
            enddo
	    xtemp = x(is)
	    x(is) = flux_pt(1,ij)

            if (xi0(ijd).eq.xi1(ijd) .and.
     &          xi0(ijd).eq.flux_pt(1,ij)) then
               iface(1) = 1
               iface(2) = 1
               xint(1) = xi0(ijd)
               xint(2) = xi1(ijd)
               if (yi0(ijd).lt.yi1(ijd)) then
                  yint(1) = yi0(ijd)
                  yint(2) = yi1(ijd)
               else
                  yint(1) = yi1(ijd)
                  yint(2) = yi0(ijd)
               end if
               go to 5
            end if
               
            ifc = 3 + (1 - isgn)/2
            xint(ifc) = flux_pt(1,ij)
            yint(ifc) = yi0(ijd) + (yi1(ijd)-yi0(ijd))*
     &                (xint(ifc)-xi0(ijd))/(xi1(ijd)-xi0(ijd)+eps)
            iface(ifc) = 1
            if ((xi0(ijd)-x(i))*(flux_pt(1,ij)-xi0(ijd)).ge.0.0) then
              xint(iface0(ijd)) = xi0(ijd)
              yint(iface0(ijd)) = yi0(ijd)
              iface(iface0(ijd)) = 1
            elseif ((xi1(ijd)-x(i))*(flux_pt(1,ij)-xi1(ijd)).ge.0.0) 
     &       then
              xint(iface1(ijd)) = xi1(ijd)
              yint(iface1(ijd)) = yi1(ijd)
              iface(iface1(ijd)) = 1
            endif

 5          continue
            istrt=iface(1)+2*iface(2)*(1-iface(1))
     &             +3*iface(3)*(1-iface(1))*(1-iface(2))
            xs(1) = xint(istrt)
            ys(1) = yint(istrt)
            npt = 2

c           ------------------------------------------------------
            vol = volume(id,j,npt,grad_vof(1,ijd),grad_vof(2,ijd),
     &            x,y,xs,ys,xint,yint,iface,istrt,cyl,ibar2,jbar2)
c           ------------------------------------------------------
            x(is) = xtemp
            dvol(ij) = min(abs(dvol(ij)),max(vol,0.0d0))
            dvol(ij) = isgn*dvol(ij)

          elseif ((flux_pt(1,ij)-x(i))*(xi0(ijd)-flux_pt(1,ij)) .ge. 0.
     &      .and. (flux_pt(1,ij)-x(i))*(xi1(ijd)-flux_pt(1,ij)) .ge. 0.
     &      ) then

            if (grad_vof(1,ijd)*u(ij)*ar(ij).lt. 0.0) dvol(ij)=0.0

          elseif ((xi0(ijd)-x(i))*(flux_pt(1,ij)-xi0(ijd)) .ge. 0.
     &      .and. (xi1(ijd)-x(i))*(flux_pt(1,ij)-xi1(ijd)) .ge. 0.
     &      ) then

            if (grad_vof(1,ijd)*u(ij)*ar(ij).gt. 0.0) then

              dvol(ij) = isgn*f(ijd)*cell_vol(ijd)

            elseif (grad_vof(1,ijd)*u(ij)*ar(ij).lt. 0.0) then

	      do n=1,4
	        xint(n) = 0.0
	        yint(n) = 0.0
	        iface(n) = 0.0
              enddo

     	      xtemp = x(is)
	      x(is) = flux_pt(1,ij)
              xint(iface0(ijd)) = xi0(ijd)
              yint(iface0(ijd)) = yi0(ijd)
              iface(iface0(ijd)) = 1
              xint(iface1(ijd)) = xi1(ijd)
              yint(iface1(ijd)) = yi1(ijd)
              iface(iface1(ijd)) = 1
              istrt=iface(1)+2*iface(2)*(1-iface(1))
     &               +3*iface(3)*(1-iface(1))*(1-iface(2))
              xs(1) = xint(istrt)
              ys(1) = yint(istrt)
              npt = 2
c             ------------------------------------------------------
              vol = volume(id,j,npt,grad_vof(1,ijd),grad_vof(2,ijd),
     &              x,y,xs,ys,xint,yint,iface,istrt,cyl,ibar2,jbar2)
c             ------------------------------------------------------
              x(is) = xtemp
              dvol(ij) = min(abs(dvol(ij)),max(vol,0.0d0))
              dvol(ij) = isgn*dvol(ij)

            endif
          endif

C.......special treatment to prevent the instability
C.......under very special case when wave breaks; force fluid to move
C.......when the free surface reconstruction criterion says no;
C.......(refer to May 27, 1997 note)
	  if (nf(ij-1).eq.6.and.nf(ij+1).eq.6.and.nf(ij).eq.1.and.
     &		dvol(ij).le.1.0e-10) 
     &          dvol(ij) = f(ij)*dvol_xy*(1.0 + cyl*(r_factor - 1.0))

   10 continue

C.......treatment for open (3), applied pressure (5), and specific (6) BC
C.......as well as periodic BC (by lin on 12/4/98)
      do j=1,jm1
        il=1
        ir=im1
        ijl=(j-1)*imax+il
        ijr=(j-1)*imax+ir

        if (kl.eq.4) dvol(ijl)=dvol(ijr-1)
        if (u(ijl).gt.0.0.and.f(ijl).lt.1.0-cutvof.and.
     &		(kl.eq.3.or.kl.eq.5.or.kl.eq.6))
     &          dvol(ijl)=u(ijl)*delt*f(ijl)*dely(j)

        if (kr.eq.4) dvol(ijr)=dvol(ijl+1)
        if (u(ijr).lt.0.0.and.f(ijr).lt.1.0-cutvof.and.
     &		(kr.eq.3.or.kr.eq.5.or.kr.eq.6))
     &          dvol(ijr)=u(ijr)*delt*f(ijr+1)*dely(j)
      end do

c     ...Update the volume fractions

	xms=0.0
      do j = 1,jm1
         do i = 1,im1

           ij = (j-1)*imax + i
           imj = ij - 1
	ipj=ij+1
	ijp=ij+imax
	ijm=ij-imax

	   if (i.gt.1) then

           	ftmp(ij) = (ftmp(ij)*pre_vol(ij) - 
     &             (dvol(ij) - dvol(imj)))/advect_vol(ij)

	xms=xms+f(ij)

	   else
		ftmp(ij) = ftmp(ij)*pre_vol(ij)/advect_vol(ij)
	   endif
         end do
      end do

      return
      end


      SUBROUTINE ADVECT_YVOF

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

      include  "comdk2.h" 

      dimension xint(4), yint(4), xs(10), ys(10), pre_vol(nxy)
      integer iface(4)

c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

      do ij = 1,nxy
	 dvol(ij) = 0.0
      end do

      do j = 1,jm1
         do i = 1,im1
            ij = (j-1)*imax + i
            pre_vol(ij) = advect_vol(ij)
         end do 
      end do

      do 10 j = 1,jm1
        do 10 i = 1,im1

          ij = (j-1)*imax + i

          jsgn = int(sign(1.0d0,v(ij)))
	  ijd = ij + imax*(1 - jsgn)/2
          jd = j + (1 - jsgn)/2
          jdm = max(1,jd - 1)
          js = max(1,j - jsgn)

          dvol(ij) = 0.0

          if ((f(ijd).lt.cutvof).or.(v(ij).eq.0.0)) go to 10

          dvol_xy = (x(i) - x(i-1))*(y(j) - flux_pt(2,ij))
          r_factor = xi(i)
          dvol(ij) = dvol_xy*(1.0 + cyl*(r_factor - 1.0))

          if (f(ijd).gt. 1.0 - cutvof) go to 10

          if ((flux_pt(2,ij)-yi0(ijd))*
     &            (yi1(ijd)-flux_pt(2,ij)).ge.0.0) then

	    do n=1,4
	      xint(n) = 0.0
	      yint(n) = 0.0
	      iface(n) = 0.0
            enddo
	    ytemp = y(js)
	    y(js) = flux_pt(2,ij)

            if (yi0(ijd).eq.yi1(ijd) .and.
     &          yi0(ijd).eq.flux_pt(2,ij)) then
               iface(3) = 1
               iface(4) = 1
               yint(3) = yi0(ijd)
               yint(4) = yi1(ijd)
               if (xi0(ijd).lt.xi1(ijd)) then
                  xint(3) = xi0(ijd)
                  xint(4) = xi1(ijd)
               else
                  xint(3) = xi1(ijd)
                  xint(4) = xi0(ijd)
               end if
               go to 5
            end if
               
            ifc = 1 + (1 - jsgn)/2
            yint(ifc) = flux_pt(2,ij)
            xint(ifc) = xi0(ijd) + (xi1(ijd)-xi0(ijd))*
     &                (yint(ifc)-yi0(ijd))/(yi1(ijd)-yi0(ijd)+eps)
            iface(ifc) = 1
            if ((yi0(ijd)-y(j))*(flux_pt(2,ij)-yi0(ijd)).ge.0.0) then
              xint(iface0(ijd)) = xi0(ijd)
              yint(iface0(ijd)) = yi0(ijd)
              iface(iface0(ijd)) = 1
            elseif ((yi1(ijd)-y(j))*(flux_pt(2,ij)-yi1(ijd)).ge.0.0) 
     &       then
              xint(iface1(ijd)) = xi1(ijd)
              yint(iface1(ijd)) = yi1(ijd)
              iface(iface1(ijd)) = 1
            endif

 5          continue
            istrt=iface(1)+2*iface(2)*(1-iface(1))
     &             +3*iface(3)*(1-iface(1))*(1-iface(2))
            xs(1) = xint(istrt)
            ys(1) = yint(istrt)
            npt = 2
c           ------------------------------------------------------
            vol = volume(i,jd,npt,grad_vof(1,ijd),grad_vof(2,ijd),
     &            x,y,xs,ys,xint,yint,iface,istrt,cyl,ibar2,jbar2)
c           ------------------------------------------------------
            y(js) = ytemp
            dvol(ij) = min(abs(dvol(ij)),max(vol,0.0d0))
            dvol(ij) = jsgn*dvol(ij)

          elseif ((flux_pt(2,ij)-y(j))*(yi0(ijd)-flux_pt(2,ij)) .ge. 0.
     &      .and. (flux_pt(2,ij)-y(j))*(yi1(ijd)-flux_pt(2,ij)) .ge. 0.) 
     &     then

            if (grad_vof(2,ijd)*v(ij)*at(ij).lt. 0.0) dvol(ij)=0.0

          elseif ((yi0(ijd)-y(j))*(flux_pt(2,ij)-yi0(ijd)) .ge. 0.
     &      .and. (yi1(ijd)-y(j))*(flux_pt(2,ij)-yi1(ijd)) .ge. 0.
     &     ) then

            if (grad_vof(2,ijd)*v(ij)*at(ij).gt. 0.0) then

              dvol(ij) = jsgn*f(ijd)*cell_vol(ijd)

            elseif (grad_vof(2,ijd)*v(ij)*at(ij).lt. 0.0) then

	      do n=1,4
	        xint(n) = 0.0
	        yint(n) = 0.0
	        iface(n) = 0.0
              enddo

     	      ytemp = y(js)
	      y(js) = flux_pt(2,ij)
              xint(iface0(ijd)) = xi0(ijd)
              yint(iface0(ijd)) = yi0(ijd)
              iface(iface0(ijd)) = 1
              xint(iface1(ijd)) = xi1(ijd)
              yint(iface1(ijd)) = yi1(ijd)
              iface(iface1(ijd)) = 1
              istrt=iface(1)+2*iface(2)*(1-iface(1))
     &               +3*iface(3)*(1-iface(1))*(1-iface(2))
              xs(1) = xint(istrt)
              ys(1) = yint(istrt)
              npt = 2
c             ------------------------------------------------------
              vol = volume(i,jd,npt,grad_vof(1,ijd),grad_vof(2,ijd),
     &              x,y,xs,ys,xint,yint,iface,istrt,cyl,ibar2,jbar2)
c             ------------------------------------------------------
              y(js) = ytemp
              dvol(ij) = min(abs(dvol(ij)),max(vol,0.0d0))
              dvol(ij) = jsgn*dvol(ij)

            endif

	  endif

C.......special treatment to prevent the instability   
C.......under very special case when wave breaks; force fluid to move
C.......when the free surface reconstruction criterion says no;
C.......(refer to May 27, 1997 note)
          if (nf(ij-imax).eq.6.and.nf(ij+imax).eq.6.and.nf(ij).eq.3.and.
     &          dvol(ij).le.1.0e-10)
     &          dvol(ij) = f(ij)*dvol_xy*(1.0 + cyl*(r_factor - 1.0))

   10 continue

	xms=0.0
c     ...Update the volume fractions

      do 20 j = 1,jm1
        do 20 i = 1,im1

          ij = (j-1)*imax + i
          ijm = ij - imax
	ijp=ij+imax

	  if (j.gt.1) then
          	ftmp(ij) = (ftmp(ij)*pre_vol(ij) - 
     &             (dvol(ij) - dvol(ijm)))/advect_vol(ij)
		xms=xms+f(ij)
	  else
		ftmp(ij) = ftmp(ij)*pre_vol(ij)/advect_vol(ij)
	  endif

   20 continue

      return
      end

      SUBROUTINE INPUT

      implicit real*8 (a-h,o-z)

      include  "comdk2.h" 

      namelist /vof/ tol, itmax0,
     &               periodic, cutvof, four_sweep,
     &               io_interval, slic

      data one /1.0d0/, zero /0.0d0/

	cutvof=1.0d-06
	eps=1.0d-20
	tol=1.0d-12
	epsmach=1.0d-15
	four_sweep=.false.
	itmax0=50
	periodic=.false.
	io_interval=2500
	preset=-123456789.0
	tracer=.false.
	slic=.false.

c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

c     -------------------------------------------------
C      open (4, file = 'Youngs.inp', status = 'old')
C      open (7, file = 'circle.small.20.dmp',
C     &      status = 'old', form='unformatted')
C      open (8, file = 'Youngs.out', status = 'unknown')
      open (19, file = 'Youngs.err', status = 'unknown')
c     -------------------------------------------------

C      read (4,vof)
C      read (7) f,x,y,xi,yj

c     ...Initialize relevant scalars and arrays

      do ij = 1,nxy
        grad_vof(1,ij) = 0.0
        grad_vof(2,ij) = 0.0
        dvol(ij) = 0.0
        dvol_tot(1,ij) = 0.0
        dvol_tot(2,ij) = 0.0
        iface0(ij) = 0
        xi0(ij) = 0.0
        yi0(ij) = 0.0
        iface1(ij) = 0
        xi1(ij) = 0.0
        yi1(ij) = 0.0
        cell_vol(ij) = 0.0
        advect_vol(ij) = 0.0
      end do

c     ...Compute flux vols, etc.

      total_volume = 0.0
      if (four_sweep) then
         frac = 0.5
         isweeps = 4
      else
         frac = 1.0
         isweeps = 2
      endif

	idir=isweeps

      do j = 2,jm1
        do i = 2,im1
          ij = (j-1)*imax + i
          ravg = cyl*0.5*(x(i)+x(i-1))+1.-cyl
          cell_vol(ij) = ravg*(x(i)-x(i-1))*(y(j)-y(j-1))
          total_volume = total_volume + f(ij)*cell_vol(ij)
        end do
      end do

      write (19,*) ' Initial fluid volume = ', total_volume

      return
      end


      SUBROUTINE INTERSECT(i,j,xs,ys,x,y,t1x,t1y,t2x,
     &           t2y,xint,yint,iface,istrt,ilook,ibar2,jbar2)

      implicit real*8 (a-h,o-z)

      dimension xs(10),ys(10),x(ibar2),y(jbar2),xint(4),yint(4),
     &          iface(4),ilook(4)

      do 100 n=1,4
        iface(n)=0
  100 continue
c
        if (j.gt.1) then
                ytmp=y(j-1)
        else
                ytmp=y(1)-(y(2)-y(1))
        endif
        if (i.gt.1) then
                xtmp=x(i-1)
        else
                xtmp=x(1)-(x(2)-x(1))
        endif
c
c.... Bottom Face
c
      if (ilook(1).eq.0) go to 25
      xintb1=xs(1)-(ys(1)-ytmp)*t1x/t1y
      if ((xintb1-xtmp)*(x(i)-xintb1).ge.0.0) then
        xint(1)=xintb1
	yint(1)=ytmp
        iface(1)=1
        go to 30
      endif
c
   25 xintb2=xs(1)-(ys(1)-ytmp)*t2x/t2y

      if ((xintb2-xtmp)*(x(i)-xintb2).ge.0.0) then
        xint(1)=xintb2
	yint(1)=ytmp
        iface(1)=1
      endif
c
c.... Top Face
c
   30 if (ilook(2).eq.0) go to 35
      xintt1=xs(1)+(y(j)-ys(1))*t1x/t1y
      if ((xintt1-xtmp)*(x(i)-xintt1).ge.0.0) then
        xint(2)=xintt1
        yint(2)=y(j)
        iface(2)=1
        go to 40
      endif
c
   35 xintt2=xs(1)+(y(j)-ys(1))*t2x/t2y
      if ((xintt2-xtmp)*(x(i)-xintt2).ge.0.0) then
        xint(2)=xintt2
        yint(2)=y(j)
        iface(2)=1
      endif
c
c.... Right Face
c
   40 if (ilook(4).eq.0) go to 5
      yintr1=ys(1)+(x(i)-xs(1))*t1y/t1x
      if ((yintr1-ytmp)*(y(j)-yintr1).ge.0.0) then
        xint(4)=x(i)
        yint(4)=yintr1
        iface(4)=1
        go to 10
      endif
c
    5 yintr2=ys(1)+(x(i)-xs(1))*t2y/t2x
      if ((yintr2-ytmp)*(y(j)-yintr2).ge.0.0) then
        xint(4)=x(i)
        yint(4)=yintr2
        iface(4)=1
      endif
c
c.... Left Face
c
   10 if (ilook(3).eq.0) go to 15
      yintl1=ys(1)-(xs(1)-xtmp)*t1y/t1x
      if ((yintl1-ytmp)*(y(j)-yintl1).ge.0.0) then
        xint(3)=xtmp
        yint(3)=yintl1
        iface(3)=1
        go to 20
      endif
c
   15 yintl2=ys(1)-(xs(1)-xtmp)*t2y/t2x
      if ((yintl2-ytmp)*(y(j)-yintl2).ge.0.0) then
        xint(3)=xtmp
        yint(3)=yintl2
        iface(3)=1
      endif
c
   20 xs(1)=iface(1)*xint(1)+iface(2)*(1-iface(1))*xint(2)
     &          +iface(3)*(1-iface(1))*(1-iface(2))*xint(3)
      ys(1)=iface(1)*yint(1)+iface(2)*(1-iface(1))*yint(2)
     &          +iface(3)*(1-iface(1))*(1-iface(2))*yint(3)
      istrt=iface(1)+2*iface(2)*(1-iface(1))
     &         +3*iface(3)*(1-iface(1))*(1-iface(2))

      return
      end


      SUBROUTINE NORMAL(ib,ie,jb,je,iskip,x,y,f,eps,cutvof,grad_vof,
     &	ibar2,jbar2,nxy)

      implicit real*8 (a-h,o-z)

      dimension x(ibar2), y(jbar2), f(nxy), grad_vof(2,nxy)
      
      do 10 j = jb, je
        do 10 i = ib, ie

          ij = (j-1)*iskip + i

          grad_vof(1,ij) = 0.0d0
          grad_vof(2,ij)   = 0.0d0

          if ((f(ij).lt.cutvof).or.(f(ij).gt. 1.0 - cutvof)) go to 10

          ipj = ij+1
	  ijp = ij+iskip
          ipjp = ipj+iskip
          ip = i + 1
          jp = j + 1
          ipjm = ipj-iskip
          imj = ij-1
          imjp = imj+iskip
          imjm = imj-iskip
          ijm = ij-iskip
          im = i - 1
          jm = j - 1
          imm = im - 1
          jmm = jm - 1
	  if (i.eq.1.and.j.gt.1) then
			imj = ij
                	imjp = ij+iskip
			im = i
			imm = im - 1
			imjm=ij-iskip
	  endif
	  if (j.eq.1.and.i.gt.1) then
	                ijm = ij
			ipjm = ij+1	
			jm = j
			jmm = jm -1
			imjm=ij-1
	  endif
	  if (i.eq.1.and.j.eq.1) then
			imjm=ij
			ijm=ij
			imj=ij
			ipjm=ij+1
			imjp=ij+iskip
			im=i
			imm=im-1
			jm=j
			jmm=jm-1
	  endif
	  if (i.gt.1) then
	        dx_i = 0.5*(x(i)-x(im))
	  else
		dx_i=0.5*(x(2)-x(1))
	  endif
          if (imm .gt. 0) then
             dx_im = 0.5*(x(im)-x(imm))
          else
             dx_im = dx_i
          endif
          dx_ip = 0.5*(x(ip)-x(i))

	  if (j.gt.1) then
          	dy_j = 0.5*(y(j)-y(jm))
	  else
		dy_j=0.5*(y(2)-y(1))
	  endif
          if (jmm .gt. 0) then
             dy_jm = 0.5*(y(jm)-y(jmm))
          else
             dy_jm = dy_j
          endif

          dy_jp = 0.5*(y(jp)-y(j))

          voftl=(dx_i*f(imjp)+dx_im*f(ijp))/(dx_i+dx_im)
          vofml=(dx_i*f(imj)+dx_im*f(ij))/(dx_i+dx_im)
          vofbl=(dx_i*f(imjm)+dx_im*f(ijm))/(dx_i+dx_im)
          voftr=(dx_i*f(ipjp)+dx_ip*f(ijp))/(dx_i+dx_ip)
          vofmr=(dx_i*f(ipj)+dx_ip*f(ij))/(dx_i+dx_ip)
          vofbr=(dx_i*f(ipjm)+dx_ip*f(ijm))/(dx_i+dx_ip)

          voftrc=(dy_j*voftr+dy_jp*vofmr)/(dy_j+dy_jp)
          vofbrc=(dy_j*vofbr+dy_jm*vofmr)/(dy_j+dy_jm)
          voftlc=(dy_j*voftl+dy_jp*vofml)/(dy_j+dy_jp)
          vofblc=(dy_j*vofbl+dy_jm*vofml)/(dy_j+dy_jm)

          grad_vof(1,ij) = 0.25*(voftrc+vofbrc-voftlc-vofblc)/dx_i
          grad_vof(2,ij) = 0.25*(voftrc+voftlc-vofblc-vofbrc)/dy_j

          fmag = sqrt(grad_vof(1,ij)**2+grad_vof(2,ij)**2)

          if (fmag.gt.eps) then
             grad_vof(2,ij) = grad_vof(2,ij)/fmag
             grad_vof(1,ij) = grad_vof(1,ij)/fmag
             if (dabs(grad_vof(2,ij)).lt.eps) grad_vof(2,ij) = 0.0d0
             if (dabs(grad_vof(1,ij)).lt.eps) grad_vof(1,ij) = 0.0d0
          endif

   10 continue

      return
      end


      SUBROUTINE RECONSTRUCT

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

      include  "comdk2.h" 

      dimension xint(4), yint(4), xs(10), ys(10)
      integer iface(4)

c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

      nmix = 0
      avgiter = 0.0

C.....initiate the xa,yx,xint,yint
	do i=1,4
	  xint(i)=0.0
	  yint(i)=0.0
	  ilook(i)=1
	end do
	do i=1,10
	  xs(i)=0.0
	  ys(i)=0.0
	end do
      do 50 j = 1,jm1
        do 50 i = 1,im1

          ij = (j-1)*imax + i

          if ((f(ij).lt.cutvof).or.(f(ij).gt. 1.0 - cutvof)) go to 50
	  if (i.eq.1.and.j.eq.1) then
		f(ij)=1.0
	  	goto 50
	  endif
            nmix = nmix + 1
            ijmix(nmix) = ij

            if (slic) then
               call slic_recon
               go to 50
            end if

            t1x = grad_vof(2,ij)  +  eps
            t2x = -t1x
            t1y = -grad_vof(1,ij)  +  eps
            t2y = -t1y
	    if (i.gt.1) then
            	dx = x(i)-x(i-1)
	    else
		dx=x(2)-x(1)
	    endif
	    if (j.gt.1) then
            	dy = y(j)-y(j-1)
	    else
		dy=y(2)-y(1)
	    endif
            iter1 = 0
            xi0(ij) = 0.0
            xi1(ij) = 0.0
            yi0(ij) = 0.0
            yi1(ij) = 0.0

            if (grad_vof(1,ij).ge.0.0.and.grad_vof(2,ij).ge.0.0) then
              x1 = x(i)
              y1 = y(j)
	      if (i.gt.1) then
              	x2 = x(i-1)
	      else
	      	x2=x(1)-dx
	      endif
	      if (j.gt.1) then
              	y2 = y(j-1)
	      else
		y2=y(1)-dy
	      endif
              xprime = xi(i) + dx*(0.5-f(ij))
              yprime = yj(j) + dy*(0.5-f(ij))
              ilook(2) = 0
              ilook(3) = 0
            elseif (grad_vof(1,ij).ge.0.0.and.grad_vof(2,ij).lt.0.0) 
     &       then
              x1 = x(i)
	      if (j.gt.1) then
              	y1 = y(j-1)
	      else
		y1=y(1)-dy
	      endif
	      if (i.gt.1) then
              	x2 = x(i-1)
	      else
		x2=x(1)-dx
	      endif
              y2 = y(j)
              xprime = xi(i) + dx*(0.5-f(ij))
              yprime = yj(j)-dy*(0.5-f(ij))
              ilook(2) = 0
              ilook(4) = 0
            elseif (grad_vof(1,ij).lt.0.0.and.grad_vof(2,ij).ge.0.0) 
     &       then
	      if (i.gt.1) then
              	x1 = x(i-1)
	      else
		x1=x(1)-dx
	      endif
              y1 = y(j)
              x2 = x(i)
	      if (j.gt.1) then
              	y2 = y(j-1)
	      else
		y2=y(1)-dy
	      endif
              xprime = xi(i)-dx*(0.5-f(ij))
              yprime = yj(j) + dy*(0.5-f(ij))
              ilook(1) = 0
              ilook(3) = 0
            else
	      if (i.gt.1) then
              	x1 = x(i-1)
	      else
		x1=x(1)-dx
	      endif
	      if (j.gt.1) then
              	y1 = y(j-1)
	      else
		y1=y(1)-dy
	      endif
              x2 = x(i)
              y2 = y(j)
              xprime = xi(i)-dx*(0.5-f(ij))
              yprime = yj(j)-dy*(0.5-f(ij))
              ilook(1) = 0
              ilook(4) = 0
            endif

            func1 = f(ij)
            func2 = f(ij)-1.0
            func = func1
            xcapp = 0.0
            ycapp = 0.0
            capq = 0.0
c
   10       if (iter1.gt.itmax0) then
CC              write (iotty,200) itmax0,ncyc,ij,f(ij),vof
              write (19,200) itmax0,ncyc,ij,f(ij),vof
              go to 50
            endif

            if (dabs(func).lt.tol + epsmach) go to 20
c
            iter1 = iter1 + 1
            xmid = 0.5*(x1 + x2)
            ymid = 0.5*(y1 + y2)
            xprime = xprime + capq*xcapp/(capq + eps)**2
            yprime = yprime + capq*ycapp/(capq + eps)**2
            if ((x1-xprime)*(x2-xprime).gt.0.0) xprime = xmid
            if ((y1-yprime)*(y2-yprime).gt.0.0) yprime = ymid
            xs(1) = xprime
            ys(1) = yprime
c
c           ------------------------------------------------
            call intersect(i,j,xs,ys,x,y,t1x,t1y,t2x,t2y,
     &                     xint,yint,iface,istrt,ilook,ibar2,jbar2)
c           ------------------------------------------------
            npt = 2
c
c           ----------------------------------------------------------
            vof = volume(i,j,npt,grad_vof(1,ij),grad_vof(2,ij),x,y,xs,
     &         ys,xint,yint,iface,istrt,cyl,ibar2,jbar2)/cell_vol(ij)
c           ----------------------------------------------------------

        if (vof.gt.-100.and.vof.lt.100) goto 555

555	continue
            func = f(ij)-vof

C            if (func.eq.0.0) go to 20
            if (abs(func).le.1.0e-10) go to 20

            capr = func/func1
            caps = func/func2
            capt = func2/func1
            capq = (capt-1.)*(capr-1.)*(caps-1.)
            xcapp = caps*(capt*(capr-capt)*(x1-xprime)-
     &             (1.-capr)*(xprime-x2))
            ycapp = caps*(capt*(capr-capt)*(y1-yprime)-
     &             (1.-capr)*(yprime-y2))

            if (func.gt.0.0) then
              x1 = xprime
              y1 = yprime
              func1 = func
              go to 10
            else
              x2 = xprime
              y2 = yprime
              func2 = func
              go to 10
            endif

   20       continue

            avgiter = avgiter + iter1

            xi0(ij) = xint(istrt)
            yi0(ij) = yint(istrt)
            iface0(ij) = istrt
            do 15 n = 1,4
              ilook(n) = 1
              if ((n .ne. istrt).and.(iface(n) .ne. 0)) then
                xi1(ij) = xint(n)
                yi1(ij) = yint(n)
                iface1(ij) = n
              endif
   15       continue

   50 continue

      if (nmix .gt. 0) then
         avgiter = avgiter/nmix
         write(19,300) avgiter
 300     format(' Interface reconstruction averaged ',1pe13.5,
     &          ' iterations')
      endif

      return
  200 format(' WARNING: Interface reconstruction convergence failure',
     &       ' after ',i4,' iterations',/,10x,'cycle ',i5,', cell ',i5,
     &       ' VOF_actual = ',1pe13.5,' VOF_computed = ',1pe13.5)
      end


      SUBROUTINE SLIC_RECON

c-----------------------------------------------------------------------
c     ...Reconstruct the interface with a modified-SLIC formulation
c-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

      include  "comdk2.h" 

c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

c     ...Get cell indices
      
      ij   = (j - 1)*imax + i
      ipj  = ij + 1
      imj  = ij - 1
      ijp  = ij + imax
      ijm  = ij - imax
      ipjp = ipj + imax
      imjp = imj + imax
      ipjm = ipj - imax
      imjm = imj - imax

c     ...Compute the center of mass coords
c        based on surrounding volume data

      sum_vol = f(ij)*cell_vol(ij) +
     &          f(ipj)*cell_vol(ipj) +
     &          f(imj)*cell_vol(imj) +
     &          f(ijp)*cell_vol(ijp) +
     &          f(ijm)*cell_vol(ijm) +
     &          f(ipjp)*cell_vol(ipjp) +
     &          f(ipjm)*cell_vol(ipjm) +
     &          f(imjp)*cell_vol(imjp) +
     &          f(imjm)*cell_vol(imjm)

      x_com = xi(i)*f(ij)*cell_vol(ij) +
     &        xi(i+1)*f(ipj)*cell_vol(ipj) +
     &        xi(i-1)*f(imj)*cell_vol(imj) +
     &        xi(i)*f(ijp)*cell_vol(ijp) +
     &        xi(i)*f(ijm)*cell_vol(ijm) +
     &        xi(i+1)*f(ipjp)*cell_vol(ipjp) +
     &        xi(i+1)*f(ipjm)*cell_vol(ipjm) +
     &        xi(i-1)*f(imjp)*cell_vol(imjp) +
     &        xi(i-1)*f(imjm)*cell_vol(imjm)
      x_com = x_com/(sum_vol + eps)

      y_com = yj(j)*f(ij)*cell_vol(ij) +
     &        yj(j)*f(ipj)*cell_vol(ipj) +
     &        yj(j)*f(imj)*cell_vol(imj) +
     &        yj(j+1)*f(ijp)*cell_vol(ijp) +
     &        yj(j-1)*f(ijm)*cell_vol(ijm) +
     &        yj(j+1)*f(ipjp)*cell_vol(ipjp) +
     &        yj(j-1)*f(ipjm)*cell_vol(ipjm) +
     &        yj(j+1)*f(imjp)*cell_vol(imjp) +
     &        yj(j-1)*f(imjm)*cell_vol(imjm)
      y_com = y_com/(sum_vol + eps)

c     ...Now go through each case

      if ((abs(grad_vof(1,ij)).lt.eps) .and.
     &    (abs(grad_vof(2,ij)).lt.eps)) then      ! Degenerate interface
         write (19,*) 'FATAL: Degenerate interface conditions',
     &               ' in cell (',i,',',j,')'
         stop
      else if (abs(grad_vof(1,ij)).lt.eps) then   ! Horizontal interface
         iface0(ij) = 4
         iface1(ij) = 3
         xi0(ij) = x(i)
         xi1(ij) = x(i-1)
         if (y_com.ge.yj(j) .and. grad_vof(2,ij).gt.0.0) then
            yint = y(j) - f(ij)*(y(j) - y(j-1))
         else if (y_com.lt.yj(j) .and. grad_vof(2,ij).lt.0.0) then
            yint = y(j-1) + f(ij)*(y(j) - y(j-1))
         else
            write (iotty,*) 'FATAL: Could not find y-coord interface',
     &               ' position in cell (',i,',',j,')'
            stop
         end if
         yi0(ij) = yint
         yi1(ij) = yint
         if (f(ijp).lt.eps .and. f(ijm).lt.eps) then
            write (iotty,*) 'WARNING: Double horizontal interface',
     &               ' in cell (',i,',',j,') approximated as 1.0'
         end if
      else if (abs(grad_vof(2,ij)).lt.eps) then   ! Vertical interface
         iface0(ij) = 2
         iface1(ij) = 1
         yi0(ij) = y(j)
         yi1(ij) = y(j-1)
         if (x_com.ge.xi(i) .and. grad_vof(1,ij).gt.0.0) then
            xint = x(i) - f(ij)*(x(i) - x(i-1))
         else if (x_com.lt.xi(i) .and. grad_vof(1,ij).lt.0.0) then
            xint = x(i-1) + f(ij)*(x(i) - x(i-1))
         else
            write (iotty,*) 'FATAL: Could not find x-coord interface',
     &               ' position in cell (',i,',',j,')'
            stop
         end if
         xi0(ij) = xint
         xi1(ij) = xint
         if (f(ipj).lt.eps .and. f(imj).lt.eps) then
            write (iotty,*) 'WARNING: Double vertical interface',
     &               ' in cell (',i,',',j,') approximated as 1.0'
         end if
      else if (abs(grad_vof(1,ij)/grad_vof(2,ij)).le.1.0) then   ! Horizontal interface
         iface0(ij) = 4
         iface1(ij) = 3
         xi0(ij) = x(i)
         xi1(ij) = x(i-1)
         if (y_com.ge.yj(j) .and. grad_vof(2,ij).gt.0.0) then
            yint = y(j) - f(ij)*(y(j) - y(j-1))
         else if (y_com.lt.yj(j) .and. grad_vof(2,ij).lt.0.0) then
            yint = y(j-1) + f(ij)*(y(j) - y(j-1))
         else
            write (iotty,*) 'FATAL: Could not find y-coord interface',
     &               ' position in cell (',i,',',j,')'
            stop
         end if
         yi0(ij) = yint
         yi1(ij) = yint
         if (f(ijp).lt.eps .and. f(ijm).lt.eps) then
            write (iotty,*) 'WARNING: Double horizontal interface',
     &               ' in cell (',i,',',j,') approximated as 1.0'
         end if
      else if (abs(grad_vof(1,ij)/grad_vof(2,ij)).gt.1.0) then   ! Vertical interface
         iface0(ij) = 2
         iface1(ij) = 1
         yi0(ij) = y(j)
         yi1(ij) = y(j-1)
         if (x_com.ge.xi(i) .and. grad_vof(1,ij).gt.0.0) then
            xint = x(i) - f(ij)*(x(i) - x(i-1))
         else if (x_com.lt.xi(i) .and. grad_vof(1,ij).lt.0.0) then
            xint = x(i-1) + f(ij)*(x(i) - x(i-1))
         else
            write (iotty,*) 'FATAL: Could not find x-coord interface',
     &               ' position in cell (',i,',',j,')'
            stop
         end if
         xi0(ij) = xint
         xi1(ij) = xint
         if (f(ipj).lt.eps .and. f(imj).lt.eps) then
            write (iotty,*) 'WARNING: Double vertical interface',
     &               ' in cell (',i,',',j,') approximated as 1.0'
         end if
      else
         write (iotty,*) 'FATAL: Unknown SLIC interface reconstruction',
     &               ' case encountered in cell (',i,',',j,')'
         stop
      end if   

      return
      end


      FUNCTION VOLUME(i,j,npt,fx,fy,x,y,xs,ys,xint,yint,iface,istrt,
     &	cyl,ibar2,jbar2)

      implicit real*8 (a-h,o-z)

      real*8 x(ibar2),y(jbar2),xs(10),ys(10),xint(4),yint(4),fx,fy,cyl
      integer*4 iface(4), i, j, npt, istrt

	if (i.gt.1) then
		xtmp=x(i-1)
	else
		xtmp=x(1)-(x(2)-x(1))
	endif
	if (j.gt.1) then
		ytmp=y(j-1)
	else
		ytmp=y(1)-(y(2)-y(1))
	endif
      go to (10,20,30) istrt
        if (istrt.gt.3.or.istrt.eq.0) write(19,*)'5',i,j,istrt,iface,
     &	npt,xs,ys,xint,yint
	close(19)
      stop 'Invalid value of ISTRT in routine VOLUME'
c
c.... Bottom intersection
c
   10 if (fx.gt.0.0) go to 100
c
c.... Clockwise Traverse
c
      xs(npt) = xtmp
      ys(npt) = ytmp
c
      npt = npt+1
      xs(npt) = (1-iface(3))*xtmp+iface(3)*xint(3)
      ys(npt) = (1-iface(3))*y(j)+iface(3)*yint(3)
      if (iface(3).eq.1) go to 9999
c
      npt = npt+1
      xs(npt) = (1-iface(2))*x(i)+iface(2)*xint(2)
      ys(npt) = (1-iface(2))*y(j)+iface(2)*yint(2)
      if (iface(2).eq.1) go to 9999
c
      npt = npt+1
      xs(npt) = xint(4)
      ys(npt) = yint(4)
      go to 9999
c
c.... Counterclockwise Traverse
c
  100 xs(npt) = x(i)
      ys(npt) = ytmp
c
      npt = npt+1
      xs(npt) = (1-iface(4))*x(i)+iface(4)*xint(4)
      ys(npt) = (1-iface(4))*y(j)+iface(4)*yint(4)
      if (iface(4).eq.1) go to 9999
c
      npt = npt+1
      xs(npt) = (1-iface(2))*xtmp+iface(2)*xint(2)
      ys(npt) = (1-iface(2))*y(j)+iface(2)*yint(2)
      if (iface(2).eq.1) go to 9999
c
      npt = npt+1
      xs(npt) = xint(3)
      ys(npt) = yint(3)
      go to 9999
c
c.... Top intersection
c
   20 if (fx.gt.0.0) go to 200
c
c.... Counterclockwise Traverse
c
      xs(npt) = xtmp
      ys(npt) = y(j)
c
      npt = npt+1
      xs(npt) = (1-iface(3))*xtmp+iface(3)*xint(3)
      ys(npt) = (1-iface(3))*ytmp+iface(3)*yint(3)
      if (iface(3).eq.1) go to 9999
c
      npt = npt+1
      xs(npt) = (1-iface(1))*x(i)+iface(1)*xint(1)
      ys(npt) = (1-iface(1))*ytmp+iface(1)*yint(1)
      if (iface(1).eq.1) go to 9999
c
      npt = npt+1
      xs(npt) = xint(4)
      ys(npt) = yint(4)
      go to 9999
c
c.... Clockwise Traverse
c
  200 xs(npt) = x(i)
      ys(npt) = y(j)
c
      npt = npt+1
      xs(npt) = (1-iface(4))*x(i)+iface(4)*xint(4)
      ys(npt) = (1-iface(4))*ytmp+iface(4)*yint(4)
      if (iface(4).eq.1) go to 9999
c
      npt = npt+1
      xs(npt) = (1-iface(1))*xtmp+iface(1)*xint(1)
      ys(npt) = (1-iface(1))*ytmp+iface(1)*yint(1)
      if (iface(1).eq.1) go to 9999
c
      npt = npt+1
      xs(npt) = xint(3)
      ys(npt) = yint(3)
      go to 9999
c
c.... Left intersection
c
   30 if (fy.gt.0.0) go to 300
c
c.... Counterclockwise Traverse
c
      xs(npt) = xtmp
      ys(npt) = ytmp
c
      npt = npt+1
      xs(npt) = (1-iface(1))*x(i)+iface(1)*xint(1)
      ys(npt) = (1-iface(1))*ytmp+iface(1)*yint(1)
      if (iface(1).eq.1) go to 9999
c
      npt = npt+1
      xs(npt) = (1-iface(4))*x(i)+iface(4)*xint(4)
      ys(npt) = (1-iface(4))*y(j)+iface(4)*yint(4)
      if (iface(4).eq.1) go to 9999
c
      npt = npt+1
      xs(npt) = xint(2)
      ys(npt) = yint(2)
      go to 9999
c
c.... Clockwise Traverse
c
  300 continue
	xs(npt) = xtmp
	ys(npt) = y(j)
c
      npt = npt+1
      xs(npt) = (1-iface(2))*x(i)+iface(2)*xint(2)
      ys(npt) = (1-iface(2))*y(j)+iface(2)*yint(2)
      if (iface(2).eq.1) go to 9999
c
      npt = npt+1
      xs(npt) = (1-iface(4))*x(i)+iface(4)*xint(4)
      ys(npt) = (1-iface(4))*ytmp+iface(4)*yint(4)
      if (iface(4).eq.1) go to 9999
c
      npt = npt+1
      xs(npt) = xint(1)
      ys(npt) = yint(1)
c
 9999 npt = npt+1
      xs(npt) = xs(1)
      ys(npt) = ys(1)
c
      vol = 0.0
      do 40 ip = 1,npt-1
        vol_xy = 0.5*(xs(ip)*ys(ip+1)-xs(ip+1)*ys(ip))
        vol = vol + vol_xy*(1.+cyl*((xs(ip)+xs(ip+1))/3.-1.))
   40 continue
      volume = abs(vol)
c
99    continue
      return
      end

      subroutine prtplt2 (n)

c

c ======================================================================

c

c   Purpose -

c     plot and provide formatted writes to paper and film,

c     where    = 1: initial mesh and problem data

c            n = 2: time step and cycle info

c              = 3: field variables

c

c   PRTPLT is called by -

c

c          name      name      name      name      name      name

c        --------  --------  --------  --------  --------  --------

c          RIPPLE    NEWCYC

c

c

c   PRTPLT calls the following subroutines and functions -

c

c          name      source    name      source    name      source

c        --------  --------  --------  --------  --------  --------

c          GLOBAL    ripple

c

c ======================================================================

c

c##############################################################

       implicit real*8 (a-h,o-z)

c       include "32bit.h"

c##############################################################

c

c############

       include "comdk2.h"
        character*4 fname
        character*1 ntype
        data kounter/0/

c      call global

      go to (20,120,140), n


   20 continue


c.... write out mesh information

c
        open(2,file=fdir//'data_status.dat')
        open(3,file=fdir//'data_xi.dat')
        open(4,file=fdir//'data_yj.dat')
        open(23,file=fdir//'data_ar.dat')
        open(24,file=fdir//'data_at.dat')
        open(43,file=fdir//'surface.dat')


        write(2,*) imax,jmax,prtdt,1
      write (3,99) (x(i),i=1,imax)
      write (4,99) (y(j),j=1,jmax)

        do  j=1,jmax
        write(23,99) (ar(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
        end do

        do  j=1,jmax
        write(24,99) (at(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
        end do

 99     format(3x,3000(f8.3))
        close(2)
        close(3)
        close(4)
        close(23)
        close(24)

      go to 9999

c.... prtplt (2)  write time step, cycle information

c

  120 continue

c       write (13,330) iter,t,idt,itc,jtc,delt,ncyc,vchgt
c     write (13,260) fvol,vvol,xmv,ymv,tke

        write(*,331)iter,t,delt,ncyc

 331  format(3x,5hiter=,i6,3x,5htime=,e15.6,3x,5hdelt=,e12.6,3x

     1  ,6hcycle=,i9)

      go to 9999

c
c.... prtplt (3)  write field variables to paper

  140 continue

c       output F function and velocity vectors  Qun

        kounter=kounter+1

        n_first=mod(kounter/1000,10)
        n_second=mod(kounter/100,10)
        n_third=mod(kounter/10,10)
        n_fourth=mod(kounter,10)


        write(fname(1:1),'(I1)')n_first
        write(fname(2:2),'(I1)')n_second
        write(fname(3:3),'(I1)')n_third
        write(fname(4:4),'(I1)')n_fourth

        open(25,file=fdir//'data_f.'//fname)
        open(26,file=fdir//'data_u.'//fname)
        open(27,file=fdir//'data_v.'//fname)
        open(28,file=fdir//'data_p.'//fname)
c        open(29,file=fdir//'data_h.'//fname)
        open(30,file=fdir//'data_prod.'//fname)
        open(31,file=fdir//'data_k.'//fname)

        print*,'printing...',fname

        do j=1,jmax
        write(25,299) (f(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
        write(26,299) (u(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
        write(27,299) (v(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
        write(28,299) (p(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
        write(31,299) (sqrt(2.0*xk(ij)),ij=(j-1)*imax+1,(j-1)*imax+imax)
        write(30,299) (xprod(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
c        write(31,99) (ken(1,ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
        end do

299     format(3x,3000e14.6)

c        write(29,99) (etah(i),i=1,imax)

        close(25)
        close(26)
        close(27)
        close(28)
        close(29)
        close(30)
        close(31)

c       open(19,file='tmp.dat')
c       do j=1,jmax
c       write(19,99) (rkc(1,ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
c       end do
c       close(19)

244     format(I3)

  190 continue

 9999 return

  260 format (10x," fluid volume ",1pe14.6," void volume ",1pe14.6,/,

     &        10x," x momentum ",1pe14.6," y momentum ",1pe14.6,/,

     &        10x," total ke ",1pe14.6)

  270 format (1h1,1x,43hnf field (incl. fictitious cells) for cycle,i6,

     1         2x,2ht=,1pe14.7,2x,5hdelt=,e12.5//5x,32i4)

  275 format (1h1,1x,43hnw field (incl. fictitious cells) for cycle,i6,

     1         2x,2ht=,1pe14.7,2x,5hdelt=,e12.5//5x,32i4)

  280 format (1x,i3,1x,63i2/(5x,63i2))

  290 format (1h1)

  300 format (a80)

  310 format (4x,1hi,5x,1hj,9x,1hu,14x,1hv,15x,1hp,15x,1hd,12x,

     1         1hf,11x,2hnf)

  320 format (2x,i3,3x,i3,6(3x,1pe12.5),3x,i3)

  330 format (1x," iter = ",i5," time = ",1pe12.5," dt",a2,"(",

     &        i2,",",i2,") = ",1pe12.5," cycle = ",i6," vchgt = ",

     &        1pe12.5)

  350 format (1h0)

  360 format (1h ,18x,a80,1x,a10,2(1x,a8))

  370 format (2x,5hnkx= ,i4)

  380 format (2x,8hmesh-x= ,i4,3x,4hxl= ,1pe12.5,3x,4hxc= ,e12.5,3x,

     1       4hxr= ,e12.5,3x,5hnxl= ,i4,3x,5hnxr= ,i4,3x,6hdxmn= ,e12.5)

  390 format (2x,8hmesh-y= ,i4,3x,4hyl= ,1pe12.5,3x,4hyc= ,e12.5,3x,

     1       4hyr= ,e12.5,3x,5hnyl= ,i4,3x,5hnyr= ,i4,3x,6hdymn= ,e12.5)

  400 format (2x,5hnky= ,i4)

  430 format (13x,i3,2x,1pe12.5,3x,e12.5)

         return
         end
