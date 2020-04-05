c----------------------------------------------
c  This subroutine is based on Buscaglia et al. 2002
c  bubble entrainment is presribed 
c  bubble breakup and coalescence are based on 
c  Prince and Blanch (1990), Carrica et al.(1999), Luo (1993)
c  Luo and Svendsen (1996)

        subroutine init_bubble
         implicit real*8 (a-h,o-z)
c############
       include "comdk2.h"
       include "bubble1.h"
       real*8 rlogR(mbub),rlogN(mbub),rN2(mbub),rlogN2(mbub)
c############
        do 10 k_bub=1,mbub
        do 10 i=2, im1
        do 10 j=2, jm1
        ij=(j-1)*imax+i
        cbub(k_bub,ij)=0.
        rnum_b(k_bub,ij)=0.
        bcprod(k_bub,ij)=0.
        bcdedu(k_bub,ij)=0.
        bnprod(k_bub,ij)=0.
        bndedu(k_bub,ij)=0.
        bnutt(k_bub,ij)=0.0005
10      continue

        do i=2, im1
        do j=2, jm1
        ij=(j-1)*imax+i
        alpha_g(ij)=0.
        enddo
        enddo

        Temp=20.+273.2

c --- Gaussian
        xcen=x(im1/2+1)
        ycen=y(jm1/10)
        aa0=0.00000005

c        do 20 k_bub=1,mbub
c        do 20 i=2, im1
c        do 20 j=2, jm1
c        ij=(j-1)*imax+i
c        if(f(ij).gt.emf.and.ac(ij).gt.emf) then
c        xr=x(i)-xcen
c        yr=y(j)-ycen
c        cbub(k_bub,ij)=aa0*exp(-20.0*(xr*xr+yr*yr))
c        endif
c20      continue

c ---   initial sizes of bubbles
        rlogR0=-1.
        rlogRn=1.
        do k_bub=1,mbub
          rlogR(k_bub)=rlogR0+(rlogRn-rlogR0)*(k_bub-1)/(mbub-1)
          rb(k_bub)=10**(rlogR(k_bub))*0.001
        enddo
c ---   bubble mid production numbers

        alpha_b=-3./2.
        beta_b=-10./3.
c   N0 = 2x10e4
        rlogN0=4.3
c  a_b = 1000. - a
        a_b=15000000.
        r_we=1.*0.001

        rlogN(1)=rlogN0
        do k_bub=2,mbub
         if(rb(k_bub).le.r_we)then
          rlogN(k_bub)=rlogN(k_bub-1)
     &       +alpha_b*(rlogR(k_bub)-rlogR(k_bub-1))
         else
         rlogN(k_bub)=rlogN(k_bub-1)
     &       +beta_b*(rlogR(k_bub)-rlogR(k_bub-1))
         endif
        enddo

        do k_bub=1,mbub
          rN2(k_bub)=10.**(rlogN(k_bub))
        enddo

        do k_bub=1,mbub-1
          prod_b_per_xprod_dt(k_bub)=
     &      a_b*rN2(k_bub)*(rb(k_bub+1)-rb(k_bub))
        enddo
          prod_b_per_xprod_dt(mbub)=
     &      a_b*rN2(mbub)*(rb(mbub)-rb(mbub-1))

c        print*,rlogR
c        print*,rb
c        print*,rN2
c        print*,prod_b_per_xprod_dt
c         stop

        call slip_velocity

c        call init_bubble_number

        return
        end

       subroutine bubble_all
       implicit real*8 (a-h,o-z)
c############
       include "comdk2.h"
       include "bubble1.h"

       call radius_update
       call c_b_transfer

c slip_velocity already initialized
c       call slip_velocity

c ---  call all bubble propagation and transformation
        do k_bub=1,mbub
        call bubble_tracking(k_bub)
        enddo

c --- calculate alpha_g

        call calculate_alpha_g

        return
        end

       subroutine slip_velocity
       implicit real*8 (a-h,o-z)
c############
       include "comdk2.h"
       include "bubble1.h"

c --- slip velocity
        do 20 k_bub=1,mbub
         wb(k_bub)=0.23
         if(rb(k_bub).le.0.0007)
     &      wb(k_bub)=4474*rb(k_bub)**(1.357)
         if(rb(k_bub).gt.0.0051)
     &      wb(k_bub)=4.202*rb(k_bub)**(0.547)
20      continue

        return
        end


       subroutine calculate_alpha_g
       implicit real*8 (a-h,o-z)
c############
       include "comdk2.h"
       include "bubble1.h"

c --- calculate alpha_g
        do 20 i=2, im1
        do 20 j=2, jm1
        ij=(j-1)*imax+i
        if(f(ij).gt.emf.and.ac(ij).gt.emf) then
        vn_total=0.
        do k_bub=1,mbub
        vn_total=vn_total+rnum_b(k_bub,ij)*4.*3.1415926/3.
     &           *rb(k_bub)**3
        enddo

        if(vn_total.ge.0.70)then
        alpha_g(ij)=0.70
        else
        alpha_g(ij)=vn_total
        endif
        
        else
        alpha_g(ij)=0.
        endif
20      continue

        return
        end


      subroutine radius_update
       implicit real*8 (a-h,o-z)
c############
       include "comdk2.h"
       include "bubble1.h"

c --- calculate radius
        do 20 ibub=1,mbub
        do 20 i=2, im1
        do 20 j=2, jm1
        ij=(j-1)*imax+i
        if(f(ij).gt.emf.and.ac(ij).gt.emf) then
        if(rnum_b(ibub,ij).ge.1)then
c --- p0=98000.
        vnu_k=Rgas*Temp*cbub(ibub,ij)/(p(ij)+9.8*rhof*10.)
     &        /rnum_b(ibub,ij)
        rb_update(ibub,ij)=rb(ibub)
        endif        
        endif
20      continue

       return
       end

      subroutine c_b_transfer
      implicit real*8 (a-h,o-z)
      include "comdk2.h"
      include "bubble1.h"

c --- we assume that T and P don't change much, so not 
c     inter-group transfer due to radius change

c --- initialize c and b transfer

      do 10 ibub=1,mbub
      do 10 i=2,im1
      do 10 j=2,jm1
      ij=(j-1)*imax+i
      cbub_transf(ibub,ij)=0.      
      rnum_transf(ibub,ij)=0.
10    continue      

c --- coalescence
      stension=0.0728

      do 20 i=2,im1
      do 20 j=2,jm1
      ij=(j-1)*imax+i
      if(f(ij).gt.emf.and.ac(ij).gt.emf) then

        do 30 ibub_k=1,mbub-1
        do 30 ibub_l=ibub_k,mbub
          v_kl=rb(ibub_k)**3+rb(ibub_l)**3
c --- coalescence rate
          ub_k2=2.*xnut(ij)**(2./3.)*(2.*rb(ibub_k))**(2./3.)
          ub_l2=2.*xnut(ij)**(2./3.)*(2.*rb(ibub_l))**(2./3.)
          ub_kl=sqrt(ub_k2+ub_l2)
          T_kl=3.14159/4.*(2.*rb(ibub_k)+2.*rb(ibub_l))**2.
     &         *ub_kl*rnum_b(ibub_k,ij)*rnum_b(ibub_l,ij)
          xi_kl=rb(ibub_k)/rb(ibub_l)
          tkl=0.5*rhof*ub_kl*(2.*rb(ibub_k))**2*(2.*rb(ibub_l))**2
     &        /(2.*rb(ibub_k)+2.*rb(ibub_l))**2/stension
          tau_kl=(1.+xi_kl)*sqrt((1./rhof+0.5)*rhof*(2.*rb(ibub_k))**3
     &            /3./(1.+xi_kl**2)/(1.+xi_kl**3)/stension)
          z_kl=exp(-tkl/tau_kl)
          theta_kl=z_kl*T_kl
c --- note that the formulus given by Prince and Blanch, 1990 may not be
c     good for N larger than 10^7. So I make a limit of 10% here
           if(theta_kl.gt.
     &       0.1*min(rnum_b(ibub_k,ij),rnum_b(ibub_l,ij)))then
              theta_kl=0.1*min(rnum_b(ibub_k,ij),rnum_b(ibub_l,ij))
           endif
          rnum_transf(ibub_k,ij)=rnum_transf(ibub_k,ij)-theta_kl
          rnum_transf(ibub_l,ij)=rnum_transf(ibub_l,ij)-theta_kl
          do ibub_i=2,mbub-1
             v_i=rb(ibub_i)**3
             v_i_n=rb(ibub_i-1)**3
             v_i_p=rb(ibub_i+1)**3
           if(v_kl.gt.v_i_n.and.v_kl.le.v_i)then
            X_ikl=(v_kl-v_i_n)/(v_i-v_i_n)
           rnum_transf(ibub_i,ij)=rnum_transf(ibub_i,ij)+theta_kl*X_ikl
           endif
           if(v_kl.gt.v_i.and.v_kl.le.v_i_p) then
            X_ikl=(v_i_p-v_kl)/(v_i_p-v_i)
           rnum_transf(ibub_i,ij)=rnum_transf(ibub_i,ij)+theta_kl*X_ikl
           endif
          enddo
          ibub_i=mbub
          v_i=rb(ibub_i)**3
          if(v_kl.gt.v_i)then
           rnum_transf(ibub_i,ij)=rnum_transf(ibub_i,ij)+theta_kl
          endif
30      continue
      endif
20    continue

      goto 4111
c --- breakup
      c_B=0.923
      beta_b=2.04
c --- here I adopted f_BV=0.5, Cf=f_BV^2/3 +(1-F_BV)^2/3-1
      Cf_b=0.2599
      do 50 i=2,im1
      do 50 j=2,jm1
      ij=(j-1)*imax+i
      if(f(ij).gt.emf.and.ac(ij).gt.emf) then
      do 60 ibub_i=1,mbub
c lambda could be calculated using van den Hengel et al 2005
c here I use lambda_m= 11.4 lambda_d, Lou and Svendsen 1996, p.1229 
c test       xnut(ij)=0.5
c test      alpha_g(ij)=0.1
c test      rnum_b(ibub_i,ij)=100000000.

       rlambda_m=11.4*(xnu**3/max(0.0001,xnut(ij)))**(1./4.)
       xi_m=rlambda_m/rb(ibub_i)
       if(xi_m.gt.1.)xi_m=1.
       rintegral=0.
       do kk=1,100
       xib=xi_m+(kk-1.)/99.*(-xi_m+1.)
       dxib=(1.-xi_m)/99.
       rintegral=rintegral+(1.+xib)**2/xib**(11./3.)*dxib*
     &            exp(-12.*cf_b*stension
     &           /beta_b/rhof*max(0.0001,xnut(ij))**(-2/3)
     &           *(2.*rb(ibub_i))**(-5./3.)*xi_m**(-11./3.))
       enddo
       b_k=c_B*(1.-alpha_g(ij))*rnum_b(ibub_i,ij)*(xnut(ij)
     &            /4./rb(ibub_i)**2)**(1./3.)*rintegral
       rnum_transf(ibub_i,ij)=rnum_transf(ibub_i,ij)-b_k
       v_i_2=(rb(ibub_i)/2.)**3
       do ibub_k=2,ibub_i-1
         v_k_p=rb(ibub_k+1)**3
         v_k_n=rb(ibub_k-1)**3
         v_k=rb(ibub_k)**3
         if(v_i_2.gt.v_k_n.and.v_i_2.le.v_k)then
         X_ik=2.*(v_i_2-v_k_n)/(v_k-v_k_n)
         rnum_transf(ibub_k,ij)=rnum_transf(ibub_k,ij)+b_k*X_ik
         endif
         if(v_i_2.gt.v_k.and.v_i_2.le.v_k_p)then
         X_ik=2.*(v_k_p-v_i_2)/(v_k_p-v_k)
         rnum_transf(ibub_k,ij)=rnum_transf(ibub_k,ij)+b_k*X_ik
         endif         
       enddo
         v_k_n=rb(1)**3
         v_k_p=rb(ibub_i)**3
         if(v_i_2.le.v_k_n)then
         rnum_transf(ibub_k,ij)=rnum_transf(ibub_k,ij)+b_k*2.
         endif
         if(v_i_2.ge.v_k_p)then
         rnum_transf(ibub_k,ij)=rnum_transf(ibub_k,ij)+b_k*2.
         endif
60    continue
      endif
50    continue

4111  continue

c --- production
      do 33 ibub=1,mbub
      do 33 i=2,im1
      do 33 j=2,jm1
      ij=(j-1)*imax+i
      if(f(ij).gt.emf.and.ac(ij).gt.emf) then
       bcprod(ibub,ij)=cbub_transf(ibub,ij)
       bnprod(ibub,ij)=rnum_transf(ibub,ij)
      else
c       bcprod(ibub,ij)=0.
       bnprod(ibub,ij)=0.
      endif
33    continue
      end

      subroutine init_bubble_number
       implicit real*8 (a-h,o-z)
c############
       include "comdk2.h"
       include "bubble1.h"

c --- calculate numbers
        do 20 ibub=1,mbub
        do 20 i=2, im1
        do 20 j=2, jm1
        ij=(j-1)*imax+i
        if(f(ij).gt.emf.and.ac(ij).gt.emf) then
c -- p0=98000
        alpha_k=Rgas*Temp*cbub(ibub,ij)/(p(ij)+100.)

        rnum_b(ibub,ij)=alpha_k*3./4./pi/(rb(ibub))**3
        endif
20      continue

       return
       end

       subroutine buoyancy
       implicit real*8 (a-h,o-z)
c############
       include "comdk2.h"
       include "bubble1.h"
c############

c buoyancy at v point

        do 20 i=2, im1
        do 20 j=2, jm1
        ij=(j-1)*imax+i
        ijp=ij+imax
        ijm=ij-imax
        if(f(ij).gt.emf.and.ac(ij).gt.emf) then
        buoy(ij)=(alpha_g(ijp)*dely(j-1)+alpha_g(ij)*dely(j))
     &           /(dely(j-1)+dely(j))
        if(buoy(ij).gt.0.5)then
        print*,'problem with buoyancy ...',ij,buoy(ij)
        stop
        endif
        endif
20      continue


         return
         end

         subroutine bubble_tracking(m)
         implicit real*8 (a-h,o-z)
c############
       include "comdk2.h"
       include "bubble1.h"
c############
        do 10 i=2, im1
        do 10 j=2, jm1

          ij=(j-1)*imax+i
          ipj=ij+1
          imj=ij-1
          ijp=ij+imax
          ijm=ij-imax
          ipjp=ipj+imax
          ipjm=ipj-imax
          imjp=imj+imax
          imjm=imj-imax

        uij=0.
        vij=0.

        if(f(ij).gt.emf.and.ac(ij).gt.emf) then
        call dccal(m,ij)
        call rncal(m,ij)

        uij=0.5*(un(ij)+un(imj))
        if(uij.ge.0) fcx=dcl*uij
        if(uij.lt.0) fcx=dcr*uij
        if(uij.ge.0) fnx=rncl*uij
        if(uij.lt.0) fnx=rncr*uij

        vij=0.5*(vn(ij)+vn(ijm))+wb(m)
        if(vij.ge.0) fcy=dcb*vij
        if(vij.lt.0) fcy=dct*vij
        if(vij.ge.0) fny=rncb*vij
        if(vij.lt.0) fny=rnct*vij


        xnuttr=(bnutt(m,ipj)*delx(i)+bnutt(m,ij)*delx(i+1))
     &      /(delx(i)+delx(i+1))

        xnuttl=(bnutt(m,imj)*delx(i)+bnutt(m,ij)*delx(i-1))
     &      /(delx(i)+delx(i-1))
        viskx=(xnuttr*(cbub(m,ipj)-cbub(m,ij))/(xi(i+1)-xi(i))
     &       -xnuttl*(cbub(m,ij)-cbub(m,imj))/(xi(i)-xi(i-1)))/delx(i)
        vnx=(xnuttr*(rnum_b(m,ipj)-rnum_b(m,ij))/(xi(i+1)-xi(i))
     &    -xnuttl*(rnum_b(m,ij)-rnum_b(m,imj))/(xi(i)-xi(i-1)))/delx(i)

        xnuttt=(bnutt(m,ij)  *dely(j+1)+bnutt(m,ijp)*dely(j))
     &     /(dely(j)+dely(j+1))

        xnuttb=(bnutt(m,ij)  *dely(j-1)+bnutt(m,ijm)*dely(j))
     &     /(dely(j)+dely(j-1))

        visky=(xnuttt*(cbub(m,ijp)-cbub(m,ij))/(yj(j+1)-yj(j))
     &       -xnuttr*(cbub(m,ij)-cbub(m,ijm))/(yj(j)-yj(j-1)))/dely(j)

        vny=(xnuttt*(rnum_b(m,ijp)-rnum_b(m,ij))/(yj(j+1)-yj(j))
     &    -xnuttr*(rnum_b(m,ij)-rnum_b(m,ijm))/(yj(j)-yj(j-1)))/dely(j)

        cbub(m,ij)=cbub(m,ij)
     1         +delt*(-fcx-fcy+viskx+visky)

        rnum_b(m,ij)=rnum_b(m,ij)
     1         +delt*(-fnx-fny+vnx+vny)+bnprod(m,ij)

        end if
 10     continue
c

c       boundary condition for bubble
c
        do 200 i=2,im1
        do 200 j=2,jm1

           ij=(j-1)*imax+i
          ipj=ij+1
          imj=ij-1
          ijp=ij+imax
          ijm=ij-imax

c  -- no bubbles in air cells except boundary points
       if(f(ij).lt.emf.or.ac(ij).lt.emf) then
        cbub(m,ij)=0.d0
        rnum_b(m,ij)=0.
        end if

c
c       free surface bubbles one-way transfer
c


        nff=nf(ij)
        nfr=nf(ipj)
        nfl=nf(imj)
        nft=nf(ijp)
        nfb=nf(ijm)

        if(nff.ge.5) then
        if(nft.lt.5) then
         cbub(m,ij)=cbub(m,ijp)
         rnum_b(m,ij)=rnum_b(m,ijp)
        endif
        if(nfb.lt.5) then
         cbub(m,ij)=cbub(m,ijm)
         rnum_b(m,ij)=rnum_b(m,ijm)
        endif
        if(nfr.lt.5) then
         cbub(m,ij)=cbub(m,ipj)
         rnum_b(m,ij)=rnum_b(m,ipj)
        endif
        if(nfl.lt.5) then
         cbub(m,ij)=cbub(m,imj)
         rnum_b(m,ij)=rnum_b(m,imj)
        endif
          end if

c
c       wall boundary
c
 200    continue

 999    continue
        return
        end

        subroutine dccal(m,ij)
c##############################################################
       implicit real*8 (a-h,o-z)
       include "comdk2.h"
       include "bubble1.h"

          ipj=ij+1
          imj=ij-1
          ijp=ij+imax
          ijm=ij-imax
          ipjp=ipj+imax
          ipjm=ipj-imax
          imjp=imj+imax
          imjm=imj-imax

        dcl=(cbub(m,ij)-cbub(m,imj))/(xi(i)-xi(i-1))
        dcr=(cbub(m,ipj)-cbub(m,ij))/(xi(i+1)-xi(i))
        dcb=(cbub(m,ij)-cbub(m,ijm))/(yj(j)-yj(j-1))
        dct=(cbub(m,ijp)-cbub(m,ij))/(yj(j+1)-yj(j))

        if(ar(ipj).lt.em6) dcr=0.
        if(ar(imj).lt.em6) dcl=0.
        if(at(ijp).lt.em6)
     &    dct=0.
        if(at(ijm).lt.em6)
     &    dcb=0.

        return
        end

        subroutine rncal(m,ij)
c##############################################################
       implicit real*8 (a-h,o-z)
       include "comdk2.h"
       include "bubble1.h"

          ipj=ij+1
          imj=ij-1
          ijp=ij+imax
          ijm=ij-imax
          ipjp=ipj+imax
          ipjm=ipj-imax
          imjp=imj+imax
          imjm=imj-imax

        rncl=(rnum_b(m,ij)-rnum_b(m,imj))/(xi(i)-xi(i-1))
        rncr=(rnum_b(m,ipj)-rnum_b(m,ij))/(xi(i+1)-xi(i))
        rncb=(rnum_b(m,ij)-rnum_b(m,ijm))/(yj(j)-yj(j-1))
        rnct=(rnum_b(m,ijp)-rnum_b(m,ij))/(yj(j+1)-yj(j))

        if(ar(ipj).lt.em6) rncr=0.
        if(ar(imj).lt.em6) rncl=0.
        if(at(ijp).lt.em6)
     &    rnct=0.
        if(at(ijm).lt.em6)
     &    rncb=0.

        return
        end


      subroutine prtplt_any
        implicit real*8 (a-h,o-z)
        include "comdk2.h"
        include "bubble1.h"
        character*4 fname
        character*2 fsize
        character*1 ntype
        data kcount/0/

        kcount=kcount+1

        n_first=mod(kcount/1000,10)
        n_second=mod(kcount/100,10)
        n_third=mod(kcount/10,10)
        n_fourth=mod(kcount,10)


        write(fname(1:1),'(I1)')n_first
        write(fname(2:2),'(I1)')n_second
        write(fname(3:3),'(I1)')n_third
        write(fname(4:4),'(I1)')n_fourth


        open(32,file=fdir//'data_ag.'//fname)
        do j=1,jmax
        write(32,999) (alpha_g(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
        end do
        close(32)

        open(32,file=fdir//'data_pr.'//fname)
        do j=1,jmax
        write(32,99) (xprod(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
        end do
        close(32)

 99     format(3x,3000(E16.3))
999     format(3x,3000(F16.8))

        return
        end

 
       subroutine prtplt_bubble(m)
        implicit real*8 (a-h,o-z)
        include "comdk2.h"
        include "bubble1.h"
        character*4 fname
        character*2 fsize
        character*1 ntype
        integer kcounter(mbub)
c        data kcounter(1:mbub)/mbub*0/

        kcounter(m)=kcounter(m)+1

        n_first=mod(kcounter(m)/1000,10)
        n_second=mod(kcounter(m)/100,10)
        n_third=mod(kcounter(m)/10,10)
        n_fourth=mod(kcounter(m),10)


        write(fname(1:1),'(I1)')n_first
        write(fname(2:2),'(I1)')n_second
        write(fname(3:3),'(I1)')n_third
        write(fname(4:4),'(I1)')n_fourth

        n_first=mod(m/10,10)
        n_second=mod(m,10)
        write(fsize(1:1),'(I1)')n_first
        write(fsize(2:2),'(I1)')n_second

        open(32,file=fdir//'data_c'//fsize//'.'//fname)

        print*,'printing bubble ...',fsize, '  at  ', fname

        do j=1,jmax
        write(32,99) (cbub(m,ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
        end do
        close(32)

     
c        open(32,file=fdir//'data_ag'//fsize//'.'//fname)

c        print*,'printing bubble ...',fsize, '  at  ', fname

c        do j=1,jmax
c        write(32,99)(alpha_g(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
c        end do
c        close(32)


        open(32,file=fdir//'data_nb'//fsize//'.'//fname)

        print*,'printing bubble ...',fsize, '  at  ', fname

        do j=1,jmax
        write(32,99) (rnum_b(m,ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
        end do
        close(32)

        open(32,file=fdir//'data_ntransf'//fsize//'.'//fname)

        print*,'printing bubble ...',fsize, '  at  ', fname

        do j=1,jmax
        write(32,99) (rnum_transf(m,ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
        end do
        close(32)


 99     format(3x,3000(E16.3))

        return
        end

       subroutine bubble_source
         implicit real*8 (a-h,o-z)
c############
       include "comdk2.h"
       include "bubble1.h"
       real*8 rnum_incr(mbub,nxy)
c############
c sand source

c gaussian

        xcen=x(im1/2+1)
        ycen=y(jm1/10)

        do 30 k_bub=1,mbub
c from 3 to im1-1 consider the thickness of surface layer, surface shear 
c equal to zero in the original code -fyshi
        do 30 i=3, im1-1
        do 30 j=3, jm1-1
        ij=(j-1)*imax+i

        if(f(ij).gt.emf) then

c generation box
         do ii=i-1,i+1
         do jj=j-1,j+1
          iijj=(jj-1)*imax+ii
          imj=iijj-1
          ipj=iijj+1
          ijp=iijj+imax
          ijm=iijj-imax
          ipjp=ipj+imax
          ipjm=ipj-imax
          imjp=imj+imax
          imjm=imj-imax
          if(nf(iijj).ne.0.or.
     &       nf(imj).ne.0.or.
     &       nf(ipj).ne.0.or.
     &       nf(ijp).ne.0.or.
     &       nf(ijm).ne.0.or.
     &       nf(ipjp).ne.0.or.
     &       nf(ipjm).ne.0.or.
     &       nf(imjp).ne.0.or.
     &       nf(imjm).ne.0)then
           isurf_layer=1
           goto 2323
          endif
         enddo
         enddo
         isurf_layer=0
2323     continue
         if(isurf_layer.eq.1)then
c ---   prod_min 0.01, max 1. log(prod) from 0 to 1 ~ 2, 100~ 4
          prod_min=max(0.001,xprod(ij))
c there is a bug pointed by Ma that log should be log10
          prod_log=2.+log10(prod_min)          
          prod_per_s=prod_b_per_xprod_dt(k_bub)*max(0.,prod_log)
          rnum_incr(k_bub,ij)=delt*prod_per_s

c          ssr=aa0*exp(-20.0*(xr*xr+yr*yr))
c          rnum_incr(k_bub,ij)=delt*ssr
c        if(j.eq.2.or.
c     &     (ac(ij).lt.1.0.and.at(ij).ge.1.))then
c        rnum_incr(k_bub,ij)=0.
c        rnum_incr(k_bub,ijp)=0.
c        endif
 
c       if(ac(ij).lt.1.0)rnum_incr(k_bub,ij)=0.

        endif
        endif
30      continue

c ---   calculate rnum and c

        do 40 k_bub=1,mbub
        do 40 i=2,im1
        do 40 j=2,jm1
        ij=(j-1)*imax+i
        if(f(ij).gt.emf) then
c p = p(ij)+rho*g*10.(m)
        rnum_b(k_bub,ij)=rnum_b(k_bub,ij)+rnum_incr(k_bub,ij)
        cbub(k_bub,ij)=cbub(k_bub,ij)+
     &     4.*3.1415926/3.*rnum_incr(k_bub,ij)*(p(ij)+rhof*9.8*10.)
     &    /Rgas/Temp*rb(k_bub)**3
        else
        rnum_b(k_bub,ij)=0.
        cbub(k_bub,ij)=0.
        endif
40      continue

       return
       end

