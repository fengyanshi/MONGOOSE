
      subroutine surface
c
c ======================================================================
c
c   Purpose - integrate surface elevation
c
c ======================================================================
c
c##############################################################
        implicit real*8 (a-h,o-z)
c        include "32bit.h"
c##############################################################
c
c############
        include "comdk2.h" 
c############
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
	if(ncyc.eq.0) then
	do i=1,imax
	etah(i)=flht
	end do
	goto 99
	end if
	do 200 i=1,imax
	do 100 j=jm1,2,-1
      ij=(j-1)*imax+i
      ijm=ij-imax
	if(f(ij).le.0.5.and.f(ijm).gt.0.5) then
	iy=j-1
	goto 150
	end if
100	continue
      ij2=imax+i
	if (f(ij2).lt.emf) etah(i)=0.
      ijm1=(jm1-1)*imax+i
	if (f(ijm1).ge.0.5) etah(i)=y(jm1)+0.5*dely(jm1)*f(ijm1)
	goto 200
150	continue
	ijiy=(iy-1)*imax+i
	ijiyp=ijiy+imax 
	yeta=dely(iy)*(f(ijiy)-0.5)+dely(iy+1)*f(ijiyp)
	etah(i)=yj(iy)+dmax1(0.d0,yeta)
200	continue
99	continue
c print
c        write(43,190)t,(etah(i),i=1,imax)
190     format(100f10.5)
	return
	end


        subroutine prtplt_surface
        include "comdk2.h"
	
	print*,imax,etah(1)
	
c         write(43,*)(etah(i),i=1,imax)
 100     format(100f10.5)
        
        end
