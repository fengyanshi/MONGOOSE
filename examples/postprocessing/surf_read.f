         parameter(m=1001)
         real h(m),ngage(15),xg(15)

	 
	 iorgin=100

         open(1,file='./Results/surface.dat')
         open(2,file='ele2.dat')
         open(3,file='sur.dat')
         do k=1,999999999
         read(1,100,end=200)t,(h(i),i=1,m)
         write(2,300)t,h(200),h(491),h(337),h(254),h(193),h(131)
c         write(3,300)(h(i),i=1,m,4)
         enddo
200      continue
         print*,k,t
100      format(100f10.5)
300      format(500f10.5)

         end
