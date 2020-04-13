       parameter(nloc_max=100,jmax=80)
c note jmax is just grid number in y
       real f(jmax),u(jmax),v(jmax),p(jmax),xxx(jmax)
       character*1 fn

       print*,'gage number=5'
       num_g=5
       do 500 ng=1,num_g
c       do 500 ng=1,1
       write(fn(1:1),'(I1)')ng
       open(1,file='./Results/gage'//fn)
       open(2,file='p'//fn//'.dat')
       open(3,file='h'//fn//'.dat')
       open(4,file='u'//fn//'.dat')
       open(7,file='f'//fn//'.dat')
       do 100 kt=1,999999999
        read(1,200,end=300)t,etah,(f(kk),kk=1,jmax),
     &                     (u(kk),kk=1,jmax),
     &                     (v(kk),kk=1,jmax),
     &                     (p(kk),kk=1,jmax),
     &                     (xxx(kk),kk=1,jmax),
     &                     (xxx(kk),kk=1,jmax),
     &                     (xxx(kk),kk=1,jmax)
        write(2,200)t,(p(kk),kk=1,jmax)
	write(3,200)t,etah
	write(4,200)t,(u(kk),kk=1,jmax)
	write(7,200)t,(f(kk),kk=1,jmax)
100    continue       
300    continue
       print*,'gage ',ng,'time=',t
       close(1)
       close(2)
       close(3)
       close(4)
       close(7)
500    continue
200    format(100e14.6)
      
      end
