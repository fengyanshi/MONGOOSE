        parameter(max_nwave=4000)
        real aawave(max_nwave),twave(max_nwave),phase(max_nwave)

        open(10,file='spectral')
        read (10,*) nwave, (aawave(nw),nw=1,nwave),
     &       (twave(nw),nw=1,nwave), (phase(nw),nw=1,nwave)

        print*,twave
        open(2,file='series.dat')
      pi = 3.1415926
      dt=0.1
      total_t=1200.
      max_t=total_t/dt
      
         t=0.
         do nt=1,max_t
         t=t+dt
         ele=0.
         do nw=1,nwave
           ele=ele+aawave(nw)*sin(2.0*pi/twave(nw)*t
     &    +phase(nw))
         enddo
         write(2,*)t,ele
         enddo

	 end
