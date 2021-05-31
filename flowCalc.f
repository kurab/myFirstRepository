c
c    differential calculus Central Method
c
        parameter(nx=120)
        dimension f(nx),fn(nx)
        data xl,ni,dt,c/100.,101,0.1,0.5/
        data xp,cp/5.,0.5/
        data etime,tint/200.,5./
c
        open(30,file='ans.d',status='unknown')
        write(30,*) xl,ni
c
        dx=xl/float(ni-1)
        do i=1,ni
         x=dx*(i-1)
         if(x.lt.xp) then
          f(i)=cp/xp*x
         else if(x.lt.xp*2.) then
          f(i)=-cp/xp*x+2.*cp
         else
          f(i)=0.
         end if
        end do
c
        time=0.
        tt=tint
  100 continue
        if(tt.ge.tint) then
         tt=0.
         write(30,*) time
         do i=1,ni
          write(30,*) f(i)
         end do
        end if
c
        do i=2,ni-1
         fn(i)=f(i)-(f(i+1)-f(i-1))*c*dt/2./dx
         fn2(i)=f(i)-(f(i+1)-f(i-1))*c*dt/2./dx**2.0d0
        end do
c
        do i=1,ni
         f(i)=fn(i)
        end do
        time=time+dt
        tt=tt+dt
        if(time.gt.etime) goto 999
        goto 100
  999 continue
        end
