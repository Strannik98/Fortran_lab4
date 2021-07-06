 program main

implicit none

integer n,k,ex

real*8 a,b,gaus5,trap,simp,lop,count(100000000),h, vib

ex=9

do while(ex.ne.0)

print *, '1-trap, 2-simp, 3-gaus, 0-exit'

read *, vib

if (vib.eq.0) then

goto 21

end if

print *, '------------------------------------------------------'

print *, 'cl-vo otrezkov'

read *, k

print *, 'start, end'

read *, a,b

print *, 'degree of polynomial'

read *, n

print *,'-------------------------------------------------------'

open (21, file='D:\Ih.txt')

do while(ex.ne.0)

h=(b-a)/k

call tabl(a,h,k,count(1))

if (vib.eq.1) then

lop = trap(n,k,h,count(1))

end if

if (vib.eq.2) then

lop = simp(n,k,h,count(1))

end if

if (vib.eq.3) then

lop = gaus5(n,k,h,count(1))

end if

print *, h, lop

write (21,*)lop

k=k*2

read *, ex

end do

print *, '-----------------------------------------------------'

print *, 'next-1, exit-0'

print *, '-----------------------------------------------------'

read *, ex

end do

21 end

real*8 function f(x,n)

implicit none

real*8 x,au,i

integer flag, n

flag = 1

c flag = 0

if(flag.ne.0) then

au = (n+1)*x**n

else

au = -10*sin(10*x)

end if

f=au

return

end

real*8 function gaus5(n,k,h,cont)

implicit none

integer n,k,i

real*8 f,au,h,q1,q2,q3,x2,x3,cont(10000)

au =0

q1=128d-0/225d-0

q2=(322d-0+13d-0*dsqrt(70d-0))/900d-0

q3=(322d-0-13d-0*dsqrt(70d-0))/900d-0

x2=1d-0/3d-0*dsqrt(5d-0-2d-0*dsqrt(10d-0/7d-0))

x3=1d-0/3d-0*dsqrt(5d-0+2d-0*dsqrt(10d-0/7d-0))

au=0

do i=1,k

au=au+h/2d-0*(q1*f((cont(i)+cont(i+1))/2d-0,n)+

*q2*f(((cont(i)+cont(i+1))/2d-0+x2*h/2d-0),n)+

*q2*f(((cont(i)+cont(i+1))/2d-0-x2*h/2d-0),n)+

*q3*f(((cont(i)+cont(i+1))/2d-0+x3*h/2d-0),n)+

*q3*f(((cont(i)+cont(i+1))/2d-0-x3*h/2d-0),n))

!print *, f(a+(i-1)*h,n), f(a+i*h,n)

end do

gaus5 = au

! print *, au

return

end

real*8 function trap(n,k,h,cont)

implicit none

integer n,k,i

real*8 h,f,au,cont(1000000000)

au =0

! print *,h

do i=1,k

au=au+h*(f(cont(i),n)+f(cont(i+1),n))/2d-0

!print *, f(a+(i-1)*h,n), f(a+i*h,n)

end do

trap = au

return

end

real*8 function simp(n,k,h,cont)

implicit none

integer n,k,i

real*8 h,f,au,cont(1000000000)

au =0

do i=1,k

au=au+h*(f(cont(i),n)+f(cont(i+1),n)+

*4d-0*f((cont(i)+cont(i+1))/2d-0,n))/6d-0

end do

simp = au

return

end

subroutine tabl(a,h,k, cont)

implicit none

integer i,k

real*8 a,h,cont(1000000000)

do i=1,k+1

cont(i)=a+h*(i-1)

! print *, h, i, cont(i)

end do

end
