PROGRAM Qst
!*********************************************************
!Qst v2.0
!This program is used to calculate adsoprtion heat by 
!pressure at three different temperature
!*********************************************************
implicit none
real(8),dimension(10000) ::q,p1,p2,p3,y1,y2,y3
real(8),dimension(3) ::q1,q2,b1,b2,n1,n2,T,x
real(8) bisection,b,m,r,qs,qe,it
integer n,i,j
real(8) ::e=3,sumx=0,sumx2=0,sumxy=0,sumy=0,sumy2=0
!*********************************************************
!'input.txt'means:
!T(1) T(2) T(3)
!the start and end of q,then the interval
!(at 3 temperature)q1 q2
!b1 b2
!n1 n2
!...........
!********************************************************
open(10,file='input.txt')
read(10,*) (T(i),i=1,3)
read(10,*) qs,qe,it
read(10,*) q1(1),q2(1)
read(10,*) b1(1),b2(1)
read(10,*) n1(1),n2(1)
read(10,*) q1(2),q2(2)
read(10,*) b1(2),b2(2)
read(10,*) n1(2),n2(2)
read(10,*) q1(3),q2(3)
read(10,*) b1(3),b2(3)
read(10,*) n1(3),n2(3)
close(10)

n=(qe-qs)/it+1
q(1)=qs
do j=2,n
q(j)=q(j-1)+it
end do
do i=1,n
p1(i)=bisection(q1(1),q2(1),b1(1),b2(1),n1(1),n2(1),q(i))
p2(i)=bisection(q1(2),q2(2),b1(2),b2(2),n1(2),n2(2),q(i))
p3(i)=bisection(q1(3),q2(3),b1(3),b2(3),n1(3),n2(3),q(i))
y1(i)=log(p1(i))
y2(i)=log(p2(i))
y3(i)=log(p3(i))
end do
do i=1,3
x(i)=1/T(i)
end do

!**************************************************************
!'output1.txt'means:
!0000  T1   T2   T3
!q(1) p(1) p(2) p(3)
!q(2) p(1) p(2) p(3)
!..................
!q(n) p(1) p(2) p(3)
!************************
!'output2,txt'means:
!q Qst intercept r
!..................
!**************************************************************
do i=1,n
sumx=x(1)+x(2)+x(3)
sumx2=x(1)**2+x(2)**2+x(3)**2
sumxy=x(1)*y1(i)+x(2)*y2(i)+x(3)*y3(i)
sumy=y1(i)+y2(i)+y3(i)
sumy2=y1(i)**2+y2(i)**2+y3(i)**2
m=(e*sumxy-sumx*sumy)/(e*sumx2-sumx**2)                         
b=(sumy*sumx2-sumx*sumxy)/(e*sumx2-sumx**2)                    
r=(sumxy-sumx*sumy/e)/sqrt((sumx2-sumx**2/e)*(sumy2 -sumy**2/e))
open(30,file='output2.txt',status='unknown')
write(30,100) q(i),m*(-8.314),b,r*(-100)
end do
open(20,file='output1.txt',status='unknown')
write(20,100) n,T(1),T(2),T(3)
do i=1,n 
write(20,100) q(i),p1(i),p2(i),p3(i)
end do
100 format(5f16.4)
end
!***************************************************************

!***************************************************************
!this function is used to calculate the p specicfic to the 
!given q by bisection method
!***************************************************************
function bisection(q1,q2,b1,b2,n1,n2,q)
implicit none
real(8) q1,q2,b1,b2,n1,n2,q,lt,rt,mi,f,bisection
lt=0
rt=1500
do while (abs(rt-lt)>1E-3)
    mi=(lt+rt)/2
    if(f(lt,q1,q2,b1,b2,n1,n2,q)*f(mi,q1,q2,b1,b2,n1,n2,q)<0)then
        rt=mi
    else
        lt=mi
    end if 
end do
bisection=mi
end function bisection
!********************************************************************

!********************************************************************
function f(p,q1,q2,b1,b2,n1,n2,q)
implicit none
real(8) q1,q2,b1,b2,n1,n2
real(8) q
real(8) f,p
f=q1*b1*p**(1/n1)/(1+b1*p**(1/n1))+q2*b2*p**(1/n2)/(1+b2*p**(1/n2))-q
end function f
!*******************************************************************
