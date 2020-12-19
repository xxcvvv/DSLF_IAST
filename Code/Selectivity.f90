PROGRAM IAST 
!******************************************************************************
! IAST v2.10 is used to predict adsorption loading of binary mixture from single 
! component isotherm by IAST-DSLF model.
! This program was last updated on 12/17/2019 
!******************************************************************************

IMPLICIT NONE
REAL(8),DIMENSION(2) :: y
REAL(8),DIMENSION(2) ::a,b,c,r,s,t
REAL(8) bisection
INTEGER i, j, n,Ps,Pe,In
REAL(8),DIMENSION(100) :: P1, P2, x1, x2, N1, N2,Se
REAL(8),DIMENSION(100) ::P

!Name of the input file is 'input.txt'
!Format of file is as followed:
!pressure ratio(use decimals):y1 y2
!the start and end of p,interval
!q1 q1
!b1 b2
!n1 n2
!q1 q2(another)
!b1 b2
!n1 n2

OPEN (UNIT=10,FILE="input.txt")
READ(10,*) (y(i),i=1,2)
READ(10,*) Ps,Pe,In
READ(10,*) (a(i),i=1,2)
READ(10,*) (b(i),i=1,2)
READ(10,*) (c(i),i=1,2)
READ(10,*) (r(i),i=1,2)
READ(10,*) (s(i),i=1,2)
READ(10,*) (t(i),i=1,2)
CLOSE(10)

n=(Pe-Ps)/In+1
DO i=1,n
    P(i)=Ps+(i-1)*In
END DO
DO i=1,n
   x1(i)=bisection(P(i),a,b,c,r,s,t,y)
   x2(i)=1-x1(i)
END DO

DO i=1,n
   P1(i)=P(i)*y(1)/x1(i)
   P2(i)=P(i)*y(2)/x2(i)
   N1(i)=a(1)*b(1)*P1(i)**c(1)/(1+b(1)*P1(i)**c(1))+a(2)*b(2)*P1(i)**(1/c(2))/(1+b(2)*P1(i)**(1/c(2)))
   N2(i)=r(1)*s(1)*P2(i)**t(1)/(1+s(1)*P2(i)**t(1))+r(2)*s(2)*P2(i)**(1/t(2))/(1+s(2)*P2(2)**(1/t(2)))
   N1(i)=N1(i)*N2(i)/(N2(i)+N1(i)*(1/x1(i)-1))
   N2(i)=N1(i)*(1/x1(i)-1)
   Se(i)=(x1(i)*y(2))/(x2(i)*y(1))
END DO

!****************************************************************
!Name of output file is'output.txt'
!Meaning:pressure x1 x2 selectivity
!****************************************************************
OPEN (UNIT=20,FILE="output.txt", STATUS="UNKNOWN")
DO i=1,n
   WRITE(20,100) P(i),Se(i),x1(i),x2(i)
END DO
100  format(5f16.3)
CLOSE(20)
END

!******************************************************************************
! Funtion bisection is used to calculte the root of equation by the method 
! of bisection
!******************************************************************************

FUNCTION bisection(P,a,b,c,r,s,t,y)
IMPLICIT NONE
REAL(8),DIMENSION(2) :: a, b, c, y,r,s,t
REAL(8) lt, rt, mi, f, P, bisection
lt=0.0
rt=1.0
DO WHILE (ABS(rt-lt)>1E-6)
   mi=(lt+rt)/2
   IF (f(P,a,b,c,r,s,t,y,lt)*f(P,a,b,c,r,s,t,y,mi)<0) THEN
      rt=mi
   ELSE
      lt=mi
   END IF
END DO
bisection=mi
END FUNCTION bisection 
!******************************************************************************

!******************************************************************************
! Funtion is the equation from IAST theory
!******************************************************************************
FUNCTION f(P,a,b,c,r,s,t,y,x1)
IMPLICIT NONE
REAL(8),DIMENSION(2) ::y
REAL(8),DIMENSION(2) ::a,b,c,r,s,t
REAL(8) P, x1, f,f1,f2
f1=a(1)*c(1)*LOG(1+b(1)*(P*y(1)/x1)**(1/c(1)))+a(2)*c(2)*LOG(1+b(2)*(P*y(1)/x1)**(1/c(2)))
f2=r(1)*t(1)*LOG(1+s(1)*(P*y(2)/(1-x1))**(1/t(1)))+r(2)*t(2)*LOG(1+s(2)*(P*y(2)/(1-x1)**(1/t(2))))
f=f1-f2
END FUNCTION f
!******************************************************************************