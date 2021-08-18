
!********************************************************************************
    subroutine fnset(Nparm,Numgr,Pttbl,Iptb,Indm,Param,Ipt,Indfn,     &
        & Icntyp,Confun)
!
implicit none

real*8 a , b , c , Confun , d , eps , one , p , Param , Pttbl ,   &
& q , x , y , zero
integer Icntyp , ii , Indfn , Indm , Ipt , Iptb , irem , Nparm ,  &
& Numgr
!
dimension Pttbl(Iptb,Indm) , Param(Nparm) , Icntyp(Numgr) ,       &
 & Confun(Numgr,Nparm+1)
!
! THIS IS THE SUBROUTINE FNSET FOR THE EXAMPLE DISCUSSED IN THE CONMAX
! USERS GUIDE.
!
! SET PRECISION DEPENDENT CONSTANTS.
one = 1.0d0
zero = one - one
!
! WE BREAK FNSET INTO SECTIONS BASED ON THE VALUE OF IPT, THAT IS, ON
! WHICH CONSTRAINT IS BEING SET.
if ( Ipt<=5 ) then
!
! HERE IPT .LE. 5 AND WE SET A CONSTRAINT OF THE FORM
! ABS(F(X,Y) - (AX+BY+C)/(DX+1)) .LE. W.
! NOTE THAT SINCE THIS IS A TYPE 2 CONSTRAINT WE DO NOT NEED TO DEAL
! WITH THE ABSOLUTE VALUE OR THE F(X,Y) HERE.
Icntyp(Ipt) = 2
a = Param(1)
b = Param(2)
c = Param(3)
d = Param(4)
x = Pttbl(Ipt,1)
y = Pttbl(Ipt,2)
p = a*x + b*y + c
q = d*x + one
Confun(Ipt,1) = p/q
if ( Indfn>0 ) then
!
! HERE IPT .LE. 5 AND INDFN=1, AND WE SET THE PARTIAL DERIVATIVES.
Confun(Ipt,2) = x/q
Confun(Ipt,3) = y/q
Confun(Ipt,4) = one/q
Confun(Ipt,5) = -p*x/(q*q)
return
endif
elseif ( Ipt<=15 ) then
!
! HERE 6 .LE. IPT .LE. 15 AND IF IPT IS EVEN WE SET THE CONSTRAINT
! (AX+BY+C)/(DX+1) .LE. W, WHICH IS HALF OF THE CONSTRAINT
! ABS((AX+BY+C)/(DX+1)) .LE. W, WHILE IF IPT IS ODD WE SET THE CONSTRAINT
! -(AX+BY+C)/(DX+1) .LE. W, WHICH IS THE OTHER HALF OF THE CONSTRAINT
! ABS((AX+BY+C)/(DX+1)) .LE. W.
Icntyp(Ipt) = 1
ii = (Ipt-4)/2
a = Param(1)
b = Param(2)
c = Param(3)
d = Param(4)
x = Pttbl(ii,1)
y = Pttbl(ii,2)
p = a*x + b*y + c
q = d*x + one
irem = Ipt - 4 - 2*ii
if ( irem<=0 ) then
!
! HERE 6 .LE. IPT .LE. 15 AND IPT IS EVEN.
Confun(Ipt,1) = p/q
if ( Indfn>0 ) then
  Confun(Ipt,2) = x/q
  Confun(Ipt,3) = y/q
  Confun(Ipt,4) = one/q
  Confun(Ipt,5) = -p*x/(q*q)
  return
endif
else
!
! HERE 6 .LE. IPT .LE. 15 AND IPT IS ODD.
Confun(Ipt,1) = -p/q
if ( Indfn>0 ) then
  Confun(Ipt,2) = -x/q
  Confun(Ipt,3) = -y/q
  Confun(Ipt,4) = -one/q
  Confun(Ipt,5) = p*x/(q*q)
  return
endif
endif
elseif ( Ipt<=20 ) then
!
! HERE 16 .LE. IPT .LE. 20 AND WE SET A CONSTRAINT OF THE FORM
! -DX - 1.0 + EPS .LE. 0.0
Icntyp(Ipt) = -1
d = Param(4)
eps = Pttbl(6,1)
ii = Ipt - 15
x = Pttbl(ii,1)
Confun(Ipt,1) = -d*x - one + eps
if ( Indfn>0 ) then
Confun(Ipt,2) = zero
Confun(Ipt,3) = zero
Confun(Ipt,4) = zero
Confun(Ipt,5) = -x
return
endif
else
!
! HERE IPT=21 AND WE SET THE CONSTRAINT
! (PARTIAL DERIVATIVE OF (AX+BY+C)/(DX+1) WITH RESPECT TO X AT
! (X,Y) = (0.0,0.0)) .LE. 0.0,
! I.E. A - CD .LE. 0.0
Icntyp(Ipt) = -2
a = Param(1)
c = Param(3)
d = Param(4)
Confun(Ipt,1) = a - c*d
if ( Indfn>0 ) then
Confun(Ipt,2) = one
Confun(Ipt,3) = zero
Confun(Ipt,4) = -d
Confun(Ipt,5) = -c
return
endif
endif

end subroutine fnset