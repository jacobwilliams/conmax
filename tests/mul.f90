
! (A) MULLERS METHOD DERIVATIVE FREE REAL ROOT FINDING
!
! (SUBPROGRAMS INVOLVED:  MULLER, FNSET (USER SUPPLIED), ILOC, D1MACH,
! ERCMP1)
!
! GIVEN A FUNCTION F OF ONE VARIABLE (WHERE F(X) IS COMPUTED IN SUBROUTINE
! FNSET AS CONFUN(1,1), WITH X = PARAM(1)), A NONNEGATIVE TOLERANCE TOLCON,
! TWO POINTS (P1,F1) AND (PROCOR,EMIN) WITH F1 = F(P1), EMIN = F(PROCOR),
! P1 .LT. PROCOR, F1 .GT. TOLCON, AND EMIN .LT. -TOLCON, THE PROGRAM WILL
! ATTEMPT TO LOCATE A NEW PROCOR WITH NEW EMIN = F(PROCOR) AND ABS(EMIN)
! .LE. TOLCON.  IF IT FAILS TO DO THIS, THE PROGRAM WILL RETURN WITH
! PROCOR = THE LEFTMOST ABSCISSA FOUND WITH EMIN = F(PROCOR) .LT. -TOLCON.
! NOTE:  IF INSTEAD OF F1 .GT. TOLCON AND EMIN .LT. -TOLCON WE START WITH
! F1 .LT. -TOLCON AND EMIN .GT. TOLCON, THIS PROGRAM CAN STILL BE USED
! BY REPLACING F BY -F BEFORE RUNNING THE PROGRAM.
! THE SOLUTION PROCEDURE IS A MODIFICATION OF THE FOLLOWING:  BISECT THE
! INTERVAL [P1,PROCOR] TO GET A THIRD POINT, PASS A QUADRATIC POLYNOMIAL
! THROUGH THE THREE POINTS AND USE ITS UNIQUE ZERO IN [P1,PROCOR] TO
! REPLACE P1 OR PROCOR, MAINTAINING THE CONDITIONS F(LEFT ENDPOINT) .GT.
! TOLCON AND F(RIGHT ENDPOINT) .LT. -TOLCON, AND CONTINUE UNTIL A SOLUTION
! IS FOUND OR F HAS BEEN COMPUTED 5 TIMES OR THE INTERVAL LENGTH FALLS
! BELOW 100.0*B**(-ITT).  THIS PROCEDURE MAY BE ESPECIALLY USEFUL IN CASES
! WHERE F IS EXPENSIVE TO COMPUTE SINCE IT MAINTAINS A SHRINKING INTERVAL
! ABOUT THE SOLUTION, HAS A HIGHER ORDER OF CONVERGENCE THAN THE REGULA
! FALSI METHOD, AND REQUIRES NO DERIVATIVES.  THE FOLLOWING SAMPLE DRIVER
! PROGRAM AND SUBROUTINE FNSET ARE SET UP TO FIND PROCOR IN [-4.0D0,2.0D0]
! WITH ABS(F(PROCOR)) .LE. 0.001D0, WHERE F(X) = 2.0D0**(-X) - 0.5D0
! (THE EXACT SOLUTION IS PROCOR = 1.0D0, EMIN = 0.0D0).
!
! SAMPLE DRIVER AND FNSET FOR (A) MULLERS METHOD
 
      IMPLICIT NONE

      REAL*8 dvec , emin , err1 , f1 , fun , p1 , parwrk , procor ,     &
           & pttbl , tolcon , work , zwork
      INTEGER iwork
      DIMENSION dvec(1) , fun(1) , pttbl(1,1) , zwork(1) , err1(4) ,    &
              & parwrk(1) , iwork(17) , work(6)
      !OPEN (6,FILE='MULOUT')
      dvec(1) = 1.0D0
      zwork(1) = 0.0D0
      iwork(16) = -2
!*****BEGIN USER SETTABLE STATEMENTS 1 OF 2
      tolcon = 1.0D-3
      p1 = -2.0D0
      f1 = 3.5D0
      procor = 2.0D0
      emin = -0.25D0
!*****END USER SETTABLE STATEMENTS 1 OF 2
      WRITE (6,99001) tolcon , p1 , f1 , procor , emin
99001 FORMAT (/' TOLCON IS',E22.13//' INITIALLY P1 IS',E22.13,'  F1 IS',&
            & E22.13//' PROCOR IS',E22.13,'  EMIN IS',E22.13)
      CALL MULLER(0,1,1,dvec,fun,1,pttbl,1,1,zwork,tolcon,0,iwork,17,   &
                & work,6,parwrk,err1,p1,f1,procor,emin)
      WRITE (6,99002) procor , emin
99002 FORMAT (/' AFTER MULLER PROCOR IS',E22.13//' EMIN IS',E22.13)
      END
!*==FNSET.spg  processed by SPAG 6.72Dc at 18:11 on 18 Aug 2021
 
      SUBROUTINE FNSET(Nparm,Numgr,Pttbl,Iptb,Indm,Param,Ipt,Indfn,     &
                     & Icntyp,Confun)
      IMPLICIT NONE

      REAL*8 Confun , Param , Pttbl
      INTEGER Icntyp , Indfn , Indm , Ipt , Iptb , Nparm , Numgr
      DIMENSION Pttbl(Iptb,Indm) , Param(Nparm) , Icntyp(Numgr) ,       &
              & Confun(Numgr,Nparm+1)
!*****BEGIN USER SETTABLE STATEMENTS 2 OF 2
      Confun(1,1) = 2.0D0**(-Param(1)) - 0.5D0
!*****END USER SETTABLE STATEMENTS 2 OF 2
      END
 
 
! OUTPUT FOR (A) MULLERS METHOD
!
! TOLCON IS   0.1000000000000E-02
!
! INITIALLY P1 IS  -0.2000000000000E+01  F1 IS   0.3500000000000E+01
!
! PROCOR IS   0.2000000000000E+01  EMIN IS  -0.2500000000000E+00
!
! AFTER MULLER PROCOR IS   0.9993746852789E+00
!
! EMIN IS   0.2167645412207E-03