!********************************************************************************
!>
! (A) MULLERS METHOD DERIVATIVE FREE REAL ROOT FINDING
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

     module muller_test_module

      use conmax_module, only: conmax_solver
      use iso_fortran_env, only: wp => real64

      implicit none

      private

      type,extends(conmax_solver) :: my_solver
            contains
            procedure :: fnset => my_fnset
      end type my_solver

      public :: muller_test

      contains

      subroutine muller_test()

      implicit none

      real(wp) dvec , emin , err1 , f1 , fun , p1 , parwrk , procor ,     &
               pttbl , tolcon , work , zwork
      integer iwork
      dimension dvec(1) , fun(1) , pttbl(1,1) , zwork(1) , err1(4) ,    &
                parwrk(1) , iwork(17) , work(6)

      type(my_solver) :: solver

      !open (6,file='mulout')

      dvec(1)   = 1.0_wp
      zwork(1)  = 0.0_wp
      iwork(16) = -2
      tolcon    = 1.0e-3_wp
      p1        = -2.0_wp
      f1        = 3.5_wp
      procor    = 2.0_wp
      emin      = -0.25_wp

      write (6,99001) tolcon , p1 , f1 , procor , emin
99001 format (/' tolcon is',e22.13//' initially p1 is',e22.13,'  f1 is',&
            & e22.13//' procor is',e22.13,'  emin is',e22.13)

      call solver%muller(0,1,1,dvec,fun,1,pttbl,1,1,zwork,tolcon,0,iwork,17,   &
                         work,6,parwrk,err1,p1,f1,procor,emin)

      write (6,99002) procor , emin
99002 format (/' after muller procor is',e22.13//' emin is',e22.13)

      end subroutine muller_test

      subroutine my_fnset(me,nparm,numgr,pttbl,iptb,indm,param,ipt,indfn, &
                          icntyp,confun)
      implicit none

      class(my_solver),intent(inout) :: me
      integer  :: nparm
      integer  :: numgr
      integer  :: iptb
      integer  :: indm
      real(wp) :: pttbl(iptb,indm)
      real(wp) :: param(nparm)
      integer  :: ipt
      integer  :: indfn
      integer  :: icntyp(numgr)
      real(wp) :: confun(numgr,nparm+1)

      confun(1,1) = 2.0_wp**(-param(1)) - 0.5_wp

      end subroutine my_fnset

      end module muller_test_module

      program main
            use muller_test_module
            implicit none
            call muller_test()
      end program main


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
