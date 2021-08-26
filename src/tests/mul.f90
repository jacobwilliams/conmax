!********************************************************************************
!>
!  Mullers method derivative free real root finding
!
!  Given a function f of one variable (where f(x) is computed in subroutine
!  fnset as confun(1,1), with:
!
!  * `x = param(1))`
!  * a nonnegative tolerance `tolcon`
!  * two points `(p1,f1)` and `(procor,emin)` with:
!    * `f1 = f(p1)`
!    * `emin = f(procor)`
!    * `p1 < procor`
!    * `f1 > tolcon`
!    * `emin < -tolcon`
!
!  The program will attempt to locate a new `procor` with new `emin = f(procor)`
!  and `abs(emin) <= tolcon`. if it fails to do this, the program will return with
!  `procor` = the leftmost abscissa found with `emin = f(procor) < -tolcon`.
!
!  Note:  if instead of `f1 > tolcon` and `emin < -tolcon` we start with
!  `f1 < -tolcon` and `emin > tolcon`, this program can still be used
!  by replacing `f` by `-f` before running the program.
!
!  The solution procedure is a modification of the following:  bisect the
!  interval [p1,procor] to get a third point, pass a quadratic polynomial
!  through the three points and use its unique zero in [p1,procor] to
!  replace p1 or procor, maintaining the conditions f(left endpoint) >
!  tolcon and f(right endpoint) < -tolcon, and continue until a solution
!  is found or f has been computed 5 times or the interval length falls
!  below 100.0*b**(-itt).  this procedure may be especially useful in cases
!  where f is expensive to compute since it maintains a shrinking interval
!  about the solution, has a higher order of convergence than the regula
!  falsi method, and requires no derivatives.
!
!  The following sample driver
!  program and subroutine fnset are set up to find `procor` in [-4.0,2.0]
!  with abs(f(procor)) <= 0.001, where f(x) = 2.0**(-x) - 0.5
!  (the exact solution is procor = 1.0, emin = 0.0).

module muller_test_module

    use conmax_module, only: conmax_solver
    use iso_fortran_env, only: wp => real64

    implicit none

    private

    type, extends(conmax_solver) :: my_solver
    contains
        procedure :: fnset => my_fnset
    end type my_solver

    public :: muller_test

contains

    subroutine muller_test()

        implicit none

        real(wp) :: dvec(1), emin, err1(4), f1, fun(1), p1, parwrk(1), procor, &
                    pttbl(1, 1), tolcon, work(6), zwork(1)
        integer :: iwork(17)

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

        write (6, 99001) tolcon, p1, f1, procor, emin
99001   format(/' tolcon is', e22.13//' initially p1 is', e22.13, '  f1 is',&
              & e22.13//' procor is', e22.13, '  emin is', e22.13)

        call solver%muller(0, 1, 1, dvec, fun, 1, pttbl, 1, 1, zwork, tolcon, 0, iwork, 17, &
                           work, 6, parwrk, err1, p1, f1, procor, emin)

        write (6, 99002) procor, emin
99002   format(/' after muller procor is', e22.13//' emin is', e22.13)

    end subroutine muller_test

    subroutine my_fnset(me, nparm, numgr, pttbl, iptb, indm, param, ipt, indfn, &
                        icntyp, confun)
        implicit none

        class(my_solver), intent(inout) :: me
        integer, intent(in)                             :: Nparm
        integer, intent(in)                             :: Numgr
        integer, intent(in)                             :: Iptb
        integer, intent(in)                             :: Indm
        real(wp), dimension(Iptb, Indm), intent(in)     :: Pttbl
        real(wp), dimension(Nparm), intent(in)          :: Param
        integer, intent(in)                             :: Ipt
        integer, intent(in)                             :: Indfn
        integer, dimension(Numgr), intent(out)          :: Icntyp
        real(wp), dimension(Numgr, Nparm + 1), intent(out) :: Confun

        confun(1, 1) = 2.0_wp**(-param(1)) - 0.5_wp

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
