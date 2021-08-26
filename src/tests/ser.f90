!********************************************************************************
!>
!  One-dimensional derivative-free quadratic search for a positive
!  local minimum
!
!  Given a function f of one variable (where f(x) is computed in subroutine
!  fnset as confun(1,1), with x = param(1)) and a positive number projct,
!  the program will attempt to locate a new projct with emin = f(projct),
!  with the new projct approximately giving a local minimum of f in
!  [(old projct)/1024, 1024*(old projct)].  if it fails to do this, the
!  program will return with projct = the abscissa found with smallest
!  emin = f(projct).  on output, nsrch will be the number of evaluations
!  of f that were done.
!
!  The solution procedure is a modification of the
!  following:  compute f(projct/2), f(projct), and f(2*projct), if the
!  conditions f(middle point) .le. f(left point) and f(middle point) .le.
!  f(right point) are not both satisfied then try to get this by computing
!  f at smaller (or larger) points at most 3 more times;  once the
!  conditions are satisfied, assuming the points are not too close to being
!  collinear, pass a quadratic polynomial through the three points and use
!  its unique minimum in the interval to replace one of the endpoints while
!  maintaining the two conditions, continuing until f has been computed 4
!  more times or the interval length falls below 100.0*b**(-itt) or the
!  points become nearly collinear.
!
!  The following sample driver and
!  subroutine fnset are set up to approximate a solution of the line search
!  problem of minimizing g((6.0d0,2.0d0) + projct*(-2.0d0,-1.0d0)), where
!  g(u,v) = 3.0d0*abs(u) + 2.0d0*abs(v).  we start with projct = 1.0d0,
!  then run the program again starting with the result of the first run.
!  (the exact solution is projct = 3.0d0, emin = 2.0d0;  convergence is
!  rather slow, mainly because f is not differentiable at the minimum.)

module search_test_module

    use conmax_module, only: conmax_solver
    use iso_fortran_env, only: wp => real64

    implicit none

    private

    type, extends(conmax_solver) :: my_solver
    contains
        procedure :: fnset => my_fnset
    end type my_solver

    public :: search_test

contains

    subroutine search_test()

        implicit none

        real(wp) :: emin, emin1, err1(4), error(4), fun(1), param(1),       &
                    parprj(1), parser(1), prjlim, projct, pttbl(1, 1), tol1, &
                    tolcon, work(42), x(2)
        integer :: iact(1), iwork(17), nsrch
        type(my_solver) :: solver

        real(wp), parameter :: spcmn = real(radix(1.0_wp), wp)**(-digits(1.0_wp))
            !! `d1mach3`: the smallest relative spacing

        !open (6,file='serout')

        tol1 = 100.0_wp*spcmn
        tolcon = sqrt(spcmn)
        iact(1) = 1
        iwork(7) = 1
        prjlim = 1.0_wp/spcmn
        param(1) = 0.0_wp
        x(1) = 1.0_wp
        projct = 1.0_wp

        write (6, 99001) projct
99001   format(/' initially projct is', e22.13)
        call solver%searsl(0, 1, 1, prjlim, tol1, x, fun, 1, pttbl, 1, 1, param, error, 2.0_wp, &
                           1, iact, 0, 1.0_wp, tolcon, 2.0_wp, 0, 0, iwork, 17, work, 42, err1, &
                           parprj, projct, emin, emin1, parser, nsrch)
        write (6, 99002) projct, emin, nsrch
99002   format(/' after searsl projct is', e22.13//' emin is', e22.13, &
               &'  nsrch is', i4)

    end subroutine search_test

    subroutine my_fnset(me, nparm, numgr, pttbl, iptb, indm, param, ipt, indfn, icntyp, confun)

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

        real(wp) :: u, v

        u = 6.0_wp + param(1)*(-2.0_wp)
        v = 2.0_wp + param(1)*(-1.0_wp)
        confun(1, 1) = 3.0_wp*abs(u) + 2.0_wp*abs(v)

    end subroutine my_fnset

end module search_test_module

program main
    use search_test_module
    implicit none
    call search_test()
end program main

! OUTPUT 1 FOR (B) ONE-DIMENSIONAL SEARCH
!
! INITIALLY PROJCT IS   0.1000000000000E+01
!
! AFTER SEARSL PROJCT IS   0.2956250000000E+01
!
! EMIN IS   0.2175000000000E+01  NSRCH IS   8
!
! OUTPUT 2 FOR (B) ONE-DIMENSIONAL SEARCH
!
! INITIALLY PROJCT IS   0.2956250000000E+01
!
! AFTER SEARSL PROJCT IS   0.3010732751581E+01
!
! EMIN IS   0.2085862012645E+01  NSRCH IS   7
