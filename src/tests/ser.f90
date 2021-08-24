
! (B) ONE-DIMENSIONAL DERIVATIVE-FREE QUADRATIC SEARCH FOR A POSITIVE
!     LOCAL MINIMUM
!
! (SUBPROGRAMS INVOLVED:  SEARSL, FNSET (USER SUPPLIED), ILOC,
! ERCMP1, RCHMOD, CORRCT, SEARCR, MULLER, WOLFE, CONENR, HOUSE,
! DOTPRD, REFWL;  ONLY THE FIRST FIVE OF THESE ARE ACTUALLY USED)
!
! GIVEN A FUNCTION F OF ONE VARIABLE (WHERE F(X) IS COMPUTED IN SUBROUTINE
! FNSET AS CONFUN(1,1), WITH X = PARAM(1)) AND A POSITIVE NUMBER PROJCT,
! THE PROGRAM WILL ATTEMPT TO LOCATE A NEW PROJCT WITH EMIN = F(PROJCT),
! WITH THE NEW PROJCT APPROXIMATELY GIVING A LOCAL MINIMUM OF F IN
! [(OLD PROJCT)/1024, 1024*(OLD PROJCT)].  IF IT FAILS TO DO THIS, THE
! PROGRAM WILL RETURN WITH PROJCT = THE ABSCISSA FOUND WITH SMALLEST
! EMIN = F(PROJCT).  ON OUTPUT, NSRCH WILL BE THE NUMBER OF EVALUATIONS
! OF F THAT WERE DONE.  THE SOLUTION PROCEDURE IS A MODIFICATION OF THE
! FOLLOWING:  COMPUTE F(PROJCT/2), F(PROJCT), AND F(2*PROJCT), IF THE
! CONDITIONS F(MIDDLE POINT) .LE. F(LEFT POINT) AND F(MIDDLE POINT) .LE.
! F(RIGHT POINT) ARE NOT BOTH SATISFIED THEN TRY TO GET THIS BY COMPUTING
! F AT SMALLER (OR LARGER) POINTS AT MOST 3 MORE TIMES;  ONCE THE
! CONDITIONS ARE SATISFIED, ASSUMING THE POINTS ARE NOT TOO CLOSE TO BEING
! COLLINEAR, PASS A QUADRATIC POLYNOMIAL THROUGH THE THREE POINTS AND USE
! ITS UNIQUE MINIMUM IN THE INTERVAL TO REPLACE ONE OF THE ENDPOINTS WHILE
! MAINTAINING THE TWO CONDITIONS, CONTINUING UNTIL F HAS BEEN COMPUTED 4
! MORE TIMES OR THE INTERVAL LENGTH FALLS BELOW 100.0*B**(-ITT) OR THE
! POINTS BECOME NEARLY COLLINEAR.  THE FOLLOWING SAMPLE DRIVER AND
! SUBROUTINE FNSET ARE SET UP TO APPROXIMATE A SOLUTION OF THE LINE SEARCH
! PROBLEM OF MINIMIZING G((6.0D0,2.0D0) + PROJCT*(-2.0D0,-1.0D0)), WHERE
! G(U,V) = 3.0D0*ABS(U) + 2.0D0*ABS(V).  WE START WITH PROJCT = 1.0D0,
! THEN RUN THE PROGRAM AGAIN STARTING WITH THE RESULT OF THE FIRST RUN.
! (THE EXACT SOLUTION IS PROJCT = 3.0D0, EMIN = 2.0D0;  CONVERGENCE IS
! RATHER SLOW, MAINLY BECAUSE F IS NOT DIFFERENTIABLE AT THE MINIMUM.)
!
! SAMPLE DRIVER AND FNSET FOR (B) ONE-DIMENSIONAL SEARCH
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
