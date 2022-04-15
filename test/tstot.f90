!********************************************************************************
!>
! THE EXAMPLE IS TO CHOOSE PARAMETERS A, B, C, AND D TO MINIMIZE
!
! MAX{ MAX (|F(X,Y) - (AX+BY+C)/(DX+1)|, |(AX+BY+C)/(DX+1)|) :
! (X,Y) IN Z}
!
! SUBJECT TO THE CONSTRAINTS
!
! DX+1 >= EPS FOR (X,Y) IN Z,
!
! AND THE FIRST PARTIAL DERIVATIVE OF (AX+BY+C)/(DX+1) WITH RESPECT TO
! X IS <= 0.0 AT (X,Y) = (0.0,0.0).
!
! HERE WE ARE TAKING Z = {(0.0,0.0), (0.0,1.0), (-2.0/3.0,1.0/3.0),
! (1.0,1.0), (1.0,2.0)}, EPS = .0001, F(0.0,0.0) = .5, F(0.0,1.0) = 1.0,
! F(-2.0/3.0,1.0/3.0) = -1.0, F(1.0,1.0) =1.5, F(1.0,2.0) =-1.0,
! AND TAKING THE INITIAL GUESSES FOR THE PARAMETERS TO BE
! A = B = C = D = 0.0
!
! TO USE CONMAX, WE WRITE THIS PROBLEM AS THE OPTIMIZATION PROBLEM
!
! MINIMIZE W, SUBJECT TO
!
! |F(X,Y) - (AX+BY+C)/(DX+1)| <= W, (X,Y) IN Z,    (TYPE 2)
!
!            (AX+BY+C)/(DX+1) <= W, (X,Y) IN Z,    (TYPE 1)
!
!           -(AX+BY+C)/(DX+1) <= W, (X,Y) IN Z,    (TYPE 1)
!
!               -DX - 1 + EPS <= 0, (X,Y) IN Z,    (TYPE -1)
!
!                      A - CD <= 0.                (TYPE -2)
!
! ONE CAN PROVE THAT THE UNIQUE BEST VALUES ARE A = -B = C = D = 1.0,
! WITH W (= ENORM = ERROR(NUMGR+1)) = 1.0.  ONE CAN ALSO PROVE THAT
! THERE ARE LOCAL SOLUTIONS WHICH ARE NOT GLOBAL SOLUTIONS, THAT IS,
! SOLUTIONS FOR WHICH W CANNOT BE REDUCED BY SMALL CHANGES IN A, B, C,
! AND D, BUT FOR WHICH A, B, C, D CAN BE FOUND SATISFYING THE CONSTRAINTS
! AND GIVING SMALLER W.  SOME SUCH SOLUTIONS ARE GIVEN BY A = B = D = 0.0,
! C = 0.25, WITH W = 1.25, AND OTHER CHOICES WHERE THE RATIONAL FUNCTION
! REDUCES TO THE CONSTANT 0.25 AND THE COEFFICIENTS SATISFY THE CONSTRAINTS.
! WHEN THE PROGRAM IS RUN, ANY OF THESE SOLUTIONS MAY BE FOUND (UP TO A
! SMALL DISCREPANCY DUE IN PART TO ROUNDOFF), DEPENDING ON THE ACCURACY OF
! THE COMPUTER BEING USED.
!
! THE OUTPUT FOR THE SAMPLE DRIVER AND FNSET FOLLOWS.  THIS WAS RUN ON
! A MAC SE (WITH ONE MEGABYTE OF RAM) WITH D1MACH(3) SET TO 16.0D0**(-14)
! USING THE ABSOFT MACFORTRAN COMPILER (VERSION 2.4);  RUNS WITH A
! DIFFERENT MACHINE AND/OR A DIFFERENT D1MACH(3) AND/OR A DIFFERENT
! COMPILER COULD PRODUCE DIFFERENT RESULTS, ESPECIALLY CONSIDERING THE
! POSSIBILITY OF LOCAL SOLUTIONS WHICH ARE NOT GLOBALY BEST.
!
! IOPTN IS    0  NPARM IS   4  NUMGR IS  21
!
! ITLIM IS  100  IFUN IS  5  IPTB IS  6  INDM IS  2
!
! THE FUNCTION VALUES ARE
!    0.50000E+00    0.10000E+01   -0.10000E+01    0.15000E+01
!   -0.10000E+01
!
! THE POINTS ARE
!    0.00000E+00    0.00000E+00
!    0.00000E+00    0.10000E+01
!   -0.66667E+00    0.33333E+00
!    0.10000E+01    0.10000E+01
!    0.10000E+01    0.20000E+01
!
! EPS IS    0.10000E-03
!
! THE INITIAL PARAMETERS ARE
!
!   0.00000000000000E+00   0.00000000000000E+00   0.00000000000000E+00
!
!   0.00000000000000E+00
!
! *****AFTER CONMAX ITER IS  13  LIWRK IS  178  LWRK IS   720
!
! THE FINAL PARAMETERS ARE
!
!  -0.90123453614326E-02  -0.12378422071992E-12   0.25000000000006E+00
!
!  -0.36049381446237E-01
!
! THE ERROR NORMS ARE
!
!    0.1250000000000E+01   -0.9638506185538E+00    0.1288396472843E-12
!
! THE ERRORS ARE
!
!   0.24999999999994E+00   0.75000000000006E+00  -0.12499999999999E+01
!   0.12499999999999E+01  -0.12499999999999E+01   0.25000000000006E+00
!  -0.25000000000006E+00   0.24999999999994E+00  -0.24999999999994E+00
!   0.24999999999994E+00  -0.24999999999994E+00   0.25000000000006E+00
!  -0.25000000000006E+00   0.24999999999994E+00  -0.24999999999994E+00
!  -0.99990000000000E+00  -0.99990000000000E+00  -0.10239329209642E+01
!  -0.96385061855376E+00  -0.96385061855376E+00   0.12883964728427E-12
!

!
! THIS TEST DRIVER PROGRAM AND FNSET ARE SET UP TO CHOOSE REAL (DOUBLE
! PRECISION) PARAMETERS A, B, C, AND D TO MINIMIZE
!
! MAX{ MAX (|F(X,Y) - (AX+BY+C)/(DX+1)|, |(AX+BY+C)/(DX+1)|) : (X,Y) IN Z}
!
! SUBJECT TO THE CONSTRAINTS
!
! DX+1 >= EPS FOR (X,Y) IN Z,
!
! AND THE FIRST PARTIAL DERIVATIVE OF (AX+BY+C)/(DX+1) WITH RESPECT TO
! X IS <= 0.0 AT (X,Y) = (0.0,0.0).
!
! HERE WE ARE TAKING:
!
!  * Z = {(0.0,0.0), (0.0,1.0), (-2.0/3.0,1.0/3.0), (1.0,1.0), (1.0,2.0)},
!  * EPS = 0.0001
!  * F(0.0,0.0)          =  0.5
!  * F(0.0,1.0)          =  1.0
!  * F(-2.0/3.0,1.0/3.0) = -1.0
!  * F(1.0,1.0)          =  1.5
!  * F(1.0,2.0)          = -1.0
!
! AND TAKING THE INITIAL GUESSES FOR THE PARAMETERS TO BE:
!
!  * A = B = C = D = 0.0

module tstot_test_module

    use conmax_module, only: conmax_solver
    use iso_fortran_env, only: wp => real64, output_unit

    implicit none

    private

    type, extends(conmax_solver) :: my_solver
    contains
        procedure :: fnset => my_fnset
    end type my_solver

    public :: tstot_test

contains

    subroutine tstot_test()

        implicit none

        real(wp), parameter :: zero  = 0.0_wp
        real(wp), parameter :: one   = 1.0_wp
        real(wp), parameter :: two   = 2.0_wp
        real(wp), parameter :: three = 3.0_wp
        integer, parameter  :: nwrit = output_unit

        real(wp) :: eps, error(24), fun(5), param(4), pttbl(6, 2), &
                    work(720)
        integer :: i, ifun, indm, ioptn, iptb, iter, itlim, &
                   iwork(178), j, liwrk, lwrk, nparm, numgr
        type(my_solver) :: solver

        !OPEN(6,FILE='TSTOT')

        ! SET PARAMETERS FOR CONMAX.
        !
        ! SET IOPTN=0 SINCE NO EXTRA OPTIONS ARE TO BE USED.
        ioptn = 0
        !
        ! SET NPARM=4 SINCE THERE ARE 4 PARAMETERS (VARIABLES) A, B, C, AND D.
        nparm = 4
        !
        ! SET NUMGR=21 SINCE THERE ARE 21 CONSTRAINTS (5 OF TYPE 2, 10 OF TYPE 1,
        ! 5 OF TYPE -1, AND ONE OF TYPE -2).
        numgr = 21
        !
        ! SET ITLIM=100, OR WHATEVER LIMIT IS DESIRED ON THE NUMBER OF ITERATIONS.
        itlim = 100
        !
        ! SET IFUN=5 (OR GREATER) SINCE THERE ARE 5 TYPE 2 CONSTRAINTS, AND USE
        ! THIS NUMBER AS THE DIMENSION OF FUN ABOVE.
        ifun = 5
        !
        ! SET IPTB=6 (OR GREATER) AND SET INDM=2 (OR GREATER), SINCE THE GREATEST
        ! FIRST SUBSCRIPT WE WILL USE IN PTTBL IS 6, AND THE GREATEST SECOND
        ! SUBSCRIPT WE WILL USE IN PTTBL IS 2.  WE USE THESE NUMBERS TO DIMENSION
        ! PTTBL ABOVE.  NOTE THAT IT IS ESSENTIAL THAT THE FIRST DIMENSION OF
        ! PTTBL BE EXACTLY IPTB.
        iptb = 6
        indm = 2
        !
        ! SET LIWRK=178 (OR GREATER), AND USE THIS NUMBER TO DIMENSION LIWRK
        ! ABOVE, BECAUSE OF THE COMPUTATION 7*21 + 7*4 + 3 = 178 (SEE CONMAX
        ! USERS GUIDE FOR AN EXPLANATION OF THIS AND OF LWRK BELOW).
        liwrk = 178
        !
        ! SET LWRK=720 (OR GREATER), AND USE THIS NUMBER TO DIMENSION LWRK ABOVE,
        ! BECAUSE OF THE COMPUTATION 2*4**2 + 4*21*4 + 11*21 + 27*4 + 13 = 720.
        lwrk = 720
        !
        ! THE DIMENSION OF PARAM ABOVE MUST BE NPARM (I.E. 4) OR GREATER, AND THE
        ! DIMENSION OF ERROR ABOVE MUST BE NUMGR+3 (I.E. 24) OR GREATER.
        !
        ! SET THE VALUES OF THE FUNCTION F.
        fun(1) = one/two
        fun(2) = one
        fun(3) = -one
        fun(4) = three/two
        fun(5) = -one
        !
        ! SET THE COORDINATES OF THE FIVE POINTS.
        pttbl(1, 1) = zero
        pttbl(1, 2) = zero
        pttbl(2, 1) = zero
        pttbl(2, 2) = one
        pttbl(3, 1) = -two/three
        pttbl(3, 2) = one/three
        pttbl(4, 1) = one
        pttbl(4, 2) = one
        pttbl(5, 1) = one
        pttbl(5, 2) = two
        !
        ! PUT EPS IN PTTBL(6,1) FOR TRANSMITTAL TO FNSET.  THIS IS WHY WE NEEDED
        ! IPTB TO BE AT LEAST 6.
        eps = (5*two)**(-4)
        pttbl(6, 1) = eps
        !
        ! SET THE INITIAL GUESSES FOR THE PARAMETERS.
        param(1) = zero
        param(2) = zero
        param(3) = zero
        param(4) = zero
        !
        ! WRITE THE INITIAL DATA.
        write (nwrit, 99001) ioptn, nparm, numgr, itlim, ifun, iptb, &
            indm
99001   format(/' IOPTN IS', i5, '  NPARM IS', i4, '  NUMGR IS', &
                i4//' ITLIM IS', i5, '  IFUN IS', i3, '  IPTB IS', i3, &
                '  INDM IS', i3)
        write (nwrit, 99002) (fun(i), i=1, 5)
99002   format(/' THE FUNCTION VALUES ARE'/(4e15.5))
        write (nwrit, 99003)
99003   format(/' THE POINTS ARE')
        do i = 1, 5
            write (nwrit, 99004) (pttbl(i, j), j=1, indm)
99004       format(2e15.5)
        end do
        write (nwrit, '(/A,e15.5)') ' EPS IS', eps
        write (nwrit, '(/A/,*(/3e23.14))') ' THE INITIAL PARAMETERS ARE', param(1:nparm)

        ! NOW CALL CONMAX.
        call solver%solve(ioptn, nparm, numgr, itlim, fun, ifun, pttbl, iptb, indm, &
                          iwork, liwrk, work, lwrk, iter, param, error)

        ! WRITE THE OUTPUT.
        ! NOTE THAT WE HAVE DEFERRED WRITING LIWRK AND LWRK UNTIL AFTER CALLING
        ! CONMAX SINCE CONMAX WILL CHANGE THEM TO THE NEGATIVE OF THE SMALLEST
        ! ALLOWABLE VALUES AND RETURN IF THEY WERE TOO SMALL.
        write (nwrit, '(/A,I4,A,I5,A,I6)') ' *****AFTER CONMAX ITER IS', iter, &
            '  LIWRK IS', liwrk, '  LWRK IS', lwrk
        write (nwrit, '(/A/,*(/3e23.14))') ' THE FINAL PARAMETERS ARE', param(1:nparm)
        write (nwrit, '(/A//,3e23.13//,a//,*(3e23.14/))') ' THE ERROR NORMS ARE', &
            error(numgr + 1), error(numgr + 2), error(numgr + 3), &
            ' THE ERRORS ARE', error(1:numgr)
! 99009 format (/' THE ERROR NORMS ARE'//3e23.13//' THE ERRORS ARE'//(3e23.14))

    end subroutine tstot_test
!********************************************************************************

!********************************************************************************
!>
!  Function for the unit test.

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

        real(wp) :: a, b, c, d, eps, p, q, x, y
        integer :: ii, irem

        real(wp), parameter :: one = 1.0_wp
        real(wp), parameter :: zero = 0.0_wp

        ! we break fnset into sections based on the value of ipt, that is, on
        ! which constraint is being set.
        if (ipt <= 5) then
            ! here ipt <= 5 and we set a constraint of the form
            ! abs(f(x,y) - (ax+by+c)/(dx+1)) <= w.
            ! note that since this is a type 2 constraint we do not need to deal
            ! with the absolute value or the f(x,y) here.
            icntyp(ipt) = 2
            a = param(1)
            b = param(2)
            c = param(3)
            d = param(4)
            x = pttbl(ipt, 1)
            y = pttbl(ipt, 2)
            p = a*x + b*y + c
            q = d*x + one
            confun(ipt, 1) = p/q
            if (indfn > 0) then
                ! here ipt <= 5 and indfn=1, and we set the partial derivatives.
                confun(ipt, 2) = x/q
                confun(ipt, 3) = y/q
                confun(ipt, 4) = one/q
                confun(ipt, 5) = -p*x/(q*q)
                return
            end if
        elseif (ipt <= 15) then
            ! here 6 <= ipt <= 15 and if ipt is even we set the constraint
            ! (ax+by+c)/(dx+1) <= w, which is half of the constraint
            ! abs((ax+by+c)/(dx+1)) <= w, while if ipt is odd we set the constraint
            ! -(ax+by+c)/(dx+1) <= w, which is the other half of the constraint
            ! abs((ax+by+c)/(dx+1)) <= w.
            icntyp(ipt) = 1
            ii = (ipt - 4)/2
            a = param(1)
            b = param(2)
            c = param(3)
            d = param(4)
            x = pttbl(ii, 1)
            y = pttbl(ii, 2)
            p = a*x + b*y + c
            q = d*x + one
            irem = ipt - 4 - 2*ii
            if (irem <= 0) then
                ! here 6 <= ipt <= 15 and ipt is even.
                confun(ipt, 1) = p/q
                if (indfn > 0) then
                    confun(ipt, 2) = x/q
                    confun(ipt, 3) = y/q
                    confun(ipt, 4) = one/q
                    confun(ipt, 5) = -p*x/(q*q)
                    return
                end if
            else
                !
                ! here 6 <= ipt <= 15 and ipt is odd.
                confun(ipt, 1) = -p/q
                if (indfn > 0) then
                    confun(ipt, 2) = -x/q
                    confun(ipt, 3) = -y/q
                    confun(ipt, 4) = -one/q
                    confun(ipt, 5) = p*x/(q*q)
                    return
                end if
            end if
        elseif (ipt <= 20) then
            ! here 16 <= ipt <= 20 and we set a constraint of the form
            ! -dx - 1.0 + eps <= 0.0
            icntyp(ipt) = -1
            d = param(4)
            eps = pttbl(6, 1)
            ii = ipt - 15
            x = pttbl(ii, 1)
            confun(ipt, 1) = -d*x - one + eps
            if (indfn > 0) then
                confun(ipt, 2) = zero
                confun(ipt, 3) = zero
                confun(ipt, 4) = zero
                confun(ipt, 5) = -x
                return
            end if
        else
            ! here ipt=21 and we set the constraint
            ! (partial derivative of (ax+by+c)/(dx+1) with respect to x at
            ! (x,y) = (0.0,0.0)) <= 0.0,
            ! i.e. a - cd <= 0.0
            icntyp(ipt) = -2
            a = param(1)
            c = param(3)
            d = param(4)
            confun(ipt, 1) = a - c*d
            if (indfn > 0) then
                confun(ipt, 2) = one
                confun(ipt, 3) = zero
                confun(ipt, 4) = -d
                confun(ipt, 5) = -c
                return
            end if
        end if

    end subroutine my_fnset

end module tstot_test_module

program main
    use tstot_test_module
    implicit none
    call tstot_test()
end program main
