
! (D) LEAST-DISTANCE QUADRATIC PROGRAMMING
!
! (SUBPROGRAMS INVOLVED:  WOLFE, ILOC, CONENR, HOUSE, DOTPRD,
! REFWL)
!
! GIVEN POSITIVE INTEGERS M AND N, AND A MATRIX PMAT, THIS PROGRAM WILL
! ATTEMPT TO LOCATE AN N-DIMENSIONAL POINT WPT IN THE POLYHEDRON
! DETERMINED BY THE INEQUALITIES
! PMAT(1,J)*WPT(1)+...+PMAT(N,J)*WPT(N)+PMAT(N+1,J) .LE. 0.0 FOR J=1,...,M
! WHOSE DISTANCE WDIST FROM THE ORIGIN IS MINIMIZED.  ON OUTPUT, JFLAG
! WILL BE AN ERROR FLAG WHOSE VALUE WILL BE 0 FOR A NORMAL SOLUTION AND
! POSITIVE FOR A LIKELY FAILURE.  THE METHOD IS AN ENHANCED VERSION OF
! THAT IN
! WOLFE, PHILIP, FINDING THE NEAREST POINT IN A POLYTOPE, MATHEMATICAL
! PROGRAMMING 11 (1976), 128-149.
! THE FOLLOWING SAMPLE DRIVER IS SET UP TO FIND THE NEAREST POINT TO
! THE ORIGIN IN THE POLYHEDRON DEFINED BY X + Y + Z .GE. 2.0D0 (I.E.
! -X - Y - Z + 2.0D0 .LE. 0.0D0) AND Z .GE. 1.0D0 (I.E. -Z + 1.0D0 .LE.
! 0.0D0).  (THE EXACT SOLUTION IS (X,Y,Z) = (0.5D0,0.5D0,1.0D0) WITH
! DISTANCE SQRT(1.5D0) FROM THE ORIGIN.)
!
! SAMPLE DRIVER FOR (D) LEAST-DISTANCE QUADRATIC PROGRAMMING
!
! THE MINIMUM DIMENSIONS ARE PMAT(NPARM+1,NUMGR), PMAT1(NPARM+1,NUMGR),
! WPT(NPARM), ICOR(NPARM+1), R(NPARM+1), PTNR(NPARM+1), COEF(NUMGR),
! WCOEF(NUMGR), IWORK(4*NUMGR+5*NPARM+3), WORK(2*NPARM**2+4*NUMGR*NPARM
! +9*NUMGR+22*NPARM+10), WHERE NPARM .GE. N AND NUMGR .GE. M.  THE FIRST
! DIMENSIONS OF PMAT AND PMAT1 MUST BE EXACTLY NPARM+1.

program main

    use conmax_module
    use iso_fortran_env, only: wp => real64

    IMPLICIT NONE

    real(wp) :: coef(2), pmat(4, 2), pmat1(4, 2), ptnr(4), r(4), s, &
                wcoef(2), wdist, work(136), wpt(3)
    INTEGER :: i, icor(4), iwork(26), j, jflag, liwrk, lwrk, m, mp1, &
               n, ncor, nmaj, nmin, np1, nparm, numgr

    n = 3
    m = 2
    pmat(1, 1) = -1.0_wp
    pmat(2, 1) = -1.0_wp
    pmat(3, 1) = -1.0_wp
    pmat(4, 1) = 2.0_wp
    pmat(1, 2) = 0.0_wp
    pmat(2, 2) = 0.0_wp
    pmat(3, 2) = -1.0_wp
    pmat(4, 2) = 1.0_wp
    nparm = n
    numgr = m
    liwrk = 4*numgr + 5*nparm + 3
    lwrk = 2*nparm**2 + 4*numgr*nparm + 9*numgr + 22*nparm + 10
    mp1 = m + 1
    np1 = n + 1

    !OPEN (6,FILE='QPOUT')
    WRITE (6, 99001) m, n
99001 FORMAT(/' THERE ARE', I5, '  BOUNDARIES AND', I5, '  DIMENSIONS'// &
              ' THE BOUNDARY COEFFICIENTS ARE')
    DO j = 1, m
        WRITE (6, 99002) (pmat(i, j), i=1, np1)
99002   FORMAT(/(3E22.13))
    END DO
    CALL WOLFE(n, m, pmat, 0, s, ncor, icor, iwork, liwrk, work, lwrk, r, coef, &
               ptnr, pmat1, nparm, numgr, wcoef, wpt, wdist, nmaj, nmin, jflag)
    WRITE (6, 99003) jflag, wdist, (wpt(i), i=1, n)
99003 FORMAT(/' AFTER WOLFE THE ERROR FLAG IS', &
         &I4//' THE DISTANCE FROM THE ORIGIN IS', &
         &E22.13//' THE POINT HAS COORDINATES'//(3E22.13))

END

! OUTPUT FOR (D) LEAST-DISTANCE QUADRATIC PROGRAMMING
!
! THERE ARE    2  BOUNDARIES AND    3  DIMENSIONS
!
! THE BOUNDARY COEFFICIENTS ARE
!
!  -0.1000000000000E+01  -0.1000000000000E+01  -0.1000000000000E+01
!   0.2000000000000E+01
!
!   0.0000000000000E+00   0.0000000000000E+00  -0.1000000000000E+01
!   0.1000000000000E+01
!
! AFTER WOLFE THE ERROR FLAG IS   0
!
! THE DISTANCE FROM THE ORIGIN IS   0.1224744871392E+01
!
! THE POINT HAS COORDINATES
!
!   0.5000000000000E+00   0.5000000000000E+00   0.1000000000000E+01
!
