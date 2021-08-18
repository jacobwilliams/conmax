!
! THE EXAMPLE IS TO CHOOSE (DOUBLE PRECISION) PARAMETERS A, B, C, AND D
! TO MINIMIZE
!
! MAX{ MAX (|F(X,Y) - (AX+BY+C)/(DX+1)|, |(AX+BY+C)/(DX+1)|) :
! (X,Y) IN Z}
!
! SUBJECT TO THE CONSTRAINTS
!
! DX+1 .GE. EPS FOR (X,Y) IN Z,
!
! AND THE FIRST PARTIAL DERIVATIVE OF (AX+BY+C)/(DX+1) WITH RESPECT TO
! X IS .LE. 0.0 AT (X,Y) = (0.0,0.0).
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
! |F(X,Y) - (AX+BY+C)/(DX+1)| .LE. W, (X,Y) IN Z,    (TYPE 2)
!
!            (AX+BY+C)/(DX+1) .LE. W, (X,Y) IN Z,    (TYPE 1)
!
!           -(AX+BY+C)/(DX+1) .LE. W, (X,Y) IN Z,    (TYPE 1)
!
!               -DX - 1 + EPS .LE. 0, (X,Y) IN Z,    (TYPE -1)
!
!                      A - CD .LE. 0.                (TYPE -2)
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


! THIS IS A TEST DRIVER PROGRAM FOR CONMAX.  FOR A DESCRIPTION OF
! CONMAX, PLEASE SEE THE CONMAX USERS GUIDE, WHICH APPEARS AT THE
! BEGINNING OF THIS PACKAGE.  FOR MORE INFORMATION ABOUT THE
! EXAMPLE WHICH IS SET UP IN THIS DRIVER PROGRAM AND IN SUBROUTINE
! FNSET, PLEASE SEE THE COMMENTS IN THESE TWO ROUTINES AS WELL AS
! THE COMMENTS IN THE AFOREMENTIONED USERS GUIDE.
!
! THIS TEST DRIVER PROGRAM AND FNSET ARE SET UP TO CHOOSE REAL (DOUBLE
! PRECISION) PARAMETERS A, B, C, AND D TO MINIMIZE
!
! MAX{ MAX (|F(X,Y) - (AX+BY+C)/(DX+1)|, |(AX+BY+C)/(DX+1)|) :
! (X,Y) IN Z}
!
! SUBJECT TO THE CONSTRAINTS
!
! DX+1 .GE. EPS FOR (X,Y) IN Z,
!
! AND THE FIRST PARTIAL DERIVATIVE OF (AX+BY+C)/(DX+1) WITH RESPECT TO
! X IS .LE. 0.0 AT (X,Y) = (0.0,0.0).
!
! HERE WE ARE TAKING Z = {(0.0,0.0), (0.0,1.0), (-2.0/3.0,1.0/3.0),
! (1.0,1.0), (1.0,2.0)}, EPS = .0001, F(0.0,0.0) = .5, F(0.0,1.0) = 1.0,
! F(-2.0/3.0,1.0/3.0) = -1.0, F(1.0,1.0) =1.5, F(1.0,2.0) =-ONE,
! AND TAKING THE INITIAL GUESSES FOR THE PARAMETERS TO BE
! A = B = C = D = 0.0

    program tstot

      implicit none

                  real*8 eps , error , fun , one , param , pttbl , three , two ,    &
                 & work , zero
            integer i , i1mach , ifun , indm , ioptn , iptb , iter , itlim ,  &
                  & iwork , j , liwrk , lwrk , nparm , nread , numgr , nwrit
            !
            dimension fun(5) , pttbl(6,2) , iwork(178) , work(720) , param(4) &
                    & , error(24)
      !
      !*****MAC INSERT
      !     OPEN(6,FILE='TSTOT')
      !*****END MAC INSERT
      !
      ! SET MACHINE AND PRECISION DEPENDENT CONSTANTS.
            one = 1.0d0
            zero = one - one
            two = one + one
            three = one + two
            nread = i1mach(1)
            nwrit = i1mach(2)
      !
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
            pttbl(1,1) = zero
            pttbl(1,2) = zero
            pttbl(2,1) = zero
            pttbl(2,2) = one
            pttbl(3,1) = -two/three
            pttbl(3,2) = one/three
            pttbl(4,1) = one
            pttbl(4,2) = one
            pttbl(5,1) = one
            pttbl(5,2) = two
      !
      ! PUT EPS IN PTTBL(6,1) FOR TRANSMITTAL TO FNSET.  THIS IS WHY WE NEEDED
      ! IPTB TO BE AT LEAST 6.
            eps = (5*two)**(-4)
            pttbl(6,1) = eps
      !
      ! SET THE INITIAL GUESSES FOR THE PARAMETERS.
            param(1) = zero
            param(2) = zero
            param(3) = zero
            param(4) = zero
      !
      ! WRITE THE INITIAL DATA.
            write (nwrit,99001) ioptn , nparm , numgr , itlim , ifun , iptb , &
                              & indm
      99001 format (/' IOPTN IS',i5,'  NPARM IS',i4,'  NUMGR IS',             &
                   &i4//' ITLIM IS',i5,'  IFUN IS',i3,'  IPTB IS',i3,         &
                   &'  INDM IS',i3)
            write (nwrit,99002) (fun(i),i=1,5)
      99002 format (/' THE FUNCTION VALUES ARE'/(4e15.5))
            write (nwrit,99003)
      99003 format (/' THE POINTS ARE')
            do i = 1 , 5
               write (nwrit,99004) (pttbl(i,j),j=1,indm)
      99004    format (2e15.5)
            enddo
            write (nwrit,99005) eps
      99005 format (/' EPS IS',e15.5)
            write (nwrit,99006) (param(j),j=1,nparm)
      99006 format (/' THE INITIAL PARAMETERS ARE'/(/3e23.14))
      !
      ! NOW CALL CONMAX.
            call conmax(ioptn,nparm,numgr,itlim,fun,ifun,pttbl,iptb,indm,     &
                      & iwork,liwrk,work,lwrk,iter,param,error)
      !
      ! WRITE THE OUTPUT.
      ! NOTE THAT WE HAVE DEFERRED WRITING LIWRK AND LWRK UNTIL AFTER CALLING
      ! CONMAX SINCE CONMAX WILL CHANGE THEM TO THE NEGATIVE OF THE SMALLEST
      ! ALLOWABLE VALUES AND RETURN IF THEY WERE TOO SMALL.
            write (nwrit,99007) iter , liwrk , lwrk
      99007 format (/' *****AFTER CONMAX ITER IS',i4,'  LIWRK IS',i5,         &
                   &'  LWRK IS',i6)
            write (nwrit,99008) (param(j),j=1,nparm)
      99008 format (/' THE FINAL PARAMETERS ARE'/(/3e23.14))
            write (nwrit,99009) error(numgr+1) , error(numgr+2) ,             &
                              & error(numgr+3) , (error(i),i=1,numgr)
      99009 format (/' THE ERROR NORMS ARE'//3e23.13//' THE ERRORS ARE'//     &
                  & (3e23.14))

    end program tstot
