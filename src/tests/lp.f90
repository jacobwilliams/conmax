
! (C) FREE-VARIABLE INEQUALITY-CONSTRAINED LINEAR PROGRAMMING
!
! (SUBPROGRAMS INVOLVED:  SLNPRO, SJELIM)
!
! GIVEN POSITIVE INTEGERS M AND N WITH M .GE. N, AND A MATRIX V, THIS
! PROGRAM WILL ATTEMPT TO SOLVE THE LINEAR PROGRAMMING PROBLEM
! MAXIMIZE  -V(M+1,1)*X(1)-...-V(M+1,N)*X(N)
! SUBJECT TO  V(I,1)*X(1)+...+V(I,N)*X(N) .LE. V(I,N+1) FOR I=1,...,M.
! ON OUTPUT, INDIC WILL BE AN ERROR FLAG WHOSE VALUE WILL BE 0 FOR A
! NORMAL SOLUTION, NEGATIVE FOR A POSSIBLY INACCURATE SOLUTION, AND
! POSITIVE FOR A PROBABLE FAILURE.  THE METHOD IS AN ENHANCED VERSION OF
! THAT IN
! AVDEYEVA, L. I. AND ZUKHOVITSKIY, S. I., LINEAR AND CONVEX PROGRAMMING,
! SAUNDERS, PHILADELPHIA, 1966.
! THE FOLLOWING SAMPLE DRIVER IS SET UP TO MAXIMIZE -Y, SUBJECT TO
! X + Y .GE. 2.0D0 (I.E. -X - Y .LE. -2.0D0) AND X - Y .LE. 4.0DO.
! (THE EXACT SOLUTION IS X = 3.0D0, Y = -1.0D0, WITH MAXIMUM OBJECTIVE
! FUNCTION VALUE 1.0D0.)
!
! SAMPLE DRIVER FOR (C) LINEAR PROGRAMMING

    program main

      use conmax_module
      use iso_fortran_env, only: wp => real64

      IMPLICIT NONE

      INTEGER i , indic , ixrct , iycct , iyrct , j , m , mp1 , n ,     &
            & np1 , nparm , numgr
      real(wp) v , x , y
!  THE MINIMUM DIMENSIONS ARE V(NUMGR+2*NPARM+1,NPARM+2), X(NPARM+1),
!  IYCCT(NPARM+1), Y(NUMGR+2*NPARM), IXRCT(NUMGR+2*NPARM),
!  IYRCT(NUMGR+2*NPARM), WHERE NPARM .GE. N AND NUMGR .GE. M. THE FIRST
!  DIMENSION OF V MUST BE EXACTLY NUMGR+2*NPARM+1.
      DIMENSION v(7,4) , x(3) , iycct(3) , y(6) , ixrct(6) , iyrct(6)

      n = 2
      m = 2
      v(1,1) = -1.0D0
      v(1,2) = -1.0D0
      v(1,3) = -2.0D0
      v(2,1) = 1.0D0
      v(2,2) = -1.0D0
      v(2,3) = 4.0D0
      v(3,1) = 0.0D0
      v(3,2) = 1.0D0
      nparm = n
      numgr = m
      mp1 = m + 1
      np1 = n + 1
      v(mp1,np1) = 0.0D0
      iyrct(1) = -1

      !OPEN (6,FILE='LPOUT')

      WRITE (6,99001) m , n
99001 FORMAT (/' THERE ARE',I5,'  CONSTRAINTS AND',I5,                  &
             &'  VARIABLES'//' THE CONSTRAINT COEFFICIENTS AND',        &
             &' RIGHT SIDES ARE')
      DO i = 1 , m
         WRITE (6,99002) (v(i,j),j=1,np1)
99002    FORMAT (/(3E22.13))
      ENDDO
      WRITE (6,99003) (v(mp1,j),j=1,n)
99003 FORMAT (/' THE NEGATIVES OF THE OBJECTIVE FUNCTION', &
             &' COEFFICIENTS ARE'//(3E22.13))
      CALL SLNPRO(v,m,n,iyrct,y,ixrct,iycct,nparm,numgr,x,indic)
      WRITE (6,99004) indic , (x(i),i=1,n)
99004 FORMAT (/' AFTER SLNPRO THE ERROR FLAG IS',I4, &
             &'  THE VARIABLES ARE'//(3E22.13))
      WRITE (6,99005) v(mp1,np1)
99005 FORMAT (/' THE OBJECTIVE FUNCTION VALUE IS',E22.13)

    END program main

! OUTPUT FOR (C) LINEAR PROGRAMMING
!
! THERE ARE    2  CONSTRAINTS AND    2  VARIABLES
!
! THE CONSTRAINT COEFFICIENTS AND RIGHT SIDES ARE
!
!  -0.1000000000000E+01  -0.1000000000000E+01  -0.2000000000000E+01
!
!   0.1000000000000E+01  -0.1000000000000E+01   0.4000000000000E+01
!
! THE NEGATIVES OF THE OBJECTIVE FUNCTION COEFFICIENTS ARE
!
!   0.0000000000000E+00   0.1000000000000E+01
!
! AFTER SLNPRO THE ERROR FLAG IS   0  THE VARIABLES ARE
!
!   0.3000000000000E+01  -0.1000000000000E+01
!
! THE OBJECTIVE FUNCTION VALUE IS   0.1000000000000E+01
