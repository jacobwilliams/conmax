
!********************************************************************************
!>
!  Module for CONMAX.
!
!### Notes
!  * We note here that the only active write statements in this package
!   are in the sample driver program, but some of those which have been
!   commented out (along with the statements  nwrit=output_unit) elsewhere
!   in the program could prove useful, especially the statement
!   write(nwrit,1100)... in subroutine conmax.)
!
!### Reference
!  * E. H. Kaufman Jr., D. J. Leeming & G. D. Taylor,
!    "An ODE-based approach to nonlinearly constrained minimax problems",
!    Numerical Algorithms, Volume 9, pages 25-37 (1995).
!    Note that since that paper was written, `tolcon` was deleted from the
!    argument list of [[conmax]].
!
!### History
!  * CONMAX, Version 1.7. This package was last revised on December 4, 1996.
!    from: http://www.netlib.org/opt/conmax.f
!  * Jacob Williams, August 2021 : Refactored into modern Fortran.

    module conmax_module

    use iso_fortran_env, only: wp => real64, output_unit

    implicit none

    private

    integer,parameter :: nwrit  = output_unit
    real(wp),parameter :: zero  = 0.0_wp
    real(wp),parameter :: one   = 1.0_wp
    real(wp),parameter :: two   = 2.0_wp
    real(wp),parameter :: three = 3.0_wp
    real(wp),parameter :: four  = 4.0_wp
    real(wp),parameter :: ten   = 10.0_wp

    real(wp),parameter :: spcmn = real(radix(1.0_wp),wp)**(-digits(1.0_wp))
                                        !! `d1mach3`: the smallest relative spacing
    real(wp),parameter :: big = one/spcmn

    type,abstract,public :: conmax_solver
        !! main CONMAX solver class. This class must be
        !! extended to define the users's `fnset` function.
        private
    contains
        private
        procedure,public :: solve => conmax  !! main solver routine
        procedure,public :: muller
        procedure,public :: searsl
        procedure(func),deferred,public :: fnset
        procedure :: ercmp1
        procedure :: rkcon
        procedure :: slpcon
        procedure :: corrct
        procedure :: rkpar
        procedure :: searcr
        procedure :: derst
        procedure :: setu1
        procedure :: pmtst
    end type conmax_solver

    abstract interface
        subroutine func(me,Nparm,Numgr,Pttbl,Iptb,Indm,Param,Ipt,Indfn,Icntyp,Confun)
        !! interface for the `fnset` function.
        !!
        !! THE FIRST EIGHT VARIABLES IN THE CALLING SEQUENCE FOR FNSET ARE FOR
        !! INPUT TO FNSET, WITH THE FIRST FIVE VARIABLES BEING EXACTLY AS THE
        !! USER SET THEM IN THE DRIVER PROGRAM.  IF THE TEN THOUSANDS DIGIT OF
        !! IOPTN WAS SET TO 0, FNSET SHOULD BE WRITTEN TO PLACE THE APPROPRIATE
        !! VALUES IN ICNTYP AND CONFUN USING THE PARAMETERS IN PARAM, AS FOLLOWS.
        !!
        !! ICNTYP(IPT) = THE TYPE OF THE IPT(TH) CONSTRAINT (I.E. 2, 1, -1,
        !!   OR -2), OR THE USER CAN SET ICNTYP(IPT)=0 AS A SIGNAL TO IGNORE
        !!   CONSTRAINT IPT.
        !!
        !! CONFUN(IPT,1) = THE APPROPRIATE VALUE AS DISCUSSED ABOVE.  (THIS CAN
        !!   BE LEFT UNDEFINED IF ICNTYP(IPT)=0.)
        !!
        !! IF INDFN=1 (WHICH IS THE ONLY POSSIBILITY OTHER THAN INDFN=0) THEN IN
        !! ADDITION TO THE ABOVE, FOR J=1,...,NPARM, FNSET SHOULD COMPUTE
        !!
        !! CONFUN(IPT,J+1) = THE VALUE OF THE PARTIAL DERIVATIVE WITH RESPECT
        !!   TO PARAM(J) OF THE FUNCTION WHOSE VALUE WAS COMPUTED IN
        !!   CONFUN(IPT,1) (UNLESS ICNTYP(IPT)=0, IN WHICH CASE THESE VALUES
        !!   NEED NOT BE COMPUTED).
        import :: wp,conmax_solver
        implicit none
        class(conmax_solver),intent(inout) :: me
        integer  :: Nparm
        integer  :: Numgr
        integer  :: Iptb
        integer  :: Indm
        real(wp) :: Pttbl(Iptb,Indm)
        real(wp) :: Param(Nparm)
        integer  :: Ipt
        integer  :: Indfn
        integer  :: Icntyp(Numgr)
        real(wp) :: Confun(Numgr,Nparm+1)
        end subroutine func
    end interface

    ! maybe just move these into the solver
    public :: wolfe
    public :: slnpro

    contains
!********************************************************************************

!********************************************************************************
!>
! CONMAX consists of two programs for solving the problem
!
!```
!     minimize w
!
!     subject to
!
!     abs(fun(i) - confun(i,1)) <= w      if icntyp(i)=2,
!
!     confun(i,1)               <= w      if icntyp(i)=1,
!
!     confun(i,1)               <= 0.0    if icntyp(i)=-1 or -2,
!```
!
! Where:
!
!  * if `icntyp(i)=-1` the constraint is linear (i.e. the left side
!    consists of a linear combination of the parameters in the vector array
!    param plus a constant)
!  * if `icntyp(i)=-2` the constraint may be nonlinear.
!
! Thus we are choosing the parameters to minimize the left sides of
! the type 2 and type 1 (i.e. primary) constraints subject to the
! type -1 and type -2 (i.e. standard) constraints.
!
! If there are no primary constraints the program will attempt to
! find a feasible point (that is, a vector param for which the type -1
! and type -2 constraints are satisfied within a tolerance tolcon
! described below) which is relatively close to the user supplied initial
! approximation, then return.
!
!### Notes
!  * Since `iwork` and `work` are changed by the program,
!    the user will need to reset these values each time conmax is called.

    subroutine conmax(me,Ioptn,Nparm,Numgr,Itlim,Fun,Ifun,Pttbl,Iptb, &
                      Indm,Iwork,Liwrk,Work,Lwrk,Iter,Param,Error)

    implicit none

    class(conmax_solver),intent(inout) :: me
    integer,  intent(in)    :: ifun   !! THIS IS THE DIMENSION OF FUN IN THE DRIVER PROGRAM.
                                      !! IT MUST BE >= THE LARGEST INDEX I FOR WHICH FUN(I) IS USED
                                      !! UNLESS NO FUN(I) IS USED, IN WHICH CASE IT MUST BE >= 1.
    integer,  intent(in)    :: Indm   !! THIS IS THE SECOND DIMENSION OF PTTBL IN THE DRIVER
                                      !! PROGRAM.  IT MUST BE >= THE LARGEST SUBSCRIPT J FOR WHICH A
                                      !! VALUE PTTBL(I,J) IS USED, AND MUST BE >= 1 IF NO SUCH VALUES
                                      !! ARE USED.
    integer,  intent(in)    :: ioptn    !! THIS IS THE OPTION SWITCH, WHICH SHOULD BE SET TO
                                        !! 0 UNLESS ONE OR MORE OF THE EXTRA OPTIONS DESCRIBED BELOW
                                        !! IS USED. THE USER HAS SEVERAL EXTRA OPTIONS WHICH ARE ACTIVATED BY SETTING
                                        !! IOPTN TO A VALUE OTHER THAN 0; MORE THAN ONE AT A TIME CAN BE USED.
                                        !! IN PARTICULAR:
                                        !!
                                        !! * IF THE ONES DIGIT OF IOPTN IS 0, THEN THE USER SHOULD GIVE FORMULAS
                                        !!   IN FNSET FOR COMPUTING THE PARTIAL DERIVATIVES WHEN INDFN=1 AS DESCRIBED
                                        !!   ABOVE.
                                        !!
                                        !! * IF THE USER SETS THE ONES DIGIT OF IOPTN TO 1, THEN INDFN WILL ALWAYS
                                        !!   BE 0 WHEN FNSET IS CALLED, AND THE PROGRAM WILL AUTOMATICALLY
                                        !!   APPROXIMATE THE PARTIAL DERIVATIVES WHEN REQUIRED USING THE FOLLOWING
                                        !!   CENTERED DIFFERENCE FORMULA:
                                        !!   PARTIAL DERIVATIVE WITH RESPECT TO THE JTH VARIABLE (WHERE 1 <= J
                                        !!   <= NPARM) OF THE FUNCTION WHOSE VALUE IS COMPUTED IN CONFUN(IPT,1)
                                        !!   (WHERE 1 <= IPT <= NUMGR) IS APPROXIMATELY THE VALUE OF THIS
                                        !!   FUNCTION WHEN THE JTH COMPONENT OF PARAM IS INCREASED BY DELT, MINUS
                                        !!   THE VALUE OF THIS FUNCTION WHEN THE JTH COMPONENT OF PARAM IS
                                        !!   DECREASED BY DELT, ALL DIVIDED BY 2.0*DELT, WHERE DELT = SQRT(B**(-ITT)),
                                        !!   WHERE B IS THE BASE FOR FOR FLOATING POINT NUMBERS AND ITT IS THE NUMBER
                                        !!   OF BASE B DIGITS IN THE MANTISSA.
                                        !!
                                        !! * IF THE TENS DIGIT OF IOPTN IS 0, THE PROGRAM WILL NOT GIVE UP
                                        !!   UNTIL EITHER AN ITERATION FAILS TO PRODUCE A REDUCTION ABS(ENCHG) OF
                                        !!   MORE THAN 100.0*B**(-ITT) IN THE PRINCIPAL ERROR NORM, OR ITLIM
                                        !!   ITERATIONS HAVE BEEN USED.
                                        !!
                                        !! * IF THE TENS DIGIT OF IOPTN IS 1, THE PROGRAM WILL ALSO GIVE UP IF
                                        !!   ABS(ENCHG) < ENCSM FOR LIMSM CONSECUTIVE STEPS IN THE MAIN PART OF
                                        !!   THE PROGRAM (I.E. AFTER TYPE -1 AND TYPE -2 FEASIBILITY HAVE BOTH BEEN
                                        !!   ACHIEVED).  THIS OPTION MAY BE USEFUL IN SAVING SOME TIME BY
                                        !!   FORESTALLING A LONG STRING OF ITERATIONS AT THE END OF A RUN WITH ONLY
                                        !!   A TINY IMPROVEMENT IN EACH ONE.  ENCSM AND LIMSM ARE TRANSMITTED TO
                                        !!   CONMAX IN WORK(1) AND IWORK(1) RESPECTIVELY.  WORK(1) AND IWORK(1) DO
                                        !!   NOT NEED TO BE ASSIGNED VALUES IF THE TENS DIGIT OF IOPTN IS 0.
                                        !!
                                        !! * IF THE HUNDREDS DIGIT OF IOPTN IS 0 OR 2, THEN THE INTERNAL PARAMETER
                                        !!   NSTEP WILL BE GIVEN THE DEFAULT VALUE 1.  NSTEP IS THE NUMBER OF
                                        !!   RUNGE-KUTTA STEPS USED IN EACH RK ITERATION.
                                        !!
                                        !! * IF THE HUNDREDS DIGIT OF IOPTN IS 1 OR 3, THEN THE VALUE OF NSTEP USED
                                        !!   WILL BE THE POSITIVE INTEGER VALUE PLACED IN IWORK(2) BY THE USER IN THE
                                        !!   DRIVER PROGRAM.  SETTING NSTEP LARGER THAN 1 MAY ALLOW THE PROGRAM TO
                                        !!   SUCCEED ON DIFFICULT PROBLEMS WHERE THE CONVERGENCE WOULD BE EXTREMELY
                                        !!   SLOW WITH NSTEP=1, BUT IT WILL GENERALLY CAUSE MORE COMPUTER TIME TO BE
                                        !!   USED IN EACH RK ITERATION.  IWORK(2) DOES NOT NEED TO BE ASSIGNED A
                                        !!   VALUE IF THE HUNDREDS DIGIT OF IOPTN IS 0 OR 2.  (NSTEP IS SOMETIMES
                                        !!   CALLED THE OVERDRIVE PARAMETER.)
                                        !!
                                        !! * IF THE HUNDREDS DIGIT OF IOPTN IS 0 OR 1, THEN THE INTERNAL PARAMETER
                                        !!   TOLCON WILL BE GIVEN THE DEFAULT VALUE SQRT(B**(-ITT)), WHERE B IS THE
                                        !!   BASE FOR FLOATING POINT NUMBERS AND ITT IS THE NUMBER OF BASE B DIGITS
                                        !!   IN THE MANTISSA.
                                        !!
                                        !! * IF THE HUNDREDS DIGIT OF IOPTN IS 2 OR 3, THEN THE VALUE OF TOLCON USED
                                        !!   WILL BE THE VALUE PLACED IN WORK(2) BY THE USER IN THE DRIVER PROGRAM.
                                        !!   THIS VALUE SHOULD ALWAYS BE NONNEGATIVE.  IF THERE ARE NO TYPE -2 OR -1
                                        !!   CONSTRAINTS THEN THE SETTING OF TOLCON WILL HAVE NO EFFECT, BUT IF
                                        !!   THERE ARE TYPE -2 OR -1 CONSTRAINTS THEN IN GENERAL SMALLER VALUES OF
                                        !!   TOLCON MAY GIVE MORE ACCURACY IN THE FINAL ANSWER, BUT MAY SLOW DOWN
                                        !!   OR PREVENT CONVERGENCE, WHILE LARGER VALUES OF TOLCON MAY HAVE THE
                                        !!   REVERSE EFFECT.  IN PARTICULAR, IF THE TYPE -2 AND -1 CONSTRAINTS
                                        !!   CANNOT BE SATISFIED SUMULTANEOUSLY WITH STRICT INEQUALITY (THIS CASE
                                        !!   COULD OCCUR, FOR EXAMPLE, IF AN EQUALITY CONSTRAINT G = 0.0 WAS
                                        !!   ENTERED AS THE TWO INEQUALITY CONSTRAINTS G <= 0.0 AND
                                        !!   -G <= 0.0), THEN SETTING TOLCON=0.0 WILL ALMOST CERTAINLY CAUSE THE
                                        !!   PROGRAM TO FAIL, SINCE AT EACH ITERATION THE PROGRAM WILL NOT ACCEPT
                                        !!   PROSPECTIVE NEW VALUES OF THE PARAMETERS UNLESS IT CAN CORRECT THEM
                                        !!   BACK INTO THE RELAXED FEASIBLE REGION WHERE CONFUN(IPT,1) <= TOLCON
                                        !!   FOR ALL THE TYPE -2 AND -1 CONSTRAINTS.
                                        !!
                                        !! * IF THE THOUSANDS DIGIT OF IOPTN IS 0, THE PROGRAM WILL USE THE RK METHOD
                                        !!   (WHICH INVOLVES APPLYING FOURTH ORDER RUNGE-KUTTA TO A SYSTEM OF
                                        !!   DIFFERENTIAL EQUATIONS), THEN IF THIS FAILS IT WILL TRY TO REDUCE
                                        !!   W WITH AN SLP STEP (I.E. SUCCESSIVE LINEAR PROGRAMMING WITH A TRUST
                                        !!   REGION), THEN GO BACK TO RK IF THE SLP STEP REDUCES W.
                                        !!
                                        !! * IF THE THOUSANDS DIGIT OF IOPTN IS 1, THE PROGRAM WILL USE SLP STEPS ONLY.
                                        !!   IN GENERAL, IN SOME PROBLEMS SLP WORKS VERY WELL, AND IN THOSE
                                        !!   PROBLEMS IT WILL USUALLY BE FASTER THAN RK BECAUSE IT REQUIRES FEWER
                                        !!   GRADIENT EVALUATIONS THAN RK, BUT IN OTHER PROBLEMS THE CONVERGENCE
                                        !!   OF SLP MAY BE EXCRUCIATINGLY SLOW, AND THE MORE POWERFUL RK METHOD
                                        !!   MAY BE MUCH FASTER.
                                        !!
                                        !! * IF THE THOUSANDS DIGIT OF IOPTN IS 2 THE PROGRAM WILL USE THE RK METHOD
                                        !!   ONLY, QUITTING WHEN RK CAN NO LONGER PRODUCE AN IMPROVEMENT.  THIS
                                        !!   MAY GIVE A LITTLE LESS ACCURACY THAN SETTING THE THOUSANDS DIGIT TO 0,
                                        !!   BUT MAY SAVE SIGNIFICANT COMPUTER TIME IN SOME CASES.
                                        !!
                                        !! * IF THE TEN THOUSANDS DIGIT OF IOPTN IS 0, THEN FNSET SHOULD BE WRITTEN AS
                                        !!   DESCRIBED ABOVE.
                                        !!
                                        !! * IF THE USER SETS THE TEN THOUSANDS DIGIT OF IOPTN TO 1, THEN FNSET SHOULD BE
                                        !!   WRITTEN AS DESCRIBED ABOVE EXCEPT THAT THE COMPUTATIONS SHOULD BE DONE
                                        !!   FOR ALL IPT=1,..,NUMGR INSTEAD OF FOR A SINGLE GIVEN VALUE OF IPT.
                                        !!   THIS OPTION MAY SAVE COMPUTER TIME IN PROBLEMS WHERE A LARGE PART OF
                                        !!   THE COMPUTATION IS THE SAME FOR DIFFERENT VALUES OF IPT, SINCE IT
                                        !!   AVOIDS UNNECESSARY REPITITION OF THIS COMMON COMPUTATION BY HAVING
                                        !!   THE LOOP OVER IPT IN FNSET INSTEAD OF OUTSIDE FNSET.
                                        !!   IF THE TEN THOUSANDS DIGIT OF IOPTN IS 1, EVEN MORE TIME CAN OFTEN BE
                                        !!   SAVED IF FNSET IS WRITTEN SO THAT ALL CONSTRAINTS ARE COMPUTED IF IPT=0,
                                        !!   BUT ONLY THE STANDARD (I.E. TYPE -1 OR -2) CONSTRAINTS ARE COMPUTED IF
                                        !!   IPT=-1.  NOTE THAT IF THE TEN THOUSANDS DIGIT OF IOPTN IS 0, THEN IPT
                                        !!   WILL BE POSITIVE WHENEVER FNSET IS CALLED, INDICATING THAT ONLY
                                        !!   CONSTRAINT IPT SHOULD BE COMPUTED.
                                        !!   THE DRAWBACK OF USING THIS OPTION IS THAT IN GENERAL SOME CONSTRAINT
                                        !!   VALUES AND DERIVATIVES WILL BE COMPUTED UNNECESSARILY, WHICH COULD COST
                                        !!   TIME IN SOME PROBLEMS.
    integer,  intent(in)    :: iptb   !! THIS IS THE FIRST DIMENSION OF PTTBL IN THE DRIVER
                                      !! PROGRAM.  IT MUST BE >= THE LARGEST SUBSCRIPT I FOR WHICH A
                                      !! VALUE PTTBL(I,J) IS USED, AND MUST BE >= 1 IF NO SUCH VALUES
                                      !! ARE USED.
    integer,  intent(in)    :: itlim  !! THIS IS THE LIMIT ON THE NUMBER OF ITERATIONS, I.E.
                                      !! THE LIMIT ON THE NUMBER OF TIMES THE PROGRAM REDUCES W.  IF
                                      !! ITLIM IS SET TO 0 THE PROGRAM WILL COMPUTE THE ERRORS FOR
                                      !! THE INITIAL APPROXIMATION AND STOP WITHOUT CHECKING
                                      !! FEASIBILITY.
    integer,  intent(in)    :: liwrk  !! THIS IS THE DIMENSION OF IWORK IN THE DRIVER PROGRAM.
                                      !! IT MUST BE AT LEAST 7*NUMGR + 7*NPARM + 3.  IF NOT, CONMAX WILL
                                      !! RETURN WITH THIS MINIMUM VALUE MULTIPLIED BY -1 AS A WARNING.
    integer,  intent(in)    :: lwrk   !! THIS IS THE DIMENSION OF WORK IN THE DRIVER PROGRAM.
                                      !! IT MUST BE AT LEAST 2*NPARM**2 + 4*NUMGR*NPARM + 11*NUMGR +
                                      !! 27*NPARM + 13.  IF NOT, CONMAX WILL RETURN WITH THIS MINIMUM
                                      !! VALUE MULTIPLIED BY -1 AS A WARNING.
    integer,  intent(in)    :: nparm  !! THIS IS THE NUMBER OF PARAMETERS IN THE PROBLEM.
                                      !! (THEY ARE STORED IN `PARAM` -- SEE BELOW.)
    integer,  intent(in)    :: numgr  !! THIS IS THE TOTAL NUMBER OF CONSTRAINTS.
    real(wp), intent(in)    :: fun(Ifun)    !! (VECTOR ARRAY OF DIMENSION IFUN)  THIS IS
                                            !! A VECTOR ARRAY IN WHICH DATA OR FUNCTION VALUES IN TYPE 2
                                            !! CONSTRAINTS (SEE ABOVE) CAN BE STORED.  FUN(I) NEED NOT BE
                                            !! ASSIGNED A VALUE IF IT IS NOT GOING TO BE USED.
    real(wp), intent(in)    :: pttbl(Iptb,Indm)  !! (MATRIX ARRAY OF DIMENSION (IPTB,INDM))
                                                 !! ROW I OF THIS ARRAY NORMALLY CONTAINS A POINT USED IN THE ITH
                                                 !! CONSTRAINT.  THE ENTRIES IN ROW I NEED NOT BE ASSIGNED VALUES IF
                                                 !! SUCH A POINT IS NOT USED IN THE ITH CONSTRAINT.
                                                 !! (EXAMPLE:  IF THE LEFT SIDE OF CONSTRAINT I IS A POLYNOMIAL IN
                                                 !! ONE INDEPENDENT VARIABLE, THEN THE VALUE OF THE INDEPENDENT
                                                 !! VARIABLE SHOULD BE IN PTTBL(I,1), AND THE COEFFICIENTS SHOULD BE
                                                 !! IN PARAM.)
                                                 !! PTTBL CAN ALSO BE USED TO PASS OTHER INFORMATION FROM THE DRIVER
                                                 !! PROGRAM TO SUBROUTINE FNSET.
    integer,  intent(inout) :: iwork(Liwrk)     !! (VECTOR ARRAY OF DIMENSION LIWRK)  THIS IS AN INTEGER
                                                !! WORK ARRAY.  THE USER NEED NOT PLACE ANY VALUES IN IT, EXCEPT
                                                !! POSSIBLY CERTAIN OPTIONAL INFORMATION AS DESCRIBED BELOW.
    real(wp), intent(inout) :: work(Lwrk)   !! (VECTOR ARRAY OF DIMENSION LWRK)  THIS IS A FLOATING
                                            !! POINT WORK ARRAY.  THE USER NEED NOT PLACE ANY VALUES IN IT,
                                            !! EXCEPT POSSIBLY CERTAIN OPTIONAL INFORMATION AS DESCRIBED BELOW.
    real(wp), intent(inout) :: param(Nparm) !! (VECTOR ARRAY OF DIMENSION AT LEAST NPARM
                                            !! IN THE DRIVER PROGRAM)  THE USER SHOULD PLACE AN INITIAL GUESS
                                            !! FOR THE PARAMETERS IN PARAM, AND ON OUTPUT PARAM WILL CONTAIN
                                            !! THE BEST PARAMETERS CONMAX HAS BEEN ABLE TO FIND.  IF THE
                                            !! INITIAL PARAM IS NOT FEASIBLE THE PROGRAM WILL FIRST TRY TO
                                            !! FIND A FEASIBLE PARAM.
    real(wp), intent(out)   :: error(Numgr+3)   !! (VECTOR ARRAY OF DIMENSION AT LEAST NUMGR + 3 IN THE
                                                !! DRIVER PROGRAM)  FOR I=1,...,NUMGR, CONMAX WILL PLACE IN
                                                !! ERROR(I) THE ERROR IN CONSTRAINT I (DEFINED TO BE THE VALUE
                                                !! OF THE LEFT SIDE OF CONSTRAINT I, EXCEPT WITHOUT THE ABSOLUTE
                                                !! VALUE IN TYPE 2 CONSTRAINTS).  FURTHER,
                                                !!
                                                !!  * ERROR(NUMGR+1) WILL BE THE (PRINCIPAL) ERROR NORM, THAT IS, THE
                                                !!    MAXIMUM VALUE OF THE LEFT SIDES OF THE TYPE 2 (INCLUDING THE
                                                !!    ABSOLUTE VALUE) AND TYPE 1 CONSTRAINTS.
                                                !!
                                                !!  * ERROR(NUMGR+2) WILL BE THE MAXIMUM VALUE OF THE LEFT SIDES OF THE
                                                !!    TYPE -1 CONSTRAINTS, OR 0.0 IF THERE ARE NO TYPE -1
                                                !!    CONSTRAINTS.  EXCEPT FOR ROUNDOFF ERROR AND SMALL TOLERANCES
                                                !!    IN SOME SUBROUTINES THIS VALUE WILL NORMALLY BE <= 0.0, AND
                                                !!    IT WILL NOT BE ALLOWED TO BE > TOLCON IN THE MAIN PART OF
                                                !!    THE PROGRAM.
                                                !!
                                                !!  * ERROR(NUMGR+3) WILL BE THE MAXIMUM VALUE OF THE LEFT SIDES OF THE
                                                !!    TYPE -2 CONSTRAINTS, OR 0.0 IF THERE ARE NO TYPE -2
                                                !!    CONSTRAINTS.  THIS VALUE SHOULD BE <= TOLCON, SINCE THE
                                                !!    PROGRAM WILL NOT EVEN ATTEMPT TO COMPUTE VALUES FOR THE
                                                !!    TYPE 2 AND TYPE 1 CONSTRAINTS OTHERWISE (EXCEPT FOR VALUES
                                                !!    CORRESPONDING TO THE INITIAL PARAMETERS PLACED IN PARAM BY
                                                !!    THE USER).  THE USER CAN USE THIS FEATURE TO INSERT TYPE -2
                                                !!    OR -1 CONSTRAINTS TO KEEP THE PARAMETERS AWAY FROM VALUES
                                                !!    WHERE A TYPE 2 OR TYPE 1 CONSTRAINT IS UNDEFINED.
    integer,  intent(out)   :: iter !! THIS IS THE NUMBER OF ITERATIONS PERFORMED BY CONMAX,
                                    !! INCLUDING THOSE USED IN ATTEMPTING TO GAIN FEASIBILITY,
                                    !! UNTIL EITHER IT CAN NO LONGER IMPROVE THE SITUATION OR THE
                                    !! ITERATION LIMIT IS REACHED.  IF ITER=ITLIM IT IS POSSIBLE
                                    !! THAT THE PROGRAM COULD FURTHER REDUCE W IF RESTARTED
                                    !! (POSSIBLY WITH THE NEW PARAMETERS).
                                    !!
                                    !!  * ITER=-1 IS A SIGNAL THAT TYPE -1 FEASIBILITY COULD NOT BE
                                    !!    ACHIEVED, IN THIS CASE ERROR WILL CONTAIN THE VALUES COMPUTED
                                    !!    USING THE USER SUPPLIED INITIAL APPROXIMATION.
                                    !!  * ITER=-2 IS A SIGNAL THAT TYPE -1 FEASIBILITY WAS ACHIEVED
                                    !!    BUT TYPE -2 FEASIBILITY COULD NOT BE ACHIEVED,
                                    !!    IN THIS CASE VALUES IN ERROR CORRESPONDING TO TYPE 1 OR
                                    !!    TYPE 2 CONSTRAINTS WILL BE ZERO.

    real(wp) :: enc1 , enchg , encsm , enor2 , enor3 , enorm , &
                prjslp , projct , &
                rchdnk , rchdwn , rchin , s
    real(wp) :: tolcon , tollin , wdist
    integer :: i , i1 , ii , ilc02 , ilc06 , ilc08 , ilc11 , &
                ilc12 , ilc13 , ilc14 , ilc15 , ilc17 , ilc20 , ilc21 , &
                ilc22 , ilc24 , ilc25 , ilc26 , ilc27
    integer :: ilc29 , ilc30 , ilc31 , ilc33 , ilc35 , ilc40 , ilc42 , &
                ilc44 , ilc46 , iophun , iopten , ioptho , &
                ioptth , iphse , ipmax , ipt
    integer :: irk , ismax , isucc , itersl , itlim1 , &
                ityp1 , ityp1k , ityp2 , ityp2k , itypm1 , itypm2 , &
                j , j1 , j2 , jflag , jiwrk , jj
    integer :: jwrk , kntsm , l , l1 , l2 , limsm , m , &
                mact1 , ncor , nmaj , nmin , npar1 , nstep , &
                numlim

    ! check to see if the dimensions liwrk and lwrk are large enough.
    jiwrk = 7*Numgr + 7*Nparm + 3
    jwrk = 2*Nparm**2 + 4*Numgr*Nparm + 11*Numgr + 27*Nparm + 13
    if (Liwrk<jiwrk .or. Lwrk<jwrk) then

        iter = -999  ! error flag
        return

    else

        ! INITIALIZE SOME OTHER PARAMETERS.
        npar1 = Nparm + 1
        isucc = 0
        Iter = 0
        itersl = 0
        itlim1 = Itlim
        enchg = zero
        ilc02 = iloc(2,Nparm,Numgr)
        ilc06 = iloc(6,Nparm,Numgr)
        ilc08 = iloc(8,Nparm,Numgr)
        ilc11 = iloc(11,Nparm,Numgr)
        ilc12 = iloc(12,Nparm,Numgr)
        ilc13 = iloc(13,Nparm,Numgr)
        ilc14 = iloc(14,Nparm,Numgr)
        ilc15 = iloc(15,Nparm,Numgr)
        ilc17 = iloc(17,Nparm,Numgr)
        ilc20 = iloc(20,Nparm,Numgr)
        ilc21 = iloc(21,Nparm,Numgr)
        ilc22 = iloc(22,Nparm,Numgr)
        ilc24 = iloc(24,Nparm,Numgr)
        ilc25 = iloc(25,Nparm,Numgr)
        ilc26 = iloc(26,Nparm,Numgr)
        ilc27 = iloc(27,Nparm,Numgr)
        ilc29 = iloc(29,Nparm,Numgr)
        ilc30 = iloc(30,Nparm,Numgr)
        ilc31 = iloc(31,Nparm,Numgr)
        ilc33 = iloc(33,Nparm,Numgr)
        ilc35 = iloc(35,Nparm,Numgr)
        ilc40 = iloc(40,Nparm,Numgr)
        ilc42 = iloc(42,Nparm,Numgr)
        ilc44 = iloc(44,Nparm,Numgr)
        ilc46 = iloc(46,Nparm,Numgr)

        ! IF THE TENS DIGIT OF IOPTN IS 1, SET KNTSM TO 0 AND GET ENCSM
        ! FROM WORK(1) AND LIMSM FROM IWORK(1).
        iopten = (Ioptn-(Ioptn/100)*100)/10
        if ( iopten>0 ) then
            kntsm = 0
            encsm = Work(1)
            limsm = Iwork(1)
        end if

        ! IF THE HUNDREDS DIGIT OF IOPTN IS 1 OR 3, SET NSTEP = IWORK(2),
        ! AND OTHERWISE SET NSTEP TO ITS DEFAULT VALUE OF 1.
        iophun = (Ioptn-(Ioptn/1000)*1000)/100
        if ( iophun<=(iophun/2)*2 ) then
            nstep = 1
        else
            nstep = Iwork(2)
        end if

        ! IF THE HUNDREDS DIGIT OF IOPTN IS 2 OR 3, SET TOLCON = WORK(2),
        ! AND OTHERWISE SET TOLCON TO ITS DEFAULT VALUE OF SQRT(SPCMN).
        if ( iophun<2 ) then
            tolcon = sqrt(spcmn)
        else
            tolcon = Work(2)
        end if

        ! IN THIS VERSION OF CONMAX WE SET THE LINEAR CONSTRAINT TOLERANCE
        ! EQUAL TO THE NONLINEAR CONSTRAINT TOLERANCE.
        tollin = tolcon

        ! SET IRK=1 IF THE THOUSANDS DIGIT OF IOPTN IS 0 AND OTHERWISE SET IRK=0.
        ioptho = (Ioptn-(Ioptn/10000)*10000)/1000
        if ( ioptho<=0 ) then
            irk = 1
        else
            irk = 0
        end if

        ! COMPUTE THE TEN THOUSANDS DIGIT OF IOPTN FOR LATER USE.
        ioptth = (Ioptn-(Ioptn/100000)*100000)/10000

        ! SET IPHSE=-1 TO INDICATE WE HAVE NOT CHECKED TYPE -1 FEASIBILITY YET.
        iphse = -1
        ! SET RCHDWN = THE NUMBER OF LENGTHS OF PROJCT IN RKSACT (OR NUMBER OF
        ! LENGTHS OF BNDLGT IN SETU1) WE WILL GO BELOW ERROR(NUMGR+1) TO DECLARE
        ! A PRIMARY CONSTRAINT TO BE ACTIVE.
        rchdwn = two
        rchdnk = rchdwn
        ! SET RCHIN = THE NUMBER OF LENGTHS OF PROJCT (OR BNDLGT) WE WILL GO
        ! BELOW 0.0 TO DECLARE A TYPE -2 CONSTRAINT TO BE ACTIVE.
        rchin = two
        ! SET A NORMAL VALUE FOR NUMLIM FOR USE IN SLPCON.
        numlim = 11
    end if
    ! END OF PRELIMINARY SECTION.  THE STATEMENTS ABOVE THIS POINT WILL NOT
    ! BE EXECUTED AGAIN IN THIS CALL TO CONMAX.

! CALL ERCMP1 WITH ICNUSE=0 TO COMPUTE THE ERRORS, ERROR NORMS, AND ICNTYP.
! WE TAKE IPHSE AS 0 SO ALL CONSTRAINTS WILL BE COMPUTED BY FNSET IN CASE
! THE TEN THOUSANDS DIGIT OF IOPTN IS 1.
! THIS IS ONE OF ONLY TWO PLACES IN THE PROGRAM WHERE WE CALL ERCMP1 WITH
! ICNUSE=0, THE OTHER BEING STATEMENT 1415 BELOW..
 100  call me%ercmp1(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,Param,0,0, &
                  Iwork,Liwrk,Work(ilc08),Iwork(ilc17),ipmax,ismax,     &
                  Error)
! IF ITLIM=0 WE RETURN.
      if ( Itlim<=0 ) then
         return
      else
!
! COMPUTE ITYP2, ITYP1, ITYPM1, AND ITYPM2 AS THE NUMBER OF CONSTRAINTS  OF
! TYPE 2 (I.E. PRIMARY, ABS(FUN(I)-CONFUN(I,1)) <= W) OR 1 (I.E. PRIMARY,
! CONFUN(I,1) <= W) OR -1 (I.E. STANDARD LINEAR, CONFUN(I,1) <= 0.0)
! OR -2 (I.E. STANDARD NONLINEAR) RESPECTIVELY.
         ityp2 = 0
         ityp1 = 0
         itypm1 = 0
         itypm2 = 0
! NOTE THAT ARRAYS NOT IN THE CALLING SEQUENCE FOR CONMAX ARE ACCESSED
! THROUGH THEIR LOCATION IN IWORK OR WORK.  CONMAX IS THE ONLY
! SUBROUTINE IN WHICH THIS IS NECESSARY.
         do i = 1 , Numgr
            ii = ilc17 - 1 + i
! HERE IWORK(II)=ICNTYP(I).
            if ( Iwork(ii)<0 ) then
               if ( Iwork(ii)+1<0 ) then
                  itypm2 = itypm2 + 1
               else
                  itypm1 = itypm1 + 1
               end if
            else if ( Iwork(ii)/=0 ) then
               if ( Iwork(ii)<=1 ) then
                  ityp1 = ityp1 + 1
               else
                  ityp2 = ityp2 + 1
               end if
            end if
         end do
      end if
!
! COMPUTE THE ERROR NORMS.  ENORM IS THE PRINCIPAL ERROR NORM.
 200  enorm = Error(Numgr+1)
      enor2 = Error(Numgr+2)
      enor3 = Error(Numgr+3)
!
! WRITE ITER, ISUCC, IRK, ENCHG, AND THE ERROR NORMS.
!1050 WRITE(NWRIT,1100)ITER,ISUCC,IRK,ENCHG,ENORM,ENOR2,ENOR3
!     WRITE(9,1100)ITER,ISUCC,IRK,ENCHG,ENORM,ENOR2,ENOR3
!1100 FORMAT(/8H ITER IS,I5,10H  ISUCC IS,I4,8H  IRK IS,I4,
!    *10H  ENCHG IS,E24.14/9H ENORM IS,E24.14,10H  ENOR2 IS,E24.14/
!    *9H ENOR3 IS,E24.14)
!
!
! THE NEXT SECTION DETERMINES WHETHER WE WILL TERMINATE DUE TO ITERATION
! COUNT, AND IF SO FOR OUTPUT PURPOSES IT MODIFIES ITER (OR TWO OF THE
! ERROR NORMS IF THE FAILURE IS DUE TO INABILITY TO GAIN TYPE -2
! FEASIBILITY).
!
! IF IOPTEN=1 AND WE HAVE DONE AT LEAST ONE ITERATION IN THE MAIN PART
! OF CONMAX, WE WILL GIVE UP IF ABS(ENCHG) HAS BEEN LESS THAN ENCSM FOR
! LIMSM CONSECUTIVE MAIN ITERATIONS (INCLUDING THIS ONE).
 300  if ( iopten==1 ) then
         if ( iphse==0 ) then
            if ( Iter>0 ) then
               if ( -enchg<encsm ) then
                  kntsm = kntsm + 1
                  if ( kntsm>=limsm ) then
!
! FOR OUTPUT PURPOSES REPLACE ITER BY ITER + ITLIM - ITLIM1, THE TRUE
! NUMBER OF ITERATIONS COUNTING INITIALIZATION.  ITLIM - ITLIM1 WILL BE
! THE NUMBER OF ITERATIONS NEEDED TO GAIN TYPE -2 FEASIBILITY.  WORK
! DONE TO GAIN TYPE -1 FEASIBILITY IS NOT COUNTED AS AN ITERATION.
                     Iter = Iter + Itlim - itlim1
                     goto 500
                  end if
               else
                  kntsm = 0
               end if
            end if
         end if
      end if
!
      if ( Iter<itlim1 ) then
!
! HERE ITER < ITLIM1.  IF IPHSE = 0 OR -2 HERE WE GO INTO THE
! ITERATIVE PHASE OF CONMAX.
         if ( iphse+1/=0 ) goto 900
!
!
! HERE IPHSE=-1 AND WE CHECK TYPE -1 FEASIBILITY, TRY TO REGAIN IT IF
! WE DONT HAVE IT, CHECK TYPE -2 FEASIBILITY, AND SET UP FOR TYPE -2
! FEASIBILITY ITERATIONS IF WE DONT HAVE IT.  THE STATEMENTS FROM HERE
! DOWN TO THE TRIPLE BLANK LINE WILL BE EXECUTED AT MOST ONCE.
!
! NOTE THAT ENOR2=0.0 IF THERE ARE NO TYPE -1 CONSTRAINTS.
         if ( enor2<=tollin ) then
!
! HERE WE HAD TYPE -1 FEASIBILITY INITIALLY.
            if ( enor3>tolcon ) goto 700
            goto 800
         else
!
! HERE WE DO NOT HAVE TYPE -1 FEASIBILITY SO WE TRY TO GET IT.
! WE WILL NEED TO TELL DERST TO COMPUTE THE VALUES OF THE LEFT SIDES
! OF THE TYPE -1 CONSTRAINTS WITH THE VARIABLES EQUAL TO ZERO (I.E.
! THE CONSTANT TERMS IN THE CONSTRAINTS), SO WE SET PARWRK TO THE
! ZERO VECTOR TO CARRY THE MESSAGE.
            do j = 1 , Nparm
               jj = ilc27 - 1 + j
! HERE WORK(JJ) = PARWRK(J).
               Work(jj) = zero
            end do
            if ( ioptth>0 ) then
! HERE IOPTTH=1 AND WE CALL DERST WITH IPT=-1 TO PUT ALL THE STANDARD
! CONSTRAINT AND DERIVATIVE VALUES IN CONFUN.
! WE SET IPT=-1 TO TELL DERST IT NEED ONLY COMPUTE STANDARD CONSTRAINTS.
               ipt = -1
               call me%derst(Ioptn,Nparm,Numgr,Pttbl,Iptb,Indm,Work(ilc27),&
                             ipt,Work(ilc24),Work(ilc35),Iwork(ilc22),     &
                             Work(ilc08))
            end if
!
            m = 0
            do i = 1 , Numgr
               ii = ilc17 - 1 + i
! HERE WE CONSIDER ONLY TYPE -1 CONSTRAINTS.  THERE MUST BE AT LEAST
! ONE OF THESE, SINCE OTHERWISE WE WOULD NOT BE HERE ATTEMPTING TO
! GAIN TYPE -1 FEASIBILITY.
! HERE IWORK(II)=ICNTYP(I).
               if ( Iwork(ii)+1==0 ) then
                  m = m + 1
                  if ( ioptth<=0 ) then
! HERE IOPTTH=0 AND WE HAVE NOT YET CALLED DERST TO PUT CONSTRAINT I
! AND ITS DERIVATIVES IN CONFUN, SO WE DO IT NOW.
                     ipt = i
                     call me%derst(Ioptn,Nparm,Numgr,Pttbl,Iptb,Indm,      &
                                   Work(ilc27),ipt,Work(ilc24),Work(ilc35),&
                                   Iwork(ilc22),Work(ilc08))
                  end if
! COPY THE DERIVATIVES INTO PMAT FOR USE BY WOLFE.
                  do l = 1 , Nparm
                     l1 = ilc29 - 1 + l + (m-1)*npar1
                     l2 = ilc08 - 1 + i + l*Numgr
! HERE WORK(L1)=PMAT(L,M) AND WORK(L2)=CONFUN(I,L+1).
                     Work(l1) = Work(l2)
                  end do
!
! NOW THE ITH CONSTRAINT (WHICH IS ALSO THE MTH TYPE -1 CONSTRAINT) HAS
! THE FORM PMAT(1,M)*Z1+...+PMAT(NPARM,M)*ZNPARM + CONFUN(I,1)  <=
! 0.0.  WE MAKE THE CHANGE OF VARIABLES ZZ = Z - PARAM TO TRANSLATE THE
! ORIGIN TO PARAM.  THE ITH CONSTRAINT WILL THEN HAVE THE FORM
! PMAT(1,M)*ZZ1+...+PMAT(NPARM,M)*ZZNPARM + (CONFUN(I,1) + PMAT(1,M)*
! PARAM(1)+...+PMAT(NPARM,M)*PARAM(NPARM)) <= 0.0.  AFTER WOLFE FINDS
! THE CLOSEST POINT TO THE ORIGIN IN THE POLYHEDRON DEFINED BY THE NEW
! CONSTRAINTS, WE WILL ADD PARAM TO TRANSLATE BACK TO THE POINT WE WANT.
                  l1 = ilc29 - 1 + npar1 + (m-1)*npar1
                  l2 = ilc08 - 1 + i
! HERE WORK(L1)=PMAT(NPAR1,1) AND WORK(L2)=CONFUN(I,1).
                  Work(l1) = Work(l2)
                  do l = 1 , Nparm
                     l2 = ilc29 - 1 + l + (m-1)*npar1
! HERE WORK(L1)=PMAT(NPAR1,1) AND WORK(L2)=PMAT(L,M).
                     Work(l1) = Work(l1) + Work(l2)*Param(l)
                  end do
               end if
            end do
! CALL WOLFE WITH ISTRT=0 TO COMPUTE THE SOLUTION IN THE ZZ COORDINATE
! SYSTEM FROM SCRATCH.
            call wolfe(Nparm,m,Work(ilc29),0,s,ncor,Iwork(ilc15),Iwork, &
                       Liwrk,Work,Lwrk,Work(ilc33),Work(ilc06),         &
                       Work(ilc31),Work(ilc30),Nparm,Numgr,Work(ilc40), &
                       Work(ilc42),wdist,nmaj,nmin,jflag)
            if ( jflag>0 ) goto 600
!
! HERE JFLAG <= 0 AND WE PUT PARAM+WPT IN PARWRK TO CHECK WHETHER
! THE TYPE -1 CONSTRAINTS ARE NOW FEASIBLE WITHIN TOLLIN.
            do j = 1 , Nparm
               j1 = ilc27 - 1 + j
               j2 = ilc42 - 1 + j
! HERE WORK(J1)=PARWRK(J) AND WORK(J2)=WPT(J).
               Work(j1) = Param(j) + Work(j2)
            end do
! FOR USE IN ERCMP1 WE SET JCNTYP(I)=-1 IF ICNTYP(I)=-1 AND SET
! JCNTYP(I)=0 OTHERWISE.
            do i = 1 , Numgr
               ii = ilc17 - 1 + i
               jj = ilc21 - 1 + i
! HERE IWORK(II)=ICNTYP(I) AND IWORK(JJ)=JCNTYP(I).
               if ( Iwork(ii)+1/=0 ) then
                  Iwork(jj) = 0
               else
                  Iwork(jj) = -1
               end if
            end do
! CALL ERCMP1 WITH ICNUSE=1.
            call me%ercmp1(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,     &
                        Work(ilc27),1,iphse,Iwork,Liwrk,Work(ilc08),    &
                        Iwork(ilc21),ipmax,ismax,Work(ilc11))
            i1 = ilc11 - 1 + (Numgr+2)
! HERE WORK(I1)=ERR1(NUMGR+2).
            if ( Work(i1)>tollin ) goto 600
!
! HERE WE HAVE ACHIEVED TYPE -1 FEASIBILITY.  WE REPLACE PARAM WITH
! PARWRK.
            do j = 1 , Nparm
               jj = ilc27 - 1 + j
! HERE WORK(JJ)=PARWRK(J).
               Param(j) = Work(jj)
            end do
            ii = ilc11 - 1 + Numgr + 2
! HERE WORK(II)=ERR1(NUMGR+2).
!     WRITE(NWRIT,1397)WORK(II),(PARAM(J),J=1,NPARM)
!1397 FORMAT(48H TYPE -1 FEASIBILITY ACHIEVED.  ERR1(NUMGR+2) IS,
!    *E15.5,10H  PARAM IS/(4E20.12))
!
! IF THERE ARE TYPE -2 CONSTRAINTS, SET JCNTYP AS ICNTYP WITH ALL BUT -2
! VALUES ZEROED OUT AND CALL ERCMP1 WITH ICNUSE=1 TO CHECK TYPE -2
! FEASIBILITY.  WE CANNOT SIMPLY CHECK THE OLD ENOR3 HERE SINCE PARAM HAS
! BEEN CHANGED.  IF THERE ARE NO TYPE -2 CONSTRAINTS WE WILL AUTOMATICALLY
! HAVE TYPE -2 FEASIBILITY.
            if ( itypm2<=0 ) then
!
! HERE WE HAVE BOTH TYPE -1 AND TYPE -2 FEASIBILITY, BUT PARAM WAS
! CHANGED IN GETTING TYPE -1 FEASIBILITY, SO WE CALL ERCMP1
! WITH ICNUSE=0 (ICNUSE=1 WOULD WORK ALSO SINCE ICNTYP HAS NOT BEEN
! CHANGED HERE) TO GET THE NEW ERROR VECTOR.
               call me%ercmp1(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,  &
                           Param,0,iphse,Iwork,Liwrk,Work(ilc08),       &
                           Iwork(ilc17),ipmax,ismax,Error)
               goto 800
            else
               do i = 1 , Numgr
                  ii = ilc17 - 1 + i
                  jj = ilc21 - 1 + i
! HERE IWORK(II)=ICNTYP(I) AND IWORK(JJ)=JCNTYP(I).
                  if ( Iwork(ii)+1<0 ) then
                     Iwork(jj) = -2
                  else
                     Iwork(jj) = 0
                  end if
               end do
               call me%ercmp1(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,  &
                           Param,1,iphse,Iwork,Liwrk,Work(ilc08),       &
                           Iwork(ilc21),ipmax,ismax,Work(ilc11))
               ii = ilc11 - 1 + Numgr + 3
! HERE WORK(II)=ERR1(NUMGR+3).
               if ( Work(ii)>tolcon ) goto 700
               call me%ercmp1(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,  &
                           Param,0,iphse,Iwork,Liwrk,Work(ilc08),       &
                           Iwork(ilc17),ipmax,ismax,Error)
               goto 800
            end if
         end if
!
! HERE ITER = ITLIM1, SO WE RETURN.
      else if ( iphse>=0 ) then
         Iter = Iter + Itlim - itlim1
         goto 500
      end if
!
! HERE WE HAVE FAILED TO ACHIEVE TYPE -2 FEASIBILITY AND WE SET ITER=-2
! AS A WARNING, PUT ERROR(NUMGR+1) IN ITS PROPER LOCATION, SET
! ERROR(NUMGR+1) = 0.0 SINCE THE PRIMARY CONSTRAINTS WERE NOT COMPUTED,
! AND RETURN.  NOTE THAT WE CANNOT HAVE IPHSE=-1 HERE SINCE THAT WOULD
! IMPLY ITER=0, THUS ITLIM=ITLIM1=0, IN WHICH CASE WE WOULD HAVE
! TERMINATED EARLIER.
 400  Iter = -2
      Error(Numgr+3) = Error(Numgr+1)
      Error(Numgr+1) = zero
!     WRITE(6,1150)
!1150 FORMAT(43H ***WARNING  NONLINEAR STANDARD FEASIBILITY,
!    *16H NOT ACHIEVED***)
      return
!
 500  return
!
! HERE WE HAVE FAILED TO ACHIEVE TYPE -1 FEASIBILITY.  WE SET ITER=-1
! AS A WARNING AND RETURN.
 600  Iter = -1
!     WRITE(NWRIT,1360)
!1360 FORMAT(40H ***WARNING  LINEAR STANDARD FEASIBILITY,
!    *16H NOT ACHIEVED***)
      return
!
! HERE WE HAVE TYPE -1 FEASIBILITY BUT NOT TYPE -2 FEASIBILITY.  WE SET
! UP FOR THE TYPE -2 FEASIBILITY ITERATIONS, IN WHICH TYPE 1 AND TYPE
! 2 CONSTRAINTS ARE IGNORED AND TYPE -2 CONSTRAINTS ARE TREATED AS
! TYPE 1 CONSTRAINTS, EXCEPT WE WILL SWITCH OVER TO NORMAL ITERATIONS
! ONCE WE CAN FORCE W <= TOLCON.  THUS WE SET THE INDICATOR IPHSE TO
! -2, RESET ICNTYP(I) TO 1 IF IT WAS -2, LEAVE IT AT -1 IF IT WAS -1,
! AND SET IT TO 0 OTHERWISE, RESET ITYP2, ITYP1, AND ITYPM2, AND CALL
! ERCMP1 WITH ICNUSE=1 TO PUT THE PROPER VALUES IN ERROR.
 700  iphse = -2
      do i = 1 , Numgr
         ii = ilc17 - 1 + i
! HERE IWORK(II)=ICNTYP(I).
         if ( Iwork(ii)+1<0 ) then
            Iwork(ii) = 1
         else if ( Iwork(ii)+1/=0 ) then
            Iwork(ii) = 0
         end if
      end do
! SAVE ITYP2 AND ITYP1.
      ityp2k = ityp2
      ityp1k = ityp1
      ityp2 = 0
      ityp1 = itypm2
      itypm2 = 0
      call me%ercmp1(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,Param,1,   &
                  iphse,Iwork,Liwrk,Work(ilc08),Iwork(ilc17),ipmax,     &
                  ismax,Error)
      goto 900
!
! HERE WE HAVE BOTH TYPE -1 AND TYPE -2 FEASIBILITY, AND WE
! SET IPHSE=0 AND GO INTO THE MAIN PART OF CONMAX (UNLESS THERE WERE
! NO TYPE 1 OR TYPE 2 CONSTRAINTS, IN WHICH CASE WE RETURN).
 800  iphse = 0
      if ( ityp1+ityp2<=0 ) goto 500
!
! END OF INITIAL FEASIBILITY CHECKING, TYPE -1 FEASIBILITY WORK, AND
! TYPE -2 SETUP.  THE BLOCK OF STATEMENTS FROM HERE UP TO THE
! PRECEDING DOUBLE BLANK LINE WILL NOT BE EXECUTED AGAIN.
!
!
!
 900  if ( irk<=0 ) then
!
! HERE IRK IS 0 OR -1 AND WE DO AN SLP STEP.  IF SLPCON CANNOT REDUCE THE
! PRINCIPAL ERROR NORM ENORM = ERROR(NUMGR+1) BY MORE THAN 100.0*B**(-ITT)
! THEN IT WILL LEAVE PARAM AND ERROR UNCHANGED.
         call me%slpcon(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,tolcon, &
                     rchin,irk,itypm1,itypm2,Iwork(ilc17),rchdwn,numlim,&
                     itersl,prjslp,Work(ilc12),Iwork(ilc20),Work(ilc44),&
                     mact1,Iwork(ilc14),Iwork(ilc21),iphse,enchg,Iwork, &
                     Liwrk,Work,Lwrk,Work(ilc26),isucc,Param,Error)
      else
!
! HERE IRK IS 1 OR 2 AND WE DO AN RK STEP.  IF RKCON CANNOT REDUCE THE
! PRINCIPAL ERROR NORM ENORM = ERROR(NUMGR+1) BY MORE THAN 100.0*B**(-ITT)
! THEN IT WILL LEAVE PARAM AND ERROR UNCHANGED.
         call me%rkcon(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,tolcon,  &
                    rchin,Iter,irk,ityp2,ityp1,itypm1,itypm2,           &
                    Iwork(ilc17),projct,rchdwn,nstep,iphse,enchg,enc1,  &
                    Work(ilc29),Work(ilc12),Iwork,Liwrk,Work,Lwrk,      &
                    Iwork(ilc13),Work(ilc02),Work(ilc25),Work(ilc26),   &
                    Work(ilc46),Work(ilc11),Work(ilc08),isucc,Param,    &
                    Error)
      end if
!
      if ( isucc<=0 ) then
! HERE THE RK OR SLP STEP REDUCED ERROR(NUMGR+1) BY MORE THAN
! 100.0*B**(-ITT), AND WE INCREMENT ITER.
         Iter = Iter + 1
!
! IF EITHER IPHSE=0, OR IPHSE=-2 AND ERROR(NUMGR+1) > TOLCON, WE GO
! ON AS USUAL TO SET UP ANOTHER STEP WITH THE SAME IPHSE.
         if ( iphse<0 ) then
            if ( Error(Numgr+1)<=tolcon ) then
!
! HERE IPHSE=-2 AND ERROR(NUMGR+1) <= TOLCON, SO WE HAVE JUST ACHIEVED
! TYPE -2 FEASIBILITY.  WE WILL SET IPHSE=0, AND IF THERE ARE ANY
! PRIMARY CONSTRAINTS WE WILL RESET ITER, ITERSL, AND ITLIM1 (SINCE
! ITER=0 AND ITERSL=0 HAVE MEANINGS TO RKCON AND SLPCON RESPECTIVELY),
! RESET RCHIN AND RCHDWN, AND GO BACK TO THE FIRST ERCMP1 CALL TO
! RESTORE ERROR AND ICNTYP (ITYP1, ITYP2, ITYPM1, AND ITYPM2 WILL ALSO
! BE RESTORED).
               iphse = 0
               if ( ityp1k+ityp2k<=0 ) goto 500
               itlim1 = Itlim - Iter
               Iter = 0
               itersl = 0
               rchin = rchdwn
               rchdwn = rchdnk
               goto 100
            end if
         end if
!
         if ( irk<0 ) then
!
! HERE WE HAD AN SLP SUCCESS AND WE ARE GOING TO TRY RK AGAIN, SO WE SET
! IRK=2 TO WARN RKCON THAT THE SUCCESS CAME FROM SLP.
            irk = 2
         else if ( irk/=0 ) then
!
! HERE IRK IS 1 OR 2, SO WE JUST HAD AN RK SUCCESS.  WE RESET IRK AND
! ITERSL.
            irk = 1
            itersl = 0
            goto 200
         end if
! HERE WE HAD AN SLP SUCCESS AND WE INCREMENT ITERSL = THE NUMBER OF SLP
! SUCCESSES SINCE THE LAST SUCCESSFUL RK STEP (IF ANY).  ITERSL IS NEEDED
! IN SUBROUTINE BNDSET (CALLED BY SLPCON).
         itersl = itersl + 1
         goto 200
      else
!
! HERE RKCON OR SLPCON FAILED TO SIGNIFICANTLY REDUCE THE PRINCIPAL ERROR
! NORM.  IF WE JUST TRIED SLP WE QUIT, AND IF WE JUST TRIED RK WE ATTEMPT
! AN SLP STEP UNLESS IOPTHO = 2, IN WHICH CASE WE QUIT.
         if ( irk>0 ) then
            if ( ioptho/=2 ) then
               irk = -1
               goto 300
            end if
         end if
!
! IF IPHSE=-2 HERE WE WILL SET ITER=-2 AS A WARNING AND CHANGE
! ERROR(NUMGR+1) AND ERROR(NUMGR+3) BEFORE RETURNING.  OTHERWISE WE WILL
! HAVE IPHSE=0 AND WE WILL ADJUST ITER BEFORE RETURNING.
         if ( iphse<0 ) goto 400
         Iter = Iter + Itlim - itlim1
         goto 500
      end if

    end subroutine conmax
!********************************************************************************

!********************************************************************************
!>
!  This function subprogram returns the subscript of the first element of
!  array `iarr` relative to `iwork` (if the array is integer, i.e. `13 <=
!  iarr <= 23`) or relative to work (if the array is floating point, i.e.
!  `1 <= iarr <= 12` or `24 <= iarr <= 48`).

    pure function iloc(Iarr,Nparm,Numgr)

    implicit none

    integer,intent(in) :: Iarr
    integer,intent(in) :: Nparm
    integer,intent(in) :: Numgr
    integer :: iloc

    select case (Iarr)
    case (2)
        ! 2  ACTDIF(NUMGR)
        iloc = 1
    case (3)
        ! 3  B(NPARM+1)  (OPPOSITE V, Y;  FOLLOWS AA)
        iloc = Nparm**2 + 3*Numgr*Nparm + 6*Numgr + 13*Nparm + 9
    case (4)
        ! 4  BETA(NPARM+1)  (OPPOSITE V, Y;  FOLLOWS B)
        iloc = Nparm**2 + 3*Numgr*Nparm + 6*Numgr + 14*Nparm + 10
    case (5)
        ! 5  BNDKP(NPARM) (FOLLOWS ACTDIF)
        iloc = Numgr + 1
    case (6)
        ! 6  COEF(NUMGR)
        iloc = Numgr + Nparm + 1
    case (7)
        ! 7  COFBND(NPARM)
        iloc = 2*Numgr + Nparm + 1
    case (8)
        ! 8  CONFUN(NUMGR,NPARM+1)  (OPPOSITE PMAT1)
        iloc = 2*Numgr + 2*Nparm + 1
    case (9)
        ! 9  D(NPARM+1)  (OPPOSITE V, Y;  FOLLOWS BETA)
        iloc = Nparm**2 + 3*Numgr*Nparm + 6*Numgr + 15*Nparm + 11
    case (10)
        !  10  DVEC(NPARM) (FOLLOWS CONFUN)
        iloc = Numgr*Nparm + 3*Numgr + 2*Nparm + 1
    case (11)
        !  11  ERR1(NUMGR+3)
        iloc = Numgr*Nparm + 3*Numgr + 3*Nparm + 1
    case (12)
        !  12  FUNTBL(NUMGR,NPARM+1)
        iloc = Numgr*Nparm + 4*Numgr + 3*Nparm + 4
    case (13)
        !  13  IACT(NUMGR)
        iloc = 1
    case (14)
        !  14  IACT1(NUMGR)
        iloc = Numgr + 1
    case (15)
        !  15  ICOR(NPARM+1)
        iloc = 2*Numgr + 1
    case (16)
        !  16  ICOR1(NPARM+1)  (DOES NOT APPEAR IN PROGRAM BY NAME)
        iloc = 2*Numgr + Nparm + 2
    case (17)
        !  17  ICNTYP(NUMGR)
        iloc = 2*Numgr + 2*Nparm + 3
    case (18)
        !  18  IXRCT(NUMGR+2*NPARM)
        iloc = 3*Numgr + 2*Nparm + 3
    case (19)
        !  19  IYCCT(NPARM+1) (OPPOSITE KPIVOT)
        iloc = 4*Numgr + 4*Nparm + 3
    case (20)
        !  20  IYRCT(NUMGR+2*NPARM)
        iloc = 4*Numgr + 5*Nparm + 4
    case (21)
        !  21  JCNTYP(NUMGR)
        iloc = 5*Numgr + 7*Nparm + 4
    case (22)
        !  22  KCNTYP(NUMGR)
        iloc = 6*Numgr + 7*Nparm + 4
    case (23)
        !  23  KPIVOT(NPARM+1)  (OPPOSITE IYCCT)
        iloc = 4*Numgr + 4*Nparm + 3
    case (24)
        !  24  PARAM1(NPARM) (FOLLOWS FUNTBL)
        iloc = 2*Numgr*Nparm + 5*Numgr + 3*Nparm + 4
    case (25)
        !  25  PARPRJ(NPARM)
        iloc = 2*Numgr*Nparm + 5*Numgr + 4*Nparm + 4
    case (26)
        !  26  PARSER(NPARM)
        iloc = 2*Numgr*Nparm + 5*Numgr + 5*Nparm + 4
    case (27)
        !  27  PARWRK(NPARM)
        iloc = 2*Numgr*Nparm + 5*Numgr + 6*Nparm + 4
    case (28)
        !  28  PICOR(NPARM+1,NPARM+1)  (OPPOSITE V, Y;  FOLLOWS D)
        iloc = Nparm**2 + 3*Numgr*Nparm + 6*Numgr + 16*Nparm + 12
    case (29)
        !  29  PMAT(NPARM+1,NUMGR) (FOLLOWS PARWRK)
        iloc = 2*Numgr*Nparm + 5*Numgr + 7*Nparm + 4
    case (30)
        !  30  PMAT1(NPARM+1,NUMGR)  (OPPOSITE CONFUN)
        iloc = 2*Numgr + 2*Nparm + 1
    case (31)
        !  31  PTNR(NPARM+1) (FOLLOWS PMAT)
        iloc = 3*Numgr*Nparm + 6*Numgr + 7*Nparm + 4
    case (32)
        !  32  PTNRR(NPARM+1)
        iloc = 3*Numgr*Nparm + 6*Numgr + 8*Nparm + 5
    case (33)
        !  33  R(NPARM+1)
        iloc = 3*Numgr*Nparm + 6*Numgr + 9*Nparm + 6
    case (34)
        !  34  SAVE(NPARM+1)
        iloc = 3*Numgr*Nparm + 6*Numgr + 10*Nparm + 7
    case (35)
        !  35  V(NUMGR+2*NPARM+1,NPARM+2)  (WITH Y, OPPOSITE AA, B, BETA, D,
        !    PICOR, ZWORK)
        iloc = 3*Numgr*Nparm + 6*Numgr + 11*Nparm + 8
    case (36)
        !  36  VDER(NPARM) (FOLLOWS Y)
        iloc = 2*Nparm**2 + 4*Numgr*Nparm + 9*Numgr + 18*Nparm + 10
    case (37)
        !  37  VDERN(NPARM)
        iloc = 2*Nparm**2 + 4*Numgr*Nparm + 9*Numgr + 19*Nparm + 10
    case (38)
        !  38  VDERS(NPARM)
        iloc = 2*Nparm**2 + 4*Numgr*Nparm + 9*Numgr + 20*Nparm + 10
    case (39)
        !  39  VEC(NPARM+1)
        iloc = 2*Nparm**2 + 4*Numgr*Nparm + 9*Numgr + 21*Nparm + 10
    case (40)
        !  40  WCOEF(NUMGR)
        iloc = 2*Nparm**2 + 4*Numgr*Nparm + 9*Numgr + 22*Nparm + 11
    case (41)
        !  41  WCOEF1(NUMGR)  (DOES NOT APPEAR IN THE PROGRAM BY NAME)
        iloc = 2*Nparm**2 + 4*Numgr*Nparm + 10*Numgr + 22*Nparm + 11
    case (42)
        !  42  WPT(NPARM)
        iloc = 2*Nparm**2 + 4*Numgr*Nparm + 11*Numgr + 22*Nparm + 11
    case (43)
        !  43  WVEC(NPARM)
        iloc = 2*Nparm**2 + 4*Numgr*Nparm + 11*Numgr + 23*Nparm + 11
    case (44)
        !  44  X(NPARM+1)
        iloc = 2*Nparm**2 + 4*Numgr*Nparm + 11*Numgr + 24*Nparm + 11
    case (45)
        !  45  XKEEP(NPARM+1)
        iloc = 2*Nparm**2 + 4*Numgr*Nparm + 11*Numgr + 25*Nparm + 12
    case (46)
        !  46  XRK(NPARM+1)
        iloc = 2*Nparm**2 + 4*Numgr*Nparm + 11*Numgr + 26*Nparm + 13
    case (47)
        !  47  Y(NUMGR+2*NPARM)  (WITH V, OPPOSITE AA, B, BETA, D, PICOR, ZWORK;  FOLLOWS V)
        iloc = 2*Nparm**2 + 4*Numgr*Nparm + 8*Numgr + 16*Nparm + 10
    case (48)
        !  48  ZWORK(NPARM)  (OPPOSITE V, Y;  FOLLOWS PICOR)
        iloc = 2*Nparm**2 + 3*Numgr*Nparm + 6*Numgr + 18*Nparm + 13
    case default
        ! 1  AA(NPARM+1,NPARM+1)  (OPPOSITE V, Y; STARTS AT V STARTING POINT)
        iloc = 3*Numgr*Nparm + 6*Numgr + 11*Nparm + 8
    end select

    end function iloc
!********************************************************************************

!********************************************************************************
!>
! This subroutine uses `fnset` to compute `confun(i,1)` and the partial
! derivatives of the function whose value is in `confun(i,1)` for
! certain value(s) of `i`.  note that we do not want the `icntyp` computed
! by `fnset` to override the `icntyp` (or `jcntyp`) carried into this
! subroutine in `icntyp`, so we use `kcntyp` when we call `fnset`.  (the
! `icntyp` computed by `fnset` was stored earlier through a call to [[ercmp1]]
! from [[conmax]].)

    subroutine derst(me,Ioptn,Nparm,Numgr,Pttbl,Iptb,Indm,Param,&
                     Ipt,Param1,v,Kcntyp,Confun)

    implicit none

    class(conmax_solver),intent(inout) :: me
    integer,intent(in)  :: Ioptn
    integer  :: Nparm
    integer  :: Numgr
    integer  :: Iptb
    integer  :: Indm
    integer  :: Ipt
    real(wp),dimension(Iptb,Indm)               :: Pttbl
    real(wp),dimension(Nparm)                   :: Param1
    real(wp),dimension(Nparm)                   :: Param
    real(wp),dimension(Numgr+2*Nparm+1,Nparm+2) :: v
    integer ,dimension(Numgr)                   :: Kcntyp
    real(wp),dimension(Numgr,Nparm+1)           :: Confun

    real(wp) :: delt , delt2, up
    integer :: iopone , ioptth , iptkp , j , k , l , npar1

    ! IF THE ONES DIGIT OF IOPTN IS 0, WE CALL FNSET WITH INDFN=1 TO DO THE
    ! COMPUTATIONS DIRECTLY USING FORMULAS SUPPLIED BY THE USER.
    iopone = Ioptn - (Ioptn/10)*10
    if ( iopone<=0 ) then
        call me%fnset(Nparm,Numgr,Pttbl,Iptb,Indm,Param,Ipt,1,Kcntyp, &
        Confun)
        return
    else
        ! HERE THE ONES DIGIT OF IOPTN IS 1, AND WE APPROXIMATE THE PARTIAL
        ! DERIVATIVES USING CENTERED DIFFERENCE APPROXIMATIONS.
        ioptth = (Ioptn-(Ioptn/100000)*100000)/10000
        ! SET PRECISION DEPENDENT CONSTANTS.
        delt = sqrt(spcmn)
        delt2 = delt + delt
        if ( ioptth<=0 ) then
            ! HERE IOPONE=1 AND IOPTTH=0, AND WE WORK ONLY WITH CONSTRAINT IPT,
            ! WHERE IPT WILL BE AN INTEGER BETWEEN 1 AND NUMGR.
            ! L WILL BE THE INDEX OF THE VARIABLE WITH RESPECT TO WHICH WE ARE
            ! COMPUTING THE PARTIAL DERIVATIVE.
            do l = 1 , Nparm
                ! SET PARAM1 EQUAL TO PARAM, ECXEPT WITH ITS LTH COMPONENT INCREASED
                ! BY DELT.
                do j = 1 , Nparm
                    Param1(j) = Param(j)
                end do
                Param1(l) = Param(l) + delt
                ! NOW CALL FNSET WITH INDFN=0 TO PLACE THE FUNCTION IN CONSTRAINT
                ! IPT EVALUATED AT POINT PARAM1 IN CONFUN(IPT,1).
                call me%fnset(Nparm,Numgr,Pttbl,Iptb,Indm,Param1,Ipt,0,Kcntyp,Confun)
                up = Confun(Ipt,1)
                ! SET PARAM1 EQUAL TO PARAM, ECXEPT WITH ITS LTH COMOPONENT DECREASED
                ! BY DELT, AND CALL FNSET AGAIN.
                Param1(l) = Param(l) - delt
                call me%fnset(Nparm,Numgr,Pttbl,Iptb,Indm,Param1,Ipt,0,Kcntyp,Confun)
                ! NOW WE CAN COMPUTE THE CENTERED-DIFFERENCE APPROXIMATION TO THE PARTIAL
                ! DERIVATIVE OF THE FUNCTION IN CONSTRAINT IPT WITH RESPECT TO THE LTH
                ! VARIABLE AT THE POINT PARAM.  THIS BELONGS IN CONFUN(IPT,L+1), AND
                ! WE COULD PUT IT THERE NOW IF THE USER FOLLOWED DIRECTIONS AND DID NOT
                ! CHANGE CONFUN(IPT,L+1) (SINCE INDFN=0) IN LATER FNSET CALLS, BUT TO
                ! BE SAFE WE TEMPORARILY STORE IT IN V(L,1).
                ! NOTE THAT V IS USED ELSEWHERE IN THE PROGRAM, BUT HERE IT IS JUST A
                ! WORK ARRAY, WHILE THE WORK ARRAY PARAM1 IS NOT USED ELSEWHERE IN
                ! THE PROGRAM.
                v(l,1) = (up-Confun(Ipt,1))/delt2
            end do
            ! NOW COMPUTE THE VALUE OF THE FUNCTION AT PARAM, AND THEN PUT THE
            ! EARLIER-COMPUTED PARTIAL DERIVATIVES INTO CONFUN.
            call me%fnset(Nparm,Numgr,Pttbl,Iptb,Indm,Param,Ipt,0,Kcntyp, &
            Confun)
            do l = 1 , Nparm
                Confun(Ipt,l+1) = v(l,1)
            end do
            return
        else
            ! HERE IOPONE=1 AND IOPTTH=1, AND EACH TIME FNSET IS CALLED IT WILL
            ! COMPUTE VALUES FOR THE FUNCTIONS IN THE LEFT SIDES OF ALL CONSTRAINTS
            ! (EXCEPT THOSE WHERE FNSET SETS ICNTYP(I)=0) IF IPT=0, AND WILL COMPUTE
            ! VALUES FOR THE FUNCTIONS IN THE LEFT SIDES OF ALL STANDARD (I.E. TYPE
            ! -1 OR -2) CONSTRAINTS IF IPT=-1.
            ! WE FIRST SAVE IPT IN CASE THE USER CHANGES IT IN A FNSET CALL;  WE WILL
            ! RESTORE IT AFTER EACH FNSET CALL.
            iptkp = Ipt
            npar1 = Nparm + 1
            ! WE WILL COMPUTE APPROXIMATIONS TO PARTIAL DERIVATIVES FOR THOSE
            ! CONSTRAINTS WHICH FNSET IS ASKED BY IPT TO COMPUTE.  TO DETERMINE WHICH
            ! THESE ARE WE ZERO OUT KCNTYP;  AFTER A FNSET CALL, THE DESIRED
            ! CONSTRAINTS WILL BE THE CONSTRAINTS K WITH KCNTYP(K) /= 0 IF IPT=0,
            ! OR THE CONSTRAINTS K WITH KCNTYP(K) < 0 IF IPT=-1.
            do k = 1 , Numgr
                Kcntyp(k) = 0
            end do
            ! NOW FOLLOW BASICALLY THE SAME PROCEDURES AS IN THE IOPTTH=0 CASE DONE
            ! ABOVE.
            do l = 1 , Nparm
                do j = 1 , Nparm
                    Param1(j) = Param(j)
                end do
                Param1(l) = Param(l) + delt
                call me%fnset(Nparm,Numgr,Pttbl,Iptb,Indm,Param1,Ipt,0, &
                Kcntyp,Confun)
                Ipt = iptkp
                do k = 1 , Numgr
                    if ( Ipt<0 ) then
                        if ( Kcntyp(k)>=0 ) cycle
                    else if ( Kcntyp(k)==0 ) then
                        cycle
                    end if
                    ! SAVE THE UPPER NUMBERS IN COLUMN NPARM+1 OF V.
                    v(k,npar1) = Confun(k,1)
                end do
                ! REVISE PARAM1 AND CALL FNSET AGAIN.
                Param1(l) = Param(l) - delt
                call me%fnset(Nparm,Numgr,Pttbl,Iptb,Indm,Param1,Ipt,0, &
                Kcntyp,Confun)
                Ipt = iptkp
                do k = 1 , Numgr
                    if ( Ipt<0 ) then
                        if ( Kcntyp(k)>=0 ) cycle
                    else if ( Kcntyp(k)==0 ) then
                        cycle
                    end if
                    ! STORE THE APPROXIMATE PARTIAL DERIVATIVES WITH RESPECT TO THE LTH
                    ! VARIABLE IN THE LTH COLUMN OF V.
                    v(k,l) = (v(k,npar1)-Confun(k,1))/delt2
                end do
            end do
            ! CALL FNSET AGAIN TO COMPUTE THE VALUES OF THE FUNCTIONS AT POINT
            ! PARAM, AND THEN PUT THE EARLIER-COMPUTED PARTIAL DERIVATIVES INTO
            ! CONFUN.
            call me%fnset(Nparm,Numgr,Pttbl,Iptb,Indm,Param,Ipt,0,Kcntyp, &
            Confun)
            do k = 1 , Numgr
                if ( Ipt<0 ) then
                    if ( Kcntyp(k)>=0 ) cycle
                else if ( Kcntyp(k)==0 ) then
                    cycle
                end if
                do l = 1 , Nparm
                    Confun(k,l+1) = v(k,l)
                end do
            end do
        end if
    end if

    end subroutine derst
!********************************************************************************

!********************************************************************************
!>
!
    subroutine slpcon(me,Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,  &
                      Tolcon,Rchin,Irk,Itypm1,Itypm2,Icntyp,Rchdwn,   &
                      Numlim,Itersl,Prjslp,Funtbl,Iyrct,x,Mact1,Iact1,&
                      Jcntyp,Iphse,Enchg,Iwork,Liwrk,Work,Lwrk,Parser,&
                      Isucc,Param,Error)

    implicit none

    class(conmax_solver),intent(inout) :: me
    integer  :: Ioptn
    integer  :: Nparm
    integer  :: Numgr
    integer  :: Ifun
    integer  :: Iptb
    integer  :: Indm
    integer  :: Irk
    integer  :: Itypm1
    integer  :: Itypm2
    integer  :: Numlim
    integer  :: Itersl
    integer  :: Mact1
    integer  :: Iphse
    integer  :: Liwrk
    integer  :: Lwrk
    integer  :: Isucc
    real(wp) :: Tolcon
    real(wp) :: Rchin
    real(wp) :: Rchdwn
    real(wp) :: Prjslp
    real(wp) :: Enchg
    real(wp) :: Fun(Ifun)
    real(wp) :: Pttbl(Iptb,Indm)
    integer  :: Icntyp(Numgr)
    integer  :: Iyrct(Numgr+2*Nparm)
    integer  :: Iact1(Numgr)
    integer  :: Jcntyp(Numgr)
    integer  :: Iwork(Liwrk)
    real(wp) :: Funtbl(Numgr,Nparm+1)
    real(wp) :: x(Nparm+1)
    real(wp) :: Work(Lwrk)
    real(wp) :: Parser(Nparm)
    real(wp) :: Param(Nparm)
    real(wp) :: Error(Numgr+3)

    real(wp) :: bndlgt , emin , emin1 , enorm , prjlim , quots , &
                ss , tol1 , tol2 , unit
    integer :: i , ilc05 , ilc07 , ilc08 , ilc11 , ilc13 , ilc18 , &
               ilc19 , ilc25 , ilc35 , ilc45 , ilc47 , indic , &
               ipmax , ismax , j , m , ng3 , npar1 , nsrch, numin

    ! SET MACHINE AND PRECISION DEPENDENT CONSTANTS.
    tol1 = ten*ten*spcmn
    tol2 = ten*spcmn
    ilc05 = iloc(5,Nparm,Numgr)
    ilc07 = iloc(7,Nparm,Numgr)
    ilc08 = iloc(8,Nparm,Numgr)
    ilc11 = iloc(11,Nparm,Numgr)
    ilc13 = iloc(13,Nparm,Numgr)
    ilc18 = iloc(18,Nparm,Numgr)
    ilc19 = iloc(19,Nparm,Numgr)
    ilc25 = iloc(25,Nparm,Numgr)
    ilc35 = iloc(35,Nparm,Numgr)
    ilc45 = iloc(45,Nparm,Numgr)
    ilc47 = iloc(47,Nparm,Numgr)
    numin = 0
    Isucc = 0
    enorm = Error(Numgr+1)
    npar1 = Nparm + 1
    ng3 = Numgr + 3

    ! IF ITERSL=0, SET IYRCT(1)=-1 FOR USE IN SETU1 AND TO TELL SLNPRO NOT
    ! TO TRY TO USE INFORMATION FROM A PREVIOUS VERTEX.
    if ( Itersl<=0 ) Iyrct(1) = -1

    do

        ! CALL BNDSET TO SET (OR RESET) THE COEFFICIENT CHANGE BOUNDS.
        call bndset(Nparm,x,Itersl,numin,Prjslp,Work(ilc07),Work(ilc45), &
                    Work(ilc05))

        ! CALL SETU1 TO SET UP FOR SLNPRO AND, IF NUMIN=0, TO DETERMINE
        ! WHICH CONSTRAINTS ARE ACTIVE AND STORE FUNCTION AND GRADIENT VALUES
        ! FOR THEM IN FUNTBL.
        call me%setu1(Ioptn,Numgr,Nparm,numin,Rchin,Pttbl,Iptb,Indm,Fun,Ifun,&
                    Funtbl,Work(ilc07),Param,Icntyp,Rchdwn,Error,Mact1,      &
                    Iact1,bndlgt,Iyrct,Iphse,Iwork,Liwrk,Work,Lwrk,          &
                    Work(ilc08),Iwork(ilc13),Work(ilc35),m)

        ! SET UNIT (FOR USE IN RCHMOD) EQUAL TO THE VALUE OF BNDLGT AFTER
        ! SETU1 IS CALLED WITH NUMIN=0.
        if ( numin<=0 ) unit = bndlgt

        ! CALL SLNPRO TO COMPUTE A SEARCH DIRECTION X.
        call slnpro(Work(ilc35),m,npar1,Iyrct,Work(ilc47),Iwork(ilc18), &
                    Iwork(ilc19),Nparm,Numgr,x,indic)

        ! IF INDIC > 0 THEN SLNPRO FAILED TO PRODUCE AN X, AND IF WE HAVE
        ! REACHED THE SLPCON ITERATION LIMIT WE RETURN WITH THE WARNING
        ! ISUCC=1.
        if ( indic<=0 ) then

            ! HERE SLNPRO SUCCEEDED AND WE SET PRJSLP=1.0 INITIALLY FOR SEARSL.
            Prjslp = one

            ! WE NOW WISH TO DETERMINE PRJLIM = THE SMALLER OF 1.0/SPCMN AND
            ! THE LARGEST VALUE OF PRJSLP FOR WHICH THE LINEAR STANDARD CONSTRAINTS
            ! ARE SATISFIED FOR THE PARAMETER VECTOR PARAM+PRJSLP*X.  THIS
            ! WILL GIVE AN UPPER BOUND FOR LINE SEARCHING.  NOTE THAT IN
            ! THEORY WE SHOULD HAVE PRJLIM >= 1.0 SINCE THE LINEAR STANDARD
            ! CONSTRAINTS SHOULD BE SATISFIED FOR PRJSLP=0.0 AND PRJSLP=1.0, BUT
            ! ROUNDOFF ERROR COULD AFFECT THIS A LITTLE.  IF THERE ARE NO
            ! LINEAR STANDARD CONSTRAINTS, WE SET PRJLIM=1.0/SPCMN.
            prjlim = big
            !*****INSERT TO MAKE SEARCHING LESS VIOLENT.
            !     PRJLIM=TWO
            !*****END INSERT
            if ( Itypm1>0 ) then
                do i = 1 , Numgr
                    if ( Icntyp(i)+1==0 ) then
                        ! WE WISH TO HAVE SUMMATION (FUNTBL(I,J+1)*(PARAM(J)+PRJSLP*X(J)))
                        ! + C(I) <= 0.0 FOR I=1,...,NUMGR, ICNTYP(I) = -1,
                        ! WHERE THE ITH CONSTRAINT APPLIED TO PARAM SAYS
                        ! SUMMATION (FUNTBL(I,J+1)*PARAM(J)) + C(I) <= 0.0, SO C(I) IS THE
                        ! CONSTANT TERM ON THE LEFT SIDE OF LINEAR CONSTRANT I.
                        ! THUS FOR I=1,...,NUMGR, ICNTYP(I) = -1, WE WANT PRJLIM*SS <= SSS,
                        ! WHERE SS = SUMMATION (FUNTBL(I,J+1)*X(J)) AND SSS = -C(I) -
                        ! SUMMATION (FUNTBL(I,J+1)*PARAM(J)) = -FUNTBL(I,1).
                        ss = zero
                        do j = 1 , Nparm
                            ss = ss + Funtbl(i,j+1)*x(j)
                        end do
                        ! IF SS < 10.0*SPCMN THIS CONSTRAINT WILL NOT PUT A SIGNIFICANT
                        ! RESTRICTION ON PRJSLP.
                        if ( ss>=tol2 ) then
                            ! HERE SS >= 10.0*SPCMN AND WE COMPARE SSS/SS AGIANST PRJLIM.
                            quots = -Funtbl(i,1)/ss
                            if ( prjlim>quots ) prjlim = quots
                        end if
                    end if
                end do
            end if
            ! DO NOT ALLOW A PRJSLP SMALLER THAN TOL1.
            if ( Prjslp<tol1 ) Prjslp = tol1
            ! CALL SEARSL TO DO A LINE SEARCH IN DIRECTION X.
            call me%searsl(Ioptn,Numgr,Nparm,prjlim,tol1,x,Fun,Ifun,Pttbl, &
                        Iptb,Indm,Param,Error,Rchdwn,Mact1,Iact1,Iphse,    &
                        unit,Tolcon,Rchin,Itypm1,Itypm2,Iwork,Liwrk,Work,  &
                        Lwrk,Work(ilc11),Work(ilc25),Prjslp,emin,emin1,    &
                        Parser,nsrch)

            ! COMPUTE THE ERROR NORM CHANGE ENCHG.
            Enchg = emin - enorm

            ! IF WE HAVE AN IMPROVEMENT IN THE ERROR NORM ENORM OF MORE THAN TOL1
            ! WE UPDATE PARAM AND ERROR AND RETURN WITH ISUCC=0, INDICATING SUCCESS.
            ! OTHERWISE WE CHECK TO SEE IF WE HAVE REACHED THE SLPCON ITERATION
            ! LIMIT, AND IF SO WE RETURN WITH ISUCC=1, INDICATING FAILURE.
            if ( Enchg+tol1<0 ) then
                ! HERE WE HAD AN IMPROVEMENT IN THE ERROR NORM ENORM OF MORE THAN TOL1.
                do j = 1 , Nparm
                    Param(j) = Parser(j)
                end do
                call me%ercmp1(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,   &
                               Param,1,Iphse,Iwork,Liwrk,Work(ilc08),Icntyp, &
                               ipmax,ismax,Error)
                return
            end if
        end if

        ! HERE WE DID NOT OBTAIN AN IMPROVED ERROR NORM SO WE RETURN WITH THE
        ! WARNING ISUCC=1 IF WE HAVE DONE NUMLIN ITERATIONS IN SLPCON.
        if ( numin<Numlim ) then
            ! HERE WE DID NOT OBTAIN AN IMPROVED ERROR NORM BUT WE HAVE NOT YET DONE
            ! NUMLIM ITERATIONS IN SLPCON SO WE INCREMENT NUMIN, SET IYRCT(1)=-1 TO
            ! TELL SLNPRO NOT TO TRY TO USE INFORMATION FROM THE PREVIOUS FAILED
            ! VERTEX, AND GO BACK TO CALL BNDSET AND TRY ANOTHER ITERATION WITH
            ! A DIFFERENT TRUST REGION.
            numin = numin + 1
            Iyrct(1) = -1
            cycle
        end if
        exit
    end do
    Isucc = 1

    end subroutine slpcon
!********************************************************************************

!********************************************************************************
!>
!  This subroutine sets the bounds on the coefficient changes in
!  [[slnpro]].

    subroutine bndset(Nparm,x,Itersl,Numin,Prjslp,Cofbnd,Xkeep,Bndkp)

    implicit none

    integer,intent(in)                        :: Nparm
    real(wp),dimension(Nparm+1),intent(in)    :: x
    integer,intent(in)                        :: Itersl
    integer,intent(in)                        :: Numin
    real(wp),intent(in)                       :: Prjslp
    real(wp),dimension(Nparm),intent(inout)   :: Cofbnd
    real(wp),dimension(Nparm+1),intent(inout) :: Xkeep
    real(wp),dimension(Nparm),intent(inout)   :: Bndkp

    ! set initial parameters.  fact1, fact3a, fact3b, chlm1, and chlm2
    ! should be between 0.0 and 1.0, while fact2 should be > 1.0.
    real(wp),parameter :: fact1 = (one+two)/four
    real(wp),parameter :: fact2 = two
    real(wp),parameter :: fact3a = one/ten
    real(wp),parameter :: fact3b = one/(ten*ten)
    real(wp),parameter :: fact4 = two/ten
    real(wp),parameter :: chlm1 = one/ten
    real(wp),parameter :: chlm2 = (four+four)/ten
    real(wp),parameter :: tstprj = one/two - one/(ten*ten*ten)
    real(wp),parameter :: epsil = ten*ten*spcmn
    real(wp),parameter :: epsil1 = (one+one/(ten*ten*ten))*epsil
    real(wp),parameter :: bnd = two/(ten*ten) !! The initial bound on all coefficient changes.

    real(wp) :: bsave , fact3
    integer :: itight , j

    if ( Numin<1 ) then

        if ( Itersl<1 ) then

        else if ( Itersl==1 ) then
            ! HERE NUMIN=0 AND ITERSL=1, SO THE LAST BNDSET CALL RESULTED IN
            ! THE FIRST SUCCESSFUL PRINCIPAL ERROR NORM IMPROVEMENT,
            ! AND SO WE SAVE COFBND IN BNDKP AND X IN XKEEP.  WE WILL NOT
            ! CHANGE COFBND HERE.
            do j = 1 , Nparm
                Xkeep(j) = x(j)
                Bndkp(j) = Cofbnd(j)
            end do
            return
        else
            ! HERE NUMIN=0 AND ITERSL >= 2, SO WE HAVE HAD AT LEAST 2 SUCCESSES,
            ! WITH THE COEFFICIENTS AND BOUNDS FOR THE LAST ONE IN X AND
            ! COFBND RESPECTIVELY, AND THE COEFFICIENTS AND BOUNDS FOR THE
            ! PREVIOUS ONE IN XKEEP AND BNDKP RESPECTIVELY.  WE WILL FORM A
            ! NEW COFBND, AND SHIFT THE OLD COFBND INTO BNDKP AND X INTO XKEEP.
            do j = 1 , Nparm
                Cofbnd_block : block
                    ! SAVE THE OLD COFBND(J) IN BSAVE.
                    bsave = Cofbnd(j)
                    ! IF AT BOTH THE LAST AND PREVIOUS SUCCESSFUL ITERATION THE CHANGES
                    ! IN A COEFFICIENT RELATIVE TO ITS BOUND WERE >= CHLM2 IN ABSOLUTE
                    ! VALUE AND IN THE SAME DIRECTION, WE LOOSEN THE BOUND BY A FACTOR
                    ! OF FACT2.  IF THE RELATIVE CHANGES WERE >= CHLM1 IN ABSOLUTE
                    ! VALUE AND IN OPPOSITE DIRECTIONS, WE TIGHTEN THE BOUND BY A FACTOR
                    ! OF FACT1 BECAUSE OF SUSPECTED OSCILLATION.  WE ALSO TIGHTEN THE
                    ! BOUND IF BOTH RELATIVE CHANGES WERE LESS THAN CHLM1 IN ABSOLUTE
                    ! VALUE IN ORDER TO PREVENT A LONG SEQUENCE OF OSCILLATIONS OF THE
                    ! SAME SMALL ORDER.  OTHERWISE WE LEAVE THE BOUND ALONE.
                    ! THE NEXT FOUR IF STATEMENTS CHECK TO SEE IF THE BOUND SHOULD BE
                    ! LOOSENED.
                    if ( x(j)<chlm2*Cofbnd(j) ) then
                        if ( x(j)+chlm2*Cofbnd(j)<=0 ) then
                            if ( Xkeep(j)+chlm2*Bndkp(j)<=0 ) then
                                ! LOOSEN THE BOUND.
                                Cofbnd(j) = fact2*Cofbnd(j)
                                exit Cofbnd_block
                            end if
                        end if
                    else if ( Xkeep(j)>=chlm2*Bndkp(j) ) then
                        Cofbnd(j) = fact2*Cofbnd(j)
                        exit Cofbnd_block
                    end if
                    ! HERE THE BOUND SHOULD NOT BE LOOSENED.  THE NEXT FIVE IF
                    ! STATEMTENTS CHECK TO SEE IF IT SHOULD BE TIGHTENED.
                    if ( x(j)<chlm1*Cofbnd(j) ) then
                        if ( x(j)+chlm1*Cofbnd(j)<=0 ) then
                            if ( Xkeep(j)<chlm1*Bndkp(j) ) exit Cofbnd_block
                            ! HERE WE HAVE ABS(X(J)) < CHLM1*COFBND(J).
                        else if ( abs(Xkeep(j))>=chlm1*Bndkp(j) ) then
                            exit Cofbnd_block
                        end if
                    else if ( Xkeep(j)+chlm1*Bndkp(j)>0 ) then
                        exit Cofbnd_block
                    end if
                    ! TIGHTEN THE BOUND.
                    Cofbnd(j) = fact1*Cofbnd(j)
                    ! DO NOT ALLOW THE BOUND TO DROP BELOW EPSIL.
                    if ( Cofbnd(j)<epsil ) Cofbnd(j) = epsil
                end block Cofbnd_block
                ! SAVE X(J) AND THE OLD COFBND(J).
                Bndkp(j) = bsave
                Xkeep(j) = x(j)
            end do
            ! IF THE LAST PROJECTION FACTOR IS SMALLER THAN .499, WE TIGHTEN THE
            ! BOUNDS BY A FACTOR OF 0.2, WITH THE RESTRICTION THAT WE DO NOT
            ! ALLOW THE BOUNDS TO DROP BELOW EPSIL.
            if ( Prjslp<tstprj ) then
                do j = 1 , Nparm
                    Cofbnd(j) = fact4*Cofbnd(j)
                    if ( Cofbnd(j)<epsil ) Cofbnd(j) = epsil
                end do
            end if
            return
        end if

        ! HERE NUMIN=0 AND ITERSL=0, SO WE ARE IN THE FIRST BNDSET CALL SINCE THE
        ! LAST RK SUCCESS (IF ANY), SO WE SET INITIAL BOUNDS.
        Cofbnd = bnd
        return

    else if ( Numin==1 ) then
        ! HERE NUMIN=1 SO THE LAST BNDSET CALL RESULTED IN FAILURE TO
        ! IMPROVE THE PRINCIPAL ERROR NORM, AND WE SET FACT3=
        ! FACT3A AND TIGHTEN THE BOUNDS.
        fact3 = fact3a
    else
        ! HERE NUMIN > 1 SO THERE HAVE BEEN AT LEAST 2 SUCCESSIVE
        ! FAILURES, AND WE SET FACT3=FACT3B AND TIGHTEN THE BOUNDS.
        fact3 = fact3b
    end if

    ! TIGHTEN THE BOUNDS BY A FACTOR OF FACT3.
    itight = 1
    do j = 1 , Nparm
        bsave = Cofbnd(j)
        Cofbnd(j) = fact3*bsave
        ! WE DO NOT ALLOW A BOUND TO DROP BELOW EPSIL.
        if ( Cofbnd(j)<epsil ) then
            ! IF THE BOUND WAS ALREADY (ESSENTIALLY) AT EPSIL, KEEP TRACK OF
            ! THIS BY NOT SETTING ITIGHT=0.
            if ( bsave>epsil1 ) itight = 0
            Cofbnd(j) = epsil
        else
            itight = 0
        end if
    end do

    ! IF ALL THE BOUNDS WERE ALREADY (ESSENTIALLY) AT EPSIL, WE TRY
    ! RESETTING THE BOUNDS TO THEIR ORIGINAL VALUES.
    !2800 WRITE(NWRIT,2900)
    !2900 FORMAT(/52H *****RESETTING BOUNDS TO THEIR ORIGINAL VALUES*****)
    if ( itight>0 ) Cofbnd = bnd

    end subroutine bndset
!********************************************************************************

!********************************************************************************
!>
! THIS SUBROUTINE SETS UP V FOR SLNPRO TO SOLVE A MODIFIED LINEARIZED
! (ABOUT THE OLD PARAMETERS IN PARAM) VERSION OF OUR PROBLEM.

    subroutine setu1(me,Ioptn,Numgr,Nparm,Numin,Rchin,Pttbl,Iptb,Indm,&
                     Fun,Ifun,Funtbl,Cofbnd,Param,Icntyp,Rchdwn,Error,&
                     Mact1,Iact1,Bndlgt,Iyrct,Iphse,Iwork,Liwrk,Work, &
                     Lwrk,Confun,Iact,v,m)

      implicit none

      class(conmax_solver),intent(inout) :: me
      real(wp) :: actlim , bndfud , Bndlgt , Cofbnd , Confun , enorm , &
                  Error , Fun , Funtbl , grdlgt , Param , &
                  Pttbl , Rchdwn , Rchin , rchind , rt , stfudg , sum
      real(wp) :: v , Work
      integer :: i , Iact , Iact1 , Icntyp , Ifun , ii , ilc22 , ilc24 , &
                 Indm , Ioptn , ioptth , Iphse , ipt , Iptb , &
                 Iwork , Iyrct , j , jj , k
      integer :: kk , l , Liwrk , Lwrk , m , mact , Mact1 , mm1 , mp1 , &
                 npar1 , npar2 , Nparm , Numgr , Numin

      dimension Pttbl(Iptb,Indm) , Fun(Ifun) , Funtbl(Numgr,Nparm+1) ,   &
                Cofbnd(Nparm) , Param(Nparm) , Error(Numgr+3) ,          &
                v(Numgr+2*Nparm+1,Nparm+2) , Iact(Numgr) , Iact1(Numgr), &
                Iyrct(Numgr+2*Nparm) , Icntyp(Numgr) ,                   &
                Confun(Numgr,Nparm+1) , Iwork(Liwrk) , Work(Lwrk)

      ilc22 = iloc(22,Nparm,Numgr)
      ilc24 = iloc(24,Nparm,Numgr)
      npar1 = Nparm + 1
      npar2 = Nparm + 2
      ioptth = (Ioptn-(Ioptn/100000)*100000)/10000
!
! THE LINEARIZED PROBLEM REPLACES THE APPROXIMATING FUNCTION BY ITS
! FIRST ORDER TAYLOR SERIES, SO FUN(I)-(APPROXIMATING FUNCTION)(I) IS
! REPLACED BY ERROR(I)-(SUMMATION OF COEFFICIENT CHANGES TIMES PARTIAL
! DERIVATIVES OF THE APPROXIMATING FUNCTION WITH RESPECT TO THOSE
! COEFFICIENTS) IF ICNTYP(I)=2, AND IF ICNTYP(I)=1 WE REPLACE THE LEFT
! SIDE OF CONSTRAINT I BY ERROR(I)+(SUMMATION OF COEFFICIENT CHANGES TIMES
! PARTIAL DERIVATIVES OF THE LEFT SIDE OF CONSTRAINT I).
! V AND M ARE THE OUTPUT QUANTITIES.  M WILL KEEP TRACK OF THE NUMBER
! OF CONSTRAINTS IN THE LP PROBLEM TO BE SOLVED BY SLNPRO.
      m = 0
      enorm = Error(Numgr+1)
      stfudg = one/ten
!
! COMPUTE THE LENGTH OF THE LONGEST X VECTOR SATISFYING THE COEFFICIENT
! CHANGE BOUNDS.
      sum = zero
      do j = 1 , Nparm
         sum = sum + (Cofbnd(j))**2
      end do
      Bndlgt = sqrt(sum)
      bndfud = stfudg*Bndlgt
!
! WE WILL SAY A PRIMARY CONSTRAINT IS ACTIVE IF ERROR(I) (OR ABS(ERROR(I
! IF ICNTYP(I)=2) >= ENORM-RCHDWN*BNDLGT.
      actlim = enorm - Rchdwn*Bndlgt
!
! WE WILL SAY A TYPE -2 CONSTRAINT IS ACTIVE IF ERROR(I) >= -RCHIND.
      rchind = Rchin*Bndlgt
!
      if ( Numin<=0 ) then
!
! HERE NUMIN=0, SO WE WILL FIRST COMPUTE A NEW SET OF ACTIVE INDICES,
! THEN PUT THE FUNCTION VALUES AND GRADIENTS FOR THESE INDICES IN
! FUNTBL, WHERE THEY WILL REMAIN THROUGHOUT THIS CALL TO SLPCON.
         do i = 1 , Numgr
            if ( Icntyp(i)<0 ) then
!
! HERE ICNTYP(I) < 0 AND WE WILL DECLARE THE CONSTRAINT TO BE ACTIVE IF
! AND ONLY IF ICNTYP(I)=-1, OR ICNTYP(I)=-2 AND ERROR(I) >= -RCHIND.
               if ( Icntyp(i)+1<0 ) then
                  if ( Error(i)+rchind<0 ) goto 50
               end if
            else if ( Icntyp(i)==0 ) then
               goto 50
            else if ( Icntyp(i)<=1 ) then
!
! HERE ICNTYP(I)=1 AND WE WILL DECLARE THE CONSTRAINT TO BE +ACTIVE IF AND
! ONLY IF ERROR(I) >= ACTLIM.
               if ( Error(i)<actlim ) goto 50
!
! HERE ICNTYP(I)=2 AND WE WILL DECLARE THE CONSTRAINT TO BE +ACTIVE IF AND
! ONLY IF ERROR(I) >= ACTLIM OR -ACTIVE IF AND ONLY IF ERROR(I)  <=
! -ACTLIM.
            else if ( Error(i)<actlim ) then
               if ( Error(i)+actlim<=0 ) then
!
! DECLARE CONSTRAINT I TO BE -ACTIVE.
                  m = m + 1
                  Iact(m) = -i
               end if
               goto 50
            end if
!
! DECLARE CONSTRAINT I TO BE (+)ACTIVE.
            m = m + 1
            Iact(m) = i
 50      end do
         mact = m
!
! NOW PUT ACTIVE VALUES AND GRADIENTS IN FUNTBL.
         if ( ioptth<=0 ) then
! HERE IOPTTH=0 AND WE CALL DERST FOR EACH ACTIVE CONSTRAINT.
            do l = 1 , mact
               i = abs(Iact(l))
               ipt = i
! CALL DERST TO COMPUTE BOTH FUNCTION AND GRADIENT VALUES.
               call me%derst(Ioptn,Nparm,Numgr,Pttbl,Iptb,Indm,Param,ipt,  &
                             Work(ilc24),v,Iwork(ilc22),Confun)
! COPY THE VALUES FOR CONSTRAINT I INTO FUNTBL.
               do j = 1 , npar1
                  Funtbl(i,j) = Confun(i,j)
               end do
            end do
            goto 300
         else
!
! HERE IOPTTH=1 AND ONLY ONE DERST CALL IS NEEDED.
! IF IPHSE < 0 OR NO ICNTYP(L) IS POSITIVE, SET IPT=-1 TO TELL DERST
! TO COMPUTE STANDARD CONSTRAINTS ONLY, WHILE OTHERWISE SET IPT=0 TO
! TELL DERST TO COMPUTE ALL CONSTRAINTS.
            if ( Iphse>=0 ) then
               do l = 1 , Numgr
                  if ( Icntyp(l)>0 ) goto 100
               end do
            end if
            ipt = -1
            goto 200
         end if
 100     ipt = 0
      else
! HERE NUMIN IS NOT 0, AND WE WILL KEEP THE OLD ACTIVE CONSTRAINT SET
! AND FOREGO RECOMPUTING FUNCTION VALUES AND GRADIENTS.
         mact = Mact1
         m = mact
         do l = 1 , mact
            Iact(l) = Iact1(l)
         end do
         goto 300
      end if
 200  call me%derst(Ioptn,Nparm,Numgr,Pttbl,Iptb,Indm,Param,ipt,Work(ilc24),&
                    v,Iwork(ilc22),Confun)
! COPY THE ACTIVE FUNCTION AND GRADIENT VALUES INTO FUNTBL.
      do l = 1 , mact
         i = abs(Iact(l))
         do j = 1 , npar1
            Funtbl(i,j) = Confun(i,j)
         end do
      end do
!
! NOW SET UP THE ACTIVE CONSTRAINTS IN V FOR SLNPRO.
 300  do l = 1 , mact
         i = abs(Iact(l))
         if ( Icntyp(i)<0 ) then
!
            if ( Icntyp(i)+1<0 ) then
!
! HERE ICNTYP(I)=-2 AND WE FIRST COMPUTE THE LENGTH OF THE GRADIENT.
               sum = zero
               do j = 1 , Nparm
                  sum = sum + (Funtbl(i,j+1))**2
               end do
               grdlgt = sqrt(sum)
! NOW SET UP A CONSTRAINT OF THE FORM GRADIENT.CHANGE  <=
! -MIN(1.0,CONSTRAINT VALUE)*BNDFUD*GRDLGT, SO IF GRDLGT > 0.0 WE
! HAVE (-GRADIENT/GRDLGT).(CHANGE/BNDLGT) >= STFUDG*MIN(1.0,
! CONSTRAINT VALUE).
               do j = 1 , Nparm
                  v(l,j) = Funtbl(i,j+1)
               end do
               v(l,npar1) = zero
               rt = Error(i)
               if ( rt>one ) rt = one
               v(l,npar2) = -rt*bndfud*grdlgt
            else
!
! HERE ICNTYP(I)=-1 AND WE SET UP A CONSTRAINT OF THE FORM
! GRADIENT.CHANGE <= -CONSTRAINT VALUE.
               do j = 1 , Nparm
                  v(l,j) = Funtbl(i,j+1)
               end do
               v(l,npar1) = zero
               v(l,npar2) = -Error(i)
            end if
         else if ( Icntyp(i)/=0 ) then
            if ( Icntyp(i)<=1 ) then
!
! HERE ICNTYP(I)=1 AND WE SET UP A CONSTRAINT OF THE FORM
! GRADIENT.CHANGE - W <= -CONSTRAINT VALUE.
               do j = 1 , Nparm
                  v(l,j) = Funtbl(i,j+1)
               end do
               v(l,npar1) = -one
               v(l,npar2) = -Error(i)
            else if ( Iact(l)<=0 ) then
!
! HERE ICNTYP(I)=2 AND IACT(L) < 0, AND WE SET UP A CONSTRAINT OF THE
! FORM GRADIENT.CHANGE - W <= FUN - CONSTRAINT VALUE.
               do j = 1 , Nparm
                  v(l,j) = Funtbl(i,j+1)
               end do
               v(l,npar1) = -one
               v(l,npar2) = Error(i)
            else
!
! HERE ICNTYP(I)=2 AND IACT(L) > 0, AND WE SET UP A CONSTRAINT OF THE
! FORM -GRADIENT.CHANGE - W <= -(FUN - CONSTRAINT VALUE).
               do j = 1 , Nparm
                  v(l,j) = -Funtbl(i,j+1)
               end do
               v(l,npar1) = -one
               v(l,npar2) = -Error(i)
            end if
         end if
      end do
!
! SET THE CONSTRAINTS OF THE FORM -X(J) <= COFBND(J) AND
! X(J) <= COFBND(J).
      do j = 1 , Nparm
         m = m + 2
         mm1 = m - 1
         do k = 1 , npar1
            v(mm1,k) = zero
            v(m,k) = zero
         end do
         v(mm1,j) = -one
         v(m,j) = one
         v(mm1,npar2) = Cofbnd(j)
         v(m,npar2) = Cofbnd(j)
      end do
!
! NOW SET THE BOTTOM ROW.  TO MINIMIZE W = X(NPARM+1) WE MAXIMIZE -W.
      mp1 = m + 1
      do j = 1 , npar2
         v(mp1,j) = zero
      end do
      v(mp1,npar1) = one
!
! THIS SECTION ADJUSTS IYRCT TO EITHER TELL SLNPRO TO DO THE INITIAL
! EXCHANGES STRICTLY ACCORDING TO A PIVOTING STRATEGY (BY SETTING
! IYRCT(1)=-1) OR TO SPECIFY AN INITIAL VERTEX FOR SLNPRO, NAMELY THE
! VERTEX CORRESPONDING TO THE LAST LINEAR PROGRAMMING SOLUTION.
! IF IYRCT(1) IS -1 ALREADY WE DO NOT ATTEMPT TO SPECIFY A VERTEX, BUT
! WE STORE MACT IN MACT1 AND IACT IN IACT1 FOR POSSIBLE LATER USE.
      if ( Iyrct(1)>=0 ) then
! HERE IYRCT(1) /= -1, AND WE CONSIDER THE PRESENT ENTRIES IN IYRCT
! ONE BY ONE.
         do j = 1 , npar1
            jj = Iyrct(j)
            if ( jj<=Mact1 ) then
! HERE ENTRY J OF IYRCT CORRESPONDS TO A FORMER ACTIVE CONSTRAINT AT
! SOME POINT abs(KK), WHERE THE SIGN OF KK WILL INDICATE WHETHER THE
! CONSTRAINT WAS +ACTIVE OR -ACTIVE.
               kk = Iact1(jj)
! WE NOW CHECK TO SEE IF THIS FORMER ACTIVE CONSTRAINT IS STILL
! ACTIVE WITH THE SAME SIGN.  IF SO, WE RESET IYRCT(J) TO THE PRESENT
! NUMBER OF THIS CONSTRAINT, AND IF NOT (WHICH WILL OCCUR IFF THE K
! LOOP BELOW IS COMPLETED), WE WILL NOT TRY TO DETERMINE A VERTEX, SO
! WE WILL SET IYRCT(1)=-1 AND LEAVE THE J LOOP.
               do k = 1 , mact
                  if ( kk==Iact(k) ) then
                     Iyrct(j) = k
                     goto 350
                  end if
               end do
               Iyrct(1) = -1
               goto 500
            else
! HERE ENTRY J OF IYRCT CORRESPONDS TO A CONSTRAINT BEYOND THE ACTIVE
! POINT CONSTRAINTS, AND WE ADJUST IYRCT(J) BY THE DIFFERENCE OF THE
! PRESENT AND FORMER NUMBER OF ACTIVE CONSTRAINTS.
               Iyrct(j) = Iyrct(j) + mact - Mact1
            end if
 350     end do
! WE HAVE NOW FILLED IN IYRCT(1),...,IYRCT(NPARM+1) WITH DISTINCT
! POSITIVE INTEGERS BETWEEN 1 AND M, AND WE FILL IN THE REST OF IYRCT
! SO THAT IYRCT WILL CONTAIN A PERMUTATION OF 1,...,M.  TO BE CONSISTENT
! WITH SLNPRO WE PUT IYRCT(NPARM+2),...,IYRCT(M) IN DECREASING ORDER.
         l = npar1
         do i = 1 , m
            ii = m - i + 1
! SKIP II IF IT IS ALREADY IN IYRCT.
            do j = 1 , npar1
               if ( ii==Iyrct(j) ) goto 400
            end do
            l = l + 1
            Iyrct(l) = ii
 400     end do
      end if
!
! SAVE MACT IN MACT1 AND IACT IN IACT1 AND RETURN.
 500  Mact1 = mact
      do j = 1 , mact
         Iact1(j) = Iact(j)
      end do
    end subroutine setu1
!********************************************************************************

!********************************************************************************
!>
!  THIS SUBROUTINE SOLVES THE LINEAR PROGRAMMING PROBLEM
!    MAXIMIZE Z = -V(M+1,1)*X(1)-...-V(M+1,N)*X(N)
!    WHERE X(1),...,X(N) ARE FREE VARIABLES, SUBJECT TO
!    V(I,1)*X(1)+...+V(I,N)*X(N) <= V(I,N+1), FOR I=1,..,M,
!    WHERE M >= N.
!  (INFORMATION CONCERNING TOLERANCES AND BASIC VARIABLES
!  IS ALSO TRANSMITTED USING M, N, AND IYRCT.)
!
! GIVEN INTEGERS M AND N (WITH M >= N) AND A MATRIX V,
! THIS SUBROUTINE SOLVES THE LINEAR PROGRAMMING PROBLEM
!    MAXIMIZE Z=-V(M+1,1)X(1)-...-V(M+1,N)X(N)+V(M+1,N+1)
! SUBJECT TO THE CONSTRAINTS
!    V(I,1)X(1)+...+V(I,N)X(N) <= V(I,N+1), I=1,...,M
! USING ESSENTIALLY THE METHOD IN THE BOOK BY AVDEYEVA AND
! ZUKHOVITSKIY.  Y(I)=-V(I,1)X(1)-...-V(I,N)X(N)+V(I,N+1),
! I=1,...,M ARE SLACK VARIABLES.  THE METHOD HAS 4 PHASES.
!
! FIRST, XS ARE EXCHANGED FOR YS TO GET A PROBLEM
! INVOLVING ONLY NONNEGATIVE VARIABLES, THE YS BEING
! SELECTED IN THE ORDER DETERMINED BY IYRCT AND A PIVOTING
! STRATEGY.  AT THE BEGINNING OF THIS ROUTINE WE MUST HAVE
! IYRCT(1) NONPOSITIVE, OR IYRCT MUST CONTAIN SOME
! PERMUTATION OF THE INTEGERS 1,...,M (SEE BELOW).
! SECOND, THE SLACK CONSTANTS OF THE DUAL PROBLEM ARE MADE
! (SIGNIFICANTLY) NONNEGATIVE.
! THIRD, THE COST COEFFICIENTS OF THE DUAL PROBLEM ARE MADE
! (SIGNIFICANTLY) NONNEGATIVE.
! FINALLY, THE SOLUTION VECTOR IS COMPUTED.
!
! THE VARIABLE INDIC WILL BE GIVEN VALUE
! 0, IF A SOLUTION WAS FOUND NORMALLY
! 1, IF THERE WAS TROUBLE IN PHASE 1
! 2, IF THERE WAS TROUBLE IN PHASE 2 (EITHER ROUND OFF
!   ERROR, OR THE ORIGINAL PROBLEM WAS INCONSISTENT OR
!   UNBOUNDED)
! 3, IF THERE WAS TROUBLE IN PHASE 3 (EITHER ROUND OFF
!   ERROR, OR THE ORIGINAL CONSTRAINTS WERE INCONSISTENT)
! 4, IF LIMJOR MODIFIED JORDAN ELIMINATIONS WERE USED IN
!   PHASES 2 AND 3
! -1, IF A SOLUTION WAS FOUND BUT IN ORDER TO OVERCOME
!   TROUBLE IN PHASE 2 OR 3 IT WAS NECESSARY TO TEMPORARILY
!   RELAX THE RESTRICTION ON DIVISION BY NUMBERS WITH SMALL
!   ABSOLUTE VALUE.  THUS THERE IS AN INCREASED CHANCE OF
!   SERIOUS ROUNDOFF ERROR IN THE RESULTS.
! -2, IF A SOLUTION WAS FOUND NORMALLY, EXCEPT THAT
!   THE PARAMETERS REA AND REA1 WERE INCREASED BY A SIGNAL
!   FROM THE CALLING PROGRAM (NAMELY, M=-M).  THE INCREASED
!   TOLERANCES MAY HAVE ALLOWED MORE ERROR THAN USUAL.
! -3, IF IN ORDER TO COMPLETE PHASE 1 IT WAS NECESSARY TO
!   TEMPORARILY RELAX THE RESTRICTION ON DIVISION BY NUMBERS
!   WITH SMALL ABSOLUTE VALUE.  THUS THERE IS AN INCREASED
!   CHANCE OF SERIOUS ROUNDOFF ERROR IN THE RESULTS.
! -4, IF A SOLUTION WAS FOUND NORMALLY, EXCEPT THAT REA AND REA1
!   WERE DECREASED BY A SIGNAL FROM THE CALLING PROGRAM (NAMELY
!   N=-N) IN ORDER TO TRY FOR MORE ACCURACY.  THIS INCREASES THE
!   CHANCES OF SERIOUS ROUNDOFF ERROR IN THE RESULTS.
! NOTE THAT INDIC=-3 WILL OVERWRITE (AND THUS CONCEAL) INDIC=-2
!   OR INDIC=-4, AND INDIC=-1 WILL OVERWRITE INDIC=-2, -3, OR -4
!
!###REFERENCE
!  * AVDEYEVA, L. I. AND ZUKHOVITSKIY, S. I.,
!    LINEAR AND CONVEX PROGRAMMING,
!    SAUNDERS, PHILADELPHIA, 1966.

    subroutine slnpro(v,m,n,Iyrct,y,Ixrct,Iycct,Nparm,Numgr,x,Indic)

      implicit none

      real(wp) :: absv , amax , amin , ampr2 , amprv , bmpr2 , bmprv , &
                 dist , dist1 , rea , rea1 , rea2 ,   &
                 rea3 , reakp , rowq , rtcol
      real(wp) :: temp , v , x , y
      integer :: i , i1 , i10 , i2 , i20 , iback , id , ifail , ii ,     &
                 inamp , Indic , indst , irlax , irow , itemp , Ixrct ,  &
                 ixrj , Iycct , iycj , Iyrct
      integer :: iyri , iytmp , j , jj , k , keep , keep1 , kkk , kpmp1 ,&
                 kpmp2 , ktjor , l , limjor , ll , lrknt , m , mp1 ,     &
                 mxrkn , n , np1
      integer :: Nparm , Numgr

      dimension v(Numgr+2*Nparm+1,Nparm+2) , Iyrct(Numgr+2*Nparm) ,    &
                x(Nparm+1) , y(Numgr+2*Nparm) , Ixrct(Numgr+2*Nparm) , &
                Iycct(Nparm+1)

! SET MACHINE DEPENDENT PARAMETERS FOR SUBROUTINE SLNPRO.
! SET SPCMN=B**(-ITT), WHERE B IS THE BASE AND ITT IS THE NUMBER
! OF BASE B DIGITS IN THE MANTISSA.  SPCMN IS THE MINIMUM
! RELATIVE SPACING ABS((X1-X2)/X2) BETWEEN TWO SUCCESSIVE
! FLOATING POINT NUMBERS, SO IT IS THE SPACING BETWEEN TWO
! SUCCESSIVE FLOATING POINT NUMBERS IN THE CLOSED INTERVAL
! (0.1,1.0).  WE ALSO HAVE SPCMN=10.0**(-ITT*(LOG10)(B))=
! 10.0**(-TNMAN), WHERE TNMAN IS THE BASE 10 EQUIVALENT OF
! THE LENGTH OF THE MANTISSA.

! SET REA (ROUND OFF ERROR ADJUSTMENT) =
! MAX(10.0**(-8),100.0*SPCMN).  THUS REA=10.0**(-EXREA),
! WHERE EXREA=MIN(8,TNMAN-2).
! DIVISION BY NUMBERS <= REA IN ABSOLUTE VALUE WILL NOT BE
! PERMITTED.
      rea = ten*ten*spcmn
      if ( rea<ten**(-8) ) rea = ten**(-8)
! SET REA1=10.0*SPCMN.  THUS REA1=10.0**(-(TNMAN-1)).
! NUMBERS IN ROW M+1 OR COLUMN N+1 WHICH ARE <= REA1 IN
! ABSOLUTE VALUE WILL BE TREATED AS ZEROES.  SLNPRO ASSUMES
! THAT 0.0 < REA1 <= REA.
      rea1 = ten*spcmn
! END OF INITIAL SETTING OF MACHINE DEPENDENT PARAMETERS FOR
! SLNPRO.  THESE PARAMETERS MAY BE ADJUSTED BY A COMMAND FROM
! THE CALLING PROGRAM.
!
      Indic = 0
      limjor = 300
! M BEING NEGATIVE IS A SIGNAL TO INCREASE REA AND REA1,
! THUS TREATING MORE NUMBERS WITH SMALL ABSOLUTE VALUES AS
! ZEROES.  THIS MAY GIVE THIS ROUTINE A BETTER CHANCE TO
! SUCCEED, BUT MAY ALSO CAUSE LARGER ERRORS.
      if ( m<0 ) then
! RESET M.
         m = -m
         rea = sqrt(rea)
         rea1 = sqrt(rea1)
         Indic = -2
      end if
! N BEING NEGATIVE IS A SIGNAL TO DECREASE REA AND REA1 TO TRY
! FOR MORE ACCURACY.  AMONG OTHER THINGS, THIS MAKES IT MORE
! LIKELY THAT THE PREVIOUS VERTEX WILL BE RETAINED IN PHASE 1
! BELOW, BUT IT ALSO COULD INCREASE ROUND OFF ERROR.
      if ( n<0 ) then
! RESET N.
         n = -n
         rea = rea1
         rea1 = rea1/(ten*ten)
         Indic = -4
      end if
! PRESERVE REA IN CASE IT MUST BE TEMPORARILY RELAXED.
! IRLAX=0 INDICATES REA IS NOT RELAXED AT THIS STAGE.
      reakp = rea
      irlax = 0
! IN COLUMN N+1, NUMBERS <= REA2 IN ABSOLUTE VALUE WILL BE
! TREATED AS ZEROES.
      rea2 = rea1
      np1 = n + 1
      mp1 = m + 1
      ktjor = 0
      iback = 0
! SET V(MP1,NP1)=0.0 SO THE DESCRIPTIONS IN AND FOLLOWING THE
! PROLOGUE WILL AGREE.
      v(mp1,np1) = zero
! THE ONLY REASON FOR THE FOLLOWING THREE STATEMENTS IS TO
! AVOID THE ERROR MESSAGE ON SOME MACHINES THAT THESE
! VARIABLES HAVE NOT BEEN ASSIGNED A VALUE.  THEY WILL BE
! REASSIGNED A VALUE BEFORE THE PROGRAM REACHES A SPOT WHERE
! THEY WILL ACTUALLY BE USED.
      dist = one
      amprv = one
      ampr2 = one
! SET IXRCT.  IXRCT(I)=0 MEANS SOME Y IS IN ROW I, WHILE
! IXRCT(I)=K/=0 MEANS X(K) IS IN ROW I.
      do i = 1 , m
         Ixrct(i) = 0
      end do
!
! EXCHANGE THE XS AT THE TOP OF THE TABLE FOR YS.
! IF IYRCT(1) IS NONPOSITIVE, WE SET IYRCT AND CHOOSE THE
! LARGEST POSSIBLE RESOLVENTS FOR THE EXCHANGES.  IF
! IYRCT(1) IS POSITIVE, IYRCT WILL HAVE BEEN PREVIOUSLY SET
! AND WE TRY TO EXCHANGE IN ROWS IYRCT(1),...,IYRCT(N),
! STILL EMPLOYING A PIVOTING STRATEGY, BUT IF WE CANNOT, WE
! EXCHANGE IN ROWS IYRCT(N+1),...,IYRCT(M).
      if ( Iyrct(1)<=0 ) then
         i10 = 1
         i20 = m
! IF WE HAVE NO INFORMATION FROM A PREVIOUS VERTEX, WE GIVE
! UP A LITTLE ACCURACY IN COLUMN N+1 TO HAVE A BETTER CHANCE
! OF SUCCESS.
         rea2 = rea
! THIS ROUTINE HAS A BACKTRACKING OPTION WHICH SOMETIMES
! INCREASES ACCURACY BUT SOMETIMES LEADS TO FAILURE DUE TO
! CYCLING.  IT IS SUGGESTED THAT THIS OPTION BE EMPLOYED IF
! INFORMATION ABOUT A STARTING VERTEX IS AVAILABLE, AND
! OTHERWISE BE DISABLED BY SETTING IBACK=1.
         iback = 1
         do i = 1 , m
            Iyrct(i) = i
         end do
      else
         i10 = 1
         i20 = n
      end if
      j = 0
! SET THE LOWER BOUND ON THE ABSOLUTE VALUE OF A RESOLVENT IN
! PHASE 1.  ALSO SET IFAIL=0 TO INDICATE THE RESOLVENT SEARCH
! IN THIS COLUMN HAS NOT FAILED.
      rea3 = rea
      ifail = 0
 100  j = j + 1
      if ( j>n ) then
!
! REARRANGE THE ROWS OF V SO THAT X(1),...,X(N) COME FIRST
! IN THAT ORDER.  REDEFINE IYRCT SO THAT AFTER THE
! REARRANGEMENT IS DONE, IYRCT(I)=K WILL MEAN Y(K) IS IN
! ROW I (FOR I GREATER THAN N).
         do i = 1 , m
            Iyrct(i) = i
         end do
         irow = 0
         goto 400
      end if
! SET I1, I2 ACCORDING TO THE STRATEGY WE ARE USING.
 200  i1 = i10
      i2 = i20
      amax = zero
! SEARCH FOR A RESOLVENT IN ROWS IYRCT(I1),...,IYRCT(I2).
 300  do i = i1 , i2
         iytmp = Iyrct(i)
         if ( Ixrct(iytmp)==0 ) then
            absv = abs(v(iytmp,j))
            if ( absv>amax ) then
               iyri = iytmp
               amax = absv
            end if
         end if
      end do
! CHECK TO SEE IF THE PROSPECTIVE RESOLVENT IS LARGE ENOUGH
! IN ABSOLUTE VALUE.
      if ( amax>rea3 ) then
! EXCHANGE X(J) FOR Y(IYRI).
         call sjelim(mp1,1,np1,iyri,j,Nparm,Numgr,v)
         Ixrct(iyri) = j
         Iycct(j) = iyri
! IYCCT(J)=IYRI MEANS Y(IYRI) IS IN COLUMN J.
! RESET REA3 AND IFAIL SINCE WE SUCCESSFULLY FOUND A RESOLVENT IN
! THIS COLUMN, AND THE RESOLVENT SEARCH IN THE NEXT COLUMN HAS
! NOT FAILED.
         rea3 = rea
         ifail = 0
         goto 100
! WE HAVE NOT FOUND A SUITABLE RESOLVENT IN ROWS IYRCT(I1),
! ...IYRCT(I2).  IF I2 < M WE SEARCH THE REST OF COLUMN J.
      else if ( i2<m ) then
         i1 = i2 + 1
         i2 = m
         goto 300
! HERE WE FAILED TO FIND A RESOLVENT IN COLUMN J WITH ABSOLUTE
! VALUE > REA3.  IF IFAIL=0 WE SET INDIC=-3 AND TRY AGAIN
! WITH REA3 REDUCED.  IF THIS HAS ALREADY BEEN TRIED WE SET
! INDIC=1 AND RETURN.
      else if ( ifail<=0 ) then
         ifail = 1
         Indic = -3
         rea3 = rea1
         goto 200
      else
!
         Indic = 1
         return
      end if
 400  irow = irow + 1
      if ( irow<=m ) then
         if ( Ixrct(irow)==0 ) goto 400
         if ( Ixrct(irow)==irow ) goto 400
      else
! NOW IXRCT IS NO LONGER NEEDED, SO STORE THE PRESENT IYCCT
! IN IT.
         do i = 1 , n
            Ixrct(i) = Iycct(i)
         end do
! END OF PHASE 1.
!
! THE FIRST N ROWS OF V GIVE THE XS IN TERMS OF CERTAIN
! YS.  THESE ROWS WILL NOT BE CHANGED BY LATER OPERATIONS.
!
! WE NOW ATTACK THE MAXIMIZATION PROBLEM BY LOOKING AT THE
! DUAL PROBLEM OF MINIMIZING A FORM GIVEN BY THE
! COEFFICIENTS IN V(N+1,N+1) THROUGH V(M,N+1) WITH SLACK
! TERMS IN THE BOTTOM ROW OF V.
! SEARCH ROW M+1 FOR A SIGNIFICANTLY NEGATIVE ELEMENT.  IF
! THERE ARE NONE, PROCEED TO THE ACTUAL MINIMIZATION
! PROBLEM.  STICK WITH COLUMN JJ UNTIL V(M+1,JJ) >= -REA1.
         jj = 0
         goto 600
      end if
! NOW X(L) IS IN ROW IROW, BUT WE WANT IT IN ROW L.
 500  l = Ixrct(irow)
      ll = Ixrct(l)
      if ( ll/=0 ) then
! X(L) IS IN ROW IROW, WHILE X(LL) IS IN ROW L.
         Ixrct(irow) = ll
         Ixrct(l) = l
      else
! X(L) IS IN ROW IROW, WHILE Y(IYRCT(L)) IS IN ROW L.
         Ixrct(irow) = 0
         Iyrct(irow) = Iyrct(l)
         Ixrct(l) = l
      end if
! NOW EXCHANGE THE CONTENTS OF ROWS IROW AND L.
      do j = 1 , np1
         temp = v(irow,j)
         v(irow,j) = v(l,j)
         v(l,j) = temp
      end do
      if ( Ixrct(irow)==0 ) goto 400
      if ( Ixrct(irow)==irow ) goto 400
      goto 500
 600  jj = jj + 1
      if ( jj>n ) then
!
! IN THE UNLIKELY EVENT THAT SOME V(M+1,J) IS STILL VERY
! SIGNIFICANTLY NEGATIVE WE BACKTRACK TO COLUMN J.  THIS
! COULD NOT HAPPEN IF THERE WERE NO ROUNDOFF ERROR AND WE
! COULD ALLOW DIVISION BY NUMBERS WITH VERY SMALL ABSOLUTE
! VALUE.  OMIT BACKTRACKING IF IBACK=1.
         if ( iback<=0 ) then
            do j = 1 , n
               if ( v(mp1,j)+rea<=0 ) then
                  jj = j
                  goto 700
               end if
            end do
         end if
         goto 900
      else if ( v(mp1,jj)+rea1>=0 ) then
         goto 600
      end if
!
! WE HAVE V(M+1,JJ) SIGNIFICANTLY NEGATIVE.  SEARCH COLUMN
! JJ FOR A POSITIVE ELEMENT, TREATING A VERY SMALL V(I,J)
! AS A ZERO.  IF THERE ARE NO POSITIVE ELEMENTS THE DUAL
! CONSTRAINTS WERE INCONSISTENT, SO THE ORIGINAL PROBLEM WAS
! INCONSISTENT OR UNBOUNDED.
 700  i = n
      inamp = 0
 800  i = i + 1
      if ( i<=m ) then
         if ( v(i,jj)<=rea ) goto 800
!
! NOW V(I,JJ) > REA.  WE SEARCH ROW I FOR INDICES K SUCH
! THAT V(M+1,K) >= 0.0.OR.K < JJ, AND V(I,K) < -REA, AND
! FIND THE MAXIMUM RATIO (I.E. THE RATIO WITH SMALLEST
! ABSOLUTE VALUE, IF V(M+1,K) >= 0.0) V(M+1,K)/V(I,K).  IF
! THERE IS NO SUCH K WE LOOK AT POSITIVE V(I,K) BELOW.
         indst = 0
         do j = 1 , n
            if ( v(mp1,j)<0 ) then
               if ( j>=jj ) goto 850
            end if
            if ( v(i,j)+rea<0 ) then
               dist1 = v(mp1,j)/v(i,j)
               if ( indst>0 ) then
                  if ( dist1<=dist ) goto 850
               end if
               dist = dist1
               indst = 1
               k = j
            end if
 850     end do
         if ( indst<=0 ) then
!
! IF THERE WAS NO INDEX K SUCH THAT V(M+1,K) >= 0.0.OR.K <
! JJ, AND V(I,K) < -REA, WE LOOK FOR THE SMALLEST (I.E.
! LARGEST IN ABSOLUTE VALUE) RATIO V(M+1,K)/V(I,K) FOR
! V(I,K) > REA AND V(M+1,K) < 0.0, AND PERFORM AN
! ELIMINATION WITH RESOLVENT V(I,K).  THERE IS AT LEAST ONE
! SUCH K, NAMELY JJ.
! THIS WILL FINISH PHASE 2 UNLESS BACKTRACKING IS NECESSARY.
            dist = one
            do j = 1 , n
               if ( v(mp1,j)<0 ) then
                  if ( v(i,j)>rea ) then
                     dist1 = v(mp1,j)/v(i,j)
                     if ( dist1<dist ) then
                        dist = dist1
                        k = j
                     end if
                  end if
               end if
            end do
         else
!
! WE NOW COMPUTE V(I,JJ)*DIST AND GO ON TO LOOK AT OTHER
! ROWS TO MINIMIZE THIS QUANTITY (I.E. TO MAXIMIZE ITS
! ABSOLUTE VALUE, IF V(M+1,K) >= 0.0).  THIS IS THE NEGATIVE
! OF THE CHANGE IN V(M+1,JJ).
            bmprv = v(i,jj)*dist
            if ( inamp>0 ) then
               if ( bmprv>=amprv ) goto 800
            end if
            amprv = bmprv
            inamp = 1
            kpmp1 = i
            kpmp2 = k
! (KPMP1,KPMP2) GIVES THE POSITION OF THE BEST PROSPECTIVE
! RESOLVENT FOUND SO FAR.
            goto 800
         end if
!
      else if ( inamp<=0 ) then
! AT THIS POINT INAMP IS POSITIVE IFF THERE WAS AT LEAST ONE
! ELEMENT > REA IN COLUMN JJ.  IF THERE WERE NONE, WE
! TEMPORARILY RELAX REA AND TRY AGAIN.
         if ( irlax<=0 ) then
            irlax = 1
            Indic = -1
            rea = rea1
            goto 700
         else
!
            Indic = 2
            return
         end if
      else
!
! CHECK TO SEE IF V(MP1,KPMP2) IS VERY SMALL IN ABSOLUTE
! VALUE OR NEGATIVE.  THIS INDICATES DEGENERACY.
         if ( v(mp1,kpmp2)<=rea ) then
!
! WE ARE NOW STUCK IN DEGENERATE COLUMN KPMP2.  WE SEARCH
! EACH DEGENERATE COLUMN IN WHICH WE ARE STUCK FOR A
! RESOLVENT WHICH WILL KEEP US FROM GETTING STUCK IN THIS
! COLUMN NEXT TIME, AND TO REDUCE THE ROUND-OFF ERROR WE
! TAKE THE SMALLEST OF THESE (I.E. LARGEST IN ABSOLUTE
! VALUE) AS OUR ACTUAL RESOLVENT.
            amin = one
            do j = 1 , n
! COLUMN J MAY BE DEGENERATE IF 0.0 <= V(M+1,J) <= REA,
! OR V(M+1,J) < 0.0.AND.J < JJ.
               if ( v(mp1,j)<0 ) then
                  if ( j>=jj ) goto 860
               else if ( v(mp1,j)>rea ) then
                  goto 860
               end if
! WE WILL BE STUCK IN COLUMN J IFF THERE IS AN INDEX ID FOR
! WHICH V(ID,JJ) > REA AND V(ID,J) < -REA.  IF THIS IS THE
! CASE, CHOOSING SUCH AN ID SO THAT V(ID,JJ)/V(ID,J) IS
! MINIMIZED (I.E. MAXIMIZED IN ABSOLUTE VALUE) AND TAKING
! V(ID,J) AS THE RESOLVENT WILL INSURE THAT WE DONT GET
! STUCK IN COLUMN J NEXT TIME.
               dist = one
               do i = np1 , m
                  if ( v(i,jj)>rea ) then
                     if ( v(i,j)+rea<0 ) then
                        dist1 = v(i,jj)/v(i,j)
                        if ( dist1<dist ) then
                           dist = dist1
                           id = i
                        end if
                     end if
                  end if
               end do
               if ( dist<one/two ) then
! WE HAVE NOW DETERMINED THAT WE ARE STUCK IN COLUMN J.
! IF V(ID,J) < AMIN THEN V(ID,J) IS THE BEST RESOLVENT
! FOUND SO FAR.
                  if ( v(id,j)<amin ) then
                     amin = v(id,j)
                     kpmp1 = id
                     kpmp2 = j
                  end if
               end if
! THE BEST RESOLVENT IS V(KPMP1,KPMP2), SO WE DO AN
! ELIMINATION.
 860        end do
         end if
! DO AN ELIMINATION WITH RESOLVENT V(KPMP1,KPMP2).
         i = kpmp1
         k = kpmp2
      end if
!
      ktjor = ktjor + 1
      if ( ktjor>limjor ) goto 1200
      call sjelim(mp1,np1,np1,i,k,Nparm,Numgr,v)
      itemp = Iyrct(i)
      Iyrct(i) = Iycct(k)
      Iycct(k) = itemp
! RESET REA AND IRLAX.
      rea = reakp
      irlax = 0
! IF NOW V(M+1,JJ) HAS BEEN MADE NOT SIGNIFICANTLY NEGATIVE,
! WE GO TO THE NEXT COLUMN.  OTHERWISE WE TRY AGAIN IN
! COLUMN JJ.
      if ( v(mp1,jj)+rea1>=0 ) goto 600
      goto 700
! END OF PHASE 2.
!
 900  i = n
      kkk = 0
!
! SEARCH FOR A SIGNIFICANTLY NEGATIVE ELEMENT BETWEEN
! V(N+1,N+1) AND V(N+1,M).  IF THERE ARE NONE WE HAVE THE
! MINIMAL POINT OF THE DUAL PROBLEM (AND THUS THE MAXIMAL
! POINT OF THE DIRECT PROBLEM) ALREADY.
 1000 i = i + 1
      if ( i<=m ) then
         if ( v(i,np1)+rea2>=0 ) goto 1000
!
! SEARCH FOR A NEGATIVE ELEMENT IN ROW I, TREATING A NUMBER
! WHICH IS VERY SMALL IN ABSOLUTE VALUE AS A ZERO.  IF THERE
! ARE NO NEGATIVE ELEMENTS THE DUAL PROBLEM WAS UNBOUNDED
! BELOW, SO THE ORIGINAL CONSTRAINTS WERE INCONSISTENT.
! FIND THE INDEX K OF THE LARGEST (I.E. SMALLEST IN ABSOLUTE
! VALUE, IF V(M+1,K) >= 0.0) RATIO V(M+1,K)/V(I,K) WITH
! V(I,K) < -REA.
 1050    indst = 0
         do j = 1 , n
            if ( v(i,j)+rea<0 ) then
               dist1 = v(mp1,j)/v(i,j)
               if ( indst>0 ) then
                  if ( dist1<=dist ) goto 1100
               end if
               k = j
               indst = 1
               dist = dist1
            end if
 1100    end do
         if ( indst>0 ) then
!
! COMPUTE THE IMPROVEMENT DIST*V(I,N+1) IN THE VALUE OF THE
! FORM USING V(I,K) AS THE RESOLVENT.  SET KKK=1 TO INDICATE
! A SIGNIFICANTLY NEGATIVE V(I,N+1) WAS FOUND, AND LOOK AT
! THE OTHER ROWS TO FIND THE RESOLVENT GIVING THE LARGEST
! IMPROVEMENT.
            bmpr2 = dist*v(i,np1)
! RESET IRLAX SO THAT THE NEXT ROW WHICH NEEDS RELAXING DOES
! NOT TERMINATE THE ROUTINE.  REA WILL REMAIN RELAXED UNTIL
! AFTER THE NEXT ELIMINATION.
            irlax = 0
            if ( kkk>0 ) then
               if ( bmpr2<=ampr2 ) goto 1000
            end if
            kkk = 1
            keep = i
            keep1 = k
            ampr2 = bmpr2
            goto 1000
! RELAX REA AND LOOK FOR NEGATIVE ELEMENTS WITH SMALLER
! ABSOLUTE VALUE.
         else if ( irlax<=0 ) then
            irlax = 1
            Indic = -1
            rea = rea1
            goto 1050
         else
!
            Indic = 3
            return
         end if
! KKK=0 HERE IFF NONE OF THE COST COEFFICIENTS ARE
! SIGNIFICANTLY NEGATIVE.
      else if ( kkk/=0 ) then
! CHECK TO SEE IF V(MP1,KEEP1) IS VERY SMALL IN ABSOLUTE
! VALUE OR NEGATIVE.  THIS INDICATES DEGENERACY.
         if ( v(mp1,keep1)<=rea ) then
!
! WE ARE NOW STUCK IN DEGENERATE COLUMN KEEP1.  WE SEARCH
! EACH DEGENERATE COLUMN IN WHICH WE ARE STUCK FOR A
! RESOLVENT WHICH WILL KEEP US FROM GETTING STUCK IN THIS
! COLUMN NEXT TIME.  IF WE ARE NOT USING THE OPTION
! DESCRIBED IN THE COMMENTS PRECEDING STATEMENT 1055, WE
! TAKE THE SMALLEST OF THESE (I.E. THE LARGEST IN ABSOLUTE
! VALUE) AS OUR ACTUAL RESOLVENT IN ORDER TO REDUCE THE
! GROWTH OF ROUND-OFF ERROR.
            amin = one
            mxrkn = np1
            do j = 1 , n
! COLUMN J MAY BE DEGENERATE IF V(M+1,J) <= REA.
               if ( v(mp1,j)<=rea ) then
! WE WILL BE STUCK IN COLUMN J IFF THERE IS AN INDEX ID FOR
! WHICH V(ID,N+1) < -REA2 AND V(ID,J) < -REA.  IF THIS
! IS THE CASE, CHOOSING SUCH AN ID SO THAT V(ID,N+1)/V(ID,J)
! IS MAXIMIZED AND TAKING V(ID,J) AS THE RESOLVENT WILL
! INSURE THAT WE DONT GET STUCK IN COLUMN J NEXT TIME.
                  dist = -one
                  do i = np1 , m
                     if ( v(i,np1)+rea2<0 ) then
                        if ( v(i,j)+rea<0 ) then
                           dist1 = v(i,np1)/v(i,j)
                           if ( dist1>dist ) then
                              dist = dist1
                              id = i
                           end if
                        end if
                     end if
                  end do
                  if ( dist+one/two<=0 ) goto 1120
!
! WE HAVE NOW DETERMINED THAT WE ARE STUCK IN COLUMN J.
! THE FOLLOWING STATEMENTS ATTEMPT TO BREAK DEGENERACY
! FASTER BY LOOKING ONE ITERATION INTO THE FUTURE, THAT IS,
! BY CHOOSING FROM THE PROSPECTIVE RESOLVENTS FOUND ABOVE
! THAT ONE WHICH MINIMIZES THE MINIMUM NUMBER OF STICKING
! PLACES IN ANY ROW AT THE NEXT STAGE.
! BECAUSE OF COMPUTER TIME AND THE POSSIBLE LOSS OF ACCURACY
! DUE TO LESSENED PIVOTING (EVEN THOUGH TIES ARE ALWAYS
! BROKEN IN FAVOR OF THE RESOLVENT WITH GREATEST ABSOLUTE
! VALUE), IT IS SUGGESTED THAT THIS OPTION BE OMITTED IF
! INFORMATION WAS AVAILABLE FROM A PREVIOUS VERTEX.  THIS
! WILL BE THE CASE IFF THE BACKTRACKING OPTION WAS USED,
! THAT IS, IFF IBACK=0.
                  if ( iback>0 ) then
! COMPUTE WHAT THE NEW BOTTOM ROW WOULD BE (EXCEPT FOR
! POSITION J) IF V(ID,J) WERE USED AS THE RESOLVENT, AND
! PUT THE RESULTS INTO Y.
                     rowq = v(mp1,j)/v(id,j)
                     do l = 1 , n
                        if ( l/=j ) y(l) = v(mp1,l) - v(id,l)*rowq
                     end do
                     lrknt = -1
! WE LOOK FOR A ROW WHICH WILL HAVE A SIGNIFICANTLY NEGATIVE
! LAST ELEMENT BUT A MINIMUM NUMBER OF PLACES WHERE WE WILL
! BE STUCK IN DEGENERATE COLUMNS.  LRKNT=-1 MEANS WE HAVE
! NOT YET FOUND A ROW WHICH WILL HAVE A SIGNIFICANTLY
! NEGATIVE LAST ELEMENT.
                     do ii = np1 , m
                        if ( ii/=id ) then
                           rowq = v(ii,j)/v(id,j)
                           rtcol = v(ii,np1) - v(id,np1)*rowq
                           if ( rtcol+rea2<0 ) then
! IF WE HAVE ALREADY LOCATED A RESOLVENT WHICH WILL FINISH
! THE ROUTINE, BUT THE PRESENT PROSPECTIVE RESOLVENT WOULD
! GIVE A ROW WITH A SIGNIFICANTLY NEGATIVE LAST ELEMENT, WE
! LOOK AT THE NEXT PROSPECTIVE RESOLVENT FOR PIVOTING
! PURPOSES.
                              if ( mxrkn+1==0 ) goto 1120
                              lrknt = 0
! NOW COUNT THE NUMBER (LRKNT) OF STICKING PLACES IN ROW II
! AT THE NEXT ITERATION.
                              do jj = 1 , n
                                 if ( jj/=j ) then
                                    if ( y(jj)<=rea ) then
                                       if ( v(ii,jj)-v(id,jj)*rowq+rea<0 ) then
                                         lrknt = lrknt + 1
                                         if ( lrknt>mxrkn ) goto 1102
                                       end if
                                    end if
                                 end if
                              end do
                              if ( lrknt<mxrkn ) then
                              else if ( lrknt==mxrkn ) then
                                 if ( v(id,j)>=amin ) goto 1102
                              else
                                 goto 1102
                              end if
                              mxrkn = lrknt
                              amin = v(id,j)
                              keep = id
                              keep1 = j
                           end if
                        end if
 1102                end do
! LRKNT=-1 HERE WOULD MEAN THIS RESOLVENT WOULD FINISH THE
! ROUTINE.  IF LRKNT >= 0 THEN MXRKN >= 0 ALSO, SO WE WILL
! NOT HAVE EARLIER FOUND A RESOLVENT WHICH WILL FINISH THE
! ROUTINE.
                     if ( lrknt+1/=0 ) goto 1120
                     if ( mxrkn+1/=0 ) goto 1105
                  end if
                  if ( v(id,j)>=amin ) goto 1120
 1105             mxrkn = -1
                  amin = v(id,j)
                  keep = id
                  keep1 = j
               end if
! THE BEST RESOLVENT IS V(KEEP,KEEP1), SO WE DO AN
! ELIMINATION.
 1120       end do
         end if
! DO AN ELIMINATION WITH RESOLVENT V(KEEP,KEEP1).
         i = keep
         k = keep1
!
         ktjor = ktjor + 1
         if ( ktjor<=limjor ) then
            call sjelim(mp1,np1,np1,i,k,Nparm,Numgr,v)
            itemp = Iyrct(i)
            Iyrct(i) = Iycct(k)
            Iycct(k) = itemp
! RESET REA AND IRLAX.
            rea = reakp
            irlax = 0
            goto 900
         end if
      else
! CHECK TO SEE IF ANY OF THE NUMBERS IN THE BOTTOM ROW HAVE
! BECOME VERY SIGNIFICANTLY NEGATIVE.  IF SO, WE MUST
! BACKTRACK TO PHASE 2 (SEE COMMENT ABOVE STATEMENT 1035).
! OMIT BACKTRACKING IF IBACK=1.
         if ( iback<=0 ) then
            do j = 1 , n
               if ( v(mp1,j)+rea<=0 ) then
                  jj = j
                  goto 700
               end if
            end do
         end if
! END OF PHASE 3.  WE NOW HAVE THE VERTEX WE ARE SEEKING.
!
! READ OFF THE Y VALUES FOR THIS VERTEX.
         do j = 1 , n
            iycj = Iycct(j)
            y(iycj) = zero
         end do
         do i = np1 , m
            iyri = Iyrct(i)
            y(iyri) = v(i,np1)
         end do
! COMPUTE THE XS FROM THE YS.  RECALL THAT IXRCT CONTAINS
! THE FORMER IYCCT.
         do i = 1 , n
            x(i) = v(i,np1)
            do j = 1 , n
               ixrj = Ixrct(j)
               x(i) = x(i) - v(i,j)*y(ixrj)
            end do
         end do
!
! NOW PUT THE VALUES IN IYCCT INTO THE FIRST N POSITIONS OF
! IYRCT IN DECREASING ORDER.
! TO ACCOMPLISH THIS, MAKE IXRCT(I)=-1 IF IYCCT(J)=I FOR
! SOME J, THEN SCAN IXRCT BACKWARDS.
         do j = 1 , n
            iycj = Iycct(j)
            Ixrct(iycj) = -1
         end do
         k = 1
         i = mp1
 1150    i = i - 1
         if ( i<=0 ) then
! NOW FILL IN THE REST OF IYRCT BY SCANNING IXRCT AGAIN.
            i = mp1
 1160       i = i - 1
            if ( i<=0 ) return
            if ( Ixrct(i)>=0 ) then
               Iyrct(k) = i
               k = k + 1
            end if
            goto 1160
         else
            if ( Ixrct(i)+1==0 ) then
               Iyrct(k) = i
               k = k + 1
            end if
            goto 1150
         end if
      end if
!
 1200 Indic = 4

    end subroutine slnpro
!********************************************************************************

!********************************************************************************
!>
!  This subroutine performs a modified jordan
!  elimination on the l-ll+1 by k matrix
!  consisting of rows ll through l of v and
!  columns 1 through k of v.  The resolvent
!  is v(ir,is).

    subroutine sjelim(l,Ll,k,Ir,Is,Nparm,Numgr,v)

    implicit none

    integer,intent(in) :: l
    integer,intent(in) :: Ll
    integer,intent(in) :: k
    integer,intent(in) :: Ir
    integer,intent(in) :: Is
    integer,intent(in) :: Nparm
    integer,intent(in) :: Numgr
    real(wp),intent(inout) :: v(Numgr+2*Nparm+1,Nparm+2)

    integer :: i , j
    real(wp) :: fact , resol

    ! DIVIDE THE ENTRIES IN THE RESOLVENT ROW (EXCEPT FOR THE
    ! RESOLVENT) BY THE RESOLVENT.
    resol = v(Ir,Is)
    do j = 1 , k
        if ( j/=Is ) v(Ir,j) = v(Ir,j)/resol
    end do
    ! SWEEP OUT IN ALL BUT ROW IR AND COLUMN IS.
    do i = Ll , l
        if ( i/=Ir ) then
            fact = -v(i,Is)
            do j = 1 , k
                if ( j/=Is ) v(i,j) = v(i,j) + v(Ir,j)*fact
            end do
        end if
    end do
    ! DIVIDE THE ENTRIES IN THE RESOLVENT COLUMN (EXCEPT FOR THE
    ! RESOLVENT) BY THE NEGATIVE OF THE RESOLVENT.
    do i = Ll , l
        if ( i/=Ir ) v(i,Is) = -v(i,Is)/resol
    end do
    ! REPLACE THE RESOLVENT BY ITS RECIPROCAL.
    v(Ir,Is) = one/resol

    end subroutine sjelim
!********************************************************************************

!********************************************************************************
!>
!  This subroutine uses a modified quadratic fitting process to
!  search for the minimum of a function f.  it requres an initial
!  guess in projct, a tolerance tol1 on the search interval length,
!  an upper bound prjlim on the minimizing point (which should be set
!  very large if no upper bound is desired), and a way to compute f(x)
!  for a given x.  the subroutine will print a warning and return if
!  it would need to compute f more than initlm times in the initialization
!  or more than nadd additional times in the main part of the program.
!  when the subroutine returns, it will have put the minimum location in
!  projct, the minimum f value in emin, the f value for the initial
!  projct in emin1, and the number of times it computed f in nsrch.

    subroutine searsl(me,Ioptn,Numgr,Nparm,Prjlim,Tol1,x,Fun,Ifun,Pttbl, &
                      Iptb,Indm,Param,Error,Rchdwn,Mact,Iact,Iphse,   &
                      Unit,Tolcon,Rchin,Itypm1,Itypm2,Iwork,Liwrk,    &
                      Work,Lwrk,Err1,Parprj,Projct,Emin,Emin1,Parser, &
                      Nsrch)

      implicit none

      class(conmax_solver),intent(inout) :: me
      real(wp) :: Emin , Emin1 , Err1 , &
                  Error , f1 , f2 , f3 , f4 , Fun , fval , fvlkp , &
                  p1 , p2 , p3
      real(wp) :: p4 , Param , Parprj , Parser , Prjlim , Projct , Pttbl ,   &
                  pval , Rchdwn , Rchin , rlf , rrt , s1 , s2 , tol4, &
                  Tol1 , Tolcon
      real(wp) :: Unit , Work , x
      integer :: Iact , icorct , Ifun , ilc08 , ilc10 , ilc17 , ilc21 , &
                 ilc27 , ilc29 , ilc48 , ilf , Indm , initlm , &
                 Ioptn , Iphse , ipmax , Iptb , irt , isave
      integer :: ismax , Itypm1 , Itypm2 , iupbar , Iwork , j , lims1 , &
                 Liwrk , lll , Lwrk , Mact , nadd , Nparm , Nsrch , Numgr

      dimension Fun(Ifun) , Pttbl(Iptb,Indm) , Param(Nparm) ,  &
                Err1(Numgr+3) , Parprj(Nparm) , x(Nparm+1) ,   &
                Error(Numgr+3) , Iact(Numgr) , Parser(Nparm) , &
                Iwork(Liwrk) , Work(Lwrk)

    real(wp),parameter :: tolden = ten*spcmn
    real(wp),parameter :: balfct = ten
    real(wp),parameter :: baladj = (ten-one)/ten

      tol4 = Tol1/four
      ilc08 = iloc(8,Nparm,Numgr)
      ilc10 = iloc(10,Nparm,Numgr)
      ilc17 = iloc(17,Nparm,Numgr)
      ilc21 = iloc(21,Nparm,Numgr)
      ilc27 = iloc(27,Nparm,Numgr)
      ilc29 = iloc(29,Nparm,Numgr)
      ilc48 = iloc(48,Nparm,Numgr)

! THE INITIAL PROJCT CAN BE INCREASED (OR DECREASED) BY A FACTOR OF
! 2.0**((INITLM-1)*INITLM-2)/2) (ASSUMING WE TAKE INITLM >= 3, AS
! WE SHOULD).  WE TAKE INITLM=6 SINCE A FACTOR OF 1024 SEEMS SUFFICIENT.
      initlm = 6
! NADD=4 SEEMS TO BE SUFFICIENT SINCE THIS NUMBER OF ITERATIONS PAST THE
! INITIALIZATION SEEMS TO ONLY RARELY BE EXCEEDED.
      nadd = 4
      Nsrch = 0
      ilf = 0
      irt = 0
      iupbar = 0
      isave = 0
! INITIALLY PUT PARAM IN PARSER SO THERE WILL BE SOMETHING THERE IF
! WE NEVER GET A CORRECTIBLE PARPRJ.
      do j = 1 , Nparm
         Parser(j) = Param(j)
      end do
! WE NOW TRY TO COMPUTE VALUES AT POINTS P2=PROJCT, P1=P2/2.0, AND
! P3=2.0*P2 (BUT P3 CANNOT EXCEED PRJLIM).
      p2 = Projct
! SET LLL=2 AS THE THREAD THROUGH THE MINOTAURS CAVERN AND JUMP
! DOWN TO PUT F(P2) IN F2.  WE WILL JUMP BACK AFTER ALL SUCH JUMPS
      lll = 2
      pval = p2
      goto 1200
! HERE SET LLL=3 AND PUT F(P3) IN F3.
 100  lll = 3
      pval = p3
      goto 1200
!
! WE NOW HAVE FOUND P1, P2, AND P3 WITH CORRESPONDING VALUES
! F1, F2, AND F3.  WE EXPAND THE INTERVAL IF NECESSARY TO TRY
! TO FIND NEW VALUES WITH F2 <= MIN(F1,F3).
 200  if ( f2<=f1 ) then
!
! HERE F2 <= F1.  IF F2 <= F3 AND WE HAVE NOT HAD ALL FAILURES OF
! THE F COMPUTATION, WE ARE DONE INITIALIZING.
         if ( f2>f3 ) goto 500
         goto 400
      else if ( f1>f3 ) then
         goto 500
      end if
!
! HERE WE WILL EXPAND THE INTERVAL TO THE LEFT, PROVIDING THAT
! NSRCH < INITLM AND P1-P1/2.0**(NSRCH-1) >= TOL4.
 300  if ( Nsrch<initlm ) then
         if ( p1-p1/two**(Nsrch-1)>=tol4 ) then
!
! EXPAND LEFT.
            p3 = p2
            f3 = f2
            p2 = p1
            f2 = f1
            p1 = p1/two**(Nsrch-1)
! SET LLL=5 AND PUT F(P1) IN F1.
            lll = 5
            pval = p1
            goto 1200
         end if
      end if
!
! HERE WE CANNOT EXPAND LEFT AND WE RETURN WITH THE BEST VALUES
! FOUND SO FAR.
      Projct = p1
      Emin = f1
      return
!
! HERE WE CHECK TO SEE IF THE F COMPUTATION HAS FAILED EVERY TIME
! (INDICATED BY F1=F2=F3=BIG), AND IF SO WE TRY TO EXPAND LEFT.
! IF NOT, WE ARE DONE WITH THE INITIALIZATION.
 400  if ( f1<big ) goto 700
      if ( f2<big ) goto 700
      if ( f3>=big ) goto 300
      goto 700
!
! HERE F3 < MIN(F1,F2) AND WE EXPAND THE INTERVAL TO THE RIGHT IF
! NSRCH < INITLM AND IUPBAR=0.
 500  if ( Nsrch<initlm ) then
         if ( iupbar<=0 ) then
!
! EXPAND RIGHT.
            p1 = p2
            f1 = f2
            p2 = p3
            f2 = f3
            p3 = two**(Nsrch-1)*p2
! IF P3 > PRJLIM, SET IUPBAR=1 AS AN INDICATOR WE CANNOT LATER
! EXPAND THE INTERVAL TO THE RIGHT.  THEN IF PRJLIM >= P2+TOL4
! REPLACE P3 BY PRJLIM, AND OTHERWISE RETURN WITH THE BEST VALUES
! FOUND SO FAR.
            if ( p3<=Prjlim ) goto 600
            iupbar = 1
            if ( Prjlim-p2<tol4 ) then
               Projct = p2
               Emin = f2
               return
            else
!
               p3 = Prjlim
               goto 600
            end if
         end if
      end if
!
! HERE WE CANNOT EXPAND RIGHT AND WE RETURN WITH THE BEST VALUES
! FOUND SO FAR.
      Projct = p3
      Emin = f3
      return
!
! SET LLL=6 AND PUT F(P3) IN F3.
 600  lll = 6
      pval = p3
      goto 1200
! END OF INITIALIZATION.
!
! ASSUMING THAT P3-P1 >= TOL1, WE NOW HAVE POINTS P1, P2, P3 WITH
! P1 <= P2-TOL4, P2 <= P3-TOL4, F1=F(P1) >= F2=F(P2), AND F3=F(P3)
! >= F2.  THESE CONDITIONS WILL BE MAINTAINED THROUGHOUT THE PROGRAM.
! SET LLL=7, WHERE IT WILL REMAIN FROM NOW ON.
 700  lll = 7
!
! RESET LIMS1 SO THAT AT MOST NADD MORE COMPUTATIONS OF F WILL BE DONE.
      lims1 = Nsrch + nadd
!
! IF WE HAVE COMPUTED F LIMS1 TIMES, WE PUT P2 IN PROJCT, PUT F2 IN
! EMIN, AND RETURN.
 800  if ( Nsrch<lims1 ) then
!
! IF THE SEARCH INTERVAL LENGTH IS LESS THAN TOL1 WE PUT P2 IN
! PROJCT, PUT F2 IN EMIN, AND RETURN.
         if ( p3-p1>=Tol1 ) then
!
! COMPUTE S1 = THE ABSOLUTE VALUE OF THE SLOPE OF THE LINE THROUGH
! (P1,F1) AND (P2,F2), AND S2 = THE (ABSOLUTE VALUE OF THE) SLOPE
! OF THE LINE THROUGH (P2,F2) AND (P3,F3).
!***MOD  CONSIDER INCREASING TOL1 TO 10**4*SPCMN
            s1 = (f1-f2)/(p2-p1)
            s2 = (f3-f2)/(p3-p2)
! IF S1+S2 IS VERY SMALL WE RETURN WITH THE BEST VALUES FOUND SO FAR.
            if ( s1+s2>=tolden ) then
!
               rlf = s2/(s1+s2)
               rrt = one - rlf
! THE MINIMUM OF THE QUADRATIC POLYNOMIAL PASSING THROUGH
! (P1,F1), (P2,F2), AND (P3,F3) WILL OCCUR AT (RLF*P1+
! RRT*P3+P2)/2.0.  NOTE THAT THE THREE POINTS CANNOT BE
! COLLNEAR, ELSE WE WOULD HAVE TERMINATED ABOVE.  SINCE THE
! MINIMUM OCCURS AT THE AVERAGE OF P2 AND A CONVEX COMBINATION
! OF P1 AND P3, IT WILL BE AT LEAST AS CLOSE TO P2 AS TO THE
! ENDPOINT ON THE SAME SIDE.
               if ( ilf>1 ) then
!
! HERE THE LEFT ENDPOINT WAS DROPPED AT THE LAST ILF > 1
! ITERATIONS, SO TO PREVENT A LONG STRING OF SUCH OCCURRENCES
! WITH LITTLE REDUCTION OF P3-P1 WE WILL SHIFT THE NEW POINT
! TO THE RIGHT BY DECREASING RLF RELATIVE TO RRT.
                  rlf = rlf/two**(ilf-1)
                  rrt = one - rlf
               else if ( irt>1 ) then
!
! HERE THE RIGHT ENDPOINT WAS DROPPED AT THE LAST IRT > 1
! ITERATIONS, AND WE WILL SHIFT THE NEW POINT TO THE LEFT.
                  rrt = rrt/two**(irt-1)
                  rlf = one - rrt
!
! HERE WE HAVE NOT JUST HAD A STRING OF TWO OR MORE MOVES IN
! THE SAME DIRECTION, BUT IF THE SUBINTERVALS ARE OUT OF
! BALANCE BY MORE THAN A FACTOR OF BALFCT, WE SHIFT THE NEW
! POINT SLIGHTLY IN THE DIRECTION OF THE LONGER INTERVAL.  THE
! IDEA HERE IS THAT THE TWO CLOSE POINTS ARE PROBABLY NEAR THE
! SOLUTION, AND IF WE CAN BRACKET THE SOLUTION WE MAY BE ABLE TO
! CUT OFF THE MAJOR PORTION OF THE LONGER SUBINTERVAL.
               else if ( p2-p1>balfct*(p3-p2) ) then
!
! HERE THE LEFT SUBINTERVAL IS MORE THAN BALFCT TIMES LONGER THAN
! THE RIGHT SUBINTERVAL, SO WE DECREASE RRT RRELATIVE TO RLF.
                  rrt = baladj*rrt
                  rlf = one - rrt
               else if ( p3-p2>balfct*(p2-p1) ) then
!
! HERE THE RIGHT SUBINTERVAL IS MORE THAN BALFCT TIMES LONGER
! THAN THE LEFT SUBINTERVAL, SO WE DECREASE RLF RELATIVE TO RRT.
                  rlf = baladj*rlf
                  rrt = one - rlf
               end if
!
! COMPUTE THE (POSSIBLY MODIFIED) MINIMUM OF THE QUADRATIC FIT.
               p4 = (rlf*p1+rrt*p3+p2)/two
!
! THE NEXT SECTION (FROM HERE TO STATEMENT 2800) MODIFIES P4 IF NECESSARY
! TO GET P1+TOL4 <= P2,P4 <= P3-TOL4 AND ABS(P4-P2) >= TOL4.  IN
! THE UNLIKELY EVENT THIS IS NOT POSSIBLE WE SET PROJCT=P2, EMIN=F2
! AND RETURN.
!
! IF ABS(P4-P2) < TOL4 WE REDEFINE P4 BY MOVING TOL4 FROM
! P2 INTO THE LONGER SUBINTERVAL.  NOTE THAT THE LENGTH OF THIS
! SUBINTERVAL MUST BE AT LEAST TOL1/2.0 = 2.0*TOL4, ELSE WE
! WOULD HAVE TERMINATED EARLIER.
               if ( abs(p4-p2)<tol4 ) then
                  if ( p3-p2>(p2-p1) ) goto 1000
                  goto 1100
! HERE WE HAD ABS(P4-P2) >= TOL4 AND WE MAKE SURE THAT P1+TOL4
! <= P4 <= P3-TOL4.
               else if ( p4<=(p3-tol4) ) then
                  if ( p4<(p1+tol4) ) then
! HERE P4 < P1+TOL4 AND WE SET P4=P1+TOL4 IF P2-P1 >= TOL1/2.0
! AND OTHERWISE WE SET P4=P2+TOL4.
                     if ( p2-p1<Tol1/two ) goto 1000
                     p4 = p1 + tol4
!
! NOW JUMP DOWN TO PUT F(P4) IN F4.
                     pval = p4
                     goto 1200
                  else
                     pval = p4
                     goto 1200
                  end if
               else
! HERE P4 > P3-TOL4 AND WE SET P4=P3-TOL4 IF P3-P2 >= TOL1/2.0,
! AND OTHERWISE WE SET P4=P2-TOL4.
                  if ( p3-p2<Tol1/two ) goto 1100
                  p4 = p3 - tol4
                  pval = p4
                  goto 1200
               end if
            end if
         end if
      end if
!
 900  Projct = p2
      Emin = f2
      return
 1000 p4 = p2 + tol4
! IF TOL4 WAS SMALL ENOUGH RELATIVE TO P2 THAT THE MACHINE THINKS P4
! STILL EQUALS P2, WHICH IS MORE LIKELY IF P2 IS LARGE, THIS COULD RESULT
! IN A DIVIDE FAULT LATER.  TO AVOID THIS, WE REDEFINE P4 AS THE AVERAGE
! OF P2 AND P3 IF NECESSARY.  IF WE STILL DONT HAVE P4 STRICTLY BETWEEN
! P2 AND P3, WE TERMINATE THE SEARCH.
      if ( p4<=p2 ) then
         p4 = (p2+p3)/two
         if ( p4<=p2 ) goto 900
      end if
      if ( p4>=p3 ) goto 900
      pval = p4
      goto 1200
 1100 p4 = p2 - tol4
! IF TOL4 WAS SMALL ENOUGH RELATIVE TO P2 THAT THE MACHINE THINKS P4
! STILL EQUALS P2, WHICH IS MORE LIKELY IF P2 IS LARGE, THIS COULD RESULT
! IN A DIVIDE FAULT LATER.  TO AVOID THIS, WE REDEFINE P4 AS THE AVERAGE
! OF P1 AND P2 IF NECESSARY.  IF WE STILL DONT HAVE P4 STRICTLY BETWEEN
! P1 AND P2, WE TERMINATE THE SEARCH.
      if ( p4>=p2 ) then
         p4 = (p1+p2)/two
         if ( p4>=p2 ) goto 900
      end if
      if ( p4<=p1 ) goto 900
      pval = p4
!
! INCREMENT NSRCH SINCE WE ARE ABOUT TO COMPUTE F.
 1200 Nsrch = Nsrch + 1
!
!
! THIS IS WHERE THE USER MUST SUPPLY A ROUTINE TO COMPUTE
! FVAL=F(PVAL).  IF THE PROCEDURE FAILS, SET FVAL=BIG.
!
      fval = big
! PROJECT X TO GET PARPRJ.
      do j = 1 , Nparm
         Parprj(j) = Param(j) + pval*x(j)
      end do
!
! CALL CORRCT TO RETURN PARPRJ TO FEASIBILITY IF NECESSARY IF ITYPM1
! OR ITYPM2 IS POSITIVE.
      if ( Itypm1+Itypm2>0 ) then
         call me%corrct(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,        &
                     Iwork(ilc17),Unit,Tolcon,Rchin,Error,Mact,Iact,    &
                     Projct,Iphse,Iwork,Liwrk,Work,Lwrk,Work(ilc27),    &
                     Err1,Work(ilc10),Work(ilc29),Work(ilc08),          &
                     Work(ilc48),Iwork(ilc21),Parprj,icorct)
         if ( icorct>0 ) goto 1300
      end if
      call me%ercmp1(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,Parprj,1,  &
                  Iphse,Iwork,Liwrk,Work(ilc08),Iwork(ilc17),ipmax,     &
                  ismax,Err1)
      fval = Err1(Numgr+1)
!
! IF NSRCH=1, INDICATING THAT WE ARE COMPUTING F WITH THE INITIAL PROJCT,
! CALL RCHMOD WITH IRCH=1 TO INCREASE RCHDWN FOR THE NEXT SETU1 OR
! RKSACT CALL IF NECESSARY.
      if ( Nsrch<=1 ) call rchmod(Numgr,Error,Err1,Iwork(ilc17),Mact,   &
                                  Iact,ipmax,ismax,Unit,1,Rchdwn,Rchin)
! WE WILL SAVE THE BEST PARPRJ FOUND IN THIS SEARSL CALL IN PARSER.
      if ( isave<=0 ) then
         isave = 1
      else if ( fval>=fvlkp ) then
         goto 1300
      end if
      do j = 1 , Nparm
         Parser(j) = Parprj(j)
      end do
      fvlkp = fval
! IF IPHSE < 0 AND FVAL <= TOLCON WE RETURN WITH THE BEST VALUES
! FOUND SO FAR.
      if ( Iphse<0 ) then
         if ( fval<=Tolcon ) then
            Projct = pval
            Emin = fval
            return
         end if
      end if
!
! CARRY THE COMPUTED F VALUE BACK TO THE APPROPRIATE PART OF THE PROGRAM.
 1300 select case (lll)
      case (1)
!
         f1 = fval
         p3 = two*p2
! IF P3 > PRJLIM, SET IUPBAR=1 AS AN INDICATOR WE CANNOT LATER
! EXPAND THE INTERVAL TO THE RIGHT.  THEN IF PRJLIM >= P2+TOL4
! REPLACE P3 BY PRJLIM, AND OTHERWISE EXPAND THE INTERVAL TO THE
! LEFT TO GET THE DESIRED THIRD POINT.
         if ( p3<=Prjlim ) goto 100
         iupbar = 1
         if ( Prjlim-p2<tol4 ) then
!
! EXPAND LEFT TO GET THE INITIAL THIRD POINT SINCE THERE IS NO ROOM
! TO EXPAND RIGHT.
            p3 = p2
            f3 = f2
            p2 = p1
            f2 = f1
            p1 = p1/two
! SET LLL=4 AND PUT F(P1) IN F1.
            lll = 4
            pval = p1
            goto 1200
         else
            p3 = Prjlim
            goto 100
         end if
      case (2)
!
         f2 = fval
! SET EMIN1 = THE VALUE OF F USING THE GIVEN PROJECTION FACTOR PROJCT.
         Emin1 = fval
         p1 = p2/two
! SET LLL=1 AND PUT F(P1) IN F1.
         lll = 1
         pval = p1
         goto 1200
      case (3)
!
         f3 = fval
         goto 200
      case (4)
!
         f1 = fval
         goto 200
      case (5)
         f1 = fval
!
! HERE F2 <= F3 AND WE HAVE JUST EXPANDED LEFT.  IF F2 > F1 WE
! TRY TO EXPAND LEFT AGAIN, OTHERWISE WE CHECK TO SEE IF WE ARE DONE
! INITIALIZING.
         if ( f2>f1 ) goto 300
         goto 400
      case (6)
         f3 = fval
!
! HERE F2 < F1 AND WE HAVE JUST EXPANDED RIGHT.  IF F2 <= F3
! WE ARE DONE INITIALIZING, OTHERWISE WE TRY TO EXPAND RIGHT AGAIN.
         if ( f2>f3 ) goto 500
         goto 700
      case (7)
!
         f4 = fval
!
! NOW WE DROP EITHER P1 OR P3 AND RELABEL THE REMAINING POINTS SO
! THAT F(P2) <= F(P1) AND F(P2) <= F(P3).
!
! IF NOW THE LEFTMOST OF THE TWO MIDDLE POINTS IS LOWER THAN THE
! RIGHTMOST OF THE TWO MIDDLE POINTS WE DROP P3, AND SET ILF=0
! AND INCREMENT IRT TO INDICATE THE RIGHT END POINT HAS BEEN DROPPED.
! OTHERWISE WE DROP P1, SET IRT=0 AND INCREMENT ILF.  IN ALL CASES
! WE THEN RESHUFFLE THE VALUES INTO P1, P2, P3, F1, F2, F3 AND TRY
! TO DO ANOTHER ITERATION.
         if ( p4<p2 ) then
!
! HERE P4 < P2.
            if ( f4<f2 ) then
               p3 = p2
               f3 = f2
               p2 = p4
               f2 = f4
               ilf = 0
               irt = irt + 1
            else
               p1 = p4
               f1 = f4
               ilf = ilf + 1
               irt = 0
            end if
!
! HERE P4 > P2.
         else if ( f2<f4 ) then
            p3 = p4
            f3 = f4
            ilf = 0
            irt = irt + 1
         else
            p1 = p2
            f1 = f2
            p2 = p4
            f2 = f4
            ilf = ilf + 1
            irt = 0
         end if
         goto 800
      case default
      end select
    end subroutine searsl
!********************************************************************************

!********************************************************************************
!>
!
    subroutine ercmp1(me,Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm, &
                      Param,Icnuse,Iphse,Iwork,Liwrk,Confun,Icntyp, &
                      Ipmax,Ismax,Error)

    implicit none

    class(conmax_solver),intent(inout) :: me
    integer,intent(in) :: Ioptn
    integer,intent(in) :: Nparm
    integer,intent(in) :: Numgr
    integer,intent(in) :: Ifun
    integer,intent(in) :: Iptb
    integer,intent(in) :: Indm
    integer,intent(in) :: Icnuse
    integer,intent(in) :: Iphse
    integer,intent(in) :: Liwrk
    integer,intent(out) :: Ipmax
    integer,intent(out) :: Ismax
    real(wp),intent(in) :: Fun(Ifun)
    real(wp),intent(in) :: Pttbl(Iptb,Indm)
    real(wp),intent(in) :: Param(Nparm)
    integer  :: Iwork(Liwrk)
    integer  :: Icntyp(Numgr)
    real(wp) :: Confun(Numgr,Nparm+1)
    real(wp) :: Error(Numgr+3)

    real(wp) :: ei , enor2 , enor3 , enorm
    integer  :: i , ilc22 , im1 , im2 , ioptth , ipt , l

    ! SET MACHINE AND PRECISION DEPENDENT CONSTANTS.
    ilc22 = iloc(22,Nparm,Numgr)
    ioptth = (Ioptn-(Ioptn/100000)*100000)/10000
    if ( ioptth<=0 ) then

        ! HERE IOPTTH=0, AND EACH CALL TO FNSET WILL COMPUTE FUNCTION VALUES
        ! FOR ONLY ONE CONSTRAINT.
        do i = 1 , Numgr
            ipt = i
            if ( Icnuse<=0 ) then
                ! HERE ICNUSE=0 SO WE WILL ACCEPT AND USE THE ICNTYP(I) COMPUTED BY
                ! FNSET.
                ! CALL FNSET WITH INDFN=0 TO COMPUTE CONFUN(I,1) AND ICNTYP(I).
                call me%fnset(Nparm,Numgr,Pttbl,Iptb,Indm,Param,ipt,0,Icntyp,Confun)
                ! HERE ICNUSE=1 AND THE ICNTYP CARRIED INTO ERCMP1 WILL OVERRIDE THAT
                ! COMPUTED BY FNSET.  THIS WILL ALSO BE TRUE IN ALL SUBROUTINES OTHER
                ! THAN CONMAX.  IF ICNTYP(I)=0 WE WILL SET ERROR(I)=0.0 AND WILL NOT
                ! NEED TO CALL FNSET.
            else if ( Icntyp(i)/=0 ) then
                ! CALL FNSET WITH INDFN=0 TO COMPUTE CONFUN(I,1).  THE COMPUTED KCNTYP
                ! WILL NOT BE USED.
                call me%fnset(Nparm,Numgr,Pttbl,Iptb,Indm,Param,ipt,0,Iwork(ilc22),Confun)
            else
                Error(i) = zero
                cycle
            end if
            call compute_error(i)
        end do

    else

        fnset: block
            ! HERE IOPTTH=1 AND A SINGLE CALL TO FNSET WITH INDFN=0 WILL COMPUTE
            ! CONFUN(.,1) AND (IF ICNUSE=0) ICNTYP(.).
            if ( Icnuse<=0 ) then
                ! HERE IOPTTH=1 AND ICNUSE=0, AND WE SET IPT=0 TO TELL FNSET TO COMPUTE
                ! ALL CONSTRAINTS (SINCE WE WANT TO BE SURE THAT ALL OF ICNTYP IS
                ! COMPUTED).  NOTE THAT IF INSTEAD WE HAD IOPTTH=0, THEN IPT WOULD
                ! BE POSITIVE AT EACH FNSET CALL, TELLING FNSET TO COMPUTE CONSTRAINT
                ! IPT ONLY.
                ipt = 0
                call me%fnset(Nparm,Numgr,Pttbl,Iptb,Indm,Param,ipt,0,Icntyp,Confun)
            else
                ! HERE IOPTTH=1 AND ICNUSE=1, AND IF IPHSE IS NEGATIVE WE SET IPT=-1
                ! TO TELL FNSET THAT ONLY STANDARD CONSTRAINTS NEED TO BE COMPUTED.
                ! IF IPHSE=0 HERE WE CHECK TO SEE IF ANY ICNTYP(L) IS POSITIVE FOR
                ! L=1,...,NUMGR, AND IF SO WE SET IPT=0 TO TELL FNSET TO COMPUTE ALL
                ! CONSTRAINTS, WHILE OTHERWISE WE SET IPT=-1.
                if ( Iphse>=0 ) then
                    do l = 1 , Numgr
                        if ( Icntyp(l)>0 ) then
                            ipt = 0
                            call me%fnset(Nparm,Numgr,Pttbl,Iptb,Indm,Param,ipt,0,Iwork(ilc22),Confun)
                            exit fnset
                        end if
                    end do
                end if
                ipt = -1
                call me%fnset(Nparm,Numgr,Pttbl,Iptb,Indm,Param,ipt,0,Iwork(ilc22),Confun)
            end if
        end block fnset

        ! Compute error as above.
        do i = 1 , Numgr
            call compute_error(i)
        end do

    end if

    ! HAVING FINISHED COMPUTING ERROR(I) AND (IF ICNUSE=0) ICNTYP(I) FOR
    ! I=1,...,NUMGR WE NOW COMPUTE THE ERROR NORMS.
    ! WE ALSO COMPUTE THE INDEX IPMAX OF THE CONSTRAINT WHERE THE PRIMARY
    ! (I.E. TYPE 1 OR TYPE 2) ERROR NORM OCCURS AND THE INDEX ISMAX OF THE
    ! CONSTRAINT WHERE THE STANDARD (I.E. TYPE -1 OR TYPE -2) ERROR NORM
    ! OCCURS.
    ! FIRST INITIALIZE THE INDICATORS AND ERROR NORMS.
    im1 = 0
    im2 = 0
    Ipmax = 0
    Ismax = 0
    enorm = zero
    enor2 = zero
    enor3 = zero

    do i = 1 , Numgr
        ei = Error(i)
        if ( Icntyp(i)<0 ) then
            if ( Icntyp(i)+1<0 ) then
                ! HERE ICNTYP(I)=-2 AND WE DO AS ABOVE EXCEPT WITH IM2 AND ENOR3.
                if ( im2>0 ) then
                    if ( ei<=enor3 ) cycle
                end if
                im2 = i
                enor3 = ei
            else
                ! HERE ICNTYP(I)=-1 AND WE DO AS ABOVE EXCEPT WITH IM1 AND ENOR2.
                if ( im1>0 ) then
                    if ( ei<=enor2 ) cycle
                end if
                im1 = i
                enor2 = ei
            end if
        else if ( Icntyp(i)/=0 ) then
            ! HERE ICNTYP(I) > 0.  IF ICNTYP(I)=2 REPLACE EI BY ABS(EI).  IF THIS
            ! IS THE FIRST I FOUND WITH ICNTYP(I) > 0 WE RESET IPMAX TO I AND PUT
            ! EI IN ENORM, AND OTHERWISE RESET IPMAX AND PUT EI IN ENORM IF AND ONLY
            ! IF EI IS BIGGER THAN THE VALUES FOUND SO FAR.
            if ( Icntyp(i)>1 ) ei = abs(ei)
            if ( Ipmax>0 ) then
                if ( ei<=enorm ) cycle
            end if
            Ipmax = i
            enorm = ei
        end if
    end do

    ! NOW RESET ISMAX IF THERE ARE ANY STANDARD CONSTRAINTS.
    if ( im1<=0 ) then
        ! HERE THERE ARE STANDARD NONLINEAR CONSTRAINTS BUT NO STANDARD LINEAR
        ! CONSTRAINTS.
        if ( im2>0 ) Ismax = im2
    else if ( im2<=0 ) then
        ! HERE THERE ARE STANDARD LINEAR CONSTRAINTS BUT NO STANDARD NONLINEAR
        ! CONSTRAINTS.
        Ismax = im1
        ! HERE THERE ARE BOTH STANDARD LINEAR CONSTRAINTS AND STANDARD NONLINEAR
        ! CONSTRAINTS.
    else if ( enor3<enor2 ) then
        Ismax = im1
    else
        Ismax = im2
    end if

    ! SET ERROR(NUMGR+1) THROUGH ERROR(NUMGR+3).
    Error(Numgr+1) = enorm
    Error(Numgr+2) = enor2
    Error(Numgr+3) = enor3

    contains

    subroutine compute_error(i)
    !! Set error(i)=0.0, or confun(i,1), or fun(i) - confun(i,1) according as
    !! icntyp(i) is 0, or -2, -1, 1, or 2.
    implicit none
    integer,intent(in) :: i
    if ( Icntyp(i)==0 ) then
        Error(i) = zero
    else if ( Icntyp(i)>1 ) then
        Error(i) = Fun(i) - Confun(i,1)
    else
        Error(i) = Confun(i,1)
    end if
    end subroutine compute_error

    end subroutine ercmp1
!********************************************************************************

!********************************************************************************
!>
!
    subroutine rkcon(me,Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,   &
                     Tolcon,Rchin,Iter,Irk,Ityp2,Ityp1,Itypm1,Itypm2, &
                     Icntyp,Projct,Rchdwn,Nstep,Iphse,Enchg,Enc1,Pmat,&
                     Funtbl,Iwork,Liwrk,Work,Lwrk,Iact,Actdif,Parprj, &
                     Parser,Xrk,Err1,Confun,Isucc,Param,Error)

    implicit none

    class(conmax_solver),intent(inout) :: me
    real(wp) :: Actdif , Confun , conup , emin , emin1 , Enc1 ,   &
                Enchg , enorm , Err1 , Error , Fun , Funtbl , &
                Param , Parprj , Parser , pe , Pmat
    real(wp) :: prjbig , prjlim , Projct , prosea , Pttbl , qt ,   &
                quots , Rchdwn , Rchin , s , ss ,    &
                steplm , Tolcon
    real(wp)  :: unit , wdist , Work , Xrk
    integer i , Iact , Icntyp , icorct , ifrkpr , Ifun , ilc06 ,      &
            ilc10 , ilc15 , ilc21 , ilc22 , ilc24 , ilc27 , ilc30 ,   &
            ilc31 , ilc33 , ilc35 , ilc36 , ilc37 , ilc38
    integer ilc40 , ilc43 , ilc48 , Indm , Ioptn , ioptth ,    &
            Iphse , ipmax , ipt , Iptb , Irk , ismax , Isucc , Iter , &
            Ityp1 , Ityp2 , Itypm1 , Itypm2 , iwarn
    integer Iwork , j , jflag , l , limfl , Liwrk , Lwrk , mactrk ,   &
            ncor , nfail , nmaj , nmin , npar1 , Nparm , nsrch ,      &
            Nstep , Numgr

    dimension Fun(Ifun) , Pttbl(Iptb,Indm) , Icntyp(Numgr) ,          &
              Param(Nparm) , Error(Numgr+3) , Pmat(Nparm+1,Numgr) ,   &
              Funtbl(Numgr,Nparm+1) , Iwork(Liwrk) , Work(Lwrk) ,     &
              Iact(Numgr) , Actdif(Numgr) , Parprj(Nparm) ,           &
              Parser(Nparm) , Xrk(Nparm+1) , Err1(Numgr+3) ,          &
              Confun(Numgr,Nparm+1)

    real(wp),parameter :: qthi = (one+two)/four
    real(wp),parameter :: qtlo = one/four
    real(wp),parameter :: tol1 = ten*ten*spcmn
    real(wp),parameter :: tol2 = ten*spcmn
    real(wp),parameter :: prden = sqrt(sqrt(spcmn))

    ioptth = (Ioptn-(Ioptn/100000)*100000)/10000
    steplm = Tolcon/ten
    ilc06 = iloc(6,Nparm,Numgr)
    ilc10 = iloc(10,Nparm,Numgr)
    ilc15 = iloc(15,Nparm,Numgr)
    ilc21 = iloc(21,Nparm,Numgr)
    ilc22 = iloc(22,Nparm,Numgr)
    ilc24 = iloc(24,Nparm,Numgr)
    ilc27 = iloc(27,Nparm,Numgr)
    ilc30 = iloc(30,Nparm,Numgr)
    ilc31 = iloc(31,Nparm,Numgr)
    ilc33 = iloc(33,Nparm,Numgr)
    ilc35 = iloc(35,Nparm,Numgr)
    ilc36 = iloc(36,Nparm,Numgr)
    ilc37 = iloc(37,Nparm,Numgr)
    ilc38 = iloc(38,Nparm,Numgr)
    ilc40 = iloc(40,Nparm,Numgr)
    ilc43 = iloc(43,Nparm,Numgr)
    ilc48 = iloc(48,Nparm,Numgr)
    Isucc = 0
    iwarn = 0
    nfail = 0
    conup = one

    ! LIMFL IS A SAFETY VALVE TO CATCH BLUNDERS; WE SET IT HIGH ENOUGH
    ! THAT IT WILL NOMALLY NOT COME INTO PLAY.
    limfl = 20
    enorm = Error(Numgr+1)
    npar1 = Nparm + 1
    prjbig = one/spcmn
    if ( Ityp2>0 ) prjbig = enorm

    ! THE NEXT GROUP OF STATEMENTS SETS AN INITIAL VALUE FOR PROJCT.

    if ( Iter>0 ) then
        if ( Irk<=1 ) then
            ! HERE ITER > 0 AND IRK=1, AND WE BUILD ON THE PREVIOUS SUCCESSFUL
            ! RK ITERATION, WHICH REDUCED THE ERROR NORM.  COMPUTE THE RATIO QT,
            ! WHICH WOULD BE 1.0 IF RUNGE-KUTTA WERE EXACT AND NO CORRECTION STEP
            ! WERE NEEDED.
            qt = -Enc1/Projct
            if ( qt>=qthi ) then
                ! HERE QT >= QTHI, SO WE INCREASE PROJCT BY A FACTOR OF 2.0.
                Projct = two*Projct
                ! IF QTLO < QT < QTHI WE LEAVE PROJCT THE SAME, WHILE IF QT <=
                ! QTLO WE DIVIDE PROJCT BY 4.0.
            else if ( qt<=qtlo ) then
                Projct = Projct/four
            end if
            goto 100
        end if
    end if

    ! HERE ITER=0, OR ELSE ITER > 0 AND IRK=2, AND WE INITIALIZE PROJCT.
    if ( Iphse+1<0 ) then
        ! HERE ITER=0 OR IRK=2, AND IPHSE=-2, SO WE ARE ATTEMPTING TO GAIN TYPE -2
        ! FEASIBILITY, AND WE SET THE INITIAL PROJCT TO ENOR3,
        ! WHICH WILL BE > TOLCON.  NOTE THAT ENOR3 IS NOW IN ERROR(NUMGR+1).
        Projct = enorm
    else if ( Iphse+1/=0 ) then

        ! HERE ITER=0 OR IRK=2, AND IPHSE=0, SO WE ARE IN THE MAIN ITERATIONS,
        ! AND WE FIRST TRY PROJCT=1.0.
        Projct = one

        ! CHECK TO SEE WHETHER ABS(ENORM) IS VERY
        ! LARGE RELATIVE TO THE INITIAL PROJCT.  IF ABS(ENORM) >
        ! PROJCT/PRDEN, WE REPLACE THE INITIAL PROJCT BY PRDEN*ABS(ENORM)
        ! SO THAT IF WE ARE SUCCESSFUL IN REDUCING ENORM TO ENORM - PROJCT,
        ! THIS QUANTITY WILL DIFFER FROM ENORM IN AT LEAST SOME SIGNIFICANT
        ! DIGITS AND THE PROGRAM WILL HAVE A CHANCE TO CONTINUE.
        pe = prden*abs(enorm)
        if ( pe>Projct ) Projct = pe
        ! IF ITYP2 > 0 WE REDUCE THE INITIAL PROJCT TO ENORM (IF NECESSARY),
        ! WHICH WILL BE THE GREATEST DECREASE IN ENORM WE CAN HOPE FOR SINCE
        ! THERE WILL BE TYPE 2 CONSTRAINTS.
        if ( Ityp2>0 ) then
            if ( enorm<Projct ) Projct = enorm
        end if
    end if

    ! WE DO NOT WISH FOR PROJCT TO BE SET BELOW 100.0*SPCMN
    if ( Projct<ten*ten*spcmn ) Projct = ten*ten*spcmn

    ! WE DO NOT WANT PROJCT TO BE BIGGER THAN PRJBIG OR SMALLER THAN
    ! STEPLM.
100 if ( Projct>prjbig ) Projct = prjbig
    if ( Projct<steplm ) then
        iwarn = 1
        Projct = steplm
    end if

    ! CALL RKSACT TO PUT THE (SIGNED) INDICES OF THE ACTIVE CONSTRAINTS IN
    ! IACT AND THEIR NUMBER IN MACTRK.
    call rksact(Ioptn,Numgr,Icntyp,Rchdwn,Rchin,conup,Projct,Error, &
                mactrk,Actdif,Iact)

    ! SET UNIT FOR USE IN RCHMOD.  UNIT WILL BE THE VALUE OF PROJCT WHEN
    ! RKSACT WAS LAST CALLED.
    unit = Projct

    ! CALL PMTST TO SET UP PMAT.
    call me%pmtst(Ioptn,Numgr,Nparm,Param,Icntyp,mactrk,Iact,Pttbl,Iptb, &
                  Indm,Actdif,Iphse,Iwork,Liwrk,Work,Lwrk,Confun,Pmat)

    ! COPY PMAT TRANSPOSE INTO FUNTBL FOR POSSIBLE LATER USE.
    do j = 1 , npar1
        do i = 1 , mactrk
            Funtbl(i,j) = Pmat(j,i)
        end do
    end do

    ! STATEMENTS ABOVE THIS POINT WILL NOT BE EXECUTED AGAIN IN THIS CALL
    ! TO RKCON.

    ! INCREMENT NFAIL, WHICH COUNTS THE NUMBER OF WOLFE CALLS IN THIS CALL TO
    ! RKCON.
200 nfail = nfail + 1

    ! CALL WOLFE WITH ISTRT=0 TO SOLVE THE LEAST DISTANCE QP PROBLEM FROM
    ! SCRATCH.
    call wolfe(Nparm,mactrk,Pmat,0,s,ncor,Iwork(ilc15),Iwork,Liwrk,   &
               Work,Lwrk,Work(ilc33),Work(ilc06),Work(ilc31),         &
               Work(ilc30),Nparm,Numgr,Work(ilc40),Work(ilc36),wdist, &
               nmaj,nmin,jflag)

    ! IF WOLFE FAILS, WE MAY TRY AGAIN WITH A SMALLER PROJCT.

    if ( jflag<=0 ) then
        ! END OF GROUP OF STATEMENTS TO REDUCE PROJCT (IF POSSIBLE) TO HANDLE
        ! A FAILURE OF SOME KIND.

        ! DO AN RK STEP.
        call me%rkpar(Ioptn,Numgr,Nparm,Icntyp,mactrk,Iact,Actdif,Projct, &
                      Param,Fun,Ifun,Pttbl,Iptb,Indm,Work(ilc36),Pmat,    &
                      ncor,s,Itypm1,Itypm2,unit,Tolcon,Rchin,Nstep,Error, &
                      Iphse,Iwork,Liwrk,Work,Lwrk,Confun,Work(ilc37),     &
                      Work(ilc38),Work(ilc43),Parprj,ifrkpr)
        if ( ifrkpr<=0 ) then
            ! HERE RKPAR SUCCEEDED.  IF THERE ARE ANY STANDARD CONSTRAINTS WE CALL
            ! CORRCT TO MOVE BACK INTO THE FEASIBLE REGION IF NECESSARY.
            if ( Itypm1+Itypm2<=0 ) goto 500
            call me%corrct(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,     &
                           Icntyp,unit,Tolcon,Rchin,Error,mactrk,Iact,     &
                           Projct,Iphse,Iwork,Liwrk,Work,Lwrk,Work(ilc27), &
                           Err1,Work(ilc10),Pmat,Confun,Work(ilc48),       &
                           Iwork(ilc21),Parprj,icorct)
            if ( icorct<=0 ) goto 500
        end if
    end if

    ! THE NEXT GROUP OF STATEMENTS IS TO REDUCE PROJCT (IF POSSIBLE) IN CASE
    ! OF A FAILURE OF SOME KIND.

300 if ( nfail<limfl ) then
        ! PREPARE TO TRY ANOTHER ITERATION IN RKCON BY
        ! REDUCING PROJCT, AND MAKING SURE PROJCT IS NOT TOO SMALL.
        Projct = Projct/(four+four)
        if ( Projct>=steplm ) goto 400
        if ( iwarn<=0 ) then
            iwarn = 1
            Projct = steplm
            goto 400
        end if
    end if

    ! HERE RKCON COULD NOT REDUCE THE ERROR NORM AND WE RETURN WITH THE
    ! WARNING ISUCC=1.
    Isucc = 1
    return

    ! NOW RESET ACTDIF FOR THIS PROJCT.
400 do l = 1 , mactrk
        i = abs(Iact(l))
        if ( Icntyp(i)<0 ) then
            if ( Icntyp(i)+1<0 ) then
                ! HERE WE HAVE AN ACTIVE TYPE -2 CONSTRAINT, AND WE SET ACTDIF(L)=
                ! MIN (CONUP, ERROR(I)/PROJCT).
                Actdif(l) = Error(i)/Projct
                if ( Actdif(l)>conup ) Actdif(l) = conup
            else
                ! HERE WE HAVE AN ACTIVE TYPE -1 CONSTRAINT.
                Actdif(l) = Error(i)/Projct
            end if
        else if ( Icntyp(i)==0 ) then
            ! ICNTYP(I)=0 SHOULD NOT OCCUR HERE SINCE CONSTRAINT I WAS DECLARED
            ! TO BE ACTIVE IN RKSACT, BUT WE ACCOUNT FOR IT ANYWAY AS A PRECAUTION.
            Actdif(i) = zero
        else if ( Icntyp(i)<=1 ) then
            ! HERE WE HAVE AN ACTIVE TYPE 1 CONSTRAINT.
            Actdif(l) = one + (Error(i)-enorm)/Projct
        else
            ! HERE WE HAVE AN ACTIVE TYPE 2 CONSTRAINT.
            Actdif(l) = one + (abs(Error(i))-enorm)/Projct
        end if
    end do

    ! COPY THE FIRST NPARM ROWS OF PMAT FROM OLD PMAT TRANSPOSE STORED
    ! IN FUNTBL, THEN APPEND ACTDIF AS THE LAST ROW.
    do j = 1 , mactrk
        do i = 1 , Nparm
            Pmat(i,j) = Funtbl(j,i)
        end do
        Pmat(npar1,j) = Actdif(j)
    end do
    goto 200

    ! PUT THE SEARCH DIRECTION VECTOR PARPRJ - PARAM INTO XRK.
500 do j = 1 , Nparm
        Xrk(j) = Parprj(j) - Param(j)
    end do

    ! CALL SEARSL TO DO A LINE SEARCH IN DIRECTION XRK AND PUT THE RESULTING
    ! VECTOR IN PARSER.  START WITH A PROJECTION FACTOR PROSEA=1.0.
    ! PARPRJ WILL BE USED TEMPORARILY AS A WORK VECTOR IN SEARSL.
    prosea = one

    ! WE NOW WISH TO DETERMINE PRJLIM = THE SMALLER OF 1.0/SPCMN AND
    ! THE LARGEST VALUE OF PROSEA FOR WHICH THE LINEAR STANDARD CONSTRAINTS
    ! ARE SATISFIED FOR THE PARAMETER VECTOR PARAM+PROSEA*XRK.  THIS
    ! WILL GIVE AN UPPER BOUND FOR LINE SEARCHING.  NOTE THAT IN
    ! THEORY WE SHOULD HAVE PRJLIM >= 1.0 SINCE THE LINEAR STANDARD
    ! CONSTRAINTS SHOULD BE SATISFIED FOR PROSEA=0.0 AND PROSEA=1.0, BUT
    ! ROUNDOFF ERROR COULD AFFECT THIS A LITTLE.  IF THERE ARE NO
    ! LINEAR STANDARD CONSTRAINTS, WE SET PRJLIM=1.0/SPCMN.
    prjlim = one/spcmn
    !*****INSERT TO MAKE SEARCHING LESS VIOLENT.
    !     PRJLIM=TWO
    !*****END INSERT
    if ( Itypm1>0 ) then
        ! HERE WE HAVE AT LEAST ONE TYPE -1 CONSTRAINT, AND IF IOPTTH=1 WE
        ! CALL DERST TO PUT ALL THE STANDARD CONSTRAINT VALUES AND GRADIENTS
        ! INTO CONFUN(.,.).
        if ( ioptth>0 ) then
            ! WE SET IPT=-1 TO TELL DERST TO COMPUTE STANDARD CONSTRAINTS ONLY.
            ipt = -1
            call me%derst(Ioptn,Nparm,Numgr,Pttbl,Iptb,Indm,Param,ipt, &
                          Work(ilc24),Work(ilc35),Iwork(ilc22),Confun)
        end if
        do i = 1 , Numgr
            if ( Icntyp(i)+1==0 ) then
                ipt = i
                ! HERE WE HAVE A TYPE -1 CONSTRAINT AND IF IOPTTH=0 WE CALL DERST
                ! TO PUT THE CONSTRAINT VALUE AND GRADIENT INTO CONFUN(IPT,.).
                if ( ioptth<=0 ) call me%derst(Ioptn,Nparm,Numgr,Pttbl,Iptb,&
                                               Indm,Param,ipt,Work(ilc24),Work(ilc35),&
                                               Iwork(ilc22),Confun)

                ! WE WISH TO HAVE SUMMATION (CONFUN(IPT,J+1)*(PARAM(J)+PROSEA*XRK(J)))
                ! + C(IPT) <= 0.0 FOR IPT=1,...,NUMGR, ICNTYP(IPT) = -1,
                ! WHERE THE IPTTH CONSTRAINT APPLIED TO PARAM SAYS
                ! SUMMATION (CONFUN(IPT,J+1)*PARAM(J)) + C(IPT) <= 0.0, SO C(IPT) IS
                ! THE CONSTANT TERM IN THE LEFT SIDE OF LINEAR CONSTRAINT IPT.
                ! THUS FOR I=1PT,...,NUMGR, ICNTYP(IPT) = -1, WE WANT PRJLIM*SS <=
                ! SSS, WHERE SS = SUMMATION (CONFUN(IPT,J+1)*XRK(J)) AND SSS = -C(IPT) -
                ! SUMMATION (CONFUN(IPT,J+1)*PARAM(J)) = -CONFUN(IPT,1).
                ss = zero
                do j = 1 , Nparm
                    ss = ss + Confun(i,j+1)*Xrk(j)
                end do
                ! IF SS < 10.0*SPCMN THIS CONSTRAINT WILL NOT PUT A SIGNIFICANT
                ! RESTRICTION ON PROSEA.
                if ( ss>=tol2 ) then
                    ! HERE SS >= 10.0*SPCMN AND WE COMPARE SSS/SS AGIANST PRJLIM.
                    quots = -Confun(i,1)/ss
                    if ( prjlim>quots ) prjlim = quots
                end if
            end if
        end do
    end if

    call me%searsl(Ioptn,Numgr,Nparm,prjlim,tol1,Xrk,Fun,Ifun,Pttbl,Iptb,&
                   Indm,Param,Error,Rchdwn,mactrk,Iact,Iphse,unit,Tolcon,&
                   Rchin,Itypm1,Itypm2,Iwork,Liwrk,Work,Lwrk,Err1,Parprj,&
                   prosea,emin,emin1,Parser,nsrch)

    ! COMPUTE THE PRINCIPAL ERROR NORM CHANGE ENCHG.  ALSO COMPUTE ENC1, THE
    ! CHANGE IN THE PRINCIPAL ERROR NORM WITHOUT THE LINE SEARCH.
    Enchg = emin - enorm
    Enc1 = emin1 - enorm

    ! IF WE OBTAINED MORE THAN A TOL1 REDUCTION IN ENORM WE UPDATE
    ! PARAM AND CALL ERCMP1 TO UPDATE ERROR, AND RETURN WITH ISUCC=0
    ! INDICATING SUCCESS.
    ! OTHERWISE WE CHECK TO SEE IF WE HAVE REACHED THE RKCON ITERATION
    ! LIMIT, AND IF SO WE RETURN WITH ISUCC=1, INDICATING FAILURE.
    if ( Enchg+tol1>=0 ) goto 300

    do j = 1 , Nparm
        Param(j) = Parser(j)
    end do
    call me%ercmp1(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,Param,1, &
                   Iphse,Iwork,Liwrk,Confun,Icntyp,ipmax,ismax,Error)

    end subroutine rkcon
!********************************************************************************

!********************************************************************************
!>
! This subroutine puts the (signed) indices of the mactrk
! active constraints in iact.  it also sets the right side vector
! actdif for the wolfe subproblem.

    subroutine rksact(Ioptn,Numgr,Icntyp,Rchdwn,Rchin,Conup,Projct, &
                      Error,Mactrk,Actdif,Iact)

    implicit none

    integer,intent(in) :: Ioptn
    integer,intent(in) :: Numgr
    real(wp),intent(in) :: Rchdwn
    real(wp),intent(in) :: Rchin
    real(wp),intent(in) :: Conup
    real(wp),intent(in) :: Projct
    real(wp),intent(in) :: Error(Numgr+3)
    real(wp) :: Actdif(Numgr)
    integer :: Iact(Numgr)
    integer,intent(in)  :: Icntyp(Numgr)
    integer,intent(out) :: Mactrk

    real(wp) :: elow , enorm , rchind
    integer :: i , l

    ! SET MACHINE AND PRECISION DEPENDENT CONSTANTS FOR RKSACT.
    enorm = Error(Numgr+1)
    elow = enorm - Rchdwn*Projct
    rchind = Rchin*Projct

    ! DETERMINE THE NUMBER MACTRK OF ACTIVE CONSTRAINTS, THEIR INDICATOR
    ! IACT, AND THE VECTOR ACTDIF OF RIGHT SIDES FOR THE WOLFE SUBPROBLEM.
    l = 0
    do i = 1 , Numgr
        if ( Icntyp(i)<0 ) then
            if ( Icntyp(i)+1>=0 ) then
                ! HERE WE HAVE A TYPE -1 CONSTRAINT, WHICH WILL AUTOMATICALLY BE
                ! DECLARED TO BE ACTIVE.
                l = l + 1
                Iact(l) = i
                Actdif(l) = Error(i)/Projct
                ! HERE WE HAVE A TYPE -2 CONSTRAINT, WHICH WILL BE DECLARED TO BE
                ! ACTIVE IFF ERROR(I) >= -RCHIND.
            else if ( Error(i)+rchind>=0 ) then
                ! HERE WE HAVE AN ACTIVE TYPE -2 CONSTRAINT, AND WE SET ACTDIF(L)=
                ! MIN (CONUP, ERROR(I)/PROJCT).
                l = l + 1
                Iact(l) = i
                Actdif(l) = Error(i)/Projct
                if ( Actdif(l)>Conup ) Actdif(l) = Conup
            end if
        else if ( Icntyp(i)/=0 ) then
            if ( Icntyp(i)>1 ) then
                ! HERE WE HAVE A TYPE 2 CONSTRAINT.
                if ( Error(i)<0 ) then
                    if ( -Error(i)>=elow ) then
                        ! HERE WE HAVE A -ACTIVE TYPE 2 CONSTRAINT.
                        l = l + 1
                        Iact(l) = -i
                        Actdif(l) = one + (-Error(i)-enorm)/Projct
                    end if
                    cycle
                end if
            end if
            ! HERE WE HAVE A TYPE 1 CONSTRAINT, OR A TYPE 2 CONSTRAINT WITH
            ! ERROR(I) >= 0.0.
            if ( Error(i)>=elow ) then
                ! HERE WE HAVE AN ACTIVE TYPE 1 CONSTRAINT OR A +ACTIVE TYPE 2 CONSTRAINT.
                l = l + 1
                Iact(l) = i
                Actdif(l) = one + (Error(i)-enorm)/Projct
            end if
        end if
    end do

    Mactrk = l

    end subroutine rksact
!********************************************************************************

!********************************************************************************
!>
! This subroutine sets up the (nparm+1) by mactrk matrix pmat.
! for 1 <= j <= mactrk, the top nparm elements of column j of pmat
! will contain the negative of the gradient of active constraint j (if
! constraint j is of type 2, i.e. of the form abs(f(x) - f(parwrk,x))
! <= w, the left side will be treated as f(x) - f(parwrk,x) if this
! quantity is nonnegative and will be treated as f(parwrk,x) - f(x)
! otherwise). the (nparm+1)st row of pmat will contain actdif, the
! right side of the inequalities gradient.vector >= actdif.

    subroutine pmtst(me,Ioptn,Numgr,Nparm,Param,Icntyp,Mactrk,Iact,Pttbl, &
                     Iptb,Indm,Actdif,Iphse,Iwork,Liwrk,Work,Lwrk, &
                     Confun,Pmat)

    implicit none

    class(conmax_solver),intent(inout) :: me
    integer,intent(in) :: Indm
    integer,intent(in) :: Ioptn
    integer,intent(in) :: Iphse
    integer,intent(in) :: Iptb
    integer,intent(in) :: Liwrk
    integer,intent(in) :: Lwrk
    integer,intent(in) :: Mactrk
    integer,intent(in) :: Nparm
    integer,intent(in) :: Numgr
    integer :: Iwork(Liwrk)
    integer :: Iact(Numgr)
    integer :: Icntyp(Numgr)
    real(wp),intent(in) :: Actdif(Numgr)
    real(wp) :: Confun(Numgr,Nparm+1)
    real(wp) :: Param(Nparm)
    real(wp) :: Pmat(Nparm+1,Numgr)
    real(wp) :: Pttbl(Iptb,Indm)
    real(wp) :: Work(Lwrk)

    integer :: i , ii , ilc22 , ilc24 , ilc35 , ioptth , ipt , j , l
    integer :: npar1

    ! SET MACHINE AND PRECISION DEPENDENT CONSTANTS FOR PMTST.
    ilc22 = iloc(22,Nparm,Numgr)
    ilc24 = iloc(24,Nparm,Numgr)
    ilc35 = iloc(35,Nparm,Numgr)
    ioptth = (Ioptn-(Ioptn/100000)*100000)/10000
    npar1 = Nparm + 1

    if ( ioptth>0 ) then
        ! HERE IOPTTH=1 AND WE CALL DERST TO PUT GRADIENT VALUES INTO CONFUN.
        ! IF IPHSE < 0 OR NO ICNTYP(L) IS POSITIVE, SET IPT=-1 TO TELL DERST
        ! TO COMPUTE STANDARD CONSTRAINTS ONLY, WHILE OTHERWISE SET IPT=0 TO
        ! TELL DERST TO COMPUTE ALL CONSTRAINTS.
        ipt = -1
        if ( Iphse>=0 ) then
            do l = 1 , Numgr
                if ( Icntyp(l)>0 ) then
                    ipt = 0
                    exit
                end if
            end do
        end if
        call me%derst(Ioptn,Nparm,Numgr,Pttbl,Iptb,Indm,Param,ipt,Work(ilc24),&
                      Work(ilc35),Iwork(ilc22),Confun)
    end if

    do i = 1 , Mactrk
        ii = Iact(i)
        ipt = abs(ii)

        ! HERE IOPTTH=0 AND WE HAVE NOT YET PLACED THE GRADIENT IN CONFUN, SO WE
        ! CALL DERST TO DO SO NOW.  DERST WILL ALSO COMPUTE THE
        ! CONSTRAINT VALUES, WHICH WILL NOT BE NEEDED HERE, BUT EXPECTING USERS TO
        ! WRITE FNSET SO THAT GRADIENT CALCULATIONS WILL NOT NEED FUNCTION VALUE
        ! CALCULATION RESULTS WOULD BE TOO MUCH OF A PROGRAMMING TRAP.
        if ( ioptth<=0 ) call me%derst(Ioptn,Nparm,Numgr,Pttbl,Iptb,Indm, &
                                       Param,ipt,Work(ilc24),Work(ilc35), &
                                       Iwork(ilc22),Confun)

        ! NOW THE GRADIENT FOR CONSTRAINT IPT IS IN CONFUN(IPT,.), AND WE PUT IT
        ! OR ITS NEGATIVE INTO PMAT.
        ! IF ICNTYP(IPT) <= 1 WE PROCEED AS IF WE HAD A -ACTIVE CONSTRAINT IN
        ! THE ICNTYP(IPT)=2 CASE.  IN ALL CASES WE PUT THE NEGATIVE OF THE
        ! CONSTRAINT GRADIENT INTO COLUMN I OF PMAT.
        if ( Icntyp(ipt)>1 ) then
            ! HERE ICNTYP(IPT)=2.
            if ( ii>0 ) then
                ! HERE WE HAVE A +ACTIVE CONSTRAINT AT POINT IPT.
                ! THE CONSTRAINT GRADIENT IS IN -CONFUN(IPT,.) SINCE THE LEFT SIDE OF
                ! CONSTRAINT I IS F(X)-F(PARWRK,X) AND DERST COMPUTES THE
                ! GRADIENT OF F(PARWRK,X).  THUS WE PUT CONFUN(IPT,.) IN COLUMN I OF PMAT.
                do j = 1 , Nparm
                    Pmat(j,i) = Confun(ipt,j+1)
                end do
                cycle
            end if
        end if

        ! HERE WE HAVE A -ACTIVE TYPE 2 CONSTRAINT AT POINT -II OR AN ACTIVE
        ! CONSTRAINT OF TYPE -2, -1, OR 1 AT POINT II.
        do j = 1 , Nparm
            Pmat(j,i) = -Confun(ipt,j+1)
        end do
    end do

    ! PUT ACTDIF IN THE LAST ROW OF PMAT.
    do i = 1 , Mactrk
        Pmat(npar1,i) = Actdif(i)
    end do

    end subroutine pmtst
!********************************************************************************

!********************************************************************************
!>
!  This subroutine computes a parameter vector parprj using fourth
!  order runge kutta with h=-projct.  h is negative since we want
!  to approximate the parameters resulting from decreasing w by
!  abs(h).  if we do nstep steps then h=-projct/nstep.

    subroutine rkpar(me,Ioptn,Numgr,Nparm,Icntyp,Mactrk,Iact,Actdif,  &
                     Projct,Param,Fun,Ifun,Pttbl,Iptb,Indm,Vder,Pmat, &
                     Ncor,s,Itypm1,Itypm2,Unit,Tolcon,Rchin,Nstep,    &
                     Error,Iphse,Iwork,Liwrk,Work,Lwrk,Confun,Vdern,  &
                     Vders,Wvec,Parprj,Ifrkpr)

    implicit none

    class(conmax_solver),intent(inout) :: me

    integer,intent(in)  :: Ifun
    integer,intent(in)  :: Indm
    integer,intent(in)  :: Iptb
    integer,intent(in)  :: Liwrk
    integer ,intent(in) :: Lwrk
    integer,intent(in)  :: Nparm
    integer,intent(in)  :: Numgr
    integer  :: Ifrkpr
    integer  :: Ioptn
    integer  :: Iphse
    integer  :: Itypm1
    integer  :: Itypm2
    integer  :: Mactrk
    integer  :: Ncor
    integer  :: Nstep
    real(wp) :: Projct
    real(wp) :: Rchin
    real(wp) :: s
    real(wp) :: Tolcon
    real(wp) :: Unit
    real(wp) :: Actdif(Numgr)
    real(wp) :: Confun(Numgr,Nparm+1)
    real(wp) :: Error(Numgr+3)
    real(wp) :: Fun(Ifun)
    integer  :: Iact(Numgr)
    integer  :: Icntyp(Numgr)
    integer  :: Iwork(Liwrk)
    real(wp) :: Param(Nparm)
    real(wp) :: Parprj(Nparm)
    real(wp) :: Pmat(Nparm+1,Numgr)
    real(wp) :: Pttbl(Iptb,Indm)
    real(wp) :: Vder(Nparm)
    real(wp) :: Vdern(Nparm)
    real(wp) :: Vders(Nparm)
    real(wp) :: Work(Lwrk)
    real(wp) :: Wvec(Nparm)

    real(wp) :: p6, wdist, proj1
    integer :: icorct , ilc06 , ilc10 , ilc11 , ilc15 , ilc21 , ilc27 , &
               ilc30 , ilc31 , ilc33 , ilc40 , ilc48 , j , jflag , nmaj , &
               nmin , npar1 , nstcnt

    ! SET MACHINE AND PRECISION DEPENDENT CONSTANTS FOR RKPAR.
    ilc06 = iloc(6,Nparm,Numgr)
    ilc10 = iloc(10,Nparm,Numgr)
    ilc11 = iloc(11,Nparm,Numgr)
    ilc15 = iloc(15,Nparm,Numgr)
    ilc21 = iloc(21,Nparm,Numgr)
    ilc27 = iloc(27,Nparm,Numgr)
    ilc30 = iloc(30,Nparm,Numgr)
    ilc31 = iloc(31,Nparm,Numgr)
    ilc33 = iloc(33,Nparm,Numgr)
    ilc40 = iloc(40,Nparm,Numgr)
    ilc48 = iloc(48,Nparm,Numgr)

    ! IFRKPR=0 IS A SIGNAL THAT THE SUBROUTINE OPERATED NORMALLY.
    Ifrkpr = 0
    proj1 = Projct/Nstep
    p6 = proj1/(two+two+two)
    npar1 = Nparm + 1
    nstcnt = 1
    ! PARPRJ WILL BE USED AS THE BASE POINT FOR THE NEXT RK STEP DURING THE
    ! OPERATION OF THIS SUBROUTINE.
    do j = 1 , Nparm
        Parprj(j) = Param(j)
        Vdern(j) = Vder(j)
    end do

    main : block

        do
            ! NOTE THAT HERE H*VDERN IS THE K1 OF THE USUAL RUNGE-KUTTA FORMULAE.
            ! SET THE WORK VECTOR WVEC = PARPRJ-PROJ1*VDERN/2.0, THEN CALL PMTST
            ! AND WOLFE TO GET THE VECTOR (AGAIN CALLED VDERN) OF DERIVATIVE VALUES.
            ! THEN H*VDERN WILL BE THE K2 OF THE USUAL RUNGE-KUTTA FORMULAE.
            ! WE WILL ACCUMULATE K1/H + 2.0*K2/H + 2.0*K3/H IN VDERS, AND ADD IN
            ! K4/H AT THE END.
            do j = 1 , Nparm
                Vders(j) = Vdern(j)
                Wvec(j) = Parprj(j) - proj1*Vdern(j)/two
            end do
            ! IF THERE ARE ANY STANDARD CONSTRAINTS, WE CORRECT BACK INTO THE
            ! FEASIBLE REGION IF POSSIBLE BEFORE CALLING PMTST.
            if ( Itypm1+Itypm2>0 ) then
                call me%corrct(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,Icntyp, &
                               Unit,Tolcon,Rchin,Error,Mactrk,Iact,Projct,Iphse,  &
                               Iwork,Liwrk,Work,Lwrk,Work(ilc27),Work(ilc11),     &
                               Work(ilc10),Pmat,Confun,Work(ilc48),Iwork(ilc21),  &
                               Wvec,icorct)
                if ( icorct>0 ) exit main ! failure
            end if
            call me%pmtst(Ioptn,Numgr,Nparm,Wvec,Icntyp,Mactrk,Iact,Pttbl,Iptb, &
                          Indm,Actdif,Iphse,Iwork,Liwrk,Work,Lwrk,Confun,Pmat)
            call wolfe(Nparm,Mactrk,Pmat,1,s,Ncor,Iwork(ilc15),Iwork,Liwrk,  &
                       Work,Lwrk,Work(ilc33),Work(ilc06),Work(ilc31),        &
                       Work(ilc30),Nparm,Numgr,Work(ilc40),Vdern,wdist,nmaj, &
                       nmin,jflag)
            ! IF WOLFE FAILED, SO WILL THIS SUBROUTINE.
            if ( jflag<=0 ) then
                ! NOW VDERN REPRESENTS K2/H.  SET WVEC = PARPRJ-PROJ1*VDERN/2.0 AND
                ! COMPUTE THE NEW VDERN, WHICH WILL REPRESENT K3/H.
                do j = 1 , Nparm
                    Vders(j) = Vders(j) + two*Vdern(j)
                    Wvec(j) = Parprj(j) - proj1*Vdern(j)/two
                end do
                ! IF THERE ARE ANY STANDARD CONSTRAINTS, WE CORRECT BACK INTO THE
                ! FEASIBLE REGION IF POSSIBLE BEFORE CALLING PMTST.
                if ( Itypm1+Itypm2>0 ) exit main ! failure
                call me%corrct(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,Icntyp, &
                               Unit,Tolcon,Rchin,Error,Mactrk,Iact,Projct,Iphse,  &
                               Iwork,Liwrk,Work,Lwrk,Work(ilc27),Work(ilc11),     &
                               Work(ilc10),Pmat,Confun,Work(ilc48),Iwork(ilc21),  &
                               Wvec,icorct)
                if ( icorct>0 ) exit main ! failure
            end if

            call me%pmtst(Ioptn,Numgr,Nparm,Wvec,Icntyp,Mactrk,Iact,Pttbl,Iptb,  &
                          Indm,Actdif,Iphse,Iwork,Liwrk,Work,Lwrk,Confun,Pmat)
            call wolfe(Nparm,Mactrk,Pmat,1,s,Ncor,Iwork(ilc15),Iwork,Liwrk,   &
                       Work,Lwrk,Work(ilc33),Work(ilc06),Work(ilc31),         &
                       Work(ilc30),Nparm,Numgr,Work(ilc40),Vdern,wdist,nmaj,  &
                       nmin,jflag)
            if ( jflag>0 ) exit main ! failure

            ! NOW VDERN REPRESENTS K3/H.  SET WVEC = PARPRJ-PROJ1*VDERN AND
            ! COMPUTE THE NEW VDERN, WHICH WILL REPRESENT K4/H.
            do j = 1 , Nparm
                Vders(j) = Vders(j) + two*Vdern(j)
                Wvec(j) = Parprj(j) - proj1*Vdern(j)
            end do
            ! IF THERE ARE ANY STANDARD CONSTRAINTS, WE CORRECT BACK INTO THE
            ! FEASIBLE REGION IF POSSIBLE BEFORE CALLING PMTST.
            if ( Itypm1+Itypm2>0 ) then
                call me%corrct(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,Icntyp, &
                               Unit,Tolcon,Rchin,Error,Mactrk,Iact,Projct,Iphse,  &
                               Iwork,Liwrk,Work,Lwrk,Work(ilc27),Work(ilc11),     &
                               Work(ilc10),Pmat,Confun,Work(ilc48),Iwork(ilc21),  &
                               Wvec,icorct)
                if ( icorct>0 ) exit main ! failure
            end if
            call me%pmtst(Ioptn,Numgr,Nparm,Wvec,Icntyp,Mactrk,Iact,Pttbl,Iptb,  &
                          Indm,Actdif,Iphse,Iwork,Liwrk,Work,Lwrk,Confun,Pmat)
            call wolfe(Nparm,Mactrk,Pmat,1,s,Ncor,Iwork(ilc15),Iwork,Liwrk,   &
                       Work,Lwrk,Work(ilc33),Work(ilc06),Work(ilc31),         &
                       Work(ilc30),Nparm,Numgr,Work(ilc40),Vdern,wdist,nmaj,  &
                       nmin,jflag)
            if ( jflag>0 ) exit main ! failure

            ! NOW VDERN REPRESENTS K4/H, SO VDERS + VDERN WILL REPRESENT (K1 +
            ! 2.0*K2 + 2.0*K3 + K4)/H.  PUT THE NEW PARAMETER VECTOR IN PARPRJ.
            do j = 1 , Nparm
                Parprj(j) = Parprj(j) - p6*(Vders(j)+Vdern(j))
            end do
            if ( nstcnt<Nstep ) then
                ! HERE NSTCNT < NSTEP AND WE SET UP FOR THE NEXT RK STEP.
                ! AFTER WE HAVE DONE THIS STEP, VDERN WILL REPRESENT THE VDER1 FOR THE
                ! NEXT STEP.  PARPRJ ALREADY IS THE BASE POINT FOR THE NEXT STEP.
                nstcnt = nstcnt + 1
                ! IF THERE ARE ANY STANDARD CONSTRAINTS, WE CORRECT BACK INTO THE
                ! FEASIBLE REGION IF POSSIBLE BEFORE CALLING PMTST.
                if ( Itypm1+Itypm2>0 ) then
                    call me%corrct(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,     &
                                   Icntyp,Unit,Tolcon,Rchin,Error,Mactrk,Iact,     &
                                   Projct,Iphse,Iwork,Liwrk,Work,Lwrk,Work(ilc27), &
                                   Work(ilc11),Work(ilc10),Pmat,Confun,Work(ilc48),&
                                   Iwork(ilc21),Parprj,icorct)
                    if ( icorct>0 ) exit main ! failure
                end if
                call me%pmtst(Ioptn,Numgr,Nparm,Parprj,Icntyp,Mactrk,Iact,Pttbl,  &
                              Iptb,Indm,Actdif,Iphse,Iwork,Liwrk,Work,Lwrk,Confun,&
                              Pmat)
                call wolfe(Nparm,Mactrk,Pmat,1,s,Ncor,Iwork(ilc15),Iwork,Liwrk,&
                           Work,Lwrk,Work(ilc33),Work(ilc06),Work(ilc31),      &
                           Work(ilc30),Nparm,Numgr,Work(ilc40),Vdern,wdist,    &
                           nmaj,nmin,jflag)
                if ( jflag>0 ) exit main ! failure
            else
                exit ! done
            end if

        end do

        return ! success
    end block main

    ! failure
    Ifrkpr = 1
    ! WRITE(NWRIT,'(A)') '*****RKPAR HAS FAILED'

    end subroutine rkpar
!********************************************************************************

!********************************************************************************
!>
!  This subroutine determines whether parprj violates any type -2
!  or type -1 (i.e. standard) constraints by more than tolcon, and if
!  so it attempts to correct back to the feasible region.  if it is
!  successful it sets icorct=0 and replaces parprj by the corrected
!  vector.  if it is not successful it sets icorct=1 and leaves parprj
!  unchanged.  if no correction was needed it sets icorct=-1 and leaves
!  parprj unchanged.

    subroutine corrct(me,Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,  &
                      Icntyp,Unit,Tolcon,Rchin,Error,Mact,Iact,Projct,&
                      Iphse,Iwork,Liwrk,Work,Lwrk,Parwrk,Err1,Dvec,   &
                      Pmat,Confun,Zwork,Jcntyp,Parprj,Icorct)

    implicit none

    class(conmax_solver),intent(inout) :: me
    real(wp) :: Confun , Dvec , emin , eold , Err1 , Error , f1 ,   &
                Fun , gain , p1 , Parprj , Parwrk , Pmat , procor ,  &
                Projct , Pttbl , rchdwn , Rchin
    real(wp) :: s , Tolcon , Unit , wdist , Work , Zwork
    integer :: i , Iact , Icntyp , Icorct , Ifun , ilc06 , ilc16 ,       &
               ilc22 , ilc24 , ilc30 , ilc31 , ilc33 , ilc35 , ilc41 ,   &
               Indm , Ioptn , ioptth , Iphse , ipmax
    integer :: ipt , Iptb , ismax , isrcr , Iwork , j , Jcntyp , jflag , &
               k , l , Liwrk , Lwrk , Mact , ncor , newtit , newtlm ,    &
               nmaj , nmin , npar1 , Nparm
    integer :: Numgr

    dimension Fun(Ifun) , Pttbl(Iptb,Indm) , Icntyp(Numgr) ,          &
              Parprj(Nparm) , Parwrk(Nparm) , Err1(Numgr+3) ,         &
              Dvec(Nparm) , Pmat(Nparm+1,Numgr) , Jcntyp(Numgr) ,     &
              Confun(Numgr,Nparm+1) , Zwork(Nparm) , Error(Numgr+3) , &
              Iact(Numgr) , Iwork(Liwrk) , Work(Lwrk)

    ! SET MACHINE AND PRECISION DEPENDENT CONSTANTS.
    ilc06 = iloc(6,Nparm,Numgr)
    ilc16 = iloc(16,Nparm,Numgr)
    ilc22 = iloc(22,Nparm,Numgr)
    ilc24 = iloc(24,Nparm,Numgr)
    ilc30 = iloc(30,Nparm,Numgr)
    ilc31 = iloc(31,Nparm,Numgr)
    ilc33 = iloc(33,Nparm,Numgr)
    ilc35 = iloc(35,Nparm,Numgr)
    ilc41 = iloc(41,Nparm,Numgr)
    ioptth = (Ioptn-(Ioptn/100000)*100000)/10000
    npar1 = Nparm + 1
    newtit = 0
    ! SET THE LIMIT NEWTLM ON THE NUMBER OF QUASI-NEWTON STEPS (I.E. CALLS
    ! TO SEARCR), AND IF NEWTLM > 1 SET THE PARAMETER GAIN SUCH THAT NO
    ! FURTHER NEWTON STEPS WILL BE TRIED UNLESS THE LAST STEP REDUCED THE
    ! MAXIMUM STANDARD ERROR BY A FACTOR OF GAIN OR BETTER.
    newtlm = 3
    gain = one/(ten*ten)
    ! FOR NOW, SET JCNTYP(I)=0 IF ICNTYP(I) > 0 AND SET JCNTYP(I)
    ! =ICNTYP(I) OTHERWISE TO DIRECT ERCMP1 TO COMPUTE THE ERRORS FOR THE
    ! STANDARD CONSTRAINTS ONLY.
    do i = 1 , Numgr
        if ( Icntyp(i)<=0 ) then
            Jcntyp(i) = Icntyp(i)
        else
            Jcntyp(i) = 0
        end if
    end do
    ! PUT PARPRJ IN PARWRK FOR USE IN ERCMP1 AND FNSET.
    do j = 1 , Nparm
        Parwrk(j) = Parprj(j)
    end do
    ! CALL ERCMP1 WITH ICNUSE=1 TO COMPUTE THE STANDARD ERRORS.
    ! WE TAKE IPHSE=-3 AS A KLUDGE TO TELL ERCMP1 TO COMPUTE ONLY STANDARD
    ! ERRORS IF THE TEN THOUSANDS DIGIT OF IOPTN IS 1, THUS SAVING ERCMP1
    ! THE WORK OF SCANNING ICNTYP.
    call me%ercmp1(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,Parwrk,1,  &
                   -3,Iwork,Liwrk,Confun,Jcntyp,ipmax,ismax,Err1)

    ! IF THE TYPE -2 AND TYPE -1 ERROR NORMS ARE BOTH <= TOLCON
    ! WE RETURN WITH ICORCT=-1.
    ! NOTE THAT IN THEORY THE TYPE -1 CONSTRAINTS SHOULD BE NO PROBLEM,
    ! BUT OCCASIONALLY THEY ARE VIOLATED DUE TO ROUNDOFF ERROR OR
    ! PROBLEMS IN WOLFE, SO WE CHECK THEM TO BE SAFE.
    if ( Err1(Numgr+3)>Tolcon ) then
        ! HERE THE TYPE -2 ERROR NORM IS > TOLCON AND WE CALL RCHMOD WITH
        ! IRCH=-1 TO SEE IF RCHIN SHOULD BE INCREASED.
        call rchmod(Numgr,Error,Err1,Icntyp,Mact,Iact,ipmax,ismax,Unit,&
                    -1,rchdwn,Rchin)
    else if ( Err1(Numgr+2)<=Tolcon ) then
        Icorct = -1
        return
    end if

    ! PUT PARPRJ INTO THE WORK VECTOR ZWORK SO PARPRJ ITSELF WILL REMAIN
    ! UNCHANGED UNLESS CORRCT IS SUCCESSFUL IN CORRECTING BACK INTO THE
    ! FEASIBLE REGION.
    do j = 1 , Nparm
        Zwork(j) = Parprj(j)
    end do
    ! COMPUTE EOLD = MAX(ERR1(NUMGR+2),ERR1(NUMGR+3)).  NOTE THAT EOLD IS
    ! POSITIVE SINCE OTHERWISE WE WOULD HAVE RETURNED ABOVE (ASSUMING
    ! TOLCON >= 0.0).  THUS IF ONLY ONE TYPE OF STANDARD CONSTRAINT IS
    ! PRESENT, THE FACT THAT ERR1(NUMGR+2) OR ERR1(NUMGR+3) IS ZERO WILL
    ! DO NO HARM.
    eold = Err1(Numgr+3)
    if ( Err1(Numgr+2)>eold ) eold = Err1(Numgr+2)

    ! STATEMENTS ABOVE THIS POINT WILL NOT BE EXECUTED AGAIN IN THIS CALL
    ! TO CORRCT.

    main_loop: do

        ! NOW WE SET UP PMAT FOR USE IN WOLFE TO TRY TO COMPUTE A VECTOR DVEC
        ! POINTING BACK INTO THE FEASIBLE REGION.
        ! IF IOPTTH=1 WE CALL DERST ONCE TO PUT THE STANDARD
        ! GRADIENTS IN CONFUN.
        if ( ioptth>0 ) then
            ! WE SET IPT=-1 TO TELL DERST TO COMPUTE STANDARD CONSTRAINTS ONLY.
            ipt = -1
            call me%derst(Ioptn,Nparm,Numgr,Pttbl,Iptb,Indm,Parwrk,ipt, &
                          Work(ilc24),Work(ilc35),Iwork(ilc22),Confun)
        end if

        l = 0
        do i = 1 , Numgr
            if ( Icntyp(i)+1<0 ) then
                ! HERE ICNTYP(I)=-2 AND WE WILL INCLUDE CONSTRAINT I IF AND ONLY IF
                ! ERR1(I) >= -RCHIN*PROJCT.  WHEN ICNTYP(I)=-1 WE HAVE A LINEAR
                ! STANDARD CONSTRAINT AND IT WILL ALWAYS BE INCLUDED.
                if ( Err1(i)+Rchin*Projct<0 ) cycle
            else if ( Icntyp(i)+1/=0 ) then
                cycle
            end if

            if ( ioptth<=0 ) then
                ! HERE IOPTTH=0 AND WE HAVE NOT YET PLACED THE GRADIENT OF THE LEFT
                ! SIDE OF CONSTRAINT I IN CONFUN(I,.) SO WE DO IT NOW.
                ipt = i
                call me%derst(Ioptn,Nparm,Numgr,Pttbl,Iptb,Indm,Parwrk,ipt, &
                              Work(ilc24),Work(ilc35),Iwork(ilc22),Confun)
            end if

            l = l + 1
            ! PUT THE GRADIENT OF THE LEFT SIDE OF CONSTRAINT I IN PMAT(1,L),...,
            ! PMAT(NPARM,L).
            do k = 1 , Nparm
                Pmat(k,l) = Confun(i,k+1)
            end do

            ! SET ROW NPARM+1 OF PMAT.  WE WILL USUALLY SET PMAT(NPARM+1,L)=
            ! ERR1(I), SO THE WOLFE CONSTRAINT WILL BE GRADIENT(I).DVEC + ERR1(I)
            ! <= 0.0, I.E. (-GRADIENT(I)).DVEC >= ERR1(I).  THE EXCEPTION
            ! OCCURS WHEN ICNTYP(I)=-1 AND ERR1(I) < 0.0, IN WHICH CASE WE
            ! REPLACE ERR1(I) BY ERR1(I)/2.0, IN ORDER TO INSURE THAT EVEN IF PROCOR
            ! TAKES ON ITS MAXIMUM ALLOWED VALUE OF 2.0, NO LINEAR STANDARD
            ! CONSTRAINT WITH NEGATIVE VALUE WILL BECOME POSITIVE VALUED (IGNORING
            ! ROUNDOFF ERROR).  NOTE THAT IF WE DENOTE CONSTRAINT I BY G(I)  <=
            ! 0.0, THEN OUR INEQUALITIES BECOME (GRAD G)(I).DVEC <= -G(I) (OR
            ! -G(I)/2.0), SO A SOLUTION DVEC IS A SOLUTION OF (GRAD G)(I).DVEC =
            ! -G(I) - EPS(I) WHERE EPS(I) = -(GRAD G)(I).DVEC - G(I) = -(GRAD G)(I).
            ! DVEC - G(I)/2.0 - G(I)/2.0 >= 0.0.  NOW WITH H(I) = G(I) + EPS(I)
            ! WE HAVE (GRAD H)(I).DVEC = -H(I), SO IF THIS SYSTEM IS SQUARE THEN
            ! PROCOR=1.0 GIVES A NEWTON STEP FOR SOLVING H(I)=0.0, I.E. G(I) =
            ! -EPS(I) <= 0.0.  THUS WE HAVE A KIND OF GENERALIZED NEWTON METHOD.
            if ( Icntyp(i)+1==0 ) then
                if ( Err1(i)<0 ) then
                    Pmat(npar1,l) = Err1(i)/two
                    cycle
                end if
            end if
            Pmat(npar1,l) = Err1(i)
        end do

        ! CALL WOLFE WITH ISTRT=0 TO COMPUTE DVEC FROM SCRATCH.
        call wolfe(Nparm,l,Pmat,0,s,ncor,Iwork(ilc16),Iwork,Liwrk,Work,   &
                   Lwrk,Work(ilc33),Work(ilc06),Work(ilc31),Work(ilc30),  &
                   Nparm,Numgr,Work(ilc41),Dvec,wdist,nmaj,nmin,jflag)
        if ( jflag<=0 ) then

            ! IN SEARCR AND MULLER WE WILL COMPUTE THE ERROR NORM FOR TYPE -2 AND
            ! TYPE -1 CONSTRAINTS, SO WE LUMP THESE TOGETHER BY SETTING
            ! JCNTYP(I)=-2 IF IT WAS -1.
            do i = 1 , Numgr
                if ( Jcntyp(i)+1==0 ) Jcntyp(i) = -2
            end do
            ! CALL SEARCR TO TRY TO FIND PROCOR SO THAT WITH PARAMETER VECTOR
            ! (OLD) ZWORK + PROCOR*DVEC WE WILL HAVE EMIN = MAX(ERR1(NUMGR+2),
            ! ERR1(NUMGR+3)) <= TOLCON.  IF SEARCR SUCCEEDS IT WILL RETURN WITH
            ! ISRCR=0, WHILE IF IT FAILS IT WILL RETURN WITH ISRCR=1.  IN BOTH
            ! CASES ZWORK WILL BE THE SAME AS BEFORE THE CALL TO SEARCR.
            call me%searcr(Ioptn,Nparm,Numgr,Dvec,Fun,Ifun,Pttbl,Iptb,Indm,   &
                           Zwork,Tolcon,Iphse,Iwork,Liwrk,Work,Lwrk,Parwrk,   &
                           Err1,p1,f1,procor,emin,isrcr)
            if ( isrcr<=0 ) then
                ! HERE THE MAXIMUM STANDARD CONSTRAINT ERROR IS SMALLER
                ! THAN -TOLCON.  SINCE OVERCORRECTION MAY ADVERSELY AFFECT CONVERGENCE,
                ! WE CALL MULLER TO TRY TO GET THE MAXIMUM STANDARD CONSTRAINT
                ! ERROR INTO THE CLOSED INTERVAL (-TOLCON, TOLCON) BY FURTHER
                ! MODIFYING PROCOR.
                if ( emin+Tolcon<0 ) call me%muller(Ioptn,Nparm,Numgr,Dvec,Fun,&
                                                    Ifun,Pttbl,Iptb,Indm,Zwork,Tolcon,Iphse,Iwork,Liwrk, &
                                                    Work,Lwrk,Parwrk,Err1,p1,f1,procor,emin)

                ! NOW COMPUTE PARPRJ = ZWORK + PROCOR*DVEC, SET ICORCT=0, AND RETURN.
                do j = 1 , Nparm
                    Parprj(j) = Zwork(j) + procor*Dvec(j)
                end do
                Icorct = 0
                return
            else
                newtit = newtit + 1
                if ( newtit<newtlm ) then
                    if ( emin<=gain*eold ) then
                        ! HERE WE UPDATE ZWORK, EOLD, PARWRK, AND ERR1, AND TRY ANOTHER NEWTON
                        ! STEP WITH SEARCR.
                        eold = emin
                        do j = 1 , Nparm
                            Zwork(j) = Zwork(j) + procor*Dvec(j)
                            Parwrk(j) = Zwork(j)
                        end do
                        ! WE TAKE IPHSE=-3 AS A KLUDGE TO TELL ERCMP1 TO COMPUTE ONLY STANDARD
                        ! ERRORS IF THE TEN THOUSANDS DIGIT OF IOPTN IS 1, THUS SAVING ERCMP1 THE
                        ! WORK OF SCANNING ICNTYP.
                        call me%ercmp1(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb, &
                                       Indm,Parwrk,1,-3,Iwork,Liwrk,Confun,   &
                                       Jcntyp,ipmax,ismax,Err1)
                        cycle main_loop
                    end if
                end if
            end if
        end if

        exit ! done
    end do main_loop

    ! HERE WE WERE UNABLE TO OBTAIN A FEASIBLE PARPRJ AND WE RETURN WITH
    ! THE WARNING ICORCT=1.
    Icorct = 1

    end subroutine corrct
!********************************************************************************

!********************************************************************************
!>
!  THIS SUBROUTINE USES A MODIFIED QUADRATIC FITTING PROCESS TO SEARCH
!  FOR A PROJECTION FACTOR PROCOR FOR WHICH THE MAXIMUM OF THE LEFT
!  SIDES OF THE TYPE -2 AND -1 CONSTRAINTS EVALUATED AT ZWORK + PROCOR*DVEC
!  IS <= TOLCON.  NOTE THAT WHEN CORRCT CALLS THIS SUBROUTINE IT WILL
!  HAVE LUMPED THE TYPE -1 CONSTRAINTS IN WITH THE TYPE -2 CONSTRAINTS
!  USING JCNTYP, WHICH IS CARRIED THROUGH THIS SUBROUTINE INTO SUBROUTINE
!  ERCMP1 IN IWORK.  IF SEARCR IS ABLE TO FORCE THIS MAXIMUM <= TOLCON
!  IT WILL RETURN WITH ISRCR=0, WITH THE MINIMUM VALUE FOUND FOR THE
!  MAXIMUM IN EMIN, WITH THE CORRESPONDING PROJECTION FACTOR IN PROCOR,
!  WITH THE NUMBER OF TIMES THE MAXIMUM WAS COMPUTED IN NSRCH, AND WITH THE
!  CLOSEST POINT FOUND TO THE LEFT WITH THE MAXIMUM > TOLCON IN (P1,F1).
!  THE SUBROUTINE WILL BEGIN BY COMPUTING THE MAXIMA FOR PROCOR = 1.0,
!  0.5, AND 2.0, AND IF NONE OF THESE MAXIMA IS <= TOLCON AND IT IS
!  NOT THE CASE THAT THE MAXIMUM AT 1.0 IS <= THE OTHER TWO MAXIMA
!  THE SUBROUTINE WILL RETURN WITH THE WARNING ISRCR=1. THE SUBROUTINE
!  WILL ALSO RETURN WITH ISRCR=1 IF IT WOULD NEED TO COMPUTE F MORE THAN
!  LIMSCR TIMES, OR THE SEARCH INTERVAL LENGTH DROPS BELOW TOL1, OR THE
!  QUADRATIC FIT BECOMES TOO FLAT.  EVEN IN THE EVENT OF A RETURN WITH
!  ISRCR=1, EMIN, PROCOR, AND NSRCH WILL BE AS ABOVE.

      subroutine searcr(me,Ioptn,Nparm,Numgr,Dvec,Fun,Ifun,Pttbl,Iptb,Indm,&
                        Zwork,Tolcon,Iphse,Iwork,Liwrk,Work,Lwrk,Parwrk,&
                        Err1,p1,f1,Procor,Emin,Isrcr)

      implicit none

      class(conmax_solver),intent(inout) :: me
      real(wp) :: baladj , balfct , Dvec , Emin , Err1 , f1 , f1kp ,&
                  f2 , f3 , f4 , Fun , fval , p1 , p2 , p3 , &
                  p4 , Parwrk
      real(wp) :: Procor , progr , Pttbl , pval , rlf , rrt , s1 , s2 , &
                  tol1 , tol4 , Tolcon , tolden , Work , &
                  Zwork
      integer :: iaddl , iext , Ifun , ilc08 , ilc21 , ilf , Indm , &
                 Ioptn , Iphse , ipmax , Iptb , irt , ismax , Isrcr , &
                 Iwork , j , limscr , Liwrk , lll
      integer :: Lwrk , Nparm , nsrch , Numgr
!
      dimension Fun(Ifun) , Pttbl(Iptb,Indm) , Parwrk(Nparm) , &
                Err1(Numgr+3) , Zwork(Nparm) , Dvec(Nparm) , &
                Iwork(Liwrk) , Work(Lwrk)

! SET MACHINE AND PRECISION DEPENDENT CONSTANTS.
      tolden = ten*spcmn
      tol1 = ten*ten*spcmn
      tol4 = tol1/four
      balfct = ten
      baladj = (ten-one)/ten
      ilc08 = iloc(8,Nparm,Numgr)
      ilc21 = iloc(21,Nparm,Numgr)
      limscr = 6
      Procor = one
      p1 = zero
      f1 = Err1(Numgr+3)
      f1kp = f1
      Isrcr = 0
      nsrch = 0
      ilf = 0
      irt = 0
! IF AFTER LIMSCR ITERATIONS HAVE BEEN DONE (WHERE LIMSCR >= 4) THE
! BEST VALUE FOUND IS <= PROGR WE WILL (ONCE ONLY) BUMP LIMSCR UP BY
! IADDL, SINCE THERE WOULD SEEM TO BE A GOOD CHANCE THAT THIS WILL
! PRODUCE SUCCESS.
! SETTING IEXT=1 HERE WILL DISABLE THE BUMPING PROCEDURE.
      iext = 0
      iaddl = 6
      progr = ten*ten*ten*Tolcon
! WE NOW TRY TO COMPUTE VALUES AT POINTS P2=PROCOR, P1=P2/2.0, AND
! P3=2.0*P2.
      p2 = Procor
! SET LLL=2 AS THE THREAD THROUGH THE MINOTAURS CAVERN AND JUMP
! DOWN TO PUT F(P2) IN F2.  WE WILL JUMP BACK AFTER ALL SUCH JUMPS
! UNLESS LIMSCR WOULD BE EXCEEDED.
      lll = 2
      pval = p2
      goto 400
 100  Emin = f3
      Procor = two
      goto 800
!
! IF THE SEARCH INTERVAL LENGTH IS LESS THAN TOL1 WE HAVE FAILED.
 200  if ( p3-p1>=tol1 ) then
!
! COMPUTE S1 = THE ABSOLUTE VALUE OF THE SLOPE OF THE LINE THROUGH
! (P1,F1) AND (P2,F2), AND S2 = THE (ABSOLUTE VALUE OF THE) SLOPE
! OF THE LINE THROUGH (P2,F2) AND (P3,F3).
         s1 = (f1-f2)/(p2-p1)
         s2 = (f3-f2)/(p3-p2)
! IF S1+S2 IS VERY SMALL WE HAVE FAILED.
         if ( s1+s2>=tolden ) then
!
            rlf = s2/(s1+s2)
            rrt = one - rlf
! THE MINIMUM OF THE QUADRATIC POLYNOMIAL PASSING THROUGH
! (P1,F1), (P2,F2), AND (P3,F3) WILL OCCUR AT (RLF*P1+
! RRT*P3+P2)/2.0.  NOTE THAT THE THREE POINTS CANNOT BE
! COLLNEAR, ELSE WE WOULD HAVE TERMINATED ABOVE.  SINCE THE
! MINIMUM OCCURS AT THE AVERAGE OF P2 AND A CONVEX COMBINATION
! OF P1 AND P3, IT WILL BE AT LEAST AS CLOSE TO P2 AS TO THE
! ENDPOINT ON THE SAME SIDE.
            if ( ilf>1 ) then
! HERE THE LEFT ENDPOINT WAS DROPPED AT THE LAST ILF > 1
! ITERATIONS, SO TO PREVENT A LONG STRING OF SUCH OCCURRENCES
! WITH LITTLE REDUCTION OF P3-P1 WE WILL SHIFT THE NEW POINT
! TO THE RIGHT BY DECREASING RLF RELATIVE TO RRT.
               rlf = rlf/two**(ilf-1)
               rrt = one - rlf
            else if ( irt>1 ) then
! HERE THE RIGHT ENDPOINT WAS DROPPED AT THE LAST IRT > 1
! ITERATIONS, AND WE WILL SHIFT THE NEW POINT TO THE LEFT.
               rrt = rrt/two**(irt-1)
               rlf = one - rrt
! HERE WE HAVE NOT JUST HAD A STRING OF TWO OR MORE MOVES IN
! THE SAME DIRECTION, BUT IF THE SUBINTERVALS ARE OUT OF
! BALANCE BY MORE THAN A FACTOR OF BALFCT, WE SHIFT THE NEW
! POINT SLIGHTLY IN THE DIRECTION OF THE LONGER INTERVAL.  THE
! IDEA HERE IS THAT THE TWO CLOSE POINTS ARE PROBABLY NEAR THE
! SOLUTION, AND IF WE CAN BRACKET THE SOLUTION WE MAY BE ABLE TO
! CUT OFF THE MAJOR PORTION OF THE LONGER SUBINTERVAL.
            else if ( p2-p1>balfct*(p3-p2) ) then
! HERE THE LEFT SUBINTERVAL IS MORE THAN BALFCT TIMES LONGER THAN
! THE RIGHT SUBINTERVAL, SO WE DECREASE RRT RRELATIVE TO RLF.
               rrt = baladj*rrt
               rlf = one - rrt
            else if ( p3-p2>balfct*(p2-p1) ) then
! HERE THE RIGHT SUBINTERVAL IS MORE THAN BALFCT TIMES LONGER
! THAN THE LEFT SUBINTERVAL, SO WE DECREASE RLF RELATIVE TO RRT.
               rlf = baladj*rlf
               rrt = one - rlf
            end if
!
! COMPUTE THE (POSSIBLY MODIFIED) MINIMUM OF THE QUADRATIC FIT.
            p4 = (rlf*p1+rrt*p3+p2)/two
!
! THE NEXT SECTION (FROM HERE TO STATEMENT 2800) MODIFIES P4, IF
! NECESSARY, TO GET P1+TOL4 <= P2,P4 <= P3-TOL4 AND ABS(P4-P2) >=
! TOL4.  THIS SECTION IS LESS COMPLICATED THAN THE CORRESPONDING SECTION
! IN SEARSL BECAUSE ALL PS LIE BETWEEN 0.5 AND 2.0, SO WEIRD ROUNDOFF
! EFFECTS ARE LESS LIKELY.
! IF ABS(P4-P2) < TOL4 WE REDEFINE P4 BY MOVING TOL4 FROM
! P2 INTO THE LONGER SUBINTERVAL.  NOTE THAT THE LENGTH OF THIS
! SUBINTERVAL MUST BE AT LEAST TOL1/2.0 = 2.0*TOL4, ELSE WE
! WOULD HAVE TERMINATED EARLIER.
            if ( abs(p4-p2)<tol4 ) then
               if ( p3-p2<=(p2-p1) ) then
                  p4 = p2 - tol4
!
! NOW JUMP DOWN TO PUT F(P4) IN F4.
                  pval = p4
               else
                  p4 = p2 + tol4
                  pval = p4
               end if
! HERE WE HAD ABS(P4-P2) >= TOL4 AND WE MAKE SURE THAT P1+TOL4
! <= P4 <= P3-TOL4.
            else if ( p4<=(p3-tol4) ) then
               if ( p4>=(p1+tol4) ) then
                  pval = p4
! HERE P4 < P1+TOL4 AND WE SET P4=P1+TOL4 IF P2-P1 >= TOL1/2.0
! AND OTHERWISE WE SET P4=P2+TOL4.
               else if ( p2-p1<tol1/two ) then
                  p4 = p2 + tol4
                  pval = p4
               else
                  p4 = p1 + tol4
                  pval = p4
               end if
! HERE P4 > P3-TOL4 AND WE SET P4=P3-TOL4 IF P3-P2 >= TOL1/2.0,
! AND OTHERWISE WE SET P4=P2-TOL4.
            else if ( p3-p2<tol1/two ) then
               p4 = p2 - tol4
               pval = p4
            else
               p4 = p3 - tol4
               pval = p4
            end if
            goto 400
         end if
      end if
!
 300  Emin = f2
      Procor = p2
      goto 800
!
! WE INCREMENT NSRCH SINCE WE ARE ABOUT TO COMPUTE F.
 400  nsrch = nsrch + 1
!
!
! THIS IS WHERE WE COMPUTE THE MAXIMUM FVAL = F(PVAL) OF THE LEFT SIDES
! OF THE TYPE -2 AND TYPE -1 CONSTRAINTS.
!
! PROJECT DVEC TO GET PARWRK.
      do j = 1 , Nparm
         Parwrk(j) = Zwork(j) + pval*Dvec(j)
      end do
! WE TAKE IPHSE=-3 AS A KLUDGE TO TELL ERCMP1 TO COMPUTE ONLY STANDARD
! ERRORS IF THE TEN THOUSANDS DIGIT OF IOPTN IS 1, THUS SAVING ERCMP1
! THE WORK OF SCANNING ICNTYP.
      call me%ercmp1(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,Parwrk,1,  &
                  -3,Iwork,Liwrk,Work(ilc08),Iwork(ilc21),ipmax,ismax,  &
                  Err1)
      fval = Err1(Numgr+3)
!
      if ( fval<=Tolcon ) then
!
! HERE FVAL <= TOLCON AND WE RETURN AFTER SETTING PROCOR, EMIN, P1,
! AND F1.
         Procor = pval
         Emin = fval
! IF LLL=1 TAKE P1=0.0 AND F1=F1KP, IF LLL=2 LEAVE P1 AND F1 ALONE (THEY
! WILL BE 0.0 AND FIKP RESPECTIVELY), IF LLL=3 TAKE P1=P2 AND F1=F2,
! IF LLL=4 AND P2 < P4 TAKE P1=P2 AND F1=F2, AND IF LLL=4 AND
! P2 >= P4 LEAVE P1 AND F1 ALONE.  IN ALL CASES (P1,F1) WILL BE THE
! POINT WITH P1 THE NEAREST VALUE LEFT OF PROCOR CONSIDERED AND WE WILL
! HAVE F1 > TOLCON.
         select case (lll)
         case (2)
         case (3)
            goto 500
         case (4)
            if ( p2<p4 ) goto 500
         case default
            p1 = zero
            f1 = f1kp
         end select
         return
      else
!
! HERE FVAL > TOLCON AND WE SEE IF LIMSCR ITERATIONS IN SEARCR HAVE
! BEEN DONE.  IF SO WE SET THE FAILURE WARNING ISRCR=1 AND RETURN
! UNLESS WE CHOOSE TO INCREASE LIMSCR.
         if ( nsrch<limscr ) goto 900
!
! HERE WE HAVE DONE LIMSCR ITERATIONS.
         if ( iext>0 ) goto 700
         if ( fval<=progr ) goto 600
         if ( f2>progr ) goto 700
         goto 600
      end if
 500  p1 = p2
      f1 = f2
      return
! HERE WE HAVE NOT BUMPED LIMSCR EARLIER, LIMSCR >= 4, AND
! MIN(FVAL,F2) <= PROGR, SO WE BUMP LIMSCR.
 600  iext = 1
      limscr = limscr + iaddl
      goto 900
!
!11000 WRITE(NWRIT,11500)
!11500 FORMAT(30H *****WARNING*****WARNING*****,
!    *30H TOO MANY ITERATIONS IN SEARCR)
!
! HERE WE HAVE FAILED AND WE SET EMIN AND PROCOR FOR OUTPUT, SET ISRCR=1,
! AND RETURN.
 700  if ( fval>f2 ) goto 300
      Emin = fval
      Procor = pval
!
 800  Isrcr = 1
      return
!
! HERE WE WILL CARRY THE COMPUTED F VALUE BACK TO THE APPROPRIATE PART
! OF THE PROGRAM.
 900  select case (lll)
      case (1)
!
         f1 = fval
         p3 = two*p2
! HERE SET LLL=3 AND PUT F(P3) IN F3.
         lll = 3
         pval = p3
         goto 400
      case (2)
!
         f2 = fval
         p1 = p2/two
! SET LLL=1 AND PUT F(P1) IN F1.
         lll = 1
         pval = p1
         goto 400
      case (3)
!
         f3 = fval
!
! WE NOW HAVE FOUND P1, P2, AND P3 WITH CORRESPONDING VALUES
! F1, F2, AND F3, AND WE CHECK WHETHER F2 <= MIN(F1,F3).
         if ( f2<=f1 ) then
!
! HERE F2 <= F1.  IF F2 <= F3 WE ARE DONE INITIALIZING.
            if ( f2>f3 ) goto 100
! END OF INITIALIZATION.
!
! ASSUMING THAT P3-P1 >= TOL1, WE NOW HAVE POINTS P1, P2, P3 WITH
! P1 <= P2-TOL4, P2 <= P3-TOL4, F1=F(P1) >= F2=F(P2), AND F3=F(P3)
! >= F2.  THESE CONDITIONS WILL BE MAINTAINED THROUGHOUT THE PROGRAM.
! SET LLL=4, WHERE IT WILL REMAIN FROM NOW ON.
            lll = 4
            goto 200
         else
            if ( f1>f3 ) goto 100
            Emin = f1
            Procor = one/two
            goto 800
         end if
      case (4)
!
         f4 = fval
!
! NOW WE DROP EITHER P1 OR P3 AND RELABEL THE REMAINING POINTS (IF
! NECESSARY) SO THAT F(P2) <= F(P1) AND F(P2) <= F(P3).
! IF NOW THE LEFTMOST OF THE TWO MIDDLE POINTS IS LOWER THAN THE
! RIGHTMOST OF THE TWO MIDDLE POINTS WE DROP P3, AND SET ILF=0
! AND INCREMENT IRT TO INDICATE THE RIGHT END POINT HAS BEEN DROPPED.
! OTHERWISE WE DROP P1, SET IRT=0 AND INCREMENT ILF.  IN ALL CASES
! WE THEN RESHUFFLE THE VALUES INTO P1, P2, P3, F1, F2, F3 AND TRY
! TO DO ANOTHER ITERATION.
         if ( p4<p2 ) then
! HERE P4 < P2.
            if ( f4<f2 ) then
               p3 = p2
               f3 = f2
               p2 = p4
               f2 = f4
               ilf = 0
               irt = irt + 1
            else
               p1 = p4
               f1 = f4
               ilf = ilf + 1
               irt = 0
            end if
! HERE P4 > P2.
         else if ( f2<f4 ) then
            p3 = p4
            f3 = f4
            ilf = 0
            irt = irt + 1
         else
            p1 = p2
            f1 = f2
            p2 = p4
            f2 = f4
            ilf = ilf + 1
            irt = 0
         end if
         goto 200
      case default
      end select
    end subroutine searcr
!********************************************************************************

!********************************************************************************
!>
!  In this subroutine we are given a base vector zwork, a direction
!  vector dvec, a scalar procor with emin = f(procor) = (the maximum type
!  -2 and -1 error with parameters zwork + procor*dvec) < -tolcon, and
!  a scalar p1 with p1 < procor and f1 = f(p1) > tolcon.  we do
!  a revised mullers method approach (with a solution contained in a
!  shrinking interval) to attempt to adjust procor so that -tolcon  <=
!  f(procor) <= tolcon, but if we are not successful we return with the
!  leftmost procor found satisfying emin = f(procor) < -tolcon on the
!  theory that overcorrection is better than no correction.  note that when
!  corrct calls this subroutine it will have lumped the type -1 constraints
!  in with the type -2 constraints using jcntyp, which is carried through
!  this subroutine into subroutine ercmp1 in iwork.

    subroutine muller(me,Ioptn,Nparm,Numgr,Dvec,Fun,Ifun,Pttbl,Iptb,Indm,&
                      Zwork,Tolcon,Iphse,Iwork,Liwrk,Work,Lwrk,Parwrk,&
                      Err1,p1,f1,Procor,Emin)

    implicit none

    class(conmax_solver),intent(inout) :: me
    integer,intent(in) :: Ifun
    integer,intent(in) :: Indm
    integer,intent(in) :: Iptb
    integer,intent(in) :: Lwrk
    integer,intent(in) :: Nparm
    integer,intent(in) :: Numgr
    integer :: Ioptn
    integer :: Iphse
    integer :: Liwrk
    real(wp) :: Emin
    real(wp) :: f1
    real(wp) :: p1
    real(wp) :: Procor
    real(wp) :: Tolcon
    integer :: Iwork(Liwrk)
    real(wp) :: Dvec(Nparm)
    real(wp) :: Err1(Numgr+3)
    real(wp) :: Fun(Ifun)
    real(wp) :: Parwrk(Nparm)
    real(wp) :: Pttbl(Iptb,Indm)
    real(wp) :: Work(Lwrk)
    real(wp) :: Zwork(Nparm)

    real(wp) :: acof , bcof , ccof , den , discr , f2 , f3 , &
                f4 , fval , p2 , p3 , p4 , pval , temp
    integer ::  ilc08 , ilc21 , imain , ipmax , ismax , j , limmul , &
                lll , nsrch

    real(wp),parameter :: tol1 = ten*ten*spcmn
    real(wp),parameter :: tol4 = tol1/four
    real(wp),parameter :: tolden = ten*spcmn

    ! SET MACHINE AND PRECISION DEPENDENT CONSTANTS.
    ilc08 = iloc(8,Nparm,Numgr)
    ilc21 = iloc(21,Nparm,Numgr)
    limmul = 5
    nsrch = 0
    imain = 0
    p3 = Procor
    f3 = Emin

    ! WE DO NOT ALLOW THE LENGTH OF THE INTERVAL (P1,P3) TO FALL BELOW
    ! TOL1.
100 if ( p3-p1<tol1 ) then
        return
    else
        ! COMPUTE P2 = (P1+P3)/2.0 AND F(P2).
        p2 = (p1+p3)/two
        pval = p2
        ! SET LLL AS THE THREAD THROUGH THE MINOTAURS CAVERN AND JUMP DOWN TO
        ! COMPUTE F(PVAL)=F(P2).  WE WILL JUMP BACK AFTER ALL SUCH JUMPS.
        lll = 1
        goto 1000
    end if

    ! HERE -TOLCON <= F2 <= TOLCON AND WE RETURN WITH PROCOR=P2 AND
    ! EMIN=F2.
200 Procor = p2
    Emin = f2
    return

    ! HERE WE HAVE NOT ACHIEVED SUCCESS YET AND WE SEE IF THE ITERATION
    ! LIMIT HAS BEEN REACHED.
300 if ( nsrch<limmul ) then

        ! HERE WE HAVE NOT REACHED THE ITERATION LIMIT SO WE TRY AGAIN.
        ! IF IMAIN=0 HERE WE WILL HAVE NO P4 TO SHUFFLE IN, AND WE WILL HAVE
        ! ALREADY CHECKED P3-P1 >= TOL1, SO WE RESET IMAIN TO 1 AND DO A FIT.
        if ( imain<=0 ) goto 900

        ! HERE WE HAVE POINTS P1, P2, P3, P4 WITH P1+TOL1/4.0 <= P2 <=
        ! P3-TOL1/4.0, P1+TOL1/4.0 <= P4 <= P3-TOL1/4.0, ABS(P4-P2) >=
        ! TOL1/4.0, F(P1) > TOLCON, F(P3) < -TOLCON, ABS(F(P2)) >
        ! TOLCON, AND ABS(F(P4)) > TOLCON.  WE WILL NOW DISCARD EITHER
        ! P1 OR P3 AND RELABEL TO GET NEW POINTS P1, P2, P3, EXCEPT IN ONE
        ! CASE WHERE TWO POINTS WILL BE DISCARDED AND WE WILL RELABEL TO GET
        ! NEW POINTS P1, P3.
        ! IF P2 > P4 HERE WE WILL, IN THE INTEREST OF A MORE READABLE
        ! PROGRAM, INTERCHANGE P2 AND P4 (AND F2 AND F4) SO WE WILL BE ABLE
        ! TO ASSUME P2 <= P4.
        if ( p2>p4 ) then
            temp = p2
            p2 = p4
            p4 = temp
            temp = f2
            f2 = f4
            f4 = temp
        end if
        if ( f2<=0 ) then
            ! HERE F2 < 0.0.
            if ( f4<=0 ) goto 700
            ! HERE F2 < 0.0 AND F4 > 0.0, AND IN THIS SAWTOOTH PATTERN WE
            ! DISCARD BOTH P4 AND P3, SET IMAIN=0, AND GO BACK TO THE BEGINNING
            ! (EXCEPT NSRCH CONTINUES TO INCREASE, INSURING EVENTUAL TERMINATION).
            imain = 0
            p3 = p2
            f3 = f2
            Procor = p3
            Emin = f3
            goto 100
        else
            ! HERE F2 > 0.0.
            if ( f4<=0 ) then
                ! HERE F2 > 0.0 AND F4 < 0.0.
                if ( p2-p1<=(p3-p4) ) goto 700
            end if
            ! HERE EITHER F2 > 0.0 AND F4 > 0.0, OR ELSE F2 > 0.0,
            ! F4 < 0.0, AND P2-P1 > P3-P4.  WE DISCARD P1, SINCE IN THE
            ! FORMER CASE THE FIRST THREE F VALUES ARE ALL POSITIVE, AND IN THE
            ! LATTER CASE ONLY THE FIRST TWO F VALUES ARE POSITIVE, BUT BY DROPPING
            ! P1 WE CAN GET MAXIMUM SHRINKAGE OF P3-P1.
            p1 = p2
            f1 = f2
            p2 = p4
            f2 = f4
            goto 800
        end if
    else

        ! HERE WE HAVE REACHED THE ITERATION LIMIT WITHOUT SUCCESS.  WE RETURN
        ! WITH PROCOR = THE LEFTMOST OF THE THREE POINTS P2, P4, AND P3 WHICH
        ! HAS NEGATIVE F VALUE (UNLESS IMAIN=0, IN WHICH CASE WE IGNORE P4).
        ! 600 WRITE(NWRIT,700)
        ! 700 FORMAT(45H ***WARNING  TOO MANY ITERATIONS IN MULLER***)
        if ( imain<=0 ) goto 600
        if ( p2<=p4 ) then
            ! HERE P2 < P4.
            if ( f2<0 ) goto 200
            if ( f4>=0 ) goto 500
            ! HERE P4 < P2.
        else if ( f4>=0 ) then
            goto 600
        end if
    end if
400 Procor = p4
    Emin = f4
    return

500 Procor = p3
    Emin = f3
    return

600 if ( f2>=0 ) goto 500
    goto 200

    ! HERE EITHER F2 < 0.0 AND F4 < 0.0, OR ELSE F2 > 0.0,
    ! F4 < 0.0, AND P2-P1 <= P3-P4.  WE DISCARD P3, SINCE IN THE
    ! FORMER CASE THE LAST THREE F VALUES ARE NEGATIVE, AND IN THE LATTER
    ! CASE ONLY THE LAST TWO F VALUES ARE NEGATIVE, BUT BY DROPPING P3 WE
    ! GET MAXIMUM SHRINKAGE OF P3-P1.
    700  p3 = p4
    f3 = f4

    ! HERE WE HAVE THREE POINTS.  IF P3-P1 < TOL1 WE WILL RETURN AFTER
    ! SETTING PROCOR AND EMIN.
800 if ( p3-p1<tol1 ) goto 600

    ! HERE WE RESET IMAIN TO 1 AND COMPUTE P4, THE UNIQUE ZERO IN THE
    ! INTERVAL (P1,P3) OF THE QUADRATIC POLYNOMIAL WHICH PASSES THROUGH
    ! (P1,F1), (P2,F2), AND (P3,F3).  RECALL THAT F1 > 0.0,
    ! F3 < 0.0, AND P1+TOL1/4.0 <= P2 <= P3-TOL1/4.0.
900 imain = 1

    ! COMPUTE THE COEFFICIENTS ACOF, BCOF, AND CCOF OF OUR POLYNOMIAL
    ! ACOF*X**2 + BCOF*X + CCOF.
    acof = ((f3-f2)*(p2-p1)-(f2-f1)*(p3-p2))/((p2-p1)*(p3-p2)*(p3-p1))
    bcof = (f3-f1)/(p3-p1) - acof*(p1+p3)
    ccof = f2 - p2*(acof*p2+bcof)
    discr = bcof**2 - four*acof*ccof
    ! IN THEORY THE DISCRIMINANT SHOULD BE POSITIVE HERE, BUT TO BE SAFE WE
    ! CHECK IT IN CASE ROUNDOFF ERROR HAS MADE IT NEGATIVE.
    if ( discr<0 ) goto 600
    if ( bcof<=0 ) then
        ! HERE BCOF <= 0.0 AND WE USE THE ALTERNATE FORM OF THE QUADRATIC
        ! FORMULA.  NOTE THAT THE DENOMINATOR CANNOT BE ZERO SINCE THAT
        ! WOULD IMPLY BOTH BCOF=0.0 AND SQRT(...)=0.0, SO ALSO EITHER
        ! ACOF=0.0 OR CCOF=0.0, BUT THIS CONTRADICTS THE FACT THAT F1  >
        ! 0.0 AND F3 < 0.0.
        ! STILL, TO BE SAFE, WE CHECK THE SIZE OF THE DENOMINATOR.
        den = -bcof + sqrt(discr)
        if ( den<tolden ) goto 600
        p4 = two*ccof/den
    else
        ! HERE BCOF > 0.0 AND WE USE THE USUAL FORM OF THE QUADRATIC
        ! FORMULA TO TRY TO REDUCE PROBLEMS WITH SUBTRACTION AND SMALL
        ! DENOMINATORS.  THE MINUS SIGN IS USED IN FRONT OF THE SQUARE ROOT
        ! BECAUSE IF ACOF > 0.0 THEN THE POLYNOMIAL IS CONCAVE UP, WHICH
        ! IMPLIES P1 MUST BE ON THE LEFT BRANCH (SINCE F1 > F3), WHICH
        ! IMPLIES WE WANT THE LEFT (I.E. SMALLER) ZERO, AGREEING WITH
        ! -SQRT(...)/ACOF <= 0.0.  IF ON THE OTHER HAND ACOF < 0.0 THEN
        ! THE POLYNOMIAL IS CONCAVE DOWN, WHICH IMPLIES P3 MUST BE ON THE
        ! RIGHT BRANCH (SINCE F1 > F3), WHICH IMPLIES WE WANT THE RIGHT
        ! (I.E. LARGER) ZERO, AGREEING WITH -SQRT(...)/ACOF >= 0.0.
        ! NOTE THAT ACOF=0.0 CANNOT OCCUR HERE SINCE IF IT DID THE POLYNOMIAL
        ! WOULD BE LINEAR, AND BCOF > 0.0 WOULD THEN CONTRADICT F1 > F3.
        ! STILL, TO BE SAFE, WE CHECK THE SIZE OF THE DENOMINATOR.
        den = two*acof
        if ( abs(den)<tolden ) goto 600
        p4 = (-bcof-sqrt(discr))/den
    end if

    ! THE NEXT SECTION (FROM HERE TO STATEMENT 3200) MODIFIES P4, IF
    ! NECESSARY, TO GET P1+TOL4 <= P2,P4 <= P3-TOL4 AND ABS(P4-P2) >=
    ! TOL4.

    ! IF ABS(P4-P2) < TOL1/4.0 WE REDEFINE P4 BY MOVING IT A DISTANCE
    ! TOL1/4.0 FROM P2 INTO THE LONGER SUBINTERVAL.  NOTE THAT THE LENGTH
    ! OF THIS SUBINTERVAL MUST BE AT LEAST TOL1/2.0 SINCE P3-P1 >= TOL1.
    if ( abs(p4-p2)<tol4 ) then
        if ( p3-p2<=(p2-p1) ) then
            p4 = p2 - tol4
        else
            p4 = p2 + tol4
        end if
        ! HERE WE HAD ABS(P4-P2) >= TOL4 AND WE MAKE SURE THAT P1+TOL4
        ! <= P4 <= P3-TOL4.
    else if ( p4<=(p3-tol4) ) then
        if ( p4<(p1+tol4) ) then
            ! HERE P4 < P1+TOL4 AND WE SET P4=P1+TOL4 IF P2-P1 >= TOL1/2.0
            ! AND OTHERWISE WE SET P4=P2+TOL4.
            if ( p2-p1<tol1/two ) then
                p4 = p2 + tol4
            else
                p4 = p1 + tol4
            end if
        end if
        ! HERE P4 > P3-TOL4 AND WE SET P4=P3-TOL4 IF P3-P2 >= TOL1/2.0,
        ! AND OTHERWISE WE SET P4=P2-TOL4.
    else if ( p3-p2<tol1/two ) then
        p4 = p2 - tol4
    else
        p4 = p3 - tol4
    end if

    ! COMPUTE F4=F(P4).
    pval = p4
    lll = 2

    ! NOW INCREMENT NSRCH SINCE WE ARE ABOUT TO COMPUTE F.
1000 nsrch = nsrch + 1

    ! HERE IS WHERE WE MUST SUPPLY A ROUTINE TO COMPUTE FVAL = F(PVAL) =
    ! THE MAXIMUM OF THE LEFT SIDES OF THE TYPE -2 AND -1 CONSTRAINTS.

    ! PROJECT DVEC TO GET PARWRK FOR USE IN ERCMP1.
    do j = 1 , Nparm
        Parwrk(j) = Zwork(j) + pval*Dvec(j)
    end do
    ! WE TAKE IPHSE=-3 AS A KLUDGE TO TELL ERCMP1 TO COMPUTE ONLY STANDARD
    ! ERRORS IF THE TEN THOUSANDS DIGIT OF IOPTN IS 1, THUS SAVING ERCMP1
    ! THE WORK OF SCANNING ICNTYP.
    call me%ercmp1(Ioptn,Nparm,Numgr,Fun,Ifun,Pttbl,Iptb,Indm,Parwrk,1,  &
                   -3,Iwork,Liwrk,Work(ilc08),Iwork(ilc21),ipmax,ismax,  &
                   Err1)
    fval = Err1(Numgr+3)

    ! CARRY THE COMPUTED F VALUE BACK TO THE APPROPRIATE PART OF THE PROGRAM.
    if ( lll==1 ) then
        f2 = fval
        if ( f2>Tolcon ) goto 300
        if ( f2+Tolcon>=0 ) goto 200
        goto 300
    else if ( lll==2 ) then
        f4 = fval
        ! IF -TOLCON <= F4 <= TOLCON WE RETURN WITH PROCOR=P4 AND EMIN
        ! =F4, AND OTHERWISE WE GO BACK UP TO SEE IF WE HAVE REACHED THE LIMIT
        ! ON THE NUMBER OF STEPS.
        if ( f4>Tolcon ) goto 300
        if ( f4+Tolcon>=0 ) goto 400
        goto 300
    end if

    end subroutine muller
!********************************************************************************

!********************************************************************************
!>
!  This subroutine increases rchdwn or rchin if it appears some
!  constraints which should have been declared active were not so
!  declared.

    subroutine rchmod(Numgr,Error,Err1,Icntyp,Mact,Iact,Ipmax,Ismax, &
                      Unit,Irch,Rchdwn,Rchin)

    implicit none

    integer,intent(in)     :: Numgr
    integer,intent(in)     :: Iact(Numgr)
    integer,intent(in)     :: Icntyp(Numgr)
    integer,intent(in)     :: Ipmax
    integer,intent(in)     :: Irch
    integer,intent(in)     :: Ismax
    integer,intent(in)     :: Mact
    real(wp),intent(in)    :: Err1(Numgr+3)
    real(wp),intent(in)    :: Error(Numgr+3)
    real(wp),intent(inout) :: Rchdwn
    real(wp),intent(inout) :: Rchin
    real(wp),intent(in)    :: Unit

    real(wp) :: ei , eipmax , enorm , epact , rch1 , rchd1
    integer  :: i , ipact , l

    real(wp),parameter :: fudge = one + one/ten
    real(wp),parameter :: rchtop = one/spcmn

    ! SET MACHINE AND PRECISION DEPENDENT CONSTANTS.
    enorm = Error(Numgr+1)

    if ( Irch<0 ) then

        ! HERE IRCH=-1 AND WE CONSIDER CHANGING RCHIN.

        ! SEE IF CONSTRAINT ISMAX IS IN THE ACTIVE SET, AND RETURN IF IT IS.
        ! NOTE THAT ISMAX > 0 SINCE WE WOULD NOT HAVE CALLED RCHMOD WITH
        ! IRCH=-1 IF THERE WERE NO STANDARD CONSTRAINTS.
        do l = 1 , Mact
            i = abs(Iact(l))
            if ( i==Ismax ) return
        end do

        ! RETURN IF RCHIN >= RCHTOP.
        if ( Rchin<rchtop ) then

            ! SET THE PROSPECTIVE NEW RCHIN.  NOTE THAT WITHOUT THE FUDGE FACTOR,
            ! RCH1 WOULD HAVE BEEN JUST BARELY LARGE ENOUGH TO HAVE CAUSED
            ! CONSTRAINT ISMAX TO BE DECLARED ACTIVE WHEN THE OLD ACTIVE SET WAS
            ! DETERMINED.  (NOTE THAT RCHIN MAY HAVE ALREADY BEEN INCREASED SINCE
            ! THEN.  NOTE ALSO THAT ERROR(ISMAX) < 0.0, ELSE CONSTRAINT ISMAX
            ! WOULD HAVE BEEN DECLARED ACTIVE.)
            rch1 = fudge*(-Error(Ismax))/Unit

            ! IF RCH1 > RCHIN WE REPLACE RCHIN BY MIN(RICH1,RCHTOP).
            if ( rch1>Rchin ) then
                Rchin = rch1
                if ( Rchin>rchtop ) Rchin = rchtop
            end if

        end if

        return

    else

        ! HERE IRCH=1 AND WE CONSIDER CHANGING RCHDWN.

        ! SEE IF CONSTRAINT IPMAX IS IN THE ACTIVE SET, AND RETURN IF IT IS.
        ! NOTE THAT IPMAX > 0 SINCE THERE WILL BE AT LEAST ONE PRIMARY
        ! CONSTRAINT AT THIS STAGE (EVEN IF THERE WERE NONE IN THE ORIGINAL
        ! PROBLEM).
        do l = 1 , Mact
            i = abs(Iact(l))
            if ( i==Ipmax ) return
        end do

        ! RETURN IF RCHDWN >= RCHTOP.
        if ( Rchdwn>=rchtop ) return

        ! WE WILL CONSIDER CHANGING RCHDWN IF THE NEW PRIMARY ERROR NORM WITH
        ! ONLY THE OLD ACTIVE CONSTRAINTS CONSIDERED IS LESS THAN THE OLD
        ! PRIMARY ERROR NORM, AND THIS WILL CERTAINLY BE THE CASE IF THE NEW
        ! PRIMARY ERROR NORM IS LESS THAN THE OLD PRIMARY ERROR NORM.
        if ( Err1(Numgr+1)>=enorm ) then
            ! COMPUTE EPACT, THE NEW PRIMARY ERROR NORM WITH ONLY THE OLD ACTIVE
            ! CONSTRAINTS CONSIDERED.
            ipact = 0
            do l = 1 , Mact
                i = abs(Iact(l))
                if ( Icntyp(i)<1 ) cycle
                if ( Icntyp(i)==1 ) then
                    ! HERE CONSTRAINT I WAS A PRIMARY ACTIVE CONSTRAINT.
                    ei = Err1(i)
                else
                    ei = abs(Err1(i))
                end if
                if ( ipact<=0 ) then
                    ipact = 1
                    epact = ei
                else if ( ei>epact ) then
                    epact = ei
                end if
            end do
            ! WE WILL RETURN IF EPACT IS >= THE OLD PRIMARY ERROR NORM, WHICH
            ! WOULD INDICATE THAT THE STEP WAS TOO INACCURATE TO BE TRUSTED TO
            ! USE IN MODIFYING RCHDWN.
        end if

        ! COMPUTE EIPMAX AS THE OLD ERROR AT CONSTRAINT IPMAX (IF ICNTYP(IPMAX)
        ! =1) OR THE OLD ABSOLUTE ERROR AT CONSTRAINT IPMAX (IF ICNTYP(IPMAX)
        ! =2).  NOTE THAT HERE ICNTYP(IPMAX) MUST BE 1 OR 2 SINCE ERCMP1
        ! COMPUTED IPMAX AS THE INDEX OF THE PRIMARY CONSTRAINT WHERE THE
        ! MAXIMUM PRIMARY CONSTRAINT ERROR (I.E. VALUE) WAS ACHIEVED.
        if ( Icntyp(Ipmax)<=1 ) then
            eipmax = Error(Ipmax)
        else
            eipmax = abs(Error(Ipmax))
        end if

        ! SET THE PROSPECTIVE NEW RCHDWN.  NOTE THAT WITHOUT THE FUDGE FACTOR,
        ! RCHD1 WOULD HAVE JUST BARELY BEEN LARGE ENOUGH TO HAVE CAUSED
        ! CONSTRAINT IPMAX TO BE DECLARED ACTIVE WHEN THE OLD ACTIVE SET WAS
        ! DETERMINED.  (NOTE THAT RCHDWN MAY HAVE ALREADY BEEN INCREASED
        ! SINCE THEN.)
        rchd1 = fudge*(enorm-eipmax)/Unit

        ! IF RCHD1 > RCHDWN WE REPLACE RCHDWN BY MIN (RCHD1, RCHTOP).
        if ( rchd1<=Rchdwn ) return
        Rchdwn = rchd1
        if ( Rchdwn>rchtop ) Rchdwn = rchtop
    end if

    ! WRITE(NWRIT,'(A,E24.14)') '***RCHDWN INCREASED TO', RCHDWN
    ! WRITE(NWRIT,'(A,E24.14)') '***RCHIN INCREASED TO ', RCHIN

    end subroutine rchmod
!********************************************************************************

!********************************************************************************
!>
!  given m inequalities of the form a(k).x + b(k) <= 0.0 for k=1,
!  ...,m, where a(k) and x are ndm dimensional vectors and b(k)
!  are numbers, this subroutine returns the nearest point to the
!  origin in the polytope defined by these inequalities (unless
!  jflag > 0, which indicates failure).  the user should put
!  the mdm+1 dimensional vectors (a(k),b(k)) in the columns of pmat.
!  the solution point will be returned in wpt, and will also be a
!  linear combination of the a(k) vectors with (nonpositive)
!  coefficients in the m dimensional vector wcoef.  wcoef may not be
!  accurate if refwl was used to refine wpt, which rarely happens. the
!  number of vectors in the final corral will be returned in ncor with
!  their indices in icor, and all entries of wcoef not corresponding to
!  indices in icor will be zero.  the distance will be returned in
!  wdist, and the numbers of major and minor cycles in the cone
!  subproblem will be returned in nmaj and nmin respectively.
!  if the user sets istrt=0 the program will start from scratch, but
!  the user can set istrt=1 (hot start) and specify ncor, icor, wcoef,
!  and the factor s.  (see later comments; set s=1.0 if no better value is
!  available.  set wcoef(j)=0.0 if icor(i) /= j for i=1,...,ncor.)  (if
!  inaccurate wcoef or s is used in a hot start attempt little will be
!  lost, since ncor and icor are more important for a successful hot start
!  than wcoef and s.)  we must always have ncor <= ndm+1 in theory since
!  the ncor ndm+1 dimensional vectors in a corral should be linearly
!  independent, and in practice we will always require ncor <= ndm+1.
!  if the user sets istrt=1 but the program fails, it will
!  automatically try from scratch before giving up.
!
!### Running by itself rather than as a part of conmax
!  To run the program, write a driver
!  program which dimensions the arrays in the calling sequence for wolfe
!  and sets the input variables as specified in the list below, then call
!  subroutine wolfe.
!
!### Reference
!  This program was developed by ed kaufman, david leeming, and jerry
!  taylor.  the method used is an enhanced version of the method described
!  in (wolfe, philip, finding the nearest point in a polytope, mathematical
!  programming 11 (1976), 128-149).

    subroutine wolfe(Ndm,m,Pmat,Istrt,s,Ncor,Icor,Iwork,Liwrk,Work, &
                     Lwrk,r,Coef,Ptnr,Pmat1,Nparm,Numgr,Wcoef,Wpt, &
                     Wdist,Nmaj,Nmin,Jflag)

    implicit none

    integer,intent(in) :: Ndm  !! The number of variables. it must be less than or equal to `nparm`.
    integer,intent(in) :: m  !! The number of inequalities defining the polytope.  it
                            !! must be less than or equal to `numgr`.
    integer,intent(in) :: Istrt   !! Set this equal to zero unless a hot start is desired--
                                !! see next paragraph of comments for more details.  if istrt is set
                                !! equal to 1, then s, wcoef, ncor, and icor must also be assigned
                                !! values initially.
    integer,intent(out) :: Ncor  !! This is the number of vectors (i.e. colunns of pmat) in
                                !! the final corral.
    integer,intent(in) :: Liwrk   !! This is the dimension of `iwork`.  It must be at least 7*nparm + 7*numgr + 3.
    integer,intent(in) :: Lwrk    !! This is the dimension of `work`.  it must be at least
                                !! 2*nparm**2 + 4*numgr*nparm + 11*numgr + 27*nparm + 13.
                                !! note that some storage could be saved by rewriting function
                                !! subprogram iloc to take out all but the arrays needed (namely 1, 3,
                                !! 4, 9, 28, 32, 34, 39 for work, 18, 23 for iwork) and scrunching
                                !! work and iwork in iloc so the remaining arrays follow one after
                                !! another.
    integer,intent(in)  :: Nparm  !! This is basically a dimension parameter here.  it must
                                !! be greater than or equal to `ndm`.
    integer,intent(in)  :: Numgr  !! This is basically a dimension parameter here.  it must
                                !! be greater than or equal to `m`.
    integer,intent(out)  :: Nmaj !! This will be the number of major cycles used in `wolfe`.
    integer,intent(out)  :: Nmin !! This will be the number of major cycles used in `wolfe`.
    integer,intent(out)  :: Jflag !! This is a flag variable which is 0 in case of a normal
                                !! solution and is positive otherwise (in which case the returned
                                !! solution may be no good).
    real(wp),intent(out) :: s !! You may ignore this scale factor unless you want to use the hot start option.
    real(wp),intent(out) :: Wdist !! This will be the (minimized) euclidean distance of `wpt` from the origin.
    integer,intent(out) :: Icor(Nparm+1)  !! This array contains the ncor indices of the vectors in
                                        !! the final corral.  its dimension in the driver program must be at
                                        !! least nparm+1.
    real(wp),intent(in) :: Pmat(Nparm+1,Numgr)  !!  This is an array whose kth column contains the vector
                                                !!  (a(k),b(k)) for k=1,...,m, where the m inequalities a(k).x + b(k)
                                                !!  <= 0.0 define the polytope whose nearest point to the origin we
                                                !!  seek.  the first dimension of pmat in the driver program must be
                                                !!  exactly nparm+1, while the second dimension of pmat in the driver
                                                !!  program must be at least numgr.
                                                !!  if we actually want the nearest point in the polytope to some point
                                                !!  y other than the origin, we translate y to the origin before calling
                                                !!  wolfe, that is, call wolfe to find the nearest point z to the origin
                                                !!  in the polytope defined by a(k).z + (b(k) + a(k).y) <= 0.0, then
                                                !!  compute x = y + z.
    integer  :: Iwork(Liwrk) !! Work array
    real(wp) :: Work(Lwrk) !! Work array
    real(wp) :: r(Nparm+1) !! Work array
    real(wp) :: Coef(Numgr) !! Work array
    real(wp) :: Ptnr(Nparm+1) !! Work array
    real(wp) :: Pmat1(Nparm+1,Numgr) !! Work array
    real(wp),intent(out) :: Wcoef(Numgr)  !! This will give the coefficients of the vectors a(k)
                                        !! needed to form a linear combination equal to the solution in wpt.
                                        !! its dimension in the driver program must be at least numgr.
                                        !! wcoef may not be accurate if it was necessary to call refwl to
                                        !! refine wpt, which rarely happens.
    real(wp),intent(out) :: Wpt(Nparm)    !! This will give the coordinates of the point we are seeking,
                                        !! namely the nearest point in the polytope to the origin.  its dimension
                                        !! in the driver program must be at least nparm.

    real(wp) :: ab , bk , dist , fackp , facsc , &
                fact , quot , s1 , s1hi , s1low , &
                scl , scl1 , scl1a , tol , tol1 , &
                tols , v1 , violm , vmax
    integer :: i , ilc18 , ilc28 , ilc32 , ilc34 , ilc39 , &
                ind , iref , istrt1 , itcon , iup , j , &
                jmax , k , l, lmcon , n

    ! set machine and precision dependent constants for wolfe.
    tol = ten*ten*spcmn
    tol1 = (ten**4)*spcmn
    tols = sqrt(spcmn)
    iref = 0
    violm = one/two
    lmcon = 3
    itcon = 0
    iup = 0
    s1low = ten*ten*ten*spcmn
    s1hi = one - s1low

    ! make sure s1low <= one third and s1hi >= two thirds to avoid
    ! squeezing the allowable region for s1 too tightly (or even making it
    ! empty).
    if ( s1low>one/three ) s1low = one/three
    if ( s1hi<two/three ) s1hi = two/three
    facsc = ten*ten*ten*ten
    fackp = facsc

    ilc18 = iloc(18,nparm,numgr)
    ilc28 = iloc(28,nparm,numgr)
    ilc32 = iloc(32,nparm,numgr)
    ilc34 = iloc(34,nparm,numgr)
    ilc39 = iloc(39,nparm,numgr)
    n = ndm + 1
    istrt1 = istrt
    do i = 1 , ndm
        r(i) = zero
    end do
    r(n) = one

    ! now compute the scale factor scl, whose main purpose is to avoid
    ! having all vectors in pmat with positive last component form an angle
    ! close to 90 degrees with r = (0...0 1), which can cause numerical
    ! problems.  we will compute scl = min(max(abs(a(i,k)): 1 <= i <=
    ! ndm)/b(k), b(k) >= tols, 1 <= k <= m) unless no b(k) is >=
    ! tols, in which case we set scl=1.0, or some b(k) is >= tols but
    ! scl would be < tol, in which case we set scl = tol.
100 scl = one
    ind = 0
    do k = 1 , m
        bk = pmat(n,k)
        if ( bk>=tols ) then
            quot = zero
            do i = 1 , ndm
                ab = abs(pmat(i,k))
                if ( ab>quot ) quot = ab
            end do
            quot = quot/bk
            if ( ind>0 ) then
                if ( quot>=scl ) cycle
            end if
            ind = 1
            scl = quot
        end if
    end do
300 if ( scl<tol ) scl = tol
    ! put scaled pmat into pmat1 for use in conenr.  pmat itself will remain
    ! unchanged.
400 do j = 1 , m
        do i = 1 , ndm
            pmat1(i,j) = pmat(i,j)/scl
        end do
        pmat1(n,j) = pmat(n,j)
    end do
    ! now do a normal scaling on each column of pmat1 which has an element
    ! with absolute value >= tol1.
    do j = 1 , m
        scl1 = zero
        do i = 1 , n
            ab = abs(pmat1(i,j))
            if ( ab>scl1 ) scl1 = ab
        end do
        if ( scl1>=tol1 ) then
            do i = 1 , n
                pmat1(i,j) = pmat1(i,j)/scl1
            end do
            if ( istrt1>0 ) coef(j) = wcoef(j)*scl1
            ! also put a scaled version of wcoef into coef if istrt1=1.
        else if ( istrt1>0 ) then
            coef(j) = wcoef(j)
        end if
    end do

    ! if istrt1=1, for use in conenr set coef = (-s1*scl**2)*coef, where
    ! s1 = s/(s + (1.0-s)*scl**2) is the s value in the scaled situation.
    ! note that a partly scaled version of wcoef (see loop ending with the
    ! statement numbered 190 above) is already in coef if istrt1=1.
    if ( istrt1>0 ) then
        ! if we had ncor > n, reset ncor to n.
        if ( ncor>n ) ncor = n
        fact = -(s/(s+(one-s)*scl**2))*scl**2
        do j = 1 , m
            coef(j) = fact*coef(j)
        end do
    end if

    ! call conenr to compute the nearest point to r in the cone of
    ! nonnegative linear combinations of columns of pmat1.
    call conenr(n,m,pmat1,r,istrt1,ncor,icor,tol,iwork,liwrk,work,    &
                lwrk,work(ilc39),work(ilc32),work(ilc28),nparm,numgr, &
                coef,ptnr,dist,nmaj,nmin,jflag)

    ! if jflag=3 then conenr has failed, possibly because scl was too large.
    if ( jflag/=3 ) then
        ! here jflag /= 3 and we compute s1 = 1.0 - ptnr(n).
        s1 = one - ptnr(n)
        if ( s1>=s1low ) then

            ! here jflag /= 3 and s1 >= s1low, so if also s1 <= s1hi we accept
            ! the result from conenr and move on.
            if ( s1<=s1hi ) goto 700

            ! here jflag /= 3 and s1 > s1hi, so if itcon < lmcon we try
            ! again with larger scl.
            ! if here jflag=0 and ncor=0 the nearest point to the origin in the
            ! polytope appears to be the origin so we forego adjusting scl.
            if ( jflag/=0 ) goto 600
            if ( ncor>0 ) goto 600
            goto 700
        end if
    end if
    ! here jflag=3 or s1 < s1low, so if itcon < lmcon we try again with
    ! smaller scl.
    if ( itcon<lmcon ) then

        ! here we increment itcon and if scl was not already very small we
        ! decrease it and try conenr again.
        itcon = itcon + 1
        if ( iup<0 ) then

        else if ( iup==0 ) then
            ! here iup=0 so either we are just starting (in which case we set iup=-1
            ! to indicate we are in a phase of decreasing scl) or we are oscillating.
            if ( itcon<=1 ) then
                iup = -1
            else
                facsc = sqrt(facsc)
            end if
        else
            ! here iup=1 and we have oscillation in the search for a usable scl so
            ! we replace the correction factor by its square root and reset iup to
            ! 0 to indicate oscillation.
            iup = 0
            facsc = sqrt(facsc)
        end if
        ! here we decrease scl if it was not already very small.
        if ( scl>=(one+one/ten)*tol ) then
            scl = scl/facsc
            goto 300
        end if
    end if

    ! here we were unable to get an acceptable s1 from conenr so we set
    ! jflag=4 as a warning and return.  first try again from scratch if this
    ! has not been done.
500 if ( istrt1<=0 ) then
        jflag = 4
        return
    else
        istrt1 = 0
        itcon = 0
        iref = 0
        iup = 0
        facsc = fackp
        goto 100
    end if
600 if ( itcon>=lmcon ) goto 500
    itcon = itcon + 1
    if ( iup<0 ) then
        ! here iup=-1 and we have oscillation in the search for a usable scl so
        ! we replace the correction factor by its square root and set iup=0
        ! to indicate oscillation.
        iup = 0
        facsc = sqrt(facsc)
        scl = scl*facsc
    else if ( iup==0 ) then
        ! here iup=0 so either we are just starting (in which case we set iup=1
        ! to indicate we are in a phase of increasing scl) or we are oscillating.
        if ( itcon<=1 ) then
            iup = 1
            scl = scl*facsc
        else
            facsc = sqrt(facsc)
            scl = scl*facsc
        end if
    else
        scl = scl*facsc
    end if
    goto 400

    ! here conenr may have succeeded and we compute the nearest point
    ! (wpt,s1)=r-ptnr to r from the dual of the cone described earlier.
    ! this new cone is the set of (x,t) such that (a(k)/scl,b(k)).(x,t) <=
    ! 0.0 for k=1,...,m.
700 do i = 1 , ndm
        wpt(i) = -ptnr(i)
    end do
    ! divide wpt by s1*scl.
    do i = 1 , ndm
        wpt(i) = wpt(i)/(s1*scl)
    end do
    ! compute the maximum wolfe constraint violation as a check.
800 do j = 1 , m
        v1 = pmat(n,j)
        do i = 1 , ndm
            v1 = v1 + pmat(i,j)*wpt(i)
        end do
        if ( j>1 ) then
            if ( v1<=vmax ) goto 900
        end if
        jmax = j
        vmax = v1
900 end do
    ! if vmax <= violm we reset jflag to 0 and accept the result.
    if ( vmax<=violm ) then
        jflag = 0

        ! divide the coefficients by -s1*scl**2.
        do j = 1 , m
            wcoef(j) = -coef(j)/(s1*scl**2)
        end do

        ! we now reconstruct the normal scaling factors computed in the loop
        ! ending with the statement labelled 190 in this subroutine.  in a later
        ! version of this subroutine an array may be created to store these in
        ! that loop, but for now we avoid the extra storage and programming work
        ! of fiddling with the variable dimensioning.  to recreate the factor
        ! scl1 corresponding to column j, we compute the maximum absolute value
        ! of the first ndm elements of pmat in this column, divide it by scl, take
        ! the maximum of this and abs(pmat(ndm+1,j)), and take scl1 to be this
        ! value unless it is less than tol1, in which we (in effect) take scl1=1.0.
        ! finally, since wcoef(j) was computed with the jth column of pmat divided
        ! by scl1 it contains a hidden factor of scl1, which we divide out.
        do j = 1 , m
            scl1a = zero
            do i = 1 , ndm
                ab = abs(pmat(i,j))
                if ( ab>scl1a ) scl1a = ab
            end do
            scl1 = scl1a/scl
            ab = abs(pmat(ndm+1,j))
            if ( ab>scl1 ) scl1 = ab
            if ( scl1>=tol1 ) wcoef(j) = wcoef(j)/scl1
        end do

        ! compute the s value for the unscaled situation.
        s = s1/(s1+(one-s1)/scl**2)
        ! copy wpt into ptnr to get the right dimension for dotprd and compute
        ! the distance.
        do i = 1 , ndm
            ptnr(i) = wpt(i)
        end do
        ptnr(n) = zero
        wdist = sqrt(dotprd(ndm,ptnr,ptnr,nparm))
        return
    else

        ! here vmax is too large.
        if ( iref>0 ) then
            ! here we have unsuccessfully tried to refine wpt with refwl at least
            ! once.  if ncor < ndm and the worst violation occurred outside
            ! icor we will put it in icor and try refwl again, otherwise we will
            ! set jflag=7 and return (first trying from scratch if this has not
            ! been done).
            if ( ncor>=ndm ) goto 1000
            if ( ncor>0 ) then
                do l = 1 , ncor
                    if ( jmax==icor(l) ) goto 1000
                end do
            end if
            ncor = ncor + 1
            icor(ncor) = jmax
        end if

        ! increment iref and call refwl to attempt to refine wpt, then go back
        ! and recheck the maximum constraint violation.
        iref = iref + 1
        call refwl(ndm,ncor,icor,pmat,pmat1,nparm,numgr,iwork(ilc18),  &
                   work(ilc34),wpt)
        goto 800
    end if

1000 if ( istrt1>0 ) then
        istrt1 = 0
        itcon = 0
        iref = 0
        iup = 0
        facsc = fackp
        goto 100
    end if

    jflag = 7

    end subroutine wolfe
!********************************************************************************

!********************************************************************************
!>
!  Given m n-dimensional vectors p(j) as the first m columns
!  of the matrix pmat1 and an n-vector r, this subroutine returns in
!  ptnr the nearest point to r in the cone of points summation(
!  coef(j)*p(j)), where coef(j) >= 0.0 for j=1,...,m (unless jflag > 0,
!  which indicates failure).  the number of vectors p(j) in
!  the final corral is returned in ncor with their indices in icor,
!  the distance is returned in dist, the number of major cycles (i.e.
!  adding a vector) is returned in nmaj, and the number of minor cycles
!  (i.e. removing a vector) is returned in nmin.  if the user sets
!  istrt1=0 the subroutne starts from scratch, but the user can set
!  istrt1=1 and initially specify ncor, icor, and coef (noting that ncor
!  must be <= n, and if j does not occur in icor, then coef(j) should
!  be set to 0.0.)

    subroutine conenr(n,m,Pmat1,r,Istrt1,Ncor,Icor,Tol,Iwork,Liwrk,   &
                      Work,Lwrk,Vec,Ptnrr,Picor,Nparm,Numgr,Coef,Ptnr,&
                      Dist,Nmaj,Nmin,Jflag)

    implicit none

    real(wp) :: amax , amin , cjj , Coef , diff , Dist , dmax ,   &
                dp , dsq , omt , pdotj , Picor ,     &
                Pmat1 , Ptnr , Ptnrr , quot
    real(wp) :: r , theta , Tol , tst , Vec ,  &
                Work
    integer :: i , Icor , icoro , ihouse , ii , ilc01 , ilc03 , ilc04 ,  &
               ilc09 , ilc23 , ilc34 , Istrt1 , itst1 , Iwork ,   &
               j , Jflag , jj , jmax , jmin
    integer :: kntsl , l , limsl , Liwrk , ll , Lwrk , m , mincf , mp1 , &
               n , Ncor , ncoro , ndm , Nmaj , Nmin , Nparm , Numgr

    dimension Pmat1(Nparm+1,Numgr) , r(Nparm+1) , Icor(Nparm+1) ,     &
              Coef(Numgr) , Ptnr(Nparm+1) , Vec(Nparm+1) ,            &
              Ptnrr(Nparm+1) , Picor(Nparm+1,Nparm+1) , Iwork(Liwrk) ,&
              Work(Lwrk)

    real(wp),parameter :: tolel = ten*ten*spcmn
    real(wp),parameter :: z1 = ten*tolel
    real(wp),parameter :: z2 = ten*z1
    real(wp),parameter :: z3 = ten*z1

    ilc01 = iloc(1,Nparm,Numgr)
    ilc03 = iloc(3,Nparm,Numgr)
    ilc04 = iloc(4,Nparm,Numgr)
    ilc09 = iloc(9,Nparm,Numgr)
    ilc23 = iloc(23,Nparm,Numgr)
    ilc34 = iloc(34,Nparm,Numgr)
    kntsl = 0
    limsl = 100
    mp1 = m + 1
    ndm = n - 1
    Nmaj = 0
    Nmin = 0
    Jflag = 0
    itst1 = 0
    ncoro = -1

    if ( Istrt1>0 ) goto 200

    ! HERE ISTRT1=0 SO WE START FROM SCRATCH.  FIND THE INDEX JMAX FOR
    ! WHICH (P(J).R)/SQRT(P(J).P(J)) IS MAXIMIZED FOR P(J).P(J) > Z1.
100 amax = zero
    jmax = 0
    do j = 1 , m
        do i = 1 , n
            Vec(i) = Pmat1(i,j)
        end do
        pdotj = dotprd(n,Vec,Vec,Nparm)
        if ( pdotj>z1 ) then
            quot = dotprd(n,Vec,r,Nparm)/sqrt(pdotj)
            if ( quot>amax ) then
                amax = quot
                jmax = j
            end if
        end if
    end do
    if ( jmax>0 ) then
        ! IF AMAX IS NOT SIGINFICANTLY POSITIVE WE PROCEED AS IF IT WERE ZERO.
        if ( amax*sqrt(ndm+one)>tolel ) then
            ! HERE WE FOUND THE RAY CLOSEST TO R AND WE COMPLETE THE
            ! INITIALIZATION BY SETTING NCOR=1, ICOR(1)=JMAX, AND COEF(JMAX)=1.0
            ! (WITH ALL OTHER ENTRIES OF COEF EQUAL TO ZERO).
            Ncor = 1
            Icor(1) = jmax
            do i = 1 , m
                Coef(i) = zero
            end do
            Coef(jmax) = one
            goto 200
        end if
    end if

    ! HERE THERE WERE NO VECTORS P(J) WHICH HAVE BOTH LENGTH SQUARED
    ! GREATER THAN Z1 AND ANGLE WITH R SIGNIFICANTLY LESS THAN 90 DEGREES,
    ! AND WE SET NCOR=0, PTNR=THE ZERO VECTOR, COEF=THE ZERO VECTOR, DIST=
    ! THE LENGTH OF R, AND WE RETURN.
    Ncor = 0
    do i = 1 , n
        Ptnr(i) = zero
    end do
    do j = 1 , m
        Coef(j) = zero
    end do
    Dist = sqrt(dotprd(n,r,r,Nparm))
    return

    ! SET PTNR TO THE CURRENT NEAREST POINT.  FIRST ZERO IT OUT.
200 do i = 1 , n
        Ptnr(i) = zero
    end do
    if ( Ncor>0 ) then
        ! HERE NCOR > 0 AND WE SET PTNR=SUMMATION(COEF(J)*P(J)).
        do j = 1 , Ncor
            jj = Icor(j)
            cjj = Coef(jj)
            do i = 1 , n
                Ptnr(i) = Ptnr(i) + cjj*Pmat1(i,jj)
            end do
        end do
    end if

    ! PUT PTNR-R INTO PTNRR AND COMPUTE THE DISTANCE FROM PTNR TO R.
    do i = 1 , n
        Ptnrr(i) = Ptnr(i) - r(i)
    end do
    dsq = dotprd(n,Ptnrr,Ptnrr,Nparm)
    Dist = sqrt(dsq)

    ! NOW CHECK OPTIMALITY.
    ! FIRST SEE WHETHER THE HYPERPLANE THROUGH PTNR PERPENDICULAR TO
    ! PTNR-R PASSES THROUGH THE ORIGIN.  IF NCOR=0 THIS WILL
    ! AUTOMATICALLY BE TRUE SINCE THEN PTNR IS THE ORIGIN.  IF IT IS NOT
    ! TRUE WE GO DOWN TO SOLVE FOR A NEW NEAREST POINT IN THE SUBSPACE
    ! DETERMINED BY THE CURRENT ICOR.
    if ( Ncor>0 ) then
        tst = dotprd(n,Ptnr,Ptnrr,Nparm)
        if ( abs(tst)>=z1 ) goto 300
    end if
    ! HERE THE HYPERPLANE ROUGHLY PASSES THROUGH THE ORIGIN, AND WE
    ! CHECK WHETHER ALL P(J) VECTORS ARE ROUGHLY SEPARATED FROM R BY IT.
    ! PUT THE MINIMUM OF (PTNR-R).(P(J)-R) IN AMIN AND THE INDEX WHERE IT
    ! IS ACHIEVED IN JMIN.
    do i = 1 , n
        Vec(i) = Pmat1(i,1) - r(i)
    end do
    jmin = 1
    amin = dotprd(n,Ptnrr,Vec,Nparm)
    if ( m>1 ) then
        do j = 2 , m
            do i = 1 , n
                Vec(i) = Pmat1(i,j) - r(i)
            end do
            dp = dotprd(n,Ptnrr,Vec,Nparm)
            if ( dp<amin ) then
                amin = dp
                jmin = j
            end if
        end do
    end if

    ! FOR TESTING PURPOSES COMPUTE THE MAXIMUM OF THE SQUARES OF THE
    ! LENGTHS OF THE DISTANCES CONSIDERED.
    do i = 1 , n
        Vec(i) = Pmat1(i,jmin) - r(i)
    end do
    dmax = dotprd(n,Vec,Vec,Nparm)
    if ( Ncor>0 ) then
        do j = 1 , Ncor
            jj = Icor(j)
            do i = 1 , n
                Vec(i) = Pmat1(i,jj) - r(i)
            end do
            dp = dotprd(n,Vec,Vec,Nparm)
            if ( dp>dmax ) dmax = dp
        end do
    end if
    ! DO THE TEST.  IF IT IS SUCCESSFUL, THEN WE HAVE (APPROXIMATE)
    ! OPTIMALITY AND WE RETURN.
    if ( amin-dsq+z1*dmax<0 ) then
        ! HERE PTNR IS NOT OPTIMAL.  AS A CHECK AGAINST BLUNDERS WE MAKE SURE
        ! NCOR < N AND JMIN IS NOT IN ICOR.
        if ( Ncor>0 ) then
            if ( Ncor<n ) then
                do l = 1 , Ncor
                    if ( jmin==Icor(l) ) goto 220
                end do
                goto 250
            end if
            ! HERE WE HAVE BLUNDERED SO WE SET JFLAG=1 AS A WARNING, COMPUTE DIST,
            ! AND RETURN.  FIRST TRY FROM SCRATCH IF THIS HAS NOT BEEN DONE.
220         if ( Istrt1+Jflag<=0 ) then
                Jflag = 1
                ! WRITE(6,'(A)') '*****JFLAG IS 1 IN CONENR'
                return
            else
                Jflag = -1
                kntsl = 0
                goto 100
            end if
        end if
        ! HERE PTNR IS NOT OPTIMAL, NCOR < N, AND JMIN IS NOT IN ICOR.
        ! WE INCREMENT THE MAJOR CYCLE COUNTER AND ADD P(JMIN).
250     Nmaj = Nmaj + 1
        Ncor = Ncor + 1
        Icor(Ncor) = jmin
        Coef(jmin) = zero
    else
        return
    end if

    ! CHECK TO SEE IF WE HAVE SOLVED THE SYSTEM BELOW LIMSL TIMES ALREADY,
    ! AND IF SO, SET JFLAG=6 AS A WARNING AND RETURN (BUT
    ! TRY FROM SCRATCH BEFORE GIVING UP IF THIS HAS NOT ALREADY BEEN DONE).
300 if ( kntsl<limsl ) then

        ! CHECK TO SEE IF NCOR AND THE LAST ELEMENT IN ICOR ARE UNCHANGED FROM THE
        ! PREVIOUS HOUSE CALL (HA HA), WHICH INDICATES FAILURE.  NOTE THAT HERE WE
        ! MUST HAVE NCOR > 0.
        if ( Ncor/=ncoro ) then
            ncoro = Ncor
        else if ( Icor(Ncor)==icoro ) then
            ! HERE WE HAVE CYCLING AND WE SET JFLAG=2 AS A WARNING AND RETURN.  FIRST
            ! TRY FROM SCRATCH IF THIS HAS NOT BEEN DONE.
            if ( Istrt1+Jflag<=0 ) then
                Jflag = 2
                return
            else
                Jflag = -1
                kntsl = 0
                goto 100
            end if
        end if
        icoro = Icor(Ncor)
        kntsl = kntsl + 1

        ! NOW WE SOLVE THE SYSTEM PICOR*VEC = R IN THE LEAST SQUARES
        ! SENSE FOR THE COEFFICIENT VECTOR VEC (RELATIVE TO
        ! ICOR) OF THE NEAREST POINT TO R IN THE SUBSPACE SPANNED BY
        ! P(ICOR(1)),...,P(ICOR(NCOR)), WHERE P(ICOR) IS THE N X NCOR MATRIX
        ! WHOSE COLUMNS ARE THESE VECTORS.
        ! NOW FILL IN PICOR AND CALL HOUSE TO COMPUTE VEC.
        do j = 1 , Ncor
            jj = Icor(j)
            do i = 1 , n
                Picor(i,j) = Pmat1(i,jj)
            end do
        end do

        call house(n,Ncor,Picor,r,Iwork(ilc23),Nparm,Work(ilc01),      &
                   Work(ilc04),Work(ilc09),Work(ilc34),Work(ilc03),Vec,&
                   ihouse)

        ! IF HOUSE FAILS (INDICATED BY IHOUSE=1) WE SET JFLAG=3 AS A
        ! WARNING AND RETURN.  FIRST TRY FROM SCRATCH IF THIS HAS NOT BEEN DONE.
        if ( ihouse<=0 ) then
            ! CHECK TO SEE IF ALL THE COEFFICIENTS IN VEC ARE > Z2, AND IF SO,
            ! PUT VEC INTO COEF AND GO BACK TO COMPUTE PTNR.  THE COEFFICIENTS IN
            ! COEF NOT CORRESPONDING TO THOSE IN VEC WILL REMAIN EQUAL TO ZERO.
            do i = 1 , Ncor
                if ( Vec(i)<=z2 ) goto 350
            end do
            do i = 1 , Ncor
                ii = Icor(i)
                Coef(ii) = Vec(i)
            end do
            goto 200
        else if ( Istrt1+Jflag<=0 ) then
            Jflag = 3
            return
        else
            Jflag = -1
            kntsl = 0
            goto 100
        end if

        ! HERE SOME ELEMENT OF VEC IS <= Z2.  COMPUTE THETA=MIN(1.0, MIN(
        ! COEF(ICOR(I))/(COEF(ICOR(I))-VEC(I)): COEF(ICOR(I))-VEC(I) > Z3)).
350     theta = one
        do l = 1 , Ncor
            ll = Icor(l)
            diff = Coef(ll) - Vec(l)
            if ( diff>z3 ) then
                quot = Coef(ll)/diff
                if ( quot<theta ) theta = quot
            end if
        end do
        ! COMPUTE THE NEW COEF AS (1.0-THETA)*COEF+THETA*VEC.
        omt = one - theta
        do l = 1 , Ncor
            ll = Icor(l)
            Coef(ll) = omt*Coef(ll) + theta*Vec(l)
        end do
        ! COMPUTE THE INDEX MINCF (RELATIVE TO ICOR) OF THE SMALLEST ELEMENT OF
        ! COEF AND SET ALL ELEMENTS OF COEF WHICH ARE <= Z2 TO ZERO.
        mincf = 0
        amin = z2
        do i = 1 , Ncor
            ii = Icor(i)
            if ( Coef(ii)<=z2 ) then
                if ( Coef(ii)<=amin ) then
                    amin = Coef(ii)
                    mincf = i
                end if
                Coef(ii) = zero
            end if
        end do

        if ( mincf<=0 ) then
            ! HERE MINCF=0 AND AN UNLIKELY BLUNDER HAS OCCURRED.  THIS MUST BE DUE TO
            ! ROUNDOFF ERROR SINCE IN THEORY (NEW) COEF(ICOR(I)) MUST BE <= Z2
            ! FOR SOME I=1,...,NCOR, WHICH MAKES MINCF > 0 IN THE LAST LOOP.
            ! TO SEE THIS, FIRST NOTE THAT FOR SOME IBAR=1,...,NCOR, VEC(IBAR) MUST
            ! BE <= Z2 SINCE OTHERWISE WE WOULD NOT BE HERE.  BY ITS DEFINITION,
            ! THETA MUST BE <= 1.0.  IF THETA = 1.0, THEN (NEW) COEF(ICOR(IBAR))
            ! = (1.0 - THETA)*(OLD) COEF(ICOR(IBAR)) + THETA*VEC(IBAR) = VEC(IBAR)
            ! <= Z2.  IF ON THE OTHER HAND THETA < 1.0, THEN FOR SOME ISTAR=1,
            ! ...,ICOR WE HAVE (OLD) COEF(ICOR(ISTAR)) - VEC(ISTAR) >= Z3 AND
            ! THETA = (OLD) COEF(ICOR(ISTAR))/((OLD) COEF(ICOR(ISTAR)) - VEC(ISTAR)),
            ! SO (NEW) COEF(ICOR(ISTAR)) = (1.0 - THETA)*(OLD) COEF(ICOR(ISTAR)) +
            ! THETA*VEC(ISTAR) = (-VEC(ISTAR)*(OLD) COEF(ICOR(ISTAR)) + (OLD)
            ! COEF(ICOR(ISTAR))*VEC(ISTAR))/((OLD) COEF(ICOR(ISTAR)) - VEC(ISTAR))
            ! = 0.0.  NOTE THAT WE HAVE Z2 >= 0.0 AND Z3 > 0.0.
            ! TO CORRECT THIS BLUNDER WE SET MINCF = AN INDEX I FOR WHICH (NEW)
            ! COEF(ICOR(I)) IS MINIMIZED AND SET COEF(ICOR(I)) = 0.0.
            do i = 1 , Ncor
                ii = Icor(i)
                if ( i>1 ) then
                    if ( Coef(ii)>=amin ) cycle
                end if
                amin = Coef(ii)
                mincf = i
            end do
            ii = Icor(mincf)
            Coef(ii) = zero
        end if
        ! INCREMENT THE MINOR ITERATION COUNTER NMIN, REMOVE ICOR(MINCF),
        ! AND DECREMENT NCOR.
        Nmin = Nmin + 1
        do l = 1 , Ncor
            if ( l>mincf ) Icor(l-1) = Icor(l)
        end do
        Ncor = Ncor - 1
        ! GO BACK TO COMPUTE PTNR.
        goto 200
    else if ( Istrt1+Jflag<=0 ) then
        Jflag = 6
        ! WRITE(6,'(A)') '*****JFLAG IS 6 IN CONENR'
        return
    else
        Jflag = -1
        kntsl = 0
        goto 100
    end if

    end subroutine conenr
!********************************************************************************

!********************************************************************************
!>
!  Given ncor n dimensional vectors as columns of the n by ncor
!  matrix picor and an n dimensional vector r, this subroutine uses
!  householder transformations to find the best least squares solution
!  vec to the linear system of equations picor*vec = r, where vec
!  is an ncor dimensional vector.  if the rank of picor is
!  (computationally) 0, the subroutine will return with the failure
!  warning ihouse=1, otherwise it will return with ihouse=0.  if the
!  rank is > 0 but < ncor, then (ncor - rank) of the components
!  of vec will be set to 0.0.  the arays picor and r will not be
!  changed by this subroutine.  the subroutine will attempt up to
!  numref iterative refinements of the solution, where the user can
!  set numref as any nonnegative integer, but to get the most out of
!  the iterative refinement process, the computation of the residual
!  summ near the end of this subroutine should be done in higher
!  precision than the other computations in the subroutine.

    subroutine house(n,Ncor,Picor,r,Kpivot,Nparm,Aa,Beta,d,Save,b,Vec,Ihouse)

    implicit none

    integer,intent(in)  :: Nparm
    integer,intent(out) :: Ihouse
    integer,intent(in)  :: n
    integer,intent(in)  :: Ncor
    integer  :: Kpivot(Nparm+1)
    real(wp) :: Aa(Nparm+1,Nparm+1)
    real(wp) :: b(Nparm+1)
    real(wp) :: Beta(Nparm+1)
    real(wp) :: d(Nparm+1)
    real(wp) :: Picor(Nparm+1,Nparm+1)
    real(wp) :: r(Nparm+1)
    real(wp) :: Save(Nparm+1)
    real(wp) :: Vec(Nparm+1)

    real(wp) :: aakk , amax , sqdk , store , sum , summ , test , testt
    integer :: i , ia , icount , ii , j , jj , k , kchnge , kk , &
               kp , krank , kt

    real(wp),parameter :: tolsq = (ten*ten*spcmn)**2
    integer,parameter :: numref = 1 !! Set numref = the limit on the number of iterative refinement steps.
    integer,parameter :: nmref1 = numref + 1
    integer,parameter :: nmref2 = numref + 2

    ! COMPUTE MACHINE AND PRECISION DEPENDENT CONSTANTS.
    Ihouse = 0

    ! SET KRANK = MIN(N,NCOR).  THIS MAY BE REDUCED LATER.
    krank = Ncor
    if ( n<Ncor ) krank = n
    ! INITIALLY SET KPIVOT.  AFTER ALL COLUMN INTERCHANGES ARE DONE
    ! KPIVOT(J) WILL BE THE ORIGINAL POSITION OF THE COLUMN WHERE THE
    ! JTH PIVOT WAS DONE.  THIS COLUMN WILL BE MOVED TO COLUMN J.
    do j = 1 , Ncor
        Kpivot(j) = j
    end do
    ! COPY R INTO B AND PICOR INTO AA, BUT IN THE PROCESS REPLACE ANY NUMBERS
    ! WITH ABSOLUTE VALUE LESS THAN SPCMN BY ZERO TO AVOID UNDERFLOWS.
    do i = 1 , n
        if ( abs(r(i))<spcmn ) then
            b(i) = zero
        else
            b(i) = r(i)
        end if
    end do
    do j = 1 , Ncor
        do i = 1 , n
            if ( abs(Picor(i,j))<spcmn ) then
               Aa(i,j) = zero
            else
               Aa(i,j) = Picor(i,j)
            end if
        end do
    end do
    do k = 1 , Ncor
        if ( k>n ) exit
        d(k) = zero
        kchnge = k
        do jj = k , Ncor
            sum = zero
            do ia = k , n
                if ( abs(Aa(ia,jj))>spcmn ) sum = sum + Aa(ia,jj)*Aa(ia,jj)
            end do
            if ( d(k)<sum ) then
                kchnge = jj
                d(k) = sum
            end if
        end do

        !  KCHNGE CONTAINS THE INDEX OF THE COLUMN OF GREATEST
        !  LENGTH BETWEEN K AND NCOR (FROM POSITION K TO THE BOTTOM).
        ! IF K=1 AND D(K) < TOLSQ WE RETURN WITH THE FAILURE WARNING
        ! IHOUSE=1.
        if ( k<=1 ) then
            if ( d(k)<tolsq ) then
                Ihouse = 1
                return
            end if
        end if

        if ( kchnge/=k ) then
            !  START COLUMN INTERCHANGE.
            do i = 1 , n
                store = Aa(i,kchnge)
                Aa(i,kchnge) = Aa(i,k)
                Aa(i,k) = store
            end do
            kk = Kpivot(k)
            Kpivot(k) = Kpivot(kchnge)
            Kpivot(kchnge) = kk
        end if
        if ( k/=1 ) then
            amax = abs(d(1))
            test = (real(n-k+1,wp)*(ten*ten*spcmn)**2)*(amax*amax)
            if ( abs(d(k))<=test ) then
                ! HERE THE LENGTH OF THE BEST OF COLUMNS K THROUGH NCOR (FROM K DOWN)
                ! WAS TOO SMALL, AND WE REDUCE KRANK TO K-1 AND LEAVE THIS LOOP.
                d(k) = sqrt(d(k))
                krank = k - 1
                exit
            end if
        end if

        ! NOW COMPUTE THE SCALAR BETA(K) AND THE N-K+1 DIMENSIONAL VECTOR
        ! GNU(K) (TO BE PLACED IN AA(K,K),...,AA(N,K)) FOR I(K) - BETA(K)*
        ! GNU(K)*(GNU(K) TRANSPOSE), WHICH IS THE ACTIVE PART OF THE
        ! HOUSEHOLDER TRANSFORMATION PH(K) = DIAG(I(K-1), ACTIVE PART).  THIS
        ! IS A SYMMETRIC ORTHOGONAL MATRIX WHICH WHEN MULTIPLIED TIMES AA WILL
        ! ZERO OUT AA(K+1,K),...,AA(N,K) AND CHANGE AA(K,K) TO -SGN(OLD
        ! AA(K,K))*SQDK, WHERE SQDK = LENGTH OF OLD (AA(K,K),...,AA(N,K)) AND
        ! WE REDEFINE THE SGN FUNCTION TO HAVE VALUE 1.0 IF ITS ARGUMENT IS
        ! 0.0.  WE WILL HAVE BETA(K) = 1.0/(SQDK**2 + ABS(OLD AA(K,K))*SQDK)
        ! AND GNU(K) = (OLD AA(K,K) + SGN(OLD AA(K,K))*SQDK, OLD AA(K+1,K),...,
        ! OLD AA(N,K)).  WE WILL ALSO REPLACE D(K) BY THE NEW AA(K,K) (WHICH
        ! WILL NOT ACTUALLY BE WRITTEN INTO AA) FOR LATER USE.
        aakk = Aa(k,k)
        sqdk = sqrt(d(k))
        if ( aakk<zero ) then
            Beta(k) = one/(d(k)-aakk*sqdk)
            Aa(k,k) = -sqdk + aakk
            d(k) = sqdk
        else
            Beta(k) = one/(d(k)+aakk*sqdk)
            Aa(k,k) = sqdk + aakk
            d(k) = -sqdk
        end if
        kt = k + 1
        if ( k/=Ncor ) then
            ! HERE K < NCOR AND WE MULTIPLY COLUMNS K+1,...,NCOR OF AA BY THE
            ! HOUSEHOLDER TRANSFORMATION PH(K), WHICH WILL CHANGE ONLY POSITIONS
            ! K THROUGH THE BOTTOM OF THESE COLUMNS.  THIS IS DONE BY, FOR J =
            ! K+1,...,NCOR, REPLACING COLUMN J (FROM K DOWN) BY COLUMN J (FROM K DOWN)
            ! - GNU(K)*(GNU(K).COLUMN J (FROM K DOWN))*BETA(K).
            do j = kt , Ncor
                Save(j) = zero
                do ia = k , n
                    Save(j) = Save(j) + Aa(ia,k)*Aa(ia,j)
                end do
                do i = k , n
                    Aa(i,j) = Aa(i,j) - Aa(i,k)*Save(j)*Beta(k)
                end do
            end do
        end if
    end do

    do i = 1 , krank
        ! IF I <= MIN(KRANK,NCOR-1), DIVIDE ROW I OF AA FROM COLUMN I+1
        ! THROUGH COLUMN NCOR BY THE NEW AA(I,I) (WHICH IS NOT ACTUALLY
        ! WRITTEN INTO AA(I,I), BUT IS STORED IN D(I)).
        ii = i + 1
        if ( i==Ncor ) exit
        do j = ii , Ncor
            Aa(i,j) = Aa(i,j)/d(i)
        end do
    end do

    ! NOW ALL THE DIAGONAL ELEMENTS OF AA (ALTHOUGH NOT WRITTEN IN)
    ! ARE 1.0 AND ALL OFF DIAGONAL ELEMENTS OF AA ARE LESS THAN OR
    ! EQUAL TO 1.0.

    ! INITIALIZE THE ITERATIVE REFINEMENT COUNTER ICOUNT AND ZERO OUT VEC
    ! INITIALLY.  THE VEC VALUES NOT CORRESPONDING TO THE FIRST KRANK
    ! COLUMNS (MODULO EARLIER COLUMN INTERCHANGES) WILL REMAIN AT 0.0.
    icount = 1
    do i = 1 , Ncor
        Vec(i) = zero
    end do

    iteration :  do

        ! PREMULTIPLY B BY THE HOUSEHOLDER TRANSFORMATIONS PH(1),...,
        ! PH(KRANK).  RECALL THAT GNU(I) IS STILL IN AA(I,I),...,AA(N,I)
        ! FOR I=1,...,KRANK.
        do i = 1 , krank
            sum = zero
            do ia = i , n
                sum = sum + Aa(ia,i)*b(ia)
            end do
            sum = sum*Beta(i)
            do j = i , n
                b(j) = b(j) - Aa(j,i)*sum
            end do
        end do

        ! NOW ONLY USE THE FIRST KRANK TERMS OF B, AS WE CANT DO ANYTHING ABOUT
        ! THE OTHERS, WHOSE SQUARE ROOT OF SUM OF SQUARES WILL GIVE THE LEAST
        ! SQUARES DISTANCE.
        ! DIVIDE B(I) BY D(I) FOR I=1,...,KRANK AS WE DID THIS TO ROW I OF AA.
        do i = 1 , krank
            b(i) = b(i)/d(i)
        end do

        ! THE PROBLEM HAS NOW BEEN REDUCED TO SOLVING (UPPER LEFT KRANK BY
        ! KRANK PART OF AA)*(FIRST KRANK TERMS OF VEC, MODULO COLUMN
        ! INTERCHANGE UNSCRAMBLING) = (FIRST KRANK TERMS OF B).  ALTHOUGH THE
        ! DIAGONAL AND BELOW DIAGONAL TERMS OF THE COEFFICIENT MATRIX HAVE NOT
        ! BEEN WRITTEN IN, THE SYSTEM IS UPPER TRIANGULAR WITH DIAGONAL ELEMENTS
        ! ALL EQUAL TO 1.0, SO WE SOLVE BY BACK SUBSTITUTION.  WE FIRST PUT
        ! THE SOLUTION TO THIS SYSTEM IN B(1),...,B(KRANK) AND SORT IT OUT
        ! LATER.  IF ICOUNT > 1 THE SOLUTION IS AN ITERATIVE CORRECTION TO
        ! VEC RATHER THAN VEC ITSELF.
        do ii = 1 , krank
            i = krank + 1 - ii
            kk = i - 1
            if ( i/=1 ) then
                ! HERE WE ALREADY HAVE B(I) (WHERE I  > 1) AND WE SUBTRACT AA(J,I)*
                ! B(I) FROM B(J) FOR J = 1,...,I-1.
                do j = 1 , kk
                    b(j) = b(j) - Aa(j,i)*b(i)
                end do
            end if
        end do

        !  TEST FOR CONVERGENCE.
        !  FIRST TEST, TOO MANY ITERATIONS.
        !  SECOND TEST, SEE IF VEC IS DECREASING.

        ! COMPUTE THE LENGTH SQUARED OF THE FIRST TOP 1 THROUGH KRANK PART OF
        ! B, WHICH WILL BE THE RESIDUAL VECTOR IF ICOUNT > 1.
        sum = zero
        do i = 1 , krank
            if ( abs(b(i))>spcmn ) sum = sum + b(i)*b(i)
        end do
        if ( icount==1 ) then
            testt = sum
        else if ( sum>test/two ) then
            icount = nmref2
        end if
        test = sum

        ! COMPUTE THE VEC VALUES, WHICH WILL BE ACTUAL VEC VALUES IF ICOUNT=1
        ! AND CORRECTIONS TO VEC VALUES IF ICOUNT > 1.  WE GET THESE BY
        ! UNSCRAMBLING THE B VALUES AND ADDING THEM TO THE APPROPRIATE OLD VEC
        ! VALUES (WHICH WILL BE 0.0 IF ICOUNT=1).
        do i = 1 , krank
            kp = Kpivot(i)
            Vec(kp) = b(i) + Vec(kp)
        end do

        ! CALCULATE THE RESIDUAL R - ACOEF*VEC.  RECALL THAT ACOEF AND R
        ! CONTAIN THE ORIGINAL COEFFICIENT AND RIGHT SIDE ARRAYS RESPECTIVELY.
        ! TO GET THE MOST OUT OF ITERATIVE REFINEMENT THIS COMPUTATION SHOULD
        ! PROBABLY BE DONE IN HIGHER PRECISION, IN WHICH CASE IT MAY BE
        ! FRUITFUL TO ALSO SET NUMREF LARGER AT THE BEGINNING OF THIS
        ! SUBROUTINE.
        do i = 1 , n
            summ = zero
            do j = 1 , Ncor
                if ( abs(Picor(i,j))>=spcmn ) summ = summ + Picor(i,j)*Vec(j)
            end do
            b(i) = r(i) - summ
        end do

        !  THIRD TEST, WAS THE CORRECTION SIGNIFICANT.
        if ( test>=spcmn*testt ) then
            if ( icount/=nmref1 ) then
                if ( icount<nmref2 ) then
                    icount = icount + 1
                    cycle iteration
                end if
            end if
        end if

        exit iteration  ! done

    end do iteration

    end subroutine house
!********************************************************************************

!********************************************************************************
!>
!  This subprogram computes the dot product of vectors vec1
!  and vec2 of length lgth.
!  vec1 and vec2 do not appear in function iloc since they are used only
!  as input names for this subprogram, and so they don't need to have
!  space reserved for them in the array work.

    pure function dotprd(Lgth,Vec1,Vec2,Nparm) result(dd)

    implicit none

    integer,intent(in) :: Lgth
    integer,intent(in) :: Nparm
    real(wp),dimension(Nparm+1),intent(in) :: Vec1
    real(wp),dimension(Nparm+1),intent(in) :: Vec2
    real(wp) :: dd

    integer :: j

    dd = Vec1(1)*Vec2(1)
    if ( Lgth>1 ) then
        do j = 2 , Lgth
            dd = dd + Vec1(j)*Vec2(j)
        end do
    end if

    end function dotprd
!********************************************************************************

!********************************************************************************
!>
!  This subroutine attempts to refine the ndm dimensional vector wpt
!  produced by wolfe by directly solving the system
!  summation(pmat(i,j)*wpt(i), i=1,...,ndm) = -pmat(ndm+1,j) for j =
!  icor(l), l=1,...,ncor.
!  nresl resolvents are chosen by total pivoting.  if nresl < ndm then
!  the remaining ndm-nresl elements of wpt are kept form the old wpt.
!  itrlm steps of iterative refinement are attempted at the end.

    subroutine refwl(Ndm,Ncor,Icor,Pmat,Pmat1,Nparm,Numgr,Ixrct,Save,Wpt)

    implicit none

    integer,intent(in) :: Nparm
    integer,intent(in) :: Numgr
    integer,intent(in) :: Ncor
    integer,intent(in) :: Ndm
    integer :: Icor(Nparm+1)
    integer :: Ixrct(2*Nparm)
    real(wp) :: Pmat1(Nparm+1,Numgr)
    real(wp) :: Pmat(Nparm+1,Numgr)
    real(wp) :: Save(Nparm)
    real(wp) :: Wpt(Nparm+1)

    real(wp) :: aa , amax , fact , wrst , wrsto
    integer :: i , imax , itrct , itrlm , j , jmax , jstrt , &
               k , kcol , kk , kp1 , l , maxrs , n , nresl

    real(wp),parameter :: tole = spcmn

    ! COMPUTE MACHINE AND PRECISION DEPENDENT CONSTANTS.
    itrlm = 2
    itrct = 0
    nresl = 0
    n = Ndm + 1

    ! IF NCOR=0 WE HAVE NOTHING TO DO SO WE RETURN.
    if ( Ncor>0 ) then

        ! COPY COLUMN ICOR(L) OF PMAT WITH THE SIGN OF THE LAST ELEMENT REVERSED
        ! INTO COLUMN L OF THE WORK MATRIX PMAT1 FOR L=1,...,NCOR.
        do l = 1 , Ncor
            j = Icor(l)
            do i = 1 , Ndm
                Pmat1(i,l) = Pmat(i,j)
            end do
            Pmat1(n,l) = -Pmat(n,j)
        end do

        ! NOW COLUMN REDUCE PMAT1.  NOTE THAT PMAT1 IS THE TRANSPOSE OF THE USUAL
        ! AUGMENTED MATRIX FOR SOLVING A LINEAR SYSTEM OF EQUATONS.
        ! THERE WILL BE AT MOST MAXRS = MIN(NDM,NCOR) RESOLVENTS.
        maxrs = Ncor
        if ( Ndm<maxrs ) maxrs = Ndm
        do k = 1 , maxrs
            ! SEARCH FOR THE INDICES IMAX AND JMAX WITH 1 <= IMAX <= NDM, 1 <=
            ! JMAX <= NCOR, PMAT1(IMAX,JMAX) IS NOT IN THE ROW OR COLUMN OF ANY
            ! OTHER RESOLVENT (I.E. PIVOT), AND ABS(PMAT1(IMAX,JMAX)) IS MAXIMIZED.
            ! WE USE THE VECTOR IXRCT TO SAVE THE RESOLVENT POSITIONS TO SAVE SPACE.
            jstrt = 0
            do j = 1 , Ncor
                if ( nresl>0 ) then
                    do l = 1 , nresl
                        if ( j==Ixrct(2*l) ) goto 20
                    end do
                end if
                ! HERE THERE IS NO EARLIER RESOLVENT IN COLUMN J.
                do i = 1 , Ndm
                    if ( nresl>0 ) then
                        do l = 1 , nresl
                            if ( i==Ixrct(2*l-1) ) goto 10
                        end do
                    end if
                    ! HERE THERE IS NO EARLIER RESOLVENT IN ROW I.
                    aa = abs(Pmat1(i,j))
                    if ( jstrt<=0 ) then
                        jstrt = 1
                    else if ( aa<=amax ) then
                        goto 10
                    end if
                    amax = aa
                    imax = i
                    jmax = j
10              end do
20          end do
            ! IF THE ABSOLUTE VALUE OF THIS RESOLVENT IS VERY SMALL WE DO NOT ATTEMPT
            ! ANY FURTHER COLUMN OPERATIONS.
            if ( amax<tole ) goto 50
            ! INCREMENT NRESL AND PUT THE LOCATION OF THE NRESLTH RESOLVENT IN
            ! (IXRCT(2*L-1),IXRCT(2*L)).
            nresl = nresl + 1
            Ixrct(2*nresl-1) = imax
            Ixrct(2*nresl) = jmax

            ! NOW ELIMINATE WPT(IMAX) FROM THOSE COLUMNS WHICH DO NOT CONTAIN ANY OF
            ! THE RESOLVENTS FOUND SO FAR (INCLUDING THE PRESENT RESOLVENT).
            do j = 1 , Ncor
                do l = 1 , nresl
                    if ( j==Ixrct(2*l) ) goto 40
                end do
                ! HERE COLUMN J DOES NOT CONTAIN ANY OF THE RESOLVENTS FOUND SO FAR, AND
                ! WE COMPUTE THE FACTOR FOR THE COLUMN OPERATION NEEDED TO ZERO OUT
                ! PMAT1(IMAX,J) (ALTHOUGH WE DO NOT ACTUALLY WRITE IN THE ZERO).
                fact = Pmat1(imax,j)/Pmat1(imax,jmax)
                ! NOW DO THE OPERATION IN COLUMN J FOR ALL ROWS NOT CONTAINING A
                ! RESOLVENT.  THE ELEMENTS IN THIS COLUMN IN THE ROWS WHICH CONTAIN AN
                ! EARLIER (OR PRESENT) RESOLVENT WILL NOT BE NEEDED LATER.
                do i = 1 , n
                    do l = 1 , nresl
                        if ( i==Ixrct(2*l-1) ) goto 30
                    end do
                    Pmat1(i,j) = Pmat1(i,j) - fact*Pmat1(i,jmax)
30              end do
40          end do
        end do
        ! END OF COLUMN REDUCTION OF PMAT1.

        ! IF NRESL=0 THEN ALL THE ELEMENTS IN PMAT1 FOR 1 <= I <= NDM AND
        ! 1 <= J <= NCOR WERE VERY SMALL IN ABSOLUTE VALUE, AND THERE IS
        ! NOTHING WE CAN DO, SO WE RETURN.
50      if ( nresl>0 ) goto 200
    end if

    return

    ! NOW DO BACK SUBSTITUTION TO COMPUTE, FOR K=NRESL,...,1,
    ! WPT(IXRCT(2*K-1)) = (PMAT1(NDM+1,IXRCT(2*K)) - SUMMATION(
    ! PMAT1(I,IXRCT(2*K))*WPT(I), FOR I = 1,...,NDM, I /= IXRCT(2*L-1)
    ! FOR ANY L=1,...,K))/PMAT1(IXRCT(2*K-1),IXRCT(2*K)).  IF WE ARE IN AN
    ! ITERATIVE REFINEMENT STEP WE WISH TO CONSIDER WPT(I) (WHICH IS THEN
    ! JUST A CORRECTION TO WPT(I)) = 0.0 IF I CORRESPONDS TO NO RESOLVENT
    ! (SINCE THE VALUE OF SUCH WPT(I) IN SAVE SHOULD NOT CHANGE) SO WE OMIT
    ! THE CORRESPONDING TERMS IN THE SUMMATION ABOVE.
200 do
        do kk = 1 , nresl
            k = nresl - kk + 1
            imax = Ixrct(2*k-1)
            jmax = Ixrct(2*k)
            Wpt(imax) = Pmat1(n,jmax)
            do i = 1 , Ndm
                do l = 1 , k
                    if ( i==Ixrct(2*l-1) ) goto 250
                end do
                ! HERE ROW I CONTAINS NO EARLIER (OR PRESENT) RESOLVENTS.
                if ( itrct>0 ) then
                    if ( k<nresl ) then
                        ! HERE WE ARE DOING ITERATIVE REFINEMENT, K < NRESL, AND I /=
                        ! IXRCT(2*L-1) FOR L=1,...,K.  WE WILL USE THE TERM CORRESPONDING TO
                        ! WPT(I) IFF I = IXRCT(2*L-1) FOR SOME L = K+1,...,NRESL.
                        kp1 = k + 1
                        do l = kp1 , nresl
                            if ( i==Ixrct(2*l-1) ) goto 220
                        end do
                    end if
                    goto 250
                end if
    220         Wpt(imax) = Wpt(imax) - Pmat1(i,jmax)*Wpt(i)
    250     end do
            Wpt(imax) = Wpt(imax)/Pmat1(imax,jmax)
        end do
        ! END OF BACK SUBSTITUTION.

        ! IF ITRCT IS POSITIVE THEN WPT WILL CONTAIN ONLY AN ITERATIVE
        ! REFINEMENT CORRECTION IN THOSE POSITIONS CORRESPONDING TO RESOLVENTS
        ! AND WE ADD THIS TO SAVE TO GET THE TRUE WPT.
        if ( itrct>0 ) then
            do i = 1 , Ndm
                do l = 1 , nresl
                    if ( i==Ixrct(2*l-1) ) then
                        Wpt(i) = Wpt(i) + Save(i)
                        exit
                    end if
                end do
            end do
        end if

        ! NOW COMPUTE THE RESIDUAL AND PUT IT INTO PMAT1(NDM+1,.).
        do k = 1 , Ncor
            ! COMPUTE THE COLUMN INDEX KCOL IN PMAT CORRESPONDING TO COLUMN K IN
            ! PMAT1.
            kcol = Icor(k)
            Pmat1(n,k) = -Pmat(n,kcol)
            do i = 1 , Ndm
                Pmat1(n,k) = Pmat1(n,k) - Pmat(i,kcol)*Wpt(i)
            end do
        end do

        ! COMPUTE THE WORST ABSOLUTE VALUE OF THE RESIDUAL ELEMENTS.
        do k = 1 , Ncor
            aa = abs(Pmat1(n,k))
            if ( k>1 ) then
                if ( aa<=wrst ) cycle
            end if
            wrst = aa
        end do

        if ( itrct<=0 ) then
            wrsto = wrst
        else if ( wrst>wrsto ) then
            ! HERE ITRCT > 0 AND WRST > WRSTO, SO WE GO BACK TO THE PREVIOUS
            ! WPT AND RETURN.
            wrst = wrsto
            do i = 1 , Ndm
                Wpt(i) = Save(i)
            end do
            return
        end if

        if ( itrct>=itrlm ) return

        ! HERE ITRCT < ITRLM AND WE INCREMENT ITRCT AND SET UP FOR THE ITRCTTH
        ! ITERATIVE REFINEMENT STEP.
        itrct = itrct + 1
        ! COPY WPT INTO SAVE.
        do i = 1 , Ndm
            Save(i) = Wpt(i)
        end do

    end do

    end subroutine refwl
!********************************************************************************

!********************************************************************************
    end module conmax_module
!********************************************************************************
