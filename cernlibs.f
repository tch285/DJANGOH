C-----------------------------------------------------------------------
*
* $Id: abend.F,v 1.1.1.1 1996/02/15 17:50:37 mclareni Exp $
*
* $Log: abend.F,v $
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
*
      SUBROUTINE ABEND
C
C CERN PROGLIB# Z035    ABEND           .VERSION KERNFOR  4.31  911111
C ORIG.  8/02/88  JZ
C

      STOP  7
      END
C-----------------------------------------------------------------------
*
* $Id: timex.F,v 1.1.1.1 1996/02/15 17:50:38 mclareni Exp $
*
* $Log: timex.F,v $
* Revision 1.1.1.1  1996/02/15 17:50:38  mclareni
* Kernlib
*
*
      SUBROUTINE TIMEX (T)
C
C CERN PROGLIB# Z007    TIMEX   DUMMY   .VERSION KERNFOR  4.05  821202
C
C-    DUMMY FOR NON-ESSENTIAL ROUTINE STILL MISSING ON YOUR MACHINE
      CALL CPU_TIME(T)
CCC      T = 9.
      RETURN
      END
C-----------------------------------------------------------------------
*
* $Id: divdif.F,v 1.1.1.1 1996/02/15 17:48:36 mclareni Exp $
*
* $Log: divdif.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:36  mclareni
* Kernlib
*
*
      FUNCTION DIVDIF(F,A,NN,X,MM)
      DIMENSION A(NN),F(NN),T(20),D(20)
      LOGICAL EXTRA
      LOGICAL MFLAG,RFLAG
      DATA MMAX/10/
C
C  TABULAR INTERPOLATION USING SYMMETRICALLY PLACED ARGUMENT POINTS.
C
C  START.  FIND SUBSCRIPT IX OF X IN ARRAY A.
      IF( (NN.LT.2) .OR. (MM.LT.1) ) GO TO 20
      N=NN
      M=MIN0(MM,MMAX,N-1)
      MPLUS=M+1
      IX=0
      IY=N+1
      IF(A(1).GT.A(N)) GO TO 4
C     (SEARCH INCREASING ARGUMENTS.)
    1    MID=(IX+IY)/2
         IF(X.GE.A(MID)) GO TO 2
            IY=MID
            GO TO 3
C        (IF TRUE.)
    2       IX=MID
    3    IF(IY-IX.GT.1) GO TO 1
         GO TO 7
C     (SEARCH DECREASING ARGUMENTS.)
    4    MID=(IX+IY)/2
         IF(X.LE.A(MID)) GO TO 5
            IY=MID
            GO TO 6
C        (IF TRUE.)
    5       IX=MID
    6    IF(IY-IX.GT.1) GO TO 4
C
C  COPY REORDERED INTERPOLATION POINTS INTO (T(I),D(I)), SETTING
C  *EXTRA* TO TRUE IF M+2 POINTS TO BE USED.
    7 NPTS=M+2-MOD(M,2)
      IP=0
      L=0
      GO TO 9
    8    L=-L
         IF(L.GE.0) L=L+1
    9    ISUB=IX+L
         IF((1.LE.ISUB).AND.(ISUB.LE.N)) GO TO 10
C        (SKIP POINT.)
            NPTS=MPLUS
            GO TO 11
C        (INSERT POINT.)
   10       IP=IP+1
            T(IP)=A(ISUB)
            D(IP)=F(ISUB)
   11    IF(IP.LT.NPTS) GO TO 8
      EXTRA=NPTS.NE.MPLUS
C
C  REPLACE D BY THE LEADING DIAGONAL OF A DIVIDED-DIFFERENCE TABLE, SUP-
C  PLEMENTED BY AN EXTRA LINE IF *EXTRA* IS TRUE.
      DO 14 L=1,M
         IF(.NOT.EXTRA) GO TO 12
            ISUB=MPLUS-L
            D(M+2)=(D(M+2)-D(M))/(T(M+2)-T(ISUB))
   12    I=MPLUS
         DO 13 J=L,M
            ISUB=I-L
            D(I)=(D(I)-D(I-1))/(T(I)-T(ISUB))
            I=I-1
   13    CONTINUE
   14 CONTINUE
C
C  EVALUATE THE NEWTON INTERPOLATION FORMULA AT X, AVERAGING TWO VALUES
C  OF LAST DIFFERENCE IF *EXTRA* IS TRUE.
      SUM=D(MPLUS)
      IF(EXTRA) SUM=0.5*(SUM+D(M+2))
      J=M
      DO 15 L=1,M
         SUM=D(J)+(X-T(J))*SUM
         J=J-1
   15 CONTINUE
      DIVDIF=SUM
      RETURN
C
   20 print *,' DIVDIF - CALL OUT OF LIMITS NM,MM ',Nm, MM
      DIVDIF=0
      IF(MFLAG) THEN
         IF(LGFILE.EQ.0) THEN
            IF(MM.LT.1) WRITE(*,101) MM
            IF(NN.LT.2) WRITE(*,102) NN
         ELSE
            IF(MM.LT.1) WRITE(LGFILE,101) MM
            IF(NN.LT.2) WRITE(LGFILE,102) NN
         ENDIF
      ENDIF
      IF(.NOT.RFLAG) CALL ABEND
      RETURN
  101 FORMAT( 7X, 'FUNCTION DIVDIF ... M =',I6,' IS LESS THAN 1')
  102 FORMAT( 7X, 'FUNCTION DIVDIF ... N =',I6,' IS LESS THAN 2')
      END
C-----------------------------------------------------------------------
*
* $Id: fint.F,v 1.1.1.1 1996/02/15 17:48:36 mclareni Exp $
*
* $Log: fint.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:36  mclareni
* Kernlib
*
*
          FUNCTION FINT(NARG,ARG,NENT,ENT,TABLE)
C
C   INTERPOLATION ROUTINE. AUTHOR C. LETERTRE.
C   MODIFIED BY B. SCHORR, 1.07.1982.
C
          INTEGER   NENT(9)
          REAL      ARG(9),   ENT(9),   TABLE(9)
          INTEGER   INDEX(32)
          REAL      WEIGHT(32)
          LOGICAL   MFLAG,    RFLAG
          FINT  =  0.
          IF(NARG .LT. 1  .OR.  NARG .GT. 5)  GOTO 300
          LMAX      =  0
          ISTEP     =  1
          KNOTS     =  1
          INDEX(1)  =  1
          WEIGHT(1) =  1.
          DO 100    N  =  1, NARG
             X     =  ARG(N)
             NDIM  =  NENT(N)
             LOCA  =  LMAX
             LMIN  =  LMAX + 1
             LMAX  =  LMAX + NDIM
             IF(NDIM .GT. 2)  GOTO 10
             IF(NDIM .EQ. 1)  GOTO 100
             H  =  X - ENT(LMIN)
             IF(H .EQ. 0.)  GOTO 90
             ISHIFT  =  ISTEP
             IF(X-ENT(LMIN+1) .EQ. 0.)  GOTO 21
             ISHIFT  =  0
             ETA     =  H / (ENT(LMIN+1) - ENT(LMIN))
             GOTO 30
  10         LOCB  =  LMAX + 1
  11         LOCC  =  (LOCA+LOCB) / 2
             IF(X-ENT(LOCC))  12, 20, 13
  12         LOCB  =  LOCC
             GOTO 14
  13         LOCA  =  LOCC
  14         IF(LOCB-LOCA .GT. 1)  GOTO 11
             LOCA    =  MIN0( MAX0(LOCA,LMIN), LMAX-1 )
             ISHIFT  =  (LOCA - LMIN) * ISTEP
             ETA     =  (X - ENT(LOCA)) / (ENT(LOCA+1) - ENT(LOCA))
             GOTO 30
  20         ISHIFT  =  (LOCC - LMIN) * ISTEP
  21         DO 22  K  =  1, KNOTS
                INDEX(K)  =  INDEX(K) + ISHIFT
  22            CONTINUE
             GOTO 90
  30         DO 31  K  =  1, KNOTS
                INDEX(K)         =  INDEX(K) + ISHIFT
                INDEX(K+KNOTS)   =  INDEX(K) + ISTEP
                WEIGHT(K+KNOTS)  =  WEIGHT(K) * ETA
                WEIGHT(K)        =  WEIGHT(K) - WEIGHT(K+KNOTS)
  31            CONTINUE
             KNOTS  =  2*KNOTS
  90         ISTEP  =  ISTEP * NDIM
 100         CONTINUE
          DO 200    K  =  1, KNOTS
             I  =  INDEX(K)
             FINT  =  FINT + WEIGHT(K) * TABLE(I)
 200         CONTINUE
          RETURN
 300      print *,' FINT - CALL OUT OF LIMITS NARG ',NARG 
          IF(MFLAG) THEN
             IF(LGFILE .EQ. 0) THEN
                WRITE(*,1000) NARG
             ELSE
                WRITE(LGFILE,1000) NARG
             ENDIF
          ENDIF
          IF(.NOT. RFLAG) CALL ABEND
          RETURN
1000      FORMAT( 7X, 24HFUNCTION FINT ... NARG =,I6,
     +              17H NOT WITHIN RANGE)
          END
C-----------------------------------------------------------------------*
* $Id: rnormx.F,v 1.1.1.1 1996/04/01 15:02:55 mclareni Exp $
*
* $Log: rnormx.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:55  mclareni
* Mathlib gen
*
*
      SUBROUTINE RNORMX(DEVIAS,NDEV,ROUTIN)
C        Generator of a vector of independent Gaussian-distributed 
C        (pseudo-)random numbers, of mean zero and variance one,
C        making use of a uniform pseudo-random generator (RANMAR).
C        The algorithm for converting uniform numbers to Gaussian
C        is that of "Ratio of Uniforms with Quadratic Bounds."  The
C        method is in principle exact (apart from rounding errors),
C        and is based on the variant published by Joseph Leva in
C        ACM TOMS vol. 18(1992), page 449 for the method and 454 for
C        the Fortran algorithm (ACM No. 712).
C        It requires at least 2 and on average 2.74 uniform deviates
C        per Gaussian (normal) deviate.
C   WARNING -- The uniform generator should not produce exact zeroes,
C   since the pair (0.0, 0.5) provokes a floating point exception.
      SAVE  S, T, A, B, R1, R2
      DIMENSION U(2), DEVIAS(*)
      EXTERNAL ROUTIN
      DATA  S, T, A, B / 0.449871, -0.386595, 0.19600, 0.25472/
      DATA  R1, R2/ 0.27597, 0.27846/
C         generate pair of uniform deviates
      DO 200 IDEV = 1, NDEV
   50 CALL ROUTIN(U,2)
      V = 1.7156 * (U(2) - 0.5)
      X = U(1) - S
      Y = ABS(V) - T
      Q = X**2 + Y*(A*Y - B*X)
C           accept P if inside inner ellipse
      IF (Q .LT. R1)  GO TO 100
C           reject P if outside outer ellipse
      IF (Q .GT. R2)  GO TO 50
C           reject P if outside acceptance region
      IF (V**2 .GT. -4.0 *ALOG(U(1)) *U(1)**2)  GO TO 50
C           ratio of P's coordinates is normal deviate
  100 DEVIAT = V/U(1)
  200 DEVIAS(IDEV) = DEVIAT
      RETURN
      END
C-----------------------------------------------------------------------
*
* $Id: ranlux.F,v 1.2 1997/09/22 13:45:47 mclareni Exp $
*
* $Log: ranlux.F,v $
* Revision 1.2  1997/09/22 13:45:47  mclareni
* Correct error in initializing RANLUX by using RLUXIN with the output of
* RLUXUT from a previous run.
*
* Revision 1.1.1.1  1996/04/01 15:02:55  mclareni
* Mathlib gen
*
*
      SUBROUTINE RANLUX(RVEC,LENV)
C         Subtract-and-borrow random number generator proposed by
C         Marsaglia and Zaman, implemented by F. James with the name
C         RCARRY in 1991, and later improved by Martin Luescher
C         in 1993 to produce "Luxury Pseudorandom Numbers".
C     Fortran 77 coded by F. James, 1993
C
C   LUXURY LEVELS.
C   ------ ------      The available luxury levels are:
C
C  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
C           and Zaman, very long period, but fails many tests.
C  level 1  (p=48): considerable improvement in quality over level 0,
C           now passes the gap test, but still fails spectral test.
C  level 2  (p=97): passes all known tests, but theoretically still
C           defective.
C  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
C           correlations have very small chance of being observed.
C  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
C
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RANLUX:                                  ++
C!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero (not included) and one (also not incl.). ++
C!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
C!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
C!!!               which is integer between zero and MAXLEV, or if   ++
C!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
C!!!               should be set to zero unless restarting at a break++ 
C!!!               point given by output of RLUXAT (see RLUXAT).     ++
C!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
C!!!               which can be used to restart the RANLUX generator ++
C!!!               at the current point by calling RLUXGO.  K1 and K2++
C!!!               specify how many numbers were generated since the ++
C!!!               initialization with LUX and INT.  The restarting  ++
C!!!               skips over  K1+K2*E9   numbers, so it can be long.++
C!!!   A more efficient but less convenient way of restarting is by: ++
C!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
C!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
C!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
C!!!                 32-bit integer seeds, to be used for restarting ++
C!!!      ISVEC must be dimensioned 25 in the calling program        ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION RVEC(LENV)
      DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
      PARAMETER (MAXLEV=4, LXDFLT=3)
      DIMENSION NDSKIP(0:MAXLEV)
      DIMENSION NEXT(24)
      PARAMETER (TWOP12=4096., IGIGA=1000000000,JSDFLT=314159265)
      PARAMETER (ITWO24=2**24, ICONS=2147483563)
      SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24, TWOM12, LUXLEV
      SAVE NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED
      INTEGER LUXLEV
      LOGICAL NOTYET
      DATA NOTYET, LUXLEV, IN24, KOUNT, MKOUNT /.TRUE., LXDFLT, 0,0,0/
      DATA I24,J24,CARRY/24,10,0./
      integer*8 initseed,idum1,idum2
      common/lseed/initseed,idum1,idum2
C                               default
C  Luxury Level   0     1     2   *3*    4
      DATA NDSKIP/0,   24,   73,  199,  365 /
Corresponds to p=24    48    97   223   389
C     time factor 1     2     3     6    10   on slow workstation
C                 1    1.5    2     3     5   on fast mainframe
C
C  NOTYET is .TRUE. if no initialization has been performed yet.
C              Default Initialization by Multiplicative Congruential
      IF (NOTYET) THEN
         NOTYET = .FALSE.
CAB         JSEED = JSDFLT  
         JSEED = INITSEED  
         INSEED = JSEED
         WRITE(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED
         LUXLEV = LXDFLT
         NSKIP = NDSKIP(LUXLEV)
         LP = NSKIP + 24
         IN24 = 0
         KOUNT = 0
         MKOUNT = 0
         WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
     +        LUXLEV,'      p =',LP
            TWOM24 = 1.
         DO 25 I= 1, 24
            TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
   25    CONTINUE
         TWOM12 = TWOM24 * 4096.
         DO 50 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
   50    CONTINUE
         NEXT(1) = 24
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
      ENDIF
C
C          The Generator proper: "Subtract-with-borrow",
C          as proposed by Marsaglia and Zaman,
C          Florida State University, March, 1989
C
      DO 100 IVEC= 1, LENV
      UNI = SEEDS(J24) - SEEDS(I24) - CARRY 
      IF (UNI .LT. 0.)  THEN
         UNI = UNI + 1.0
         CARRY = TWOM24
      ELSE
         CARRY = 0.
      ENDIF
      SEEDS(I24) = UNI
      I24 = NEXT(I24)
      J24 = NEXT(J24)
      RVEC(IVEC) = UNI
C  small numbers (with less than 12 "significant" bits) are "padded".
      IF (UNI .LT. TWOM12)  THEN
         RVEC(IVEC) = RVEC(IVEC) + TWOM24*SEEDS(J24)
C        and zero is forbidden in case someone takes a logarithm
         IF (RVEC(IVEC) .EQ. 0.)  RVEC(IVEC) = TWOM24*TWOM24
      ENDIF
C        Skipping to luxury.  As proposed by Martin Luscher.
      IN24 = IN24 + 1
      IF (IN24 .EQ. 24)  THEN
         IN24 = 0
         KOUNT = KOUNT + NSKIP
         DO 90 ISK= 1, NSKIP
         UNI = SEEDS(J24) - SEEDS(I24) - CARRY
         IF (UNI .LT. 0.)  THEN
            UNI = UNI + 1.0
            CARRY = TWOM24
         ELSE
            CARRY = 0.
         ENDIF
         SEEDS(I24) = UNI
         I24 = NEXT(I24)
         J24 = NEXT(J24)
   90    CONTINUE
      ENDIF
  100 CONTINUE
      KOUNT = KOUNT + LENV
      IF (KOUNT .GE. IGIGA)  THEN
         MKOUNT = MKOUNT + 1
         KOUNT = KOUNT - IGIGA
      ENDIF
      RETURN
C
C           Entry to input and float integer seeds from previous run
      ENTRY RLUXIN(ISDEXT)
         NOTYET = .FALSE.
         TWOM24 = 1.
         DO 195 I= 1, 24
         NEXT(I) = I-1
  195    TWOM24 = TWOM24 * 0.5
         NEXT(1) = 24
         TWOM12 = TWOM24 * 4096.
      WRITE(6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
      WRITE(6,'(5X,5I12)') ISDEXT
      DO 200 I= 1, 24
      SEEDS(I) = REAL(ISDEXT(I))*TWOM24
  200 CONTINUE
      CARRY = 0.
      IF (ISDEXT(25) .LT. 0)  CARRY = TWOM24
      ISD = IABS(ISDEXT(25))
      I24 = MOD(ISD,100)
      ISD = ISD/100
      J24 = MOD(ISD,100)
      ISD = ISD/100
      IN24 = MOD(ISD,100)
      ISD = ISD/100
      LUXLEV = ISD
        IF (LUXLEV .LE. MAXLEV) THEN
          NSKIP = NDSKIP(LUXLEV)
          WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',
     +                         LUXLEV
        ELSE  IF (LUXLEV .GE. 24) THEN
          NSKIP = LUXLEV - 24
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV
        ELSE
          NSKIP = NDSKIP(MAXLEV)
          WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV
          LUXLEV = MAXLEV
        ENDIF
      INSEED = -1
      RETURN
C
C                    Entry to ouput seeds as integers
      ENTRY RLUXUT(ISDEXT)
      DO 300 I= 1, 24
         ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
  300 CONTINUE
      ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV
      IF (CARRY .GT. 0.)  ISDEXT(25) = -ISDEXT(25)
      RETURN
C
C                    Entry to output the "convenient" restart point
      ENTRY RLUXAT(LOUT,INOUT,K1,K2)
      LOUT = LUXLEV
      INOUT = INSEED
      K1 = KOUNT
      K2 = MKOUNT
      RETURN
C
C                    Entry to initialize from one or three integers
      ENTRY RLUXGO(LUX,INS,K1,K2)
         IF (LUX .LT. 0) THEN
            LUXLEV = LXDFLT
         ELSE IF (LUX .LE. MAXLEV) THEN
            LUXLEV = LUX
         ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN
            LUXLEV = MAXLEV
            WRITE (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',LUX
         ELSE
            LUXLEV = LUX
            DO 310 ILX= 0, MAXLEV
              IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX
  310       CONTINUE
         ENDIF
      IF (LUXLEV .LE. MAXLEV)  THEN
         NSKIP = NDSKIP(LUXLEV)
         WRITE(6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',
     +        LUXLEV,'     P=', NSKIP+24
      ELSE
          NSKIP = LUXLEV - 24
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV
      ENDIF
      IN24 = 0
      IF (INS .LT. 0)  WRITE (6,'(A)')   
     +   ' Illegal initialization by RLUXGO, negative input seed'
      IF (INS .GT. 0)  THEN
        JSEED = INS
        WRITE(6,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
     +      JSEED, K1,K2
      ELSE
        JSEED = JSDFLT
        WRITE(6,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
      ENDIF
      INSEED = JSEED
      NOTYET = .FALSE.
      TWOM24 = 1.
         DO 325 I= 1, 24
           TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
  325    CONTINUE
      TWOM12 = TWOM24 * 4096.
         DO 350 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
  350    CONTINUE
      NEXT(1) = 24
      I24 = 24
      J24 = 10
      CARRY = 0.
      IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
C        If restarting at a break point, skip K1 + IGIGA*K2
C        Note that this is the number of numbers delivered to
C        the user PLUS the number skipped (if luxury .GT. 0).
      KOUNT = K1
      MKOUNT = K2
      IF (K1+K2 .NE. 0)  THEN
        DO 500 IOUTER= 1, K2+1
          INNER = IGIGA
          IF (IOUTER .EQ. K2+1)  INNER = K1
          DO 450 ISK= 1, INNER
            UNI = SEEDS(J24) - SEEDS(I24) - CARRY 
            IF (UNI .LT. 0.)  THEN
               UNI = UNI + 1.0
               CARRY = TWOM24
            ELSE
               CARRY = 0.
            ENDIF
            SEEDS(I24) = UNI
            I24 = NEXT(I24)
            J24 = NEXT(J24)
  450     CONTINUE
  500   CONTINUE
C         Get the right value of IN24 by direct calculation
        IN24 = MOD(KOUNT, NSKIP+24)
        IF (MKOUNT .GT. 0)  THEN
           IZIP = MOD(IGIGA, NSKIP+24)
           IZIP2 = MKOUNT*IZIP + IN24
           IN24 = MOD(IZIP2, NSKIP+24)
        ENDIF
C       Now IN24 had better be between zero and 23 inclusive
        IF (IN24 .GT. 23) THEN
           WRITE (6,'(A/A,3I11,A,I5)')  
     +    '  Error in RESTARTING with RLUXGO:','  The values', INS,
     +     K1, K2, ' cannot occur at luxury level', LUXLEV
           IN24 = 0
        ENDIF
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
*
* $Id: cross.F,v 1.1.1.1 1996/02/15 17:49:51 mclareni Exp $
*
* $Log: cross.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:51  mclareni
* Kernlib
*
*
      SUBROUTINE CROSS(X,Y,Z)
C
C CERN PROGLIB# F117    CROSS           .VERSION KERNFOR  2.02  730125
C ORIG. 01/01/64
C
      COMMON /SLATE/Z1,Z2,XQX(38)
      DIMENSION X(3),Y(3),Z(3)
C
C
      Z1=X(2)*Y(3)-X(3)*Y(2)
      Z2=X(3)*Y(1)-X(1)*Y(3)
      Z(3)=X(1)*Y(2)-X(2)*Y(1)
      Z(1)=Z1
      Z(2)=Z2
      RETURN
      END
*
* $Id: vdot.F,v 1.1.1.1 1996/02/15 17:50:16 mclareni Exp $
*
* $Log: vdot.F,v $
* Revision 1.1.1.1  1996/02/15 17:50:16  mclareni
* Kernlib
*
*
      FUNCTION VDOT (X,Y,N)
C
C CERN PROGLIB# F121    VDOT            .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION X(*),Y(*)
C
      XX= 0.
      IF (N.LE.0)  GO TO 100
      DO 9 I= 1,N
    9 XX= XX + X(I)*Y(I)
C
  100 VDOT= XX
      RETURN
      END
*
* $Id: vdotn.F,v 1.1.1.1 1996/02/15 17:50:16 mclareni Exp $
*
* $Log: vdotn.F,v $
* Revision 1.1.1.1  1996/02/15 17:50:16  mclareni
* Kernlib
*
*
      FUNCTION VDOTN (X,Y,N)
C
C CERN PROGLIB# F121    VDOTN           .VERSION KERNFOR  2.09  751101
C ORIG. 01/07/71
C
      DIMENSION X(*), Y(*)
C
      IF (N.LE.0) GO TO 100
      XX= 0.
      XY= 0.
      YY= 0.
C
      DO 9 J=1,N
      XX = XX + X(J)*X(J)
      XY = XY + X(J)*Y(J)
    9 YY = YY + Y(J)*Y(J)
C
      VDOTN= XY / SQRT(XX*YY)
      IF (ABS (VDOTN).LT.1.) RETURN
      VDOTN= SIGN (1.,VDOTN)
      RETURN
  100 VDOTN = 0.
      RETURN
      END
