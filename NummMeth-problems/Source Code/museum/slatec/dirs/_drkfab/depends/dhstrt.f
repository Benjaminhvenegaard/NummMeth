      SUBROUTINE DHSTRT (DF, NEQ, A, B, Y, YPRIME, ETOL, MORDER, SMALL,
     +   BIG, SPY, PV, YP, SF, RPAR, IPAR, H)
C
      INTEGER IPAR, J, K, LK, MORDER, NEQ
      DOUBLE PRECISION A, ABSDX, B, BIG, DA, DELF, DELY,
     1      DFDUB, DFDXB, DHVNRM,
     2      DX, DY, ETOL, FBND, H, PV, RELPER, RPAR, SF, SMALL, SPY,
     3      SRYDPB, TOLEXP, TOLMIN, TOLP, TOLSUM, Y, YDPB, YP, YPRIME
      DIMENSION Y(*),YPRIME(*),ETOL(*),SPY(*),PV(*),YP(*),
     1          SF(*),RPAR(*),IPAR(*)
      EXTERNAL DF
C
C     ..................................................................
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 160
C***FIRST EXECUTABLE STATEMENT  DHSTRT
         DX = B - A
         ABSDX = ABS(DX)
         RELPER = SMALL**0.375D0
C
C        ...............................................................
C
C             COMPUTE AN APPROXIMATE BOUND (DFDXB) ON THE PARTIAL
C             DERIVATIVE OF THE EQUATION WITH RESPECT TO THE
C             INDEPENDENT VARIABLE. PROTECT AGAINST AN OVERFLOW.
C             ALSO COMPUTE A BOUND (FBND) ON THE FIRST DERIVATIVE
C             LOCALLY.
C
         DA = SIGN(MAX(MIN(RELPER*ABS(A),ABSDX),
     1                    100.0D0*SMALL*ABS(A)),DX)
         IF (DA .EQ. 0.0D0) DA = RELPER*DX
         CALL DF(A+DA,Y,SF,RPAR,IPAR)
         DO 10 J = 1, NEQ
            YP(J) = SF(J) - YPRIME(J)
   10    CONTINUE
         DELF = DHVNRM(YP,NEQ)
         DFDXB = BIG
         IF (DELF .LT. BIG*ABS(DA)) DFDXB = DELF/ABS(DA)
         FBND = DHVNRM(SF,NEQ)
C
C        ...............................................................
C
C             COMPUTE AN ESTIMATE (DFDUB) OF THE LOCAL LIPSCHITZ
C             CONSTANT FOR THE SYSTEM OF DIFFERENTIAL EQUATIONS. THIS
C             ALSO REPRESENTS AN ESTIMATE OF THE NORM OF THE JACOBIAN
C             LOCALLY.  THREE ITERATIONS (TWO WHEN NEQ=1) ARE USED TO
C             ESTIMATE THE LIPSCHITZ CONSTANT BY NUMERICAL DIFFERENCES.
C             THE FIRST PERTURBATION VECTOR IS BASED ON THE INITIAL
C             DERIVATIVES AND DIRECTION OF INTEGRATION. THE SECOND
C             PERTURBATION VECTOR IS FORMED USING ANOTHER EVALUATION OF
C             THE DIFFERENTIAL EQUATION.  THE THIRD PERTURBATION VECTOR
C             IS FORMED USING PERTURBATIONS BASED ONLY ON THE INITIAL
C             VALUES. COMPONENTS THAT ARE ZERO ARE ALWAYS CHANGED TO
C             NON-ZERO VALUES (EXCEPT ON THE FIRST ITERATION). WHEN
C             INFORMATION IS AVAILABLE, CARE IS TAKEN TO ENSURE THAT
C             COMPONENTS OF THE PERTURBATION VECTOR HAVE SIGNS WHICH ARE
C             CONSISTENT WITH THE SLOPES OF LOCAL SOLUTION CURVES.
C             ALSO CHOOSE THE LARGEST BOUND (FBND) FOR THE FIRST
C             DERIVATIVE.
C
C                               PERTURBATION VECTOR SIZE IS HELD
C                               CONSTANT FOR ALL ITERATIONS. COMPUTE
C                               THIS CHANGE FROM THE
C                                       SIZE OF THE VECTOR OF INITIAL
C                                       VALUES.
         DELY = RELPER*DHVNRM(Y,NEQ)
         IF (DELY .EQ. 0.0D0) DELY = RELPER
         DELY = SIGN(DELY,DX)
         DELF = DHVNRM(YPRIME,NEQ)
         FBND = MAX(FBND,DELF)
         IF (DELF .EQ. 0.0D0) GO TO 30
C           USE INITIAL DERIVATIVES FOR FIRST PERTURBATION
            DO 20 J = 1, NEQ
               SPY(J) = YPRIME(J)
               YP(J) = YPRIME(J)
   20       CONTINUE
         GO TO 50
   30    CONTINUE
C           CANNOT HAVE A NULL PERTURBATION VECTOR
            DO 40 J = 1, NEQ
               SPY(J) = 0.0D0
               YP(J) = 1.0D0
   40       CONTINUE
            DELF = DHVNRM(YP,NEQ)
   50    CONTINUE
C
         DFDUB = 0.0D0
         LK = MIN(NEQ+1,3)
         DO 140 K = 1, LK
C           DEFINE PERTURBED VECTOR OF INITIAL VALUES
            DO 60 J = 1, NEQ
               PV(J) = Y(J) + DELY*(YP(J)/DELF)
   60       CONTINUE
            IF (K .EQ. 2) GO TO 80
C              EVALUATE DERIVATIVES ASSOCIATED WITH PERTURBED
C              VECTOR  AND  COMPUTE CORRESPONDING DIFFERENCES
               CALL DF(A,PV,YP,RPAR,IPAR)
               DO 70 J = 1, NEQ
                  PV(J) = YP(J) - YPRIME(J)
   70          CONTINUE
            GO TO 100
   80       CONTINUE
C              USE A SHIFTED VALUE OF THE INDEPENDENT VARIABLE
C                                    IN COMPUTING ONE ESTIMATE
               CALL DF(A+DA,PV,YP,RPAR,IPAR)
               DO 90 J = 1, NEQ
                  PV(J) = YP(J) - SF(J)
   90          CONTINUE
  100       CONTINUE
C           CHOOSE LARGEST BOUNDS ON THE FIRST DERIVATIVE
C                          AND A LOCAL LIPSCHITZ CONSTANT
            FBND = MAX(FBND,DHVNRM(YP,NEQ))
            DELF = DHVNRM(PV,NEQ)
C        ...EXIT
            IF (DELF .GE. BIG*ABS(DELY)) GO TO 150
            DFDUB = MAX(DFDUB,DELF/ABS(DELY))
C     ......EXIT
            IF (K .EQ. LK) GO TO 160
C           CHOOSE NEXT PERTURBATION VECTOR
            IF (DELF .EQ. 0.0D0) DELF = 1.0D0
            DO 130 J = 1, NEQ
               IF (K .EQ. 2) GO TO 110
                  DY = ABS(PV(J))
                  IF (DY .EQ. 0.0D0) DY = DELF
               GO TO 120
  110          CONTINUE
                  DY = Y(J)
                  IF (DY .EQ. 0.0D0) DY = DELY/RELPER
  120          CONTINUE
               IF (SPY(J) .EQ. 0.0D0) SPY(J) = YP(J)
               IF (SPY(J) .NE. 0.0D0) DY = SIGN(DY,SPY(J))
               YP(J) = DY
  130       CONTINUE
            DELF = DHVNRM(YP,NEQ)
  140    CONTINUE
  150    CONTINUE
C
C        PROTECT AGAINST AN OVERFLOW
         DFDUB = BIG
  160 CONTINUE
C
C     ..................................................................
C
C          COMPUTE A BOUND (YDPB) ON THE NORM OF THE SECOND DERIVATIVE
C
      YDPB = DFDXB + DFDUB*FBND
C
C     ..................................................................
C
C          DEFINE THE TOLERANCE PARAMETER UPON WHICH THE STARTING STEP
C          SIZE IS TO BE BASED.  A VALUE IN THE MIDDLE OF THE ERROR
C          TOLERANCE RANGE IS SELECTED.
C
      TOLMIN = BIG
      TOLSUM = 0.0D0
      DO 170 K = 1, NEQ
         TOLEXP = LOG10(ETOL(K))
         TOLMIN = MIN(TOLMIN,TOLEXP)
         TOLSUM = TOLSUM + TOLEXP
  170 CONTINUE
      TOLP = 10.0D0**(0.5D0*(TOLSUM/NEQ + TOLMIN)/(MORDER+1))
C
C     ..................................................................
C
C          COMPUTE A STARTING STEP SIZE BASED ON THE ABOVE FIRST AND
C          SECOND DERIVATIVE INFORMATION
C
C                            RESTRICT THE STEP LENGTH TO BE NOT BIGGER
C                            THAN ABS(B-A).   (UNLESS  B  IS TOO CLOSE
C                            TO  A)
      H = ABSDX
C
      IF (YDPB .NE. 0.0D0 .OR. FBND .NE. 0.0D0) GO TO 180
C
C        BOTH FIRST DERIVATIVE TERM (FBND) AND SECOND
C                     DERIVATIVE TERM (YDPB) ARE ZERO
         IF (TOLP .LT. 1.0D0) H = ABSDX*TOLP
      GO TO 200
  180 CONTINUE
C
      IF (YDPB .NE. 0.0D0) GO TO 190
C
C        ONLY SECOND DERIVATIVE TERM (YDPB) IS ZERO
         IF (TOLP .LT. FBND*ABSDX) H = TOLP/FBND
      GO TO 200
  190 CONTINUE
C
C        SECOND DERIVATIVE TERM (YDPB) IS NON-ZERO
         SRYDPB = SQRT(0.5D0*YDPB)
         IF (TOLP .LT. SRYDPB*ABSDX) H = TOLP/SRYDPB
  200 CONTINUE
C
C     FURTHER RESTRICT THE STEP LENGTH TO BE NOT
C                               BIGGER THAN  1/DFDUB
      IF (H*DFDUB .GT. 1.0D0) H = 1.0D0/DFDUB
C
C     FINALLY, RESTRICT THE STEP LENGTH TO BE NOT
C     SMALLER THAN  100*SMALL*ABS(A).  HOWEVER, IF
C     A=0. AND THE COMPUTED H UNDERFLOWED TO ZERO,
C     THE ALGORITHM RETURNS  SMALL*ABS(B)  FOR THE
C                                     STEP LENGTH.
      H = MAX(H,100.0D0*SMALL*ABS(A))
      IF (H .EQ. 0.0D0) H = SMALL*ABS(B)
C
C     NOW SET DIRECTION OF INTEGRATION
      H = SIGN(H,DX)
C
      RETURN
      END
