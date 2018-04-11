      DOUBLE PRECISION FUNCTION DBVALU (T, A, N, K, IDERIV, X, INBV,
     +   WORK)
C
      INTEGER I,IDERIV,IDERP1,IHI,IHMKMJ,ILO,IMK,IMKPJ, INBV, IPJ,
     1 IP1, IP1MJ, J, JJ, J1, J2, K, KMIDER, KMJ, KM1, KPK, MFLAG, N
      DOUBLE PRECISION A, FKMJ, T, WORK, X
      DIMENSION T(*), A(*), WORK(*)
C***FIRST EXECUTABLE STATEMENT  DBVALU
      DBVALU = 0.0D0
      IF(K.LT.1) GO TO 102
      IF(N.LT.K) GO TO 101
      IF(IDERIV.LT.0 .OR. IDERIV.GE.K) GO TO 110
      KMIDER = K - IDERIV
C
C *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
C     (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
      KM1 = K - 1
      CALL DINTRV(T, N+1, X, INBV, I, MFLAG)
      IF (X.LT.T(K)) GO TO 120
      IF (MFLAG.EQ.0) GO TO 20
      IF (X.GT.T(I)) GO TO 130
   10 IF (I.EQ.K) GO TO 140
      I = I - 1
      IF (X.EQ.T(I)) GO TO 10
C
C *** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES
C     WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K
C
   20 IMK = I - K
      DO 30 J=1,K
        IMKPJ = IMK + J
        WORK(J) = A(IMKPJ)
   30 CONTINUE
      IF (IDERIV.EQ.0) GO TO 60
      DO 50 J=1,IDERIV
        KMJ = K - J
        FKMJ = KMJ
        DO 40 JJ=1,KMJ
          IHI = I + JJ
          IHMKMJ = IHI - KMJ
          WORK(JJ) = (WORK(JJ+1)-WORK(JJ))/(T(IHI)-T(IHMKMJ))*FKMJ
   40   CONTINUE
   50 CONTINUE
C
C *** COMPUTE VALUE AT *X* IN (T(I),(T(I+1)) OF IDERIV-TH DERIVATIVE,
C     GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV).
   60 IF (IDERIV.EQ.KM1) GO TO 100
      IP1 = I + 1
      KPK = K + K
      J1 = K + 1
      J2 = KPK + 1
      DO 70 J=1,KMIDER
        IPJ = I + J
        WORK(J1) = T(IPJ) - X
        IP1MJ = IP1 - J
        WORK(J2) = X - T(IP1MJ)
        J1 = J1 + 1
        J2 = J2 + 1
   70 CONTINUE
      IDERP1 = IDERIV + 1
      DO 90 J=IDERP1,KM1
        KMJ = K - J
        ILO = KMJ
        DO 80 JJ=1,KMJ
          WORK(JJ) = (WORK(JJ+1)*WORK(KPK+ILO)+WORK(JJ)
     1              *WORK(K+JJ))/(WORK(KPK+ILO)+WORK(K+JJ))
          ILO = ILO - 1
   80   CONTINUE
   90 CONTINUE
  100 DBVALU = WORK(1)
      RETURN
C
C
  101 CONTINUE
      CALL XERMSG ('SLATEC', 'DBVALU', 'N DOES NOT SATISFY N.GE.K', 2,
     +   1)
      RETURN
  102 CONTINUE
      CALL XERMSG ('SLATEC', 'DBVALU', 'K DOES NOT SATISFY K.GE.1', 2,
     +   1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('SLATEC', 'DBVALU',
     +   'IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K', 2, 1)
      RETURN
  120 CONTINUE
      CALL XERMSG ('SLATEC', 'DBVALU',
     +   'X IS N0T GREATER THAN OR EQUAL TO T(K)', 2, 1)
      RETURN
  130 CONTINUE
      CALL XERMSG ('SLATEC', 'DBVALU',
     +   'X IS NOT LESS THAN OR EQUAL TO T(N+1)', 2, 1)
      RETURN
  140 CONTINUE
      CALL XERMSG ('SLATEC', 'DBVALU',
     +   'A LEFT LIMITING VALUE CANNOT BE OBTAINED AT T(K)', 2, 1)
      RETURN
      END
