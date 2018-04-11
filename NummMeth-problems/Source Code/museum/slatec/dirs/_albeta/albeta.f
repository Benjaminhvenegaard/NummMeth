      FUNCTION ALBETA (A, B)
      EXTERNAL GAMMA
      SAVE SQ2PIL
      DATA SQ2PIL / 0.9189385332 0467274 E0 /
C***FIRST EXECUTABLE STATEMENT  ALBETA
      P = MIN (A, B)
      Q = MAX (A, B)
C
      IF (P .LE. 0.0) CALL XERMSG ('SLATEC', 'ALBETA',
     +   'BOTH ARGUMENTS MUST BE GT ZERO', 1, 2)
      IF (P.GE.10.0) GO TO 30
      IF (Q.GE.10.0) GO TO 20
C
C P AND Q ARE SMALL.
C
      ALBETA = LOG(GAMMA(P) * (GAMMA(Q)/GAMMA(P+Q)) )
      RETURN
C
C P IS SMALL, BUT Q IS BIG.
C
 20   CORR = R9LGMC(Q) - R9LGMC(P+Q)
      ALBETA = ALNGAM(P) + CORR + P - P*LOG(P+Q) +
     1  (Q-0.5)*ALNREL(-P/(P+Q))
      RETURN
C
C P AND Q ARE BIG.
C
 30   CORR = R9LGMC(P) + R9LGMC(Q) - R9LGMC(P+Q)
      ALBETA = -0.5*LOG(Q) + SQ2PIL + CORR + (P-0.5)*LOG(P/(P+Q))
     1  + Q*ALNREL(-P/(P+Q))
      RETURN
C
      END