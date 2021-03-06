      FUNCTION ALNGAM (X)
      LOGICAL FIRST
      EXTERNAL GAMMA
      SAVE SQ2PIL, SQPI2L, PI, XMAX, DXREL, FIRST
      DATA SQ2PIL / 0.9189385332 0467274E0/
      DATA SQPI2L / 0.2257913526 4472743E0/
      DATA PI     / 3.1415926535 8979324E0/
      DATA FIRST  /.TRUE./
C***FIRST EXECUTABLE STATEMENT  ALNGAM
      IF (FIRST) THEN
         XMAX = R1MACH(2)/LOG(R1MACH(2))
         DXREL = SQRT (R1MACH(4))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.10.0) GO TO 20
C
C LOG (ABS (GAMMA(X))) FOR  ABS(X) .LE. 10.0
C
      ALNGAM = LOG (ABS (GAMMA(X)))
      RETURN
C
C LOG (ABS (GAMMA(X))) FOR ABS(X) .GT. 10.0
C
 20   IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'ALNGAM',
     +   'ABS(X) SO BIG ALNGAM OVERFLOWS', 2, 2)
C
      IF (X.GT.0.) ALNGAM = SQ2PIL + (X-0.5)*LOG(X) - X + R9LGMC(Y)
      IF (X.GT.0.) RETURN
C
      SINPIY = ABS (SIN(PI*Y))
      IF (SINPIY .EQ. 0.) CALL XERMSG ('SLATEC', 'ALNGAM',
     +   'X IS A NEGATIVE INTEGER', 3, 2)
C
      IF (ABS((X-AINT(X-0.5))/X) .LT. DXREL) CALL XERMSG ('SLATEC',
     +   'ALNGAM', 'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR ' //
     +   'NEGATIVE INTEGER', 1, 1)
C
      ALNGAM = SQPI2L + (X-0.5)*LOG(Y) - X - LOG(SINPIY) - R9LGMC(Y)
      RETURN
C
      END
