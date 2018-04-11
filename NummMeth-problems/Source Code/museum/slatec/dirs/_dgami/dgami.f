      DOUBLE PRECISION FUNCTION DGAMI (A, X)
      DOUBLE PRECISION A, X, FACTOR, DLNGAM, DGAMIT
C***FIRST EXECUTABLE STATEMENT  DGAMI
      IF (A .LE. 0.D0) CALL XERMSG ('SLATEC', 'DGAMI',
     +   'A MUST BE GT ZERO', 1, 2)
      IF (X .LT. 0.D0) CALL XERMSG ('SLATEC', 'DGAMI',
     +   'X MUST BE GE ZERO', 2, 2)
C
      DGAMI = 0.D0
      IF (X.EQ.0.0D0) RETURN
C
C THE ONLY ERROR POSSIBLE IN THE EXPRESSION BELOW IS A FATAL OVERFLOW.
      FACTOR = EXP (DLNGAM(A) + A*LOG(X))
C
      DGAMI = FACTOR * DGAMIT (A, X)
C
      RETURN
      END
