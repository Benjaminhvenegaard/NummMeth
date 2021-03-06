      FUNCTION EXPREL (X)
      LOGICAL FIRST
      SAVE NTERMS, XBND, FIRST
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  EXPREL
      IF (FIRST) THEN
         ALNEPS = LOG(R1MACH(3))
         XN = 3.72 - 0.3*ALNEPS
         XLN = LOG((XN+1.0)/1.36)
         NTERMS = XN - (XN*XLN+ALNEPS)/(XLN+1.36) + 1.5
         XBND = R1MACH(3)
      ENDIF
      FIRST = .FALSE.
C
      ABSX = ABS(X)
      IF (ABSX.GT.0.5) EXPREL = (EXP(X) - 1.0) / X
      IF (ABSX.GT.0.5) RETURN
C
      EXPREL = 1.0
      IF (ABSX.LT.XBND) RETURN
C
      EXPREL = 0.0
      DO 20 I=1,NTERMS
        EXPREL = 1.0 + EXPREL*X/(NTERMS+2-I)
 20   CONTINUE
C
      RETURN
      END
