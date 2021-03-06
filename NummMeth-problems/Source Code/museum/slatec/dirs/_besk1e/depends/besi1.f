      FUNCTION BESI1 (X)
      DIMENSION BI1CS(11)
      LOGICAL FIRST
      SAVE BI1CS, NTI1, XMIN, XSML, XMAX, FIRST
      DATA BI1CS( 1) /   -.0019717132 61099859E0 /
      DATA BI1CS( 2) /    .4073488766 7546481E0 /
      DATA BI1CS( 3) /    .0348389942 99959456E0 /
      DATA BI1CS( 4) /    .0015453945 56300123E0 /
      DATA BI1CS( 5) /    .0000418885 21098377E0 /
      DATA BI1CS( 6) /    .0000007649 02676483E0 /
      DATA BI1CS( 7) /    .0000000100 42493924E0 /
      DATA BI1CS( 8) /    .0000000000 99322077E0 /
      DATA BI1CS( 9) /    .0000000000 00766380E0 /
      DATA BI1CS(10) /    .0000000000 00004741E0 /
      DATA BI1CS(11) /    .0000000000 00000024E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  BESI1
      IF (FIRST) THEN
         NTI1 = INITS (BI1CS, 11, 0.1*R1MACH(3))
         XMIN = 2.0*R1MACH(1)
         XSML = SQRT (4.5*R1MACH(3))
         XMAX = LOG (R1MACH(2))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.3.0) GO TO 20
C
      BESI1 = 0.0
      IF (Y.EQ.0.0)  RETURN
C
      IF (Y .LE. XMIN) CALL XERMSG ('SLATEC', 'BESI1',
     +   'ABS(X) SO SMALL I1 UNDERFLOWS', 1, 1)
      IF (Y.GT.XMIN) BESI1 = 0.5*X
      IF (Y.GT.XSML) BESI1 = X * (.875 + CSEVL(Y*Y/4.5-1., BI1CS, NTI1))
      RETURN
C
 20   IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'BESI1',
     +   'ABS(X) SO BIG I1 OVERFLOWS', 2, 2)
C
      BESI1 = EXP(Y) * BESI1E(X)
C
      RETURN
      END
