      FUNCTION CBRT (X)
      DIMENSION CBRT2(5)
      SAVE CBRT2, NITER
      DATA CBRT2(1) / 0.6299605249 4743658E0 /
      DATA CBRT2(2) / 0.7937005259 8409974E0 /
      DATA CBRT2(3) / 1.0E0 /
      DATA CBRT2(4) / 1.2599210498 9487316E0 /
      DATA CBRT2(5) / 1.5874010519 6819947E0 /
      DATA NITER / 0 /
C***FIRST EXECUTABLE STATEMENT  CBRT
      IF (NITER.EQ.0) NITER = 1.443*LOG(-.106*LOG(0.1*R1MACH(3))) + 1.
C
      CBRT = 0.0
      IF (X.EQ.0.) RETURN
C
      CALL R9UPAK (ABS(X), Y, N)
      IXPNT = N/3
      IREM = N - 3*IXPNT + 3
C
C THE APPROXIMATION BELOW IS A GENERALIZED CHEBYSHEV SERIES CONVERTED
C TO POLYNOMIAL FORM.  THE APPROX IS NEARLY BEST IN THE SENSE OF
C RELATIVE ERROR WITH 4.085 DIGITS ACCURACY.
C
      CBRT = .439581E0 + Y*(.928549E0 + Y*(-.512653E0 + Y*.144586E0))
C
      DO 10 ITER=1,NITER
        CBRTSQ = CBRT*CBRT
        CBRT = CBRT + (Y-CBRT*CBRTSQ)/(3.0*CBRTSQ)
 10   CONTINUE
C
      CBRT = R9PAK (CBRT2(IREM)*SIGN(CBRT,X), IXPNT)
      RETURN
C
      END