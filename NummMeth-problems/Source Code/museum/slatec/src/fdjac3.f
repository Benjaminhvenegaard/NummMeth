      SUBROUTINE FDJAC3 (FCN, M, N, X, FVEC, FJAC, LDFJAC, IFLAG,
     +   EPSFCN, WA)
      INTEGER M,N,LDFJAC,IFLAG
      REAL EPSFCN
      REAL X(*),FVEC(*),FJAC(LDFJAC,*),WA(*)
      INTEGER I,J
      REAL EPS,EPSMCH,H,TEMP,ZERO
      REAL R1MACH
      SAVE ZERO
      DATA ZERO /0.0E0/
C***FIRST EXECUTABLE STATEMENT  FDJAC3
      EPSMCH = R1MACH(4)
C
      EPS = SQRT(MAX(EPSFCN,EPSMCH))
C      SET IFLAG=1 TO INDICATE THAT FUNCTION VALUES
C           ARE TO BE RETURNED BY FCN.
      IFLAG = 1
      DO 20 J = 1, N
         TEMP = X(J)
         H = EPS*ABS(TEMP)
         IF (H .EQ. ZERO) H = EPS
         X(J) = TEMP + H
         CALL FCN(IFLAG,M,N,X,WA,FJAC,LDFJAC)
         IF (IFLAG .LT. 0) GO TO 30
         X(J) = TEMP
         DO 10 I = 1, M
            FJAC(I,J) = (WA(I) - FVEC(I))/H
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE FDJAC3.
C
      END