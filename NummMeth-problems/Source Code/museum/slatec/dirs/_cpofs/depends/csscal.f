      SUBROUTINE CSSCAL (N, SA, CX, INCX)
      COMPLEX CX(*)
      REAL SA
      INTEGER I, INCX, IX, N
C***FIRST EXECUTABLE STATEMENT  CSSCAL
      IF (N .LE. 0) RETURN
C
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        CX(IX) = SA*CX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
   20 DO 30 I = 1,N
        CX(I) = SA*CX(I)
   30 CONTINUE
      RETURN
      END
