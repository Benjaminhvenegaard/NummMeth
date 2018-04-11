      SUBROUTINE IVOUT (N, IX, IFMT, IDIGIT)
      DIMENSION IX(*)
      CHARACTER IFMT*(*)
C
C     GET THE UNIT NUMBER WHERE OUTPUT WILL BE WRITTEN.
C***FIRST EXECUTABLE STATEMENT  IVOUT
      J=2
      LOUT=I1MACH(J)
      WRITE(LOUT,IFMT)
      IF(N.LE.0) RETURN
      NDIGIT = IDIGIT
      IF(IDIGIT.EQ.0) NDIGIT = 4
      IF(IDIGIT.GE.0) GO TO 80
C
      NDIGIT = -IDIGIT
      IF(NDIGIT.GT.4) GO TO 20
C
      DO 10 K1=1,N,10
      K2 = MIN(N,K1+9)
      WRITE(LOUT,1000) K1,K2,(IX(I),I=K1,K2)
   10 CONTINUE
      RETURN
C
   20 CONTINUE
      IF(NDIGIT.GT.6) GO TO 40
C
      DO 30 K1=1,N,7
      K2 = MIN(N,K1+6)
      WRITE(LOUT,1001) K1,K2,(IX(I),I=K1,K2)
   30 CONTINUE
      RETURN
C
   40 CONTINUE
      IF(NDIGIT.GT.10) GO TO 60
C
      DO 50 K1=1,N,5
      K2=MIN(N,K1+4)
      WRITE(LOUT,1002) K1,K2,(IX(I),I=K1,K2)
   50 CONTINUE
      RETURN
C
   60 CONTINUE
      DO 70 K1=1,N,3
      K2 = MIN(N,K1+2)
      WRITE(LOUT,1003) K1,K2,(IX(I),I=K1,K2)
   70 CONTINUE
      RETURN
C
   80 CONTINUE
      IF(NDIGIT.GT.4) GO TO 100
C
      DO 90 K1=1,N,20
      K2 = MIN(N,K1+19)
      WRITE(LOUT,1000) K1,K2,(IX(I),I=K1,K2)
   90 CONTINUE
      RETURN
C
  100 CONTINUE
      IF(NDIGIT.GT.6) GO TO 120
C
      DO 110 K1=1,N,15
      K2 = MIN(N,K1+14)
      WRITE(LOUT,1001) K1,K2,(IX(I),I=K1,K2)
  110 CONTINUE
      RETURN
C
  120 CONTINUE
      IF(NDIGIT.GT.10) GO TO 140
C
      DO 130 K1=1,N,10
      K2 = MIN(N,K1+9)
      WRITE(LOUT,1002) K1,K2,(IX(I),I=K1,K2)
  130 CONTINUE
      RETURN
C
  140 CONTINUE
      DO 150 K1=1,N,7
      K2 = MIN(N,K1+6)
      WRITE(LOUT,1003) K1,K2,(IX(I),I=K1,K2)
  150 CONTINUE
      RETURN
 1000 FORMAT(1X,I4,' - ',I4,20(1X,I5))
 1001 FORMAT(1X,I4,' - ',I4,15(1X,I7))
 1002 FORMAT(1X,I4,' - ',I4,10(1X,I11))
 1003 FORMAT(1X,I4,' - ',I4,7(1X,I15))
      END