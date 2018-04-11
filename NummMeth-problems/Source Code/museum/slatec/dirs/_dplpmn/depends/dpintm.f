      SUBROUTINE DPINTM (M, N, SX, IX, LMX, IPAGEF)
      DOUBLE PRECISION SX(*),ZERO,ONE
      DIMENSION IX(*)
      SAVE ZERO, ONE
      DATA ZERO,ONE /0.D0,1.D0/
C***FIRST EXECUTABLE STATEMENT  DPINTM
      IOPT=1
C
C     CHECK FOR INPUT ERRORS.
C
      IF (.NOT.(M.LE.0 .OR. N.LE.0)) GO TO 20002
      NERR=55
      CALL XERMSG ('SLATEC', 'DPINTM',
     +   'MATRIX DIMENSION M OR N .LE. 0', NERR, IOPT)
C
C     VERIFY IF VALUE OF LMX IS LARGE ENOUGH.
C
20002 IF (.NOT.(LMX.LT.N+7)) GO TO 20005
      NERR=55
      CALL XERMSG ('SLATEC', 'DPINTM',
     +   'THE VALUE OF LMX IS TOO SMALL', NERR, IOPT)
C
C     INITIALIZE DATA STRUCTURE INDEPENDENT VALUES.
C
20005 SX(1)=ZERO
      SX(2)=ZERO
      SX(3)=IPAGEF
      IX(1)=LMX
      IX(2)=M
      IX(3)=N
      IX(4)=0
      SX(LMX-1)=ZERO
      SX(LMX)=-ONE
      IX(LMX-1)=-1
      LP4=N+4
C
C     INITIALIZE DATA STRUCTURE DEPENDENT VALUES.
C
      I=4
      N20008=LP4
      GO TO 20009
20008 I=I+1
20009 IF ((N20008-I).LT.0) GO TO 20010
      SX(I)=ZERO
      GO TO 20008
20010 I=5
      N20012=LP4
      GO TO 20013
20012 I=I+1
20013 IF ((N20012-I).LT.0) GO TO 20014
      IX(I)=LP4
      GO TO 20012
20014 SX(N+5)=ZERO
      IX(N+5)=0
      IX(LMX)=0
C
C     INITIALIZATION COMPLETE.
C
      RETURN
      END
