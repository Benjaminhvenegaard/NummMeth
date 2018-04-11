      SUBROUTINE DQAWS (F, A, B, ALFA, BETA, INTEGR, EPSABS, EPSREL,
     +   RESULT, ABSERR, NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
C
      DOUBLE PRECISION A,ABSERR,ALFA,B,BETA,EPSABS,EPSREL,F,RESULT,WORK
      INTEGER IER,INTEGR,IWORK,LAST,LENW,LIMIT,LVL,L1,L2,L3,NEVAL
C
      DIMENSION IWORK(*),WORK(*)
C
      EXTERNAL F
C
C         CHECK VALIDITY OF LIMIT AND LENW.
C
C***FIRST EXECUTABLE STATEMENT  DQAWS
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF(LIMIT.LT.2.OR.LENW.LT.LIMIT*4) GO TO 10
C
C         PREPARE CALL FOR DQAWSE.
C
      L1 = LIMIT+1
      L2 = LIMIT+L1
      L3 = LIMIT+L2
C
      CALL DQAWSE(F,A,B,ALFA,BETA,INTEGR,EPSABS,EPSREL,LIMIT,RESULT,
     1  ABSERR,NEVAL,IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),IWORK,LAST)
C
C         CALL ERROR HANDLER IF NECESSARY.
C
      LVL = 0
10    IF(IER.EQ.6) LVL = 1
      IF (IER .NE. 0) CALL XERMSG ('SLATEC', 'DQAWS',
     +   'ABNORMAL RETURN', IER, LVL)
      RETURN
      END
