      SUBROUTINE DQAWF (F, A, OMEGA, INTEGR, EPSABS, RESULT, ABSERR,
     +   NEVAL, IER, LIMLST, LST, LENIW, MAXP1, LENW, IWORK, WORK)
C
       DOUBLE PRECISION A,ABSERR,EPSABS,F,OMEGA,RESULT,WORK
       INTEGER IER,INTEGR,IWORK,LENIW,LENW,LIMIT,LIMLST,LL2,LVL,
     1  LST,L1,L2,L3,L4,L5,L6,MAXP1,NEVAL
C
       DIMENSION IWORK(*),WORK(*)
C
       EXTERNAL F
C
C         CHECK VALIDITY OF LIMLST, LENIW, MAXP1 AND LENW.
C
C***FIRST EXECUTABLE STATEMENT  DQAWF
      IER = 6
      NEVAL = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF(LIMLST.LT.3.OR.LENIW.LT.(LIMLST+2).OR.MAXP1.LT.1.OR.LENW.LT.
     1   (LENIW*2+MAXP1*25)) GO TO 10
C
C         PREPARE CALL FOR DQAWFE
C
      LIMIT = (LENIW-LIMLST)/2
      L1 = LIMLST+1
      L2 = LIMLST+L1
      L3 = LIMIT+L2
      L4 = LIMIT+L3
      L5 = LIMIT+L4
      L6 = LIMIT+L5
      LL2 = LIMIT+L1
      CALL DQAWFE(F,A,OMEGA,INTEGR,EPSABS,LIMLST,LIMIT,MAXP1,RESULT,
     1  ABSERR,NEVAL,IER,WORK(1),WORK(L1),IWORK(1),LST,WORK(L2),
     2  WORK(L3),WORK(L4),WORK(L5),IWORK(L1),IWORK(LL2),WORK(L6))
C
C         CALL ERROR HANDLER IF NECESSARY
C
      LVL = 0
10    IF(IER.EQ.6) LVL = 1
      IF (IER .NE. 0) CALL XERMSG ('SLATEC', 'DQAWF',
     +   'ABNORMAL RETURN', IER, LVL)
      RETURN
      END
