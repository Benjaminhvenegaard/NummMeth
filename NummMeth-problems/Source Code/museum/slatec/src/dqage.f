      SUBROUTINE DQAGE (F, A, B, EPSABS, EPSREL, KEY, LIMIT, RESULT,
     +   ABSERR, NEVAL, IER, ALIST, BLIST, RLIST, ELIST, IORD, LAST)
C
      DOUBLE PRECISION A,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,A2,B,
     1  BLIST,B1,B2,DEFABS,DEFAB1,DEFAB2,D1MACH,ELIST,EPMACH,
     2  EPSABS,EPSREL,ERRBND,ERRMAX,ERROR1,ERROR2,ERRO12,ERRSUM,F,
     3  RESABS,RESULT,RLIST,UFLOW
      INTEGER IER,IORD,IROFF1,IROFF2,K,KEY,KEYF,LAST,LIMIT,MAXERR,NEVAL,
     1  NRMAX
C
      DIMENSION ALIST(*),BLIST(*),ELIST(*),IORD(*),
     1  RLIST(*)
C
      EXTERNAL F
C
C            LIST OF MAJOR VARIABLES
C            -----------------------
C
C           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
C                      (ALIST(I),BLIST(I))
C           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
C           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST
C                       ERROR ESTIMATE
C           ERRMAX    - ELIST(MAXERR)
C           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
C           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
C           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
C                       ABS(RESULT))
C           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
C           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
C           LAST      - INDEX FOR SUBDIVISION
C
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH  IS THE LARGEST RELATIVE SPACING.
C           UFLOW  IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQAGE
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
C
C           TEST ON VALIDITY OF PARAMETERS
C           ------------------------------
C
      IER = 0
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = 0.0D+00
      ELIST(1) = 0.0D+00
      IORD(1) = 0
      IF(EPSABS.LE.0.0D+00.AND.
     1  EPSREL.LT.MAX(0.5D+02*EPMACH,0.5D-28)) IER = 6
      IF(IER.EQ.6) GO TO 999
C
C           FIRST APPROXIMATION TO THE INTEGRAL
C           -----------------------------------
C
      KEYF = KEY
      IF(KEY.LE.0) KEYF = 1
      IF(KEY.GE.7) KEYF = 6
      NEVAL = 0
      IF(KEYF.EQ.1) CALL DQK15(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF(KEYF.EQ.2) CALL DQK21(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF(KEYF.EQ.3) CALL DQK31(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF(KEYF.EQ.4) CALL DQK41(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF(KEYF.EQ.5) CALL DQK51(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF(KEYF.EQ.6) CALL DQK61(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
C
C           TEST ON ACCURACY.
C
      ERRBND = MAX(EPSABS,EPSREL*ABS(RESULT))
      IF(ABSERR.LE.0.5D+02*EPMACH*DEFABS.AND.ABSERR.GT.ERRBND) IER = 2
      IF(LIMIT.EQ.1) IER = 1
      IF(IER.NE.0.OR.(ABSERR.LE.ERRBND.AND.ABSERR.NE.RESABS)
     1  .OR.ABSERR.EQ.0.0D+00) GO TO 60
C
C           INITIALIZATION
C           --------------
C
C
      ERRMAX = ABSERR
      MAXERR = 1
      AREA = RESULT
      ERRSUM = ABSERR
      NRMAX = 1
      IROFF1 = 0
      IROFF2 = 0
C
C           MAIN DO-LOOP
C           ------------
C
      DO 30 LAST = 2,LIMIT
C
C           BISECT THE SUBINTERVAL WITH THE LARGEST ERROR ESTIMATE.
C
        A1 = ALIST(MAXERR)
        B1 = 0.5D+00*(ALIST(MAXERR)+BLIST(MAXERR))
        A2 = B1
        B2 = BLIST(MAXERR)
        IF(KEYF.EQ.1) CALL DQK15(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.2) CALL DQK21(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.3) CALL DQK31(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.4) CALL DQK41(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.5) CALL DQK51(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.6) CALL DQK61(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.1) CALL DQK15(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
        IF(KEYF.EQ.2) CALL DQK21(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
        IF(KEYF.EQ.3) CALL DQK31(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
        IF(KEYF.EQ.4) CALL DQK41(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
        IF(KEYF.EQ.5) CALL DQK51(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
        IF(KEYF.EQ.6) CALL DQK61(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
C
C           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
C           AND ERROR AND TEST FOR ACCURACY.
C
        NEVAL = NEVAL+1
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2) GO TO 5
        IF(ABS(RLIST(MAXERR)-AREA12).LE.0.1D-04*ABS(AREA12)
     1  .AND.ERRO12.GE.0.99D+00*ERRMAX) IROFF1 = IROFF1+1
        IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF2 = IROFF2+1
    5   RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
        IF(ERRSUM.LE.ERRBND) GO TO 8
C
C           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
C
        IF(IROFF1.GE.6.OR.IROFF2.GE.20) IER = 2
C
C           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF SUBINTERVALS
C           EQUALS LIMIT.
C
        IF(LAST.EQ.LIMIT) IER = 1
C
C           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
C           AT A POINT OF THE INTEGRATION RANGE.
C
        IF(MAX(ABS(A1),ABS(B2)).LE.(0.1D+01+0.1D+03*
     1  EPMACH)*(ABS(A2)+0.1D+04*UFLOW)) IER = 3
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
    8   IF(ERROR2.GT.ERROR1) GO TO 10
        ALIST(LAST) = A2
        BLIST(MAXERR) = B1
        BLIST(LAST) = B2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 20
   10   ALIST(MAXERR) = A2
        ALIST(LAST) = A1
        BLIST(LAST) = B1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
C           WITH THE LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
C
   20   CALL DQPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
C ***JUMP OUT OF DO-LOOP
        IF(IER.NE.0.OR.ERRSUM.LE.ERRBND) GO TO 40
   30 CONTINUE
C
C           COMPUTE FINAL RESULT.
C           ---------------------
C
   40 RESULT = 0.0D+00
      DO 50 K=1,LAST
        RESULT = RESULT+RLIST(K)
   50 CONTINUE
      ABSERR = ERRSUM
   60 IF(KEYF.NE.1) NEVAL = (10*KEYF+1)*(2*NEVAL+1)
      IF(KEYF.EQ.1) NEVAL = 30*NEVAL+15
  999 RETURN
      END
