      SUBROUTINE QAGPE (F, A, B, NPTS2, POINTS, EPSABS, EPSREL, LIMIT,
     +   RESULT, ABSERR, NEVAL, IER, ALIST, BLIST, RLIST, ELIST, PTS,
     +   IORD, LEVEL, NDIN, LAST)
      REAL A,ABSEPS,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,
     1  A2,B,BLIST,B1,B2,CORREC,DEFABS,DEFAB1,DEFAB2,
     2  DRES,R1MACH,ELIST,EPMACH,EPSABS,EPSREL,ERLARG,ERLAST,ERRBND,
     3  ERRMAX,ERROR1,ERRO12,ERROR2,ERRSUM,ERTEST,F,OFLOW,POINTS,PTS,
     4  RESA,RESABS,RESEPS,RESULT,RES3LA,RLIST,RLIST2,SIGN,TEMP,
     5  UFLOW
      INTEGER I,ID,IER,IERRO,IND1,IND2,IORD,IP1,IROFF1,IROFF2,
     1  IROFF3,J,JLOW,JUPBND,K,KSGN,KTMIN,LAST,LEVCUR,LEVEL,LEVMAX,
     2  LIMIT,MAXERR,NDIN,NEVAL,NINT,NINTP1,NPTS,NPTS2,NRES,
     3  NRMAX,NUMRL2
      LOGICAL EXTRAP,NOEXT
C
C
      DIMENSION ALIST(*),BLIST(*),ELIST(*),IORD(*),
     1  LEVEL(*),NDIN(*),POINTS(*),PTS(*),RES3LA(3),
     2  RLIST(*),RLIST2(52)
C
      EXTERNAL F
C
C            THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF
C            LIMEXP IN SUBROUTINE EPSALG (RLIST2 SHOULD BE OF DIMENSION
C            (LIMEXP+2) AT LEAST).
C
C
C            LIST OF MAJOR VARIABLES
C            -----------------------
C
C           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
C                       (ALIST(I),BLIST(I))
C           RLIST2    - ARRAY OF DIMENSION AT LEAST LIMEXP+2
C                       CONTAINING THE PART OF THE EPSILON TABLE WHICH
C                       IS STILL NEEDED FOR FURTHER COMPUTATIONS
C           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
C           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST ERROR
C                       ESTIMATE
C           ERRMAX    - ELIST(MAXERR)
C           ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
C                       (BEFORE THAT SUBDIVISION HAS TAKEN PLACE)
C           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
C           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
C           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
C                       ABS(RESULT))
C           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
C           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
C           LAST      - INDEX FOR SUBDIVISION
C           NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE
C           NUMRL2    - NUMBER OF ELEMENTS IN RLIST2. IF AN
C                       APPROPRIATE APPROXIMATION TO THE COMPOUNDED
C                       INTEGRAL HAS BEEN OBTAINED, IT IS PUT IN
C                       RLIST2(NUMRL2) AFTER NUMRL2 HAS BEEN INCREASED
C                       BY ONE.
C           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
C                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
C           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE
C                       IS ATTEMPTING TO PERFORM EXTRAPOLATION. I.E.
C                       BEFORE SUBDIVIDING THE SMALLEST INTERVAL WE
C                       TRY TO DECREASE THE VALUE OF ERLARG.
C           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION IS
C                       NO LONGER ALLOWED (TRUE-VALUE)
C
C            MACHINE DEPENDENT CONSTANTS
C            ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  QAGPE
      EPMACH = R1MACH(4)
C
C            TEST ON VALIDITY OF PARAMETERS
C            -----------------------------
C
      IER = 0
      NEVAL = 0
      LAST = 0
      RESULT = 0.0E+00
      ABSERR = 0.0E+00
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = 0.0E+00
      ELIST(1) = 0.0E+00
      IORD(1) = 0
      LEVEL(1) = 0
      NPTS = NPTS2-2
      IF(NPTS2.LT.2.OR.LIMIT.LE.NPTS.OR.(EPSABS.LE.0.0E+00.AND.
     1  EPSREL.LT.MAX(0.5E+02*EPMACH,0.5E-14))) IER = 6
      IF(IER.EQ.6) GO TO 999
C
C            IF ANY BREAK POINTS ARE PROVIDED, SORT THEM INTO AN
C            ASCENDING SEQUENCE.
C
      SIGN = 1.0E+00
      IF(A.GT.B) SIGN = -1.0E+00
      PTS(1) = MIN(A,B)
      IF(NPTS.EQ.0) GO TO 15
      DO 10 I = 1,NPTS
        PTS(I+1) = POINTS(I)
   10 CONTINUE
   15 PTS(NPTS+2) = MAX(A,B)
      NINT = NPTS+1
      A1 = PTS(1)
      IF(NPTS.EQ.0) GO TO 40
      NINTP1 = NINT+1
      DO 20 I = 1,NINT
        IP1 = I+1
        DO 20 J = IP1,NINTP1
          IF(PTS(I).LE.PTS(J)) GO TO 20
          TEMP = PTS(I)
          PTS(I) = PTS(J)
          PTS(J) = TEMP
   20 CONTINUE
      IF(PTS(1).NE.MIN(A,B).OR.PTS(NINTP1).NE.
     1  MAX(A,B)) IER = 6
      IF(IER.EQ.6) GO TO 999
C
C            COMPUTE FIRST INTEGRAL AND ERROR APPROXIMATIONS.
C            ------------------------------------------------
C
   40 RESABS = 0.0E+00
      DO 50 I = 1,NINT
        B1 = PTS(I+1)
        CALL QK21(F,A1,B1,AREA1,ERROR1,DEFABS,RESA)
        ABSERR = ABSERR+ERROR1
        RESULT = RESULT+AREA1
        NDIN(I) = 0
        IF(ERROR1.EQ.RESA.AND.ERROR1.NE.0.0E+00) NDIN(I) = 1
        RESABS = RESABS+DEFABS
        LEVEL(I) = 0
        ELIST(I) = ERROR1
        ALIST(I) = A1
        BLIST(I) = B1
        RLIST(I) = AREA1
        IORD(I) = I
        A1 = B1
   50 CONTINUE
      ERRSUM = 0.0E+00
      DO 55 I = 1,NINT
        IF(NDIN(I).EQ.1) ELIST(I) = ABSERR
        ERRSUM = ERRSUM+ELIST(I)
   55 CONTINUE
C
C           TEST ON ACCURACY.
C
      LAST = NINT
      NEVAL = 21*NINT
      DRES = ABS(RESULT)
      ERRBND = MAX(EPSABS,EPSREL*DRES)
      IF(ABSERR.LE.0.1E+03*EPMACH*RESABS.AND.ABSERR.GT.
     1  ERRBND) IER = 2
      IF(NINT.EQ.1) GO TO 80
      DO 70 I = 1,NPTS
        JLOW = I+1
        IND1 = IORD(I)
        DO 60 J = JLOW,NINT
          IND2 = IORD(J)
          IF(ELIST(IND1).GT.ELIST(IND2)) GO TO 60
          IND1 = IND2
          K = J
   60   CONTINUE
        IF(IND1.EQ.IORD(I)) GO TO 70
        IORD(K) = IORD(I)
        IORD(I) = IND1
   70 CONTINUE
      IF(LIMIT.LT.NPTS2) IER = 1
   80 IF(IER.NE.0.OR.ABSERR.LE.ERRBND) GO TO 999
C
C           INITIALIZATION
C           --------------
C
      RLIST2(1) = RESULT
      MAXERR = IORD(1)
      ERRMAX = ELIST(MAXERR)
      AREA = RESULT
      NRMAX = 1
      NRES = 0
      NUMRL2 = 1
      KTMIN = 0
      EXTRAP = .FALSE.
      NOEXT = .FALSE.
      ERLARG = ERRSUM
      ERTEST = ERRBND
      LEVMAX = 1
      IROFF1 = 0
      IROFF2 = 0
      IROFF3 = 0
      IERRO = 0
      UFLOW = R1MACH(1)
      OFLOW = R1MACH(2)
      ABSERR = OFLOW
      KSGN = -1
      IF(DRES.GE.(0.1E+01-0.5E+02*EPMACH)*RESABS) KSGN = 1
C
C           MAIN DO-LOOP
C           ------------
C
      DO 160 LAST = NPTS2,LIMIT
C
C           BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST
C           ERROR ESTIMATE.
C
        LEVCUR = LEVEL(MAXERR)+1
        A1 = ALIST(MAXERR)
        B1 = 0.5E+00*(ALIST(MAXERR)+BLIST(MAXERR))
        A2 = B1
        B2 = BLIST(MAXERR)
        ERLAST = ERRMAX
        CALL QK21(F,A1,B1,AREA1,ERROR1,RESA,DEFAB1)
        CALL QK21(F,A2,B2,AREA2,ERROR2,RESA,DEFAB2)
C
C           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
C           AND ERROR AND TEST FOR ACCURACY.
C
        NEVAL = NEVAL+42
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2) GO TO 95
        IF(ABS(RLIST(MAXERR)-AREA12).GT.0.1E-04*ABS(AREA12)
     1  .OR.ERRO12.LT.0.99E+00*ERRMAX) GO TO 90
        IF(EXTRAP) IROFF2 = IROFF2+1
        IF(.NOT.EXTRAP) IROFF1 = IROFF1+1
   90   IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF3 = IROFF3+1
   95   LEVEL(MAXERR) = LEVCUR
        LEVEL(LAST) = LEVCUR
        RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
C
C           TEST FOR ROUNDOFF ERROR AND EVENTUALLY
C           SET ERROR FLAG.
C
        IF(IROFF1+IROFF2.GE.10.OR.IROFF3.GE.20) IER = 2
        IF(IROFF2.GE.5) IERRO = 3
C
C           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF
C           SUBINTERVALS EQUALS LIMIT.
C
        IF(LAST.EQ.LIMIT) IER = 1
C
C           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
C           AT A POINT OF THE INTEGRATION RANGE
C
        IF(MAX(ABS(A1),ABS(B2)).LE.(0.1E+01+0.1E+03*EPMACH)*
     1  (ABS(A2)+0.1E+04*UFLOW)) IER = 4
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
        IF(ERROR2.GT.ERROR1) GO TO 100
        ALIST(LAST) = A2
        BLIST(MAXERR) = B1
        BLIST(LAST) = B2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 110
  100   ALIST(MAXERR) = A2
        ALIST(LAST) = A1
        BLIST(LAST) = B1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE QPSRT TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE
C           SUBINTERVAL WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE
C           BISECTED NEXT).
C
  110   CALL QPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
C ***JUMP OUT OF DO-LOOP
        IF(ERRSUM.LE.ERRBND) GO TO 190
C ***JUMP OUT OF DO-LOOP
        IF(IER.NE.0) GO TO 170
        IF(NOEXT) GO TO 160
        ERLARG = ERLARG-ERLAST
        IF(LEVCUR+1.LE.LEVMAX) ERLARG = ERLARG+ERRO12
        IF(EXTRAP) GO TO 120
C
C           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE
C           SMALLEST INTERVAL.
C
        IF(LEVEL(MAXERR)+1.LE.LEVMAX) GO TO 160
        EXTRAP = .TRUE.
        NRMAX = 2
  120   IF(IERRO.EQ.3.OR.ERLARG.LE.ERTEST) GO TO 140
C
C           THE SMALLEST INTERVAL HAS THE LARGEST ERROR.
C           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS
C           OVER THE LARGER INTERVALS (ERLARG) AND PERFORM
C           EXTRAPOLATION.
C
        ID = NRMAX
        JUPBND = LAST
        IF(LAST.GT.(2+LIMIT/2)) JUPBND = LIMIT+3-LAST
        DO 130 K = ID,JUPBND
          MAXERR = IORD(NRMAX)
          ERRMAX = ELIST(MAXERR)
C ***JUMP OUT OF DO-LOOP
          IF(LEVEL(MAXERR)+1.LE.LEVMAX) GO TO 160
          NRMAX = NRMAX+1
  130   CONTINUE
C
C           PERFORM EXTRAPOLATION.
C
  140   NUMRL2 = NUMRL2+1
        RLIST2(NUMRL2) = AREA
        IF(NUMRL2.LE.2) GO TO 155
        CALL QELG(NUMRL2,RLIST2,RESEPS,ABSEPS,RES3LA,NRES)
        KTMIN = KTMIN+1
        IF(KTMIN.GT.5.AND.ABSERR.LT.0.1E-02*ERRSUM) IER = 5
        IF(ABSEPS.GE.ABSERR) GO TO 150
        KTMIN = 0
        ABSERR = ABSEPS
        RESULT = RESEPS
        CORREC = ERLARG
        ERTEST = MAX(EPSABS,EPSREL*ABS(RESEPS))
C ***JUMP OUT OF DO-LOOP
        IF(ABSERR.LT.ERTEST) GO TO 170
C
C           PREPARE BISECTION OF THE SMALLEST INTERVAL.
C
  150   IF(NUMRL2.EQ.1) NOEXT = .TRUE.
        IF(IER.GE.5) GO TO 170
  155   MAXERR = IORD(1)
        ERRMAX = ELIST(MAXERR)
        NRMAX = 1
        EXTRAP = .FALSE.
        LEVMAX = LEVMAX+1
        ERLARG = ERRSUM
  160 CONTINUE
C
C           SET THE FINAL RESULT.
C           ---------------------
C
C
  170 IF(ABSERR.EQ.OFLOW) GO TO 190
      IF((IER+IERRO).EQ.0) GO TO 180
      IF(IERRO.EQ.3) ABSERR = ABSERR+CORREC
      IF(IER.EQ.0) IER = 3
      IF(RESULT.NE.0.0E+00.AND.AREA.NE.0.0E+00)GO TO 175
      IF(ABSERR.GT.ERRSUM)GO TO 190
      IF(AREA.EQ.0.0E+00) GO TO 210
      GO TO 180
  175 IF(ABSERR/ABS(RESULT).GT.ERRSUM/ABS(AREA))GO TO 190
C
C           TEST ON DIVERGENCE.
C
  180 IF(KSGN.EQ.(-1).AND.MAX(ABS(RESULT),ABS(AREA)).LE.
     1  DEFABS*0.1E-01) GO TO 210
      IF(0.1E-01.GT.(RESULT/AREA).OR.(RESULT/AREA).GT.0.1E+03.OR.
     1  ERRSUM.GT.ABS(AREA)) IER = 6
      GO TO 210
C
C           COMPUTE GLOBAL INTEGRAL SUM.
C
  190 RESULT = 0.0E+00
      DO 200 K = 1,LAST
        RESULT = RESULT+RLIST(K)
  200 CONTINUE
      ABSERR = ERRSUM
  210 IF(IER.GT.2) IER = IER - 1
      RESULT = RESULT*SIGN
 999  RETURN
      END
