      SUBROUTINE SDAJAC (NEQ, X, Y, YPRIME, DELTA, CJ, H, IER, WT, E,
     *   WM, IWM, RES, IRES, UROUND, JAC, RPAR, IPAR, NTEMP)
C
      INTEGER  NEQ, IER, IWM(*), IRES, IPAR(*), NTEMP
      REAL  X, Y(*), YPRIME(*), DELTA(*), CJ, H, WT(*), E(*), WM(*),
     *   UROUND, RPAR(*)
      EXTERNAL  RES, JAC
C
      EXTERNAL  SGBFA, SGEFA
C
      INTEGER  I, I1, I2, II, IPSAVE, ISAVE, J, K, L, LENPD, LIPVT,
     *   LML, LMTYPE, LMU, MBA, MBAND, MEB1, MEBAND, MSAVE, MTYPE, N,
     *   NPD, NPDM1, NROW
      REAL  DEL, DELINV, SQUR, YPSAVE, YSAVE
C
      PARAMETER (NPD=1)
      PARAMETER (LML=1)
      PARAMETER (LMU=2)
      PARAMETER (LMTYPE=4)
      PARAMETER (LIPVT=21)
C
C***FIRST EXECUTABLE STATEMENT  SDAJAC
      IER = 0
      NPDM1=NPD-1
      MTYPE=IWM(LMTYPE)
      GO TO (100,200,300,400,500),MTYPE
C
C
C     DENSE USER-SUPPLIED MATRIX
100   LENPD=NEQ*NEQ
      DO 110 I=1,LENPD
110      WM(NPDM1+I)=0.0E0
      CALL JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
      GO TO 230
C
C
C     DENSE FINITE-DIFFERENCE-GENERATED MATRIX
200   IRES=0
      NROW=NPDM1
      SQUR = SQRT(UROUND)
      DO 210 I=1,NEQ
         DEL=SQUR*MAX(ABS(Y(I)),ABS(H*YPRIME(I)),ABS(WT(I)))
         DEL=SIGN(DEL,H*YPRIME(I))
         DEL=(Y(I)+DEL)-Y(I)
         YSAVE=Y(I)
         YPSAVE=YPRIME(I)
         Y(I)=Y(I)+DEL
         YPRIME(I)=YPRIME(I)+CJ*DEL
         CALL RES(X,Y,YPRIME,E,IRES,RPAR,IPAR)
         IF (IRES .LT. 0) RETURN
         DELINV=1.0E0/DEL
         DO 220 L=1,NEQ
220      WM(NROW+L)=(E(L)-DELTA(L))*DELINV
      NROW=NROW+NEQ
      Y(I)=YSAVE
      YPRIME(I)=YPSAVE
210   CONTINUE
C
C
C     DO DENSE-MATRIX LU DECOMPOSITION ON PD
230      CALL SGEFA(WM(NPD),NEQ,NEQ,IWM(LIPVT),IER)
      RETURN
C
C
C     DUMMY SECTION FOR IWM(MTYPE)=3
300   RETURN
C
C
C     BANDED USER-SUPPLIED MATRIX
400   LENPD=(2*IWM(LML)+IWM(LMU)+1)*NEQ
      DO 410 I=1,LENPD
410      WM(NPDM1+I)=0.0E0
      CALL JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
      MEBAND=2*IWM(LML)+IWM(LMU)+1
      GO TO 550
C
C
C     BANDED FINITE-DIFFERENCE-GENERATED MATRIX
500   MBAND=IWM(LML)+IWM(LMU)+1
      MBA=MIN(MBAND,NEQ)
      MEBAND=MBAND+IWM(LML)
      MEB1=MEBAND-1
      MSAVE=(NEQ/MBAND)+1
      ISAVE=NTEMP-1
      IPSAVE=ISAVE+MSAVE
      IRES=0
      SQUR=SQRT(UROUND)
      DO 540 J=1,MBA
         DO 510 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          WM(ISAVE+K)=Y(N)
          WM(IPSAVE+K)=YPRIME(N)
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),ABS(WT(N)))
          DEL=SIGN(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          Y(N)=Y(N)+DEL
510       YPRIME(N)=YPRIME(N)+CJ*DEL
      CALL RES(X,Y,YPRIME,E,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) RETURN
      DO 530 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          Y(N)=WM(ISAVE+K)
          YPRIME(N)=WM(IPSAVE+K)
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),ABS(WT(N)))
          DEL=SIGN(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          DELINV=1.0E0/DEL
          I1=MAX(1,(N-IWM(LMU)))
          I2=MIN(NEQ,(N+IWM(LML)))
          II=N*MEB1-IWM(LML)+NPDM1
          DO 520 I=I1,I2
520         WM(II+I)=(E(I)-DELTA(I))*DELINV
530      CONTINUE
540   CONTINUE
C
C
C     DO LU DECOMPOSITION OF BANDED PD
550   CALL SGBFA(WM(NPD),MEBAND,NEQ,
     *    IWM(LML),IWM(LMU),IWM(LIPVT),IER)
      RETURN
C------END OF SUBROUTINE SDAJAC------
      END
