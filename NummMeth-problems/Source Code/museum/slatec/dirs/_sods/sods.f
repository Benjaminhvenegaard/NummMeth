      SUBROUTINE SODS (A, X, B, NEQ, NUK, NRDA, IFLAG, WORK, IWORK)
      DIMENSION A(NRDA,*),X(*),B(*),WORK(*),IWORK(*)
C
C***FIRST EXECUTABLE STATEMENT  SODS
      ITER=0
      IS=2
      IP=3
      KS=2
      KD=3
      KZ=KD+NUK
      KV=KZ+NUK
      KT=KV+NUK
      KC=KT+NUK
C
      CALL LSSODS(A,X,B,NEQ,NUK,NRDA,IFLAG,IWORK(1),IWORK(IS),A,
     1            WORK(KD),IWORK(IP),ITER,WORK(1),WORK(KS),
     2            WORK(KZ),B,WORK(KV),WORK(KT),WORK(KC))
C
      RETURN
      END
