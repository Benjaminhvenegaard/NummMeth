      SUBROUTINE POLCOF (XX, N, X, C, D, WORK)
C
      DIMENSION X(*), C(*), D(*), WORK(*)
C***FIRST EXECUTABLE STATEMENT  POLCOF
      DO 10010 K=1,N
      D(K)=C(K)
10010 CONTINUE
      IF (N.EQ.1) RETURN
      WORK(1)=1.0
      PONE=C(1)
      NM1=N-1
      DO 10020 K=2,N
      KM1=K-1
      NPKM1=N+K-1
      WORK(NPKM1)=XX-X(KM1)
      WORK(K)=WORK(NPKM1)*WORK(KM1)
      PTWO=PONE+WORK(K)*C(K)
      PONE=PTWO
10020 CONTINUE
      D(1)=PTWO
      IF (N.EQ.2) RETURN
      DO 10030 K=2,NM1
      KM1=K-1
      KM2N=K-2+N
      NMKP1=N-K+1
      DO 10030 I=2,NMKP1
      KM2NPI=KM2N+I
      IM1=I-1
      KM1PI=KM1+I
      WORK(I)=WORK(KM2NPI)*WORK(IM1)+WORK(I)
      D(K)=D(K)+WORK(I)*D(KM1PI)
10030 CONTINUE
      RETURN
      END
