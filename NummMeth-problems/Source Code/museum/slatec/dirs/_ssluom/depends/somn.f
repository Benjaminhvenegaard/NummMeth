      SUBROUTINE SOMN (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,
     +   NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, AP,
     +   EMAP, DZ, CSAV, RWORK, IWORK)
C     .. Scalar Arguments ..
      REAL ERR, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, N, NELT, NSAVE
C     .. Array Arguments ..
      REAL A(NELT), AP(N,0:NSAVE), B(N), CSAV(NSAVE), DZ(N),
     +     EMAP(N,0:NSAVE), P(N,0:NSAVE), R(N), RWORK(*), X(N), Z(N)
      INTEGER IA(NELT), IWORK(*), JA(NELT)
C     .. Subroutine Arguments ..
      EXTERNAL MATVEC, MSOLVE
C     .. Local Scalars ..
      REAL AK, AKDEN, AKNUM, BKL, BNRM, FUZZ, SOLNRM
      INTEGER I, IP, IPO, K, L, LMAX
C     .. External Functions ..
      REAL R1MACH, SDOT
      INTEGER ISSOMN
      EXTERNAL R1MACH, SDOT, ISSOMN
C     .. External Subroutines ..
      EXTERNAL SAXPY, SCOPY
C     .. Intrinsic Functions ..
      INTRINSIC ABS, MIN, MOD
C***FIRST EXECUTABLE STATEMENT  SOMN
C
C         Check some of the input data.
C
      ITER = 0
      IERR = 0
      IF( N.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      FUZZ = R1MACH(3)
      IF( TOL.LT.500*FUZZ ) THEN
         TOL = 500*FUZZ
         IERR = 4
      ENDIF
      FUZZ = FUZZ*FUZZ
C
C         Calculate initial residual and pseudo-residual, and check
C         stopping criterion.
      CALL MATVEC(N, X, R, NELT, IA, JA, A, ISYM)
      DO 10 I = 1, N
         R(I)  = B(I) - R(I)
 10   CONTINUE
      CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C
      IF( ISSOMN(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, NSAVE,
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $     R, Z, P, AP, EMAP, DZ, CSAV,
     $     RWORK, IWORK, AK, BNRM, SOLNRM) .NE. 0 ) GO TO 200
      IF( IERR.NE.0 ) RETURN
C
C
C         ***** iteration loop *****
C
CVD$R NOVECTOR
CVD$R NOCONCUR
      DO 100 K = 1, ITMAX
         ITER = K
         IP = MOD( ITER-1, NSAVE+1 )
C
C         calculate direction vector p, a*p, and (m-inv)*a*p,
C         and save if desired.
         CALL SCOPY(N, Z, 1, P(1,IP), 1)
         CALL MATVEC(N, P(1,IP), AP(1,IP), NELT, IA, JA, A, ISYM)
         CALL MSOLVE(N, AP(1,IP), EMAP(1,IP), NELT, IA, JA, A, ISYM,
     $        RWORK, IWORK)
         IF( NSAVE.EQ.0 ) THEN
            AKDEN = SDOT(N, EMAP, 1, EMAP, 1)
         ELSE
            IF( ITER.GT.1 ) THEN
               LMAX = MIN( NSAVE, ITER-1 )
               DO 20 L = 1, LMAX
                  IPO = MOD(IP+(NSAVE+1-L),NSAVE+1)
                  BKL = SDOT(N, EMAP(1,IP), 1, EMAP(1,IPO), 1)
                  BKL = BKL*CSAV(L)
                  CALL SAXPY(N, -BKL,    P(1,IPO), 1,    P(1,IP), 1)
                  CALL SAXPY(N, -BKL,   AP(1,IPO), 1,   AP(1,IP), 1)
                  CALL SAXPY(N, -BKL, EMAP(1,IPO), 1, EMAP(1,IP), 1)
 20            CONTINUE
               IF( NSAVE.GT.1 ) THEN
                  DO 30 L = NSAVE-1, 1, -1
                     CSAV(L+1) = CSAV(L)
 30               CONTINUE
               ENDIF
            ENDIF
            AKDEN = SDOT(N, EMAP(1,IP), 1, EMAP(1,IP), 1)
            IF( ABS(AKDEN).LT.FUZZ ) THEN
               IERR = 6
               RETURN
            ENDIF
            CSAV(1) = 1.0E0/AKDEN
C
C         calculate coefficient ak, new iterate x, new residual r, and
C         new pseudo-residual z.
         ENDIF
         AKNUM = SDOT(N, Z, 1, EMAP(1,IP), 1)
         AK = AKNUM/AKDEN
         CALL SAXPY(N,  AK,    P(1,IP), 1, X, 1)
         CALL SAXPY(N, -AK,   AP(1,IP), 1, R, 1)
         CALL SAXPY(N, -AK, EMAP(1,IP), 1, Z, 1)
C
C         check stopping criterion.
         IF( ISSOMN(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, NSAVE,
     $        ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $        R, Z, P, AP, EMAP, DZ, CSAV,
     $        RWORK, IWORK, AK, BNRM, SOLNRM) .NE. 0 ) GO TO 200
C
 100  CONTINUE
C
C         *****   end of loop  *****
C
C         Stopping criterion not satisfied.
      ITER = ITMAX + 1
      IERR = 2
C
 200  RETURN
C------------- LAST LINE OF SOMN FOLLOWS ----------------------------
      END