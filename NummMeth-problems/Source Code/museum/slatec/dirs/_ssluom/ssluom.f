      SUBROUTINE SSLUOM (N, B, X, NELT, IA, JA, A, ISYM, NSAVE, ITOL,
     +   TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW)
C     .. Parameters ..
      INTEGER LOCRB, LOCIB
      PARAMETER (LOCRB=1, LOCIB=11)
C     .. Scalar Arguments ..
      REAL ERR, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, LENIW, LENW, N,
     +        NELT, NSAVE
C     .. Array Arguments ..
      REAL A(N), B(N), RWORK(LENW), X(N)
      INTEGER IA(NELT), IWORK(LENIW), JA(NELT)
C     .. Local Scalars ..
      INTEGER ICOL, J, JBGN, JEND, LOCAP, LOCCSA, LOCDIN, LOCDZ, LOCEMA,
     +        LOCIL, LOCIU, LOCIW, LOCJL, LOCJU, LOCL, LOCNC, LOCNR,
     +        LOCP, LOCR, LOCU, LOCW, LOCZ, NL, NU
C     .. External Subroutines ..
      EXTERNAL SCHKW, SOMN, SS2Y, SSILUS, SSLUI, SSMV
C***FIRST EXECUTABLE STATEMENT  SSLUOM
C
      IERR = 0
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
C
C         Change the SLAP input matrix IA, JA, A to SLAP-Column format.
      CALL SS2Y( N, NELT, IA, JA, A, ISYM )
C
C         Count number of Non-Zero elements preconditioner ILU matrix.
C         Then set up the work arrays.
      NL = 0
      NU = 0
      DO 20 ICOL = 1, N
C         Don't count diagonal.
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
CVD$ NOVECTOR
            DO 10 J = JBGN, JEND
               IF( IA(J).GT.ICOL ) THEN
                  NL = NL + 1
                  IF( ISYM.NE.0 ) NU = NU + 1
               ELSE
                  NU = NU + 1
               ENDIF
 10         CONTINUE
         ENDIF
 20   CONTINUE
C
      LOCIL = LOCIB
      LOCJL = LOCIL + N+1
      LOCIU = LOCJL + NL
      LOCJU = LOCIU + NU
      LOCNR = LOCJU + N+1
      LOCNC = LOCNR + N
      LOCIW = LOCNC + N
C
      LOCL   = LOCRB
      LOCDIN = LOCL + NL
      LOCU   = LOCDIN + N
      LOCR   = LOCU + NU
      LOCZ   = LOCR + N
      LOCP   = LOCZ + N
      LOCAP  = LOCP + N*(NSAVE+1)
      LOCEMA = LOCAP + N*(NSAVE+1)
      LOCDZ  = LOCEMA + N*(NSAVE+1)
      LOCCSA = LOCDZ + N
      LOCW   = LOCCSA + NSAVE
C
C         Check the workspace allocations.
      CALL SCHKW( 'SSLUOM', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
      IF( IERR.NE.0 ) RETURN
C
      IWORK(1) = LOCIL
      IWORK(2) = LOCJL
      IWORK(3) = LOCIU
      IWORK(4) = LOCJU
      IWORK(5) = LOCL
      IWORK(6) = LOCDIN
      IWORK(7) = LOCU
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
C
C         Compute the Incomplete LU decomposition.
      CALL SSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL),
     $     IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU),
     $     IWORK(LOCJU), RWORK(LOCU), IWORK(LOCNR), IWORK(LOCNC) )
C
C         Perform the incomplete LU preconditioned OrthoMin algorithm.
      CALL SOMN(N, B, X, NELT, IA, JA, A, ISYM, SSMV,
     $     SSLUI, NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $     RWORK(LOCR), RWORK(LOCZ), RWORK(LOCP), RWORK(LOCAP),
     $     RWORK(LOCEMA), RWORK(LOCDZ), RWORK(LOCCSA),
     $     RWORK, IWORK )
      RETURN
      END