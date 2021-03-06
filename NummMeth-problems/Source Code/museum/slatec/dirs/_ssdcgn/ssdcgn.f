      SUBROUTINE SSDCGN (N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,
     +   ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW)
C     .. Parameters ..
      INTEGER LOCRB, LOCIB
      PARAMETER (LOCRB=1, LOCIB=11)
C     .. Scalar Arguments ..
      REAL ERR, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, LENIW, LENW, N, NELT
C     .. Array Arguments ..
      REAL A(NELT), B(N), RWORK(LENW), X(N)
      INTEGER IA(NELT), IWORK(LENIW), JA(NELT)
C     .. Local Scalars ..
      INTEGER LOCATD, LOCATP, LOCATZ, LOCD, LOCDZ, LOCIW, LOCP, LOCR,
     +        LOCW, LOCZ
C     .. External Subroutines ..
      EXTERNAL SCGN, SCHKW, SS2Y, SSD2S, SSDI, SSMTV, SSMV
C***FIRST EXECUTABLE STATEMENT  SSDCGN
C
      IERR = 0
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
C
C         Modify the SLAP matrix data structure to YSMP-Column.
      CALL SS2Y( N, NELT, IA, JA, A, ISYM )
C
C         Set up the work arrays.
      LOCIW = LOCIB
C
      LOCD = LOCRB
      LOCR = LOCD + N
      LOCZ = LOCR + N
      LOCP = LOCZ + N
      LOCATP = LOCP + N
      LOCATZ = LOCATP + N
      LOCDZ = LOCATZ + N
      LOCATD = LOCDZ + N
      LOCW = LOCATD + N
C
C         Check the workspace allocations.
      CALL SCHKW( 'SSDCGN', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
      IF( IERR.NE.0 ) RETURN
C
      IWORK(4) = LOCD
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
C
C         Compute the inverse of the diagonal of AA'.  This will be
C         used as the preconditioner.
      CALL SSD2S(N, NELT, IA, JA, A, ISYM, RWORK(1))
C
C         Perform Conjugate Gradient algorithm on the normal equations.
      CALL SCGN( N, B, X, NELT, IA, JA, A, ISYM, SSMV, SSMTV, SSDI,
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK(LOCR),
     $     RWORK(LOCZ), RWORK(LOCP), RWORK(LOCATP), RWORK(LOCATZ),
     $     RWORK(LOCDZ), RWORK(LOCATD), RWORK, IWORK )
C
      IF( ITER.GT.ITMAX ) IERR = 2
      RETURN
C------------- LAST LINE OF SSDCGN FOLLOWS ----------------------------
      END
