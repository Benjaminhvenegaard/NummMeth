      SUBROUTINE CTBSV (UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX)
C     .. Scalar Arguments ..
      INTEGER            INCX, K, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
C     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * )
C     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
C     .. Local Scalars ..
      COMPLEX            TEMP
      INTEGER            I, INFO, IX, J, JX, KPLUS1, KX, L
      LOGICAL            NOCONJ, NOUNIT
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, MIN
C***FIRST EXECUTABLE STATEMENT  CTBSV
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( K.LT.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.( K + 1 ) )THEN
         INFO = 7
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CTBSV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      NOCONJ = LSAME( TRANS, 'T' )
      NOUNIT = LSAME( DIAG , 'N' )
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed by sequentially with one pass through A.
C
      IF( LSAME( TRANS, 'N' ) )THEN
C
C        Form  x := inv( A )*x.
C
         IF( LSAME( UPLO, 'U' ) )THEN
            KPLUS1 = K + 1
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     L = KPLUS1 - J
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( KPLUS1, J )
                     TEMP = X( J )
                     DO 10, I = J - 1, MAX( 1, J - K ), -1
                        X( I ) = X( I ) - TEMP*A( L + I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 40, J = N, 1, -1
                  KX = KX - INCX
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     L  = KPLUS1 - J
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( KPLUS1, J )
                     TEMP = X( JX )
                     DO 30, I = J - 1, MAX( 1, J - K ), -1
                        X( IX ) = X( IX ) - TEMP*A( L + I, J )
                        IX      = IX      - INCX
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     L = 1 - J
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( 1, J )
                     TEMP = X( J )
                     DO 50, I = J + 1, MIN( N, J + K )
                        X( I ) = X( I ) - TEMP*A( L + I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  KX = KX + INCX
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     L  = 1  - J
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( 1, J )
                     TEMP = X( JX )
                     DO 70, I = J + 1, MIN( N, J + K )
                        X( IX ) = X( IX ) - TEMP*A( L + I, J )
                        IX      = IX      + INCX
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
C
C        Form  x := inv( A' )*x  or  x := inv( conjg( A') )*x.
C
         IF( LSAME( UPLO, 'U' ) )THEN
            KPLUS1 = K + 1
            IF( INCX.EQ.1 )THEN
               DO 110, J = 1, N
                  TEMP = X( J )
                  L    = KPLUS1 - J
                  IF( NOCONJ )THEN
                     DO 90, I = MAX( 1, J - K ), J - 1
                        TEMP = TEMP - A( L + I, J )*X( I )
   90                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( KPLUS1, J )
                  ELSE
                     DO 100, I = MAX( 1, J - K ), J - 1
                        TEMP = TEMP - CONJG( A( L + I, J ) )*X( I )
  100                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/CONJG( A( KPLUS1, J ) )
                  END IF
                  X( J ) = TEMP
  110          CONTINUE
            ELSE
               JX = KX
               DO 140, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  L    = KPLUS1  - J
                  IF( NOCONJ )THEN
                     DO 120, I = MAX( 1, J - K ), J - 1
                        TEMP = TEMP - A( L + I, J )*X( IX )
                        IX   = IX   + INCX
  120                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( KPLUS1, J )
                  ELSE
                     DO 130, I = MAX( 1, J - K ), J - 1
                        TEMP = TEMP - CONJG( A( L + I, J ) )*X( IX )
                        IX   = IX   + INCX
  130                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/CONJG( A( KPLUS1, J ) )
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   + INCX
                  IF( J.GT.K )
     $               KX = KX + INCX
  140          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 170, J = N, 1, -1
                  TEMP = X( J )
                  L    = 1      - J
                  IF( NOCONJ )THEN
                     DO 150, I = MIN( N, J + K ), J + 1, -1
                        TEMP = TEMP - A( L + I, J )*X( I )
  150                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( 1, J )
                  ELSE
                     DO 160, I = MIN( N, J + K ), J + 1, -1
                        TEMP = TEMP - CONJG( A( L + I, J ) )*X( I )
  160                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/CONJG( A( 1, J ) )
                  END IF
                  X( J ) = TEMP
  170          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 200, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  L    = 1       - J
                  IF( NOCONJ )THEN
                     DO 180, I = MIN( N, J + K ), J + 1, -1
                        TEMP = TEMP - A( L + I, J )*X( IX )
                        IX   = IX   - INCX
  180                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( 1, J )
                  ELSE
                     DO 190, I = MIN( N, J + K ), J + 1, -1
                        TEMP = TEMP - CONJG( A( L + I, J ) )*X( IX )
                        IX   = IX   - INCX
  190                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/CONJG( A( 1, J ) )
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   - INCX
                  IF( ( N - J ).GE.K )
     $               KX = KX - INCX
  200          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of CTBSV .
C
      END
