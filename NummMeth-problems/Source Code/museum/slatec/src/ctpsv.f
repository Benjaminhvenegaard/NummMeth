      SUBROUTINE CTPSV (UPLO, TRANS, DIAG, N, AP, X, INCX)
C     .. Scalar Arguments ..
      INTEGER            INCX, N
      CHARACTER*1        DIAG, TRANS, UPLO
C     .. Array Arguments ..
      COMPLEX            AP( * ), X( * )
C     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
C     .. Local Scalars ..
      COMPLEX            TEMP
      INTEGER            I, INFO, IX, J, JX, K, KK, KX
      LOGICAL            NOCONJ, NOUNIT
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          CONJG
C***FIRST EXECUTABLE STATEMENT  CTPSV
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
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CTPSV ', INFO )
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
C     Start the operations. In this version the elements of AP are
C     accessed sequentially with one pass through AP.
C
      IF( LSAME( TRANS, 'N' ) )THEN
C
C        Form  x := inv( A )*x.
C
         IF( LSAME( UPLO, 'U' ) )THEN
            KK = ( N*( N + 1 ) )/2
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/AP( KK )
                     TEMP = X( J )
                     K    = KK     - 1
                     DO 10, I = J - 1, 1, -1
                        X( I ) = X( I ) - TEMP*AP( K )
                        K      = K      - 1
   10                CONTINUE
                  END IF
                  KK = KK - J
   20          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 40, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/AP( KK )
                     TEMP = X( JX )
                     IX   = JX
                     DO 30, K = KK - 1, KK - J + 1, -1
                        IX      = IX      - INCX
                        X( IX ) = X( IX ) - TEMP*AP( K )
   30                CONTINUE
                  END IF
                  JX = JX - INCX
                  KK = KK - J
   40          CONTINUE
            END IF
         ELSE
            KK = 1
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/AP( KK )
                     TEMP = X( J )
                     K    = KK     + 1
                     DO 50, I = J + 1, N
                        X( I ) = X( I ) - TEMP*AP( K )
                        K      = K      + 1
   50                CONTINUE
                  END IF
                  KK = KK + ( N - J + 1 )
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/AP( KK )
                     TEMP = X( JX )
                     IX   = JX
                     DO 70, K = KK + 1, KK + N - J
                        IX      = IX      + INCX
                        X( IX ) = X( IX ) - TEMP*AP( K )
   70                CONTINUE
                  END IF
                  JX = JX + INCX
                  KK = KK + ( N - J + 1 )
   80          CONTINUE
            END IF
         END IF
      ELSE
C
C        Form  x := inv( A' )*x  or  x := inv( conjg( A' ) )*x.
C
         IF( LSAME( UPLO, 'U' ) )THEN
            KK = 1
            IF( INCX.EQ.1 )THEN
               DO 110, J = 1, N
                  TEMP = X( J )
                  K    = KK
                  IF( NOCONJ )THEN
                     DO 90, I = 1, J - 1
                        TEMP = TEMP - AP( K )*X( I )
                        K    = K    + 1
   90                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/AP( KK + J - 1 )
                  ELSE
                     DO 100, I = 1, J - 1
                        TEMP = TEMP - CONJG( AP( K ) )*X( I )
                        K    = K    + 1
  100                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/CONJG( AP( KK + J - 1 ) )
                  END IF
                  X( J ) = TEMP
                  KK     = KK   + J
  110          CONTINUE
            ELSE
               JX = KX
               DO 140, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  IF( NOCONJ )THEN
                     DO 120, K = KK, KK + J - 2
                        TEMP = TEMP - AP( K )*X( IX )
                        IX   = IX   + INCX
  120                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/AP( KK + J - 1 )
                  ELSE
                     DO 130, K = KK, KK + J - 2
                        TEMP = TEMP - CONJG( AP( K ) )*X( IX )
                        IX   = IX   + INCX
  130                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/CONJG( AP( KK + J - 1 ) )
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   + INCX
                  KK      = KK   + J
  140          CONTINUE
            END IF
         ELSE
            KK = ( N*( N + 1 ) )/2
            IF( INCX.EQ.1 )THEN
               DO 170, J = N, 1, -1
                  TEMP = X( J )
                  K    = KK
                  IF( NOCONJ )THEN
                     DO 150, I = N, J + 1, -1
                        TEMP = TEMP - AP( K )*X( I )
                        K    = K    - 1
  150                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/AP( KK - N + J )
                  ELSE
                     DO 160, I = N, J + 1, -1
                        TEMP = TEMP - CONJG( AP( K ) )*X( I )
                        K    = K    - 1
  160                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/CONJG( AP( KK - N + J ) )
                  END IF
                  X( J ) = TEMP
                  KK     = KK   - ( N - J + 1 )
  170          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 200, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  IF( NOCONJ )THEN
                     DO 180, K = KK, KK - ( N - ( J + 1 ) ), -1
                        TEMP = TEMP - AP( K )*X( IX )
                        IX   = IX   - INCX
  180                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/AP( KK - N + J )
                  ELSE
                     DO 190, K = KK, KK - ( N - ( J + 1 ) ), -1
                        TEMP = TEMP - CONJG( AP( K ) )*X( IX )
                        IX   = IX   - INCX
  190                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/CONJG( AP( KK - N + J ) )
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   - INCX
                  KK      = KK   - ( N - J + 1 )
  200          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of CTPSV .
C
      END
