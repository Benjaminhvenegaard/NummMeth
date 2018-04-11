      SUBROUTINE CHPR2 (UPLO, N, ALPHA, X, INCX, Y, INCY, AP)
C     .. Scalar Arguments ..
      COMPLEX            ALPHA
      INTEGER            INCX, INCY, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      COMPLEX            AP( * ), X( * ), Y( * )
C     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
C     .. Local Scalars ..
      COMPLEX            TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY
C     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
C     .. External Subroutines ..
      EXTERNAL           XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC          CONJG, REAL
C***FIRST EXECUTABLE STATEMENT  CHPR2
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CHPR2 ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Set up the start points in X and Y if the increments are not both
C     unity.
C
      IF( ( INCX.NE.1 ).OR.( INCY.NE.1 ) )THEN
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( N - 1 )*INCY
         END IF
         JX = KX
         JY = KY
      END IF
C
C     Start the operations. In this version the elements of the array AP
C     are accessed sequentially with one pass through AP.
C
      KK = 1
      IF( LSAME( UPLO, 'U' ) )THEN
C
C        Form  A  when upper triangle is stored in AP.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 20, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*CONJG( Y( J ) )
                  TEMP2 = CONJG( ALPHA*X( J ) )
                  K     = KK
                  DO 10, I = 1, J - 1
                     AP( K ) = AP( K ) + X( I )*TEMP1 + Y( I )*TEMP2
                     K       = K       + 1
   10             CONTINUE
                  AP( KK + J - 1 ) = REAL( AP( KK + J - 1 ) ) +
     $                               REAL( X( J )*TEMP1 + Y( J )*TEMP2 )
               ELSE
                  AP( KK + J - 1 ) = REAL( AP( KK + J - 1 ) )
               END IF
               KK = KK + J
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*CONJG( Y( JY ) )
                  TEMP2 = CONJG( ALPHA*X( JX ) )
                  IX    = KX
                  IY    = KY
                  DO 30, K = KK, KK + J - 2
                     AP( K ) = AP( K ) + X( IX )*TEMP1 + Y( IY )*TEMP2
                     IX      = IX      + INCX
                     IY      = IY      + INCY
   30             CONTINUE
                  AP( KK + J - 1 ) = REAL( AP( KK + J - 1 ) ) +
     $                               REAL( X( JX )*TEMP1 +
     $                                     Y( JY )*TEMP2 )
               ELSE
                  AP( KK + J - 1 ) = REAL( AP( KK + J - 1 ) )
               END IF
               JX = JX + INCX
               JY = JY + INCY
               KK = KK + J
   40       CONTINUE
         END IF
      ELSE
C
C        Form  A  when lower triangle is stored in AP.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1   = ALPHA*CONJG( Y( J ) )
                  TEMP2   = CONJG( ALPHA*X( J ) )
                  AP( KK ) = REAL( AP( KK ) ) +
     $                       REAL( X( J )*TEMP1 + Y( J )*TEMP2 )
                  K        = KK               + 1
                  DO 50, I = J + 1, N
                     AP( K ) = AP( K ) + X( I )*TEMP1 + Y( I )*TEMP2
                     K       = K       + 1
   50             CONTINUE
               ELSE
                  AP( KK ) = REAL( AP( KK ) )
               END IF
               KK = KK + N - J + 1
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1    = ALPHA*CONJG( Y( JY ) )
                  TEMP2    = CONJG( ALPHA*X( JX ) )
                  AP( KK ) = REAL( AP( KK ) ) +
     $                       REAL( X( JX )*TEMP1 + Y( JY )*TEMP2 )
                  IX       = JX
                  IY       = JY
                  DO 70, K = KK + 1, KK + N - J
                     IX      = IX      + INCX
                     IY      = IY      + INCY
                     AP( K ) = AP( K ) + X( IX )*TEMP1 + Y( IY )*TEMP2
   70             CONTINUE
               ELSE
                  AP( KK ) = REAL( AP( KK ) )
               END IF
               JX = JX + INCX
               JY = JY + INCY
               KK = KK + N - J + 1
   80       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of CHPR2 .
C
      END
