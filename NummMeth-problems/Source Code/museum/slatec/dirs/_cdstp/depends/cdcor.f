      SUBROUTINE CDCOR (DFDY, EL, FA, H, IERROR, IMPL, IPVT, MATDIM,
     8   MITER, ML, MU, N, NDE, NQ, T, USERS, Y, YH, YWT, EVALFA, SAVE1,
     8   SAVE2, A, D, JSTATE)
      INTEGER I, IERROR, IFLAG, IMPL, J, JSTATE, MATDIM, MITER, ML, MU,
     8        MW, N, NDE, NQ
      COMPLEX A(MATDIM,*), DFDY(MATDIM,*), SAVE1(*), SAVE2(*), Y(*),
     8        YH(N,*), YWT(*)
      REAL D, EL(13,12), H, SCNRM2, T
      INTEGER IPVT(*)
      LOGICAL EVALFA
C***FIRST EXECUTABLE STATEMENT  CDCOR
      IF (MITER .EQ. 0) THEN
        IF (IERROR .EQ. 1 .OR. IERROR .EQ. 5) THEN
          DO 100 I = 1,N
 100        SAVE1(I) = (H*SAVE2(I) - YH(I,2) - SAVE1(I))/YWT(I)
        ELSE
          DO 102 I = 1,N
            SAVE1(I) = (H*SAVE2(I) - YH(I,2) - SAVE1(I))/
     8      MAX(ABS(Y(I)), ABS(YWT(I)))
 102        CONTINUE
        END IF
        D = SCNRM2(N, SAVE1, 1)/SQRT(REAL(N))
        DO 105 I = 1,N
 105      SAVE1(I) = H*SAVE2(I) - YH(I,2)
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IF (IMPL .EQ. 0) THEN
          DO 130 I = 1,N
 130        SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)
        ELSE IF (IMPL .EQ. 1) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 150 I = 1,N
 150        SAVE2(I) = H*SAVE2(I)
          DO 160 J = 1,N
            DO 160 I = 1,N
 160          SAVE2(I) = SAVE2(I) - A(I,J)*(YH(J,2) + SAVE1(J))
        ELSE IF (IMPL .EQ. 2) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 180 I = 1,N
 180        SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))
        ELSE IF (IMPL .EQ. 3) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 140 I = 1,N
 140        SAVE2(I) = H*SAVE2(I)
          DO 170 J = 1,NDE
            DO 170 I = 1,NDE
 170          SAVE2(I) = SAVE2(I) - A(I,J)*(YH(J,2) + SAVE1(J))
        END IF
        CALL CGESL (DFDY, MATDIM, N, IPVT, SAVE2, 0)
        IF (IERROR .EQ. 1 .OR. IERROR .EQ. 5) THEN
          DO 200 I = 1,N
            SAVE1(I) = SAVE1(I) + SAVE2(I)
 200        SAVE2(I) = SAVE2(I)/YWT(I)
        ELSE
          DO 205 I = 1,N
            SAVE1(I) = SAVE1(I) + SAVE2(I)
 205        SAVE2(I) = SAVE2(I)/MAX(ABS(Y(I)), ABS(YWT(I)))
        END IF
        D = SCNRM2(N, SAVE2, 1)/SQRT(REAL(N))
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        IF (IMPL .EQ. 0) THEN
          DO 230 I = 1,N
 230        SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)
        ELSE IF (IMPL .EQ. 1) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 250 I = 1,N
 250        SAVE2(I) = H*SAVE2(I)
          MW = ML + 1 + MU
          DO 260 J = 1,N
            DO 260 I = MAX(ML+1, MW+1-J), MIN(MW+N-J, MW+ML)
              SAVE2(I+J-MW) = SAVE2(I+J-MW)
     8                        - A(I,J)*(YH(J,2) + SAVE1(J))
 260        CONTINUE
        ELSE IF (IMPL .EQ. 2) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 280 I = 1,N
 280        SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))
        ELSE IF (IMPL .EQ. 3) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 270 I = 1,N
 270        SAVE2(I) = H*SAVE2(I)
          MW = ML + 1 + MU
          DO 290 J = 1,NDE
            DO 290 I = MAX(ML+1, MW+1-J), MIN(MW+NDE-J, MW+ML)
              SAVE2(I+J-MW) = SAVE2(I+J-MW)
     8                        - A(I,J)*(YH(J,2) + SAVE1(J))
 290        CONTINUE
        END IF
        CALL CGBSL (DFDY, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
        IF (IERROR .EQ. 1 .OR. IERROR .EQ. 5) THEN
          DO 300 I = 1,N
            SAVE1(I) = SAVE1(I) + SAVE2(I)
 300        SAVE2(I) = SAVE2(I)/YWT(I)
        ELSE
          DO 305 I = 1,N
            SAVE1(I) = SAVE1(I) + SAVE2(I)
 305        SAVE2(I) = SAVE2(I)/MAX(ABS(Y(I)), ABS(YWT(I)))
        END IF
        D = SCNRM2(N, SAVE2, 1)/SQRT(REAL(N))
      ELSE IF (MITER .EQ. 3) THEN
        IFLAG = 2
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)
        IF (N .EQ. 0) THEN
          JSTATE = 10
          RETURN
        END IF
        IF (IERROR .EQ. 1 .OR. IERROR .EQ. 5) THEN
          DO 320 I = 1,N
            SAVE1(I) = SAVE1(I) + SAVE2(I)
 320        SAVE2(I) = SAVE2(I)/YWT(I)
        ELSE
          DO 325 I = 1,N
            SAVE1(I) = SAVE1(I) + SAVE2(I)
 325        SAVE2(I) = SAVE2(I)/MAX(ABS(Y(I)), ABS(YWT(I)))
        END IF
        D = SCNRM2(N, SAVE2, 1)/SQRT(REAL(N))
      END IF
      RETURN
      END
