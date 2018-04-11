      DOUBLE PRECISION FUNCTION DPRVEC (M, U, V)
C
      DOUBLE PRECISION DDOT
      INTEGER M, N, NP
      DOUBLE PRECISION U(*), V(*), VP
C***FIRST EXECUTABLE STATEMENT  DPRVEC
      N = M/2
      NP = N + 1
      VP = DDOT(N,U(1),1,V(NP),1)
      DPRVEC = DDOT(N,U(NP),1,V(1),1) - VP
      RETURN
      END
