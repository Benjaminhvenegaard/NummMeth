      DOUBLE PRECISION FUNCTION DHVNRM (V, NCOMP)
C
      INTEGER K, NCOMP
      DOUBLE PRECISION V
      DIMENSION V(*)
C***FIRST EXECUTABLE STATEMENT  DHVNRM
      DHVNRM = 0.0D0
      DO 10 K = 1, NCOMP
         DHVNRM = MAX(DHVNRM,ABS(V(K)))
   10 CONTINUE
      RETURN
      END
