      COMPLEX FUNCTION CATANH (Z)
      COMPLEX Z, CI, CATAN
      SAVE CI
      DATA CI /(0.,1.)/
C***FIRST EXECUTABLE STATEMENT  CATANH
      CATANH = -CI*CATAN(CI*Z)
C
      RETURN
      END
