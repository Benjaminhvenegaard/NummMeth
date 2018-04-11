DECLARE SUB AVEVAR (DATQ!(), N!, AVE!, VAR!)
DECLARE FUNCTION BETAI! (A!, B!, X!)

SUB TPTEST (DATA1(), DATA2(), N, T, PROB)
CALL AVEVAR(DATA1(), N, AVE1, VAR1)
CALL AVEVAR(DATA2(), N, AVE2, VAR2)
COV = 0!
FOR J = 1 TO N
  COV = COV + (DATA1(J) - AVE1) * (DATA2(J) - AVE2)
NEXT J
DF = N - 1
COV = COV / DF
SD = SQR((VAR1 + VAR2 - 2! * COV) / CSNG(N))
T = (AVE1 - AVE2) / SD
PROB = BETAI(.5 * DF, .5, DF / (DF + T ^ 2))
END SUB
