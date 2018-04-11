DECLARE SUB AVEVAR (DATQ!(), N!, AVE!, VAR!)
DECLARE FUNCTION BETAI! (A!, B!, X!)

SUB TUTEST (DATA1(), N1, DATA2(), N2, T, PROB)
CALL AVEVAR(DATA1(), N1, AVE1, VAR1)
CALL AVEVAR(DATA2(), N2, AVE2, VAR2)
T = (AVE1 - AVE2) / SQR(VAR1 / N1 + VAR2 / N2)
DUM = (VAR1 / N1 + VAR2 / N2) ^ 2
DF = DUM / ((VAR1 / N1) ^ 2 / (N1 - 1) + (VAR2 / N2) ^ 2 / (N2 - 1))
PROB = BETAI(.5 * DF, .5, DF / (DF + T ^ 2))
END SUB

