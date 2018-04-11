      SUBROUTINE DBKIAS (X, N, KTRMS, T, ANS, IND, MS, GMRN, H, IERR)
      INTEGER I, II, IND, J, JMI, JN, K, KK, KM, KTRMS, MM, MP, MS, N,
     * IERR
      DOUBLE PRECISION ANS, B, BND, DEN1, DEN2, DEN3, ER, ERR, FJ, FK,
     * FLN, FM1, GMRN, G1, GS, H, HN, HRTPI, RAT, RG1, RXP, RZ, RZX, S,
     * SS, SUMI, SUMJ, T, TOL, V, W, X, XP, Z
      DOUBLE PRECISION DGAMRN, D1MACH
      DIMENSION B(120), XP(16), S(31), H(*), V(52), W(52), T(50),
     * BND(15)
      SAVE B, BND, HRTPI
C-----------------------------------------------------------------------
C             COEFFICIENTS OF POLYNOMIAL P(J-1,X), J=1,15
C-----------------------------------------------------------------------
      DATA B(1), B(2), B(3), B(4), B(5), B(6), B(7), B(8), B(9), B(10),
     * B(11), B(12), B(13), B(14), B(15), B(16), B(17), B(18), B(19),
     * B(20), B(21), B(22), B(23), B(24) /1.00000000000000000D+00,
     * 1.00000000000000000D+00,-2.00000000000000000D+00,
     * 1.00000000000000000D+00,-8.00000000000000000D+00,
     * 6.00000000000000000D+00,1.00000000000000000D+00,
     * -2.20000000000000000D+01,5.80000000000000000D+01,
     * -2.40000000000000000D+01,1.00000000000000000D+00,
     * -5.20000000000000000D+01,3.28000000000000000D+02,
     * -4.44000000000000000D+02,1.20000000000000000D+02,
     * 1.00000000000000000D+00,-1.14000000000000000D+02,
     * 1.45200000000000000D+03,-4.40000000000000000D+03,
     * 3.70800000000000000D+03,-7.20000000000000000D+02,
     * 1.00000000000000000D+00,-2.40000000000000000D+02,
     * 5.61000000000000000D+03/
      DATA B(25), B(26), B(27), B(28), B(29), B(30), B(31), B(32),
     * B(33), B(34), B(35), B(36), B(37), B(38), B(39), B(40), B(41),
     * B(42), B(43), B(44), B(45), B(46), B(47), B(48)
     * /-3.21200000000000000D+04,5.81400000000000000D+04,
     * -3.39840000000000000D+04,5.04000000000000000D+03,
     * 1.00000000000000000D+00,-4.94000000000000000D+02,
     * 1.99500000000000000D+04,-1.95800000000000000D+05,
     * 6.44020000000000000D+05,-7.85304000000000000D+05,
     * 3.41136000000000000D+05,-4.03200000000000000D+04,
     * 1.00000000000000000D+00,-1.00400000000000000D+03,
     * 6.72600000000000000D+04,-1.06250000000000000D+06,
     * 5.76550000000000000D+06,-1.24400640000000000D+07,
     * 1.10262960000000000D+07,-3.73392000000000000D+06,
     * 3.62880000000000000D+05,1.00000000000000000D+00,
     * -2.02600000000000000D+03,2.18848000000000000D+05/
      DATA B(49), B(50), B(51), B(52), B(53), B(54), B(55), B(56),
     * B(57), B(58), B(59), B(60), B(61), B(62), B(63), B(64), B(65),
     * B(66), B(67), B(68), B(69), B(70), B(71), B(72)
     * /-5.32616000000000000D+06,4.47650000000000000D+07,
     * -1.55357384000000000D+08,2.38904904000000000D+08,
     * -1.62186912000000000D+08,4.43390400000000000D+07,
     * -3.62880000000000000D+06,1.00000000000000000D+00,
     * -4.07200000000000000D+03,6.95038000000000000D+05,
     * -2.52439040000000000D+07,3.14369720000000000D+08,
     * -1.64838430400000000D+09,4.00269508800000000D+09,
     * -4.64216395200000000D+09,2.50748121600000000D+09,
     * -5.68356480000000000D+08,3.99168000000000000D+07,
     * 1.00000000000000000D+00,-8.16600000000000000D+03,
     * 2.17062600000000000D+06,-1.14876376000000000D+08,
     * 2.05148277600000000D+09,-1.55489607840000000D+10/
      DATA B(73), B(74), B(75), B(76), B(77), B(78), B(79), B(80),
     * B(81), B(82), B(83), B(84), B(85), B(86), B(87), B(88), B(89),
     * B(90), B(91), B(92), B(93), B(94), B(95), B(96)
     * /5.60413987840000000D+10,-1.01180433024000000D+11,
     * 9.21997902240000000D+10,-4.07883018240000000D+10,
     * 7.82771904000000000D+09,-4.79001600000000000D+08,
     * 1.00000000000000000D+00,-1.63560000000000000D+04,
     * 6.69969600000000000D+06,-5.07259276000000000D+08,
     * 1.26698177760000000D+10,-1.34323420224000000D+11,
     * 6.87720046384000000D+11,-1.81818864230400000D+12,
     * 2.54986547342400000D+12,-1.88307966182400000D+12,
     * 6.97929436800000000D+11,-1.15336085760000000D+11,
     * 6.22702080000000000D+09,1.00000000000000000D+00,
     * -3.27380000000000000D+04,2.05079880000000000D+07,
     * -2.18982980800000000D+09,7.50160522280000000D+10/
      DATA B(97), B(98), B(99), B(100), B(101), B(102), B(103), B(104),
     * B(105), B(106), B(107), B(108), B(109), B(110), B(111), B(112),
     * B(113), B(114), B(115), B(116), B(117), B(118)
     * /-1.08467651241600000D+12,7.63483214939200000D+12,
     * -2.82999100661120000D+13,5.74943734645920000D+13,
     * -6.47283751398720000D+13,3.96895780558080000D+13,
     * -1.25509040179200000D+13,1.81099255680000000D+12,
     * -8.71782912000000000D+10,1.00000000000000000D+00,
     * -6.55040000000000000D+04,6.24078900000000000D+07,
     * -9.29252692000000000D+09,4.29826006340000000D+11,
     * -8.30844432796800000D+12,7.83913848313120000D+13,
     * -3.94365587815520000D+14,1.11174747256968000D+15,
     * -1.79717122069056000D+15,1.66642448627145600D+15,
     * -8.65023253219584000D+14,2.36908271543040000D+14/
      DATA B(119), B(120) /-3.01963769856000000D+13,
     * 1.30767436800000000D+12/
C-----------------------------------------------------------------------
C             BOUNDS B(M,K) , K=M-3
C-----------------------------------------------------------------------
      DATA BND(1), BND(2), BND(3), BND(4), BND(5), BND(6), BND(7),
     * BND(8), BND(9), BND(10), BND(11), BND(12), BND(13), BND(14),
     * BND(15) /1.0D0,1.0D0,1.0D0,1.0D0,3.10D0,5.18D0,11.7D0,29.8D0,
     * 90.4D0,297.0D0,1070.0D0,4290.0D0,18100.0D0,84700.0D0,408000.0D0/
      DATA HRTPI /8.86226925452758014D-01/
C
C***FIRST EXECUTABLE STATEMENT  DBKIAS
      IERR=0
      TOL = MAX(D1MACH(4),1.0D-18)
      FLN = N
      RZ = 1.0D0/(X+FLN)
      RZX = X*RZ
      Z = 0.5D0*(X+FLN)
      IF (IND.GT.1) GO TO 10
      GMRN = DGAMRN(Z)
   10 CONTINUE
      GS = HRTPI*GMRN
      G1 = GS + GS
      RG1 = 1.0D0/G1
      GMRN = (RZ+RZ)/GMRN
      IF (IND.GT.1) GO TO 70
C-----------------------------------------------------------------------
C     EVALUATE ERROR FOR M=MS
C-----------------------------------------------------------------------
      HN = 0.5D0*FLN
      DEN2 = KTRMS + KTRMS + N
      DEN3 = DEN2 - 2.0D0
      DEN1 = X + DEN2
      ERR = RG1*(X+X)/(DEN1-1.0D0)
      IF (N.EQ.0) GO TO 20
      RAT = 1.0D0/(FLN*FLN)
   20 CONTINUE
      IF (KTRMS.EQ.0) GO TO 30
      FJ = KTRMS
      RAT = 0.25D0/(HRTPI*DEN3*SQRT(FJ))
   30 CONTINUE
      ERR = ERR*RAT
      FJ = -3.0D0
      DO 50 J=1,15
        IF (J.LE.5) ERR = ERR/DEN1
        FM1 = MAX(1.0D0,FJ)
        FJ = FJ + 1.0D0
        ER = BND(J)*ERR
        IF (KTRMS.EQ.0) GO TO 40
        ER = ER/FM1
        IF (ER.LT.TOL) GO TO 60
        IF (J.GE.5) ERR = ERR/DEN3
        GO TO 50
   40   CONTINUE
        ER = ER*(1.0D0+HN/FM1)
        IF (ER.LT.TOL) GO TO 60
        IF (J.GE.5) ERR = ERR/FLN
   50 CONTINUE
      GO TO 200
   60 CONTINUE
      MS = J
   70 CONTINUE
      MM = MS + MS
      MP = MM + 1
C-----------------------------------------------------------------------
C     H(K)=(-Z)**(K)*(PSI(K-1,Z)-PSI(K-1,Z+0.5))/GAMMA(K) , K=1,2,...,MM
C-----------------------------------------------------------------------
      IF (IND.GT.1) GO TO 80
      CALL DHKSEQ(Z, MM, H, IERR)
      GO TO 100
   80 CONTINUE
      RAT = Z/(Z-0.5D0)
      RXP = RAT
      DO 90 I=1,MM
        H(I) = RXP*(1.0D0-H(I))
        RXP = RXP*RAT
   90 CONTINUE
  100 CONTINUE
C-----------------------------------------------------------------------
C     SCALED S SEQUENCE
C-----------------------------------------------------------------------
      S(1) = 1.0D0
      FK = 1.0D0
      DO 120 K=2,MP
        SS = 0.0D0
        KM = K - 1
        I = KM
        DO 110 J=1,KM
          SS = SS + S(J)*H(I)
          I = I - 1
  110   CONTINUE
        S(K) = SS/FK
        FK = FK + 1.0D0
  120 CONTINUE
C-----------------------------------------------------------------------
C     SCALED S-TILDA SEQUENCE
C-----------------------------------------------------------------------
      IF (KTRMS.EQ.0) GO TO 160
      FK = 0.0D0
      SS = 0.0D0
      RG1 = RG1/Z
      DO 130 K=1,KTRMS
        V(K) = Z/(Z+FK)
        W(K) = T(K)*V(K)
        SS = SS + W(K)
        FK = FK + 1.0D0
  130 CONTINUE
      S(1) = S(1) - SS*RG1
      DO 150 I=2,MP
        SS = 0.0D0
        DO 140 K=1,KTRMS
          W(K) = W(K)*V(K)
          SS = SS + W(K)
  140   CONTINUE
        S(I) = S(I) - SS*RG1
  150 CONTINUE
  160 CONTINUE
C-----------------------------------------------------------------------
C     SUM ON J
C-----------------------------------------------------------------------
      SUMJ = 0.0D0
      JN = 1
      RXP = 1.0D0
      XP(1) = 1.0D0
      DO 190 J=1,MS
        JN = JN + J - 1
        XP(J+1) = XP(J)*RZX
        RXP = RXP*RZ
C-----------------------------------------------------------------------
C     SUM ON I
C-----------------------------------------------------------------------
        SUMI = 0.0D0
        II = JN
        DO 180 I=1,J
          JMI = J - I + 1
          KK = J + I + 1
          DO 170 K=1,JMI
            V(K) = S(KK)*XP(K)
            KK = KK + 1
  170     CONTINUE
          CALL DBDIFF(JMI, V)
          SUMI = SUMI + B(II)*V(JMI)*XP(I+1)
          II = II + 1
  180   CONTINUE
        SUMJ = SUMJ + SUMI*RXP
  190 CONTINUE
      ANS = GS*(S(1)-SUMJ)
      RETURN
  200 CONTINUE
      IERR=2
      RETURN
      END
