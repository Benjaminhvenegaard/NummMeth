      SUBROUTINE DQK31 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
      DOUBLE PRECISION A,ABSC,ABSERR,B,CENTR,DHLGTH,
     1  D1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC,
     2  RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(15),FV2(15),XGK(16),WGK(16),WG(8)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 31-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 15-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 15-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 31-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 15-POINT GAUSS RULE
C
C
C GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
C AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
C BELL LABS, NOV. 1981.
C
      SAVE WG, XGK, WGK
      DATA WG  (  1) / 0.0307532419 9611726835 4628393577 204 D0 /
      DATA WG  (  2) / 0.0703660474 8810812470 9267416450 667 D0 /
      DATA WG  (  3) / 0.1071592204 6717193501 1869546685 869 D0 /
      DATA WG  (  4) / 0.1395706779 2615431444 7804794511 028 D0 /
      DATA WG  (  5) / 0.1662692058 1699393355 3200860481 209 D0 /
      DATA WG  (  6) / 0.1861610000 1556221102 6800561866 423 D0 /
      DATA WG  (  7) / 0.1984314853 2711157645 6118326443 839 D0 /
      DATA WG  (  8) / 0.2025782419 2556127288 0620199967 519 D0 /
C
      DATA XGK (  1) / 0.9980022986 9339706028 5172840152 271 D0 /
      DATA XGK (  2) / 0.9879925180 2048542848 9565718586 613 D0 /
      DATA XGK (  3) / 0.9677390756 7913913425 7347978784 337 D0 /
      DATA XGK (  4) / 0.9372733924 0070590430 7758947710 209 D0 /
      DATA XGK (  5) / 0.8972645323 4408190088 2509656454 496 D0 /
      DATA XGK (  6) / 0.8482065834 1042721620 0648320774 217 D0 /
      DATA XGK (  7) / 0.7904185014 4246593296 7649294817 947 D0 /
      DATA XGK (  8) / 0.7244177313 6017004741 6186054613 938 D0 /
      DATA XGK (  9) / 0.6509967412 9741697053 3735895313 275 D0 /
      DATA XGK ( 10) / 0.5709721726 0853884753 7226737253 911 D0 /
      DATA XGK ( 11) / 0.4850818636 4023968069 3655740232 351 D0 /
      DATA XGK ( 12) / 0.3941513470 7756336989 7207370981 045 D0 /
      DATA XGK ( 13) / 0.2991800071 5316881216 6780024266 389 D0 /
      DATA XGK ( 14) / 0.2011940939 9743452230 0628303394 596 D0 /
      DATA XGK ( 15) / 0.1011420669 1871749902 7074231447 392 D0 /
      DATA XGK ( 16) / 0.0000000000 0000000000 0000000000 000 D0 /
C
      DATA WGK (  1) / 0.0053774798 7292334898 7792051430 128 D0 /
      DATA WGK (  2) / 0.0150079473 2931612253 8374763075 807 D0 /
      DATA WGK (  3) / 0.0254608473 2671532018 6874001019 653 D0 /
      DATA WGK (  4) / 0.0353463607 9137584622 2037948478 360 D0 /
      DATA WGK (  5) / 0.0445897513 2476487660 8227299373 280 D0 /
      DATA WGK (  6) / 0.0534815246 9092808726 5343147239 430 D0 /
      DATA WGK (  7) / 0.0620095678 0067064028 5139230960 803 D0 /
      DATA WGK (  8) / 0.0698541213 1872825870 9520077099 147 D0 /
      DATA WGK (  9) / 0.0768496807 5772037889 4432777482 659 D0 /
      DATA WGK ( 10) / 0.0830805028 2313302103 8289247286 104 D0 /
      DATA WGK ( 11) / 0.0885644430 5621177064 7275443693 774 D0 /
      DATA WGK ( 12) / 0.0931265981 7082532122 5486872747 346 D0 /
      DATA WGK ( 13) / 0.0966427269 8362367850 5179907627 589 D0 /
      DATA WGK ( 14) / 0.0991735987 2179195933 2393173484 603 D0 /
      DATA WGK ( 15) / 0.1007698455 2387559504 4946662617 570 D0 /
      DATA WGK ( 16) / 0.1013300070 1479154901 7374792767 493 D0 /
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 15-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 31-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C***FIRST EXECUTABLE STATEMENT  DQK31
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
C
      CENTR = 0.5D+00*(A+B)
      HLGTH = 0.5D+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 31-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      FC = F(CENTR)
      RESG = WG(8)*FC
      RESK = WGK(16)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,7
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,8
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5D+00
      RESASC = WGK(16)*ABS(FC-RESKH)
      DO 20 J=1,15
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00)
     1  ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
      IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX
     1  ((EPMACH*0.5D+02)*RESABS,ABSERR)
      RETURN
      END
