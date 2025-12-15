PROGRAM KELOMPOK_4
      IMPLICIT NONE

      REAL :: S_IN, T_IN, P_IN, LAT_IN, PR_IN, CND_IN
      INTEGER :: M_IN
      REAL :: V_OUT_CND, V_OUT_SVAN, V_OUT_DEPTH, V_OUT_TF
      REAL :: V_OUT_CPSW, V_OUT_ATG, V_OUT_THETA, V_OUT_SVEL
      REAL :: SIGMA_RES
      
      REAL, EXTERNAL :: SAL78, SVAN, DEPTH, TF, CPSW, ATG, THETA, SVEL

      PRINT *, "=== KALKULATOR OSEANOGRAFI (UNESCO 44) ==="
      PRINT *, ""

      ! 1. PILIH MODE INPUT
      PRINT *, "PILIH MODE INPUT:"
      PRINT *, " 0 = Input via KONDUKTIVITAS (Hitung Salinitas)"
      PRINT *, " 1 = Input via SALINITAS (Hitung Konduktivitas)"
      READ *, M_IN
      PRINT *, ""

      ! 2. INPUT DATA UMUM
      PRINT *, "INPUT: Temperatur (C):"
      READ *, T_IN
      PRINT *, "INPUT: Tekanan (Decibars):"
      READ *, P_IN

      ! 3. LOGIKA SALINITAS/KONDUKTIVITAS
      IF (M_IN .EQ. 0) THEN
          PRINT *, "INPUT: Rasio Konduktivitas (R):"
          READ *, CND_IN
          ! Hitung Salinitas dari Konduktivitas
          S_IN = SAL78(CND_IN, T_IN, P_IN, 0)
          PRINT *, ">> HASIL SALINITAS (PSS-78)   : ", S_IN
      ELSE
          PRINT *, "INPUT: Salinitas (PSS-78):"
          READ *, S_IN
          ! Hitung Konduktivitas dari Salinitas
          V_OUT_CND = SAL78(S_IN, T_IN, P_IN, 1)
          PRINT *, ">> HASIL KONDUKTIVITAS (R)    : ", V_OUT_CND
      END IF
      PRINT *, ""

      ! 4. HITUNG PARAMETER LAIN (MENGGUNAKAN S_IN)
      PRINT *, "--- PARAMETER TURUNAN ---"
      
      V_OUT_TF = TF(S_IN, P_IN)
      PRINT *, "1. Freezing Point (TF)        : ", V_OUT_TF

      V_OUT_SVEL = SVEL(S_IN, T_IN, P_IN)
      PRINT *, "2. Sound Speed (SVEL)         : ", V_OUT_SVEL

      V_OUT_CPSW = CPSW(S_IN, T_IN, P_IN)
      PRINT *, "3. Specific Heat (CPSW)       : ", V_OUT_CPSW

      V_OUT_ATG = ATG(S_IN, T_IN, P_IN)
      PRINT *, "4. Adiabatic Lapse Rate (ATG) : ", V_OUT_ATG

      V_OUT_SVAN = SVAN(S_IN, T_IN, P_IN, SIGMA_RES)
      PRINT *, "5. Spec. Vol. Anomaly (SVAN)  : ", V_OUT_SVAN
      PRINT *, "   Density Anomaly (SIGMA)    : ", SIGMA_RES
      PRINT *, ""

      ! 5. PARAMETER DENGAN INPUT TAMBAHAN
      PRINT *, "--- CALCULATIONS REQUIRING EXTRA INPUT ---"
      
      PRINT *, "INPUT: Lintang (Degrees) for DEPTH:"
      READ *, LAT_IN
      V_OUT_DEPTH = DEPTH(P_IN, LAT_IN)
      PRINT *, ">> Depth (M)                  : ", V_OUT_DEPTH
      PRINT *, ""

      PRINT *, "INPUT: Tekanan Referensi (db) for THETA:"
      READ *, PR_IN
      V_OUT_THETA = THETA(S_IN, T_IN, P_IN, PR_IN)
      PRINT *, ">> Potential Temp (THETA)     : ", V_OUT_THETA
      PRINT *, ""

      PRINT *, "Selesai. Tekan Enter untuk keluar..."
      READ *, S_IN
      END PROGRAM KELOMPOK_4

! ------------------------------------------------------------------
! 1. KONVERSI KONDUKTIVITAS <-> SALINITAS (SAL78)
! ------------------------------------------------------------------
      REAL FUNCTION SAL78(CND,T,P,M)

      SAL(XR,XT) = ((((2.7081*XR-7.0261)*XR+14.0941)*XR+25.3851)*XR     &
     & -0.1692)*XR + 0.0080                                             &
     & + (XT/(1.0+0.0162*XT))*(((((-0.0144*XR+                          &
     & 0.0636)*XR-0.0375)*XR-0.0066)*XR-0.0056)*XR+0.0005)

      DSAL(XR,XT) = ((((13.5405*XR-28.1044)*XR+42.2823)*XR+50.7702)*XR  &
     & -0.1692) + (XT/(1.0+0.0162*XT))*((((-0.072*XR+0.2544)*XR         &
     & -0.1125)*XR-0.0132)*XR-0.0056)

      RT35(XT) = (((1.0031E-9*XT-6.9698E-7)*XT+1.104259E-4)*XT          &
     & + 2.00564E-2)*XT + 0.6766097

      C_POLY(XP) = ((3.989E-15*XP-6.370E-10)*XP+2.070E-5)*XP
      B_POLY(XT) = (4.464E-4*XT+3.426E-2)*XT + 1.0
      A_POLY(XT) = -3.107E-3*XT + 0.4215

      SAL78 = 0.0
      IF((M.EQ.0).AND.(CND.LE.5E-4)) RETURN
      IF((M.EQ.1).AND.(CND.LE.0.02)) RETURN

      DT = T - 15.0
      IF(M.EQ.1) GO TO 10

      R = CND
      RT = R/(RT35(T)*(1.0 + C_POLY(P)/(B_POLY(T) + A_POLY(T)*R)))
      RT = SQRT(ABS(RT))
      SAL78 = SAL(RT, DT)
      RETURN

   10 RT = SQRT(CND/35.0)
      SI = SAL(RT,DT)
      N = 0

   15 RT = RT + (CND - SI)/DSAL(RT,DT)
      SI = SAL(RT,DT)
      N = N + 1
      DELS = ABS(SI - CND)
      IF((DELS.GT.1.0E-4).AND.(N.LT.10)) GO TO 15

      RTT = RT35(T)*RT*RT
      AT = A_POLY(T)
      BT = B_POLY(T)
      CP = C_POLY(P)
      CP = RTT*(CP+BT)
      BT = BT - RTT*AT

      R = (SQRT(ABS(BT*BT + 4.0*AT*CP)) - BT)/(2.0*AT)
      SAL78 = R
      RETURN
      END

! ------------------------------------------------------------------
! 3. ANOMALI VOLUME SPESIFIK & DENSITAS (SVAN)
! ------------------------------------------------------------------
      REAL FUNCTION SVAN(S,T,P0,SIGMA)
      
      REAL P,T,S,SIG,SR,R1,R2,R3,R4
      REAL A,B,C,D,E,A1,B1,AW,BW,K,KO,KW,K35
      EQUIVALENCE (E,D,B1), (BW,B,R3), (C,A1,R2)
      EQUIVALENCE (AW,A,R1), (KW,KO,K)

      DATA R3500, R4/1028.1063, 4.8314E-4/
      DATA DR350/28.106331/

      P = P0/10.
      SR = SQRT(ABS(S))

      R1 = ((((6.536332E-9*T-1.120083E-6)*T+1.001685E-4)*T              &
     & -9.095290E-3)*T+6.793952E-2)*T-28.263737
      R2 = (((5.3875E-9*T-8.2467E-7)*T+7.6438E-5)*T-4.0899E-3)*T        &
     & + 8.24493E-1
      R3 = (-1.6546E-6*T+1.0227E-4)*T-5.72466E-3
      SIG = (R4*S + R3*SR + R2)*S + R1

      V350P = 1.0/R3500
      SVA = -SIG*V350P/(R3500+SIG)
      SIGMA = SIG+DR350

      SVAN = SVA*1.0E+8
      IF(P.EQ.0.0) RETURN

      E = (9.1697E-10*T+2.0816E-8)*T-9.9348E-7
      BW = (5.2787E-8*T-6.12293E-6)*T+3.47718E-5
      B = BW + E*S

      D = 1.91075E-4
      C = (-1.6078E-6*T-1.0981E-5)*T+2.2838E-3
      AW = ((-5.77905E-7*T+1.16092E-4)*T+1.43713E-3)*T - 0.1194975
      A = AW + (C + D*SR)*S

      B1 = (-5.3009E-4*T+1.6483E-2)*T+7.944E-2
      A1 = ((-6.1670E-5*T+1.09987E-2)*T-0.603459)*T+54.6746
      KW = ((-5.155288E-5*T+1.360477E-2)*T-2.327105)*T                  &
     & + 148.4206
      KW = KW*T - 1930.06
      KO = KW + (A1 + B1*SR)*S

      DK = (B*P + A)*P + KO
      K35 = (5.03217E-5*P+3.359406)*P+21582.27
      GAM = P/K35
      PK = 1.0 - GAM
      SVA = SVA*PK + (V350P+SVA)*P*DK/(K35*(K35+DK))

      SVAN = SVA*1.0E+8
      V350P = V350P*PK

      DR35P = GAM/V350P
      DVAN = SVA/(V350P*(V350P+SVA))
      SIGMA = DR350 + DR35P - DVAN
      RETURN
      END

! ------------------------------------------------------------------
! 4. KONVERSI TEKANAN KE KEDALAMAN (DEPTH)
! ------------------------------------------------------------------
      REAL FUNCTION DEPTH(P, LAT)
      REAL LAT
      X = SIN(LAT/57.29578)
      X = X*X
      
      GR = 9.780318*(1.0+(5.2788E-3+2.36E-5*X)*X) + 1.092E-6*P
      DEPTH = (((-1.82E-15*P+2.279E-10)*P-2.2512E-5)*P+9.72659)*P
      DEPTH = DEPTH/GR
      RETURN
      END

! ------------------------------------------------------------------
! 5. TEMPERATUR TITIK BEKU (TF)
! ------------------------------------------------------------------
      REAL FUNCTION TF(S,P)
      TF = (-.0575+1.710523E-3*SQRT(ABS(S))-2.154996E-4*S)*S - 7.53E-4*P
      RETURN
      END

! ------------------------------------------------------------------
! 6. PANAS JENIS (CPSW)
! ------------------------------------------------------------------
      REAL FUNCTION CPSW(S,T,P0)
      P = P0/10.
      SR = SQRT(ABS(S))

      A = (-1.38385E-3*T+0.1072763)*T-7.643575
      B = (5.148E-5*T-4.07718E-3)*T+0.1770383
      C = (((2.093236E-5*T-2.654387E-3)*T+0.1412855)*T                  &
     & -3.720283)*T+4217.4
      CP0 = (B*SR + A)*S + C

      A = (((1.7168E-8*T+2.0357E-6)*T-3.13885E-4)*T+1.45747E-2)*T       &
     & -0.49592
      B = (((2.2956E-11*T-4.0027E-9)*T+2.87533E-7)*T-1.08645E-5)*T      &
     & + 2.4931E-4
      C = ((6.136E-13*T-6.5637E-11)*T+2.6380E-9)*T-5.422E-8
      CP1 = ((C*P+B)*P+A)*P

      A = (((-2.9179E-10*T+2.5941E-8)*T+9.802E-7)*T-1.28315E-4)*T       &
     & + 4.9247E-3
      B = (3.122E-8*T-1.517E-6)*T-1.2331E-4
      A = (A+B*SR)*S

      B = ((1.8448E-11*T-2.3905E-9)*T+1.17054E-7)*T-2.9558E-6
      B = (B+9.971E-8*SR)*S

      C = (3.513E-13*T-1.7682E-11)*T+5.540E-10
      C = (C-1.4300E-12*T*SR)*S
      CP2 = ((C*P+B)*P+A)*P

      CPSW = CP0 + CP1 + CP2
      RETURN
      END

! ------------------------------------------------------------------
! 7. ADIABATIC LAPSE RATE (ATG)
! ------------------------------------------------------------------
      REAL FUNCTION ATG(S,T,P)
      DS = S - 35.0
      ATG = (((-2.1687E-16*T+1.8676E-14)*T-4.6206E-13)*P                &
     & + ((2.7759E-12*T-1.1351E-10)*DS+((-5.4481E-14*T                  &
     & + 8.733E-12)*T-6.7795E-10)*T+1.8741E-8))*P                       &
     & + (-4.2393E-8*T+1.8932E-6)*DS                                    &
     & + ((6.6228E-10*T-6.836E-8)*T+8.5258E-6)*T+3.5803E-5
      RETURN
      END

! ------------------------------------------------------------------
! 8. TEMPERATUR POTENSIAL (THETA)
! ------------------------------------------------------------------
      REAL FUNCTION THETA(S,TO,PO,PR)
      P = PO
      T = TO
      H = PR - P

      XK = H*ATG(S,T,P)
      T = T + 0.5*XK
      Q = XK
      P = P + 0.5*H

      XK = H*ATG(S,T,P)
      T = T + 0.29289322*(XK-Q)
      Q = 0.58578644*XK + 0.121320344*Q

      XK = H*ATG(S,T,P)
      T = T + 1.707106781*(XK-Q)
      Q = 3.0*XK - 4.121320344*Q
      P = P + 0.5*H

      XK = H*ATG(S,T,P)
      THETA = T + (XK-2.0*Q)/6.0
      RETURN
      END

! ------------------------------------------------------------------
! 9. KECEPATAN SUARA (SVEL)
! ------------------------------------------------------------------
      REAL FUNCTION SVEL(S,T,P0)
      EQUIVALENCE (A0,B0,C0), (A1,B1,C1), (A2,C2), (A3,C3)

      P = P0/10.
      SR = SQRT(ABS(S))
      D = 1.727E-3 - 7.9836E-6*P

      B1 = 7.3637E-5 + 1.7945E-7*T
      B0 = -1.922E-2 - 4.42E-5*T
      B = B0 + B1*P

      A3 = (-3.389E-13*T+6.649E-12)*T+1.100E-10
      A2 = ((7.988E-12*T-1.6002E-10)*T+9.1041E-9)*T-3.9064E-7
      A1 = (((-2.0122E-10*T+1.0507E-6)*T-6.4885E-8)*T-1.2580E-5)*T      &
     & + 9.4742E-5
      A0 = (((-3.21E-8*T+2.006E-6)*T+7.164E-5)*T-1.262E-2)*T+1.389
      A = ((A3*P+A2)*P+A1)*P+A0

      C3 = (-2.3643E-12*T+3.8504E-10)*T-9.7729E-9
      C2 = (((1.0405E-12*T-2.5335E-10)*T+2.5974E-8)*T-1.7107E-6)*T      &
     & + 3.1260E-5
      C1 = (((-6.1185E-10*T+1.3621E-7)*T-8.1788E-6)*T+6.8982E-4)*T      &
     & + 0.153563
      C0 = ((((3.1464E-9*T-1.47800E-6)*T+3.3420E-4)*T-5.80852E-2)*T     &
     & + 5.03711)*T+1402.388
      C = ((C3*P+C2)*P+C1)*P+C0

      SVEL = C + (A + B*SR + D*S)*S
      RETURN
      END