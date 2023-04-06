RateSolver == SummationDecapode(parse_decapode(quote
    ρ_Ap == ρ_N2p + ρ_Np + ρ_O2p + ρ_Op + ρ_O4p + ρ_NOp + ρ_Yp

    ρ_Bm == ρ_Om + ρ_O2m + ρ_O3m + ρ_O4m + ρ_OHm + ρ_CO3m + ρ_CO4m + ρ_NO2m +
        ρ_NO3m + ρ_O2mNO + ρ_HCO3m

    # N2(A)
    P_N2A == r10 + r119 + r120
    L_N2A == r114 + r115 + r116 + r117 + r118

    # N2(B)
    P_N2B == r11 + r122 + r125
    L_N2B == r119 + r120 + r121

    # N2(a')
    P_N2a == r12 + r126
    L_N2a == r122 + r123 + r124
        
    # N2(C)
    P_N2C == r13
    L_N2C == r125 + r126 + r127

    # N(4S)
    P_N4S == r3 + 2 .*r14 + r15 + r85 + r88 + r90 + r101 + r102 + r104 + r106 +
        r107 + r124 + r153 + r154 + r155

    L_N4S == r32 + r95 + r118 + r131 + r139 + r140 + r141 + r143 + r148 + r149 +
        r150 + r151

    # N*
    P_N2D == r4 + r15 + r93 + r101 + r102 + 2 .*r103 + r105 + r108 + r109 + r117

    L_N2D == r142 + r144 + r153 + r154 + r155

    # O2(a)
    P_O2a == r17 + r114 + r135 + r136 + r137 + r157

    L_O2a == r34 + r37 + r44 + r61 + r128 + r129 + r130 + r131 + r132 + r133

    # O2(b)
    P_O2b == r18 + r115
    L_O2b == r134 + r135 + r136 + r137 + r138

    # O  see bottom for details
    P_O == 0   
    L_O == 0

    # NO
    P_NO == r32 + r77 + r79 + r87 + r92 + r94 + r117 + r131 + r139 + r140 +
        r141 + r142 + 2 .*r151 + r152

    L_NO == r9 + r36 + r49 + r59 + r64 + r71 + r72 + r76 + r83 + r96 + r124 +
        r143 + r144 + r145 + r146 + r147 

    # NO2
    P_NO2 == r36 + r43 + r82 + r83 + r145 + r146 + r147

    L_NO2 == r50 + r55 + r60 + r66 + r79 + r148 + r149 + r150 + r151 + r152

    # e-
    P_e == r1 + r2 + r5 + r6 + r8 + r9 + r16 + r22 + r23 + r30 + r31 + r32 +
        r33 + r34 + r35 + r36 + r37 + r38 + r39 + r40 + r41 + r42 + r43

    L_e == r25 + r26 + r27 + r28 + r29 + r101 + r102 + r103 + r104 + r105 +
        r106 + r107 + r108 + r109 + r110 + r111 + r112

    # O-
    P_Om == r24 + r27 + r28 + r51
    L_Om == r30 + r31 + r32 + r33 + r34 + r35 + r36 + r44 + r45 + r46 + r47 +
        r48 + r49 + r50 + k113 .* ρ_Ap .* ρ_Om

    # O2-
    P_O2m == r25 + r26 + r29 + r44 + r56 + r61 + r68

    L_O2m == r37 + r38 + r51 + r52 + r53 + r54 + r55 + k113 .* ρ_Ap .* ρ_O2m

    # O3-
    P_O3m == r45 + r46 + r52 + r62 + r65 + r69 + r73

    L_O3m == r39 + r40 + r56 + r57 + r58 + r59 + r60 + k113 .* ρ_Ap .* ρ_O3m

    # O4-
    P_O4m == r53
    L_O4m == r61 + r62 + r63 + r64 + k113 .* ρ_Ap .* ρ_O4m

    # OH-
    P_OHm == r47 + r57 + r70 + r77
    L_OHm == r41 + r42 + r65 + r66 + r67 + k113 .* ρ_Ap .* ρ_OHm

    # CO3-
    P_CO3m == r48 + r58 + r74 + r75 + r82
    L_CO3m == r68 + r69 + r70 + r71 + r72 + k113 .* ρ_Ap .* ρ_CO3m

    # CO4-
    P_CO4m == r54 + r63
    L_CO4m == r73 + r74 + r75 + r76 + k113 .* ρ_Ap .* ρ_CO4m

    # NO2-
    P_NO2m == r49 + r50 + r55 + r66 + r71 + r80 + r81 + r83 + r84

    L_NO2m == r43 + r77 + r78 + r79 + k113 .* ρ_Ap .* ρ_NO2m

    # NO3-
    P_NO3m == r59 + r60 + r72 + r78 + r79
    L_NO3m == r80 + r81 + k113 .* ρ_Ap .* ρ_NO3m

    # O2-.NO
    P_O2mNO == r64 + r76
    L_O2mNO == r82 + r83 + r84 + k113 .* ρ_Ap .* ρ_O2mNO

    # HCO3-
    P_HCO3m == r67
    L_HCO3m == k113 .* ρ_Ap .* ρ_HCO3m

    # N2+
    P_N2p == r1 + r16
    L_N2p == r91 + r92 + r93 + r101 + r102 + r103 + k113 .* ρ_N2p .* ρ_Bm

    # N+
    P_Np == r2
    L_Np == r85 + r86 + r87 + k113 .* ρ_Np .* ρ_Bm

    # O2+
    P_O2p== r5 + r8 + r22 + r85 + r89 + r91 + r98
    L_O2p == r94 + r95 + r96 + r97 + r110 + k113 .* ρ_O2p .* ρ_Bm

    # O+
    P_Op == r6 + r23 + r24 + r87
    L_Op == r88 + r89 + r90 + k113 .* ρ_Op .* ρ_Bm

    # O4+
    P_O4p == r97
    L_O4p == r98 + r99 + r111 + k113 .* ρ_O4p .* ρ_Bm

    # NO+
    P_NOp == r9 + r86 + r88 + r90 + r92 + r93 + r94 + r95 + r96
    L_NOp == r100 + r104 + r105 + r106 + r107 + r108 + r109 + k113 .* ρ_NOp .*
        ρ_Bm

    # Y+
    P_Yp == r99 + r100
    L_Yp == r112 + k113 .* ρ_Yp .* ρ_Bm




    # ATOMIC OXYGEN production and loss rates depend heavily upon the odd
    # hydrogen family, which are not provided for by the NRLMSISE modeL_  Given
    # that its production via electric field heating is small, the density of 
    # atomic oxygen may safetly be considered to be constant throughout the 
    # simulation.  NOTE: This approximation must be verified, and can be done
    # so most readily by comparison with electron-impact dissociation of N2 .

    P_O == r7 + 2 .*r19 + 2 .*r20 + 2 .*r21 + r23 + r27 + r29 + r30 + r44 +
        r46 + r50 + r59 + r86 + r89 + r95 + r104 + r105 + r106 + r107 + r108 +
        r109 + 2 .*r110 + 2 .*r116 + 2 .*r121 + 2 .*r123 + r124 + 2 .*r127 +
        r131 + r138 + r141 + r142 + r143 + r144 + r148 + 2 .*r149 + r163 + r165
    
    L_O == r35 + r39 + r41 + r51 + r56 + r62 + r68 + r74 + r80 + r93 +
        r117 + r139 + r145 + r146 + r152 + r156 + r157 + r158 + r159 + r160 +
        r161 + r162



    # POSITIVE NITROGEN ION, N+, densities depend only upon production via GCR,
    # and loss by charge exchange with O2 .  Given that the production and loss
    # rates are entirely unaffected by electric field heating (dissociative
    # ionization of N2 is not included), the density of N+ remains at its 
    # equilibrium value throughout the simulation.  Due to its volatile nature 
    # (large time-rate of change relative to its low density), can be taken as
    # constant for heating simulations.
    
    # P_Np == r2
    # 
    # L_Np == r85 + r86 + r87 + k113 .* ρ_Np .* ρ_Bm
end))