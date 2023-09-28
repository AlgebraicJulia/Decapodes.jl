Chemistry = SummationDecapode(parse_decapode(quote
    Tn::Form0{X}
    θ::Form0{X}
    ρ_e::Form0{X}
    ρ_gas::Form0{X}
    kB::Constant{X}

    ρ_gas == ρ_N2 + ρ_O2

    Te_Less_65::Form0
    Te_Less_65 == .≤(θ, 65)
    Te_Greater_65::Form0
    Te_Greater_65 == .>(θ, 65)

    # See Kotovsky pp. 92 Eq. 5-1
    # Note that the constant (3/2) is already divided through.
    Te_unthresholded::Form0
    Te_unthresholded == (Te_Less_65 * 1.54e-2 * θ) +
        (Te_Greater_65 * 2 * (θ ./ 65).^(2.6) ./ (1 .+ (θ ./65) .^2))

    Te_Greater_kB_times_Tn::Form0
    Te_Greater_kB_times_Tn == .>(Te_unthresholded, (kB * Tn))
    Te_Less_kB_times_Tn::Form0
    Te_Less_kB_times_Tn == .≤(Te_unthresholded, (kB * Tn))
    # Effective electron temperature [eV]
    Te::Form0
    Te == (Te_Greater_kB_times_Tn .* Te_unthresholded) +
        (Te_Less_kB_times_Tn .* (kB * Tn))

    # Electron temperature [K]
    TeK::Form0
    TeK == Te ./ kB

    # Electron impact excitation, dissociation, and ionization
    # TODO: Annotate each equation with the paper it came from. See Kotovsky Table B-1 references
    # Note that some papers, like Kossyi et al. (92) have some equations cited from sources that Kotovsky did not cite, like for the equation of ion-ion recombination. (It may be of interest to expand these reference dependencies and cite them here.)
    # TODO: Check whether this should be .> or .<
    k10_through_24_mask == .>(Te, 0.1)
    # TODO: Check if θ can ever be zero.
    k10 == k10_through_24_mask .* (10 .^ (-1 * (8.4 .+ 140 ./ θ)))
    k11 == k10_through_24_mask .* (10 .^ (-1 * (8.2 .+ 148 ./ θ)) + 10 .^ (-1 * (8.3 .+ 154 ./ θ)) + 10 .^ (-1 * (8.7 .+ 168 ./ θ)))
    k12 == k10_through_24_mask .* (10 .^ (-1 * (8.8 .+ 167 ./ θ)) + 10 .^ (-1 * (8.5 .+ 174 ./ θ)) + 10 .^ (-1 * (8.7 .+ 175 ./ θ)))
    k13 == k10_through_24_mask .* (10 .^ (-1 * (8.2 .+ 211 ./ θ)) + 10 .^ (-1 * (10.1 .+ 254 ./ θ)) + 10 .^ (-1 * (9.2 .+ 262 ./ θ)))

    k14_mask == .>(θ, 76)
    k14 == k10_through_24_mask .* k14_mask .* (2e-10 .* (3.096 .- 6.71e-2.*θ + 3e-4.*(θ .^2) + 1.59e-6.*(θ .^3) .- 1.57e-9.*(θ .^4)))

    k15 == k10_through_24_mask .* (k14 ./2)
    k16 == k10_through_24_mask .* (10 .^ (-1 * (8.3 + 365 ./θ)))

    k17_mask == .>(θ, 40)
    k17 == k10_through_24_mask .* (
        k17_mask .* (10 .^ (-1 * (10.2 + 3.5 ./θ))) +
        invert_mask(k17_mask) .* (10 .^ (-1 * (9 + 52 ./θ))))

    k18_mask == .>(θ, 30)
    k18 == k10_through_24_mask .* (
        k18_mask .* ( 10 .^ (- 1 * (11.2 + 7.2./θ))) +
        invert_mask(k18_mask) .* (10 .^ (-1 * (9.5 + 60 ./θ))))

    k19 == k10_through_24_mask .* (10 .^ (-1 * (7.9 + 134 ./θ)))
    k20 == k10_through_24_mask .* (10 .^ (-1 * (8 + 169 ./θ)))
    k21 == k10_through_24_mask .* (10 .^ (-1 * (8.8 + 119 ./θ)))
    k22 == k10_through_24_mask .* (10 .^ (-1 * (8.8 + 281 ./θ)))
    k23 == k10_through_24_mask .* (3.18e-14 * (TeK ./300) .* exp(-206055 ./TeK))
    k24 == k10_through_24_mask .* (7.10e-11 * Te .^(0.5) .* exp(-17 ./Te))

    # Electron attachment
    k25 == 1.07e-31 .* (300 ./TeK).^2 .* exp(-70 ./TeK) .* exp(1500 .*(TeK .-Tn)./TeK./Tn)
    k26 == 1.4e-29 .* (300 ./TeK) .* exp(-600 ./TeK) .* exp(700 .*(TeK .-Tn)./TeK./Tn)
    k27 == 1.07e-9 .* Te .^ (-1.39) .* exp(-6.26 ./Te)

    k28_mask == .>(Te, 0.13)
    k28 == k28_mask .* (2.12e-9 .* Te .^ (-1.06) .* exp(-0.93 ./Te)) +
        invert_mask(k28_mask) .* 1e-11

    k29_mask == .>(Te, 0.14)
    k29 == k29_mask .* (9.76e-8 .* Te .^ (-1.309) .* exp(-1.007 ./Te)) +
        invert_mask(k29_mask) .* 1e-9

    # Electron detachment
    k30 == 5.47e-8 * Te .^ (0.324) .* exp(-2.98 ./Te)

    k31_mask == .≥(Te, 0.1)
    k31 == k31_mask .* (10 .^ (-15.66 + 2.97 .* log10(θ) - 0.58 .* (log10(θ) .^ 2)))

    # Electron-ion recombination
    k101 == 0.37 * 2.2e-7 * (300 ./Tn).^(0.2) .* (TeK./Tn).^(-0.39)
    k102 == 0.11 * 2.2e-7 * (300 ./Tn).^(0.2) .* (TeK./Tn).^(-0.39)
    k103 == 0.52 * 2.2e-7 * (300 ./Tn).^(0.2) .* (TeK./Tn).^(-0.39)

    # Observe: k104 and k105 are set to 0 where Te > 1.25, but k106, k107, k108,
    # and k109 are set to 0 where Te ≤ 1.25.
    k104_through_105_mask == .≤(Te, 1.25)
    k106_through_109_mask == .>(Te, 1.25)
    
    k104 == k104_through_105_mask .* (0.05 * 3.5e-7 .* (300 ./Tn) .* (TeK./Tn).^(-0.69))

    k105 == k104_through_105_mask .* (0.95 * 3.5e-7 .* (300 ./Tn) .* (TeK./Tn).^(-0.69))

    k106 == k106_through_109_mask .* (0.10 * 3.5e-7 .* (300 ./Tn) .* (TeK./Tn).^(-0.69))

    k107 == k106_through_109_mask .* (0.10 * 3.5e-7 .* (300 ./Tn) .* (TeK./Tn).^(-0.69))

    k108 == k106_through_109_mask .* (0.70 * 3.5e-7 .* (300 ./Tn) .* (TeK./Tn).^(-0.69))

    k109 == k106_through_109_mask .* (0.10 * 3.5e-7 .* (300 ./Tn) .* (TeK./Tn).^(-0.69))

    k110 == 1.95e-7 * (300 ./Tn).^(0.7) .* (TeK./Tn).^(-0.7)
    k111 == 4.2e-6 * (TeK./Tn).^(-0.48)
    k112 == 6.5e-6 * (300 ./TeK).^(1/2) # Recombination rate for cluster order n = 3

    # Total gas density (dominated by N2 and O2)
    ρ_M == ρ_N2 + ρ_O2

    # TODO: Maybe break down each class of rate computation into a composition pattern.
    # Cosmic ray and photo-ionization
    r1 == k1 .* ρ_N2
    r2 == k2 .* ρ_N2
    r3 == k3 .* ρ_N2
    r4 == k4 .* ρ_N2
    r5 == k5 .* ρ_O2
    r6 == k6 .* ρ_O2
    r7 == k7 .* ρ_O2
    r8 == k8 .* ρ_O2
    r9 == k9 .* ρ_NO

    # Electron impact excitation, dissociation, and ionization
    r10 == k10 .* ρ_e .* ρ_N2
    r11 == k11 .* ρ_e .* ρ_N2
    r12 == k12 .* ρ_e .* ρ_N2
    r13 == k13 .* ρ_e .* ρ_N2
    r14 == k14 .* ρ_e .* ρ_N2
    r15 == k15 .* ρ_e .* ρ_N2
    r16 == k16 .* ρ_e .* ρ_N2
    r17 == k17 .* ρ_e .* ρ_O2
    r18 == k18 .* ρ_e .* ρ_O2
    r19 == k19 .* ρ_e .* ρ_O2
    r20 == k20 .* ρ_e .* ρ_O2
    r21 == k21 .* ρ_e .* ρ_O2
    r22 == k22 .* ρ_e .* ρ_O2
    r23 == k23 .* ρ_e .* ρ_O2
    r24 == k24 .* ρ_e .* ρ_O2

    # Electron attachment
    r25 == k25 .*ρ_e .* ρ_O2 .* ρ_N2
    r26 == k26 .*ρ_e .* ρ_O2 .* ρ_O2
    r27 == k27 .*ρ_e .* ρ_O2
    r28 == k28 .*ρ_e .* ρ_O3
    r29 == k29 .*ρ_e .* ρ_O3

    # Electron detachment
    r30 == k30 .* ρ_Om .* ρ_e
    r31 == k31 .* ρ_Om .* ρ_N2
    r32 == k32 .* ρ_Om .* ρ_N4S
    r33 == k33 .* ρ_Om .* ρ_O2
    r34 == k34 .* ρ_Om .* ρ_O2a
    r35 == k35 .* ρ_Om .* ρ_O
    r36 == k36 .* ρ_Om .* ρ_NO
    r37 == k37 .* ρ_O2m .* ρ_O2a
    r38 == k38 .* ρ_O2m .* ρ_O
    r39 == k39 .* ρ_O3m .* ρ_O
    r40 == k40 .* ρ_O3m .* ρ_O3
    r41 == k41 .* ρ_OHm .* ρ_O
    r42 == k42 .* ρ_OHm .* ρ_H
    r43 == k43 .* ρ_NO2m .* ρ_O

    # Negative ion conversion
    r44 == k44 .* ρ_Om .* ρ_O2a
    r45 == k45 .* ρ_Om .* ρ_O2 .* ρ_M
    r46 == k46 .* ρ_Om .* ρ_O3
    r47 == k47 .* ρ_Om .* ρ_H2O
    r48 == k48 .* ρ_Om .* ρ_CO2 .* ρ_M
    r49 == k49 .* ρ_Om .* ρ_NO .* ρ_M
    r50 == k50 .* ρ_Om .* ρ_NO2
    r51 == k51 .* ρ_O2m .* ρ_O
    r52 == k52 .* ρ_O2m .* ρ_O3
    r53 == k53 .* ρ_O2m .* ρ_O2 .* ρ_M
    r54 == k54 .* ρ_O2m .* ρ_CO2 .* ρ_M
    r55 == k55 .* ρ_O2m .* ρ_NO2
    r56 == k56 .* ρ_O3m .* ρ_O
    r57 == k57 .* ρ_O3m .* ρ_H
    r58 == k58 .* ρ_O3m .* ρ_CO2
    r59 == k59 .* ρ_O3m .* ρ_NO
    r60 == k50 .* ρ_O3m .* ρ_NO2
    r61 == k61 .* ρ_O4m .* ρ_O2a
    r62 == k62 .* ρ_O4m .* ρ_O
    r63 == k63 .* ρ_O4m .* ρ_CO2
    r64 == k64 .* ρ_O4m .* ρ_NO
    r65 == k65 .* ρ_OHm .* ρ_O3
    r66 == k66 .* ρ_OHm .* ρ_NO2
    r67 == k67 .* ρ_OHm .* ρ_CO2 .* ρ_M
    r68 == k68 .* ρ_CO3m .* ρ_O
    r69 == k69 .* ρ_CO3m .* ρ_O2
    r70 == k70 .* ρ_CO3m .* ρ_H
    r71 == k71 .* ρ_CO3m .* ρ_NO
    r72 == k72 .* ρ_CO3m .* ρ_NO2
    r73 == k73 .* ρ_CO4m .* ρ_O3
    r74 == k74 .* ρ_CO4m .* ρ_O
    r75 == k75 .* ρ_CO4m .* ρ_H
    r76 == k76 .* ρ_CO4m .* ρ_NO
    r77 == k77 .* ρ_NO2m .* ρ_H
    r78 == k78 .* ρ_NO2m .* ρ_O3
    r79 == k79 .* ρ_NO2m .* ρ_NO2
    r80 == k80 .* ρ_NO3m .* ρ_O
    r81 == k81 .* ρ_NO3m .* ρ_O3
    r82 == k82 .* ρ_O2mNO .* ρ_CO2
    r83 == k83 .* ρ_O2mNO .* ρ_NO
    r84 == k84 .* ρ_O2mNO .* ρ_H

    # Positive ion conversion
    r85 == k85 .* ρ_Np .* ρ_O2
    r86 == k86 .* ρ_Np .* ρ_O2
    r87 == k87 .* ρ_Np .* ρ_O2
    r88 == k88 .* ρ_Op .* ρ_N2
    r89 == k89 .* ρ_Op .* ρ_O2
    r90 == k90 .* ρ_Op .* ρ_N2 .* ρ_M
    r91 == k91 .* ρ_N2p .* ρ_O2
    r92 == k92 .* ρ_N2p .* ρ_O2
    r93 == k93 .* ρ_N2p .* ρ_O
    r94 == k94 .* ρ_O2p .* ρ_N2
    r95 == k95 .* ρ_O2p .* ρ_N4S
    r96 == k96 .* ρ_O2p .* ρ_NO
    r97 == k97 .* ρ_O2p .* ρ_O2 .* ρ_M
    r98 == k98 .* ρ_O4p .* ρ_O
    r99 == k99 .* ρ_O4p .* ρ_H2O
    r100 == k100 .* ρ_NOp .* ρ_M .* ρ_M

    # Electrion-ion recombination
    r101 == k101 .* ρ_N2p .* ρ_e
    r102 == k102 .* ρ_N2p .* ρ_e
    r103 == k103 .* ρ_N2p .* ρ_e
    r104 == k104 .* ρ_NOp .* ρ_e
    r105 == k105 .* ρ_NOp .* ρ_e
    r106 == k106 .* ρ_NOp .* ρ_e
    r107 == k107 .* ρ_NOp .* ρ_e
    r108 == k108 .* ρ_NOp .* ρ_e
    r109 == k109 .* ρ_NOp .* ρ_e
    r110 == k110 .* ρ_O2p .* ρ_e
    r111 == k111 .* ρ_O4p .* ρ_e
    r112 == k112 .* ρ_Yp .* ρ_e

    # Ion-ion recombination
    # r113 determined in rateSolver.m

    # Active-state N2 chemistry
    r114 == k114 .* ρ_N2A .* ρ_O2
    r115 == k115 .* ρ_N2A .* ρ_O2
    r116 == k116 .* ρ_N2A .* ρ_O2
    r117 == k117 .* ρ_N2A .* ρ_O
    r118 == k118 .* ρ_N2A .* ρ_N4S
    r119 == k119 .* ρ_N2B
    r120 == k120 .* ρ_N2B .* ρ_N2
    r121 == k121 .* ρ_N2B .* ρ_O2
    r122 == k122 .* ρ_N2a .* ρ_N2
    r123 == k123 .* ρ_N2a .* ρ_O2
    r124 == k124 .* ρ_N2a .* ρ_NO
    r125 == k125 .* ρ_N2C
    r126 == k126 .* ρ_N2C .* ρ_N2
    r127 == k127 .* ρ_N2C .* ρ_O2

    # Active-state O2 chemistry
    r128 == k128 .* ρ_O2a
    r129 == k129 .* ρ_O2a .* ρ_N2
    r130 == k130 .* ρ_O2a .* ρ_O2
    r131 == k131 .* ρ_O2a .* ρ_N4S
    r132 == k132 .* ρ_O2a .* ρ_O
    r133 == k133 .* ρ_O2a .* ρ_NO
    r134 == k134 .* ρ_O2b
    r135 == k135 .* ρ_O2b .* ρ_N2
    r136 == k136 .* ρ_O2b .* ρ_O2
    r137 == k137 .* ρ_O2b .* ρ_O
    r138 == k138 .* ρ_O2b .* ρ_O3

    # Odd nitrogen chemistry
    r139 == k139 .* ρ_N4S .* ρ_O .* ρ_M
    r140 == k140 .* ρ_N4S .* ρ_O3
    r141 == k141 .* ρ_N4S .* ρ_O2
    r142 == k142 .* ρ_N2D .* ρ_O2
    r143 == k143 .* ρ_N4S .* ρ_NO
    r144 == k144 .* ρ_N2D .* ρ_NO
    r145 == k145 .* ρ_NO .* ρ_O
    r146 == k146 .* ρ_NO .* ρ_O .* ρ_M
    r147 == k147 .* ρ_NO .* ρ_O3
    r148 == k148 .* ρ_N4S .* ρ_NO2
    r149 == k149 .* ρ_N4S .* ρ_NO2
    r150 == k150 .* ρ_N4S .* ρ_NO2
    r151 == k151 .* ρ_N4S .* ρ_NO2
    r152 == k152 .* ρ_O .* ρ_NO2
    r153 == k153 .* ρ_N2D
    r154 == k154 .* ρ_N2D .* ρ_N2
    r155 == k155 .* ρ_N2D .* ρ_O

    # Odd oxygen chemistry
    r156 == k156 .* ρ_O .* ρ_O .* ρ_M
    r157 == k157 .* ρ_O .* ρ_O .* ρ_M
    r158 == k158 .* ρ_O .* ρ_O2 .* ρ_N2
    r159 == k159 .* ρ_O .* ρ_O2 .* ρ_O2
    r160 == k160 .* ρ_O .* ρ_O3
    r161 == k161 .* ρ_O .* ρ_OH 
    r162 == k162 .* ρ_O .* ρ_HO2
    r163 == k163 .* ρ_OH .* ρ_OH
    r164 == k164 .* ρ_OH .* ρ_O3
    r165 == k165 .* ρ_H .* ρ_HO2
    r166 == k166 .* ρ_H .* ρ_O3
    r167 == k167 .* ρ_HO2 .* ρ_O3

    # Positively charged ions are grouped.
    ρ_Ap == ρ_N2p + ρ_Np + ρ_O2p + ρ_Op + ρ_O4p + ρ_NOp + ρ_Yp

    # Negatively charged ions are grouped.
    ρ_Bm == ρ_Om + ρ_O2m + ρ_O3m + ρ_O4m + ρ_OHm + ρ_CO3m + ρ_CO4m + ρ_NO2m +
        ρ_NO3m + ρ_O2mNO + ρ_HCO3m

    # N2(A)
    P_N2A == r10 + r119 + r120
    L_N2A == r114 + r115 + r116 + r117 + r118
    ∂ₜ(ρ_N2A) == maskalt(P_N2A - L_N2A)
    # TODO: Either hard code the above line for the other 38 species or use a
    # composition pattern.

    # N2(B)
    P_N2B == r11 + r122 + r125
    L_N2B == r119 + r120 + r121
    ∂ₜ(ρ_N2B) == maskalt(P_N2B - L_N2B)

    # N2(a')
    P_N2a == r12 + r126
    L_N2a == r122 + r123 + r124
    ∂ₜ(ρ_N2a) == maskalt(P_N2a - L_N2a)
        
    # N2(C)
    P_N2C == id(r13)
    L_N2C == r125 + r126 + r127
    ∂ₜ(ρ_N2C) == maskalt(P_N2C - L_N2C)

    # N(4S)
    P_N4S == r3 + 2 .*r14 + r15 + r85 + r88 + r90 + r101 + r102 + r104 + r106 +
        r107 + r124 + r153 + r154 + r155

    L_N4S == r32 + r95 + r118 + r131 + r139 + r140 + r141 + r143 + r148 + r149 +
        r150 + r151
    ∂ₜ(ρ_N4S) == maskalt(P_N4S - L_N4S)

    # N*
    P_N2D == r4 + r15 + r93 + r101 + r102 + 2 .*r103 + r105 + r108 + r109 + r117

    L_N2D == r142 + r144 + r153 + r154 + r155
    ∂ₜ(ρ_N2D) == maskalt(P_N2D - L_N2D)

    # O2(a)
    P_O2a == r17 + r114 + r135 + r136 + r137 + r157

    L_O2a == r34 + r37 + r44 + r61 + r128 + r129 + r130 + r131 + r132 + r133
    ∂ₜ(ρ_O2a) == maskalt(P_O2a - L_O2a)

    # O2(b)
    P_O2b == r18 + r115
    L_O2b == r134 + r135 + r136 + r137 + r138
    ∂ₜ(ρ_O2b) == maskalt(P_O2b - L_O2b)

    # O  see bottom for details
    #P_O == 0   
    #L_O == 0

    # NO
    P_NO == r32 + r77 + r79 + r87 + r92 + r94 + r117 + r131 + r139 + r140 +
        r141 + r142 + 2 .*r151 + r152

    L_NO == r9 + r36 + r49 + r59 + r64 + r71 + r72 + r76 + r83 + r96 + r124 +
        r143 + r144 + r145 + r146 + r147 
    ∂ₜ(ρ_NO) == maskalt(P_NO - L_NO)

    # NO2
    P_NO2 == r36 + r43 + r82 + r83 + r145 + r146 + r147

    L_NO2 == r50 + r55 + r60 + r66 + r79 + r148 + r149 + r150 + r151 + r152
    ∂ₜ(ρ_NO2) == maskalt(P_NO2 - L_NO2)

    # e-
    P_e == r1 + r2 + r5 + r6 + r8 + r9 + r16 + r22 + r23 + r30 + r31 + r32 +
        r33 + r34 + r35 + r36 + r37 + r38 + r39 + r40 + r41 + r42 + r43

    L_e == r25 + r26 + r27 + r28 + r29 + r101 + r102 + r103 + r104 + r105 +
        r106 + r107 + r108 + r109 + r110 + r111 + r112
    ∂ₜ(ρ_e) == maskalt(P_e - L_e)

    # O-
    P_Om == r24 + r27 + r28 + r51
    L_Om == r30 + r31 + r32 + r33 + r34 + r35 + r36 + r44 + r45 + r46 + r47 +
        r48 + r49 + r50 + k113 .* ρ_Ap .* ρ_Om
    ∂ₜ(ρ_Om) == maskalt(P_Om - L_Om)

    # O2-
    P_O2m == r25 + r26 + r29 + r44 + r56 + r61 + r68

    L_O2m == r37 + r38 + r51 + r52 + r53 + r54 + r55 + k113 .* ρ_Ap .* ρ_O2m
    ∂ₜ(ρ_O2m) == maskalt(P_O2m - L_O2m)

    # O3-
    P_O3m == r45 + r46 + r52 + r62 + r65 + r69 + r73

    L_O3m == r39 + r40 + r56 + r57 + r58 + r59 + r60 + k113 .* ρ_Ap .* ρ_O3m
    ∂ₜ(ρ_O3m) == maskalt(P_O3m - L_O3m)

    # O4-
    P_O4m == id(r53)
    L_O4m == r61 + r62 + r63 + r64 + k113 .* ρ_Ap .* ρ_O4m
    ∂ₜ(ρ_O4m) == maskalt(P_O4m - L_O4m)

    # OH-
    P_OHm == r47 + r57 + r70 + r77
    L_OHm == r41 + r42 + r65 + r66 + r67 + k113 .* ρ_Ap .* ρ_OHm
    ∂ₜ(ρ_OHm) == maskalt(P_OHm - L_OHm)

    # CO3-
    P_CO3m == r48 + r58 + r74 + r75 + r82
    L_CO3m == r68 + r69 + r70 + r71 + r72 + k113 .* ρ_Ap .* ρ_CO3m
    ∂ₜ(ρ_CO3m) == maskalt(P_CO3m - L_CO3m)

    # CO4-
    P_CO4m == r54 + r63
    L_CO4m == r73 + r74 + r75 + r76 + k113 .* ρ_Ap .* ρ_CO4m
    ∂ₜ(ρ_CO4m) == maskalt(P_CO4m - L_CO4m)

    # NO2-
    P_NO2m == r49 + r50 + r55 + r66 + r71 + r80 + r81 + r83 + r84

    L_NO2m == r43 + r77 + r78 + r79 + k113 .* ρ_Ap .* ρ_NO2m
    ∂ₜ(ρ_NO2m) == maskalt(P_NO2m - L_NO2m)

    # NO3-
    P_NO3m == r59 + r60 + r72 + r78 + r79
    L_NO3m == r80 + r81 + k113 .* ρ_Ap .* ρ_NO3m
    ∂ₜ(ρ_NO3m) == maskalt(P_NO3m - L_NO3m)

    # O2-.NO
    P_O2mNO == r64 + r76
    L_O2mNO == r82 + r83 + r84 + k113 .* ρ_Ap .* ρ_O2mNO
    ∂ₜ(ρ_O2mNO) == maskalt(P_O2mNO - L_O2mNO)

    # HCO3-
    P_HCO3m == id(r67)
    L_HCO3m == k113 .* ρ_Ap .* ρ_HCO3m
    ∂ₜ(ρ_HCO3m) == maskalt(P_HCO3m - L_HCO3m)

    # N2+
    P_N2p == r1 + r16
    L_N2p == r91 + r92 + r93 + r101 + r102 + r103 + k113 .* ρ_N2p .* ρ_Bm
    ∂ₜ(ρ_N2p) == maskalt(P_N2p - L_N2p)

    # N+
    P_Np == id(r2)
    L_Np == r85 + r86 + r87 + k113 .* ρ_Np .* ρ_Bm
    ∂ₜ(ρ_Np) == maskalt(P_Np - L_Np)

    # O2+
    P_O2p== r5 + r8 + r22 + r85 + r89 + r91 + r98
    L_O2p == r94 + r95 + r96 + r97 + r110 + k113 .* ρ_O2p .* ρ_Bm
    ∂ₜ(ρ_O2p) == maskalt(P_O2p - L_O2p)

    # O+
    P_Op == r6 + r23 + r24 + r87
    L_Op == r88 + r89 + r90 + k113 .* ρ_Op .* ρ_Bm
    ∂ₜ(ρ_Op) == maskalt(P_Op - L_Op)

    # O4+
    P_O4p == id(r97)
    L_O4p == r98 + r99 + r111 + k113 .* ρ_O4p .* ρ_Bm
    ∂ₜ(ρ_O4p) == maskalt(P_O4p - L_O4p)

    # NO+
    P_NOp == r9 + r86 + r88 + r90 + r92 + r93 + r94 + r95 + r96
    L_NOp == r100 + r104 + r105 + r106 + r107 + r108 + r109 + k113 .* ρ_NOp .*
        ρ_Bm
    ∂ₜ(ρ_NOp) == maskalt(P_NOp - L_NOp)

    # Y+
    P_Yp == r99 + r100
    L_Yp == r112 + k113 .* ρ_Yp .* ρ_Bm
    ∂ₜ(ρ_Yp) == maskalt(P_Yp - L_Yp)




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
    ∂ₜ(ρ_O) == maskalt(P_O - L_O)



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

function rateCoefficents_dynam(rateCoef, theta, Tn, Te, useTe)

    # If useTe == 0, then the ambient solver is calling this script.
    # TODO: Look at this again after investigating the ambient solver.
    # For now, just call this as a parameter.
    if(useTe == 0)
        # TODO: use NaN instead of zeroes
        theta = zeros(size(Te))
        for index in 1:size((Te, 1))
            for j in 1:size((Te, 2))
                if(Te[index, j]) < 1
                    theta[index, j] = 65 .* Te[index, j]
                else
                    eqn = x -> 2*(x)^(2.6) - Te[index, j] * (x)^(2) - Te[index, j]
                    # TODO: Verify this works the same as matlab's fzero
                    theta[index, j] = 65 * find_zero(eqn, Te[index, j])
                end
            end
        end
    end

    # Electron temperature (eV)
    kB = 8.617e-5 # Boltzmann constant

    # If useTe == 1, then LIMA_Danny is calling this script.
    if(useTe == 1)
        index1 = (theta .<= 65)
        #Te .= 1.54e-2 * (theta .* index1)
        Te[index1] .= 1.54e-2 * (theta[index1])

        index2 = (theta .> 65)
        #Te .= 2 * ((theta .* index2)./65).^(2.6) ./ (1 .+ ((theta .* index2) ./65) .^2)
        Te[index2] .= 2 * ((theta[index2])./65).^(2.6) ./ (1 .+ ((theta[index2]) ./65) .^2)
        Te = max.(Te, kB * Tn)
    end
    
    # Electron temperature in Kelvin (K)
    TeK = Te ./ kB

    # Electron impact excitation, dissociation, and ionization

    rateCoef[10] = 10 .^ -(8.4 .+ 140 ./ theta)
    rateCoef[11] = 10 .^ -(8.2 .+ 148 ./ theta) + 10 .^ -(8.3 .+ 154 ./ theta) + 10 .^ -(8.7 .+ 168 ./ theta);
    rateCoef[12] = 10 .^ -(8.8 .+ 167 ./ theta) + 10 .^ -(8.5 .+ 174 ./ theta) + 10 .^ -(8.7 .+ 175 ./ theta);
    rateCoef[13] = 10 .^ -(8.2 .+ 211 ./ theta) + 10 .^ -(10.1 .+ 254 ./ theta) + 10 .^ -(9.2 .+ 262 ./ theta);

    rateCoef[14] = zeros(size(Te))
    index = theta .> 76
    rateCoef[14][index] = 2e-10 .* (3.096 - 6.71e-2.*(theta[index]) + 3e-4.*(theta[index].^2) + 1.59e-6.*(theta[index].^3) - 1.57e-9.*(theta[index].^4));


    rateCoef[15] = rateCoef[14] ./2;
    rateCoef[16] = 10 .^-(8.3 + 365 ./theta);

    rateCoef[17] = 10 .^-(9 + 52 ./theta);
    index = (theta .> 40);
    rateCoef[17][index] = 10 .^-(10.2 + 3.5 ./theta[index]);

    rateCoef[18] = 10 .^-(9.5 + 60 ./theta);
    index = (theta .> 30);
    rateCoef[18][index] = 10 .^-(11.2 + 7.2./theta[index]);

    rateCoef[19] = 10 .^-(7.9 + 134 ./theta);
    rateCoef[20] = 10 .^-(8 + 169 ./theta);
    rateCoef[21] = 10 .^-(8.8 + 119 ./theta);
    rateCoef[22] = 10 .^-(8.8 + 281 ./theta);
    rateCoef[23] = 3.18e-14 * (TeK ./300) .* exp(-206055 ./TeK);
    rateCoef[24] = 7.10e-11 * Te .^(0.5) .* exp(-17 ./Te); # TODO: Order of ops here?

    index = (Te .< 0.1);
    rateCoef[10][index] = 0; # TODO: These might have to be .=
    rateCoef[11][index] = 0;
    rateCoef[12][index] = 0;
    rateCoef[13][index] = 0;
    rateCoef[14][index] = 0;
    rateCoef[15][index] = 0;
    rateCoef[16][index] = 0;
    rateCoef[17][index] = 0;
    rateCoef[18][index] = 0;
    rateCoef[19][index] = 0;
    rateCoef[20][index] = 0;
    rateCoef[21][index] = 0;
    rateCoef[22][index] = 0;
    rateCoef[23][index] = 0;
    rateCoef[24][index] = 0;

    # Electron attachment
    rateCoef[25] = 1.07e-31 .* (300 ./TeK).^2 .* exp(-70 ./TeK) .* exp(1500 .*(TeK-Tn)./TeK./Tn);
    rateCoef[26] = 1.4e-29 .* (300 ./TeK) .* exp(-600 ./TeK) .* exp(700 .*(TeK-Tn)./TeK./Tn);
    rateCoef[27] = 1.07e-9 .* Te .^(-1.39) .* exp(-6.26./Te);

    rateCoef[28] = 1e-11 .* ones(size(Tn));
    index = (Te .> 0.13);
    rateCoef[28][index] = 2.12e-9 .* Te[index].^(-1.06) .* exp(-0.93./Te[index]);

    rateCoef[29] = 1e-9 .* ones(size(Tn));
    index = (Te .> 0.14);
    rateCoef[29][index] = 9.76e-8 .* Te[index].^(-1.309) .* exp(-1.007./Te[index]);

    # Electron detachment
    rateCoef[30] = 5.47e-8 * Te.^(0.324) .* exp(-2.98./Te);

    rateCoef[31] = 10 .^(-15.66 + 2.97.*log10.(theta) - 0.58.*(log10.(theta).^2));
    index = (Te .< 0.1);
    rateCoef[31][index] = 0;

    # Electron-ion recombination
    rateCoef[101] = 0.37 * 2.2e-7 * (300 ./Tn).^(0.2) .* (TeK./Tn).^(-0.39);
    rateCoef[102] = 0.11 * 2.2e-7 * (300 ./Tn).^(0.2) .* (TeK./Tn).^(-0.39);
    rateCoef[103] = 0.52 * 2.2e-7 * (300 ./Tn).^(0.2) .* (TeK./Tn).^(-0.39);

    rateCoef[104] = 0.05 * 3.5e-7 .* (300 ./Tn) .* (TeK./Tn).^(-0.69);
    # This call to index is the same for 104 through 109
    index = (Te .> 1.25);
    rateCoef[104][index] = 0;

    rateCoef[105] = 0.95 * 3.5e-7 .* (300 ./Tn) .* (TeK./Tn).^(-0.69);
    #index = (Te .> 1.25);

    rateCoef[106] = zeros(size(Tn));
    rateCoef[107] = zeros(size(Tn));
    rateCoef[108] = zeros(size(Tn));
    rateCoef[109] = zeros(size(Tn));

    rateCoef[106][index] = 0.10 * 3.5e-7 .* (300 ./Tn(index)) .* (TeK(index)./Tn(index)).^(-0.69);
    rateCoef[107][index] = 0.10 * 3.5e-7 .* (300 ./Tn(index)) .* (TeK(index)./Tn(index)).^(-0.69);
    rateCoef[108][index] = 0.70 * 3.5e-7 .* (300 ./Tn(index)) .* (TeK(index)./Tn(index)).^(-0.69);
    rateCoef[109][index] = 0.10 * 3.5e-7 .* (300 ./Tn(index)) .* (TeK(index)./Tn(index)).^(-0.69);

    rateCoef[105][index] = 0;
    rateCoef[106][index] = 0;
    rateCoef[107][index] = 0;
    rateCoef[108][index] = 0;
    rateCoef[109][index] = 0;

    rateCoef[110] = 1.95e-7 * (300 ./Tn).^(0.7) .* (TeK./Tn).^(-0.7);
    rateCoef[111] = 4.2e-6 * (TeK./Tn).^(-0.48);
    rateCoef[112] = 6.5e-6 * (300 ./TeK).^(1/2);   # Recombination rate for cluster order n = 3

    # Note: red_vec is always empty
    # Note: This never runs
    #d_array = [10:31, 101:112]
    #d_array = reduce(vcat, collect.(d_array))
    #fields = fieldnames(rateCoef);
    #if(!isempty(red_vec))
    #    for i = 1:length(red_vec)
    #        reduce = find(red_vec(i) == d_array);
    #        red_ind = 130+reduce; # 130 is the number of static rate coefficients
    #        rateCoef.(fields{red_ind}) = zeros(size(rateCoef.(fields{red_ind})));
    #    end
    #end

    # Re-pasted back everything below here.

    # Total gas density (dominated by N2 and O2)
    ρ_M = ρ_N2 + ρ_O2

    # TODO: Maybe break down each class of rate computation into a composition pattern.
    # Cosmic ray and photo-ionization
    rate[1] = rateCoef[1] .* ρ_N2
    rate[2] = rateCoef[2] .* ρ_N2
    rate[3] = rateCoef[3] .* ρ_N2
    rate[4] = rateCoef[4] .* ρ_N2
    rate[5] = rateCoef[5] .* ρ_O2
    rate[6] = rateCoef[6] .* ρ_O2
    rate[7] = rateCoef[7] .* ρ_O2
    rate[8] = rateCoef[8] .* ρ_O2
    rate[9] = rateCoef[9] .* ρ_NO

    # Electron impact excitation, dissociation, and ionization
    rate[10] = rateCoef[10] .* ρ_e .* ρ_N2
    rate[11] = rateCoef[11] .* ρ_e .* ρ_N2
    rate[12] = rateCoef[12] .* ρ_e .* ρ_N2
    rate[13] = rateCoef[13] .* ρ_e .* ρ_N2
    rate[14] = rateCoef[14] .* ρ_e .* ρ_N2
    rate[15] = rateCoef[15] .* ρ_e .* ρ_N2
    rate[16] = rateCoef[16] .* ρ_e .* ρ_N2
    rate[17] = rateCoef[17] .* ρ_e .* ρ_O2
    rate[18] = rateCoef[18] .* ρ_e .* ρ_O2
    rate[19] = rateCoef[19] .* ρ_e .* ρ_O2
    rate[20] = rateCoef[20] .* ρ_e .* ρ_O2
    rate[21] = rateCoef[21] .* ρ_e .* ρ_O2
    rate[22] = rateCoef[22] .* ρ_e .* ρ_O2
    rate[23] = rateCoef[23] .* ρ_e .* ρ_O2
    rate[24] = rateCoef[24] .* ρ_e .* ρ_O2

    # Electron attachment
    rate[25] = rateCoef[25] .*ρ_e .* ρ_O2 .* ρ_N2
    rate[26] = rateCoef[26] .*ρ_e .* ρ_O2 .* ρ_O2
    rate[27] = rateCoef[27] .*ρ_e .* ρ_O2
    rate[28] = rateCoef[28] .*ρ_e .* ρ_O3
    rate[29] = rateCoef[29] .*ρ_e .* ρ_O3

    # Electron detachment
    rate[30] = rateCoef[30] .* ρ_Om .* ρ_e
    rate[31] = rateCoef[31] .* ρ_Om .* ρ_N2
    rate[32] = rateCoef[32] .* ρ_Om .* ρ_N4S
    rate[33] = rateCoef[33] .* ρ_Om .* ρ_O2
    rate[34] = rateCoef[34] .* ρ_Om .* ρ_O2a
    rate[35] = rateCoef[35] .* ρ_Om .* ρ_O
    rate[36] = rateCoef[36] .* ρ_Om .* ρ_NO
    rate[37] = rateCoef[37] .* ρ_O2m .* ρ_O2a
    rate[38] = rateCoef[38] .* ρ_O2m .* ρ_O
    rate[39] = rateCoef[39] .* ρ_O3m .* ρ_O
    rate[40] = rateCoef[40] .* ρ_O3m .* ρ_O3
    rate[41] = rateCoef[41] .* ρ_OHm .* ρ_O
    rate[42] = rateCoef[42] .* ρ_OHm .* ρ_H
    rate[43] = rateCoef[43] .* ρ_NO2m .* ρ_O

    # Negative ion conversion
    rate[44] = rateCoef[44] .* ρ_Om .* ρ_O2a
    rate[45] = rateCoef[45] .* ρ_Om .* ρ_O2 .* ρ_M
    rate[46] = rateCoef[46] .* ρ_Om .* ρ_O3
    rate[47] = rateCoef[47] .* ρ_Om .* ρ_H2O
    rate[48] = rateCoef[48] .* ρ_Om .* ρ_CO2 .* ρ_M
    rate[49] = rateCoef[49] .* ρ_Om .* ρ_NO .* ρ_M
    rate[50] = rateCoef[50] .* ρ_Om .* ρ_NO2
    rate[51] = rateCoef[51] .* ρ_O2m .* ρ_O
    rate[52] = rateCoef[52] .* ρ_O2m .* ρ_O3
    rate[53] = rateCoef[53] .* ρ_O2m .* ρ_O2 .* ρ_M
    rate[54] = rateCoef[54] .* ρ_O2m .* ρ_CO2 .* ρ_M
    rate[55] = rateCoef[55] .* ρ_O2m .* ρ_NO2
    rate[56] = rateCoef[56] .* ρ_O3m .* ρ_O
    rate[57] = rateCoef[57] .* ρ_O3m .* ρ_H
    rate[58] = rateCoef[58] .* ρ_O3m .* ρ_CO2
    rate[59] = rateCoef[59] .* ρ_O3m .* ρ_NO
    rate[60] = rateCoef[50] .* ρ_O3m .* ρ_NO2
    rate[61] = rateCoef[61] .* ρ_O4m .* ρ_O2a
    rate[62] = rateCoef[62] .* ρ_O4m .* ρ_O
    rate[63] = rateCoef[63] .* ρ_O4m .* ρ_CO2
    rate[64] = rateCoef[64] .* ρ_O4m .* ρ_NO
    rate[65] = rateCoef[65] .* ρ_OHm .* ρ_O3
    rate[66] = rateCoef[66] .* ρ_OHm .* ρ_NO2
    rate[67] = rateCoef[67] .* ρ_OHm .* ρ_CO2 .* ρ_M
    rate[68] = rateCoef[68] .* ρ_CO3m .* ρ_O
    rate[69] = rateCoef[69] .* ρ_CO3m .* ρ_O2
    rate[70] = rateCoef[70] .* ρ_CO3m .* ρ_H
    rate[71] = rateCoef[71] .* ρ_CO3m .* ρ_NO
    rate[72] = rateCoef[72] .* ρ_CO3m .* ρ_NO2
    rate[73] = rateCoef[73] .* ρ_CO4m .* ρ_O3
    rate[74] = rateCoef[74] .* ρ_CO4m .* ρ_O
    rate[75] = rateCoef[75] .* ρ_CO4m .* ρ_H
    rate[76] = rateCoef[76] .* ρ_CO4m .* ρ_NO
    rate[77] = rateCoef[77] .* ρ_NO2m .* ρ_H
    rate[78] = rateCoef[78] .* ρ_NO2m .* ρ_O3
    rate[79] = rateCoef[79] .* ρ_NO2m .* ρ_NO2
    rate[80] = rateCoef[80] .* ρ_NO3m .* ρ_O
    rate[81] = rateCoef[81] .* ρ_NO3m .* ρ_O3
    rate[82] = rateCoef[82] .* ρ_O2mNO .* ρ_CO2
    rate[83] = rateCoef[83] .* ρ_O2mNO .* ρ_NO
    rate[84] = rateCoef[84] .* ρ_O2mNO .* ρ_H

    # Positive ion conversion
    rate[85] = rateCoef[85] .* ρ_Np .* ρ_O2
    rate[86] = rateCoef[86] .* ρ_Np .* ρ_O2
    rate[87] = rateCoef[87] .* ρ_Np .* ρ_O2
    rate[88] = rateCoef[88] .* ρ_Op .* ρ_N2
    rate[89] = rateCoef[89] .* ρ_Op .* ρ_O2
    rate[90] = rateCoef[90] .* ρ_Op .* ρ_N2 .* ρ_M
    rate[91] = rateCoef[91] .* ρ_N2p .* ρ_O2
    rate[92] = rateCoef[92] .* ρ_N2p .* ρ_O2
    rate[93] = rateCoef[93] .* ρ_N2p .* ρ_O
    rate[94] = rateCoef[94] .* ρ_O2p .* ρ_N2
    rate[95] = rateCoef[95] .* ρ_O2p .* ρ_N4S
    rate[96] = rateCoef[96] .* ρ_O2p .* ρ_NO
    rate[97] = rateCoef[97] .* ρ_O2p .* ρ_O2 .* ρ_M
    rate[98] = rateCoef[98] .* ρ_O4p .* ρ_O
    rate[99] = rateCoef[99] .* ρ_O4p .* ρ_H2O
    rate[100] = rateCoef[100] .* ρ_NOp .* ρ_M .* ρ_M

    # Electrion-ion recombination
    rate[101] = rateCoef[101] .* ρ_N2p .* ρ_e
    rate[102] = rateCoef[102] .* ρ_N2p .* ρ_e
    rate[103] = rateCoef[103] .* ρ_N2p .* ρ_e
    rate[104] = rateCoef[104] .* ρ_NOp .* ρ_e
    rate[105] = rateCoef[105] .* ρ_NOp .* ρ_e
    rate[106] = rateCoef[106] .* ρ_NOp .* ρ_e
    rate[107] = rateCoef[107] .* ρ_NOp .* ρ_e
    rate[108] = rateCoef[108] .* ρ_NOp .* ρ_e
    rate[109] = rateCoef[109] .* ρ_NOp .* ρ_e
    rate[110] = rateCoef[110] .* ρ_O2p .* ρ_e
    rate[111] = rateCoef[111] .* ρ_O4p .* ρ_e
    rate[112] = rateCoef[112] .* ρ_Yp .* ρ_e

    # Ion-ion recombination
    # r113 determined in rateSolver.m

    # Active-state N2 chemistry
    rate[114] = rateCoef[114] .* ρ_N2A .* ρ_O2
    rate[115] = rateCoef[115] .* ρ_N2A .* ρ_O2
    rate[116] = rateCoef[116] .* ρ_N2A .* ρ_O2
    rate[117] = rateCoef[117] .* ρ_N2A .* ρ_O
    rate[118] = rateCoef[118] .* ρ_N2A .* ρ_N4S
    rate[119] = rateCoef[119] .* ρ_N2B
    rate[120] = rateCoef[120] .* ρ_N2B .* ρ_N2
    rate[121] = rateCoef[121] .* ρ_N2B .* ρ_O2
    rate[122] = rateCoef[122] .* ρ_N2a .* ρ_N2
    rate[123] = rateCoef[123] .* ρ_N2a .* ρ_O2
    rate[124] = rateCoef[124] .* ρ_N2a .* ρ_NO
    rate[125] = rateCoef[125] .* ρ_N2C
    rate[126] = rateCoef[126] .* ρ_N2C .* ρ_N2
    rate[127] = rateCoef[127] .* ρ_N2C .* ρ_O2

    # Active-state O2 chemistry
    rate[128] = rateCoef[128] .* ρ_O2a
    rate[129] = rateCoef[129] .* ρ_O2a .* ρ_N2
    rate[130] = rateCoef[130] .* ρ_O2a .* ρ_O2
    rate[131] = rateCoef[131] .* ρ_O2a .* ρ_N4S
    rate[132] = rateCoef[132] .* ρ_O2a .* ρ_O
    rate[133] = rateCoef[133] .* ρ_O2a .* ρ_NO
    rate[134] = rateCoef[134] .* ρ_O2b
    rate[135] = rateCoef[135] .* ρ_O2b .* ρ_N2
    rate[136] = rateCoef[136] .* ρ_O2b .* ρ_O2
    rate[137] = rateCoef[137] .* ρ_O2b .* ρ_O
    rate[138] = rateCoef[138] .* ρ_O2b .* ρ_O3

    # Odd nitrogen chemistry
    rate[139] = rateCoef[139] .* ρ_N4S .* ρ_O .* ρ_M
    rate[140] = rateCoef[140] .* ρ_N4S .* ρ_O3
    rate[141] = rateCoef[141] .* ρ_N4S .* ρ_O2
    rate[142] = rateCoef[142] .* ρ_N2D .* ρ_O2
    rate[143] = rateCoef[143] .* ρ_N4S .* ρ_NO
    rate[144] = rateCoef[144] .* ρ_N2D .* ρ_NO
    rate[145] = rateCoef[145] .* ρ_NO .* ρ_O
    rate[146] = rateCoef[146] .* ρ_NO .* ρ_O .* ρ_M
    rate[147] = rateCoef[147] .* ρ_NO .* ρ_O3
    rate[148] = rateCoef[148] .* ρ_N4S .* ρ_NO2
    rate[149] = rateCoef[149] .* ρ_N4S .* ρ_NO2
    rate[150] = rateCoef[150] .* ρ_N4S .* ρ_NO2
    rate[151] = rateCoef[151] .* ρ_N4S .* ρ_NO2
    rate[152] = rateCoef[152] .* ρ_O .* ρ_NO2
    rate[153] = rateCoef[153] .* ρ_N2D
    rate[154] = rateCoef[154] .* ρ_N2D .* ρ_N2
    rate[155] = rateCoef[155] .* ρ_N2D .* ρ_O

    # Odd oxygen chemistry
    rate[156] = rateCoef[156] .* ρ_O .* ρ_O .* ρ_M
    rate[157] = rateCoef[157] .* ρ_O .* ρ_O .* ρ_M
    rate[158] = rateCoef[158] .* ρ_O .* ρ_O2 .* ρ_N2
    rate[159] = rateCoef[159] .* ρ_O .* ρ_O2 .* ρ_O2
    rate[160] = rateCoef[160] .* ρ_O .* ρ_O3
    rate[161] = rateCoef[161] .* ρ_O .* ρ_OH 
    rate[162] = rateCoef[162] .* ρ_O .* ρ_HO2
    rate[163] = rateCoef[163] .* ρ_OH .* ρ_OH
    rate[164] = rateCoef[164] .* ρ_OH .* ρ_O3
    rate[165] = rateCoef[165] .* ρ_H .* ρ_HO2
    rate[166] = rateCoef[166] .* ρ_H .* ρ_O3
    rate[167] = rateCoef[167] .* ρ_HO2 .* ρ_O3

    # Positively charged ions are grouped.
    ρ_Ap = ρ_N2p + ρ_Np + ρ_O2p + ρ_Op + ρ_O4p + ρ_NOp + ρ_Yp

    # Negatively charged ions are grouped.
    ρ_Bm = ρ_Om + ρ_O2m + ρ_O3m + ρ_O4m + ρ_OHm + ρ_CO3m + ρ_CO4m + ρ_NO2m +
        ρ_NO3m + ρ_O2mNO + ρ_HCO3m

    # N2(A)
    P_N2A = r10 + r119 + r120
    L_N2A = r114 + r115 + r116 + r117 + r118
    ∂ₜ_ρ_N2A = P_N2A - L_N2A
    # TODO: Either hard code the above line for the other 38 species or use a
    # composition pattern.

    # N2(B)
    P_N2B = r11 + r122 + r125
    L_N2B = r119 + r120 + r121
    ∂ₜ_ρ_N2B = P_N2B - L_N2B

    # N2(a')
    P_N2a = r12 + r126
    L_N2a = r122 + r123 + r124
    ∂ₜ_ρ_N2a = P_N2a - L_N2a
        
    # N2(C)
    P_N2C = (r13)
    L_N2C = r125 + r126 + r127
    ∂ₜ_ρ_N2C = P_N2C - L_N2C

    # N(4S)
    P_N4S = r3 + 2 .*r14 + r15 + r85 + r88 + r90 + r101 + r102 + r104 + r106 +
        r107 + r124 + r153 + r154 + r155

    L_N4S = r32 + r95 + r118 + r131 + r139 + r140 + r141 + r143 + r148 + r149 +
        r150 + r151
    ∂ₜ_ρ_N4S = P_N4S - L_N4S

    # N*
    P_N2D = r4 + r15 + r93 + r101 + r102 + 2 .*r103 + r105 + r108 + r109 + r117

    L_N2D = r142 + r144 + r153 + r154 + r155
    ∂ₜ_ρ_N2D = P_N2D - L_N2D

    # O2(a)
    P_O2a = r17 + r114 + r135 + r136 + r137 + r157

    L_O2a = r34 + r37 + r44 + r61 + r128 + r129 + r130 + r131 + r132 + r133
    ∂ₜ_ρ_O2a = P_O2a - L_O2a

    # O2(b)
    P_O2b = r18 + r115
    L_O2b = r134 + r135 + r136 + r137 + r138
    ∂ₜ_ρ_O2b = P_O2b - L_O2b

    # O  see bottom for details
    #P_O == 0   
    #L_O == 0

    # NO
    P_NO = r32 + r77 + r79 + r87 + r92 + r94 + r117 + r131 + r139 + r140 +
        r141 + r142 + 2 .*r151 + r152

    L_NO = r9 + r36 + r49 + r59 + r64 + r71 + r72 + r76 + r83 + r96 + r124 +
        r143 + r144 + r145 + r146 + r147 
    ∂ₜ_ρ_NO = P_NO - L_NO

    # NO2
    P_NO2 = r36 + r43 + r82 + r83 + r145 + r146 + r147

    L_NO2 = r50 + r55 + r60 + r66 + r79 + r148 + r149 + r150 + r151 + r152
    ∂ₜ_ρ_NO2 = P_NO2 - L_NO2

    # e-
    P_e = r1 + r2 + r5 + r6 + r8 + r9 + r16 + r22 + r23 + r30 + r31 + r32 +
        r33 + r34 + r35 + r36 + r37 + r38 + r39 + r40 + r41 + r42 + r43

    L_e = r25 + r26 + r27 + r28 + r29 + r101 + r102 + r103 + r104 + r105 +
        r106 + r107 + r108 + r109 + r110 + r111 + r112
    ∂ₜ_ρ_e = P_e - L_e

    # O-
    P_Om = r24 + r27 + r28 + r51
    L_Om = r30 + r31 + r32 + r33 + r34 + r35 + r36 + r44 + r45 + r46 + r47 +
        r48 + r49 + r50 + rateCoef[113] .* ρ_Ap .* ρ_Om
    ∂ₜ_ρ_Om = P_Om - L_Om

    # O2-
    P_O2m = r25 + r26 + r29 + r44 + r56 + r61 + r68

    L_O2m = r37 + r38 + r51 + r52 + r53 + r54 + r55 + rateCoef[113] .* ρ_Ap .* ρ_O2m
    ∂ₜ_ρ_O2m = P_O2m - L_O2m

    # O3-
    P_O3m = r45 + r46 + r52 + r62 + r65 + r69 + r73

    L_O3m = r39 + r40 + r56 + r57 + r58 + r59 + r60 + rateCoef[113] .* ρ_Ap .* ρ_O3m
    ∂ₜ_ρ_O3m = P_O3m - L_O3m

    # O4-
    P_O4m = (r53)
    L_O4m = r61 + r62 + r63 + r64 + rateCoef[113] .* ρ_Ap .* ρ_O4m
    ∂ₜ_ρ_O4m = P_O4m - L_O4m

    # OH-
    P_OHm = r47 + r57 + r70 + r77
    L_OHm = r41 + r42 + r65 + r66 + r67 + rateCoef[113] .* ρ_Ap .* ρ_OHm
    ∂ₜ_ρ_OHm = P_OHm - L_OHm

    # CO3-
    P_CO3m = r48 + r58 + r74 + r75 + r82
    L_CO3m = r68 + r69 + r70 + r71 + r72 + rateCoef[113] .* ρ_Ap .* ρ_CO3m
    ∂ₜ_ρ_CO3m = P_CO3m - L_CO3m

    # CO4-
    P_CO4m = r54 + r63
    L_CO4m = r73 + r74 + r75 + r76 + rateCoef[113] .* ρ_Ap .* ρ_CO4m
    ∂ₜ_ρ_CO4m = P_CO4m - L_CO4m

    # NO2-
    P_NO2m = r49 + r50 + r55 + r66 + r71 + r80 + r81 + r83 + r84

    L_NO2m = r43 + r77 + r78 + r79 + rateCoef[113] .* ρ_Ap .* ρ_NO2m
    ∂ₜ_ρ_NO2m = P_NO2m - L_NO2m

    # NO3-
    P_NO3m = r59 + r60 + r72 + r78 + r79
    L_NO3m = r80 + r81 + rateCoef[113] .* ρ_Ap .* ρ_NO3m
    ∂ₜ_ρ_NO3m = P_NO3m - L_NO3m

    # O2-.NO
    P_O2mNO = r64 + r76
    L_O2mNO = r82 + r83 + r84 + rateCoef[113] .* ρ_Ap .* ρ_O2mNO
    ∂ₜ_ρ_O2mNO = P_O2mNO - L_O2mNO

    # HCO3-
    P_HCO3m = (r67)
    L_HCO3m = rateCoef[113] .* ρ_Ap .* ρ_HCO3m
    ∂ₜ_ρ_HCO3m = P_HCO3m - L_HCO3m

    # N2+
    P_N2p = r1 + r16
    L_N2p = r91 + r92 + r93 + r101 + r102 + r103 + rateCoef[113] .* ρ_N2p .* ρ_Bm
    ∂ₜ_ρ_N2p = P_N2p - L_N2p

    # N+
    P_Np = (r2)
    L_Np = r85 + r86 + r87 + rateCoef[113] .* ρ_Np .* ρ_Bm
    ∂ₜ_ρ_Np = P_Np - L_Np

    # O2+
    P_O2p= r5 + r8 + r22 + r85 + r89 + r91 + r98
    L_O2p = r94 + r95 + r96 + r97 + r110 + rateCoef[113] .* ρ_O2p .* ρ_Bm
    ∂ₜ_ρ_O2p = P_O2p - L_O2p

    # O+
    P_Op = r6 + r23 + r24 + r87
    L_Op = r88 + r89 + r90 + rateCoef[113] .* ρ_Op .* ρ_Bm
    ∂ₜ_ρ_Op = P_Op - L_Op

    # O4+
    P_O4p = (r97)
    L_O4p = r98 + r99 + r111 + rateCoef[113] .* ρ_O4p .* ρ_Bm
    ∂ₜ_ρ_O4p = P_O4p - L_O4p

    # NO+
    P_NOp = r9 + r86 + r88 + r90 + r92 + r93 + r94 + r95 + r96
    L_NOp = r100 + r104 + r105 + r106 + r107 + r108 + r109 + rateCoef[113] .* ρ_NOp .*
        ρ_Bm
    ∂ₜ_ρ_NOp = P_NOp - L_NOp

    # Y+
    P_Yp = r99 + r100
    L_Yp = r112 + rateCoef[113] .* ρ_Yp .* ρ_Bm
    ∂ₜ_ρ_Yp = P_Yp - L_Yp




    # ATOMIC OXYGEN production and loss rates depend heavily upon the odd
    # hydrogen family, which are not provided for by the NRLMSISE modeL_  Given
    # that its production via electric field heating is small, the density of 
    # atomic oxygen may safetly be considered to be constant throughout the 
    # simulation.  NOTE: This approximation must be verified, and can be done
    # so most readily by comparison with electron-impact dissociation of N2 .

    P_O = r7 + 2 .*r19 + 2 .*r20 + 2 .*r21 + r23 + r27 + r29 + r30 + r44 +
        r46 + r50 + r59 + r86 + r89 + r95 + r104 + r105 + r106 + r107 + r108 +
        r109 + 2 .*r110 + 2 .*r116 + 2 .*r121 + 2 .*r123 + r124 + 2 .*r127 +
        r131 + r138 + r141 + r142 + r143 + r144 + r148 + 2 .*r149 + r163 + r165
    
    L_O = r35 + r39 + r41 + r51 + r56 + r62 + r68 + r74 + r80 + r93 +
        r117 + r139 + r145 + r146 + r152 + r156 + r157 + r158 + r159 + r160 +
        r161 + r162
    ∂ₜ_ρ_O = P_O - L_O



    # POSITIVE NITROGEN ION, N+, densities depend only upon production via GCR,
    # and loss by charge exchange with O2 .  Given that the production and loss
    # rates are entirely unaffected by electric field heating (dissociative
    # ionization of N2 is not included), the density of N+ remains at its 
    # equilibrium value throughout the simulation.  Due to its volatile nature 
    # (large time-rate of change relative to its low density), can be tarateCoef[en as
    # constant for heating simulations.
    
    # P_Np = r2
    # 
    # L_Np = r85 + r86 + r87 + rateCoef[113] .* ρ_Np .* ρ_Bm
end
