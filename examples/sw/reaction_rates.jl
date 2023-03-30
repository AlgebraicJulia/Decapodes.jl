
"""
    These reaction rates are used to define production and loss rates for each species.
"""
ReactionRates = SummationDecapode(parse_decapode(quote
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
end))

# All rates are defined as 0-Forms.
ReactionRates[:type] .= :Form0

