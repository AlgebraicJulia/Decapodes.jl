# TODO: Write this file as a Decapode

rateCoef = Any[i for i in 1:167]

RateCoefficientsDynamic = SummationDecapode(parse_decapode(quote
    kB::Constant{X}

    Te_Less_65 == θ .≤ 65
    Te_Greater_65 == θ .> 65

    # See Kotovsky pp. 92 Eq. 5-1
    Te_unthresholded == (Te_Less_65 * 1.54e-2 * θ) +
        (Te_Greater_65 * 2 * (θ ./ 65).^(2.6) ./ (1 .+ (θ ./65) .^2))

    Te_Greater_kB_times_Tn == Te_unthresholded .> (kB * Tn)
    Te_Less_kB_times_Tn == Te_unthresholded .≤ (kB * Tn)
    Te == (Te_Greater_kB_times_Tn * Te_unthresholded) +
        (Te_Less_kB_times_Tn * (kB * Tn))

    # Electron temperature [K]
    TeK == Te ./ kB

    # Electron impact excitation, dissociation, and ionization
    # TODO: Annotate each equation with the paper it came from. See Kotovsky Table B-1 references
    # Note that some papers, like Kossyi et al. (92) have some equations cited from sources that Kotovsky did not cite, like for the equation of ion-ion recombination.
    # TODO: Check whether this should be .> or .<
    k10_through_24_mask == Te .> 0.1
    # TODO: Check if θ can ever be zero.
    k10 == k10_through_24_mask * (10 .^ -(8.4 .+ 140 ./ θ))
    k11 == k10_through_24_mask * (10 .^ -(8.2 .+ 148 ./ θ) + 10 .^ -(8.3 .+ 154 ./ θ) + 10 .^ -(8.7 .+ 168 ./ θ))
    k12 == k10_through_24_mask * (10 .^ -(8.8 .+ 167 ./ θ) + 10 .^ -(8.5 .+ 174 ./ θ) + 10 .^ -(8.7 .+ 175 ./ θ))
    k13 == k10_through_24_mask * (10 .^ -(8.2 .+ 211 ./ θ) + 10 .^ -(10.1 .+ 254 ./ θ) + 10 .^ -(9.2 .+ 262 ./ θ))

    k14_mask == θ .> 76
    k14 == k10_through_24_mask * k14_mask * (2e-10 .* (3.096 .- 6.71e-2.*θ + 3e-4.*(θ .^2) + 1.59e-6.*(θ .^3) .- 1.57e-9.*(θ .^4)))

    k15 == k10_through_24_mask * (k14 ./2)
    k16 == k10_through_24_mask * (10 .^ -(8.3 + 365 ./θ))

    k17_mask == θ .> 40
    k17 == k10_through_24_mask * (
        k17_mask * (10 .^ -(10.2 + 3.5 ./θ)) +
        invert_mask(k17_mask) * (10 .^ -(9 + 52 ./θ)))

    k18_mask == θ .> 30
    k18 == k10_through_24_mask * (
        k18_mask * ( 10 .^ -(11.2 + 7.2./θ)) +
        invert_mask(k18_mask) * (10 .^ -(9.5 + 60 ./θ)))

    k19 == k10_through_24_mask * (10 .^ -(7.9 + 134 ./θ))
    k20 == k10_through_24_mask * (10 .^ -(8 + 169 ./θ))
    k21 == k10_through_24_mask * (10 .^ -(8.8 + 119 ./θ))
    k22 == k10_through_24_mask * (10 .^ -(8.8 + 281 ./θ))
    k23 == k10_through_24_mask * (3.18e-14 * (TeK ./300) .* exp(-206055 ./TeK))
    k24 == k10_through_24_mask * (7.10e-11 * Te .^(0.5) .* exp(-17 ./Te))

    # Electron attachment
    k25 == 1.07e-31 .* (300 ./TeK).^2 .* exp(-70 ./TeK) .* exp(1500 .*(TeK .-Tn)./TeK./Tn)
    k26 == 1.4e-29 .* (300 ./TeK) .* exp(-600 ./TeK) .* exp(700 .*(TeK .-Tn)./TeK./Tn)
    k27 == 1.07e-9 .* Te .^ (-1.39) .* exp(-6.26 ./Te)

    k28_mask == Te .> 0.13
    k28 == k28_mask * (2.12e-9 .* Te .^ (-1.06) .* exp(-0.93 ./Te)) +
        invert_mask(k28_mask) * 1e-11

    k29_mask == Te .> 0.14
    k29 == k29_mask * (9.76e-8 .* Te .^ (-1.309) .* exp(-1.007 ./Te)) +
        invert_mask(k29_mask) * 1e-9

    # Electron detachment
    k30 == 5.47e-8 * Te .^ (0.324) .* exp(-2.98 ./Te)

    k31_mask == Te .≥ 0.1
    k31 == k31_mask * (10 .^ (-15.66 + 2.97 .* log10.(θ) - 0.58 .* (log10.(θ) .^ 2)))

    # Electron-ion recombination
    k101 == 0.37 * 2.2e-7 * (300 ./Tn).^(0.2) .* (TeK./Tn).^(-0.39)
    k102 == 0.11 * 2.2e-7 * (300 ./Tn).^(0.2) .* (TeK./Tn).^(-0.39)
    k103 == 0.52 * 2.2e-7 * (300 ./Tn).^(0.2) .* (TeK./Tn).^(-0.39)

    # Observe: k104 and k105 are set to 0 where Te > 1.25, but k106, k107, k108,
    # and k109 are set to 0 where Te ≤ 1.25.
    k104_through_105_mask == Te .≤ 1.25
    k106_through_109_mask == Te .> 1.25
    
    k104 == k104_through_105_mask * (0.05 * 3.5e-7 .* (300 ./Tn) .* (TeK./Tn).^(-0.69))

    k105 == k104_through_105_mask * (0.95 * 3.5e-7 .* (300 ./Tn) .* (TeK./Tn).^(-0.69))

    k106 == k106_through_109_mask * (0.10 * 3.5e-7 .* (300 ./Tn) .* (TeK./Tn).^(-0.69))

    k107 == k106_through_109_mask * (0.10 * 3.5e-7 .* (300 ./Tn) .* (TeK./Tn).^(-0.69))

    k108 == k106_through_109_mask * (0.70 * 3.5e-7 .* (300 ./Tn) .* (TeK./Tn).^(-0.69))

    k109 == k106_through_109_mask * (0.10 * 3.5e-7 .* (300 ./Tn) .* (TeK./Tn).^(-0.69))

    k110 == 1.95e-7 * (300 ./Tn).^(0.7) .* (TeK./Tn).^(-0.7)
    k111 == 4.2e-6 * (TeK./Tn).^(-0.48)
    k112 == 6.5e-6 * (300 ./TeK).^(1/2) # Recombination rate for cluster order n = 3
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
end
