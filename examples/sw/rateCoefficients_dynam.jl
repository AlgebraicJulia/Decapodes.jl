rateCoef = Any[i for i in 1:167]

function rateCoefficents_dynam(rateCoef, theta, Tn, Te, useTe)

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
    index = (Te .> 1.25);
    rateCoef[104][index] = 0;

    rateCoef[105] = 0.95 * 3.5e-7 .* (300 ./Tn) .* (TeK./Tn).^(-0.69);
    # This call to index is the same for 104 through 109
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

    d_array = [10:31, 101:112]
    d_array = reduce(vcat, collect.(d_array))
    #fields = fieldnames(rateCoef);
    # TODO: red_vec is always empty
    # TODO: This never runs
    if(!isempty(red_vec))
        for i = 1:length(red_vec)
            reduce = find(red_vec(i) == d_array);
            red_ind = 130+reduce; # 130 is the number of static rate coefficients
            rateCoef.(fields{red_ind}) = zeros(size(rateCoef.(fields{red_ind})));
        end
    end
end
