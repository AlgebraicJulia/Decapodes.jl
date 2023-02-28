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
                    egn = x -> 2*(x)^(2.6) - Te[index, j] * (x)^(2) - Te[index, j]
                    # TODO: Solve for the zeroes when Te[index, j] is plugged In

                end
            end
        end
    end

    # Boltzmann constant
    kB = 8.617e-5

    if(useTe == 1)
        index1 = (theta .<= 65)
        Te .= 1.54e-2 * (theta .* index1)

        index2 = (theta .> 65)
        Te .= 2 * ((theta .* index2)./65).^(2.6) ./ (1 .+ ((theta .* index2) ./65) .^2)
        Te = max.(Te, kB * Tn)
    end
    
    TeK = Te ./ kB

    rateCoef[10] = 10 .^ -(8.4 .+ 140 ./ theta)
    rateCoef[11] = 10 .^ -(8.2 .+ 148 ./ theta) + 10 .^ -(8.3 .+ 154 ./ theta) + 10 .^ -(8.7 .+ 168 ./ theta);
    rateCoef[12] = 10 .^ -(8.8 .+ 167 ./ theta) + 10 .^ -(8.5 .+ 174 ./ theta) + 10 .^ -(8.7 .+ 175 ./ theta);
    rateCoef[13] = 10 .^ -(8.2 .+ 211 ./ theta) + 10 .^ -(10.1 .+ 254 ./ theta) + 10 .^ -(9.2 .+ 262 ./ theta);

end