using Interpolations

function formatAtmosphere(filePath, dz, imax, jmax, sd)

# TODO: Define densities for primals, and just average to get the values on duals.

#----------------------------------------

# TODO: Use appropriate includes.

zLow = 60 + 0.5 * dz/1e3;             # Bottom height, 60 km + dz/2
zHigh = jmax*dz/1e3 - 0.5 * dz/1e3;   # Top height, jmax*dz - dz/2

indLow  = findfirst(==(zLow), z);
indHigh = findfirst(==(zHigh), z);
# indLow : 0.2 : indHigh
ind = indLow : round( (dz/1e3)/(z[2]-z[1]) ) : indHigh;

# Note: Make sure that you do not try to iterate through :gas.
species = [:N2 ,:O2 ,:O3 ,:H2O ,:H ,:OH ,:HO2 ,:CO2 ,:N2A ,:N2B ,:N2a ,:N2C ,:N4S ,:N2D ,:O2a ,:O2b ,:O ,:NO ,:NO2 ,:e ,:Om ,:O2m ,:O3m ,:O4m ,:OHm ,:CO3m ,:CO4m ,:NO2m ,:NO3m ,:O2mNO ,:HCO3m ,:N2p ,:Np ,:O2p ,:Op ,:O4p ,:NOp ,:Yp, :gas]
# TODO: Put this in a loop that goes over each species.
#for curr_species ∈ setdiff(species, [:gas])
 
# TODO: Could be a bug if dz does not evenly partition the range aLow:zHigh
# TODO: Change this range here to include upper limit in sucha case
interp_linear_N2 = linear_interpolation(zLow:dz:zHigh, density_N2)
output_density_N2 = VForm(map(sd[:point]) do p
  # TODO: Use an interface that picks out 3 just by giving it symbol z
  interp_linear_N2(p[2])
end)
#end

#GLMakie.mesh(sd, color=output_density_N2, colorrange=(2, 100))

# TODO: Convert cartesian -> cylindrical coordinates
# z is z
# rho is sqrt(x^2+y^2), y==0 => x
# B_ϕ is orthogonal to surface.
# First vector is defined with a corresponding vector of zs from 60.1 to 99.90 . Then we want to interpolate values between those to define the rest of the mesh.

outputTemp = Tn[ind]

outputDensity.gas = outputDensity.N2 + outputDensity.O2;

outputRateCoef = rateCoef;

return outputDensity, outputTemp, outputRateCoef

end
