"""
Define the initial conditions for the density fields.

We use here example data that has been collected from the NRLMSIS model.

The available data is defined for elevations from 60 to 110 km, and so we assume that densities below 60 km are 0.0.
"""
function formatAtmosphere(mat_filename, sd)

# TODO: Define densities for primals, and just average to get the values on duals.

#----------------------------------------

# TODO: Integrate with the NRLMSIS api to get this data on-demand.
vars = matread(mat_filename)
chi_densities = vars["density"]
chi_rateCoefs = vars["rateCoef"]
# TODO: This gives values for elevations up to 110 km, but we set the mesh to up to 100 km.
chi_zs = vars["z"] # 60.0, 60.1, 60.2, ..., 110.0

# TODO: Use appropriate includes.

#zLow = 60 + 0.5 * dz/1e3;             # Bottom height, 60 km + dz/2
#zHigh = jmax*dz/1e3 - 0.5 * dz/1e3;   # Top height, jmax*dz - dz/2

#indLow  = findfirst(==(zLow), z);
#indHigh = findfirst(==(zHigh), z);
# indLow : 0.2 : indHigh
#ind = indLow : round( (dz/1e3)/(z[2]-z[1]) ) : indHigh;

# Note: Make sure that you do not try to iterate through :gas.
species = [:N2 ,:O2 ,:O3 ,:H2O ,:H ,:OH ,:HO2 ,:CO2 ,:N2A ,:N2B ,:N2a ,:N2C ,:N4S ,:N2D ,:O2a ,:O2b ,:O ,:NO ,:NO2 ,:e ,:Om ,:O2m ,:O3m ,:O4m ,:OHm ,:CO3m ,:CO4m ,:NO2m ,:NO3m ,:O2mNO ,:HCO3m ,:N2p ,:Np ,:O2p ,:Op ,:O4p ,:NOp ,:Yp, :gas]
densities = Dict{Symbol, VForm{Float64, Vector{Float64}}}()
for curr_species ∈ setdiff(species, [:gas])
  # TODO: Could be a bug if dz does not evenly partition the range aLow:zHigh
  # TODO: Change this range here to include upper limit in such a case
  #interp_linear_N2 = linear_interpolation(zLow:dz:zHigh, density_N2)
  #output_density_N2 = VForm(map(sd[:point]) do p
  #  # TODO: Use an interface that picks out 3 just by giving it symbol z
  #  interp_linear_N2(p[2])
  #end)

  #interp_linear = linear_interpolation(zLow:dz:zHigh, chi_densities[String(species)])
  interp_linear = linear_interpolation(chi_zs[:], chi_densities[String(curr_species)][:])
  output_density = VForm(map(sd[:point]) do p
    # TODO: Use an interface that picks out 3 just by giving it symbol z
    z = p[2] * 1.0e-3
    # The default chi data we have only starts at 60 km, so we set values below this altitude to 0.
    z < minimum(chi_zs[:]) ? interp_linear(minimum(chi_zs[:])) : interp_linear(z)
    #z < minimum(chi_zs[:]) ? 0.0 : interp_linear(z)
  end)
  densities[curr_species] = output_density
end

#GLMakie.mesh(s, color=densities[:N2])
#GLMakie.mesh(s, color=densities[:O2])
#GLMakie.mesh(s, color=densities[:Yp])

# TODO: Convert cartesian -> cylindrical coordinates
# z is z
# rho is sqrt(x^2+y^2), y==0 => x
# B_ϕ is orthogonal to surface.
# First vector is defined with a corresponding vector of zs from 60.1 to 99.90 . Then we want to interpolate values between those to define the rest of the mesh.

#outputTemp = Tn[ind]
interp_linear = linear_interpolation(chi_zs[:], vars["Tn"][:])
outputTemp = VForm(map(sd[:point]) do p
  # TODO: Use an interface that picks out 3 just by giving it symbol z
  z = p[2] * 1e-3
  # The default chi data we have only starts at 60 km, so we set values below this altitude to the value at 60 km.
  # TODO: Maybe just go ahead and set this to 0.0. below 60 km.
  z < minimum(chi_zs[:]) ? interp_linear(minimum(chi_zs[:])) : interp_linear(z)
end)

#outputDensity.gas = outputDensity.N2 + outputDensity.O2;
densities[:gas] = VForm(densities[:N2].data + densities[:O2].data)
#GLMakie.mesh(s, color=densities[:gas])

# These are the rate coefficients k1, k2, ...
rateCoef_names = Symbol.(keys(vars["rateCoef"]))
rateCoefs = Dict{Symbol, VForm{Float64, Vector{Float64}}}()
for rateCoef ∈ rateCoef_names
  interp_linear = linear_interpolation(chi_zs[:], chi_rateCoefs[String(rateCoef)][:])
  output_rateCoef = VForm(map(sd[:point]) do p
    # TODO: Use an interface that picks out 3 just by giving it symbol z
    z = p[2] * 1.0e-3
    # The default chi data we have only starts at 60 km, so we set values below this altitude to 0.
    # If we do not want any reactions to happen below 60 km, we should exponentially decay to 0.
    #z < minimum(chi_zs[:]) ? interp_linear(minimum(chi_zs[:])) : interp_linear(z)
    z < minimum(chi_zs[:]) ? 0.0 : interp_linear(z)
  end)
  rateCoefs[rateCoef] = output_rateCoef
end

return species, densities, outputTemp, rateCoef_names, rateCoefs

end
