#import the functions from the Pkg manager
import Pkg

#activate the project
Pkg.activate(".")

#install the packages 
Pkg.instantiate()

#check the pkg status
Pkg.status()

#Packages I am  using 
using BioEnergeticFoodWebs, EcologicalNetworks, CSV, Random, Plots, DataFrames, Statistics, Pkg, VegaLite

## Utilities for Scaling and Temperature ##
include("common_utils.jl")

#Set the random seed 
Random.seed!(21)

#Set the temperature range
T_range = [0:1:40;] .+ 273.15# temperature ranging from 0 to 40C

Z_range = 10.0 .^ [1:1:5;] # consumer-resource mass ratio ranging from 10 to 100000

A1 = [0 1 0 ; 0 0 1 ; 0 0 0] #linear food chain

A = BioEnergeticFoodWebs.nichemodel(20,0.15)

df = DataFrame(T = [], Z = [], biomass = [], persistence = [], diversity = [], stability = []) # initialize outputs data frame
df_deactivation = DataFrame(T = [], Z = [], biomass = [], persistence = [], diversity = [], stability = [])

year = Int(60*60*24*364.25)


for z in Z_range 
	for t in T_range
		println("T = $t and log(z) = $(log10(z))")
		p_tmp = model_parameters(A, T = t, Z = z, h = 2.0, functional_response = :classical)
		p_deactivated = model_parameters(A, T=t, Z=z, h = 2.0, functional_response = :classical, growthrate = ExtendedBA(:growthrate, parameters_tuple = (T_opt = 300.15, )),# change the parameters values for the allometric exponent using a named tuple
								metabolicrate = ExtendedBA(:metabolicrate, parameters_tuple = (deactivation_energy_vertebrate = 1.02, T_opt_invertebrate = 293.15)),
								attackrate = ExtendedBA(:attackrate, parameters_tuple = (deactivation_energy_vertebrate = 1.02, T_opt_invertebrate = 293.15))) 
        ScaleRates!(p_tmp, 10.0)
		ScaleRates!(p_deactivated, 10.0)
		out = simulate(p_tmp, rand(20), stop = year*3000, interval_tkeep = year)
		out2 = simulate(p_deactivated, rand(20), stop = year*3000, interval_tkeep = year)
		tmp = (T = t, Z = z, biomass = total_biomass(out, last = 500), persistence = species_persistence(out, last = 500), diversity = foodweb_evenness(out, last = 1000), stability = population_stability(out, last = 1000))
		deactivated = (T = t, Z = z, biomass = total_biomass(out2, last = 500), persistence = species_persistence(out2, last = 500), diversity = foodweb_evenness(out2, last = 500), stability = population_stability(out2, last = 500))
		push!(df, tmp)
		push!(df_deactivation, deactivated)
	end
end

TZ_array_biomass = fill(NaN, length(Z_range), length(T_range)) #preallocate a 2d array for biomass value
TZ_array_biomass_deac = fill(NaN, length(Z_range), length(T_range))

TZ_array_persistence = fill(NaN, length(Z_range), length(T_range)) #preallocate a 2d array for persistence values
TZ_array_persistence_deac = fill(NaN, length(Z_range), length(T_range))

TZ_array_diversity = fill(NaN, length(Z_range), length(T_range)) #preallocate a 2d array for diversity value
TZ_array_diversity_deac = fill(NaN, length(Z_range), length(T_range))

TZ_array_stability = fill(NaN, length(Z_range), length(T_range)) #preallocate a 2d array for stability value
TZ_array_stability_deac = fill(NaN, length(Z_range), length(T_range))

for (i,z) in enumerate(Z_range) 
	for (j,t) in enumerate(T_range)
		
		tmp = df[(df.Z .== z) .& (df.T .== t),:]
		deac = df_deactivation[(df_deactivation.Z .== z) .& (df_deactivation.T .== t),:]
		
		TZ_array_biomass[i,j] = tmp.biomass[1]
		TZ_array_persistence[i,j] = tmp.persistence[1]
		TZ_array_diversity[i,j] = tmp.diversity[1]
		TZ_array_stability[i,j] = tmp.stability[1]
		
		TZ_array_biomass_deac[i,j] = deac.biomass[1]
		TZ_array_persistence_deac[i,j] = deac.persistence[1]
		TZ_array_diversity_deac[i,j] = deac.diversity[1]
		TZ_array_stability_deac[i,j] = deac.stability[1]

	end
end

#HEATMAPS
#----------------

p1 = heatmap(string.(T_range .- 273.15), string.(log10.(Z_range)), log10.(TZ_array_biomass), title  = "log10(Total biomass) - Heatmap", xlabel = "Temperature (C)", ylabel = "log10(Z)") 
savefig(p1 , "/Users/emilyhatchwell/Desktop/Biomass Heatmap(1b)")

p1_deac= heatmap(string.(T_range .- 273.15), string.(log10.(Z_range)), log10.(TZ_array_biomass_deac), title  = "log10(Total biomass) - Heatmap with deactivation energy", xlabel = "Temperature (C)", ylabel = "log10(Z)") 
savefig(p1_deac , "/Users/emilyhatchwell/Desktop/Biomass Heatmap Deactivation(1b)")

#---------------

p2 = heatmap(string.(T_range .- 273.15), string.(log10.(Z_range)), TZ_array_persistence, title  = "Persistence - Heatmap", xlabel = "Temperature (C)", ylabel = "log10(Z)")
savefig(p2 , "/Users/emilyhatchwell/Desktop/Persistence Heatmap(1b)")

p2_deac = heatmap(string.(T_range .- 273.15), string.(log10.(Z_range)), TZ_array_persistence_deac, title  = "Persistence - Heatmap with deactivation energy", xlabel = "Temperature (C)", ylabel = "log10(Z)")
savefig(p2_deac , "/Users/emilyhatchwell/Desktop/Persistence Heatmap Deactivation(1b)")

#----------------

p3 = heatmap(string.(T_range .- 273.15), string.(log10.(Z_range)), TZ_array_diversity, title  = "Diversity - Heatmap", xlabel = "Temperature (C)", ylabel = "log10(Z)") 
savefig(p3 , "/Users/emilyhatchwell/Desktop/Diversity Heatmap(1b)")

p3_deac= heatmap(string.(T_range .- 273.15), string.(log10.(Z_range)), TZ_array_diversity_deac, title  = "Diversity - Heatmap with deactivation energy", xlabel = "Temperature (C)", ylabel = "log10(Z)") 
savefig(p3_deac , "/Users/emilyhatchwell/Desktop/Diversity Heatmap Deactivation(1b)")

#----------------

p4 = heatmap(string.(T_range .- 273.15), string.(log10.(Z_range)), TZ_array_stability, title  = "Stability - Heatmap", xlabel = "Temperature (C)", ylabel = "log10(Z)") 
savefig(p4 , "/Users/emilyhatchwell/Desktop/Stability Heatmap(1b)")

p4_deac= heatmap(string.(T_range .- 273.15), string.(log10.(Z_range)), TZ_array_stability_deac, title  = "Stability - Heatmap with deactivation energy", xlabel = "Temperature (C)", ylabel = "log10(Z)") 
savefig(p4_deac , "/Users/emilyhatchwell/Desktop/Stability Heatmap Deactivation(1b)")

#-----------------