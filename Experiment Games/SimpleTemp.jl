using BioEnergeticFoodWebs, Plots, Random, Pkg
using DataFrames, CSV

include("utils.jl")

## The Starting Matrix ## 
# need to generate reps....
Random.seed!(21)
A = BioEnergeticFoodWebs.nichemodel(20,0.15)
bm = rand(size(A,1))

## The Grid of Values ##
T = 293.15
# T = [293.15-10, 293.15-5, 293.15, 293.15+5, 293.15+10]

parms = model_parameters(A, Z = 10.0, T = 293.15, functional_response = :classical)
ScaleRates!(parms)

sims = simulate(parms, bm, start = 0, stop = 60*60*24*365*10, interval_tkeep = 10000)

persistence_T = species_persistence(sims, last = 1000)
biomass_T = population_biomass(sims, last = 1000)
stability_T = population_stability(sims, last = 1000)
diversity_T = foodweb_evenness(sims, last = 1000)

plot(sims[:B])