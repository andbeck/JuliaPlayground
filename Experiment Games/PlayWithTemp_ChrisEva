using BioEnergeticFoodWebs, Plots, Random, Pkg
using DataFrames, CSV

## Utilities for Scaling and Temperature ##
include("utils.jl")

## The Starting Matrix ## 
# need to generate reps....
Random.seed!(21)
A = BioEnergeticFoodWebs.nichemodel(20,0.15)
bm = rand(size(A,1))

## The Grid of Values ##

Z = [10.0]
h = [2.0]
K = exp10.(range(-1,1,length = 10))[8]
T = [293.15, 293.15+3]

## The Collector ##

outdf_no_rewire = DataFrame(Z = [], h = [], K = [], T = [], persistence = [], biomass = [], stability = [], diversity = [])
outdf_rewire = DataFrame(Z = [], h = [], K = [], T = [], persistence = [], biomass = [], stability = [], diversity = [])

#= 
This is the master loop
Use functional_response=:classical and ScaleRates! and don't forget we are in real time
So that stop is big and interval_tkeep to keep it fast.
=#

for i in 1:length(Z)
    for j in 1:length(h)
        for k in 1:length(K)
            for l in 1:length(T)
                println(i, j, k, l)
                p_nr = model_parameters(A, Z = Z[i], h = h[j], K = [K[k]], T = T[l], rewire_method = :none, functional_response = :classical)
                p_r = model_parameters(A, Z = Z[i], h = h[j], K = [K[k]], T = T[l], rewire_method = :ADBM, functional_response = :classical)
                ScaleRates!(p_nr)
                ScaleRates!(p_r)
                # 10 years
                out_nr = simulate(p_nr, bm, start = 0, stop = 60*60*24*365*10, interval_tkeep = 10000)
                out_r = simulate(p_r, bm, start = 0, stop = 60*60*24*365*10, interval_tkeep = 10000)
                
                persistence_nr = species_persistence(out_nr, last = 1000)
                biomass_nr = population_biomass(out_nr, last = 1000)
                stability_nr = population_stability(out_nr, last = 1000)
                diversity_nr = foodweb_evenness(out_nr, last = 1000)

                persistence_r = species_persistence(out_r, last = 1000)
                biomass_r = population_biomass(out_r, last = 1000)
                stability_r = population_stability(out_r, last = 1000)
                diversity_r = foodweb_evenness(out_r, last = 1000)

                push!(outdf_no_rewire, [Z[i], h[j], K[k], T[l], persistence_nr, biomass_nr, stability_nr, diversity_nr])
                push!(outdf_rewire, [Z[i], h[j], K[k], T[l], persistence_r, biomass_r, stability_r, diversity_r])
            end
        end
    end
end


outdf_no_rewire
outdf_rewire

CSV.write("outdf_nr.csv", outdf_no_rewire)
CSV.write("outdf_r.csv", outdf_rewire)

plot(outdf_no_rewire[1,6])
plot!(outdf_no_rewire[2,6])

(3*3*10*5)

p = model_parameters(A, Z = Z[1], h = h[1], K = [K[1]], T = T[1], functional_response = :classical)
ScaleRates!(p)
p[:r] # should be a set of smalls and 0.

# simulate for x years 60*60*24*365*x in seconds, keep relatively few
out = simulate(p, bm, start = 0, stop = 60*60*24*365*10, interval_tkeep = 10000) 
plot(out[:B], legend = false)