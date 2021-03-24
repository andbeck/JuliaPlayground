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

Z = [1.0,10.0,100.0]
h = [1.0,1.2, 2.0]
K = exp10.(range(-1,1,length = 10))
T = [293.15-10, 293.15-5, 293.15, 293.15+5, 293.15+10]

## The Collector ##

outdf = DataFrame(Z = [], h = [], K = [], T = [],
                    persistence = [], biomass = [], stability = [])

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
                p = model_parameters(A, Z = Z[i], h = h[j], K = [K[k]], T = T[l], rewire_method = :none, functional_response = :classical)
                ScaleRates!(p)
                out = simulate(p, bm, start = 0, stop = 2000000000, interval_tkeep = 10000)
                persistence = species_persistence(out, last = 1000)
                biomass = population_biomass(out, last = 1000)
                stability = population_stability(out, last = 1000)
                diversity = foodweb_evenness(out, last = 1000)

                push!(outdf, [Z[i], h[j], K[k], T[l], persistence, biomass, stability])
            end
        end
    end
end


outdf

CSV.write("outdf.csv", outdf)


(3*3*10*5)

p = model_parameters(A, Z = Z[1], h = h[1], K = [K[1]], T = T[1], functional_response = :classical)
ScaleRates!(p)
out = simulate(p, bm, start = 0, stop = 2000000000, interval_tkeep = 10000) 
plot(out[:B], legend = false)