using BioEnergeticFoodWebs, EcologicalNetworks, Plots, Random, Pkg
using DataFrames

Random.seed!(21)
A = BioEnergeticFoodWebs.nichemodel(20,0.15)
bm = rand(size(A,1))

Z = [1.0,10.0]
h = [1,2]
K = [1.0, 10.0]
T = [273.15-10, 273.15, 273.15+10]


p = model_parameters(A, Z = Z[1], h = h[1], K = [K[1]], T = T[2])
out = simulate(p, bm, start = 0, stop = 2000, interval_tkeep = 1) 

outdf = DataFrame(Z = [], h = [], K = [], T = [],
                    persistence = [], biomass = [], stability = [])


for i in 1:length(Z)
    for j in 1:length(h)
        for k in 1:length(K)
            for l in 1:length(T)
                println(i, j, k, l)
                p = model_parameters(A, Z = Z[i], h = h[j], K = [K[k]], T = T[l])
                out = simulate(p, bm, start = 0, stop = 2000, interval_tkeep = 1)
                persistence = species_persistence(out, last = 1000)
                biomass = population_biomass(out, last = 1000)
                stability = population_stability(out, last = 1000)
                diversity = foodweb_evenness(out, last = 1000)

                push!(outdf, [Z[i], h[j], K[k], T[l], persistence, biomass, stability])
            end
        end
    end
end






