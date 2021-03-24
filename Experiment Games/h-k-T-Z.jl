using BioEnergeticFoodWebs, Plots, Random, Pkg
using DataFrames, CSV

Random.seed!(21)
A = BioEnergeticFoodWebs.nichemodel(20,0.15)
bm = rand(size(A,1))

Z = [1.0,10.0,100.0]
h = [1.0,1.2, 2.0]
K = exp10.(range(-1,1,length = 10))
T = [293.15-10, 293.15-5, 293.15, 293.15-5, 293.15+10]

(3*3*10*5)

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


outdf

plot(outdf[!,"T"], outdf[!,"persistence"])

CSV.write("outdf.csv", outdf)


