## Working with BioEnergeticFoodWebs

using Plots, BioEnergeticFoodWebs, Random, DataFrames

# gr
gr()

# basic simulation with BEFW from MEE paper
# generate random food web

Random.seed!(42) # replaces srand()

# 9 x 9 x 2 design
a = range(0.5, 1, length = 9)
K = exp10.(range(-1,1,length = 9))

diversity = zeros((9,9,2))
stability = zeros((9,9,2))
tot_biomass = zeros((9,9,2))

for k in 1:2

    global(diversity)
    global(stability)
    global(total_biomass)

    A = nichemodel(20,0.15, tolerance = 0.01)

    for j in 1:length(a)  # change from MEE with J 1.0
        for i in 1:length(K)  # changed from MEE paper with J 1.0 + typo in paper

            p = model_parameters(A, α = a[j], K = K[i], productivity = :competitive)

            # We start each simulation with
            # random biomasses in ]0;1[
            bm = rand(size(A, 1))

            # Andfinally,wesimulate.
            out = simulate(p, bm, start = 0, stop = 2000) # omit method from MEE

            # And measure the output
            diversity[i,j,k] = foodweb_evenness(out, last = 1000) # omit eps() from MEE
            stability[i,j,k] = population_stability(out, last = 1000)
            tot_biomass[i,j,k] = total_biomass(out, last = 1000)
        end
    end
end

diversity
stability
tot_biomass

# facets
plot(K, diversity, layout = (3,1))

# layout grid.arrange
p1 = plot(K, diversity, ylabel = "Evenness")
p2 = plot(K, stability, ylabel = "Stability")
p3 = plot(K, tot_biomass, ylabel = "Biomass",xlabel = "K")

plot(p1, p2, p3, layout = (3, 1), legend = true)

heatmap(K, a, diversity)
heatmap(K, a, stability)
heatmap(K, a, tot_biomass)
