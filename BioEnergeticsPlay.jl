## Working with BioEnergeticFoodWebs

using Plots, BioEnergeticFoodWebs, Random

# gr
gr()

# basic simulation with BEFW from MEE paper
# generate random food web
Random.seed!(123)
A = nichemodel(20,0.15, toltype = :abs)
a = range(0.92, stop = 1.08, length = 3)
K = exp10.(range(-1,1,length = 9))

diversity = zeros((9,3))

for j in 1:length(a)  # change from MEE with J 1.0
    for i in 1:length(K)  # changed from MEE paper with J 1.0 + typo in paper

        global(diversity)

        p = model_parameters(A, α = a[j], K = K[i], productivity = :competitive)

        # We start each simulation with
        # random biomasses in ]0;1[
        bm = rand(size(A, 1))

        # Andfinally,wesimulate.
        out = simulate(p, bm, start = 0, stop = 2000) # omit method from MEE

        # And measure the output
        diversity[i,j] = foodweb_evenness(out, last = 1000) # omit eps() from MEE
    end
end

diversity
plot(K, diversity[1:9,1])
plot!(K, diversity[1:9,2])
plot!(K, diversity[1:9,3])

# plot(K, diversity)
