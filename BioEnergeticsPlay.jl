## Working with BioEnergeticFoodWebs

using Plots
using BioEnergeticFoodWebs

# gr
gr()

# basic simulation with BEFW from MEE paper
# generate random food web

A = nichemodel(20,0.15, toltype = :abs)

for a in range(0.92, stop = 1.08, length = 3) # change from MEE with J 1.0
    for K in exp10.(range(-1,1,length = 9)) # changed from MEE paper with J 1.0 + typo in paper
        p = model_parameters(A, α = a, K = K, productivity = :competitive)

        # We start each simulation with
        # random biomasses in ]0;1[
        bm = rand(size(A, 1))

        # Andfinally,wesimulate.
        out = simulate(p, bm, start = 0, stop = 2000) # omit method from MEE

        # And measure the output
        diversity = foodweb_evenness(out, last = 1000) # omit eps() from MEE
    end
end
