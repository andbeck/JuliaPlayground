using BioEnergeticFoodWebs EcologicalNetworks

Random.seed!(21)
A = Int.(adjacency(EcologicalNetworks.nichemodel(20,0.15)))
p = model_parameters(A, h= 1.0, K = 1.0)
bm = rand(size(A,1))
out = simulate(p, bm, start = 0, stop = 2000)