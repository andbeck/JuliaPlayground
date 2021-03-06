 
Random.seed!(21)
A_en = Int.(adjacency(EcologicalNetworks.nichemodel(20,0.15)))

Random.seed!(21)
A_befw = BioEnergeticFoodWebs.nichemodel(20,0.15)

A_en == A_befw

#Pkg.status()

Random.seed!(21)
A = Int.(adjacency(EcologicalNetworks.nichemodel(20,0.15)))
p = model_parameters(A, h= 1.0, K = [1.0])
bm = rand(size(A,1))
out = simulate(p, bm, start = 0, stop = 2000)

Random.seed!(21)
A = BioEnergeticFoodWebs.nichemodel(20,0.15)
p = model_parameters(A, h= 1.0, K = [1.0])
bm = rand(size(A,1))
out2 = simulate(p, bm, start = 0, stop = 2000)