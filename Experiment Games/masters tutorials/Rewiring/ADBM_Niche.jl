using BioEnergeticFoodWebs
using EcologicalNetworks
include("Rewiring/utils.jl") 
import Random.seed!
seed!(22)

# Generate a food web from ADBM and Niche

S = 30
A, adbm_p = ADBM_foodweb(S)

C_niche = sum(A)/S^2
A2 = BioEnergeticFoodWebs.nichemodel(S, C_niche)

Z = 10 #consumers are in average 10 times bigger than their prey
M = Z .^ (trophic_rank(A) .- 1)
#now for visualisation purposes let's order the matrix by increasing body mass
M_order = sortperm(M)
M = M[M_order]
A2 = A2[M_order, M_order]
#and we can visualize the interaction matrix 


sum(A2)/S^2 #interesting.

#You can check the food web with webplot
pltA = webplot(A, false)
pltA2 = webplot(A2, false)
plts=[pltA, pltA2]
plot(plts..., layout = grid(2,1))

# not sure this is right
p_ADBM = model_parameters_modif(A, adbm_p, :none) # default is h=2.0
p_niche = model_parameters(A2, h = 2.0, Z = 10.0, T = 293.15)
ScaleRates!(p_niche, 10.0)

ScaleRates!(p_niche,)

tstop = Int.(60*60*24*364.25*50) #simulation time = 50 years here
Δt = Int.(60*60*24*30*12) #sample every year
ϵ = 1e-10 #extinction threshold

b0 = rand(S)

# this isn't working - different time scales.
s_niche = simulate(p_niche, b0, stop = tstop, extinction_threshold = ϵ)
s_adbm = simulate(p_adbm, b0, stop = tstop, extinction_threshold = ϵ)

plt_niche = plot(s_niche[:B])
plt_adbm = plot(s_adbm[:B])

plts_trace = [plt_niche, plt_adbm]

plot(plts_trace..., layout = grid(2,1))




# Generate a food web

S = 30
C = 0.25
A = EcologicalNetworks.nichemodel(S, C).edges 
A = Int.(Array(A))
#we can calculate species mass based on their trophic rank and Z, the consumer resource mass ratio
Z = 10 #consumers are in average 10 times bigger than their prey
M = Z .^ (trophic_rank(A) .- 1)
#now for visualisation purposes let's order the matrix by increasing body mass
M_order = sortperm(M)
M = M[M_order]
A = A[M_order, M_order]
#and we can visualize the interaction matrix 
pltA = webplot(A, true) #true means that we want to have consumers plotted as rows(i) and resources as columns (j)
#the resulting plot shows the interaction matrix A[i,j] = 1 is represented by a black dot and 
#means that i eats j. Dots on the diagonal indicate canibalism (i eats i).
# niche model food webs will look messy because their generation is not based on mass

# Generate the sets of parameters

#no rewiring is the default, technically you don't need to specify rewire_method in that case
p_none = model_parameters(A, bodymass = M, h = 2.0, rewire_method = :none) 
#DO stands for diet overlap, that's from Staniczencko's paper
p_do = model_parameters(A, bodymass = M, h = 2.0, rewire_method = :DO)
#DS stands for diet similarity, from Gilljam's paper 
p_ds = model_parameters(A, bodymass = M, h = 2.0, rewire_method = :DS)