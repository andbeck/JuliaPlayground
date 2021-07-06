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

sum(A2)/S^2 #interesting.

#You can check the food web with webplot
pltA = webplot(A, false)
pltA2 = webplot(A2, false)
plts=[pltA, pltA2]
plot(plts..., layout = grid(2,1))

# not sure this is right
p_ADBM = model_parameters_modif(A, adbm_p, :none)
p_niche = model_parameters(A)

tstop = Int.(60*60*24*364.25*50) #simulation time = 50 years here
Δt = Int.(60*60*24*30) #sample every month
ϵ = 1e-10 #extinction threshold

b0 = rand(S)

# this isn't working - different time scales.
s_niche = simulate(p_niche, b0, stop = tstop, interval_tkeep = Δt, extinction_threshold = ϵ)
s_adbm = simulate(p_adbm, b0, stop = tstop, interval_tkeep = Δt, extinction_threshold = ϵ)

plt_niche = plot(s_niche[:B])
plt_adbm = plot(s_adbm[:B])

plts_trace = [plt_niche, plt_adbm]

plot(plts_trace..., layout = grid(2,1))