# Rhi Test
# https://poisotlab.github.io/BioEnergeticFoodWebs.jl/latest/man/temperature/#Temperature-dependence-for-body-sizes-1

using BioEnergeticFoodWebs 
using Plots
using Random

A = [0 1 1 ; 0 0 1 ; 0 0 0] #omnivory motif
T_use = 290.0
Z_use = 10.0

p_basic = model_parameters(A, T = T_use, Z = Z_use)
p_aqua = model_parameters(A, T = T_use, Z = Z_use, TSR_type = :mean_aquatic)
p_terr = model_parameters(A, T = T_use, Z = Z_use, TSR_type = :mean_terrestrial)
p_max = model_parameters(A, T = T_use, Z = Z_use, TSR_type = :maximum)

Random.seed!(123)
b0 = rand(3)
#b0 = [26.3, 15.2, 4.3]
#b0 = [1.8, 0.7, 0.2]

out_basic = simulate(p_basic, b0)
out_TSR_aq = simulate(p_aqua, b0)
out_TSR_terr = simulate(p_terr, b0)
out_TSR_max = simulate(p_max, b0)

p1 = plot(out_basic[:t], out_basic[:B])
p2 = plot(out_TSR_aq[:t], out_TSR_aq[:B])
p3 = plot(out_TSR_terr[:t], out_TSR_terr[:B])
p4 = plot(out_TSR_max[:t], out_TSR_max[:B])

# All the Same!? 
plot(p1, p2, p3, p4,  layout = (2,2))
