using BioEnergeticFoodWebs, EcologicalNetworks, Random, Pkg
using Plots

gr()

# set up initial network and burn in sims
Random.seed!(21)
S = 20
A_init = BioEnergeticFoodWebs.nichemodel(S,0.15)
p_init = model_parameters(A_init, h= 2.0, K = [10.0])
bm_init = rand(size(A_init,1))
out_burn = simulate(p_init, bm_init, start = 0, stop = 2000, 
    interval_tkeep = 1, use = :stiff)

# Intersting to see that there are micro-fluctuations.
# affected by use = :stiff though
p0 = plot(log10.(out_burn[:B][1700:2000,:]), legend = false)
p1 = plot(out_burn[:B][1700:2000,[1]], legend = false)
p2 = plot(out_burn[:B][1700:2000,[5]], legend = false)
p3 = plot(out_burn[:B][1700:2000,[10]], legend = false)
p4 = plot(out_burn[:B][1700:2000,[20]], legend = false)
l = @layout [a ; b c; d e]
plot(p0, p1, p2,p3,p4, layout=l)

# identify extinctions
extinct_1 = p_init[:extinctions]
length(extinct_1)
non_extinct = trues(S) 
non_extinct[extinct_1] .= false 
non_extinct

# subset the matrix
A_next = A_init[non_extinct, non_extinct]
S-size(A_next)[1] == length(extinct_1)

# subset the bodymasses
bodymass_next = p_init[:bodymass][non_extinct]
# subset and collect the biomass (last value approach vs. mean?)
# whether the mean or the last value is taken doesn't seem to matter
biomass_next = population_biomass(out_burn, last = 1000)[non_extinct] 

# reset the model parameters with the new A-matrix
p_next = model_parameters(A_next, h = 2.0, K = [10.0])
# reset the bodymass vector with subsetted bodymass
p_next[:bodymass] = bodymass_next

# resimulate
out_forward = simulate(p_next, biomass_next, stop = 2000, 
    interval_tkeep = 1, use = :stiff)

# resimulate again to test for whether the extinctions cause it or not (they do)
biomass_next_2 = population_biomass(out_forward, last = 1000)
out_forward_2 = simulate(p_next, biomass_next_2, stop = 2000, 
    interval_tkeep = 1, use = :stiff)


l = @layout [a ; b c]

# burn in vs. first restart
plot_init = plot(log10.(out_burn[:B][1900:2000,:]), 
    ylims=(-1.5,1), legend = false)
plot_forward = plot(log10.(out_forward[:B][1:100,:]), 
    ylims=(-1.5,1), legend = false)
plot(plot_init, plot_forward, layout = (1,2))


# first restart vs. second restart
plot_forward = plot(log10.(out_forward[:B][1900:2000,:]), 
    ylims=(-1.5,1), legend = false)
plot_forward_2 = plot(log10.(out_forward_2[:B][1:100,:]), 
    ylims=(-1.5,1), legend = false)
plot(plot_forward, plot_forward_2, layout = (1,2))




# names for plotting and for named array if we could
sp_original = 1:1.0:20
sp_after = sp_original[non_extinct]

