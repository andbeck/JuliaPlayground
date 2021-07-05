#=
How to use the rewirinh methods in BioEnergeticFoodWebs
author: Eva Delmas
date: April 25th, 2021
=#

# Set up

#we will nee a different version of BEFWM with the newest features 
#so start my removing the package with rm and reinstalling with 
#add BioEnergeticFoodWebs#dev-2.0.0
using BioEnergeticFoodWebs
using EcologicalNetworks
include("utils.jl") 
import Random.seed!
seed!(22)

# Generate a food web
#=
Because we are going to use the ADBM as one of the rewiring models, 
we are going to also use the ADBM as a generative model (instead of
the niche model). This is because the ADBM rewire on a more global
scale than the other two, as a consequence, if we start with a niche 
model, the first rewiring event will reorganize the whole food web, 
and the ADBM will work with a different food web than the Diet 
Overlap and DIet Generality models. 
=#

S = 30
#we don't need to specify connectance anymore - connectance is an 
#output of the network generation process with ADBM
#C = 0.25 
A, adbm_p = ADBM_foodweb(S)
#You can check the food web with webplot
pltA = webplot(A, true) #true means that we want to have consumers plotted as rows(i) and resources as columns (j)
#the resulting plot shows the interaction matrix A[i,j] = 1 is represented by a black dot and 
#means that i eats j. Dots on the diagonal indicate canibalism (i eats i).

# Generate the sets of parameters for the BEFWM

#no rewiring is the default, technically you don't need to specify rewire_method in that case
p_none = model_parameters_modif(A, adbm_p, :none)
#DO stands for diet overlap, that's from Staniczencko's paper
p_do = model_parameters_modif(A, adbm_p, :DO)
#DS stands for diet similarity, from Gilljam's paper 
p_ds = model_parameters_modif(A, adbm_p, :DS)
#ADBM stands for Allometric Diet Breadth Model, from Petchey's paper
p_adbm = model_parameters_modif(A, adbm_p, :ADBM)

# Set up simulations

# /!\ You'll need toc hange the simulation time because we have changed the units of the biological rates
tstop = Int.(60*60*24*364.25*500) #simulation time = 500 years
Δt = Int.(60*60*24*30) #sample every month
ϵ = 1e-10 #extinction threshold
#when a species biomass reaches the extinction threshold, it's considered as
#extinct. It's important to know the extinction threshold when using rewiring because 
#rewiring will be triggered at each extinction event during the simulations. 
# (We'll see later how this can be different for the ADBM)
b0 = rand(S)

# Simulate biomass dynamics for each

#no rewiring
s_none = simulate(p_none, b0, stop = tstop, interval_tkeep = Δt, extinction_threshold = ϵ)
#diet overlap
s_do = simulate(p_do, b0, stop = tstop, interval_tkeep = Δt, extinction_threshold = ϵ)
#diet similarity 
s_ds = simulate(p_ds, b0, stop = tstop, interval_tkeep = Δt, extinction_threshold = ϵ)
#ADBM
s_adbm = simulate(p_adbm, b0, stop = tstop, interval_tkeep = Δt, extinction_threshold = ϵ)

# Visualize the results
plt_dyn_none = plot(s_none[:B], leg = false, c = :black, ylims = (0,8), xlabel = "time", ylabel = "biomass")
plt_dyn_do = plot(s_do[:B], leg = false, c = :black, ylims = (0,8),xlabel = "time", ylabel = "biomass")
plt_dyn_ds = plot(s_ds[:B], leg = false, c = :black, ylims = (0,8),xlabel = "time", ylabel = "biomass")
plt_dyn_adbm = plot(s_adbm[:B], leg = false, c = :black, ylims = (0,8),xlabel = "time", ylabel = "biomass")
plt_mat_none = webplot(p_none[:A], true)
plt_mat_do = webplot(p_do[:A], true)
plt_mat_ds = webplot(p_ds[:A], true)
plt_mat_adbm = webplot(p_adbm[:A], true)

plt = [plt_dyn_none, plt_dyn_do, plt_dyn_ds, plt_dyn_adbm
     , plt_mat_none, plt_mat_do, plt_mat_ds, plt_mat_adbm]

plot(plt..., layout = grid(2,4), size = (1000,400))
