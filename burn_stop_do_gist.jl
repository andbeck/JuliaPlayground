#=
Modifying biomass values "during" simulation
This Gist show how to
- initialize the BEFWm with a food web and your choice of parameters
- run a burn-in to remove structurally doomed species (to avoid detecting them as false secondary extinctions)
- target a species for primary extinction or perturbation
- run the model
- detect secondary extinctions / measure things
=#

using BioEnergeticFoodWebs, Distributions, Plots
import Random.seed!
#[9b49b652] BioEnergeticFoodWebs v1.2.0
#[31c24e10] Distributions v0.24.15

seed!(5436) # setting a seed makes it possible to replicate exactly the same results, even when there is randomness involved

#=
Step 1: Import or generate a food web
You can use an empirical food web or use a model (niche model, ADBM...) to generate realistic food webs
We start by creating a food web using the niche model
- we initialize with 100 species and a connectance of 0.2
=#
S = 100
A = nichemodel(S,0.2)

#=
Step 2: Set the model parameters
Here we are running very basix simulations, with:
- h (hill exponent) = 2
- K (carrying capacity for basal species / producers) = 10
=#
p = model_parameters(A, h = 2.0, K = 10.0)

#=
Step 3: Run a burn-in
We run a burn-in phase to remove some structurally doomed species
=#
b0 = rand(S) #initial biomass
s_full = simulate(p, b0, stop = 2000)
p1=plot(s_full[:B])
#=
Step 4: Discard extinct species / Record equilibrium dynamics
we see that even without primary extinction some species move toward extinction
we want to remove these species because otherwise we'll detect them as secondary extinctio later
(and they're not, they were just doomed from the start)
=#
extinct = p[:extinctions] #id of extinct species
non_extinct = trues(S) #initialize a vector for non extinct species
non_extinct[extinct] .= false #set extinct species to false - don't forget the broadcast operator "." because we are modifying multiple entries
A = A[non_extinct, non_extinct] #we keep only non extinct species
b1 = population_biomass(s_full, last = 500)[non_extinct] #average biomass of remaining species during the last 500 timesteps (equilibrium dynamics)
# the equilibrium dynamics will serve as a starting point for later simulations

#=
Step 5: Primary extinction or perturbation
=#
# set a rule for targetting a species, here the consumer with the highest trophic rank
tl = trophic_rank(A) #calculate trophic rank
id_primext = findfirst(tl .== maximum(tl)) #first species that fit our criterion

#=
IF EXTINCTION
=#
b1[id_primext] = 0.0 #set biomass of primary extinct species to 0
p = model_parameters(A, h = 2.0, K = 10.0) #recalculate parameters wit the new food web
s2 = simulate(p, b1, stop = 2000) #run the simulation
sec_ext = p[:extinctions]

p2=plot(s2[:B])
plot(p1, p2, layout = (2,1))

#=
What percentage are extinct?
=#

1-length(extinct)/S


#=
IF PERTURBATION
=#
# set a rule
b_loss = 0.8 # here for example, the biomass of the target species will be reduced to 80% of its original equilibrium value
b1[id_primext] = b1[id_primext] * b_loss #set biomass of primary extinct species to 0
p = model_parameters(A, h = 2.0, K = [10.0]) #recalculate parameters wit the new food web
s = simulate(p, b1, stop = 2000) #run the simulation
sec_ext = p[:extinctions]

#=
Step 4,5...: Do it again
If you want to continue to remove or perturb species, you can repeat steps 4 and 5
#=
