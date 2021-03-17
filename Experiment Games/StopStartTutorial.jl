#=
Modifying biomass values "during" simulation 
This Gist show how to 
- initialize the BEFWm with a food web and your choice of parameters
- run a burn-in to remove structurally doomed species and get to the equilibrium
- community (to avoid detecting them as false secondary extinctions)
- then target a species for primary extinction or harvesting or other perturbation
- re-start the model 
- detect secondary extinctions / measure things that happen because of the extinction 
- or harvesting
=#

# set up the toolbox
using BioEnergeticFoodWebs
using Distributions
import Random.seed!

# setting a seed makes it possible to replicate exactly the same results, 
# even when there is randomness involved
seed!(5436) 

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

p = model_parameters(A, h = 2.0, K = [10.0])

#=
Step 3: Run a burn-in
We run a burn-in phase to remove some structurally doomed species 
# note too that we've set bodymass as random....
=#

#initial biomass distribution.
b0 = rand(S) 

# intitial 'burn-in' simulation
s = simulate(p, b0, stop = 2000)

#=
Step 4: 
1. identify extinct species
2. Discard these extinct species and reorganise the matrix of interactions
and the distribution of bodymasses
    - we see that even without primary extinction some species move toward extinction
    we want to remove these species because otherwise we'll detect them as secondary extinctio later 
    (and they're not, they were just doomed from the start)
3. Record equilibrium dynamics from the initial 2000 steps
=#

extinct = p[:extinctions] #id of extinct species 
non_extinct = trues(S) # create vector of TRUE for all species.

# Now set extinct species to false - don't forget the broadcast operator "." 
# because we are modifying multiple entries

non_extinct[extinct] .= false 

# Now - subset to keep only non extinct species
# there were 70 extinctions.... down to a 30 x 30 array
A = A[non_extinct, non_extinct]

# Now subset body mass information to keep only non extinct species
bodymass_vector = p[:bodymass][non_extinct]  

# now collect the equilibrium dynamics of these non-extinct:
# These serve as a starting point for next simulation
# average biomass of remaining species during the last 500 timesteps 
# (equilibrium dynamics)
b1 = population_biomass(s, last = 500)[non_extinct] 

#= 
 --- REVIEW ----
At this stage, all you've done is what we call the burn-in.  
We've initialised a network, run it for 2000 timesteps
Collected the network at the end of it
Reduced the bodymass distribution to the same species
Estimated the biomass of each of these species.

These last three are the STARTING point for the experiments
=#

#=
Step 5: Primary extinction or perturbation experiment
=#

# set a rule for targetting a species for extinction OR harvesting
# here we target the consumer with the highest trophic rank
tl = trophic_rank(A) #calculate trophic rank
id_primext = findfirst(tl .== maximum(tl)) #first species that fit our criterion

#=
Step 5.a -> IF EXTINCTION
=#

b1[id_primext] = 0.0 #set biomass of primary extinct species to 0
p = model_parameters(A, h = 2.0, K = [10.0]) #recalculate parameters with the new food web

# re-define body masses so that it's only the non extinct species from the burnin.
# If you don't do this then the model_parameters() function will reset body masses 
# based on Z and the body mass distribution of the network will change
p[:bodymass] = bodymass_vector

#run the simulation for another 2000 time steps!
s = simulate(p, b1, stop = 2000) 

# this is the list of species that went extinct after removing 
# the target 1˚ extinction
sec_ext = p[:extinctions]

# what would you do now to estimate biomass or stability after the perturbation?

#=
Step 5.b -> IF PERTURBATION
=#

# set a rule for the % loss of a species
# here for example, the biomass of the target species will be reduced to 
# 80% of its original equilibrium value
b_loss = 0.8 

# apply loss of biomass to the target species
b1[id_primext] = b1[id_primext] * b_loss 

#recalculate parameters with the new food web 
p = model_parameters(A, h = 2.0, K = [10.0]) 

# re-define body masses so that it's only the non extinct species from the burn in.
# If you don't do this then the model_parameters() function will reset body masses 
# based on Z and the body mass distribution of the network will change
p[:bodymass] = bodymass_vector  

#run the simulation for another 2000 time steps
s = simulate(p, b1, stop = 2000) 

# These are the species that went extinct as a result of the harvesting
sec_ext = p[:extinctions] 

# what would you do now to estimate biomass or stability after the perturbation?

#=
Step 4,5...: Do an order of extinctions again and again?
Step 4,5...: Apply harvesting to more than 1 species?   
If you want to continue to remove or perturb species, you can repeat steps 4 and 5
(e.g. make loops!)
=#
