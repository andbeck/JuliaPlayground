#=
Extinction cascade tutorial
This script shows how to 
- initialize the BEFWm with a food web and your choice of parameters 
- run a burn-in to remove structurally doomed species (to avoid detecting them as false secondary extinctions) and achieve equilibrium dynamics
- target a species for primary extinction
- run the model 
- detect secondary extinctions / measure things
=#

# Packages
import Pkg
using BioEnergeticFoodWebs
using Distributions
using Plots
using DataFrames

# Make sure you are using the right BEFW!
Pkg.status()

# Set seed
import Random.seed!
seed!(5436) # setting a seed makes it possible to replicate exactly the same results, even when there is randomness involved

#=
Step 1: Import or generate a food web
You can use an empirical food web or use a model (niche model, ADBM...) to generate realistic food webs
We start by creating a food web using the niche model
- we initialize with 100 species and a connectance of 0.2 
=#
S = 100 # number of species
con = 0.2 # connectance
A = nichemodel(S,con) # use niche model to generate network 

#=
APB EDITS HERE
Step 2: Set the model parameters 
Here we are running a very basic simulation, with:
- h (hill exponent) = 2 (Type III functional response - stabilising)
- K (carrying capacity for basal species / producers) = 10 
- Z = (consumer-resource size ratio) = 10.0 (consumers are 10x bigger than resources)
=#
p = model_parameters(A, h = 2.0, K = [10.0], Z = 10.0)

#=
Step 3: Run a burn-in
Run a burn-in phase to remove some structurally doomed species and achieve equilibrium dynamics
We've also provided the code to save some information about each species (TL, in and out degree) and the community (total biomass, persistence) in two seperate dataframes
=#
b0 = rand(S) # set some initial biomasses at random
s = simulate(p, b0, stop = 2000) # simulate 
s = simulate(p, b0, stop = 2000) # simulate #APB: should we call this out_burn, rather than s?  I fear s and S will get confused.
plot(s[:B], legend = false) # plot

# Below, we've created a dataframe that stores species level information for all species in the network, including their TL, indegree (number species that eat them) and outdegree (number of species that they eat), as well as information about when some of those species went extinct (which timestep) and during which simulation (0 = burn-in and 1 = second simulation post perturbation)
# To provide some biological context, a specialist species will have a small outdegree value (the consumer has few prey species), whereas a generalist species will have a large outdegree value (they consume many prey species). Also, a basal species might be expected to have a large indegree value (a lot of species consume them), whereas a large predator near the top of the food web will have a small indegree value (few species consume them) 

# we first create a dataframe of S rows (100) by preallocating the column names and size of the dataframe, and filling it with NaNs (same as NAs in R)
species_data = DataFrame(fill(NaN,S,6), [:ID, :TL, :indegree, :outdegree, :extinction_time, :extinction_sim])

# APB EDITS HERE
# Here we fill in some of the details from the end of the burnin

# species names
species_data.ID = collect(1:S) # <- should be S, not 100?
# species Trophic Level (their rank)
species_data.TL .= p[:trophic_rank]
# speciues indegree
species_data.indegree .= vec(sum(A, dims = 1))
# species outdegree
species_data.outdegree .= vec(sum(A, dims = 2))

# APB EDITS HERE
# species extinsions_sim
species_data[p[:extinctions],:extinction_sim] .= Int.(zeros(length((p[:extinctions]))))

# This one line is really valuable: 
 # view the dataframe for species 
 # that went extinct during the burn-in (38 species went extinct)
species_data[p[:extinctions], :]

# NOTE THAT species extinsions_time is filled in AFTER everthing is done (see very end)

# APB EDITS HERE - Should we call the network_data to match the species_data above rather than df?
# here we construct a dataframe for community level metrics by preallocating the type and name of each column
df = DataFrame([Int64, Float64, Float64, Float64, Float64],[:sim, :total_biomass, :persistence, :richness, :R_x]) 


# populate df using push!
push!(df, [0, total_biomass(s, last = 500), species_persistence(s, last = 500), species_richness(s, last = 500), species_richness(s, last = 500)/S]) 
# again, we've used 0 for this first simulation (sim column) which is the burn in phase
# we could also have calculated stability or diversity, or some other out of the box metrics like network height (max TL) or size structure (the realised Z value of the community)

#=
Step 4: Save the biomass of each species at equilibrium 
Basically, extract biomass of each species at the end of the first simulation
=#
b1 = s[:B][end,:] # b1 will serve as a starting point for our next simulation

#=
Step 5: Using a while loop, lets loop through primary extinction events until 50% (R50) of all species have gone extinct
Here, we will remove species in decending order from the largest to the smallest based on TL
The catch here is we can only remove species that are persistent in the community and have to take this into account when looping
After each simulation/loop our databases (species_data and df) will be updated and our i counter will be updated based on the species richness of the community
=#

# we need to estimate species richness at the end of the burnin as our reference point for R50.
# we calcuate that with the species_richness function applied to the burnin simulations

global i = species_richness(s, last = 500) #we need the global macro to be able to use that variable in the while loop
global j = 0.0

while i >= S/2
    println("i = $i and j = $j") # keep track of loop
    
    #=
    Step 5a: Primary extinction
    =#
    
    # There are many options here: random extinctions, Large-Small, Top to Bottom.....

    # set a rule for targetting a species, here the consumer with the highest trophic rank
    # this is quite complicated but we've scattered it onto multiple lines to provide clarity (or sanity checks)
    # Basically, we had to make sure that the first species with maximum TL wasn't already extinct and therefore that the removal of species was having the correct effect in our simulations
    
    is_extinct_0 = falses(S) # make a vector of 0's (falses) the size of S
    is_extinct_0[p[:extinctions]] .= true # set extinction species to 1 (true)
    tl = p[:trophic_rank] # calculate trophic rank of each species
    id_sorted = sortperm(tl, rev = true) # create a vector of species identity sorted by decreasing trophic level
    is_extinct_sorted = is_extinct_0[id_sorted] # reorder the extinct vector by species id
    id_primext = id_sorted[.!is_extinct_sorted][1] # select first species that fit our criterion that isn't already extinct

    #=
    Step 5b: Simulate forward (e.g. from where the burnin finished, but now with making an extinction!)
    =#
    b1[id_primext] = 0.0 #set biomass of primary extinct species to 0
    s1 = simulate(p, b1, stop = 2000) #run the simulation APB EDIT/SUGGEST: s1 or out_SOMETHING?

    #=
    Step 5c: Update community dataframe
    =#
    # update community dataframe - easy using push!
    # this calcuates stuff caused by the extinction above
    push!(df, [j, total_biomass(s1, last = 500), species_persistence(s1, last = 500), species_richness(s1, last = 500),
        species_richness(s1, last = 500)/S])

    #=
    Step 5d: Update species dataframe
    =#
    # identify any newly extinct species
    global is_extinct_1 = falses(S)
    global is_extinct_1[p[:extinctions]] .= true
    new_extinct = findall(is_extinct_1 .!= is_extinct_0) 
    # for newly extinct species, specify that they went extinct during the jth run by setting :extinction_sim in species_data as j
    # add relevant information to the species_data data frame
    species_data[new_extinct,:extinction_sim] .= fill(j,length(new_extinct)) 

    #=
    Step 5e: Update biomasses again for the next extinction event
    =#
    b1 = s1[:B][end,:]

    #=
    Step 5f: Update counters
    =#
    global i = minimum(df.richness)
    global j = j + 1.0

end

# Check out df data DataFrame
df
# look at richess declining with events
plot(df[!,"richness"], df[!,"sim"])

# Following our loop, for all extinct species, we can specify the time step at which they went extinct
# this can be done at the end of the loop because p[:extinctionstime] is not overwritten and keeps all records of extinction times 

ext_time = p[:extinctionstime] #this vector is ordered by time of extinction and not species id so we need to reorder it
ext_time_sp = [i[2] for i in ext_time] #extract species identity from this object
sort_extinct = sortperm(ext_time_sp) #create a vector to order by identity
ext_time_ts = [i[1] for i in ext_time[sort_extinct]] #use it to reorder extinction times
species_data[is_extinct_1,:extinction_time] .= ext_time_ts #pass extinction times to the data frame

# look at the species data!
species_data
# look at the distribution of extinction times
histogram(species_data[!,"extinction_time"])

# THE END #








