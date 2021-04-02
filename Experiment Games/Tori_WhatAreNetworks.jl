#= This demonstrates what nichemodel networks are
- They are constructed by the nichemodel (Williams and Martinez, 2000 Nature)
- You provide S = Species Richness and C = Connectance.
- For each of these numbers it makes a network of species interactions
- Each time you do this, the network is different (it is a probabilistic model)
=#

import(Random)
using BioEnergeticFoodWebs

Random.seed!(21)

# set up the collecting zone for the networks
global networks = []
# define how many you want
reps = 3 # how many networks
# set up a counter for the loop to make the networks
global l = length(networks)

# define the Species Richness and Connectance to make them
S = 10 # species richness
C = 0.15 # connectance

# create the neworks.
while l < reps
    # make a network with the nichemodel
    A = nichemodel(S,0.15)
    
    # this checks to make sure the connectance is 0.15
    co = sum(A)/size(A,1)^2
    
    # if it is definitely 0.15, add it to the list
    if co == 0.15
        push!(networks, A)
    end
    global l = length(networks)
end

# look at all of the networks
networks # all three networks

networks[1] # the first one
networks[2]

# are the first and second the same?  No.
# and you can see that above.
networks[1]==networks[2]