#setup
import(Pkg)
Pkg.activate(".")

# packages
using BioEnergeticFoodWebs
using EcologicalNetworks
using Random
using Plots

# utils
include("common_utils.jl")


# set seed
Random.seed!(12345)

# make networks

# set S values and reps
global S=[20,50,100]
reps = 5

# define master collection
global networks = []

# list to store networks
for i in 1:3
    # define temporary collector for while loop
    global nets_collect = []
    # monitoring variable
    global l = length(nets_collect)
    # while loop
    while l < reps
        # generate network
        A_bool = EcologicalNetworks.nichemodel(S[i],0.15)
        # convert the UnipartiteNetwork object into a matrix of 1s and 0s
        Ad = adjacency(A_bool)
        A = Int.(Ad)
        # calculate connectance
        co = sum(A)/(size(A,1)^2)
        # ensure that connectance = 0.15
        if co == 0.15
            push!(nets_collect, A)
            # save network is co = 0.15
        end
        # update counter for while loop
        global l = length(nets_collect)
    end
    # once you have 5, while stops and we push the temporary
    # to the master collection
    push!(networks, nets_collect)
end

# master collection has three components (S = 20, 50 and 100)
length(networks)

# to see/use them, notice the double [][]
# [1][1] is S = 20, first of the 5
# [1][5] is S = 20, fifth of the 5
# [3][1] is S = 100, 1 of the five
networks[1][1]
networks[1][5]
networks[2][1]
networks[3][1]



# simulate
