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

S=[20,50,100]

begin
	# list to store networks
	global networks = []
	# monitoring variable 
	global l = length(networks)
    for i in enumerate(S)
        # while loop
        while l < reps
            # generate network
            A_bool = EcologicalNetworks.nichemodel(S[i],0.15) 
            # convert the UnipartiteNetwork object into a matrix of 1s and 0s
            Ad = adjacency(A_bool)
            A = Int.(Ad)
            # ensure size of network is correct
            println(size(A,1))
            # calculate connectance
            co = sum(A)/(size(A,1)^2)
            # ensure that connectance = 0.15
            if co == 0.15
                push!(networks, A)
                # save network is co = 0.15
            end
            global l = length(networks)
        end
    end
end
# simulate