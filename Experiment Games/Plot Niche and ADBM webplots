# Just making plots of niche and adbms.
using DataFrames, BioEnergeticFoodWebs, ProgressMeter
import Random
Random.seed!(22)

include("src/common_utils.jl") #get all the functions!

global niche = []
global l = length(niche)

## creates reps networks
while l < 10
    A_bool = EcologicalNetworks.nichemodel(20,0.15)
    A = adjacency(A_bool)
    A = Int.(A)
    co = sum(A)/size(A,1)^2
    if co == 0.15
        push!(niche, A)
    end
    global l = length(niche)
end


global adbm = []
global l = length(adbm)

## creates reps networks
while l < 10
    A= ADBM_foodweb(S, parameters = :random)[1]
    co = sum(A)/size(A,1)^2
    if co == 0.15
        push!(adbm, A)
    end
    global l = length(adbm)
end


pts = [webplot(adbm[1]);webplot(adbm[2]);webplot(adbm[3]);webplot(adbm[4]);
webplot(niche[1]);webplot(niche[2]);webplot(niche[3]); webplot(niche[4])]
plot(pts..., layout = (2,4))


