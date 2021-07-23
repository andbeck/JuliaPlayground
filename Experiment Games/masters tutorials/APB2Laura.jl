using BioEnergeticFoodWebs: connectance
using BioEnergeticFoodWebs, EcologicalNetworks,DataFrames, Plots, CSV

import Random.seed!
seed!(5436)


reps = 1000 #number of repetitions (used to create loop)

df2 = DataFrame(network = [], connectance = [])

# loop to create and store 30 networks with co = 0.15
global networks1 = [] 
for h in 1:reps
    A_bool = EcologicalNetworks.nichemodel(20,0.15)
    A = Int.(A_bool.A)
    co = sum(A)/(size(A,1)^2) # calculate connectance
    push!(networks1, A)
    push!(df2, [h, co]) #h counter of networks co = connectance ## dataframe with network and connectance value
    println(h)
end

df2

# subsets using the data frame package
df2[(df2[:connectance].<0.15),:] # 12
df2[(df2[:connectance].>0.15),:] # 16
df2[(df2[:connectance].==0.15),:] # 2

# use the .| operator for OR
df2[(df2[:connectance].<0.15) .| (df2[:connectance].>0.15) ,:]

# number of observations
nrow(df2[(df2[:connectance].<0.15),:])

# histogram
histogram(df2[:connectance], bins = 10)
vline!([0.15])

df2[(df2[:connectance].<=0.1),:] #3
nrow(df2[(df2[:connectance].<=0.1),:])/reps

# to work in R..... which is more better for graphing and data manipulation
CSV.write("df2.csv", df2)
