#= The ICQ =#

import Random
using Plots, Distributions, Statistics

Random.seed!(123)

# define the mean and sd of the distribution
mean_sd = Normal(10,3) # mean 10, sd 3
# create the distribution - this could be the body sizes
bodymass = rand(mean_sd, 1000) # sample random normal 100 times

# find the Interquartile Range - it's the 25th and 75th quantile
IQR = quantile(bodymass, [0.25, 0.75])
IQR

# plot the bodysizes
histogram(bodymass)
# add the IQR to the picture to see
vline!(IQR, linewidth = 4)
xlabel!("Body Masses")
