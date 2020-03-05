# Andrew Learns More Julia
# cmd-enter binding is for running code by line/block
# cmd-shift-enter runs entire script
# ] to add packages
# ctrl-c to get out of packages

# --------------------------------

# using is used to load packages after installation with ]
using Plots
using Distributions
using Statistics
using DataFrames
using GLM, Random

# specify gr module in Plot
gr()

# create some data

x = [1,2,3]
y = [10,20,55]

# plot the data as a line and add the points
plot(x,y)
scatter!(x, y)

# sequences of numbers are from, by, to
xx = collect(0:0.1:1)
xx

# multiply vectors
xx*2

# plot the multiplied vector
plot(xx,xx*2)
scatter!(xx,xx*2)

# randn is in base J.  Only 0,1 mean sd
# need Distributions to do more
xx = randn(1000)
histogram(xx, bins= 5, alpha = 0.5)
histogram!(xx, bins= 20, alpha = 0.5, colour = "red")

# with Distributions and rand and Normal
nn = rand(Normal(10,1),1000)
nn
histogram(nn, bins = 20, colour = "goldenrod", alpha = 0.4)

# simple ANOVA?
Random.seed!(1)
data = DataFrame(y = rand(100), x = categorical(repeat([1, 2, 3, 4], 25)))
data
lm(@formula(y ~ x), data)
lm(@formula(y ~ x), data, contrasts = Dict(:x => DummyCoding()))
