import Pkg

using BioEnergeticFoodWebs, DataFrames

import Random.seed!

Pkg.status()

# Utils functions
# Totally not working
inlcude("utils_new.jl") # from Chris' temperature work
include("TempScaleRates_APB.jl") # andrew isolates the temp scaling functions only

# make a network
A = [0 1 0 ; 0 0 1 ; 0 0 0]

#### TL;dr
# Setting T in the model_parameters and using our utils are the same.
# Seems to default to 20C and produce scaling in line with utils
# ANY use of built in functions produces different output
# possible that

# model parameters with no T specified
p1 = model_parameters(A)

# model with T specific - should trigger a change
# appears to use 293.15 = 20C

p2 = model_parameters(A, T = 293.15)

# model with T specified and formal calls to ExponentialBA() Temp Dependence
# According to the Documentation, T0 defaults to 298.15K = 25C
p3 = model_parameters(A, T = 293.15,
    growthrate = ExponentialBA(:growth),
    metabolicrate = ExponentialBA(:metabolism),
    handlingtime = ExponentialBA(:handlingtime),
    attackrate = ExponentialBA(:attackrate)
)

# model with T specified and formal calls to ExtendedEpply() Temp Dependence
# According to the Documentation, T0 defaults to 298.15K = 25C

# DOES NOT WORK
# p4 = ScaleRates!(p2)

# Manually create P4 and P4 with Utils Functions
# p4 is at 20 reference
# p5 is at 25 reference (as in ExponentialBA in the documentation function)

p4 = p2
p5 = p2

# set referece T
T = p4[:T] #293.15 = 20C
M = p4[:bodymass]

# use bits from utils
# skipping metabolism and HT
ri = ScaleGrowth(M, T)
ar = ScaleAttack(M, T)
p4[:r] = ri .* p4[:is_producer]
p4[:ar] = ar


ri_25 = ScaleGrowth2(M, T)
ar_25 = ScaleAttack2(M, T)
p5[:r] = ri_25 .* p5[:is_producer]
p5[:ar] = ar_25


[p1[:r], p1[:ar]]
[p2[:r], p2[:ar]]
[p3[:r], p3[:ar]]
[p4[:r], p4[:ar]]
[p5[:r], p5[:ar]]
