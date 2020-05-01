using BioEnergeticFoodWebs,Plots
Plots.gr()

srand(101)
# Using this interaction matrix
A = nichemodel(20,0.3)
#Extinctions are detected at first
p = model_parameters(A, rewire_method = :ADBM,Nmethod = :biomass,Z=10.0)
b = rand(20)

A = BioEnergeticFoodWebs.ADBM(20,p,b)
p = model_parameters(A, rewire_method = :ADBM,Nmethod = :biomass,Z=10.0)

s = simulate(p,b,stop = 1000)

plot(s[:B])

s[:p][:Arrays]



for i = 1:length(s[:p][:Arrays])
    file = string("array",i,".csv")
    writecsv(joinpath("./../Results/ArrayImage",file),s[:p][:Arrays][i])
end

writecsv(joinpath("./../Results/ArrayImage","Biomass.csv"),s[:B])
