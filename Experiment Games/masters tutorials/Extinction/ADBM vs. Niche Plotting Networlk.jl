using(BioEnergeticFoodWebs)
using(EcologicalNetworks)
using(EcologicalNetworksPlots)

import Random.seed!

# 1 is strange
# 87, 12, 870

seedN = 870

seed!(seedN)
A = EcologicalNetworks.nichemodel(20,0.1)
I = initial(FoodwebInitialLayout, A)
plot(I, A)
scatter!(I,A)


seed!(seedN)
N = BioEnergeticFoodWebs.nichemodel(20, 0.1)
A2 = EcologicalNetworks.UnipartiteNetwork(EcologicalNetworks.Bool.(N)) #transform the interaction matrix into an EcologicalNetworks object
I2 = initial(FoodwebInitialLayout, A2)
plot(I2, A2)
scatter!(I2,A2)
sum(A2)

seed!(seedN)
diffN = true

while diffN
    A3 = ADBM_foodweb(20)[1]
    print(sum(A3));
    diffN = sum(A3) != 50
end

sum(A3)
A3 = EcologicalNetworks.UnipartiteNetwork(EcologicalNetworks.Bool.(A3))
I3 = initial(FoodwebInitialLayout, A3)
plot(I3, A3)
scatter!(I3,A3)

plts = [plot(I,A), plot(I2,A2), plot(I3,A3)]
plot(plts..., layout = (1,3))