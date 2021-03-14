using EcologicalNetworks, Plots, EcologicalNetworksPlots

Fweb = simplify(nz_stream_foodweb()[5])

# first by trophic level with random postion within
I = initial(FoodwebInitialLayout, Fweb)
plot(I, Fweb)
scatter!(I, Fweb)

# generate ForceDirectedLayout
for step in 1:4000
  position!(ForceDirectedLayout(true, false, 2.5), I, Fweb)
end

# with ForceDirectedLayout
plot(I, Fweb)
scatter!(I, Fweb)