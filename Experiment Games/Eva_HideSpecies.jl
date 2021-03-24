using BioEnergeticFoodWebs, Plots

#[9b49b652] BioEnergeticFoodWebs v1.2.0 #dev-1.3.0 (https://github.com/PoisotLab/BioEnergeticFoodWebs.jl.git)

include("utils.jl")

# Simple example

#======================== PART 1 : Issues not solved =================================#
A = [0 1 0 ; 0 0 1 ; 0 0 0] #linear food chain

p = model_parameters(A, T = 293.15, h = 2.0, Z = 10.0)
ScaleRates!(p)
b0 = B0(A, p[:K]) #biomass starting values are K for plants and 1/8 mean K for consumers
b0[2] = 0.0 #force extinction of primary consumers

#simulate for just a few time steps -> we don't want top consumer to go extinct yet
s = simulate(p, b0, stop = Int(60*60*24), interval_tkeep = Int(60*60)) 
plot(s[:B], labels = ["Top pred" "Prim. cons." "Prod."], 
    c = :black, linestyle = [:dash :dot :solid], leg = :right) #nothing has really happened yet, rates are small => change is slow...

id_extinct = p[:extinctions] #sp 2 is extinct, as expected

    #Now we can update the parameter object and see if the top predator then becomes a plant
nonextinct = trues(3)
nonextinct[id_extinct] .= false
newA = A[nonextinct, nonextinct]
M = p[:bodymass][nonextinct]
newp = model_parameters(newA, bodymass = M, h = 2.0, T = 293.15)
newp[:is_producer] #non-green thing is now green thing...

#======================== PART 2 : Test solution ======================================#
# We can try conserving species identity by keeping the original matrix 
# and just turning off extinct species... 
A = [0 1 0 ; 0 0 1 ; 0 0 0] #linear food chain
p = model_parameters(A, T = 293.15, h = 2.0, Z = 10.0)
ScaleRates!(p)
b0 = B0(A, p[:K]) 
b0[2] = 0.0 #force extinction of primary consumers

s = simulate(p, b0, stop = Int(60*60*24), interval_tkeep = Int(60*60)) 
plt1 = plot(s[:B], labels = ["Top pred" "Prim. cons." "Prod."], 
    c = :black, linestyle = [:dash :dot :solid], leg = :right) 
xlabel!("time (mn)")
ylabel!("biomass")
b1 = s[:B][end,:] #final biomass

p[:extinctions]

#we can check that extinct species have a final biomass of 0.0 (and are not zombie species that still feed the dynamics)
@assert all(b1[p[:extinctions]] .== [0.0])
s1 = simulate(p, b1, stop = Int(60*60*24*364.25*50), interval_tkeep = Int(60*60*24*364.25))
plt2 = plot(s1[:B], c = :black, linestyle = [:dash :dot :solid], leg = false) 
xlabel!("time (year)")
plot(plt1, plt2, layout = grid(1,2))

p[:extinctions]