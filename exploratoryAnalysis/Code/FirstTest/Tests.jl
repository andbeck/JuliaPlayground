using BioEnergeticFoodWebs,Plots
Plots.gr()

srand(1)
# Using this interaction matrix
A = [0 0 0 0 0;
     1 0 0 0 0;
     1 0 0 0 0;
     0 1 1 0 0;
     0 0 0 1 0]

#Extinctions are detected at first
p = model_parameters(A, rewire_method = :Gilljam,Z=10.0)
b = [1.0,1.0,10.0,10.0,10.0]
s = simulate(p,b,stop = 1000)
plot(s[:B])


function simulate_rewire()
    A = nichemodel(20,0.3)
    p = model_parameters(A,Z=10.0)
    b = rand(20)

    sim = simulate(p,b,stop = 500)
    # plot(sim[:B])

    n = sum(sim[:B][end,:] .> 0.0)
    bio = sim[:B][end,:]
    bio[(find(rand(bio[bio .> 0.0])))] = 0.0

    p2 = model_parameters(A,Z=10.0)
    sim2 = simulate(p2,bio)
    n2 = sum(sim2[:B][end,:] .> 0.0)
    # plot(vcat(sim[:B],sim2[:B]))

    p3 = model_parameters(A,Z=10.0,rewire_method = :Gilljam)
    sim3 = simulate(p3,bio)
    n3 = sum(sim3[:B][end,:] .> 0.0)
    # plot(vcat(sim[:B],sim3[:B]))

    p4 = model_parameters(A,Z=10.0,rewire_method = :ADBM,Nmethod = :biomass)
    sim4 = simulate(p4,bio)
    n4 = sum(sim4[:B][end,:] .> 0.0)
    #  plot(vcat(sim[:B],sim4[:B]))

    return([n2/n,n3/n,n4/n])
end

simulate_rewire()

function f()
    results = Array{Float64,2}(100,3)
    for i = 1:100
        results[i,:] = simulate_rewire()
        println(i)
    end
    return(results)
end


a = f()

writecsv("../../Results/Test/testGilljam.csv",a)
histogram(a[:,1])
histogram(a[:,2])
