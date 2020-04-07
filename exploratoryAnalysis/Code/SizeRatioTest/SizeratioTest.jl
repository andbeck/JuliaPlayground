using BioEnergeticFoodWebs,Plots
Plots.gr()

srand(1)

function simulate_rewire()
    #generate inital web to get bodysizes
    A = nichemodel(20,0.3)
    biomass = rand(20)

    #results
    results = Vector{Float64}(10)

    for i = 1:10
        println("sim: ",i)
        #each loop we generate the initiial ADBM web with a set b
        p = model_parameters(A,Z=10.0,rewire_method = :ADBM, Nmethod = :biomass,b = linspace(0.01,10.0,10)[i]) #gives params for ADBM
        Adbm = BioEnergeticFoodWebs.ADBM(20,p,biomass)

        #We then get the simulation parameters using the ADBM web
        p = model_parameters(Adbm,bodymass = p[:bodymass], b = linspace(0.01,10.0,10)[i],rewire_method = :ADBM, Nmethod = :biomass)

        #We burn in for 1000 steps
        s = simulate(p,biomass,stop = 500)
        n = sum(s[:B][end,:] .> 0) - 1 #take the number of SP before

        #continue the simulation with a forced extinction
        b = s[:B][end,:]
        b[(find(rand(b[b .> 0.0])))] = 0.0

        #simulate for 1000 steps
        s2 = simulate(p,b,stop = 1000)

        results[i] = sum(s2[:B][end,:] .> 0)/n
    end

    return(results)
end

simulate_rewire()

function f()
    results = Array{Float64,2}(100,10)
    for i = 1:100
        println("rep: " , i)
        results[i,:] = simulate_rewire()

    end
    return(results)
end


a = f()

writecsv("../../Results/SizeRatioTest/test.csv",a)
histogram(a[:,1])
histogram(a[:,2])
