## Working with BioEnergeticFoodWebs

using Plots, BioEnergeticFoodWebs, Random, DataFrames, CSV

# gr
gr()

# basic simulation with BEFW from MEE paper
# generate random food web

Random.seed!(42) # replaces srand()

# 9 x 9 x rep design
a = range(0.5, 1, length=9)
K = exp10.(range(-1, 1, length=9))
reps = 10

# set up DF to write into 
df = DataFrame(a=Float64[], K=Float64[], replicate=Float64[],
    Diversity=Float64[], Stability=Float64[], Biomass=Float64[])

# rep loop over niche models (probably should have made them all first)
for k in 1:reps

    global (diversity)
    global (stability)
    global (total_biomass)

    A = nichemodel(20, 0.15, tolerance=0.01)
    
    # loop over a and K values
    for j in 1:length(a)  # change from MEE with J 1.0
        for i in 1:length(K)  # changed from MEE paper with J 1.0 + typo in paper
            
            # run the model    
            p = model_parameters(A, α=a[j], K=K[i], productivity=:competitive)

            # We start each simulation with
            # random biomasses in ]0;1[
            bm = rand(size(A, 1))

            # And finally,we simulate.
            out = simulate(p, bm, start=0, stop=2000) # omit method from MEE

            # And measure the output
            diversity = foodweb_evenness(out, last=1000) # omit eps() from MEE
            stability = population_stability(out, last=1000)
            tot_biomass = total_biomass(out, last=1000)

            # write each line to the df
            push!(df, [a[j] K[i] k diversity  stability  tot_biomass]) 
        end
    end
end


# looking at it
df
last(df, 6)
head(df)
df[[11,15,20],:]

# write it to file
CSV.write("myDF.csv", df)