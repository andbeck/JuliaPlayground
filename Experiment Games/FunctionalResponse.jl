using BioEnergeticFoodWebs, Distributions, Plots, EcologicalNetworks, Statistics, 
CSV, DataFrames, Random, Gadfly, Cairo, Fontconfig

## Parameters
α=[0.92, 1.0, 1.08]
k = exp10.(range(-1,1,length = 10))

## replicates (networks)
reps = 2

## Data Frame set up to collect
df = DataFrame(α = [], K = [], network = [], diversity = [], stability = [], biomass = [])

# Make N Networks

global networks = []
global l = length(networks)

## creatnes reps networks 
while l < reps
    A_bool = EcologicalNetworks.nichemodel(20,0.15)
    A = adjacency(A_bool)
    A = Int.(A)
    co = sum(A)/size(A,1)^2
    if co == 0.15
        push!(networks, A)
    end
    global l = length(networks)
end

# =  check network = #
networks[1]

for h in 1:reps
    A = networks[h]
    for i in 1:length(α)
        for j in 1:length(k)
            
            println(h,i,j)
            
            p = model_parameters(A, α = α[i], K = [k[j]])
            bm = rand(size(A,1))
            out = simulate(p, bm, start = 0, stop = 2000)

            diversity = foodweb_evenness(out, last = 1000)
            stability = population_stability(out, last = 1000)
            biomass = total_biomass(out, last = 1000)

            push!(df, [α[i], k[j], h, diversity, stability, biomass])

            println("moving to next iteration")
        end
    end
end

# are there any NaN - yes.
# get rid of the NaN
bb=df[!,"biomass"]
dd=df[!, "diversity"]
ss=df[!,"stability"]

filter(isnan,bb2)
filter(isnan,dd2)
filter(isnan,ss2)

# get them gone
df2 = df[.!isnan.(df[!, :biomass]),:]
df2 = df[.!isnan.(df[!, :diversity]),:]
df2 = df[.!isnan.(df[!, :stability]),:]

df2

df_biomass = combine(groupby(df2, [:α, :K]), :biomass .=> mean)
df_stability = combine(groupby(df2, [:α, :K]), :stability .=> mean)
df_diversity = combine(groupby(df2, [:α, :K]), :diversity .=> mean)

p1=Gadfly.plot(df_biomass, xgroup = :α, 
    x = :K, y = :biomass_mean,
    Geom.subplot_grid(Geom.point))

p2=Gadfly.plot(df_stability, xgroup = :α, 
    x = :K, y = :stability_mean,
    Geom.subplot_grid(Geom.point))

p3=Gadfly.plot(df_diversity, xgroup = :α, 
    x = :K, y = :diversity_mean,
    Geom.subplot_grid(Geom.point))

plotStack = hstack(p1, p2, p3)

# does not work
draw(PNG("test.png"))

@rlibrary ggplot2
R"ggplot2::ggplot($df_biomass, aes(x = K, y = biomass_mean, colour = α))+
geom_point()"


gasoline = dataset("Ecdat", "Gasoline");