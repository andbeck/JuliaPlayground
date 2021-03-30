using BioEnergeticFoodWebs, Distributions, Plots, EcologicalNetworks
using Statistics, CSV, DataFrames, Random, Cairo, Fontconfig
using Gadfly

## Parameters
# FR parameter from Williams and Martinez 2004
q = [0.0, 0.1, 0.25, 1.0]

# paramters for BEFW here.
# note .+ for adding scalar to vector
h = q .+ 1

#α=[0.92, 1.0, 1.08]
k = exp10.(range(-1,1,length = 10))

## replicates (networks)
reps = 5

## Data Frame set up to collect df = DataFrame(α = [], K = [], network = [], diversity = [], stability = [], biomass = [])
df2 = DataFrame(h = [], K = [], network = [], diversity = [], stability = [], biomass = [])

# Make N Networks
Random.seed!(21)
global networks = []
global l = length(networks)

## creates reps networks
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

for n in 1:reps
    A = networks[n]
    ## choose h or α
    for i in 1:length(h)
    # for i in 1:length(α)
        for j in 1:length(k)

            println(n,i,j)

            ## choose h or α
            p = model_parameters(A, h = h[i], K = [k[j]])
            # p = model_parameters(A, α = α[i], K = [k[j]])
            bm = rand(size(A,1))
            out = simulate(p, bm, start = 0, stop = 2000)

            diversity = foodweb_evenness(out, last = 1000)
            stability = population_stability(out, last = 1000)
            biomass = total_biomass(out, last = 1000)

            ## use for alpha
            #push!(df, [α[i], k[j], n, diversity, stability, biomass])
            ## use for h
            push!(df2, [h[i], k[j], n, diversity, stability, biomass])

            println("moving to next iteration")
        end
    end
end

df # alpha
df2 # h


## are there any NaN - yes.
bb=df2[!,"biomass"]
dd=df[!, "diversity"]
ss=df[!,"stability"]

filter(isnan,bb)
filter(isnan,dd2)
filter(isnan,ss2)

# get them gone
df_clean = df[.!isnan.(df[!, :biomass]),:]
df_clean = df[.!isnan.(df[!, :diversity]),:]
df_clean = df[.!isnan.(df[!, :stability]),:]

df2_clean = df2[.!isnan.(df2[!, :biomass]),:]
df2_clean = df2[.!isnan.(df2[!, :diversity]),:]
df2_clean = df2[.!isnan.(df2[!, :stability]),:]

# check
df_clean #alpha
df2_clean # h

# generate summaries
df_biomass = combine(groupby(df2, [:α, :K]), :biomass .=> mean)
df_stability = combine(groupby(df2, [:α, :K]), :stability .=> mean)
df_diversity = combine(groupby(df2, [:α, :K]), :diversity .=> mean)

# summaries h
df2_biomass = combine(groupby(df2, [:h, :K]), :biomass .=> mean)
df2_stability = combine(groupby(df2, [:h, :K]), :stability .=> mean)
df2_diversity = combine(groupby(df2, [:h, :K]), :diversity .=> mean)


# make plots of summaries using Gadfly
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


p1=Gadfly.plot(df2_biomass, x = :K, y = :biomass_mean, colour = :h, Geom.line,
Geom.point, Scale.color_discrete_manual("black","red","purple","green"))

p2=Gadfly.plot(df2_stability, x = :K, y = :stability_mean, colour = :h, Geom.line,
Geom.point, Scale.color_discrete_manual("black","red","purple","green"))

p3=Gadfly.plot(df2_diversity, x = :K, y = :diversity_mean, colour = :h, Geom.line,
Geom.point, Scale.color_discrete_manual("black","red","purple","green"))

plotStack = hstack(p1, p2, p3)
draw(PDF("test.pdf", 10inch, 4inch), plotStack)




# does the BEFW work for h = 2?

Random.seed!(21)
A = Int.(adjacency(EcologicalNetworks.nichemodel(20,0.15)))
p = model_parameters(A, h= 1.0, K = [1.0])
bm = rand(size(A,1))
out = simulate(p, bm, start = 0, stop = 2000)

Random.seed!(21)
A = Int.(adjacency(EcologicalNetworks.nichemodel(20,0.15)))
p = model_parameters(A, h= 1.25, K = [1.0])
bm = rand(size(A,1))
out = simulate(p, bm, start = 0, stop = 2000)

Random.seed!(21)
A = Int.(adjacency(EcologicalNetworks.nichemodel(20,0.15)))
p = model_parameters(A, h= 2.0, K = [1.0])
bm = rand(size(A,1))
out = simulate(p, bm, start = 0, stop = 2000)



pl = Plots.plot([NaN], [NaN],
                label = "",
                ylims = (0,1.1),
                leg = :bottomright,
                foreground_colour_legend = nothing,
                xticks = (log10.(k), string.(round.(k, digits = 1))),
                xlabel = "Carrying Capacity",
                ylabel = "")
#Set marker shapes.
shp = [:square, :diamond, :utriangle]
#Set line types.
ls = [:solid, :dash, :dot]
#Set colours.
clr = [RGB(174/255, 139/255, 194/255), RGB(188/255, 188/255, 188/255), RGB(124/255, 189/255, 122/255)]
#When we define colours in JUlia they are printed.
#Set legend labels.
lbl = ["1.0", "1.1", "1.25", "2.0"]

pl
#Make the plot.
for (i, h) in enumerate(h)
    #Subset.
    tmp = df2[df2.h .== h, :]
    #Remove NaN values.
    tmp = tmp[.!(isnan.(tmp.diversity)), :]
    #Calculate mean across reps.
    meandf = combine(groupby(tmp, :K), :diversity => mean)
    #Command to avoid printing legends multiple times.
    l = i == 1 ? lbl[i] : ""
    #Add to pl.
    plot!(pl, log10.(meandf.K), meandf.diversity_mean,
    mc = :white,
    msw = 3,
    lw = 2,
    label = lbl[i],
    seriestype = [:line :scatter])
end
