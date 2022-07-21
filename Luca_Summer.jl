#setting up the packages
using BioEnergeticFoodWebs
using EcologicalNetworks
using JLD2
using Statistics
using Plots
using CSV
using DataFrames
using Random
Random.seed!(21)

#generate initial network
A_bool = EcologicalNetworks.nichemodel(20,0.15)
#convert into matrix
Ad = adjacency(A_bool)
A = Int.(Ad)
#calculate connectance
co = sum(A)/(size(A,1)^2)
#parameters, use default values
p = model_parameters(A)
# assign biomasses
bm = rand(size(A,1))
#simulate outcome
out = simulate(p, bm, start=0, stop=2000)
#experiment
#define experiment by creating vectors of variables, fix repetition
metabolicrate = [0.1, 1, 5]

#number of reps
reps = 5
#DataFrame
df= DataFrame(metabolicrate = [], network = [], diversity = [], stability = [], biomass = [])
#while loop
begin
    global networks = []
    global l = length(networks)
    while l < reps
        for i in 1:length(metabolicrate)
        A_bool = EcologicalNetworks.nichemodel(20,0.15)
        Ad = adjacency(A_bool)
        A = Int.(Ad)
        co = sum(A)/(size(A,1)^2)
        if co == 0.15
            push!(networks, A)
        end
        global l = length(networks)
    end
end
end
#for loop over networks
for h in 1:reps
    A = networks[h]
            for i in 1:length(metabolicrate)
                p = model_parameters(A, metabolicrate = metabolicrate[i], productivity=:competitive)
                bm = rand(size(A,1))
                out = simulate(p, bm, start=0, stop=2000)
                metabolicrate_num = metabolicrate[i]
                diversity = foodweb_evenness(out, last = 1000)
                stability = population_stability(out, last = 1000)
                biomass = total_biomass(out, last = 1000)
                println(("metabolicrate = $metabolicrate_num", "network = $h"))
            end
        end

#plots
pl = Plots.plot([NaN], [NaN],
                label = "",
                ylims = (0,1.1),
                leg = :bottomright,
                foreground_colour_legend = nothing,
                xticks = (log10.(k), string.(round.(k, digits = 1))),
                xlabel = "Metabolic Rate",
                ylabel = "Food web diversity (evenness)")

shp = [:square, :diamond, :utriangle]
ls = [:solid, :dash, :dot]
clr = [RGB(174/255, 139/255, 194/255), RGB(188/255, 188/255, 188/255), RGB(124/255, 189/255, 122/255)]
lbl = ["coexistence", "Neutral", "Exclusion"]

for (i, x) in enumerate(x)
    tmp = df[df.x .== x, :]
    tmp = tmp[.!(isnan.(tmp.diversity)), :]
    gdf = groupby(tmp, :K)
    meandf = combine(gdf, :diversity => mean)
    l = i == 1 ? lbl[i] : ""
    plot!(pl, log10.(meandf.K), meandf.diversity_mean,
    msc = clr[i],
    mc = :white,
    msw = 3,
    markershape = shp[i],
    linestyle = ls[i],
    lc = clr[i],
    lw = 2,
    label = lbl[i],
    seriestype = [:line :scatter])
end