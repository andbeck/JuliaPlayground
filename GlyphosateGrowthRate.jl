# packages
using BioEnergeticFoodWebs , EcologicalNetworks , JLD2 , Statistics , Plots , CSV , DataFrames , Random

# random seed
Random.seed!(21)

# vector of growth rate ,r, (-50%, -20%, -10%, 0, +5%, +20%)
r = [0.50, 0.80, 0.90, 1.0, 1.05, 1.20]

# vector of carrying capacity
k = [0.50, 0.80, 0.90, 1.0]

# set repititions and DataFrame
reps = 5
df = DataFrame(r = [], K = [], network = [], diversity = [], stability = [], biomass = [], biomass2 = [], time = [])

# make while loop
begin
	# list to store networks
	global networks = []
	# monitoring variable
	global l = length(networks)
	# while loop
	while l < reps
	    # generate network
	    A_bool = EcologicalNetworks.nichemodel(20,0.15)
	    # convert the UN object into a matrix of 1s and 0s
		Ad = adjacency(A_bool)
	    A = Int.(Ad)
	    # calculate connectance
	    co = sum(A)/(size(A,1)^2)
	    # ensure that connectance = 0.15
	    if co == 0.15
	        push!(networks, A)
	        # save network if co = 0.15
	    end
	    global l = length(networks)
	end
end

# loop over networks

for h in 1:reps

    A = networks[h]

    # loop over r
    for i in 1:length(r)

        # loop over K
        for j in 1:length(k)

        # create model parameters
        p = model_parameters(A, r = [r[i]], K = [k[j]], productivity=:competitive)

        # assign biomass
        bm = rand(size(A,1))

        # simulate
        out = simulate(p, bm, start=0, stop=2000)

        # calculate output metrics
        diversity = foodweb_evenness(out, last = 1000)
        stability = population_stability(out, last = 1000)
        biomass = total_biomass(out, last = 1000)
        biomass2 = out[:B]
        time = out[:t]

        # dummy naming variables
        r_num = r[i]
        K_num = k[j]

        # push to df
        push!(df, [r[i], k[j], h, diversity, stability, biomass, biomass2, time])

        # print some stuff
        println(("r=$r_num", "network = $h"))
        end
    end
end



# see output
df
first(df,6)
last(df,6)

# make plot
Plots.plot(out[:t], out[:B], legend = true, ylabel = "Biomass", xlabel = "Time")

Plots.plot(r, biomass, legends = true, ylabel = "Biomass", xlabel = "growth rate")

# make plot 2: Biomass vs Growth Rate
# initialise empty plot
pl = Plots.plot([NaN], [NaN],
label = "",
ylims = (0,1.1),
xlims = (0, 1.5),
leg = :bottomright,
foreground_colour_legend = nothing,
xlabel = "Growth Rate",
ylabel = "Biomass",
title = "How biomass changes at different growth rates")

# marker and line type
shp = [:circle]
ls = [:solid]

for (i, r) in enumerate(r)
    # subset
    tmp = df[df.r .== r, :]
    # remove NaN values
    tmp = tmp[.!(isnan.(tmp.biomass)), :]
    # calculate mean across reps
    gdf = groupby(tmp, :r)
	meandf = combine(gdf, :biomass => mean)
    # command to avoid printing legends multiple times
    l = i == 1 ? lbl[i] : ""
    # add to pl
    plot!(pl, log10.(meandf.r), meandf.biomass_mean,
              mc = :white,
              msw = 3,
              markershape = shp[i],
              linestyle = ls[i],
              lw = 2,
              seriestype = [:line :scatter])
end






# CODE BELOW HERE NOT USING AT THE MOMENT

# initialise empty plot
pl = Plots.plot([NaN], [NaN],
label = "",
ylims = (0,1.1),
xlims = (0, 2000),
leg = :bottomright,
foreground_colour_legend = nothing,
xlabel = "Time",
ylabel = "Biomass",
title = "How biomass changes at different growth rates")

# set colours (red, orange, yellow, blue, light blue, green)
clr = [RGB(200/255, 0/255, 0/255), RGB(230/255, 140/255, 22/255), RGB(210/255, 220/255, 65/255), RGB(0/255, 0/255, 255/255), RGB(22/255, 223/255, 230/255), RGB(0/255, 205/255, 0/255)]

# set legend labels
lbl = ["-50%","-20%", "-10%", "0%", "+5%", "+20%"]

for (i, r) in enumerate(r)
    # subset
    tmp = df[df.α .== α, :]
    # remove NaN values
    tmp = tmp[.!(isnan.(tmp.diversity)), :]
    # calculate mean across reps
    gdf = groupby(tmp, :K)
	meandf = combine(gdf, :diversity => mean)
    # command to avoid printing legends multiple times
    l = i == 1 ? lbl[i] : ""
    # add to pl
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
