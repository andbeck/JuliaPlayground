using DataFrames, Plots, StatsPlots
gr()

df = DataFrame(A = 1:2:1000, B = repeat(1:10, inner=50), C = 1:500)
df

@df df plot(:B,:A)


