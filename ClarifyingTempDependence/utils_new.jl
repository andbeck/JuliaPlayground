#=
Useful functions
=#
using Statistics
import BioEnergeticFoodWebs.trophic_rank
using LinearAlgebra

# Functions:
## Generation rate:
# For each unique model_params object calculate the median generation rate of autrophs and output as gen object
#function get_gen_rate(i_b, r_vec, siz)
#    b = i_b * 2 # Calculate doubled biomass
#    gen_time = zeros(siz) # make gen_time object
#    for j in 1:siz#
#        # B(t) = B(0)exp(rt) - exponential growth equation
#        gen_time[j] = (1/r_vec[j]) * log(b[j]/i_b[j]) # Solve exponential growth equation for t - i.e. what is the doubling time (generation rate) of each species?
#    end
#    gen = median(gen_time) # Calculate median generation time
#    # Next step - multiple median generation time of autrophs by 2000 and run the model for that number - so 2000 generations of the median autroph in the system.
#    return(gen)
#end

function get_temp_df(gen, m, a, c, red_r, blue_r)
    # gen = generation rate of the median autroph - will vary by model
    # m = mean temperature - 20 degrees
    # a = amplitude of seasonal wave - 2 degrees (plus or minus)
    # c = number of waves in gen*2000 (gets put into b = 1/c)
    # red_r = red noise temporal autocorrelation value - always positive
    # blue_r = blue noise temporal autocorrelation value - always negative
    t = gen*100 # gives us the equivalent of 2000 generations
    x = collect(1:1:trunc(Int,t)) # seq from 1 - t
    temp_df = DataFrame(sc1 = 1:t, sc2 = 1:t, sc3 = 1:t, sc4 = 1:t, sc5 = 1:t, sc6 = 1:t, sc7 = 1:t, sc8 = 1:t, sc9 = 1:t, sc10 = 1:t, sc11 = 1:t, sc12 = 1:t, sc13 = 1:t, sc14 = 1:t, sc15 = 1:t, sc16 = 1:t)
    temp_df.sc1 = fill(20.0, trunc(Int, t)) # sc1 = No increase in temperature (fixed at 20 degrees) - Control
    temp_df.sc2 = collect(range(0, 40, length = trunc(Int, t))) # sc.2 = General increase (GI) in temperature (range - 0-40 degrees)
    # Generate seasonality (sine wave):
    b = 1/c # 1/number of units per wave/cycle
    for i in 1:trunc(Int, t)
        temp_df.sc3[i] = a*sin(2*pi*b*x[i]) + m # sc.3 = Seasonality (fixed at 20 degrees)
    end
    for i in 1:trunc(Int, t)
        temp_df.sc4[i] = a*sin(2*pi*b*x[i]) + temp_df.sc2[i] # sc.4 = Seasonality + GI in temperature (range - 0-40 degrees)
    end
    # Generate noise:
    white = randn(trunc(Int, t))
    red = zeros(trunc(Int, t))
    red[1] = white[1]
    for j in 1:(trunc(Int, t)-1)
        red[j+1] = red_r*red[j] + (1 - red_r^2)^0.5 * white[j+1]
    end
    blu = zeros(trunc(Int, t))
    blu[1] = white[1]
    for j in 1:(trunc(Int, t)-1)
        blu[j+1] = blue_r*blu[j] + (1 - blue_r^2)^0.5 * white[j+1]
    end
    # White noise scenarios:
    for i in 1:trunc(Int, t)
        temp_df.sc5[i] = white[i] + 20 # sc.5w = Stochasticity - white (fixed at 20 degrees)
    end
    for i in 1:trunc(Int, t)
        temp_df.sc6[i] = white[i] + temp_df.sc2[i] # sc.6w = Stochasticity (white) + GI in temperature (range - 0-40 degrees)
    end
    for i in 1:trunc(Int, t)
        temp_df.sc7[i] = white[i] + temp_df.sc3[i] # sc.7w = Stochasticity (white) + Seasonality (fixed at 20 degrees)
    end
    for i in 1:trunc(Int, t)
        temp_df.sc8[i] = white[i] + temp_df.sc4[i] # sc.8w = Stochasticity (white) + Seasonality + GI in temperature (range - 0-40 degrees)
    end
    # Red noise scenarios:
    for i in 1:trunc(Int, t)
        temp_df.sc9[i] = red[i] + 20 # sc.5r = Stochasticity - red (fixed at 20 degrees)
    end
    for i in 1:trunc(Int, t)
        temp_df.sc10[i] = red[i] + temp_df.sc2[i] # sc.6r = Stochasticity (red) + GI in temperature (range - 0-40 degrees)
    end
    for i in 1:trunc(Int, t)
        temp_df.sc11[i] = red[i] + temp_df.sc3[i] # sc.7r = Stochasticity (red) + Seasonality (fixed at 20 degrees)
    end
    for i in 1:trunc(Int, t)
        temp_df.sc12[i] = red[i] + temp_df.sc4[i] # sc.8r = Stochasticity (red) + Seasonality + GI in temperature (range - 0-40 degrees)
    end
    # Blue noise scenarios:
    for i in 1:trunc(Int, t)
        temp_df.sc13[i] = blu[i] + 20 # sc.5b = Stochasticity - blue (fixed at 20 degrees)
    end
    for i in 1:trunc(Int, t)
        temp_df.sc14[i] = blu[i] + temp_df.sc2[i] # sc.6b = Stochasticity (blue) + GI in temperature (range - 0-40 degrees)
    end
    for i in 1:trunc(Int, t)
        temp_df.sc15[i] = blu[i] + temp_df.sc3[i] # sc.7b = Stochasticity (blue) + Seasonality (fixed at 20 degrees)
    end
    for i in 1:trunc(Int, t)
        temp_df.sc16[i] = blu[i] + temp_df.sc4[i] # sc.8b = Stochasticity (blue) + Seasonality + GI in temperature (range - 0-40 degrees)
    end
    # Covert degrees to kelvin:
    temp_df = temp_df .+ 273.15

    return(temp_df) # output temp_df object
end

## Construct niche models:
# Will have to construct it by siz and con index
function niche_func(siz, con)
    mod = Any[] # List of lists
    for i in range(1, stop=8000)
        A_bool = EcologicalNetworks.nichemodel(siz, con) # use the niche model from Ecological Network to generate food webs
        A = Int.(A_bool.A) # convert the UnipartiteNetwork object to Matrix of 1s and 0s
        co = sum(A)/(size(A,1)^2) # calculate connectance
        if co == con
            push!(mod, A)
        end
    end
    mod_out = mod[1:100] #select first 100 with connectance = con
    return mod_out # output
end

function niche_func2(siz, con, rep)
    mod = Any[] # List of lists
    global l = length(mod)
    while l < rep
        A_bool = EcologicalNetworks.nichemodel(siz, con) # use the niche model from Ecological Network to generate food webs
        A = Int.(A_bool.A) # convert the UnipartiteNetwork object to Matrix of 1s and 0s
        co = sum(A)/(size(A,1)^2) # calculate connectance
        if co == con
            push!(mod, A)
        end
        global l = length(mod)
    end
    return mod # output
end

## Need two functions to make sure the network is stable and has enough species in it.
# (1) This function returns the difference between the log(biomass) of each species between two consecutive time steps
function diff_biomass(out)
    B = Array(out[:B])
    logB = log.(complex(vec(B)))
    S = size(B, 2)
    len = size(B,1)
    dlogB = zeros(len-1, S)
    for k in 1:S
        k_hold = @fastmath diff(log.(vec(B[:,k])))
        dlogB[1:end,k] = k_hold
    end
    return dlogB
end

#delta_log = diff_biomass(out)
# (2) This function identifies whether a species has reached a steady state and when (iteration number) steady state is reached
function is_steady(out; threshold::Float64 = 5e-4)
    dlog = diff_biomass(out)
    p = out[:p]
    S = size(dlog,2)
    is_ss = falses(S)
    first_iter_ss = zeros(S)
    for s in 1:S
        if in(s, p[:extinctions])
            is_ss[s] = true
            first_iter_ss[s] = findfirst(isnan.(dlog[:,s]))
        else
            try
                first_iter_ss[s] = findfirst(abs.(dlog[:,s]) .< threshold)
                is_ss[s] = true
            catch
                first_iter_ss[s] = NaN
                is_ss[s] = false
            end
        end
    end
    return is_ss, first_iter_ss
end
#id_steady, iteration_steady = is_steady(out)

# Crate new dir
function makeoutputdir(path)
    if !isdir(path)
        mkdir(path)
        mkdir(path * "Initial matrices")
        mkdir(path * "Non time varying test")
        mkdir(path * "Non time varying test/out_objects")
        #mkdir(path * "Non time varying test/biomasses")
        #mkdir(path * "Non time varying test/time")
        #mkdir(path * "Non time varying test/extinctions")
        #mkdir(path * "Non time varying test/species_traits")
        mkdir(path * "Time varying test")
        mkdir(path * "Time varying test/out_objects")
        #mkdir(path * "Time varying test/biomasses")
        #mkdir(path * "Time varying test/time")
        #mkdir(path * "Time varying test/extinctions")
        #mkdir(path * "Time varying test/species_traits")
    end
end

function normalize_matrix(A)
    A2 = transpose(A)
    colsum = sum(A2, dims = 1)
    colsum[colsum .== 0] .= 1
    normA = (A2'./vec(colsum))'
    return normA
end

function trophic_position(A)
    S = size(A,1)
    if S < 3
        return trophic_rank(A)
    else
        Mt = normalize_matrix(A)
        m = Int.(zeros(30,30))
        [m[i,i] = 1 for i in 1:S] #fill diag with 1
        detM = det(m .- Mt')
        if detM != 0
            tp = \(m .- Mt', repeat([1], S))
        else
            tmp = m 
            for i in 1:9
                tmp = tmp * Mt' .+ m
            end
            tp = tmp * repeat([1], S)
        end
        return tp
    end
end

function bodymass_calc(Z, A)
    trophiclevel = trophic_position(A)
    #dist_eps = Normal(0, 1)
    #eps_L = rand(dist_eps, size(A,1)) #we don't want all the species of a trophic level having exactly the same mass
    #m = 0.01 .* Z .^ (trophiclevel .- 1 .+ eps_L)
    m = 0.01 .* (Z .^ (trophiclevel .- 1))
    return m
end

function bodymass_calc2(Z, A)
    trophiclevel = trophic_rank(A)
    dist_eps = Normal(0, 1)
    eps_L = rand(dist_eps, size(A,1)) #we don't want all the species of a trophic level having exactly the same mass
    m = 0.01 .* (Z .^ ((trophiclevel .- 1) .+ eps_L))
    return m
end


#METABOLIC RATES

#Growth rate (producers)
function ScaleGrowth(M, T)
    r0 = exp(-15.68)
    βr = -0.25
    Er = -0.84
    T0 = 293.15 #20 celsius
    k = 8.617e-5
    return r0 .* (M .^ βr) .* exp(Er .* ((T0 .- T) ./ (k .* T .* T0)))
end

#Metabolic rate
function ScaleMetabolism(M, T)
    x0 = exp(-16.54)
    sx = -0.31
    Ex = -0.69
    T0 = 293.15
    k = 8.617e-5
    return x0 .* (M .^ sx) .* exp(Ex .* ((T0 .- T) ./ (k .* T .* T0)))
end

#Handling time
function ScaleHandling(m, T)
    h0 = exp(9.66)
    βres = -0.45
    βcons = 0.47
    Eh = 0.26
    T0 = 293.15
    k = 8.617E-5
    boltz = exp(Eh * ((T0-T)/(k*T*T0)))
    hij = zeros(length(m), length(m))
    for i in eachindex(m) #i = rows => consumers
      for j in eachindex(m) #j = cols => resources
        mcons = m[i] ^ βcons #mass scaled for cons
        mres = m[j] ^ βres #mass scaled for res
        hij[i,j] = h0 * mres * mcons * boltz
      end
    end
    return hij
end

#Attack rate
function ScaleAttack(m, T)
    a0 = exp(-13.1)
    βres = 0.25 #resource
    βcons = -0.8 #consumer
    Ea = -0.38
    T0 = 293.15
    k = 8.617E-5
    boltz = exp(Ea * ((T0-T)/(k*T*T0)))
    aij = zeros(length(m), length(m))
    for i in eachindex(m) #i = rows => consumers
      for j in eachindex(m) #j = cols => resources
        mcons = m[i] ^ βcons #mass scaled for cons
        mres = m[j] ^ βres #mass scaled for res
        aij[i,j] = a0 * mres * mcons * boltz
      end
    end
    return aij
end

"""
ScaleRates(p, k0 = 10, prod_metab = false)

Update the biological rates of a BioEnergeticFoodWebs model_parameters object (p) 

- p : an object created by BioEnergeticFoodWebs.model_parameters
- k0 : (Int or Float) scaling coefficient for the carrying capacity (default = 10)
- prod_metab : (bool) do basal species have metabolic losses (default = false) 
"""
function ScaleRates!(p; k0 = 10, prod_metab = false)
    T = p[:T]
    M = p[:bodymass]
    ri = ScaleGrowth(M, T)
    r = ri .* p[:is_producer]
    xi = ScaleMetabolism(M, T)
    x = prod_metab ? xi : xi .* .!p[:is_producer]
    ar = ScaleAttack(M, T)
    ht = ScaleHandling(M, T)
    ki = carrying(M, k0, T)
    k = ki .* p[:is_producer]
    p[:r] = r
    p[:x] = x
    p[:ar] = ar
    p[:ht] = ht
end

#Maximum consumption
function MaxCons(ht, x)
    y = 1 ./ ht
    y_norm = y ./ x
    y_norm[y_norm .== Inf] .= 0.0
    return y_norm
end

#Half saturation density
function HalfSaturation(ht, ar)
    hs = 1 ./ (ht .* ar)
    hs[hs .== Inf] .= 0.0
    return hs
end

#Carrying capacity
function carrying(m, k0, T)
    βk = 0.28
    Ek = 0.71 
    return k0 .* (m .^ βk) .* exp.(Ek .* (293.15 .- T ) ./ (8.617e-5 .* T .* 293.15))
end

#initial biomass
function B0(A, kcap)
    b0 = zeros(size(A,1))
    is_producer = sum(A, dims = 2) .== 0
    carrying_prod = kcap[vec(is_producer)]
    kmean = mean(carrying_prod)
    b0[vec(is_producer)] .= kcap[vec(is_producer)]
    b0[vec(.!is_producer)] .= kmean/8
    return b0
end

#=
Plotting + Updatin the matrix
=#

using Plots

"""
This function is used internally by the webplot function, it just transform 
an interaction matrix into a list of interaction.
"""
function matrix_to_list(A)
    idx = findall(A .== 1)
    from_sp = (i->i[1]).(idx)
    to_sp = (i->i[2]).(idx)
    return hcat(from_sp, to_sp)
end

"""
This function is also used internally by the webplot function. 
"""
function invert(v, S)
    1 .+ (S .- v)
end

"""
Plot the interaction matrix :tada: with consumer as either rows or columns :tada: 
whichever you prefer! 

    `plot(A, consasrow)`

- use `plot(A, true)` to have consumers in rows
- use `plot(A, false)` to have consumers in columns
"""
function webplot(A::Array{Int64,2}; consasrow::Bool = true)
    S = size(A,1)
    if !consasrow
        floatA = Float64.(A)
        listA = matrix_to_list(floatA')
        plt = scatter(listA[:,2], invert(listA[:,1], S)
            , ms = 3, c = :black
            , size = (300,300), leg = false
            , framestyle = :box, ticks = ([1.5:1:S;], repeat([""], S)), foreground_color_axis = :white
            , xlims = (0.5,S+0.5), ylims = (0.5,S+0.5))    
        plot!([1,S], [S,1], c = :black, linestyle = :dash)
        xlabel!("consumers")
        ylabel!("resources")
    else
        floatA = Float64.(A)
        listA = matrix_to_list(floatA)
        plt = scatter(listA[:,2], invert(listA[:,1],S)
        , ms = 3, c = :black
        , size = (300,300), leg = false
        , framestyle = :box, ticks = ([1.5:1:S;], repeat([""], S)), foreground_color_axis = :white
        , xlims = (0.5,S+0.5), ylims = (0.5,S+0.5))    
        plot!([1,S], [S,1], c = :black, linestyle = :dash)
        ylabel!("consumers")
        xlabel!("resources")
    end
    return plt
end

"""
Plot the interaction matrix :tada: with consumer as either rows or columns :tada: 
whichever you prefer! 

    `plot(out, consasrow, colorextinct)`

- `out` is the BEFWM simulation output
- `consasrow` (Bool) specifies whether we want consumers as rows (true) or column (false)
- `colorextinct` (Bool) specifies whether we want to shade the extinct species
"""
function webplot(out::Dict{Symbol,Any}; consasrow::Bool = true, colorextinct::Bool = false)
    A = out[:p][:A]
    S = size(A,1)
    id_alive = trues(S)
    id_alive[out[:p][:extinctions]] .= false
    wp = webplot(A; consasrow = consasrow)
    if colorextinct
        ii = invert(i,S)
        for i in out[:p][:extinctions]
            if consasrow
                shp1 = Shape([0.5, S+0.5, S+0.5, 0.5], [ii-0.5, ii-0.5, ii+0.5, ii+0.5])
                shp2 = Shape([i-0.5, i+0.5, i+0.5, i-0.5], [0.5, 0.5, S+0.5, S+0.5])    
            else
                shp1 = Shape([i-0.5, i-0.5, i+0.5, i+0.5], [0.5, S+0.5, S+0.5, 0.5])
                shp2 = Shape([0.5, 0.5, S+0.5, S+0.5], [ii-0.5, ii+0.5, ii+0.5, ii-0.5])    
            end
            plot!(shp1, c = :grey, opacity=.3, lw = 0)
            plot!(shp2, c = :grey, opacity=.3, lw = 0)
        end
    end
    return wp
end

"""
Updates the interaction matrix based on the vector of extinct species 

    `updateA(out)`

- `out` is the BEFWM simulation output

Returns a matrix of dimension S where S is the number of persistent species. This function also
checks for disconnected consumers. If it detect any, it will print a message and return a vector
with the list of disconnected species instead of the updated matrix.
"""
function updateA(out)
    A = out[:p][:A]
    S = size(A,1)
    id_alive = trues(S)
    id_alive[out[:p][:extinctions]] .= false
    Anew = A[id_alive, id_alive]
    #check for status change (consumers can't become producers)
    id_th_prod = sum(A, dims = 2) .== 0
    alive_prod = findall((id_th_prod) .& (id_alive))
    alive_prod = [i[1] for i in alive_prod]
    discosp = []
    for i in alive_prod
        if i ∉ findall(out[:p][:is_producer])
            println("/!\\ disconnected consumer $i identified as producer")
            push!(discosp, i)
        end
    end
    toreturn = length(discosp) == 0 ? Anew : discosp
    return toreturn
end