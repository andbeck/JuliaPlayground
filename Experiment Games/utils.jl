using Statistics
import BioEnergeticFoodWebs.trophic_rank
using LinearAlgebra

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
