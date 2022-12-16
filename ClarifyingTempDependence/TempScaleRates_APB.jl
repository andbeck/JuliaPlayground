#Growth rate (producers)
function ScaleGrowth(M, T)
    r0 = exp(-15.68)
    βr = -0.25
    Er = -0.84
    T0 = 293.15 #20 celsius
    k = 8.617e-5
    return r0 .* (M .^ βr) .* exp(Er .* ((T0 .- T) ./ (k .* T .* T0)))
end

function ScaleGrowth2(M, T)
    r0 = exp(-15.68)
    βr = -0.25
    Er = -0.84
    T0 = 298.15 #25 celsius
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

function ScaleAttack2(m, T)
    a0 = exp(-13.1)
    βres = 0.25 #resource
    βcons = -0.8 #consumer
    Ea = -0.38
    T0 = 298.15
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

#Carrying capacity
function carrying(m, k0, T)
    βk = 0.28
    Ek = 0.71
    return k0 .* (m .^ βk) .* exp.(Ek .* (293.15 .- T ) ./ (8.617e-5 .* T .* 293.15))
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
