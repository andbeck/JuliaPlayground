using Plots, DiffEqSensitivity, DifferentialEquations, ParameterizedFunctions

# plot set
gr()

# model
function f(du,u,p,t)
  du[1] = p[1]*u[1] - p[2]*u[1]*u[2]
  du[2] = -3*u[2] + u[1]*u[2]
end
u0 = [1.0;1.0]
tspan = (0.0,10.0)
p = [1.5,1.0]
prob = ODEProblem(f,u0,tspan,p)
t = collect(range(0, stop=10, length=200))

# Morris

m = DiffEqSensitivity.morris_sensitivity(prob,Tsit5(),t,[[1,5],[0.5,5]],[10,10],
len_trajectory=1500,
total_num_trajectory=1000,
num_trajectory=150)

m.means
m.variances

stdv1 = sqrt.(m.variances[1])
p = plot(m.means[1]', yerror=stdv1)

s0 = sobol_sensitivity(prob,Tsit5(),t,[[1,5],[0.5,5]],N,0)
Out[8]: 2-element Array{Array{Float64,2},1}:
 [NaN 0.507831 … 1.00731 1.00436; NaN 1.92336 … 0.732384 0.730945]
 [NaN 0.47214 … 0.676224 0.681525; NaN -1.68656 … 0.879557 0.877603]

s1 = sobol_sensitivity(prob,Tsit5(),t,[[1,5],[0.5,5]],N,1)
Out[9]: 2-element Array{Array{Float64,2},1}:
 [NaN 0.39537 … 0.341697 0.343645; NaN -2.06101 … 0.10922 0.106976]
 [NaN 0.652815 … 0.00910675 0.00815206; NaN 5.24832 … 0.296978 0.296639]

s2 = sobol_sensitivity(prob,Tsit5(),t,[[1,5],[0.5,5]],N,2)
Out[10]: 1-element Array{Array{Float64,2},1}:
 [NaN -0.0596478 … 0.652303 0.657847; NaN -1.84504 … 0.645139 0.620036]

p1 = bar(["a","b"],[s0[1][end-2],s0[2][end-2]],color=[:red,:blue],title="Total Order Indices at t=9.949748743718592",legend=false)
p2 = bar(["a","b"],[s1[1][end-2],s1[2][end-2]],color=[:red,:blue],title="First Order Indices at t=9.949748743718592",legend=false)
p3 = bar(["a","b"],[s0[1][3],s0[2][3]],color=[:red,:blue],title="Total Order Indices at t=0.05025125628140704",legend=false)
p4 = bar(["a","b"],[s1[1][3],s1[2][3]],color=[:red,:blue],title="First Order Indices at t=0.05025125628140704",legend=false)
plo = plot(p1,p2,p3,p4,layout=(4,1),size=(600,500))
