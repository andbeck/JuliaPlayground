using Plots, DiffEqSensitivity
using DifferentialEquations, ParameterizedFunctions
using Statistics

# plot set
gr()

function f(du,u,p,t)
  du[1] = p[1]*u[1] - p[2]*u[1]*u[2] #prey
  du[2] = -p[3]*u[2] + p[4]*u[1]*u[2] #predator
end

u0 = [1.0;1.0]
tspan = (0.0,10.0)
p = [1.5,1.0,3.0,1.0]
prob = ODEProblem(f,u0,tspan,p)
t = collect(range(0, stop=10, length=200))

#define and solve problem with defaults
sol = solve(prob)

# time series plot
plot(sol)
# phase plane plot
plot(sol, vars = (1,2))

f1 = function (p)
  prob1 = remake(prob;p=p)
  sol = solve(prob1,Tsit5();saveat=t)
  [mean(sol[1,:]), maximum(sol[2,:])]
end


m = gsa(f1,Morris(total_num_trajectory=1000,num_trajectory=150),[[1,5],[1,5],[1,5],[1,5]])
m_sobol = gsa(f1,Sobol(),N = 1,[[1,5],[1,5],[1,5],[1,5]])


m.means
scatter(m.means[1,:], m.variances[1,:],series_annotations=[:a,:b,:c,:d],color=:gray)
scatter(m.means[2,:], m.variances[2,:],series_annotations=[:a,:b,:c,:d],color=:gray)
