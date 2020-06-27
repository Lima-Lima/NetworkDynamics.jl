using LightGraphs
using OrdinaryDiffEq
using NetworkDynamics
#using Plots, LaTeXStrings

N = 16 # number of nodes
k = 3  # average degree
g = barabasi_albert(N, k)

#=
  Berner, Rico, Eckehard Schöll, and Serhiy Yanchuk.
  "Multiclusters in Networks of Adaptively Coupled Phase Oscillators."
  SIAM Journal on Applied Dynamical Systems 18.4 (2019): 2227-2266.
=#


@inline function kuramoto_plastic_edge!(de, e, ex, v_s, v_d, p, t)
    # Source to Destination coupling

    ϵ = 0.1
    α = .2π
    β = -.95π

    # The coupling function is modeled by a differential algebraic equation with mass matrix 0
    # 0 * de[1] = e[2] * sin(v_s[1] - v_d[1] + α) / N - e[1] is equivalent to e[1] = e[2] * sin(v_s[1] - v_d[1] + α) / N

    ex[1] =  e[1] * sin(v_s[1] - v_d[1] + α) / N  # Not updating the dynamic state vector
    de[1] = - ϵ * (sin(v_s[1] - v_d[1] + β) + e[1]) # Updating the dynamics state vector

    # Destination to source coupling
    # since the coupling function is not symmetric we have to compute the other direction as well

    ex[2] =  e[2] * sin(v_d[1] - v_s[1] + α) / N
    de[2] = - ϵ * (sin(v_d[1] - v_s[1] + β) + e[2])
    nothing
end

@inline function kuramoto_plastic_vertex!(dv, v, e_s, e_d, p, t)
    dv .= 0
    for e in e_s
        dv .-= e[1]
    end
    for e in e_d # other direction is stored at other index
        dv .-= e[2]
    end
end

function run()
    # Parameter definitions
# NetworkDynamics Setup
plasticvertex = ODEVertex(f! = kuramoto_plastic_vertex!, dim =1)

plasticedge = ODEEdge(f! = kuramoto_plastic_edge!, dim=2, expr_dim=2, sym=[:ks,:kd], expr_sym=[:es,:ed])
kuramoto_plastic! = network_dynamics(plasticvertex, plasticedge, g)

# ODE Setup & Solution
x0_plastic        = vcat(randn(N), ones(2ne(g)))
tspan_plastic     = (0., 200.)
params_plastic    = (nothing, nothing)
prob_plastic      = ODEProblem(kuramoto_plastic!, x0_plastic, tspan_plastic, params_plastic)
sol_plastic       = solve(prob_plastic, Rosenbrock23(), abstol = 1e-3, reltol = 1e-3)

# Plotting
v_idx = idx_containing(kuramoto_plastic!, :v)

#plot(sol_plastic, vars=v_idx, legend=false, ylabel=L"\theta")

# Shows Coupling terms
e_idx = idx_containing(kuramoto_plastic!, :k)

return sol_plastic
#plot!(sol_plastic, vars=e_idx, legend=false, color=:black, linewidth=0.001)
end
