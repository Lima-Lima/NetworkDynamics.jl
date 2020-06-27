using LightGraphs
using OrdinaryDiffEq
using NetworkDynamics
using Debugger

N = 16 # number of nodes
k = 3  # average degree
g = barabasi_albert(N, k)

#=
  Berner, Rico, Eckehard Schöll, and Serhiy Yanchuk.
  "Multiclusters in Networks of Adaptively Coupled Phase Oscillators."
  SIAM Journal on Applied Dynamical Systems 18.4 (2019): 2227-2266.
=#


@inline function kuramoto_plastic_edge!(de, e, v_s, v_d, p, t)
    # Source to Destination coupling

    # Simplified for a symmetric example, as a proof of concept.
    
    # Here is the differential equation. It is numerically integrated with the
    # ode solver. Just like in regular networkdynamics, the differential de[]
    # states are appended to the end of the vertex state vector and this is
    # passed to the integrator.
    @bp
    de[1] = -ϵ * (sin(v_s[1] + v_d[1]) + e[1])
   
    # Now, we want to offer some algebraic expression to the verticies. In this
    # POC, the algebraic expression is returned. See the update to the main
    # network loop over edges, where this then populates the new (non-state)
    # array expr
    return e[1] * sin(v_s[1] - v_d[1])
end

@inline function kuramoto_plastic_vertex!(dv, v, e_s, e_d, p, t)
    dv .= 0
    for e in e_s
        dv .-= e[1]
    end
    for e in e_d # other direction is stored at other index
        dv .-= e[1]
    end
end


    # Parameter definitions
ϵ = 0.1
α = .2π
β = -.95π

# NetworkDynamics Setup
plasticvertex = ODEVertex(f! = kuramoto_plastic_vertex!, dim =1)
plasticedge = ODEEdge(f! = kuramoto_plastic_edge!, dim=1, sym=[:k])
kuramoto_plastic! = network_dynamics(plasticvertex, plasticedge, g)

# ODE Setup & Solution
x0_plastic        = vcat(randn(N), ones(ne(g)))
tspan_plastic     = (0., 200.)
params_plastic    = (nothing, nothing)
prob_plastic      = ODEProblem(kuramoto_plastic!, x0_plastic, tspan_plastic, params_plastic)
sol_plastic       = solve(prob_plastic, Rosenbrock23(), abstol = 1e-3, reltol = 1e-3)

# Plotting
#v_idx = idx_containing(kuramoto_plastic!, :v)

#plot(sol_plastic, vars=v_idx, legend=false, ylabel=L"\theta")

# Shows Coupling terms
e_idx = idx_containing(kuramoto_plastic!, :k)
# plot!(sol_plastic, vars=e_idx, legend=false, color=:black, linewidth=0.001)
