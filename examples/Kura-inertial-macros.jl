using NetworkDynamics
using LightGraphs
using OrdinaryDiffEq
# using DiffEqOperators
using LinearAlgebra

N = 10
g = barabasi_albert(N, 5)

ω = Array{Float64}(1:N) ./ N

@Edge kuraEdge(k) begin 
    # Prep work
end [[F]] begin
    # Dynamics
F = k * sin(v_s[1] - v_d[1]) #TODO:
#F = k * sin(v_s[:ϕ] - v_d[:ϕ]) #TODO: symbolic indexing
end

@Vertex kuraVertex(SetPower,D) begin
    I # MassMatrix
end begin
end [[ϕ,ω],[ω,dω]] begin
    dϕ = ω
    dω = SetPower - D*ω
    for e in e_s
        dω -= e[1] #TODO: this should
        # allow symbolic indexing (e[:F])
    end
    for e in e_d
        dω += e[1]
    end
end
tspan = (0.0,1.0)
vertex_list = [construct(kuraVertex(1.0,0.5)) for v in vertices(g)]
edge_list = [construct(kuraEdge(1.0)) for e in edges(g)]

p = (ω, nothing)

kur_network_nd = network_dynamics(vertex_list,edge_list, g)

x0_nd = 0.1 .* Array{Float64}(1:2N)

prob_nd = ODEProblem(kur_network_nd, x0_nd, tspan, nothing)

sol_nd = solve(prob_nd, Tsit5())
sol_nd = solve(prob_nd, Tsit5())
@time sol_nd = solve(prob_nd, Tsit5())



