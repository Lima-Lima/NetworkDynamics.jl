using MacroTools

function vertex_creator!(constructcall,name,params,massmatrix,states,body)
    rhssig = :(rhs!(dx,x,e_s,e_d,p,t))
    #rhssig = Meta.parse("(this::$name)(dx,x,e_s,e_d,p,t)")
    rhsbody = quote end
    rhsbody.args[1] = body.args[1]
    # Prepare the parts
    state_vars = [:($sym = x[$index]) for (index,sym) in 
                  enumerate(states.statevec)]
    dstate_update = [:($sym = dx[$index]) for (index,sym) in 
                  enumerate(states.dstatevec)]
    append!(rhsbody.args, state_vars)
    append!(rhsbody.args, body.args)
    append!(rhsbody.args, dstate_update)
    rhsfunc = Expr(:function, rhssig, rhsbody)
    ode_constructor = :(ODEVertex(f! = rhs!, 
                                  dim = $(states.statevec),
                                  mass_matrix=$massmatrix,
                                  sym = $states.statevec))
    append!(constructcall.args[2].args, [rhsfunc, ode_constructor])
    nothing
end

function getstatevector(states)
    return (statevec  = [state.args[1] for state in states.args],
            dstatevec = [state.args[2] for state in states.args])
end

function definestruct(name, parameters)
    return struct_def = Expr(
                             :struct, false, :($name),
                      Expr(:block, parameters...)
                     )
end

macro Vertex(typedef,massmatrix,prep,states,body)
    return create(typedef, massmatrix, prep, states, body)
    #return create(typedef, massmatrix, prep, states, bodys..)
end

function create(typedef, massmatrix, prep, states, body)
    @capture(typedef,name_(parameters__))
    struct_expr = definestruct(name,parameters)
    constructcall = :(construct(vtx::$(name)))
    #constructcall = Meta.parse("(this::$name)(dx,x,e_s,e_d,p,t)")
    paramlocal = map(sym -> :($sym = vtx.$sym), parameters)
    constructbody = quote end
    append!(constructbody.args,paramlocal)
    append!(constructbody.args,prep.args)
    constructfunction = Expr(:function, constructcall, constructbody)
    expression = vertex_creator!(constructfunction,name, massmatrix, prep,getstatevector(states),body)
    ex = quote
        $(struct_expr)
        $(constructfunction)
    end
    return esc(ex)
end


#macro Vertex(typedef,massmatrix,prep,states,body)
#    @capture(typedef, name_(parameters__))
#    states_tuple = getstatevector(states)
#    struct_def = definestruct(name,parameters)
#    massmatrix = massmatrix == nothind ? I : massmatrix
#    return vertex_creator(params,massmatrix,states,body)
#end

#@DynamicNode SwingEq(H, P, D, Ω) begin
#    @assert D >= 0 "damping (D) should be >=0"
#    @assert H > 0 "inertia (H) should be >0"
#    Ω_H = Ω / (2*H)
#end [[ω, dω]] begin
#    p = real(u * conj(i))
#    dϕ = ω # dϕ is only a temp variable that Julia should optimize out
#    du = u * im * dϕ
#    dω = (P - D*ω - p)*Ω_H
#end
#@ODEVertex Test(Params...) begin
#    #preparatory stuff
#end begin DynamicsLevel 
#    @states [state1, state2] 
#    @algebraic [algebraic1, algebraic2]
#    #dynamics function
#end
#
#@ODEMultilayerVertex Test(Params...) begin
#    #preparatory stuff
#end begin DynamicsLevel 
#    @states [state1, state2] 
#    @algebraic [algebraic1, algebraic2]
#    #dynamics function
#end begin StochasticsLevel [state1, state2]
#    #dynamics function
#end begin CommunicationsLevel [state1, state2]
#    #comms functions
#end begin
#    # combination level
#end
