using MacroTools
using LinearAlgebra

function definestruct(name, parameters)
    return struct_def = Expr(:struct, false, :($name),
                      Expr(:block, parameters...)
                     )
end

macro Vertex(typedef,statements...)
    @capture(typedef, name_(parameters__))
    
    # Get satements and find block ranges 
    stmts = MacroTools.prewalk(rmlines,statements[1])
    bidx = [x for x in enumerate(stmts.args) if begin 
            typeof(x[2])==Expr && x[2].head == :macrocall end]
    body_statements_range = [(bidx[n][2],bidx[n][1]+1,bidx[n+1][1]-1) for n in 1:length(bidx)-1]
    append!(body_statements_range,[(bidx[end][2],bidx[end][1]+1,length(stmts.args))])
    
    # Build struct definition expression
    struct_expr = definestruct(name,parameters)
     
    # The construct function signature dispatch
    constructcall = :(construct(vtx::$(name)))
    # The construct function body
    paramlocal = map(sym -> :($sym = vtx.$sym), parameters)
    constructbody = quote end
    append!(constructbody.args,paramlocal)
    

    # Construct Function (completed later)
    constructfunction = Expr(:function, constructcall, constructbody)
    
    #Vertex creation
    
    # Decompose the statments
    states_stmts = [:(nothing)]
    prepare_stmts = [:(nothing)]
    massmatrix_stmts = [:(nothing)]
    dynamics_stmts = [:(nothing)]
    for section in body_statements_range
        if section[1].args[1] == Symbol("@states")
            states_stmts = stmts.args[section[2]:section[3]]
        elseif section[1].args[1] == Symbol("@prepare")
            prepare_stmts = stmts.args[section[2]:section[3]]
        elseif section[1].args[1] == Symbol("@massmatrix")
            massmatrix_stmts = stmts.args[section[2]:section[3]]
        elseif section[1].args[1] == Symbol("@dynamics")
            dynamics_stmts = stmts.args[section[2]:section[3]]
        else
            # TODO: Throw error on unknown section header
        end
    end
    
    # states
    state_vars = []
    state_vars = [ex.args[1] for ex in states_stmts if 
                  typeof(ex)==Expr] #TODO: length check
    dstate_vars = []
    dstate_vars = [ex.args[2] for ex in states_stmts if 
                   typeof(ex)==Expr]
    algebraic_vars = [] 
    algebraic_vars = [sym for sym in states_stmts if 
                      typeof(sym)==Symbol]
    
    # mass matrix 
    if massmatrix_stmts == [:(nothing)]
        mm = ones(length(state_vars) + length(algebraic_vars))
        mm[length(state_vars)+1:end] .= 0
        massmatrix = Diagonal(mm)
    end
    
    # prepare statements
    append!(constructbody.args,prepare_stmts)
    
    # RHS function
    if isempty(state_vars) #then it's a static vertex
        rhssig = :(rhs!(x,e_s,e_d,p,t))
    else 
        rhssig = :(rhs!(dx,x,e_s,e_d,p,t))
    end
    rhsbody = quote end
    
    # state vector assignments (state1= x[0], dstate1 = dxh 
    state_stmts = [:($sym = x[$index]) for (index,sym) in 
                  enumerate([state_vars; algebraic_vars])]
    if isempty(state_vars) # then it's a static vertex 
        dstate_stmts = [:(x[$index] = $sym) for (index,sym) in 
                   enumerate([dstate_vars; algebraic_vars])]
    else 
        dstate_stmts = [:(dx[$index] = $sym) for (index,sym) in 
                   enumerate([dstate_vars; algebraic_vars])]
    end 
    append!(rhsbody.args, state_stmts)  
    append!(rhsbody.args, dynamics_stmts)  
    append!(rhsbody.args, dstate_stmts) #TODO: Catch if there is no assignment to a 
                                        # dstate in body
    append!(rhsbody.args, [:(nothing)])
    rhsfunc = Expr(:function,rhssig,rhsbody)

    # construction call
    if isempty(dstate_vars) # then it's not an ode
        nd_constructor = :(StaticVertex(f! =rhs!,
                            dim=$(length(algebraic_vars)),
                            sym=$(algebraic_vars)))
    else 
        nd_constructor = :(ODEVertex(f! =rhs!,
                            dim=$(length(state_vars)+length(algebraic_vars)),
                            mass_matrix=$massmatrix,
                            sym=$([state_vars; algebraic_vars])))
    end
    append!(constructfunction.args[2].args,[rhsfunc,nd_constructor])
    
    ex = quote
        $(struct_expr)
        $(constructfunction)
    end
    @show ex
    return esc(ex)
end


macro Edge(typedef,statements...)
    @capture(typedef, name_(parameters__))
    # Read statements and determine where keywords define blocks 
    stmts = MacroTools.prewalk(rmlines,statements[1])
    bidx = [x for x in enumerate(stmts.args) if begin 
            typeof(x[2])==Expr && x[2].head == :macrocall end]
    body_statements_range = [(bidx[n][2],bidx[n][1]+1,bidx[n+1][1]-1) for n in 1:length(bidx)-1]
    append!(body_statements_range,[(bidx[end][2],bidx[end][1]+1,length(stmts.args))])
    
    # The expression to create the struct
    struct_expr = definestruct(name,parameters)
    
    # The construct function signature dispatch
    constructcall = :(construct(vtx::$(name)))
    # The construct function body
    paramlocal = map(sym -> :($sym = vtx.$sym), parameters)
    constructbody = quote end
    append!(constructbody.args,paramlocal)
    
    # Construct Function (completed later)
    constructfunction = Expr(:function, constructcall, constructbody)
    
    #Edge creation
    
    # Decompose the statments
    states_stmts = [:(nothing)]
    prepare_stmts = [:(nothing)]
    massmatrix_stmts = [:(nothing)]
    dynamics_stmts = [:(nothing)]
    for section in body_statements_range
        if section[1].args[1] == Symbol("@states")
            states_stmts = stmts.args[section[2]:section[3]]
        elseif section[1].args[1] == Symbol("@prepare")
            prepare_stmts = stmts.args[section[2]:section[3]]
        elseif section[1].args[1] == Symbol("@massmatrix")
            massmatrix_stmts = stmts.args[section[2]:section[3]]
        elseif section[1].args[1] == Symbol("@dynamics")
            dynamics_stmts = stmts.args[section[2]:section[3]]
        else
            # TODO: Throw error on unknown section header
        end
    end
    
    # states
    state_vars = []
    state_vars = [ex.args[1] for ex in states_stmts if 
                  typeof(ex)==Expr] #TODO: length check
    dstate_vars = [] 
    dstate_vars = [ex.args[2] for ex in states_stmts if 
                   typeof(ex)==Expr]
    algebraic_vars = [] 
    algebraic_vars = [sym for sym in states_stmts if 
                      typeof(sym)==Symbol]
    
    # mass matrix 
    if massmatrix_stmts == [:(nothing)]
        mm = ones(length(state_vars) + length(algebraic_vars))
        mm[length(state_vars)+1:end] .= 0
        massmatrix = Diagonal(mm)
    end
    
    # prepare statements
    append!(constructbody.args,prepare_stmts)
    
    # RHS function
    if isempty(state_vars) # it's a static edge
        rhssig = :(rhs!(e,v_s,v_d,p,t))
    else 
        rhssig = :(rhs!(de,e,v_s,v_d,p,t))
    end 
    rhsbody = quote end
    
    # state vector assignments (state1= x[0], dstate1 = dxh 
    state_stmts = [:($sym = e[$index]) for (index,sym) in 
                  enumerate([state_vars; algebraic_vars])]
    if isempty(state_vars) # it's a static edge
        dstate_stmts = [:(e[$index] = $sym) for (index,sym) in 
                   enumerate(algebraic_vars)]
    else 
        dstate_stmts = [:(de[$index] = $sym) for (index,sym) in 
                   enumerate([dstate_vars; algebraic_vars])]
    end 
    append!(rhsbody.args, state_stmts)  
    append!(rhsbody.args, dynamics_stmts)  
    append!(rhsbody.args, dstate_stmts) #TODO: Catch if there is no assignment to a 
                                        # dstate in body
    append!(rhsbody.args, [:(nothing)])
    rhsfunc = Expr(:function,rhssig,rhsbody)

    # static edge constructor call
    if isempty(dstate_vars) # then it's not an ode
        nd_constructor = :(StaticEdge(f! = rhs!,
                            dim=$(length(algebraic_vars)),
                            sym=$(algebraic_vars)))
    # ODE construction call
    else 
        nd_constructor = :(ODEEdge(f! = rhs!,
                            dim=$(length(state_vars)+length(algebraic_vars)),
                            mass_matrix=$massmatrix,
                            sym=$([state_vars; algebraic_vars])))
    end 
    
    append!(constructfunction.args[2].args,[rhsfunc,nd_constructor])
    
    ex = quote
        $(struct_expr)
        $(constructfunction)
    end
    @show ex
    return esc(ex)
end


#@Vertex TestV(Params...) begin
#       @prepare
#       
#       d = e + f
#       @states
#       x,dx
#       y,dy
#       z,dz
#       @dynamics
#       a = x+dx
#       dz = y + z
# end
#@Edge TestE(Params...) begin
#       @prepare
#       
#       d = e + f
#       @states
#       x,dx
#       y,dy
#       z,dz
#       @dynamics
#       a = x+dx
#       dz = y + z
# end
