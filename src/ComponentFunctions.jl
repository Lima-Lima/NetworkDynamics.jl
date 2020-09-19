"""
Together with Constructors this module forms the backbone of the core API.
It provide the basic types to construct Arrays of VertexFunction and
EdgeFunction which can be handled by network_dynamics.
"""
module ComponentFunctions

using LinearAlgebra
# using SparseArrays

import Base.convert
import Base.promote_rule


export VertexFunction
export EdgeFunction
export StaticVertex
export StaticEdge
export ODEVertex
export ODEEdge
export ODEEAEdge
export DDEVertex
export StaticDelayEdge
#export DDEEdge
"""
Abstract supertype for all vertex functions.
"""
abstract type VertexFunction end
"""
Abstract supertype for all edge functions.
"""
abstract type EdgeFunction end

"""
    StaticVertex(f!, dim, sym)

Wrapper that ensures compatibility of a **mutating** function **`f!`** with
the key constructor `network_dynamics`.

**`f!`**  describes the local behaviour at a static node and has to respect
the following calling syntax

```julia
f!(v, e_s, e_t, p, t) -> nothing
```

Here  `v`, `p` and `t` are the usual arguments, while
`e_s` and `e_d` are arrays containing the edges for which the
described vertex is the source or the destination respectively.

**`dim`** is the number of independent variables in the vertex equations and
**`sym`** is an array of symbols for these variables.

For more details see the documentation.
"""
@Base.kwdef struct StaticVertex{T} <: VertexFunction
    f!::T # (v, e_s, e_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    sym=[:v for i in 1:dim] # Symbols for the dimensions
end

"""
    StaticEdge(f!, dim, sym)

Wrapper that ensures compatibility of a **mutating** function **`f!`** with
the key constructor `network_dynamics`.

**`f!`**  describes the local behaviour at a static edge and has to respect
the following calling syntax

```julia
f!(e, v_s, v_t, p, t) -> nothing
```

Here  `e`, `p` and `t` are the usual arguments, while
`v_s` and `v_d` are arrays containing the vertices which are
the source and destination of the described edge.

**`dim`** is the number of independent variables in the edge equations and
**`sym`** is an array of symbols for these variables.

For more details see the documentation.
"""
# The StaticEdge function returns an ODEEAEdge
function StaticEdge{T}(f!::T,dim::Int;sym=[:e for i in 1:dim])
    fea! = (ded,ed,eS,v_s,v_d,p,t) -> f!(ded,ed,v_s,v_d,p,t)
    return ODEEAEdge(fea!,ed_dim=0,eS_dim=dim,eS_sym=sym)
end

"""
    ODEVertex(f!, dim, mass_matrix, sym)

Wrapper that ensures compatibility of a **mutating** function **`f!`** with
the key constructor `network_dynamics`.

**`f!`**  describes the local behaviour at a dynamic node and has to respect
the following calling syntax

```julia
f!(dv, v, e_s, e_t, p, t) -> nothing
```

Here `dv`, `v`, `p` and `t` are the usual ODE arguments, while
`e_s` and `e_d` are arrays containing the edges for which the
described vertex is the source or the destination respectively.

**`dim`** is the number of independent variables in the vertex equations and
**`sym`** is an array of symbols for these variables.
**`mass_matrix`** is an optional argument that defaults to the identity
matrix `I`. If a mass matrix M is given the system `M * dv = f!` will be
solved.

For more details see the documentation.
"""
@Base.kwdef struct ODEVertex{T} <: VertexFunction
    f!::T # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    mass_matrix=I # Mass matrix for the equation
    sym=[:v for i in 1:dim] # Symbols for the dimensions
end
"""
    ODEEAVertex(f!, dim, mass_matrix, sym)

    and f! must be
    
    `f!(dv,v,e_s_v,e_d_v,eS_s_v,eS_d_v,p,t)`
    with differential vertex states `dv` and `v`, 
    `ed_s_v` differential edge states source
    `ed_d_v` differential edge states dest
    `eS_s_v` explicit edge states source
    `eS_d_v` explicit edge states dest
    `p` params
    `t` time

"""
@Base.kwdef struct ODEEAVertex{T} <: VertexFunction
    f!::T # The function (dx,x,v,eS_s_v,eS_d,v,p,t) 
    dim::Int # number of dimensions of x
    mass_matrix=I # Mass matrix for the equation
    sym=[:v for i in 1:dim] # Symbols for the dimensions
end

#TODO: Pick up here
"""
    ODEEdge(f!, dim, mass_matrix, sym)

Wrapper that ensures compatibility of a **mutating** function **`f!`** with
the key constructor `network_dynamics`.

**`f!`**  describes the local behaviour at a dynamic edge and has to respect
the following calling syntax

```julia
f!(de, e, v_s, v_t, p, t) -> nothing
```

Here  `de`, `e`, `p` and `t` are the usual arguments, while
`v_s` and `v_d` are arrays containing the vertices which are
the source and destination of the described edge.

**`dim`** is the number of independent variables in the edge equations and
**`sym`** is an array of symbols for these variables. For more details see
the documentation.
**`mass_matrix`** is an optional argument that defaults to the identity
matrix `I`. If a mass matrix M is given the system `M * de = f!` will be
solved.

For more details see the documentation.
"""
# The ODEEdge contructor creates an ODEEAEdge
function ODEEdge{T}(f!::T,dim::Int;mass_matrix=I,sym=[:e for i in 1:ed_dim])
    fea! = (ded,ed,eS,v_s,v_d,p,t) -> f!(ded,ed,v_s,v_d,p,t)
    return ODEEAEdge(fea!,
                     ed_dim=dim,
                     eS_dim=0,
                     mass_matrix=mass_matrix,
                     ed_sym=sym)
end

"""
    ODEEAEdge(f!, dyn_dim, ex_dim, mass_matrix, dyn_sym, ex_sym)

This is an ODE edge that also contains an explicit algebraic equation. 

Inputs are 
    `f!` rhs function (described below)
    `dyn_dim` dynamic dimension of the edge
    `ex_dim` explicit algebraic dimension of the edge
    `mass_matrix` mass matrix elements for the dynamic portion
    `dyn_sym` symbols associated with `dyn_dim`s
    `ex_sym` symbols associated with  `ex_dim`s. 

**f!** describes the dynamics behavior of the function,  both with explicit algebraic and differential states. 

```julia
f!(ded,ed,eS,v_s,v_d,p,t) -> nothing
```
with `ded` edge' differential array (output)
     `ed` edge differential array (input)
     `eS` edge explicit algebraic static terms (output)
     `v_s`, `v_d` the source and destination of the edge
     `p` parameters
     `t` time

"""
@Base.kwdef struct ODEEAEdge{T} <: EdgeFunction
    f!::T # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    ed_dim::Int # number of dimensions of x
    eS_dim::Int # number of non-dynamic dimensions of x
    mass_matrix=I # Mass matrix for the equation
    ed_sym=[:e for i in 1:ed_dim] # Symbols for the dimensions
    eS_sym=[:e for i in 1:eS_dim] # Symbols for the dimensions
end



"""
    DDEVertex(f!, dim, mass_matrix, sym)

Wrapper that ensures compatibility of a **mutating** function **`f!`** with
the key constructor `network_dynamics`.

**`f!`**  describes the local behaviour at a dynamic node and has to respect
the following calling syntax

```julia
f!(dv, v, e_s, e_t, h, p, t) -> nothing
```

Here `dv`, `v`, `p` and `t` are the usual ODE arguments, while
`e_s` and `e_d` are arrays containing the edges for which the
described vertex is the source or the destination respectively. `h` is the history array for `v`.

**`dim`** is the number of independent variables in the vertex equations and
**`sym`** is an array of symbols for these variables.
**`mass_matrix`** is an optional argument that defaults to the identity
matrix `I`. If a mass matrix M is given the system `M * dv = f!` will be
solved.

For more details see the documentation.
"""
@Base.kwdef struct DDEVertex{T} <: VertexFunction
    f!::T # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    ed_dim::Int # number of dimensions of x
    eS_dim::Int # number of non-dynamic dimensions of x
    mass_matrix=I # Mass matrix for the equation
    ed_sym=[:e for i in 1:dim] # Symbols for the dimensions
    eS_sym=[:ex for i in 1:dim] # Symbols for the dimensions
end

function DDEVertex(ov::ODEVertex)
    f! = (dv, v, e_s, e_d, h_v, p, t) -> ov.f!(dv, v, e_s, e_d, p, t)
    DDEVertex(f!, ov.dim, ov.mass_matrix, ov.sym)
end

# function DDEVertex(sv::StaticVertex)
#     f! = (dv, v, e_s, e_d, h_v, p, t) -> sv.f!(v, e_s, e_d, p, t)
#     DDEVertex(f! = f!, dim = sv.dim, sym= sv.sym)
# end
function DDEVertex(sv::StaticVertex)
    f! = (dv, v, e_s, e_d, h_v, p, t) -> ODE_from_Static(sv.f!)(dv, v, e_s, e_d, p, t)
    DDEVertex(f!, sv.dim, 0., sv.sym)
end

# Promotion rules [eventually there might be too many types to hardcode everyhting]

convert(::Type{DDEVertex}, x::StaticVertex) = DDEVertex(x)
promote_rule(::Type{DDEVertex}, ::Type{StaticVertex}) = DDEVertex


convert(::Type{DDEVertex}, x::ODEVertex) = DDEVertex(x)
promote_rule(::Type{DDEVertex}, ::Type{ODEVertex}) = DDEVertex


"""
Like a static edge but with extra arguments for the history of the source and destination vertices. This is NOT a DDEEdge.
"""
 @Base.kwdef struct StaticDelayEdge{T} <: EdgeFunction
    f!::T # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    eS_dim::Int # number of non-dynamic dimensions of x
    eS_sym=[:ex for i in 1:dim] # Symbols for the dimensions
end

function StaticDelayEdge(se::ODEEAEdge)
    @assert(se.ed_dim=0) # for now, this can't have differential terms.
    f! = (e, v_s, v_d, h_v_s, h_v_d, p, t) -> se.f!(
                                ded=nothing,
                                ed=nothing, 
                                eS=e,
                                v_s=v_s, 
                                v_d=v_d, 
                                p=p, 
                                t=t)
    StaticDelayEdge(f!, se.dim, se.sym)
    # TODO: This approach makes a lot of small structs,
    # as each StaticEdge that gets promoted now exists as a StaticEdge
    # and a StaticDelayEdge. In a large system of mostly static edges,
    # this is not insignificant as they all get doubled. A better style
    # should be realizable.
    #
    # TODO: This could also be a delay ODE/EA edge, with ded!= nothing.
    # That resolves the assert at the start of this function.
end
# Return a ODEEAEdge

# Promotion rules

convert(::Type{StaticDelayEdge}, x::ODEEAEdge) = StaticDelayEdge(x)
promote_rule(::Type{StaticDelayEdge}, ::Type{ODEEAEdge}) = StaticDelayEdge


convert(::Type{ODEVertex}, x::StaticVertex) = ODEVertex(x)
promote_rule(::Type{ODEVertex}, ::Type{StaticVertex}) = ODEVertex

# Not sure if the next line does something?
promote_rule(::Type{ODEVertex{T}}, ::Type{ODEVertex{U}}) where {T, U} = ODEVertex


