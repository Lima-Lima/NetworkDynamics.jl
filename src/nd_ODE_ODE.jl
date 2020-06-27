module nd_ODE_ODE_mod

using ..NetworkStructures
using ..NDFunctions
using ..Utilities
export nd_ODE_ODE

#= The arguments of the vertex functions must be of the form (dv,v,e_s,e_d,p,t),
where dv is the vertex variable derivative, v the vertex variable and e_s and e_d Arrays of edge variables that
have the vertex as source and destination respectively. p and t are as usual.
The arguments of the edge functions must be of the form (de,e,v_s,v_d,p,t),
where de is the derivative of the edge variable, e the edge variable, v_s and v_d the vertex variables of the vertices
the edge has as source and destination respectively. This works for arbitrary dimensional vertex and edge functions, they only
need to fit, i.e. don't do something like edges! = v_s - v_d when v_s and v_d have not the same dimension. =#


# In order to match the type, we need to pass both, a view that matches the type
# to be constructed, and the original array we want to construct a GD on top of.
@inline function prep_gd(dy::T, y::T, x, gd::GraphData{T, T}, gs) where T
    # println("Type match")
    gd.v_array = view(x, 1:gs.dim_v)
    gd.e_array = view(x, gs.dim_v+1:gs.dim_v+gs.dim_e)
    gd
end

@inline function prep_gd(dy, y, x, gd, gs)
    # println("Type mismatch")
    v_array = view(x, 1:gs.dim_v)
    e_array = view(x, gs.dim_v+1:gs.dim_v+gs.dim_e)
    gd 
    #GraphData(v_array, e_array, gs)
end


@Base.kwdef struct nd_ODE_ODE{G, T, T1, T2}
    vertices!::T1
    edges!::T2
    graph::G
    graph_structure::GraphStruct
    graph_data::GraphData{T, T}
    parallel::Bool
end

function (d::nd_ODE_ODE)(dx, x, p, t)
    gd = prep_gd(view(dx, 1:2), view(x, 1:2), x, d.graph_data, d.graph_structure)
    @nd_threads d.parallel for i in 1:d.graph_structure.num_e
            d.graph_data.expr_array[i] = maybe_idx(d.edges!, i).f!(
                view(dx,d.graph_structure.e_idx[i] .+ d.graph_structure.dim_v), 
                gd.e[i],
                gd.v_s_e[i], gd.v_d_e[i],
                p_e_idx(p, i),
                t)
            # I want edge.f! to update the state variable in [x], just as it
            # does now. But then it should also update a expression array.  So,
            # there will be [v1,v2,v3,v4,e1,e2]::float as state vector x, and
            # [expr1,expr2]::float as the edge command vector.  The edge command
            # vector is then drawn directly into the vertex aggregation.
            #
            # So let's keep edge.f! mutating the array for the differential
            # variables, but it can return the result of the explicit edge
            # function. 
            #
            # If a system was something like:
            #   \dot Y = P - \sum K_{ij} * sin(y_d + y_s) 
            #   \dot K = K - (Yi+Yj) 
            #
            #   Then the edge functions would look like
            #function f(.....)
            #   de[1] = e[1] - (V_d[1] + V_s[1]}
            #   return e[1] * V_d[1] + V_s[1]
            #end
            # 
            # A nicer alternative might be to pass an optional expr_array to the
            # rhs functions.
            #
            # The fundamental thing is that we're trying to express two
            # different ideas in the edge function: how states on the edges
            # evolve over time AND what expression is behind the sum in the node
            # equation. They shouldn't both get the same treatment.
            #
    end

    @nd_threads d.parallel for i in 1:d.graph_structure.num_v
            maybe_idx(d.vertices!,i).f!(view(dx,d.graph_structure.v_idx[i]),
                                        gd.v[i], 
                                        gd.expr_s_v[i], 
                                        gd.expr_d_v[i],
                                        p_v_idx(p, i), 
                                        t)
            # This loop now takes expr_s_v, which are views into the underlying
            # gd.expr_array at the correct location for source and destination
            # edges.
        end
    
    nothing
end


function (d::nd_ODE_ODE)(x, p, t, ::Type{GetGD})
    prep_gd(view(x, 1:2), view(x, 1:2), x, d.graph_data, d.graph_structure)
end

function (d::nd_ODE_ODE)(::Type{GetGS})
    d.graph_structure
end


end #module
