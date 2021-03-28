## Build symbolic model

function symbolic_model(domain::Domain{N,T}, sys) where {N,T}
    println("symbolic_model! started")
    idxn = Indexing(domain)
    graph = Graph(1:get_ncells(domain))
    r = domain.grid.h/2
    _H_ = SMatrix{N,N}(I).*r
    _ONE_ = ones(SVector{N})
    nedges = 0
    Fe = sys.error_map(norm(r, Inf))
    Fr = SVector{N,T}(r .+ Fe)

    for pos in enum_pos(domain)
        source = get_index_by_pos(idxn, pos)
        x = get_coord_by_pos(domain.grid, pos)
        Fx, DFx = sys.linsys_map(x, _H_)
        A = inv(DFx)
        b = abs.(A)*Fr .+ 1
        HP = CenteredPolyhedron(A, b)
        # TODO: can we improve abs.(DFx)*_ONE_?
        rad = abs.(DFx)*_ONE_ .+ Fe
        rng = get_pos_lims_outer(domain.grid, HyperRectangle(Fx - rad, Fx + rad))
        # HyperRectangle(Fx - rad, Fx + rad) can be smaller than HP. Therefore,
        # in the plots, we may have cells not in idxn while the 1st-order approx
        # cover them.
        for fpos in enum_elems(rng)
            dx = get_coord_by_pos(domain.grid, fpos) - Fx
            !(dx in HP) && continue
            if fpos in domain
                target = get_index_by_pos(idxn, fpos)
                add_edge!(graph, source, target)
                nedges += 1
            end
        end
    end

    println("symbolic_model! terminated with success: $(nedges) edges created")
    return graph, idxn
end

# Essential graph

function viable_states!(statelist, graph::Graph{S}, viablelist) where S
    println("viable_states! started")
    viableset = Set(viablelist)
    npre = Dict((state, 0) for state in enum_states(graph))
    npost = Dict((state, 0) for state in enum_states(graph))
    for edge in enum_edges(graph)
        if edge.source ∈ viableset && edge.target ∈ viableset
            npre[edge.target] += 1
            npost[edge.source] += 1
        end
    end

    valid_edgeset = Set(enum_edges(graph))
    wrong_edgeset = Set{Edge{S}}()
    while true
        empty!(wrong_edgeset)
        for edge in valid_edgeset
            if npre[edge.source] < 1 || npost[edge.target] < 1
                push!(wrong_edgeset, edge)
                npre[edge.target] -= 1
                npost[edge.source] -= 1
            end
        end
        if isempty(wrong_edgeset)
            break
        end
        setdiff!(valid_edgeset, wrong_edgeset)
    end

    nstates = 0
    for state in keys(npre)
        if npre[state] > 0 && npost[state] > 0
            push!(statelist, state)
            nstates += 1
        end
    end

    println("viable_states! terminated with success: $(nstates) viable states added")
end
