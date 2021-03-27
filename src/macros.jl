## Build symbolic model

function symbolic_model(domain::Domain{N,T}, sys) where {N,T}
    println("symbolic_model! started")
    symb = Symbolic(domain)
    ncell = get_ncells(domain)
    graph = Graph(ncell)
    r = domain.grid.h/2
    _H_ = SMatrix{N,N}(I).*r
    _ONE_ = ones(SVector{N})
    ntrans = 0
    Fe = sys.error_map(norm(r, Inf))
    Fr = SVector{N,T}(r .+ Fe)

    for pos in enum_pos(domain)
        source = get_state_by_pos(symb, pos)
        x = get_coord_by_pos(domain.grid, pos)
        Fx, DFx = sys.linsys_map(x, _H_)
        A = inv(DFx)
        b = abs.(A)*Fr .+ 1
        HP = CenteredPolyhedron(A, b)
        # TODO: can we improve abs.(DFx)*_ONE_?
        rad = abs.(DFx)*_ONE_ .+ Fe
        rectI = get_pos_lims_outer(domain.grid, HyperRectangle(Fx - rad, Fx + rad))
        # HyperRectangle(Fx - rad, Fx + rad) can be smaller than HP. Therefore,
        # in the plots, we may have cells not in symb while the 1st-order approx
        # cover them.
        fpos_iter = Iterators.product(_ranges(rectI)...)
        for fpos in fpos_iter
            dx = get_coord_by_pos(domain.grid, fpos) - Fx
            !(dx in HP) && continue
            if fpos in domain
                target = get_state_by_pos(symb, fpos)
                add_transition!(graph, source, target)
                ntrans += 1
            end
        end
    end

    println("symbolic_model! terminated with success: $(ntrans) transitions created")
    return graph, symb
end

# Essential graph

function viable_states!(statelist, graph::Graph{S}, viablelist) where S
    println("viable_states! started")
    viableset = Set(viablelist)
    npre = Dict((state, 0) for state in enum_states(graph))
    npost = Dict((state, 0) for state in enum_states(graph))
    for trans in enum_transitions(graph)
        if trans.source ∈ viableset && trans.target ∈ viableset
            npre[trans.target] += 1
            npost[trans.source] += 1
        end
    end

    transset = Set(enum_transitions(graph))
    leavingtransset = Set{Transition{S}}()
    while true
        empty!(leavingtransset)
        for trans in transset
            if npre[trans.source] < 1 || npost[trans.target] < 1
                push!(leavingtransset, trans)
                npre[trans.target] -= 1
                npost[trans.source] -= 1
            end
        end
        if isempty(leavingtransset)
            break
        end
        setdiff!(transset, leavingtransset)
    end

    nviable = 0
    for state in keys(npre)
        if npre[state] > 0 && npost[state] > 0
            push!(statelist, state)
            nviable += 1
        end
    end

    println("viable_states! terminated with success: $(nviable) viable states added")
end
