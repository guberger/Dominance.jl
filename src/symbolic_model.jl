## Build symbolic model

function symbolic_model(domain::Domain{N,T}, sys) where {N,T}
    println("symbolic_model! started")
    idxn = Indexing(enum_pos(domain))
    graph = Graph(get_ncells(domain))
    r = domain.grid.h/2
    _H_ = SMatrix{N,N}(I).*r
    _ONE_ = ones(SVector{N})
    nedges = 0
    Fe = sys.error_map(norm(r, Inf))
    Fr = SVector{N,T}(r .+ Fe)

    for pos in enum_pos(domain)
        source = get_index_by_elem(idxn, pos)
        x = get_coord_by_pos(domain.grid, pos)
        Fx, DFx = sys.linsys_map(x, _H_)
        A = inv(DFx)
        b = abs.(A)*Fr .+ 1
        HP = CenteredPolyhedron(A, b)
        # TODO: can we improve abs.(DFx)*_ONE_?
        rad = abs.(DFx)*_ONE_ .+ Fe
        fpos_iter = get_pos_lims_outer(domain.grid, HyperRectangle(Fx - rad, Fx + rad))
        # HyperRectangle(Fx - rad, Fx + rad) can be smaller than HP. Therefore,
        # in the plots, we may have cells not in idxn while the 1st-order approx
        # cover them.
        for fpos in fpos_iter
            dx = get_coord_by_pos(domain.grid, fpos) - Fx
            !(dx in HP) && continue
            if fpos in domain
                target = get_index_by_elem(idxn, fpos)
                add_edge!(graph, source, target)
                nedges += 1
            end
        end
    end

    println("symbolic_model! terminated: $(nedges) edges created")
    return graph, idxn
end
