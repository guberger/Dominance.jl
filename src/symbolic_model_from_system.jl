## Build symbolic model

function symbolic_model_from_system(domain::Domain{N,T}, sys, nsub) where {N,T}
    println("symbolic_model_from_system started")
    grid = domain.grid
    pos2ind, ind2pos = indexing(NTuple{N,Int}, enum_pos(domain))
    graph = Graph(get_ncells(domain))

    subpos_lims = ((nsub .- 1)./2)./nsub
    subpos_iter = hyper_range((-1).*subpos_lims, subpos_lims, nsub)
    hsub = domain.grid.h./nsub
    r = hsub/2
    _H_ = SMatrix{N,N}(I).*r
    _ONE_ = ones(SVector{N})
    nedges = 0
    Fe = sys.error_map(norm(r, Inf))
    Fr = SVector{N,T}(domain.grid.h/2 .+ Fe)
    edgetupleset = Set{Tuple{Int,Int}}()

    for pos in enum_pos(domain)
        empty!(edgetupleset)
        source = pos2ind[pos]
        for subpos in subpos_iter
            pos2 = pos .+ subpos
            x = get_coord_by_pos(grid, pos2)
            Fx, DFx = sys.linsys_map(x, _H_)
            A = inv(DFx)
            b = abs.(A)*Fr .+ 1
            HP = CenteredPolyhedron(A, b)
            # TODO: can we improve abs.(DFx)*_ONE_?
            rad = abs.(DFx)*_ONE_ .+ Fe
            fpos_iter = get_pos_lims_outer(grid, HyperRectangle(Fx - rad, Fx + rad))
            # HyperRectangle(Fx - rad, Fx + rad) can be smaller than HP. Therefore,
            # in the plots, we may have cells not in symmod while the 1st-order approx
            # cover them.
            for fpos in fpos_iter
                dx = get_coord_by_pos(grid, fpos) - Fx
                !(dx in HP) && continue
                if fpos in domain
                    target = pos2ind[fpos]
                    push!(edgetupleset, (source, target))
                end
            end
        end
        nedges += length(edgetupleset)
        for edgetuple in edgetupleset
            add_edge!(graph, edgetuple[1], edgetuple[2])
        end
    end

    println("symbolic_model_from_system terminated: $(nedges) edges created")
    return SymbolicModel(grid, graph, pos2ind, ind2pos)
end

## Build containing the Jacobian at the center of the cells
# TODO: add subdivision with `nsub` ???
function sensitivity_matrices(domain::Domain{N,T}, sys) where {N,T}
    A_field = Dict{NTuple{N,Int},SMatrix{N,N,T}}()
    grid = domain.grid
    _H_ = SMatrix{N,N}(I)
    for pos in enum_pos(domain)
        x = get_coord_by_pos(grid, pos)
        ~, DFx = sys.linsys_map(x, _H_)
        A_field[pos] = DFx
    end
    return A_field
end
