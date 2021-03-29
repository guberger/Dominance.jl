## Build containing the Jacobian at the center of the cells

function matrix_field(domain::Domain{N,T}, sys, idxn, statelist) where {N,T}
    A_field = Dict{Int,SMatrix{N,N,T}}()
    r = domain.grid.h/2
    _H_ = SMatrix{N,N}(I)
    for state in statelist
        pos = get_elem_by_index(idxn, state)
        x = get_coord_by_pos(domain.grid, pos)
        ~, DFx = sys.linsys_map(x, _H_)
        A_field[state] = DFx
    end
    return A_field
end
