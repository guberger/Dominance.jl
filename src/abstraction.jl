## Grid

struct Grid{N,T}
    orig::SVector{N,T}
    h::SVector{N,T}
end

function get_pos_by_coord(grid::Grid{N}, x) where N
    return ntuple(i -> round(Int, (x[i] - grid.orig[i])/grid.h[i]), Val(N))
end

function get_coord_by_pos(grid, pos)
    return grid.orig + pos.*grid.h
end

function get_pos_lims_inner(grid::Grid{N}, rect) where N
    lbI = ntuple(i -> ceil(Int, (rect.lb[i] - grid.orig[i])/grid.h[i] + 0.5), Val(N))
    ubI = ntuple(i -> floor(Int, (rect.ub[i] - grid.orig[i])/grid.h[i] - 0.5), Val(N))
    return hyper_range(lbI, ubI)
end

function get_pos_lims_outer(grid::Grid{N}, rect) where N
    lbI = ntuple(i -> ceil(Int, (rect.lb[i] - grid.orig[i])/grid.h[i] - 0.5), Val(N))
    ubI = ntuple(i -> floor(Int, (rect.ub[i] - grid.orig[i])/grid.h[i] + 0.5), Val(N))
    return hyper_range(lbI, ubI)
end

function get_pos_lims(grid, rect, incl_mode::INCL_MODE)
    if incl_mode == INNER
        return get_pos_lims_inner(grid, rect)
    else
        return get_pos_lims_outer(grid, rect)
    end
end

## Domain

struct Domain{N,T}
    grid::Grid{N,T}
    elems::Set{NTuple{N,Int}}
end

function Domain(grid::Grid{N}) where N
    return Domain(grid, Set{NTuple{N,Int}}())
end

function add_pos!(domain, pos)
    push!(domain.elems, pos)
end

function add_coord!(domain, x)
    add_pos!(domain, get_pos_by_coord(domain.grid, x))
end

function add_set!(domain, rect::HyperRectangle, incl_mode::INCL_MODE)
    pos_iter = get_pos_lims(domain.grid, rect, incl_mode)
    for pos in pos_iter
        add_pos!(domain, pos)
    end
end

function add_subset!(domain1, domain2, rect::HyperRectangle, incl_mode::INCL_MODE)
    pos_iter = get_pos_lims(domain1.grid, rect, incl_mode)
    if length(pos_iter) < get_ncells(domain2)
        for pos in pos_iter
            if pos ∈ domain2
                add_pos!(domain1, pos)
            end
        end
    else
        for pos in enum_pos(domain2)
            if pos ∈ pos_iter
                add_pos!(domain1, pos)
            end
        end
    end
end

function remove_pos!(domain, pos)
    delete!(domain.elems, pos)
end

function remove_coord!(domain, x)
    remove_pos!(domain, get_pos_by_coord(domain.grid, x))
end

function remove_set!(domain, rect::HyperRectangle, incl_mode::INCL_MODE)
    pos_iter = get_pos_lims(domain.grid, rect, incl_mode)
    if length(pos_iter) < get_ncells(domain)
        for pos in pos_iter
            remove_pos!(domain, pos)
        end
    else
        for pos in enum_pos(domain)
            if pos ∈ pos_iter
                remove_pos!(domain, pos)
            end
        end
    end
end

function Base.union!(domain1::Domain, domain2::Domain)
    union!(domain1.elems, domain2.elems)
end

function Base.setdiff!(domain1::Domain, domain2::Domain)
    setdiff!(domain1.elems, domain2.elems)
end

function Base.empty!(domain::Domain)
    empty!(domain.elems)
end

function Base.in(pos, domain::Domain)
    return in(pos, domain.elems)
end

function Base.isempty(domain::Domain)
    return isempty(domain.elems)
end

function Base.issubset(domain1::Domain, domain2::Domain)
    return issubset(domain1.elems, domain2.elems)
end

function get_ncells(domain::Domain)
    return length(domain.elems)
end

function get_somepos(domain::Domain)
    return first(domain.elems)
end

function enum_pos(domain::Domain)
    return domain.elems
end
