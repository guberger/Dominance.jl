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

## Symbolic model

function indexing(::Type{T}, elemlist) where T
    elemset = Set(elemlist)
    ind2elem = Vector{T}(undef, length(elemset))
    elem2ind = Dict{T,Int}()
    sizehint!(elem2ind, length(elemset))
    for (i, elem) in enumerate(elemset)
        ind2elem[i] = elem
        elem2ind[elem] = i
    end
    return elem2ind, ind2elem
end

struct SymbolicModel{N,T}
    grid::Grid{N,T}
    graph::Graph
    pos2ind::Dict{NTuple{N,Int},Int}
    ind2pos::Vector{NTuple{N,Int}}
end

function get_state_by_pos(symmod, pos)
    return get(symmod.pos2ind, pos, 0)
end

function get_pos_by_state(symmod, state)
    return symmod.ind2pos[state]
end

function get_nstates(symmod::SymbolicModel)
    return get_nstates(symmod.graph)
end

## Macros

# Create a domain from a symbolic model and a list of states
function support_domain(symmod, statelist)
    domain = Domain(symmod.grid)
    for state in statelist
        pos = get_pos_by_state(symmod, state)
        add_pos!(domain, pos)
    end
    return domain
end

# Refine discretizetaion
function refine_domain(domain, nsub)
    println("refine_domain started")
    orig = domain.grid.orig
    h = domain.grid.h
    h2 = h./nsub
    orig2 = orig - h/2 + h2/2
    grid2 = Grid(orig2, h2)
    domain2 = Domain(grid2)
    sub_iter = hyper_range(nsub.*0, nsub .- 1)
    for pos in enum_pos(domain)
        for sub in sub_iter
            add_pos!(domain2, pos.*nsub .+ sub)
        end
    end
    ncells = get_ncells(domain2)
    println("refine_domain terminated: $(ncells) cells created")
    return domain2
end

# Remove edges and reindex

function trim_symbolic_model(symmod, statelist)
    stateset = Set(statelist)
    nstates2 = length(stateset)
    graph2 = Graph(nstates2)
    elem2ind, ind2elem = indexing(Int, stateset)
    ind2pos = [get_pos_by_state(symmod, ind2elem[i]) for i = 1:nstates2]
    pos2ind = Dict(ind2pos[i] => i for i = 1:nstates2)
    for edge in enum_edges(symmod.graph)
        if edge.source ∈ stateset && edge.target ∈ stateset
            add_edge!(graph2, elem2ind[edge.source], elem2ind[edge.target])
        end
    end
    return SymbolicModel(symmod.grid, graph2, pos2ind, ind2pos)
end
