# Some higher-level macro for building and making operations on the graph related
# to the abstraction

################################################################################
# Computes the max norm_Inf_2_matp of DDs on every cell in idx_list, with nsub
# points given in theAbs.nsubopt_list.
function compute_dslipInf!(theAbs, idx_list)
    idx_list, DDs_list = _compute_dslip_(theAbs, idx_list, Inf)
    theAbs.dslipInf_list[idx_list] = DDs_list
    return DDs_list
end

function compute_dslip2!(theAbs, idx_list)
    idx_list, DDs_list = _compute_dslip_!(theAbs, idx_list, 2)
    theAbs.dslip2_list[idx_list] = DDs_list
    return DDs_list
end

function compute_dslip(theAbs, idx_list, p)
    return _compute_dslip_(theAbs, idx_list, p)[2]
end

function _compute_dslip_(theAbs, idx_list, p)
    idx_list = preprocess_idxlist(theAbs, idx_list)
    check_nsublist(theAbs.nsubopt_list[idx_list], theAbs.dim)
    nidx = length(idx_list)
    DDs_list = Vector{Float64}(undef, nidx)
    F0 = DDsNorm_function(theAbs, p)

    for (i, idx) in enumerate(idx_list)
        dslip = minimize_function_cell(theAbs, idx, F0)
        @printf("Iter %d (/%d) (idx: %d) -> dslip (p: %s): %f\n",
            i, nidx, idx, p, dslip)
        DDs_list[i] = dslip
    end

    return (idx_list, DDs_list)
end

################################################################################
function DDs_normInf2matp_rect(theAbs, p, x0, h, nsub)
    # Compute the max norm_Inf_2_matp of DDs on a rectangle given by x0 and h and
    # from different initial points given by nsub.
    @assert length(x0) == theAbs.dim
    F0 = DDsNorm_function(theAbs, p)
    res = minimum_grid(F0, x0, h, nsub, theAbs.dim, true)
    return (-res[1], res[2])
end

################################################################################
function compute_nsubim!(theAbs, idx_list, relmargin, absmargin)
    # Computes the nusbim_list for all indexes in idx_list, such that 0.5*dslipInf*httm^2
    # (where httm = max(h.*Lev./nsub)) is smaller than min(absmargin, relmargin*min(h)).
    # Requires dslipInf to be >= 0.0.
    idx_list = preprocess_idxlist(theAbs, idx_list)
    @assert all(theAbs.dslipInf_list[idx_list] .>= 0.0)
    nidx = length(idx_list)
    hmin = minimum(theAbs.h)
    margin = min(absmargin, relmargin*hmin)
    @assert margin > 0.0
    tol = 1e-5

    for (i, idx) in enumerate(idx_list)
        theCell = theAbs.cell_list[idx]
        ht = theCell.Lev.*theAbs.h
        dslipInf = theAbs.dslipInf_list[idx]
        coeff = sqrt(0.5*dslipInf/margin)
        nsubim = ceil.(Int, coeff*ht .+ tol)
        theAbs.nsubim_list[idx] = nsubim
        @printf("iter %d (/%d) (idx: %s) -> dslipInf: %f, nsub: %s\n",
            i, nidx, idx, dslipInf, nsubim)
    end

    return deepcopy(theAbs.nsubim_list[idx_list])
end

################################################################################
# Extends the set of cells starting from the cells in idx_list. Also returns the
# list of indexes of all new cells, and a map (true or false) if the cell has been
# visited (i.e., its edges computed). The nsubim of new cells are computed
# accordingly to relmargin and absmargin.
#
# Lazy -> dslipInf of new cells is given by dslipInfnew.
# Sharp -> dslipInf of new cells is computed using DDx_normInf2matp_cell and with
# nsubopt of new cells given by nsuboptnew.

function DXLIP_lazy!(theAbs, idx_list, dslipInfnew)
    theAbs.dslipInf_list[idx_list] = fill(dslipInfnew, length(idx_list))
end

function extend_invariant_lazy!(theAbs, idx_list, maxIter, dslipInfnew, relmargin, absmargin)
    @inline DXLIP!(theAbs, idx_list) = DXLIP_lazy!(theAbs, idx_list, dslipInfnew)
    return _extend_invariant_!(theAbs, idx_list, maxIter, DXLIP!,
        relmargin, absmargin, fill(-1, theAbs.dim))
end

function extend_invariant_sharp!(theAbs, idx_list, maxIter, nsuboptnew, relmargin, absmargin)
    @inline DXLIP!(theAbs, idx_list) = compute_dslipInf!(theAbs, idx_list)
    return _extend_invariant_!(theAbs, idx_list, maxIter, DXLIP!,
        relmargin, absmargin, nsuboptnew)
end

function _extend_cells_!(theAbs, idx_list, batch)
    ncellold = length(theAbs.cell_list)
    ncell = ncellold
    @inline ONE_() = fill(1, theAbs.dim) # "copy"

    for (i, idx) in enumerate(idx_list)
        idxim_list, coordsim_list = indexlist_imagecell(theAbs, idx)
        for dest_idx in idxim_list
            push!(theAbs.edge_list, [idx, dest_idx])
        end
        for coords in coordsim_list
            ncell += 1
            push!(theAbs.edge_list, [idx, ncell])
            theCell = CellType(coords, ONE_())
            push!(theAbs.cell_list, theCell)
            @printf("Batch: %d, Iter: %d (/%d) (idx: %d) -> New cell: %d: %s\n",
                batch, i, length(idx_list), idx, ncell, coords)
        end
    end

    idxnew_list = Vector{Int}(ncellold+1:ncell)
    nidxnew = length(idxnew_list)
    append!(theAbs.nsubim_list, fill(-ONE_(), nidxnew))
    append!(theAbs.dslipInf_list, fill(-1.0, nidxnew))
    append!(theAbs.nsubopt_list, fill(-ONE_(), nidxnew))

    return idxnew_list
end

function _extend_invariant_!(theAbs, idx_list, maxIter, DXLIP!,
        relmargin, absmargin, nsuboptnew)
    #---------------------------------------------------------------------------
    # Core function of extend_invariant!
    idx_list = preprocess_idxlist(theAbs, idx_list)
    @assert all(theAbs.dslipInf_list[idx_list] .>= 0.0)
    check_nsublist(theAbs.nsubim_list[idx_list], theAbs.dim)
    @assert length(nsuboptnew) == theAbs.dim

    batch = 0
    nidx_extended = length(idx_list)
    ncellold = length(theAbs.cell_list)
    visited_bool = falses(ncellold)

    while nidx_extended <= maxIter && ~isempty(idx_list)
        batch += 1
        idxnew_list = _extend_cells_!(theAbs, idx_list, batch)
        nidxnew = length(idxnew_list)
        theAbs.nsubopt_list[idxnew_list] = fill(nsuboptnew, nidxnew)
        DXLIP!(theAbs, idxnew_list)
        compute_nsubim!(theAbs, idxnew_list, relmargin, absmargin)
        append!(visited_bool, falses(nidxnew))
        visited_bool[idx_list] = trues(length(idx_list))
        idx_list = idxnew_list
        nidx_extended += nidxnew
    end

    if ~isempty(idx_list)
        @printf("Warning: extend_invariant! has reached maxIter (%d)\n", maxIter)
    end

    ncell = length(theAbs.cell_list)
    idxnew_list = Vector{Int}(ncellold+1:ncell)
    nidxnew = length(idxnew_list)
    append!(theAbs.DS_list, fill(Matrix{Float64}[], nidxnew))
    append!(theAbs.dslip2_list, fill(-1.0, nidxnew))
    append!(theAbs.nsubds_list, fill(fill(-1, theAbs.dim), nidxnew))
    @printf("extend_invariant! finished -> new number of cells: %d (old: %d)\n",
        nidxnew, ncellold)

    return (idxnew_list, visited_bool)
end

################################################################################
function remove_cells!(theAbs, idx_list)
    # Removes the cells with index in idx_list, updates the edge_list, nsubim_list,
    # dslipInf_list and nsubopt_list accordingly.
    # Also returns an array for the mapping between the old and new indexes (a "-1"
    # means that the cell has been removed).
    ncell = length(theAbs.cell_list)
    if isempty(idx_list)
        return Vector{Int}(1:ncell)
    end
    @assert all(idx_list .> 0) && all(idx_list .<= ncell)
    idx_list = unique(idx_list)
    nidx = length(idx_list)
    idx_bool = trues(ncell)
    idx_bool[idx_list] = falses(nidx)
    newidx_list = cumsum(idx_bool)
    theAbs.cell_list = theAbs.cell_list[idx_bool]
    theAbs.nsubim_list = theAbs.nsubim_list[idx_bool]
    theAbs.nsubopt_list = theAbs.nsubopt_list[idx_bool]
    theAbs.dslipInf_list = theAbs.dslipInf_list[idx_bool]
    theAbs.DS_list = theAbs.DS_list[idx_bool]
    theAbs.dslip2_list = theAbs.dslip2_list[idx_bool]
    theAbs.nsubds_list = theAbs.nsubds_list[idx_bool]
    idx_map = fill(-1, ncell)
    idx_map[idx_bool] = newidx_list[idx_bool]

    edge_list = theAbs.edge_list
    nedge0 = length(edge_list)
    nedge = 0

    for i = 1:nedge0
        if idx_bool[edge_list[i][1]] && idx_bool[edge_list[i][2]]
            nedge += 1
            edge_list[nedge][1] = newidx_list[edge_list[i][1]]
            edge_list[nedge][2] = newidx_list[edge_list[i][2]]
        end
    end

    resize!(edge_list, nedge)
    @printf("remove_cells! finished -> removed %d cells (/%d),\n",
        nidx, ncell)
    @printf("                       -> number of remaining edges %d (/%d)\n",
        nedge, nedge0)

    return idx_map
end

################################################################################
@inline function onestepblocking!(edge_list, nonblocking_list, out_list, in_list)
    # Puts a true in out_list[idx] iff there is an edge idx->i* to a cell i* that
    # is non-blocking (i.e., blocking_list[i*] = true).
    for i = 1:length(edge_list)
        if nonblocking_list[edge_list[i][2]]
            out_list[edge_list[i][1]] = true
        end
        if nonblocking_list[edge_list[i][1]]
            in_list[edge_list[i][2]] = true
        end
    end
end

function find_blockingcells(theAbs)
    # Finds the index of all cells that are not in a nontrivial cycle.
    ncell = length(theAbs.cell_list)
    # All cells that are not considered as blocking at step k.
    # Initially, they are all non-blocking.
    nonblocking_list = trues(ncell)
    # Used to identify the cells that have an edge FROM (in) or TO (out)
    # a non-blocking cell:
    out_list = falses(ncell)
    in_list = falses(ncell)
    countnew = ncell
    countprev = countnew + 1

    while countnew < countprev
        countprev = countnew
        onestepblocking!(theAbs.edge_list, nonblocking_list, out_list, in_list)
        # The new non-blocking cells are the ones with an edge FROM and TO a previous
        # non-blocking cell.
        nonblocking_list = out_list .& in_list
        countnew = sum(nonblocking_list)
        fill!(out_list, false)
        fill!(in_list, false)
    end

    idx_list = findall(.~nonblocking_list)

    @printf("find_blockingcells! finished -> number of non-blocking cells %d (/%d)\n",
        sum(nonblocking_list), ncell)

    return idx_list
end

################################################################################
function compute_edges!(theAbs, idx_list)
    # Builds the edges starting from the cells with index in idx_list.
    idx_list = preprocess_idxlist(theAbs, idx_list)
    @assert all(theAbs.dslipInf_list[idx_list] .>= 0.0)
    check_nsublist(theAbs.nsubim_list[idx_list], theAbs.dim)
    nedge = 0

    for (i, idx) in enumerate(idx_list)
        idxim_list, = indexlist_imagecell(theAbs, idx)
        nedge += length(idxim_list)
        for dest_idx in idxim_list
            push!(theAbs.edge_list, [idx, dest_idx])
        end
    end

    @printf("compute_edges! finished -> %d edges computed\n", nedge)
end

################################################################################
function refine_abstraction!(theAbs, idx_list, coeffs)
    idx_list = preprocess_idxlist(theAbs, idx_list)
    check_nsublist([coeffs], theAbs.dim)
    theAbs.h = theAbs.h./coeffs
    cell_list = theAbs.cell_list
    ncellold = length(cell_list)
    resize!(theAbs.edge_list, 0)

    for i = 1:ncellold
        theCell = cell_list[i]
        cell_list[i] = CellType(theCell.coords.*coeffs, theCell.Lev.*coeffs)
    end

    subs = [0:coeffs[i]-1 for i = 1:theAbs.dim]
    sub_iter_full = Iterators.product(subs...)
    sub_iter = Iterators.drop(sub_iter_full, 1)
    nnsub = length(sub_iter)
    nidx = length(idx_list)

    for arr in [cell_list, theAbs.nsubim_list, theAbs.dslipInf_list,
            theAbs.nsubopt_list, theAbs.dslip2_list, theAbs.nsubds_list]
        resize!(arr, ncellold + nidx*(nnsub))
    end

    idx_map = Vector{Int}(undef, ncellold + nidx*(nnsub))
    idx_map[1:ncellold] = 1:ncellold
    run = ncellold

    for idx in idx_list
        theCell = cell_list[idx]
        coords = theCell.coords
        Lev = theCell.Lev./coeffs
        cell_list[idx] = CellType(coords, Lev)

        for it in sub_iter
            run += 1
            coordst = coords .+ it
            cell_list[run] = CellType(coordst, Lev)
            idx_map[run] = idx
            theAbs.nsubim_list[run] = theAbs.nsubim_list[idx]
            theAbs.dslipInf_list[run] = theAbs.dslipInf_list[idx]
            theAbs.nsubopt_list[run] = theAbs.nsubopt_list[idx]
            theAbs.dslip2_list[run] = theAbs.dslip2_list[idx]
            theAbs.nsubds_list[run] = theAbs.nsubds_list[idx]
        end
    end

    ncell = length(cell_list)
    theAbs.DS_list = fill(Matrix{Float64}[], ncell)

    @printf("refine_abstraction! finished -> number cells: %d (old: %d)\n",
        ncell, ncellold)

    return (Vector{Int}(ncellold+1:ncell), idx_map)
end
