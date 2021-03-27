# Analysis of p-dominance of the system

import JuMP
const JP = JuMP
using MosekTools

################################################################################
function DSlist_cell(theAbs, idx)
    # Computes the matrices of the linearized systems at nsubds points in the cell
    # with index given by idx.
    dim = theAbs.dim
    x_or = theAbs.x_or
    h = theAbs.h
    nsub = theAbs.nsubds_list[idx]
    # "t" refers to the cell cell_list[idx].
    coordst = theAbs.cell_list[idx].coords
    Levt = theAbs.cell_list[idx].Lev
    ht = Levt.*h
    xt = x_or + coordst.*h
    htt = ht./nsub

    subs = [1:nsub[i] for i = 1:dim]
    sub_iter = Iterators.product(subs...)
    DS_list = Vector{Matrix{Float64}}(undef, length(sub_iter))

    for (i, it) in enumerate(sub_iter)
        # Center of current subrectangle
        # (Tuple .* Array is well-defined and outputs an array)
        xttc = xt + it.*htt - 0.5*htt
        DS_list[i] = theAbs.linsys_map(xttc)[2]
    end

    return DS_list
end

function compute_DSlist!(theAbs, idx_list)
    # Computes the matrices of the linearized system for all cells with index
    # in idx_list.
    idx_list = preprocess_idxlist(theAbs, idx_list)
    check_nsublist(theAbs.nsubds_list[idx_list], theAbs.dim)
    nidx = length(idx_list)

    for idx in idx_list
        theAbs.DS_list[idx] = DSlist_cell(theAbs, idx)
    end
end

################################################################################
function compute_nsubds!(theAbs, idx_list, delta)
    # Computes the nusbds_list for all indexes in idx_list, such that 0.5*httm
    # (where httm = max(h.*Lev./nsub)) is smaller than delta.
    # Requires dslip2 to be >= 0.0.
    idx_list = preprocess_idxlist(theAbs, idx_list)
    @assert all(theAbs.dslip2_list[idx_list] .>= 0.0)
    nidx = length(idx_list)
    @assert delta > 0.0
    tol = 1e-5

    for (i, idx) in enumerate(idx_list)
        theCell = theAbs.cell_list[idx]
        ht = theCell.Lev.*theAbs.h
        dslip2 = theAbs.dslip2_list[idx]
        nsubds = ceil.(Int, 0.5*ht*dslip2/delta .+ tol)
        theAbs.nsubds_list[idx] = nsubds
        @printf("iter %d (/%d) (idx: %s) -> dslip2: %f, nsub: %s\n",
            i, nidx, idx, dslip2, nsubds)
    end

    return deepcopy(theAbs.nsubds_list[idx_list])
end

################################################################################
function compute_cones(theAbs, gamma)
    DS_list = theAbs.DS_list
    @assert all(x -> ~isempty(x), DS_list)
    edge_list = theAbs.edge_list
    ncell = length(theAbs.cell_list)
    dim = theAbs.dim
    G2 = gamma^2

    model = JP.Model(JP.with_optimizer(Mosek.Optimizer))
    var_list = [JP.@variable(model, [1:dim, 1:dim], Symmetric,
        base_name = string("P", i)) for i = 1:ncell]
    ee = JP.@variable(model, base_name = "ee")
    EYE_ = Matrix{Float64}(I, dim, dim)

    for edge in edge_list
        idx1 = edge[1] # source
        idx2 = edge[2] # sink
        DSl = DS_list[idx1]
        nmat = length(DSl)
        JP.@constraint(model, [k = 1:nmat], Symmetric(G2*var_list[idx1]
            - DSl[k]'*var_list[idx2]*DSl[k]
            - ee.*EYE_) in JP.PSDCone())
    end

    JP.@constraint(model, [i = 1:ncell], Symmetric(EYE_ - var_list[i]) in JP.PSDCone())
    JP.@constraint(model, [i = 1:ncell], Symmetric(EYE_ + var_list[i]) in JP.PSDCone())
    JP.@objective(model, Max, ee)

    display(model)
    JP.optimize!(model)
    display(JP.primal_status(model))
    display(JP.dual_status(model))
    @printf("value of ee: %f\n", JP.value(ee))

    eignegmax = -Inf
    eigposmin = Inf
    detmax = -Inf
    P_list = Vector{Symmetric{Float64,Array{Float64, 2}}}(undef, ncell)

    for (i, var) in enumerate(var_list)
        P = JP.value.(var)
        eigs = eigvals(P)
        eignegmax = max(eignegmax, minimum(eigs))
        eigposmin = min(eigposmin, maximum(eigs))
        detmax = max(detmax, det(P))
        P_list[i] = Symmetric(P)
    end

    println((eignegmax, eigposmin, detmax))
    return P_list
end
