# Compute the images of the cells and whose cell they intersect

#=
In computing, inline expansion, or inlining, is a manual or compiler optimization
that replaces a function call site with the body of the called function. Inline
expansion is similar to macro expansion, but occurs during compilation, without
changing the source code (the text), while macro expansion occurs prior to compilation,
and results in different text that is then processed by the compiler.
=#

################################################################################
@inline function cross_product_1d(A, dim)
    return [1.0]
end

@inline function cross_product_2d(A, dim)
    return [A[2, 1], -A[1, 1]]
end

@inline function cross_product_3d(A, dim)
    return [A[2, 1]*A[3, 2] - A[3, 1]*A[2, 2],
        -(A[1, 1]*A[3, 2] - A[3, 1]*A[1, 2]),
        A[1, 1]*A[2, 2] - A[2, 1]*A[1, 2]]
end

@inline function cross_product_nd(A, dim)
    C = Vector{Float64}(undef, dim)
    idx = fill(true, dim)
    rho = 1
    for i = 1:dim
        idx[i] = false
        C[i] = rho*det(A[idx, :])
        rho = -rho
        idx[i] = true
    end
    return C
end

function augmented_linearbox!(A_lineq, A1, A2, colidx_iter, vert_iter, dim, reltol, abstol)
    # Computes the sum of A1*B0 + A2*B0 where B0 is the hypercube [-1,1]^dim
    # and returns a matrix A such that the sum is given by {x : |A*x| <= 1}.
    A = hcat(A1, A2)
    if dim == 1
        cross_prod = cross_product_1d
    elseif dim == 2
        cross_prod = cross_product_2d
    elseif dim == 3
        cross_prod = cross_product_3d
    else
        cross_prod = cross_product_nd
    end
    nA = [norm(A[:, i]) for i = 1:2*dim]
    # colidx_iter is an iterator over all (dim-1)-subsets of {1,...,2*dim}.
    # vert_iter is an iterator on the vertices of [-1,1]^dim.
    for (i, idx) in enumerate(colidx_iter)
        a = cross_prod(A[:, idx], dim)
        # Test whether a is nondegenerate.
        if norm(a) <= max(abstol, reltol*prod(nA[idx]))
            continue
        end
        normalize!(a)
        # dot(Array, Tuple) is well-defined and outputs an array.
        aA1 = a'*A1
        b1 = maximum(t -> dot(aA1, t), vert_iter)
        aA2 = a'*A2
        b2 = maximum(t -> dot(aA2, t), vert_iter)
        A_lineq[i, :] = a/(b1 + b2)
    end
end

#=
Old version of augmented_linearbox (may be faster for high dimensional systems)
=#
function augmented_linearbox_old!(A_lineq, A1, A2, colidx_iter, vert_iter, dim)
    # Computes the sum of A1*B0 + A2*B0 where B0 is the hypercube [-1,1]^dim
    # and returns a matrix A such that the sum is given by {x : |A*x| <= 1}.
    AA = hcat(A1, A2)
    # colidx_iter is an iterator over all (dim-1)-subsets of {1,...,2*dim}
    # vert_iter is an iterator on the vertices of [-1,1]^dim
    for (i, idx) in enumerate(colidx_iter)
        B = AA[:, idx]
        # B is a matrix and it must be, otherwise nullspace throws an error...
        a = nullspace(B')'
        if size(a, 1) == 1
            aA1 = a*A1
            # dot(Array, Tuple) is well-defined and outputs an array
            f1 = t -> dot(aA1, t)
            b1 = maximum(f1, vert_iter)
            aA2 = a*A2
            f2 = t -> dot(aA2, t)
            b2 = maximum(f2, vert_iter)
            A_lineq[i, :] = a/(b1 + b2)
        end
    end
end

@inline function coordslist_convexset!(coords_list, A, xc, x_or, h, coords_iter)
    # Given an iterator of coords, keeps only the coords in coords_iter that are
    # in {x+xc : |A*x| <= 1}, and stores them in coords_list.
    for coords in coords_iter
        cc = collect(coords)
        x = x_or + h.*cc - xc
        v = abs.(A*x)
        if all(p -> (p <= 1.0), v)
            push!(coords_list, cc)
        end
    end
end

function coordslist_imagecell(theAbs, idx)
    # Gives the lists of coords such that the image of cell_list[idx] intersects
    # the cell with coords and lev = [1,...,1].
    # The image is computed by dividing the cell in prod(nsub) subrectangles and
    # by computing the first-order approximation of each subrectangle.
    dim = theAbs.dim
    reltol = theAbs.set_reltol
    abstol = theAbs.set_abstol
    x_or = theAbs.x_or
    h = theAbs.h
    nsub = theAbs.nsubim_list[idx]
    dslipInf = theAbs.dslipInf_list[idx]
    # "t" refers to the cell cell_list[idx].
    coordst = theAbs.cell_list[idx].coords
    Levt = theAbs.cell_list[idx].Lev
    ht = Levt.*h
    xt = x_or + coordst.*h
    # Center of the cell
    xtc = xt + 0.5*ht
    Fxtc, DFxtc = theAbs.linsys_map(xtc)
    # The cell can be enclosed into a box [-htm,htm]^dim.
    htm = 0.5*maximum(ht)
    # The 0.5 (in 0.5*dslipInf) comes from the integration int_0^1 t dt = 0.5.
    rad = (opnorm(DFxtc, Inf)*htm + 0.5*dslipInf*htm^2)*1.05
    cmin = floor.(Int, (Fxtc .- rad - x_or)./h)
    cmax = floor.(Int, (Fxtc .+ rad - x_or)./h)

    # "tt" refers to the subrectangles.
    htt = ht./nsub
    # Each subrectangle can be enclosed into a box [-httm,httm]^dim.
    httm = 0.5*maximum(htt)
    # A cell with margin M = 0.5*dslipInf*httm^2 is given by prod_{i=1}^dim [-M,h[i]+M]
    # and thus is the translated (by h/2) image of [-1,1]^dim by A2 below.
    A2 = 0.5*(h .+ dslipInf*httm^2).*Matrix{Float64}(I, dim, dim)

    subs = [1:nsub[i] for i = 1:dim]
    sub_iter = Iterators.product(subs...)
    single_coords = [cmin[i]:cmax[i] for i = 1:dim]
    coords_iter = Iterators.product(single_coords...)
    coords_list = Vector{Int}[]
    sizehint!(coords_list, length(sub_iter)*length(coords_iter))
    colidx_iter = combinations(1:2*dim, dim-1)
    vert_iter = Iterators.product(fill(-1:2:1, dim)...)
    A_lineq = zeros(length(colidx_iter), dim)

    for it in sub_iter
        # Center of current subrectangle
        # (Tuple .* Array is well-defined and outputs an array)
        xttc = xt + it.*htt - 0.5*htt
        Sxttc, DSxttc = theAbs.linsys_map(xttc)
        A1 = 0.5*DSxttc.*htt'
        fill!(A_lineq, 0.0)
        augmented_linearbox!(A_lineq, A1, A2, colidx_iter, vert_iter, dim, reltol, abstol)
        # Translated image of [-1,1]^dim by A2 (see above)
        xc = Sxttc - 0.5*h
        coordslist_convexset!(coords_list, A_lineq, xc, x_or, h, coords_iter)
    end

    return unique(coords_list)
end

@inline function coords_isin_cell(theCell, coords)
    # Returns true iff the standard cell (i.e., with Lev = [1,...,1]) and with
    # coordinates coords is in theCell.
    return all(coords .>= theCell.coords) && all(coords .< theCell.coords + theCell.Lev)
end

function indexmapping_coordslist(cell_list, coords_list)
    # For each coords in coords_list, gives the index of the cell in cell_list
    # to which it belongs, or returns -1 if there is no such cell.
    idx_mapping = fill(-1, length(coords_list))
    for (idx, theCell) in enumerate(cell_list)
        for (i, coords) in enumerate(coords_list)
            if coords_isin_cell(theCell, coords)
                idx_mapping[i] = idx
            end
        end
    end
    return idx_mapping
end

function indexlist_imagecell(theAbs, idx)
    # Returns the indexes of the cells that intersect the image of cell_list[idx]
    # and returns also the coords of the standard cells (i.e., with Lev = [1,...,1])
    # that intersect the image of cell_list[idx] but do not belong to any cell
    # in cell_list.
    coords_list = coordslist_imagecell(theAbs, idx)
    idx_mapping = indexmapping_coordslist(theAbs.cell_list, coords_list)
    return (filter(x -> (x >= 0), idx_mapping), coords_list[idx_mapping .< 0])
end
