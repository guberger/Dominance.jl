# All plotting functions

function set_axes(fig, theAbs)
    if theAbs.dim < 3
        ax = fig.gca()
    elseif theAbs.dim == 3
        ax = fig.gca(projection = "3d")
    else
        error("Dimension must be at most 3")
    end
    return ax
end

function vector_2_colormap(v_list, colormap_name)
    v_list = convert(Vector{Float64}, v_list)
    if isempty(v_list)
        v_list = [0.0]
    end
    vmin = minimum(v_list)
    vmax = maximum(v_list)
    norm = matplotlib.colors.Normalize(vmin = vmin, vmax = vmax)
    cmap = matplotlib.cm.get_cmap(colormap_name)
    # scal = matplotlib.cm.ScalarMappable(norm = norm, cmap = cmap)
    scal = matplotlib.cm.ScalarMappable(cmap = cmap)
    scal.set_array(v_list)
    fc_mat = scal.to_rgba(v_list)
    fc_list = [Tuple(fc_mat[i, :]) for i = 1:size(fc_mat, 1)]
    return (scal, fc_list)
end

function colorbar_from_scal(fig, scal)
    fig.colorbar(scal)
end

function preprocess_colors(ec, ea, fc_list, fa_list)
    @inline FC(c, a) =  matplotlib.colors.colorConverter.to_rgba(c, alpha = a)
    ec = FC(ec, ea)
    fca_list = map(FC, fc_list, fa_list)
    return (ec, fca_list)
end

@inline function array_to_matrix(Arr)
    return [Arr[i][j] for i = 1:length(Arr), j = 1:length(Arr[1])]
end

@inline function matrix_to_array(Mat)
    return [Mat[i, :] for i = 1:size(Mat, 1)]
end

function rectangle_vertlist(x0, h, nsub)
    # Give the contour of the rectangle with lower-left corner x0 and sides h
    # and with nsub[i] points on each side i
    n1 = nsub[1]
    n2 = nsub[2]
    x1 = range(x0[1], x0[1]+h[1], length = n1)
    x2 = range(x0[2], x0[2]+h[2], length = n2)
    x1List = vcat(x1, fill(x0[1]+h[1], n2-1), x1[n1-1:-1:1], fill(x0[1], n2-1))
    x2List = vcat(fill(x0[2], n1), x2[2:n2], fill(x0[2]+h[2], n1-1), x2[n2-1:-1:1])
    return [[x1List[i], x2List[i]] for i = 1:length(x1List)]
end

function cube_facemeshlist(x0, h, nsub)
    # Gives a mesh of each face of a cube with lower-left corner x0 and sides h
    # and with nsub[i] points on each side i
    dim = length(x0)
    xsub = [range(x0[i], x0[i]+h[i], length = nsub[i]) for i = 1:dim]
    facemesh_list = Vector{Vector{Array{Float64, dim-1}}}(undef, 2*dim)
    edgepoint_list = Vector{Vector{Vector{Float64}}}(undef, (2^(dim-1))*dim)
    fixedvars_iter = Iterators.product(fill([0, 1], dim-1)...)
    runFace = 0
    runEdge = 0
    idx0 = 1:dim

    for i = 1:dim
        idx = circshift(idx0, -i)
        idxnot = idx[1:dim-1]
        nsubnot = nsub[idxnot]
        hnot = h[idxnot]
        x0not = x0[idxnot]
        xsubnot = xsub[idxnot]
        xsubnot_iter = Iterators.product(xsubnot...)
        for j = 1:2
            runFace += 1
            theFace = Vector{Array{Float64, dim-1}}(undef, dim)
            theFace[i] = fill(x0[i] + h[i]*(j-1), nsubnot...)
            for k = 1:dim-1
                theFace[idx[k]] = map(x -> x[k], xsubnot_iter)
            end
            facemesh_list[runFace] = theFace
        end
        for fixedvars in fixedvars_iter
            runEdge += 1
            fixedx = x0not + fixedvars.*hnot
            fpoint = j -> circshift([fixedx..., xsub[i][j]], i)
            edgepoint_list[runEdge] = [fpoint(j) for j = 1:nsub[i]]
        end
    end

    return (facemesh_list, edgepoint_list)
end

################################################################################
# Trajectory
function plottrajectory!(ax, theAbs, x0, nstep;
        linecolor = "red", linewidth = 1.5,
        markercolor = "black", markersize = 5.0,
        density = 100)
    #---------------------------------------------------------------------------
    dim = theAbs.dim
    @assert dim <= 3
    @assert length(x0) == dim
    tstep = theAbs.tstep
    u0 = vcat(x0, zeros(dim*dim))
    prob = ODEProblem(theAbs.linsys_vector!, u0, (0.0, nstep*tstep))
    sol = solve(prob, Tsit5(), reltol = theAbs.ode_reltol, abstol = theAbs.ode_abstol)

    t1 = range(0.0, stop = nstep*tstep, length = nstep*density + 1)
    t2 = 0:tstep:nstep*tstep
    sol1 = sol(t1)
    sol2 = sol(t2)

    points_mat = array_to_matrix(sol1)
    points_col = [points_mat[:, i] for i = 1:dim]
    ax.plot(points_col..., linewidth = linewidth, color = linecolor)
    points_mat = array_to_matrix(sol2)
    points_col = [points_mat[:, i] for i = 1:dim]
    ax.plot(points_col..., linestyle = "None", marker = ".",
        markersize = markersize, markerfacecolor = markercolor,
        markeredgecolor = markercolor)
end

################################################################################
# Cells
function plotcells!(ax, theAbs, idx_list;
        facecolor_list = ["red"], facealpha_list = [0.5],
        edgecolor = "black", edgealpha = 1.0, edgewidth = 1.5,
        show_numbers = false)
    #---------------------------------------------------------------------------
    # Plots cells whose index is in idx_list
    if isempty(idx_list)
        return
    end
    idx_list = preprocess_idxlist(theAbs, idx_list, soft = true)
    fc_list = preprocess_anylist(idx_list, facecolor_list)
    fa_list = preprocess_anylist(idx_list, facealpha_list)
    ec, fca_list = preprocess_colors(edgecolor, edgealpha, fc_list, fa_list)

    if theAbs.dim == 2
        _plotcells2d_!(ax, theAbs, idx_list, fca_list, ec, edgewidth, show_numbers)
    elseif theAbs.dim == 3
        _plotcells3d_!(ax, theAbs, idx_list, fca_list, ec, edgewidth, show_numbers)
    else
        error("Dimension must be either 2 or 3")
    end
end

function _plotcells2d_!(ax, theAbs, idx_list, fca_list, ec, edgewidth, show_numbers)
    x_or = theAbs.x_or
    h = theAbs.h
    nidx = length(idx_list)
    vertlist_list = Vector{Vector{Vector{Float64}}}(undef, nidx)
    numcoords_list = Vector{Vector{Float64}}(undef, nidx)
    numtxt_list = Vector{String}(undef, nidx)

    for (i, idx) in enumerate(idx_list)
        theCell = theAbs.cell_list[idx]
        xt = x_or + theCell.coords.*h
        ht = theCell.Lev.*h
        vertlist_list[i] = rectangle_vertlist(xt, ht, [2, 2])
        numcoords_list[i] = xt + ht*0.5
        numtxt_list[i] = string(idx)
    end

    poly_list = matplotlib.collections.PolyCollection(vertlist_list)
    poly_list.set_facecolor(fca_list)
    poly_list.set_edgecolor(ec)
    poly_list.set_linewidth(edgewidth)
    ax.add_collection(poly_list)

    if show_numbers
        for i = 1:nidx
            ax.text(numcoords_list[i]..., numtxt_list[i],
            ha = "center", va = "center", fontsize = 10)
        end
    end
end

function _plotcells3d_!(ax, theAbs, idx_list, fca_list, ec, edgewidth, show_numbers)
    x_or = theAbs.x_or
    h = theAbs.h
    nidx = length(idx_list)
    numcoords_list = Vector{Vector{Float64}}(undef, nidx)
    numtxt_list = Vector{String}(undef, nidx)

    for (i, idx) in enumerate(idx_list)
        theCell = theAbs.cell_list[idx]
        xt = x_or + theCell.coords.*h
        ht = theCell.Lev.*h
        mesh_list, = cube_facemeshlist(xt, ht, fill(2, theAbs.dim))
        for j = 1:length(mesh_list)
            surf = ax.plot_surface(mesh_list[j]...)
            surf.set_facecolor(fca_list[i])
            surf.set_edgecolor(ec)
            surf.set_linewidth(edgewidth)
        end
        numcoords_list[i] = xt + ht*0.5
        numtxt_list[i] = string(idx)
    end

    if show_numbers
        for i = 1:nidx
            ax.text(numcoords_list[i]..., numtxt_list[i],
            ha = "center", va = "center", fontsize = 10)
        end
    end
end

################################################################################
# Images
function plotimages!(ax, theAbs, idx_list; nsub_list = [fill(5, theAbs.dim)],
        facecolor_list = ["blue"], facealpha_list = [0.5],
        edgecolor = "darkblue", edgealpha = 1.0, edgewidth = 1.5,
        wirecolor_list = facecolor_list, wirealpha_list = facealpha_list,
        wirewidth = 0.1)
    #---------------------------------------------------------------------------
    # Plots images of the cells whose index is in idx_list
    if isempty(idx_list)
        return
    end
    check_nsublist(nsub_list, theAbs.dim)
    idx_list = preprocess_idxlist(theAbs, idx_list, soft = true)
    nsub_list = preprocess_anylist(idx_list, nsub_list)
    fc_list = preprocess_anylist(idx_list, facecolor_list)
    fa_list = preprocess_anylist(idx_list, facealpha_list)
    wc_list = preprocess_anylist(idx_list, wirecolor_list)
    wa_list = preprocess_anylist(idx_list, wirealpha_list)
    ec, fca_list = preprocess_colors(edgecolor, edgealpha, fc_list, fa_list)
    void, wca_list = preprocess_colors(edgecolor, edgealpha, wc_list, wa_list)

    if theAbs.dim == 2
        _plotimages2d_!(ax, theAbs, idx_list, nsub_list, fca_list, ec, edgewidth)
    elseif theAbs.dim == 3
        _plotimages3d_!(ax, theAbs, idx_list, nsub_list, fca_list, ec, edgewidth,
            wca_list, wirewidth)
    else
        error("Dimension must be either 2 or 3")
    end
end

function _plotimages2d_!(ax, theAbs, idx_list, nsub_list, fca_list, ec, edgewidth)
    x_or = theAbs.x_or
    h = theAbs.h
    nidx = length(idx_list)
    vertlist_list = Vector{Vector{Vector{Float64}}}(undef, nidx)
    TWO_ = fill(2, theAbs.dim)

    for (i, idx) in enumerate(idx_list)
        theCell = theAbs.cell_list[idx]
        xt = x_or + theCell.coords.*h
        ht = theCell.Lev.*h
        nsub = max.(nsub_list[i] .* theCell.Lev .- theCell.Lev .+ 1, TWO_)
        vertlist = rectangle_vertlist(xt, ht, nsub)
        vertlist_list[i] = map(x -> theAbs.linsys_map(x)[1], vertlist)
    end

    poly_list = matplotlib.collections.PolyCollection(vertlist_list)
    poly_list.set_facecolor(fca_list)
    poly_list.set_edgecolor(ec)
    poly_list.set_linewidth(edgewidth)
    ax.add_collection(poly_list)
end

function _plotimages3d_!(ax, theAbs, idx_list, nsub_list, fca_list, ec, edgewidth,
        wca_list, wirewidth)
    #---------------------------------------------------------------------------
    x_or = theAbs.x_or
    h = theAbs.h
    nidx = length(idx_list)
    TWO_ = fill(2, theAbs.dim)
    fmap = (x...) -> theAbs.linsys_map([x...])[1]

    for (i, idx) in enumerate(idx_list)
        theCell = theAbs.cell_list[idx]
        xt = x_or + theCell.coords.*h
        ht = theCell.Lev.*h
        nsub = max.(nsub_list[i].*theCell.Lev .- theCell.Lev .+ 1, TWO_)
        mesh_list, edge_list = cube_facemeshlist(xt, ht, nsub)
        for j = 1:length(mesh_list)
            mesh_vec = map(fmap, mesh_list[j]...)
            mesh_sliced = [[x[i] for x in mesh_vec] for i = 1:theAbs.dim]
            surf = ax.plot_surface(mesh_sliced...)
            surf.set_facecolor(fca_list[i])
            surf.set_edgecolor(wca_list[i])
            surf.set_linewidth(wirewidth)
        end
        for j = 1:length(edge_list)
            points_vec = map(x -> theAbs.linsys_map(x)[1], edge_list[j])
            points_sliced = [[x[i] for x in points_vec] for i = 1:theAbs.dim]
            line = ax.plot(points_sliced...)
            line[1].set_color(ec)
            line[1].set_linewidth(edgewidth)
        end
    end
end

################################################################################
# Forst-order approximation
function plotapprox!(ax, theAbs, idx_list;
        edgeincolor = "black", edgeinalpha = 1.0, edgeinwidth = 0.5,
        faceoutcolor_list = ["yellow"], faceoutalpha_list = [0.5],
        edgeoutcolor_list = faceoutcolor_list, edgeoutalpha_list = faceoutalpha_list,
        edgeoutwidth = 0.5)
    #---------------------------------------------------------------------------
    # Plots first-order outer approximation of the cells whose index is in idx_list
    if isempty(idx_list)
        return
    end
    idx_list = preprocess_idxlist(theAbs, idx_list, soft = true)
    check_nsublist(theAbs.nsubim_list[idx_list], theAbs.dim)
    fc_list = preprocess_anylist(idx_list, faceoutcolor_list)
    fa_list = preprocess_anylist(idx_list, faceoutalpha_list)
    ec_list = preprocess_anylist(idx_list, edgeoutcolor_list)
    ea_list = preprocess_anylist(idx_list, edgeoutalpha_list)
    einc, fca_list = preprocess_colors(edgeincolor, edgeinalpha, fc_list, fa_list)
    void, eoutca_list = preprocess_colors(edgeincolor, edgeinalpha, ec_list, ea_list)

    if theAbs.dim == 2
        _plotapprox2d_!(ax, theAbs, idx_list, einc, edgeinwidth, fca_list,
            eoutca_list, edgeoutwidth)
    elseif theAbs.dim == 3
        _plotapprox3d_!(ax, theAbs, idx_list, einc, edgeinwidth, fca_list,
            eoutca_list, edgeoutwidth)
    else
        error("Dimension must be either 2 or 3")
    end
end

function _plotapprox2d_!(ax, theAbs, idx_list, einc, edgeinwidth, fca_list,
        eoutca_list, edgeoutwidth)
    #---------------------------------------------------------------------------
    for (i, idx) in enumerate(idx_list)
        vertinl_list, vertoutl_list = _approx_imagecell2d_!(theAbs, idx)
        polyin_list = matplotlib.collections.PolyCollection(vertinl_list)
        polyin_list.set_facecolor("none")
        polyin_list.set_edgecolor(einc)
        polyin_list.set_linewidth(edgeinwidth)
        ax.add_collection(polyin_list)
        polyout_list = matplotlib.collections.PolyCollection(vertoutl_list)
        polyout_list.set_facecolor(fca_list[i])
        polyout_list.set_edgecolor(eoutca_list[i])
        polyout_list.set_linewidth(edgeoutwidth)
        ax.add_collection(polyout_list)
    end
end

function _plotapprox3d_!(ax, theAbs, idx_list, einc, edgeinwidth, fca_list,
        eoutca_list, edgeoutwidth)
    #---------------------------------------------------------------------------
    for (i, idx) in enumerate(idx_list)
        meshin_list, vertoutl_list = _approx_imagecell3d_!(theAbs, idx)
        for j = 1:length(meshin_list)
            surf = ax.plot_surface(meshin_list[j]...)
            surf.set_facecolor(fca_list[i])
            surf.set_edgecolor(einc)
            surf.set_linewidth(edgeinwidth)
        end
        polyout_list = art3d.Poly3DCollection(vertoutl_list)
        polyout_list.set_facecolor(fca_list[i])
        polyout_list.set_edgecolor(eoutca_list[i])
        polyout_list.set_linewidth(edgeoutwidth)
        ax.add_collection3d(polyout_list)
    end
end

# Computes the first-order approximation and its outer buffer for each subcell
# the cell with index idx.
# See function "coordslist_imagecell"

function _approx_imagecell2d_!(theAbs, idx)
    dim = theAbs.dim
    x_or = theAbs.x_or
    h = theAbs.h
    nsub = theAbs.nsubim_list[idx]
    dslipInf = theAbs.dslipInf_list[idx]
    coordst = theAbs.cell_list[idx].coords
    Levt = theAbs.cell_list[idx].Lev
    ht = Levt.*h
    xt = x_or + coordst.*h

    htt = ht./nsub
    httm = 0.5*maximum(htt)
    A2 = 0.5*dslipInf*(httm^2)*Matrix{Float64}(I, dim, dim)

    subs = [1:nsub[i] for i = 1:dim]
    sub_iter = Iterators.product(subs...)
    nnsub = prod(nsub)
    vertinlist_list = Vector{Vector{Float64}}[]
    vertoutlist_list = Vector{Vector{Float64}}[]

    for it in sub_iter
        xttc = xt + it.*htt - 0.5*htt
        Sxttc, DSxttc = theAbs.linsys_map(xttc)
        A1 = 0.5*DSxttc.*htt'
        vertin_list = rectangle_vertlist(fill(-1.0, dim), fill(2.0, dim), fill(2, dim))
        push!(vertinlist_list, map(x -> A1*x + Sxttc, vertin_list))
        vertout_iter = Iterators.product(fill(-1:2:1, 2*dim)...)
        fout = x -> A1*collect(x[1:dim]) + A2*collect(x[dim+1:2*dim]) + Sxttc
        pointout_list = vec(map(fout, vertout_iter))
        hull = spatial.ConvexHull(array_to_matrix(pointout_list))
        vertices = hull.vertices .+ 1
        push!(vertoutlist_list, pointout_list[vertices])
    end

    return (vertinlist_list, vertoutlist_list)
end

function _approx_imagecell3d_!(theAbs, idx)
    dim = theAbs.dim
    x_or = theAbs.x_or
    h = theAbs.h
    nsub = theAbs.nsubim_list[idx]
    dslipInf = theAbs.dslipInf_list[idx]
    coordst = theAbs.cell_list[idx].coords
    Levt = theAbs.cell_list[idx].Lev
    ht = Levt.*h
    xt = x_or + coordst.*h

    htt = ht./nsub
    httm = 0.5*maximum(htt)
    A2 = 0.5*dslipInf*(httm^2)*Matrix{Float64}(I, dim, dim)

    subs = [1:nsub[i] for i = 1:dim]
    sub_iter = Iterators.product(subs...)
    nnsub = prod(nsub)
    faceinmesh_list = Vector{Array{Float64, dim-1}}[]
    vertoutlist_list = Vector{Vector{Float64}}[]

    for it in sub_iter
        xttc = xt + it.*htt - 0.5*htt
        Sxttc, DSxttc = theAbs.linsys_map(xttc)
        A1 = 0.5*DSxttc.*htt'
        mesh_list, = cube_facemeshlist(fill(-1.0, dim), fill(2.0, dim), fill(2, dim))
        fmap = (x...) -> A1*[x...] + Sxttc
        for j = 1:length(mesh_list)
            mesh_vec = map(fmap, mesh_list[j]...)
            mesh_sliced = [[x[i] for x in mesh_vec] for i = 1:theAbs.dim]
            push!(faceinmesh_list, mesh_sliced)
        end
        vertout_iter = Iterators.product(fill(-1:2:1, 2*dim)...)
        fout = x -> A1*collect(x[1:dim]) + A2*collect(x[dim+1:2*dim]) + Sxttc
        pointout_list = vec(map(fout, vertout_iter))
        hull = spatial.ConvexHull(array_to_matrix(pointout_list))
        simplices = matrix_to_array(hull.simplices .+ 1)
        for simplex in simplices
            push!(vertoutlist_list, pointout_list[simplex])
        end
    end

    return (faceinmesh_list, vertoutlist_list)
end

#=
function matrix_to_cone(P::Symmetric, rad, np)
    EV = eigen(P, 1:2)
    @assert EV.values[1] < 0.0 && EV.values[2] > 0.0
    U = EV.vectors
    d1 = EV.values[1]
    d2 = EV.values[2]
    beta = atan(sqrt(-d1/d2))
    PL0 = map(x -> [x[1], x[2]], Plots.partialcircle(-beta, beta, np, rad))
    Point_list = vcat([[0.0, 0.0]], PL0, [[0.0, 0.0]], reverse(map(x -> -x, PL0)), [[0.0, 0.0]])
    map!(x -> U * x, Point_list, Point_list)
    return Point_list
end

function plotcones(plt, PList, XList, YList, rad, np;
    linecolor = :black, linewidth = 0.5,
    fillcolor = :green, fillalpha = 0.5)

    (fillcolor, fillalpha, linecolor, linewidth) =
        preprocessargs(fillcolor, fillalpha, linecolor, linewidth)

    nidx = length(PList)
    npTot = 2*np + 4
    Cone_list = Array{Float64, 2}(undef, nidx*npTot, 2)

    for k = 1:nidx
        Point_list = matrix_to_cone(PList[k], rad, np)
        xListBis = map(p -> p[1] + XList[k], Point_list)
        yListBis = map(p -> p[2] + YList[k], Point_list)
        l = npTot * (k-1)
        Cone_list[l+1:l+npTot, :] = vcat(hcat(xListBis, yListBis), [NaN NaN])
    end

    plot!(plt, Cone_list[:, 1], Cone_list[:, 2], seriestype = :shape,
        fillcolor = fillcolor, fillalpha = fillalpha,
        linecolor = linecolor, linewidth = linewidth)
end
=#
