# All plotting functions

module Plot

import ..Dominance
DO = Dominance

using LinearAlgebra
using StaticArrays
using LazySets
using Polyhedra
using CDDLib
using PyPlot
using PyCall

art3d = PyObject(PyPlot.art3D)
ColorConv(c, a) =  matplotlib.colors.colorConverter.to_rgba(c, alpha = a)
const _colors = repeat(matplotlib.rcParams["axes.prop_cycle"].by_key()["color"], 5)
const _mult = (SVector(-1,-1), SVector(-1,1), SVector(1,1), SVector(1,-1))

## Abstraction

function verts_rect(c, h)
    return ntuple(i -> c + h.*_mult[i], 4)
end

function project(x, vars)
    return SVector(x[vars[1]], x[vars[2]])
end

# Cells
function domain!(ax, vars, domain::DO.Domain{N,T};
        fc = "red", fa = 0.5, ec = "black", ea = 1.0, ew = 1.5) where {N,T}
    grid = domain.grid
    @assert length(vars) == 2 && N >= 2
    fca = ColorConv(fc, fa)
    eca = ColorConv(ec, ea)
    h = project(grid.h, vars)

    vertslist = NTuple{4,SVector{2,T}}[]

    for pos in unique(x -> x[vars], DO.enum_pos(domain))
        c = project(DO.get_coord_by_pos(grid, pos), vars)
        push!(vertslist, verts_rect(c, h/2.0))
    end

    polylist = matplotlib.collections.PolyCollection(vertslist)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(eca)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
end

# Sets
function set!(ax, vars, rect::DO.HyperRectangle;
        fc = "green", fa = 0.5, ec = "black", ea = 1.0, ew = 1.5)
    @assert length(vars) == 2 && length(rect.lb) == length(rect.ub) >= 2
    c = (rect.lb + rect.ub)/2.0
    h = (rect.ub - rect.lb)/2.0
    polylist = matplotlib.collections.PolyCollection((verts_rect(c[vars], h[vars]),))
    fca = ColorConv(fc, fa)
    eca = ColorConv(ec, ea)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(eca)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
end

# Trajectory open loop
function trajectory!(ax, vars, sys, x0, nstep;
        lc = "red", lw = 1.5, mc = "black", ms = 5.0)
    @assert length(vars) == length(x0) == 2
    X1list = Vector{eltype(x0)}(undef, nstep + 1)
    X2list = Vector{eltype(x0)}(undef, nstep + 1)
    X1list[1] = x0[vars[1]]
    X2list[1] = x0[vars[2]]
    x = x0

    for i = 1:nstep
        x = sys.sys_map(x)
        X1list[i + 1] = x[vars[1]]
        X2list[i + 1] = x[vars[2]]
    end

    ax.plot(X1list, X2list, lw = lw, c = lc,
        marker = ".", ms = ms, mfc = mc, mec = mc)
end

# Images
function cell_image!(ax, vars, domain::DO.Domain{N,T}, sys;
        nsub = ntuple(i -> 5, Val(N)),
        fc = "blue", fa = 0.5, ec = "darkblue", ea = 1.0, ew = 1.5) where {N,T}
    @assert length(vars) == 2 && N >= 2
    fca = ColorConv(fc, fa)
    eca = ColorConv(ec, ea)
    vertslist = Vector{SVector{2,T}}[]
    ns = nsub .- 1

    for pos in DO.enum_pos(domain)
        x = DO.get_coord_by_pos(domain.grid, pos)
        subpos_axes = ((0:ns[i])./ns[i] .- 0.5 for i in eachindex(ns))
        subpos_iter = Iterators.product(subpos_axes...)
        x_iter = (x + subpos.*domain.grid.h for subpos in subpos_iter)
        Fx_iter = (sys.sys_map(x) for x in x_iter)
        push!(vertslist, convex_hull([project(Fx, vars) for Fx in Fx_iter][:]))
    end

    polylist = matplotlib.collections.PolyCollection(vertslist)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(eca)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
end

function cell_approx!(ax, vars, domain::DO.Domain{N,T}, sys;
        fc = "yellow", fa = 0.5, ec = "gold", ea = 1.0, ew = 0.5) where {N,T}
    @assert length(vars) == 2 && N >= 2
    fca = ColorConv(fc, fa)
    eca = ColorConv(ec, ea)
    vertslist = Vector{SVector{2,T}}[]
    r = domain.grid.h/2
    _H_ = SMatrix{N,N}(I).*r
    _ONE_ = ones(SVector{N})
    Fe = sys.error_map(norm(r, Inf))
    Fr = typeof(r)(_ONE_*Fe)

    for pos in DO.enum_pos(domain)
        x = DO.get_coord_by_pos(domain.grid, pos)
        Fx, DFx = sys.linsys_map(x, _H_)
        A = inv(DFx)
        b = abs.(A)*Fr .+ 1
        HP = HPolytope([A; -A], vcat(b, b))
        verts1 = vertices_list(HP, backend = CDDLib.Library())
        verts2 = [project(Fx + v, vars) for v in verts1]
        push!(vertslist, convex_hull(verts2))
    end

    polylist = matplotlib.collections.PolyCollection(vertslist)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(eca)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
end

## Cones

function make_collection(verts::Vector{Vector{SVector{N,T}}};
        fc = "blue", fa = 0.5, ec = "black", ea = 1.0, lw = 1.0) where {N,T}
    if fc == "none"
        fa = 0.0
    end
    fca = ColorConv(fc, fa)
    eca = ColorConv(ec, ea)
    if ec == "none"
        eca = "none"
        lw = 0.0
    end
    poly_coll = N == 2 ? matplotlib.collections.PolyCollection(verts) :
        N == 3 ? art3d.Poly3DCollection(verts) :
        error("dimension must be 2 or 3")
    poly_coll.set_facecolor(fca)
    poly_coll.set_edgecolor(eca)
    poly_coll.set_linewidth(lw)
    return poly_coll
end

function matrix_to_cone2d(P::SMatrix{2,2}, rad, np)
    EV = eigen(Symmetric(P))
    evals = EV.values
    @assert evals[1] < 0 && evals[2] > 0
    U = EV.vectors
    beta = atan(sqrt(-evals[1]/evals[2]))
    ang = range(-beta, beta, length = np)
    arc = map(x -> rad*SVector(cos(x), sin(x)), ang)
    arcrev = reverse(map(x -> -x, arc))
    verts1 = vcat([SVector(0.0, 0.0)], arc, [SVector(0.0, 0.0)])
    map!(x -> U*x, verts1, verts1)
    verts2 = map(x -> -x, verts1)
    return [verts1, verts2]
end

_mirror(vert) = SVector(-1, 1, 1).*vert

function draw_icecream3d(rad, np, ang_start, ang_stop)
    ang = range(ang_start, ang_stop, length = np + 1)
    arcpos = map(x -> rad*SVector(1.0, cos(x), sin(x)), ang)
    verts_side = [[SVector(0.0, 0.0, 0.0), arcpos[i], arcpos[i+1]] for i = 1:np]
    append!(verts_side, map(verts -> _mirror.(verts), verts_side))
    arcneg = map(x -> rad*SVector(-1.0, cos(x), sin(x)), ang)
    verts_top = [arcpos, arcneg]
    return verts_side, verts_top
end

function matrix_to_cone3d(P::SMatrix{3,3}, rad, np;
        ang_start = 0.0, ang_stop = 2*pi)
    EV = eigen(Symmetric(P))
    evals = EV.values
    @assert evals[1] < 0 && evals[2] > 0 && evals[3] > 0
    evals = abs.(evals)./(-evals[1])
    verts_side, verts_top = draw_icecream3d(rad, np, ang_start, ang_stop)
    U = EV.vectors./sqrt.(evals')
    map!(x -> map(y -> U*y, x), verts_side, verts_side)
    map!(x -> map(y -> U*y, x), verts_top, verts_top)
    return verts_side, verts_top
end

function cones!(ax, grid::DO.Grid{2}, sys, x0, idxn, P_field, nsteps, rad, np;
        fc1 = "blue", fa1 = 0.5, ec1 = "blue", ea1 = 1.0, lw1 = 0.5,
        fc2 = "green", fa2 = 0.5, ec2 = "green", ea2 = 1.0, lw2 = 0.5,
        fact = 1.1)
    @assert length(x0) == 2
    x = x0
    _H_ = SMatrix{2,2}(I)
    for i = 0:nsteps
        Fx, DFx = sys.linsys_map(x, _H_)
        DFx = DFx/opnorm(DFx)
        pos = DO.get_pos_by_coord(grid, x)
        index = DO.get_index_by_elem(idxn, pos)
        P = P_field[index]
        verts = Plot.matrix_to_cone2d(P, rad, np)
        f_shift1(vert) = vert + x
        map!(verts -> f_shift1.(verts), verts, verts)
        poly_list = Plot.make_collection(verts,
            fc = fc1, fa = fa1, ec = ec1, ea = ea1, lw = lw1)
        ax.add_collection(poly_list)
        i == nsteps && break
        FP = DFx'\P/DFx
        verts = Plot.matrix_to_cone2d(FP, rad*fact, np)
        f_shift2(vert) = vert + Fx
        map!(verts -> f_shift2.(verts), verts, verts)
        poly_list = Plot.make_collection(verts,
            fc = fc2, fa = fa2, ec = ec2, ea = ea2, lw = lw2)
        ax.add_collection(poly_list)
        x = Fx
    end
end

end  # Plot
