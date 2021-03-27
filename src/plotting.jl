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

FC(c, a) =  matplotlib.colors.colorConverter.to_rgba(c, alpha = a)
const _mult = (SVector(-1,-1), SVector(-1,1), SVector(1,1), SVector(1,-1))

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
    fca = FC(fc, fa)
    eca = FC(ec, ea)
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
    fca = FC(fc, fa)
    eca = FC(ec, ea)
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

    for i in 1:nstep
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
    fca = FC(fc, fa)
    eca = FC(ec, ea)
    vertslist = Vector{SVector{2,T}}[]
    ns = nsub .- 1

    for pos in DO.enum_pos(domain)
        x = DO.get_coord_by_pos(domain.grid, pos)
        subpos_axes = ((0:ns[i])./ns[i] .- 0.5 for i in 1:length(ns))
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
    fca = FC(fc, fa)
    eca = FC(ec, ea)
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

end  # Plot
