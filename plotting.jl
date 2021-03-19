using LinearAlgebra
using PyPlot
using PyCall
art3d = PyObject(PyPlot.art3D)
CConv = matplotlib.colors.colorConverter

# function proj_transform(vec, axM)
#     v = [vec..., 1.0]
#     w = axM*v
#     return w[1:3]/w[4]
# end

function matrix_to_cone2d(P, rad, np)
    EV = eigen(Symmetric(P), 1:2)
    evals = EV.values
    @assert evals[1] < 0.0 && evals[2] > 0.0
    U = EV.vectors
    beta = atan(sqrt(-evals[1]/evals[2]))
    ang = range(-beta, beta, length = np)
    arc = map(x -> rad*[cos(x), sin(x)], ang)
    arcrev = reverse(map(x -> -x, arc))
    verts1 = vcat([[0.0, 0.0]], arc, [[0.0, 0.0]])
    map!(x -> U * x, verts1, verts1)
    verts2 = map(x -> -x, verts1)
    return [verts1, verts2]
end

function make_collection(verts; facecolor = "blue", facealpha = 0.5,
    edgecolor = "black", edgealpha = 1.0, linewidth = 1.0)
    if facecolor == "none"
        facealpha = 0.0
    end
    fca = CConv.to_rgba(facecolor, alpha = facealpha)
    eca = CConv.to_rgba(edgecolor, alpha = edgealpha)
    if edgecolor == "none"
        eca = "none"
        linewidth = 0.0
    end
    poly_coll = if length(verts[1][1]) == 2
        matplotlib.collections.PolyCollection(verts)
    elseif length(verts[1][1]) == 3
        art3d.Poly3DCollection(verts)
    else
        error("dimension must be 2 or 3")
    end
    poly_coll.set_facecolor(fca)
    poly_coll.set_edgecolor(eca)
    poly_coll.set_linewidth(linewidth)
    return poly_coll
end

#=
function _make_collection3d(verts; facecolor = "blue", facealpha = 0.5,
    edgecolor = "black", edgealpha = 1.0, linewidth = 1.0)
    if facecolor == "none"
        facealpha = 0.0
    end
    fca = CConv.to_rgba(facecolor, alpha = facealpha)
    eca = CConv.to_rgba(edgecolor, alpha = edgealpha)
    if edgecolor == "none"
        eca = "none"
        linewidth = 0.0
    end
    poly_coll = art3d.Poly3DCollection(verts)
    poly_coll.set_facecolor(fca)
    poly_coll.set_edgecolor(eca)
    poly_coll.set_linewidth(linewidth)
    return poly_coll
end
=#

function draw_icecream3d(rad, np)
    ang = range(0.0, 2.0*pi, length = np + 1)
    arc = map(x -> rad*[1.0, cos(x), sin(x)], ang)
    verts_side = [[[0.0, 0.0, 0.0], arc[i], arc[i+1]] for i = 1:np]
    append!(verts_side, map(x -> map(y -> -y, x), verts_side))
    arcmirror = map(x -> rad*[-1.0, cos(x), sin(x)], ang)
    verts_top = [arc, arcmirror]
    return (verts_side, verts_top)
end

function normalize_cone(verts, rad)
    fnorm(x) = norm(x) > 1e-8 ? rad*normalize(x) : 0.0*x
    return map(x -> map(fnorm, x), verts)
end

#=
function draw_diabolo3d(rad, np)
    ang = range(0.0, 2.0*pi, length = np + 1)
    arc = map(x -> rad*[1.0, cos(x), sin(x)], ang)
    verts = [[[0.0, 0.0, 0.0], arc[i], arc[i+1]] for i = 1:np]
    append!(verts, map(x -> map(y -> -y, x), verts))
    arcmirror = map(x -> rad*[-1.0, cos(x), sin(x)], ang)
    strip = [[arc[i], arc[i+1], arcmirror[i+1], arcmirror[i]] for i = 1:np]
    append!(verts, strip)
    return verts
end
=#

function matrix_to_cone3d(P, rad, np)
    EV = eigen(Symmetric(P), 1:3)
    evals = EV.values
    @assert evals[1] < 0.0 && evals[2] > 0.0 && evals[3] > 0.0
    evals = abs.(evals)./(-evals[1])
    verts_side, verts_top = draw_icecream3d(rad, np)
    U = EV.vectors./sqrt.(evals')
    map!(x -> map(y -> U*y, x), verts_side, verts_side)
    map!(x -> map(y -> U*y, x), verts_top, verts_top)
    return verts_side, verts_top
end

#=
function matrix_to_ellipse2d(P, rad, np)
    EV = eigen(Symmetric(P), 1:2)
    evals = EV.values
    @assert evals[1] > 0.0 && evals[2] > 0.0
    evals = evals./evals[1]
    U = EV.vectors./sqrt.(evals')
    ang = range(0.0, 2.0*pi, length = np + 1)
    arc = map(x -> rad*[cos(x), sin(x)], ang)
    points = map(x -> U*x, arc)
    return points
end
=#
