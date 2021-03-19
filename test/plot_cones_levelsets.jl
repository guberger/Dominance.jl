module TestMain

include("../plotting.jl")

using LinearAlgebra
using Printf
using PyPlot
using PyCall
art3d = PyObject(PyPlot.art3D)
_cols = repeat(matplotlib.rcParams["axes.prop_cycle"].by_key()["color"], 10, 1)
CConv = matplotlib.colors.colorConverter
axes_grid1 = pyimport("mpl_toolkits.axes_grid1")

P = [-0.5 0; 0 0.25]
V(x, y) = [x y]*P*[x; y]

np = 100
xl = 1.0
yl = 1.0
x = range(-xl, xl, length = np)
y = range(-yl, yl, length = np)
IT = collect(Iterators.product(x, y))
X = map(x -> x[1], IT)
Y = map(x -> x[2], IT)
Z = map((x, y) -> V(x, y)[1], X, Y)

fig = PyPlot.figure(figsize = (9.4, 4.8))
gs = matplotlib.gridspec.GridSpec(1, 2, figure = fig,
    width_ratios = (1, 1.5), wspace = 0.2)
ax = fig.add_subplot(get(gs, 0), aspect = "equal")

CS = ax.contour(X, Y, Z)
ax.clabel(CS, CS.levels, inline = true, fmt = "%2.1f", fontsize = 14)

points = matrix_to_cone2d(P, 2.0, 10)
poly_list = matplotlib.collections.PolyCollection(points, alpha = 0.25)
poly_list.set_facecolor("blue")
ax.add_collection(poly_list)
# points = matrix_to_cone2d(-P, 2.0, 10)
# poly_list = matplotlib.collections.PolyCollection(points, alpha = 0.25)
# poly_list.set_facecolor("red")
# ax.add_collection(poly_list)

ax.set_xlim(-1.0, 1.0)
ax.set_ylim(-1.0, 1.0)
ax.set_xlabel("x")
ax.set_ylabel("y")

P = [1.0 0 0; 0 1.0 0; 0 0 -2.0]
np = 50
rad = sqrt(2)

ang = range(0.5*pi, 2.0*pi, length = np + 1)
arc = map(x -> rad*[cos(x), sin(x), sqrt(0.5)], ang)
points = [[[0.0, 0.0, 0.0], arc[i], arc[i+1]] for i = 1:np]
append!(points, map(x -> map(y -> [y[1], y[2], -y[3]], x), points))

ax = fig.add_subplot(get(gs, 1), projection = "3d")
pcoll = make_collection(points, facecolor = "blue", facealpha = 0.4,
    edgecolor = "none")
ax.add_collection3d(pcoll)

function levelcurves(a, b, np, lim1, lim2, levels)
    points = Vector{Vector{Vector{Float64}}}(undef, length(levels))
    for (i, lev) in enumerate(levels)
        if lev >= 0.0
            ps = _lcurves(b, a, np, lim2, lim1, lev)
            points[i] = [ps[2], ps[1]]
        else
            points[i] = _lcurves(a, b, np, lim1, lim2, -lev)
        end
    end
    return points
end

function _lcurves(a, b, np, lim1, lim2, lev)
    @assert lev >= 0.0
    if lev > b*lim2^2
        return [Float64[] for i = 1:2]
    end
    l1 = min(lim1, sqrt((b*lim2^2 - lev)/a))
    np = ceil(Int, (l1/lim1 + 1)*np)
    xL = collect(range(0.0, l1, length = np))
    yL = map(x -> sqrt((a*x^2 + lev)/b), xL)
    return [xL, yL]
end

levels = -2.0:0.4:2.0
lim1 = sqrt(2)
lim2 = 1.0
cmap = matplotlib.cm.get_cmap("viridis")
scal = matplotlib.cm.ScalarMappable(cmap = cmap)
scal.set_array(-levels)
fc_mat = scal.to_rgba(-levels)
fc_list = [Tuple(fc_mat[i, :]) for i = 1:size(fc_mat, 1)]

points = levelcurves(1.0, 2.0, 20, lim1, lim2, levels)

for (i, ps) in enumerate(points)
    ax.plot(ps[1], zeros(length(ps[1])), ps[2], c = fc_list[i], lw = 2.5)
    ax.plot(ps[1], zeros(length(ps[1])), -ps[2], c = fc_list[i], lw = 2.5)
    ax.plot(zeros(length(ps[1])), ps[1], ps[2], c = fc_list[i], lw = 2.5)
    ax.plot(zeros(length(ps[1])), ps[1], -ps[2], c = fc_list[i], lw = 2.5)
end

fig.colorbar(scal, shrink = 0.7)

ax.plot([0.0, lim1, lim1, 0.0, 0.0], zeros(5),
    [-lim2, -lim2, lim2, lim2, -lim2], c = "grey", lw = 1.0)
ax.plot(zeros(5), [0.0, lim1, lim1, 0.0, 0.0],
    [-lim2, -lim2, lim2, lim2, -lim2], c = "grey", lw = 1.0)
ax.plot(map(x -> x[1], arc), map(x -> x[2], arc), map(x -> x[3], arc),
    c = "grey", lw = 1.0)
ax.plot(map(x -> x[1], arc), map(x -> x[2], arc), map(x -> -x[3], arc),
    c = "grey", lw = 1.0)

ax.set_xlim(-1.2, 1.2)
ax.set_ylim(-1.2, 1.2)
ax.set_zlim(-1.2, 1.2)
ax.xaxis.pane.fill = false
ax.yaxis.pane.fill = false
ax.zaxis.pane.fill = false
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.view_init(elev = 14.0, azim = 52.0)
fig.savefig("./figures/fig_cones_levelsets.png", transparent = false,
    bbox_inches = matplotlib.transforms.Bbox([[0.45, 0.6], [8.45, 3.8]]))
end
