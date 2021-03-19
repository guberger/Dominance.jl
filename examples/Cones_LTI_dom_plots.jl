module TestMain

include("../macros.jl")
include("../plotting.jl")

using LinearAlgebra
using Printf
using PyPlot
using PyCall
using Random
art3d = PyObject(PyPlot.art3D)
_cols = repeat(matplotlib.rcParams["axes.prop_cycle"].by_key()["color"], 10, 1)
CConv = matplotlib.colors.colorConverter
axes_grid1 = pyimport("mpl_toolkits.axes_grid1")

fig = PyPlot.figure(figsize = (9.0, 10.5))

## Cones

zlevel1 = 0.65

ax = fig.add_axes((0.1, zlevel1, 0.35, 0.35), aspect = "equal")

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
ax.set_xticks(Tuple(-1.0:0.5:1.0))
ax.set_yticks(Tuple(-1.0:0.5:1.0))

P = [1.0 0 0; 0 1.0 0; 0 0 -2.0]
np = 50
rad = sqrt(2)

ang = range(0.5*pi, 2.0*pi, length = np + 1)
arc = map(x -> rad*[cos(x), sin(x), sqrt(0.5)], ang)
points = [[[0.0, 0.0, 0.0], arc[i], arc[i+1]] for i = 1:np]
append!(points, map(x -> map(y -> [y[1], y[2], -y[3]], x), points))

ax = fig.add_axes((0.5, zlevel1, 0.45, 0.35), projection = "3d")

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
ax.set_xticks(Tuple(-1.0:0.5:1.0))
ax.set_yticks(Tuple(-1.0:0.5:1.0))
ax.set_zticks(Tuple(-1.0:0.5:1.0))

## LTI 1Dom

zlevel2 = zlevel1 - 0.38

ax = fig.add_axes((0.1, zlevel2, 0.35, 0.35), aspect = "equal")

α = 0.2 # limit: (6.0 - sqrt(32.0))/2.0 = 0.1715728752538097
A = [α α-1.0; 0.0 1.0]

np = 50
t = range(0.0, 2.0*pi, length = np)

ev = eigvals(A)
γ = 0.6

ax.plot(γ*cos.(t), γ*sin.(t), ls = "--", c = "green", lw = 3.0,
    label = L"$\gamma = 0.6$")
ax.plot(real.(ev), imag.(ev), marker = ".", ms = 18, ls = "none",
    mfc = _cols[1], mec = "black", mew = 1.5)

ax.legend(fontsize = 14, ncol = 2)
ax.grid(true)
ax.set_xlabel("Re(λ)")
ax.set_ylabel("Im(λ)")
ax.set_xlim(-1.2, 1.2)
ax.set_ylim(-1.2, 1.2)
ax.text(α-0.05, -0.4, "α", fontsize = 14)
ax.plot([α, α], [-0.05, -0.28], c = "black")

A_list = [A]
edge_list = [(1, 1, 1, 1)]
rates_list = ndgrid_array([γ])
display(rates_list)

println("")
P_max, ee_max, rates_max = solve_lmi_disc_path(A_list, edge_list, rates_list)
println("")

np = 50
rad = 1.0
fact = 1.1

ax = fig.add_axes((0.55, zlevel2, 0.35, 0.35), aspect = "equal")

verts = matrix_to_cone2d(P_max[1], rad, np)
pcoll = make_collection(verts, facecolor = _cols[1], facealpha = 0.75,
    edgecolor = _cols[1])
ax.add_collection(pcoll)

verts = matrix_to_cone2d(A'\P_max[1]/A, rad*fact, np)
pcoll = make_collection(verts, facecolor = _cols[2], facealpha = 0.75,
    edgecolor = _cols[2])
ax.add_collection(pcoll)
ax.set_xlim(-1.3, 1.3)
ax.set_ylim(-1.3, 1.3)
ax.set_yticks(-1.0:0.5:1.0)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.grid(true)

circ = [matplotlib.patches.Patch(fc = _cols[i], ec = _cols[i], alpha = 0.8)
    for i = 1:2]

LAB = [L"$\mathcal{K}(P)$", L"$A\mathcal{K}(P)$"]
ax.legend(circ, LAB, fontsize = 14)

## LTI 2Dom

zlevel3 = zlevel2 - 0.38

ax = fig.add_axes((0.1, zlevel3, 0.35, 0.35), aspect = "equal")

A = [1 0.5 0; -0.5 0.75 0.5; -0.5 0 1.0]

np = 50
t = range(0.0, 2.0*pi, length = np)

ev = eigvals(A)
γ = 0.9

ax.plot(γ*cos.(t), γ*sin.(t), ls = "--", c = "green", lw = 3.0,
    label = "\$\\gamma = $(γ)\$")
ax.plot(real.(ev), imag.(ev), marker = ".", ms = 18, ls = "none",
    mfc = _cols[1], mec = "black", mew = 1.5)

ax.legend(fontsize = 14, ncol = 2)
ax.grid(true)
ax.set_xlabel("Re(λ)")
ax.set_ylabel("Im(λ)")
ax.set_xlim(-1.75, 1.75)
ax.set_ylim(-1.75, 1.75)
ax.set_xticks(Tuple(-1.5:0.5:1.5))
ax.set_yticks(Tuple(-1.5:0.5:1.5))

A_list = [A]
edge_list = [(1, 1, 1, 1)]
rates_list = ndgrid_array([γ])
display(rates_list)

println("")
P_max, ee_max, rates_max = solve_lmi_disc_path(A_list, edge_list, rates_list)
println("")

np = 50
rad = 1.0

ax = fig.add_axes((0.5, zlevel3, 0.42, 0.35), projection = "3d")

verts_side, verts_top = matrix_to_cone3d(-P_max[1], rad, np)
pcoll = make_collection(verts_side, facecolor = _cols[1], facealpha = 0.35,
    edgecolor = "none")
ax.add_collection3d(pcoll)

verts_side = map(x -> map(y -> A*y, x), verts_side)
pcoll = make_collection(verts_side, facecolor = _cols[2], facealpha = 0.35,
    edgecolor = "none")
ax.add_collection3d(pcoll)
ax.set_xlim(-1.2, 1.2)
ax.set_ylim(-1.2, 1.2)
ax.set_zlim(-1.2, 1.2)
ax.xaxis.pane.fill = false
ax.yaxis.pane.fill = false
ax.zaxis.pane.fill = false
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.view_init(elev = 20.0, azim = 60.0)

circ = [matplotlib.patches.Patch(fc = _cols[i], ec = _cols[i], alpha = 0.8)
    for i = 1:2]

LAB = [L"$\mathcal{K}(P)$", L"$A\mathcal{K}(P)$"]
ax.legend(circ, LAB, fontsize = 14)

fig.text(0.03, zlevel1 + 0.32, "a", weight = "bold", fontsize = 15)
fig.text(0.48, zlevel1 + 0.32, "b", weight = "bold", fontsize = 15)
fig.text(0.03, zlevel2 + 0.32, "c", weight = "bold", fontsize = 15)
fig.text(0.48, zlevel2 + 0.32, "d", weight = "bold", fontsize = 15)
fig.text(0.03, zlevel3 + 0.32, "e", weight = "bold", fontsize = 15)
fig.text(0.48, zlevel3 + 0.32, "f", weight = "bold", fontsize = 15)

fig.savefig("./figures/fig_cones_LTI_dom.png",
    transparent = false, bbox_inches = "tight")

end
