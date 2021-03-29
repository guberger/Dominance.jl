include("../src/Dominance.jl")

module ExampleMain

using LinearAlgebra
using StaticArrays
using PyPlot
using PyCall
using Main.Dominance
DO = Main.Dominance

include("../src/plotting.jl")

sleep(0.1) # used for good printing
println("Plot cones levelsets")

matplotlib.rc("legend", fontsize = 15)
matplotlib.rc("axes", labelsize = 15)
matplotlib.rc("xtick", labelsize = 11)
matplotlib.rc("ytick", labelsize = 11)

P = @SMatrix [-0.5 0; 0 0.25]
VP(x, y) = P[1, 1]*x*x + (P[1, 2] + P[2, 1])*x*y + P[2, 2]*y*y

np = 100
xl = 1.0
yl = 1.0
x_vec = range(-xl, xl, length = np)
y_vec = range(-yl, yl, length = np)
IT = collect(Iterators.product(x_vec, y_vec))
X = map(x -> x[1], IT)
Y = map(x -> x[2], IT)
Z = map(VP, X, Y)

fig = PyPlot.figure(figsize = (9.8, 4.8))
gs = matplotlib.gridspec.GridSpec(1, 2, figure = fig,
    width_ratios = (1, 1.5), wspace = 0.25)
ax = fig.add_subplot(get(gs, 0), aspect = "equal")

CS = ax.contour(X, Y, Z)
ax.clabel(CS, CS.levels, inline = true, fmt = "%2.1f", fontsize = 14)

verts = Plot.matrix_to_cone2d(P, 2.0, 10)
poly_list = Plot.make_collection(verts, fc = "blue", fa = 0.25, ec = "none")
ax.add_collection(poly_list)

ax.set_xlim(-1.0, 1.0)
ax.set_ylim(-1.0, 1.0)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xticks(-1.0:0.5:1.0)
ax.set_yticks(-1.0:0.5:1.0)

P = @SMatrix [1.0 0 0; 0 1.0 0; 0 0 -2.0]
np = 50
rad = sqrt(1.0)

ax = fig.add_subplot(get(gs, 1), projection = "3d")
verts, ~ = Plot.matrix_to_cone3d(P, rad, np, ang_start = 0.5*pi)
poly_list = Plot.make_collection(verts, fc = "blue", fa = 0.4, ec = "none")
ax.add_collection3d(poly_list)

function levelcurves(a, b, np, lim1, lim2, levels)
    verts = Vector{NTuple{2,Vector{Float64}}}(undef, length(levels))
    for (i, lev) in enumerate(levels)
        verts[i] = lev >= 0.0 ? reverse(_lcurves(b, a, np, lim2, lim1, lev)) :
            _lcurves(a, b, np, lim1, lim2, -lev)
    end
    return verts
end

function _lcurves(a, b, np, lim1, lim2, lev)
    @assert lev >= 0.0
    if lev > b*lim2^2
        return Float64[], Float64[]
    end
    l1 = min(lim1, sqrt((b*lim2^2 - lev)/a))
    np = ceil(Int, (l1/lim1 + 1)*np)
    xL = collect(range(0.0, l1, length = np))
    yL = map(x -> sqrt((a*x^2 + lev)/b), xL)
    return (xL, yL)
end

levels = -2.0:0.4:2.0
lim1 = sqrt(2)
lim2 = 1.0
cmap = matplotlib.cm.get_cmap("viridis")
scal = matplotlib.cm.ScalarMappable(cmap = cmap)
scal.set_array(-levels)
fc_mat = scal.to_rgba(-levels)
fc_list = [Tuple(fc_mat[i, :]) for i = 1:size(fc_mat, 1)]

verts = levelcurves(1.0, 2.0, 20, lim1, lim2, levels)

for (i, ps) in enumerate(verts)
    ax.plot(ps[1], zeros(length(ps[1])), ps[2], c = fc_list[i], lw = 2.5)
    ax.plot(ps[1], zeros(length(ps[1])), -ps[2], c = fc_list[i], lw = 2.5)
    ax.plot(zeros(length(ps[1])), ps[1], ps[2], c = fc_list[i], lw = 2.5)
    ax.plot(zeros(length(ps[1])), ps[1], -ps[2], c = fc_list[i], lw = 2.5)
end

fig.colorbar(scal, shrink = 0.7)

# arc = map(x -> rad*SVector(cos(x), sin(x), sqrt(0.5)), ang)
# ax.plot([0.0, lim1, lim1, 0.0, 0.0], zeros(5),
#     [-lim2, -lim2, lim2, lim2, -lim2], c = "grey", lw = 1.0)
# ax.plot(zeros(5), [0.0, lim1, lim1, 0.0, 0.0],
#     [-lim2, -lim2, lim2, lim2, -lim2], c = "grey", lw = 1.0)
# ax.plot(map(x -> x[1], arc), map(x -> x[2], arc), map(x -> x[3], arc),
#     c = "grey", lw = 1.0)
# ax.plot(map(x -> x[1], arc), map(x -> x[2], arc), map(x -> -x[3], arc),
#     c = "grey", lw = 1.0)

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

# zlevel = 0.79
# fig.text(0.05, zlevel, "a", weight = "bold", fontsize = 15)
# fig.text(0.46, zlevel, "b", weight = "bold", fontsize = 15)

fig.savefig("./figures/fig_cones_levelsets.png", transparent = false,
    bbox_inches = matplotlib.transforms.Bbox(((0.5, 0.53), (8.85, 3.8))))
end  # module ExampleMain
