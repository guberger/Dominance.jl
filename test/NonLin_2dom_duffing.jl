module TestMain

include("../macros.jl")
include("../plotting.jl")

using LinearAlgebra
using Printf
using PyPlot
using PyCall
using Random
art3d = PyObject(PyPlot.art3D)
_cols = repeat(matplotlib.rcParams["axes.prop_cycle"].by_key()["color"], 100, 1)
CConv = matplotlib.colors.colorConverter
axes_grid1 = pyimport("mpl_toolkits.axes_grid1")

## Definitions

function duffing_map!(x, y)
    y[1] = x[1] + 0.3*x[2]
    y[2] = -0.15*x[1] + 0.3*sin(x[1]) + 0.7*x[2] + 0.03*x[3]
    y[3] = -1.5*x[1] + 0.925*x[3]
end

function duffing_deriv(cosx)
    return [1.0 0.3 0.0;
        (-0.15 + 0.3*cosx) 0.7 0.03;
        -1.5 0.0 0.925]
end

fig = PyPlot.figure(figsize = [9.1, 8.3])

## Cones
nmat = 4
dcosx = 1.0/nmat
ccosx_list = -1.0+dcosx:2.0*dcosx:1.0-dcosx
Ac_list = [duffing_deriv(ccosx) for ccosx in ccosx_list]
DELTA = duffing_deriv(dcosx) - duffing_deriv(0.0)
Ad_list = [ss*DELTA for ss in [-1.0, 1.0]]

ax = fig.add_axes([0.6, 0.55, 0.35, 0.35], aspect = "equal")
γ_min = 0.0
γ_max = Inf

Ap_list = [duffing_deriv(ccosx) for ccosx in -1:0.4:1.0]

for (i, A) in enumerate(Ap_list)
    ev = eigvals(A)
    EG = pth_eigen(A, 2, 1e-9)
    global γ_max = min(γ_max, EG[1])
    global γ_min = max(γ_min, EG[2])
    ax.plot(real.(ev), imag.(ev), marker = ".", ms = 10, ls = "none",
        mfc = _cols[1], mec = "k", mew = 0.5)
end

np = 100
tL = range(0.0, 2.0*pi, length = np)
verts1 = map(t -> [γ_min*cos.(t), γ_min*sin.(t)], tL)
verts2 = map(t -> [γ_max*cos.(t), γ_max*sin.(t)], reverse(tL))
verts = vcat(verts1, verts2, [[γ_min, 0.0]])
poly_list = matplotlib.collections.PolyCollection([verts], alpha = 0.25)
poly_list.set_facecolor("green")
ax.add_collection(poly_list)
ax.plot(0.0, 0.0, marker = "x", c = "k", ms = 8, mew = 3.0)

edge_list = [(i, j, i, [1, 2], [1, 1])
    for (i, j) in Iterators.product(1:nmat, 1:nmat)]
rates_list = ndgrid_array([range(γ_min, γ_max, length = 9)])

println("")
P_max, ee_max, rates_max = solve_lmi_disc_path_convex(Ac_list, Ad_list, edge_list, rates_list)
println("")

for i = 1:length(P_max)
    display(P_max[i])
    display(eigvals(P_max[i])')
end

display(rates_max)

γ = rates_max[1]
ax.plot(γ*cos.(tL), γ*sin.(tL), ls = "--", c = "green", lw = 2.0)
ax.set_xlabel("Re(λ)")
ax.set_ylabel("Im(λ)")
ax.set_xlim(-1.2, 1.2)
ax.set_ylim(-1.2, 1.2)
ax.grid(true)

## Simulations
nstep = 1000

x01 = [-1.5, -2.5, 40.0]
x02 = [-0.7, 1.7, -5.0]
x_list = [[Vector{Float64}(undef, 3) for i = 1:nstep] for k = 1:2]
x_list[1][1] = x01
x_list[2][1] = x02

for k = 1:2, i = 1:nstep-1
    duffing_map!(x_list[k][i], x_list[k][i+1])
end

ax = fig.add_axes([0.5, 0.0, 0.5, 0.55], projection = "3d")
ax.patch.set_facecolor("none")

X_plot = [[map(x -> x[i], x_list[k]) for i = 1:3] for k = 1:2]
ax.plot(X_plot[1]..., marker = ".")
ax.plot(X_plot[2]..., marker = ".")
ax.xaxis.pane.fill = false
ax.yaxis.pane.fill = false
ax.zaxis.pane.fill = false
ax.set_xlabel(L"$x_1$")
ax.set_ylabel(L"$x_2$")
ax.set_zlabel(L"$x_3$")
ax.view_init(elev = 7.5, azim = -127.25)

display(collect(ccosx_list))
xx = x_list[1]
i1 = 52
xx1 = xx[i1]
display(cos(xx1[1]))
i2 = 66
xx2 = xx[i2]
display(cos(xx2[1]))

np = 50
rad = 7.0
verts_side, verts_top = matrix_to_cone3d(-P_max[2], rad, np)
map!(x -> map(y -> y + xx1, x), verts_side, verts_side)
pcoll = make_collection(verts_side, facecolor = "green", facealpha = 0.3,
    edgecolor = "none")
ax.add_collection3d(pcoll)
verts_side, verts_top = matrix_to_cone3d(-P_max[2], rad, np)
map!(x -> map(y -> y + xx2, x), verts_side, verts_side)
pcoll = make_collection(verts_side, facecolor = "green", facealpha = 0.3,
    edgecolor = "none")
ax.add_collection3d(pcoll)

## Cone plot
x0 = 0.02
y0 = 0.04
ax1 = fig.add_axes([x0 + 0.0, y0 + 0.22, 0.22, 0.22], projection = "3d")
ax2 = fig.add_axes([x0 + 0.22, y0 + 0.22, 0.22, 0.22], projection = "3d")
ax3 = fig.add_axes([x0 + 0.0, y0 + 0.0, 0.22, 0.22], projection = "3d")
ax4 = fig.add_axes([x0 + 0.22, y0 + 0.0, 0.22, 0.22], projection = "3d")
AX_ = [ax1, ax2, ax3, ax4]

np = 50
rad = 1.0
circ = Vector{Any}(undef, 4)
LS = 0.6
LAB = [L"$P_{\mathtt{a}}$", L"$P_{\mathtt{b}}$", L"$P_{\mathtt{c}}$", L"$P_{\mathtt{d}}$"]

for i = 1:4
    ax = AX_[i]
    verts_side, verts_top = matrix_to_cone3d(-P_max[i], rad, np)
    pcoll = make_collection(verts_side, facecolor = _cols[i], facealpha = 0.3,
        edgecolor = "none")
    ax.add_collection3d(pcoll)
    circ[i] = matplotlib.patches.Patch(fc = _cols[i], ec = "none", alpha = 0.3)
    ax.patch.set_facecolor("none")
    ax.quiver3D(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, arrow_length_ratio = 0.3, length = 0.5)
    ax.quiver3D(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, arrow_length_ratio = 0.3, length = 0.5)
    ax.quiver3D(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, arrow_length_ratio = 0.3, length = 0.5)
    ax.set_xlim(-LS, LS)
    ax.set_ylim(-LS, LS)
    ax.set_zlim(-LS, LS)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    ax.text3D(-0.2*LS, -1.2*LS, -1.4*LS, L"$x_1$")
    ax.text3D(1.1*LS, 0.0, -1.4*LS, L"$x_2$")
    ax.text3D(1.2*LS, 1.1*LS, 0.0, L"$x_3$")
    ax.xaxis.pane.fill = false
    ax.yaxis.pane.fill = false
    ax.zaxis.pane.fill = false
    # ax.axis("off")
    ax.text3D(-0.5*LS, -LS, LS, LAB[i], color = _cols[i],
        bbox = Dict([("facecolor", "white"), ("alpha", 0.5)]),
        fontsize = 15)
end

ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
ax.patch.set_facecolor("none")
ax.text(0.04, 0.47, "c", weight = "bold", fontsize = 15)
ax.text(0.52, 0.47, "d", weight = "bold", fontsize = 15)
ax.text(0.52, 0.9, "b", weight = "bold", fontsize = 15)
ax.text(0.04, 0.9, "a", weight = "bold", fontsize = 15)

fig.savefig("test/figures/fig_NonLin_2dom_duffing_allllllll.png",
    transparent = false,
    bbox_inches = matplotlib.transforms.Bbox([[0.3, 0.3], [8.9, 7.7]]))
end
