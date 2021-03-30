include("../src/Dominance.jl")

module ExampleMain

using LinearAlgebra
using StaticArrays
using PyPlot
using PyCall
using JuMP
using MosekTools
using Random
using Main.Dominance
DO = Main.Dominance

include("../src/plotting.jl")

sleep(0.1) # used for good printing
println("Plot SLS 1-dom simple")

matplotlib.rc("legend", fontsize = 15)
matplotlib.rc("axes", labelsize = 15)
matplotlib.rc("xtick", labelsize = 11)
matplotlib.rc("ytick", labelsize = 11)

## Definitions

function duffing_map(x)
    return SVector(x[1] + 0.3*x[2],
        -0.15*x[1] + 0.3*sin(x[1]) + 0.7*x[2] + 0.03*x[3],
        -1.5*x[1] + 0.925*x[3])
end

function duffing_deriv(cosx)
    return @SMatrix [1.0 0.3 0.0;
        (-0.15 + 0.3*cosx) 0.7 0.03;
        -1.5 0.0 0.925]
end

## Cones
nregions = 4
dcosx = 1.0/nregions
ccosx_list = range(-1.0+dcosx, 1.0-dcosx, length = nregions)
Ac_list = [duffing_deriv(ccosx) for ccosx in ccosx_list]
Ad = duffing_deriv(dcosx) - duffing_deriv(0.0)

graph = DO.Graph(nregions)
ASri_tmp = Any[]
for i = 1:nregions, j = 1:nregions
    DO.add_edge!(graph, i, j)
    push!(ASri_tmp, DO.Edge(i, j) => [(DO.MatrixSet(Ac_list[i], [Ad, -Ad]), 1)])
end
ASri_field = Dict(ASri_tmp)

fig = PyPlot.figure(figsize = (9.1, 8.3))
ax = fig.add_axes((0.6, 0.55, 0.35, 0.35), aspect = "equal")

Aev_list = [duffing_deriv(ccosx) for ccosx in -1.0:0.4:1.0]
γ_min = 0.0
γ_max = Inf
for A in Aev_list
    ev = eigvals(Matrix(A))
    EG = DO.pth_eigval(Matrix(A), 2, 1e-9)
    global γ_max = min(γ_max, EG[1])
    global γ_min = max(γ_min, EG[2])
    ax.plot(real.(ev), imag.(ev), marker = ".", ms = 10, ls = "none",
        mfc = Plot._colors[1], mec = "k", mew = 0.5)
end
nr = 9
rate_tuple_iter = DO.hyper_range((γ_min,), (γ_max,), nr)

np = 100
ang = range(0.0, 2.0*pi, length = np)
verts1 = map(t -> SVector(γ_min*cos(t), γ_min*sin(t)), ang)
verts2 = map(t -> SVector(γ_max*cos(t), γ_max*sin(t)), reverse(ang))
verts = vcat(verts1, verts2, [SVector(γ_min, 0.0)])
poly_list = Plot.make_collection([verts], fc = "green", fa = 0.25, ec = "none")
ax.add_collection(poly_list)
ax.plot(0.0, 0.0, marker = "x", c = "k", ms = 8, mew = 3.0)

optim_solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
P_opt, δ_opt, rates_opt = DO.cone_optim(graph, ASri_field, rate_tuple_iter, optim_solver)

for i in eachindex(P_opt)
    println(eigvals(P_opt[i]))
end
println(δ_opt)
println(rates_opt)

γ = rates_opt[1]
ax.plot(γ*cos.(ang), γ*sin.(ang), ls = "--", c = "green", lw = 2.0)
ax.set_xlabel("Re(λ)")
ax.set_ylabel("Im(λ)")
ax.set_xlim(-1.2, 1.2)
ax.set_ylim(-1.2, 1.2)
ax.grid(true)

## Trajectories
nstep = 1000

x01 = SVector(-1.5, -2.5, 40.0)
x02 = SVector(-0.7, 1.7, -5.0)
x_list = [Vector{SVector{3,Float64}}(undef, nstep) for k = 1:2]
x_list[1][1] = x01
x_list[2][1] = x02

for k = 1:2, i = 1:nstep-1
    x_list[k][i+1] = duffing_map(x_list[k][i])
end

ax = fig.add_axes((0.5, 0.0, 0.5, 0.55), projection = "3d")
ax.patch.set_facecolor("none")

X_plot = [[getindex.(x_list[k], i) for i = 1:3] for k = 1:2]
ax.plot(X_plot[1]..., marker = ".")
ax.plot(X_plot[2]..., marker = ".")
ax.xaxis.pane.fill = false
ax.yaxis.pane.fill = false
ax.zaxis.pane.fill = false
ax.set_xlabel(L"$x_1$")
ax.set_ylabel(L"$x_2$")
ax.set_zlabel(L"$x_3$")
ax.view_init(elev = 7.5, azim = -127.25)

println(collect(ccosx_list))
xx = x_list[1]
i1 = 52
xx1 = xx[i1]
println(cos(xx1[1]))
i2 = 66
xx2 = xx[i2]
println(cos(xx2[1]))

np = 50
rad = 3.0
verts_side, ~ = Plot.matrix_to_cone3d(-P_opt[2], rad, np)
f_shift(x) = x + xx1
map!(verts -> f_shift.(verts), verts_side, verts_side)
poly_list = Plot.make_collection(verts_side, fc = "green", fa = 0.3, ec = "none")
ax.add_collection3d(poly_list)
verts_side, ~ = Plot.matrix_to_cone3d(-P_opt[2], rad, np)
f_shift(x) = x + xx2
map!(verts -> f_shift.(verts), verts_side, verts_side)
poly_list = Plot.make_collection(verts_side, fc = "green", fa = 0.3, ec = "none")
ax.add_collection3d(poly_list)

## Cone plot
x0 = 0.02
y0 = 0.04
ax1 = fig.add_axes((x0 + 0.0, y0 + 0.20, 0.22, 0.24), projection = "3d")
ax2 = fig.add_axes((x0 + 0.22, y0 + 0.20, 0.22, 0.24), projection = "3d")
ax3 = fig.add_axes((x0 + 0.0, y0 + 0.0, 0.22, 0.24), projection = "3d")
ax4 = fig.add_axes((x0 + 0.22, y0 + 0.0, 0.22, 0.24), projection = "3d")
AX_ = (ax1, ax2, ax3, ax4)

np = 50
rad = 0.5
circ = Vector{Any}(undef, 4)
LS = 0.6
LAB = [L"$P_{\mathtt{a}}$", L"$P_{\mathtt{b}}$", L"$P_{\mathtt{c}}$", L"$P_{\mathtt{d}}$"]

for i = 1:4
    ax = AX_[i]
    ax.patch.set_facecolor("none")
    verts_side, ~ = Plot.matrix_to_cone3d(-P_opt[i], rad, np)
    poly_list = Plot.make_collection(verts_side, fc = Plot._colors[i], fa = 0.3, ec = "none")
    ax.add_collection3d(poly_list)
    circ[i] = matplotlib.patches.Patch(fc = Plot._colors[i], ec = "none", alpha = 0.3)
    ax.quiver3D(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, arrow_length_ratio = 0.3, length = 0.5)
    ax.quiver3D(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, arrow_length_ratio = 0.3, length = 0.5)
    ax.quiver3D(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, arrow_length_ratio = 0.3, length = 0.5)
    ax.set_xlim(-LS, LS)
    ax.set_ylim(-LS, LS)
    ax.set_zlim(-LS, LS)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    ax.text3D(-0.2*LS, -1.2*LS, -1.3*LS, L"$x_1$")
    ax.text3D(-1.1*LS, 0.5*LS, -1.4*LS, L"$x_2$")
    ax.text3D(-1.5*LS, 1.1*LS, 0.0, L"$x_3$")
    ax.xaxis.pane.fill = false
    ax.yaxis.pane.fill = false
    ax.zaxis.pane.fill = false
    ax.view_init(elev = 7.5, azim = -127.25)
    # ax.axis("off")
    ax.text3D(-0.6*LS, LS, 0.6*LS, LAB[i], color = Plot._colors[i],
        bbox = Dict([("fc", "white"), ("alpha", 0.5)]),
        fontsize = 15)
end

ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
ax.patch.set_facecolor("none")
ax.text(0.04, 0.47, "c", weight = "bold", fontsize = 15)
ax.text(0.52, 0.47, "d", weight = "bold", fontsize = 15)
ax.text(0.52, 0.9, "b", weight = "bold", fontsize = 15)
ax.text(0.04, 0.9, "a", weight = "bold", fontsize = 15)

# fig.savefig("test/figures/fig_NonLin_2dom_duffing_allllllll.png",
#     transparent = false,
#     bbox_inches = matplotlib.transforms.Bbox([[0.3, 0.3], [8.9, 7.7]]))
end
