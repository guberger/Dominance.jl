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

matplotlib.rc("legend", fontsize = 15)
matplotlib.rc("axes", labelsize = 15)
matplotlib.rc("xtick", labelsize = 11)
matplotlib.rc("ytick", labelsize = 11)

α = 0.2 # limit: (6.0 - sqrt(32.0))/2.0 = 0.1715728752538097
A = [α α-1.0; 0.0 1.0]

fig = PyPlot.figure(figsize = (9.8, 4.8))
gs = matplotlib.gridspec.GridSpec(1, 2, figure = fig, wspace = 0.35)
ax = fig.add_subplot(get(gs, 0), aspect = "equal")

np = 50
t = range(0.0, 2.0*pi, length = np)

ev = eigvals(A)
γ = 0.6

ax.plot(γ*cos.(t), γ*sin.(t), ls = "--", c = "green", lw = 3.0,
    label = L"$\gamma = 0.6$")
ax.plot(real.(ev), imag.(ev), marker = ".", ms = 18, ls = "none",
    mfc = _cols[1], mec = "black", mew = 1.5)

ax.legend(ncol = 2)
ax.grid(true)
ax.set_xlabel("Re(λ)")
ax.set_ylabel("Im(λ)")
ax.set_xlim(-1.2, 1.2)
ax.set_ylim(-1.2, 1.2)
ax.text(α-0.05, -0.4, "α", fontsize = 15)
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

ax = fig.add_subplot(get(gs, 1), aspect = "equal")

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
ax.legend(circ, LAB)

zlevel = 0.79
fig.text(0.05, zlevel, "a", weight = "bold", fontsize = 18)
fig.text(0.5, zlevel, "b", weight = "bold", fontsize = 18)

fig.savefig("./figures/fig_LTI_1dom_simple_eigs_cone.png", transparent = false,
    bbox_inches = matplotlib.transforms.Bbox([[0.47, 0.24], [8.86, 4.0]]))

end
