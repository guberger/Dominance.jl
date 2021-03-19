module TestMain

include("../macros.jl")
include("../plotting.jl")

using LinearAlgebra
using Printf
using PyPlot
using PyCall
art3d = PyObject(PyPlot.art3D)
_cols = repeat(matplotlib.rcParams["axes.prop_cycle"].by_key()["color"], 10, 1)
CConv = matplotlib.colors.colorConverter
axes_grid1 = pyimport("mpl_toolkits.axes_grid1")

A = [1 0.5 0; -0.5 0.75 0.5; -0.5 0 1.0]

fig = PyPlot.figure(figsize = (9.8, 4.8))
gs = matplotlib.gridspec.GridSpec(1, 2, figure = fig,
    width_ratios = (1, 1.5), wspace = 0.05)
ax = fig.add_subplot(get(gs, 0), aspect = "equal")
np = 50
t = range(0.0, 2.0*pi, length = np)

ev = eigvals(A)
γ = 0.9

ax.plot(γ*cos.(t), γ*sin.(t), ls = "--", c = "green", lw = 3.0,
    label = "\$\\gamma = $(γ)\$")
ax.plot(real.(ev), imag.(ev), marker = ".", ms = 18, ls = "none", mfc = _cols[1],
    mec = "black", mew = 1.5)

ax.legend(fontsize = 14, ncol = 2)
ax.grid(true)
ax.set_xlabel("Re(λ)")
ax.set_ylabel("Im(λ)")
ax.set_xlim(-1.75, 1.75)
ax.set_ylim(-1.75, 1.75)

A_list = [A]
edge_list = [(1, 1, 1, 1)]
rates_list = ndgrid_array([γ])
display(rates_list)

println("")
P_max, ee_max, rates_max = solve_lmi_disc_path(A_list, edge_list, rates_list)
println("")

np = 50
rad = 1.0

ax = fig.add_subplot(get(gs, 1), projection = "3d")

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

fig.savefig("./figures/fig_LTI_2dom_simple_eigs_cone.png",
    transparent = false, bbox_inches = "tight")
end
