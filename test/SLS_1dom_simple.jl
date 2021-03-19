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

α = 0.1 # limit: (6.0 - sqrt(32.0))/2.0 = 0.1715728752538097
A1 = [1.0 0.0; 1.0-α α]
A2 = [α α-1.0; 0.0 1.0]

ev1 = eigvals(A1)
γ1 = sqrt(prod(norm.(ev1)))
ev2 = eigvals(A1*A2)
γ2 = sqrt(prod(norm.(ev2)))
display([γ1, γ2])

## Eigenspaces

fig = PyPlot.figure(figsize = [9.8, 4.8])
gs = matplotlib.gridspec.GridSpec(1, 2, figure = fig, wspace = 0.35)
ax = fig.add_subplot(get(gs, 0), aspect = "equal")
np = 50
t = range(0.0, 2.0*pi, length = np)

ax.plot(γ1*cos.(t), γ1*sin.(t), ls = "--", c = _cols[1], lw = 3.0, label = L"$\bar{\gamma}$")
ax.plot(γ2*cos.(t), γ2*sin.(t), ls = "--", c = _cols[2], lw = 3.0, label = L"$\bar{\gamma}^2$")

ax.plot(real.(ev1), imag.(ev1), clip_on = false, marker = ".", ms = 18, ls = "none",
    mfc = _cols[1], mec = "black", mew = 1.5, label = L"$A_1$, $A_2$", zorder = 100)
ax.plot(real.(ev2), imag.(ev2), marker = ".", ms = 18, ls = "none",
    mfc = _cols[2], mec = "black", mew = 1.5, label = L"$A_1A_2$, $A_2A_1$")

ax.legend(fontsize = 14, ncol = 2, loc = "upper center")
ax.grid(true)
ax.set_xlabel("Re(λ)")
ax.set_ylabel("Im(λ)")
ax.set_xlim(-1.0, 1.0)
ax.set_ylim(-1.0, 1.0)
ax.set_xticks(-1.0:0.5:1.0)
ax.set_yticks(-1.0:0.5:1.0)
ax.text(α-0.05, -0.4, "α", fontsize = 14)
ax.plot([α, α], [-0.05, -0.28], c = "black")
ax.text(-0.15, 1.0, "a", weight = "bold", fontsize = 14,
    ha = "right", va = "center", transform = ax.transAxes)

ax = fig.add_subplot(get(gs, 1), aspect = "equal")
LG = 0.9
ax.plot(LG*[-1.0, 1.0], LG*[-1.0, 1.0], ls = "solid", c = "blue", lw = 2.0, label = L"$A_1$")
ax.plot([0.0, 0.0], [-1.0, 1.0], ls = "dashed", c = "blue", lw = 2.0)
ax.plot(LG*[-1.0, 1.0], -LG*[-1.0, 1.0], ls = "solid", c = "red", lw = 2.0, label = L"$A_2$")
ax.plot([-1.0, 1.0], [0.0, 0.0], ls = "dashed", c = "red", lw = 2.0)
ax.set_xlim(-1.0, 1.0)
ax.set_ylim(-1.0, 1.0)
ax.text(0.05, 0.7, L"$\lambda=\alpha$", fontsize=12, rotation = -90)
ax.text(0.5, 0.65, L"$\lambda=1$", fontsize=12, rotation = 45)
ax.text(0.7, 0.05, L"$\lambda=\alpha$", fontsize=12)
ax.text(0.61, -0.73, L"$\lambda=1$", fontsize=12, rotation = -45)
circ1 = matplotlib.lines.Line2D([0.0], [0.0], c = "blue", lw = 2.0, label = L"$A_1$")
circ2 = matplotlib.lines.Line2D([0.0], [0.0], c = "red", lw = 2.0, label = L"$A_2$")
ax.legend(handles = [circ1, circ2], fontsize = 14, loc = "lower left", framealpha = 1.0)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_yticks(-1.0:0.5:1.0)
ax.grid(true)
ax.text(-0.15, 1.0, "b", weight = "bold", fontsize = 14,
    ha = "right", va = "center", transform = ax.transAxes)

fig.savefig("test/figures/fig_SLS_1dom_simple_eigs.png", transparent = false,
    bbox_inches = "tight")

## Cones

fig = PyPlot.figure(figsize = [6.4, 4.8])
ax = fig.gca(aspect = "equal")

A_list = [A1, A2]
edge_list = [(1, 1, 1, 1),
    (1, 2, 2, 2),
    (2, 1, 1, 2),
    (2, 2, 2, 1)]
rates_list = ndgrid_array([[γ1], [sqrt(γ2)]])
display(rates_list)

println("")
P_max, ee_max, rates_max = solve_lmi_disc_path(A_list, edge_list, rates_list)
println("")

np = 50
rad = 1.0
fact = 1.1

verts = matrix_to_cone2d(P_max[1], rad, np)
pcoll = make_collection(verts, facecolor = _cols[1], facealpha = 0.75,
    edgecolor = _cols[1])
ax.add_collection(pcoll)
verts = matrix_to_cone2d(P_max[2], rad, np)
pcoll = make_collection(verts, facecolor = _cols[2], facealpha = 0.75,
    edgecolor = _cols[2])
ax.add_collection(pcoll)

verts = matrix_to_cone2d(A1'\P_max[1]/A1, rad*fact, np)
pcoll = make_collection(verts, facecolor = _cols[3], facealpha = 0.75,
    edgecolor = _cols[3])
ax.add_collection(pcoll)
verts = matrix_to_cone2d(A1'\P_max[2]/A1, rad*fact, np)
pcoll = make_collection(verts, facecolor = _cols[4], facealpha = 0.75,
    edgecolor = _cols[4])
ax.add_collection(pcoll)
verts = matrix_to_cone2d(A2'\P_max[1]/A2, rad*fact, np)
pcoll = make_collection(verts, facecolor = _cols[5], facealpha = 0.75,
    edgecolor = _cols[5])
ax.add_collection(pcoll)
verts = matrix_to_cone2d(A2'\P_max[2]/A2, rad*fact, np)
pcoll = make_collection(verts, facecolor = _cols[6], facealpha = 0.75,
    edgecolor = _cols[6])
ax.add_collection(pcoll)
ax.set_xlim(-1.15, 1.15)
ax.set_ylim(-1.15, 1.15)
ax.set_yticks(-1.0:0.5:1.0)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.grid(true)

circ = Vector{Any}(undef, 6)

for i = 1:6
    circ[i] = matplotlib.patches.Patch(fc = _cols[i], ec = _cols[i], alpha = 0.7)
end

LAB = [L"$\mathcal{K}(P_{\mathtt{a}})$", L"$\mathcal{K}(P_{\mathtt{b}})$",
    L"$A_1(\mathcal{K}(P_{\mathtt{a}}))$", L"$A_1(\mathcal{K}(P_{\mathtt{b}}))$",
    L"$A_2(\mathcal{K}(P_{\mathtt{a}}))$", L"$A_2(\mathcal{K}(P_{\mathtt{b}}))$"]
ax.legend(circ, LAB, loc = "upper left", bbox_to_anchor = (1.0, 1.0), fontsize = 14)

fig.savefig("test/figures/fig_SLS_1dom_simple_cones.png",
    transparent = false, bbox_inches = "tight")

## Trajectories

Random.seed!(1)
np = 40
tmax = 15
x_list = Vector{Vector{Vector{Float64}}}(undef, tmax)
x_list[1] = [rand(2) .- 0.5 for n = 1: np]
fnorm = x -> x/norm(x)
map!(fnorm, x_list[1], x_list[1])
seq = rand([1, length(A_list)], tmax)
display(seq)

fig = PyPlot.figure(figsize = [11.3, 4.0])
ax = fig.gca()

for t = 1:tmax-1
    x_list[t+1] = map(x -> fnorm(A_list[seq[t]]*x), x_list[t])
end

X_plot = [[x_list[t][k][i] for t = 1:tmax] for k = 1:np, i = 1:2]

for i = 1:2, k = 1:np
    ax.plot(X_plot[k,i], marker = ".", ms = 10, c = _cols[i])
end

circ = Vector{Any}(undef, 2)

for i = 1:2
    circ[i] = matplotlib.lines.Line2D([0.0], [0.0], ls = "solid", c = _cols[i],
        marker = ".", ms = 10, mfc = _cols[i], mec = _cols[i])
end

LAB = [L"$x_1$", L"$x_2$"]
ax.legend(circ, LAB, ncol = 2, fontsize = 14)

ax.set_xlabel(L"$t$", fontsize = 14)
ax.set_ylabel(L"$x_i(t)$", fontsize = 14)
ax.set_xlim(-0.5, tmax - 0.5)
ax.set_ylim(-1.2, 1.2)

fig.savefig("test/figures/fig_SLS_1dom_simple_traj.png",
    transparent = false, bbox_inches = "tight")

## Frgigile cone
α = 0.1715 # limit: (6.0 - sqrt(32.0))/2.0 = 0.1715728752538097
A1 = [1.0 0.0; 1.0-α α]
A2 = [α α-1.0; 0.0 1.0]

ev1 = eigvals(A1)
γ1 = sqrt(prod(norm.(ev1)))
ev2 = eigvals(A1*A2)
γ2 = sqrt(prod(norm.(ev2)))
display([γ1, γ2])

fig = PyPlot.figure(figsize = [9.8, 4.8])
gs = matplotlib.gridspec.GridSpec(1, 2, figure = fig, wspace = 0.35)
ax = fig.add_subplot(get(gs, 0), aspect = "equal")
np = 50
t = range(0.0, 2.0*pi, length = np)

ax.plot(γ1*cos.(t), γ1*sin.(t), ls = "--", c = _cols[1], lw = 3.0, label = L"$\bar{\gamma}$")
ax.plot(γ2*cos.(t), γ2*sin.(t), ls = "--", c = _cols[2], lw = 3.0, label = L"$\bar{\gamma}^2$")

ax.plot(real.(ev1), imag.(ev1), clip_on = false, marker = ".", ms = 18, ls = "none",
    mfc = _cols[1], mec = "black", mew = 1.5, label = L"$A_1$, $A_2$", zorder = 100)
ax.plot(real.(ev2), imag.(ev2), marker = ".", ms = 18, ls = "none",
    mfc = _cols[2], mec = "black", mew = 1.5, label = L"$A_1A_2$, $A_2A_1$")

ax.legend(fontsize = 14, ncol = 2, loc = "upper center")
ax.grid(true)
ax.set_xlabel("Re(λ)")
ax.set_ylabel("Im(λ)")
ax.set_xlim(-1.0, 1.0)
ax.set_ylim(-1.0, 1.0)
ax.set_xticks(-1.0:0.5:1.0)
ax.set_yticks(-1.0:0.5:1.0)
ax.text(α-0.05, -0.25, "α", fontsize = 14)
ax.plot([α, α], [-0.05, -0.13], c = "black")
ax.text(-0.15, 1.0, "c", weight = "bold", fontsize = 14,
    ha = "right", va = "center", transform = ax.transAxes)

A_list = [A1, A2]
edge_list = [(1, 1, 1, 1),
    (1, 2, 2, 2),
    (2, 1, 1, 2),
    (2, 2, 2, 1)]
rates_list = ndgrid_array([[γ1], [sqrt(γ2)]])
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
verts = matrix_to_cone2d(P_max[2], rad, np)
pcoll = make_collection(verts, facecolor = _cols[2], facealpha = 0.75,
    edgecolor = _cols[2])
ax.add_collection(pcoll)

verts = matrix_to_cone2d(A1'\P_max[1]/A1, rad*fact, np)
pcoll = make_collection(verts, facecolor = _cols[3], facealpha = 0.75,
    edgecolor = _cols[3])
ax.add_collection(pcoll)
verts = matrix_to_cone2d(A1'\P_max[2]/A1, rad*fact, np)
pcoll = make_collection(verts, facecolor = _cols[4], facealpha = 0.75,
    edgecolor = _cols[4])
ax.add_collection(pcoll)
verts = matrix_to_cone2d(A2'\P_max[1]/A2, rad*fact, np)
pcoll = make_collection(verts, facecolor = _cols[5], facealpha = 0.75,
    edgecolor = _cols[5])
ax.add_collection(pcoll)
verts = matrix_to_cone2d(A2'\P_max[2]/A2, rad*fact, np)
pcoll = make_collection(verts, facecolor = _cols[6], facealpha = 0.75,
    edgecolor = _cols[6])
ax.add_collection(pcoll)
ax.set_xlim(-1.15, 1.15)
ax.set_ylim(-1.15, 1.15)
ax.set_yticks(-1.0:0.5:1.0)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.grid(true)
ax.text(-0.15, 1.0, "d", weight = "bold", fontsize = 14,
    ha = "right", va = "center", transform = ax.transAxes)
ax.text(1.0, -0.03, "1.0", ha = "center", va = "center", c = "w",
    transform = ax.transAxes)

fig.savefig("test/figures/fig_SLS_1dom_simple_fragile_eigs_cones.png",
    transparent = false, bbox_inches = "tight")

end
