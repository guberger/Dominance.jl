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

## Cones

α = 0.1 # limit: (6.0 - sqrt(32.0))/2.0 = 0.1715728752538097
A1 = @SMatrix [1.0 0.0; 1.0-α α]
A2 = @SMatrix [α α-1.0; 0.0 1.0]

ev1 = eigvals(Matrix(A1))
γ1 = sqrt(prod(norm.(ev1)))
ev2 = eigvals(Matrix(A1*A2))
γ2 = sqrt(prod(norm.(ev2)))

graph = DO.Graph(2)
DO.add_edge!(graph, 1, 1)
DO.add_edge!(graph, 1, 2)
DO.add_edge!(graph, 2, 1)
DO.add_edge!(graph, 2, 2)
ASri_lab = Dict([DO.Edge(1, 1) => [(DO.MatrixSet(A1), 1)],
    DO.Edge(1, 2) => [(DO.MatrixSet(A2), 2)],
    DO.Edge(2, 1) => [(DO.MatrixSet(A1), 2)],
    DO.Edge(2, 2) => [(DO.MatrixSet(A2), 1)]])
rate_tuple_iter = Iterators.product((γ1,), (sqrt(γ2),))

optim_solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
P_opt, ~, ~ = DO.cone_optim(graph, ASri_lab, rate_tuple_iter, optim_solver)

np = 50
rad = 1.0
fact = 1.1

fig = PyPlot.figure(figsize = (6.4, 4.8))
ax = fig.add_subplot(aspect = "equal")

verts = Plot.matrix_to_cone2d(P_opt[1], rad, np)
poly_list = Plot.make_collection(verts, fc = Plot._colors[1], fa = 0.75, ec = Plot._colors[1])
ax.add_collection(poly_list)
verts = Plot.matrix_to_cone2d(P_opt[2], rad, np)
poly_list = Plot.make_collection(verts, fc = Plot._colors[2], fa = 0.75, ec = Plot._colors[2])
ax.add_collection(poly_list)
verts = Plot.matrix_to_cone2d(A1'\P_opt[1]/A1, rad*fact, np)
poly_list = Plot.make_collection(verts, fc = Plot._colors[3], fa = 0.75, ec = Plot._colors[3])
ax.add_collection(poly_list)
verts = Plot.matrix_to_cone2d(A1'\P_opt[2]/A1, rad*fact, np)
poly_list = Plot.make_collection(verts, fc = Plot._colors[4], fa = 0.75, ec = Plot._colors[4])
ax.add_collection(poly_list)
verts = Plot.matrix_to_cone2d(A2'\P_opt[1]/A2, rad*fact, np)
poly_list = Plot.make_collection(verts, fc = Plot._colors[5], fa = 0.75, ec = Plot._colors[5])
ax.add_collection(poly_list)
verts = Plot.matrix_to_cone2d(A2'\P_opt[2]/A2, rad*fact, np)
poly_list = Plot.make_collection(verts, fc = Plot._colors[6], fa = 0.75, ec = Plot._colors[6])
ax.add_collection(poly_list)

ax.set_xlim(-1.15, 1.15)
ax.set_ylim(-1.15, 1.15)
ax.set_yticks(-1.0:0.5:1.0)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.grid(true)

circ = [matplotlib.patches.Patch(
    fc = Plot._colors[i], ec = Plot._colors[i], alpha = 0.8) for i = 1:6]
LAB = (L"$\mathcal{K}(P_{\mathtt{a}})$", L"$\mathcal{K}(P_{\mathtt{b}})$",
    L"$A_1(\mathcal{K}(P_{\mathtt{a}}))$", L"$A_1(\mathcal{K}(P_{\mathtt{b}}))$",
    L"$A_2(\mathcal{K}(P_{\mathtt{a}}))$", L"$A_2(\mathcal{K}(P_{\mathtt{b}}))$")
ax.legend(circ, LAB, loc = "upper left", bbox_to_anchor = (1.0, 1.0))

fig.savefig("./figures/fig_SLS_1dom_simple_cones.png", transparent = false,
    bbox_inches = matplotlib.transforms.Bbox(((0.72, 0.05), (6.99, 4.24))))

## Trajectories

Random.seed!(5)
A_list = (A1, A2)
np = 40
tmax = 15
x_list = Vector{Vector{SVector{2,Float64}}}(undef, tmax)
x_list[1] = [SVector{2}(rand() - 0.5 for i = 1:2) for n = 1:np]
fnorm = x -> x/norm(x)
map!(fnorm, x_list[1], x_list[1])
seq = rand((1, length(A_list)), tmax)

fig = PyPlot.figure(figsize = (11.3, 4.0))
ax = fig.add_subplot()

for t = 1:tmax-1
    x_list[t+1] = map(x -> fnorm(A_list[seq[t]]*x), x_list[t])
end

X_plot = [[x_list[t][k][i] for t = 1:tmax] for k = 1:np, i = 1:2]

for i = 1:2, k = 1:np
    ax.plot(X_plot[k,i], marker = ".", ms = 10, c = Plot._colors[i])
end

circ = [matplotlib.lines.Line2D(
    (0.0,), (0.0,), ls = "solid", c = Plot._colors[i], marker = ".", ms = 10,
    mfc = Plot._colors[i], mec = Plot._colors[i]) for i = 1:2]
LAB = (L"$\xi^{(1)}$", L"$\xi^{(2)}$")
ax.legend(circ, LAB, ncol = 2)

ax.set_xlabel(L"$t$")
ax.set_ylabel(L"$\xi^{(i)}(t)$")
ax.set_xlim(-0.5, tmax - 0.5)
ax.set_ylim(-1.2, 1.3)

fig.savefig("./figures/fig_SLS_1dom_simple_traj.png", transparent = false,
    bbox_inches = matplotlib.transforms.Bbox(((0.6, -0.05), (10.2, 3.55))))

## Fragile cone

fig = PyPlot.figure(figsize = (9.8, 9.8))
gs = matplotlib.gridspec.GridSpec(2, 2, figure = fig, wspace = 0.35, hspace = 0.15)
ax = fig.add_subplot(get(gs, 0), aspect = "equal")

np = 50
t = range(0.0, 2.0*pi, length = np)

ax.plot(γ1*cos.(t), γ1*sin.(t), ls = "--", c = Plot._colors[1], lw = 3.0, label = L"$\bar{\gamma}$")
ax.plot(γ2*cos.(t), γ2*sin.(t), ls = "--", c = Plot._colors[2], lw = 3.0, label = L"$\bar{\gamma}^2$")
ax.plot(real.(ev1), imag.(ev1), clip_on = false, marker = ".", ms = 18,
    ls = "none", mfc = Plot._colors[1], mec = "black", mew = 1.5,
    label = L"$A_1$, $A_2$", zorder = 100)
ax.plot(real.(ev2), imag.(ev2), marker = ".", ms = 18, ls = "none",
    mfc = Plot._colors[2], mec = "black", mew = 1.5, label = L"$A_1A_2$, $A_2A_1$")

ax.legend(ncol = 2, loc = "upper center")
ax.grid(true)
ax.set_xlabel("Re(λ)")
ax.set_ylabel("Im(λ)")
ax.set_xlim(-1.0, 1.0)
ax.set_ylim(-1.0, 1.0)
ax.set_xticks(-1.0:0.5:1.0)
ax.set_yticks(-1.0:0.5:1.0)
ax.text(α-0.05, -0.4, "α", fontsize = 15)
ax.plot((α, α), (-0.05, -0.28), c = "black")
ax.text(-0.17, 1.0, "a", weight = "bold", fontsize = 18,
    ha = "right", va = "center", transform = ax.transAxes)

ax = fig.add_subplot(get(gs, 1), aspect = "equal")
LG = 0.9
ax.plot(LG.*(-1.0, 1.0), LG.*(-1.0, 1.0), ls = "solid", c = "blue", lw = 2.0, label = L"$A_1$")
ax.plot((0.0, 0.0), (-1.0, 1.0), ls = "dashed", c = "blue", lw = 2.0)
ax.plot(LG.*(-1.0, 1.0), -LG.*(-1.0, 1.0), ls = "solid", c = "red", lw = 2.0, label = L"$A_2$")
ax.plot((-1.0, 1.0), (0.0, 0.0), ls = "dashed", c = "red", lw = 2.0)
ax.set_xlim(-1.0, 1.0)
ax.set_ylim(-1.0, 1.0)
ax.text(0.05, 0.65, L"$\lambda=\alpha$", fontsize = 15, rotation = -90)
ax.text(0.48, 0.65, L"$\lambda=1$", fontsize = 15, rotation = 45)
ax.text(0.65, 0.075, L"$\lambda=\alpha$", fontsize = 15)
ax.text(0.61, -0.73, L"$\lambda=1$", fontsize = 15, rotation = -45)
circ1 = matplotlib.lines.Line2D((0.0,), (0.0,), c = "blue", lw = 2.0, label = L"$A_1$")
circ2 = matplotlib.lines.Line2D((0.0,), (0.0,), c = "red", lw = 2.0, label = L"$A_2$")
ax.legend(handles = (circ1, circ2), loc = "lower left", framealpha = 1.0)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_yticks(-1.0:0.5:1.0)
ax.grid(true)
ax.text(-0.17, 1.0, "b", weight = "bold", fontsize = 18,
    ha = "right", va = "center", transform = ax.transAxes)

α = 0.1715 # limit: (6.0 - sqrt(32.0))/2.0 = 0.1715728752538097
A1 = @SMatrix [1.0 0.0; 1.0-α α]
A2 = @SMatrix [α α-1.0; 0.0 1.0]

ev1 = eigvals(Matrix(A1))
γ1 = sqrt(prod(norm.(ev1)))
ev2 = eigvals(Matrix(A1*A2))
γ2 = sqrt(prod(norm.(ev2)))

ax = fig.add_subplot(get(gs, 2), aspect = "equal")

ax.plot(γ1*cos.(t), γ1*sin.(t), ls = "--", c = Plot._colors[1], lw = 3.0, label = L"$\bar{\gamma}$")
ax.plot(γ2*cos.(t), γ2*sin.(t), ls = "--", c = Plot._colors[2], lw = 3.0, label = L"$\bar{\gamma}^2$")
ax.plot(real.(ev1), imag.(ev1), clip_on = false, marker = ".", ms = 18,
    ls = "none", mfc = Plot._colors[1], mec = "black", mew = 1.5,
    label = L"$A_1$, $A_2$", zorder = 100)
ax.plot(real.(ev2), imag.(ev2), marker = ".", ms = 18, ls = "none",
    mfc = Plot._colors[2], mec = "black", mew = 1.5, label = L"$A_1A_2$, $A_2A_1$")

ax.legend(ncol = 2, loc = "upper center")
ax.grid(true)
ax.set_xlabel("Re(λ)")
ax.set_ylabel("Im(λ)")
ax.set_xlim(-1.0, 1.0)
ax.set_ylim(-1.0, 1.0)
ax.set_xticks(-1.0:0.5:1.0)
ax.set_yticks(-1.0:0.5:1.0)
ax.text(α-0.05, -0.25, "α", fontsize = 18)
ax.plot((α, α), (-0.05, -0.13), c = "black")
ax.text(-0.17, 1.0, "c", weight = "bold", fontsize = 18,
    ha = "right", va = "center", transform = ax.transAxes)

rate_tuple_iter = Iterators.product((γ1,), (sqrt(γ2),))

optim_solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
P_opt, ~, ~ = DO.cone_optim(graph, ASri_lab, rate_tuple_iter, optim_solver)

np = 50
rad = 1.0
fact = 1.1

ax = fig.add_subplot(get(gs, 3), aspect = "equal")

verts = Plot.matrix_to_cone2d(P_opt[1], rad, np)
poly_list = Plot.make_collection(verts, fc = Plot._colors[1], fa = 0.75, ec = Plot._colors[1])
ax.add_collection(poly_list)
verts = Plot.matrix_to_cone2d(P_opt[2], rad, np)
poly_list = Plot.make_collection(verts, fc = Plot._colors[2], fa = 0.75, ec = Plot._colors[2])
ax.add_collection(poly_list)
verts = Plot.matrix_to_cone2d(A1'\P_opt[1]/A1, rad*fact, np)
poly_list = Plot.make_collection(verts, fc = Plot._colors[3], fa = 0.75, ec = Plot._colors[3])
ax.add_collection(poly_list)
verts = Plot.matrix_to_cone2d(A1'\P_opt[2]/A1, rad*fact, np)
poly_list = Plot.make_collection(verts, fc = Plot._colors[4], fa = 0.75, ec = Plot._colors[4])
ax.add_collection(poly_list)
verts = Plot.matrix_to_cone2d(A2'\P_opt[1]/A2, rad*fact, np)
poly_list = Plot.make_collection(verts, fc = Plot._colors[5], fa = 0.75, ec = Plot._colors[5])
ax.add_collection(poly_list)
verts = Plot.matrix_to_cone2d(A2'\P_opt[2]/A2, rad*fact, np)
poly_list = Plot.make_collection(verts, fc = Plot._colors[6], fa = 0.75, ec = Plot._colors[6])
ax.add_collection(poly_list)

ax.set_xlim(-1.15, 1.15)
ax.set_ylim(-1.15, 1.15)
ax.set_yticks(-1.0:0.5:1.0)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.grid(true)
ax.text(-0.17, 1.0, "d", weight = "bold", fontsize = 18,
    ha = "right", va = "center", transform = ax.transAxes)
ax.text(1.0, -0.03, "1.0", ha = "center", va = "center", c = "w",
    transform = ax.transAxes)

fig.savefig("./figures/fig_SLS_1dom_simple_fragile_eigs_cones.png", transparent = false,
    bbox_inches = matplotlib.transforms.Bbox(((0.5, 0.7), (8.95, 8.6))))
end  # module ExampleMain
