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
println("Plot SLS 1-dom Leslie")

matplotlib.rc("legend", fontsize = 15)
matplotlib.rc("axes", labelsize = 15)
matplotlib.rc("xtick", labelsize = 11)
matplotlib.rc("ytick", labelsize = 11)

Ga = @SMatrix [0.1 0.2 0.2; 0.95 0 0; 0 0.9 0]
Gb = @SMatrix [0.3 0.9 0.7; 0.9 0 0; 0 0.85 0]
alpha = 0.3
B1 = @SMatrix [1 alpha; 0 1-alpha]
B2 = B1[SVector(2, 1), SVector(2, 1)]
A1tmp = [B1[1, 1] * Ga B1[1, 2] * Gb;
    B1[2, 1] * Ga B1[2, 2] * Gb]
setindex!(A1tmp, zeros(3), 1, 4:6)
A1 = SMatrix{6,6}(A1tmp)
A2tmp = [B2[1, 1] * Ga B2[1, 2] * Gb;
    B2[2, 1] * Ga B2[2, 2] * Gb]
setindex!(A2tmp, zeros(3), 4, 1:3)
A2 = SMatrix{6,6}(A2tmp)

## Eigenvalues

EG1 = DO.pth_eigval(A1tmp, 1, 1e-9)
EG2 = DO.pth_eigval(A2tmp, 1, 1e-9)
γ1 = sqrt(prod(norm.(EG1)))
γ2 = sqrt(prod(norm.(EG2)))

ev1 = eigvals(A1tmp)
ev2 = eigvals(A2tmp)

fig = PyPlot.figure(figsize = (6.0, 4.8))
ax = fig.gca(aspect = "equal")

np = 50
t = range(0.0, 2.0*pi, length = np)

ax.plot(γ1*cos.(t), γ1*sin.(t), ls = "--", c = Plot._colors[1], lw = 3.0,
    label = L"$\gamma_{\mathtt{a}1\mathtt{a}}$")
ax.plot(γ2*cos.(t), γ2*sin.(t), ls = "--", c = Plot._colors[2], lw = 3.0,
    label = L"$\gamma_{\mathtt{a}2\mathtt{a}}$")
ax.plot(0.0, 0.0, marker = "x", c = "k", ms = 8, mew = 3.0)
ax.plot(real.(ev1), imag.(ev1), marker = ".", ms = 18, ls = "none",
    mfc = Plot._colors[1], mec = "black", mew = 1.5, label = L"$A_1$")
ax.plot(real.(ev2), imag.(ev2), marker = ".", ms = 18, ls = "none",
    mfc = Plot._colors[2], mec = "black", mew = 1.5, label = L"$A_2$")

Leg = ax.legend(loc = "upper left", bbox_to_anchor = (1.01, 1.0))
ax.grid(true)
ax.set_xlabel("Re(λ)")
ax.set_ylabel("Im(λ)")
ax.set_xlim(-1.2, 1.5)
ax.set_ylim(-1.3, 1.3)

## Cones

graph = DO.Graph(1)
DO.add_edge!(graph, 1, 1)
ASri_lab = Dict([DO.Edge(1, 1) => [(DO.MatrixSet(A1), 1), (DO.MatrixSet(A2), 2)]])
rate_tuple_iter = Iterators.product((γ1,), (γ2,))

optim_solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
P_opt, δ_opt, rates_opt = DO.cone_optim(graph, ASri_lab, rate_tuple_iter, optim_solver)

for i in eachindex(P_opt)
    println(eigvals(P_opt[i]))
end
println(δ_opt)
println(rates_opt)

## Trajectories

Random.seed!(1)
A_list = (A1, A2)
np = 10
tmax = 15
x_list = Vector{Vector{SVector{6,Float64}}}(undef, tmax)
x_list[1] = [SVector{6}(rand() for i = 1:6) for n = 1:np]
fnorm = x -> x/norm(x, 1)
map!(fnorm, x_list[1], x_list[1])
seq = rand((1, length(A_list)), tmax)

fig = PyPlot.figure(figsize = (9.8, 4.8))
ax = fig.add_subplot()

for t = 1:tmax-1
    x_list[t+1] = map(x -> fnorm(A_list[seq[t]]*x), x_list[t])
end

X_plot = [[x_list[t][k][i] for t = 1:tmax] for k = 1:np, i = 1:6]

for i = 1:6, k = 1:np
    ax.plot(X_plot[k,i], marker = ".", ms = 10, c = Plot._colors[i])
end

circ = [matplotlib.lines.Line2D(
    (0.0,), (0.0,), ls = "solid", c = Plot._colors[i], marker = ".", ms = 10,
    mfc = Plot._colors[i], mec = Plot._colors[i]) for i = 1:6]
LAB = (L"$\xi^{(1)}$", L"$\xi^{(2)}$", L"$\xi^{(3)}$", L"$\xi^{(4)}$",
    L"$\xi^{(5)}$", L"$\xi^{(6)}$")
ax.legend(circ, LAB, ncol = 6, loc = "upper right", columnspacing = 1.4)

ax.set_xlabel(L"$t$")
ax.set_ylabel(L"$\xi^{(i)}(t)$")
ax.set_xlim(-0.3, tmax - 0.7)
ax.set_ylim(-0.02, 0.55)

fig.savefig("./figures/fig_SLS_1dom_Leslie_traj.png", transparent = false,
    bbox_inches = matplotlib.transforms.Bbox(((0.54, 0.05), (8.85, 4.25))))
end  # module ExampleMain
