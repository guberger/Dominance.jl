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

Ga = [0.1 0.2 0.2; 0.95 0 0; 0 0.9 0]
Gb = [0.3 0.9 0.7; 0.9 0 0; 0 0.85 0]
alpha = 0.3
B1 = [1 alpha; 0 1-alpha]
B2 = B1[[2, 1], [2, 1]]
A1 = [B1[1, 1] * Ga B1[1, 2] * Gb;
    B1[2, 1] * Ga B1[2, 2] * Gb]
setindex!(A1, zeros(3), 1, 4:6)
A2 = [B2[1, 1] * Ga B2[1, 2] * Gb;
    B2[2, 1] * Ga B2[2, 2] * Gb]
setindex!(A2, zeros(3), 4, 1:3)
A_list = [A1, A2]
EG1 = pth_eigen(A1, 1, 1e-9)
EG2 = pth_eigen(A2, 1, 1e-9)
display([EG1, EG2])
γ1 = sqrt(prod(norm.(EG1)))
γ2 = sqrt(prod(norm.(EG2)))
display([γ1, γ2])

ev1 = eigvals(A1)
ev2 = eigvals(A2)

fig = PyPlot.figure(figsize = [6.0, 4.8])
ax = fig.gca(aspect = "equal")
np = 50
t = range(0.0, 2.0*pi, length = np)

ax.plot(γ1*cos.(t), γ1*sin.(t), ls = "--", c = _cols[1], lw = 3.0,
    label = L"$\gamma_{\mathtt{a}1\mathtt{a}}$")
ax.plot(γ2*cos.(t), γ2*sin.(t), ls = "--", c = _cols[2], lw = 3.0,
    label = L"$\gamma_{\mathtt{a}2\mathtt{a}}$")
ax.plot(0.0, 0.0, marker = "x", c = "k", ms = 8, mew = 3.0)

ax.plot(real.(ev1), imag.(ev1), marker = ".", ms = 18, ls = "none",
    mfc = _cols[1], mec = "black", mew = 1.5, label = L"$A_1$")
ax.plot(real.(ev2), imag.(ev2), marker = ".", ms = 18, ls = "none",
    mfc = _cols[2], mec = "black", mew = 1.5, label = L"$A_2$")

Leg = ax.legend(fontsize = 14, loc = "upper left", bbox_to_anchor = (1.01, 1.0))
ax.grid(true)
ax.set_xlabel("Re(λ)")
ax.set_ylabel("Im(λ)")
ax.set_xlim(-1.2, 1.5)
ax.set_ylim(-1.3, 1.3)
fig.savefig("test/figures/fig_SLS_1dom_Leslie_eigs.png",
    transparent = false, bbox_inches = "tight")

## Cones

edge_list = [(1, 1, 1, 1),
    (1, 1, 2, 2,)]
rates_list = ndgrid_array([[γ1], [sqrt(γ2)]])
display(rates_list)

println("")
P_max, ee_max, rates_max = solve_lmi_disc_path(A_list, edge_list, rates_list)
println("")

for i = 1:length(P_max)
    display(P_max[i])
    display(eigvals(P_max[i]))
end

@printf("\n")
display(ee_max)
display(rates_max)
@printf("\n")

## Asymptotic behavior

Random.seed!(1)
np = 10
tmax = 15
x_list = Vector{Vector{Vector{Float64}}}(undef, tmax)
x_list[1] = [rand(6) for n = 1: np]
fnorm = x -> x/norm(x, 1)
map!(fnorm, x_list[1], x_list[1])
seq = rand([1, length(A_list)], tmax)
display(seq)

fig = PyPlot.figure(figsize = [9.8, 4.8])
ax = fig.gca()

for t = 1:tmax-1
    x_list[t+1] = map(x -> fnorm(A_list[seq[t]]*x), x_list[t])
end

X_plot = [[x_list[t][k][i] for t = 1:tmax] for k = 1:np, i = 1:6]

for i = 1:6, k = 1:np
    ax.plot(X_plot[k,i], marker = ".", ms = 10, c = _cols[i])
end

circ = Vector{Any}(undef, 6)

for i = 1:6
    circ[i] = matplotlib.lines.Line2D([0.0], [0.0], ls = "solid", c = _cols[i],
        marker = ".", ms = 10, mfc = _cols[i], mec = _cols[i])
end

LAB = [L"$x_1$", L"$x_2$", L"$x_3$", L"$x_4$", L"$x_5$", L"$x_6$",]
ax.legend(circ, LAB, ncol = 6, fontsize = 12, loc = "upper right")

ax.set_xlabel(L"$t$")
ax.set_ylabel(L"$x_i(t)$")
ax.set_xlim(-0.3, tmax - 0.7)
ax.set_ylim(-0.02, 0.55)

fig.savefig("test/figures/fig_SLS_1dom_Leslie_traj.png",
    transparent = false, bbox_inches = "tight")

end
