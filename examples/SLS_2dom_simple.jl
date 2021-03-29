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
println("Plot SLS 2-dom simple")

matplotlib.rc("legend", fontsize = 15)
matplotlib.rc("axes", labelsize = 15)
matplotlib.rc("xtick", labelsize = 11)
matplotlib.rc("ytick", labelsize = 11)

alpha_list = -1.0:0.25:0.0
A_list = [(@SMatrix [1 0.5 0; α 0.75 0.5; -0.5 0 1.0]) for α in alpha_list]

fig = PyPlot.figure(figsize = (10.05, 5.7))
gs = matplotlib.gridspec.GridSpec(2, 3, figure = fig, wspace = 0.06, hspace = 0.15)

np_s = 20
rad_s = 1.0
u = range(0.0, 2.0*pi, length = np_s)
v = range(0.0, pi, length = np_s)
XS = rad_s*cos.(u)'.*sin.(v)
YS = rad_s*sin.(u)'.*sin.(v)
ZS = rad_s*ones(size(u))'.*cos.(v)

np = 20
rad = 1.0
u = range(0.0, 2.0*pi, length = np)
v = range(0.0, pi, length = np)
X = rad*cos.(u)'.*sin.(v)
Y = rad*sin.(u)'.*sin.(v)
Z = rad*ones(size(u))'.*cos.(v)
tmax = 50
x_list = Vector{Vector{SVector{3, Float64}}}(undef, tmax)
x_list[1] = [SVector(X[i], Y[i], Z[i]) for i in eachindex(X)]
Random.seed!(4)
seq = rand((1, length(A_list)), tmax)
seq[5] = 5

for t = 1:tmax-1
    x_list[t+1] = map(x -> normalize(A_list[seq[t]]*x), x_list[t])
end

for i = 1:6
    ax = fig.add_subplot(get(gs, i-1), projection = "3d")
    ax.plot_surface(XS, YS, ZS, color = "yellow", alpha = 0.3)
    x = x_list[i]
    x1 = map(x -> x[1], x)
    x2 = map(x -> x[2], x)
    x3 = map(x -> x[3], x)
    ax.plot3D(x1, x2, x3, ls = "none", marker = ".")
    ax.view_init(elev = 20.0, azim = 40.0)
    ax.xaxis.pane.fill = false
    ax.yaxis.pane.fill = false
    ax.zaxis.pane.fill = false
    ax.set_title(string("t = ", i - 1), fontsize = 12)
    ax.tick_params(axis = "both", labelsize = 8)
end

fig.savefig("./figures/fig_SLS_2dom_simple_traj.png", transparent = false,
    bbox_inches = matplotlib.transforms.Bbox(((1.15, 0.5), (9.05, 5.1))))
end  # module ExampleMain
