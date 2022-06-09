module ExampleMain

using LinearAlgebra
using StaticArrays
using PyPlot
using JuMP
using MosekTools
include("../src/Dominance.jl")
DO = Dominance
include("../src/plotting.jl")

sleep(0.1) # used for good printing
println("Plot SLS 1-dom entropy")

matplotlib.rc("legend", fontsize = 15)
matplotlib.rc("axes", labelsize = 15)
matplotlib.rc("xtick", labelsize = 11)
matplotlib.rc("ytick", labelsize = 11)

α = 0.3
A0 = @SMatrix [2 α α; 0 1-α 0; 0 0 1-α]
A_list_tmp = (A0, circshift(A0, (1, 1)), circshift(A0, (2, 2)))
M = @SMatrix [-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
A_list = map(A -> M*A*M, A_list_tmp)

EG1 = DO.pth_eigval(Matrix(A_list[1]), 1, 1e-9)
EG2 = DO.pth_eigval(Matrix(A_list[1]*A_list[2]*A_list[3]), 1, 1e-9)

graph = DO.Graph(3)
ASri_tmp = Any[]
for i = 1:3, j = 1:3
    DO.add_edge!(graph, i, j)
    ri = i == j ? 1 : 2
    push!(ASri_tmp, DO.Edge(i, j) => [(DO.MatrixSet(A_list[j]), ri)])
end
ASri_lab = Dict(ASri_tmp)
nr = 9
rate_tuple_iter = DO.hyper_range((EG1[2], EG2[2]), (EG1[1], EG2[1]), (nr, nr))

optim_solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
P_opt, ~, ~ = DO.cone_optim(graph, ASri_lab, rate_tuple_iter, optim_solver)

fig = PyPlot.figure(figsize = (9.8, 9.0))
gs = matplotlib.gridspec.GridSpec(2, 2, figure = fig, wspace = 0.1, hspace = 0.0)
AX_ = [fig.add_subplot(get(gs, i), projection = "3d") for i = 1:3]

np = 50
rad = 0.5
circ = Vector{Any}(undef, 4)
LS = 0.6
LAB = (L"$P_{\mathtt{a}}$", L"$P_{\mathtt{b}}$", L"$P_{\mathtt{c}}$")

for i = 1:3
    ax = AX_[i]
    verts, ~ = Plot.matrix_to_cone3d(M*P_opt[i]*M, rad, np)
    poly_list = Plot.make_collection(verts, fc = Plot._colors[i], fa = 0.3, ec = "none")
    ax.add_collection3d(poly_list)
    circ[i] = matplotlib.patches.Patch(fc = Plot._colors[i], ec = "none", alpha = 0.3)
    ax.quiver3D(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, arrow_length_ratio = 0.3, length = 0.5)
    ax.quiver3D(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, arrow_length_ratio = 0.3, length = 0.5)
    ax.quiver3D(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, arrow_length_ratio = 0.3, length = 0.5)
    ax.set_xlim(-LS, LS)
    ax.set_ylim(-LS, LS)
    ax.set_zlim(-LS, LS)
    ax.set_xticklabels(())
    ax.set_yticklabels(())
    ax.set_zticklabels(())
    ax.text3D(-0.2*LS, -1.2*LS, -1.4*LS, L"$x$", fontsize = 18)
    ax.text3D(1.1*LS, 0.0, -1.4*LS, L"$y$", fontsize = 18)
    ax.text3D(1.2*LS, 1.1*LS, 0.0, L"$z$", fontsize = 18)
    ax.xaxis.pane.fill = false
    ax.yaxis.pane.fill = false
    ax.zaxis.pane.fill = false
    # ax.axis("off")
    ax.text3D(-0.5*LS, -LS, LS, LAB[i], color = Plot._colors[i],
        bbox = Dict((("fc", "white"), ("alpha", 0.5))), fontsize = 18)
    ax.view_init(elev = 30.0, azim = -55.0)
end

fig.savefig("./figures/fig_SLS_1dom_entropy.png", transparent = false,
    bbox_inches = matplotlib.transforms.Bbox(((1.6, 1.15), (8.78, 7.73))))
end  # module ExampleMain
