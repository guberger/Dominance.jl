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

α = 0.3
A0 = [2 α α; 0 1-α 0; 0 0 1-α]
A_list = [A0, circshift(A0, (1, 1)), circshift(A0, (2, 2))]
M = Matrix{Float64}(I, 3, 3)
M[1, 1] = -1.0
map!(A -> M*A*M, A_list, A_list)
display(A_list)

EG1 = pth_eigen(A_list[1], 1, 1e-9)
EG2 = pth_eigen(A_list[1]*A_list[2]*A_list[3], 1, 1e-9)

the_rate(i, j) = i == j ? 1 : 2
edge_list = [(i, j, j, the_rate(i, j)) for i = 1:3, j = 1:3][:]
nr = 9
rates_list = ndgrid_array([range(EG1[2], EG1[1], length = nr),
    range(EG2[2], EG2[1], length = nr)])

println("")
P_max, ee_max, rates_max = solve_lmi_disc_path(A_list, edge_list, rates_list)
println("")

fig = PyPlot.figure(figsize = (9.8, 9.0))
gs = matplotlib.gridspec.GridSpec(2, 2, figure = fig, wspace = 0.1,
    hspace = 0.0)
AX_ = [fig.add_subplot(get(gs, i), projection = "3d") for i = 1:3]

np = 50
rad = 0.5
circ = Vector{Any}(undef, 4)
LS = 0.6
LAB = [L"$P_{\mathtt{a}}$", L"$P_{\mathtt{b}}$", L"$P_{\mathtt{c}}$"]

for i = 1:3
    ax = AX_[i]
    verts_side, verts_top = matrix_to_cone3d(M*P_max[i]*M, rad, np)
    pcoll = make_collection(verts_side, facecolor = _cols[i], facealpha = 0.3,
        edgecolor = "none")
    ax.add_collection3d(pcoll)
    circ[i] = matplotlib.patches.Patch(fc = _cols[i], ec = "none", alpha = 0.3)
    ax.patch.set_facecolor("none")
    ax.quiver3D(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, arrow_length_ratio = 0.3,
        length = 0.5)
    ax.quiver3D(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, arrow_length_ratio = 0.3,
        length = 0.5)
    ax.quiver3D(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, arrow_length_ratio = 0.3,
        length = 0.5)
    ax.set_xlim(-LS, LS)
    ax.set_ylim(-LS, LS)
    ax.set_zlim(-LS, LS)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    ax.text3D(-0.2*LS, -1.2*LS, -1.4*LS, L"$x$", fontsize = 18)
    ax.text3D(1.1*LS, 0.0, -1.4*LS, L"$y$", fontsize = 18)
    ax.text3D(1.2*LS, 1.1*LS, 0.0, L"$z$", fontsize = 18)
    ax.xaxis.pane.fill = false
    ax.yaxis.pane.fill = false
    ax.zaxis.pane.fill = false
    # ax.axis("off")
    ax.text3D(-0.5*LS, -LS, LS, LAB[i], color = _cols[i],
        bbox = Dict([("facecolor", "white"), ("alpha", 0.5)]),
        fontsize = 18)
    ax.view_init(elev = 30.0, azim = -55.0)
end

fig.savefig("./figures/fig_SLS_1dom_entropy.png", transparent = false,
    bbox_inches = matplotlib.transforms.Bbox([[1.6, 1.15], [8.78, 7.73]]))
end
