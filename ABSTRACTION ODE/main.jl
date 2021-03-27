# Main dominance abstraction

module TestMain

include("abstraction.jl")
MABS = Abstraction
include("system_models.jl")
using LinearAlgebra
using Printf
using PyPlot
using PyCall
const art3d = PyObject(PyPlot.art3D)

println("Started main_2")

AbsGenerator = newDuffingControlled3D()

theAbs = AbsGenerator()

# (DDxNorm, xopt) = MABS.DDxNorm_vecInf_2_matp_global(theAbs, Inf,
#     [-5.0, -3.0], [10.0, 6.0], [1, 1])
# @printf("DF Lipschitz Inf: %f\n", DDxNorm)

nsubopt = fill(1, theAbs.dim)
theAbs.nsubopt_list = fill(nsubopt, length(theAbs.cell_list))
DS_list = MABS.compute_dslipInf!(theAbs, [-1])
nsubim_list = MABS.compute_nsubim!(theAbs, [-1], 0.5, Inf)
nsuboptnew = fill(1, theAbs.dim)
MABS.extend_invariant_sharp!(theAbs, [-1], 200, nsuboptnew, 0.5, Inf)
# MABS.extend_invariant_lazy!(theAbs, [-1], 800, 15.0, 0.1, Inf)

# theAbs.dslip2_list = sqrt(theAbs.dim)*theAbs.dslipInf_list

nloop = 1

for i = 1:nloop
    # MABS.refine_abstraction!(theAbs, [-1], fill(2, theAbs.dim))
    # MABS.compute_edges!(theAbs, [-1])
    idxb_list = MABS.find_blockingcells(theAbs)
    MABS.remove_cells!(theAbs, idxb_list)
end
# =#

# theAbs = Main.theAbs0

println("Start plotting")

fig = PyPlot.figure() # figsize = (20, 10)
ax = MABS.set_axes(fig, theAbs)

#=
h = theAbs.h
x = theAbs.x_or + theAbs.cell_list[1].coords.*h
sub_iter = Iterators.product([0, 1], [0, 1])

cols = [:red, :yellow, :green, :blue]

for (i, it) in enumerate(sub_iter)
    x0 = x + h.*it
    nstep = 10
    MABS.plottrajectory!(ax, theAbs, x0, nstep, markersize = 12, linecolor = cols[i])
end

x0 = x + 0.5*h
A = theAbs.linsys_map(x0)[2]
display(A)
=#

MABS.plottrajectory!(ax, theAbs, theAbs.x_or, 10, markersize = 12)

scal, fc_list = MABS.vector_2_colormap(theAbs.dslipInf_list, "inferno")
MABS.colorbar_from_scal(fig, scal)

MABS.plotcells!(ax, theAbs, [-1], facecolor_list = fc_list, facealpha_list = [0.3],
    edgewidth = 0.5, edgealpha = 0.25, show_numbers = false)

MABS.plotimages!(ax, theAbs, [-1], nsub_list = [fill(8, theAbs.dim)],
    facealpha_list = [0.2], facecolor_list = ["blue"],
    edgecolor = "black", edgewidth = 0.5, edgealpha = 0.5,
    wirecolor_list = ["darkblue"], wirewidth = 0.1, wirealpha_list = [0.1])

nidx = length(theAbs.cell_list)
MABS.plotapprox!(ax, theAbs, [-1], faceoutcolor_list = ["green"],
    faceoutalpha_list = [0.1], edgeoutalpha_list = [0.2])
ax.set_xlim([-4.0, 4.0])
ax.set_ylim([-2.3, 2.3])
if theAbs.dim == 3
    ax.set_zlim([-2.8, 2.8])
end
ax.set_aspect("equal")

#=
idx_list = [-1]
TWO_ = fill(2, theAbs.dim)

while ~isempty(idx_list)
    global idx_list
    MABS.compute_nsubds!(theAbs, idx_list, 0.5)
    idx_list = findall(x -> any(x .> 1), theAbs.nsubds_list)
    display(idx_list)
    a = readline()
    @assert a != "5"
    idxnew_list = MABS.refine_abstraction!(theAbs, idx_list, TWO_)[1]
    append!(idx_list, idxnew_list)
end

MABS.compute_edges!(theAbs, [-1])
idxb_list = MABS.find_blockingcells(theAbs)
MABS.remove_cells!(theAbs, idxb_list)
MABS.compute_DSlist!(theAbs, [-1])
display(maximum(x -> maximum(x), theAbs.nsubds_list))
display(sum(length, theAbs.DS_list))
display(sum(x -> length(theAbs.DS_list[x[1]]), theAbs.edge_list) + 2*length(theAbs.cell_list))
# MABS.compute_cones(theAbs, 0.9)
=#

println("Finished")

end # module TestMain
