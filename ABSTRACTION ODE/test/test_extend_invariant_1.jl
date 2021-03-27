include("../abstraction.jl")
MABS = Abstraction
include("../system_models.jl")
using LinearAlgebra
using Printf
using PyPlot
using PyCall
const art3d = PyObject(PyPlot.art3D)

sleep(0.1) # used for good printing
println("Started test")

AbsGenerator = newNonlinearOscillator2D(5, 5)

theAbs = AbsGenerator()

nsubopt = fill(1, theAbs.dim)
theAbs.nsubopt_list = fill(nsubopt, length(theAbs.cell_list))
DS_list = MABS.compute_dslipInf!(theAbs, [-1])
nsubim_list = MABS.compute_nsubim!(theAbs, [-1], 0.5, Inf)

################################################################################
println("Started plotting")

fig = PyPlot.figure() # figsize = (20, 10)
ax = MABS.set_axes(fig, theAbs)

MABS.plottrajectory!(ax, theAbs, theAbs.x_or, 10, markersize = 12)

scal, fc_list = MABS.vector_2_colormap(theAbs.dslipInf_list, "inferno")
MABS.colorbar_from_scal(fig, scal)

MABS.plotcells!(ax, theAbs, [-1], facecolor_list = fc_list, facealpha_list = [0.3],
    edgewidth = 0.5, edgealpha = 0.25, show_numbers = false)

MABS.plotimages!(ax, theAbs, [-1], nsub_list = [fill(8, theAbs.dim)],
    facealpha_list = [0.2], facecolor_list = ["blue"],
    edgecolor = "black", edgewidth = 0.5, edgealpha = 0.5,
    wirecolor_list = ["darkblue"], wirewidth = 0.1, wirealpha_list = [0.1])

################################################################################
println("Started extend_invariant")

nsuboptnew = fill(1, theAbs.dim)
MABS.extend_invariant_sharp!(theAbs, [-1], 200, nsuboptnew, 0.5, Inf)
# MABS.extend_invariant_lazy!(theAbs, [-1], 800, 15.0, 0.1, Inf)

################################################################################
println("Started plotting")

fig = PyPlot.figure() # figsize = (20, 10)
ax = MABS.set_axes(fig, theAbs)

MABS.plottrajectory!(ax, theAbs, theAbs.x_or, 10, markersize = 12)

scal, fc_list = MABS.vector_2_colormap(theAbs.dslipInf_list, "inferno")
MABS.colorbar_from_scal(fig, scal)

MABS.plotcells!(ax, theAbs, [-1], facecolor_list = fc_list, facealpha_list = [0.3],
    edgewidth = 0.5, edgealpha = 0.25, show_numbers = false)

MABS.plotimages!(ax, theAbs, [-1], nsub_list = [fill(8, theAbs.dim)],
    facealpha_list = [0.2], facecolor_list = ["blue"],
    edgecolor = "black", edgewidth = 0.5, edgealpha = 0.5,
    wirecolor_list = ["darkblue"], wirewidth = 0.1, wirealpha_list = [0.1])

ax.set_xlim([-3.0, 3.0])
ax.set_ylim([-3.0, 3.0])

println("Finished")
