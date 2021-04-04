module TestMain

using Test
using LinearAlgebra
using StaticArrays
@static if isdefined(Main, :TestLocal)
    include("../src/Dominance.jl")
else
    using Dominance
end
DO = Dominance

sleep(0.1) # used for good printing
println("Started test")

@testset "System" begin
lb = SVector(0.0, 0.0)
ub = SVector(10.0, 11.0)
x0 = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)
grid = DO.Grid(x0, h)
domain = DO.Domain(grid)
@test DO.get_ncells(domain) == 0
DO.add_set!(domain, DO.HyperRectangle(lb, ub), DO.OUTER)
@test DO.get_ncells(domain) == 77

tstep = 1.0
nsys = 3
F_sys(x) = SVector(0.5, -cos(x[1]))
DF_sys(x) = SMatrix{2,2}(0.0, sin(x[1]), 0.0, 0.0)
bound_DF = 1.0
bound_DDF = 1.0

sys = DO.ContSystemRK4(tstep, F_sys, DF_sys, bound_DF, bound_DDF, nsys)

pos = (1, 1)
x = DO.get_coord_by_pos(grid, pos)
dom1 = DO.Domain(grid)
DO.add_pos!(dom1, pos)
@test DO.get_ncells(dom1) == 1

@static if get(ENV, "CI", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-10.0, 10.0))
    ax.set_ylim((-10.0, 10.0))
    Plot.domain!(ax, 1:2, dom1)
    Plot.trajectory!(ax, 1:2, sys, x, 50)
    Plot.cell_image!(ax, 1:2, dom1, sys)
    Plot.cell_approx!(ax, 1:2, dom1, sys, (10, 10))
end

θ = π/5.0
U = 2*SMatrix{2,2}(cos(θ), -sin(θ), sin(θ), cos(θ))
F_sys(x) = U*SVector(atan(x[1]), atan(x[2]))
DF_sys(x) = U*SMatrix{2,2}(1/(1 + x[1]^2), 0, 0, 1/(1 + x[2]^2))
bound_DDF = opnorm(U, Inf)*3*sqrt(3)/8

sys = DO.DiscSystem(F_sys, DF_sys, bound_DDF)

pos = (1, 0)
x = DO.get_coord_by_pos(grid, pos)
dom1 = DO.Domain(grid)
DO.add_pos!(dom1, pos)
@test DO.get_ncells(dom1) == 1

@static if get(ENV, "CI", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-10.0, 10.0))
    ax.set_ylim((-10.0, 10.0))
    Plot.domain!(ax, 1:2, dom1)
    Plot.trajectory!(ax, 1:2, sys, x, 50)
    Plot.cell_image!(ax, 1:2, dom1, sys)
    Plot.cell_approx!(ax, 1:2, dom1, sys, (10, 10))
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
