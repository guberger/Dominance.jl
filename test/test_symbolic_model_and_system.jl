include("../src/Dominance.jl")

module TestMain

using Test
using LinearAlgebra
using StaticArrays
using Main.Dominance
DO = Main.Dominance

sleep(0.1) # used for good printing
println("Started test")

@testset "symbolic_model_from_system_from_system" begin
lb = SVector(0.0, 0.0)
ub = SVector(10.0, 11.0)
x0 = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)
grid = DO.Grid(x0, h)
domain = DO.Domain(grid)
DO.add_set!(domain, DO.HyperRectangle(lb, ub), DO.OUTER)

tstep = 0.5
nsys = 3
F_sys(x) = SVector(0.5, -cos(x[1]))
DF_sys(x) = SMatrix{2,2}(0.0, sin(x[1]), 0.0, 0.0)
bound_DF = 1.0
bound_DDF = 1.0

sys = DO.ContSystemRK4(tstep, F_sys, DF_sys, bound_DF, bound_DDF, nsys)
symmod = DO.symbolic_model_from_system(domain, sys)
@test DO.get_nedges(symmod.graph) == 589

pos = (1, 2)
x = DO.get_coord_by_pos(grid, pos)
source = DO.get_state_by_pos(symmod, pos)

dom1 = DO.Domain(grid)
DO.add_pos!(dom1, pos)
dom2 = DO.Domain(grid)
edgelist = DO.Edge[]
DO.compute_post!(edgelist, symmod.graph, source)
for edge in edgelist
    DO.add_pos!(dom2, DO.get_pos_by_state(symmod, edge.target))
end

@static if get(ENV, "CI", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-1.0, 11.0))
    ax.set_ylim((-2.0, 14.0))
    Plot.domain!(ax, 1:2, domain, fa = 0.1)
    Plot.domain!(ax, 1:2, dom1)
    Plot.domain!(ax, 1:2, dom2)
    Plot.trajectory!(ax, 1:2, sys, x, 50)
    Plot.cell_image!(ax, 1:2, dom1, sys)
    Plot.cell_approx!(ax, 1:2, dom1, sys)
end

lb = SVector(-7.0, -7.0)
ub = SVector(7.0, 7.0)
x0 = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)
grid = DO.Grid(x0, h)
domain = DO.Domain(grid)
DO.add_set!(domain, DO.HyperRectangle(lb, ub), DO.OUTER)

θ = π/5.0
U = 2*SMatrix{2,2}(cos(θ), -sin(θ), sin(θ), cos(θ))
F_sys(x) = U*SVector(atan(x[1]), atan(x[2]))
DF_sys(x) = U*SMatrix{2,2}(1/(1 + x[1]^2), 0, 0, 1/(1 + x[2]^2))
bound_DDF = opnorm(U, Inf)*3*sqrt(3)/8

sys = DO.DiscSystem(F_sys, DF_sys, bound_DDF)
symmod = DO.symbolic_model_from_system(domain, sys)
@test DO.get_nedges(symmod.graph) == 1065

pos = (1, 2)
x = DO.get_coord_by_pos(grid, pos)
source = DO.get_state_by_pos(symmod, pos)

dom1 = DO.Domain(grid)
DO.add_pos!(dom1, pos)
dom2 = DO.Domain(grid)
edgelist = DO.Edge[]
DO.compute_post!(edgelist, symmod.graph, source)
for edge in edgelist
    DO.add_pos!(dom2, DO.get_pos_by_state(symmod, edge.target))
end

@static if get(ENV, "CI", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-8.0, 8.0))
    ax.set_ylim((-9.5, 9.5))
    Plot.domain!(ax, 1:2, domain, fa = 0.1)
    Plot.domain!(ax, 1:2, dom1)
    Plot.domain!(ax, 1:2, dom2)
    Plot.trajectory!(ax, 1:2, sys, x, 50)
    Plot.cell_image!(ax, 1:2, dom1, sys)
    Plot.cell_approx!(ax, 1:2, dom1, sys)
end
end

@testset "Sensitivity matrices" begin
lb = SVector(-7.0, -7.0)
ub = SVector(7.0, 7.0)
x0 = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)/10
grid = DO.Grid(x0, h)
domain = DO.Domain(grid)
DO.add_set!(domain, DO.HyperRectangle(lb, ub), DO.OUTER)
DO.remove_set!(domain, DO.HyperRectangle(lb/5, ub/5), DO.OUTER)

θ = π/5.0
U = 2*SMatrix{2,2}(cos(θ), -sin(θ), sin(θ), cos(θ))
F_sys(x) = U*SVector(atan(x[1]), atan(x[2]))
DF_sys(x) = U*SMatrix{2,2}(1/(1 + x[1]^2), 0, 0, 1/(1 + x[2]^2))
bound_DDF = opnorm(U, Inf)*3*sqrt(3)/8

sys = DO.DiscSystem(F_sys, DF_sys, bound_DDF)
A_field = DO.sensitivity_matrices(domain, sys)
@test length(A_field) == DO.get_ncells(domain)
run = 0

for pos in DO.enum_pos(domain)
    run += 1
    run > 100 && break
    x = DO.get_coord_by_pos(domain.grid, pos)
    @test A_field[pos] == U*SMatrix{2,2}(1/(1 + x[1]^2), 0, 0, 1/(1 + x[2]^2))
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
