include("../src/Dominance.jl")

module TestMain

using Test
using LinearAlgebra
using StaticArrays
using Main.Dominance
DO = Main.Dominance

sleep(0.1) # used for good printing
println("Started test")

@testset "Macros: symbolic_model" begin
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
graph, idxn = DO.symbolic_model(domain, sys)
@test DO.get_nedges(graph) == 589

pos = (1, 2)
x = DO.get_coord_by_pos(grid, pos)
source = DO.get_index_by_pos(idxn, pos)

dom1 = DO.Domain(grid)
DO.add_pos!(dom1, pos)
dom2 = DO.Domain(grid)
edgelist = DO.Edge{Int}[]
DO.compute_post!(edgelist, graph, source)
for edge in edgelist
    DO.add_pos!(dom2, DO.get_pos_by_index(idxn, edge.target))
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
bound_DDF = norm(U, Inf)*3*sqrt(3)/8

sys = DO.DiscSystem(F_sys, DF_sys, bound_DDF)
graph, idxn = DO.symbolic_model(domain, sys)
@test DO.get_nedges(graph) == 717

pos = (1, 2)
x = DO.get_coord_by_pos(grid, pos)
source = DO.get_index_by_pos(idxn, pos)

dom1 = DO.Domain(grid)
DO.add_pos!(dom1, pos)
dom2 = DO.Domain(grid)
edgelist = DO.Edge{Int}[]
DO.compute_post!(edgelist, graph, source)
for edge in edgelist
    DO.add_pos!(dom2, DO.get_pos_by_index(idxn, edge.target))
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

@testset "Macros: viable_states!" begin
nstates = 10
graph = DO.Graph(1:nstates)

DO.add_edge!(graph, 5, 9)
DO.add_edge!(graph, 5, 8)
DO.add_edge!(graph, 5, 3)
DO.add_edge!(graph, 8, 3)
DO.add_edge!(graph, 5, 5)
DO.add_edge!(graph, 8, 5)
DO.add_edge!(graph, 1, 2)
DO.add_edge!(graph, 2, 4)
DO.add_edge!(graph, 4, 6)
DO.add_edge!(graph, 6, 7)
DO.add_edge!(graph, 7, 8)
DO.add_edge!(graph, 9, 10)
@test DO.get_nedges(graph) == 12

viablelist = 1:DO.get_nstates(graph)
statelist = Int[]
DO.viable_states!(statelist, graph, viablelist)
@test Set((statelist)) == Set([5, 8])
end

@testset "Macros: symbolic_model + viabel_statelistoller" begin
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
bound_DDF = norm(U, Inf)*3*sqrt(3)/8

sys = DO.DiscSystem(F_sys, DF_sys, bound_DDF)
graph, idxn = DO.symbolic_model(domain, sys)
DO.symbolic_model(domain, sys)
@test DO.get_nedges(graph) == 23468

viablelist = Int[]
for pos in DO.enum_pos(domain)
    push!(viablelist, DO.get_index_by_pos(idxn, pos))
end

statelist = Int[]
DO.viable_states!(statelist, graph, viablelist)
@test length(statelist) == 292

pos = (1, 2)
x = DO.get_coord_by_pos(grid, pos)

dom1 = DO.Domain(grid)
for state in statelist
    DO.add_pos!(dom1, DO.get_pos_by_index(idxn, state))
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
    Plot.trajectory!(ax, 1:2, sys, x, 50)
    Plot.cell_image!(ax, 1:2, dom1, sys)
    Plot.cell_approx!(ax, 1:2, dom1, sys)
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
