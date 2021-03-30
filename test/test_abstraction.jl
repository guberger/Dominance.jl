include("../src/Dominance.jl")

module TestMain

using Test
using StaticArrays
using Main.Dominance
DO = Main.Dominance

sleep(0.1) # used for good printing
println("Started test")

@testset "Grid + Domain" begin
orig = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)
grid = DO.Grid(orig, h)
domain1 = DO.Domain(grid)

DO.add_coord!(domain1, SVector(1.2, 3.5))
@test DO.get_ncells(domain1) == 1
DO.add_coord!(domain1, SVector(-0.5002, -1.0))
@test DO.get_ncells(domain1) == 2

DO.add_set!(domain1, DO.HyperRectangle(SVector(1.0, 0.0), SVector(11.0, 10.0)), DO.OUTER)
@test DO.get_ncells(domain1) == 67

DO.remove_coord!(domain1, SVector(2.0, 2.0))
@test DO.get_ncells(domain1) == 66
DO.remove_set!(domain1, DO.HyperRectangle(SVector(5.0, 5.0), SVector(10000.0, 10000.0)), DO.INNER)
@test DO.get_ncells(domain1) == 48

pos_iter = DO.enum_pos(domain1)
@test length(pos_iter) == 48

domain2 = DO.Domain(grid)
union!(domain2, domain1)
@test DO.get_ncells(domain2) == 48
DO.remove_set!(domain2, DO.HyperRectangle(SVector(1.0, 1.0), SVector(2.0, 2.0)), DO.OUTER)
@test DO.get_ncells(domain2) == 45
DO.add_subset!(domain2, domain1, DO.HyperRectangle(SVector(0.0, 0.0), SVector(5.0, 5.0)), DO.INNER)
@test DO.get_ncells(domain2) == 46

@static if get(ENV, "CI", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-2.0, 14.0))
    ax.set_ylim((-2.0, 14.0))
    Plot.domain!(ax, 1:2, domain2, fa = 0.3)
    Plot.set!(ax, 1:2, DO.HyperRectangle(SVector(1.0, 0.0), SVector(8.0, 10.0)))
    Plot.set!(ax, 1:2, DO.HyperRectangle(SVector(5.0, 5.0), SVector(10000.0, 10000.0)))
    Plot.set!(ax, 1:2, DO.HyperRectangle(SVector(1.0, 1.0), SVector(2.0, 2.0)))
    Plot.set!(ax, 1:2, DO.HyperRectangle(SVector(0.0, 0.0), SVector(5.0, 5.0)))
end
end

@testset "Symoblic model" begin
grid = DO.Grid(SVector(0, 0), SVector(0, 0))
graph = DO.Graph(1)
pos2ind, ind2pos = DO.indexing(Tuple{Int,Int}, [(1, 1), (2, 2)])
symmod = DO.SymbolicModel(grid, graph, pos2ind, ind2pos)

statelist = Int[]
push!(statelist, DO.get_state_by_pos(symmod, (1, 1)))
push!(statelist, DO.get_state_by_pos(symmod, (2, 2)))
sort!(statelist)
@test all(statelist .== [1, 2])

poslist = Tuple{Int, Int}[]
push!(poslist, DO.get_pos_by_state(symmod, 1))
push!(poslist, DO.get_pos_by_state(symmod, 2))
sort!(poslist)
@test all(poslist .== [(1, 1), (2, 2)])
end

@testset "Refine domain" begin
lb = SVector(0.0, 0.0)
ub = SVector(10.0, 11.0)
x0 = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)
grid = DO.Grid(x0, h)
domain = DO.Domain(grid)
DO.add_set!(domain, DO.HyperRectangle(lb, ub), DO.OUTER)
ind2pos = collect(DO.enum_pos(domain))
statelist = 1:4:DO.get_ncells(domain)
domain2 = DO.Domain(grid)
for state in statelist
    pos = ind2pos[state]
    DO.add_pos!(domain2, pos)
end
@test DO.get_ncells(domain2) == ((DO.get_ncells(domain) - 1) ÷ 4) + 1

domain3 = DO.refine_domain(domain2, (2, 3))
@test DO.get_ncells(domain3) === DO.get_ncells(domain2)*6
for pos in DO.enum_pos(domain2)
    x = DO.get_coord_by_pos(domain2.grid, pos)
    pos3 = DO.get_pos_by_coord(domain3.grid, x)
    @test pos3 ∈ domain3
end
for pos in DO.enum_pos(domain3)
    x = DO.get_coord_by_pos(domain3.grid, pos)
    pos2 = DO.get_pos_by_coord(domain2.grid, x)
    @test pos2 ∈ domain2
end

@static if get(ENV, "CI", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-1.0, 11.0))
    ax.set_ylim((-2.0, 14.0))
    Plot.domain!(ax, 1:2, domain, fa = 0.1)
    Plot.domain!(ax, 1:2, domain2)
    Plot.domain!(ax, 1:2, domain3, fc = "none", fa = 0.5)
end
end

@testset "Support domain" begin
grid = DO.Grid(SVector(0, 0), SVector(0, 0))
graph = DO.Graph(1)
pos2ind, ind2pos = DO.indexing(Tuple{Int,Int}, [(1, 1), (2, 2), (4, 6), (7, 8)])
symmod = DO.SymbolicModel(grid, graph, pos2ind, ind2pos)

domain = DO.support_domain(symmod, 1:3)
@test Set(DO.enum_pos(domain)) == Set(ind2pos[i] for i = 1:3)
domain = DO.support_domain(symmod, 4:4)
@test Set(DO.enum_pos(domain)) == Set((ind2pos[4],))
end

@testset "Trim symbolic model" begin
grid = DO.Grid(SVector(0), SVector(0))
nstates = 10
graph = DO.Graph(nstates)
pos2ind, ind2pos = DO.indexing(Tuple{Int}, (i,) for i = 1:nstates)
symmod = DO.SymbolicModel(grid, graph, pos2ind, ind2pos)

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

statelist = 1:2:DO.get_nstates(graph)
symmod2 = DO.trim_symbolic_model(symmod, statelist)
@test DO.get_nedges(symmod2.graph) == 3
for edge in DO.enum_edges(symmod2.graph)
    @test DO.Edge(
            DO.get_state_by_pos(symmod, DO.get_pos_by_state(symmod2, edge.source)),
            DO.get_state_by_pos(symmod, DO.get_pos_by_state(symmod2, edge.target))
        ) ∈ DO.enum_edges(graph)
end
for edge in DO.enum_edges(graph)
    @test iseven(edge.source) || iseven(edge.target) ||
        DO.Edge(
            DO.get_state_by_pos(symmod2, DO.get_pos_by_state(symmod, edge.source)),
            DO.get_state_by_pos(symmod2, DO.get_pos_by_state(symmod, edge.target))
        ) ∈ DO.enum_edges(symmod2.graph)
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
