module TestMain

using Test
@static if isdefined(Main, :TestLocal)
    include("../src/Dominance.jl")
else
    using Dominance
end
DO = Dominance

sleep(0.1) # used for good printing
println("Started test")

@testset "Graph" begin
nstates = 10
graph = DO.Graph(nstates)

DO.add_edge!(graph, 5, 9)
DO.add_edge!(graph, 5, 8)
DO.add_edge!(graph, 5, 3)
DO.add_edge!(graph, 8, 3)
DO.add_edge!(graph, 5, 5)
DO.add_edge!(graph, 8, 3)
@test DO.get_nedges(graph) == 5
@test length(collect(DO.enum_edges(graph))) == 5
@test length(collect(DO.enum_states(graph))) == 10
edgelist = DO.Edge[]
DO.compute_post!(edgelist, graph, 5)
@test Set(edgelist) == Set(map(x -> DO.Edge(x...), [(5, 9), (5, 8), (5, 3), (5, 5)]))
empty!(edgelist)
DO.compute_pre!(edgelist, graph, 3)
@test Set(edgelist) == Set(map(x -> DO.Edge(x...), [(5, 3), (8, 3), (8, 3)]))
end

@testset "viable_states" begin
nstates = 10
graph = DO.Graph(nstates)

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

statelist = 1:DO.get_nstates(graph)
viablelist = DO.viable_states(graph, statelist)
@test Set(viablelist) == Set([5, 8])
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
