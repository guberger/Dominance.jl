include("../src/Dominance.jl")

module TestMain

using Test
using Main.Dominance
DO = Main.Dominance

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

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
