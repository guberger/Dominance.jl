include("../src/Dominance.jl")

module TestMain

using Test
using LinearAlgebra
using StaticArrays
using Main.Dominance
DO = Main.Dominance

sleep(0.1) # used for good printing
println("Started test")

@testset "viable_states!" begin
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

statelist = 1:2:DO.get_nstates(graph)
subgraph, idxn = DO.sub_graph(graph, statelist)
@test DO.get_nedges(subgraph) == 3
for edge in DO.enum_edges(subgraph)
    @test DO.Edge(DO.get_elem_by_index(idxn, edge.source),
        DO.get_elem_by_index(idxn, edge.target)) ∈ DO.enum_edges(graph)
end
for edge in DO.enum_edges(graph)
    @test iseven(edge.source) || iseven(edge.target) ||
        DO.Edge(DO.get_index_by_elem(idxn, edge.source),
            DO.get_index_by_elem(idxn, edge.target)) ∈ DO.enum_edges(subgraph)

end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
