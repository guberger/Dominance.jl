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

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
