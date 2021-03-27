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

DO.add_transition!(graph, 5, 9)
DO.add_transition!(graph, 5, 8)
DO.add_transition!(graph, 5, 3)
DO.add_transition!(graph, 8, 3)
DO.add_transition!(graph, 5, 5)
DO.add_transition!(graph, 8, 3)
@test DO.get_ntransitions(graph) == 5
@test length(collect(DO.enum_transitions(graph))) == 5
translist = DO.Transition{Int}[]
DO.compute_post!(translist, graph, 5)
@test Set(translist) == Set(map(x -> DO.Transition(x...), [(5, 9), (5, 8), (5, 3), (5, 5)]))
empty!(translist)
DO.compute_pre!(translist, graph, 3)
@test Set(translist) == Set(map(x -> DO.Transition(x...), [(5, 3), (8, 3), (8, 3)]))
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
