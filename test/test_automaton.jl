include("../src/dominance.jl")

module TestMain

using Test
using Main.Dominance
DO = Main.Dominance

sleep(0.1) # used for good printing
println("Started test")

@testset "Symobic Model" begin
nstates = 10
nsymbols = 5
autom = DO.Automaton(nstates, nsymbols)

DO.add_transition!(autom, (5, 1, 9))
DO.add_transition!(autom, (5, 1, 8))
DO.add_transition!(autom, (5, 1, 3))
DO.add_transition!(autom, (8, 1, 3))
DO.add_transition!(autom, (5, 3, 5))
DO.add_transition!(autom, (8, 1, 3))
@test DO.get_ntransitions(autom) == 5
@test length(collect(DO.enum_transitions(autom))) == 5
translist = Tuple{Int,Int,Int}[]
DO.compute_post!(translist, autom, 5)
@test Set(translist) == Set([(5, 1, 9),(5, 1, 8),(5, 1, 3),(5, 3, 5)])
empty!(translist)
DO.compute_pre!(translist, autom, 3)
@test Set(translist) == Set([(5, 1, 3),(8, 1, 3),(8, 1, 3)])
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
