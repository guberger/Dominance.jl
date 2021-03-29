include("../src/Dominance.jl")

module TestMain

using Test
using StaticArrays
using Main.Dominance
DO = Main.Dominance

sleep(0.1) # used for good printing
println("Started test")

@testset "Mapping: Indexing" begin
elemlist = Set([(1, 1), (2, 2)])
idxn = DO.Indexing(elemlist)

indexlist = Int[]
push!(indexlist, DO.get_index_by_elem(idxn, (1, 1)))
push!(indexlist, DO.get_index_by_elem(idxn, (2, 2)))
sort!(indexlist)
@test all(indexlist .== [1, 2])

elemlist = Tuple{Int, Int}[]
push!(elemlist, DO.get_elem_by_index(idxn, 1))
push!(elemlist, DO.get_elem_by_index(idxn, 2))
sort!(elemlist)
@test all(elemlist .== [(1, 1), (2, 2)])
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
