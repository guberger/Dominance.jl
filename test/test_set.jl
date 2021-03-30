include("../src/Dominance.jl")

module TestMain

using Test
using StaticArrays
using Main.Dominance
DO = Main.Dominance

sleep(0.1) # used for good printing
println("Started test")

@testset "Hyper range" begin
iter = DO.hyper_range((1,2,3), (1,3,6))
@test length(iter) == 8
@test (1,2,6) ∈ iter
@test (1,2.2,6) ∉ iter
iter = DO.hyper_range((1,2,3), (1,3,6), 5)
@test length(iter) == 125
@test (1.0,2,6) ∈ iter
@test (1,2,4) ∉ iter
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
