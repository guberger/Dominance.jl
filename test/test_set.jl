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

@testset "Hyper range" begin
iter = DO.hyper_range((1,2,3), (1,3,6))
@test length(iter) == 8
@test (1,2,6) ∈ iter
@test (1,2.2,6) ∉ iter
iter = DO.hyper_range((1,2,3), (1,3,6), (5, 5, 2))
@test length(iter) == 50
@test (1.0,2,6) ∈ iter
@test (1,2,4) ∉ iter
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
