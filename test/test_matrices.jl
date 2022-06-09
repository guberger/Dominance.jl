module TestMain

using Test
using StaticArrays
@static if isdefined(Main, :TestLocal)
    include("../src/Dominance.jl")
else
    using Dominance
end
DO = Dominance

sleep(0.1) # used for good printing
println("Started test")

@testset "MatrixSet" begin
matset = DO.MatrixSet([1 2; 4 5])
@test matset.center == [1 2; 4 5]
@test iszero(matset.hull) && eltype(matset.hull) == Matrix{Int}
@test iszero(matset.radius) && typeof(matset.radius) == Int

matset = DO.MatrixSet((@SMatrix [1 2; 4.0 5]), [@SMatrix [1 2.0; 4 5]])
@test matset.center == @SMatrix [1 2; 4 5]
@test eltype(matset.hull) <: SMatrix{2,2,Float64}
@test iszero(matset.radius) && typeof(matset.radius) == Float64

matset = DO.MatrixSet([1 2; 4 5], 1)
@test matset.center == [1 2; 4 5]
@test iszero(matset.hull) && eltype(matset.hull) == Matrix{Int}
@test typeof(matset.radius) == Int
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
