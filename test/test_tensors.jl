include("../src/Dominance.jl")

module TestMain

using Test
using LinearAlgebra
using StaticArrays
using Main.Dominance
DO = Main.Dominance

sleep(0.1) # used for good printing
println("Started test")

@testset "tensors: norm" begin
θ = π/5.0
U = 2*SMatrix{2,2}(cos(θ), -sin(θ), sin(θ), cos(θ))
F_sys(x) = U*SVector(atan(x[1]), atan(x[2]))
DF_sys(x) = U*SMatrix{2,2}(1/(1 + x[1]^2), 0, 0, 1/(1 + x[2]^2))
T_tmp = Array{Float64}(undef, 2, 2, 2)
T_tmp[:, :, 1] = U*[1.0 0.0; 0.0 0.0]
T_tmp[:, :, 2] = U*[0.0 0.0; 0.0 -1.0]
T = SArray{Tuple{2,2,2}}(T_tmp...)
@test DO.tensor3d_normInf2matp(T, Inf) == opnorm(U, Inf)
@test DO.tensor3d_normInf2matp(T, 2) == opnorm(U, 2)
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
