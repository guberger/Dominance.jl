include("../src/Dominance.jl")

module TestMain

using Test
using LinearAlgebra
using StaticArrays
using JuMP
using SDPA
using Main.Dominance
DO = Main.Dominance

sleep(0.1) # used for good printing
println("Started test")

@testset "cone_optim_single" begin
α = 0.1 # limit: (6.0 - sqrt(32.0))/2.0 = 0.1715728752538097
A1 = @SMatrix [1.0 0.0; 1.0-α α]
A2 = @SMatrix [α α-1.0; 0.0 1.0]

ev1 = eigvals(Matrix(A1))
γ1 = sqrt(prod(norm.(ev1)))
ev2 = eigvals(Matrix(A1*A2))
γ2 = sqrt(prod(norm.(ev2)))

graph = DO.Graph(2)
DO.add_edge!(graph, 1, 1)
DO.add_edge!(graph, 1, 2)
DO.add_edge!(graph, 2, 1)
DO.add_edge!(graph, 2, 2)

Ari_field = Dict([DO.Edge(1, 1) => [(A1, 1)],
    DO.Edge(1, 2) => [(A2, 2)],
    DO.Edge(2, 1) => [(A1, 2)],
    DO.Edge(2, 2) => [(A2, 1)]])
rate_tuple_iter = Iterators.product((γ1,), (sqrt(γ2),))

optim_solver = optimizer_with_attributes(SDPA.Optimizer)
~, δ_opt, rates_opt = DO.cone_optim_single(graph, Ari_field, rate_tuple_iter, optim_solver)
@test δ_opt > 0.0225941*0.999
@test rates_opt == (γ1, sqrt(γ2))

@static if get(ENV, "CI", "false") == "false"
    using MosekTools
    optim_solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    ~, δ_opt, rates_opt = DO.cone_optim_single(graph, Ari_field, rate_tuple_iter, optim_solver)
    @test δ_opt > 0.0225946*0.999
    @test rates_opt == (γ1, sqrt(γ2))
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
