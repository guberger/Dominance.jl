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

@testset "cone_optim single" begin
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

ASri_field = Dict([DO.Edge(1, 1) => [(DO.MatrixSet(A1), 1)],
    DO.Edge(1, 2) => [(DO.MatrixSet(A2), 2)],
    DO.Edge(2, 1) => [(DO.MatrixSet(A1), 2)],
    DO.Edge(2, 2) => [(DO.MatrixSet(A2), 1)]])
rate_tuple_iter = Iterators.product((γ1,), (sqrt(γ2),))

optim_solver = optimizer_with_attributes(SDPA.Optimizer)
~, δ_opt, rates_opt = DO.cone_optim(graph, ASri_field, rate_tuple_iter, optim_solver)
@test δ_opt > 0.0225941*0.999
@test rates_opt == (γ1, sqrt(γ2))

@static if get(ENV, "CI", "false") == "false"
    using MosekTools
    optim_solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    ~, δ_opt, rates_opt = DO.cone_optim(graph, ASri_field, rate_tuple_iter, optim_solver)
    @test δ_opt > 0.0225946*0.999
    @test rates_opt == (γ1, sqrt(γ2))
end
end

@testset "cone_optim set simple" begin
α = 0.15 # limit: (6.0 - sqrt(32.0))/2.0 = 0.1715728752538097
A1 = @SMatrix [1.0 0.0; 1.0-α α]
A2 = @SMatrix [α α-1.0; 0.0 1.0]
Ad = (@SMatrix [1.0 1.0; 1.0 1.0])*1e-6
radius = 1e-6

ev1 = eigvals(Matrix(A1))
γ1 = sqrt(prod(norm.(ev1)))
ev2 = eigvals(Matrix(A1*A2))
γ2 = sqrt(prod(norm.(ev2)))

graph = DO.Graph(2)
DO.add_edge!(graph, 1, 1)
DO.add_edge!(graph, 1, 2)
DO.add_edge!(graph, 2, 1)
DO.add_edge!(graph, 2, 2)

ASri_field_1 = Dict([DO.Edge(1, 1) => [(DO.MatrixSet(A1), 1)],
    DO.Edge(1, 2) => [(DO.MatrixSet(A2), 2)],
    DO.Edge(2, 1) => [(DO.MatrixSet(A1), 2)],
    DO.Edge(2, 2) => [(DO.MatrixSet(A2), 1)]])
ASri_field_2 = Dict([DO.Edge(1, 1) => [(DO.MatrixSet(A1), 1)],
    DO.Edge(1, 2) => [(DO.MatrixSet(A2, [-Ad, Ad], radius), 2)],
    DO.Edge(2, 1) => [(DO.MatrixSet(A1, [-Ad, Ad], radius), 2)],
    DO.Edge(2, 2) => [(DO.MatrixSet(A2, [-Ad, Ad], radius), 1)]])
rate_tuple_iter = Iterators.product((γ1,), (sqrt(γ2),))

optim_solver = optimizer_with_attributes(SDPA.Optimizer)
~, δ_opt_1, rates_opt_1 = DO.cone_optim(graph, ASri_field_1, rate_tuple_iter, optim_solver)
~, δ_opt_2, rates_opt_2 = DO.cone_optim(graph, ASri_field_2, rate_tuple_iter, optim_solver)
@test δ_opt_1*0.999 < δ_opt_2 < δ_opt_1*1.001
@test all(rates_opt_1 .≈ rates_opt_2)
@test δ_opt_2 > 0.0117813*0.999

@static if get(ENV, "CI", "false") == "false"
    using MosekTools
    optim_solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    ~, δ_opt_1, rates_opt_1 = DO.cone_optim(graph, ASri_field_1, rate_tuple_iter, optim_solver)
    ~, δ_opt_2, rates_opt_2 = DO.cone_optim(graph, ASri_field_2, rate_tuple_iter, optim_solver)
    @test δ_opt_1*0.999 < δ_opt_2 < δ_opt_1*1.001
    @test all(rates_opt_1 .≈ rates_opt_2)
    @test δ_opt_2 > 0.0117819*0.999
end
end

@testset "cone_optim set duffing" begin
Ac_list = ((@SMatrix [1.0 0.3 0.0; -0.375 0.7 0.03; -1.5 0.0 0.925]),
(@SMatrix [1.0 0.3 0.0; -0.225 0.7 0.03; -1.5 0.0 0.925]),
(@SMatrix [1.0 0.3 0.0; -0.075 0.7 0.03; -1.5 0.0 0.925]),
(@SMatrix [1.0 0.3 0.0; 0.075 0.7 0.03; -1.5 0.0 0.925]))
Ad = @SMatrix [0.0 0.0 0.0; 0.075 0.0 0.0; 0.0 0.0 0.0]

graph = DO.Graph(4)
ASri_tmp = Any[]
for i = 1:4, j = 1:4
    DO.add_edge!(graph, i, j)
    push!(ASri_tmp, DO.Edge(i, j) => [(DO.MatrixSet(Ac_list[i], [Ad, -Ad]), 1)])
end
ASri_field = Dict(ASri_tmp)
γ_min = 0.0
γ_max = Inf
for Ac in Ac_list
    for Avert in (-Ad, 0.0, Ad)
        EG = DO.pth_eigval(Matrix(Ac .+ Avert), 2, 1e-9)
        γ_max = min(γ_max, EG[1])
        γ_min = max(γ_min, EG[2])
    end
end
nr = 9
rate_tuple_iter = DO.hyper_range((γ_min,), (γ_max,), nr)

optim_solver = optimizer_with_attributes(SDPA.Optimizer)
~, δ_opt, rates_opt = DO.cone_optim(graph, ASri_field, rate_tuple_iter, optim_solver)
@test δ_opt > 0.000198015*0.999
@test all(rates_opt .≈ (0.8272967578345953,))

@static if get(ENV, "CI", "false") == "false"
    using MosekTools
    optim_solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    ~, δ_opt, rates_opt = DO.cone_optim(graph, ASri_field, rate_tuple_iter, optim_solver)
    @test δ_opt > 0.000198138*0.999
    @test all(rates_opt .≈ (0.8272967578345953,))
end
end

@testset "cone_optim set radius" begin
Ac_list = ((@SMatrix [0.0 2.0; 1.0 0.0]),
    (@SMatrix [0.0 1.0; 2.0 0.0]))
Ad = @SMatrix [0.1 0.1; 0.1 0.1]
radius = 0.1

graph = DO.Graph(2)
DO.add_edge!(graph, 1, 2)
DO.add_edge!(graph, 2, 1)
ASri_field = Dict([DO.Edge(1, 2) => [(DO.MatrixSet(Ac_list[1], [Ad, -Ad], radius), 1)],
    DO.Edge(2, 1) => [(DO.MatrixSet(Ac_list[2], [Ad, -Ad], radius), 1)]])
γ_min = 0.0
γ_max = Inf
for Ac in Ac_list
    for Avert in (-Ad, 0.0, Ad)
        EG = DO.pth_eigval(Matrix(Ac .+ Avert), 1, 1e-9)
        γ_max = min(γ_max, EG[1])
        γ_min = max(γ_min, EG[2])
    end
end
nr = 9
rate_tuple_iter = DO.hyper_range((γ_min,), (γ_max,), nr)

optim_solver = optimizer_with_attributes(SDPA.Optimizer)
~, δ_opt, rates_opt = DO.cone_optim(graph, ASri_field, rate_tuple_iter, optim_solver)
@test δ_opt > 0.579326*0.999
@test all(rates_opt .≈ (1.4198684167769349,))

@static if get(ENV, "CI", "false") == "false"
    using MosekTools
    optim_solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    ~, δ_opt, rates_opt = DO.cone_optim(graph, ASri_field, rate_tuple_iter, optim_solver)
    @test δ_opt > 0.579327*0.999
    @test all(rates_opt .≈ (1.4198684167769349,))
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
