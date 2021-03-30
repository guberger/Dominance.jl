include("../src/Dominance.jl")

module TestMain

using Test
using LinearAlgebra
using StaticArrays
using Main.Dominance
DO = Main.Dominance

sleep(0.1) # used for good printing
println("Started test")

@testset "symbolic_model + viabel_states" begin
lb = SVector(-7.0, -7.0)
ub = SVector(7.0, 7.0)
x0 = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)/10
grid = DO.Grid(x0, h)
domain = DO.Domain(grid)
DO.add_set!(domain, DO.HyperRectangle(lb, ub), DO.OUTER)
DO.remove_set!(domain, DO.HyperRectangle(lb/5, ub/5), DO.OUTER)

θ = π/5.0
U = 2*SMatrix{2,2}(cos(θ), -sin(θ), sin(θ), cos(θ))
F_sys(x) = U*SVector(atan(x[1]), atan(x[2]))
DF_sys(x) = U*SMatrix{2,2}(1/(1 + x[1]^2), 0, 0, 1/(1 + x[2]^2))
bound_DDF = opnorm(U, Inf)*3*sqrt(3)/8

sys = DO.DiscSystem(F_sys, DF_sys, bound_DDF)
graph, idxn = DO.symbolic_model(domain, sys)
@test DO.get_nedges(graph) == 25214

statelist = collect(1:2:DO.get_nstates(graph))
A_field = DO.matrix_field(domain, sys, idxn, statelist)
run = 0

for pos in DO.enum_pos(domain)
    run += 1
    run > 100 && break
    state = DO.get_index_by_elem(idxn, pos)
    @test iseven(state) ⊻ haskey(A_field, state)
    iseven(state) && continue
    x = DO.get_coord_by_pos(domain.grid, pos)
    @test A_field[state] == U*SMatrix{2,2}(1/(1 + x[1]^2), 0, 0, 1/(1 + x[2]^2))
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
