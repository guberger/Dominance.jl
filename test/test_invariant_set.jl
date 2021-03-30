include("../src/Dominance.jl")

module TestMain

using Test
using LinearAlgebra
using StaticArrays
using Main.Dominance
DO = Main.Dominance

sleep(0.1) # used for good printing
println("Started test")

@testset "Invariant set: symbolic_model_from_system + viabel_states" begin
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
symmod = DO.symbolic_model_from_system(domain, sys)
@test DO.get_nedges(symmod.graph) == 25214

statelist = 1:DO.get_ncells(domain)
viablelist = DO.viable_states(symmod.graph, statelist)
@test length(viablelist) == 294

pos = (1, 2)
x = DO.get_coord_by_pos(grid, pos)

dom1 = DO.support_domain(symmod, viablelist)

@static if get(ENV, "CI", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-8.0, 8.0))
    ax.set_ylim((-9.5, 9.5))
    Plot.domain!(ax, 1:2, domain, fa = 0.1)
    Plot.domain!(ax, 1:2, dom1)
    Plot.trajectory!(ax, 1:2, sys, x, 50)
    Plot.cell_image!(ax, 1:2, dom1, sys)
    Plot.cell_approx!(ax, 1:2, dom1, sys)
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
