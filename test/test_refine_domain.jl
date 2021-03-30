include("../src/Dominance.jl")

module TestMain

using Test
using LinearAlgebra
using StaticArrays
using Main.Dominance
DO = Main.Dominance

sleep(0.1) # used for good printing
println("Started test")

@testset "refine_domain" begin
lb = SVector(0.0, 0.0)
ub = SVector(10.0, 11.0)
x0 = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)
grid = DO.Grid(x0, h)
domain = DO.Domain(grid)
DO.add_set!(domain, DO.HyperRectangle(lb, ub), DO.OUTER)
idxn = DO.Indexing(DO.enum_pos(domain))
statelist = 1:4:DO.get_ncells(domain)
domain2 = DO.Domain(grid)
for state in statelist
    pos = DO.get_elem_by_index(idxn, state)
    DO.add_pos!(domain2, pos)
end
@test DO.get_ncells(domain2) == ((DO.get_ncells(domain) - 1) ÷ 4) + 1

domain3 = DO.refine_domain(domain2, (2, 3))
@test DO.get_ncells(domain3) === DO.get_ncells(domain2)*6
for pos in DO.enum_pos(domain2)
    x = DO.get_coord_by_pos(domain2.grid, pos)
    pos3 = DO.get_pos_by_coord(domain3.grid, x)
    @test pos3 ∈ domain3
end
for pos in DO.enum_pos(domain3)
    x = DO.get_coord_by_pos(domain3.grid, pos)
    pos2 = DO.get_pos_by_coord(domain2.grid, x)
    @test pos2 ∈ domain2
end

@static if get(ENV, "CI", "false") == "false"
    include("../src/plotting.jl")
    using PyPlot
    fig = PyPlot.figure()
    ax = fig.gca()
    ax.set_xlim((-1.0, 11.0))
    ax.set_ylim((-2.0, 14.0))
    Plot.domain!(ax, 1:2, domain, fa = 0.1)
    Plot.domain!(ax, 1:2, domain2)
    Plot.domain!(ax, 1:2, domain3, fc = "none", fa = 0.5)
end
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
