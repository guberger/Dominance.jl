include("../src/Dominance.jl")

module ExampleMain

using LinearAlgebra
using StaticArrays
using PyPlot
using PyCall
using JuMP
using MosekTools
using Main.Dominance
DO = Main.Dominance

include("../src/plotting.jl")

sleep(0.1) # used for good printing
println("Plot SLS 1-dom simple")

matplotlib.rc("legend", fontsize = 15)
matplotlib.rc("axes", labelsize = 15)
matplotlib.rc("xtick", labelsize = 11)
matplotlib.rc("ytick", labelsize = 11)

## System

lb = SVector(-7.0, -7.0)
ub = SVector(7.0, 7.0)
x0 = SVector(0.0, 0.0)
h = SVector(1.0, 2.0)/20
grid = DO.Grid(x0, h)
domain = DO.Domain(grid)
DO.add_set!(domain, DO.HyperRectangle(lb, ub), DO.OUTER)
DO.remove_set!(domain, DO.HyperRectangle(lb/5, ub/5), DO.OUTER)

θ = π/5.0
U = 2*SMatrix{2,2}(cos(θ), -sin(θ), sin(θ), cos(θ))
F_sys(x) = U*SVector(atan(x[1]), atan(x[2]))
DF_sys(x) = U*SMatrix{2,2}(1/(1 + x[1]^2), 0, 0, 1/(1 + x[2]^2))
bound_DDF = opnorm(U, Inf)*3*sqrt(3)/8
bound_DDF_2 = opnorm(U)*3*sqrt(3)/8

sys = DO.DiscSystem(F_sys, DF_sys, bound_DDF)
graph, idxn1 = DO.symbolic_model(domain, sys)

viablelist = Int[]
for pos in DO.enum_pos(domain)
    push!(viablelist, DO.get_index_by_elem(idxn1, pos))
end

statelist = Int[]
DO.viable_states!(statelist, graph, viablelist)
subgraph, idxn2 = DO.sub_graph(graph, statelist)
nstates = DO.get_nstates(subgraph)
idxn = DO.compose(idxn2, idxn1)
A_field = DO.matrix_field(domain, sys, idxn, 1:nstates)

radius = bound_DDF_2*norm(h, Inf)/2
println(radius)
ASri_tmp = Any[]
for edge in DO.enum_edges(subgraph)
    source = edge.source
    target = edge.target
    push!(ASri_tmp, edge => [(DO.MatrixSet(A_field[source], radius), 1)])
    # push!(ASri_tmp, edge => [(A_field[source], 1)])
end
ASri_field = Dict(ASri_tmp)
nr = 3
rate_tuple_iter = DO.hyper_range((0.5,), (0.7,), nr)
rate_tuple_iter = DO.hyper_range((0.6,), (0.6,), 1)

println("start optim")
optim_solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
P_opt, δ_opt, rates_opt = DO.cone_optim(subgraph, ASri_field, rate_tuple_iter, optim_solver)

for i in eachindex(P_opt)
    ev = eigvals(P_opt[i])
    if ev[1] >= 0 || ev[2] <= 0
        print(eigvals(P_opt[i]), ", ")
    end
end
println(δ_opt)
println(rates_opt)

pos = DO.get_elem_by_index(idxn1, first(statelist))
x = DO.get_coord_by_pos(grid, pos)

dom1 = DO.Domain(grid)
for state = 1:DO.get_nstates(subgraph)
    DO.add_pos!(dom1, DO.get_elem_by_index(idxn, state))
end

nsteps = 5
np = 50
rad = 0.5
fact = 1.1
fig = PyPlot.figure()
ax = fig.gca()
ax.set_xlim((-8.0, 8.0))
ax.set_ylim((-9.5, 9.5))
Plot.domain!(ax, 1:2, dom1, ew = 0.1)
Plot.trajectory!(ax, 1:2, sys, x, nsteps)
Plot.cones!(ax, grid, sys, x, idxn, P_opt, nsteps, rad, np)

end  # module ExampleMain
