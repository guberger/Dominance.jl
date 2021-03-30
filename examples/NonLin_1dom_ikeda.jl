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
include("./_ikeda_.jl")

sleep(0.1) # used for good printing
println("Plot SLS 1-dom simple")

matplotlib.rc("legend", fontsize = 15)
matplotlib.rc("axes", labelsize = 15)
matplotlib.rc("xtick", labelsize = 11)
matplotlib.rc("ytick", labelsize = 11)

## System

lb = SVector(-1.1, -1.5)
ub = SVector(3.4, 1.8)
x0 = SVector(0.0, 0.0)
h = (ub - lb)./(20, 15)
grid = DO.Grid(x0, h)
domain = DO.Domain(grid)
DO.add_set!(domain, DO.HyperRectangle(lb, ub), DO.OUTER)
DDF_norm_inf(x) = -DO.tensor3d_normInf2matp(DDIkeda(x), Inf)
DDF_norm_2(x) = -DO.tensor3d_normInf2matp(DDIkeda(x), 2)
f_opt, ~ = DO.minimum_domain(DDF_norm_inf, domain, (10, 10))
println(-f_opt)
bound_DDF_inf = -f_opt
bound_DDF_inf = 22.5

sys = DO.DiscSystem(Ikeda, DIkeda, bound_DDF_inf)

nrounds = 4
statelist = Int[]

for i = 1:nrounds
    global domain
    global grid
    graph, idxn1 = DO.symbolic_model(domain, sys)
    global graph
    global idxn1
    viablelist = 1:DO.get_ncells(domain)
    empty!(statelist)
    DO.viable_states!(statelist, graph, viablelist)
    domain = DO.Domain(grid)
    for state in statelist
        DO.add_pos!(domain, DO.get_elem_by_index(idxn1, state))
    end
    i == nrounds && break
    domain = DO.refine_domain(domain, (2, 2))
    grid = domain.grid
end

subgraph, idxn2 = DO.sub_graph(graph, statelist)
nstates = DO.get_nstates(subgraph)
idxn = DO.compose(idxn2, idxn1)


#=
println("compute bound_DDF_2")
f_opt, ~ = DO.minimum_domain(DDF_norm_2, domain, nsub)
println(-f_opt)
bound_DDF_2 = -f_opt
# bound_DDF_2 = opnorm(U)*3*sqrt(3)/8
radius = bound_DDF_2*norm(h, Inf)/2
println(radius)
A_field = DO.matrix_field(domain, sys, idxn, 1:nstates)
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
=#

pos = DO.get_elem_by_index(idxn1, first(statelist))
x = DO.get_coord_by_pos(grid, pos)

nsteps = 5
np = 50
rad = 0.5
fact = 1.1
fig = PyPlot.figure()
ax = fig.gca()
ax.set_xlim((-1.1, 3.4).*1.1)
ax.set_ylim((-1.5, 1.8).*1.1)
Plot.domain!(ax, 1:2, domain, ew = 0.1)
Plot.trajectory!(ax, 1:2, sys, x, nsteps)
# Plot.cones!(ax, grid, sys, x, idxn, P_opt, nsteps, rad, np)

end  # module ExampleMain
