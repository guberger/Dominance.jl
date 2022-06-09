module ExampleMain

using LinearAlgebra
using StaticArrays
using PyPlot
using JuMP
using MosekTools
include("../src/Dominance.jl")
DO = Dominance
include("../src/plotting.jl")
include("./_contracting_.jl")

sleep(0.1) # used for good printing
println("Plot NonLin 1-dom contracting")

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

nsub = (3, 3)
DDF_norm_inf(x) = -DO.tensor3d_normInf2matp(DDContracting(x), Inf)
DDF_norm_2(x) = -DO.tensor3d_normInf2matp(DDContracting(x), 2)
f_opt, ~ = DO.minimize_over_domain(DDF_norm_inf, domain, nsub)
println(-f_opt)
bound_DDF_inf = -f_opt
# bound_DDF_inf = opnorm(U, Inf)*3*sqrt(3)/8

sys = DO.DiscSystem(Contracting, DContracting, bound_DDF_inf)
symmod = DO.symbolic_model_from_system(domain, sys, (2, 2))

statelist = 1:DO.get_ncells(domain)
viablelist = DO.viable_states(symmod.graph, statelist)
domain = DO.support_domain(symmod, viablelist)
symmod = DO.trim_symbolic_model(symmod, viablelist)
graph = symmod.graph

println("compute bound_DDF_2")
f_opt, ~ = DO.minimize_over_domain(DDF_norm_2, domain, nsub)
println(-f_opt)
bound_DDF_2 = -f_opt
# bound_DDF_2 = opnorm(U)*3*sqrt(3)/8
radius = bound_DDF_2*norm(h, Inf)/2
println(radius)
A_field = DO.sensitivity_matrices(domain, sys)
ASri_tmp = Any[]
for edge in DO.enum_edges(graph)
    source = edge.source
    target = edge.target
    pos = DO.get_pos_by_state(symmod, source)
    push!(ASri_tmp, edge => [(DO.MatrixSet(A_field[pos], radius), 1)])
    # push!(ASri_tmp, edge => [(A_field[source], 1)])
end
ASri_lab = Dict(ASri_tmp)
nr = 3
rate_tuple_iter = DO.hyper_range((0.5,), (0.7,), (nr,))
rate_tuple_iter = DO.hyper_range((0.6,), (0.6,), (1,))

println("start optim")
optim_solver = optimizer_with_attributes(Mosek.Optimizer)
P_opt, δ_opt, rates_opt = DO.cone_optim(graph, ASri_lab, rate_tuple_iter, optim_solver)

for i in eachindex(P_opt)
    ev = eigvals(P_opt[i])
    if ev[1] >= 0 || ev[2] <= 0
        print(eigvals(P_opt[i]), ", ")
    end
end
println(δ_opt)
println(rates_opt)
P_field = Dict([DO.get_pos_by_state(symmod, i) => P_opt[i] for i in eachindex(P_opt)])

pos = DO.get_pos_by_state(symmod, 1)
x = DO.get_coord_by_pos(symmod.grid, pos)

nsteps = 5
np = 50
rad = 0.5
fact = 1.1
fig = PyPlot.figure()
ax = fig.gca()
ax.set_xlim((-8.0, 8.0))
ax.set_ylim((-9.5, 9.5))
Plot.domain!(ax, 1:2, domain, ew = 0.1)
line = Plot.trajectory!(ax, 1:2, sys, x, nsteps, lc = "black")
Plot.add_arrow!(line[1])
Plot.cones!(ax, symmod.grid, sys, x, P_field, nsteps, rad, np)

end  # module ExampleMain
