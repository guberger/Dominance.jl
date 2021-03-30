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
h = (ub - lb)./(20, 15)./10
grid = DO.Grid(x0, h)
domain = DO.Domain(grid)
DO.add_set!(domain, DO.HyperRectangle(lb, ub), DO.OUTER)
DDF_norm_inf(x) = -DO.tensor3d_normInf2matp(DDIkeda(x), Inf)
DDF_norm_2(x) = -DO.tensor3d_normInf2matp(DDIkeda(x), 2)
f_opt, ~ = DO.minimize_over_domain(DDF_norm_inf, domain, (10, 10))
println(-f_opt)
# bound_DDF_inf = -f_opt
bound_DDF_inf = 22.5

sys = DO.DiscSystem(Ikeda, DIkeda, bound_DDF_inf)

nrounds = 3

for i = 1:nrounds
    global domain
    global grid
    symmod = DO.symbolic_model_from_system(domain, sys)
    global symmod
    statelist = 1:DO.get_ncells(domain)
    viablelist = DO.viable_states(symmod.graph, statelist)
    global viablelist
    domain = DO.support_domain(symmod, viablelist)
    i == nrounds && break
    domain = DO.refine_domain(domain, (2, 2))
end

symmod = DO.trim_symbolic_model(symmod, viablelist)
graph = symmod.graph

println("compute bound_DDF_2")
f_opt, ~ = DO.minimize_over_domain(DDF_norm_2, domain, (3, 3))
println(-f_opt)
bound_DDF_2 = -f_opt
# bound_DDF_2 = opnorm(U)*3*sqrt(3)/8
radius = bound_DDF_2*norm(symmod.grid.h, Inf)/2
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
rate_tuple_iter = DO.hyper_range((1.0,), (1.0,), 1)

println("$(DO.get_nedges(graph)) edges")

# println("start optim")
# optim_solver = optimizer_with_attributes(Mosek.Optimizer)
# P_opt, δ_opt, rates_opt = DO.cone_optim(graph, ASri_lab, rate_tuple_iter, optim_solver)
#
# for i in eachindex(P_opt)
#     ev = eigvals(P_opt[i])
#     if ev[1] >= 0 || ev[2] <= 0
#         print(eigvals(P_opt[i]), ", ")
#     end
# end
# println(δ_opt)
# println(rates_opt)
# P_field = Dict([DO.get_pos_by_state(symmod, i) => P_opt[i] for i in eachindex(P_opt)])

pos = DO.get_pos_by_state(symmod, 1)
x = DO.get_coord_by_pos(symmod.grid, pos)
domain1 = DO.support_domain(symmod, 1)

nsteps = 5
np = 50
rad = 0.5
fact = 1.1
fig = PyPlot.figure()
ax = fig.gca()
ax.set_xlim((-1.1, 3.4).*1.1)
ax.set_ylim((-1.5, 1.8).*1.1)
Plot.domain!(ax, 1:2, domain, ew = 0.1)
Plot.cell_image!(ax, 1:2, domain1, sys)
Plot.cell_approx!(ax, 1:2, domain1, sys)
Plot.trajectory!(ax, 1:2, sys, x, nsteps)
Plot.cones!(ax, symmod.grid, sys, x, P_field, nsteps, rad, np)

end  # module ExampleMain
