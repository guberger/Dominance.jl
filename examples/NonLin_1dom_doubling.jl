module ExampleMain

using LinearAlgebra
using StaticArrays
using PyPlot
using JuMP
using MosekTools
include("../src/Dominance.jl")
DO = Dominance
include("../src/plotting.jl")
include("./_doubling_.jl")

sleep(0.1) # used for good printing
println("Plot NonLin 1-dom doubling")

matplotlib.rc("legend", fontsize = 15)
matplotlib.rc("axes", labelsize = 15)
matplotlib.rc("xtick", labelsize = 11)
matplotlib.rc("ytick", labelsize = 11)

## System

lim1 = 3.0
lb = SVector(-lim1, -lim1)
ub = SVector(lim1, lim1)
x0 = SVector(0.0, 0.0)
h = SVector(1.0, 1.0)/40
grid = DO.Grid(x0, h)
domain = DO.Domain(grid)
domain1 = domain
DO.add_set!(domain, DO.HyperRectangle(lb, ub), DO.OUTER)
DO.remove_set!(domain, DO.HyperRectangle(lb/10, ub/10), DO.OUTER)

nsub = (3, 3)
DDF_norm_inf(x) = -DO.tensor3d_normInf2matp(DDDoubling(x), Inf)
DDF_norm_2(x) = -DO.tensor3d_normInf2matp(DDDoubling(x), 2)
f_opt, ~ = DO.minimize_over_domain(DDF_norm_inf, domain, nsub)
println(-f_opt)
bound_DDF_inf = -f_opt
# bound_DDF_inf = opnorm(U, Inf)*3*sqrt(3)/8

sys = DO.DiscSystem(Doubling, DDoubling, bound_DDF_inf)
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
rate_tuple_iter = DO.hyper_range((1.85,), (1.85,), (1,))
# rate_tuple_iter = DO.hyper_range((2.15,), (2.15,), (1,))

println("$(DO.get_nedges(graph)) edges")

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

pos = DO.get_pos_by_state(symmod, 10)
x = DO.get_coord_by_pos(symmod.grid, pos)

nsteps = 3
np = 50
rad = 0.5
fact = 1.1
fig = PyPlot.figure(figsize = (5.8, 5.8))
ax = fig.add_subplot(aspect = "equal")
lim2 = 1.5
ax.set_xlim((-lim2, lim2))
ax.set_ylim((-lim2, lim2))
Plot.domain!(ax, 1:2, domain1, ew = 0.1)
line = Plot.trajectory!(ax, 1:2, sys, x, nsteps, lc = "black")
Plot.add_arrow!(line[1])
Plot.cones!(ax, symmod.grid, sys, x, P_field, nsteps, rad, np)
ax.set_xlabel(L"$x_1$")
ax.set_ylabel(L"$x_2$")

fig.savefig("./figures/fig_NonLin_1dom_doubling_cones.png", transparent = false,
    bbox_inches = matplotlib.transforms.Bbox(((0.01, 0.13), (5.32, 5.18))))

end  # module ExampleMain
