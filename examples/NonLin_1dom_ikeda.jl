module ExampleMain

using LinearAlgebra
using StaticArrays
using PyPlot
using JuMP
using MosekTools
include("../src/Dominance.jl")
DO = Dominance
include("../src/plotting.jl")
include("./_ikeda_.jl")

sleep(0.1) # used for good printing
println("Plot NonLin 1-dom ikeda")

matplotlib.rc("legend", fontsize = 18)
matplotlib.rc("axes", labelsize = 16)
matplotlib.rc("xtick", labelsize = 15)
matplotlib.rc("ytick", labelsize = 15)

## System

lb = SVector(-1.1, -1.5)
ub = SVector(3.4, 1.8)
h = (ub - lb)./(20, 15)
x0 = lb + h/2
grid = DO.Grid(x0, h)
domain = DO.Domain(grid)
DO.add_set!(domain, DO.HyperRectangle(lb - h/4, ub + h/4), DO.INNER)
DDF_norm_inf(x) = -DO.tensor3d_normInf2matp(DDIkeda(x), Inf)
DDF_norm_2(x) = -DO.tensor3d_normInf2matp(DDIkeda(x), 2)
f_opt, ~ = DO.minimize_over_domain(DDF_norm_inf, domain, (10, 10))
println(-f_opt)
# bound_DDF_inf = -f_opt
bound_DDF_inf = 22.5

sys = DO.DiscSystem(Ikeda, DIkeda, bound_DDF_inf)
domain_list = DO.Domain[]
h_list = typeof(h)[]

nrounds = 6

for i = 1:nrounds
    global domain
    global grid
    symmod = DO.symbolic_model_from_system(domain, sys, (6, 6))
    global symmod
    statelist = 1:DO.get_ncells(domain)
    viablelist = DO.viable_states(symmod.graph, statelist)
    global viablelist
    domain = DO.support_domain(symmod, viablelist)
    push!(domain_list, domain)
    push!(h_list, domain.grid.h)
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
end
ASri_lab = Dict(ASri_tmp)
rate_tuple_iter = DO.hyper_range((1.0,), (1.0,), (1,))

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

fig = PyPlot.figure(figsize = (9.8, 7.0))
irounds = (2, 3, 4, 5)
ncols = ceil(Int, length(irounds)/2)
gs = matplotlib.gridspec.GridSpec(2, ncols, figure = fig, wspace = 0.2, hspace = 0.3)
AX_ = [fig.add_subplot(get(gs, i - 1)) for i in eachindex(irounds)]
extend = (-1, 1)

for (i, iround) in enumerate(irounds)
    ax = AX_[i]
    ax.set_xlim((lb[1], ub[1]) .+ 0.2 .*extend)
    ax.set_ylim((lb[2], 1.25) .+ 0.12 .*extend)
    if isodd(i)
        ax.set_ylabel(L"$x_2$")
    end
    if i >= ncols*2 - 1
        ax.set_xlabel(L"$x_1$")
    end
    ax.set_title("\$h=[$(h_list[iround][1]), $(h_list[iround][2])]^\\top\$", fontsize = 16)
    ew = norm(h_list[iround]) * 5.0
    display(ew)
    Plot.domain!(ax, 1:2, domain_list[iround], ew = ew)
end

fig.savefig("./figures/fig_NonLin_1dom_ikeda_model.png", transparent = false, dpi = 400,
    bbox_inches = matplotlib.transforms.Bbox(((0.58, 0.18), (8.85, 6.44))))

pos = DO.get_pos_by_state(symmod, 1)
x = DO.get_coord_by_pos(symmod.grid, pos)
domain1 = DO.support_domain(symmod, 1)

nsteps = 5
np = 50
rad = 0.4
fact = 1.1
fig = PyPlot.figure(figsize = (9.8, 7.0))
ax = fig.add_subplot()
ax.set_xlim((lb[1], ub[1]) .+ 0.3 .*extend)
ax.set_ylim((-1.6, 0.83) .+ 0.3 .*extend)
Plot.domain!(ax, 1:2, domain_list[4], ew = 0.3)
# Plot.cell_image!(ax, 1:2, domain1, sys)
# Plot.cell_approx!(ax, 1:2, domain1, sys)
line = Plot.trajectory!(ax, 1:2, sys, x, nsteps, lc = "black")
Plot.add_arrow!(line[1])
if !isdefined(ExampleMain, :P_field)
    P_field = Dict([DO.get_pos_by_state(symmod, i) => SMatrix{2,2}(1.0, 0.0, 0.0, -1.0)
        for i = 1:DO.get_nstates(symmod)])
end
Plot.cones!(ax, symmod.grid, sys, x, P_field, nsteps, rad, np)
ax.set_xlabel(L"$x_1$", fontsize = 20)
ax.set_ylabel(L"$x_2$", fontsize = 20)
for ax in fig.get_axes()
    ax.tick_params(axis = "both", which = "major", labelsize = 19)
end

fig.savefig("./figures/fig_NonLin_1dom_ikeda_cones.png", transparent = false, dpi = 400,
    bbox_inches = matplotlib.transforms.Bbox(((0.17, 0.06), (8.85, 6.18))))

end  # module ExampleMain
