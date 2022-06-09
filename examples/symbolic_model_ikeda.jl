module ExampleMain

using Test
using LinearAlgebra
using StaticArrays
using PyPlot
using LightGraphs
using GraphPlot
using Compose
using Colors
using Cairo
using Fontconfig
include("../src/Dominance.jl")
DO = Dominance
include("../src/plotting.jl")
include("./_ikeda_.jl")

sleep(0.1) # used for good printing
println("Plot symbolic model ikeda")

matplotlib.rc("legend", fontsize = 15)
matplotlib.rc("axes", labelsize = 18)
matplotlib.rc("xtick", labelsize = 14)
matplotlib.rc("ytick", labelsize = 14)

## System

lb = SVector(-1.1, -1.5)
ub = SVector(3.4, 1.8)
h = (ub - lb)/4
x0 = lb + h/2
grid = DO.Grid(x0, h)
domain = DO.Domain(grid)
DO.add_set!(domain, DO.HyperRectangle(lb - h/4, ub + h/4), DO.INNER)
bound_DDF_inf = 22.5
# bound_DDF_inf = 0.0

sys = DO.DiscSystem(Ikeda, DIkeda, bound_DDF_inf)

x = SVector(2.0, 1.5)
pos = DO.get_pos_by_coord(grid, x)
domain1 = DO.Domain(grid)
DO.add_pos!(domain1, pos)

symmod = DO.symbolic_model_from_system(domain, sys, (10, 10))
graph = symmod.graph
domain2 = DO.Domain(grid)
source = DO.get_state_by_pos(symmod, pos)
edgelist = DO.Edge[]
DO.compute_post!(edgelist, symmod.graph, source)
for edge in edgelist
    DO.add_pos!(domain2, DO.get_pos_by_state(symmod, edge.target))
end

nsteps = 5
np = 50
rad = 0.5
fact = 1.1
fig = PyPlot.figure(figsize = (7.8, 5.8))
ax = fig.add_subplot(aspect = "equal")
extend = (-1, 1)
ax.set_xlim((lb[1], ub[1]) .+ 0.2 .*extend)
ax.set_ylim((lb[2], ub[2]) .+ 0.2 .*extend)
Plot.domain!(ax, 1:2, domain, ew = 1.5, fc = "none")
Plot.domain!(ax, 1:2, domain1, ew = 0.1)
Plot.domain!(ax, 1:2, domain2, ew = 0.1, fc = "blue", fa = 0.4)
Plot.cell_image!(ax, 1:2, domain1, sys, fa = 0.7)
for (i, pos) in enumerate(DO.enum_pos(domain))
    x = DO.get_coord_by_pos(grid, pos)
    ax.text(x..., i, fontsize = "18")
end
# Plot.cell_approx!(ax, 1:2, domain1, sys, fa = 0.3)
ax.set_xlabel(L"$x_1$")
ax.set_ylabel(L"$x_2$")

fig.savefig("./figures/fig_symbolic_model_image.png", transparent = false,
    bbox_inches = matplotlib.transforms.Bbox(((0.18, 0.07), (7.0, 5.18))))

nstates = DO.get_nstates(graph)
g = SimpleDiGraph(nstates)
for edge in DO.enum_edges(graph)
    add_edge!(g, edge.source, edge.target)
end
node_colors = [RGBA(0.0,0.0,0.0,0.0) for i = 1:nv(g)]
for edge in edgelist
    node_colors[edge.target] = RGBA(0.0,0.0,1.0,0.5)
    node_colors[edge.source] = RGBA(1.0,0.0,0.0,0.5)
end
edge_colors = [colorant"gray70" for i = 1:ne(g)]
edge_lws = [1.0 for i = 1:ne(g)]
for (i, edge) in enumerate(edges(g))
    if edge.src == source
        edge_colors[i] = colorant"salmon"
        edge_lws[i] = 2.0
    end
end
Compose.draw(PNG("./figures/fig_symbolic_model_graph.png", 16cm, 16cm, dpi = 400),
    gplot(g, layout = spring_layout, nodelabel = 1:nv(g),
        nodefillc = node_colors, nodestrokec = "black", nodestrokelw = 0.01,
        NODELABELSIZE = 10, NODESIZE = 0.1,
        edgestrokec = edge_colors, edgelinewidth = edge_lws))
end  # module ExampleMain
