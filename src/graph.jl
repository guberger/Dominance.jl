## Graph

struct Edge{S}
    source::S
    target::S
end

mutable struct Graph{S,C}
    states::C
    edges::Set{Edge{S}}
end

function Graph(states::UnitRange{S}) where S
    return Graph(states, Set{Edge{S}}())
end

function get_nstates(graph)
    return length(graph.states)
end

function get_nedges(graph)
    return length(graph.edges)
end

# In add_edge:
# Do not check that source, target are "inbounds"
# Assumes not add twice same transition...
function add_edge!(graph::Graph{S}, source::S, target::S) where S
    push!(graph.edges, Edge(source, target))
end

Base.empty!(graph::Graph) = empty!(graph.edges)

function enum_states(graph)
    return graph.states
end

function enum_edges(graph)
    return graph.edges
end

function compute_pre!(translist, graph, target)
    for trans in enum_edges(graph)
        if trans.target == target
            push!(translist, trans)
        end
    end
end

function compute_post!(translist, graph, source)
    for trans in enum_edges(graph)
        if trans.source == source
            push!(translist, trans)
        end
    end
end
