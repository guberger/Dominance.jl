## Graph

struct Edge
    source::Int
    target::Int
end

struct Graph
    nstates::Int
    edges::Set{Edge}
end

function Graph(nstates)
    return Graph(nstates, Set{Edge}())
end

function get_nstates(graph)
    return graph.nstates
end

function get_nedges(graph)
    return length(graph.edges)
end

# In add_edge:
# Do not check that source, target are "inbounds"
# Assumes not add twice same edge...
function add_edge!(graph, source, target)
    push!(graph.edges, Edge(source, target))
end

Base.empty!(graph::Graph) = empty!(graph.edges)

function enum_states(graph)
    return 1:graph.nstates
end

function enum_edges(graph)
    return graph.edges
end

function compute_pre!(edgelist, graph, target)
    for edge in enum_edges(graph)
        if edge.target == target
            push!(edgelist, edge)
        end
    end
end

function compute_post!(edgelist, graph, source)
    for edge in enum_edges(graph)
        if edge.source == source
            push!(edgelist, edge)
        end
    end
end
