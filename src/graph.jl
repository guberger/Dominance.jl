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
# Assumes not add twice same transition...
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
