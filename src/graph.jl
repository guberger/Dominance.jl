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

function get_nstates(graph::Graph)
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

## Macros

# Essential graph

function viable_states(graph, statelist)
    println("viable_states started")
    stateset = Set(statelist)
    npre = zeros(get_nstates(graph))
    npost = zeros(get_nstates(graph))
    for edge in enum_edges(graph)
        if edge.source ∈ stateset && edge.target ∈ stateset
            npre[edge.target] += 1
            npost[edge.source] += 1
        end
    end

    valid_edgeset = Set(enum_edges(graph))
    wrong_edgeset = Set{Edge}()
    while true
        empty!(wrong_edgeset)
        for edge in valid_edgeset
            if npre[edge.source] < 1 || npost[edge.target] < 1
                push!(wrong_edgeset, edge)
                npre[edge.target] -= 1
                npost[edge.source] -= 1
            end
        end
        if isempty(wrong_edgeset)
            break
        end
        setdiff!(valid_edgeset, wrong_edgeset)
    end

    nstates = 0
    statelist = Int[]
    for state in keys(npre)
        if npre[state] > 0 && npost[state] > 0
            push!(statelist, state)
            nstates += 1
        end
    end

    println("viable_states terminated: $(nstates) viable states")
    return statelist
end
