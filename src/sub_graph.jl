# Essential graph

function sub_graph(graph, statelist)
    println("sub_graph started")
    idxn = Indexing(statelist)
    subgraph = Graph(length(statelist))
    nedges = 0
    for edge in graph.edges
        if edge.source ∈ statelist && edge.target ∈ statelist
            nedges += 1
            add_edge!(subgraph,
                get_index_by_elem(idxn, edge.source),
                get_index_by_elem(idxn, edge.target))
        end
    end
    println("sub_graph terminated: $(nedges) edges created")
    return subgraph, idxn
end
