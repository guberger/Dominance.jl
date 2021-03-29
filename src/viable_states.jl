# Essential graph

function viable_states!(statelist, graph, viablelist)
    println("viable_states! started")
    viableset = Set(viablelist)
    npre = zeros(get_nstates(graph))
    npost = zeros(get_nstates(graph))
    for edge in enum_edges(graph)
        if edge.source ∈ viableset && edge.target ∈ viableset
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
    for state in keys(npre)
        if npre[state] > 0 && npost[state] > 0
            push!(statelist, state)
            nstates += 1
        end
    end

    println("viable_states! terminated with success: $(nstates) viable states added")
end
