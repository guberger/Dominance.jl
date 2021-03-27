## Graph

struct Transition{S}
    source::S
    target::S
end

mutable struct Graph{S}
    states::Set{S}
    transitions::Set{Transition{S}}
end

function Graph(nstates::Int)
    return Graph(Set(1:nstates), Set{Transition{Int}}())
end

function get_nstates(graph)
    return length(graph.states)
end

function get_ntransitions(graph)
    return length(graph.transitions)
end

# In add_transition and add_transitions:
# Do not check that source, target are "inbounds"
# Assumes not add twice same transition...
function add_transition!(graph::Graph{S}, source::S, target::S) where S
    push!(graph.transitions, Transition(source, target))
end

Base.empty!(graph::Graph) = empty!(graph.transitions)

function enum_states(graph)
    return graph.states
end

function enum_transitions(graph)
    return graph.transitions
end

function compute_pre!(translist, graph, target)
    for trans in enum_transitions(graph)
        if trans.target == target
            push!(translist, trans)
        end
    end
end

function compute_post!(translist, graph, source)
    for trans in enum_transitions(graph)
        if trans.source == source
            push!(translist, trans)
        end
    end
end
