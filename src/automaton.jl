## Automaton

# transitions: (source, symbol, target)
mutable struct Automaton
    nstates::Int
    nsymbols::Int
    transitions::Set{Tuple{Int,Int,Int}}
end

function Automaton(nstates, nsymbols)
    return Automaton(nstates, nsymbols, Set{Tuple{Int,Int,Int}}())
end

function get_ntransitions(autom)
    return length(autom.transitions)
end

# In add_transition and add_transitions:
# Do not check that source, target are "inbounds"
# Assumes not add twice same transition...
function add_transition!(autom, trans)
    push!(autom.transitions, trans)
end

Base.empty!(autom::Automaton) = empty!(autom.transitions)

function enum_transitions(autom)
    return autom.transitions
end

function compute_pre!(translist, autom, target)
    for trans in enum_transitions(autom)
        if trans[3] == target
            push!(translist, trans)
        end
    end
end

function compute_post!(translist, autom, source)
    for trans in enum_transitions(autom)
        if trans[1] == source
            push!(translist, trans)
        end
    end
end

## Controller

mutable struct Controller
    states::Set{Int}
end

Controller() = Controller(Set{Int}())

function add_state!(contr, state)
    push!(contr.states, state)
end

Base.empty!(contr::Controller) = empty!(contr.states)

function enum_states(contr)
    return contr.states
end
