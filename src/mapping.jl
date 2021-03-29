## Indexing

struct Indexing{T}
    elem2ind::Dict{T,Int}
    ind2elem::Vector{T}
end

function Indexing(elemlist)
    ind2elem = Vector{eltype(elemlist)}(undef, length(elemlist))
    elem2ind = Dict{eltype(elemlist),Int}()
    sizehint!(elem2ind, length(elemlist))
    for (i, elem) in enumerate(elemlist)
        ind2elem[i] = elem
        elem2ind[elem] = i
    end
    return Indexing(elem2ind, ind2elem)
end

function get_elem_by_index(idxn, index)
    return idxn.ind2elem[index]
end

function get_index_by_elem(idxn, elem)
    return idxn.elem2ind[elem]
end

function compose(idxn_first, idxn_second)
    return Indexing(get_elem_by_index(idxn_second,
            get_elem_by_index(idxn_first, index)
        ) for index in 1:length(idxn_first.ind2elem))
end
