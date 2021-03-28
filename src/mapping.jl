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

# Fields

struct Field{S,T}
    elem2labels::Dict{S,Vector{T}}
end

function Field(::Type{S}, ::Type{T}) where {S,T}
    return Field(Dict{S,Vector{T}}())
end

function get_labels(field::Field{S,T}, elem) where {S,T}
    return get(field.elem2labels, elem, T[])
end

function add_label(field, elem, label)
    labels = get_labels(field, elem)
    push!(labels, label)
    field.elem2labels[elem] = labels
end
