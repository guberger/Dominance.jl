## Matrix set

struct MatrixSet{T,MT<:AbstractMatrix{T}}
    center::MT
    hull::Vector{MT}
    radius::T
end

function MatrixSet(Ac::AbstractMatrix{T}) where T
    return MatrixSet(Ac, [zero(Ac)], zero(T))
end

function MatrixSet(Ac::MT, A_hull::Vector{MT}) where {T,MT<:AbstractMatrix{T}}
    return MatrixSet(Ac, A_hull, zero(T))
end

function MatrixSet(Ac::MT, rad::T) where {T,MT<:AbstractMatrix{T}}
    return MatrixSet(Ac, [zero(Ac)], rad)
end
