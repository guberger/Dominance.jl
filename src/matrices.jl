## Matrix set

struct MatrixSet{T,MT<:AbstractArray{T,2}}
    center::MT
    hull::Vector{MT}
    radius::T
end

function MatrixSet(Ac::AbstractArray{T,2}) where T
    return MatrixSet(Ac, [zero(Ac)], zero(T))
end

function MatrixSet(Ac::MT, A_hull::Vector{MT}) where {T,MT<:AbstractArray{T,2}}
    return MatrixSet(Ac, A_hull, zero(T))
end

function MatrixSet(Ac::MT, rad::T) where {T,MT<:AbstractArray{T,2}}
    return MatrixSet(Ac, [zero(Ac)], rad)
end
