## Hyper rectangle

struct HyperRectangle{VT}
    lb::VT
    ub::VT
end

function Base.in(x, rect::HyperRectangle)
    return all(rect.lb .<= x .<= rect.ub)
end
# all(x .<= y) is (surprisingly) faster than all(i -> x[i] <= y[i], eachindex(x))

function Base.isempty(rect::HyperRectangle)
    return any(rect.lb .> rect.ub)
end

function Base.intersect(a::HyperRectangle, b::HyperRectangle)
    return HyperRectangle(max.(a.lb, b.lb), min.(a.ub, b.ub))
end

function Base.issubset(a::HyperRectangle, b::HyperRectangle)
    return all(a.lb .>= b.lb) && all(a.ub .<= b.ub)
end

## Polyhedron

# Polyhedron defined by {x : |[Ax]_i| ≦ b_i ∀ i}
struct CenteredPolyhedron{MT,VT}
    A::MT
    b::VT
end

function Base.in(x, H::CenteredPolyhedron)
    return all(abs.(H.A*x) .<= H.b)
end
# all(x .<= y) is (surprisingly) faster than all(i -> x[i] <= y[i], eachindex(x))

## Hyper range

function hyper_range(lb::NTuple{N,T}, ub::NTuple{N,T}) where {N,T}
    ranges = ntuple(i -> UnitRange(lb[i], ub[i]), Val(N))
    return Iterators.product(ranges...)
end

function hyper_range(lb::NTuple{N,T}, ub::NTuple{N,T}, lengths::NTuple{N,Int}) where {N,T}
    ranges = ntuple(i -> range(lb[i], ub[i], length = lengths[i]), Val(N))
    return Iterators.product(ranges...)
end
