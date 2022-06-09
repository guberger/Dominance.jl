## Tensor norm

# max { ‖ ∑_i T[:, :, i]*h[i] ‖_p : ‖h‖_∞ = 1 }
function tensor3d_normInf2matp(T::SArray{Tuple{N1,N2,N3}}, p) where {N1,N2,N3}
    _ONE_ = ntuple(i -> (-1, 1), Val(N3))
    vert_iter = Iterators.product(_ONE_...)
    norm_vert(v) = opnorm(sum(T[:, :, i]*v[i] for i = 1:N3), p)
    return maximum(norm_vert, vert_iter)
end
