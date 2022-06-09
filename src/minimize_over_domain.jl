# Minimize a function over a domain

function minimize_over_domain(f, domain::Domain{N,T}, nsub) where {N,T}
    subpos_lims = ((nsub .- 1)./2)./nsub
    subpos_iter = hyper_range((-1).*subpos_lims, subpos_lims, nsub)
    x_opt = zero(SVector{N,T})
    f_opt = typemax(T)
    for pos in enum_pos(domain)
        for subpos in subpos_iter
            pos2 = pos .+ subpos
            x = get_coord_by_pos(domain.grid, pos2)
            fval = f(x)
            if fval < f_opt
                f_opt = fval
                x_opt = x
            end
        end
    end
    return f_opt, x_opt
end
