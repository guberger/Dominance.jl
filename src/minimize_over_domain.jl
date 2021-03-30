# Minimize a function over a domain

function minimize_over_domain(f, domain::Domain{N,T}, nsub) where {N,T}
    hsub = domain.grid.h./nsub
    nsub_bis = (nsub .- 1)./2
    sub_iter = hyper_range((-1).*nsub_bis, nsub_bis)
    x_opt = zero(SVector{N,T})
    f_opt = typemax(T)
    for pos in enum_pos(domain)
        for sub in sub_iter
            pos2 = pos.*sub
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
