# Analysis of Lipschitz constant, etc.

function tensor3d_normInf2matp(Tensor3d, p, vert_iter, dim)
    # max { || D^2f[h] \in R^{2x2} ||_{p, p} : ||h||_inf = 1 }.
    @inline f(v) = opnorm(sum([Tensor3d[:, :, i]*v[i] for i = 1:dim]), p)
    return maximum(f, vert_iter)
end

function rectangle_penalty(x, x0, h)
    # Normalized inf-distance from rectangle.
    delta = max.(0.0, x0 - x, x - x0 - h)
    return sum(delta./h)
end

function minimum_grid(f, x0, h, nsub, dim, do_print)
    # Min of f using NelderMead starting from prod(nsub) initial points in rectangle
    # defined by x0 and h.
    ht = h./nsub
    nnsub = prod(nsub)
    subs = [1:nsub[i] for i = 1:dim]
    sub_iter = Iterators.product(subs...)
    fmin = Inf
    xopt = fill(NaN, dim)

    @inline F0(x) = f(x) + 1e5 * rectangle_penalty(x, x0, h)
    Optimizer = dim > 1 ? NelderMead() : BFGS()

    for (i, it) in enumerate(sub_iter)
        # Center of current subrectangle (Tuple .* Array is well-defined)
        xt = x0 + it.*ht - 0.5*ht
        res = optimize(F0, xt, Optimizer)
        if do_print
            @printf("iter %d (/%d) (sub: %s) -> fmin: %f, xopt: %s\n",
                i, nnsub, it, res.minimum, res.minimizer)
        end
        if res.minimum < fmin
            fmin = res.minimum
            xopt = res.minimizer
        end
    end

    return (fmin, xopt)
end

function DDsNorm_function(theAbs, p)
    dim = theAbs.dim
    F! = theAbs.F!
    DF! = theAbs.DF!
    DDF! = theAbs.DDF!
    vert_iter = Iterators.product(fill(-1:2:1, dim)...)
    @inline secordsys_vector!(du, u, p_dummy, t) =
        secondordersystem_vector!(F!, DF!, DDF!, u, du, dim)
    @inline secordsys_map(x) = secondordersystem_map(secordsys_vector!, x,
            theAbs.tstep, theAbs.ode_reltol, theAbs.ode_abstol, dim)
    @inline F0(x) = -tensor3d_normInf2matp(secordsys_map(x)[3], p, vert_iter, dim)
    return F0
end

function minimize_function_cell(theAbs, idx, F0)
    # Minimize the function F0 over the cell with index idx, and with nsub given
    # by theAbs.nsubopt_list[idx]
    ht = theAbs.cell_list[idx].Lev.*theAbs.h
    xt = theAbs.x_or + theAbs.cell_list[idx].coords.*theAbs.h
    nsub = theAbs.nsubopt_list[idx]
    fmin = -minimum_grid(F0, xt, ht, nsub, theAbs.dim, false)[1]
    return fmin
end
