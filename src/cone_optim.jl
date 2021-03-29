## Solve dominance LMIs

function pth_eigval(A, p, shrink)
    e = sort(broadcast(abs, eigvals(A)), rev = true)
    return (e[p]*(1.0-shrink), e[p+1]/(1.0-shrink))
end

# Ari_field : dict (edge -> [(matrix, rate index)])
function cone_optim_single(graph, Ari_field, rate_tuple_iter, optim_solver)
    A_type = fieldtype(eltype(valtype(Ari_field)), 1)
    A_size = size(A_type)
    scalar_type = eltype(A_type)
    rate_tuple_type = eltype(rate_tuple_iter)
    ee_opt = -one(scalar_type)
    rates_opt = ntuple(i -> -one(fieldtype(rate_tuple_type, i)), ndims(rate_tuple_iter))
    nstates = get_nstates(graph)
    P_opt = Vector{A_type}(undef, nstates)
    nrate_tuples = length(rate_tuple_iter)
    _EYE_ = A_type(I)

    for (iter, rate_tuple) in enumerate(rate_tuple_iter)
        model = Model(optim_solver)

        P_list = [@variable(model, [1:A_size[1], 1:A_size[2]], Symmetric,
            base_name = string("P", q)) for q = 1:nstates]
        ee = @variable(model, base_name = "ee")

        for edge in enum_edges(graph)
            P1 = SMatrix{A_size...}(P_list[edge.source])
            P2 = SMatrix{A_size...}(P_list[edge.target])
            Ari_list = get(Ari_field, edge, Tuple{A_type,Int}[])
            for Ari in Ari_list
                A = Ari[1]
                r2 = rate_tuple[Ari[2]]^2
                @constraint(model, Symmetric(r2*P1 - A'*P2*A - ee.*_EYE_) in PSDCone())
            end
        end

        for q in 1:nstates
            P = SMatrix{A_size...}(P_list[q])
            @constraint(model, Symmetric(_EYE_ - P) in PSDCone())
            @constraint(model, Symmetric(_EYE_ + P) in PSDCone())
        end

        @objective(model, Max, ee)

        optimize!(model)

        if value(ee) > ee_opt
            ee_opt = value(ee)
            rates_opt = rate_tuple
            map!(x -> A_type(value.(x)), P_opt, P_list)
        end

        @printf("Iter: %d (/%d): rates = %s, ee: %g (ee*: %g)\n",
            iter, nrate_tuples, rate_tuple, value(ee), ee_opt)
    end

    return P_opt, ee_opt, rates_opt
end
