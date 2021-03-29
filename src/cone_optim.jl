## Solve dominance LMIs

function pth_eigval(A, p, shrink)
    e = sort(broadcast(abs, eigvals(A)), rev = true)
    return (e[p]*(1.0-shrink), e[p+1]/(1.0-shrink))
end

## Single

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

        for q = 1:nstates
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

## Set

_matrix_type(::Type{MatrixSet{T,MT}}) where {T,MT} = MT

# ASri_field : dict (edge -> [(matrix set, rate index)])
function cone_optim_set(graph, ASri_field, rate_tuple_iter, optim_solver)
    AS_type = fieldtype(eltype(valtype(ASri_field)), 1)
    A_type = _matrix_type(AS_type)
    A_size = size(A_type)
    scalar_type = eltype(A_type)
    rate_tuple_type = eltype(rate_tuple_iter)
    ee_opt = -one(scalar_type)
    rates_opt = ntuple(i -> -one(fieldtype(rate_tuple_type, i)), ndims(rate_tuple_iter))
    nstates = get_nstates(graph)
    P_opt = Vector{A_type}(undef, nstates)
    nrate_tuples = length(rate_tuple_iter)
    _EYE_ = A_type(I)
    VAR_TYPE = Symmetric{VariableRef, Matrix{VariableRef}}

    for (iter, rate_tuple) in enumerate(rate_tuple_iter)
        model = Model(optim_solver)

        P_list = [@variable(model, [1:A_size[1], 1:A_size[2]], Symmetric,
            base_name = string("P", q)) for q = 1:nstates]
        E_list = Dict{Int,VAR_TYPE}()
        ee = @variable(model, base_name = "ee")

        for edge in enum_edges(graph)
            q1 = edge.source
            q2 = edge.target
            P1 = SMatrix{A_size...}(P_list[q1])
            P2 = SMatrix{A_size...}(P_list[q2])
            ASri_list = get(ASri_field, edge, Tuple{A_type,Int}[])
            for ASri in ASri_list
                AS = ASri[1]
                r2 = rate_tuple[ASri[2]]^2
                Ac = AS.center
                AHull = AS.hull
                if isempty(AHull)
                    continue
                end
                if length(AHull) == 1
                    A = Ac + AHull[1]
                    @constraint(model, Symmetric(r2*P1 - A'*P2*A - ee.*_EYE_) in PSDCone())
                    continue
                end
                if !haskey(E_list, q2)
                    E_list[q2] = @variable(model, [1:A_size[1], 1:A_size[2]], PSD,
                        base_name = string("E", q2))
                    E2 = SMatrix{A_size...}(E_list[q2])
                    @constraint(model, Symmetric(E2 - P2) in PSDCone())
                end
                E2 = SMatrix{A_size...}(E_list[q2])
                for Avert in AHull
                    @constraint(model, Symmetric(r2*P1 - Ac'*P2*Ac + Avert'*P2*Ac
                        + Ac'*P2*Avert + Avert'*E2*Avert - ee.*_EYE_) in PSDCone())
                end
            end
        end

        for q = 1:nstates
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
