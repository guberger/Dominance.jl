## Solve dominance LMIs

function pth_eigval(A, p, shrink)
    e = sort(broadcast(abs, eigvals(A)), rev = true)
    return (e[p]*(1.0-shrink), e[p+1]/(1.0-shrink))
end

## Set

_matrix_type(::Type{MatrixSet{T,MT}}) where {T,MT} = MT

# ASri_lab : dict (edge -> [(matrix set, rate index)])
# lab for labelling
function cone_optim(graph, ASri_lab, rate_tuple_iter, optim_solver)
    ASri_type = eltype(valtype(ASri_lab))
    AS_type = fieldtype(ASri_type, 1)
    A_type = _matrix_type(AS_type)
    A_size = size(A_type)
    scalar_type = eltype(A_type)
    rate_tuple_type = eltype(rate_tuple_iter)
    δ_opt = -one(scalar_type)
    rates_opt = ntuple(i -> -one(fieldtype(rate_tuple_type, i)), ndims(rate_tuple_iter))
    nstates = get_nstates(graph)
    P_opt = Vector{A_type}(undef, nstates)
    nrate_tuples = length(rate_tuple_iter)
    _EYE_ = Symmetric(A_type(I))
    MAT_VAR = Symmetric{VariableRef, Matrix{VariableRef}}

    for (iter, rate_tuple) in enumerate(rate_tuple_iter)
        model = Model(optim_solver)

        P_list = [@variable(model, [1:A_size[1], 1:A_size[2]], Symmetric,
            base_name = string("P", q)) for q = 1:nstates]::Vector{MAT_VAR}
        E_list = Dict{Int,MAT_VAR}()
        radius_done = falses(nstates)
        ε = @variable(model, base_name = "ε")
        μ = @variable(model, base_name = "μ", lower_bound = 0.0)

        for edge in enum_edges(graph)
            q1 = edge.source
            q2 = edge.target
            P1 = P_list[q1]
            P2 = P_list[q2]
            ASri_list = get(ASri_lab, edge, ASri_type[])
            for ASri in ASri_list
                AS = ASri[1]
                r2 = rate_tuple[ASri[2]]^2
                Ac = AS.center
                AHull = AS.hull
                radius = AS.radius
                radius2 = radius^2
                if isempty(AHull)
                    continue
                end
                if radius > 0 && !radius_done[q2]
                    μbis = @variable(model)
                    for Avert in AHull
                        A = Ac + Avert
                        P2A = radius*P2*A
                        B = [μbis.*_EYE_ P2A; P2A' μbis.*_EYE_]
                        @constraint(model, Symmetric(B) ∈ PSDCone())
                    end
                    η = μ - 2*μbis
                    @constraint(model, Symmetric(η.*_EYE_ - radius2*P2) ∈ PSDCone())
                    radius_done[q2] = true
                end
                if length(AHull) == 1
                    A = Ac + AHull[1]
                    @constraint(model,
                        Symmetric(r2*P1 - A'*P2*A - ε.*_EYE_) ∈ PSDCone())
                    continue
                end
                if !haskey(E_list, q2)
                    E_list[q2] = @variable(model, [1:A_size[1], 1:A_size[2]], PSD,
                        base_name = string("E", q2))::MAT_VAR
                    E2 = E_list[q2]
                    @constraint(model, Symmetric(E2 - P2) ∈ PSDCone())
                end
                E2 = E_list[q2]
                for Avert in AHull
                    @constraint(model, Symmetric(r2*P1 - Ac'*P2*Ac - Avert'*P2*Ac
                        - Ac'*P2*Avert - Avert'*E2*Avert - ε.*_EYE_) ∈ PSDCone())
                end
            end
        end

        for q = 1:nstates
            @constraint(model, Symmetric(_EYE_ - P_list[q]) ∈ PSDCone())
            @constraint(model, Symmetric(_EYE_ + P_list[q]) ∈ PSDCone())
        end

        if !any(radius_done)
            println("Fix μ")
            fix(μ, 0.0, force = true)
        end

        δ = ε - μ
        @objective(model, Max, δ)

        println("Optimize:")
        optimize!(model)

        if value(δ) > δ_opt
            δ_opt = value(δ)
            rates_opt = rate_tuple
            map!(x -> A_type(value.(x)), P_opt, P_list)
        end

        @printf("Iter: %d (/%d): rates = %s, δ: %g (δ*: %g)\n",
            iter, nrate_tuples, rate_tuple, value(δ), δ_opt)
    end

    return P_opt, δ_opt, rates_opt
end

#=
function hyperbolic_optim(graph, A_lab, optim_solver)
    A_type = eltype(valtype(A_lab))
    A_size = size(A_type)
    scalar_type = eltype(A_type)
    δ_opt = -one(scalar_type)
    rates_opt = (1,)
    nstates = get_nstates(graph)
    P_opt = Vector{A_type}(undef, nstates)
    _EYE_ = A_type(I)

    model = Model(optim_solver)

    P_list = [@variable(model, [1:A_size[1], 1:A_size[2]], Symmetric,
        base_name = string("P", q)) for q = 1:nstates]
    ε = @variable(model, base_name = "ε")
    ε_EYE_ = ε.*_EYE_

    for edge in enum_edges(graph)
        P1 = P_list[edge.source]
        P2 = P_list[edge.target]
        A_list = get(A_lab, edge, A_type[])
        @constraint(model, [i = 1:length(A_list)], Symmetric(P1 -
            A_list[i]'*P2*A_list[i] - ε_EYE_) ∈ PSDCone())
    end

    for q = 1:nstates
        P = P_list[q]
        @constraint(model, Symmetric(_EYE_ - P) ∈ PSDCone())
        @constraint(model, Symmetric(_EYE_ + P) ∈ PSDCone())
    end

    @objective(model, Max, ε)

    optimize!(model)

    δ_opt = value(ε)
    map!(x -> A_type(value.(x)), P_opt, P_list)

    @printf("Iter: 1 (/1): rates = (1,), δ: %g (δ*: %g)\n", δ_opt, δ_opt)

    return P_opt, δ_opt, rates_opt
end
=#
