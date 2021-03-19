using JuMP
using MosekTools
using LinearAlgebra
using Printf

function solve_lmi_disc_path(A_list, edge_list, rates_list)
    # Here each edge is a tuple (q1, q2, i, gamma)
    nNode = maximum(map(x -> max(x[1], x[2]), edge_list))
    nEdge = length(edge_list)
    dim = size(A_list[1], 1)
    nRates = length(rates_list)
    EYE_ = Matrix{Float64}(I, dim, dim)

    ee_max = -1.0
    rates_max = rates_list[1]
    P_max = Vector{Matrix{Float64}}(undef, nNode)

    for (iter, rates) in enumerate(rates_list)
        model = Model(optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true))

        P_list = [@variable(model, [1:dim, 1:dim], Symmetric,
            base_name = string("P", q)) for q = 1:nNode]
        ee = @variable(model, base_name = "ee")

        for edge in edge_list
            P1 = P_list[edge[1]]
            P2 = P_list[edge[2]]
            A = A_list[edge[3]]
            r = rates[edge[4]]^2
            @constraint(model, Symmetric(r*P1 - A'*P2*A - ee.*EYE_) in PSDCone())
        end

        @constraint(model, [q = 1:nNode], Symmetric(EYE_ - P_list[q]) in PSDCone())
        @constraint(model, [q = 1:nNode], Symmetric(EYE_ + P_list[q]) in PSDCone())

        @objective(model, Max, ee)

        optimize!(model)

        if value(ee) > ee_max
            ee_max = value(ee)
            rates_max = rates
            map!(x -> value.(x), P_max, P_list)
        end

        @printf("Iter: %d (/%d): rates = %s, ee: %g (ee_max: %g)\n",
            iter, nRates, rates, value(ee), ee_max)
    end

    return (P_max, ee_max, rates_max)
end

function solve_lmi_disc_path_convex(Ac_list, Ad_list, edge_list, rates_list)
    # Here each edge is a tuple (q1, q2, ic, [id], [i_gamma])
    nNode = maximum(map(x -> max(x[1], x[2]), edge_list))
    nEdge = length(edge_list)
    dim = size(Ac_list[1], 1)
    nRates = length(rates_list)
    EYE_ = Matrix{Float64}(I, dim, dim)

    ee_max = -1.0
    rates_max = rates_list[1]
    P_max = Vector{Matrix{Float64}}(undef, nNode)

    for (iter, rates) in enumerate(rates_list)
        model = Model(optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true))

        P_list = [@variable(model, [1:dim, 1:dim], Symmetric,
            base_name = string("P", q)) for q = 1:nNode]
        E_list = [@variable(model, [1:dim, 1:dim], PSD,
            base_name = string("E", d)) for d = 1:nEdge]
        ee = @variable(model, base_name = "ee")

        for (d, edge) in enumerate(edge_list)
            P1 = P_list[edge[1]]
            P2 = P_list[edge[2]]
            E = E_list[d]
            Ac = Ac_list[edge[3]]
            @constraint(model, Symmetric(E - P2) in PSDCone())
            @constraint(model, Symmetric(EYE_ - E) in PSDCone())
            for (id, ig) in zip(edge[4], edge[5])
                Ad = Ad_list[id]
                r = rates[ig]^2
                @constraint(model, Symmetric(r*P1 - Ac'*P2*Ac - Ac'*P2*Ad
                    - Ad'*P2*Ac - Ad'*E*Ad - ee.*EYE_) in PSDCone())
            end
        end

        @objective(model, Max, ee)

        optimize!(model)

        if value(ee) > ee_max
            ee_max = value(ee)
            rates_max = rates
            map!(x -> value.(x), P_max, P_list)
        end

        @printf("Iter: %d (/%d): rates = %s, ee: %g (ee_max: %g)\n",
            iter, nRates, rates, value(ee), ee_max)
    end

    return (P_max, ee_max, rates_max)
end

function solve_lmi_cont_path(A_list, edge_list, rates_list)
    nNode = maximum(map(x -> max(x[1], x[2]), edge_list))
    nEdge = length(edge_list)
    dim = size(A_list[1], 1)
    nRates = length(rates_list)
    EYE_ = Matrix{Float64}(I, dim, dim)

    ee_max = -1.0
    rates_max = rates_list[1]
    P_max = Vector{Matrix{Float64}}(undef, nNode)

    for (iter, rates) in enumerate(rates_list)
        model = Model(optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true))

        P_list = [@variable(model, [1:dim, 1:dim], Symmetric,
            base_name = string("P", q)) for q = 1:nNode]
        ee = @variable(model, base_name = "ee")

        for edge in edge_list
            P1 = P_list[edge[1]]
            P2 = P_list[edge[2]]
            A = A_list[edge[3]]
            r = 2.0*rates[edge[4]]
            @constraint(model, Symmetric(r*P1 - A'*P2 - P2*A - ee.*EYE_) in PSDCone())
        end

        @constraint(model, [q = 1:nNode], Symmetric(EYE_ - P_list[q]) in PSDCone())
        @constraint(model, [q = 1:nNode], Symmetric(EYE_ + P_list[q]) in PSDCone())

        @objective(model, Max, ee)

        optimize!(model)

        if value(ee) > ee_max
            ee_max = value(ee)
            rates_max = rates
            map!(x -> value.(x), P_max, P_list)
        end

        @printf("Iter: %d (/%d): rates = %s, ee: %g (ee_max: %g)\n",
            iter, nRates, rates, value(ee), ee_max)
    end

    return (P_max, ee_max, rates_max)
end

function ndgrid_array(vec_array)
    vec_iter = Iterators.product(vec_array...)
    return vec(collect(vec_iter))
end

function linspace_array(bounds_list, ndisc)
    return [range(bounds_list[i][1], stop = bounds_list[i][2], length = ndisc) for i = 1:length(bounds_list)]
end

function pth_eigen(A, p, shrink)
    e = sort(broadcast(abs, eigvals(A)), rev = true)
    return (e[p]*(1.0-shrink), e[p+1]/(1.0-shrink))
end
