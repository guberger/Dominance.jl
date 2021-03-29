using JuMP
using MosekTools
using LinearAlgebra

model = Model(optimizer_with_attributes(Mosek.Optimizer))
P = @variable(model, [1:2, 1:2], PSD)
C = @variable(model, [1:2], lower_bound = 0.0, base_name = "rt")
@objective(model, Max, -P[2, 2])
optimize!(model)
