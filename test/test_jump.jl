using JuMP
using MosekTools

model = Model(optimizer_with_attributes(Mosek.Optimizer))
P = @variable(model, [1:2, 1:2], PSD)
C = @variable(model)
@objective(model, Max, -P[2, 2])
fix(C, 5.0)
optimize!(model)
