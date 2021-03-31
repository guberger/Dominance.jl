module TestMain

using Test
using LinearAlgebra
using StaticArrays
@static if isdefined(Main, :TestLocal)
    include("../src/Dominance.jl")
else
    using Dominance
end
DO = Dominance

sleep(0.1) # used for good printing
println("Started test")

@testset "Minimize over domain" begin
lb = SVector(-7.0, -7.0)
ub = SVector(7.0, 7.0)
x0 = SVector(0.01, 0.01)
h = SVector(1.0, 2.0)/10
grid = DO.Grid(x0, h)
domain = DO.Domain(grid)
DO.add_set!(domain, DO.HyperRectangle(lb, ub), DO.OUTER)

f(x) = -1/((x[1] - 1)^2 + (x[2] - 2)^2 + 1)
nsub = (3, 3)
f_opt, x_opt = DO.minimize_over_domain(f, domain, nsub)
@test f_opt < -0.99
@test norm(x_opt - SVector(1.0, 2.0), Inf) <= norm(h, Inf)/minimum(nsub)

θ = π/5.0
U = 2*SMatrix{2,2}(cos(θ), -sin(θ), sin(θ), cos(θ))
F_sys(x) = U*SVector(atan(x[1]), atan(x[2]))
DF_sys(x) = U*SMatrix{2,2}(1/(1 + x[1]^2), 0, 0, 1/(1 + x[2]^2))
DDF_sys(x) = SArray{Tuple{2,2,2}}(-3.23606797749979*x[1]/(x[1]^2 + 1)^2,
    2.35114100916989*x[1]/(x[1]^2 + 1)^2, 0, 0, 0, 0,
    -2.35114100916989*x[2]/(x[2]^2 + 1)^2,
    -3.23606797749979*x[2]/(x[2]^2 + 1)^2, )
f(x) = -DO.tensor3d_normInf2matp(DDF_sys(x), Inf)
f_opt, x_opt = DO.minimize_over_domain(f, domain, nsub)
bound_DDF = opnorm(U, Inf)*3*sqrt(3)/8
display(f_opt)
@test f_opt < -bound_DDF*0.99
f(x) = -DO.tensor3d_normInf2matp(DDF_sys(x), 2)
f_opt, x_opt = DO.minimize_over_domain(f, domain, nsub)
bound_DDF = opnorm(U, 2)*3*sqrt(3)/8
@test f_opt < -bound_DDF*0.99
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
