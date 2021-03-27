include("../abstraction.jl")
MABS = Abstraction
include("../system_models.jl")
using LinearAlgebra
using Printf

sleep(0.1) # used for good printing
println("\nStarted test")

AbsGenerator = newDamped1D()

theAbs = AbsGenerator()

(DDxNorm, xopt) = MABS.DDs_normInf2matp_rect(theAbs, Inf,
    [-5.0], [10.0], [10])
@printf("DDs Norm Inf: %f\n", DDxNorm)

AbsGenerator = newDamped2D(5, 5)

theAbs = AbsGenerator()

(DDxNorm, xopt) = MABS.DDs_normInf2matp_rect(theAbs, Inf,
    [-5.0, -5.0], [10.0, 10.0], [5, 5])
@printf("DDs Norm Inf: %f\n", DDxNorm)

AbsGenerator = newDamped3D()

theAbs = AbsGenerator()

(DDxNorm, xopt) = MABS.DDs_normInf2matp_rect(theAbs, Inf,
    [-5.0, -5.0, -5.0], [10.0, 10.0, 10.0], [5, 5, 5])
@printf("DDs Norm Inf: %f\n", DDxNorm)

# To be validated with matlab_test_DDsNorm.m
