include("../src/symfunc.jl")

# Build Abstraction
tstep = 1.0
mu = 2.0
mu = 0.1
# Symbolic computations
x = SymPy.Sym("x")
y = SymPy.Sym("y")
# f(x,y): R^2 to R^2
fx = mu*(1.0 - y^2)*x - y
fy = x
symfunc_2_juliafunc((fx, fy), (x, y), base_name = "Vanderpol")
