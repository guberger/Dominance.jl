include("../src/symfunc.jl")

x = SymPy.Sym("x")
y = SymPy.Sym("y")
sym_func = [x^2 - y^2, 2*x*y]*((x^2 + y^2)^(-3/4))
# sym_func = 2*atan.([x^2 - y^2, 2*x*y])
symfunc_2_juliafunc(sym_func, (x, y), base_name = "Doubling")
