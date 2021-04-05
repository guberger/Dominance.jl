include("../src/symfunc.jl")

x = SymPy.Sym("x")
y = SymPy.Sym("y")
θ = π/5.0
U = 2*SMatrix{2,2}(cos(θ), -sin(θ), sin(θ), cos(θ))
sym_func = U*[atan(x), atan(y)]
symfunc_2_juliafunc(sym_func, (x, y), base_name = "Contracting")
