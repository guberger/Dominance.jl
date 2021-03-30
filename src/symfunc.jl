# Symbolic manipulations

using SymPy

function symfunc_2_juliafunc(sym_func, sym_vars; base_name = "F_float")
    dim = length(sym_vars)
    xlist = [SymPy.Sym(string("x[", i, ']')) for i = 1:dim]
    subsmap = [sym_vars[i] => xlist[i] for i = 1:dim]
    F_Sym = [sym_func[i](subsmap...) for i =1:dim]
    DF_Sym = [diff(F_Sym[i], xlist[j]) for i = 1:dim, j = 1:dim]
    DDF_Sym = [diff(F_Sym[i], xlist[j], xlist[k]) for i = 1:dim, j = 1:dim, k = 1:dim]

    F_string = [string(F_Sym[i], ", ") for i = 1:dim]

    S = string("function $(base_name)(x)\n", "    ",
        "return SVector(", F_string..., ")\nend")
    println(S)
    ss = Meta.parse(S)
    eval(ss)

    DF_string = [string(DF_Sym[i, j], ", ") for i = 1:dim, j = 1:dim]

    S = string("function D$(base_name)(x)\n", "    ",
        "return SMatrix{$(dim),$(dim)}(", DF_string..., ")\nend")
    println(S)
    ss = Meta.parse(S)
    eval(ss)

    DDF_string = [string(DDF_Sym[i, j, k], ", ") for i = 1:dim, j = 1:dim, k = 1:dim]

    S = string("function DD$(base_name)(x)\n", "    ",
        "return SArray{Tuple{$(dim),$(dim),$(dim)}}(", DDF_string..., ")\nend")
    println(S)
    ss = Meta.parse(S)
    eval(ss)
end
