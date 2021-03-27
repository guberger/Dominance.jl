# Symbolic manipulations
import SymPy

function sympy_2_juliafunction(F_Sym, vars)
    # Given a symbolic function, creates Julia functions that provide the first
    # and second derivatives of F_Sym, using symbolic computations (SymPy).
    dim = length(vars)
    xlist = Vector{Any}(undef, dim)

    for i = 1:dim
        xlist[i] = SymPy.Sym(string("x[", i, ']'))
    end

    subsmap = [(vars[i]=>xlist[i]) for i = 1:dim]
    F_Sym = [F_Sym[i](subsmap...) for i = 1:dim]
    DF_Sym = [diff(F_Sym[i], u) for i = 1:dim, u in xlist]
    DDF_Sym = [diff(F_Sym[i], u1, u2) for i = 1:dim, u1 in xlist, u2 in xlist]

    S = string([string("\tFx[", i, "] = ", F_Sym[i], "\n") for i = 1:dim]...)
    S = string("@inline function F_float!(x, Fx)\n", S, "end")
    println(S)
    ss = Meta.parse(S)
    eval(ss)

    S = string([string("\tDFx[", i + (j-1)*dim, "] = ", DF_Sym[i, j], "\n")
        for i = 1:dim, j = 1:dim]...)
    S = string("@inline function DF_float!(x, DFx)\n", S, "end")
    println(S)
    ss = Meta.parse(S)
    eval(ss)

    S = string([string("\tDDFx[", i + (j-1)*dim + (k-1)*dim*dim,
        "] = ", DDF_Sym[i, j, k], "\n")
        for i = 1:dim, j = 1:dim, k = 1:dim]...)
    S = string("@inline function DDF_float!(x, DDFx)\n", S, "end")
    println(S)
    ss = Meta.parse(S)
    eval(ss)
end
