include("symbolic.jl")

function make_chessboard(cell_list, n1, n2)
    for i = 1:n1
        for j = 1:n2
            if mod(i+j, 2) == 0
                push!(cell_list, Abstraction.CellType([i-1, j-1], [1, 1]))
            end
        end
    end
end

function newDamped1D()
    tstep = 2.0
    df_lip = 0.0
    x_or = [-1.0]
    h = [0.5]
    cell_list = [Abstraction.CellType([1], [1])]

    # Symbolic computations
    x = SymPy.Sym("x")
    # f(x): R to R
    F_Sym = [1.0/(1.0 + x^2)]

    sympy_2_juliafunction(F_Sym, [x])

    function newAbstraction()
        theAbs = Abstraction.AbstractionType(tstep,
            F_float!, DF_float!, DDF_float!,
            x_or, h, cell_list)
        return theAbs
    end

    return newAbstraction
end

function newDamped2D(n1, n2)
    tstep = 2.0
    df_lip = 0.0
    x_or = [-1.0, -1.0]
    h = [0.5, 0.5]
    cell_list = Abstraction.CellType[]
    make_chessboard(cell_list, n1, n2)

    # Symbolic computations
    x = SymPy.Sym("x")
    y = SymPy.Sym("y")
    # f(x,y): R^2 to R^2
    F_Sym = [1.0/(1.0 + x^2), x^2]

    sympy_2_juliafunction(F_Sym, [x, y])

    function newAbstraction()
        theAbs = Abstraction.AbstractionType(tstep,
            F_float!, DF_float!, DDF_float!,
            x_or, h, cell_list)
        return theAbs
    end

    return newAbstraction
end

function newDamped3D()
    tstep = 2.0
    df_lip = 0.0
    x_or = [-1.0, -1.0, -1.0]
    h = [0.5, 0.5, 0.5]
    cell_list = [Abstraction.CellType([1, 1, 1], [1, 1, 1])]

    # Symbolic computations
    x = SymPy.Sym("x")
    y = SymPy.Sym("y")
    z = SymPy.Sym("z")
    # f(x,y): R^2 to R^2
    F_Sym = [1.0/(1.0 + x^2), x^2, 2.0*y]

    sympy_2_juliafunction(F_Sym, [x, y, z])

    function newAbstraction()
        theAbs = Abstraction.AbstractionType(tstep,
            F_float!, DF_float!, DDF_float!,
            x_or, h, cell_list)
        return theAbs
    end

    return newAbstraction
end

function newNonlinearOscillator2D(n1, n2)
    # Build Abstraction
    tstep = 1.0
    df_lip = 4.362169
    x_or = [-1.0, -1.0]
    h = [0.5, 0.6]
    cell_list = Abstraction.CellType[]
    make_chessboard(cell_list, n1, n2)

    c = 0.25
    u = 2.0
    omega = 1.0
    # Symbolic computations
    x = SymPy.Sym("x")
    y = SymPy.Sym("y")
    # f(x,y): R^2 to R^2
    F_Sym = [c*(u^2-x^2-y^2)*x - omega*y, omega*x]

    sympy_2_juliafunction(F_Sym, [x, y])

    function newAbstraction()
        theAbs = Abstraction.AbstractionType(tstep,
            F_float!, DF_float!, DDF_float!,
            x_or, h, cell_list)
        return theAbs
    end

    return newAbstraction
end

function newNonlinearOscillator3D()
    # Build Abstraction
    tstep = 1.0
    x_or = [-3.0, -3.0, -3.0]*2.25
    h = [0.5, 0.5, 0.5]
    cell_list =[Abstraction.CellType([1, 1, 1], [1, 1, 1])]

    c = 0.25
    u = 2.0
    omega = 1.0
    # Symbolic computations
    x = SymPy.Sym("x")
    y = SymPy.Sym("y")
    z = SymPy.Sym("z")
    # f(x,y): R^2 to R^2
    F_Sym = [c*(u^2-x^2-y^2)*x - omega*y, omega*x, -z]

    sympy_2_juliafunction(F_Sym, [x, y, z])

    function newAbstraction()
        theAbs = Abstraction.AbstractionType(tstep,
            F_float!, DF_float!, DDF_float!,
            x_or, h, cell_list)
        return theAbs
    end

    return newAbstraction
end

function newVanderpol()
    # Build Abstraction
    tstep = 1.0
    # tstep = 7.6299
    # df_lip = 548.42
    x_or = [-1.0, -1.0]
    h = [0.3, 0.3]
    cell_list = [Abstraction.CellType([0, 0], [1, 1])]

    mu = 2.0
    # Symbolic computations
    x = SymPy.Sym("x")
    y = SymPy.Sym("y")
    # f(x,y): R^2 to R^2
    F_Sym = [mu*(1.0-y^2)*x - y, x]

    sympy_2_juliafunction(F_Sym, [x, y])

    function newAbstraction()
        theAbs = Abstraction.AbstractionType(tstep,
            F_float!, DF_float!, DDF_float!,
            x_or, h, cell_list)
        return theAbs
    end

    return newAbstraction
end

function newDuffingControlled3D()
    # Build Abstraction
    tstep = 1.5
    x_or = [0.5, 0.5, 0.0]
    h = [0.5, 0.5, 0.5]
    cell_list =[Abstraction.CellType([0, 0, 0], [1, 1, 1])]

    delta = 0.3
    alpha = -1.0
    beta = 1.0
    kF = 0.2
    yref = 0.2
    kX = 1.0
    kV = 0.0
    R = 0.3
    L = 1.0
    # Symbolic computations
    x = SymPy.Sym("x")
    y = SymPy.Sym("y")
    z = SymPy.Sym("z")
    # f(x,y): R^2 to R^2
    F_Sym = [-delta*x - alpha*y - beta*y^3 + kF*z,
        x,
        (kX/L)*(yref - y) - (kV/L)*x - (R/L)*z]
    sympy_2_juliafunction(F_Sym, [x, y, z])

    function newAbstraction()
        theAbs = Abstraction.AbstractionType(tstep,
            F_float!, DF_float!, DDF_float!,
            x_or, h, cell_list)
        return theAbs
    end

    return newAbstraction
end

#=
function newDamped1D()
    tstep = 2.0
    df_lip = 0.0
    x_or = [-1.0]
    h = [0.5]
    cell_list = [Abstraction.CellType([1], [1])]

    # Symbolic computations
    x = SymPy.Sym("x")
    # f(x): R to R
    F_Sym = [1.0/(1.0 + x^2)]

    Abstraction.sympy_2_juliafunction(F_Sym, [x])

    theAbs = Abstraction.AbstractionType(tstep,
        Abstraction.F_float!, Abstraction.DF_float!, Abstraction.DDF_float!,
        df_lip, x_or, h, cell_list)

    return theAbs
end

function newDamped2D(n1, n2)
    tstep = 2.0
    df_lip = 0.0
    x_or = [-1.0, -1.0]
    h = [0.5, 0.5]
    cell_list = Abstraction.CellType[]
    make_chessboard(cell_list, n1, n2)

    # Symbolic computations
    x = SymPy.Sym("x")
    y = SymPy.Sym("y")
    # f(x,y): R^2 to R^2
    F_Sym = [1.0/(1.0 + x^2), x^2]

    Abstraction.sympy_2_juliafunction(F_Sym, [x, y])

    theAbs = Abstraction.AbstractionType(tstep,
        Abstraction.F_float!, Abstraction.DF_float!, Abstraction.DDF_float!,
        df_lip, x_or, h, cell_list)

    return theAbs
end

function newDamped3D()
    tstep = 2.0
    df_lip = 0.0
    x_or = [-1.0, -1.0, -1.0]
    h = [0.5, 0.5, 0.5]
    cell_list = [Abstraction.CellType([1, 1, 1], [1, 1, 1])]

    # Symbolic computations
    x = SymPy.Sym("x")
    y = SymPy.Sym("y")
    z = SymPy.Sym("z")
    # f(x,y): R^2 to R^2
    F_Sym = [1.0/(1.0 + x^2), x^2, 2*y]

    Abstraction.sympy_2_juliafunction(F_Sym, [x, y, z])

    theAbs = Abstraction.AbstractionType(tstep,
        Abstraction.F_float!, Abstraction.DF_float!, Abstraction.DDF_float!,
        df_lip, x_or, h, cell_list)

    return theAbs
end

function newLinearOscillator(n1, n2)
    # Build Abstraction
    tstep = 1.0
    df_lip = 0.0
    x_or = [-1.0, -1.0]
    h = [0.5, 0.6]
    cell_list = Abstraction.CellType[]
    make_chessboard(cell_list, n1, n2)

    c = 0.25
    omega = 1.0
    # Symbolic computations
    x = SymPy.Sym("x")
    y = SymPy.Sym("y")
    # f(x,y): R^2 to R^2
    F_Sym = [-c*x - omega*y, omega*x]

    Abstraction.sympy_2_juliafunction(F_Sym, [x, y])

    theAbs = Abstraction.AbstractionType(tstep,
        Abstraction.F_float!, Abstraction.DF_float!, Abstraction.DDF_float!,
        df_lip, x_or, h, cell_list)

    return theAbs
end

function newNonlinearOscillator2D(n1, n2)
    # Build Abstraction
    tstep = 1.0
    df_lip = 4.362169
    x_or = [-1.0, -1.0]
    h = [0.5, 0.6]
    cell_list = Abstraction.CellType[]
    make_chessboard(cell_list, n1, n2)

    c = 0.25
    u = 2.0
    omega = 1.0
    # Symbolic computations
    x = SymPy.Sym("x")
    y = SymPy.Sym("y")
    # f(x,y): R^2 to R^2
    F_Sym = [c*(u^2-x^2-y^2)*x - omega*y, omega*x]

    Abstraction.sympy_2_juliafunction(F_Sym, [x, y])

    theAbs = Abstraction.AbstractionType(tstep,
        Abstraction.F_float!, Abstraction.DF_float!, Abstraction.DDF_float!,
        df_lip, x_or, h, cell_list)

    return theAbs
end
=#
