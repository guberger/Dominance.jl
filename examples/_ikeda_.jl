include("../src/symfunc.jl")

# Parameters; see Osipenko (2006)
r = 2.0
a = -0.9
b = 0.9
C1 = 0.4
C3 = 6.0

# Symbolic computations
x = Sym("x")
y = Sym("y")
tau = C1 - C3 / (1.0 + x^2 + y^2)
fx = r + a * (x * cos(tau) - y * sin(tau))
fy = b * (x * sin(tau) + y * cos(tau))
symfunc_2_juliafunc((fx, fy), (x, y), base_name = "Ikeda")



function newikeda(nx, ny)
    # Build Abstraction
    df_lip = 22.5
    xmin = -1.1
    xmax = 3.4
    ymin = -1.5
    ymax = 1.8
    hx = (xmax-xmin) / nx
    hy = (ymax-ymin) / ny
    xsw = xmin
    ysw = ymin
    nidx = nx * ny
    cell_array = Vector{Abstraction.Cell}(undef, nidx)

    theAbs = Abstraction.AbstractionType(
        fx_float, fy_float, DF_MAT_float, DDF_MAT_float, df_lip,
        cell_array, xsw, ysw, hx, hy,
        Vector{Int}[], false,
        Vector{Array{Float64, 2}}[], false)

    # Initialize cell_array
    for i = 1:nx
        for j = 1:ny
            idx = i + (j-1) * nx
            theAbs.cell_array[idx] = Abstraction.Cell(i-1, j-1, 1, 1)
        end
    end

    return theAbs
end
