# ODE, linearized and second-order systems

@inline function linearizedsystem_vector!(F!, DF!, u, du, dim)
    # Vector field of the linearized system, with solution (xsi, dxsi/dx) vectorized.
    x = @view u[1:dim]
    Dx = @view u[dim+1:dim+dim^2]
    dx = @view du[1:dim]
    dDx = @view du[dim+1:dim+dim^2]
    DFx = Matrix{Float64}(undef, dim, dim)
    F!(x, dx)
    DF!(x, DFx)
    # Faster than reshape and than using three loops (tested).
    for j = 1:dim
        dDx[1+(j-1)*dim:j*dim] = DFx*Dx[1+(j-1)*dim:j*dim]
    end
end

function linearizedsystem_map(linsys_vector!, x, tstep, reltol, abstol, dim)
    # Computes the solution of the linearized system starting from (x, I) and at time tstep.
    u0 = vcat(x, vec(Matrix{Float64}(I, dim, dim)))
    prob = ODEProblem(linsys_vector!, u0, (0.0, tstep))
    sol = solve(prob, Tsit5(), reltol = reltol, abstol = abstol, saveat = [tstep])
    return (sol[1][1:dim], reshape(sol[1][dim+1:dim+dim^2], dim, dim))
end

function set_linearizedsystem!(theAbs)
    @inline linsys_vector!(du, u, p_dummy, t) =
        linearizedsystem_vector!(theAbs.F!, theAbs.DF!, u, du, theAbs.dim)
    linsys_map(x) = linearizedsystem_map(linsys_vector!, x,
        theAbs.tstep, theAbs.ode_reltol, theAbs.ode_abstol, theAbs.dim)
    theAbs.linsys_vector! = linsys_vector!
    theAbs.linsys_map = linsys_map
end

@inline function secondordersystem_vector!(F!, DF!, DDF!, u, du, dim)
    # Vector field for xsi, dxsi/dx, and d^2dxsi/dx^2.
    i0 = dim + dim^2
    x = @view u[1:dim]
    Dx = @view u[dim+1:dim+dim^2]
    DDx = @view u[i0+1:i0+dim^3]
    DxResh = reshape(Dx, dim, dim)
    DDxResh = reshape(DDx, dim, dim, dim)
    dx = @view du[1:dim]
    dDx = @view du[dim+1:dim+dim^2]
    dDDx = @view du[i0+1:i0+dim^3]
    DFx = Matrix{Float64}(undef, dim, dim)
    DDFx = Array{Float64, 3}(undef, dim, dim, dim)
    F!(x, dx)
    DF!(x, DFx)
    DDF!(x, DDFx)

    # Faster than reshape and than using three loops (tested).
    for j = 1:dim
        dDx[1+(j-1)*dim:j*dim] = DFx*Dx[1+(j-1)*dim:j*dim]
    end
    for i = 1:dim
        DDFx[i, :, :] = DxResh'*DDFx[i, :, :]*DxResh
    end
    for k = 1:dim
        DDFx[:, :, k] = DDFx[:, :, k] + DFx*DDxResh[:, :, k]
    end

    dDDx[:] = vec(DDFx)
end

function secondordersystem_map(secordsys_vector!, x, tstep, reltol, abstol, dim)
    # Computes the solution of the 2nd-order system starting from (x, I, 0) and at time tstep.
    u0 = vcat(x, vec(Matrix{Float64}(I, dim, dim)), zeros(dim^3))
    prob = ODEProblem(secordsys_vector!, u0, (0.0, tstep))
    sol = solve(prob, Tsit5(), reltol = reltol, abstol = abstol, saveat = [tstep])
    i0 = dim + dim^2
    return (sol[1][1:dim],
        reshape(sol[1][dim+1:dim^2+dim], dim, dim),
        reshape(sol[1][i0+1:i0+dim^3], dim, dim, dim))
end
