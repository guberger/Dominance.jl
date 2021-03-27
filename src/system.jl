struct System{F1<:Function,F2<:Function,F3<:Function}
    sys_map::F1
    linsys_map::F2
    error_map::F3
end

function RungeKutta4(F, x, tstep, nsub::Int)
    τ = tstep/nsub
    for i in 1:nsub
        Fx1 = F(x)
        xrk = x + Fx1*(τ/2)
        Fx2 = F(xrk)
        xrk = x + Fx2*(τ/2)
        Fx3 = F(xrk)
        xrk = x + Fx3*τ
        Fx4 = F(xrk)
        x = x + (Fx1 + Fx2*2 + Fx3*2 + Fx4)*(τ/6)
    end
    return x
end

function RungeKutta4Linearized(F, DF, x, dx, tstep, nsub::Int)
    τ = tstep/nsub
    for i in 1:nsub
        Fx1 = F(x)
        DFx1 = DF(x)*dx
        xrk = x + Fx1*(τ/2)
        dxrk = dx + DFx1*(τ/2)
        Fx2 = F(xrk)
        DFx2 = DF(xrk)*dxrk
        xrk = x + Fx2*(τ/2)
        dxrk = dx + DFx2*(τ/2)
        Fx3 = F(xrk)
        DFx3 = DF(xrk)*dxrk
        xrk = x + Fx3*τ
        dxrk = dx + DFx3*τ
        Fx4 = F(xrk)
        DFx4 = DF(xrk)*dxrk
        x += (Fx1 + Fx2*2 + Fx3*2 + Fx4)*(τ/6)
        dx += (DFx1 + DFx2*2 + DFx3*2 + DFx4)*(τ/6)
    end
    return (x, dx)
end

# Give an upper-bound on x(t) staisfying x'(t) ≦ a*x(t) + b*exp(2at)
function BoundSecondOrder(a, b, tstep)
    if abs(a) < 1e-8
        return b*tstep
    else
        ρ = exp(a*tstep)
        return (b/a)*ρ*(ρ - 1)/2
    end
end

# Bounds with respect to \infty-norm
function ContSystemRK4(tstep, F_sys, DF_sys, bound_DF, bound_DDF, nsys)
    sys_map = let tstep = tstep, nsys = nsys
        x -> RungeKutta4(F_sys, x, tstep, nsys)
    end
    linsys_map = let tstep = tstep, nsys = nsys
        (x, dx) -> RungeKutta4Linearized(F_sys, DF_sys, x, dx, tstep, nsys)
    end
    error_map = let tstep = tstep, bound_DF = bound_DF, bound_DDF = bound_DDF
        r -> BoundSecondOrder(bound_DF, bound_DDF, tstep)*(r^2)
    end
    return System(sys_map, linsys_map, error_map)
end

# Bounds with respect to \infty-norm
function DiscSystem(F_sys, DF_sys, bound_DDF)
    sys_map = F_sys
    linsys_map = (x, dx) -> (sys_map(x), DF_sys(x)*dx)
    error_map = let bound_DDF = bound_DDF
        r -> bound_DDF*(r^2)/2
    end
    return System(sys_map, linsys_map, error_map)
end
