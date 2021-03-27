module Abstraction

using LinearAlgebra
using DifferentialEquations
using Printf
using Combinatorics
using PyPlot
using PyCall
const art3d = PyObject(PyPlot.art3D)
const spatial = pyimport_conda("scipy.spatial", "scipy")
using Optim

struct CellType
    coords::Vector{Int}
    Lev::Vector{Int}
end

mutable struct AbstractionType
    dim::Int
    tstep::Float64
    ode_reltol::Float64
    ode_abstol::Float64
    set_reltol::Float64
    set_abstol::Float64
    F!::Function
    DF!::Function
    DDF!::Function
    linsys_vector!::Function
    linsys_map::Function
    x_or::Vector{Float64}
    h::Vector{Float64}
    cell_list::Vector{CellType}
    edge_list::Vector{Vector{Int}}
    nsubim_list::Vector{Vector{Int}}
    # || dxsi/dx(tstep, x) - dxsi/dx(tstep, y) ||_{inf, inf} <= df_lip || x-y ||_inf
    dslipInf_list::Vector{Float64}
    nsubopt_list::Vector{Vector{Int}}
    DS_list::Vector{Vector{Matrix{Float64}}}
    dslip2_list::Vector{Float64}
    nsubds_list::Vector{Vector{Int}}
    ############################################################################
    # Constructor
    function AbstractionType(tstep, F!, DF!, DDF!, x_or, h, cell_list;
            ode_abstol = 1e-8, ode_reltol = 1e-8,
            set_abstol = 1e-10, set_reltol = 1e-10)
        #-----------------------------------------------------------------------
        ncell = length(cell_list)
        dim = length(x_or)
        @assert length(h) == dim
        @inline nsub_list() = fill(fill(-1, dim), ncell) # "copy"
        @inline dslip_list() = fill(-1.0, ncell)
        DS_list = fill(Matrix{Float64}[], ncell)

        theAbs = new(dim, tstep, ode_reltol, ode_abstol, set_reltol, set_abstol,
            F!, DF!, DDF!, x -> x, x -> x,
            x_or, h, cell_list, Vector{Int}[], nsub_list(), dslip_list(),
            nsub_list(), DS_list, dslip_list(), nsub_list())
        set_linearizedsystem!(theAbs)

        return theAbs
    end
    ############################################################################
end

include("interface.jl")
include("differential.jl")
include("deviation.jl")
include("images.jl")
include("macros.jl")
include("plotting.jl")
include("dominance.jl")
# include("hyperbolic.jl")
# include("saveload.jl")

end  # module Abstraction
