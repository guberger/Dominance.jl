module Dominance

using LinearAlgebra
using StaticArrays
using JuMP
using Printf

@enum INCL_MODE INNER OUTER

include("mapping.jl")
include("matrices.jl")
include("tensors.jl")
include("graph.jl")
include("sub_graph.jl")
include("system.jl")
include("set.jl")
include("abstraction.jl")
include("refine_domain.jl")
include("minimum_domain.jl")
include("symbolic_model.jl")
include("viable_states.jl")
include("cone_optim.jl")
include("matrix_field.jl")
end  # module Dominance
