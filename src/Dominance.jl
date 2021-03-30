module Dominance

using LinearAlgebra
using StaticArrays
using JuMP
using Printf

@enum INCL_MODE INNER OUTER

include("matrices.jl")
include("tensors.jl")
include("set.jl")

include("system.jl")

include("graph.jl")
include("abstraction.jl")

# macros
include("symbolic_model_from_system.jl")
include("minimize_over_domain.jl")
include("cone_optim.jl")

end  # module Dominance
