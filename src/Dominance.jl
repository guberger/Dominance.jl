module Dominance

using LinearAlgebra
using StaticArrays

@enum INCL_MODE INNER OUTER

include("mapping.jl")
include("matrices.jl")
include("graph.jl")
include("system.jl")
include("set.jl")
include("abstraction.jl")
include("symbolic_model.jl")
include("viable_states.jl")

end  # module Dominance
