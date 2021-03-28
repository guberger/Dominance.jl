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
include("macros.jl")

end  # module Dominance
