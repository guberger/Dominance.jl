module Dominance

using LinearAlgebra
using StaticArrays

@enum INCL_MODE INNER OUTER

include("Graph.jl")
include("system.jl")
include("abstraction.jl")
include("macros.jl")

end  # module Dominance
