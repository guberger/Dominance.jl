using LinearAlgebra
using Combinatorics

include("../abstraction.jl")

println("Hello")

dim = 2
col_idx = combinations(1:2*dim, dim-1)
vert_list = Iterators.product(fill(-1:2:1, dim)...)
A_lineq1 = zeros(length(col_idx), dim)
A_lineq2 = zeros(length(col_idx), dim)
A1 = rand(dim, dim)
A2 = rand(dim, dim)
reltol = 1e-8
abstol = 1e-8

np = 500

function test()
    @time for i = 1:np
        Abstraction.augmented_linearbox!(A_lineq1, A1, A2, col_idx, vert_list, dim, reltol, abstol)
    end
    @time for i = 1:np
        Abstraction.augmented_linearbox_old!(A_lineq2, A1, A2, col_idx, vert_list, dim)
    end
    # -> Old is faster !!!
    display(A_lineq1)
    display(A_lineq2)
end

sleep(0.1) # used for good printing
println("\nHey")
test()
