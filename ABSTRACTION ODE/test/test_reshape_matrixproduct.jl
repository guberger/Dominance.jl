# Test perf of matrix-matrix product A*B when B is vectorized

function meth_reshape!(A, b, c, dim)
    bb = reshape(b, (dim, dim))
    c[:] = vec(A*bb)
end

function meth_loop!(A, b, c, dim)
    for i = 1:dim
        for j = 1:dim
            c[i + (j-1)*dim] = sum(A[i, :] .* b[(j-1)*dim+1:j*dim])
        end
    end
end

function meth_loopmat!(A, b, c, dim)
    for j = 1:dim
        c[(j-1)*dim+1:j*dim] = A*b[(j-1)*dim+1:j*dim]
    end
end

function test()
    dim = 7
    A = rand(dim, dim)
    b = rand(dim*dim)
    c1 = Vector{Float64}(undef, dim*dim)
    c2 = Vector{Float64}(undef, dim*dim)
    c3 = Vector{Float64}(undef, dim*dim)

    @time for i = 1:100
        meth_reshape!(A, b, c1, dim)
    end
    @time for i = 1:100
        meth_loop!(A, b, c2, dim)
    end
    @time for i = 1:100
        meth_loopmat!(A, b, c3, dim)
    end
    # -> meth_reshape faster

    # display(c1)
    # display(c2)
    # display(c3)

    println()

    x = rand(50)
    y = rand(60)

    @time v1 = [x..., y...]
    @time v2 = vcat(x, y)
    # vcat much faster
end

println("\nhey")
test()
