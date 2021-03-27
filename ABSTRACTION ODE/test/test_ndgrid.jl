function test()
    x = (1,2,3,6,6,6,6,7,6)
    @time for i = 1:100
        [x...]
    end
    @time for i = 1:100
        collect(x)
    end
    # -> Comparable

    println()

    a = [1:0.02:i for i = 1:4]
    it = Iterators.product(a..., [1, 3])

    @time vec(map(sum, it))
    @time vec([sum(x) for x in it])
    @time [sum(x) for x in Iterators.take(it, length(it))]
    # -> [f(x) for x in it] slightly faster

    println()

    np = 500
    @time [1 for i = 1:np]
    @time convert(Vector{Int}, ones(np))
    @time Int.(ones(np))
    @time fill(1, np)
    # -> first and fourth faster and identical
end

sleep(0.1)
println("\nhey")
sleep(0.1) # used for good printing
test()
