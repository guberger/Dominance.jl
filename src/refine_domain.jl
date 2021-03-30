# Refine discretizetaion

function refine_domain(domain, nsub)
    println("refine_domain started")
    orig = domain.grid.orig
    h = domain.grid.h
    h2 = h./nsub
    orig2 = orig - h/2 + h2/2
    grid2 = Grid(orig2, h2)
    domain2 = Domain(grid2)
    sub_iter = hyper_range(nsub.*0, nsub .- 1)
    for pos in enum_pos(domain)
        for sub in sub_iter
            add_pos!(domain2, pos.*nsub .+ sub)
        end
    end
    ncells = get_ncells(domain2)
    println("refine_domain terminated: $(ncells) cells created")
    return domain2
end
