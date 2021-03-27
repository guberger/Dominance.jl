## Build symbolic model

# Assumes that automaton is "empty"
function symmodel_from_system!(symmodel::SymbolicModel{N}, sys) where N
    println("symmodel_from_system! started")
    domain = symmodel.domain
    r = domain.grid.h/2
    _H_ = SMatrix{N,N}(I).*r
    _ONE_ = ones(SVector{N})
    ntrans = 0
    Fe = sys.error_map(norm(r, Inf))
    Fr = typeof(r)(r .+ Fe)

    for pos in enum_pos(domain)
        source = get_state_by_pos(symmodel, pos)
        x = get_coord_by_pos(domain.grid, pos)
        Fx, DFx = sys.linsys_map(x, _H_)
        A = inv(DFx)
        b = abs.(A)*Fr .+ 1
        HP = CenteredPolyhedron(A, b)
        # TODO: can we improve abs.(DFx)*_ONE_?
        rad = abs.(DFx)*_ONE_ .+ Fe
        rectI = get_pos_lims_outer(domain.grid, HyperRectangle(Fx - rad, Fx + rad))
        # HyperRectangle(Fx - rad, Fx + rad) can be smaller than HP. Therefore,
        # in the plots, we may have cells not in symmodel while the 1st-order approx
        # cover them.
        fpos_iter = Iterators.product(_ranges(rectI)...)
        for fpos in fpos_iter
            dx = get_coord_by_pos(domain.grid, fpos) - Fx
            !(dx in HP) && continue
            if fpos in domain
                target = get_state_by_pos(symmodel, fpos)
                add_transition!(symmodel.autom, (source, source, target))
                ntrans += 1
            end
        end
    end

    println("symmodel_from_system! terminated with success: $(ntrans) transitions created")
end

## Essential automaton

function viable_controller!(contr, autom, viablelist)
    println("invariant_controller! started")
    nstates = autom.nstates
    viableset = Set(viablelist)
    unviableset = Set{Int}()
    pairstable = [false for i in 1:nstates, j in 1:nstates]

    for trans in enum_transitions(autom)
        if trans[1] in viableset && trans[3] in viableset
            pairstable[trans[1], trans[3]] = true
        end
    end

    while true
        empty!(unviableset)
        for state in viableset
            targets = view(pairstable, state, :)
            if !any(targets)
                pairstable[:, state] .= false
                push!(unviableset, state)
            end
            sources = view(pairstable, :, state)
            if !any(sources)
                pairstable[state, :] .= false
                push!(unviableset, state)
            end
        end
        if isempty(unviableset)
            break
        end
        setdiff!(viableset, unviableset)
    end

    for state in viableset
        add_state!(contr, state)
    end

    nviable = length(viableset)

    println("viable_controller! terminated with success: $(nviable) viable states added")
end
