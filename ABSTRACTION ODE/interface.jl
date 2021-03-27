function preprocess_idxlist(theAbs, idx_list; soft = false)
    # If idx_list is a singleton with a negative value, then it returns the list
    # of all indexes of the cells in theAbs.
    if isempty(idx_list)
        return Int[]
    end
    ncell = length(theAbs.cell_list)
    if length(idx_list) == 1 && idx_list[1] < 0
        return Vector{Int}(1:ncell)
    end
    finbounds = x -> (x > 0 && x <= ncell)
    if soft
        idx_list = filter(finbounds, idx_list)
    else
        @assert all(finbounds, idx_list)
        idx_list = copy(idx_list)
    end
    return unique(idx_list)
end

function check_nsublist(nsub_list, dim)
    @assert all(x -> length(x) == dim, nsub_list)
    @assert all(x -> all(x .> 0), nsub_list)
end

function preprocess_anylist(idx_list, any_list)
    # Checks whether any_list matches with idx_list and parses single lists to be
    # a uniform list of length equal to the one of idx_list.
    nidx = length(idx_list)
    if length(any_list) == 1
        any_list = fill(any_list[1], nidx)
    else
        @assert length(any_list) == nidx
    end
    return any_list
end
