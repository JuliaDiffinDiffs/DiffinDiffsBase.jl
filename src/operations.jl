function _relabel!(a1::Array, a2::AbstractArray)
    labels, invpool, pool = _label(a2)
    _relabel!(a1, labels, length(pool))
end

_relabel!(a1::Array, a2::PooledArray) = _relabel!(a1, a2.refs, length(a2.pool))

# Obtain unique labels for all row-wise pairs of values from a1 and a2
function _relabel!(a1::Array, a2::Array, mult::Integer)
    size(a1) == size(a2) || throw(DimensionMismatch(
        "cannot match array of size $(size(a1)) with array of size $(size(a2))"))
    a1 .+= mult .* a2
end

# A variant of SplitApplyCombine.groupfind using IdDict instead of Dictionaries.Dictionary
function _groupfind(container)
    T = keytype(container)
    inds = IdDict{eltype(container), Vector{T}}()
    @inbounds for i in keys(container)
        push!(get!(Vector{T}, inds, container[i]), i)
    end
    return inds
end

# Convert the keys from refs to cells
function _cellrows(cols::SubColumns, refrows::IdDict)
    ncol = size(cols, 2)
    cellrows = IdDict{Tuple, Vector{Int}}()
    for rows in values(refrows)
        cell = ntuple(n->cols[n][rows[1]], ncol)
        cellrows[cell] = rows
    end
    return cellrows
end

"""
    findcell(cellnames, data, esample=Colon())

Group the row indices of a subsample of `data` over `esample`
so that the row-wise combinations of values from columns indexed by `cellnames`
are the same within each group.

Note that unless `esample` covers all rows of `data`,
the row indices are those for the subsample selected based on `esample`
rather than those for the full `data`.

# Returns
- `IdDict{Tuple, Vector{Int}}`: a map from row-wise combinations of values to row indices of these combinations.
"""
function findcell(cellnames, data, esample=Colon())
    cols = SubColumns(data, cellnames, esample)
    isempty(cols) && throw(ArgumentError("empty data columns"))
    ncol = size(cols, 2)
    pooled = cols[1] isa PooledArray
    refs = pooled ? cols[1].refs : _label(cols[1])[1]
    if ncol > 1
        pooled && (refs = copy(refs))
        for n in 2:ncol
            _relabel!(refs, cols[n])
        end
    end
    cellrows = _cellrows(cols, _groupfind(refs))
    return cellrows
end
